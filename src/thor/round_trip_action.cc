#include "midgard/encoded.h"
#include "midgard/logging.h"
#include "midgard/polyline2.h"
#include "midgard/tiles.h"
#include "thor/worker.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace valhalla;
using namespace valhalla::baldr;
using namespace valhalla::midgard;

namespace {
constexpr size_t kMinimumCandidateEvaluations = 2;
// The public planner considers an "around" target accurate within 10%. Do not
// spend many additional full graph searches improving an already exact result.
constexpr double kEarlyAcceptanceTargetError = .10;
// Cost multiplier for roads an earlier leg used, and the corridor around them that also counts
// as used so that a parallel carriageway or service road isn't taken back instead
constexpr double kReuseFactor = 8.;
constexpr float kOverlapRadiusMeters = 50.f;
constexpr float kCorridorSampleSpacingMeters = 25.f;
// generalizing the leg shape keeps the corridor search cheap, well below the corridor radius
constexpr double kCorridorGeneralizeMeters = 10.;

// Cost factors for the edges a later leg should avoid: the leg's own edges plus every edge
// within the overlap corridor, in both directions, and the shortcuts those are part of.
// Resolved once when the leg is routed, then appended to the costing of every later leg.
using ReusedEdges = std::vector<CostFactorEdge>;

struct CandidateResult {
  TripRoute route;
  double target_error = std::numeric_limits<double>::max();
  double overlap = 1.;
  double score = std::numeric_limits<double>::max();
  bool target_mode_satisfied = false;
  bool anchor_uturn = false;
};

// Buckets the line's segments into cells at least `radius` wide, so a point only has to be
// compared against the segments of its own cell instead of against the whole line
class SegmentGrid {
public:
  SegmentGrid(const std::vector<PointLL>& line, float radius) : line_(line) {
    double max_lat = 0.;
    for (const auto& p : line)
      max_lat = std::max(max_lat, std::abs(p.lat()));
    lat_size_ = radius / kMetersPerDegreeLat;
    lng_size_ = lat_size_ / std::cos(max_lat * kRadPerDegD);

    // split segments into cell sized pieces, then register each piece with every cell its
    // bounding box, padded by the radius, overlaps
    for (size_t i = 0; i + 1 < line.size(); ++i) {
      const auto& u = line[i];
      const auto& v = line[i + 1];
      const auto pieces = std::max(1., std::ceil(std::max(std::abs(v.lat() - u.lat()) / lat_size_,
                                                          std::abs(v.lng() - u.lng()) / lng_size_)));
      for (double piece = 0; piece < pieces; ++piece) {
        const PointLL a(u.lng() + (v.lng() - u.lng()) * piece / pieces,
                        u.lat() + (v.lat() - u.lat()) * piece / pieces);
        const PointLL b(u.lng() + (v.lng() - u.lng()) * (piece + 1) / pieces,
                        u.lat() + (v.lat() - u.lat()) * (piece + 1) / pieces);
        for (auto x = cell(std::min(a.lng(), b.lng()) - lng_size_, lng_size_);
             x <= cell(std::max(a.lng(), b.lng()) + lng_size_, lng_size_); ++x) {
          for (auto y = cell(std::min(a.lat(), b.lat()) - lat_size_, lat_size_);
               y <= cell(std::max(a.lat(), b.lat()) + lat_size_, lat_size_); ++y) {
            auto& segments = cells_[key(x, y)];
            if (segments.empty() || segments.back() != i)
              segments.push_back(i);
          }
        }
      }
    }
  }

  bool is_within(const PointLL& point, float radius) const {
    const auto found = cells_.find(key(cell(point.lng(), lng_size_), cell(point.lat(), lat_size_)));
    if (found == cells_.cend())
      return false;
    for (const auto segment : found->second) {
      // a tiny forward cutoff limits ClosestPoint to the one segment starting at the pivot
      if (std::get<1>(point.ClosestPoint(line_, segment, 1e-9)) <= radius)
        return true;
    }
    return false;
  }

private:
  static int64_t cell(double coordinate, double size) {
    return static_cast<int64_t>(std::floor(coordinate / size));
  }
  static uint64_t key(int64_t x, int64_t y) {
    return (static_cast<uint64_t>(x) << 32) ^ static_cast<uint32_t>(y);
  }

  const std::vector<PointLL>& line_;
  double lat_size_;
  double lng_size_;
  std::unordered_map<uint64_t, std::vector<size_t>> cells_;
};

// Test a candidate edge at regular intervals rather than only at its shape
// vertices. This catches a long parallel edge whose endpoints lie outside the
// corridor, while keeping the distance calculation local and inexpensive.
bool is_near_line(const std::vector<PointLL>& candidate, const SegmentGrid& line, float radius) {
  if (candidate.size() == 1)
    return line.is_within(candidate.front(), radius);

  for (size_t i = 1; i < candidate.size(); ++i) {
    const auto& from = candidate[i - 1];
    const auto& to = candidate[i];
    const auto samples =
        std::max(1, static_cast<int>(std::ceil(from.Distance(to) / kCorridorSampleSpacingMeters)));
    for (int sample = 0; sample <= samples; ++sample) {
      const auto ratio = static_cast<double>(sample) / samples;
      PointLL point(from.lng() + (to.lng() - from.lng()) * ratio,
                    from.lat() + (to.lat() - from.lat()) * ratio);
      if (line.is_within(point, radius))
        return true;
    }
  }
  return false;
}

// Find graph edges in the line's tile bins and their immediate neighbours,
// then retain only edges whose geometry enters the corridor. The radius
// has to stay well below the size of a bin for the neighbours to suffice.
std::unordered_set<GraphId>
edges_near_line(GraphReader& reader, const std::vector<PointLL>& line_shape, float radius) {
  if (line_shape.size() < 2 || radius <= 0)
    return {};

  const auto& level = TileHierarchy::levels().back();
  const auto intersected_bins = level.tiles.Intersect(line_shape);
  std::unordered_set<uint64_t> visited_bins;
  std::unordered_set<GraphId> candidates;

  auto collect_bin = [&](uint32_t tile_id, unsigned short bin_id) {
    const auto key = static_cast<uint64_t>(tile_id) | (static_cast<uint64_t>(bin_id) << 32);
    if (!visited_bins.insert(key).second)
      return;
    auto tile = reader.GetGraphTile(GraphId(tile_id, level.level, 0));
    if (!tile)
      return;
    for (const auto& edge_id : tile->GetBin(bin_id))
      candidates.insert(edge_id);
  };

  for (const auto& [tile_id, bins] : intersected_bins) {
    for (const auto bin_id : bins) {
      collect_bin(tile_id, bin_id);
      for (const auto& neighbour : level.tiles.GetNeighboringBins<false>(tile_id, bin_id))
        collect_bin(neighbour.first, neighbour.second);
    }
  }

  const SegmentGrid grid(line_shape, radius);
  std::unordered_set<GraphId> nearby;
  for (const auto& edge_id : candidates) {
    auto tile = reader.GetGraphTile(edge_id);
    if (!tile)
      continue;
    const auto shape = tile->edgeinfo(tile->directededge(edge_id)).shape();
    if (is_near_line(shape, grid, radius))
      nearby.insert(edge_id);
  }
  return nearby;
}

uint64_t canonical_edge(GraphReader& reader, uint64_t value) {
  const GraphId edge(value);
  const auto opposing = reader.GetOpposingEdgeId(edge);
  return opposing.is_valid() ? std::min(value, opposing.value) : value;
}

std::vector<uint64_t> edge_ids(const TripLeg& leg) {
  std::vector<uint64_t> ids;
  ids.reserve(leg.node_size());
  for (const auto& node : leg.node()) {
    if (node.has_edge())
      ids.push_back(node.edge().id());
  }
  return ids;
}

CostFactorEdge reuse_factor(const GraphId& edge, double start = 0., double end = 1.) {
  CostFactorEdge factor;
  factor.set_id(edge);
  factor.set_factor(kReuseFactor);
  factor.set_start(start);
  factor.set_end(end);
  return factor;
}

ReusedEdges reused_edges(GraphReader& reader, const TripLeg& leg) {
  std::unordered_set<GraphId> edges;
  auto add_both_directions = [&](const GraphId& edge) {
    edges.insert(edge);
    const auto opposing = reader.GetOpposingEdgeId(edge);
    if (opposing.is_valid())
      edges.insert(opposing);
  };

  for (const auto id : edge_ids(leg))
    add_both_directions(GraphId(id));

  auto shape = decode<std::vector<PointLL>>(leg.shape());
  if (shape.size() > 2)
    Polyline2<PointLL>::Generalize(shape, kCorridorGeneralizeMeters);
  for (const auto& edge : edges_near_line(reader, shape, kOverlapRadiusMeters))
    add_both_directions(edge);

  ReusedEdges factors;
  factors.reserve(edges.size());
  std::unordered_set<GraphId> shortcuts;
  const auto local_level = TileHierarchy::levels().back().level;
  for (const auto& edge : edges) {
    factors.push_back(reuse_factor(edge));

    // the search may take a shortcut instead of its base edges, so it has to be just as
    // expensive over the stretch those cover
    const auto* directededge = reader.directededge(edge);
    if (edge.level() == local_level || !directededge || !directededge->superseded())
      continue;
    const auto shortcut = reader.GetShortcut(edge);
    if (!shortcut.is_valid() || !shortcuts.insert(shortcut).second)
      continue;
    const double shortcut_length = reader.directededge(shortcut)->length();
    double offset = 0.;
    for (const auto& constituent : reader.RecoverShortcut(shortcut)) {
      const auto length = reader.directededge(constituent)->length();
      if (edges.count(constituent))
        factors.push_back(
            reuse_factor(shortcut, offset / shortcut_length, (offset + length) / shortcut_length));
      offset += length;
    }
  }
  return factors;
}

void add_reuse_factors(Options& options, const std::vector<ReusedEdges>& used) {
  auto& costing_options =
      *options.mutable_costings()->find(options.costing_type())->second.mutable_options();
  for (const auto& prior : used) {
    for (const auto& factor : prior)
      costing_options.add_cost_factor_edges()->CopyFrom(factor);
  }
}

double overlap_fraction(GraphReader& reader, const TripRoute& route) {
  std::unordered_set<uint64_t> seen;
  double repeated_km = 0.;
  double total_km = 0.;
  for (const auto& leg : route.legs()) {
    for (const auto& node : leg.node()) {
      if (!node.has_edge())
        continue;
      const auto length = static_cast<double>(node.edge().length_km());
      total_km += length;
      if (!seen.insert(canonical_edge(reader, node.edge().id())).second)
        repeated_km += length;
    }
  }
  return total_km > 0. ? repeated_km / total_km : 1.;
}

double route_time_seconds(const TripRoute& route) {
  double seconds = 0.;
  for (const auto& leg : route.legs()) {
    if (!leg.node().empty())
      seconds += leg.node().rbegin()->cost().elapsed_cost().seconds();
  }
  return seconds;
}

double route_length(const TripRoute& route) {
  double length = 0.;
  for (const auto& leg : route.legs()) {
    for (const auto& node : leg.node()) {
      if (node.has_edge())
        length += node.edge().length_km();
    }
  }
  return length;
}
} // namespace

namespace valhalla {
namespace thor {

// Route each leg independently and rebuild costing after every completed leg.
// This is the key distinction from an ordinary multi-location route: edges
// used by leg N are already expensive while leg N+1 is being searched.
void thor_worker_t::round_trip(Api& request) {
  const auto search_started = std::chrono::steady_clock::now();
  const auto& round_trip = request.options().round_trip_options();
  if (round_trip.anchors().empty())
    throw valhalla_exception_t{171, "no viable round-trip anchors"};

  // every leg is routed as a plain two location route without the round trip specifics
  Options leg_options = request.options();
  leg_options.clear_round_trip_options();
  leg_options.clear_locations();
  leg_options.set_alternates(0);

  // Route a single leg. When more legs follow, also resolve the edges they
  // should avoid; the closing leg back to the origin skips that work.
  auto append_leg = [&](const Location& from, const Location& to, std::vector<ReusedEdges>& used,
                        TripRoute& output, bool more_legs_follow) {
    Api leg_request;
    auto& options = *leg_request.mutable_options();
    options.CopyFrom(leg_options);
    options.add_locations()->CopyFrom(from);
    options.add_locations()->CopyFrom(to);
    options.mutable_locations(0)->set_type(Location::kBreak);
    options.mutable_locations(1)->set_type(Location::kBreak);
    add_reuse_factors(options, used);

    try {
      route(leg_request);
    } catch (const valhalla_exception_t&) { return false; }
    if (leg_request.trip().routes().empty() || leg_request.trip().routes(0).legs().empty())
      return false;

    for (const auto& leg : leg_request.trip().routes(0).legs()) {
      if (edge_ids(leg).empty())
        return false;
      if (more_legs_follow)
        used.push_back(reused_edges(*reader, leg));
      output.add_legs()->CopyFrom(leg);
    }
    return true;
  };

  // All user-supplied vias precede the generated anchor, so their legs are
  // identical for every candidate. Compute this prefix once instead of once
  // per anchor (the dominant cost for POI-shaped loops).
  const auto& options = request.options();
  const auto& origin = options.locations(0);
  const auto& last_user_location = options.locations(options.locations_size() - 1);
  TripRoute prefix_route;
  std::vector<ReusedEdges> prefix_used;
  for (int i = 0; i + 1 < options.locations_size(); ++i) {
    if (!append_leg(options.locations(i), options.locations(i + 1), prefix_used, prefix_route, true))
      throw valhalla_exception_t{442};
  }

  const auto target = round_trip.target_distance() > 0
                          ? static_cast<double>(round_trip.target_distance()) *
                                (options.units() == Options::miles ? 1.609344 : 1.)
                          : static_cast<double>(round_trip.target_time()) * 60.;
  CandidateResult best;
  bool found = false;
  size_t evaluated = 0;

  auto consider = [&](CandidateResult candidate) {
    ++evaluated;
    const auto actual = round_trip.target_distance() > 0 ? route_length(candidate.route)
                                                         : route_time_seconds(candidate.route);
    candidate.target_error = std::abs(actual - target) / target;
    candidate.overlap = overlap_fraction(*reader, candidate.route);
    const auto over_target = actual > target ? (actual - target) / target : 0.;
    const auto under_target = actual < target ? (target - actual) / target : 0.;
    candidate.target_mode_satisfied = (!round_trip.prefer_under() || actual <= target) &&
                                      (!round_trip.prefer_at_least() || actual >= target);
    candidate.score = candidate.target_error + 4. * candidate.overlap;
    if (candidate.overlap > round_trip.max_overlap())
      candidate.score += 10. * (candidate.overlap - round_trip.max_overlap());
    if (round_trip.prefer_under())
      candidate.score += 8. * over_target;
    if (round_trip.prefer_at_least())
      candidate.score += 8. * under_target;
    // A generated anchor is a shaping point, not a destination. An immediate
    // reversal there creates a conspicuous stem. Keep it as a last-resort
    // fallback, but let any continuous candidate outrank it.
    if (candidate.anchor_uturn)
      candidate.score += .25;
    LOG_DEBUG("round_trip candidate length=" + std::to_string(route_length(candidate.route)) +
              " overlap=" + std::to_string(candidate.overlap) +
              " target_error=" + std::to_string(candidate.target_error) + " anchor_uturn=" +
              std::to_string(candidate.anchor_uturn) + " score=" + std::to_string(candidate.score));

    if (!found || candidate.score < best.score) {
      best = std::move(candidate);
      found = true;
    }
  };

  // When the user's vias already shape a loop of the right size, going straight
  // back to the origin beats any detour through a generated anchor. It costs a
  // single leg and competes with the anchors on the same score.
  const auto has_vias = options.locations_size() > 1;
  if (has_vias) {
    CandidateResult direct;
    direct.route.CopyFrom(prefix_route);
    if (append_leg(last_user_location, origin, prefix_used, direct.route, false))
      consider(std::move(direct));
  }

  const auto minimum_evaluations =
      std::min(kMinimumCandidateEvaluations,
               static_cast<size_t>(round_trip.anchors_size()) + (has_vias ? 1 : 0));
  size_t examined = 0;
  for (const auto& anchor : round_trip.anchors()) {
    // Stop early once a diverse minimum set contains a route within the
    // planner's target tolerance and the requested overlap threshold.
    if (found && evaluated >= minimum_evaluations && best.target_mode_satisfied &&
        !best.anchor_uturn && best.target_error <= kEarlyAcceptanceTargetError &&
        best.overlap <= round_trip.max_overlap())
      break;
    // Bound the latency of long loops: once out of time, a usable loop now
    // beats a slightly better one later
    const std::chrono::duration<float> elapsed = std::chrono::steady_clock::now() - search_started;
    if (found && elapsed.count() > round_trip_time_budget_seconds) {
      LOG_WARN("round_trip time budget spent after " + std::to_string(examined) + " anchors");
      break;
    }

    ++examined;
    CandidateResult candidate;
    candidate.route.CopyFrom(prefix_route);
    auto used = prefix_used;
    if (!append_leg(last_user_location, anchor, used, candidate.route, true))
      continue;
    const auto arrival_edge = GraphId(edge_ids(*candidate.route.legs().rbegin()).back());
    const auto return_leg_index = candidate.route.legs_size();
    if (!append_leg(anchor, origin, used, candidate.route, false))
      continue;
    const auto departure_edges = edge_ids(candidate.route.legs(return_leg_index));
    const auto uturn_edge = reader->GetOpposingEdgeId(arrival_edge);
    candidate.anchor_uturn = !departure_edges.empty() && uturn_edge.is_valid() &&
                             departure_edges.front() == uturn_edge.value;
    consider(std::move(candidate));
  }

  if (!found)
    throw valhalla_exception_t{442};
  LOG_INFO("round_trip selected overlap=" + std::to_string(best.overlap) + " target_error=" +
           std::to_string(best.target_error) + " candidates=" + std::to_string(evaluated) + "/" +
           std::to_string(examined) + " prefix_legs=" + std::to_string(prefix_route.legs_size()));

  const auto leg_count = best.route.legs_size();
  for (int i = 0; i < leg_count; ++i) {
    best.route.mutable_legs(i)->set_leg_id(i);
    best.route.mutable_legs(i)->set_leg_count(leg_count);
  }
  request.clear_trip();
  request.mutable_trip()->add_routes()->CopyFrom(best.route);
}

} // namespace thor
} // namespace valhalla
