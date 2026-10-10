#include "midgard/encoded.h"
#include "midgard/logging.h"
#include "midgard/polyline2.h"
#include "thor/worker.h"

#include <algorithm>
#include <cmath>
#include <limits>
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

struct UsedLeg {
  std::vector<uint64_t> edge_ids;
  std::vector<PointLL> shape;
};

struct CandidateResult {
  TripRoute route;
  double target_error = std::numeric_limits<double>::max();
  double overlap = 1.;
  double score = std::numeric_limits<double>::max();
  bool target_mode_satisfied = false;
  bool anchor_uturn = false;
};

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

void add_reuse_factors(Options& options,
                       const std::vector<UsedLeg>& used,
                       const RoundTripOptions& round_trip) {
  for (const auto& prior : used) {
    auto* line = options.add_cost_factor_lines();
    line->set_cost_factor(round_trip.reuse_factor());
    line->set_bidirectional(true);
    line->set_radius(round_trip.overlap_radius());
    for (const auto id : prior.edge_ids)
      line->add_edge_ids(id);
    for (const auto& point : prior.shape) {
      auto* location = line->add_shape();
      location->mutable_ll()->set_lng(point.lng());
      location->mutable_ll()->set_lat(point.lat());
    }
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
  const auto round_trip = request.options().round_trip_options();
  if (round_trip.anchors().empty())
    throw valhalla_exception_t{171, "no viable round-trip anchors"};

  // Route a single leg and retain both its exact graph edges and geometry for
  // the reuse/corridor penalty applied to later legs.
  auto append_leg = [&](const Location& from,
                        const Location& to,
                        std::vector<UsedLeg>& used,
                        TripRoute& output) {
    Api leg_request;
    leg_request.mutable_options()->CopyFrom(request.options());
    leg_request.mutable_options()->clear_round_trip_options();
    leg_request.mutable_options()->set_alternates(0);
    leg_request.mutable_options()->mutable_locations()->Clear();
    leg_request.mutable_options()->add_locations()->CopyFrom(from);
    leg_request.mutable_options()->add_locations()->CopyFrom(to);
    leg_request.mutable_options()->mutable_locations(0)->set_type(Location::kBreak);
    leg_request.mutable_options()->mutable_locations(1)->set_type(Location::kBreak);
    add_reuse_factors(*leg_request.mutable_options(), used, round_trip);

    try {
      route(leg_request);
    } catch (const valhalla_exception_t&) { return false; }
    if (leg_request.trip().routes().empty() || leg_request.trip().routes(0).legs().empty())
      return false;

    for (const auto& leg : leg_request.trip().routes(0).legs()) {
      auto ids = edge_ids(leg);
      if (ids.empty())
        return false;
      auto shape = decode<std::vector<PointLL>>(leg.shape());
      if (round_trip.overlap_radius() > 0 && shape.size() > 2) {
        // Corridor matching compares nearby graph edges against this shape. A
        // 10 m generalization is safely below the 50 m avoidance radius and
        // prevents that comparison from becoming quadratic on long routes.
        Polyline2<PointLL>::Generalize(
            shape, std::max(2., std::min(10., round_trip.overlap_radius() * .2)));
      }
      used.push_back({std::move(ids), std::move(shape)});
      output.add_legs()->CopyFrom(leg);
    }
    return true;
  };

  // All user-supplied vias precede the generated anchor, so their legs are
  // identical for every candidate. Compute this prefix once instead of once
  // per anchor (the dominant cost for POI-shaped loops).
  TripRoute prefix_route;
  std::vector<UsedLeg> prefix_used;
  for (int i = 0; i + 1 < request.options().locations_size(); ++i) {
    if (!append_leg(request.options().locations(i), request.options().locations(i + 1),
                    prefix_used, prefix_route))
      throw valhalla_exception_t{442};
  }

  CandidateResult best;
  bool found = false;
  size_t evaluated = 0;
  size_t examined = 0;
  const auto minimum_evaluations =
      std::min(kMinimumCandidateEvaluations, static_cast<size_t>(round_trip.anchors_size()));
  for (const auto& anchor : round_trip.anchors()) {
    ++examined;
    CandidateResult candidate;
    candidate.route.CopyFrom(prefix_route);
    auto used = prefix_used;
    const auto& last_user_location =
        request.options().locations(request.options().locations_size() - 1);
    if (!append_leg(last_user_location, anchor, used, candidate.route))
      continue;
    const auto arrival_edge = GraphId(used.back().edge_ids.back());
    const auto return_leg_index = candidate.route.legs_size();
    if (!append_leg(anchor, request.options().locations(0), used, candidate.route))
      continue;
    const auto departure_edges = edge_ids(candidate.route.legs(return_leg_index));
    const auto uturn_edge = reader->GetOpposingEdgeId(arrival_edge);
    candidate.anchor_uturn = !departure_edges.empty() && uturn_edge.is_valid() &&
                             departure_edges.front() == uturn_edge.value;
    ++evaluated;

    const auto actual = round_trip.target_distance() > 0
                            ? route_length(candidate.route)
                            : route_time_seconds(candidate.route);
    const auto target = round_trip.target_distance() > 0
                            ? static_cast<double>(round_trip.target_distance()) *
                                  (request.options().units() == Options::miles ? 1.609344 : 1.)
                            : static_cast<double>(round_trip.target_time()) * 60.;
    candidate.target_error = std::abs(actual - target) / target;
    candidate.overlap = overlap_fraction(*reader, candidate.route);
    const auto over_target = actual > target ? (actual - target) / target : 0.;
    const auto under_target = actual < target ? (target - actual) / target : 0.;
    candidate.target_mode_satisfied =
        (!round_trip.prefer_under() || actual <= target) &&
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
              " target_error=" + std::to_string(candidate.target_error) +
              " anchor_uturn=" + std::to_string(candidate.anchor_uturn) +
              " score=" + std::to_string(candidate.score));

    if (!found || candidate.score < best.score) {
      best = std::move(candidate);
      found = true;
    }

    // Stop early once a diverse minimum set contains a route within the
    // planner's target tolerance and the requested overlap threshold.
    if (evaluated >= minimum_evaluations && best.target_mode_satisfied && !best.anchor_uturn &&
        best.target_error <= kEarlyAcceptanceTargetError &&
        best.overlap <= round_trip.max_overlap())
      break;
  }

  if (!found)
    throw valhalla_exception_t{442};
  LOG_INFO("round_trip selected overlap=" + std::to_string(best.overlap) +
           " target_error=" + std::to_string(best.target_error) +
           " candidates=" + std::to_string(evaluated) + "/" + std::to_string(examined) +
           " prefix_legs=" + std::to_string(prefix_route.legs_size()));

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
