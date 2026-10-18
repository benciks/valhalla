#include "loki/linear_cost_factors.h"
#include "exceptions.h"
#include "midgard/tiles.h"
#include "thor/route_matcher.h"

#include <cmath>
#include <unordered_map>
#include <unordered_set>

using namespace valhalla;
using namespace valhalla::baldr;
using namespace valhalla::midgard;
using namespace valhalla::thor;

namespace {

// first and last edges of a cost factor line covering less than this are dropped
constexpr double kMinCostFactorEdgeLength = 1.0; // meters
constexpr float kMaxLinearFeatureRadiusMeters = 200.0f;
constexpr float kCorridorSampleSpacingMeters = 25.0f;

/**
 * Adds a shortcut to the cost factor edges given one
 * of its constituents
 */
void add_shortcut(GraphReader& reader,
                  GraphId shortcut,
                  valhalla::Costing_Options* options,
                  valhalla::CostFactorEdge* cost_factor) {

  // allow applies to the whole edge
  if (cost_factor->allow_types() && cost_factor->factor() != 1.0) {
    auto* e = options->add_cost_factor_edges();
    e->set_id(shortcut.value);
    e->set_allow_types(cost_factor->allow_types());
    return;
  }

  GraphId edge = static_cast<GraphId>(cost_factor->id());
  graph_tile_ptr tile = reader.GetGraphTile(shortcut);
  // it's part of a shortcut
  auto constituents = reader.RecoverShortcut(shortcut);
  auto* shortcut_edge = tile->directededge(shortcut);

  tile = reader.GetGraphTile(edge);
  auto* current_edge = tile->directededge(edge);

  // walk the base edges until we find ours
  uint64_t accumulated_length = 0;
  for (const auto& constituent : constituents) {
    if (edge == constituent)
      break;

    tile = reader.GetGraphTile(constituent, tile);
    if (!tile)
      break;

    auto* de = tile->directededge(constituent);
    accumulated_length += de->length();
  }
  auto* e = options->add_cost_factor_edges();
  e->set_id(shortcut);
  e->set_factor(cost_factor->factor());
  e->set_start(static_cast<double>(accumulated_length + (static_cast<double>(current_edge->length()) *
                                                         cost_factor->start())) /
               static_cast<double>(shortcut_edge->length()));
  e->set_end(static_cast<double>(accumulated_length +
                                 (static_cast<double>(current_edge->length()) * cost_factor->end())) /
             static_cast<double>(shortcut_edge->length()));
}

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

} // namespace

namespace valhalla {
namespace loki {

// Find graph edges in the route's tile bins and their immediate neighbours,
// then retain only edges whose geometry enters the requested corridor. The
// supported radius is deliberately smaller than a hierarchy-level bin.
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

void add_cost_factor_edges(const sif::mode_costing_t& mode_costing,
                           const sif::TravelMode& mode,
                           GraphReader& reader,
                           valhalla::Options& options,
                           double min_allowed_factor,
                           uint64_t max_allowed_edges) {
  Costing_Options* costing_options =
      options.mutable_costings()->find(options.costing_type())->second.mutable_options();

  // keep track of how many edges we're adding
  uint64_t edge_count = 0;

  for (auto& line : *options.mutable_cost_factor_lines()) {
    const auto first_line_edge = costing_options->cost_factor_edges_size();
    if (!line.edge_ids().empty()) {
      for (const auto edge_id_value : line.edge_ids()) {
        if (edge_count >= max_allowed_edges)
          throw valhalla_exception_t{234};
        edge_count++;
        const GraphId edge_id(edge_id_value);
        auto* edge = costing_options->add_cost_factor_edges();
        edge->set_id(edge_id);
        edge->set_factor(std::max(line.cost_factor(), min_allowed_factor));
        edge->set_start(0.);
        edge->set_end(1.);
        edge->set_allow_types(line.allow_types());
        const auto shortcut = reader.GetShortcut(edge_id);
        if (shortcut.is_valid())
          add_shortcut(reader, shortcut, costing_options, edge);
      }
    } else {
      std::vector<std::vector<PathInfo>> legs;
      if (!RouteMatcher::FormPath(mode_costing, mode, reader, line, false, /* use_shortcuts=*/true,
                                  legs)) {
        throw valhalla_exception_t{233};
      }
      for (const auto& leg : legs) {
        for (size_t i = 0; i < leg.size(); ++i) {
          if (edge_count > max_allowed_edges)
            throw valhalla_exception_t{234};
          auto& path_info = leg[i];
          bool is_first = i == 0;
          bool is_last = i == leg.size() - 1;
          if (is_first && is_last) { // trivial path
            edge_count++;
            auto* e = costing_options->add_cost_factor_edges();
            e->set_id(path_info.edgeid);
            e->set_factor(line.cost_factor());
            e->set_allow_types(line.allow_types());
            for (const auto& edge : line.locations(0).correlation().edges()) {
              if (path_info.edgeid == edge.graph_id()) {
                e->set_start(edge.percent_along());
                break;
              }
            }
            for (const auto& edge : line.locations(1).correlation().edges()) {
              if (path_info.edgeid == edge.graph_id()) {
                e->set_end(edge.percent_along());
                break;
              }
            }
            auto shortcut = reader.GetShortcut(path_info.edgeid);
            if (shortcut.is_valid()) {
              add_shortcut(reader, shortcut, costing_options, e);
            }
          } else if (is_first || is_last) { // beginning or end edge
            for (const auto& edge :
                 line.locations(static_cast<size_t>(is_last)).correlation().edges()) {
              if (path_info.edgeid == edge.graph_id()) {
                double start = is_first ? edge.percent_along() : 0.;
                double end = is_last ? edge.percent_along() : 1.;
                // an endpoint just off a node also correlates to the adjacent edges, skip those small
                // parts
                const auto* de = reader.directededge(path_info.edgeid);
                if (de && de->length() * (end - start) < kMinCostFactorEdgeLength) {
                  break;
                }
                edge_count++;
                auto* e = costing_options->add_cost_factor_edges();
                e->set_id(path_info.edgeid);
                // apply the minimum allowed value specified in the config
                e->set_factor(std::max(line.cost_factor(), min_allowed_factor));
                e->set_allow_types(line.allow_types());
                e->set_start(start);
                e->set_end(end);
                auto shortcut = reader.GetShortcut(path_info.edgeid);
                if (shortcut.is_valid()) {
                  add_shortcut(reader, shortcut, costing_options, e);
                }
                break;
              }
            }
          } else { // intermediate edges
            edge_count++;
            auto* e = costing_options->add_cost_factor_edges();
            e->set_id(path_info.edgeid);
            e->set_factor(std::max(line.cost_factor(), min_allowed_factor));
            e->set_allow_types(line.allow_types());
            e->set_start(0.);
            e->set_end(1.);

            // if it's a shortcut, also add all of its constituent edges
            if (path_info.is_shortcut) {
              auto constituents = reader.RecoverShortcut(path_info.edgeid);
              for (const auto& constituent : constituents) {
                edge_count++;
                auto* e = costing_options->add_cost_factor_edges();
                e->set_id(constituent);
                e->set_factor(std::max(line.cost_factor(), min_allowed_factor));
                e->set_allow_types(line.allow_types());
                e->set_start(0);
                e->set_end(1);
              }
            } else {
              // if it's not a shortcut, it may be part of one
              // TODO: this is an expensive operation, since we need to expand the graph
              // a little, can't we persist this information somehow?
              auto shortcut = reader.GetShortcut(path_info.edgeid);
              if (shortcut.is_valid()) {
                add_shortcut(reader, shortcut, costing_options, e);
              }
            }
          }
        }
      }
    }

    // Graph edges are directed. For use cases such as round-trip route
    // generation, reusing the same road in the opposite direction must carry
    // the same factor. Add the opposing edge after the source line has been
    // edge-walked so one-way portions do not make reverse map matching fail.
    if (line.bidirectional()) {
      const auto last_line_edge = costing_options->cost_factor_edges_size();
      for (int i = first_line_edge; i < last_line_edge; ++i) {
        const auto source = costing_options->cost_factor_edges(i);
        const auto opposing = reader.GetOpposingEdgeId(static_cast<GraphId>(source.id()));
        if (!opposing.is_valid())
          continue;
        if (edge_count >= max_allowed_edges)
          throw valhalla_exception_t{234};

        edge_count++;
        auto* edge = costing_options->add_cost_factor_edges();
        edge->set_id(opposing);
        edge->set_factor(source.factor());
        edge->set_start(1. - source.end());
        edge->set_end(1. - source.start());
        edge->set_allow_types(source.allow_types());
      }
    }

    const auto radius = std::max(0.0f, std::min(line.radius(), kMaxLinearFeatureRadiusMeters));
    if (radius > 0) {
      std::unordered_set<GraphId> factored_edges;
      for (const auto& edge : costing_options->cost_factor_edges())
        factored_edges.insert(static_cast<GraphId>(edge.id()));

      auto add_full_edge_factor = [&](const GraphId& edge_id) {
        if (!edge_id.is_valid() || !factored_edges.insert(edge_id).second)
          return;
        if (edge_count >= max_allowed_edges)
          throw valhalla_exception_t{234};

        edge_count++;
        auto* edge = costing_options->add_cost_factor_edges();
        edge->set_id(edge_id);
        edge->set_factor(std::max(line.cost_factor(), min_allowed_factor));
        edge->set_start(0.);
        edge->set_end(1.);
        edge->set_allow_types(line.allow_types());

        const auto shortcut = reader.GetShortcut(edge_id);
        if (shortcut.is_valid())
          add_shortcut(reader, shortcut, costing_options, edge);
      };

      std::vector<PointLL> line_shape;
      line_shape.reserve(line.shape_size());
      for (const auto& location : line.shape())
        line_shape.emplace_back(location.ll().lng(), location.ll().lat());
      for (const auto& edge_id : edges_near_line(reader, line_shape, radius)) {
        add_full_edge_factor(edge_id);
        if (line.bidirectional())
          add_full_edge_factor(reader.GetOpposingEdgeId(edge_id));
      }
    }
  }
}

} // namespace loki
} // namespace valhalla
