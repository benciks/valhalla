#include "baldr/tilehierarchy.h"
#include "loki/search.h"
#include "loki/worker.h"

#include <algorithm>
#include <cmath>

using namespace valhalla;
using namespace valhalla::baldr;

namespace {
constexpr double kEarthRadiusMeters = 6371000.;
constexpr double kGoldenAngleDegrees = 137.50776405003785;
constexpr double kEstimatedRoadStretch[] = {1.18, 1.28, 1.10};

midgard::PointLL project(const midgard::PointLL& origin, double bearing_degrees, double meters) {
  const auto bearing = bearing_degrees * M_PI / 180.;
  const auto lat = origin.lat() * M_PI / 180.;
  const auto angular = meters / kEarthRadiusMeters;
  const auto projected_lat = std::asin(std::sin(lat) * std::cos(angular) +
                                       std::cos(lat) * std::sin(angular) * std::cos(bearing));
  const auto projected_lon = origin.lng() * M_PI / 180. +
                             std::atan2(std::sin(bearing) * std::sin(angular) * std::cos(lat),
                                        std::cos(angular) - std::sin(lat) * std::sin(projected_lat));
  return {projected_lon * 180. / M_PI, projected_lat * 180. / M_PI};
}

size_t append_round_trip_candidates(Options& options) {
  if (!options.has_round_trip_options())
    return 0;
  auto& round_trip = *options.mutable_round_trip_options();
  round_trip.clear_anchors();
  if (options.locations().empty())
    throw valhalla_exception_t{120};
  if ((round_trip.target_distance() <= 0) == (round_trip.target_time() <= 0))
    throw valhalla_exception_t{113, "exactly one of target_distance or target_time"};
  // A motorcycle minute is estimated as 0.92 km solely to place initial
  // anchors. Thor scores actual route time and distance afterwards.
  double target_meters = round_trip.target_distance() > 0 ? round_trip.target_distance() * 1000.
                                                          : round_trip.target_time() * 920.;
  if (options.units() == Options::miles && round_trip.target_distance() > 0)
    target_meters *= 1.609344;
  const auto origin =
      midgard::PointLL(options.locations(0).ll().lng(), options.locations(0).ll().lat());
  const auto last = midgard::PointLL(options.locations().rbegin()->ll().lng(),
                                     options.locations().rbegin()->ll().lat());
  double prefix_meters = 0.;
  for (int i = 1; i < options.locations_size(); ++i) {
    const auto from =
        midgard::PointLL(options.locations(i - 1).ll().lng(), options.locations(i - 1).ll().lat());
    const auto to =
        midgard::PointLL(options.locations(i).ll().lng(), options.locations(i).ll().lat());
    prefix_meters += from.Distance(to);
  }
  const auto count = std::clamp(round_trip.candidates(), 1u, 24u);
  const auto rotation = static_cast<double>(round_trip.seed() % 360);
  const auto outward_bearing = options.locations_size() > 1 ? origin.Heading(last) : rotation;
  const auto side = (round_trip.seed() & 1) == 0 ? 1. : -1.;
  for (uint32_t i = 0; i < count; ++i) {
    // Solve for an anchor whose straight-line itinerary consumes the target
    // after accounting for fixed vias. This makes the first candidates useful
    // instead of blindly trying many radii. Alternating road-stretch estimates
    // brackets regions whose road network is more or less direct.
    const auto bearing =
        options.locations_size() > 1
            ? std::fmod(outward_bearing + side * (105. + kGoldenAngleDegrees * i) + 720., 360.)
            : std::fmod(rotation + kGoldenAngleDegrees * i, 360.);
    const auto target_air_meters =
        target_meters / kEstimatedRoadStretch[i % std::size(kEstimatedRoadStretch)];
    const auto minimum_anchor_legs = last.Distance(origin) + std::max(1000., target_meters * .04);
    const auto desired_anchor_legs = std::max(minimum_anchor_legs, target_air_meters - prefix_meters);
    double low = 500.;
    double high = std::clamp(target_meters, 5000., 2500000.);
    for (size_t iteration = 0; iteration < 28; ++iteration) {
      const auto radius = (low + high) * .5;
      const auto candidate = project(origin, bearing, radius);
      if (radius + last.Distance(candidate) < desired_anchor_legs)
        low = radius;
      else
        high = radius;
    }
    const auto radius = (low + high) * .5;
    const auto ll = project(origin, bearing, radius);
    auto* candidate = options.add_locations();
    candidate->mutable_ll()->set_lng(ll.lng());
    candidate->mutable_ll()->set_lat(ll.lat());
    candidate->set_type(Location::kBreak);
    candidate->set_minimum_reachability(50);
  }
  // Keep the real origin as the final correlation so generated candidates do
  // not receive route-end reachability relaxation.
  options.add_locations()->CopyFrom(options.locations(0));
  return count;
}

void check_locations(const size_t location_count, const size_t max_locations) {
  // check that location size does not exceed max.
  if (location_count > max_locations) {
    throw valhalla_exception_t{150, std::to_string(max_locations)};
  };
}

void check_distance(const google::protobuf::RepeatedPtrField<valhalla::Location>& locations,
                    float max_distance,
                    bool all_pairs) {
  // test if total distance along a polyline formed by connecting locations exceeds the maximum
  // or if all_pairs is specified test all pairs of locations to see if any are over the threshold
  float total_path_distance = 0.0f;
  for (int i = 0; i < locations.size(); ++i) {
    for (int j = i + 1; j < locations.size(); ++j) {
      auto dist = to_ll(locations.Get(i)).Distance(to_ll(locations.Get(j)));
      total_path_distance += i + 1 == j ? dist : 0;
      if ((!all_pairs && total_path_distance > max_distance) || (all_pairs && dist > max_distance))
        throw valhalla_exception_t{154,
                                   std::to_string(static_cast<size_t>(max_distance)) + " meters"};
      if (!all_pairs)
        break;
    }
  }
}

} // namespace

namespace valhalla {
namespace loki {

void loki_worker_t::init_route(Api& request) {
  parse_locations(request.mutable_options()->mutable_locations(), request);
  // need to check location size here instead of in parse_locations because of locate action needing
  // a different size
  const auto minimum_locations = request.options().has_round_trip_options() ? 1 : 2;
  if (request.options().locations_size() < minimum_locations) {
    throw valhalla_exception_t{120};
  };
  parse_costing(request);
}

void loki_worker_t::route(Api& request) {
  // time this whole method and save that statistic
  auto _ = measure_scope_time(request);

  init_route(request);
  auto& options = *request.mutable_options();
  if (options.has_round_trip_options()) {
    const auto& round_trip = options.round_trip_options();
    if (static_cast<size_t>(options.locations_size()) > max_round_trip_locations)
      throw valhalla_exception_t{150, " for round_trip: " + std::to_string(max_round_trip_locations)};
    // A motorcycle minute is estimated as 0.92 km, as for anchor placement
    auto target_km = round_trip.target_distance() > 0 ? round_trip.target_distance()
                                                      : round_trip.target_time() * .92f;
    if (options.units() == Options::miles && round_trip.target_distance() > 0)
      target_km *= 1.609344f;
    if (target_km > max_round_trip_distance_km)
      throw valhalla_exception_t{154,
                                 " for round_trip: " +
                                     std::to_string(static_cast<int>(max_round_trip_distance_km)) +
                                     " km"};
  }
  const auto original_location_count = options.locations_size();
  const auto generated_candidate_count = append_round_trip_candidates(options);
  for (size_t i = original_location_count; i < original_location_count + generated_candidate_count;
       ++i)
    parse_location(*options.mutable_locations(i));
  const auto& costing_name = Costing_Enum_Name(options.costing_type());
  if (request.options().action() == Options::centroid) {
    check_locations(options.locations_size(), max_locations.find("centroid")->second);
    check_distance(options.locations(), max_distance.find("centroid")->second, true);
  } else {
    // Generated anchors are engine search candidates rather than mandatory
    // user stops; do not count their fan-out against route request limits.
    check_locations(generated_candidate_count ? original_location_count : options.locations_size(),
                    max_locations.find(costing_name)->second);
    if (!generated_candidate_count)
      check_distance(options.locations(), max_distance.find(costing_name)->second, false);
  }

  // check distance for hierarchy pruning
  check_hierarchy_distance(request);

  auto connectivity_level = TileHierarchy::levels().back();
  uint32_t connectivity_radius = 0;
  // Validate walking distances (make sure they are in the accepted range)
  if (costing_name == "multimodal" || costing_name == "transit" ||
      costing_name == "auto_pedestrian") {
    auto& ped_opts = *options.mutable_costings()->find(Costing::pedestrian)->second.mutable_options();

    // "transit_start_end_max_distance" is deprecated
    // but we will still allow it for some time
    if (ped_opts.has_transit_transfer_max_distance_case()) {
      ped_opts.set_multimodal_start_end_max_distance(ped_opts.transit_start_end_max_distance());
    }

    // we have renamed this parameter
    if (!ped_opts.has_multimodal_start_end_max_distance_case())
      ped_opts.set_multimodal_start_end_max_distance(min_multimodal_walking_dist);
    auto multimodal_start_end_max_distance = ped_opts.multimodal_start_end_max_distance();

    if (!ped_opts.has_transit_transfer_max_distance_case())
      ped_opts.set_transit_transfer_max_distance(min_multimodal_walking_dist);
    auto transit_transfer_max_distance = ped_opts.transit_transfer_max_distance();

    if (multimodal_start_end_max_distance < min_multimodal_walking_dist ||
        multimodal_start_end_max_distance > max_multimodal_walking_dist) {
      throw valhalla_exception_t{155, " Min: " + std::to_string(min_multimodal_walking_dist) +
                                          " Max: " + std::to_string(max_multimodal_walking_dist) +
                                          " (Meters)"};
    }
    if (transit_transfer_max_distance < min_multimodal_walking_dist ||
        transit_transfer_max_distance > max_multimodal_walking_dist) {
      throw valhalla_exception_t{156, " Min: " + std::to_string(min_multimodal_walking_dist) +
                                          " Max: " + std::to_string(max_multimodal_walking_dist) +
                                          " (Meters)"};
    }
    if (costing_name != "auto_pedestrian") {
      connectivity_level = TileHierarchy::GetTransitLevel();
    }
    connectivity_radius = ped_opts.transit_start_end_max_distance();
  }

  // correlate the various locations to the underlying graph
  std::unordered_map<size_t, size_t> color_counts;
  try {
    auto* locations = options.mutable_locations();
    locations->begin()->set_minimum_inbound_reachability(0);
    locations->rbegin()->set_minimum_outbound_reachability(0);
    auto locations_size = locations->size();

    // in case of auto_pedestrian costing, we 1) only allow two locations
    // and 2) need two different costings for the start and end location.
    // Search::search does not allow for multiple costings per location so instead
    // we bite the bullet and call search twice, merging the results.
    // TODO(chris): right now we hash location based purely on ll, but what if the user wants to start
    // and end at the same location but with different costing? What does that even mean? Find the
    // nearest parking and walk back to where I am? Seems like a plausible use case...
    if (costing_name == "auto_pedestrian") {
      if (locations_size > 2) {
        throw valhalla_exception_t{150, "for auto_pedestrian: " + std::to_string(locations->size())};
      }
      google::protobuf::RepeatedPtrField<Location> start_loc(locations->begin(),
                                                             locations->begin() + 1);
      search_.search(start_loc, mode_costing[static_cast<size_t>(TravelMode::kDrive)]);
      google::protobuf::RepeatedPtrField<Location> end_loc(locations->begin() + 1,
                                                           locations->begin() + 2);
      search_.clear();
      search_.search(end_loc, mode_costing[static_cast<size_t>(TravelMode::kPedestrian)]);
      // merge them again
      locations->at(0).CopyFrom(start_loc.at(0));
      locations->at(1).CopyFrom(end_loc.at(0));
    } else {
      search_.search(*locations, mode_costing[static_cast<size_t>(mode)]);
    }

    // throw if there's a location we did not find any candidates for. Generated
    // round-trip anchors are exempt: one projected outside the tile coverage or
    // far from any road is dropped below, and the remaining anchors still work.
    for (int i = 0; i < locations->size(); ++i) {
      const auto is_generated_anchor =
          i >= original_location_count &&
          i < original_location_count + static_cast<int>(generated_candidate_count);
      const auto& location = locations->at(i);
      if (!is_generated_anchor && location.correlation().edges().empty() &&
          location.correlation().filtered_edges().empty()) {
        throw valhalla_exception_t(171);
      }
    }

    if (connectivity_map) {
      for (int i = 0; i < locations_size; ++i) {
        auto colors =
            connectivity_map->get_colors(connectivity_level, locations->at(i), connectivity_radius);
        for (auto color : colors) {
          auto itr = color_counts.find(color);
          if (itr == color_counts.cend()) {
            color_counts[color] = 1;
          } else {
            ++itr->second;
          }
        }
      }
    }

    if (generated_candidate_count) {
      auto* round_trip = options.mutable_round_trip_options();
      for (size_t i = 0; i < generated_candidate_count; ++i) {
        auto& candidate = locations->at(original_location_count + i);
        auto* edges = candidate.mutable_correlation()->mutable_edges();
        if (edges->empty())
          edges->MergeFrom(candidate.correlation().filtered_edges());
        auto new_end = std::remove_if(edges->begin(), edges->end(), [&](const PathEdge& edge) {
          auto tile = reader->GetGraphTile(GraphId(edge.graph_id()));
          return !tile || tile->directededge(GraphId(edge.graph_id()))->use() == Use::kCuldesac;
        });
        edges->DeleteSubrange(std::distance(edges->begin(), new_end),
                              std::distance(new_end, edges->end()));
        if (!edges->empty())
          round_trip->add_anchors()->CopyFrom(candidate);
      }
      // Remove generated anchors and the duplicate origin. Thor consumes the
      // correlated anchors from RoundTripOptions instead.
      locations->DeleteSubrange(original_location_count, locations->size() - original_location_count);
      if (round_trip->anchors().empty())
        throw valhalla_exception_t{171};
      return;
    }
  } catch (const valhalla_exception_t& e) { throw e; } catch (const std::exception&) {
    throw valhalla_exception_t{171};
  }

  // are all the locations in the same color regions
  if (!connectivity_map) {
    return;
  }
  bool connected = false;
  for (const auto& c : color_counts) {
    if (c.second == static_cast<size_t>(options.locations_size())) {
      connected = true;
      break;
    }
  }
  if (!connected) {
    throw valhalla_exception_t{170};
  };
}
} // namespace loki
} // namespace valhalla
