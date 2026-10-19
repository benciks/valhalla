#include "baldr/graphreader.h"
#include "loki/search.h"
#include "loki/worker.h"
#include "midgard/constants.h"

#include <algorithm>
#include <cmath>

using namespace valhalla;
using namespace valhalla::baldr;
using namespace valhalla::midgard;

namespace {
constexpr double kEarthRadiusMeters = 6371000.;
constexpr double kGoldenAngleDegrees = 137.50776405003785;
// How much longer than the straight line the roads to an anchor are assumed to be. Thor can
// shorten a loop that comes out too long but can't stretch one that is too short, so these lean
// low and alternate to bracket more and less direct road networks.
constexpr double kEstimatedRoadStretch[] = {1.18, 1.28, 1.10};
// rough average driving speed, only to turn a time target into a distance for placing anchors
constexpr float kDrivingSpeedKph = 55.f;
constexpr uint32_t kAnchorReachability = 50;
constexpr uint32_t kMaxAnchors = 24;

PointLL project(const PointLL& origin, double bearing_degrees, double meters) {
  const auto bearing = bearing_degrees * kRadPerDegD;
  const auto lat = origin.lat() * kRadPerDegD;
  const auto angular = meters / kEarthRadiusMeters;
  const auto projected_lat = std::asin(std::sin(lat) * std::cos(angular) +
                                       std::cos(lat) * std::sin(angular) * std::cos(bearing));
  const auto projected_lon = origin.lng() * kRadPerDegD +
                             std::atan2(std::sin(bearing) * std::sin(angular) * std::cos(lat),
                                        std::cos(angular) - std::sin(lat) * std::sin(projected_lat));
  return {projected_lon * kDegPerRadD, projected_lat * kDegPerRadD};
}

float average_speed_kph(const Options& options) {
  const auto& costing = options.costings().find(options.costing_type())->second.options();
  switch (options.costing_type()) {
    case Costing::pedestrian:
      return costing.walking_speed();
    case Costing::bicycle:
      return costing.cycling_speed();
    case Costing::motor_scooter:
      return std::min(static_cast<float>(costing.top_speed()), kDrivingSpeedKph);
    default:
      return kDrivingSpeedKph;
  }
}

PointLL ll(const Location& location) {
  return {location.ll().lng(), location.ll().lat()};
}

// Places the anchors so that a straight line itinerary out through the user's vias to the anchor
// and back to the origin would cover the target
google::protobuf::RepeatedPtrField<Location>
place_anchors(const Options& options, double target_meters, uint32_t count) {
  const auto& round_trip = options.round_trip_options();
  const auto origin = ll(options.locations(0));
  const auto last = ll(*options.locations().rbegin());
  double prefix_meters = 0.;
  for (int i = 1; i < options.locations_size(); ++i)
    prefix_meters += ll(options.locations(i - 1)).Distance(ll(options.locations(i)));

  const auto has_vias = options.locations_size() > 1;
  const auto rotation = static_cast<double>(round_trip.seed() % 360);
  const auto outward_bearing = has_vias ? origin.Heading(last) : rotation;
  const auto side = (round_trip.seed() & 1) == 0 ? 1. : -1.;
  google::protobuf::RepeatedPtrField<Location> anchors;
  for (uint32_t i = 0; i < count; ++i) {
    // with vias the loop carries on sideways from the last one, otherwise anchors spread evenly
    const auto bearing =
        has_vias ? std::fmod(outward_bearing + side * (105. + kGoldenAngleDegrees * i) + 720., 360.)
                 : std::fmod(rotation + kGoldenAngleDegrees * i, 360.);
    const auto target_air_meters =
        target_meters / kEstimatedRoadStretch[i % std::size(kEstimatedRoadStretch)];
    const auto minimum_anchor_legs = last.Distance(origin) + std::max(1000., target_meters * .04);
    const auto desired_anchor_legs = std::max(minimum_anchor_legs, target_air_meters - prefix_meters);
    // the distance out to the anchor plus the one from it back to the last via
    double low = 500.;
    double high = std::clamp(target_meters, 5000., 2500000.);
    for (size_t iteration = 0; iteration < 28; ++iteration) {
      const auto radius = (low + high) * .5;
      if (radius + last.Distance(project(origin, bearing, radius)) < desired_anchor_legs)
        low = radius;
      else
        high = radius;
    }
    const auto anchor_ll = project(origin, bearing, (low + high) * .5);
    auto* anchor = anchors.Add();
    anchor->mutable_ll()->set_lng(anchor_ll.lng());
    anchor->mutable_ll()->set_lat(anchor_ll.lat());
    anchor->set_type(Location::kBreak);
    anchor->set_minimum_reachability(kAnchorReachability);
  }
  return anchors;
}

} // namespace

namespace valhalla {
namespace loki {

void loki_worker_t::round_trip(Api& request) {
  auto& options = *request.mutable_options();
  if ((options.round_trip_options().target_distance() > 0) ==
      (options.round_trip_options().target_time() > 0))
    throw valhalla_exception_t{116};
  if (options.costing_type() == Costing::multimodal || options.costing_type() == Costing::transit ||
      options.costing_type() == Costing::bikeshare ||
      options.costing_type() == Costing::auto_pedestrian)
    throw valhalla_exception_t{140, Options_Action_Enum_Name(options.action())};
  if (static_cast<size_t>(options.locations_size()) > max_round_trip_locations)
    throw valhalla_exception_t{150, " for round_trip: " + std::to_string(max_round_trip_locations)};

  // the origin and the vias are located like those of any route
  route(request);

  auto& round_trip = *options.mutable_round_trip_options();
  const auto target_km =
      round_trip.target_distance() > 0
          ? round_trip.target_distance() * (options.units() == Options::miles ? kKmPerMile : 1.f)
          : round_trip.target_time() / 60.f * average_speed_kph(options);
  if (target_km > max_round_trip_distance_km)
    throw valhalla_exception_t{154, " for round_trip: " +
                                        std::to_string(static_cast<int>(max_round_trip_distance_km)) +
                                        " km"};

  // every anchor costs two graph searches over the whole loop, so short loops can afford more
  // of them and a stricter overlap limit while very long ones get fewer
  const auto count = std::min(target_km < 100 ? 5u : target_km >= 500 ? 2u : 3u, kMaxAnchors);
  round_trip.set_max_overlap(target_km < 100 ? .005f : .02f);

  auto anchors = place_anchors(options, target_km * 1000., count);
  for (auto& anchor : anchors)
    parse_location(anchor);
  search_.search(anchors, mode_costing[static_cast<size_t>(mode)]);
  search_.clear();

  // an anchor off the tiles or far from any road is just dropped, as is one on a cul-de-sac
  round_trip.clear_anchors();
  for (auto& anchor : anchors) {
    auto* edges = anchor.mutable_correlation()->mutable_edges();
    if (edges->empty())
      edges->MergeFrom(anchor.correlation().filtered_edges());
    edges->erase(std::remove_if(edges->begin(), edges->end(),
                                [&](const PathEdge& edge) {
                                  const auto* directededge =
                                      reader->directededge(GraphId(edge.graph_id()));
                                  return !directededge || directededge->use() == Use::kCuldesac;
                                }),
                 edges->end());
    if (!edges->empty())
      round_trip.add_anchors()->Swap(&anchor);
  }
  if (round_trip.anchors().empty())
    throw valhalla_exception_t{171, " for any round trip anchor"};
}

} // namespace loki
} // namespace valhalla
