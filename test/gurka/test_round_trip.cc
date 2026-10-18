#include "gurka.h"

#include <gtest/gtest.h>

#include <format>

using namespace valhalla;

TEST(Standalone, GraphNativeRoundTrip) {
  constexpr double grid_size_metres = 100.;
  const std::string ascii_map = R"(
    A---B---C
    |       |
    H       D
    |       |
    G---F---E
  )";
  const gurka::ways ways = {
      {"ABC", {{"highway", "residential"}}},
      {"CDE", {{"highway", "residential"}}},
      {"EFG", {{"highway", "residential"}}},
      {"GHA", {{"highway", "residential"}}},
  };
  const auto layout = gurka::detail::map_to_coordinates(ascii_map, grid_size_metres);
  const auto map = gurka::buildtiles(layout, ways, {}, {}, VALHALLA_BUILD_DIR "test/data/round_trip");
  const auto request = std::format(
      R"({{"locations":[{{"lat":{:.6f},"lon":{:.6f}}}],"target_distance":0.8,"costing":"motorcycle","round_trip_options":{{"candidates":8,"seed":3,"max_overlap":0.02}}}})",
      layout.at("A").lat(), layout.at("A").lng());

  std::string response;
  const auto result = gurka::do_action(Options::route, map, request, nullptr, &response);
  ASSERT_EQ(result.trip().routes_size(), 1);
  ASSERT_EQ(result.trip().routes(0).legs_size(), 2);
  ASSERT_EQ(result.options().locations_size(), 1);
  double length = 0.;
  for (const auto& leg : result.directions().routes(0).legs())
    length += leg.summary().length();
  EXPECT_NEAR(length, 0.8, 0.35);
  EXPECT_NE(response.find("Found route between points"), std::string::npos);

  // An explicit via is a fixed prefix. Thor should preserve it in the result
  // while reusing its already-routed leg across anchor candidates.
  const auto via_request = std::format(
      R"({{"locations":[{{"lat":{:.6f},"lon":{:.6f}}},{{"lat":{:.6f},"lon":{:.6f}}}],"target_distance":1.0,"costing":"motorcycle","round_trip_options":{{"candidates":8,"seed":3,"max_overlap":0.1}}}})",
      layout.at("A").lat(), layout.at("A").lng(), layout.at("C").lat(), layout.at("C").lng());
  const auto via_result = gurka::do_action(Options::route, map, via_request);
  ASSERT_EQ(via_result.trip().routes_size(), 1);
  ASSERT_EQ(via_result.trip().routes(0).legs_size(), 3);
  ASSERT_EQ(via_result.options().locations_size(), 2);
}
