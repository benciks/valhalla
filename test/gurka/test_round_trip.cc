#include "exceptions.h"
#include "gurka.h"

#include <gtest/gtest.h>

#include <format>

using namespace valhalla;

class RoundTrip : public ::testing::Test {
protected:
  static gurka::map map;
  static gurka::nodelayout layout;

  static void SetUpTestSuite() {
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
    layout = gurka::detail::map_to_coordinates(ascii_map, grid_size_metres);
    map = gurka::buildtiles(layout, ways, {}, {}, VALHALLA_BUILD_DIR "test/data/round_trip");
  }

  static std::string location(const std::string& node) {
    return std::format(R"({{"lat":{:.6f},"lon":{:.6f}}})", layout.at(node).lat(),
                       layout.at(node).lng());
  }
};

gurka::map RoundTrip::map = {};
gurka::nodelayout RoundTrip::layout = {};

TEST_F(RoundTrip, LoopFromOrigin) {
  const auto request = std::format(
      R"({{"locations":[{}],"target_distance":2.4,"costing":"motorcycle","round_trip_options":{{"seed":3}}}})",
      location("A"));
  std::string response;
  const auto result = gurka::do_action(Options::round_trip, map, request, nullptr, &response);
  ASSERT_EQ(result.trip().routes_size(), 1);
  ASSERT_EQ(result.trip().routes(0).legs_size(), 2);
  ASSERT_EQ(result.options().locations_size(), 1);
  double length = 0.;
  for (const auto& leg : result.directions().routes(0).legs())
    length += leg.summary().length();
  // the only loop through the anchor is once around the block
  EXPECT_NEAR(length, 2.4, 0.1);
  EXPECT_NE(response.find("Found route between points"), std::string::npos);
}

TEST_F(RoundTrip, OutAndBackWhenNoLoopFits) {
  // the block is 2.4 km around, so 0.8 km can only be out and back
  const auto request = std::format(
      R"({{"locations":[{}],"target_distance":0.8,"costing":"motorcycle","round_trip_options":{{"seed":3}}}})",
      location("A"));
  const auto result = gurka::do_action(Options::round_trip, map, request);
  ASSERT_EQ(result.trip().routes(0).legs_size(), 2);
  double length = 0.;
  for (const auto& leg : result.directions().routes(0).legs())
    length += leg.summary().length();
  EXPECT_NEAR(length, 0.8, 0.35);
}

TEST_F(RoundTrip, ViasAreKept) {
  const auto request = std::format(
      R"({{"locations":[{},{}],"target_distance":1.0,"costing":"motorcycle","round_trip_options":{{"seed":3}}}})",
      location("A"), location("C"));
  const auto result = gurka::do_action(Options::round_trip, map, request);
  ASSERT_EQ(result.trip().routes_size(), 1);
  ASSERT_EQ(result.trip().routes(0).legs_size(), 3);
  ASSERT_EQ(result.options().locations_size(), 2);
}

TEST_F(RoundTrip, RouteIgnoresTargets) {
  const auto request =
      std::format(R"({{"locations":[{},{}],"target_distance":5,"costing":"motorcycle"}})",
                  location("A"), location("C"));
  const auto result = gurka::do_action(Options::route, map, request);
  ASSERT_EQ(result.trip().routes(0).legs_size(), 1);
  EXPECT_FALSE(result.options().has_round_trip_options());
}

TEST_F(RoundTrip, ExactlyOneTarget) {
  for (const auto* targets : {R"("target_distance":1,"target_time":10,)", ""}) {
    const auto request =
        std::format(R"({{"locations":[{}],{}"costing":"motorcycle"}})", location("A"), targets);
    try {
      gurka::do_action(Options::round_trip, map, request);
      FAIL() << "expected an exception";
    } catch (const valhalla_exception_t& e) { EXPECT_EQ(e.code, 116); }
  }
}
