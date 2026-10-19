#include "gurka.h"

#include <gtest/gtest.h>

using namespace valhalla;

class UseCurvatureTest : public ::testing::Test {
protected:
  static gurka::map map;

  static void SetUpTestSuite() {
    constexpr double gridsize = 50;

    // the zig zag way around is 60% longer than the straight one
    const std::string ascii_map = R"(
      A----------B
      |          |
      C a c e g  D
       b d f h i
    )";

    const gurka::ways ways = {{"AB", {{"highway", "residential"}}},
                              {"ACbadcfehgiDB", {{"highway", "residential"}}}};

    const auto layout = gurka::detail::map_to_coordinates(ascii_map, gridsize);
    map = gurka::buildtiles(layout, ways, {}, {}, "test/data/use_curvature");
  }
};

gurka::map UseCurvatureTest::map = {};

TEST_F(UseCurvatureTest, StraightByDefault) {
  const auto result = gurka::do_action(Options::route, map, {"A", "B"}, "motorcycle");
  gurka::assert::raw::expect_path(result, {"AB"});
}

TEST_F(UseCurvatureTest, CurvyWhenPreferred) {
  const auto result = gurka::do_action(Options::route, map, {"A", "B"}, "motorcycle",
                                       {{"/costing_options/motorcycle/use_curvature", "1"}});
  gurka::assert::raw::expect_path(result, {"ACbadcfehgiDB"});
}
