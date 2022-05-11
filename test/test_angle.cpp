#include "cit.hpp"

#include "gtest/gtest.h"

TEST(AngleTest, AbsDifference) {
  static const double s = M_PI / 180.;
  static const double e = 1.e-10;
  EXPECT_LT(cit::angleAbsDifference(170 * s, -170 * s), 20 * s + e);
  EXPECT_LT(cit::angleAbsDifference(-170 * s, -10 * s), 160 * s + e);
}
