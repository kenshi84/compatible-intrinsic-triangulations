#include "cit.hpp"

#include "gtest/gtest.h"

#include <boost/rational.hpp>
#include <boost/multiprecision/cpp_int.hpp>

using rational = boost::rational<boost::multiprecision::int512_t>;
using Vector2r = Eigen::Matrix<rational, 2, 1>;
using Matrix2r = Eigen::Matrix<rational, 2, 2>;

inline rational cross(const Vector2r& u, const Vector2r& v) { return u.x() * v.y() - u.y() * v.x(); }
inline std::array<rational, 3> getLengthSquared(const std::array<Vector2r,3>& P) {
  return {
    (P[1] - P[0]).squaredNorm(),
    (P[2] - P[1]).squaredNorm(),
    (P[0] - P[2]).squaredNorm()
  };
}

rational computeJacobianSquaredNorm1(const std::array<Vector2r,3>& A, const std::array<Vector2r,3>& B) {
  Matrix2r A_D;
  Matrix2r B_D;
  A_D << A[1] - A[0], A[2] - A[0];
  B_D << B[1] - B[0], B[2] - B[0];
  Matrix2r J = B_D * A_D.inverse();
  return J.squaredNorm();
}

rational computeJacobianSquaredNorm2(const std::array<rational,3>& from_lensq, const std::array<rational,3>& to_lensq, rational from_area) {
  std::array<rational,3> s;
  for (int i = 0; i < 3; ++i) {
    s[i] = from_lensq[i] * (-to_lensq[i] + to_lensq[(i+1) % 3] + to_lensq[(i+2) % 3]);
  }
  return (s[0] + s[1] + s[2]) / (8 * from_area * from_area);
}

TEST(EnergyTest, CompareUsingRational) {
  int cnt = 0;
  for (; cnt < 1000;) {
    std::array<Vector2r, 3> A, B;
    for (int i = 0; i < 3; ++i) {
      A[i].setRandom();
      B[i].setRandom();
    }
    rational A_area = cross(A[1] - A[0], A[2] - A[0]) / 2;
    rational B_area = cross(B[1] - B[0], B[2] - B[0]) / 2;
    if (A_area < 0) continue;
    if (B_area < 0) continue;
    rational r1 = computeJacobianSquaredNorm1(A, B);
    rational r2 = computeJacobianSquaredNorm2(getLengthSquared(A), getLengthSquared(B), A_area);
    rational e = r2 - r1;
    EXPECT_EQ(e, 0);
    ++cnt;
  }
}
