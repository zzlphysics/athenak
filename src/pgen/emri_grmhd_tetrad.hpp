#ifndef PGEN_EMRI_GRMHD_TETRAD_HPP_
#define PGEN_EMRI_GRMHD_TETRAD_HPP_
//========================================================================================
// AthenaK astrophysical plasma code
// Copyright(C) 2026 AthenaK collaboration
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file emri_grmhd_tetrad.hpp
//! \brief Device-portable Eulerian orthonormal frame aligned with a Cartesian face.

#include <cmath>

#include "athena.hpp"

namespace emri_grmhd {

KOKKOS_INLINE_FUNCTION
void ExpandSpatialMetric(const Real compact[NSPMETRIC], Real metric[3][3]) {
  metric[0][0] = compact[S11];
  metric[0][1] = compact[S12];
  metric[0][2] = compact[S13];
  metric[1][0] = compact[S12];
  metric[1][1] = compact[S22];
  metric[1][2] = compact[S23];
  metric[2][0] = compact[S13];
  metric[2][1] = compact[S23];
  metric[2][2] = compact[S33];
}

//! Build e_(a)^i and theta^(a)_i with e_(0) the outward unit face normal.
//!
//! The first tangent follows the next cyclic coordinate direction.  The second is
//! Gram-Schmidt orthogonalized and signed so (normal,tangent1,tangent2) is right handed.
//! Coordinate vectors transform as V^(a)=theta^(a)_i V^i and
//! V^i=e_(a)^i V^(a).
KOKKOS_INLINE_FUNCTION
bool BuildFaceFrame(const Real compact[NSPMETRIC], int normal_axis, int normal_sign,
                    Real basis[3][3], Real dual[3][3], Real &sqrt_determinant,
                    Real &sqrt_inverse_normal_metric) {
  if (normal_axis < 0 || normal_axis > 2 ||
      (normal_sign != -1 && normal_sign != 1)) return false;
  Real metric[3][3];
  ExpandSpatialMetric(compact, metric);
  const Real determinant =
      metric[0][0]*(metric[1][1]*metric[2][2] - metric[1][2]*metric[2][1])
      - metric[0][1]*(metric[1][0]*metric[2][2] - metric[1][2]*metric[2][0])
      + metric[0][2]*(metric[1][0]*metric[2][1] - metric[1][1]*metric[2][0]);
  if (!(determinant > 0.0) || !isfinite(determinant)) return false;
  const Real inverse_determinant = 1.0/determinant;
  Real inverse[3][3];
  inverse[0][0] = (metric[1][1]*metric[2][2] - metric[1][2]*metric[2][1])
                  *inverse_determinant;
  inverse[0][1] = (metric[0][2]*metric[2][1] - metric[0][1]*metric[2][2])
                  *inverse_determinant;
  inverse[0][2] = (metric[0][1]*metric[1][2] - metric[0][2]*metric[1][1])
                  *inverse_determinant;
  inverse[1][0] = inverse[0][1];
  inverse[1][1] = (metric[0][0]*metric[2][2] - metric[0][2]*metric[2][0])
                  *inverse_determinant;
  inverse[1][2] = (metric[0][2]*metric[1][0] - metric[0][0]*metric[1][2])
                  *inverse_determinant;
  inverse[2][0] = inverse[0][2];
  inverse[2][1] = inverse[1][2];
  inverse[2][2] = (metric[0][0]*metric[1][1] - metric[0][1]*metric[1][0])
                  *inverse_determinant;
  if (!(inverse[normal_axis][normal_axis] > 0.0)) return false;
  sqrt_determinant = sqrt(determinant);
  sqrt_inverse_normal_metric = sqrt(inverse[normal_axis][normal_axis]);

  for (int component = 0; component < 3; ++component) {
    basis[0][component] = normal_sign*inverse[component][normal_axis]
                          /sqrt_inverse_normal_metric;
    basis[1][component] = 0.0;
    basis[2][component] = 0.0;
  }
  const int tangent1 = (normal_axis + 1)%3;
  const int tangent2 = (normal_axis + 2)%3;
  if (!(metric[tangent1][tangent1] > 0.0)) return false;
  basis[1][tangent1] = 1.0/sqrt(metric[tangent1][tangent1]);
  basis[2][tangent2] = normal_sign;
  Real projection = 0.0;
  for (int first = 0; first < 3; ++first) {
    for (int second = 0; second < 3; ++second) {
      projection += basis[2][first]*metric[first][second]*basis[1][second];
    }
  }
  for (int component = 0; component < 3; ++component) {
    basis[2][component] -= projection*basis[1][component];
  }
  Real tangent2_norm_squared = 0.0;
  for (int first = 0; first < 3; ++first) {
    for (int second = 0; second < 3; ++second) {
      tangent2_norm_squared +=
          basis[2][first]*metric[first][second]*basis[2][second];
    }
  }
  if (!(tangent2_norm_squared > 0.0) || !isfinite(tangent2_norm_squared)) return false;
  const Real inverse_tangent2_norm = 1.0/sqrt(tangent2_norm_squared);
  for (int component = 0; component < 3; ++component) {
    basis[2][component] *= inverse_tangent2_norm;
  }
  for (int frame = 0; frame < 3; ++frame) {
    for (int coordinate = 0; coordinate < 3; ++coordinate) {
      dual[frame][coordinate] = 0.0;
      for (int component = 0; component < 3; ++component) {
        dual[frame][coordinate] += metric[coordinate][component]
                                   *basis[frame][component];
      }
      if (!isfinite(basis[frame][coordinate]) ||
          !isfinite(dual[frame][coordinate])) return false;
    }
  }
  return true;
}

KOKKOS_INLINE_FUNCTION
void CoordinateToFrame(const Real dual[3][3], const Real coordinate[3], Real frame[3]) {
  for (int local = 0; local < 3; ++local) {
    frame[local] = 0.0;
    for (int component = 0; component < 3; ++component) {
      frame[local] += dual[local][component]*coordinate[component];
    }
  }
}

KOKKOS_INLINE_FUNCTION
void FrameToCoordinate(const Real basis[3][3], const Real frame[3], Real coordinate[3]) {
  for (int component = 0; component < 3; ++component) {
    coordinate[component] = 0.0;
    for (int local = 0; local < 3; ++local) {
      coordinate[component] += basis[local][component]*frame[local];
    }
  }
}

//! Outward physical velocity of a fixed-coordinate face relative to Eulerian observers.
KOKKOS_INLINE_FUNCTION
Real OutwardGridSpeed(Real lapse, Real normal_shift, int normal_sign,
                      Real sqrt_inverse_normal_metric) {
  return normal_sign*normal_shift/(lapse*sqrt_inverse_normal_metric);
}

}  // namespace emri_grmhd

#endif  // PGEN_EMRI_GRMHD_TETRAD_HPP_
