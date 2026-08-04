#ifndef COORDINATES_SUPERPOSED_BBH_HPP_
#define COORDINATES_SUPERPOSED_BBH_HPP_

//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file superposed_bbh.hpp
//! \brief Boosted, superposed Kerr-Schild binary-black-hole metric.
//!
//! This implements Eq. (11) of Combi & Ressler, "A binary black hole metric
//! approximation from inspiral to merger", arXiv:2403.13308.  Orbital dynamics and the
//! transition to a remnant are deliberately represented by the supplied trajectory state;
//! the metric evaluation itself is independent of how that state was generated.

#include "athena.hpp"

namespace binary_bh {

enum TrajectoryIndex {
  X1, Y1, Z1, X2, Y2, Z2,
  VX1, VY1, VZ1, VX2, VY2, VZ2,
  AX1, AY1, AZ1, AX2, AY2, AZ2,
  M1T, M2T, NTRAJ
};

//! Per-hole mass scaling and parameters for the regularized region.
struct MetricParameters {
  Real mass_scale1;
  Real mass_scale2;
  Real spin_buffer1;
  Real spin_buffer2;
  Real singularity_floor;
};

KOKKOS_INLINE_FUNCTION
void RestFramePosition(const Real x, const Real y, const Real z,
                       const Real position[3], const Real velocity[3], Real x_bh[3]) {
  const Real v2 = SQR(velocity[0]) + SQR(velocity[1]) + SQR(velocity[2]);
  const Real gamma = 1.0/sqrt(1.0 - v2);
  const Real dx[3] = {x-position[0], y-position[1], z-position[2]};
  const Real vdotx = velocity[0]*dx[0] + velocity[1]*dx[1] + velocity[2]*dx[2];
  // (gamma-1)/v^2 = gamma^2/(gamma+1) is regular as v -> 0.
  const Real boost_coefficient = gamma*gamma/(gamma + 1.0);
  for (int i=0; i<3; ++i) {
    x_bh[i] = dx[i] + boost_coefficient*vdotx*velocity[i];
  }
}

KOKKOS_INLINE_FUNCTION
Real KerrRadiusSquared(const Real x_bh[3], const Real spin[3]) {
  const Real spin2 = SQR(spin[0]) + SQR(spin[1]) + SQR(spin[2]);
  const Real radius2 = SQR(x_bh[0]) + SQR(x_bh[1]) + SQR(x_bh[2]);
  const Real adotx = spin[0]*x_bh[0] + spin[1]*x_bh[1] + spin[2]*x_bh[2];
  const Real radial_term = radius2 - spin2;
  return 0.5*(radial_term + sqrt(SQR(radial_term) + 4.0*SQR(adotx)));
}

KOKKOS_INLINE_FUNCTION
Real RegularizationRadius(const Real spin[3], const Real mass, const Real spin_buffer,
                          const Real singularity_floor) {
  const Real spin2 = SQR(spin[0]) + SQR(spin[1]) + SQR(spin[2]);
  const Real spin_mag = sqrt(spin2);
  // For a subextremal component, taper only inside its Kerr horizon.  During the
  // canonical merger interpolation an individual half-mass term can be superextremal,
  // even though the two identical terms sum to a regular remnant; the ring-based bound
  // then supplies the required interior cutoff.
  const Real horizon_radius = mass + sqrt(fmax(SQR(mass)-spin2, 0.0));
  // The auxiliary ring floor scales with the component mass, so an extreme mass ratio
  // cannot move the regularization into the causal exterior of the smaller hole.
  const Real ring_radius = spin_mag*spin_buffer + mass*singularity_floor;
  return fmax(horizon_radius, ring_radius);
}

//----------------------------------------------------------------------------------------
//! Add one boosted Kerr-Schild term to a covariant spacetime metric.
//!
//! Position and velocity are given in global Cartesian coordinates.  The Kerr spin
//! components use the instantaneous moving-hole axes of the paper; here those axes are
//! chosen parallel to the global axes before the pure boost, with no additional spatial
//! rotation.  The spin has dimensions of length (a^i = S^i/M), rather than being the
//! dimensionless spin chi^i.  Inside the component's horizon-scale excision region, the
//! Kerr-Schild perturbation is tapered to zero with a compact C2 window before reaching
//! the degenerate r=0 Kerr disk.  The exact metric and its null vector are unchanged
//! outside that region; unlike the coordinate jump in the reference implementation, the
//! taper is continuous across the spin equator and its boundary.

KOKKOS_INLINE_FUNCTION
void AddBoostedKerrSchildTerm(const Real x, const Real y, const Real z,
                             const Real position[3], const Real velocity[3],
                             const Real spin[3], const Real mass,
                             const Real spin_buffer, const Real singularity_floor,
                             Real gcov[4][4]) {
  if (mass == 0.0) return;

  const Real v2 = SQR(velocity[0]) + SQR(velocity[1]) + SQR(velocity[2]);
  if (!(mass > 0.0) || !(v2 < 1.0)) {
    Kokkos::abort("Invalid mass or superluminal velocity in binary-BH trajectory");
  }

  const Real gamma = 1.0/sqrt(1.0 - v2);
  const Real boost_coefficient = gamma*gamma/(gamma + 1.0);
  Real x_bh[3];
  RestFramePosition(x, y, z, position, velocity, x_bh);

  const Real spin2 = SQR(spin[0]) + SQR(spin[1]) + SQR(spin[2]);
  // Arbitrarily oriented Cartesian Kerr-Schild metric, Eqs. (1)--(4).
  const Real adotx = spin[0]*x_bh[0] + spin[1]*x_bh[1] + spin[2]*x_bh[2];
  const Real r2 = KerrRadiusSquared(x_bh, spin);
  const Real regularization_radius = RegularizationRadius(
      spin, mass, spin_buffer, singularity_floor);
  const Real regularization_radius2 = SQR(regularization_radius);
  const Real inner_radius2 = 0.25*regularization_radius2;
  if (r2 <= inner_radius2) return;
  Real taper = 1.0;
  if (r2 < regularization_radius2) {
    const Real u = (r2-inner_radius2)/(regularization_radius2-inner_radius2);
    const Real u2 = u*u;
    const Real u3 = u2*u;
    // Quintic smootherstep: value and first two derivatives match at both ends.
    taper = u3*(10.0 - 15.0*u + 6.0*u2);
  }
  const Real r = sqrt(r2);
  const Real denominator = SQR(r2) + SQR(adotx);
  if (!(r > 0.0) || !(denominator > 0.0)) {
    Kokkos::abort("Binary-BH Kerr-Schild singularity was not regularized");
  }

  const Real hfun = taper*mass*r2*r/denominator;
  const Real spatial_denominator = r2 + spin2;
  // x_bh cross spin implements -epsilon^i_jk a^j X^k in the paper.
  const Real x_cross_a[3] = {
    x_bh[1]*spin[2] - x_bh[2]*spin[1],
    x_bh[2]*spin[0] - x_bh[0]*spin[2],
    x_bh[0]*spin[1] - x_bh[1]*spin[0]
  };
  Real null_covector[4] = {1.0, 0.0, 0.0, 0.0};
  for (int i=0; i<3; ++i) {
    null_covector[i+1] =
        (r*x_bh[i] + x_cross_a[i] + (adotx/r)*spin[i])/spatial_denominator;
  }

  // Instantaneous Lorentz Jacobian, Eq. (6).  It is symmetric, so the same array can be
  // used for Lambda^c_a and Lambda^d_b in the covariant transformation.
  Real lorentz[4][4] = {
    {gamma, 0.0, 0.0, 0.0},
    {0.0, 1.0, 0.0, 0.0},
    {0.0, 0.0, 1.0, 0.0},
    {0.0, 0.0, 0.0, 1.0}
  };
  for (int i=0; i<3; ++i) {
    lorentz[0][i+1] = -gamma*velocity[i];
    lorentz[i+1][0] = lorentz[0][i+1];
  }
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      lorentz[i+1][j+1] += boost_coefficient*velocity[i]*velocity[j];
    }
  }

  Real boosted_null[4] = {0.0, 0.0, 0.0, 0.0};
  for (int a=0; a<4; ++a) {
    for (int c=0; c<4; ++c) boosted_null[a] += lorentz[c][a]*null_covector[c];
  }
  for (int a=0; a<4; ++a) {
    for (int b=a; b<4; ++b) {
      const Real term = 2.0*hfun*boosted_null[a]*boosted_null[b];
      gcov[a][b] += term;
      if (a != b) gcov[b][a] += term;
    }
  }
}

//----------------------------------------------------------------------------------------
//! Evaluate the full superposed binary metric in global Cartesian coordinates.

KOKKOS_INLINE_FUNCTION
void ComputeMetric(const Real x, const Real y, const Real z,
                   const Real trajectory[NTRAJ], const MetricParameters &parameters,
                   Real gcov[4][4]) {
  for (int a=0; a<4; ++a) {
    for (int b=0; b<4; ++b) gcov[a][b] = (a == b) ? 1.0 : 0.0;
  }
  gcov[0][0] = -1.0;

  const Real position1[3] = {trajectory[X1], trajectory[Y1], trajectory[Z1]};
  const Real position2[3] = {trajectory[X2], trajectory[Y2], trajectory[Z2]};
  const Real velocity1[3] = {trajectory[VX1], trajectory[VY1], trajectory[VZ1]};
  const Real velocity2[3] = {trajectory[VX2], trajectory[VY2], trajectory[VZ2]};
  const Real spin1[3] = {
    trajectory[AX1]*parameters.mass_scale1,
    trajectory[AY1]*parameters.mass_scale1,
    trajectory[AZ1]*parameters.mass_scale1
  };
  const Real spin2[3] = {
    trajectory[AX2]*parameters.mass_scale2,
    trajectory[AY2]*parameters.mass_scale2,
    trajectory[AZ2]*parameters.mass_scale2
  };

  AddBoostedKerrSchildTerm(x, y, z, position1, velocity1, spin1,
      trajectory[M1T]*parameters.mass_scale1, parameters.spin_buffer1,
      parameters.singularity_floor, gcov);
  AddBoostedKerrSchildTerm(x, y, z, position2, velocity2, spin2,
      trajectory[M2T]*parameters.mass_scale2, parameters.spin_buffer2,
      parameters.singularity_floor, gcov);
}

} // namespace binary_bh

#endif // COORDINATES_SUPERPOSED_BBH_HPP_
