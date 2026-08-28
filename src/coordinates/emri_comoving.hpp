#ifndef COORDINATES_EMRI_COMOVING_HPP_
#define COORDINATES_EMRI_COMOVING_HPP_

//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file emri_comoving.hpp
//! \brief Helically symmetric local coordinates centered on an EMRI secondary.
//!
//! The external spacetime is a stationary primary Kerr-Schild term pulled back into
//! coordinates that rotate about the primary and keep the secondary at x^i=0.  By
//! default, the secondary Kerr perturbation is constructed in an orthonormal tangent
//! frame of that external metric at the orbit and mapped back as a covariant tensor.
//! This avoids interpreting a nearly constant primary potential as a global
//! special-relativistic boost of the small hole.  A legacy global-boost embedding is
//! retained as a controlled diagnostic.  Modes can omit the primary term and/or the
//! orbital pullback, and spins are restricted to the orbital z axis in this first model.
//!
//! If X^i are asymptotically inertial Cartesian coordinates and x^i are local
//! coordinates, the map at t=0 is
//!
//!   X^i = (R_0 + x, y, z),
//!   dX^i/dt = Omega (-y, R_0+x, 0),
//!
//! with the spatial axes rotating at Omega.  Axisymmetry and aligned spins make the
//! pulled-back metric independent of local time.  Evaluating at t=0 avoids subtracting
//! the large orbital radius from the secondary position, which is important at extreme
//! mass ratios.  The tangent embedding is a local effective metric, not a matched
//! asymptotic or constraint-solved binary spacetime.

#include "athena.hpp"
#include "coordinates/superposed_bbh.hpp"

namespace emri_comoving {

struct MetricParameters {
  Real primary_mass;
  Real secondary_mass;
  Real primary_spin;       // Kerr a=S/M, with dimensions of length
  Real secondary_spin;     // Kerr a=S/m, with dimensions of length
  Real orbital_radius;     // primary Kerr-Schild/Boyer-Lindquist r
  Real coordinate_radius;  // equatorial Cartesian KS radius sqrt(r^2+a^2)
  Real omega;
  Real spin_buffer_primary;
  Real spin_buffer_secondary;
  Real singularity_floor;
  Real secondary_tetrad_covector[4][4]; // source-frame co-basis at the orbit
  bool include_primary;       // include the primary Kerr-Schild contribution
  bool include_orbital_frame; // apply the co-rotating/translated coordinate pullback
  bool embed_secondary_in_tetrad;
};

//! Test-particle circular equatorial Kerr angular frequency.  direction=+1 and -1 select
//! positive and negative coordinate angular velocity, respectively.
KOKKOS_INLINE_FUNCTION
Real CircularKerrOmega(const Real primary_mass, const Real primary_spin,
                       const Real orbital_radius, const Real direction) {
  const Real sqrt_mass = sqrt(primary_mass);
  return direction*sqrt_mass
      /(orbital_radius*sqrt(orbital_radius)
        + direction*primary_spin*sqrt_mass);
}

KOKKOS_INLINE_FUNCTION
void PullBackMetric(const Real global_metric[4][4], const Real x, const Real y,
                    const MetricParameters &parameters, Real local_metric[4][4]) {
  const Real frame_omega = parameters.include_orbital_frame ? parameters.omega : 0.0;
  const Real grid_velocity[3] = {
    -frame_omega*y,
    frame_omega*(parameters.coordinate_radius+x),
    0.0
  };
  const Real jacobian[4][4] = {
    {1.0, 0.0, 0.0, 0.0},
    {grid_velocity[0], 1.0, 0.0, 0.0},
    {grid_velocity[1], 0.0, 1.0, 0.0},
    {grid_velocity[2], 0.0, 0.0, 1.0}
  };
  for (int a=0; a<4; ++a) {
    for (int b=a; b<4; ++b) {
      Real value = 0.0;
      for (int global_a=0; global_a<4; ++global_a) {
        for (int global_b=0; global_b<4; ++global_b) {
          value += global_metric[global_a][global_b]
                 *jacobian[global_a][a]*jacobian[global_b][b];
        }
      }
      local_metric[a][b] = value;
      local_metric[b][a] = value;
    }
  }
}

//! Evaluate Minkowski plus the primary Kerr term in the local rotating chart.
KOKKOS_INLINE_FUNCTION
void ComputeExternalMetric(const Real x, const Real y, const Real z,
                           const MetricParameters &parameters,
                           Real local_metric[4][4]) {
  Real global_metric[4][4] = {
    {-1.0, 0.0, 0.0, 0.0},
    { 0.0, 1.0, 0.0, 0.0},
    { 0.0, 0.0, 1.0, 0.0},
    { 0.0, 0.0, 0.0, 1.0}
  };

  if (parameters.include_primary) {
    const Real zero[3] = {0.0, 0.0, 0.0};
    const Real primary_spin[3] = {0.0, 0.0, parameters.primary_spin};
    const Real primary_point_x = parameters.coordinate_radius+x;
    binary_bh::AddBoostedKerrSchildTerm(
        primary_point_x, y, z, zero, zero, primary_spin, parameters.primary_mass,
        parameters.spin_buffer_primary, parameters.singularity_floor, global_metric);
  }
  PullBackMetric(global_metric, x, y, parameters, local_metric);
}

//! Evaluate the secondary perturbation at a chosen displacement while holding the local
//! rotating-coordinate basis fixed.  This separation is needed when differentiating
//! with respect to the source position rather than the field-point coordinates.
KOKKOS_INLINE_FUNCTION
void ComputeSecondaryMetricPerturbationAtDisplacement(
    const Real displacement_x, const Real displacement_y, const Real displacement_z,
    const Real basis_x, const Real basis_y, const MetricParameters &parameters,
    Real local_metric[4][4]) {
  if (parameters.embed_secondary_in_tetrad) {
    Real tetrad_position[3] = {0.0, 0.0, 0.0};
    const Real displacement[4] = {
      0.0, displacement_x, displacement_y, displacement_z
    };
    for (int axis=0; axis<3; ++axis) {
      for (int component=0; component<4; ++component) {
        tetrad_position[axis] +=
            parameters.secondary_tetrad_covector[axis+1][component]
            *displacement[component];
      }
    }
    Real tetrad_metric[4][4] = {{0.0}};
    const Real zero[3] = {0.0, 0.0, 0.0};
    const Real secondary_spin[3] = {0.0, 0.0, parameters.secondary_spin};
    binary_bh::AddBoostedKerrSchildTerm(
        tetrad_position[0], tetrad_position[1], tetrad_position[2], zero, zero,
        secondary_spin, parameters.secondary_mass, parameters.spin_buffer_secondary,
        parameters.singularity_floor, tetrad_metric);
    for (int mu=0; mu<4; ++mu) {
      for (int nu=mu; nu<4; ++nu) {
        Real value = 0.0;
        for (int leg_a=0; leg_a<4; ++leg_a) {
          for (int leg_b=0; leg_b<4; ++leg_b) {
            value += tetrad_metric[leg_a][leg_b]
                *parameters.secondary_tetrad_covector[leg_a][mu]
                *parameters.secondary_tetrad_covector[leg_b][nu];
          }
        }
        local_metric[mu][nu] = value;
        local_metric[nu][mu] = value;
      }
    }
    return;
  }

  Real global_metric[4][4] = {{0.0}};
  const Real zero[3] = {0.0, 0.0, 0.0};
  // Translation invariance lets the secondary term be evaluated directly from its
  // displacement (x,y,z), avoiding a loss of precision from (R_0+x)-R_0 when q << 1.
  const Real orbital_speed = parameters.include_orbital_frame
      ? parameters.omega*parameters.coordinate_radius : 0.0;
  const Real secondary_velocity[3] = {0.0, orbital_speed, 0.0};
  const Real secondary_spin[3] = {0.0, 0.0, parameters.secondary_spin};
  binary_bh::AddBoostedKerrSchildTerm(
      displacement_x, displacement_y, displacement_z, zero, secondary_velocity,
      secondary_spin, parameters.secondary_mass, parameters.spin_buffer_secondary,
      parameters.singularity_floor, global_metric);
  PullBackMetric(global_metric, basis_x, basis_y, parameters, local_metric);
}

//! Evaluate only the secondary Kerr-Schild perturbation in the local rotating chart.
KOKKOS_INLINE_FUNCTION
void ComputeSecondaryMetricPerturbation(const Real x, const Real y, const Real z,
                                        const MetricParameters &parameters,
                                        Real local_metric[4][4]) {
  ComputeSecondaryMetricPerturbationAtDisplacement(
      x, y, z, x, y, parameters, local_metric);
}

//! Evaluate the complete covariant effective metric in the local rotating chart.
KOKKOS_INLINE_FUNCTION
void ComputeMetric(const Real x, const Real y, const Real z,
                   const MetricParameters &parameters, Real local_metric[4][4]) {
  Real external_metric[4][4];
  Real secondary_perturbation[4][4];
  ComputeExternalMetric(x, y, z, parameters, external_metric);
  ComputeSecondaryMetricPerturbation(x, y, z, parameters, secondary_perturbation);
  for (int a=0; a<4; ++a) {
    for (int b=0; b<4; ++b) {
      local_metric[a][b] = external_metric[a][b]+secondary_perturbation[a][b];
    }
  }
}

//! Small-hole Kerr radius on a local t=constant slice, evaluated either through the
//! source tangent co-basis or through the legacy global orbital boost.
KOKKOS_INLINE_FUNCTION
Real SecondaryKerrRadiusSquared(const Real x, const Real y, const Real z,
                                const MetricParameters &parameters) {
  if (parameters.embed_secondary_in_tetrad) {
    const Real displacement[4] = {0.0, x, y, z};
    Real tetrad_position[3] = {0.0, 0.0, 0.0};
    for (int axis=0; axis<3; ++axis) {
      for (int component=0; component<4; ++component) {
        tetrad_position[axis] +=
            parameters.secondary_tetrad_covector[axis+1][component]
            *displacement[component];
      }
    }
    const Real spin[3] = {0.0, 0.0, parameters.secondary_spin};
    return binary_bh::KerrRadiusSquared(tetrad_position, spin);
  }
  const Real position[3] = {0.0, 0.0, 0.0};
  const Real orbital_speed = parameters.include_orbital_frame
      ? parameters.omega*parameters.coordinate_radius : 0.0;
  const Real velocity[3] = {0.0, orbital_speed, 0.0};
  const Real spin[3] = {0.0, 0.0, parameters.secondary_spin};
  Real rest_position[3];
  binary_bh::RestFramePosition(x, y, z, position, velocity, rest_position);
  return binary_bh::KerrRadiusSquared(rest_position, spin);
}

KOKKOS_INLINE_FUNCTION
void SecondaryRestFramePosition(const Real x, const Real y, const Real z,
                                const MetricParameters &parameters,
                                Real rest_position[3]) {
  if (parameters.embed_secondary_in_tetrad) {
    const Real displacement[4] = {0.0, x, y, z};
    for (int axis=0; axis<3; ++axis) {
      rest_position[axis] = 0.0;
      for (int component=0; component<4; ++component) {
        rest_position[axis] +=
            parameters.secondary_tetrad_covector[axis+1][component]
            *displacement[component];
      }
    }
    return;
  }
  const Real position[3] = {0.0, 0.0, 0.0};
  const Real orbital_speed = parameters.include_orbital_frame
      ? parameters.omega*parameters.coordinate_radius : 0.0;
  const Real velocity[3] = {0.0, orbital_speed, 0.0};
  binary_bh::RestFramePosition(x, y, z, position, velocity, rest_position);
}

KOKKOS_INLINE_FUNCTION
Real SecondaryRegularizationRadius(const MetricParameters &parameters) {
  const Real spin[3] = {0.0, 0.0, parameters.secondary_spin};
  return binary_bh::RegularizationRadius(
      spin, parameters.secondary_mass, parameters.spin_buffer_secondary,
      parameters.singularity_floor);
}

} // namespace emri_comoving

#endif // COORDINATES_EMRI_COMOVING_HPP_
