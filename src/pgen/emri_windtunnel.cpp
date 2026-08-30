//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file emri_windtunnel.cpp
//! \brief Local GRMHD wind tunnel centered on a circular equatorial EMRI secondary.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "coordinates/adm.hpp"
#include "coordinates/cartesian_ks.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/emri_comoving.hpp"
#include "dyn_grmhd/dyn_grmhd.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "mhd/mhd.hpp"
#include "outputs/outputs.hpp"
#include "parameter_input.hpp"
#include "pgen/emri_inner_worldtube.hpp"

namespace {

struct WindTunnelParameters {
  emri_comoving::MetricParameters metric;
  Real metric_fd_step;
  Real external_metric_fd_step;
  Real rho;
  Real pressure;
  Real wind_u[3];
  Real magnetic_field[3];
  Real log_density_gradient[3];
  Real log_pressure_gradient[3];
  Real velocity_gradient[3][3];
  Real magnetic_gradient_source[3][3];
  Real densitized_magnetic_gradient[3][3];
  Real max_log_contrast;
  Real wind_eulerian_u[3];
  Real densitized_magnetic_field[3];
  Real source_tetrad[4][4];
  Real force_projection[3][3];
  Real source_spatial_inverse[3][3];
  Real source_surface_covector[3][3];
  Real source_spatial_determinant;
  Real source_dt_dtau;
  Real source_coordinate_stretch;
  Real source_tetrad_stretch;
  Real secondary_coordinate_stretch;
  Real secondary_rest_stretch;
  Real refinement_radius;
  Real refinement_radius_ratio;
  Real refinement_hysteresis;
  Real refinement_horizon_factor;
  Real force_surface_radius;
  Real force_outer_radius[3];
  Real adiabatic_index;
  int force_surface_nlevel;
  bool force_subtract_background;
  bool reinitialize_wind_on_restart;
  bool wind_is_source_tetrad;
  bool force_is_source_tetrad;
};

struct WindProfileSample {
  Real time;
  Real rho;
  Real pressure;
  Real wind_u[3];
  Real magnetic_field[3];
  Real log_density_gradient[3];
  Real log_pressure_gradient[3];
  Real velocity_gradient[3][3];
  Real magnetic_gradient_source[3][3];
};

struct MetricWithDerivatives {
  Real g[4][4];
  Real dg[4][4][4];  // dg[coordinate][first metric index][second metric index]
};

struct ADMPoint {
  Real gamma[3][3];
  Real curvature[3][3];
  Real alpha;
  Real beta[3];
  Real psi4;
};

WindTunnelParameters wind_tunnel;
std::vector<Real> refinement_shell_radii;
std::vector<WindProfileSample> wind_profile_table;
std::uint64_t wind_profile_hash = 0;
Real wind_profile_time_offset = 0.0;
Real wind_profile_current_table_time = std::numeric_limits<Real>::quiet_NaN();
Real wind_profile_source_coordinate_radius = std::numeric_limits<Real>::quiet_NaN();
Real wind_profile_source_omega = std::numeric_limits<Real>::quiet_NaN();
Real wind_profile_orbit_tolerance = std::numeric_limits<Real>::quiet_NaN();
bool wind_profile_enabled = false;
bool wind_profile_hold_outside = false;
bool omega_is_geodesic = false;
Real primary_geodesic_residual = 0.0;

[[noreturn]] void Fatal(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << "\n" << message << std::endl;
  std::exit(EXIT_FAILURE);
}

void SetADMVariablesToEMRI(MeshBlockPack *pmbp);
void AugmentEMRIExcisionMasks(MeshBlockPack *pmbp);
void RefineSecondary(MeshBlockPack *pmbp);
void EMRIHistory(HistoryData *pdata, Mesh *);
void EMRIWindBoundary(Mesh *pm);
bool EMRIOutputRegionCenter(const std::string &name, Real time, Real center[3]);
void InitializeAnalyticMagneticField(MeshBlockPack *pmbp);

KOKKOS_INLINE_FUNCTION
void EvaluateMetric(const Real x, const Real y, const Real z,
                    const WindTunnelParameters &parameters, Real gcov[4][4]) {
  emri_comoving::ComputeMetric(x, y, z, parameters.metric, gcov);
}

KOKKOS_INLINE_FUNCTION
Real MetricInnerProduct(const Real metric[4][4], const Real left[4],
                        const Real right[4]) {
  Real product = 0.0;
  for (int a=0; a<4; ++a) {
    for (int b=0; b<4; ++b) product += metric[a][b]*left[a]*right[b];
  }
  return product;
}

//! Construct an orthonormal tetrad at the secondary's orbit using only the external
//! background.  Its timelike leg follows the stationary source worldline x^i=0; the
//! spatial legs are Gram-Schmidt aligned with local radial, tangential, and vertical
//! axes.
KOKKOS_INLINE_FUNCTION
void BuildSourceTetrad(const emri_comoving::MetricParameters &parameters,
                       Real tetrad[4][4]) {
  Real external_metric[4][4];
  emri_comoving::ComputeExternalMetric(0.0, 0.0, 0.0, parameters,
                                       external_metric);
  for (int leg=0; leg<4; ++leg) {
    for (int component=0; component<4; ++component) tetrad[leg][component] = 0.0;
  }
  const Real time_norm2 = -external_metric[0][0];
  if (!(time_norm2 > 0.0) || !isfinite(time_norm2)) {
    Kokkos::abort("EMRI source worldline is not timelike in the external metric");
  }
  tetrad[0][0] = 1.0/sqrt(time_norm2);

  for (int spatial_leg=0; spatial_leg<3; ++spatial_leg) {
    Real candidate[4] = {0.0, 0.0, 0.0, 0.0};
    candidate[spatial_leg+1] = 1.0;
    const Real time_overlap =
        MetricInnerProduct(external_metric, candidate, tetrad[0]);
    for (int component=0; component<4; ++component) {
      candidate[component] += time_overlap*tetrad[0][component];
    }
    for (int previous=0; previous<spatial_leg; ++previous) {
      const Real spatial_overlap = MetricInnerProduct(
          external_metric, candidate, tetrad[previous+1]);
      for (int component=0; component<4; ++component) {
        candidate[component] -= spatial_overlap*tetrad[previous+1][component];
      }
    }
    const Real norm2 = MetricInnerProduct(external_metric, candidate, candidate);
    if (!(norm2 > 0.0) || !isfinite(norm2)) {
      Kokkos::abort("failed to construct the EMRI source spatial tetrad");
    }
    for (int component=0; component<4; ++component) {
      tetrad[spatial_leg+1][component] = candidate[component]/sqrt(norm2);
    }
  }
}

//! Build a contravariant orthonormal triad aligned with the coordinate spatial axes.
KOKKOS_INLINE_FUNCTION
void BuildSpatialTriad(const Real metric[4][4], Real triad[3][3]) {
  for (int leg=0; leg<3; ++leg) {
    Real candidate[3] = {0.0, 0.0, 0.0};
    candidate[leg] = 1.0;
    for (int previous=0; previous<leg; ++previous) {
      Real overlap = 0.0;
      for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
          overlap += metric[i+1][j+1]*candidate[i]*triad[previous][j];
        }
      }
      for (int i=0; i<3; ++i) candidate[i] -= overlap*triad[previous][i];
    }
    Real norm2 = 0.0;
    for (int i=0; i<3; ++i) {
      for (int j=0; j<3; ++j) {
        norm2 += metric[i+1][j+1]*candidate[i]*candidate[j];
      }
    }
    if (!(norm2 > 0.0) || !isfinite(norm2)) {
      Kokkos::abort("failed to construct an EMRI spatial triad");
    }
    for (int i=0; i<3; ++i) triad[leg][i] = candidate[i]/sqrt(norm2);
  }
}

KOKKOS_INLINE_FUNCTION
void PrimitiveToFourVelocity(const Real metric[4][4], const Real primitive[3],
                             Real four_velocity[4]) {
  const Real determinant = adm::SpatialDet(
      metric[1][1], metric[1][2], metric[1][3],
      metric[2][2], metric[2][3], metric[3][3]);
  Real inverse_gamma[3][3];
  adm::SpatialInv(1.0/determinant,
      metric[1][1], metric[1][2], metric[1][3],
      metric[2][2], metric[2][3], metric[3][3],
      &inverse_gamma[0][0], &inverse_gamma[0][1], &inverse_gamma[0][2],
      &inverse_gamma[1][1], &inverse_gamma[1][2], &inverse_gamma[2][2]);
  inverse_gamma[1][0] = inverse_gamma[0][1];
  inverse_gamma[2][0] = inverse_gamma[0][2];
  inverse_gamma[2][1] = inverse_gamma[1][2];
  Real beta[3] = {0.0, 0.0, 0.0};
  Real beta_squared = 0.0;
  Real lorentz2 = 1.0;
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      beta[i] += inverse_gamma[i][j]*metric[0][j+1];
      lorentz2 += metric[i+1][j+1]*primitive[i]*primitive[j];
    }
    beta_squared += metric[0][i+1]*beta[i];
  }
  const Real alpha = sqrt(beta_squared-metric[0][0]);
  four_velocity[0] = sqrt(lorentz2)/alpha;
  for (int i=0; i<3; ++i) {
    four_velocity[i+1] = primitive[i]-four_velocity[0]*beta[i];
  }
}

KOKKOS_INLINE_FUNCTION
void FourVelocityToPrimitive(const Real metric[4][4], const Real four_velocity[4],
                             Real primitive[3]) {
  const Real determinant = adm::SpatialDet(
      metric[1][1], metric[1][2], metric[1][3],
      metric[2][2], metric[2][3], metric[3][3]);
  Real inverse_gamma[3][3];
  adm::SpatialInv(1.0/determinant,
      metric[1][1], metric[1][2], metric[1][3],
      metric[2][2], metric[2][3], metric[3][3],
      &inverse_gamma[0][0], &inverse_gamma[0][1], &inverse_gamma[0][2],
      &inverse_gamma[1][1], &inverse_gamma[1][2], &inverse_gamma[2][2]);
  inverse_gamma[1][0] = inverse_gamma[0][1];
  inverse_gamma[2][0] = inverse_gamma[0][2];
  inverse_gamma[2][1] = inverse_gamma[1][2];
  Real beta[3] = {0.0, 0.0, 0.0};
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) beta[i] += inverse_gamma[i][j]*metric[0][j+1];
    primitive[i] = four_velocity[i+1]+four_velocity[0]*beta[i];
  }
}

//! Map a coordinate displacement at fixed numerical time into the spatial source
//! tangent chart.  The axes are radial, prograde tangential, and vertical.
KOKKOS_INLINE_FUNCTION
void ComputeSourceDisplacement(const Real x, const Real y, const Real z,
                               const WindTunnelParameters &parameters,
                               Real source_position[3]) {
  const Real coordinate[3] = {x, y, z};
  for (int source_axis=0; source_axis<3; ++source_axis) {
    source_position[source_axis] = 0.0;
    for (int coordinate_axis=0; coordinate_axis<3; ++coordinate_axis) {
      source_position[source_axis] +=
          parameters.metric.secondary_tetrad_covector[source_axis+1]
                                                        [coordinate_axis+1]
          *coordinate[coordinate_axis];
    }
  }
}

KOKKOS_INLINE_FUNCTION
void ComputeWindThermodynamics(const Real x, const Real y, const Real z,
                               const WindTunnelParameters &parameters,
                               Real &density, Real &pressure) {
  Real source_position[3];
  ComputeSourceDisplacement(x, y, z, parameters, source_position);
  Real log_density_ratio = 0.0;
  Real log_pressure_ratio = 0.0;
  for (int axis=0; axis<3; ++axis) {
    log_density_ratio +=
        parameters.log_density_gradient[axis]*source_position[axis];
    log_pressure_ratio +=
        parameters.log_pressure_gradient[axis]*source_position[axis];
  }
  density = parameters.rho*exp(log_density_ratio);
  pressure = parameters.pressure*exp(log_pressure_ratio);
}

//! Evaluate the analytic densitized magnetic field represented by the discrete curl
//! initial data.  A trace-free gradient makes its coordinate divergence identically
//! zero.  This evaluator is also used to extend the analytic state into user ghosts.
KOKKOS_INLINE_FUNCTION
void ComputeDensitizedMagneticField(const Real x, const Real y, const Real z,
                                    const WindTunnelParameters &parameters,
                                    Real field[3]) {
  const Real coordinate[3] = {x, y, z};
  for (int component=0; component<3; ++component) {
    field[component] = parameters.densitized_magnetic_field[component];
    for (int direction=0; direction<3; ++direction) {
      field[component] +=
          parameters.densitized_magnetic_gradient[component][direction]
          *coordinate[direction];
    }
  }
}

//! Vector potential for B^i_dens = B0^i + G^i_j x^j with tr(G)=0:
//! A = (B0 x x)/2 - x x (G x)/3.  Since AthenaK's CT variable is the densitized
//! magnetic field, this ordinary coordinate curl is the required face flux density.
KOKKOS_INLINE_FUNCTION
void ComputeMagneticVectorPotential(const Real x, const Real y, const Real z,
                                    const WindTunnelParameters &parameters,
                                    Real potential[3]) {
  const Real coordinate[3] = {x, y, z};
  Real gradient_field[3] = {0.0, 0.0, 0.0};
  for (int component=0; component<3; ++component) {
    for (int direction=0; direction<3; ++direction) {
      gradient_field[component] +=
          parameters.densitized_magnetic_gradient[component][direction]
          *coordinate[direction];
    }
  }
  const Real base_cross[3] = {
    parameters.densitized_magnetic_field[1]*coordinate[2]
        - parameters.densitized_magnetic_field[2]*coordinate[1],
    parameters.densitized_magnetic_field[2]*coordinate[0]
        - parameters.densitized_magnetic_field[0]*coordinate[2],
    parameters.densitized_magnetic_field[0]*coordinate[1]
        - parameters.densitized_magnetic_field[1]*coordinate[0]
  };
  const Real position_cross_gradient[3] = {
    coordinate[1]*gradient_field[2]-coordinate[2]*gradient_field[1],
    coordinate[2]*gradient_field[0]-coordinate[0]*gradient_field[2],
    coordinate[0]*gradient_field[1]-coordinate[1]*gradient_field[0]
  };
  for (int component=0; component<3; ++component) {
    potential[component] =
        0.5*base_cross[component]-position_cross_gradient[component]/3.0;
  }
}

//! Convert the common source-frame wind to the normal-frame spatial four-velocity used
//! by dynamical GRMHD.  Source-frame matching is exact at the orbital anchor.  Its
//! Eulerian orthonormal components are held fixed in the source tangent chart and then
//! transformed into the numerical slicing.  This makes rotating-Minkowski and isolated
//! controls differ only by the intended non-inertial gradients.
KOKKOS_INLINE_FUNCTION
void ComputeWindPrimitive(const Real x, const Real y, const Real z,
                          const WindTunnelParameters &parameters,
                          Real primitive_velocity[3]) {
  Real source_position[3];
  ComputeSourceDisplacement(x, y, z, parameters, source_position);
  Real local_wind_u[3];
  for (int component=0; component<3; ++component) {
    local_wind_u[component] = parameters.wind_u[component];
    for (int direction=0; direction<3; ++direction) {
      local_wind_u[component] +=
          parameters.velocity_gradient[component][direction]
          *source_position[direction];
    }
  }
  if (!parameters.wind_is_source_tetrad) {
    for (int i=0; i<3; ++i) primitive_velocity[i] = local_wind_u[i];
    return;
  }

  Real metric[4][4];
  EvaluateMetric(x, y, z, parameters, metric);
  Real tangent_metric[4][4];
  for (int leg_a=0; leg_a<4; ++leg_a) {
    for (int leg_b=0; leg_b<4; ++leg_b) {
      tangent_metric[leg_a][leg_b] = 0.0;
      for (int mu=0; mu<4; ++mu) {
        for (int nu=0; nu<4; ++nu) {
          tangent_metric[leg_a][leg_b] += metric[mu][nu]
              *parameters.source_tetrad[leg_a][mu]
              *parameters.source_tetrad[leg_b][nu];
        }
      }
    }
  }
  Real tangent_triad[3][3];
  BuildSpatialTriad(tangent_metric, tangent_triad);
  Real tangent_primitive[3] = {0.0, 0.0, 0.0};
  for (int axis=0; axis<3; ++axis) {
    for (int leg=0; leg<3; ++leg) {
      tangent_primitive[axis] += local_wind_u[leg]*tangent_triad[leg][axis];
    }
  }
  Real tangent_velocity[4];
  PrimitiveToFourVelocity(tangent_metric, tangent_primitive, tangent_velocity);
  Real coordinate_velocity[4] = {0.0, 0.0, 0.0, 0.0};
  for (int mu=0; mu<4; ++mu) {
    for (int leg=0; leg<4; ++leg) {
      coordinate_velocity[mu] +=
          parameters.source_tetrad[leg][mu]*tangent_velocity[leg];
    }
  }
  FourVelocityToPrimitive(metric, coordinate_velocity, primitive_velocity);
}

KOKKOS_INLINE_FUNCTION
void DifferentiateMetric(const Real x, const Real y, const Real z,
                         const WindTunnelParameters &parameters,
                         MetricWithDerivatives &metric) {
  EvaluateMetric(x, y, z, parameters, metric.g);
  for (int a=0; a<4; ++a) {
    for (int b=0; b<4; ++b) metric.dg[0][a][b] = 0.0;
  }

  Real external_minus[4][4];
  Real external_plus[4][4];
  Real secondary_minus[4][4];
  Real secondary_plus[4][4];
  const Real coordinate[3] = {x, y, z};
  for (int direction=0; direction<3; ++direction) {
    Real external_lower[3] = {coordinate[0], coordinate[1], coordinate[2]};
    Real external_upper[3] = {coordinate[0], coordinate[1], coordinate[2]};
    Real secondary_lower[3] = {coordinate[0], coordinate[1], coordinate[2]};
    Real secondary_upper[3] = {coordinate[0], coordinate[1], coordinate[2]};
    external_lower[direction] -= parameters.external_metric_fd_step;
    external_upper[direction] += parameters.external_metric_fd_step;
    secondary_lower[direction] -= parameters.metric_fd_step;
    secondary_upper[direction] += parameters.metric_fd_step;
    if (!(external_upper[direction] > external_lower[direction])
        || !(secondary_upper[direction] > secondary_lower[direction])) {
      Kokkos::abort("metric finite-difference step is too small for the precision");
    }
    emri_comoving::ComputeExternalMetric(
        external_lower[0], external_lower[1], external_lower[2], parameters.metric,
        external_minus);
    emri_comoving::ComputeExternalMetric(
        external_upper[0], external_upper[1], external_upper[2], parameters.metric,
        external_plus);
    emri_comoving::ComputeSecondaryMetricPerturbation(
        secondary_lower[0], secondary_lower[1], secondary_lower[2], parameters.metric,
        secondary_minus);
    emri_comoving::ComputeSecondaryMetricPerturbation(
        secondary_upper[0], secondary_upper[1], secondary_upper[2], parameters.metric,
        secondary_plus);
    const Real inverse_external_width =
        1.0/(external_upper[direction]-external_lower[direction]);
    const Real inverse_secondary_width =
        1.0/(secondary_upper[direction]-secondary_lower[direction]);
    for (int a=0; a<4; ++a) {
      for (int b=0; b<4; ++b) {
        metric.dg[direction+1][a][b] =
            (external_plus[a][b]-external_minus[a][b])*inverse_external_width
            +(secondary_plus[a][b]-secondary_minus[a][b])*inverse_secondary_width;
      }
    }
  }
}

//! Differentiate the secondary perturbation with respect to field-point displacement
//! while holding the rotating-coordinate Jacobian fixed.  The negative of this
//! derivative is the metric derivative with respect to the secondary's source position.
KOKKOS_INLINE_FUNCTION
void DifferentiateSecondaryDisplacement(
    const Real x, const Real y, const Real z,
    const emri_comoving::MetricParameters &metric_parameters,
    const Real fd_step, Real derivative[3][4][4]) {
  const Real coordinate[3] = {x, y, z};
  for (int direction=0; direction<3; ++direction) {
    Real lower[3] = {coordinate[0], coordinate[1], coordinate[2]};
    Real upper[3] = {coordinate[0], coordinate[1], coordinate[2]};
    lower[direction] -= fd_step;
    upper[direction] += fd_step;
    Real metric_lower[4][4];
    Real metric_upper[4][4];
    emri_comoving::ComputeSecondaryMetricPerturbationAtDisplacement(
        lower[0], lower[1], lower[2], x, y, metric_parameters, metric_lower);
    emri_comoving::ComputeSecondaryMetricPerturbationAtDisplacement(
        upper[0], upper[1], upper[2], x, y, metric_parameters, metric_upper);
    const Real inverse_width = 1.0/(upper[direction]-lower[direction]);
    for (int a=0; a<4; ++a) {
      for (int b=0; b<4; ++b) {
        derivative[direction][a][b] =
            (metric_upper[a][b]-metric_lower[a][b])*inverse_width;
      }
    }
  }
}

//! Reconstruct the covariant four-metric from its ADM decomposition.
KOKKOS_INLINE_FUNCTION
void ADMToFourMetric(const Real gamma[3][3], const Real alpha,
                     const Real beta[3], Real metric[4][4]) {
  Real beta_lower[3] = {0.0, 0.0, 0.0};
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      metric[i+1][j+1] = gamma[i][j];
      beta_lower[i] += gamma[i][j]*beta[j];
    }
    metric[0][i+1] = beta_lower[i];
    metric[i+1][0] = beta_lower[i];
  }
  metric[0][0] = -SQR(alpha);
  for (int i=0; i<3; ++i) metric[0][0] += beta_lower[i]*beta[i];
}

//! Reconstruct T^{mu nu} from dynamical-GRMHD primitives and a supplied metric.  The
//! velocity primitive is W v^i measured by the Eulerian observer, while the CT magnetic
//! field is densitized by sqrt(gamma).
KOKKOS_INLINE_FUNCTION
void ComputeStressEnergyFromMetric(
    const Real metric[4][4],
    const WindTunnelParameters &parameters, const Real density, const Real pressure,
    const Real velocity[3], const Real densitized_field[3],
    Real four_velocity[4], Real stress[4][4], Real &sqrt_gamma,
    Real &sqrt_minus_g) {
  const Real determinant = adm::SpatialDet(
      metric[1][1], metric[1][2], metric[1][3], metric[2][2], metric[2][3],
      metric[3][3]);
  if (!(determinant > 0.0) || !isfinite(determinant)) {
    Kokkos::abort("invalid metric determinant in EMRI force diagnostic");
  }
  sqrt_gamma = sqrt(determinant);

  Real inverse_gamma[3][3];
  adm::SpatialInv(1.0/determinant,
      metric[1][1], metric[1][2], metric[1][3], metric[2][2], metric[2][3],
      metric[3][3], &inverse_gamma[0][0], &inverse_gamma[0][1],
      &inverse_gamma[0][2], &inverse_gamma[1][1], &inverse_gamma[1][2],
      &inverse_gamma[2][2]);
  inverse_gamma[1][0] = inverse_gamma[0][1];
  inverse_gamma[2][0] = inverse_gamma[0][2];
  inverse_gamma[2][1] = inverse_gamma[1][2];

  Real beta[3] = {0.0, 0.0, 0.0};
  Real beta_squared = 0.0;
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) beta[i] += inverse_gamma[i][j]*metric[0][j+1];
    beta_squared += metric[0][i+1]*beta[i];
  }
  const Real alpha_squared = beta_squared-metric[0][0];
  if (!(alpha_squared > 0.0) || !isfinite(alpha_squared)) {
    Kokkos::abort("invalid lapse in EMRI force diagnostic");
  }
  const Real alpha = sqrt(alpha_squared);
  sqrt_minus_g = alpha*sqrt_gamma;

  Real inverse_metric[4][4] = {{0.0}};
  inverse_metric[0][0] = -1.0/alpha_squared;
  for (int i=0; i<3; ++i) {
    inverse_metric[0][i+1] = beta[i]/alpha_squared;
    inverse_metric[i+1][0] = inverse_metric[0][i+1];
    for (int j=0; j<3; ++j) {
      inverse_metric[i+1][j+1] = inverse_gamma[i][j]
          - beta[i]*beta[j]/alpha_squared;
    }
  }

  Real velocity_lower[3] = {0.0, 0.0, 0.0};
  Real velocity_squared = 0.0;
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      velocity_lower[i] += metric[i+1][j+1]*velocity[j];
    }
    velocity_squared += velocity[i]*velocity_lower[i];
  }
  const Real lorentz = sqrt(1.0+velocity_squared);
  four_velocity[0] = lorentz/alpha;
  for (int i=0; i<3; ++i) {
    four_velocity[i+1] = velocity[i]-lorentz*beta[i]/alpha;
  }

  Real field[3];
  Real field_velocity = 0.0;
  for (int i=0; i<3; ++i) {
    field[i] = densitized_field[i]/sqrt_gamma;
    field_velocity += field[i]*velocity_lower[i];
  }
  Real field_squared = 0.0;
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      field_squared += metric[i+1][j+1]*field[i]*field[j];
    }
  }
  Real magnetic_four[4];
  magnetic_four[0] = field_velocity/alpha;
  for (int i=0; i<3; ++i) {
    magnetic_four[i+1] =
        (field[i]+alpha*magnetic_four[0]*four_velocity[i+1])/lorentz;
  }

  const Real magnetic_squared =
      (field_squared+SQR(field_velocity))/SQR(lorentz);

  const Real enthalpy_density = density
      + parameters.adiabatic_index*pressure/(parameters.adiabatic_index-1.0);
  const Real total_pressure = pressure+0.5*magnetic_squared;
  for (int a=0; a<4; ++a) {
    for (int b=0; b<4; ++b) {
      stress[a][b] = (enthalpy_density+magnetic_squared)
          *four_velocity[a]*four_velocity[b] + total_pressure*inverse_metric[a][b]
          - magnetic_four[a]*magnetic_four[b];
    }
  }
}

KOKKOS_INLINE_FUNCTION
void ComputeStressEnergy(
    const Real x, const Real y, const Real z,
    const WindTunnelParameters &parameters, const Real density, const Real pressure,
    const Real velocity[3], const Real densitized_field[3], Real metric[4][4],
    Real four_velocity[4], Real stress[4][4], Real &sqrt_gamma,
    Real &sqrt_minus_g) {
  EvaluateMetric(x, y, z, parameters, metric);
  ComputeStressEnergyFromMetric(
      metric, parameters, density, pressure, velocity, densitized_field,
      four_velocity, stress, sqrt_gamma, sqrt_minus_g);
}

// dt(gamma_ij) = -2 alpha K_ij + D_i beta_j + D_j beta_i.
KOKKOS_INLINE_FUNCTION
void DecomposeMetric(const MetricWithDerivatives &metric, ADMPoint &result) {
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) result.gamma[i][j] = metric.g[i+1][j+1];
  }

  const Real determinant = adm::SpatialDet(
      result.gamma[0][0], result.gamma[0][1], result.gamma[0][2],
      result.gamma[1][1], result.gamma[1][2], result.gamma[2][2]);
  const Real leading_minor2 = result.gamma[0][0]*result.gamma[1][1]
                            - SQR(result.gamma[0][1]);
  if (!isfinite(determinant) || !(result.gamma[0][0] > 0.0)
      || !(leading_minor2 > 0.0) || !(determinant > 0.0)) {
    Kokkos::abort("Non-positive spatial metric in the local EMRI spacetime");
  }

  Real inverse[3][3];
  adm::SpatialInv(1.0/determinant,
      result.gamma[0][0], result.gamma[0][1], result.gamma[0][2],
      result.gamma[1][1], result.gamma[1][2], result.gamma[2][2],
      &inverse[0][0], &inverse[0][1], &inverse[0][2],
      &inverse[1][1], &inverse[1][2], &inverse[2][2]);
  inverse[1][0] = inverse[0][1];
  inverse[2][0] = inverse[0][2];
  inverse[2][1] = inverse[1][2];

  const Real beta_lower[3] = {metric.g[0][1], metric.g[0][2], metric.g[0][3]};
  Real beta_squared = 0.0;
  for (int i=0; i<3; ++i) {
    result.beta[i] = 0.0;
    for (int j=0; j<3; ++j) result.beta[i] += inverse[i][j]*beta_lower[j];
    beta_squared += beta_lower[i]*result.beta[i];
  }
  const Real alpha_squared = beta_squared-metric.g[0][0];
  if (!isfinite(alpha_squared) || !(alpha_squared > 0.0)) {
    Kokkos::abort("Non-positive lapse squared in the local EMRI spacetime");
  }
  result.alpha = sqrt(alpha_squared);
  result.psi4 = cbrt(determinant);

  for (int i=0; i<3; ++i) {
    for (int j=i; j<3; ++j) {
      Real covariant_i_beta_j = metric.dg[i+1][0][j+1];
      Real covariant_j_beta_i = metric.dg[j+1][0][i+1];
      for (int k=0; k<3; ++k) {
        Real christoffel = 0.0;
        for (int l=0; l<3; ++l) {
          christoffel += 0.5*inverse[k][l]*(
              metric.dg[i+1][l+1][j+1] + metric.dg[j+1][l+1][i+1]
              - metric.dg[l+1][i+1][j+1]);
        }
        covariant_i_beta_j -= christoffel*beta_lower[k];
        covariant_j_beta_i -= christoffel*beta_lower[k];
      }
      result.curvature[i][j] = (covariant_i_beta_j+covariant_j_beta_i
          - metric.dg[0][i+1][j+1])/(2.0*result.alpha);
      result.curvature[j][i] = result.curvature[i][j];
    }
  }
}

Real KerrISCO(const Real mass, const Real chi, const Real direction) {
  const Real effective_spin = direction*chi;
  const Real z1 = 1.0 + std::cbrt(1.0-SQR(effective_spin))
      *(std::cbrt(1.0+effective_spin)+std::cbrt(1.0-effective_spin));
  const Real z2 = std::sqrt(3.0*SQR(effective_spin)+SQR(z1));
  Real signed_root = 0.0;
  if (effective_spin != 0.0) {
    signed_root = std::copysign(
        std::sqrt((3.0-z1)*(3.0+z1+2.0*z2)), effective_spin);
  }
  return mass*(3.0+z2-signed_root);
}

Real MatrixSpectralNorm(const Real matrix[3][3]) {
  Real normal_matrix[3][3] = {{0.0}};
  for (int row=0; row<3; ++row) {
    for (int column=0; column<3; ++column) {
      for (int coordinate=0; coordinate<3; ++coordinate) {
        normal_matrix[row][column] +=
            matrix[coordinate][row]*matrix[coordinate][column];
      }
    }
  }
  int seed_axis = 0;
  if (normal_matrix[1][1] > normal_matrix[seed_axis][seed_axis]) seed_axis = 1;
  if (normal_matrix[2][2] > normal_matrix[seed_axis][seed_axis]) seed_axis = 2;
  Real vector[3] = {0.0, 0.0, 0.0};
  vector[seed_axis] = 1.0;
  Real eigenvalue = 0.0;
  for (int iteration=0; iteration<32; ++iteration) {
    Real product[3] = {0.0, 0.0, 0.0};
    for (int row=0; row<3; ++row) {
      for (int column=0; column<3; ++column) {
        product[row] += normal_matrix[row][column]*vector[column];
      }
    }
    const Real norm = std::sqrt(SQR(product[0])+SQR(product[1])+SQR(product[2]));
    if (!(norm > 0.0) || !std::isfinite(norm)) {
      Fatal("failed to determine a tangent-frame matrix norm");
    }
    for (int axis=0; axis<3; ++axis) vector[axis] = product[axis]/norm;
    eigenvalue = 0.0;
    for (int row=0; row<3; ++row) {
      for (int column=0; column<3; ++column) {
        eigenvalue += vector[row]*normal_matrix[row][column]*vector[column];
      }
    }
  }
  return std::sqrt(eigenvalue);
}

void CoframeSpatialStretches(const Real coframe[4][4],
                             Real &coordinate_stretch, Real &rest_stretch) {
  const Real a = coframe[1][1];
  const Real b = coframe[1][2];
  const Real c = coframe[1][3];
  const Real d = coframe[2][1];
  const Real e = coframe[2][2];
  const Real f = coframe[2][3];
  const Real g = coframe[3][1];
  const Real h = coframe[3][2];
  const Real i = coframe[3][3];
  const Real determinant = a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g);
  if (!(std::abs(determinant) > std::numeric_limits<Real>::epsilon()) ||
      !std::isfinite(determinant)) {
    Fatal("secondary replay coframe has a singular spatial block");
  }
  const Real inverse[3][3] = {
    {(e*i-f*h)/determinant, (c*h-b*i)/determinant,
     (b*f-c*e)/determinant},
    {(f*g-d*i)/determinant, (a*i-c*g)/determinant,
     (c*d-a*f)/determinant},
    {(d*h-e*g)/determinant, (b*g-a*h)/determinant,
     (a*e-b*d)/determinant}
  };
  const Real spatial[3][3] = {{a, b, c}, {d, e, f}, {g, h, i}};
  coordinate_stretch = MatrixSpectralNorm(inverse);
  rest_stretch = MatrixSpectralNorm(spatial);
}

void ConfigureSourceFrame() {
  auto &parameters = wind_tunnel;
  BuildSourceTetrad(parameters.metric, parameters.source_tetrad);
  parameters.source_dt_dtau = parameters.source_tetrad[0][0];
  for (int source_axis=0; source_axis<3; ++source_axis) {
    for (int coordinate_axis=0; coordinate_axis<3; ++coordinate_axis) {
      parameters.force_projection[source_axis][coordinate_axis] =
          parameters.source_dt_dtau
          *parameters.source_tetrad[source_axis+1][coordinate_axis+1];
    }
  }

  Real external_metric[4][4];
  emri_comoving::ComputeExternalMetric(0.0, 0.0, 0.0, parameters.metric,
                                       external_metric);
  // The dual co-basis theta^A_mu = eta^AA g_mu_nu e_A^nu maps the local
  // coordinate displacement into the secondary's orthonormal tangent frame.
  for (int leg=0; leg<4; ++leg) {
    const Real signature = (leg == 0) ? -1.0 : 1.0;
    for (int component=0; component<4; ++component) {
      parameters.metric.secondary_tetrad_covector[leg][component] = 0.0;
      for (int contracted=0; contracted<4; ++contracted) {
        parameters.metric.secondary_tetrad_covector[leg][component] +=
            signature*external_metric[component][contracted]
            *parameters.source_tetrad[leg][contracted];
      }
    }
  }

  const Real a = parameters.metric.secondary_tetrad_covector[1][1];
  const Real b = parameters.metric.secondary_tetrad_covector[1][2];
  const Real c = parameters.metric.secondary_tetrad_covector[1][3];
  const Real d = parameters.metric.secondary_tetrad_covector[2][1];
  const Real e = parameters.metric.secondary_tetrad_covector[2][2];
  const Real f = parameters.metric.secondary_tetrad_covector[2][3];
  const Real g = parameters.metric.secondary_tetrad_covector[3][1];
  const Real h = parameters.metric.secondary_tetrad_covector[3][2];
  const Real i = parameters.metric.secondary_tetrad_covector[3][3];
  const Real coframe_determinant = a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g);
  if (!(std::abs(coframe_determinant) > std::numeric_limits<Real>::epsilon())
      || !std::isfinite(coframe_determinant)) {
    Fatal("secondary tangent-frame spatial co-basis is singular");
  }
  const Real inverse_c[3][3] = {
    {(e*i-f*h)/coframe_determinant, (c*h-b*i)/coframe_determinant,
     (b*f-c*e)/coframe_determinant},
    {(f*g-d*i)/coframe_determinant, (a*i-c*g)/coframe_determinant,
     (c*d-a*f)/coframe_determinant},
    {(d*h-e*g)/coframe_determinant, (b*g-a*h)/coframe_determinant,
     (a*e-b*d)/coframe_determinant}
  };
  const Real inverse_determinant = 1.0/std::abs(coframe_determinant);
  parameters.source_spatial_determinant = std::abs(coframe_determinant);
  for (int coordinate=0; coordinate<3; ++coordinate) {
    for (int source_axis=0; source_axis<3; ++source_axis) {
      parameters.source_spatial_inverse[coordinate][source_axis] =
          inverse_c[coordinate][source_axis];
      parameters.source_surface_covector[source_axis][coordinate] =
          inverse_determinant
          *parameters.metric.secondary_tetrad_covector[source_axis+1][coordinate+1];
    }
  }
  const Real spatial_c[3][3] = {{a, b, c}, {d, e, f}, {g, h, i}};
  parameters.source_coordinate_stretch = MatrixSpectralNorm(inverse_c);
  parameters.source_tetrad_stretch = MatrixSpectralNorm(spatial_c);

  if (parameters.metric.embed_secondary_in_tetrad) {
    // Bound both directions of the tangent-frame map so refinement and excision
    // remain conservative in a sheared or stretched local chart.
    parameters.secondary_coordinate_stretch = parameters.source_coordinate_stretch;
    parameters.secondary_rest_stretch = parameters.source_tetrad_stretch;
  } else {
    const Real speed = parameters.metric.include_orbital_frame
        ? parameters.metric.omega*parameters.metric.coordinate_radius : 0.0;
    parameters.secondary_coordinate_stretch = 1.0/std::sqrt(1.0-SQR(speed));
    parameters.secondary_rest_stretch = parameters.secondary_coordinate_stretch;
  }
#if SINGLE_PRECISION_ENABLED
  const Real tolerance = 2.0e-5;
#else
  const Real tolerance = 2.0e-12;
#endif
  for (int leg_a=0; leg_a<4; ++leg_a) {
    for (int leg_b=0; leg_b<4; ++leg_b) {
      const Real expected = (leg_a == leg_b) ? ((leg_a == 0) ? -1.0 : 1.0) : 0.0;
      const Real product = MetricInnerProduct(
          external_metric, parameters.source_tetrad[leg_a],
          parameters.source_tetrad[leg_b]);
      if (std::abs(product-expected) > tolerance) {
        Fatal("EMRI source tetrad failed its orthonormality check");
      }
    }
  }

  if (!parameters.wind_is_source_tetrad) {
    for (int i=0; i<3; ++i) {
      parameters.wind_eulerian_u[i] = parameters.wind_u[i];
      parameters.densitized_magnetic_field[i] = parameters.magnetic_field[i];
    }
    return;
  }

  Real source_lorentz2 = 1.0;
  for (int i=0; i<3; ++i) source_lorentz2 += SQR(parameters.wind_u[i]);
  const Real source_lorentz = std::sqrt(source_lorentz2);
  Real coordinate_velocity[4] = {0.0, 0.0, 0.0, 0.0};
  for (int component=0; component<4; ++component) {
    coordinate_velocity[component] =
        source_lorentz*parameters.source_tetrad[0][component];
    for (int axis=0; axis<3; ++axis) {
      coordinate_velocity[component] +=
          parameters.wind_u[axis]*parameters.source_tetrad[axis+1][component];
    }
  }
  if (!(coordinate_velocity[0] > 0.0) || !std::isfinite(coordinate_velocity[0])) {
    Fatal("source-frame wind has a non-future-directed coordinate four-velocity");
  }
  Real magnetic_source[4];
  magnetic_source[0] = 0.0;
  for (int i=0; i<3; ++i) {
    magnetic_source[0] += parameters.magnetic_field[i]*parameters.wind_u[i];
  }
  for (int i=0; i<3; ++i) {
    magnetic_source[i+1] = (parameters.magnetic_field[i]
        + magnetic_source[0]*parameters.wind_u[i])/source_lorentz;
  }
  Real magnetic_coordinate[4] = {0.0, 0.0, 0.0, 0.0};
  for (int component=0; component<4; ++component) {
    for (int leg=0; leg<4; ++leg) {
      magnetic_coordinate[component] +=
          magnetic_source[leg]*parameters.source_tetrad[leg][component];
    }
  }

  const Real determinant = adm::SpatialDet(
      external_metric[1][1], external_metric[1][2], external_metric[1][3],
      external_metric[2][2], external_metric[2][3], external_metric[3][3]);
  Real inverse_gamma[3][3];
  adm::SpatialInv(1.0/determinant,
      external_metric[1][1], external_metric[1][2], external_metric[1][3],
      external_metric[2][2], external_metric[2][3], external_metric[3][3],
      &inverse_gamma[0][0], &inverse_gamma[0][1], &inverse_gamma[0][2],
      &inverse_gamma[1][1], &inverse_gamma[1][2], &inverse_gamma[2][2]);
  inverse_gamma[1][0] = inverse_gamma[0][1];
  inverse_gamma[2][0] = inverse_gamma[0][2];
  inverse_gamma[2][1] = inverse_gamma[1][2];
  Real beta[3] = {0.0, 0.0, 0.0};
  Real beta_squared = 0.0;
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      beta[i] += inverse_gamma[i][j]*external_metric[0][j+1];
    }
    beta_squared += external_metric[0][i+1]*beta[i];
  }
  const Real alpha = std::sqrt(beta_squared-external_metric[0][0]);
  const Real eulerian_lorentz = alpha*coordinate_velocity[0];
  const Real sqrt_gamma = std::sqrt(determinant);
  Real spatial_triad[3][3];
  BuildSpatialTriad(external_metric, spatial_triad);
  Real coordinate_primitive[3];
  for (int i=0; i<3; ++i) {
    coordinate_primitive[i] = coordinate_velocity[i+1]
                            + coordinate_velocity[0]*beta[i];
  }
  for (int leg=0; leg<3; ++leg) {
    parameters.wind_eulerian_u[leg] = 0.0;
    for (int i=0; i<3; ++i) {
      for (int j=0; j<3; ++j) {
        parameters.wind_eulerian_u[leg] += external_metric[i+1][j+1]
            *coordinate_primitive[i]*spatial_triad[leg][j];
      }
    }
  }
  Real eulerian_lorentz2_check = 1.0;
  for (int leg=0; leg<3; ++leg) {
    eulerian_lorentz2_check += SQR(parameters.wind_eulerian_u[leg]);
  }
  if (std::abs(eulerian_lorentz2_check-SQR(eulerian_lorentz))
      > tolerance*SQR(eulerian_lorentz)) {
    Fatal("source-frame wind failed its Eulerian Lorentz-factor check");
  }
  for (int i=0; i<3; ++i) {
    const Real eulerian_field = eulerian_lorentz*magnetic_coordinate[i+1]
        - alpha*magnetic_coordinate[0]*coordinate_velocity[i+1];
    parameters.densitized_magnetic_field[i] = sqrt_gamma*eulerian_field;
  }

  const Real recovered_lorentz = -MetricInnerProduct(
      external_metric, coordinate_velocity, parameters.source_tetrad[0]);
  if (std::abs(recovered_lorentz-source_lorentz) > tolerance*source_lorentz) {
    Fatal("source-frame wind failed its tetrad transformation check");
  }
  for (int axis=0; axis<3; ++axis) {
    const Real recovered_velocity = MetricInnerProduct(
        external_metric, coordinate_velocity, parameters.source_tetrad[axis+1]);
    if (std::abs(recovered_velocity-parameters.wind_u[axis])
        > tolerance*source_lorentz) {
      Fatal("source-frame wind failed its spatial tetrad transformation check");
    }
    parameters.wind_eulerian_u[axis] = parameters.wind_u[axis];
  }
}

//! Convert a source-tangent magnetic-flux gradient into the numerical coordinate
//! flux-density gradient.  For xhat=C x, a contravariant vector density transforms as
//! Bdens=det(C) C^{-1} Bhat.  Consequently G=det(C) C^{-1} Ghat C and its trace remains
//! zero.  The constant field keeps the exact legacy source-to-slicing transformation.
void ConfigureMagneticGradient() {
  auto &parameters = wind_tunnel;
  const Real determinant = parameters.source_spatial_determinant;
  for (int coordinate_component=0; coordinate_component<3;
       ++coordinate_component) {
    for (int coordinate_direction=0; coordinate_direction<3;
         ++coordinate_direction) {
      parameters.densitized_magnetic_gradient[coordinate_component]
                                                [coordinate_direction] = 0.0;
      for (int source_component=0; source_component<3; ++source_component) {
        for (int source_direction=0; source_direction<3; ++source_direction) {
          parameters.densitized_magnetic_gradient[coordinate_component]
                                                  [coordinate_direction] +=
              determinant
              *parameters.source_spatial_inverse[coordinate_component]
                                                   [source_component]
              *parameters.magnetic_gradient_source[source_component]
                                                    [source_direction]
              *parameters.metric.secondary_tetrad_covector[source_direction+1]
                                                            [coordinate_direction+1];
        }
      }
    }
  }

  Real trace = 0.0;
  Real scale = 0.0;
  for (int component=0; component<3; ++component) {
    trace += parameters.densitized_magnetic_gradient[component][component];
    for (int direction=0; direction<3; ++direction) {
      scale = std::max(scale, std::abs(
          parameters.densitized_magnetic_gradient[component][direction]));
    }
  }
  const Real tolerance = 256.0*std::numeric_limits<Real>::epsilon()
                       *std::max(scale, Real(1.0));
  if (std::abs(trace) > tolerance) {
    Fatal("source-frame trace-free magnetic gradient lost zero divergence");
  }
}

void ValidateSpatialProfiles(const RegionSize &domain) {
  const auto &parameters = wind_tunnel;
  Real maximum_density_exponent = 0.0;
  Real maximum_pressure_exponent = 0.0;
  for (int corner=0; corner<8; ++corner) {
    const Real x = (corner & 1) ? domain.x1max : domain.x1min;
    const Real y = (corner & 2) ? domain.x2max : domain.x2min;
    const Real z = (corner & 4) ? domain.x3max : domain.x3min;
    Real source_position[3];
    ComputeSourceDisplacement(x, y, z, parameters, source_position);
    Real density_exponent = 0.0;
    Real pressure_exponent = 0.0;
    Real local_wind[3] = {
      parameters.wind_u[0], parameters.wind_u[1], parameters.wind_u[2]
    };
    for (int direction=0; direction<3; ++direction) {
      density_exponent +=
          parameters.log_density_gradient[direction]*source_position[direction];
      pressure_exponent +=
          parameters.log_pressure_gradient[direction]*source_position[direction];
      for (int component=0; component<3; ++component) {
        local_wind[component] +=
            parameters.velocity_gradient[component][direction]
            *source_position[direction];
      }
    }
    maximum_density_exponent = std::max(
        maximum_density_exponent, std::abs(density_exponent));
    maximum_pressure_exponent = std::max(
        maximum_pressure_exponent, std::abs(pressure_exponent));
    for (const Real component : local_wind) {
      if (!std::isfinite(component)) {
        Fatal("velocity-gradient profile is non-finite at a domain corner");
      }
    }
    Real field[3];
    ComputeDensitizedMagneticField(x, y, z, parameters, field);
    for (const Real component : field) {
      if (!std::isfinite(component)) {
        Fatal("magnetic-gradient profile is non-finite at a domain corner");
      }
    }
  }
  if (maximum_density_exponent > parameters.max_log_contrast
      || maximum_pressure_exponent > parameters.max_log_contrast) {
    Fatal("analytic density/pressure profile exceeds max_log_contrast at a domain "
          "corner; reduce the gradient or increase the limit deliberately");
  }
}

std::string FormatFNV1a64(const std::uint64_t hash) {
  std::ostringstream formatted;
  formatted << std::hex << std::setfill('0') << std::setw(16) << hash;
  return formatted.str();
}

std::string HashStringFNV1a64(const std::string &value) {
  std::uint64_t hash = UINT64_C(14695981039346656037);
  for (const unsigned char byte : value) {
    hash ^= byte;
    hash *= UINT64_C(1099511628211);
  }
  return FormatFNV1a64(hash);
}

std::uint64_t HashWindProfileFile(const std::string &path) {
  std::ifstream stream(path, std::ios::binary);
  if (!stream) Fatal("cannot open <problem> profile_file: "+path);
  std::uint64_t hash = UINT64_C(14695981039346656037);
  char buffer[8192];
  while (stream) {
    stream.read(buffer, sizeof(buffer));
    const std::streamsize count = stream.gcount();
    for (std::streamsize index=0; index<count; ++index) {
      hash ^= static_cast<unsigned char>(buffer[index]);
      hash *= UINT64_C(1099511628211);
    }
  }
  if (!stream.eof()) Fatal("failed while hashing <problem> profile_file: "+path);
  return hash;
}

void ValidateWindProfileSample(const WindProfileSample &sample,
                               const int line_number) {
  std::vector<Real> values = {
    sample.time, sample.rho, sample.pressure,
    sample.wind_u[0], sample.wind_u[1], sample.wind_u[2],
    sample.magnetic_field[0], sample.magnetic_field[1], sample.magnetic_field[2]
  };
  for (int direction=0; direction<3; ++direction) {
    values.push_back(sample.log_density_gradient[direction]);
    values.push_back(sample.log_pressure_gradient[direction]);
    for (int component=0; component<3; ++component) {
      values.push_back(sample.velocity_gradient[component][direction]);
      values.push_back(sample.magnetic_gradient_source[component][direction]);
    }
  }
  for (const Real value : values) {
    if (!std::isfinite(value)) {
      Fatal("non-finite value in EMRI profile table line "
            +std::to_string(line_number));
    }
  }
  if (!(sample.rho > 0.0) || !(sample.pressure > 0.0)) {
    Fatal("non-positive density or pressure in EMRI profile table line "
          +std::to_string(line_number));
  }
  Real trace = 0.0;
  Real scale = 0.0;
  for (int component=0; component<3; ++component) {
    trace += sample.magnetic_gradient_source[component][component];
    for (int direction=0; direction<3; ++direction) {
      scale = std::max(scale, std::abs(
          sample.magnetic_gradient_source[component][direction]));
    }
  }
  const Real tolerance = 256.0*std::numeric_limits<Real>::epsilon()
                       *std::max(scale, Real(1.0));
  if (std::abs(trace) > tolerance) {
    Fatal("magnetic gradient is not trace-free in EMRI profile table line "
          +std::to_string(line_number));
  }
}

bool ParseWindProfileMetadata(const std::string &line, const std::string &name,
                              Real &value) {
  const std::string prefix = "# "+name+":";
  if (line.compare(0, prefix.size(), prefix) != 0) return false;
  if (std::isfinite(value)) Fatal("duplicate EMRI profile metadata for "+name);
  std::istringstream stream(line.substr(prefix.size()));
  stream >> value >> std::ws;
  if (!stream.eof() || !std::isfinite(value)) {
    Fatal("invalid EMRI profile metadata value for "+name);
  }
  return true;
}

void LoadWindProfileTable(const std::string &path) {
  std::ifstream stream(path);
  if (!stream) Fatal("cannot open <problem> profile_file: "+path);
  wind_profile_table.clear();
  wind_profile_source_coordinate_radius = std::numeric_limits<Real>::quiet_NaN();
  wind_profile_source_omega = std::numeric_limits<Real>::quiet_NaN();
  wind_profile_orbit_tolerance = std::numeric_limits<Real>::quiet_NaN();
  std::string line;
  int line_number = 0;
  bool has_table_contract = false;
  while (std::getline(stream, line)) {
    ++line_number;
    const std::size_t first = line.find_first_not_of(" \t\r");
    if (first == std::string::npos) continue;
    if (line[first] == '#') {
      const std::string metadata = line.substr(first);
      has_table_contract = has_table_contract
          || metadata == "# athenak-emri-taylor-series-v2";
      ParseWindProfileMetadata(
          metadata, "source_coordinate_radius_local_units",
          wind_profile_source_coordinate_radius);
      ParseWindProfileMetadata(
          metadata, "source_coordinate_angular_frequency_local_units",
          wind_profile_source_omega);
      ParseWindProfileMetadata(
          metadata, "orbit_tolerance", wind_profile_orbit_tolerance);
      continue;
    }
    std::istringstream row(line);
    std::vector<Real> values;
    Real value;
    while (row >> value) values.push_back(value);
    row >> std::ws;
    if (!row.eof() || values.size() != 33) {
      Fatal("EMRI profile table line "+std::to_string(line_number)
            +" must contain time plus 32 profile values");
    }
    WindProfileSample sample{};
    int column = 0;
    sample.time = values[column++];
    sample.rho = values[column++];
    sample.pressure = values[column++];
    for (int component=0; component<3; ++component) {
      sample.wind_u[component] = values[column++];
    }
    for (int component=0; component<3; ++component) {
      sample.magnetic_field[component] = values[column++];
    }
    for (int direction=0; direction<3; ++direction) {
      sample.log_density_gradient[direction] = values[column++];
    }
    for (int direction=0; direction<3; ++direction) {
      sample.log_pressure_gradient[direction] = values[column++];
    }
    for (int component=0; component<3; ++component) {
      for (int direction=0; direction<3; ++direction) {
        sample.velocity_gradient[component][direction] = values[column++];
      }
    }
    for (int component=0; component<3; ++component) {
      for (int direction=0; direction<3; ++direction) {
        sample.magnetic_gradient_source[component][direction] = values[column++];
      }
    }
    ValidateWindProfileSample(sample, line_number);
    if (!wind_profile_table.empty()
        && !(sample.time > wind_profile_table.back().time)) {
      Fatal("EMRI profile table times must increase strictly");
    }
    wind_profile_table.push_back(sample);
  }
  if (!stream.eof()) Fatal("failed while reading <problem> profile_file: "+path);
  if (!has_table_contract) {
    Fatal("EMRI profile table is missing the athenak-emri-taylor-series-v2 contract");
  }
  if (!(wind_profile_source_coordinate_radius > 0.0)
      || !std::isfinite(wind_profile_source_omega)
      || !(wind_profile_orbit_tolerance > 0.0)) {
    Fatal("EMRI profile table is missing valid orbit metadata");
  }
  if (wind_profile_table.size() < 2) {
    Fatal("EMRI profile table requires at least two data rows");
  }
  wind_profile_hash = HashWindProfileFile(path);
  wind_profile_current_table_time = std::numeric_limits<Real>::quiet_NaN();
}

Real InterpolateProfileValue(const Real lower, const Real upper,
                             const Real fraction) {
  return lower+fraction*(upper-lower);
}

void UpdateWindProfileForTime(const Real simulation_time,
                              const RegionSize &domain) {
  if (!wind_profile_enabled) return;
  Real table_time = simulation_time+wind_profile_time_offset;
  const Real first_time = wind_profile_table.front().time;
  const Real last_time = wind_profile_table.back().time;
  const Real time_scale = std::max({Real(1.0), std::abs(table_time),
                                    std::abs(first_time), std::abs(last_time)});
  const Real tolerance = 64.0*std::numeric_limits<Real>::epsilon()*time_scale;
  if (table_time < first_time-tolerance || table_time > last_time+tolerance) {
    if (!wind_profile_hold_outside) {
      Fatal("EMRI profile request lies outside the table time range: table_time="
            +std::to_string(table_time));
    }
    table_time = std::min(std::max(table_time, first_time), last_time);
  } else {
    table_time = std::min(std::max(table_time, first_time), last_time);
  }
  if (table_time == wind_profile_current_table_time) return;

  const auto upper = std::upper_bound(
      wind_profile_table.begin(), wind_profile_table.end(), table_time,
      [](const Real time, const WindProfileSample &sample) {
        return time < sample.time;
      });
  const WindProfileSample *left = nullptr;
  const WindProfileSample *right = nullptr;
  Real fraction = 0.0;
  if (upper == wind_profile_table.begin()) {
    left = right = &wind_profile_table.front();
  } else if (upper == wind_profile_table.end()) {
    left = right = &wind_profile_table.back();
  } else {
    right = &(*upper);
    left = &(*(upper-1));
    fraction = (table_time-left->time)/(right->time-left->time);
  }
  wind_tunnel.rho = std::exp(InterpolateProfileValue(
      std::log(left->rho), std::log(right->rho), fraction));
  wind_tunnel.pressure = std::exp(InterpolateProfileValue(
      std::log(left->pressure), std::log(right->pressure), fraction));
  for (int component=0; component<3; ++component) {
    wind_tunnel.wind_u[component] = InterpolateProfileValue(
        left->wind_u[component], right->wind_u[component], fraction);
    wind_tunnel.magnetic_field[component] = InterpolateProfileValue(
        left->magnetic_field[component], right->magnetic_field[component], fraction);
  }
  for (int direction=0; direction<3; ++direction) {
    wind_tunnel.log_density_gradient[direction] = InterpolateProfileValue(
        left->log_density_gradient[direction],
        right->log_density_gradient[direction], fraction);
    wind_tunnel.log_pressure_gradient[direction] = InterpolateProfileValue(
        left->log_pressure_gradient[direction],
        right->log_pressure_gradient[direction], fraction);
    for (int component=0; component<3; ++component) {
      wind_tunnel.velocity_gradient[component][direction] = InterpolateProfileValue(
          left->velocity_gradient[component][direction],
          right->velocity_gradient[component][direction], fraction);
      wind_tunnel.magnetic_gradient_source[component][direction] =
          InterpolateProfileValue(
              left->magnetic_gradient_source[component][direction],
              right->magnetic_gradient_source[component][direction], fraction);
    }
  }
  ConfigureSourceFrame();
  ConfigureMagneticGradient();
  ValidateSpatialProfiles(domain);
  wind_profile_current_table_time = table_time;
}

std::string MetricRestartContractPayload() {
  const auto &metric = wind_tunnel.metric;
  std::ostringstream contract;
  contract << std::setprecision(std::numeric_limits<Real>::max_digits10)
           << "emri-comoving-v5"
           << ";primary_mass=" << metric.primary_mass
           << ";secondary_mass=" << metric.secondary_mass
           << ";primary_spin=" << metric.primary_spin
           << ";secondary_spin=" << metric.secondary_spin
           << ";orbital_radius=" << metric.orbital_radius
           << ";coordinate_radius=" << metric.coordinate_radius
           << ";omega=" << metric.omega
           << ";spin_buffer_primary=" << metric.spin_buffer_primary
           << ";spin_buffer_secondary=" << metric.spin_buffer_secondary
           << ";singularity_floor=" << metric.singularity_floor
           << ";include_primary=" << metric.include_primary
           << ";include_orbital_frame=" << metric.include_orbital_frame
           << ";embed_secondary_in_tetrad=" << metric.embed_secondary_in_tetrad
           << ";metric_fd_step=" << wind_tunnel.metric_fd_step
           << ";external_metric_fd_step=" << wind_tunnel.external_metric_fd_step
           << ";max_log_contrast=" << wind_tunnel.max_log_contrast;
  if (wind_profile_enabled) {
    contract << ";profile_hash_fnv1a64=" << FormatFNV1a64(wind_profile_hash)
             << ";profile_rows=" << wind_profile_table.size()
             << ";profile_first_time=" << wind_profile_table.front().time
             << ";profile_last_time=" << wind_profile_table.back().time
             << ";profile_time_offset=" << wind_profile_time_offset
             << ";profile_hold_outside=" << wind_profile_hold_outside;
  } else {
    contract << ";rho=" << wind_tunnel.rho
             << ";pressure=" << wind_tunnel.pressure
             << ";wind_u1=" << wind_tunnel.wind_u[0]
             << ";wind_u2=" << wind_tunnel.wind_u[1]
             << ";wind_u3=" << wind_tunnel.wind_u[2]
             << ";magnetic_b1=" << wind_tunnel.magnetic_field[0]
             << ";magnetic_b2=" << wind_tunnel.magnetic_field[1]
             << ";magnetic_b3=" << wind_tunnel.magnetic_field[2];
    for (int direction=0; direction<3; ++direction) {
      contract << ";dlnrho_dxh" << direction+1 << "="
               << wind_tunnel.log_density_gradient[direction]
               << ";dlnpgas_dxh" << direction+1 << "="
               << wind_tunnel.log_pressure_gradient[direction];
    }
    for (int direction=0; direction<3; ++direction) {
      for (int component=0; component<3; ++component) {
        contract << ";du" << component+1 << "_dxh" << direction+1 << "="
                 << wind_tunnel.velocity_gradient[component][direction]
                 << ";db" << component+1 << "_dxh" << direction+1 << "="
                 << wind_tunnel.magnetic_gradient_source[component][direction];
      }
    }
  }
  contract
           << ";wind_is_source_tetrad=" << wind_tunnel.wind_is_source_tetrad
           << ";force_is_source_tetrad=" << wind_tunnel.force_is_source_tetrad;
  return contract.str();
}

std::string MetricRestartContract() {
  return "emri-comoving-v5-fnv1a64="
         +HashStringFNV1a64(MetricRestartContractPayload());
}

void ValidateOrStoreRestartContract(ParameterInput *pin, const bool restart) {
  constexpr const char *name = "emri_metric_restart_contract";
  const std::string current = MetricRestartContract();
  if (restart) {
    if (!pin->DoesParameterExist("problem", name)
        || pin->GetString("problem", name) != current) {
      Fatal("EMRI metric/profile parameters differ from the restart contract");
    }
  } else {
    if (pin->DoesParameterExist("problem", name)) {
      Fatal(std::string("<problem> ")+name
            + " is reserved for restart provenance and must not be supplied");
    }
    pin->SetString("problem", name, current);
  }
}

void ValidateMetricKernel() {
#if SINGLE_PRECISION_ENABLED
  const Real tolerance = 5.0e-5;
#else
  const Real tolerance = 3.0e-12;
#endif
  emri_comoving::MetricParameters isolated{};
  isolated.secondary_mass = 1.0;
  isolated.secondary_spin = 0.37;
  isolated.spin_buffer_primary = 0.05;
  isolated.spin_buffer_secondary = 0.05;
  isolated.singularity_floor = 1.0e-3;
  Real metric[4][4];
  emri_comoving::ComputeMetric(2.1, -0.7, 0.9, isolated, metric);
  Real reference[4][4];
  Real inverse[4][4];
  ComputeMetricAndInverse(2.1, -0.7, 0.9, false, 0.37, reference, inverse);
  for (int a=0; a<4; ++a) {
    for (int b=0; b<4; ++b) {
      if (std::abs(metric[a][b]-reference[a][b]) > tolerance) {
        Fatal("local EMRI metric failed its isolated Kerr limit");
      }
    }
  }
  isolated.embed_secondary_in_tetrad = true;
  for (int leg=0; leg<4; ++leg) {
    for (int component=0; component<4; ++component) {
      isolated.secondary_tetrad_covector[leg][component] =
          (leg == component) ? 1.0 : 0.0;
    }
  }
  emri_comoving::ComputeMetric(2.1, -0.7, 0.9, isolated, metric);
  for (int a=0; a<4; ++a) {
    for (int b=0; b<4; ++b) {
      if (std::abs(metric[a][b]-reference[a][b]) > tolerance) {
        Fatal("tangent-embedded EMRI metric failed its isolated Kerr limit");
      }
    }
  }

  WindTunnelParameters weak_field{};
  weak_field.metric.secondary_mass = 1.0;
  weak_field.metric.spin_buffer_primary = 0.05;
  weak_field.metric.spin_buffer_secondary = 0.05;
  weak_field.metric.singularity_floor = 1.0e-3;
#if SINGLE_PRECISION_ENABLED
  weak_field.metric_fd_step = 1.0e-2;
  const Real force_kernel_tolerance = 2.0e-4;
#else
  weak_field.metric_fd_step = 1.0e-3;
  const Real force_kernel_tolerance = 2.0e-8;
#endif
  Real displacement_derivative[3][4][4];
  DifferentiateSecondaryDisplacement(
      100.0, 0.0, 0.0, weak_field.metric, weak_field.metric_fd_step,
      displacement_derivative);
  const Real newtonian_kernel = 1.0e-4;
  const Real relativistic_kernel = -0.5*displacement_derivative[0][0][0];
  if (std::abs(relativistic_kernel/newtonian_kernel-1.0)
      > force_kernel_tolerance) {
    Fatal("EMRI source-position force kernel failed its Newtonian limit");
  }

  WindTunnelParameters flat_fluid{};
  flat_fluid.metric.spin_buffer_primary = 0.05;
  flat_fluid.metric.spin_buffer_secondary = 0.05;
  flat_fluid.metric.singularity_floor = 1.0e-3;
  flat_fluid.adiabatic_index = 4.0/3.0;
  const Real rest_velocity[3] = {0.0, 0.0, 0.0};
  const Real test_field[3] = {0.4, 0.1, -0.2};
  Real test_metric[4][4];
  Real test_velocity[4];
  Real test_stress[4][4];
  Real test_sqrt_gamma;
  Real test_sqrt_minus_g;
  ComputeStressEnergy(0.2, -0.4, 0.1, flat_fluid, 2.0, 0.3,
                      rest_velocity, test_field, test_metric, test_velocity,
                      test_stress, test_sqrt_gamma, test_sqrt_minus_g);
  const Real test_magnetic2 = SQR(test_field[0])+SQR(test_field[1])
                            +SQR(test_field[2]);
  const Real test_energy = 2.0+0.3/(flat_fluid.adiabatic_index-1.0)
                         +0.5*test_magnetic2;
  if (std::abs(test_stress[0][0]-test_energy) > tolerance
      || std::abs(test_stress[1][1]
                  -(0.3+0.5*test_magnetic2-SQR(test_field[0]))) > tolerance
      || std::abs(test_stress[1][2]+test_field[0]*test_field[1]) > tolerance
      || std::abs(test_velocity[0]-1.0) > tolerance
      || std::abs(test_sqrt_gamma-1.0) > tolerance
      || std::abs(test_sqrt_minus_g-1.0) > tolerance) {
    Fatal("EMRI force diagnostic failed its flat-space stress-energy check");
  }

  emri_comoving::MetricParameters rotating{};
  rotating.coordinate_radius = 3.0;
  rotating.omega = 0.1;
  rotating.include_orbital_frame = true;
  rotating.spin_buffer_primary = 0.05;
  rotating.spin_buffer_secondary = 0.05;
  rotating.singularity_floor = 1.0e-3;
  const Real x = 0.2;
  const Real y = -0.4;
  emri_comoving::ComputeMetric(x, y, 0.1, rotating, metric);
  const Real vx = -rotating.omega*y;
  const Real vy = rotating.omega*(rotating.coordinate_radius+x);
  if (std::abs(metric[0][0]-(-1.0+SQR(vx)+SQR(vy))) > tolerance
      || std::abs(metric[0][1]-vx) > tolerance
      || std::abs(metric[0][2]-vy) > tolerance
      || std::abs(metric[1][1]-1.0) > tolerance
      || std::abs(metric[2][2]-1.0) > tolerance
      || std::abs(metric[3][3]-1.0) > tolerance) {
    Fatal("local EMRI metric failed its rotating-Minkowski limit");
  }

  const Real primary_position[3] = {
    wind_tunnel.metric.coordinate_radius, 0.0, 0.0
  };
  const Real primary_spin[3] = {0.0, 0.0, wind_tunnel.metric.primary_spin};
  const Real recovered_radius2 =
      binary_bh::KerrRadiusSquared(primary_position, primary_spin);
  const Real scale = std::max(SQR(wind_tunnel.metric.orbital_radius), Real(1.0));
  if (std::abs(recovered_radius2-SQR(wind_tunnel.metric.orbital_radius))
      > 64.0*std::numeric_limits<Real>::epsilon()*scale) {
    Fatal("orbital_radius is inconsistent with the primary Kerr-Schild radius");
  }

  MetricWithDerivatives differentiated;
  ADMPoint decomposed;
  DifferentiateMetric(2.1*wind_tunnel.metric.secondary_mass,
                      -0.7*wind_tunnel.metric.secondary_mass,
                      0.9*wind_tunnel.metric.secondary_mass,
                      wind_tunnel, differentiated);
  DecomposeMetric(differentiated, decomposed);
  if (!std::isfinite(decomposed.alpha) || !std::isfinite(decomposed.psi4)) {
    Fatal("local EMRI ADM decomposition is not finite");
  }
  Real reconstructed_metric[4][4];
  ADMToFourMetric(
      decomposed.gamma, decomposed.alpha, decomposed.beta,
      reconstructed_metric);
  for (int a=0; a<4; ++a) {
    for (int b=0; b<4; ++b) {
      if (std::abs(reconstructed_metric[a][b]-differentiated.g[a][b])
          > 10.0*tolerance) {
        Fatal("local EMRI ADM-to-four-metric reconstruction failed");
      }
    }
  }

  // In the primary-only spacetime, the local origin must be a geodesic when Omega is
  // the circular Kerr value.  This catches a missing rotational/translation term in the
  // pullback, which a mere metric-symmetry test would not detect.
  WindTunnelParameters primary_only = wind_tunnel;
  primary_only.metric.secondary_mass = 0.0;
  primary_only.metric.secondary_spin = 0.0;
  DifferentiateMetric(0.0, 0.0, 0.0, primary_only, differentiated);
  DecomposeMetric(differentiated, decomposed);
  const Real determinant = adm::SpatialDet(
      decomposed.gamma[0][0], decomposed.gamma[0][1], decomposed.gamma[0][2],
      decomposed.gamma[1][1], decomposed.gamma[1][2], decomposed.gamma[2][2]);
  Real inverse_gamma[3][3];
  adm::SpatialInv(1.0/determinant,
      decomposed.gamma[0][0], decomposed.gamma[0][1], decomposed.gamma[0][2],
      decomposed.gamma[1][1], decomposed.gamma[1][2], decomposed.gamma[2][2],
      &inverse_gamma[0][0], &inverse_gamma[0][1], &inverse_gamma[0][2],
      &inverse_gamma[1][1], &inverse_gamma[1][2], &inverse_gamma[2][2]);
  inverse_gamma[1][0] = inverse_gamma[0][1];
  inverse_gamma[2][0] = inverse_gamma[0][2];
  inverse_gamma[2][1] = inverse_gamma[1][2];
  Real gamma_i00[3] = {0.0, 0.0, 0.0};
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      const Real inverse_spacetime_ij = inverse_gamma[i][j]
          - decomposed.beta[i]*decomposed.beta[j]/SQR(decomposed.alpha);
      gamma_i00[i] -= 0.5*inverse_spacetime_ij*differentiated.dg[j+1][0][0];
    }
  }
  Real residual2 = 0.0;
  for (int i=0; i<3; ++i) {
    for (int j=0; j<3; ++j) {
      residual2 += decomposed.gamma[i][j]*gamma_i00[i]*gamma_i00[j];
    }
  }
  primary_geodesic_residual =
      wind_tunnel.metric.primary_mass*std::sqrt(std::max(residual2, Real(0.0)));
#if SINGLE_PRECISION_ENABLED
  const Real geodesic_tolerance = 3.0e-3;
#else
  const Real geodesic_tolerance = 3.0e-7;
#endif
  if (omega_is_geodesic && wind_tunnel.metric.include_primary
      && wind_tunnel.metric.include_orbital_frame
      && primary_geodesic_residual > geodesic_tolerance) {
    Fatal("the secondary-centered origin failed the circular-Kerr geodesic check");
  }

  // Manufactured analytic-profile checks.  A quadratic vector potential must recover
  // the requested linear field exactly up to floating-point differencing, and rotating
  // x, B, and G together must commute with field evaluation.
  WindTunnelParameters analytic{};
  analytic.densitized_magnetic_field[0] = 0.4;
  analytic.densitized_magnetic_field[1] = -0.2;
  analytic.densitized_magnetic_field[2] = 0.1;
  analytic.densitized_magnetic_gradient[0][0] = 0.03;
  analytic.densitized_magnetic_gradient[0][1] = -0.02;
  analytic.densitized_magnetic_gradient[0][2] = 0.01;
  analytic.densitized_magnetic_gradient[1][0] = 0.04;
  analytic.densitized_magnetic_gradient[1][1] = -0.01;
  analytic.densitized_magnetic_gradient[1][2] = -0.03;
  analytic.densitized_magnetic_gradient[2][0] = 0.02;
  analytic.densitized_magnetic_gradient[2][1] = 0.05;
  analytic.densitized_magnetic_gradient[2][2] = -0.02;
  const Real profile_position[3] = {0.7, -0.4, 0.2};
  Real expected_field[3];
  ComputeDensitizedMagneticField(profile_position[0], profile_position[1],
                                 profile_position[2], analytic, expected_field);
#if SINGLE_PRECISION_ENABLED
  const Real curl_step = 2.0e-2;
  const Real profile_tolerance = 2.0e-5;
#else
  const Real curl_step = 1.0e-3;
  const Real profile_tolerance = 2.0e-11;
#endif
  Real derivative[3][3];
  for (int direction=0; direction<3; ++direction) {
    Real lower_position[3] = {
      profile_position[0], profile_position[1], profile_position[2]
    };
    Real upper_position[3] = {
      profile_position[0], profile_position[1], profile_position[2]
    };
    lower_position[direction] -= curl_step;
    upper_position[direction] += curl_step;
    Real lower_potential[3];
    Real upper_potential[3];
    ComputeMagneticVectorPotential(lower_position[0], lower_position[1],
                                   lower_position[2], analytic, lower_potential);
    ComputeMagneticVectorPotential(upper_position[0], upper_position[1],
                                   upper_position[2], analytic, upper_potential);
    for (int component=0; component<3; ++component) {
      derivative[component][direction] =
          (upper_potential[component]-lower_potential[component])/(2.0*curl_step);
    }
  }
  const Real recovered_field[3] = {
    derivative[2][1]-derivative[1][2],
    derivative[0][2]-derivative[2][0],
    derivative[1][0]-derivative[0][1]
  };
  for (int component=0; component<3; ++component) {
    if (std::abs(recovered_field[component]-expected_field[component])
        > profile_tolerance) {
      Fatal("analytic EMRI magnetic vector potential failed its curl check");
    }
  }

  WindTunnelParameters rotated{};
  const Real rotation[3][3] = {
    {0.0, -1.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}
  };
  Real rotated_position[3] = {0.0, 0.0, 0.0};
  for (int row=0; row<3; ++row) {
    for (int column=0; column<3; ++column) {
      rotated_position[row] += rotation[row][column]*profile_position[column];
      rotated.densitized_magnetic_field[row] +=
          rotation[row][column]*analytic.densitized_magnetic_field[column];
    }
  }
  for (int row=0; row<3; ++row) {
    for (int column=0; column<3; ++column) {
      for (int left=0; left<3; ++left) {
        for (int right=0; right<3; ++right) {
          rotated.densitized_magnetic_gradient[row][column] +=
              rotation[row][left]
              *analytic.densitized_magnetic_gradient[left][right]
              *rotation[column][right];
        }
      }
    }
  }
  Real rotated_field[3];
  ComputeDensitizedMagneticField(rotated_position[0], rotated_position[1],
                                 rotated_position[2], rotated, rotated_field);
  for (int row=0; row<3; ++row) {
    Real rotated_expected = 0.0;
    for (int column=0; column<3; ++column) {
      rotated_expected += rotation[row][column]*expected_field[column];
    }
    if (std::abs(rotated_field[row]-rotated_expected) > profile_tolerance) {
      Fatal("analytic EMRI magnetic profile failed its rotation-covariance check");
    }
  }
}

void InitializeInflowArrays(MeshBlockPack *pmbp) {
  auto &u_in = pmbp->pmhd->pbval_u->u_in;
  auto &b_in = pmbp->pmhd->pbval_b->b_in;
  const int nvariables = pmbp->pmhd->nmhd+pmbp->pmhd->nscalars;
  for (int face=0; face<6; ++face) {
    for (int n=0; n<nvariables; ++n) u_in.h_view(n, face) = 0.0;
    u_in.h_view(IDN, face) = wind_tunnel.rho;
    u_in.h_view(IVX, face) = wind_tunnel.wind_u[0];
    u_in.h_view(IVY, face) = wind_tunnel.wind_u[1];
    u_in.h_view(IVZ, face) = wind_tunnel.wind_u[2];
    u_in.h_view(IPR, face) = wind_tunnel.pressure;
    b_in.h_view(IBX, face) = wind_tunnel.densitized_magnetic_field[0];
    b_in.h_view(IBY, face) = wind_tunnel.densitized_magnetic_field[1];
    b_in.h_view(IBZ, face) = wind_tunnel.densitized_magnetic_field[2];
  }
  u_in.template modify<HostMemSpace>();
  u_in.template sync<DevExeSpace>();
  b_in.template modify<HostMemSpace>();
  b_in.template sync<DevExeSpace>();
}

//! Initialize CT face fluxes from a vector potential evaluated on grid edges.  The
//! fine-neighbor correction makes the edge integral shared by coarse and fine faces
//! identical, following the standard AthenaK vector-potential initialization pattern.
void InitializeAnalyticMagneticField(MeshBlockPack *pmbp) {
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int ncells1 = indcs.nx1+2*indcs.ng;
  const int ncells2 = (indcs.nx2 > 1) ? indcs.nx2+2*indcs.ng : 2;
  const int ncells3 = (indcs.nx3 > 1) ? indcs.nx3+2*indcs.ng : 2;
  const int nmb = pmbp->nmb_thispack;
  DvceArray4D<Real> a1;
  DvceArray4D<Real> a2;
  DvceArray4D<Real> a3;
  Kokkos::realloc(a1, nmb, ncells3, ncells2, ncells1);
  Kokkos::realloc(a2, nmb, ncells3, ncells2, ncells1);
  Kokkos::realloc(a3, nmb, ncells3, ncells2, ncells1);

  auto &size = pmbp->pmb->mb_size;
  auto &nghbr = pmbp->pmb->nghbr;
  auto &mblev = pmbp->pmb->mb_lev;
  const WindTunnelParameters parameters = wind_tunnel;
  par_for("emri_vector_potential", DevExeSpace(), 0, nmb-1, ks, ke+1,
          js, je+1, is, ie+1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x1v = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min,
                                 size.d_view(m).x1max);
    const Real x1f = LeftEdgeX(i-is, indcs.nx1, size.d_view(m).x1min,
                               size.d_view(m).x1max);
    const Real x2v = CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min,
                                 size.d_view(m).x2max);
    const Real x2f = LeftEdgeX(j-js, indcs.nx2, size.d_view(m).x2min,
                               size.d_view(m).x2max);
    const Real x3v = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min,
                                 size.d_view(m).x3max);
    const Real x3f = LeftEdgeX(k-ks, indcs.nx3, size.d_view(m).x3min,
                               size.d_view(m).x3max);
    const Real dx1 = size.d_view(m).dx1;
    const Real dx2 = size.d_view(m).dx2;
    const Real dx3 = size.d_view(m).dx3;
    Real potential[3];
    ComputeMagneticVectorPotential(x1v, x2f, x3f, parameters, potential);
    a1(m, k, j, i) = potential[0];
    ComputeMagneticVectorPotential(x1f, x2v, x3f, parameters, potential);
    a2(m, k, j, i) = potential[1];
    ComputeMagneticVectorPotential(x1f, x2f, x3v, parameters, potential);
    a3(m, k, j, i) = potential[2];

    if (EdgeTouchesFinerNeighbor(nghbr.d_view, m, mblev.d_view(m), 1,
        i, j, k, is, ie, js, je, ks, ke)) {
      Real lower[3];
      Real upper[3];
      ComputeMagneticVectorPotential(x1v-0.25*dx1, x2f, x3f,
                                     parameters, lower);
      ComputeMagneticVectorPotential(x1v+0.25*dx1, x2f, x3f,
                                     parameters, upper);
      a1(m, k, j, i) = 0.5*(lower[0]+upper[0]);
    }
    if (EdgeTouchesFinerNeighbor(nghbr.d_view, m, mblev.d_view(m), 2,
        i, j, k, is, ie, js, je, ks, ke)) {
      Real lower[3];
      Real upper[3];
      ComputeMagneticVectorPotential(x1f, x2v-0.25*dx2, x3f,
                                     parameters, lower);
      ComputeMagneticVectorPotential(x1f, x2v+0.25*dx2, x3f,
                                     parameters, upper);
      a2(m, k, j, i) = 0.5*(lower[1]+upper[1]);
    }
    if (EdgeTouchesFinerNeighbor(nghbr.d_view, m, mblev.d_view(m), 3,
        i, j, k, is, ie, js, je, ks, ke)) {
      Real lower[3];
      Real upper[3];
      ComputeMagneticVectorPotential(x1f, x2f, x3v-0.25*dx3,
                                     parameters, lower);
      ComputeMagneticVectorPotential(x1f, x2f, x3v+0.25*dx3,
                                     parameters, upper);
      a3(m, k, j, i) = 0.5*(lower[2]+upper[2]);
    }
  });

  auto &b0 = pmbp->pmhd->b0;
  par_for("emri_curl_vector_potential", DevExeSpace(), 0, nmb-1, ks, ke,
          js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real dx1 = size.d_view(m).dx1;
    const Real dx2 = size.d_view(m).dx2;
    const Real dx3 = size.d_view(m).dx3;
    b0.x1f(m, k, j, i) =
        (a3(m, k, j+1, i)-a3(m, k, j, i))/dx2
        -(a2(m, k+1, j, i)-a2(m, k, j, i))/dx3;
    b0.x2f(m, k, j, i) =
        (a1(m, k+1, j, i)-a1(m, k, j, i))/dx3
        -(a3(m, k, j, i+1)-a3(m, k, j, i))/dx1;
    b0.x3f(m, k, j, i) =
        (a2(m, k, j, i+1)-a2(m, k, j, i))/dx1
        -(a1(m, k, j+1, i)-a1(m, k, j, i))/dx2;
    if (i == ie) {
      b0.x1f(m, k, j, i+1) =
          (a3(m, k, j+1, i+1)-a3(m, k, j, i+1))/dx2
          -(a2(m, k+1, j, i+1)-a2(m, k, j, i+1))/dx3;
    }
    if (j == je) {
      b0.x2f(m, k, j+1, i) =
          (a1(m, k+1, j+1, i)-a1(m, k, j+1, i))/dx3
          -(a3(m, k, j+1, i+1)-a3(m, k, j+1, i))/dx1;
    }
    if (k == ke) {
      b0.x3f(m, k+1, j, i) =
          (a2(m, k+1, j, i+1)-a2(m, k+1, j, i))/dx1
          -(a1(m, k+1, j+1, i)-a1(m, k+1, j, i))/dx2;
    }
  });

  auto &bcc0 = pmbp->pmhd->bcc0;
  par_for("emri_cell_centered_field", DevExeSpace(), 0, nmb-1, ks, ke,
          js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    bcc0(m, IBX, k, j, i) =
        0.5*(b0.x1f(m, k, j, i)+b0.x1f(m, k, j, i+1));
    bcc0(m, IBY, k, j, i) =
        0.5*(b0.x2f(m, k, j, i)+b0.x2f(m, k, j+1, i));
    bcc0(m, IBZ, k, j, i) =
        0.5*(b0.x3f(m, k, j, i)+b0.x3f(m, k+1, j, i));
  });
}

Real DistanceSquaredToBlock(const RegionSize &block) {
  const Real dx = std::max(std::max(block.x1min, Real(0.0)), -block.x1max);
  const Real dy = std::max(std::max(block.x2min, Real(0.0)), -block.x2max);
  const Real dz = std::max(std::max(block.x3min, Real(0.0)), -block.x3max);
  return SQR(dx)+SQR(dy)+SQR(dz);
}

Real SecondaryRestEnclosingRadius() {
  const Real horizon = emri_comoving::SecondaryRegularizationRadius(wind_tunnel.metric);
  return std::sqrt(SQR(wind_tunnel.metric.secondary_spin)+SQR(horizon));
}

Real SecondaryEnclosingRadius() {
  return wind_tunnel.secondary_coordinate_stretch*SecondaryRestEnclosingRadius();
}

Real RefinementRadius(const int physical_level) {
  return std::max(refinement_shell_radii[physical_level],
      wind_tunnel.refinement_horizon_factor*SecondaryEnclosingRadius());
}

} // namespace

//----------------------------------------------------------------------------------------
//! Configure the local metric and initialize an analytic magnetized wind profile.

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  wind_profile_table.clear();
  wind_profile_hash = 0;
  wind_profile_time_offset = 0.0;
  wind_profile_current_table_time = std::numeric_limits<Real>::quiet_NaN();
  wind_profile_source_coordinate_radius = std::numeric_limits<Real>::quiet_NaN();
  wind_profile_source_omega = std::numeric_limits<Real>::quiet_NaN();
  wind_profile_orbit_tolerance = std::numeric_limits<Real>::quiet_NaN();
  wind_profile_enabled = false;
  wind_profile_hold_outside = false;
  if (!pmbp->pcoord->is_dynamical_relativistic || pmbp->padm == nullptr
      || pmbp->pdyngr == nullptr || pmbp->pmhd == nullptr) {
    Fatal("emri_windtunnel requires <adm> and <mhd> blocks (dynamical GRMHD)");
  }
  if (pin->DoesBlockExist("radiation")) {
    Fatal("emri_windtunnel does not yet support radiation tetrads");
  }
  for (const auto &input_block : pin->block) {
    if (input_block.block_name.compare(0, 6, "output") != 0) continue;
    const std::string &block_name = input_block.block_name;
    const bool active = pin->DoesParameterExist(block_name, "dcycle")
                            ? pin->GetInteger(block_name, "dcycle") != 0
                            : pin->GetReal(block_name, "dt") > 0.0;
    if (!active) continue;
    const std::string file_type = pin->GetString(block_name, "file_type");
    if (file_type == "hst" || file_type == "rst" || file_type == "log"
        || file_type == "trk") {
      continue;
    }
    bool requests_jcon = pin->GetString(block_name, "variable") == "mhd_jcon";
    if (file_type == "pdf" && pin->DoesParameterExist(block_name, "variable_2")
        && pin->GetOrAddInteger(block_name, "nbin2", 0) > 1) {
      requests_jcon = requests_jcon
          || pin->GetString(block_name, "variable_2") == "mhd_jcon";
    }
    if (requests_jcon) {
      Fatal("emri_windtunnel does not support mhd_jcon, whose tetrad assumes a "
            "single Kerr metric");
    }
  }
  if (!pmbp->pcoord->coord_data.bh_excise
      || pmbp->pcoord->coord_data.excision_scheme != ExcisionScheme::lapse) {
    Fatal("emri_windtunnel requires <coord> excise=true and excision_scheme=lapse");
  }
  if (pmy_mesh_->multilevel && pmy_mesh_->pmr->prolong_prims) {
    Fatal("emri_windtunnel AMR requires prolong_primitives=false");
  }

  pmbp->padm->SetADMVariables = &SetADMVariablesToEMRI;
  pmbp->pcoord->AugmentExcisionMasks = &AugmentEMRIExcisionMasks;
  user_ref_func = RefineSecondary;
  user_hist_func = EMRIHistory;
  user_output_region_func = EMRIOutputRegionCenter;
  user_bcs_func = EMRIWindBoundary;
  user_bcs_level_subcycling_safe = true;

  auto &metric = wind_tunnel.metric;
  metric.primary_mass = pin->GetOrAddReal("problem", "primary_mass", 1.0e5);
  metric.secondary_mass = pin->GetOrAddReal("problem", "secondary_mass", 1.0);
  const Real primary_chi = pin->GetOrAddReal("problem", "primary_chi", 0.0);
  const Real secondary_chi = pin->GetOrAddReal("problem", "secondary_chi", 0.0);
  metric.primary_spin = primary_chi*metric.primary_mass;
  metric.secondary_spin = secondary_chi*metric.secondary_mass;
  metric.orbital_radius = pin->GetOrAddReal(
      "problem", "orbital_radius", 10.0*metric.primary_mass);
  metric.coordinate_radius = std::sqrt(
      SQR(metric.orbital_radius)+SQR(metric.primary_spin));
  const std::string background_mode =
      pin->GetOrAddString("problem", "background_mode", "full");
  if (background_mode == "full") {
    metric.include_primary = true;
    metric.include_orbital_frame = true;
  } else if (background_mode == "frame_only") {
    metric.include_primary = false;
    metric.include_orbital_frame = true;
  } else if (background_mode == "isolated") {
    metric.include_primary = false;
    metric.include_orbital_frame = false;
  } else {
    Fatal("unknown <problem> background_mode: "+background_mode);
  }
  const std::string secondary_embedding =
      pin->GetOrAddString("problem", "secondary_embedding", "tangent_tetrad");
  if (secondary_embedding == "tangent_tetrad") {
    metric.embed_secondary_in_tetrad = true;
  } else if (secondary_embedding == "global_boost") {
    metric.embed_secondary_in_tetrad = false;
  } else {
    Fatal("unknown <problem> secondary_embedding: "+secondary_embedding);
  }
  const int orbit_direction = pin->GetOrAddInteger("problem", "orbit_direction", 1);
  const std::string omega_mode =
      pin->GetOrAddString("problem", "omega_mode", "kerr_geodesic");
  if (omega_mode == "kerr_geodesic") {
    omega_is_geodesic = true;
    metric.omega = emri_comoving::CircularKerrOmega(
        metric.primary_mass, metric.primary_spin, metric.orbital_radius,
        static_cast<Real>(orbit_direction));
    const Real recorded_omega =
        pin->GetOrAddReal("problem", "orbital_omega", metric.omega);
    const Real omega_scale = std::max(std::abs(metric.omega), Real(1.0));
    if (std::abs(recorded_omega-metric.omega)
        > 64.0*std::numeric_limits<Real>::epsilon()*omega_scale) {
      Fatal("orbital_omega cannot override omega_mode=kerr_geodesic; use "
            "omega_mode=custom");
    }
  } else if (omega_mode == "custom") {
    omega_is_geodesic = false;
    metric.omega = pin->GetReal("problem", "orbital_omega");
  } else {
    Fatal("unknown <problem> omega_mode: "+omega_mode);
  }
  metric.spin_buffer_primary =
      pin->GetOrAddReal("problem", "spin_buffer_primary", 0.05);
  metric.spin_buffer_secondary =
      pin->GetOrAddReal("problem", "spin_buffer_secondary", 0.05);
  metric.singularity_floor =
      pin->GetOrAddReal("problem", "singularity_floor", 1.0e-3);
#if SINGLE_PRECISION_ENABLED
  const Real default_fd_step = 5.0e-4*metric.secondary_mass;
#else
  const Real default_fd_step = 5.0e-5*metric.secondary_mass;
#endif
  wind_tunnel.metric_fd_step =
      pin->GetOrAddReal("problem", "metric_fd_step", default_fd_step);
  const Real default_external_fd_step =
      std::cbrt(std::numeric_limits<Real>::epsilon())*metric.primary_mass;
  wind_tunnel.external_metric_fd_step = pin->GetOrAddReal(
      "problem", "external_metric_fd_step", default_external_fd_step);
  wind_tunnel.rho = pin->GetOrAddReal("problem", "rho0", 1.0e-5);
  wind_tunnel.pressure = pin->GetOrAddReal("problem", "pgas0", 1.0e-7);
  wind_tunnel.wind_u[0] = pin->GetOrAddReal("problem", "u1", 0.5);
  wind_tunnel.wind_u[1] = pin->GetOrAddReal("problem", "u2", 0.0);
  wind_tunnel.wind_u[2] = pin->GetOrAddReal("problem", "u3", 0.0);
  wind_tunnel.magnetic_field[0] = pin->GetOrAddReal("problem", "b1", 0.0);
  wind_tunnel.magnetic_field[1] = pin->GetOrAddReal("problem", "b2", 0.0);
  wind_tunnel.magnetic_field[2] = pin->GetOrAddReal("problem", "b3", 1.0e-6);
  wind_tunnel.max_log_contrast =
      pin->GetOrAddReal("problem", "max_log_contrast", 20.0);
  for (int direction=0; direction<3; ++direction) {
    const std::string suffix = std::to_string(direction+1);
    wind_tunnel.log_density_gradient[direction] = pin->GetOrAddReal(
        "problem", "dlnrho_dxh"+suffix, 0.0);
    wind_tunnel.log_pressure_gradient[direction] = pin->GetOrAddReal(
        "problem", "dlnpgas_dxh"+suffix, 0.0);
  }
  for (int component=0; component<3; ++component) {
    for (int direction=0; direction<3; ++direction) {
      const std::string component_suffix = std::to_string(component+1);
      const std::string direction_suffix = std::to_string(direction+1);
      wind_tunnel.velocity_gradient[component][direction] = pin->GetOrAddReal(
          "problem", "du"+component_suffix+"_dxh"+direction_suffix, 0.0);
      const Real magnetic_default = (component == 2 && direction == 2)
          ? -wind_tunnel.magnetic_gradient_source[0][0]
            -wind_tunnel.magnetic_gradient_source[1][1]
          : 0.0;
      wind_tunnel.magnetic_gradient_source[component][direction] =
          pin->GetOrAddReal(
              "problem", "db"+component_suffix+"_dxh"+direction_suffix,
              magnetic_default);
    }
  }
  const std::string wind_frame =
      pin->GetOrAddString("problem", "wind_frame", "source_tetrad");
  if (wind_frame == "source_tetrad") {
    wind_tunnel.wind_is_source_tetrad = true;
  } else if (wind_frame == "normal_frame") {
    wind_tunnel.wind_is_source_tetrad = false;
  } else {
    Fatal("unknown <problem> wind_frame: "+wind_frame);
  }
  const std::string force_frame =
      pin->GetOrAddString("problem", "force_frame", "source_tetrad");
  if (force_frame == "source_tetrad") {
    wind_tunnel.force_is_source_tetrad = true;
  } else if (force_frame == "coordinate") {
    wind_tunnel.force_is_source_tetrad = false;
  } else {
    Fatal("unknown <problem> force_frame: "+force_frame);
  }
  const std::string profile_file =
      pin->GetOrAddString("problem", "profile_file", "");
  if (!profile_file.empty()) {
    if (!wind_tunnel.wind_is_source_tetrad) {
      Fatal("time-dependent EMRI profiles require <problem> wind_frame=source_tetrad");
    }
    if (!pmbp->padm->is_dynamic) {
      Fatal("time-dependent EMRI profiles require <adm> dynamic=true for RK-stage "
            "boundary synchronization");
    }
    if (pin->DoesParameterExist("time", "subcycling")
        && pin->GetString("time", "subcycling") == "level") {
      Fatal("time-dependent EMRI profiles do not yet support time/subcycling=level");
    }
    LoadWindProfileTable(profile_file);
    const Real radius_scale = std::max(
        std::abs(metric.coordinate_radius),
        std::abs(wind_profile_source_coordinate_radius));
    const Real omega_scale = std::max({
        std::abs(metric.omega), std::abs(wind_profile_source_omega),
        std::numeric_limits<Real>::min()});
    if (std::abs(metric.coordinate_radius-wind_profile_source_coordinate_radius)
            > wind_profile_orbit_tolerance*radius_scale
        || std::abs(metric.omega-wind_profile_source_omega)
            > wind_profile_orbit_tolerance*omega_scale) {
      Fatal("EMRI profile orbit metadata is incompatible with the local metric");
    }
    wind_profile_enabled = true;
    const std::string offset_setting = pin->GetOrAddString(
        "problem", "profile_time_offset", "auto");
    if (offset_setting == "auto") {
      wind_profile_time_offset = wind_profile_table.front().time-pmy_mesh_->time;
      pin->SetReal("problem", "profile_time_offset", wind_profile_time_offset);
    } else {
      std::istringstream offset_stream(offset_setting);
      offset_stream >> wind_profile_time_offset >> std::ws;
      if (!offset_stream.eof()) {
        Fatal("<problem> profile_time_offset must be a finite real or auto");
      }
    }
    const std::string extrapolation = pin->GetOrAddString(
        "problem", "profile_extrapolation", "error");
    if (extrapolation == "error") {
      wind_profile_hold_outside = false;
    } else if (extrapolation == "hold") {
      wind_profile_hold_outside = true;
    } else {
      Fatal("unknown <problem> profile_extrapolation: "+extrapolation);
    }
  }
  wind_tunnel.refinement_radius =
      pin->GetOrAddReal("problem", "refinement_radius", 6.0*metric.secondary_mass);
  wind_tunnel.refinement_radius_ratio =
      pin->GetOrAddReal("problem", "refinement_radius_ratio", 2.0);
  wind_tunnel.refinement_hysteresis =
      pin->GetOrAddReal("problem", "refinement_hysteresis", 1.25);
  wind_tunnel.refinement_horizon_factor =
      pin->GetOrAddReal("problem", "refinement_horizon_factor", 1.25);

  const RegionSize &domain = pmy_mesh_->mesh_size;
  const Real inscribed_domain_radius = std::min({
      -domain.x1min, domain.x1max, -domain.x2min, domain.x2max,
      -domain.x3min, domain.x3max});
  wind_tunnel.force_surface_radius = pin->GetOrAddReal(
      "problem", "force_surface_radius", 3.0*metric.secondary_mass);
  wind_tunnel.force_outer_radius[0] = pin->GetOrAddReal(
      "problem", "force_outer_radius_1", 0.5*inscribed_domain_radius);
  wind_tunnel.force_outer_radius[1] = pin->GetOrAddReal(
      "problem", "force_outer_radius_2", 0.75*inscribed_domain_radius);
  wind_tunnel.force_outer_radius[2] = pin->GetOrAddReal(
      "problem", "force_outer_radius_3", 0.95*inscribed_domain_radius);
  wind_tunnel.force_surface_nlevel =
      pin->GetOrAddInteger("problem", "force_surface_nlevel", 5);
  wind_tunnel.force_subtract_background = pin->GetOrAddBoolean(
      "problem", "force_subtract_background", true);
  wind_tunnel.reinitialize_wind_on_restart = pin->GetOrAddBoolean(
      "problem", "reinitialize_wind_on_restart", false);
  if (wind_tunnel.reinitialize_wind_on_restart && !restart) {
    Fatal("<problem> reinitialize_wind_on_restart=true is only valid with -r");
  }
  wind_tunnel.adiabatic_index = pin->GetOrAddReal("mhd", "gamma", 5.0/3.0);
  const std::string dynamic_eos =
      pin->GetOrAddString("mhd", "dyn_eos", "ideal");
  if (user_hist && dynamic_eos != "ideal") {
    Fatal("EMRI force history currently requires <mhd> dyn_eos=ideal");
  }

  std::vector<Real> finite_values = {
    metric.primary_mass, metric.secondary_mass, primary_chi, secondary_chi,
    metric.orbital_radius, metric.coordinate_radius, metric.omega,
    metric.spin_buffer_primary, metric.spin_buffer_secondary,
    metric.singularity_floor, wind_tunnel.metric_fd_step,
    wind_tunnel.external_metric_fd_step, wind_tunnel.rho, wind_tunnel.pressure,
    wind_tunnel.wind_u[0], wind_tunnel.wind_u[1],
    wind_tunnel.wind_u[2], wind_tunnel.magnetic_field[0],
    wind_tunnel.magnetic_field[1], wind_tunnel.magnetic_field[2],
    wind_tunnel.max_log_contrast,
    wind_tunnel.refinement_radius, wind_tunnel.refinement_radius_ratio,
    wind_tunnel.refinement_hysteresis, wind_tunnel.refinement_horizon_factor,
    wind_tunnel.force_surface_radius, wind_tunnel.force_outer_radius[0],
    wind_tunnel.force_outer_radius[1], wind_tunnel.force_outer_radius[2],
    wind_tunnel.adiabatic_index
  };
  if (wind_profile_enabled) finite_values.push_back(wind_profile_time_offset);
  for (int direction=0; direction<3; ++direction) {
    finite_values.push_back(wind_tunnel.log_density_gradient[direction]);
    finite_values.push_back(wind_tunnel.log_pressure_gradient[direction]);
    for (int component=0; component<3; ++component) {
      finite_values.push_back(wind_tunnel.velocity_gradient[component][direction]);
      finite_values.push_back(
          wind_tunnel.magnetic_gradient_source[component][direction]);
    }
  }
  for (const Real value : finite_values) {
    if (!std::isfinite(value)) Fatal("emri_windtunnel parameters must be finite");
  }
  if (!(metric.primary_mass > 0.0) || !(metric.secondary_mass > 0.0)
      || std::abs(primary_chi) > 1.0 || std::abs(secondary_chi) > 1.0
      || !(metric.orbital_radius > 0.0)
      || (orbit_direction != -1 && orbit_direction != 1)
      || !(metric.spin_buffer_primary >= 0.0)
      || !(metric.spin_buffer_secondary >= 0.0)
      || !(metric.singularity_floor > 0.0)
      || metric.spin_buffer_primary+metric.singularity_floor > 1.0
      || metric.spin_buffer_secondary+metric.singularity_floor > 1.0
      || !(wind_tunnel.metric_fd_step > 0.0)
      || !(wind_tunnel.external_metric_fd_step > 0.0)
      || !(wind_tunnel.rho > 0.0) || !(wind_tunnel.pressure > 0.0)
      || !(wind_tunnel.max_log_contrast > 0.0)
      || !(wind_tunnel.refinement_radius > 0.0)
      || !(wind_tunnel.refinement_radius_ratio >= 1.0)
      || !(wind_tunnel.refinement_hysteresis > 1.0)
      || !(wind_tunnel.refinement_horizon_factor > 0.0)
      || !(wind_tunnel.adiabatic_index > 1.0)) {
    Fatal("invalid emri_windtunnel parameter");
  }
  Real source_magnetic_trace = 0.0;
  Real source_magnetic_scale = 0.0;
  for (int component=0; component<3; ++component) {
    source_magnetic_trace +=
        wind_tunnel.magnetic_gradient_source[component][component];
    for (int direction=0; direction<3; ++direction) {
      source_magnetic_scale = std::max(source_magnetic_scale, std::abs(
          wind_tunnel.magnetic_gradient_source[component][direction]));
    }
  }
  const Real magnetic_trace_tolerance =
      256.0*std::numeric_limits<Real>::epsilon()
      *std::max(source_magnetic_scale, Real(1.0));
  if (std::abs(source_magnetic_trace) > magnetic_trace_tolerance) {
    Fatal("dbi_dxhj must be trace-free: set db3_dxh3=-db1_dxh1-db2_dxh2");
  }
  const Real orbital_speed = std::abs(metric.omega*metric.coordinate_radius);
  const Real active_frame_speed = metric.include_orbital_frame ? orbital_speed : 0.0;
  if (!(orbital_speed < 1.0)) {
    Fatal("the effective secondary trajectory is superluminal in global KS coordinates");
  }
  if (pin->GetOrAddBoolean("problem", "require_stable_orbit", true)) {
    const Real isco = KerrISCO(metric.primary_mass, primary_chi,
                               static_cast<Real>(orbit_direction));
    if (metric.orbital_radius < isco) {
      Fatal("orbital_radius lies inside the equatorial Kerr ISCO; set "
            "require_stable_orbit=false only for an explicitly non-equilibrium study");
    }
  }
#if SINGLE_PRECISION_ENABLED
  if (metric.primary_mass/metric.secondary_mass > 1.0e3) {
    Fatal("extreme-mass-ratio wind tunnels require a double-precision build");
  }
#endif
  ConfigureSourceFrame();
  ConfigureMagneticGradient();
  if (wind_profile_enabled) {
    UpdateWindProfileForTime(pmy_mesh_->time, domain);
  } else {
    ValidateSpatialProfiles(domain);
  }

  bool has_analytic_gradient = wind_profile_enabled;
  for (int direction=0; direction<3; ++direction) {
    has_analytic_gradient = has_analytic_gradient
        || wind_tunnel.log_density_gradient[direction] != 0.0
        || wind_tunnel.log_pressure_gradient[direction] != 0.0;
    for (int component=0; component<3; ++component) {
      has_analytic_gradient = has_analytic_gradient
          || wind_tunnel.velocity_gradient[component][direction] != 0.0
          || wind_tunnel.magnetic_gradient_source[component][direction] != 0.0;
    }
  }

  bool has_user_boundary = false;
  for (int face=0; face<6; ++face) {
    has_user_boundary = has_user_boundary
        || pmy_mesh_->mesh_bcs[face] == BoundaryFlag::user;
    if ((wind_tunnel.wind_is_source_tetrad || has_analytic_gradient)
        && pmy_mesh_->mesh_bcs[face] == BoundaryFlag::inflow) {
      Fatal("source-frame or gradient winds require spatially varying user boundaries; "
            "replace each inflow boundary flag with user");
    }
  }
  if ((wind_tunnel.wind_is_source_tetrad || has_analytic_gradient)
      && !has_user_boundary) {
    Fatal("source-frame or gradient winds require at least one user inflow boundary");
  }

  if (!(domain.x1min < 0.0 && domain.x1max > 0.0
        && domain.x2min < 0.0 && domain.x2max > 0.0
        && domain.x3min < 0.0 && domain.x3max > 0.0)) {
    Fatal("the local wind-tunnel domain must contain the secondary at the origin");
  }
  const Real source_surface_inner_bound = metric.embed_secondary_in_tetrad
      ? SecondaryRestEnclosingRadius()
      : wind_tunnel.source_tetrad_stretch*SecondaryEnclosingRadius();
  const Real force_surface_inner_bound = wind_tunnel.force_is_source_tetrad
      ? source_surface_inner_bound : SecondaryEnclosingRadius();
  const Real force_domain_bound = wind_tunnel.force_is_source_tetrad
      ? inscribed_domain_radius/wind_tunnel.source_coordinate_stretch
      : inscribed_domain_radius;
  if (user_hist
      && (!(wind_tunnel.force_surface_radius > force_surface_inner_bound)
          || !(wind_tunnel.force_surface_radius
               < wind_tunnel.force_outer_radius[0])
          || !(wind_tunnel.force_outer_radius[0]
               < wind_tunnel.force_outer_radius[1])
          || !(wind_tunnel.force_outer_radius[1]
               < wind_tunnel.force_outer_radius[2])
          || !(wind_tunnel.force_outer_radius[2] <= force_domain_bound)
          || wind_tunnel.force_surface_nlevel < 0
          || wind_tunnel.force_surface_nlevel > 12)) {
    Fatal("EMRI force radii must enclose the secondary horizon, increase strictly, "
          "and fit inside the largest origin-centered sphere in the domain");
  }
  if (user_hist && std::abs(pmbp->pcoord->coord_data.bh_spin) > 0.0) {
    Fatal("EMRI source-frame force extraction requires <coord> bh_spin=0; "
          "the secondary spin is supplied by <problem> secondary_chi");
  }
  if (metric.include_primary) {
    const Real primary_horizon = metric.primary_mass
        + std::sqrt(SQR(metric.primary_mass)-SQR(metric.primary_spin));
    const Real primary_enclosing = std::sqrt(
        SQR(metric.primary_spin)+SQR(primary_horizon));
    if (domain.x1min <= -metric.coordinate_radius+primary_enclosing) {
      Fatal("the local domain intersects the primary horizon; reduce its radial extent");
    }
  }
  const Real largest_extent = std::max({
      std::abs(domain.x1min), std::abs(domain.x1max),
      std::abs(domain.x2min), std::abs(domain.x2max),
      std::abs(domain.x3min), std::abs(domain.x3max)});
  if (global_variable::my_rank == 0
      && largest_extent > 0.2*metric.orbital_radius) {
    std::cout << "### WARNING: local EMRI box extent/orbital_radius="
              << largest_extent/metric.orbital_radius
              << "; interpretation as a local wind tunnel is becoming marginal"
              << std::endl;
  }

  const int tracked_levels = pmy_mesh_->max_level-pmy_mesh_->root_level;
  refinement_shell_radii.assign(tracked_levels+1, 0.0);
  bool explicit_shells = false;
  for (int level=1; level<=tracked_levels; ++level) {
    const std::string name = "refinement_radius_level_"+std::to_string(level);
    explicit_shells = explicit_shells || pin->DoesParameterExist("problem", name);
  }
  for (int level=1; level<=tracked_levels; ++level) {
    const std::string name = "refinement_radius_level_"+std::to_string(level);
    if (explicit_shells && !pin->DoesParameterExist("problem", name)) {
      Fatal("explicit EMRI refinement shells must specify every configured level");
    }
    const Real default_radius = wind_tunnel.refinement_radius*std::pow(
        wind_tunnel.refinement_radius_ratio, tracked_levels-level);
    refinement_shell_radii[level] = pin->GetOrAddReal(
        "problem", name, default_radius);
    if (!(refinement_shell_radii[level] > 0.0)
        || (level > 1
            && refinement_shell_radii[level] > refinement_shell_radii[level-1])) {
      Fatal("EMRI refinement radii must be positive and non-increasing at finer levels");
    }
  }

  if (user_hist) {
    const Real *surface_transform = wind_tunnel.force_is_source_tetrad
        ? &wind_tunnel.source_spatial_inverse[0][0] : nullptr;
    spherical_grids.push_back(std::make_unique<SphericalGrid>(
        pmbp, wind_tunnel.force_surface_nlevel,
        wind_tunnel.force_surface_radius, -1, surface_transform));
  }

  ValidateOrStoreRestartContract(pin, restart);
  ValidateMetricKernel();
  InitializeInflowArrays(pmbp);
  pmbp->padm->SetADMVariables(pmbp);
  pmbp->pcoord->UpdateExcisionMasks();

  if (global_variable::my_rank == 0) {
    const Real q = metric.secondary_mass/metric.primary_mass;
    const Real hill_radius = metric.orbital_radius*std::cbrt(q/3.0);
    std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10)
              << "Local EMRI wind tunnel: q=" << q
              << ", background=" << background_mode
              << ", secondary embedding=" << secondary_embedding
              << ", r/M=" << metric.orbital_radius/metric.primary_mass
              << ", Omega*M=" << metric.omega*metric.primary_mass
              << ", reference orbital KS speed=" << orbital_speed
              << ", active frame speed=" << active_frame_speed
              << ", wind frame=" << wind_frame
              << ", force frame=" << force_frame
              << ", analytic gradients=" << has_analytic_gradient
              << ", profile replay=" << wind_profile_enabled
              << ", r_Hill/m=" << hill_radius/metric.secondary_mass
              << ", box/orbit=" << largest_extent/metric.orbital_radius
              << ", metric_fd_step/m="
              << wind_tunnel.metric_fd_step/metric.secondary_mass
              << ", external_fd_step/M="
              << wind_tunnel.external_metric_fd_step/metric.primary_mass
              << ", M|Gamma^i_00|=" << primary_geodesic_residual
              << ", force_shell/m={"
              << wind_tunnel.force_surface_radius/metric.secondary_mass << ","
              << wind_tunnel.force_outer_radius[0]/metric.secondary_mass << ","
              << wind_tunnel.force_outer_radius[1]/metric.secondary_mass << ","
              << wind_tunnel.force_outer_radius[2]/metric.secondary_mass << "}"
              << ", ADM cache="
              << (pmbp->padm->is_dynamic ? "stage-refresh" : "stationary")
              << ", restart wind reset="
              << wind_tunnel.reinitialize_wind_on_restart
              << std::endl;
    if (wind_profile_enabled) {
      std::cout << "EMRI profile table: rows=" << wind_profile_table.size()
                << ", time=[" << wind_profile_table.front().time << ","
                << wind_profile_table.back().time << "]"
                << ", offset=" << wind_profile_time_offset
                << ", FNV1a64=" << FormatFNV1a64(wind_profile_hash) << std::endl;
    }
  }
  if (restart && !wind_tunnel.reinitialize_wind_on_restart) return;

  auto &indcs = pmy_mesh_->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int nmb = pmbp->nmb_thispack;
  const int nscalars = pmbp->pmhd->nscalars;
  auto &w0 = pmbp->pmhd->w0;
  auto &size = pmbp->pmb->mb_size;
  const WindTunnelParameters parameters = wind_tunnel;

  par_for("emri_analytic_wind", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min,
                               size.d_view(m).x1max);
    const Real y = CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min,
                               size.d_view(m).x2max);
    const Real z = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min,
                               size.d_view(m).x3max);
    Real density;
    Real pressure;
    ComputeWindThermodynamics(x, y, z, parameters, density, pressure);
    Real primitive_velocity[3];
    ComputeWindPrimitive(x, y, z, parameters, primitive_velocity);
    w0(m, IDN, k, j, i) = density;
    w0(m, IVX, k, j, i) = primitive_velocity[0];
    w0(m, IVY, k, j, i) = primitive_velocity[1];
    w0(m, IVZ, k, j, i) = primitive_velocity[2];
    w0(m, IPR, k, j, i) = pressure;
    for (int n=0; n<nscalars; ++n) w0(m, IYF+n, k, j, i) = 0.0;
  });
  InitializeAnalyticMagneticField(pmbp);
  pmbp->pdyngr->PrimToConInit(is, ie, js, je, ks, ke);
}

namespace {

void SetADMVariablesToEMRI(MeshBlockPack *pmbp) {
  UpdateWindProfileForTime(pmbp->pmesh->time, pmbp->pmesh->mesh_size);
  auto &adm_vars = pmbp->padm->adm;
  auto &size = pmbp->pmb->mb_size;
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int ng = indcs.ng;
  const int is = indcs.is;
  const int js = indcs.js;
  const int ks = indcs.ks;
  const int n1 = indcs.nx1+2*ng;
  const int n2 = (indcs.nx2 > 1) ? indcs.nx2+2*ng : 1;
  const int n3 = (indcs.nx3 > 1) ? indcs.nx3+2*ng : 1;
  auto active_lids = pmbp->active_lids.d_view;
  const int active_offset = pmbp->active_offset;
  const int nmb_active = pmbp->nmb_active;
  const WindTunnelParameters parameters = wind_tunnel;

  par_for_active("emri_update_adm", DevExeSpace(), active_lids, active_offset,
  nmb_active, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min,
                               size.d_view(m).x1max);
    const Real y = CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min,
                               size.d_view(m).x2max);
    const Real z = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min,
                               size.d_view(m).x3max);
    MetricWithDerivatives differentiated;
    ADMPoint point;
    DifferentiateMetric(x, y, z, parameters, differentiated);
    DecomposeMetric(differentiated, point);
    for (int a=0; a<3; ++a) {
      adm_vars.beta_u(m, a, k, j, i) = point.beta[a];
      for (int b=a; b<3; ++b) {
        adm_vars.g_dd(m, a, b, k, j, i) = point.gamma[a][b];
        adm_vars.vK_dd(m, a, b, k, j, i) = point.curvature[a][b];
      }
    }
    adm_vars.alpha(m, k, j, i) = point.alpha;
    adm_vars.psi4(m, k, j, i) = point.psi4;
  });
}

//! Resolve the fixed small-hole center for filtered local binary output.
bool EMRIOutputRegionCenter(const std::string &name, const Real time,
                            Real center[3]) {
  (void)time;
  if (name != "secondary") return false;
  for (int axis=0; axis<3; ++axis) center[axis] = 0.0;
  return true;
}

//----------------------------------------------------------------------------------------
//! Extend the analytic source-frame profile into every user ghost region.  Magnetic
//! ghost faces are sampled from the same trace-free field whose active faces came from
//! curl(A); CT continues to evolve the active face fluxes.

void EMRIWindBoundary(Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &indcs = pm->mb_indcs;
  const int ng = indcs.ng;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int n1 = indcs.nx1+2*ng;
  const int n2 = (indcs.nx2 > 1) ? indcs.nx2+2*ng : 1;
  const int n3 = (indcs.nx3 > 1) ? indcs.nx3+2*ng : 1;
  auto active_lids = pmbp->active_lids.d_view;
  const int active_offset = pmbp->active_offset;
  const int nmb_active = pmbp->nmb_active;
  const int nscalars = pmbp->pmhd->nscalars;
  auto &mb_bcs = pmbp->pmb->mb_bcs;
  auto &size = pmbp->pmb->mb_size;
  auto &w0 = pmbp->pmhd->w0;
  auto &b0 = pmbp->pmhd->b0;
  auto &bcc0 = pmbp->pmhd->bcc0;
  const WindTunnelParameters parameters = wind_tunnel;

  par_for_active("emri_user_b1", DevExeSpace(), active_lids, active_offset,
      nmb_active, 0, n3-1, 0, n2-1, 0, n1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const bool user_ghost =
        (i < is && mb_bcs.d_view(m, BoundaryFace::inner_x1) == BoundaryFlag::user)
        || (i > ie+1
            && mb_bcs.d_view(m, BoundaryFace::outer_x1) == BoundaryFlag::user)
        || (j < js
            && mb_bcs.d_view(m, BoundaryFace::inner_x2) == BoundaryFlag::user)
        || (j > je
            && mb_bcs.d_view(m, BoundaryFace::outer_x2) == BoundaryFlag::user)
        || (k < ks
            && mb_bcs.d_view(m, BoundaryFace::inner_x3) == BoundaryFlag::user)
        || (k > ke
            && mb_bcs.d_view(m, BoundaryFace::outer_x3) == BoundaryFlag::user);
    if (user_ghost) {
      const Real x = LeftEdgeX(i-is, indcs.nx1, size.d_view(m).x1min,
                               size.d_view(m).x1max);
      const Real y = CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min,
                                 size.d_view(m).x2max);
      const Real z = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min,
                                 size.d_view(m).x3max);
      Real field[3];
      ComputeDensitizedMagneticField(x, y, z, parameters, field);
      b0.x1f(m, k, j, i) = field[0];
    }
  });
  par_for_active("emri_user_b2", DevExeSpace(), active_lids, active_offset,
      nmb_active, 0, n3-1, 0, n2, 0, n1-1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const bool user_ghost =
        (i < is && mb_bcs.d_view(m, BoundaryFace::inner_x1) == BoundaryFlag::user)
        || (i > ie
            && mb_bcs.d_view(m, BoundaryFace::outer_x1) == BoundaryFlag::user)
        || (j < js
            && mb_bcs.d_view(m, BoundaryFace::inner_x2) == BoundaryFlag::user)
        || (j > je+1
            && mb_bcs.d_view(m, BoundaryFace::outer_x2) == BoundaryFlag::user)
        || (k < ks
            && mb_bcs.d_view(m, BoundaryFace::inner_x3) == BoundaryFlag::user)
        || (k > ke
            && mb_bcs.d_view(m, BoundaryFace::outer_x3) == BoundaryFlag::user);
    if (user_ghost) {
      const Real x = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min,
                                 size.d_view(m).x1max);
      const Real y = LeftEdgeX(j-js, indcs.nx2, size.d_view(m).x2min,
                               size.d_view(m).x2max);
      const Real z = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min,
                                 size.d_view(m).x3max);
      Real field[3];
      ComputeDensitizedMagneticField(x, y, z, parameters, field);
      b0.x2f(m, k, j, i) = field[1];
    }
  });
  par_for_active("emri_user_b3", DevExeSpace(), active_lids, active_offset,
      nmb_active, 0, n3, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const bool user_ghost =
        (i < is && mb_bcs.d_view(m, BoundaryFace::inner_x1) == BoundaryFlag::user)
        || (i > ie
            && mb_bcs.d_view(m, BoundaryFace::outer_x1) == BoundaryFlag::user)
        || (j < js
            && mb_bcs.d_view(m, BoundaryFace::inner_x2) == BoundaryFlag::user)
        || (j > je
            && mb_bcs.d_view(m, BoundaryFace::outer_x2) == BoundaryFlag::user)
        || (k < ks
            && mb_bcs.d_view(m, BoundaryFace::inner_x3) == BoundaryFlag::user)
        || (k > ke+1
            && mb_bcs.d_view(m, BoundaryFace::outer_x3) == BoundaryFlag::user);
    if (user_ghost) {
      const Real x = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min,
                                 size.d_view(m).x1max);
      const Real y = CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min,
                                 size.d_view(m).x2max);
      const Real z = LeftEdgeX(k-ks, indcs.nx3, size.d_view(m).x3min,
                               size.d_view(m).x3max);
      Real field[3];
      ComputeDensitizedMagneticField(x, y, z, parameters, field);
      b0.x3f(m, k, j, i) = field[2];
    }
  });

  par_for_active("emri_user_wind", DevExeSpace(), active_lids, active_offset,
      nmb_active, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const bool user_ghost =
        (i < is && mb_bcs.d_view(m, BoundaryFace::inner_x1) == BoundaryFlag::user)
        || (i > ie
            && mb_bcs.d_view(m, BoundaryFace::outer_x1) == BoundaryFlag::user)
        || (j < js
            && mb_bcs.d_view(m, BoundaryFace::inner_x2) == BoundaryFlag::user)
        || (j > je
            && mb_bcs.d_view(m, BoundaryFace::outer_x2) == BoundaryFlag::user)
        || (k < ks
            && mb_bcs.d_view(m, BoundaryFace::inner_x3) == BoundaryFlag::user)
        || (k > ke
            && mb_bcs.d_view(m, BoundaryFace::outer_x3) == BoundaryFlag::user);
    if (!user_ghost) return;

    const Real x = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min,
                               size.d_view(m).x1max);
    const Real y = CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min,
                               size.d_view(m).x2max);
    const Real z = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min,
                               size.d_view(m).x3max);
    Real density;
    Real pressure;
    ComputeWindThermodynamics(x, y, z, parameters, density, pressure);
    Real primitive_velocity[3];
    ComputeWindPrimitive(x, y, z, parameters, primitive_velocity);
    w0(m, IDN, k, j, i) = density;
    w0(m, IVX, k, j, i) = primitive_velocity[0];
    w0(m, IVY, k, j, i) = primitive_velocity[1];
    w0(m, IVZ, k, j, i) = primitive_velocity[2];
    w0(m, IPR, k, j, i) = pressure;
    for (int n=0; n<nscalars; ++n) w0(m, IYF+n, k, j, i) = 0.0;
    bcc0(m, IBX, k, j, i) =
        0.5*(b0.x1f(m, k, j, i)+b0.x1f(m, k, j, i+1));
    bcc0(m, IBY, k, j, i) =
        0.5*(b0.x2f(m, k, j, i)+b0.x2f(m, k, j+1, i));
    bcc0(m, IBZ, k, j, i) =
        0.5*(b0.x3f(m, k, j, i)+b0.x3f(m, k+1, j, i));
  });
}

void AugmentEMRIExcisionMasks(MeshBlockPack *pmbp) {
  auto &floor = pmbp->pcoord->excision_floor;
  auto &flux = pmbp->pcoord->excision_flux;
  auto &size = pmbp->pmb->mb_size;
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int ng = indcs.ng;
  const int is = indcs.is;
  const int js = indcs.js;
  const int ks = indcs.ks;
  const int n1 = indcs.nx1+2*ng;
  const int n2 = (indcs.nx2 > 1) ? indcs.nx2+2*ng : 1;
  const int n3 = (indcs.nx3 > 1) ? indcs.nx3+2*ng : 1;
  auto active_lids = pmbp->active_lids.d_view;
  const int active_offset = pmbp->active_offset;
  const int nmb_active = pmbp->nmb_active;
  const WindTunnelParameters parameters = wind_tunnel;
  EmriADMReplayGeometry replay_geometry{};
  const bool numerical_adm = EmriInnerWorldtubeADMReplayGeometry(
      pmbp->pmesh, replay_geometry);
  emri_comoving::MetricParameters metric_left = parameters.metric;
  emri_comoving::MetricParameters metric_right = parameters.metric;
  Real rest_stretch_left = parameters.secondary_rest_stretch;
  Real rest_stretch_right = parameters.secondary_rest_stretch;
  if (numerical_adm) {
    for (int a=0; a<4; ++a) {
      for (int b=0; b<4; ++b) {
        metric_left.secondary_tetrad_covector[a][b] =
            replay_geometry.coframe_left[a][b];
        metric_right.secondary_tetrad_covector[a][b] =
            replay_geometry.coframe_right[a][b];
      }
    }
    metric_left.embed_secondary_in_tetrad = true;
    metric_right.embed_secondary_in_tetrad = true;
    Real unused_coordinate_stretch;
    CoframeSpatialStretches(
        replay_geometry.coframe_left, unused_coordinate_stretch,
        rest_stretch_left);
    CoframeSpatialStretches(
        replay_geometry.coframe_right, unused_coordinate_stretch,
        rest_stretch_right);
  }

  par_for_active("emri_geometric_excision", DevExeSpace(), active_lids, active_offset,
  nmb_active, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min,
                               size.d_view(m).x1max);
    const Real y = CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min,
                               size.d_view(m).x2max);
    const Real z = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min,
                               size.d_view(m).x3max);
    const Real radius2_left = emri_comoving::SecondaryKerrRadiusSquared(
        x, y, z, metric_left);
    const Real radius2_right = numerical_adm
        ? emri_comoving::SecondaryKerrRadiusSquared(x, y, z, metric_right)
        : radius2_left;
    const Real regularization_radius =
        emri_comoving::SecondaryRegularizationRadius(metric_left);
    const Real dx2 = (indcs.nx2 > 1) ? SQR(size.d_view(m).dx2) : 0.0;
    const Real dx3 = (indcs.nx3 > 1) ? SQR(size.d_view(m).dx3) : 0.0;
    const Real padding = 0.5*sqrt(SQR(size.d_view(m).dx1)+dx2+dx3)
                       + parameters.metric_fd_step;
    const Real spin2 = SQR(parameters.metric.secondary_spin);
    const Real enclosing_radius = sqrt(spin2+SQR(regularization_radius));
    Real rest_position_left[3];
    emri_comoving::SecondaryRestFramePosition(
        x, y, z, metric_left, rest_position_left);
    const Real rest_distance_left = sqrt(
        SQR(rest_position_left[0])+SQR(rest_position_left[1])
        +SQR(rest_position_left[2]));
    Real rest_distance_right = rest_distance_left;
    if (numerical_adm) {
      Real rest_position_right[3];
      emri_comoving::SecondaryRestFramePosition(
          x, y, z, metric_right, rest_position_right);
      rest_distance_right = sqrt(
          SQR(rest_position_right[0])+SQR(rest_position_right[1])
          +SQR(rest_position_right[2]));
    }
    floor(m, k, j, i) = floor(m, k, j, i)
                     || radius2_left <= SQR(regularization_radius)
                     || radius2_right <= SQR(regularization_radius);
    flux(m, k, j, i) = flux(m, k, j, i)
                    || rest_distance_left <= enclosing_radius
                                             +rest_stretch_left*padding
                    || rest_distance_right <= enclosing_radius
                                              +rest_stretch_right*padding;
  });
}

void RefineSecondary(MeshBlockPack *pmbp) {
  Mesh *pmesh = pmbp->pmesh;
  auto &refine_flag = pmesh->pmr->refine_flag;
  auto &size = pmbp->pmb->mb_size;
  const int nmb = pmbp->nmb_thispack;
  const int meshblock_start = pmesh->gids_eachrank[global_variable::my_rank];
  for (int m=0; m<nmb; ++m) {
    const int level = pmbp->pmb->mb_lev.h_view(m);
    const int physical_level = level-pmesh->root_level;
    const Real distance2 = DistanceSquaredToBlock(size.h_view(m));
    const bool refine = level < pmesh->max_level
        && distance2 < SQR(RefinementRadius(physical_level+1));
    const bool derefine = level > pmesh->root_level && !refine
        && distance2 > SQR(wind_tunnel.refinement_hysteresis
                           *RefinementRadius(physical_level));
    int &flag = refine_flag.h_view(m+meshblock_start);
    if (refine) flag = 1;
    if (derefine && flag == 0) flag = -1;
  }
  refine_flag.template modify<HostMemSpace>();
  refine_flag.template sync<DevExeSpace>();
}

void EMRIHistory(HistoryData *pdata, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  EmriADMReplayGeometry replay_geometry{};
  const bool numerical_adm_force =
      EmriInnerWorldtubeADMReplayGeometry(pm, replay_geometry);
  if (numerical_adm_force && wind_tunnel.force_is_source_tetrad) {
    Fatal("numerical ADM force history currently requires coordinate force frame");
  }
  pdata->nhist = 20;
  const char *coordinate_labels[20] = {
    "mass_ratio", "orbit_r_M", "omega_M", "mdot",
    "Fmom_x", "Fmom_y", "Fmom_z",
    "Fnewt_x", "Fnewt_y", "Fnewt_z",
    "Frel1_x", "Frel1_y", "Frel1_z",
    "Frel2_x", "Frel2_y", "Frel2_z",
    "Frel3_x", "Frel3_y", "Frel3_z", "geo_resid"
  };
  const char *source_labels[20] = {
    "mass_ratio", "orbit_r_M", "omega_M", "mdot_hat",
    "FmomH_x", "FmomH_y", "FmomH_z",
    "FnewtH_x", "FnewtH_y", "FnewtH_z",
    "Frel1H_x", "Frel1H_y", "Frel1H_z",
    "Frel2H_x", "Frel2H_y", "Frel2H_z",
    "Frel3H_x", "Frel3H_y", "Frel3H_z", "geo_resid"
  };
  for (int n=0; n<pdata->nhist; ++n) {
    pdata->label[n] = wind_tunnel.force_is_source_tetrad
        ? source_labels[n] : coordinate_labels[n];
    pdata->hdata[n] = 0.0;
  }
  const auto &metric_parameters = wind_tunnel.metric;
  if (global_variable::my_rank == 0) {
    pdata->hdata[0] = metric_parameters.secondary_mass
                    /metric_parameters.primary_mass;
    pdata->hdata[1] = metric_parameters.orbital_radius
                    /metric_parameters.primary_mass;
    pdata->hdata[2] = metric_parameters.omega*metric_parameters.primary_mass;
    pdata->hdata[19] = primary_geodesic_residual;
  }

  if (pm->pgen->spherical_grids.size() != 1) {
    Fatal("EMRI force history requires exactly one momentum-flux sphere");
  }
  auto &grid = *pm->pgen->spherical_grids[0];
  DualArray2D<Real> interpolated_field;
  DualArray2D<Real> interpolated_adm;
  grid.InterpolateToSphere(3, pmbp->pmhd->bcc0);
  Kokkos::realloc(interpolated_field, grid.nangles, 3);
  Kokkos::deep_copy(interpolated_field, grid.interp_vals);
  interpolated_field.template modify<DevExeSpace>();
  interpolated_field.template sync<HostMemSpace>();
  if (numerical_adm_force) {
    grid.InterpolateToSphere(adm::ADM::nadm, pmbp->padm->u_adm);
    Kokkos::realloc(interpolated_adm, grid.nangles, adm::ADM::nadm);
    Kokkos::deep_copy(interpolated_adm, grid.interp_vals);
    interpolated_adm.template modify<DevExeSpace>();
    interpolated_adm.template sync<HostMemSpace>();
  }
  grid.InterpolateToSphere(0, IPR, pmbp->pmhd->w0);

  const WindTunnelParameters parameters = wind_tunnel;
  emri_comoving::MetricParameters secondary_metric_left = parameters.metric;
  emri_comoving::MetricParameters secondary_metric_right = parameters.metric;
  if (numerical_adm_force) {
    for (int a=0; a<4; ++a) {
      for (int b=0; b<4; ++b) {
        secondary_metric_left.secondary_tetrad_covector[a][b] =
            replay_geometry.coframe_left[a][b];
        secondary_metric_right.secondary_tetrad_covector[a][b] =
            replay_geometry.coframe_right[a][b];
      }
    }
    secondary_metric_left.embed_secondary_in_tetrad = true;
    secondary_metric_right.embed_secondary_in_tetrad = true;
    Real coordinate_stretch_left;
    Real coordinate_stretch_right;
    Real unused_rest_stretch;
    CoframeSpatialStretches(
        replay_geometry.coframe_left, coordinate_stretch_left,
        unused_rest_stretch);
    CoframeSpatialStretches(
        replay_geometry.coframe_right, coordinate_stretch_right,
        unused_rest_stretch);
    const Real regularization_radius =
        emri_comoving::SecondaryRegularizationRadius(secondary_metric_left);
    const Real rest_enclosing_radius = sqrt(
        SQR(secondary_metric_left.secondary_spin)+SQR(regularization_radius));
    const Real coordinate_enclosing_radius = rest_enclosing_radius*std::max(
        coordinate_stretch_left, coordinate_stretch_right);
    if (!(parameters.force_surface_radius > coordinate_enclosing_radius)) {
      Fatal("numerical ADM force surface does not enclose the replayed secondary");
    }
  }
  const Real surface_radius2 = SQR(parameters.force_surface_radius);
  for (int n=0; n<grid.nangles; ++n) {
    const Real x = grid.interp_coord.h_view(n, 0);
    const Real y = grid.interp_coord.h_view(n, 1);
    const Real z = grid.interp_coord.h_view(n, 2);
    const Real coordinate[3] = {x, y, z};
    const Real inverse_radius = 1.0/parameters.force_surface_radius;
    Real surface_covector[3];
    if (parameters.force_is_source_tetrad) {
      Real source_normal[3] = {0.0, 0.0, 0.0};
      for (int source_axis=0; source_axis<3; ++source_axis) {
        for (int coordinate_axis=0; coordinate_axis<3; ++coordinate_axis) {
          source_normal[source_axis] +=
              parameters.metric.secondary_tetrad_covector[source_axis+1]
                  [coordinate_axis+1]*coordinate[coordinate_axis]*inverse_radius;
        }
      }
      for (int coordinate_axis=0; coordinate_axis<3; ++coordinate_axis) {
        surface_covector[coordinate_axis] = 0.0;
        for (int source_axis=0; source_axis<3; ++source_axis) {
          surface_covector[coordinate_axis] +=
              parameters.source_surface_covector[source_axis][coordinate_axis]
              *source_normal[source_axis];
        }
      }
    } else {
      for (int axis=0; axis<3; ++axis) {
        surface_covector[axis] = coordinate[axis]*inverse_radius;
      }
    }
    const Real velocity[3] = {
      grid.interp_vals.h_view(n, IVX),
      grid.interp_vals.h_view(n, IVY),
      grid.interp_vals.h_view(n, IVZ)
    };
    const Real densitized_field[3] = {
      interpolated_field.h_view(n, IBX),
      interpolated_field.h_view(n, IBY),
      interpolated_field.h_view(n, IBZ)
    };
    const Real density = grid.interp_vals.h_view(n, IDN);
    const Real pressure = grid.interp_vals.h_view(n, IPR);
    if (density == 0.0 && pressure == 0.0) continue;
    Real metric[4][4];
    Real four_velocity[4];
    Real stress[4][4];
    Real sqrt_gamma;
    Real sqrt_minus_g;
    if (numerical_adm_force) {
      Real gamma[3][3];
      const int gamma_index[3][3] = {
        {adm::ADM::I_ADM_GXX, adm::ADM::I_ADM_GXY, adm::ADM::I_ADM_GXZ},
        {adm::ADM::I_ADM_GXY, adm::ADM::I_ADM_GYY, adm::ADM::I_ADM_GYZ},
        {adm::ADM::I_ADM_GXZ, adm::ADM::I_ADM_GYZ, adm::ADM::I_ADM_GZZ}
      };
      for (int i=0; i<3; ++i) {
        for (int j=0; j<3; ++j) {
          gamma[i][j] = interpolated_adm.h_view(n, gamma_index[i][j]);
        }
      }
      const Real beta[3] = {
        interpolated_adm.h_view(n, adm::ADM::I_ADM_BETAX),
        interpolated_adm.h_view(n, adm::ADM::I_ADM_BETAY),
        interpolated_adm.h_view(n, adm::ADM::I_ADM_BETAZ)
      };
      ADMToFourMetric(
          gamma, interpolated_adm.h_view(n, adm::ADM::I_ADM_ALPHA),
          beta, metric);
      ComputeStressEnergyFromMetric(
          metric, parameters, density, pressure, velocity, densitized_field,
          four_velocity, stress, sqrt_gamma, sqrt_minus_g);
    } else {
      ComputeStressEnergy(x, y, z, parameters, density, pressure, velocity,
                          densitized_field, metric, four_velocity, stress,
                          sqrt_gamma, sqrt_minus_g);
    }
    const Real area_weight = surface_radius2
        *grid.solid_angles.h_view(n)*sqrt_minus_g;
    Real radial_velocity = 0.0;
    for (int j=0; j<3; ++j) {
      radial_velocity += surface_covector[j]*four_velocity[j+1];
    }
    pdata->hdata[3] -= density*radial_velocity*area_weight;

    // Form the perturbative momentum flux on the same surface quadrature used for
    // the evolved state.  The analytic wind is evaluated at the quadrature point
    // and contracted with the same metric, so this removes the leading upstream
    // stress on the finite sphere.  A fixed-tree control measures the remaining
    // interpolation/quadrature residual.  mdot intentionally remains the total mass
    // flux rather than a perturbative diagnostic.
    Real background_stress[4][4] = {{0.0}};
    if (parameters.force_subtract_background) {
      Real background_density;
      Real background_pressure;
      ComputeWindThermodynamics(
          x, y, z, parameters, background_density, background_pressure);
      Real background_velocity[3];
      ComputeWindPrimitive(x, y, z, parameters, background_velocity);
      Real background_field[3];
      ComputeDensitizedMagneticField(x, y, z, parameters, background_field);
      Real background_four_velocity[4];
      Real background_sqrt_gamma;
      Real background_sqrt_minus_g;
      ComputeStressEnergyFromMetric(
          metric, parameters, background_density, background_pressure,
          background_velocity, background_field, background_four_velocity,
          background_stress, background_sqrt_gamma, background_sqrt_minus_g);
    }
    Real momentum_covector[4] = {0.0, 0.0, 0.0, 0.0};
    for (int covector_component=0; covector_component<4; ++covector_component) {
      for (int surface_direction=0; surface_direction<3; ++surface_direction) {
        Real mixed_stress = 0.0;
        for (int a=0; a<4; ++a) {
          mixed_stress +=
              (stress[surface_direction+1][a]
               -background_stress[surface_direction+1][a])
              *metric[a][covector_component];
        }
        momentum_covector[covector_component] +=
            surface_covector[surface_direction]*mixed_stress;
      }
    }
    if (parameters.force_is_source_tetrad) {
      for (int source_axis=0; source_axis<3; ++source_axis) {
        Real source_momentum = 0.0;
        for (int component=0; component<4; ++component) {
          source_momentum += parameters.source_tetrad[source_axis+1][component]
                           *momentum_covector[component];
        }
        pdata->hdata[4+source_axis] +=
            parameters.source_dt_dtau*source_momentum*area_weight;
      }
    } else {
      for (int coordinate_axis=0; coordinate_axis<3; ++coordinate_axis) {
        pdata->hdata[4+coordinate_axis] +=
            momentum_covector[coordinate_axis+1]*area_weight;
      }
    }
  }

  auto &w0 = pmbp->pmhd->w0;
  auto &bcc = pmbp->pmhd->bcc0;
  auto adm_variables = pmbp->padm->adm;
  auto &excision_floor = pmbp->pcoord->excision_floor;
  auto &size = pmbp->pmb->mb_size;
  auto &indcs = pm->mb_indcs;
  const int is = indcs.is;
  const int js = indcs.js;
  const int ks = indcs.ks;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;
  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;
  const int nmkji = pmbp->nmb_thispack*nkji;
  const Real outer_radius2[3] = {
    SQR(parameters.force_outer_radius[0]),
    SQR(parameters.force_outer_radius[1]),
    SQR(parameters.force_outer_radius[2])
  };

  array_sum::GlobalSum force_integrals;
  Kokkos::parallel_reduce(
    "emri_force_volume", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
    KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &sum) {
      const int m = idx/nkji;
      const int k0 = (idx-m*nkji)/nji;
      const int j0 = (idx-m*nkji-k0*nji)/nx1;
      const int i0 = idx-m*nkji-k0*nji-j0*nx1;
      const int i = i0+is;
      const int j = j0+js;
      const int k = k0+ks;
      if (excision_floor(m,k,j,i)) return;

      const Real x = CellCenterX(i0, nx1, size.d_view(m).x1min,
                                size.d_view(m).x1max);
      const Real y = CellCenterX(j0, nx2, size.d_view(m).x2min,
                                size.d_view(m).x2max);
      const Real z = CellCenterX(k0, nx3, size.d_view(m).x3min,
                                size.d_view(m).x3max);
      const Real coordinate[3] = {x, y, z};
      Real force_position[3] = {x, y, z};
      if (parameters.force_is_source_tetrad) {
        for (int source_axis=0; source_axis<3; ++source_axis) {
          force_position[source_axis] = 0.0;
          for (int coordinate_axis=0; coordinate_axis<3; ++coordinate_axis) {
            force_position[source_axis] +=
                parameters.metric.secondary_tetrad_covector[source_axis+1]
                    [coordinate_axis+1]*coordinate[coordinate_axis];
          }
        }
      }
      const Real radius2 = SQR(force_position[0])+SQR(force_position[1])
                         +SQR(force_position[2]);
      if (radius2 <= surface_radius2 || radius2 > outer_radius2[2]) return;

      const Real velocity[3] = {
        w0(m, IVX, k, j, i), w0(m, IVY, k, j, i), w0(m, IVZ, k, j, i)
      };
      const Real densitized_field[3] = {
        bcc(m, IBX, k, j, i), bcc(m, IBY, k, j, i), bcc(m, IBZ, k, j, i)
      };
      const Real density = w0(m, IDN, k, j, i);
      const Real pressure = w0(m, IPR, k, j, i);
      Real local_metric[4][4];
      Real four_velocity[4];
      Real stress[4][4];
      Real sqrt_gamma;
      Real sqrt_minus_g;
      if (numerical_adm_force) {
        Real gamma[3][3];
        Real beta[3];
        for (int a=0; a<3; ++a) {
          beta[a] = adm_variables.beta_u(m, a, k, j, i);
          for (int b=0; b<3; ++b) {
            gamma[a][b] = adm_variables.g_dd(m, a, b, k, j, i);
          }
        }
        ADMToFourMetric(
            gamma, adm_variables.alpha(m, k, j, i), beta, local_metric);
        ComputeStressEnergyFromMetric(
            local_metric, parameters, density, pressure, velocity,
            densitized_field, four_velocity, stress, sqrt_gamma, sqrt_minus_g);
      } else {
        ComputeStressEnergy(x, y, z, parameters, density, pressure, velocity,
                            densitized_field, local_metric, four_velocity, stress,
                            sqrt_gamma, sqrt_minus_g);
      }
      Real background_density = 0.0;
      Real background_stress[4][4] = {{0.0}};
      if (parameters.force_subtract_background) {
        Real background_pressure;
        ComputeWindThermodynamics(
            x, y, z, parameters, background_density, background_pressure);
        Real background_velocity[3];
        ComputeWindPrimitive(x, y, z, parameters, background_velocity);
        Real background_field[3];
        ComputeDensitizedMagneticField(x, y, z, parameters, background_field);
        Real background_four_velocity[4];
        Real background_sqrt_gamma;
        Real background_sqrt_minus_g;
        ComputeStressEnergyFromMetric(
            local_metric, parameters, background_density, background_pressure,
            background_velocity, background_field, background_four_velocity,
            background_stress, background_sqrt_gamma, background_sqrt_minus_g);
      }
      Real metric_derivative[3][4][4];
      if (numerical_adm_force) {
        Real derivative_left[3][4][4];
        Real derivative_right[3][4][4];
        DifferentiateSecondaryDisplacement(
            x, y, z, secondary_metric_left, parameters.metric_fd_step,
            derivative_left);
        DifferentiateSecondaryDisplacement(
            x, y, z, secondary_metric_right, parameters.metric_fd_step,
            derivative_right);
        for (int direction=0; direction<3; ++direction) {
          for (int a=0; a<4; ++a) {
            for (int b=0; b<4; ++b) {
              metric_derivative[direction][a][b] =
                  (1.0-replay_geometry.time_fraction)
                      *derivative_left[direction][a][b]
                  +replay_geometry.time_fraction
                      *derivative_right[direction][a][b];
            }
          }
        }
      } else {
        DifferentiateSecondaryDisplacement(
            x, y, z, parameters.metric, parameters.metric_fd_step,
            metric_derivative);
      }
      const Real coordinate_volume = size.d_view(m).dx1*size.d_view(m).dx2
                                   *size.d_view(m).dx3;
      const Real radius = sqrt(radius2);
      const Real source_density = parameters.force_subtract_background
          ? density-background_density : density;
      const Real newtonian_volume = parameters.force_is_source_tetrad
          ? parameters.source_spatial_determinant*coordinate_volume
          : sqrt_gamma*coordinate_volume;
      const Real newtonian_coefficient = parameters.metric.secondary_mass
          *source_density*newtonian_volume/(radius2*radius);

      array_sum::GlobalSum cell;
      cell.the_array[7] = newtonian_coefficient*force_position[0];
      cell.the_array[8] = newtonian_coefficient*force_position[1];
      cell.the_array[9] = newtonian_coefficient*force_position[2];
      for (int force_direction=0; force_direction<3; ++force_direction) {
        Real contraction = 0.0;
        for (int a=0; a<4; ++a) {
          for (int b=0; b<4; ++b) {
            contraction += (stress[a][b]-background_stress[a][b])
                         *metric_derivative[force_direction][a][b];
          }
        }
        const Real relativistic_force = -0.5*contraction*sqrt_minus_g
                                      *coordinate_volume;
        for (int cutoff=0; cutoff<3; ++cutoff) {
          if (radius2 <= outer_radius2[cutoff]) {
            cell.the_array[10+3*cutoff+force_direction] = relativistic_force;
          }
        }
      }
      sum += cell;
    }, Kokkos::Sum<array_sum::GlobalSum>(force_integrals));

  for (int n=7; n<=18; ++n) {
    pdata->hdata[n] += force_integrals.the_array[n];
  }

  if (wind_tunnel.force_is_source_tetrad) {
    pdata->hdata[3] *= wind_tunnel.source_dt_dtau;
    // Momentum flux and the Newtonian estimator are accumulated directly in the
    // source frame; the relativistic generalized forces remain coordinate covectors.
    const int vector_starts[3] = {10, 13, 16};
    for (const int start : vector_starts) {
      const Real coordinate_force[3] = {
        pdata->hdata[start], pdata->hdata[start+1], pdata->hdata[start+2]
      };
      for (int source_axis=0; source_axis<3; ++source_axis) {
        pdata->hdata[start+source_axis] = 0.0;
        for (int coordinate_axis=0; coordinate_axis<3; ++coordinate_axis) {
          pdata->hdata[start+source_axis] +=
              wind_tunnel.force_projection[source_axis][coordinate_axis]
              *coordinate_force[coordinate_axis];
        }
      }
    }
  }
}

} // namespace
