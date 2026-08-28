//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file emri_windtunnel.cpp
//! \brief Local GRMHD wind tunnel centered on a circular equatorial EMRI secondary.

#include <algorithm>
#include <cmath>
#include <cstdlib>
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

namespace {

struct WindTunnelParameters {
  emri_comoving::MetricParameters metric;
  Real metric_fd_step;
  Real external_metric_fd_step;
  Real rho;
  Real pressure;
  Real wind_u[3];
  Real magnetic_field[3];
  Real refinement_radius;
  Real refinement_radius_ratio;
  Real refinement_hysteresis;
  Real refinement_horizon_factor;
  Real force_surface_radius;
  Real force_outer_radius[3];
  Real adiabatic_index;
  int force_surface_nlevel;
  bool force_subtract_background;
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

KOKKOS_INLINE_FUNCTION
void EvaluateMetric(const Real x, const Real y, const Real z,
                    const WindTunnelParameters &parameters, Real gcov[4][4]) {
  emri_comoving::ComputeMetric(x, y, z, parameters.metric, gcov);
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
    const WindTunnelParameters &parameters, Real derivative[3][4][4]) {
  const Real coordinate[3] = {x, y, z};
  for (int direction=0; direction<3; ++direction) {
    Real lower[3] = {coordinate[0], coordinate[1], coordinate[2]};
    Real upper[3] = {coordinate[0], coordinate[1], coordinate[2]};
    lower[direction] -= parameters.metric_fd_step;
    upper[direction] += parameters.metric_fd_step;
    Real metric_lower[4][4];
    Real metric_upper[4][4];
    emri_comoving::ComputeSecondaryMetricPerturbationAtDisplacement(
        lower[0], lower[1], lower[2], x, y, parameters.metric, metric_lower);
    emri_comoving::ComputeSecondaryMetricPerturbationAtDisplacement(
        upper[0], upper[1], upper[2], x, y, parameters.metric, metric_upper);
    const Real inverse_width = 1.0/(upper[direction]-lower[direction]);
    for (int a=0; a<4; ++a) {
      for (int b=0; b<4; ++b) {
        derivative[direction][a][b] =
            (metric_upper[a][b]-metric_lower[a][b])*inverse_width;
      }
    }
  }
}

//! Reconstruct T^{mu nu} from dynamical-GRMHD primitives.  The velocity primitive is
//! W v^i measured by the Eulerian observer, while the CT magnetic field is densitized
//! by sqrt(gamma).
KOKKOS_INLINE_FUNCTION
void ComputeStressEnergy(
    const Real x, const Real y, const Real z,
    const WindTunnelParameters &parameters, const Real density, const Real pressure,
    const Real velocity[3], const Real densitized_field[3], Real metric[4][4],
    Real four_velocity[4], Real stress[4][4], Real &sqrt_gamma,
    Real &sqrt_minus_g) {
  EvaluateMetric(x, y, z, parameters, metric);
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

std::string MetricRestartContract() {
  const auto &metric = wind_tunnel.metric;
  std::ostringstream contract;
  contract << std::setprecision(std::numeric_limits<Real>::max_digits10)
           << "emri-comoving-v2"
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
           << ";metric_fd_step=" << wind_tunnel.metric_fd_step
           << ";external_metric_fd_step=" << wind_tunnel.external_metric_fd_step;
  return contract.str();
}

void ValidateOrStoreRestartContract(ParameterInput *pin, const bool restart) {
  constexpr const char *name = "emri_metric_restart_contract";
  const std::string current = MetricRestartContract();
  if (restart) {
    if (!pin->DoesParameterExist("problem", name)
        || pin->GetString("problem", name) != current) {
      Fatal("EMRI metric parameters differ from the restart contract");
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
      100.0, 0.0, 0.0, weak_field, displacement_derivative);
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
    b_in.h_view(IBX, face) = wind_tunnel.magnetic_field[0];
    b_in.h_view(IBY, face) = wind_tunnel.magnetic_field[1];
    b_in.h_view(IBZ, face) = wind_tunnel.magnetic_field[2];
  }
  u_in.template modify<HostMemSpace>();
  u_in.template sync<DevExeSpace>();
  b_in.template modify<HostMemSpace>();
  b_in.template sync<DevExeSpace>();
}

Real DistanceSquaredToBlock(const RegionSize &block) {
  const Real dx = std::max(std::max(block.x1min, Real(0.0)), -block.x1max);
  const Real dy = std::max(std::max(block.x2min, Real(0.0)), -block.x2max);
  const Real dz = std::max(std::max(block.x3min, Real(0.0)), -block.x3max);
  return SQR(dx)+SQR(dy)+SQR(dz);
}

Real SecondaryEnclosingRadius() {
  const Real horizon = emri_comoving::SecondaryRegularizationRadius(wind_tunnel.metric);
  const Real rest_enclosing = std::sqrt(
      SQR(wind_tunnel.metric.secondary_spin)+SQR(horizon));
  const Real speed = wind_tunnel.metric.include_orbital_frame
      ? wind_tunnel.metric.omega*wind_tunnel.metric.coordinate_radius : 0.0;
  return rest_enclosing/std::sqrt(1.0-SQR(speed));
}

Real RefinementRadius(const int physical_level) {
  return std::max(refinement_shell_radii[physical_level],
      wind_tunnel.refinement_horizon_factor*SecondaryEnclosingRadius());
}

} // namespace

//----------------------------------------------------------------------------------------
//! Configure the local metric and initialize a uniform magnetized wind.

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
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
  wind_tunnel.adiabatic_index = pin->GetOrAddReal("mhd", "gamma", 5.0/3.0);
  const std::string dynamic_eos =
      pin->GetOrAddString("mhd", "dyn_eos", "ideal");
  if (user_hist && dynamic_eos != "ideal") {
    Fatal("EMRI force history currently requires <mhd> dyn_eos=ideal");
  }

  const Real finite_values[] = {
    metric.primary_mass, metric.secondary_mass, primary_chi, secondary_chi,
    metric.orbital_radius, metric.coordinate_radius, metric.omega,
    metric.spin_buffer_primary, metric.spin_buffer_secondary,
    metric.singularity_floor, wind_tunnel.metric_fd_step,
    wind_tunnel.external_metric_fd_step, wind_tunnel.rho, wind_tunnel.pressure,
    wind_tunnel.wind_u[0], wind_tunnel.wind_u[1],
    wind_tunnel.wind_u[2], wind_tunnel.magnetic_field[0],
    wind_tunnel.magnetic_field[1], wind_tunnel.magnetic_field[2],
    wind_tunnel.refinement_radius, wind_tunnel.refinement_radius_ratio,
    wind_tunnel.refinement_hysteresis, wind_tunnel.refinement_horizon_factor,
    wind_tunnel.force_surface_radius, wind_tunnel.force_outer_radius[0],
    wind_tunnel.force_outer_radius[1], wind_tunnel.force_outer_radius[2],
    wind_tunnel.adiabatic_index
  };
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
      || !(wind_tunnel.refinement_radius > 0.0)
      || !(wind_tunnel.refinement_radius_ratio >= 1.0)
      || !(wind_tunnel.refinement_hysteresis > 1.0)
      || !(wind_tunnel.refinement_horizon_factor > 0.0)
      || !(wind_tunnel.adiabatic_index > 1.0)) {
    Fatal("invalid emri_windtunnel parameter");
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

  if (!(domain.x1min < 0.0 && domain.x1max > 0.0
        && domain.x2min < 0.0 && domain.x2max > 0.0
        && domain.x3min < 0.0 && domain.x3max > 0.0)) {
    Fatal("the local wind-tunnel domain must contain the secondary at the origin");
  }
  if (user_hist
      && (!(wind_tunnel.force_surface_radius > SecondaryEnclosingRadius())
          || !(wind_tunnel.force_surface_radius
               < wind_tunnel.force_outer_radius[0])
          || !(wind_tunnel.force_outer_radius[0]
               < wind_tunnel.force_outer_radius[1])
          || !(wind_tunnel.force_outer_radius[1]
               < wind_tunnel.force_outer_radius[2])
          || !(wind_tunnel.force_outer_radius[2] <= inscribed_domain_radius)
          || wind_tunnel.force_surface_nlevel < 0
          || wind_tunnel.force_surface_nlevel > 12)) {
    Fatal("EMRI force radii must enclose the secondary horizon, increase strictly, "
          "and fit inside the largest origin-centered sphere in the domain");
  }
  if (user_hist && std::abs(pmbp->pcoord->coord_data.bh_spin) > 0.0) {
    Fatal("EMRI force extraction uses local Cartesian spheres and requires "
          "<coord> bh_spin=0");
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
    spherical_grids.push_back(std::make_unique<SphericalGrid>(
        pmbp, wind_tunnel.force_surface_nlevel,
        wind_tunnel.force_surface_radius));
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
              << ", r/M=" << metric.orbital_radius/metric.primary_mass
              << ", Omega*M=" << metric.omega*metric.primary_mass
              << ", reference orbital KS speed=" << orbital_speed
              << ", active frame speed=" << active_frame_speed
              << ", r_H/m=" << hill_radius/metric.secondary_mass
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
              << std::endl;
  }
  if (restart) return;

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
  auto &b0 = pmbp->pmhd->b0;
  auto &bcc0 = pmbp->pmhd->bcc0;
  const WindTunnelParameters parameters = wind_tunnel;

  par_for("emri_uniform_wind", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    w0(m, IDN, k, j, i) = parameters.rho;
    w0(m, IVX, k, j, i) = parameters.wind_u[0];
    w0(m, IVY, k, j, i) = parameters.wind_u[1];
    w0(m, IVZ, k, j, i) = parameters.wind_u[2];
    w0(m, IPR, k, j, i) = parameters.pressure;
    for (int n=0; n<nscalars; ++n) w0(m, IYF+n, k, j, i) = 0.0;
    bcc0(m, IBX, k, j, i) = parameters.magnetic_field[0];
    bcc0(m, IBY, k, j, i) = parameters.magnetic_field[1];
    bcc0(m, IBZ, k, j, i) = parameters.magnetic_field[2];
  });
  Kokkos::deep_copy(DevExeSpace(), b0.x1f, wind_tunnel.magnetic_field[0]);
  Kokkos::deep_copy(DevExeSpace(), b0.x2f, wind_tunnel.magnetic_field[1]);
  Kokkos::deep_copy(DevExeSpace(), b0.x3f, wind_tunnel.magnetic_field[2]);
  pmbp->pdyngr->PrimToConInit(is, ie, js, je, ks, ke);
}

namespace {

void SetADMVariablesToEMRI(MeshBlockPack *pmbp) {
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

  par_for_active("emri_geometric_excision", DevExeSpace(), active_lids, active_offset,
  nmb_active, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min,
                               size.d_view(m).x1max);
    const Real y = CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min,
                               size.d_view(m).x2max);
    const Real z = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min,
                               size.d_view(m).x3max);
    const Real radius2 = emri_comoving::SecondaryKerrRadiusSquared(
        x, y, z, parameters.metric);
    const Real regularization_radius =
        emri_comoving::SecondaryRegularizationRadius(parameters.metric);
    const Real dx2 = (indcs.nx2 > 1) ? SQR(size.d_view(m).dx2) : 0.0;
    const Real dx3 = (indcs.nx3 > 1) ? SQR(size.d_view(m).dx3) : 0.0;
    const Real padding = 0.5*sqrt(SQR(size.d_view(m).dx1)+dx2+dx3)
                       + parameters.metric_fd_step;
    const Real spin2 = SQR(parameters.metric.secondary_spin);
    const Real enclosing_radius = sqrt(spin2+SQR(regularization_radius));
    const Real orbital_speed = parameters.metric.include_orbital_frame
        ? parameters.metric.omega*parameters.metric.coordinate_radius : 0.0;
    const Real gamma = 1.0/sqrt(1.0-SQR(orbital_speed));
    const Real velocity[3] = {0.0, orbital_speed, 0.0};
    const Real origin[3] = {0.0, 0.0, 0.0};
    Real rest_position[3];
    binary_bh::RestFramePosition(x, y, z, origin, velocity, rest_position);
    const Real rest_distance = sqrt(SQR(rest_position[0])+SQR(rest_position[1])
                                  +SQR(rest_position[2]));
    floor(m, k, j, i) = floor(m, k, j, i)
                     || radius2 <= SQR(regularization_radius);
    flux(m, k, j, i) = flux(m, k, j, i)
                    || rest_distance <= enclosing_radius+gamma*padding;
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
  pdata->nhist = 20;
  const char *labels[20] = {
    "mass_ratio", "orbit_r_M", "omega_M", "mdot",
    "Fmom_x", "Fmom_y", "Fmom_z",
    "Fnewt_x", "Fnewt_y", "Fnewt_z",
    "Frel1_x", "Frel1_y", "Frel1_z",
    "Frel2_x", "Frel2_y", "Frel2_z",
    "Frel3_x", "Frel3_y", "Frel3_z", "geo_resid"
  };
  for (int n=0; n<pdata->nhist; ++n) {
    pdata->label[n] = labels[n];
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

  MeshBlockPack *pmbp = pm->pmb_pack;
  if (pm->pgen->spherical_grids.size() != 1) {
    Fatal("EMRI force history requires exactly one momentum-flux sphere");
  }
  auto &grid = *pm->pgen->spherical_grids[0];
  DualArray2D<Real> interpolated_field;
  grid.InterpolateToSphere(3, pmbp->pmhd->bcc0);
  Kokkos::realloc(interpolated_field, grid.nangles, 3);
  Kokkos::deep_copy(interpolated_field, grid.interp_vals);
  interpolated_field.template modify<DevExeSpace>();
  interpolated_field.template sync<HostMemSpace>();
  grid.InterpolateToSphere(0, IPR, pmbp->pmhd->w0);

  const WindTunnelParameters parameters = wind_tunnel;
  const Real surface_radius2 = SQR(parameters.force_surface_radius);
  for (int n=0; n<grid.nangles; ++n) {
    const Real x = grid.interp_coord.h_view(n, 0);
    const Real y = grid.interp_coord.h_view(n, 1);
    const Real z = grid.interp_coord.h_view(n, 2);
    const Real inverse_radius = 1.0/parameters.force_surface_radius;
    const Real normal[3] = {x*inverse_radius, y*inverse_radius,
                            z*inverse_radius};
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
    Real metric[4][4];
    Real four_velocity[4];
    Real stress[4][4];
    Real sqrt_gamma;
    Real sqrt_minus_g;
    ComputeStressEnergy(x, y, z, parameters, density, pressure, velocity,
                        densitized_field, metric, four_velocity, stress,
                        sqrt_gamma, sqrt_minus_g);
    const Real area_weight = surface_radius2
        *grid.solid_angles.h_view(n)*sqrt_minus_g;
    Real radial_velocity = 0.0;
    for (int j=0; j<3; ++j) radial_velocity += normal[j]*four_velocity[j+1];
    pdata->hdata[3] -= density*radial_velocity*area_weight;

    for (int force_direction=0; force_direction<3; ++force_direction) {
      Real radial_momentum = 0.0;
      for (int surface_direction=0; surface_direction<3; ++surface_direction) {
        Real mixed_stress = 0.0;
        for (int a=0; a<4; ++a) {
          mixed_stress += stress[surface_direction+1][a]
                         *metric[a][force_direction+1];
        }
        radial_momentum += normal[surface_direction]*mixed_stress;
      }
      pdata->hdata[4+force_direction] += radial_momentum*area_weight;
    }
  }

  auto &w0 = pmbp->pmhd->w0;
  auto &bcc = pmbp->pmhd->bcc0;
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
      const Real radius2 = SQR(x)+SQR(y)+SQR(z);
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
      ComputeStressEnergy(x, y, z, parameters, density, pressure, velocity,
                          densitized_field, local_metric, four_velocity, stress,
                          sqrt_gamma, sqrt_minus_g);
      Real metric_derivative[3][4][4];
      DifferentiateSecondaryDisplacement(
          x, y, z, parameters, metric_derivative);
      const Real coordinate_volume = size.d_view(m).dx1*size.d_view(m).dx2
                                   *size.d_view(m).dx3;
      const Real radius = sqrt(radius2);
      const Real source_density = parameters.force_subtract_background
          ? density-parameters.rho : density;
      const Real newtonian_coefficient = parameters.metric.secondary_mass
          *source_density*sqrt_gamma*coordinate_volume/(radius2*radius);

      array_sum::GlobalSum cell;
      cell.the_array[7] = newtonian_coefficient*x;
      cell.the_array[8] = newtonian_coefficient*y;
      cell.the_array[9] = newtonian_coefficient*z;
      for (int force_direction=0; force_direction<3; ++force_direction) {
        Real contraction = 0.0;
        for (int a=0; a<4; ++a) {
          for (int b=0; b<4; ++b) {
            contraction += stress[a][b]
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
}

} // namespace
