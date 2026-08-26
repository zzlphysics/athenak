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
           << "emri-comoving-v1"
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

  emri_comoving::MetricParameters rotating{};
  rotating.coordinate_radius = 3.0;
  rotating.omega = 0.1;
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
  if (omega_is_geodesic && primary_geodesic_residual > geodesic_tolerance) {
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
  const Real speed = wind_tunnel.metric.omega*wind_tunnel.metric.coordinate_radius;
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
    wind_tunnel.refinement_hysteresis, wind_tunnel.refinement_horizon_factor
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
      || !(wind_tunnel.refinement_horizon_factor > 0.0)) {
    Fatal("invalid emri_windtunnel parameter");
  }
  const Real orbital_speed = std::abs(metric.omega*metric.coordinate_radius);
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

  const RegionSize &domain = pmy_mesh_->mesh_size;
  if (!(domain.x1min < 0.0 && domain.x1max > 0.0
        && domain.x2min < 0.0 && domain.x2max > 0.0
        && domain.x3min < 0.0 && domain.x3max > 0.0)) {
    Fatal("the local wind-tunnel domain must contain the secondary at the origin");
  }
  const Real primary_horizon = metric.primary_mass
      + std::sqrt(SQR(metric.primary_mass)-SQR(metric.primary_spin));
  const Real primary_enclosing = std::sqrt(
      SQR(metric.primary_spin)+SQR(primary_horizon));
  if (domain.x1min <= -metric.coordinate_radius+primary_enclosing) {
    Fatal("the local domain intersects the primary horizon; reduce its radial extent");
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
              << ", r/M=" << metric.orbital_radius/metric.primary_mass
              << ", Omega*M=" << metric.omega*metric.primary_mass
              << ", orbital KS speed=" << orbital_speed
              << ", r_H/m=" << hill_radius/metric.secondary_mass
              << ", box/orbit=" << largest_extent/metric.orbital_radius
              << ", metric_fd_step/m="
              << wind_tunnel.metric_fd_step/metric.secondary_mass
              << ", external_fd_step/M="
              << wind_tunnel.external_metric_fd_step/metric.primary_mass
              << ", M|Gamma^i_00|=" << primary_geodesic_residual
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
    const Real orbital_speed =
        parameters.metric.omega*parameters.metric.coordinate_radius;
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

void EMRIHistory(HistoryData *pdata, Mesh *) {
  pdata->nhist = 7;
  const char *labels[7] = {
    "primary_m", "secondary_m", "mass_ratio", "orbit_r_M", "omega_M", "hill_r_m",
    "geo_resid"
  };
  for (int n=0; n<7; ++n) pdata->label[n] = labels[n];
  const auto &metric = wind_tunnel.metric;
  const Real q = metric.secondary_mass/metric.primary_mass;
  const Real hill_radius = metric.orbital_radius*std::cbrt(q/3.0);
  const Real values[7] = {
    metric.primary_mass, metric.secondary_mass, q,
    metric.orbital_radius/metric.primary_mass,
    metric.omega*metric.primary_mass,
    hill_radius/metric.secondary_mass,
    primary_geodesic_residual
  };
  for (int n=0; n<7; ++n) {
    pdata->hdata[n] = (global_variable::my_rank == 0) ? values[n] : 0.0;
  }
}

} // namespace
