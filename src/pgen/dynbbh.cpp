//========================================================================================
// AthenaK astrophysical fluid dynamics and numerical relativity code
// Copyright(C) 2020 James M. Stone and the AthenaK collaboration
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file dynbbh.cpp
//! \brief GRMHD in the effective binary-black-hole spacetime of arXiv:2403.13308.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "athena.hpp"
#include "coordinates/adm.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/cartesian_ks.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/superposed_bbh.hpp"
#include "dyn_grmhd/dyn_grmhd.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "mhd/mhd.hpp"
#include "parameter_input.hpp"

namespace {

using binary_bh::AX1;
using binary_bh::AX2;
using binary_bh::AY1;
using binary_bh::AY2;
using binary_bh::AZ1;
using binary_bh::AZ2;
using binary_bh::M1T;
using binary_bh::M2T;
using binary_bh::NTRAJ;
using binary_bh::VX1;
using binary_bh::VX2;
using binary_bh::VY1;
using binary_bh::VY2;
using binary_bh::VZ1;
using binary_bh::VZ2;
using binary_bh::X1;
using binary_bh::X2;
using binary_bh::Y1;
using binary_bh::Y2;
using binary_bh::Z1;
using binary_bh::Z2;

struct BBHParameters {
  int use_table;
  Real separation;
  Real omega;
  Real mass_ratio;
  Real chi1;
  Real chi2;
  Real theta1;
  Real theta2;
  Real phi1;
  Real phi2;
  Real time_offset;
  Real fd_step;
  Real rho0;
  Real pgas0;
  Real u10;
  Real u20;
  Real u30;
  Real b10;
  Real b20;
  Real b30;
  Real alpha_threshold;
  Real refinement_radius;
  binary_bh::MetricParameters metric;
};

struct TrajectorySample {
  Real time;
  std::array<Real, NTRAJ> state;
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

struct TrajectoryStencil {
  Real minus[NTRAJ];
  Real center[NTRAJ];
  Real plus[NTRAJ];
  Real inverse_time_width;
};

BBHParameters bbh;
std::vector<TrajectorySample> trajectory_table;
Real bbh_metric_time = 0.0;

[[noreturn]] void Fatal(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << "\n" << message << std::endl;
  std::exit(EXIT_FAILURE);
}

void FindTrajectory(Real simulation_time, Real state[NTRAJ]);
void FindTrajectoryAtTableTime(Real table_time, Real state[NTRAJ]);
void TimeDerivativeEndpoints(Real center, Real &earlier, Real &later);
void LoadTrajectoryTable(const std::string &filename);
void ValidateTrajectorySignature(Real table_time_start, Real table_time_end);
void ValidateMetricKernel();
void SetADMVariablesToBBH(MeshBlockPack *pmbp);
void AugmentBBHExcisionMasks(MeshBlockPack *pmbp);
void RefineAlphaMin(MeshBlockPack *pmbp);
void RefineTracker(MeshBlockPack *pmbp);

KOKKOS_INLINE_FUNCTION
void EvaluateMetric(const Real x, const Real y, const Real z, const Real state[NTRAJ],
                    const BBHParameters &parameters, Real gcov[4][4]) {
  binary_bh::ComputeMetric(x, y, z, state, parameters.metric, gcov);
}

// The r=constant Kerr surfaces are oblate spheroids contained by a Euclidean sphere of
// radius sqrt(a^2+r^2) in the instantaneous rest frame.  A pure boost stretches any
// global spatial displacement by at most gamma, so this gives a conservative cell mask
// without assuming that Kerr r is a Euclidean distance.
KOKKOS_INLINE_FUNCTION
void CheckRegularizedHole(const Real x, const Real y, const Real z,
                          const Real position[3], const Real velocity[3],
                          const Real spin[3], const Real mass, const Real spin_buffer,
                          const Real global_padding,
                          const binary_bh::MetricParameters &metric_parameters,
                          bool &regularized, bool &near_regularized) {
  if (mass <= 0.0) return;

  Real x_bh[3];
  binary_bh::RestFramePosition(x, y, z, position, velocity, x_bh);
  const Real spin2 = SQR(spin[0]) + SQR(spin[1]) + SQR(spin[2]);
  const Real radius = binary_bh::RegularizationRadius(
      spin, mass, spin_buffer, metric_parameters.singularity_floor);
  const Real r2 = binary_bh::KerrRadiusSquared(x_bh, spin);
  regularized = regularized || (r2 <= SQR(radius));

  const Real velocity2 = SQR(velocity[0]) + SQR(velocity[1]) + SQR(velocity[2]);
  const Real gamma = 1.0/sqrt(1.0-velocity2);
  const Real rest_distance = sqrt(SQR(x_bh[0]) + SQR(x_bh[1]) + SQR(x_bh[2]));
  const Real enclosing_radius = sqrt(spin2 + SQR(radius));
  near_regularized = near_regularized
      || (rest_distance <= enclosing_radius + gamma*global_padding);
}

KOKKOS_INLINE_FUNCTION
void CheckRegularizedState(const Real x, const Real y, const Real z,
                           const Real global_padding, const Real state[NTRAJ],
                           const BBHParameters &parameters, bool &regularized,
                           bool &near_regularized) {
  const Real position1[3] = {state[X1], state[Y1], state[Z1]};
  const Real velocity1[3] = {state[VX1], state[VY1], state[VZ1]};
  const Real spin1[3] = {
    state[AX1]*parameters.metric.mass_scale1,
    state[AY1]*parameters.metric.mass_scale1,
    state[AZ1]*parameters.metric.mass_scale1
  };
  CheckRegularizedHole(x, y, z, position1, velocity1, spin1,
      state[M1T]*parameters.metric.mass_scale1, parameters.metric.spin_buffer1,
      global_padding, parameters.metric, regularized, near_regularized);

  const Real position2[3] = {state[X2], state[Y2], state[Z2]};
  const Real velocity2[3] = {state[VX2], state[VY2], state[VZ2]};
  const Real spin2[3] = {
    state[AX2]*parameters.metric.mass_scale2,
    state[AY2]*parameters.metric.mass_scale2,
    state[AZ2]*parameters.metric.mass_scale2
  };
  CheckRegularizedHole(x, y, z, position2, velocity2, spin2,
      state[M2T]*parameters.metric.mass_scale2, parameters.metric.spin_buffer2,
      global_padding, parameters.metric, regularized, near_regularized);
}

KOKKOS_INLINE_FUNCTION
void DifferentiateMetric(const Real x, const Real y, const Real z,
                         const TrajectoryStencil &trajectory,
                         const BBHParameters &parameters, MetricWithDerivatives &metric) {
  const Real step = parameters.fd_step;
  Real minus[4][4];
  Real plus[4][4];

  EvaluateMetric(x, y, z, trajectory.center, parameters, metric.g);

  EvaluateMetric(x, y, z, trajectory.minus, parameters, minus);
  EvaluateMetric(x, y, z, trajectory.plus, parameters, plus);
  for (int a = 0; a < 4; ++a) {
    for (int c = 0; c < 4; ++c) {
      metric.dg[0][a][c] =
          (plus[a][c] - minus[a][c])*trajectory.inverse_time_width;
    }
  }

  const Real xminus = x-step;
  const Real xplus = x+step;
  if (!(xplus > xminus)) {
    Kokkos::abort("metric_fd_step is too small for the x-coordinate precision");
  }
  EvaluateMetric(xminus, y, z, trajectory.center, parameters, minus);
  EvaluateMetric(xplus, y, z, trajectory.center, parameters, plus);
  for (int a = 0; a < 4; ++a) {
    for (int c = 0; c < 4; ++c) {
      metric.dg[1][a][c] = (plus[a][c] - minus[a][c])/(xplus-xminus);
    }
  }

  const Real yminus = y-step;
  const Real yplus = y+step;
  if (!(yplus > yminus)) {
    Kokkos::abort("metric_fd_step is too small for the y-coordinate precision");
  }
  EvaluateMetric(x, yminus, z, trajectory.center, parameters, minus);
  EvaluateMetric(x, yplus, z, trajectory.center, parameters, plus);
  for (int a = 0; a < 4; ++a) {
    for (int c = 0; c < 4; ++c) {
      metric.dg[2][a][c] = (plus[a][c] - minus[a][c])/(yplus-yminus);
    }
  }

  const Real zminus = z-step;
  const Real zplus = z+step;
  if (!(zplus > zminus)) {
    Kokkos::abort("metric_fd_step is too small for the z-coordinate precision");
  }
  EvaluateMetric(x, y, zminus, trajectory.center, parameters, minus);
  EvaluateMetric(x, y, zplus, trajectory.center, parameters, plus);
  for (int a = 0; a < 4; ++a) {
    for (int c = 0; c < 4; ++c) {
      metric.dg[3][a][c] = (plus[a][c] - minus[a][c])/(zplus-zminus);
    }
  }
}

// Convert g_ab and its Cartesian derivatives to the ADM variables.  The sign convention
// is dt(gamma_ij) = -2 alpha K_ij + D_i beta_j + D_j beta_i.
KOKKOS_INLINE_FUNCTION
void DecomposeMetric(const MetricWithDerivatives &metric, ADMPoint &result) {
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) result.gamma[i][j] = metric.g[i+1][j+1];
  }

  const Real determinant = adm::SpatialDet(
      result.gamma[0][0], result.gamma[0][1], result.gamma[0][2],
      result.gamma[1][1], result.gamma[1][2], result.gamma[2][2]);
  const Real leading_minor2 = result.gamma[0][0]*result.gamma[1][1]
      - SQR(result.gamma[0][1]);
  if (!isfinite(determinant) || !(result.gamma[0][0] > 0.0)
      || !(leading_minor2 > 0.0) || !(determinant > 0.0)) {
    Kokkos::abort("Non-positive-definite spatial metric in binary-BH metric");
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

  Real beta_lower[3] = {metric.g[0][1], metric.g[0][2], metric.g[0][3]};
  Real beta_squared = 0.0;
  for (int i = 0; i < 3; ++i) {
    result.beta[i] = 0.0;
    for (int j = 0; j < 3; ++j) result.beta[i] += inverse[i][j]*beta_lower[j];
    beta_squared += beta_lower[i]*result.beta[i];
  }
  const Real alpha_squared = beta_squared - metric.g[0][0];
  if (!isfinite(alpha_squared) || !(alpha_squared > 0.0)) {
    Kokkos::abort("Non-positive lapse squared in binary-BH metric");
  }
  result.alpha = sqrt(alpha_squared);
  result.psi4 = cbrt(determinant);

  for (int i = 0; i < 3; ++i) {
    for (int j = i; j < 3; ++j) {
      Real covariant_i_beta_j = metric.dg[i+1][0][j+1];
      Real covariant_j_beta_i = metric.dg[j+1][0][i+1];
      for (int k = 0; k < 3; ++k) {
        Real christoffel = 0.0;
        for (int l = 0; l < 3; ++l) {
          christoffel += 0.5*inverse[k][l]*(
              metric.dg[i+1][l+1][j+1] + metric.dg[j+1][l+1][i+1]
              - metric.dg[l+1][i+1][j+1]);
        }
        covariant_i_beta_j -= christoffel*beta_lower[k];
        covariant_j_beta_i -= christoffel*beta_lower[k];
      }
      result.curvature[i][j] = (covariant_i_beta_j + covariant_j_beta_i
          - metric.dg[0][i+1][j+1])/(2.0*result.alpha);
      result.curvature[j][i] = result.curvature[i][j];
    }
  }
}

} // namespace

//----------------------------------------------------------------------------------------
//! Initialize a uniform magnetized gas and the prescribed binary metric.

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
  if (!pmbp->pcoord->is_dynamical_relativistic || pmbp->padm == nullptr
      || pmbp->pdyngr == nullptr || pmbp->pmhd == nullptr) {
    Fatal("dynbbh requires <adm> and <mhd> blocks (dynamical GRMHD)");
  }
  if (!pmbp->padm->is_dynamic) Fatal("dynbbh requires <adm> dynamic = true");
  if (pin->DoesBlockExist("radiation")) {
    Fatal("dynbbh does not support radiation, whose tetrads assume a single Kerr metric");
  }
  for (const auto &input_block : pin->block) {
    if (input_block.block_name.compare(0, 6, "output") != 0) continue;
    for (const auto &input_line : input_block.line) {
      const bool is_variable = (input_line.param_name == "variable"
                                || input_line.param_name == "variable_2");
      if (is_variable && input_line.param_value == "mhd_jcon") {
        Fatal("dynbbh does not support mhd_jcon, which assumes a single Kerr metric");
      }
    }
  }
  if (!pmbp->pcoord->coord_data.bh_excise) {
    Fatal("dynbbh requires excision for the regularized moving black-hole interiors");
  }
  if (pmbp->pcoord->coord_data.excision_scheme != ExcisionScheme::lapse) {
    Fatal("a moving binary requires <coord> excision_scheme = lapse");
  }
  if (pmy_mesh_->multilevel && pmy_mesh_->pmr->prolong_prims) {
    Fatal("dynbbh AMR requires <mesh_refinement> prolong_primitives = false");
  }
  const std::string integrator = pin->GetOrAddString("time", "integrator", "rk2");
  if (integrator != "rk1" && integrator != "rk2" && integrator != "rk3") {
    Fatal("dynbbh currently supports the rk1, rk2, and rk3 time integrators");
  }

  std::string amr_condition = pin->GetOrAddString("problem", "amr_condition", "track");
  if (amr_condition == "alpha_min") {
    user_ref_func = RefineAlphaMin;
  } else if (amr_condition == "track") {
    user_ref_func = RefineTracker;
  } else {
    Fatal("unknown <problem> amr_condition: " + amr_condition);
  }
  pmbp->padm->SetADMVariables = &SetADMVariablesToBBH;
  pmbp->pcoord->AugmentExcisionMasks = &AugmentBBHExcisionMasks;

  const std::string trajectory_mode =
      pin->GetOrAddString("problem", "trajectory_mode", "circular");
  std::string trajectory_filename;
  if (trajectory_mode == "table") {
    bbh.use_table = 1;
    trajectory_filename = pin->GetString("problem", "trajectory_file");
  } else if (trajectory_mode == "circular") {
    bbh.use_table = 0;
    trajectory_filename.clear();
  } else {
    Fatal("unknown <problem> trajectory_mode: " + trajectory_mode);
  }

  bbh.separation = pin->GetOrAddReal("problem", "separation",
      pin->GetOrAddReal("problem", "sep", 20.0));
  bbh.mass_ratio = pin->GetOrAddReal("problem", "mass_ratio",
      pin->GetOrAddReal("problem", "q", 1.0));
  const Real default_omega = std::pow(bbh.separation, -1.5);
  bbh.omega = pin->GetOrAddReal("problem", "omega", default_omega);
  bbh.chi1 = pin->GetOrAddReal("problem", "chi1",
      pin->GetOrAddReal("problem", "a1", 0.0));
  bbh.chi2 = pin->GetOrAddReal("problem", "chi2",
      pin->GetOrAddReal("problem", "a2", 0.0));
  bbh.theta1 = pin->GetOrAddReal("problem", "theta1",
      pin->GetOrAddReal("problem", "th_a1", 0.0));
  bbh.theta2 = pin->GetOrAddReal("problem", "theta2",
      pin->GetOrAddReal("problem", "th_a2", 0.0));
  bbh.phi1 = pin->GetOrAddReal("problem", "phi1",
      pin->GetOrAddReal("problem", "ph_a1", 0.0));
  bbh.phi2 = pin->GetOrAddReal("problem", "phi2",
      pin->GetOrAddReal("problem", "ph_a2", 0.0));
  bbh.time_offset = pin->GetOrAddReal("problem", "trajectory_time_offset", 0.0);
#if SINGLE_PRECISION_ENABLED
  const Real default_fd_step = 5.0e-4;
#else
  const Real default_fd_step = 5.0e-5;
#endif
  bbh.fd_step = pin->GetOrAddReal("problem", "metric_fd_step", default_fd_step);
  bbh.rho0 = pin->GetOrAddReal("problem", "rho0", 1.0e-5);
  bbh.pgas0 = pin->GetOrAddReal("problem", "pgas0", 1.0e-7);
  bbh.u10 = pin->GetOrAddReal("problem", "u1", 0.0);
  bbh.u20 = pin->GetOrAddReal("problem", "u2", 0.0);
  bbh.u30 = pin->GetOrAddReal("problem", "u3", 0.0);
  bbh.b10 = pin->GetOrAddReal("problem", "b1", 0.0);
  bbh.b20 = pin->GetOrAddReal("problem", "b2", 0.0);
  bbh.b30 = pin->GetOrAddReal("problem", "b3", 0.0);
  bbh.alpha_threshold = pin->GetOrAddReal("problem", "alpha_threshold", 0.6);
  bbh.refinement_radius = pin->GetOrAddReal("problem", "refinement_radius", 6.0);
  bbh.metric.mass_scale1 = pin->GetOrAddReal("problem", "mass_scale1",
      pin->GetOrAddReal("problem", "adjust_mass1", 1.0));
  bbh.metric.mass_scale2 = pin->GetOrAddReal("problem", "mass_scale2",
      pin->GetOrAddReal("problem", "adjust_mass2", 1.0));
  bbh.metric.spin_buffer1 = pin->GetOrAddReal("problem", "spin_buffer1",
      pin->GetOrAddReal("problem", "a1_buffer", 0.05));
  bbh.metric.spin_buffer2 = pin->GetOrAddReal("problem", "spin_buffer2",
      pin->GetOrAddReal("problem", "a2_buffer", 0.05));
  bbh.metric.singularity_floor = pin->GetOrAddReal("problem", "singularity_floor",
      pin->GetOrAddReal("problem", "cutoff_floor", 1.0e-3));

  const Real finite_parameters[] = {
    bbh.separation, bbh.omega, bbh.mass_ratio, bbh.chi1, bbh.chi2,
    bbh.theta1, bbh.theta2, bbh.phi1, bbh.phi2, bbh.time_offset, bbh.fd_step,
    bbh.rho0, bbh.pgas0, bbh.u10, bbh.u20, bbh.u30, bbh.b10, bbh.b20, bbh.b30,
    bbh.alpha_threshold, bbh.refinement_radius, bbh.metric.mass_scale1,
    bbh.metric.mass_scale2, bbh.metric.spin_buffer1, bbh.metric.spin_buffer2,
    bbh.metric.singularity_floor
  };
  for (Real value : finite_parameters) {
    if (!std::isfinite(value)) Fatal("dynbbh parameters must be finite");
  }

  if (!(bbh.separation > 0.0) || !(bbh.mass_ratio > 0.0)
      || !(bbh.fd_step > 0.0) || !(bbh.rho0 > 0.0) || !(bbh.pgas0 > 0.0)
      || !(bbh.metric.mass_scale1 > 0.0) || !(bbh.metric.mass_scale2 > 0.0)
      || !(bbh.metric.spin_buffer1 >= 0.0) || !(bbh.metric.spin_buffer2 >= 0.0)
      || !(bbh.metric.singularity_floor > 0.0) || !(bbh.refinement_radius > 0.0)) {
    Fatal("invalid non-positive dynbbh parameter");
  }
  if (bbh.metric.spin_buffer1 + bbh.metric.singularity_floor > 1.0
      || bbh.metric.spin_buffer2 + bbh.metric.singularity_floor > 1.0) {
    Fatal("spin_buffer plus singularity_floor must not exceed one");
  }
  if (!bbh.use_table && (std::abs(bbh.chi1) > 1.0 || std::abs(bbh.chi2) > 1.0)) {
    Fatal("circular-orbit dimensionless spins chi1 and chi2 must satisfy |chi| <= 1");
  }
  if (!bbh.use_table) {
    const Real radius1 = bbh.mass_ratio*bbh.separation/(1.0 + bbh.mass_ratio);
    const Real radius2 = bbh.separation/(1.0 + bbh.mass_ratio);
    if (!(std::max(radius1, radius2)*std::abs(bbh.omega) < 1.0)) {
      Fatal("circular-orbit parameters give a superluminal black-hole velocity");
    }
  } else if (bbh.metric.mass_scale1 != bbh.metric.mass_scale2) {
    Fatal("table trajectories require equal mass scales to preserve the remnant limit");
  }

  Real required_start, unused_plus, unused_minus, required_end;
  TimeDerivativeEndpoints(pmy_mesh_->time + bbh.time_offset,
                          required_start, unused_plus);
  TimeDerivativeEndpoints(pin->GetReal("time", "tlim") + bbh.time_offset,
                          unused_minus, required_end);
  if (bbh.use_table) {
    LoadTrajectoryTable(trajectory_filename);
    if (required_start < trajectory_table.front().time
        || required_end > trajectory_table.back().time) {
      std::ostringstream message;
      message << "trajectory table covers [" << trajectory_table.front().time << ", "
              << trajectory_table.back().time
              << "] but the run and metric stencil require ["
              << required_start << ", " << required_end << "]";
      Fatal(message.str());
    }
  } else {
    trajectory_table.clear();
  }
  ValidateTrajectorySignature(required_start, required_end);

  if (global_variable::my_rank == 0) {
    std::cout << "Effective BBH metric: trajectory=" << trajectory_mode;
    if (bbh.use_table) std::cout << " file=" << trajectory_filename;
    std::cout << ", finite-difference step=" << bbh.fd_step << std::endl;
  }
  ValidateMetricKernel();

  // Metric data is not reliably reusable from a restart, so recompute it from the
  // trajectory.
  pmbp->padm->SetADMVariables(pmbp);
  if (pmbp->pcoord->coord_data.bh_excise) pmbp->pcoord->UpdateExcisionMasks();
  if (restart) return;

  auto &indcs = pmy_mesh_->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int nmb = pmbp->nmb_thispack;
  auto &w0 = pmbp->pmhd->w0;
  auto &b0 = pmbp->pmhd->b0;
  auto &bcc0 = pmbp->pmhd->bcc0;
  const int nscalars = pmbp->pmhd->nscalars;
  const BBHParameters parameters = bbh;

  par_for("dynbbh_uniform_gas", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    w0(m, IDN, k, j, i) = parameters.rho0;
    w0(m, IVX, k, j, i) = parameters.u10;
    w0(m, IVY, k, j, i) = parameters.u20;
    w0(m, IVZ, k, j, i) = parameters.u30;
    w0(m, IPR, k, j, i) = parameters.pgas0;
    for (int n = 0; n < nscalars; ++n) w0(m, IYF+n, k, j, i) = 0.0;
    bcc0(m, IBX, k, j, i) = parameters.b10;
    bcc0(m, IBY, k, j, i) = parameters.b20;
    bcc0(m, IBZ, k, j, i) = parameters.b30;
  });
  Kokkos::deep_copy(DevExeSpace(), b0.x1f, bbh.b10);
  Kokkos::deep_copy(DevExeSpace(), b0.x2f, bbh.b20);
  Kokkos::deep_copy(DevExeSpace(), b0.x3f, bbh.b30);

  pmbp->pdyngr->PrimToConInit(is, ie, js, je, ks, ke);
}

namespace {

void ValidateMetricKernel() {
#if SINGLE_PRECISION_ENABLED
  const Real tolerance = 3.0e-5;
#else
  const Real tolerance = 2.0e-12;
#endif
  const binary_bh::MetricParameters parameters = {1.0, 1.0, 0.05, 0.05, 1.0e-3};
  Real state[NTRAJ] = {0.0};
  state[AZ1] = 0.37;
  state[M1T] = 1.0;
  Real metric[4][4];
  binary_bh::ComputeMetric(2.1, -0.7, 0.9, state, parameters, metric);
  Real reference[4][4];
  Real inverse[4][4];
  ComputeMetricAndInverse(2.1, -0.7, 0.9, false, 0.37, reference, inverse);
  for (int a = 0; a < 4; ++a) {
    for (int c = 0; c < 4; ++c) {
      if (std::abs(metric[a][c]-reference[a][c]) > tolerance) {
        Fatal("single-hole limit of the binary metric failed its Kerr-Schild check");
      }
      if (std::abs(metric[a][c]-metric[c][a]) > tolerance) {
        Fatal("binary metric failed its symmetry check");
      }
    }
  }

  // The centered-derivative ADM decomposition must also recover the existing analytic
  // stationary Kerr lapse, shift, spatial metric, and extrinsic curvature.
  BBHParameters test_parameters{};
  test_parameters.fd_step = 5.0e-5;
  test_parameters.metric = parameters;
  TrajectoryStencil static_trajectory{};
  for (int n = 0; n < NTRAJ; ++n) {
    static_trajectory.minus[n] = state[n];
    static_trajectory.center[n] = state[n];
    static_trajectory.plus[n] = state[n];
  }
  static_trajectory.inverse_time_width = 0.5/test_parameters.fd_step;
  MetricWithDerivatives differentiated;
  ADMPoint decomposed;
  DifferentiateMetric(2.1, -0.7, 0.9, static_trajectory, test_parameters,
                      differentiated);
  DecomposeMetric(differentiated, decomposed);
  Real ref_alpha, ref_betax, ref_betay, ref_betaz, ref_psi4;
  Real ref_gxx, ref_gxy, ref_gxz, ref_gyy, ref_gyz, ref_gzz;
  Real ref_Kxx, ref_Kxy, ref_Kxz, ref_Kyy, ref_Kyz, ref_Kzz;
  ComputeADMDecomposition(2.1, -0.7, 0.9, false, 0.37,
      &ref_alpha, &ref_betax, &ref_betay, &ref_betaz, &ref_psi4,
      &ref_gxx, &ref_gxy, &ref_gxz, &ref_gyy, &ref_gyz, &ref_gzz,
      &ref_Kxx, &ref_Kxy, &ref_Kxz, &ref_Kyy, &ref_Kyz, &ref_Kzz);
  const Real adm_values[] = {
    decomposed.alpha, decomposed.beta[0], decomposed.beta[1], decomposed.beta[2],
    decomposed.psi4, decomposed.gamma[0][0], decomposed.gamma[0][1],
    decomposed.gamma[0][2], decomposed.gamma[1][1], decomposed.gamma[1][2],
    decomposed.gamma[2][2], decomposed.curvature[0][0], decomposed.curvature[0][1],
    decomposed.curvature[0][2], decomposed.curvature[1][1], decomposed.curvature[1][2],
    decomposed.curvature[2][2]
  };
  const Real adm_reference[] = {
    ref_alpha, ref_betax, ref_betay, ref_betaz, ref_psi4,
    ref_gxx, ref_gxy, ref_gxz, ref_gyy, ref_gyz, ref_gzz,
    ref_Kxx, ref_Kxy, ref_Kxz, ref_Kyy, ref_Kyz, ref_Kzz
  };
#if SINGLE_PRECISION_ENABLED
  const Real adm_tolerance = 3.0e-3;
#else
  const Real adm_tolerance = 2.0e-7;
#endif
  for (int n = 0; n < 17; ++n) {
    if (std::abs(adm_values[n]-adm_reference[n]) > adm_tolerance) {
      Fatal("single-hole ADM decomposition failed its analytic Kerr-Schild check");
    }
  }

  // The compact H taper must remain continuous across the spin equator at the Kerr ring.
  // This specifically guards against the sign-dependent coordinate jump used by older
  // generated implementations.
  for (int n = 0; n < NTRAJ; ++n) state[n] = 0.0;
  state[AZ1] = 0.45;
  state[M1T] = 0.5;
  Real above_ring[4][4], below_ring[4][4];
  binary_bh::ComputeMetric(0.45, 0.0, 1.0e-6, state, parameters, above_ring);
  binary_bh::ComputeMetric(0.45, 0.0, -1.0e-6, state, parameters, below_ring);
  for (int a = 0; a < 4; ++a) {
    for (int c = 0; c < 4; ++c) {
      if (std::abs(above_ring[a][c]-below_ring[a][c]) > 2.0e-2) {
        Fatal("binary metric regularization is discontinuous across the Kerr ring");
      }
    }
  }
  for (int n = 0; n < NTRAJ; ++n) {
    static_trajectory.minus[n] = state[n];
    static_trajectory.center[n] = state[n];
    static_trajectory.plus[n] = state[n];
  }
  DifferentiateMetric(0.45, 0.0, 0.0, static_trajectory, test_parameters,
                      differentiated);
  DecomposeMetric(differentiated, decomposed);
  if (!std::isfinite(decomposed.alpha) || !std::isfinite(decomposed.psi4)) {
    Fatal("binary metric regularization produced a non-finite ADM decomposition");
  }

  // Exercise the exact center, the inner taper boundary, and the transition back to the
  // unmodified Kerr metric along the spin axis.  All metric derivatives must remain
  // finite at these otherwise easy-to-miss locations.
  const Real transition_spin[3] = {state[AX1], state[AY1], state[AZ1]};
  const Real regularization_radius = binary_bh::RegularizationRadius(
      transition_spin, state[M1T], parameters.spin_buffer1,
      parameters.singularity_floor);
  const Real transition_points[] = {
    0.0, 0.5*regularization_radius, regularization_radius
  };
  for (Real z : transition_points) {
    DifferentiateMetric(0.0, 0.0, z, static_trajectory, test_parameters,
                        differentiated);
    DecomposeMetric(differentiated, decomposed);
    const Real transition_values[] = {
      decomposed.alpha, decomposed.psi4,
      decomposed.beta[0], decomposed.beta[1], decomposed.beta[2],
      decomposed.curvature[0][0], decomposed.curvature[0][1],
      decomposed.curvature[0][2], decomposed.curvature[1][1],
      decomposed.curvature[1][2], decomposed.curvature[2][2]
    };
    for (Real value : transition_values) {
      if (!std::isfinite(value)) {
        Fatal("binary metric taper produced a non-finite transition value");
      }
    }
  }

  // When both terms have identical position, velocity, and spin, linearity in H means
  // that their sum must be exactly the corresponding remnant Kerr metric.
  for (int n = 0; n < NTRAJ; ++n) state[n] = 0.0;
  state[X1] = state[X2] = 0.3;
  state[Y1] = state[Y2] = -0.2;
  state[Z1] = state[Z2] = 0.1;
  state[VX1] = state[VX2] = 0.08;
  state[VY1] = state[VY2] = -0.03;
  state[VZ1] = state[VZ2] = 0.01;
  state[AX1] = state[AX2] = 0.2;
  state[AY1] = state[AY2] = 0.1;
  state[AZ1] = state[AZ2] = -0.05;
  state[M1T] = 0.4;
  state[M2T] = 0.6;
  Real superposed[4][4];
  binary_bh::ComputeMetric(2.0, 1.0, -0.5, state, parameters, superposed);
  state[M1T] = 1.0;
  state[M2T] = 0.0;
  Real remnant[4][4];
  binary_bh::ComputeMetric(2.0, 1.0, -0.5, state, parameters, remnant);
  for (int a = 0; a < 4; ++a) {
    for (int c = 0; c < 4; ++c) {
      if (std::abs(superposed[a][c]-remnant[a][c]) > tolerance) {
        Fatal("binary metric failed its merged-remnant identity check");
      }
    }
  }
}

void LoadTrajectoryTable(const std::string &filename) {
  std::ifstream input(filename);
  if (!input.is_open()) Fatal("could not open trajectory table: " + filename);

  trajectory_table.clear();
  std::string line;
  int line_number = 0;
  while (std::getline(input, line)) {
    ++line_number;
    const std::size_t comment = line.find('#');
    if (comment != std::string::npos) line.erase(comment);
    std::istringstream row(line);
    TrajectorySample sample;
    if (!(row >> sample.time)) continue;
    for (int n = 0; n < NTRAJ; ++n) {
      if (!(row >> sample.state[n])) {
        Fatal("trajectory table line " + std::to_string(line_number)
            + " does not contain 21 numeric columns");
      }
    }
    Real extra;
    if (row >> extra) {
      Fatal("trajectory table line " + std::to_string(line_number)
          + " contains more than 21 numeric columns");
    }
    if (!std::isfinite(sample.time)) Fatal("non-finite trajectory time");
    for (int n = 0; n < NTRAJ; ++n) {
      if (!std::isfinite(sample.state[n])) Fatal("non-finite trajectory value");
    }
    if (!trajectory_table.empty() && !(sample.time > trajectory_table.back().time)) {
      Fatal("trajectory times must be strictly increasing");
    }
    if (!(sample.state[M1T] >= 0.0) || !(sample.state[M2T] >= 0.0)
        || !(sample.state[M1T] + sample.state[M2T] > 0.0)) {
      Fatal("trajectory component masses must be non-negative with positive total mass");
    }
    const Real speed1_squared = SQR(sample.state[VX1]) + SQR(sample.state[VY1])
        + SQR(sample.state[VZ1]);
    const Real speed2_squared = SQR(sample.state[VX2]) + SQR(sample.state[VY2])
        + SQR(sample.state[VZ2]);
    if (!(speed1_squared < 1.0) || !(speed2_squared < 1.0)) {
      Fatal("trajectory contains a superluminal black-hole velocity");
    }
    trajectory_table.push_back(sample);
  }
  if (trajectory_table.size() < 2) {
    Fatal("trajectory table must contain at least two rows");
  }

  // The Hermite velocity on each interval is a quadratic Bezier curve.  Its control
  // vectors are v0, 3*(p1-p0)/dt-v0-v1, and v1.  Since the open unit ball is convex,
  // requiring all three controls to be subluminal certifies every interpolated time,
  // including extrema between the representative signature-validation samples.
  const int position_index[2][3] = {{X1, Y1, Z1}, {X2, Y2, Z2}};
  const int velocity_index[2][3] = {{VX1, VY1, VZ1}, {VX2, VY2, VZ2}};
  for (std::size_t interval_index = 0;
       interval_index+1 < trajectory_table.size(); ++interval_index) {
    const TrajectorySample &lower = trajectory_table[interval_index];
    const TrajectorySample &upper = trajectory_table[interval_index+1];
    const Real interval = upper.time-lower.time;
    for (int hole = 0; hole < 2; ++hole) {
      Real control_speed_squared = 0.0;
      for (int direction = 0; direction < 3; ++direction) {
        const int position = position_index[hole][direction];
        const int velocity = velocity_index[hole][direction];
        const Real middle_control = 3.0*(upper.state[position]-lower.state[position])
            /interval-lower.state[velocity]-upper.state[velocity];
        control_speed_squared += SQR(middle_control);
      }
      if (!std::isfinite(control_speed_squared) || !(control_speed_squared < 1.0)) {
        std::ostringstream message;
        message << "trajectory interval [" << lower.time << ", " << upper.time
                << "] for hole " << hole+1
                << " cannot certify subluminal Hermite interpolation; middle velocity "
                << "control has |v|^2=" << control_speed_squared;
        Fatal(message.str());
      }
    }
  }
}

void FindTrajectory(Real simulation_time, Real state[NTRAJ]) {
  FindTrajectoryAtTableTime(simulation_time + bbh.time_offset, state);
}

void TimeDerivativeEndpoints(Real center, Real &earlier, Real &later) {
  if (!std::isfinite(center)) Fatal("metric stencil time must be finite");
  earlier = center-bbh.fd_step;
  later = center+bbh.fd_step;
  if (!(earlier < center)) {
    earlier = std::nextafter(center, -std::numeric_limits<Real>::infinity());
  }
  if (!(later > center)) {
    later = std::nextafter(center, std::numeric_limits<Real>::infinity());
  }
  if (!std::isfinite(earlier) || !std::isfinite(later) || !(earlier < later)) {
    Fatal("could not construct a finite time-derivative stencil");
  }
}

void FindTrajectoryAtTableTime(Real table_time, Real state[NTRAJ]) {
  if (!std::isfinite(table_time)) Fatal("requested trajectory time must be finite");
  if (!bbh.use_table) {
    const Real mass1 = 1.0/(1.0 + bbh.mass_ratio);
    const Real mass2 = bbh.mass_ratio*mass1;
    const Real radius1 = bbh.mass_ratio*bbh.separation/(1.0 + bbh.mass_ratio);
    const Real radius2 = -bbh.separation/(1.0 + bbh.mass_ratio);
    const Real phase = bbh.omega*table_time;
    const Real cosine = std::cos(phase);
    const Real sine = std::sin(phase);
    state[X1] = radius1*cosine;
    state[Y1] = radius1*sine;
    state[Z1] = 0.0;
    state[X2] = radius2*cosine;
    state[Y2] = radius2*sine;
    state[Z2] = 0.0;
    state[VX1] = -radius1*bbh.omega*sine;
    state[VY1] = radius1*bbh.omega*cosine;
    state[VZ1] = 0.0;
    state[VX2] = -radius2*bbh.omega*sine;
    state[VY2] = radius2*bbh.omega*cosine;
    state[VZ2] = 0.0;
    state[AX1] = bbh.chi1*mass1*std::sin(bbh.theta1)*std::cos(bbh.phi1);
    state[AY1] = bbh.chi1*mass1*std::sin(bbh.theta1)*std::sin(bbh.phi1);
    state[AZ1] = bbh.chi1*mass1*std::cos(bbh.theta1);
    state[AX2] = bbh.chi2*mass2*std::sin(bbh.theta2)*std::cos(bbh.phi2);
    state[AY2] = bbh.chi2*mass2*std::sin(bbh.theta2)*std::sin(bbh.phi2);
    state[AZ2] = bbh.chi2*mass2*std::cos(bbh.theta2);
    state[M1T] = mass1;
    state[M2T] = mass2;
    return;
  }

  const Real scale = std::max(Real(1.0), std::abs(table_time));
  const Real tolerance = 16.0*std::numeric_limits<Real>::epsilon()*scale;
  if (table_time < trajectory_table.front().time-tolerance
      || table_time > trajectory_table.back().time+tolerance) {
    Fatal("requested time lies outside the trajectory table");
  }
  if (table_time <= trajectory_table.front().time) {
    std::copy(trajectory_table.front().state.begin(),
              trajectory_table.front().state.end(), state);
    return;
  }
  if (table_time >= trajectory_table.back().time) {
    std::copy(trajectory_table.back().state.begin(),
              trajectory_table.back().state.end(), state);
    return;
  }

  const auto upper = std::lower_bound(
      trajectory_table.begin(), trajectory_table.end(), table_time,
      [](const TrajectorySample &sample, Real value) { return sample.time < value; });
  const auto lower = upper-1;
  const Real interval = upper->time-lower->time;
  const Real weight = (table_time-lower->time)/interval;
  for (int n = 0; n < NTRAJ; ++n) {
    state[n] = lower->state[n] + weight*(upper->state[n]-lower->state[n]);
  }

  // Cubic Hermite interpolation makes the supplied positions and velocities a
  // derivative-consistent pair.  Masses and rest-frame spin components remain linear,
  // matching the piecewise trajectory data without inventing their derivatives.
  const int position_index[2][3] = {{X1, Y1, Z1}, {X2, Y2, Z2}};
  const int velocity_index[2][3] = {{VX1, VY1, VZ1}, {VX2, VY2, VZ2}};
  const Real weight2 = weight*weight;
  const Real weight3 = weight2*weight;
  const Real h00 = 2.0*weight3 - 3.0*weight2 + 1.0;
  const Real h10 = weight3 - 2.0*weight2 + weight;
  const Real h01 = -2.0*weight3 + 3.0*weight2;
  const Real h11 = weight3 - weight2;
  for (int hole = 0; hole < 2; ++hole) {
    for (int direction = 0; direction < 3; ++direction) {
      const int position = position_index[hole][direction];
      const int velocity = velocity_index[hole][direction];
      const Real p0 = lower->state[position];
      const Real p1 = upper->state[position];
      const Real v0 = lower->state[velocity];
      const Real v1 = upper->state[velocity];
      state[position] = h00*p0 + h10*interval*v0
          + h01*p1 + h11*interval*v1;
      state[velocity] = ((6.0*weight2-6.0*weight)/interval)*p0
          + (3.0*weight2-4.0*weight+1.0)*v0
          + ((-6.0*weight2+6.0*weight)/interval)*p1
          + (3.0*weight2-2.0*weight)*v1;
    }
  }
  const Real speed1_squared = SQR(state[VX1]) + SQR(state[VY1]) + SQR(state[VZ1]);
  const Real speed2_squared = SQR(state[VX2]) + SQR(state[VY2]) + SQR(state[VZ2]);
  if (!(speed1_squared < 1.0) || !(speed2_squared < 1.0)) {
    Fatal("Hermite-interpolated trajectory contains a superluminal velocity");
  }
}

bool HasLorentzianADMDecomposition(const Real state[NTRAJ], const Real x,
                                   const Real y, const Real z,
                                   Real &alpha_squared) {
  Real metric[4][4];
  binary_bh::ComputeMetric(x, y, z, state, bbh.metric, metric);
  for (int a = 0; a < 4; ++a) {
    for (int c = 0; c < 4; ++c) {
      if (!std::isfinite(metric[a][c])) return false;
    }
  }

  const Real determinant = adm::SpatialDet(
      metric[1][1], metric[1][2], metric[1][3],
      metric[2][2], metric[2][3], metric[3][3]);
  const Real leading_minor2 = metric[1][1]*metric[2][2] - SQR(metric[1][2]);
  if (!(metric[1][1] > 0.0) || !(leading_minor2 > 0.0)
      || !(determinant > 0.0) || !std::isfinite(determinant)) {
    return false;
  }

  Real inverse[3][3];
  adm::SpatialInv(1.0/determinant,
      metric[1][1], metric[1][2], metric[1][3],
      metric[2][2], metric[2][3], metric[3][3],
      &inverse[0][0], &inverse[0][1], &inverse[0][2],
      &inverse[1][1], &inverse[1][2], &inverse[2][2]);
  inverse[1][0] = inverse[0][1];
  inverse[2][0] = inverse[0][2];
  inverse[2][1] = inverse[1][2];
  const Real beta_lower[3] = {metric[0][1], metric[0][2], metric[0][3]};
  Real beta_squared = 0.0;
  for (int i = 0; i < 3; ++i) {
    Real beta_upper = 0.0;
    for (int j = 0; j < 3; ++j) beta_upper += inverse[i][j]*beta_lower[j];
    beta_squared += beta_lower[i]*beta_upper;
  }
  alpha_squared = beta_squared-metric[0][0];
  const Real scale = std::max(Real(1.0),
      std::max(std::abs(beta_squared), std::abs(metric[0][0])));
  const Real margin = 64.0*std::numeric_limits<Real>::epsilon()*scale;
  return std::isfinite(alpha_squared) && alpha_squared > margin;
}

void GlobalPointAtKerrRadius(const Real position[3], const Real velocity[3],
                             const Real spin[3], const Real direction[3],
                             const Real radius, Real point[3]) {
  const Real direction_norm = sqrt(
      SQR(direction[0]) + SQR(direction[1]) + SQR(direction[2]));
  const Real unit_direction[3] = {
    direction[0]/direction_norm,
    direction[1]/direction_norm,
    direction[2]/direction_norm
  };
  const Real spin2 = SQR(spin[0]) + SQR(spin[1]) + SQR(spin[2]);
  Real mu = 0.0;
  if (spin2 > 0.0) {
    mu = (unit_direction[0]*spin[0] + unit_direction[1]*spin[1]
          + unit_direction[2]*spin[2])/sqrt(spin2);
  }
  const Real rest_length = radius*sqrt(
      (SQR(radius)+spin2)/(SQR(radius)+spin2*SQR(mu)));
  const Real rest_position[3] = {
    rest_length*unit_direction[0],
    rest_length*unit_direction[1],
    rest_length*unit_direction[2]
  };
  const Real velocity2 = SQR(velocity[0]) + SQR(velocity[1]) + SQR(velocity[2]);
  const Real gamma = 1.0/sqrt(1.0-velocity2);
  const Real vdotx = velocity[0]*rest_position[0]
      + velocity[1]*rest_position[1] + velocity[2]*rest_position[2];
  const Real inverse_coefficient = gamma/(gamma+1.0);
  for (int i = 0; i < 3; ++i) {
    point[i] = position[i] + rest_position[i]
        - inverse_coefficient*vdotx*velocity[i];
  }
}

void ValidateTrajectoryStateSignature(const Real table_time, const Real state[NTRAJ]) {
  auto CheckPoint = [&](const Real x, const Real y, const Real z) {
    Real alpha_squared = -std::numeric_limits<Real>::infinity();
    if (!HasLorentzianADMDecomposition(state, x, y, z, alpha_squared)) {
      std::ostringstream message;
      message << "effective BBH metric has no Lorentzian ADM decomposition at "
              << "trajectory time "
              << table_time << " near (" << x << ", " << y << ", " << z
              << "), alpha^2=" << alpha_squared
              << ". Start the inspiral-to-remnant transition earlier or increase the "
              << "circular smoke-test separation.";
      Fatal(message.str());
    }
  };

  const Real positions[2][3] = {
    {state[X1], state[Y1], state[Z1]},
    {state[X2], state[Y2], state[Z2]}
  };
  const Real velocities[2][3] = {
    {state[VX1], state[VY1], state[VZ1]},
    {state[VX2], state[VY2], state[VZ2]}
  };
  const Real spins[2][3] = {
    {state[AX1]*bbh.metric.mass_scale1,
     state[AY1]*bbh.metric.mass_scale1,
     state[AZ1]*bbh.metric.mass_scale1},
    {state[AX2]*bbh.metric.mass_scale2,
     state[AY2]*bbh.metric.mass_scale2,
     state[AZ2]*bbh.metric.mass_scale2}
  };
  const Real masses[2] = {
    state[M1T]*bbh.metric.mass_scale1,
    state[M2T]*bbh.metric.mass_scale2
  };

  if (masses[0] > 0.0 && masses[1] > 0.0) {
    for (int sample = 0; sample <= 64; ++sample) {
      const Real weight = static_cast<Real>(sample)/64.0;
      CheckPoint(positions[0][0] + weight*(positions[1][0]-positions[0][0]),
                 positions[0][1] + weight*(positions[1][1]-positions[0][1]),
                 positions[0][2] + weight*(positions[1][2]-positions[0][2]));
    }
  }

  for (int hole = 0; hole < 2; ++hole) {
    if (!(masses[hole] > 0.0)) continue;
    const int companion = 1-hole;
    const Real spin_buffer = (hole == 0) ?
        bbh.metric.spin_buffer1 : bbh.metric.spin_buffer2;
    const Real radius = binary_bh::RegularizationRadius(
        spins[hole], masses[hole], spin_buffer, bbh.metric.singularity_floor);
    Real rest_toward[3];
    binary_bh::RestFramePosition(
        positions[companion][0], positions[companion][1], positions[companion][2],
        positions[hole], velocities[hole], rest_toward);
    Real rest_separation = sqrt(
        SQR(rest_toward[0]) + SQR(rest_toward[1]) + SQR(rest_toward[2]));
    if (!(rest_separation > 0.0) || !(masses[companion] > 0.0)) {
      rest_toward[0] = 1.0;
      rest_toward[1] = rest_toward[2] = 0.0;
      rest_separation = 2.0*radius;
    }
    Real equatorial_toward[3] = {
      rest_toward[0], rest_toward[1], rest_toward[2]
    };
    const Real spin2 = SQR(spins[hole][0]) + SQR(spins[hole][1])
        + SQR(spins[hole][2]);
    if (spin2 > 0.0) {
      const Real projection = (rest_toward[0]*spins[hole][0]
          + rest_toward[1]*spins[hole][1]
          + rest_toward[2]*spins[hole][2])/spin2;
      for (int i = 0; i < 3; ++i) {
        equatorial_toward[i] -= projection*spins[hole][i];
      }
    }
    const Real equatorial_norm2 = SQR(equatorial_toward[0])
        + SQR(equatorial_toward[1]) + SQR(equatorial_toward[2]);
    const Real shell_factors[] = {0.5, 0.75, 0.95, 1.0, 1.05};
    constexpr int angular_samples = 256;
    const Real golden_angle = std::acos(-1.0)*(3.0-std::sqrt(5.0));
    for (Real shell_factor : shell_factors) {
      // A dense Fibonacci sphere avoids the large angular holes of a Cartesian 26-point
      // stencil for arbitrarily oriented spins.
      for (int sample = 0; sample < angular_samples; ++sample) {
        const Real z_direction = 1.0
            - 2.0*(static_cast<Real>(sample)+0.5)/angular_samples;
        const Real planar_radius = sqrt(fmax(1.0-SQR(z_direction), 0.0));
        const Real azimuth = golden_angle*sample;
        const Real direction[3] = {
          planar_radius*std::cos(azimuth),
          planar_radius*std::sin(azimuth),
          z_direction
        };
        Real point[3];
        GlobalPointAtKerrRadius(positions[hole], velocities[hole], spins[hole],
                                direction, shell_factor*radius, point);
        CheckPoint(point[0], point[1], point[2]);
      }
      const Real *special_directions[] = {rest_toward, equatorial_toward};
      const int direction_count = (equatorial_norm2 > 0.0) ? 2 : 1;
      for (int direction = 0; direction < direction_count; ++direction) {
        Real point[3];
        GlobalPointAtKerrRadius(positions[hole], velocities[hole], spins[hole],
                                special_directions[direction], shell_factor*radius,
                                point);
        CheckPoint(point[0], point[1], point[2]);
      }
    }

    Real ray_radius = 1.05*radius;
    for (int sample = 0;
         sample < 12 && ray_radius < 0.5*rest_separation; ++sample) {
      Real point[3];
      GlobalPointAtKerrRadius(positions[hole], velocities[hole], spins[hole],
                              rest_toward, ray_radius, point);
      CheckPoint(point[0], point[1], point[2]);
      ray_radius *= 2.0;
    }
  }
}

void ValidateTrajectorySignature(Real table_time_start, Real table_time_end) {
  if (!(table_time_end >= table_time_start)) {
    Fatal("metric-signature validation interval is reversed");
  }
  Real state[NTRAJ];
  if (!bbh.use_table) {
    const Real duration = table_time_end-table_time_start;
    Real validation_duration = duration;
    if (bbh.omega != 0.0) {
      validation_duration = std::min(duration,
          2.0*std::acos(-1.0)/std::abs(bbh.omega));
    }
    const int samples = (validation_duration > 0.0) ? 32 : 0;
    for (int sample = 0; sample <= samples; ++sample) {
      const Real weight = (samples > 0) ?
          static_cast<Real>(sample)/static_cast<Real>(samples) : 0.0;
      const Real time = table_time_start + weight*validation_duration;
      FindTrajectoryAtTableTime(time, state);
      ValidateTrajectoryStateSignature(time, state);
    }
    return;
  }

  FindTrajectoryAtTableTime(table_time_start, state);
  ValidateTrajectoryStateSignature(table_time_start, state);
  const auto first_upper = std::upper_bound(
      trajectory_table.begin(), trajectory_table.end(), table_time_start,
      [](Real value, const TrajectorySample &sample) { return value < sample.time; });
  std::size_t interval = (first_upper == trajectory_table.begin()) ? 0
      : static_cast<std::size_t>(first_upper-trajectory_table.begin()-1);
  while (interval+1 < trajectory_table.size()
         && trajectory_table[interval].time < table_time_end) {
    const Real segment_start = std::max(
        table_time_start, trajectory_table[interval].time);
    const Real segment_end = std::min(
        table_time_end, trajectory_table[interval+1].time);
    if (segment_end > segment_start) {
      const Real fractions[] = {0.25, 0.5, 0.75, 1.0};
      for (Real fraction : fractions) {
        const Real time = segment_start + fraction*(segment_end-segment_start);
        FindTrajectoryAtTableTime(time, state);
        ValidateTrajectoryStateSignature(time, state);
      }
    }
    ++interval;
  }
}

void SetADMVariablesToBBH(MeshBlockPack *pmbp) {
  const Real time = pmbp->pmesh->time;
  bbh_metric_time = time;
  TrajectoryStencil trajectory{};
  const Real table_time = time+bbh.time_offset;
  Real minus_time, plus_time;
  TimeDerivativeEndpoints(table_time, minus_time, plus_time);
  FindTrajectoryAtTableTime(minus_time, trajectory.minus);
  FindTrajectoryAtTableTime(table_time, trajectory.center);
  FindTrajectoryAtTableTime(plus_time, trajectory.plus);
  trajectory.inverse_time_width = 1.0/(plus_time-minus_time);

  auto &adm_vars = pmbp->padm->adm;
  auto &size = pmbp->pmb->mb_size;
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int ng = indcs.ng;
  const int is = indcs.is;
  const int js = indcs.js;
  const int ks = indcs.ks;
  const int n1 = indcs.nx1 + 2*ng;
  const int n2 = (indcs.nx2 > 1) ? indcs.nx2 + 2*ng : 1;
  const int n3 = (indcs.nx3 > 1) ? indcs.nx3 + 2*ng : 1;
  auto active_lids = pmbp->active_lids.d_view;
  const int active_offset = pmbp->active_offset;
  const int nmb_active = pmbp->nmb_active;
  const BBHParameters parameters = bbh;

  par_for_active("dynbbh_update_adm", DevExeSpace(), active_lids, active_offset,
  nmb_active, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min,
                               size.d_view(m).x1max);
    const Real y = CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min,
                               size.d_view(m).x2max);
    const Real z = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min,
                               size.d_view(m).x3max);
    MetricWithDerivatives metric;
    ADMPoint point;
    DifferentiateMetric(x, y, z, trajectory, parameters, metric);
    DecomposeMetric(metric, point);

    for (int a = 0; a < 3; ++a) {
      adm_vars.beta_u(m, a, k, j, i) = point.beta[a];
      for (int c = a; c < 3; ++c) {
        adm_vars.g_dd(m, a, c, k, j, i) = point.gamma[a][c];
        adm_vars.vK_dd(m, a, c, k, j, i) = point.curvature[a][c];
      }
    }
    adm_vars.alpha(m, k, j, i) = point.alpha;
    adm_vars.psi4(m, k, j, i) = point.psi4;
  });
}

void AugmentBBHExcisionMasks(MeshBlockPack *pmbp) {
  TrajectoryStencil trajectory{};
  const Real table_time = bbh_metric_time+bbh.time_offset;
  Real minus_time, plus_time;
  TimeDerivativeEndpoints(table_time, minus_time, plus_time);
  FindTrajectoryAtTableTime(minus_time, trajectory.minus);
  FindTrajectoryAtTableTime(table_time, trajectory.center);
  FindTrajectoryAtTableTime(plus_time, trajectory.plus);

  auto &floor = pmbp->pcoord->excision_floor;
  auto &flux = pmbp->pcoord->excision_flux;
  auto &size = pmbp->pmb->mb_size;
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int ng = indcs.ng;
  const int is = indcs.is;
  const int js = indcs.js;
  const int ks = indcs.ks;
  const int n1 = indcs.nx1 + 2*ng;
  const int n2 = (indcs.nx2 > 1) ? indcs.nx2 + 2*ng : 1;
  const int n3 = (indcs.nx3 > 1) ? indcs.nx3 + 2*ng : 1;
  auto active_lids = pmbp->active_lids.d_view;
  const int active_offset = pmbp->active_offset;
  const int nmb_active = pmbp->nmb_active;
  const BBHParameters parameters = bbh;

  par_for_active("dynbbh_geometric_excision", DevExeSpace(), active_lids, active_offset,
          nmb_active, 0, n3-1, 0, n2-1, 0, n1-1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x = CellCenterX(i-is, indcs.nx1, size.d_view(m).x1min,
                               size.d_view(m).x1max);
    const Real y = CellCenterX(j-js, indcs.nx2, size.d_view(m).x2min,
                               size.d_view(m).x2max);
    const Real z = CellCenterX(k-ks, indcs.nx3, size.d_view(m).x3min,
                               size.d_view(m).x3max);
    const Real dx2 = (indcs.nx2 > 1) ? SQR(size.d_view(m).dx2) : 0.0;
    const Real dx3 = (indcs.nx3 > 1) ? SQR(size.d_view(m).dx3) : 0.0;
    const Real global_padding = 0.5*sqrt(SQR(size.d_view(m).dx1) + dx2 + dx3)
        + parameters.fd_step;
    bool regularized = false;
    bool near_regularized = false;
    CheckRegularizedState(x, y, z, global_padding, trajectory.minus, parameters,
                          regularized, near_regularized);
    CheckRegularizedState(x, y, z, global_padding, trajectory.center, parameters,
                          regularized, near_regularized);
    CheckRegularizedState(x, y, z, global_padding, trajectory.plus, parameters,
                          regularized, near_regularized);

    floor(m, k, j, i) = floor(m, k, j, i) || regularized;
    flux(m, k, j, i) = flux(m, k, j, i) || near_regularized;
  });
}

void RefineAlphaMin(MeshBlockPack *pmbp) {
  Mesh *pmesh = pmbp->pmesh;
  const int nmb = pmbp->nmb_thispack;
  const int meshblock_start = pmesh->gids_eachrank[global_variable::my_rank];
  auto &refine_flag = pmesh->pmr->refine_flag;
  auto &indcs = pmesh->mb_indcs;
  const int is = indcs.is;
  const int js = indcs.js;
  const int ks = indcs.ks;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;
  const int cells_per_block = nx3*nx2*nx1;
  const int cells_per_plane = nx2*nx1;
  auto &alpha = pmbp->padm->adm.alpha;
  const Real threshold = bbh.alpha_threshold;

  par_for_outer("dynbbh_refine_lapse", DevExeSpace(), 0, 0, 0, nmb-1,
  KOKKOS_LAMBDA(TeamMember_t member, const int m) {
    Real minimum;
    Kokkos::parallel_reduce(Kokkos::TeamThreadRange(member, cells_per_block),
    [=](const int index, Real &value) {
      int k = index/cells_per_plane;
      int j = (index-k*cells_per_plane)/nx1;
      int i = index-k*cells_per_plane-j*nx1;
      value = fmin(value, alpha(m, k+ks, j+js, i+is));
    }, Kokkos::Min<Real>(minimum));
    if (minimum < threshold) refine_flag.d_view(m+meshblock_start) = 1;
    if (minimum > 1.25*threshold) refine_flag.d_view(m+meshblock_start) = -1;
  });
  refine_flag.template modify<DevExeSpace>();
  refine_flag.template sync<HostMemSpace>();
}

Real DistanceSquaredToBlock(Real x, Real y, Real z, const RegionSize &block) {
  const Real dx = std::max(std::max(block.x1min-x, Real(0.0)), x-block.x1max);
  const Real dy = std::max(std::max(block.x2min-y, Real(0.0)), y-block.x2max);
  const Real dz = std::max(std::max(block.x3min-z, Real(0.0)), z-block.x3max);
  return SQR(dx) + SQR(dy) + SQR(dz);
}

void RefineTracker(MeshBlockPack *pmbp) {
  Mesh *pmesh = pmbp->pmesh;
  auto &refine_flag = pmesh->pmr->refine_flag;
  auto &size = pmbp->pmb->mb_size;
  const int nmb = pmbp->nmb_thispack;
  const int meshblock_start = pmesh->gids_eachrank[global_variable::my_rank];
  Real state[NTRAJ];
  FindTrajectory(pmesh->time, state);
  for (int m = 0; m < nmb; ++m) {
    const RegionSize &block = size.h_view(m);
    Real minimum_distance = std::numeric_limits<Real>::infinity();
    if (state[M1T]*bbh.metric.mass_scale1 > 0.0) {
      minimum_distance = DistanceSquaredToBlock(
          state[X1], state[Y1], state[Z1], block);
    }
    if (state[M2T]*bbh.metric.mass_scale2 > 0.0) {
      minimum_distance = std::min(minimum_distance, DistanceSquaredToBlock(
          state[X2], state[Y2], state[Z2], block));
    }
    refine_flag.h_view(m+meshblock_start) =
        (minimum_distance < SQR(bbh.refinement_radius)) ? 1 : -1;
  }
  refine_flag.template modify<HostMemSpace>();
  refine_flag.template sync<DevExeSpace>();
}

} // namespace
