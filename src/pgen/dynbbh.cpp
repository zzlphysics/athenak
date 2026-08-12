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
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

#include "athena.hpp"
#include "coordinates/adm.hpp"
#include "coordinates/cell_locations.hpp"
#include "coordinates/cartesian_ks.hpp"
#include "coordinates/coordinates.hpp"
#include "coordinates/superposed_bbh.hpp"
#include "dyn_grmhd/dyn_grmhd.hpp"
#include "eos/eos.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "mhd/mhd.hpp"
#include "outputs/outputs.hpp"
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
  int initial_data;
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
  Real refinement_radius_ratio;
  Real refinement_hysteresis;
  Real refinement_horizon_factor;
  Real refinement_floor_radius;
  Real refinement_floor_center[3];
  int refinement_floor_level;
  Real history_inner_radius;
  Real trajectory_center_end_time;
  Real root_dt_max;
  binary_bh::MetricParameters metric;
};

enum class TorusMagneticField : int {
  none = 0,
  single_loop = 1
};

enum class TorusMagNorm : int {
  density_weighted_beta = 0,
  peak_beta = 1,
  integrated_pressure = 2,
  integrated_internal_energy = 3
};

const char *TorusMagNormName(const TorusMagNorm norm) {
  switch (norm) {
    case TorusMagNorm::density_weighted_beta: return "density_weighted_beta";
    case TorusMagNorm::peak_beta: return "peak_beta";
    case TorusMagNorm::integrated_pressure: return "integrated_pressure";
    case TorusMagNorm::integrated_internal_energy: return "integrated_internal_energy";
  }
  return "unknown";
}

struct ReferenceFMTorusParameters {
  Real gamma_adi;
  Real reference_mass;
  Real reference_center[3];
  Real reference_velocity[3];
  Real r_edge;
  Real r_peak;
  Real rho_max;
  Real rho_min;
  Real rho_pow;
  Real pgas_min;
  Real pgas_pow;
  Real rho_atm;
  Real temp_atm;
  Real pert_amp;
  Real potential_cutoff;
  Real mag_target;
  Real min_grid_peak_fraction;
  Real l_peak;
  Real log_h_edge;
  Real rho_normalization;
  Real r_outer;
  int seed;
  int min_magnetic_cells;
  bool require_full_domain;
  TorusMagneticField magnetic_field;
  TorusMagNorm mag_norm;
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

struct SignatureValidationStats {
  std::uint64_t points = 0;
  Real minimum_alpha_squared = std::numeric_limits<Real>::infinity();
  Real minimum_spatial_pivot = std::numeric_limits<Real>::infinity();
};

BBHParameters bbh;
ReferenceFMTorusParameters reference_torus;
std::vector<TrajectorySample> trajectory_table;
std::vector<Real> refinement_shell_radii;
std::string trajectory_content_fingerprint;  // NOLINT(runtime/string)
Real bbh_metric_time = 0.0;

[[noreturn]] void Fatal(const std::string &message) {
  std::cout << "### FATAL ERROR in " << __FILE__ << "\n" << message << std::endl;
  std::exit(EXIT_FAILURE);
}

void FindTrajectory(Real simulation_time, Real state[NTRAJ]);
void FindTrajectoryAtTableTime(Real table_time, Real state[NTRAJ]);
void TimeDerivativeEndpoints(Real center, Real &earlier, Real &later);
void LoadTrajectoryTable(const std::string &filename);
void ValidateOrStoreTrajectoryRestartContract(ParameterInput *pin, bool restart);
void ValidateTrajectorySignature(Real table_time_start, Real table_time_end);
void ValidateMetricKernel();
void SetADMVariablesToBBH(MeshBlockPack *pmbp);
void AugmentBBHExcisionMasks(MeshBlockPack *pmbp);
void RefineAlphaMin(MeshBlockPack *pmbp);
void RefineTracker(MeshBlockPack *pmbp);
void BBHHistory(HistoryData *pdata, Mesh *pm);
bool BBHOutputRegionCenter(const std::string &name, Real time, Real center[3]);
void InitializeReferenceFMTorus(MeshBlockPack *pmbp);

KOKKOS_INLINE_FUNCTION
Real ReferenceFMTorusLogH(const ReferenceFMTorusParameters &torus, Real r,
                          Real sin_theta);
KOKKOS_INLINE_FUNCTION
Real ReferenceFMTorusDensity(const ReferenceFMTorusParameters &torus, Real x, Real y,
                             Real z);
KOKKOS_INLINE_FUNCTION
Real ReferenceFMTorusA1(const ReferenceFMTorusParameters &torus, Real x, Real y, Real z);
KOKKOS_INLINE_FUNCTION
Real ReferenceFMTorusA2(const ReferenceFMTorusParameters &torus, Real x, Real y, Real z);

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

// The compact circumbinary-disk test in Combi & Ressler (2024) is initialized from a
// constant-angular-momentum Fishbone-Moncrief torus around a single, nonspinning hole of
// the total binary mass.  These helpers implement only that reference solution.  The
// velocity is mapped into the actual binary ADM metric below, so this is intentionally an
// approximate t=0 state rather than an equilibrium solution of the time-dependent metric.
KOKKOS_INLINE_FUNCTION
Real ReferenceFMTorusLFromRPeak(const Real reference_mass, const Real r_peak) {
  const Real dimensionless_r_peak = r_peak/reference_mass;
  return reference_mass*sqrt(dimensionless_r_peak)
      /(1.0 - 3.0/dimensionless_r_peak);
}

KOKKOS_INLINE_FUNCTION
Real ReferenceFMTorusLogH(const ReferenceFMTorusParameters &torus, const Real r,
                          const Real sin_theta) {
  const Real dimensionless_r = r/torus.reference_mass;
  if (!(dimensionless_r > 2.0) || !(fabs(sin_theta) > 0.0)) return -1.0e30;
  const Real sin2 = SQR(sin_theta);
  const Real delta = SQR(dimensionless_r) - 2.0*dimensionless_r;
  const Real sigma = SQR(dimensionless_r);
  const Real aa = SQR(SQR(dimensionless_r));
  const Real exp_2nu = sigma*delta/aa;
  const Real exp_2psi = aa/sigma*sin2;
  if (!(exp_2nu > 0.0) || !(exp_2psi > 0.0)) return -1.0e30;
  const Real dimensionless_l = torus.l_peak/torus.reference_mass;
  const Real var_a = sqrt(1.0 + 4.0*SQR(dimensionless_l)*exp_2nu/exp_2psi);
  return 0.5*log((1.0 + var_a)/exp_2nu) - 0.5*var_a;
}

KOKKOS_INLINE_FUNCTION
Real ReferenceFMTorusDensity(const ReferenceFMTorusParameters &torus, const Real x,
                             const Real y, const Real z) {
  const Real x_relative = x-torus.reference_center[0];
  const Real y_relative = y-torus.reference_center[1];
  const Real z_relative = z-torus.reference_center[2];
  const Real r = sqrt(SQR(x_relative) + SQR(y_relative) + SQR(z_relative));
  if (!(r >= torus.r_edge)) return 0.0;
  const Real cylindrical_radius = sqrt(SQR(x_relative) + SQR(y_relative));
  const Real log_h = ReferenceFMTorusLogH(torus, r, cylindrical_radius/r)
                   - torus.log_h_edge;
  if (!(log_h > 0.0)) return 0.0;
  const Real gm1 = torus.gamma_adi - 1.0;
  const Real p_over_rho = gm1/torus.gamma_adi*(exp(log_h) - 1.0);
  return pow(p_over_rho, 1.0/gm1)/torus.rho_normalization;
}

KOKKOS_INLINE_FUNCTION
Real ReferenceFMTorusAphi(const ReferenceFMTorusParameters &torus, const Real x,
                          const Real y, const Real z) {
  const Real rho = ReferenceFMTorusDensity(torus, x, y, z);
  return fmax(rho/torus.rho_max - torus.potential_cutoff, 0.0);
}

KOKKOS_INLINE_FUNCTION
Real ReferenceFMTorusA1(const ReferenceFMTorusParameters &torus, const Real x,
                        const Real y, const Real z) {
  const Real x_relative = x-torus.reference_center[0];
  const Real y_relative = y-torus.reference_center[1];
  const Real cylindrical_radius2 = SQR(x_relative) + SQR(y_relative);
  if (!(cylindrical_radius2 > 0.0)) return 0.0;
  return -y_relative/cylindrical_radius2*ReferenceFMTorusAphi(torus, x, y, z);
}

KOKKOS_INLINE_FUNCTION
Real ReferenceFMTorusA2(const ReferenceFMTorusParameters &torus, const Real x,
                        const Real y, const Real z) {
  const Real x_relative = x-torus.reference_center[0];
  const Real y_relative = y-torus.reference_center[1];
  const Real cylindrical_radius2 = SQR(x_relative) + SQR(y_relative);
  if (!(cylindrical_radius2 > 0.0)) return 0.0;
  return x_relative/cylindrical_radius2*ReferenceFMTorusAphi(torus, x, y, z);
}

KOKKOS_INLINE_FUNCTION
Real ReferenceFMTorusOmega(const ReferenceFMTorusParameters &torus, const Real r,
                           const Real sin_theta) {
  const Real dimensionless_r = r/torus.reference_mass;
  const Real dimensionless_l = torus.l_peak/torus.reference_mass;
  const Real sin2 = SQR(sin_theta);
  const Real delta = SQR(dimensionless_r) - 2.0*dimensionless_r;
  const Real sigma = SQR(dimensionless_r);
  const Real aa = SQR(SQR(dimensionless_r));
  const Real g_00 = -(1.0 - 2.0/dimensionless_r);
  const Real g_33 = sigma*sin2;
  const Real exp_2nu = sigma*delta/aa;
  const Real exp_2psi = aa/sigma*sin2;
  const Real u_phi_projected = sqrt(0.5*(-1.0 + sqrt(
      1.0 + 4.0*SQR(dimensionless_l)*exp_2nu/exp_2psi)));
  const Real u3 = sqrt(sigma/aa)/sin_theta*u_phi_projected;
  const Real u0 = -1.0/g_00*sqrt(-g_00*g_33*SQR(u3) - g_00);
  return u3/u0/torus.reference_mass;
}

KOKKOS_INLINE_FUNCTION
Real ReferenceFMTorusPerturbation(const int seed, const int gid, const int i,
                                  const int j, const int k) {
  // SplitMix64 gives a decomposition- and execution-order-independent cell perturbation.
  std::uint64_t value = static_cast<std::uint64_t>(static_cast<std::uint32_t>(seed));
  value ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(gid))
         * 0x9e3779b97f4a7c15ULL;
  value ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(i))
         * 0xbf58476d1ce4e5b9ULL;
  value ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(j))
         * 0x94d049bb133111ebULL;
  value ^= static_cast<std::uint64_t>(static_cast<std::uint32_t>(k))
         * 0xd6e8feb86659fd93ULL;
  value += 0x9e3779b97f4a7c15ULL;
  value = (value ^ (value >> 30))*0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27))*0x94d049bb133111ebULL;
  value ^= value >> 31;
  const Real unit = static_cast<Real>(value >> 11)*static_cast<Real>(0x1.0p-53);
  return 2.0*unit - 1.0;
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
    const std::string &block_name = input_block.block_name;
    const bool active = pin->DoesParameterExist(block_name, "dcycle")
                            ? pin->GetInteger(block_name, "dcycle") != 0
                            : pin->GetReal(block_name, "dt") > 0.0;
    if (!active) continue;
    const std::string file_type = pin->GetString(block_name, "file_type");
    const bool consumes_variables =
        (file_type != "hst" && file_type != "rst" && file_type != "log" &&
         file_type != "trk");
    if (!consumes_variables) continue;
    bool requests_jcon = pin->GetString(block_name, "variable") == "mhd_jcon";
    if (file_type == "pdf" && pin->DoesParameterExist(block_name, "variable_2") &&
        pin->GetOrAddInteger(block_name, "nbin2", 0) > 1) {
      requests_jcon = requests_jcon ||
                      pin->GetString(block_name, "variable_2") == "mhd_jcon";
    }
    if (requests_jcon) {
      Fatal("dynbbh does not support mhd_jcon, which assumes a single Kerr metric");
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
    if (pmy_mesh_->adaptive &&
        (pmy_mesh_->pmr->ncyc_check_amr != 1 ||
         pmy_mesh_->pmr->refinement_interval != 1)) {
      Fatal("moving-BBH track AMR requires ncycle_check=1 and refinement_interval=1");
    }
  } else {
    Fatal("unknown <problem> amr_condition: " + amr_condition);
  }
  pmbp->padm->SetADMVariables = &SetADMVariablesToBBH;
  pmbp->pcoord->AugmentExcisionMasks = &AugmentBBHExcisionMasks;
  user_hist_func = BBHHistory;
  user_output_region_func = BBHOutputRegionCenter;

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
  const std::string initial_data =
      pin->GetOrAddString("problem", "initial_data", "uniform");
  if (initial_data == "uniform") {
    bbh.initial_data = 0;
  } else if (initial_data == "reference_fm_torus") {
    bbh.initial_data = 1;
  } else {
    Fatal("unknown <problem> initial_data: " + initial_data);
  }
  const bool torus_r_edge_explicit =
      pin->DoesParameterExist("problem", "torus_r_edge");
  const bool torus_r_peak_explicit =
      pin->DoesParameterExist("problem", "torus_r_peak");
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
  bbh.refinement_radius_ratio =
      pin->GetOrAddReal("problem", "refinement_radius_ratio", 2.0);
  bbh.refinement_hysteresis =
      pin->GetOrAddReal("problem", "refinement_hysteresis", 1.25);
  bbh.refinement_horizon_factor =
      pin->GetOrAddReal("problem", "refinement_horizon_factor", 1.25);
  bbh.refinement_floor_radius =
      pin->GetOrAddReal("problem", "refinement_floor_radius", 0.0);
  bbh.refinement_floor_center[0] =
      pin->GetOrAddReal("problem", "refinement_floor_center1", 0.0);
  bbh.refinement_floor_center[1] =
      pin->GetOrAddReal("problem", "refinement_floor_center2", 0.0);
  bbh.refinement_floor_center[2] =
      pin->GetOrAddReal("problem", "refinement_floor_center3", 0.0);
  bbh.refinement_floor_level =
      pin->GetOrAddInteger("problem", "refinement_floor_level", 0);
  bbh.history_inner_radius = pin->GetOrAddReal(
      "problem", "history_inner_radius", 4.0*bbh.separation);
  bbh.root_dt_max = 0.0;
  if (pin->DoesParameterExist("time", "subcycling") &&
      pin->GetString("time", "subcycling") == "level" &&
      pin->DoesParameterExist("time", "root_dt_max")) {
    bbh.root_dt_max = pin->GetReal("time", "root_dt_max");
  }
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

  const int tracked_levels = pmy_mesh_->max_level-pmy_mesh_->root_level;
  refinement_shell_radii.clear();
  refinement_shell_radii.resize(tracked_levels+1, 0.0);
  bool has_explicit_shell_radius = false;
  for (int level=1; level<=tracked_levels; ++level) {
    const std::string name = "refinement_radius_level_" + std::to_string(level);
    has_explicit_shell_radius =
        has_explicit_shell_radius || pin->DoesParameterExist("problem", name);
  }
  for (int level=1; level<=tracked_levels; ++level) {
    const std::string name = "refinement_radius_level_" + std::to_string(level);
    if (has_explicit_shell_radius && !pin->DoesParameterExist("problem", name)) {
      Fatal("explicit track AMR requires every " + name
          + " through the configured finest level");
    }
    const Real default_radius = bbh.refinement_radius*std::pow(
        bbh.refinement_radius_ratio, tracked_levels-level);
    refinement_shell_radii[level] =
        pin->GetOrAddReal("problem", name, default_radius);
  }

  if (bbh.initial_data == 1) {
    if (pmy_mesh_->mb_indcs.nx2 <= 1 || pmy_mesh_->mb_indcs.nx3 <= 1) {
      Fatal("reference_fm_torus requires a three-dimensional mesh");
    }
    if (pin->GetString("mhd", "eos") != "ideal"
        || pin->GetString("mhd", "dyn_eos") != "ideal") {
      Fatal("reference_fm_torus currently requires mhd/eos=ideal and "
            "mhd/dyn_eos=ideal");
    }
    reference_torus.gamma_adi = pmbp->pmhd->peos->eos_data.gamma;
    reference_torus.r_edge = pin->GetOrAddReal("problem", "torus_r_edge", 18.0);
    reference_torus.r_peak = pin->GetOrAddReal("problem", "torus_r_peak", 29.0);
    reference_torus.rho_max = pin->GetOrAddReal("problem", "torus_rho_max", 1.0);
    reference_torus.rho_min = pin->GetOrAddReal("problem", "torus_rho_min", 1.0e-5);
    reference_torus.rho_pow = pin->GetOrAddReal("problem", "torus_rho_pow", -1.5);
    reference_torus.pgas_min =
        pin->GetOrAddReal("problem", "torus_pgas_min", 3.3333333333333334e-8);
    reference_torus.pgas_pow =
        pin->GetOrAddReal("problem", "torus_pgas_pow", -2.5);
    reference_torus.rho_atm =
        pin->GetOrAddReal("problem", "torus_rho_atm", 1.0e-10);
    if (pin->DoesParameterExist("problem", "torus_pgas_atm")) {
      Fatal("torus_pgas_atm has ambiguous pressure-floor semantics; use "
            "torus_temp_atm for the Valencia temperature floor");
    }
    reference_torus.temp_atm =
        pin->GetOrAddReal("problem", "torus_temp_atm", 3.3333333333333334e-13);
    reference_torus.pert_amp =
        pin->GetOrAddReal("problem", "torus_pert_amp", 0.02);
    reference_torus.seed = pin->GetOrAddInteger("problem", "torus_seed", 1);
    reference_torus.potential_cutoff =
        pin->GetOrAddReal("problem", "torus_potential_cutoff", 0.2);
    reference_torus.mag_target =
        pin->GetOrAddReal("problem", "torus_mag_target", 100.0);
    reference_torus.min_grid_peak_fraction = pin->GetOrAddReal(
        "problem", "torus_min_grid_peak_fraction", 0.9);
    reference_torus.min_magnetic_cells = pin->GetOrAddInteger(
        "problem", "torus_min_magnetic_cells", 64);
    reference_torus.require_full_domain = pin->GetOrAddBoolean(
        "problem", "torus_require_full_domain", true);

    const std::string magnetic_field =
        pin->GetOrAddString("problem", "torus_magnetic_field", "single_loop");
    if (magnetic_field == "single_loop") {
      reference_torus.magnetic_field = TorusMagneticField::single_loop;
    } else if (magnetic_field == "none") {
      reference_torus.magnetic_field = TorusMagneticField::none;
    } else {
      Fatal("unknown <problem> torus_magnetic_field: " + magnetic_field);
    }

    const std::string mag_norm = pin->GetOrAddString(
        "problem", "torus_mag_norm", "density_weighted_beta");
    if (mag_norm == "density_weighted_beta") {
      reference_torus.mag_norm = TorusMagNorm::density_weighted_beta;
    } else if (mag_norm == "peak_beta") {
      reference_torus.mag_norm = TorusMagNorm::peak_beta;
    } else if (mag_norm == "integrated_pressure") {
      reference_torus.mag_norm = TorusMagNorm::integrated_pressure;
    } else if (mag_norm == "integrated_internal_energy") {
      reference_torus.mag_norm = TorusMagNorm::integrated_internal_energy;
    } else {
      Fatal("unknown <problem> torus_mag_norm: " + mag_norm);
    }

    const Real torus_values[] = {
      reference_torus.gamma_adi, reference_torus.r_edge, reference_torus.r_peak,
      reference_torus.rho_max, reference_torus.rho_min, reference_torus.rho_pow,
      reference_torus.pgas_min, reference_torus.pgas_pow, reference_torus.rho_atm,
      reference_torus.temp_atm, reference_torus.pert_amp,
      reference_torus.potential_cutoff, reference_torus.mag_target,
      reference_torus.min_grid_peak_fraction
    };
    for (Real value : torus_values) {
      if (!std::isfinite(value)) Fatal("reference_fm_torus parameters must be finite");
    }
    if (!(reference_torus.gamma_adi > 1.0)
        || !(reference_torus.r_edge > 0.0)
        || !(reference_torus.r_peak > reference_torus.r_edge)
        || !(reference_torus.rho_max > 0.0)
        || !(reference_torus.rho_min >= 0.0)
        || !(reference_torus.pgas_min >= 0.0)
        || !(reference_torus.rho_atm > 0.0)
        || !(reference_torus.temp_atm > 0.0)
        || !(reference_torus.pert_amp >= 0.0 && reference_torus.pert_amp < 1.0)
        || !(reference_torus.potential_cutoff >= 0.0
             && reference_torus.potential_cutoff < 1.0)
        || !(reference_torus.mag_target > 0.0)
        || !(reference_torus.min_grid_peak_fraction >= 0.0
             && reference_torus.min_grid_peak_fraction <= 1.0)
        || !(reference_torus.min_magnetic_cells >= 0)) {
      Fatal("invalid reference_fm_torus parameter");
    }
    if (pin->GetReal("mhd", "dfloor") > reference_torus.rho_atm
        || pin->GetReal("mhd", "tfloor") > reference_torus.temp_atm) {
      Fatal("mhd density and temperature floors must not exceed the corresponding "
            "reference_fm_torus atmosphere floors");
    }
  }

  const Real finite_parameters[] = {
    bbh.separation, bbh.omega, bbh.mass_ratio, bbh.chi1, bbh.chi2,
    bbh.theta1, bbh.theta2, bbh.phi1, bbh.phi2, bbh.time_offset, bbh.fd_step,
    bbh.rho0, bbh.pgas0, bbh.u10, bbh.u20, bbh.u30, bbh.b10, bbh.b20, bbh.b30,
    bbh.alpha_threshold, bbh.refinement_radius, bbh.refinement_radius_ratio,
    bbh.refinement_hysteresis, bbh.refinement_horizon_factor,
    bbh.refinement_floor_radius, bbh.refinement_floor_center[0],
    bbh.refinement_floor_center[1], bbh.refinement_floor_center[2],
    bbh.history_inner_radius,
    bbh.metric.mass_scale1,
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
      || !(bbh.metric.singularity_floor > 0.0) || !(bbh.refinement_radius > 0.0)
      || !(bbh.refinement_radius_ratio >= 1.0)
      || !(bbh.refinement_hysteresis > 1.0)
      || !(bbh.refinement_horizon_factor > 0.0)
      || !(bbh.history_inner_radius > 0.0)) {
    Fatal("invalid non-positive dynbbh parameter");
  }
  if (!(bbh.refinement_floor_radius >= 0.0)
      || bbh.refinement_floor_level < 0
      || bbh.refinement_floor_level > tracked_levels
      || ((bbh.refinement_floor_radius == 0.0)
          != (bbh.refinement_floor_level == 0))) {
    Fatal("track AMR refinement floor requires radius>0 with a physical level in "
          "[1, num_levels-1], or radius=level=0 to disable it");
  }
  for (int level=1; level<=tracked_levels; ++level) {
    const Real radius = refinement_shell_radii[level];
    if (!std::isfinite(radius) || !(radius > 0.0)) {
      Fatal("track AMR refinement radii must be finite and positive");
    }
    if (level > 1 && radius > refinement_shell_radii[level-1]) {
      Fatal("track AMR refinement radii must be non-increasing toward finer levels");
    }
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
    // A segment tlim is only an execution/checkpoint boundary.  Moving-hole AMR
    // must look through that boundary when the immutable trajectory table has more
    // coverage, otherwise an ordinary segmented restart receives a light-speed
    // padding sphere instead of the actual future hole positions.  Reserve the
    // final finite-difference half stencil so every sampled endpoint is also a
    // valid metric-update time after restart.
    const Real center_end =
        trajectory_table.back().time-bbh.time_offset-bbh.fd_step;
    if (!std::isfinite(center_end)) {
      Fatal("trajectory table and time offset do not give a finite metric-time bound");
    }
    bbh.trajectory_center_end_time = std::nextafter(
        center_end, -std::numeric_limits<Real>::infinity());
  } else {
    trajectory_table.clear();
    trajectory_content_fingerprint = "analytic-circular";
    bbh.trajectory_center_end_time = std::numeric_limits<Real>::max();
  }
  ValidateOrStoreTrajectoryRestartContract(pin, restart);
  ValidateTrajectorySignature(required_start, required_end);

  Real inferred_torus_mass = 0.0;
  Real inferred_torus_center[3] = {0.0, 0.0, 0.0};
  Real inferred_torus_velocity[3] = {0.0, 0.0, 0.0};
  if (bbh.initial_data == 1) {
    Real initial_state[NTRAJ];
    FindTrajectory(pmy_mesh_->time, initial_state);
    const Real mass1 = initial_state[M1T]*bbh.metric.mass_scale1;
    const Real mass2 = initial_state[M2T]*bbh.metric.mass_scale2;
    inferred_torus_mass = mass1+mass2;
    if (!std::isfinite(inferred_torus_mass) || !(inferred_torus_mass > 0.0)) {
      Fatal("could not infer a positive reference_fm_torus mass from the trajectory");
    }
    inferred_torus_center[0] =
        (mass1*initial_state[X1] + mass2*initial_state[X2])/inferred_torus_mass;
    inferred_torus_center[1] =
        (mass1*initial_state[Y1] + mass2*initial_state[Y2])/inferred_torus_mass;
    inferred_torus_center[2] =
        (mass1*initial_state[Z1] + mass2*initial_state[Z2])/inferred_torus_mass;
    inferred_torus_velocity[0] =
        (mass1*initial_state[VX1] + mass2*initial_state[VX2])/inferred_torus_mass;
    inferred_torus_velocity[1] =
        (mass1*initial_state[VY1] + mass2*initial_state[VY2])/inferred_torus_mass;
    inferred_torus_velocity[2] =
        (mass1*initial_state[VZ1] + mass2*initial_state[VZ2])/inferred_torus_mass;

    reference_torus.reference_mass = pin->GetOrAddReal(
        "problem", "torus_reference_mass", inferred_torus_mass);
    reference_torus.reference_center[0] = pin->GetOrAddReal(
        "problem", "torus_reference_center1", inferred_torus_center[0]);
    reference_torus.reference_center[1] = pin->GetOrAddReal(
        "problem", "torus_reference_center2", inferred_torus_center[1]);
    reference_torus.reference_center[2] = pin->GetOrAddReal(
        "problem", "torus_reference_center3", inferred_torus_center[2]);
    reference_torus.reference_velocity[0] = pin->GetOrAddReal(
        "problem", "torus_reference_velocity1", inferred_torus_velocity[0]);
    reference_torus.reference_velocity[1] = pin->GetOrAddReal(
        "problem", "torus_reference_velocity2", inferred_torus_velocity[1]);
    reference_torus.reference_velocity[2] = pin->GetOrAddReal(
        "problem", "torus_reference_velocity3", inferred_torus_velocity[2]);
    if (!torus_r_edge_explicit) {
      reference_torus.r_edge = 18.0*reference_torus.reference_mass;
      pin->SetReal("problem", "torus_r_edge", reference_torus.r_edge);
    }
    if (!torus_r_peak_explicit) {
      reference_torus.r_peak = 29.0*reference_torus.reference_mass;
      pin->SetReal("problem", "torus_r_peak", reference_torus.r_peak);
    }
    const Real reference_values[] = {
      reference_torus.reference_mass,
      reference_torus.reference_center[0], reference_torus.reference_center[1],
      reference_torus.reference_center[2], reference_torus.reference_velocity[0],
      reference_torus.reference_velocity[1], reference_torus.reference_velocity[2]
    };
    for (Real value : reference_values) {
      if (!std::isfinite(value)) Fatal("reference_fm_torus frame must be finite");
    }
    const Real reference_velocity2 =
        SQR(reference_torus.reference_velocity[0])
        + SQR(reference_torus.reference_velocity[1])
        + SQR(reference_torus.reference_velocity[2]);
    if (!(reference_torus.reference_mass > 0.0)
        || !(reference_velocity2 < 1.0)
        || !(reference_torus.r_edge > 2.0*reference_torus.reference_mass)
        || !(reference_torus.r_peak > reference_torus.r_edge)) {
      Fatal("reference_fm_torus requires |v_ref|<1 and "
            "r_peak > r_edge > 2*M_ref");
    }

    reference_torus.l_peak = ReferenceFMTorusLFromRPeak(
        reference_torus.reference_mass, reference_torus.r_peak);
    reference_torus.log_h_edge = ReferenceFMTorusLogH(
        reference_torus, reference_torus.r_edge, 1.0);
    const Real log_h_peak = ReferenceFMTorusLogH(
        reference_torus, reference_torus.r_peak, 1.0) - reference_torus.log_h_edge;
    const Real gm1 = reference_torus.gamma_adi - 1.0;
    const Real p_over_rho_peak =
        gm1/reference_torus.gamma_adi*(std::exp(log_h_peak) - 1.0);
    reference_torus.rho_normalization =
        std::pow(p_over_rho_peak, 1.0/gm1)/reference_torus.rho_max;
    if (!std::isfinite(reference_torus.l_peak)
        || !std::isfinite(reference_torus.log_h_edge)
        || !(log_h_peak > 0.0)
        || !std::isfinite(reference_torus.rho_normalization)
        || !(reference_torus.rho_normalization > 0.0)) {
      Fatal("reference_fm_torus geometry does not define a finite torus");
    }

    Real outer_lower = reference_torus.r_peak;
    Real outer_upper = 2.0*reference_torus.r_peak;
    Real outer_surface = ReferenceFMTorusLogH(
        reference_torus, outer_upper, 1.0) - reference_torus.log_h_edge;
    for (int bracket=0; bracket<64 && outer_surface > 0.0; ++bracket) {
      outer_upper *= 2.0;
      outer_surface = ReferenceFMTorusLogH(
          reference_torus, outer_upper, 1.0) - reference_torus.log_h_edge;
    }
    if (!std::isfinite(outer_upper) || !(outer_surface <= 0.0)) {
      Fatal("could not bracket the reference_fm_torus outer edge");
    }
    for (int iteration=0; iteration<100; ++iteration) {
      const Real middle = 0.5*(outer_lower+outer_upper);
      const Real middle_surface = ReferenceFMTorusLogH(
          reference_torus, middle, 1.0) - reference_torus.log_h_edge;
      if (middle_surface > 0.0) {
        outer_lower = middle;
      } else {
        outer_upper = middle;
      }
    }
    reference_torus.r_outer = 0.5*(outer_lower+outer_upper);
    if (reference_torus.require_full_domain) {
      const RegionSize &domain = pmy_mesh_->mesh_size;
      if (reference_torus.reference_center[0]-reference_torus.r_outer < domain.x1min
          || reference_torus.reference_center[0]+reference_torus.r_outer > domain.x1max
          || reference_torus.reference_center[1]-reference_torus.r_outer < domain.x2min
          || reference_torus.reference_center[1]+reference_torus.r_outer > domain.x2max
          || reference_torus.reference_center[2]-reference_torus.r_outer < domain.x3min
          || reference_torus.reference_center[2]+reference_torus.r_outer > domain.x3max) {
        Fatal("reference_fm_torus is not fully contained in the mesh; enlarge the domain "
              "or explicitly set torus_require_full_domain=false");
      }
    }
  }

  if (global_variable::my_rank == 0) {
    std::cout << "Effective BBH metric: trajectory=" << trajectory_mode;
    if (bbh.use_table) {
      std::cout << " file=" << trajectory_filename
                << " fingerprint=" << trajectory_content_fingerprint;
    }
    std::cout << ", finite-difference step=" << bbh.fd_step
              << ", initial_data=" << initial_data << std::endl;
    if (!refinement_shell_radii.empty() && refinement_shell_radii.size() > 1) {
      std::cout << "BBH track AMR radii:";
      for (std::size_t level=1; level<refinement_shell_radii.size(); ++level) {
        std::cout << " L" << level << "=" << refinement_shell_radii[level];
      }
      std::cout << ", horizon factor=" << bbh.refinement_horizon_factor;
      if (bbh.root_dt_max > 0.0) {
        std::cout << ", root-step lookahead cap=" << bbh.root_dt_max;
      }
      std::cout << std::endl;
      if (bbh.refinement_floor_level > 0) {
        std::cout << "BBH AMR resolution floor: radius="
                  << bbh.refinement_floor_radius << " through L"
                  << bbh.refinement_floor_level << ", center=("
                  << bbh.refinement_floor_center[0] << ","
                  << bbh.refinement_floor_center[1] << ","
                  << bbh.refinement_floor_center[2] << ")" << std::endl;
      }
      if (user_hist) {
        std::cout << "BBH history: outside-excision global diagnostics, inner radius="
                  << bbh.history_inner_radius << std::endl;
      }
    }
    if (bbh.initial_data == 1) {
      std::cout << "Reference FM torus: M_ref=" << reference_torus.reference_mass
                << " (trajectory total=" << inferred_torus_mass << "), center=("
                << reference_torus.reference_center[0] << ","
                << reference_torus.reference_center[1] << ","
                << reference_torus.reference_center[2] << "), velocity=("
                << reference_torus.reference_velocity[0] << ","
                << reference_torus.reference_velocity[1] << ","
                << reference_torus.reference_velocity[2] << "), r_edge="
                << reference_torus.r_edge
                << ", r_peak=" << reference_torus.r_peak
                << ", r_outer=" << reference_torus.r_outer
                << ", rho_max=" << reference_torus.rho_max
                << ", magnetic norm=" << TorusMagNormName(reference_torus.mag_norm)
                << ", magnetic target=" << reference_torus.mag_target << std::endl;
    }
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

  if (bbh.initial_data == 0) {
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
  } else {
    InitializeReferenceFMTorus(pmbp);
  }

  pmbp->pdyngr->PrimToConInit(is, ie, js, je, ks, ke);
}

namespace {

void InitializeReferenceFMTorus(MeshBlockPack *pmbp) {
  auto &indcs = pmbp->pmesh->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const int nmb = pmbp->nmb_thispack;
  const int nscalars = pmbp->pmhd->nscalars;
  const int nx1 = indcs.nx1;
  const int nx2 = indcs.nx2;
  const int nx3 = indcs.nx3;
  auto &size = pmbp->pmb->mb_size;
  auto &mb_gid = pmbp->pmb->mb_gid;
  auto &adm = pmbp->padm->adm;
  auto &w0 = pmbp->pmhd->w0;
  auto &b0 = pmbp->pmhd->b0;
  auto &bcc0 = pmbp->pmhd->bcc0;
  auto &excision_floor = pmbp->pcoord->excision_floor;
  const Real dexcise = pmbp->pcoord->coord_data.dexcise;
  const Real pexcise = pmbp->pcoord->coord_data.pexcise;
  const ReferenceFMTorusParameters torus = reference_torus;

  Kokkos::deep_copy(DevExeSpace(), b0.x1f, 0.0);
  Kokkos::deep_copy(DevExeSpace(), b0.x2f, 0.0);
  Kokkos::deep_copy(DevExeSpace(), b0.x3f, 0.0);

  par_for("dynbbh_reference_fm_torus", DevExeSpace(), 0, nmb-1, ks, ke, js, je,
          is, ie, KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x = CellCenterX(i-is, nx1, size.d_view(m).x1min,
                              size.d_view(m).x1max);
    const Real y = CellCenterX(j-js, nx2, size.d_view(m).x2min,
                              size.d_view(m).x2max);
    const Real z = CellCenterX(k-ks, nx3, size.d_view(m).x3min,
                              size.d_view(m).x3max);
    const Real x_relative = x-torus.reference_center[0];
    const Real y_relative = y-torus.reference_center[1];
    const Real z_relative = z-torus.reference_center[2];
    const Real r = sqrt(SQR(x_relative) + SQR(y_relative) + SQR(z_relative));
    const Real r_background = fmax(r/torus.reference_mass, 1.0);
    const Real rho_background = fmax(torus.rho_atm,
        torus.rho_min*pow(r_background, torus.rho_pow));
    const Real pgas_background = fmax(rho_background*torus.temp_atm,
        torus.pgas_min*pow(r_background, torus.pgas_pow));
    const Real rho_torus = ReferenceFMTorusDensity(torus, x, y, z);

    const bool excised = excision_floor(m,k,j,i);
    Real rho = excised ? dexcise : rho_background;
    Real pgas = excised ? pexcise : pgas_background;
    Real utilde[3] = {0.0, 0.0, 0.0};
    if (!excised && rho_torus > 0.0) {
      const Real cylindrical_radius = sqrt(SQR(x_relative) + SQR(y_relative));
      const Real log_h = ReferenceFMTorusLogH(torus, r, cylindrical_radius/r)
                       - torus.log_h_edge;
      const Real p_over_rho = (torus.gamma_adi - 1.0)/torus.gamma_adi
                            *(exp(log_h) - 1.0);
      const Real perturbation = torus.pert_amp*ReferenceFMTorusPerturbation(
          torus.seed, mb_gid.d_view(m), i-is, j-js, k-ks);
      rho = fmax(rho_torus, rho_background);
      pgas = fmax(rho_torus*p_over_rho*(1.0 + perturbation), pgas_background);
      pgas = fmax(pgas, rho*torus.temp_atm);

      const Real omega = ReferenceFMTorusOmega(torus, r, cylindrical_radius/r);
      const Real coordinate_velocity[3] = {
        torus.reference_velocity[0]-omega*y_relative,
        torus.reference_velocity[1]+omega*x_relative,
        torus.reference_velocity[2]
      };
      Real normal_velocity[3];
      for (int a = 0; a < 3; ++a) {
        normal_velocity[a] = coordinate_velocity[a] + adm.beta_u(m,a,k,j,i);
      }
      Real normal_velocity2 = 0.0;
      for (int a = 0; a < 3; ++a) {
        for (int c = 0; c < 3; ++c) {
          normal_velocity2 += adm.g_dd(m,a,c,k,j,i)
                            *normal_velocity[a]*normal_velocity[c];
        }
      }
      const Real timelike_margin = SQR(adm.alpha(m,k,j,i)) - normal_velocity2;
      if (!isfinite(timelike_margin) || !(timelike_margin > 0.0)) {
        Kokkos::abort("reference_fm_torus velocity is not timelike in the BBH metric");
      }
      const Real u0 = 1.0/sqrt(timelike_margin);
      for (int a = 0; a < 3; ++a) utilde[a] = u0*normal_velocity[a];
    }

    w0(m, IDN, k, j, i) = rho;
    w0(m, IVX, k, j, i) = utilde[0];
    w0(m, IVY, k, j, i) = utilde[1];
    w0(m, IVZ, k, j, i) = utilde[2];
    w0(m, IPR, k, j, i) = pgas;
    for (int n = 0; n < nscalars; ++n) w0(m, IYF+n, k, j, i) = 0.0;
    bcc0(m, IBX, k, j, i) = 0.0;
    bcc0(m, IBY, k, j, i) = 0.0;
    bcc0(m, IBZ, k, j, i) = 0.0;
  });

  const int nkji = nx3*nx2*nx1;
  const int nji = nx2*nx1;
  Real grid_rho_max = 0.0;
  Real min_timelike_margin = std::numeric_limits<Real>::max();
  Kokkos::parallel_reduce("dynbbh_torus_checks",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
      KOKKOS_LAMBDA(const int idx, Real &rho_max, Real &margin_min) {
    const int m = idx/nkji;
    const int k0 = (idx-m*nkji)/nji;
    const int j0 = (idx-m*nkji-k0*nji)/nx1;
    const int i = idx-m*nkji-k0*nji-j0*nx1+is;
    const int k = k0+ks;
    const int j = j0+js;
    const Real x = CellCenterX(i-is, nx1, size.d_view(m).x1min,
                              size.d_view(m).x1max);
    const Real y = CellCenterX(j-js, nx2, size.d_view(m).x2min,
                              size.d_view(m).x2max);
    const Real z = CellCenterX(k-ks, nx3, size.d_view(m).x3min,
                              size.d_view(m).x3max);
    if (excision_floor(m,k,j,i)) return;
    const Real rho_torus = ReferenceFMTorusDensity(torus, x, y, z);
    rho_max = fmax(rho_max, rho_torus);
    if (rho_torus > 0.0) {
      const Real x_relative = x-torus.reference_center[0];
      const Real y_relative = y-torus.reference_center[1];
      const Real z_relative = z-torus.reference_center[2];
      const Real r = sqrt(SQR(x_relative) + SQR(y_relative) + SQR(z_relative));
      const Real omega = ReferenceFMTorusOmega(
          torus, r, sqrt(SQR(x_relative)+SQR(y_relative))/r);
      const Real coordinate_velocity[3] = {
        torus.reference_velocity[0]-omega*y_relative,
        torus.reference_velocity[1]+omega*x_relative,
        torus.reference_velocity[2]
      };
      Real normal_velocity[3];
      for (int a = 0; a < 3; ++a) {
        normal_velocity[a] = coordinate_velocity[a] + adm.beta_u(m,a,k,j,i);
      }
      Real normal_velocity2 = 0.0;
      for (int a = 0; a < 3; ++a) {
        for (int c = 0; c < 3; ++c) {
          normal_velocity2 += adm.g_dd(m,a,c,k,j,i)
                            *normal_velocity[a]*normal_velocity[c];
        }
      }
      margin_min = fmin(margin_min, SQR(adm.alpha(m,k,j,i))-normal_velocity2);
    }
  }, Kokkos::Max<Real>(grid_rho_max), Kokkos::Min<Real>(min_timelike_margin));

#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &grid_rho_max, 1, MPI_ATHENA_REAL, MPI_MAX,
                MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &min_timelike_margin, 1, MPI_ATHENA_REAL, MPI_MIN,
                MPI_COMM_WORLD);
#endif
  if (!(grid_rho_max > 0.0) || !std::isfinite(grid_rho_max)
      || !(min_timelike_margin > 0.0) || !std::isfinite(min_timelike_margin)) {
    Fatal("the mesh does not resolve a finite, timelike reference_fm_torus");
  }
  if (grid_rho_max < torus.min_grid_peak_fraction*torus.rho_max) {
    std::ostringstream message;
    message << "reference_fm_torus grid peak is only "
            << grid_rho_max/torus.rho_max << " of the analytic rho_max; required "
            << torus.min_grid_peak_fraction;
    Fatal(message.str());
  }

  if (torus.magnetic_field == TorusMagneticField::none) {
    if (global_variable::my_rank == 0) {
      std::cout << "Reference FM torus initialized: grid rho_max=" << grid_rho_max
                << ", minimum timelike margin=" << min_timelike_margin
                << ", magnetic field=none" << std::endl;
    }
    return;
  }

  std::int64_t magnetic_cells = 0;
  Kokkos::parallel_reduce("dynbbh_torus_magnetic_cells",
      Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
      KOKKOS_LAMBDA(const int idx, std::int64_t &count) {
    const int m = idx/nkji;
    const int k0 = (idx-m*nkji)/nji;
    const int j0 = (idx-m*nkji-k0*nji)/nx1;
    const int i = idx-m*nkji-k0*nji-j0*nx1+is;
    const int k = k0+ks;
    const int j = j0+js;
    if (excision_floor(m,k,j,i)) return;
    const Real x = CellCenterX(i-is, nx1, size.d_view(m).x1min,
                              size.d_view(m).x1max);
    const Real y = CellCenterX(j-js, nx2, size.d_view(m).x2min,
                              size.d_view(m).x2max);
    const Real z = CellCenterX(k-ks, nx3, size.d_view(m).x3min,
                              size.d_view(m).x3max);
    if (ReferenceFMTorusDensity(torus, x, y, z)
        > torus.potential_cutoff*torus.rho_max) {
      ++count;
    }
  }, Kokkos::Sum<std::int64_t>(magnetic_cells));
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &magnetic_cells, 1, MPI_INT64_T, MPI_SUM,
                MPI_COMM_WORLD);
#endif
  if (magnetic_cells < torus.min_magnetic_cells) {
    std::ostringstream message;
    message << "reference_fm_torus magnetic loop has only " << magnetic_cells
            << " active cells; required " << torus.min_magnetic_cells;
    Fatal(message.str());
  }

  const int ncells1 = nx1 + 2*indcs.ng;
  const int ncells2 = nx2 + 2*indcs.ng;
  const int ncells3 = nx3 + 2*indcs.ng;
  DvceArray4D<Real> a1("dynbbh_torus_a1", nmb, ncells3, ncells2, ncells1);
  DvceArray4D<Real> a2("dynbbh_torus_a2", nmb, ncells3, ncells2, ncells1);
  auto &nghbr = pmbp->pmb->nghbr;
  auto &mblev = pmbp->pmb->mb_lev;

  par_for("dynbbh_torus_vector_potential", DevExeSpace(), 0, nmb-1, ks, ke+1,
          js, je+1, is, ie+1, KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x1v = CellCenterX(i-is, nx1, size.d_view(m).x1min,
                                size.d_view(m).x1max);
    const Real x1f = LeftEdgeX(i-is, nx1, size.d_view(m).x1min,
                              size.d_view(m).x1max);
    const Real x2v = CellCenterX(j-js, nx2, size.d_view(m).x2min,
                                size.d_view(m).x2max);
    const Real x2f = LeftEdgeX(j-js, nx2, size.d_view(m).x2min,
                              size.d_view(m).x2max);
    const Real x3v = CellCenterX(k-ks, nx3, size.d_view(m).x3min,
                                size.d_view(m).x3max);
    const Real x3f = LeftEdgeX(k-ks, nx3, size.d_view(m).x3min,
                              size.d_view(m).x3max);
    const Real dx1 = size.d_view(m).dx1;
    const Real dx2 = size.d_view(m).dx2;

    a1(m,k,j,i) = ReferenceFMTorusA1(torus, x1v, x2f, x3f);
    a2(m,k,j,i) = ReferenceFMTorusA2(torus, x1f, x2v, x3f);
    if (EdgeTouchesFinerNeighbor(nghbr.d_view, m, mblev.d_view(m), 1,
        i, j, k, is, ie, js, je, ks, ke)) {
      a1(m,k,j,i) = 0.5*(ReferenceFMTorusA1(torus, x1v+0.25*dx1, x2f, x3f)
                         + ReferenceFMTorusA1(torus, x1v-0.25*dx1, x2f, x3f));
    }
    if (EdgeTouchesFinerNeighbor(nghbr.d_view, m, mblev.d_view(m), 2,
        i, j, k, is, ie, js, je, ks, ke)) {
      a2(m,k,j,i) = 0.5*(ReferenceFMTorusA2(torus, x1f, x2v+0.25*dx2, x3f)
                         + ReferenceFMTorusA2(torus, x1f, x2v-0.25*dx2, x3f));
    }
  });

  par_for("dynbbh_torus_curl_a", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
          KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real dx1 = size.d_view(m).dx1;
    const Real dx2 = size.d_view(m).dx2;
    const Real dx3 = size.d_view(m).dx3;
    b0.x1f(m,k,j,i) = -(a2(m,k+1,j,i)-a2(m,k,j,i))/dx3;
    b0.x2f(m,k,j,i) =  (a1(m,k+1,j,i)-a1(m,k,j,i))/dx3;
    b0.x3f(m,k,j,i) =  (a2(m,k,j,i+1)-a2(m,k,j,i))/dx1
                         -(a1(m,k,j+1,i)-a1(m,k,j,i))/dx2;
    if (i == ie) {
      b0.x1f(m,k,j,i+1) = -(a2(m,k+1,j,i+1)-a2(m,k,j,i+1))/dx3;
    }
    if (j == je) {
      b0.x2f(m,k,j+1,i) = (a1(m,k+1,j+1,i)-a1(m,k,j+1,i))/dx3;
    }
    if (k == ke) {
      b0.x3f(m,k+1,j,i) = (a2(m,k+1,j,i+1)-a2(m,k+1,j,i))/dx1
                            -(a1(m,k+1,j+1,i)-a1(m,k+1,j,i))/dx2;
    }
  });

  par_for("dynbbh_torus_cell_b", DevExeSpace(), 0, nmb-1, ks, ke, js, je, is, ie,
          KOKKOS_LAMBDA(int m, int k, int j, int i) {
    bcc0(m,IBX,k,j,i) = 0.5*(b0.x1f(m,k,j,i)+b0.x1f(m,k,j,i+1));
    bcc0(m,IBY,k,j,i) = 0.5*(b0.x2f(m,k,j,i)+b0.x2f(m,k,j+1,i));
    bcc0(m,IBZ,k,j,i) = 0.5*(b0.x3f(m,k,j,i)+b0.x3f(m,k+1,j,i));
  });

  Real gas_measure = 0.0;
  Real magnetic_measure = 0.0;
  if (torus.mag_norm == TorusMagNorm::peak_beta) {
    Kokkos::parallel_reduce("dynbbh_torus_peak_beta",
        Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb*nkji),
        KOKKOS_LAMBDA(const int idx, Real &gas_max, Real &magnetic_max) {
      const int m = idx/nkji;
      const int k0 = (idx-m*nkji)/nji;
      const int j0 = (idx-m*nkji-k0*nji)/nx1;
      const int i = idx-m*nkji-k0*nji-j0*nx1+is;
      const int k = k0+ks;
      const int j = j0+js;
      const Real x = CellCenterX(i-is, nx1, size.d_view(m).x1min,
                                size.d_view(m).x1max);
      const Real y = CellCenterX(j-js, nx2, size.d_view(m).x2min,
                                size.d_view(m).x2max);
      const Real z = CellCenterX(k-ks, nx3, size.d_view(m).x3min,
                                size.d_view(m).x3max);
      if (excision_floor(m,k,j,i)) return;
      if (!(ReferenceFMTorusDensity(torus, x, y, z) > 0.0)) return;
      const Real detg = adm::SpatialDet(
          adm.g_dd(m,0,0,k,j,i), adm.g_dd(m,0,1,k,j,i),
          adm.g_dd(m,0,2,k,j,i), adm.g_dd(m,1,1,k,j,i),
          adm.g_dd(m,1,2,k,j,i), adm.g_dd(m,2,2,k,j,i));
      const Real inv_sqrt_detg = 1.0/sqrt(detg);
      Real velocity_lower[3] = {0.0, 0.0, 0.0};
      Real velocity2 = 0.0;
      Real magnetic_lower[3] = {0.0, 0.0, 0.0};
      for (int a = 0; a < 3; ++a) {
        for (int c = 0; c < 3; ++c) {
          velocity_lower[a] += adm.g_dd(m,a,c,k,j,i)*w0(m,IVX+c,k,j,i);
          velocity2 += adm.g_dd(m,a,c,k,j,i)*w0(m,IVX+a,k,j,i)
                     *w0(m,IVX+c,k,j,i);
          magnetic_lower[a] += adm.g_dd(m,a,c,k,j,i)*bcc0(m,c,k,j,i)
                              *inv_sqrt_detg;
        }
      }
      Real magnetic_dot_velocity = 0.0;
      Real magnetic2 = 0.0;
      for (int a = 0; a < 3; ++a) {
        magnetic_dot_velocity += bcc0(m,a,k,j,i)*velocity_lower[a]*inv_sqrt_detg;
        magnetic2 += bcc0(m,a,k,j,i)*magnetic_lower[a]*inv_sqrt_detg;
      }
      const Real bsq = (magnetic2+SQR(magnetic_dot_velocity))/(1.0+velocity2);
      gas_max = fmax(gas_max, w0(m,IPR,k,j,i));
      magnetic_max = fmax(magnetic_max, 0.5*bsq);
    }, Kokkos::Max<Real>(gas_measure), Kokkos::Max<Real>(magnetic_measure));
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, &gas_measure, 1, MPI_ATHENA_REAL, MPI_MAX,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &magnetic_measure, 1, MPI_ATHENA_REAL, MPI_MAX,
                  MPI_COMM_WORLD);
#endif
  } else {
    // Preserve the beta normalization across MPI decompositions.  A single device-wide
    // floating-point reduction changes its association when the local MeshBlock count
    // changes.  Instead, form one deterministic compensated sum per global block,
    // gather those values in GID order, and add them once in that fixed order on rank 0.
    DvceArray2D<double> block_measure("dynbbh_torus_block_beta", nmb, 2);
    Kokkos::parallel_for("dynbbh_torus_integrated_beta_by_block",
        Kokkos::RangePolicy<>(DevExeSpace(), 0, nmb),
        KOKKOS_LAMBDA(const int m) {
      double gas_sum = 0.0;
      double magnetic_sum = 0.0;
      double gas_compensation = 0.0;
      double magnetic_compensation = 0.0;
      for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
          for (int i = is; i <= ie; ++i) {
            if (excision_floor(m,k,j,i)) continue;
            const Real x = CellCenterX(i-is, nx1, size.d_view(m).x1min,
                                      size.d_view(m).x1max);
            const Real y = CellCenterX(j-js, nx2, size.d_view(m).x2min,
                                      size.d_view(m).x2max);
            const Real z = CellCenterX(k-ks, nx3, size.d_view(m).x3min,
                                      size.d_view(m).x3max);
            if (!(ReferenceFMTorusDensity(torus, x, y, z) > 0.0)) continue;
            const Real detg = adm::SpatialDet(
                adm.g_dd(m,0,0,k,j,i), adm.g_dd(m,0,1,k,j,i),
                adm.g_dd(m,0,2,k,j,i), adm.g_dd(m,1,1,k,j,i),
                adm.g_dd(m,1,2,k,j,i), adm.g_dd(m,2,2,k,j,i));
            const Real sqrt_detg = sqrt(detg);
            const Real inv_sqrt_detg = 1.0/sqrt_detg;
            Real velocity_lower[3] = {0.0, 0.0, 0.0};
            Real velocity2 = 0.0;
            Real magnetic_lower[3] = {0.0, 0.0, 0.0};
            for (int a = 0; a < 3; ++a) {
              for (int c = 0; c < 3; ++c) {
                velocity_lower[a] += adm.g_dd(m,a,c,k,j,i)*w0(m,IVX+c,k,j,i);
                velocity2 += adm.g_dd(m,a,c,k,j,i)*w0(m,IVX+a,k,j,i)
                           *w0(m,IVX+c,k,j,i);
                magnetic_lower[a] += adm.g_dd(m,a,c,k,j,i)*bcc0(m,c,k,j,i)
                                    *inv_sqrt_detg;
              }
            }
            Real magnetic_dot_velocity = 0.0;
            Real magnetic2 = 0.0;
            for (int a = 0; a < 3; ++a) {
              magnetic_dot_velocity += bcc0(m,a,k,j,i)*velocity_lower[a]*inv_sqrt_detg;
              magnetic2 += bcc0(m,a,k,j,i)*magnetic_lower[a]*inv_sqrt_detg;
            }
            const Real bsq = (magnetic2+SQR(magnetic_dot_velocity))/(1.0+velocity2);
            const Real volume = size.d_view(m).dx1*size.d_view(m).dx2
                              * size.d_view(m).dx3*sqrt_detg;
            Real gas_weight = 1.0;
            if (torus.mag_norm == TorusMagNorm::density_weighted_beta) {
              gas_weight = w0(m,IDN,k,j,i);
            }
            Real gas_energy = w0(m,IPR,k,j,i);
            if (torus.mag_norm == TorusMagNorm::integrated_internal_energy) {
              gas_energy /= torus.gamma_adi-1.0;
            }
            const Real gas_cell_measure = gas_weight*gas_energy*volume;
            const double gas_value = static_cast<double>(gas_cell_measure);
            const double gas_increment = gas_value-gas_compensation;
            const double new_gas_sum = gas_sum+gas_increment;
            gas_compensation = (new_gas_sum-gas_sum)-gas_increment;
            gas_sum = new_gas_sum;
            const Real magnetic_cell_measure = gas_weight*0.5*bsq*volume;
            const double magnetic_value = static_cast<double>(magnetic_cell_measure);
            const double magnetic_increment = magnetic_value-magnetic_compensation;
            const double new_magnetic_sum = magnetic_sum+magnetic_increment;
            magnetic_compensation =
                (new_magnetic_sum-magnetic_sum)-magnetic_increment;
            magnetic_sum = new_magnetic_sum;
          }
        }
      }
      block_measure(m,0) = gas_sum;
      block_measure(m,1) = magnetic_sum;
    });

    auto block_measure_h =
        Kokkos::create_mirror_view_and_copy(HostMemSpace(), block_measure);
    std::vector<double> local_gas(nmb);
    std::vector<double> local_magnetic(nmb);
    for (int m = 0; m < nmb; ++m) {
      local_gas[m] = block_measure_h(m,0);
      local_magnetic[m] = block_measure_h(m,1);
    }

    const int nmb_total = pmbp->pmesh->nmb_total;
    int gid_order_valid =
        (nmb == pmbp->pmesh->nmb_eachrank[global_variable::my_rank]) ? 1 : 0;
    const int first_gid = pmbp->pmesh->gids_eachrank[global_variable::my_rank];
    for (int m = 0; m < nmb; ++m) {
      if (mb_gid.h_view(m) != first_gid+m) gid_order_valid = 0;
    }
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, &gid_order_valid, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#endif
    if (gid_order_valid == 0) {
      Fatal("reference_fm_torus deterministic magnetic normalization requires "
            "contiguous MeshBlock GIDs on every rank");
    }
    std::vector<double> gas_by_gid;
    std::vector<double> magnetic_by_gid;
    if (global_variable::my_rank == 0) {
      gas_by_gid.resize(nmb_total);
      magnetic_by_gid.resize(nmb_total);
    }
#if MPI_PARALLEL_ENABLED
    MPI_Gatherv(local_gas.data(), nmb, MPI_DOUBLE,
                global_variable::my_rank == 0 ? gas_by_gid.data() : nullptr,
                pmbp->pmesh->nmb_eachrank, pmbp->pmesh->gids_eachrank,
                MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gatherv(local_magnetic.data(), nmb, MPI_DOUBLE,
                global_variable::my_rank == 0 ? magnetic_by_gid.data() : nullptr,
                pmbp->pmesh->nmb_eachrank, pmbp->pmesh->gids_eachrank,
                MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
    gas_by_gid = local_gas;
    magnetic_by_gid = local_magnetic;
#endif
    if (global_variable::my_rank == 0) {
      double gas_sum = 0.0;
      double magnetic_sum = 0.0;
      double gas_compensation = 0.0;
      double magnetic_compensation = 0.0;
      for (int gid = 0; gid < nmb_total; ++gid) {
        const double gas_value = gas_by_gid[gid];
        const double new_gas_sum = gas_sum+gas_value;
        if (std::abs(gas_sum) >= std::abs(gas_value)) {
          gas_compensation += (gas_sum-new_gas_sum)+gas_value;
        } else {
          gas_compensation += (gas_value-new_gas_sum)+gas_sum;
        }
        gas_sum = new_gas_sum;
        const double magnetic_value = magnetic_by_gid[gid];
        const double new_magnetic_sum = magnetic_sum+magnetic_value;
        if (std::abs(magnetic_sum) >= std::abs(magnetic_value)) {
          magnetic_compensation +=
              (magnetic_sum-new_magnetic_sum)+magnetic_value;
        } else {
          magnetic_compensation +=
              (magnetic_value-new_magnetic_sum)+magnetic_sum;
        }
        magnetic_sum = new_magnetic_sum;
      }
      gas_measure = static_cast<Real>(gas_sum+gas_compensation);
      magnetic_measure = static_cast<Real>(magnetic_sum+magnetic_compensation);
    }
  }

  Real unnormalized_ratio = 0.0;
  Real magnetic_normalization = 0.0;
  int normalization_valid = 1;
  if (global_variable::my_rank == 0) {
    if (!std::isfinite(gas_measure) || !(gas_measure > 0.0)
        || !std::isfinite(magnetic_measure) || !(magnetic_measure > 0.0)) {
      normalization_valid = 0;
    } else {
      unnormalized_ratio = gas_measure/magnetic_measure;
      magnetic_normalization = sqrt(unnormalized_ratio/torus.mag_target);
      if (!std::isfinite(magnetic_normalization)
          || !(magnetic_normalization > 0.0)) {
        normalization_valid = 0;
      }
    }
  }
#if MPI_PARALLEL_ENABLED
  MPI_Bcast(&normalization_valid, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&magnetic_normalization, 1, MPI_ATHENA_REAL, 0, MPI_COMM_WORLD);
#endif
  if (normalization_valid == 0) {
    Fatal("reference_fm_torus magnetic loop has invalid normalization energy");
  }

  par_for("dynbbh_torus_normalize_face_b", DevExeSpace(), 0, nmb-1, ks, ke, js, je,
          is, ie, KOKKOS_LAMBDA(int m, int k, int j, int i) {
    b0.x1f(m,k,j,i) *= magnetic_normalization;
    b0.x2f(m,k,j,i) *= magnetic_normalization;
    b0.x3f(m,k,j,i) *= magnetic_normalization;
    if (i == ie) b0.x1f(m,k,j,i+1) *= magnetic_normalization;
    if (j == je) b0.x2f(m,k,j+1,i) *= magnetic_normalization;
    if (k == ke) b0.x3f(m,k+1,j,i) *= magnetic_normalization;
  });

  // Reconstruct only after every face has been normalized.  Combining these operations
  // in one kernel creates a RAW race: a cell can read its upper face while the adjacent
  // CUDA thread block is still scaling that face.  Besides making bcc nondeterministic,
  // PrimToConInit() would then construct conserved variables from a stale magnetic field.
  par_for("dynbbh_torus_normalized_cell_b", DevExeSpace(), 0, nmb-1, ks, ke, js, je,
          is, ie, KOKKOS_LAMBDA(int m, int k, int j, int i) {
    bcc0(m,IBX,k,j,i) = 0.5*(b0.x1f(m,k,j,i)+b0.x1f(m,k,j,i+1));
    bcc0(m,IBY,k,j,i) = 0.5*(b0.x2f(m,k,j,i)+b0.x2f(m,k,j+1,i));
    bcc0(m,IBZ,k,j,i) = 0.5*(b0.x3f(m,k,j,i)+b0.x3f(m,k+1,j,i));
  });

  if (global_variable::my_rank == 0) {
    std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10)
              << "Reference FM torus initialized: grid rho_max=" << grid_rho_max
              << ", minimum timelike margin=" << min_timelike_margin
              << ", magnetic cells=" << magnetic_cells
              << ", " << TorusMagNormName(torus.mag_norm)
              << " unnormalized ratio=" << unnormalized_ratio
              << ", magnetic scale=" << magnetic_normalization
              << ", normalized target=" << torus.mag_target
              << std::setprecision(6) << std::endl;
  }
}

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
    0.0, static_cast<Real>(0.5)*regularization_radius, regularization_radius
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

std::string FingerprintTrajectoryBytes(const std::string &contents) {
  std::uint64_t hash1 = 14695981039346656037ULL;
  std::uint64_t hash2 = 7809847782465536322ULL;
  for (char byte : contents) {
    const std::uint64_t value = static_cast<unsigned char>(byte);
    hash1 = (hash1 ^ value)*1099511628211ULL;
    hash2 ^= value + 0x9e3779b97f4a7c15ULL + (hash2 << 6) + (hash2 >> 2);
    hash2 *= 0xbf58476d1ce4e5b9ULL;
  }
  const std::uint64_t bytes = static_cast<std::uint64_t>(contents.size());

  std::uint64_t local[3] = {hash1, hash2, bytes};
#if MPI_PARALLEL_ENABLED
  std::uint64_t minimum[3];
  std::uint64_t maximum[3];
  MPI_Allreduce(local, minimum, 3, MPI_UINT64_T, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(local, maximum, 3, MPI_UINT64_T, MPI_MAX, MPI_COMM_WORLD);
  for (int n=0; n<3; ++n) {
    if (minimum[n] != maximum[n]) {
      Fatal("MPI ranks read different trajectory-table contents");
    }
  }
#endif

  std::ostringstream result;
  result << "content128-v1:" << std::hex << std::setfill('0')
         << std::setw(16) << local[0] << std::setw(16) << local[1]
         << std::dec << ":" << local[2];
  return result.str();
}

std::string CurrentTrajectoryRestartContract() {
  std::ostringstream contract;
  contract << std::setprecision(std::numeric_limits<Real>::max_digits10)
           << "dynbbh-trajectory-v1;mode=" << (bbh.use_table ? "table" : "circular")
           << ";content=" << trajectory_content_fingerprint
           << ";time_offset=" << bbh.time_offset
           << ";metric_fd_step=" << bbh.fd_step
           << ";mass_scale1=" << bbh.metric.mass_scale1
           << ";mass_scale2=" << bbh.metric.mass_scale2
           << ";spin_buffer1=" << bbh.metric.spin_buffer1
           << ";spin_buffer2=" << bbh.metric.spin_buffer2
           << ";singularity_floor=" << bbh.metric.singularity_floor;
  if (!bbh.use_table) {
    contract << ";separation=" << bbh.separation << ";omega=" << bbh.omega
             << ";mass_ratio=" << bbh.mass_ratio
             << ";chi1=" << bbh.chi1 << ";chi2=" << bbh.chi2
             << ";theta1=" << bbh.theta1 << ";theta2=" << bbh.theta2
             << ";phi1=" << bbh.phi1 << ";phi2=" << bbh.phi2;
  }
  return contract.str();
}

void ValidateOrStoreTrajectoryRestartContract(ParameterInput *pin, bool restart) {
  constexpr const char *contract_name = "trajectory_restart_contract";
  const std::string current = CurrentTrajectoryRestartContract();
  if (restart) {
    if (!pin->DoesParameterExist("problem", contract_name)) {
      Fatal("dynbbh restart lacks a trajectory contract; restart from a checkpoint "
            "written by the current code instead of silently accepting an unverified "
            "trajectory");
    }
    const std::string stored = pin->GetString("problem", contract_name);
    if (stored != current) {
      Fatal("dynbbh trajectory or metric parameters differ from the restart contract; "
            "the trajectory table, time offset, finite-difference step, mass scaling, "
            "and regularization parameters must remain unchanged");
    }
  } else {
    if (pin->DoesParameterExist("problem", contract_name)) {
      Fatal(std::string("<problem> ") + contract_name
          + " is reserved for restart provenance and must not be supplied by the user");
    }
    pin->SetString("problem", contract_name, current);
  }
}

void LoadTrajectoryTable(const std::string &filename) {
  // Read once, then fingerprint and parse the same immutable byte snapshot.  This keeps
  // the restart identity coupled to the states actually evolved even if an external
  // process replaces the path while startup is in progress.
  std::ifstream source(filename, std::ios::binary);
  if (!source.is_open()) Fatal("could not open trajectory table: " + filename);
  const std::string contents{
      std::istreambuf_iterator<char>(source), std::istreambuf_iterator<char>()};
  if (source.bad()) Fatal("failed while reading trajectory table: " + filename);
  trajectory_content_fingerprint = FingerprintTrajectoryBytes(contents);
  std::istringstream input(contents);

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
                                   Real &alpha_squared, Real &spatial_pivot) {
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
  spatial_pivot = metric[1][1];
  if (metric[1][1] > 0.0) {
    spatial_pivot = std::min(spatial_pivot, leading_minor2/metric[1][1]);
  }
  if (leading_minor2 > 0.0) {
    spatial_pivot = std::min(spatial_pivot, determinant/leading_minor2);
  }
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

void ValidateTrajectoryStateSignature(const Real table_time, const Real state[NTRAJ],
                                      SignatureValidationStats &stats) {
  auto CheckPoint = [&](const Real x, const Real y, const Real z) {
    Real alpha_squared = -std::numeric_limits<Real>::infinity();
    Real spatial_pivot = -std::numeric_limits<Real>::infinity();
    if (!HasLorentzianADMDecomposition(
            state, x, y, z, alpha_squared, spatial_pivot)) {
      std::ostringstream message;
      message << "effective BBH metric has no Lorentzian ADM decomposition at "
              << "trajectory time "
              << table_time << " near (" << x << ", " << y << ", " << z
              << "), alpha^2=" << alpha_squared
              << ", minimum spatial LDL pivot=" << spatial_pivot
              << ". Start the inspiral-to-remnant transition earlier or increase the "
              << "circular smoke-test separation.";
      Fatal(message.str());
    }
    ++stats.points;
    stats.minimum_alpha_squared =
        std::min(stats.minimum_alpha_squared, alpha_squared);
    stats.minimum_spatial_pivot =
        std::min(stats.minimum_spatial_pivot, spatial_pivot);
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
  SignatureValidationStats stats;
  Real state[NTRAJ];
  if (!bbh.use_table) {
    const Real duration = table_time_end-table_time_start;
    Real validation_duration = duration;
    if (bbh.omega != 0.0) {
      validation_duration = std::min(duration,
          static_cast<Real>(2.0*std::acos(-1.0))/std::abs(bbh.omega));
    }
    const int samples = (validation_duration > 0.0) ? 32 : 0;
    for (int sample = 0; sample <= samples; ++sample) {
      const Real weight = (samples > 0) ?
          static_cast<Real>(sample)/static_cast<Real>(samples) : 0.0;
      const Real time = table_time_start + weight*validation_duration;
      FindTrajectoryAtTableTime(time, state);
      ValidateTrajectoryStateSignature(time, state, stats);
    }
  } else {
    FindTrajectoryAtTableTime(table_time_start, state);
    ValidateTrajectoryStateSignature(table_time_start, state, stats);
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
          ValidateTrajectoryStateSignature(time, state, stats);
        }
      }
      ++interval;
    }
  }
  if (global_variable::my_rank == 0) {
    std::cout << std::setprecision(std::numeric_limits<Real>::max_digits10)
              << "Effective BBH signature certificate: interval=["
              << table_time_start << "," << table_time_end << "], points="
              << stats.points << ", min(alpha^2)=" << stats.minimum_alpha_squared
              << ", min(spatial LDL pivot)=" << stats.minimum_spatial_pivot
              << std::endl;
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
    int &flag = refine_flag.d_view(m+meshblock_start);
    if (minimum < threshold) flag = 1;
    if (minimum > 1.25*threshold && flag == 0) flag = -1;
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

struct TrackingTarget {
  Real position[3];
  Real velocity[3];
  Real spin[3];
  Real mass;
};

bool TrackingValuesCoincide(const Real first, const Real second) {
  const Real scale = std::max(Real(1.0), std::max(std::abs(first), std::abs(second)));
  return std::abs(first-second)
      <= 128.0*std::numeric_limits<Real>::epsilon()*scale;
}

void LoadTrackingTarget(const Real state[NTRAJ], const int hole,
                        TrackingTarget &target) {
  const int position_index[2][3] = {{X1, Y1, Z1}, {X2, Y2, Z2}};
  const int velocity_index[2][3] = {{VX1, VY1, VZ1}, {VX2, VY2, VZ2}};
  const int spin_index[2][3] = {{AX1, AY1, AZ1}, {AX2, AY2, AZ2}};
  const Real scale = (hole == 0) ? bbh.metric.mass_scale1 : bbh.metric.mass_scale2;
  for (int direction=0; direction<3; ++direction) {
    target.position[direction] = state[position_index[hole][direction]];
    target.velocity[direction] = state[velocity_index[hole][direction]];
    target.spin[direction] = state[spin_index[hole][direction]]*scale;
  }
  target.mass = state[(hole == 0) ? M1T : M2T]*scale;
}

int BuildTrackingTargets(const Real state[NTRAJ], TrackingTarget targets[2]) {
  LoadTrackingTarget(state, 0, targets[0]);
  LoadTrackingTarget(state, 1, targets[1]);
  const bool active1 = targets[0].mass > 0.0;
  const bool active2 = targets[1].mass > 0.0;
  if (!active1 && !active2) return 0;
  if (!active1) {
    targets[0] = targets[1];
    return 1;
  }
  if (!active2) return 1;

  bool coincident = true;
  for (int direction=0; direction<3; ++direction) {
    coincident = coincident
        && TrackingValuesCoincide(targets[0].position[direction],
                                  targets[1].position[direction])
        && TrackingValuesCoincide(targets[0].velocity[direction],
                                  targets[1].velocity[direction])
        && TrackingValuesCoincide(targets[0].spin[direction],
                                  targets[1].spin[direction]);
  }
  if (!coincident) return 2;

  targets[0].mass += targets[1].mass;
  for (int direction=0; direction<3; ++direction) {
    targets[0].position[direction] =
        0.5*(targets[0].position[direction]+targets[1].position[direction]);
    targets[0].velocity[direction] =
        0.5*(targets[0].velocity[direction]+targets[1].velocity[direction]);
    targets[0].spin[direction] =
        0.5*(targets[0].spin[direction]+targets[1].spin[direction]);
  }
  return 1;
}

//----------------------------------------------------------------------------------------
//! Resolve moving centers for filtered BBH binary output.

bool BBHOutputRegionCenter(const std::string &name, const Real time, Real center[3]) {
  Real state[NTRAJ];
  FindTrajectory(time, state);
  TrackingTarget holes[2];
  LoadTrackingTarget(state, 0, holes[0]);
  LoadTrackingTarget(state, 1, holes[1]);

  if (name == "bh1") {
    for (int axis=0; axis<3; ++axis) center[axis] = holes[0].position[axis];
    return true;
  }
  if (name == "bh2") {
    for (int axis=0; axis<3; ++axis) center[axis] = holes[1].position[axis];
    return true;
  }
  if (name == "bbh_com") {
    const Real total_mass = holes[0].mass+holes[1].mass;
    if (!(total_mass > 0.0)) return false;
    for (int axis=0; axis<3; ++axis) {
      center[axis] = (holes[0].mass*holes[0].position[axis]
                     +holes[1].mass*holes[1].position[axis])/total_mass;
    }
    return true;
  }
  return false;
}

//----------------------------------------------------------------------------------------
//! Global, metric-aware BBH diagnostics at synchronized root-level endpoints.
//!
//! The first ten entries are trajectory metadata.  HistoryOutput combines rank-local
//! buffers with MPI_SUM, so only rank zero contributes those replicated values.  The
//! remaining entries include coordinate-volume integrals of densitized conserved
//! variables, proper-volume integrals using sqrt(det(gamma_ij)), and two maxima.  Cells
//! in the excision-floor mask are deliberately omitted: their atmosphere state is
//! numerical regularization, not physical disk material.

void BBHHistory(HistoryData *pdata, Mesh *pm) {
  constexpr int ntrajectory = 10;
  pdata->nhist = 20;
  const char *labels[20] = {
    "bh1_x", "bh1_y", "bh1_z", "bh2_x", "bh2_y", "bh2_z",
    "bh_sep", "orb_omega", "bh1_mass", "bh2_mass",
    "baryon_m", "rho_prp", "pgas_prp", "emag_prp", "lor_D",
    "sigma_D", "angmom_z", "inner_D", "rho_max", "sigma_max"
  };
  for (int n=0; n<pdata->nhist; ++n) pdata->label[n] = labels[n];

  Real state[NTRAJ];
  FindTrajectory(pm->time, state);
  TrackingTarget holes[2];
  LoadTrackingTarget(state, 0, holes[0]);
  LoadTrackingTarget(state, 1, holes[1]);

  const Real rx = holes[0].position[0]-holes[1].position[0];
  const Real ry = holes[0].position[1]-holes[1].position[1];
  const Real rz = holes[0].position[2]-holes[1].position[2];
  const Real vx = holes[0].velocity[0]-holes[1].velocity[0];
  const Real vy = holes[0].velocity[1]-holes[1].velocity[1];
  const Real vz = holes[0].velocity[2]-holes[1].velocity[2];
  const Real separation2 = SQR(rx)+SQR(ry)+SQR(rz);
  const Real separation = sqrt(separation2);
  const Real omega_x = ry*vz-rz*vy;
  const Real omega_y = rz*vx-rx*vz;
  const Real omega_z = rx*vy-ry*vx;
  const Real orbital_omega = (separation2 > 0.0)
      ? sqrt(SQR(omega_x)+SQR(omega_y)+SQR(omega_z))/separation2 : 0.0;
  const Real trajectory_values[ntrajectory] = {
    holes[0].position[0], holes[0].position[1], holes[0].position[2],
    holes[1].position[0], holes[1].position[1], holes[1].position[2],
    separation, orbital_omega, holes[0].mass, holes[1].mass
  };
  for (int n=0; n<ntrajectory; ++n) {
    pdata->hdata[n] = (global_variable::my_rank == 0) ? trajectory_values[n] : 0.0;
  }

  const Real total_mass = holes[0].mass+holes[1].mass;
  Real center_x = 0.0;
  Real center_y = 0.0;
  Real center_z = 0.0;
  if (total_mass > 0.0) {
    center_x = (holes[0].mass*holes[0].position[0]
               +holes[1].mass*holes[1].position[0])/total_mass;
    center_y = (holes[0].mass*holes[0].position[1]
               +holes[1].mass*holes[1].position[1])/total_mass;
    center_z = (holes[0].mass*holes[0].position[2]
               +holes[1].mass*holes[1].position[2])/total_mass;
  }

  MeshBlockPack *pmbp = pm->pmb_pack;
  auto &u0 = pmbp->pmhd->u0;
  auto &w0 = pmbp->pmhd->w0;
  auto &bcc = pmbp->pmhd->bcc0;
  auto &adm = pmbp->padm->adm;
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
  const Real inner_radius2 = SQR(bbh.history_inner_radius);

  array_sum::GlobalSum integral;
  Real rho_max = 0.0;
  Real sigma_max = 0.0;
  Kokkos::parallel_reduce(
    "dynbbh_history", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
    KOKKOS_LAMBDA(const int idx, array_sum::GlobalSum &sum, Real &max_rho,
                  Real &max_sigma) {
      const int m = idx/nkji;
      const int k0 = (idx-m*nkji)/nji;
      const int j0 = (idx-m*nkji-k0*nji)/nx1;
      const int i0 = idx-m*nkji-k0*nji-j0*nx1;
      const int i = i0+is;
      const int j = j0+js;
      const int k = k0+ks;
      if (excision_floor(m,k,j,i)) return;

      const Real gxx = adm.g_dd(m,0,0,k,j,i);
      const Real gxy = adm.g_dd(m,0,1,k,j,i);
      const Real gxz = adm.g_dd(m,0,2,k,j,i);
      const Real gyy = adm.g_dd(m,1,1,k,j,i);
      const Real gyz = adm.g_dd(m,1,2,k,j,i);
      const Real gzz = adm.g_dd(m,2,2,k,j,i);
      const Real detg = adm::SpatialDet(gxx, gxy, gxz, gyy, gyz, gzz);
      const Real sqrt_detg = sqrt(detg);
      const Real inv_sqrt_detg = 1.0/sqrt_detg;

      const Real wx = w0(m,IVX,k,j,i);
      const Real wy = w0(m,IVY,k,j,i);
      const Real wz = w0(m,IVZ,k,j,i);
      const Real wx_d = gxx*wx+gxy*wy+gxz*wz;
      const Real wy_d = gxy*wx+gyy*wy+gyz*wz;
      const Real wz_d = gxz*wx+gyz*wy+gzz*wz;
      const Real lorentz2 = 1.0+wx*wx_d+wy*wy_d+wz*wz_d;
      const Real lorentz = sqrt(lorentz2);

      const Real bx = bcc(m,IBX,k,j,i)*inv_sqrt_detg;
      const Real by = bcc(m,IBY,k,j,i)*inv_sqrt_detg;
      const Real bz = bcc(m,IBZ,k,j,i)*inv_sqrt_detg;
      const Real bx_d = gxx*bx+gxy*by+gxz*bz;
      const Real by_d = gxy*bx+gyy*by+gyz*bz;
      const Real bz_d = gxz*bx+gyz*by+gzz*bz;
      const Real magnetic2 = bx*bx_d+by*by_d+bz*bz_d;
      const Real magnetic_velocity = (bx*wx_d+by*wy_d+bz*wz_d)/lorentz;
      const Real bsq = SQR(magnetic_velocity)+magnetic2/lorentz2;
      const Real rho = w0(m,IDN,k,j,i);
      const Real pgas = w0(m,IPR,k,j,i);
      const Real sigma = (rho > 0.0) ? bsq/rho : 0.0;
      const Real coordinate_volume = size.d_view(m).dx1*size.d_view(m).dx2
                                   *size.d_view(m).dx3;
      const Real proper_volume = coordinate_volume*sqrt_detg;
      const Real densitized_mass = u0(m,IDN,k,j,i);
      const Real x = CellCenterX(i0, nx1, size.d_view(m).x1min,
                                size.d_view(m).x1max);
      const Real y = CellCenterX(j0, nx2, size.d_view(m).x2min,
                                size.d_view(m).x2max);
      const Real z = CellCenterX(k0, nx3, size.d_view(m).x3min,
                                size.d_view(m).x3max);

      array_sum::GlobalSum cell;
      cell.the_array[10] = coordinate_volume*densitized_mass;
      cell.the_array[11] = proper_volume*rho;
      cell.the_array[12] = proper_volume*pgas;
      cell.the_array[13] = proper_volume*0.5*bsq;
      cell.the_array[14] = coordinate_volume*densitized_mass*lorentz;
      cell.the_array[15] = coordinate_volume*densitized_mass*sigma;
      cell.the_array[16] = coordinate_volume*((x-center_x)*u0(m,IM2,k,j,i)
                                             -(y-center_y)*u0(m,IM1,k,j,i));
      cell.the_array[17] =
          (SQR(x-center_x)+SQR(y-center_y)+SQR(z-center_z) <= inner_radius2)
          ? coordinate_volume*densitized_mass : 0.0;
      sum += cell;
      max_rho = fmax(max_rho, rho);
      max_sigma = fmax(max_sigma, sigma);
    }, Kokkos::Sum<array_sum::GlobalSum>(integral), Kokkos::Max<Real>(rho_max),
    Kokkos::Max<Real>(sigma_max));

  for (int n=10; n<=17; ++n) pdata->hdata[n] = integral.the_array[n];
#if MPI_PARALLEL_ENABLED
  if (global_variable::my_rank == 0) {
    MPI_Reduce(MPI_IN_PLACE, &rho_max, 1, MPI_ATHENA_REAL, MPI_MAX, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(MPI_IN_PLACE, &sigma_max, 1, MPI_ATHENA_REAL, MPI_MAX, 0,
               MPI_COMM_WORLD);
  } else {
    Real unused_max = 0.0;
    MPI_Reduce(&rho_max, &unused_max, 1, MPI_ATHENA_REAL, MPI_MAX, 0,
               MPI_COMM_WORLD);
    MPI_Reduce(&sigma_max, &unused_max, 1, MPI_ATHENA_REAL, MPI_MAX, 0,
               MPI_COMM_WORLD);
    rho_max = 0.0;
    sigma_max = 0.0;
  }
#endif
  pdata->hdata[18] = rho_max;
  pdata->hdata[19] = sigma_max;
}

Real HorizonEnclosingRadius(const TrackingTarget &target) {
  const Real spin2 = SQR(target.spin[0]) + SQR(target.spin[1]) + SQR(target.spin[2]);
  const Real horizon_radius = target.mass
      + sqrt(std::max(SQR(target.mass)-spin2, Real(0.0)));
  const Real rest_enclosing_radius = sqrt(spin2+SQR(horizon_radius));
  const Real velocity2 = SQR(target.velocity[0]) + SQR(target.velocity[1])
      + SQR(target.velocity[2]);
  return rest_enclosing_radius/sqrt(1.0-velocity2);
}

Real TrackingShellRadius(const TrackingTarget &target, const int physical_level) {
  // The user-specified shells already encode the nesting geometry.  The horizon
  // guard is a physical lower bound at every level, not another geometric shell:
  // scaling it by the number of parent levels would refine most of a large domain.
  const Real horizon_radius = bbh.refinement_horizon_factor
      *HorizonEnclosingRadius(target);
  return std::max(refinement_shell_radii[physical_level], horizon_radius);
}

void RefineTracker(MeshBlockPack *pmbp) {
  Mesh *pmesh = pmbp->pmesh;
  auto &refine_flag = pmesh->pmr->refine_flag;
  auto &size = pmbp->pmb->mb_size;
  const int nmb = pmbp->nmb_thispack;
  const int meshblock_start = pmesh->gids_eachrank[global_variable::my_rank];
  // AMR is evaluated at a synchronized root endpoint, before the next CFL timestep is
  // known.  Its growth limiter allows that next root step to be at most twice the last
  // completed step.  Sample the complete look-ahead interval and enlarge each sample
  // sphere by a light-speed half-sample bound, so a subluminal hole cannot leave the
  // refined region before the next synchronized regrid.  The execution tlim is a
  // checkpoint boundary, not a physical trajectory boundary: sample through it when
  // the immutable trajectory has coverage.  If the actual trajectory ends sooner,
  // truncate there; a restart cannot evolve beyond that point with the same table.
  // A segment tlim can shorten the completed step without changing the CFL/growth
  // step that a restart may take next.  Driver preserves that pre-clip reference in
  // dt_restart_growth.  Using the clipped pmesh->dt here would checkpoint a hierarchy
  // that is too small for an immediate full-size continuation.
  const Real fresh_dt = std::numeric_limits<float>::max();
  const bool growth_reference_valid =
      std::isfinite(pmesh->dt_restart_growth) &&
      pmesh->dt_restart_growth > 0.0 && pmesh->dt_restart_growth < fresh_dt;
  const bool completed_dt_valid =
      std::isfinite(pmesh->dt) && pmesh->dt > 0.0 && pmesh->dt < fresh_dt;
  const Real growth_reference = growth_reference_valid
      ? pmesh->dt_restart_growth
      : (completed_dt_valid ? pmesh->dt : 0.0);
  Real lookahead = 2.0*growth_reference;
  // The strict level scheduler may impose a tighter absolute root-step ceiling than
  // its generic factor-two growth limiter.  Honor that same reachable-step bound here;
  // otherwise moving-hole AMR creates and audits a halo for a timestep the Driver can
  // never take, which can cause a one-level startup transient and needless memory use.
  if (bbh.root_dt_max > 0.0) {
    lookahead = std::min(lookahead, bbh.root_dt_max);
  }
  const Real trajectory_horizon = std::max(
      bbh.trajectory_center_end_time-pmesh->time, Real(0.0));
  const Real sample_horizon = std::min(lookahead, trajectory_horizon);
  const Real sample_times[3] = {
    pmesh->time, pmesh->time+static_cast<Real>(0.5)*sample_horizon,
    pmesh->time+sample_horizon
  };
  std::array<std::array<Real, NTRAJ>, 3> states{};
  for (int sample=0; sample<3; ++sample) {
    FindTrajectory(sample_times[sample], states[sample].data());
  }
  const Real padding = static_cast<Real>(0.25)*sample_horizon;
  for (int m = 0; m < nmb; ++m) {
    const RegionSize &block = size.h_view(m);
    const int level = pmbp->pmb->mb_lev.h_view(m);
    const int physical_level = level-pmesh->root_level;
    bool tracker_refine = false;
    if (physical_level < bbh.refinement_floor_level) {
      tracker_refine = DistanceSquaredToBlock(
          bbh.refinement_floor_center[0], bbh.refinement_floor_center[1],
          bbh.refinement_floor_center[2], block) < SQR(bbh.refinement_floor_radius);
    }
    if (level < pmesh->max_level) {
      const int child_physical_level = physical_level+1;
      for (const auto &state : states) {
        TrackingTarget targets[2];
        const int count = BuildTrackingTargets(state.data(), targets);
        for (int target=0; target<count; ++target) {
          const Real radius =
              TrackingShellRadius(targets[target], child_physical_level)+padding;
          tracker_refine = tracker_refine ||
              DistanceSquaredToBlock(targets[target].position[0],
                                     targets[target].position[1],
                                     targets[target].position[2], block) < SQR(radius);
        }
      }
    }
    bool tracker_derefine = level > pmesh->root_level;
    if (!tracker_refine && tracker_derefine) {
      for (const auto &state : states) {
        TrackingTarget targets[2];
        const int count = BuildTrackingTargets(state.data(), targets);
        for (int target=0; target<count; ++target) {
          const Real keep_radius = bbh.refinement_hysteresis
              *TrackingShellRadius(targets[target], physical_level)+padding;
          if (DistanceSquaredToBlock(targets[target].position[0],
                                     targets[target].position[1],
                                     targets[target].position[2], block)
              <= SQR(keep_radius)) {
            tracker_derefine = false;
          }
        }
      }
      if (physical_level <= bbh.refinement_floor_level
          && DistanceSquaredToBlock(
                 bbh.refinement_floor_center[0], bbh.refinement_floor_center[1],
                 bbh.refinement_floor_center[2], block)
              <= SQR(bbh.refinement_hysteresis*bbh.refinement_floor_radius)) {
        tracker_derefine = false;
      }
    }
    int &flag = refine_flag.h_view(m+meshblock_start);
    if (tracker_refine) flag = 1;
    if (tracker_derefine && flag == 0) flag = -1;
  }
  refine_flag.template modify<HostMemSpace>();
  refine_flag.template sync<DevExeSpace>();
}

} // namespace
