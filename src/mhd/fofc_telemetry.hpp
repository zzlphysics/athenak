#ifndef MHD_FOFC_TELEMETRY_HPP_
#define MHD_FOFC_TELEMETRY_HPP_
//========================================================================================
// Athena++K astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file fofc_telemetry.hpp
//! \brief Fixed, diagnostic-only binning for dynamical-GRMHD FOFC corrections.

#include <cmath>
#include <cstddef>
#include <cstdint>

#include "athena.hpp"
#include "eos/primitive-solver/ps_error.hpp"

namespace mhd {
namespace fofc_telemetry {

// Keep this layout fixed so comment records in the event log remain machine-readable
// across MPI decompositions and AMR changes.  Logical levels 0--31 have distinct bins;
// anything outside that range is retained in the overflow bin.
constexpr int kExplicitLevelBins = 32;
constexpr int kLevelBins = kExplicitLevelBins + 1;
constexpr int kLevelOverflow = kExplicitLevelBins;

// Bin zero is reserved for restart-carried or otherwise unattributed corrections.
constexpr int kStageBins = 4;
constexpr int kStageOther = 0;

enum Reason : std::uint8_t {
  unknown = 0,
  dmp_preflag,
  scalar,
  cons_density_floor,
  cons_energy_floor,
  prim_density_floor,
  prim_temperature_floor,
  rho_too_big,
  rho_too_small,
  nans_in_cons,
  mag_too_big,
  bracketing_failed,
  no_solution,
  invalid_geometry,
  other_c2p,
  nreason
};

constexpr int kReasonBins = static_cast<int>(Reason::nreason);

// Spatial bins use cylindrical distance from a configurable fixed center, |z-center_z|,
// and lapse.  Lapse supplies a moving-hole proximity proxy without coupling this generic
// MHD diagnostic to a particular problem generator or trajectory implementation.
//
//   R_cyl edges: 2, 4, 8, 16, 32, 64
//   |z| edges:   0.5, 1, 2, 4, 8, 16
//   lapse edges: 0.2, 0.4, 0.6, 0.8, 1.0
constexpr int kRadiusBins = 7;
constexpr int kAbsZBins = 7;
constexpr int kLapseBins = 6;

constexpr std::size_t kHistogramSize =
    static_cast<std::size_t>(kLevelBins)*kStageBins*kReasonBins*kRadiusBins*
    kAbsZBins*kLapseBins;

KOKKOS_INLINE_FUNCTION
constexpr std::uint8_t ReasonFromSolver(const Primitive::Error error,
                                        const std::uint32_t events) {
  // Floors-only C2P currently emits at most one of these floor bits.  Keep an explicit,
  // deterministic priority nevertheless: if a future solver emits several, report the
  // first intervention in conserved-to-primitive order.  A magnetization adjustment may
  // precede a later failure; the failure that actually triggers FOFC remains the cause.
  if ((events & Primitive::CONS_DENSITY_FLOOR) != 0U) {
    return Reason::cons_density_floor;
  }
  if ((events & Primitive::CONS_ENERGY_FLOOR) != 0U) {
    return Reason::cons_energy_floor;
  }
  if ((events & Primitive::PRIM_DENSITY_FLOOR) != 0U) {
    return Reason::prim_density_floor;
  }
  if ((events & Primitive::PRIM_TEMPERATURE_FLOOR) != 0U) {
    return Reason::prim_temperature_floor;
  }
  switch (error) {
    case Primitive::Error::RHO_TOO_BIG:
      return Reason::rho_too_big;
    case Primitive::Error::RHO_TOO_SMALL:
      return Reason::rho_too_small;
    case Primitive::Error::NANS_IN_CONS:
      return Reason::nans_in_cons;
    case Primitive::Error::MAG_TOO_BIG:
      return Reason::mag_too_big;
    case Primitive::Error::BRACKETING_FAILED:
      return Reason::bracketing_failed;
    case Primitive::Error::NO_SOLUTION:
      return Reason::no_solution;
    default:
      return Reason::other_c2p;
  }
}

KOKKOS_INLINE_FUNCTION
constexpr int LevelBin(const int level) {
  return (level >= 0 && level < kExplicitLevelBins) ? level : kLevelOverflow;
}

KOKKOS_INLINE_FUNCTION
constexpr int StageBin(const int stage) {
  return (stage >= 1 && stage < kStageBins) ? stage : kStageOther;
}

KOKKOS_INLINE_FUNCTION
int RadiusBin(const Real x, const Real y, const Real center_x, const Real center_y) {
  const Real dx = x-center_x;
  const Real dy = y-center_y;
  const Real radius2 = dx*dx + dy*dy;
  if (radius2 < 4.0) return 0;
  if (radius2 < 16.0) return 1;
  if (radius2 < 64.0) return 2;
  if (radius2 < 256.0) return 3;
  if (radius2 < 1024.0) return 4;
  if (radius2 < 4096.0) return 5;
  return 6;
}

KOKKOS_INLINE_FUNCTION
int AbsZBin(const Real z, const Real center_z) {
  const Real abs_z = fabs(z-center_z);
  if (abs_z < 0.5) return 0;
  if (abs_z < 1.0) return 1;
  if (abs_z < 2.0) return 2;
  if (abs_z < 4.0) return 3;
  if (abs_z < 8.0) return 4;
  if (abs_z < 16.0) return 5;
  return 6;
}

KOKKOS_INLINE_FUNCTION
int LapseBin(const Real lapse) {
  if (lapse >= 0.0 && lapse < 0.2) return 0;
  if (lapse >= 0.2 && lapse < 0.4) return 1;
  if (lapse >= 0.4 && lapse < 0.6) return 2;
  if (lapse >= 0.6 && lapse < 0.8) return 3;
  if (lapse >= 0.8 && lapse < 1.0) return 4;
  return 5;
}

KOKKOS_INLINE_FUNCTION
constexpr std::size_t HistogramIndex(const int level_bin, const int stage_bin,
                                     const int reason_bin, const int radius_bin,
                                     const int abs_z_bin, const int lapse_bin) {
  return (((((static_cast<std::size_t>(level_bin)*kStageBins + stage_bin)*kReasonBins
              + reason_bin)*kRadiusBins + radius_bin)*kAbsZBins + abs_z_bin)*
              kLapseBins + lapse_bin);
}

inline const char *ReasonName(const int reason) {
  static constexpr const char *names[kReasonBins] = {
    "unknown", "dmp_preflag", "scalar", "cons_density_floor",
    "cons_energy_floor", "prim_density_floor", "prim_temperature_floor",
    "rho_too_big", "rho_too_small", "nans_in_cons", "mag_too_big",
    "bracketing_failed", "no_solution", "invalid_geometry", "other_c2p"
  };
  return (reason >= 0 && reason < kReasonBins) ? names[reason] : "invalid";
}

static_assert(ReasonFromSolver(Primitive::Error::CONS_FLOOR,
                               Primitive::CONS_DENSITY_FLOOR) ==
              Reason::cons_density_floor);
static_assert(ReasonFromSolver(Primitive::Error::PRIM_FLOOR,
                               Primitive::PRIM_TEMPERATURE_FLOOR) ==
              Reason::prim_temperature_floor);
static_assert(ReasonFromSolver(Primitive::Error::NO_SOLUTION, 0U) ==
              Reason::no_solution);
static_assert(kHistogramSize == 582120,
              "FOFC telemetry schema/memory budget changed");
static_assert(HistogramIndex(kLevelBins-1, kStageBins-1, kReasonBins-1,
                             kRadiusBins-1, kAbsZBins-1, kLapseBins-1) + 1 ==
              kHistogramSize);

} // namespace fofc_telemetry
} // namespace mhd

#endif // MHD_FOFC_TELEMETRY_HPP_
