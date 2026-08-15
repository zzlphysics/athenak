#ifndef MHD_C2P_TELEMETRY_HPP_
#define MHD_C2P_TELEMETRY_HPP_
//========================================================================================
// Athena++K astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file c2p_telemetry.hpp
//! \brief Fixed, diagnostic-only binning for dynamical-GRMHD C2P interventions.

#include <cmath>
#include <cstddef>
#include <cstdint>

#include "athena.hpp"
#include "athena_tensor.hpp"
#include "coordinates/cell_locations.hpp"
#include "eos/primitive-solver/ps_error.hpp"
#include "mesh/mesh.hpp"
#include "mhd/fofc_telemetry.hpp"

namespace mhd {
namespace c2p_telemetry {

// Keep the layout fixed so event-log comments remain machine-readable across MPI
// decompositions and AMR changes.  The spatial bins deliberately match FOFC v1.
constexpr int kLevelBins = fofc_telemetry::kLevelBins;
constexpr int kLevelOverflow = fofc_telemetry::kLevelOverflow;
constexpr int kRadiusBins = fofc_telemetry::kRadiusBins;
constexpr int kAbsZBins = fofc_telemetry::kAbsZBins;
constexpr int kLapseBins = fofc_telemetry::kLapseBins;
constexpr int kStageBins = fofc_telemetry::kStageBins;
constexpr int kStageOther = fofc_telemetry::kStageOther;

enum Intervention : std::uint8_t {
  cons_adjust = 0,
  mag_adjust,
  nintervention
};
constexpr int kInterventionBins = static_cast<int>(Intervention::nintervention);

// Density is normalized by the configured rest-mass atmosphere floor D_floor.
// Finite-bin edges: 1, 2, 4, 16, 64, 256.  The final bin is invalid/non-finite.
constexpr int kFiniteDensityRatioBins = 7;
constexpr int kDensityRatioInvalid = kFiniteDensityRatioBins;
constexpr int kDensityRatioBins = kFiniteDensityRatioBins + 1;

// Magnetization is B^2/D immediately before the magnetization response (therefore
// after any conserved-density floor) normalized by max_bsq.
// Finite-bin edges: 0.01, 0.1, 0.5, 1, 2, 10.  The final bin is invalid/non-finite.
constexpr int kFiniteMagRatioBins = 7;
constexpr int kMagRatioInvalid = kFiniteMagRatioBins;
constexpr int kMagRatioBins = kFiniteMagRatioBins + 1;

// The main joint distribution preserves the correlations required to distinguish a
// floor-dominated atmosphere/funnel intervention from one in resolved disk material.
// Stage and geometry-validity marginals are appended to the same allocation.
constexpr std::size_t kJointHistogramSize =
    static_cast<std::size_t>(kInterventionBins)*kLevelBins*kRadiusBins*kAbsZBins*
    kLapseBins*kDensityRatioBins*kMagRatioBins;
constexpr std::size_t kStageHistogramSize =
    static_cast<std::size_t>(kInterventionBins)*kStageBins;
constexpr std::size_t kGeometryHistogramSize =
    static_cast<std::size_t>(kInterventionBins)*2;
constexpr std::size_t kStageHistogramOffset = kJointHistogramSize;
constexpr std::size_t kGeometryHistogramOffset =
    kStageHistogramOffset + kStageHistogramSize;
constexpr std::size_t kPendingSize =
    kGeometryHistogramOffset + kGeometryHistogramSize;

KOKKOS_INLINE_FUNCTION
int DensityRatioBin(const Real ratio) {
  if (!Kokkos::isfinite(ratio)) return kDensityRatioInvalid;
  if (ratio < 1.0) return 0;
  if (ratio < 2.0) return 1;
  if (ratio < 4.0) return 2;
  if (ratio < 16.0) return 3;
  if (ratio < 64.0) return 4;
  if (ratio < 256.0) return 5;
  return 6;
}

KOKKOS_INLINE_FUNCTION
int MagnetizationRatioBin(const Real ratio) {
  if (!Kokkos::isfinite(ratio)) return kMagRatioInvalid;
  if (ratio < 0.01) return 0;
  if (ratio < 0.1) return 1;
  if (ratio < 0.5) return 2;
  if (ratio < 1.0) return 3;
  if (ratio < 2.0) return 4;
  if (ratio < 10.0) return 5;
  return 6;
}

KOKKOS_INLINE_FUNCTION
constexpr std::size_t JointHistogramIndex(
    const int intervention, const int level, const int radius, const int abs_z,
    const int lapse, const int density_ratio, const int mag_ratio) {
  return ((((((static_cast<std::size_t>(intervention)*kLevelBins + level)*
              kRadiusBins + radius)*kAbsZBins + abs_z)*kLapseBins + lapse)*
              kDensityRatioBins + density_ratio)*kMagRatioBins + mag_ratio);
}

KOKKOS_INLINE_FUNCTION
constexpr std::size_t StageHistogramIndex(const int intervention, const int stage) {
  return kStageHistogramOffset +
         static_cast<std::size_t>(intervention)*kStageBins + stage;
}

KOKKOS_INLINE_FUNCTION
constexpr std::size_t GeometryHistogramIndex(
    const int intervention, const bool valid) {
  return kGeometryHistogramOffset +
         static_cast<std::size_t>(intervention)*2 + static_cast<int>(valid);
}

inline const char *InterventionName(const int intervention) {
  static constexpr const char *names[kInterventionBins] = {
    "cons_adjust", "mag_adjust"
  };
  return (intervention >= 0 && intervention < kInterventionBins)
      ? names[intervention] : "invalid";
}

// Keep the disabled instantiation empty.  This prevents the default C2P device kernel
// from capturing telemetry Views or adding a per-cell runtime branch.
template<bool ENABLED>
struct Recorder {
  KOKKOS_INLINE_FUNCTION
  void Record(const bool, const bool, const int, const int, const int, const int,
              const Real, const Real, const Real) const {}
};

template<>
struct Recorder<true> {
  DvceArray1D<std::uint64_t> pending;
  DvceArray1D<int> level;
  DvceArray1D<RegionSize> size;
  AthenaTensor<Real, TensorSymm::NONE, 3, 0> lapse;
  int is, js, ks;
  int nx1, nx2, nx3;
  int stage;
  Real center1, center2, center3;
  Real density_floor;
  Real max_magnetization;

  KOKKOS_INLINE_FUNCTION
  void Record(const bool cons_adjusted, const bool mag_adjusted,
              const int m, const int k, const int j, const int i,
              const Real input_density, const Real magnetization_density,
              const Real input_bsq) const {
    if (!cons_adjusted && !mag_adjusted) return;

    const Real x1 = CellCenterX(i-is, nx1, size(m).x1min, size(m).x1max);
    const Real x2 = CellCenterX(j-js, nx2, size(m).x2min, size(m).x2max);
    const Real x3 = CellCenterX(k-ks, nx3, size(m).x3min, size(m).x3max);
    const Real alpha = lapse(m,k,j,i);
    const bool geometry_valid =
        Kokkos::isfinite(x1) && Kokkos::isfinite(x2) && Kokkos::isfinite(x3) &&
        Kokkos::isfinite(alpha) && alpha >= 0.0;
    const int level_bin = fofc_telemetry::LevelBin(level(m));
    const int radius_bin = geometry_valid
        ? fofc_telemetry::RadiusBin(x1, x2, center1, center2) : kRadiusBins-1;
    const int abs_z_bin = geometry_valid
        ? fofc_telemetry::AbsZBin(x3, center3) : kAbsZBins-1;
    const int lapse_bin = geometry_valid
        ? fofc_telemetry::LapseBin(alpha) : kLapseBins-1;
    const int stage_bin = fofc_telemetry::StageBin(stage);
    const bool density_scale_valid = density_floor > 0.0 &&
                                     Kokkos::isfinite(density_floor);
    const bool magnetization_scale_valid = magnetization_density > 0.0 &&
        max_magnetization > 0.0 && Kokkos::isfinite(magnetization_density) &&
        Kokkos::isfinite(max_magnetization);
    const int density_bin = density_scale_valid
        ? DensityRatioBin(input_density/density_floor) : kDensityRatioInvalid;
    const int mag_bin = magnetization_scale_valid
        ? MagnetizationRatioBin(
            input_bsq/magnetization_density/max_magnetization)
        : kMagRatioInvalid;

    for (int intervention=0; intervention<kInterventionBins; ++intervention) {
      const bool selected = (intervention == Intervention::cons_adjust)
          ? cons_adjusted : mag_adjusted;
      if (!selected) continue;
      const std::size_t joint_index = JointHistogramIndex(
          intervention, level_bin, radius_bin, abs_z_bin, lapse_bin,
          density_bin, mag_bin);
      Kokkos::atomic_fetch_add(&pending(joint_index), std::uint64_t{1});
      Kokkos::atomic_fetch_add(
          &pending(StageHistogramIndex(intervention, stage_bin)), std::uint64_t{1});
      Kokkos::atomic_fetch_add(
          &pending(GeometryHistogramIndex(intervention, geometry_valid)),
          std::uint64_t{1});
    }
  }
};

static_assert(kJointHistogramSize == 1241856,
              "C2P telemetry schema/memory budget changed");
static_assert(kPendingSize == 1241868,
              "C2P telemetry appended marginals changed");
static_assert(JointHistogramIndex(kInterventionBins-1, kLevelBins-1,
                                  kRadiusBins-1, kAbsZBins-1, kLapseBins-1,
                                  kDensityRatioBins-1, kMagRatioBins-1) + 1 ==
              kJointHistogramSize);
static_assert(StageHistogramIndex(kInterventionBins-1, kStageBins-1) + 1 ==
              kGeometryHistogramOffset);
static_assert(GeometryHistogramIndex(kInterventionBins-1, true) + 1 == kPendingSize);

} // namespace c2p_telemetry
} // namespace mhd

#endif // MHD_C2P_TELEMETRY_HPP_
