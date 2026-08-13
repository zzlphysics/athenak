#ifndef DYN_GRMHD_C2P_PROJECTION_STATS_HPP_
#define DYN_GRMHD_C2P_PROJECTION_STATS_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file c2p_projection_stats.hpp
//! \brief Lightweight diagnostics for conservative-to-primitive cache projections.

#include <cstdint>

#include "athena.hpp"

namespace dyngr {

struct C2PProjectionStats {
  std::uint64_t solver_failures = 0;
  std::uint64_t preserved_cons_mismatches = 0;
  std::uint64_t preserved_cons_nonfinite = 0;
  Real max_preserved_cons_relative_change = 0.0;
  Real max_preserved_cons_absolute_change = 0.0;
};

} // namespace dyngr

#endif // DYN_GRMHD_C2P_PROJECTION_STATS_HPP_
