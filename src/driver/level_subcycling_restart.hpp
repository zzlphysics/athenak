#ifndef DRIVER_LEVEL_SUBCYCLING_RESTART_HPP_
#define DRIVER_LEVEL_SUBCYCLING_RESTART_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file level_subcycling_restart.hpp
//! \brief Binary restart contract for exact strict-subcycling continuation.

#include <cstdint>

namespace level_subcycling {

// "ATKSCC01": checkpoints carrying this extension serialize the synchronized U/B
// active and ghost state consumed by the next strict-subcycling root step.
constexpr std::uint64_t kRestartCacheContractMagic =
    UINT64_C(0x41544b5343433031);
constexpr int kRestartCacheContractVersion = 1;

} // namespace level_subcycling

#endif // DRIVER_LEVEL_SUBCYCLING_RESTART_HPP_
