#ifndef EOS_PRIMITIVE_SOLVER_PS_ERROR_HPP_
#define EOS_PRIMITIVE_SOLVER_PS_ERROR_HPP_
//========================================================================================
// PrimitiveSolver equation-of-state framework
// Copyright(C) 2023 Jacob M. Fields <jmf6719@psu.edu>
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file ps_error.hpp
//  \brief defines an enumerator struct for error types.

#include <cstdint>

namespace Primitive {
enum struct Error {
  SUCCESS,
  RHO_TOO_BIG,
  RHO_TOO_SMALL,
  NANS_IN_CONS,
  MAG_TOO_BIG,
  BRACKETING_FAILED,
  NO_SOLUTION,
  CONS_FLOOR,
  PRIM_FLOOR,
  CONS_ADJUSTED,
};

enum SolverEvent : std::uint32_t {
  CONS_DENSITY_FLOOR = 1U << 0,
  CONS_ENERGY_FLOOR = 1U << 1,
  PRIM_DENSITY_FLOOR = 1U << 2,
  PRIM_TEMPERATURE_FLOOR = 1U << 3,
  MAGNETIZATION_ADJUSTED = 1U << 4,
};

struct SolverResult {
  Error error;
  int  iterations;
  bool cons_floor;
  bool prim_floor;
  bool cons_adjusted;
  std::uint32_t events = 0;
};

} // namespace Primitive

#endif  // EOS_PRIMITIVE_SOLVER_PS_ERROR_HPP_
