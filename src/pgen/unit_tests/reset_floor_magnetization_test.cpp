//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file reset_floor_magnetization_test.cpp
//! \brief Device regression for PrimitiveSolver excess-magnetization normalization.

#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "eos/primitive-solver/eos.hpp"
#include "eos/primitive-solver/idealgas.hpp"
#include "eos/primitive-solver/primitive_solver.hpp"
#include "eos/primitive-solver/reset_floor.hpp"
#include "pgen/pgen.hpp"

namespace {

// A test-only derived policy gives the metric regression a unique solver
// template instantiation, so it cannot accidentally resolve to an identical
// weak instantiation emitted by another translation unit.
class MetricResetFloor : public Primitive::ResetFloor {
 public:
  MetricResetFloor() = default;
};

}  // namespace

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  (void)pin;
  (void)restart;

  Primitive::EOS<Primitive::IdealGas, Primitive::ResetFloor> eos;
  eos.SetMaximumMagnetization(3.0);
  Primitive::PrimitiveSolver<Primitive::IdealGas, MetricResetFloor> metric_ps;
  metric_ps.GetEOSMutable().SetMaximumMagnetization(3.0);

  int failures = 0;
  Kokkos::parallel_reduce(
      "reset_floor_magnetization_test",
      Kokkos::RangePolicy<DevExeSpace>(0, 3),
      KOKKOS_LAMBDA(const int test, int &sum) {
        if (test == 2) {
          // The two states below differ only by the diagonal coordinate
          // transformation dx_phys^i = scale[i] dx^i.  A correct C2P solve
          // must therefore return the same physical state after transforming
          // the vector components back.  In particular, this exercises the
          // post-magnetization recomputation of r^i with the inverse metric.
          // Using g_ij there instead of g^ij breaks this covariance check by a
          // large margin while still passing identity-metric tests.
          const Real scale[3] = {2.0, 3.0, 4.0};
          Real g3d_cart[NSPMETRIC] = {1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
          Real g3u_cart[NSPMETRIC] = {1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
          Real g3d_scaled[NSPMETRIC] = {
              scale[0]*scale[0], 0.0, 0.0,
              scale[1]*scale[1], 0.0, scale[2]*scale[2]};
          Real g3u_scaled[NSPMETRIC] = {
              1.0/(scale[0]*scale[0]), 0.0, 0.0,
              1.0/(scale[1]*scale[1]), 0.0, 1.0/(scale[2]*scale[2])};

          Real prim_cart[NPRIM] = {0.0};
          prim_cart[PRH] = 1.0;
          prim_cart[PVX] = 0.18;
          prim_cart[PVY] = -0.11;
          prim_cart[PVZ] = 0.07;
          prim_cart[PPR] = 0.2;
          prim_cart[PTM] = 0.2;
          Real prim_scaled[NPRIM] = {0.0};
          prim_scaled[PRH] = prim_cart[PRH];
          prim_scaled[PVX] = prim_cart[PVX]/scale[0];
          prim_scaled[PVY] = prim_cart[PVY]/scale[1];
          prim_scaled[PVZ] = prim_cart[PVZ]/scale[2];
          prim_scaled[PPR] = prim_cart[PPR];
          prim_scaled[PTM] = prim_cart[PTM];

          Real b_cart[NMAG] = {2.0, -1.0, 0.5};
          Real b_scaled[NMAG] = {
              b_cart[IBX]/scale[0], b_cart[IBY]/scale[1], b_cart[IBZ]/scale[2]};
          Real cons_cart[NCONS] = {0.0};
          Real cons_scaled[NCONS] = {0.0};
          metric_ps.PrimToCon(prim_cart, cons_cart, b_cart, g3d_cart);
          metric_ps.PrimToCon(prim_scaled, cons_scaled, b_scaled, g3d_scaled);

          Real recovered_cart[NPRIM] = {0.0};
          Real recovered_scaled[NPRIM] = {0.0};
          const Primitive::SolverResult result_cart = metric_ps.ConToPrim(
              recovered_cart, cons_cart, b_cart, g3d_cart, g3u_cart);
          const Primitive::SolverResult result_scaled = metric_ps.ConToPrim(
              recovered_scaled, cons_scaled, b_scaled, g3d_scaled, g3u_scaled);

          const Real tolerance = 1.0e-11;
          bool mismatch = result_cart.error != Primitive::Error::SUCCESS ||
                          result_scaled.error != Primitive::Error::SUCCESS ||
                          !result_cart.cons_adjusted || !result_scaled.cons_adjusted ||
                          (result_cart.events & Primitive::MAGNETIZATION_ADJUSTED) == 0U ||
                          (result_scaled.events & Primitive::MAGNETIZATION_ADJUSTED) == 0U;
          for (int n = 0; n < NPRIM; ++n) {
            mismatch = mismatch || !isfinite(recovered_cart[n]) ||
                       !isfinite(recovered_scaled[n]);
          }
          for (int n = 0; n < NCONS; ++n) {
            mismatch = mismatch || !isfinite(cons_cart[n]) || !isfinite(cons_scaled[n]);
          }
          mismatch = mismatch || fabs(recovered_cart[PRH] - recovered_scaled[PRH])
                                    > tolerance;
          mismatch = mismatch || fabs(recovered_cart[PPR] - recovered_scaled[PPR])
                                    > tolerance;
          mismatch = mismatch || fabs(recovered_cart[PTM] - recovered_scaled[PTM])
                                    > tolerance;
          mismatch = mismatch || fabs(cons_cart[CDN] - cons_scaled[CDN]) > tolerance;
          mismatch = mismatch || fabs(cons_cart[CTA] - cons_scaled[CTA]) > tolerance;
          for (int n = 0; n < 3; ++n) {
            mismatch = mismatch ||
                fabs(recovered_cart[PVX + n] - scale[n]*recovered_scaled[PVX + n])
                    > tolerance;
            mismatch = mismatch ||
                fabs(cons_cart[CSX + n] - cons_scaled[CSX + n]/scale[n]) > tolerance;
          }
          if (mismatch) {
            ++sum;
          }
          return;
        }

        Real bsq = 3.0;
        Real b_u[3] = {1.0, -1.0, 1.0};
        Primitive::Error expected = Primitive::Error::SUCCESS;
        if (test == 0) {
          bsq = 12.0;
          b_u[0] = 2.0;
          b_u[1] = -2.0;
          b_u[2] = 2.0;
          expected = Primitive::Error::CONS_ADJUSTED;
        }

        const Primitive::Error actual = eos.DoMagnetizationResponse(bsq, b_u);
        const Real expected_b0 = 1.0;
        const Real expected_b1 = -1.0;
        const Real expected_b2 = 1.0;
        const Real norm = b_u[0]*b_u[0] + b_u[1]*b_u[1] + b_u[2]*b_u[2];
        if (actual != expected || bsq != 3.0 || b_u[0] != expected_b0 ||
            b_u[1] != expected_b1 || b_u[2] != expected_b2 || norm != bsq) {
          ++sum;
        }
      },
      failures);

  if (failures != 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << ": " << failures
              << " reset-floor magnetization checks failed" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::cout << "ResetFloor magnetization device test passed" << std::endl;
}
