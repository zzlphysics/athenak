//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file mhd_newdt.cpp
//! \brief function to compute MHD timestep across all MeshBlock(s) in a MeshBlockPack

#include <math.h>

#include <limits>
#include <iostream>
#include <algorithm> // min

#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "driver/driver.hpp"
#include "coordinates/adm.hpp"
#include "eos/eos.hpp"
#include "mhd.hpp"
#include "diffusion/conduction.hpp"
#include "srcterms/srcterms.hpp"

namespace mhd {

KOKKOS_INLINE_FUNCTION
Real ADMCoordinateLightSpeed(const Real metric[NSPMETRIC], const Real shift[3],
                             const Real lapse, const int direction) {
  for (int n=0; n<NSPMETRIC; ++n) {
    if (!isfinite(metric[n])) {
      Kokkos::abort("Non-finite ADM metric while computing the GRMHD timestep");
    }
  }
  for (int n=0; n<3; ++n) {
    if (!isfinite(shift[n])) {
      Kokkos::abort("Non-finite ADM shift while computing the GRMHD timestep");
    }
  }
  const Real determinant = adm::SpatialDet(
      metric[S11], metric[S12], metric[S13], metric[S22], metric[S23], metric[S33]);
  const Real leading_minor2 = metric[S11]*metric[S22] - SQR(metric[S12]);
  if (!isfinite(lapse) || !(metric[S11] > 0.0) || !(leading_minor2 > 0.0)
      || !(determinant > 0.0) || !(lapse > 0.0)) {
    Kokkos::abort("Invalid ADM metric while computing the GRMHD timestep");
  }
  Real guxx, guxy, guxz, guyy, guyz, guzz;
  adm::SpatialInv(1.0/determinant,
      metric[S11], metric[S12], metric[S13], metric[S22], metric[S23], metric[S33],
      &guxx, &guxy, &guxz, &guyy, &guyz, &guzz);
  const Real inverse_diagonal[3] = {guxx, guyy, guzz};
  if (!(inverse_diagonal[direction] > 0.0)) {
    Kokkos::abort("Non-positive inverse ADM metric while computing the GRMHD timestep");
  }
  const Real speed = fabs(shift[direction])
      + lapse*sqrt(inverse_diagonal[direction]);
  if (!isfinite(speed) || !(speed > 0.0)) {
    Kokkos::abort(
        "Invalid ADM coordinate light speed while computing the GRMHD timestep");
  }
  return speed;
}

//----------------------------------------------------------------------------------------
// \!fn void MHD::NewTimeStep()
// \brief calculate the minimum timestep within a MeshBlockPack for MHD problems

TaskStatus MHD::NewTimeStep(Driver *pdriver, int stage) {
  if (stage != (pdriver->nexp_stages)) {
    return TaskStatus::complete; // only execute last stage
  }

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  int is = indcs.is, nx1 = indcs.nx1;
  int js = indcs.js, nx2 = indcs.nx2;
  int ks = indcs.ks, nx3 = indcs.nx3;

  Real dt1 = std::numeric_limits<float>::max();
  Real dt2 = std::numeric_limits<float>::max();
  Real dt3 = std::numeric_limits<float>::max();

  // capture class variables for kernel
  auto &w0_ = w0;
  auto &eos = pmy_pack->pmhd->peos->eos_data;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto &is_special_relativistic_ = pmy_pack->pcoord->is_special_relativistic;
  auto &is_general_relativistic_ = pmy_pack->pcoord->is_general_relativistic;
  auto &is_dynamical_relativistic_ = pmy_pack->pcoord->is_dynamical_relativistic;
  adm::ADM::ADM_vars adm_vars;
  if (is_dynamical_relativistic_) adm_vars = pmy_pack->padm->adm;
  auto active_lids = pmy_pack->active_lids.d_view;
  const int active_offset = pmy_pack->active_offset;
  const int nmb_active = pmy_pack->nmb_active;
  if (nmb_active <= 0) {
    dtnew = std::numeric_limits<float>::max();
    return TaskStatus::complete;
  }
  const int nmkji = nmb_active*nx3*nx2*nx1;
  const int nkji = nx3*nx2*nx1;
  const int nji  = nx2*nx1;

  if (pdriver->time_evolution == TimeEvolution::kinematic) {
    // find smallest (dx/v) in each direction for advection problems
    Kokkos::parallel_reduce("MHDNudt1",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
    KOKKOS_LAMBDA(const int &idx, Real &min_dt1, Real &min_dt2, Real &min_dt3) {
      // compute m,k,j,i indices of thread and call function
      const int a = idx/nkji;
      const int m = active_lids(active_offset + a);
      int k = (idx - a*nkji)/nji;
      int j = (idx - a*nkji - k*nji)/nx1;
      int i = (idx - a*nkji - k*nji - j*nx1) + is;
      k += ks;
      j += js;

      min_dt1 = fmin((mbsize.d_view(m).dx1/fabs(w0_(m,IVX,k,j,i))), min_dt1);
      min_dt2 = fmin((mbsize.d_view(m).dx2/fabs(w0_(m,IVY,k,j,i))), min_dt2);
      min_dt3 = fmin((mbsize.d_view(m).dx3/fabs(w0_(m,IVZ,k,j,i))), min_dt3);
    }, Kokkos::Min<Real>(dt1), Kokkos::Min<Real>(dt2),Kokkos::Min<Real>(dt3));
  } else {
    // find smallest dx/(v +/- Cf) in each direction for mhd problems
    auto &bcc0_ = bcc0;

    Kokkos::parallel_reduce("MHDNudt2",Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
    KOKKOS_LAMBDA(const int &idx, Real &min_dt1, Real &min_dt2, Real &min_dt3) {
      // compute m,k,j,i indices of thread and call function
      const int a = idx/nkji;
      const int m = active_lids(active_offset + a);
      int k = (idx - a*nkji)/nji;
      int j = (idx - a*nkji - k*nji)/nx1;
      int i = (idx - a*nkji - k*nji - j*nx1) + is;
      k += ks;
      j += js;
      Real max_dv1 = 0.0, max_dv2 = 0.0, max_dv3 = 0.0;

      // Bound each characteristic speed by the coordinate light speed at the same two
      // face metrics used by the Riemann solver.  The old unit-speed assumption is unsafe
      // for a boosted or otherwise time-dependent metric.
      if (is_dynamical_relativistic_) {
        Real face_metric[NSPMETRIC], face_shift[3], face_lapse;
        adm::Face1Metric(m, k, j, i, adm_vars.g_dd, adm_vars.beta_u, adm_vars.alpha,
                         face_metric, face_shift, face_lapse);
        max_dv1 = ADMCoordinateLightSpeed(face_metric, face_shift, face_lapse, 0);
        adm::Face1Metric(m, k, j, i+1, adm_vars.g_dd, adm_vars.beta_u, adm_vars.alpha,
                         face_metric, face_shift, face_lapse);
        max_dv1 = fmax(max_dv1,
                       ADMCoordinateLightSpeed(face_metric, face_shift, face_lapse, 0));
        if (nx2 > 1) {
          adm::Face2Metric(m, k, j, i, adm_vars.g_dd, adm_vars.beta_u, adm_vars.alpha,
                           face_metric, face_shift, face_lapse);
          max_dv2 = ADMCoordinateLightSpeed(face_metric, face_shift, face_lapse, 1);
          adm::Face2Metric(m, k, j+1, i, adm_vars.g_dd, adm_vars.beta_u, adm_vars.alpha,
                           face_metric, face_shift, face_lapse);
          max_dv2 = fmax(max_dv2,
                         ADMCoordinateLightSpeed(face_metric, face_shift, face_lapse, 1));
        }
        if (nx3 > 1) {
          adm::Face3Metric(m, k, j, i, adm_vars.g_dd, adm_vars.beta_u, adm_vars.alpha,
                           face_metric, face_shift, face_lapse);
          max_dv3 = ADMCoordinateLightSpeed(face_metric, face_shift, face_lapse, 2);
          adm::Face3Metric(m, k+1, j, i, adm_vars.g_dd, adm_vars.beta_u, adm_vars.alpha,
                           face_metric, face_shift, face_lapse);
          max_dv3 = fmax(max_dv3,
                         ADMCoordinateLightSpeed(face_metric, face_shift, face_lapse, 2));
        }
      // timestep in stationary GR MHD (analytic Kerr-Schild path)
      } else if (is_general_relativistic_) {
        max_dv1 = 1.0;
        max_dv2 = 1.0;
        max_dv3 = 1.0;
      // timestep in SR MHD
      } else if (is_special_relativistic_) {
        Real &wd = w0_(m,IDN,k,j,i);
        Real &ux = w0_(m,IVX,k,j,i);
        Real &uy = w0_(m,IVY,k,j,i);
        Real &uz = w0_(m,IVZ,k,j,i);
        Real &bcc1 = bcc0_(m,IBX,k,j,i);
        Real &bcc2 = bcc0_(m,IBY,k,j,i);
        Real &bcc3 = bcc0_(m,IBZ,k,j,i);

        Real v2 = SQR(ux) + SQR(uy) + SQR(uz);
        Real lor = sqrt(1.0 + v2);
        // FIXME ERM: Ideal fluid for now
        Real p = eos.IdealGasPressure(w0_(m,IEN,k,j,i));
        // Calculate 4-magnetic field in left state
        Real b_0 = bcc1*ux + bcc2*uy + bcc3*uz;
        Real b_1 = (bcc1 + b_0 * ux) / lor;
        Real b_2 = (bcc2 + b_0 * uy) / lor;
        Real b_3 = (bcc3 + b_0 * uz) / lor;
        Real b_sq = -SQR(b_0) + SQR(b_1) + SQR(b_2) + SQR(b_3);

        Real lm, lp;
        eos.IdealSRMHDFastSpeeds(wd, p, ux, lor, b_sq, lp, lm);
        max_dv1 = fmax(fabs(lm), lp);

        eos.IdealSRMHDFastSpeeds(wd, p, uy, lor, b_sq, lp, lm);
        max_dv2 = fmax(fabs(lm), lp);

        eos.IdealSRMHDFastSpeeds(wd, p, uz, lor, b_sq, lp, lm);
        max_dv3 = fmax(fabs(lm), lp);
      // timestep in Newtonian MHD
      } else {
        Real &w_d = w0_(m,IDN,k,j,i);
        Real &w_bx = bcc0_(m,IBX,k,j,i);
        Real &w_by = bcc0_(m,IBY,k,j,i);
        Real &w_bz = bcc0_(m,IBZ,k,j,i);
        Real cf;
        if (eos.is_ideal) {
          Real p = eos.IdealGasPressure(w0_(m,IEN,k,j,i));
          cf = eos.IdealMHDFastSpeed(w_d, p, w_bx, w_by, w_bz);
          max_dv1 = fabs(w0_(m,IVX,k,j,i)) + cf;
          cf = eos.IdealMHDFastSpeed(w_d, p, w_by, w_bz, w_bx);
          max_dv2 = fabs(w0_(m,IVY,k,j,i)) + cf;
          cf = eos.IdealMHDFastSpeed(w_d, p, w_bz, w_bx, w_by);
          max_dv3 = fabs(w0_(m,IVZ,k,j,i)) + cf;
        } else {
          cf = eos.IdealMHDFastSpeed(w_d, w_bx, w_by, w_bz);
          max_dv1 = fabs(w0_(m,IVX,k,j,i)) + cf;
          cf = eos.IdealMHDFastSpeed(w_d, w_by, w_bz, w_bx);
          max_dv2 = fabs(w0_(m,IVY,k,j,i)) + cf;
          cf = eos.IdealMHDFastSpeed(w_d, w_bz, w_bx, w_by);
          max_dv3 = fabs(w0_(m,IVZ,k,j,i)) + cf;
        }
      }

      min_dt1 = fmin((mbsize.d_view(m).dx1/max_dv1), min_dt1);
      min_dt2 = fmin((mbsize.d_view(m).dx2/max_dv2), min_dt2);
      min_dt3 = fmin((mbsize.d_view(m).dx3/max_dv3), min_dt3);
    }, Kokkos::Min<Real>(dt1), Kokkos::Min<Real>(dt2),Kokkos::Min<Real>(dt3));
  }

  // compute minimum of dt1/dt2/dt3 for 1D/2D/3D problems
  dtnew = dt1;
  if (pmy_pack->pmesh->multi_d) { dtnew = std::min(dtnew, dt2); }
  if (pmy_pack->pmesh->three_d) { dtnew = std::min(dtnew, dt3); }

  // compute timestep for diffusion
  if (pcond != nullptr) {
    pcond->NewTimeStep(w0, peos->eos_data);
  }
  // compute source terms timestep
  if (psrc != nullptr) {
    psrc->NewTimeStep(w0, peos->eos_data);
  }

  return TaskStatus::complete;
}
} // namespace mhd
