//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file flux_correction_cc.cpp
//! \brief functions to pack/send and recv/unpack fluxes for cell-centered variables at
//! fine/coarse boundaries for the flux correction step.

#include <algorithm>
#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "bvals.hpp"

namespace {

bool IsCCFluxFace(const int n) {
  return (n < 16) || ((n >= 24) && (n < 32));
}

[[noreturn]] void FluxRegisterCCError(const char *message) {
  if (global_variable::my_rank == 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << std::endl
              << "Cell-centered level flux register: " << message << std::endl;
  }
#if MPI_PARALLEL_ENABLED
  if (global_variable::nranks > 1) {
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
#endif
  std::exit(EXIT_FAILURE);
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void MeshBoundaryValuesCC::InitializeFluxRegistersCC()
//! \brief Explicitly allocate persistent flux registers for level subcycling.
//!
//! Receive-side storage is keyed by (coarse MeshBlock, coarse-side neighbor buffer).
//! Send-side storage integrates fine fluxes whose coarse destination is on another rank.
//! Both are separate from the transient flux communication buffer and are allocated only
//! when requested, leaving the legacy global-timestep path allocation-free.

void MeshBoundaryValuesCC::InitializeFluxRegistersCC(int nvar) {
  if (nvar <= 0) {
    FluxRegisterCCError("the number of variables must be positive");
  }
  if (pmy_pack->pmesh->nmb_packs_thisrank != 1) {
    FluxRegisterCCError("level registers currently require one MeshBlockPack per rank");
  }
  if (!(pmy_pack->pmesh->multilevel) || pmy_pack->pmesh->adaptive) {
    FluxRegisterCCError("the initial implementation requires static mesh refinement");
  }
  if (flux_reg_nvar_ != 0) {
    if (flux_reg_nvar_ != nvar) {
      FluxRegisterCCError("cannot reinitialize registers with a different variable count");
    }
    return;
  }

  const int nmb = std::max(pmy_pack->nmb_thispack,
                           pmy_pack->pmesh->nmb_maxperrank);
  const int nnghbr = pmy_pack->pmb->nnghbr;
  for (int n=0; n<nnghbr; ++n) {
    if (IsCCFluxFace(n)) {
      const int recv_ndata = nvar*recvbuf[n].iflxc_ndat;
      if (recv_ndata > 0) {
        Kokkos::realloc(recvbuf[n].flux_reg, nmb, recv_ndata);
        Kokkos::deep_copy(recvbuf[n].flux_reg, 0.0);
      }
      const int send_ndata = nvar*sendbuf[n].iflxc_ndat;
      if (send_ndata > 0) {
        Kokkos::realloc(sendbuf[n].flux_reg, nmb, send_ndata);
        Kokkos::deep_copy(sendbuf[n].flux_reg, 0.0);
      }
    }
  }
  flux_reg_nvar_ = nvar;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MeshBoundaryValuesCC::ResetFluxRegistersCC()
//! \brief Clear registers for the coarse_level/coarse_level+1 interface pair.
//!
//! Receive registers live on coarse blocks.  Outbound registers live on fine blocks whose
//! coarse neighbor may be remote.  Registers owned by a coarser recursive caller are
//! deliberately preserved while this level takes either of its two child steps.

TaskStatus MeshBoundaryValuesCC::ResetFluxRegistersCC(int coarse_level) {
  if (flux_reg_nvar_ <= 0) {
    FluxRegisterCCError("ResetFluxRegistersCC called before initialization");
  }

  const int nmb = pmy_pack->nmb_thispack;
  const int nnghbr = pmy_pack->pmb->nnghbr;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mblev = pmy_pack->pmb->mb_lev;
  auto &sbuf = sendbuf;
  auto &rbuf = recvbuf;

  for (int n=0; n<nnghbr; ++n) {
    if (!IsCCFluxFace(n)) {
      continue;
    }
    const int recv_ndata = flux_reg_nvar_*rbuf[n].iflxc_ndat;
    const int send_ndata = flux_reg_nvar_*sbuf[n].iflxc_ndat;
    if (recv_ndata <= 0 && send_ndata <= 0) {
      continue;
    }
    auto recv_reg = rbuf[n].flux_reg;
    auto send_reg = sbuf[n].flux_reg;
    Kokkos::TeamPolicy<> policy(DevExeSpace(), nmb, Kokkos::AUTO);
    Kokkos::parallel_for("reset_cc_flux_register", policy,
    KOKKOS_LAMBDA(TeamMember_t tmember) {
      const int m = tmember.league_rank();
      if ((mblev.d_view(m) == coarse_level) &&
          (nghbr.d_view(m,n).gid >= 0) &&
          (nghbr.d_view(m,n).lev == coarse_level + 1)) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, recv_ndata),
        [&](const int q) {
          recv_reg(m,q) = 0.0;
        });
      } else if ((mblev.d_view(m) == coarse_level + 1) &&
                 (nghbr.d_view(m,n).gid >= 0) &&
                 (nghbr.d_view(m,n).lev == coarse_level)) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, send_ndata),
        [&](const int q) {
          send_reg(m,q) = 0.0;
        });
      }
    });
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MeshBoundaryValuesCC::AccumulateFluxRegistersCC()
//! \brief Add one RK2 stage to persistent coarse/fine flux mismatch registers.
//!
//! stage_weight includes both the RK quadrature coefficient and the active level dt.  For
//! Heun RK2 it is dt/2 for both stages (not Driver::beta).  An active coarse block adds
//! -weight*F_coarse to each finer-neighbor register.  An active fine block area-restricts
//! its flux exactly as PackAndSendFluxCC.  A local coarse destination is updated directly;
//! a remote destination is accumulated in an outbound register for exchange at the next
//! synchronization point.

TaskStatus MeshBoundaryValuesCC::AccumulateFluxRegistersCC(
    DvceFaceFld5D<Real> &flx, Real stage_weight) {
  if (flux_reg_nvar_ <= 0) {
    FluxRegisterCCError("AccumulateFluxRegistersCC called before initialization");
  }
  if (flx.x1f.extent_int(1) != flux_reg_nvar_) {
    FluxRegisterCCError("flux variable count does not match the allocated register");
  }
  if (!(stage_weight > 0.0)) {
    FluxRegisterCCError("stage weight must be positive and include the RK time weight");
  }
  if (pmy_pack->all_blocks_active || pmy_pack->active_level < 0) {
    FluxRegisterCCError("an explicit active logical level must be selected before accumulation");
  }
  if (pmy_pack->nmb_active <= 0) {
    return TaskStatus::complete;
  }

  const int nactive = pmy_pack->nmb_active;
  const int active_offset = pmy_pack->active_offset;
  const int nnghbr = pmy_pack->pmb->nnghbr;
  const int nvar = flux_reg_nvar_;
  const int cis = pmy_pack->pmesh->mb_indcs.cis;
  const int cjs = pmy_pack->pmesh->mb_indcs.cjs;
  const int cks = pmy_pack->pmesh->mb_indcs.cks;
  const int my_rank = global_variable::my_rank;
  const bool one_d = pmy_pack->pmesh->one_d;
  const bool two_d = pmy_pack->pmesh->two_d;

  auto active_lids = pmy_pack->active_lids.d_view;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mbgid = pmy_pack->pmb->mb_gid;
  auto &mblev = pmy_pack->pmb->mb_lev;
  auto &sbuf = sendbuf;
  auto &rbuf = recvbuf;

  Kokkos::TeamPolicy<> policy(DevExeSpace(), nactive*nnghbr*nvar, Kokkos::AUTO);
  Kokkos::parallel_for("accumulate_cc_flux_register", policy,
  KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int a = tmember.league_rank()/(nnghbr*nvar);
    const int n = (tmember.league_rank() - a*(nnghbr*nvar))/nvar;
    const int v = tmember.league_rank() - a*(nnghbr*nvar) - n*nvar;
    const int m = active_lids(active_offset + a);

    if ((nghbr.d_view(m,n).gid < 0) ||
        !((n < 16) || ((n >= 24) && (n < 32)))) {
      return;
    }

    const int neighbor_level = nghbr.d_view(m,n).lev;
    const int source_level = mblev.d_view(m);

    // Coarse contribution: initialize/increment this coarse-side key with -integral(Fc).
    if (neighbor_level == source_level + 1) {
      int il = rbuf[n].iflux_coar[0].bis;
      int iu = rbuf[n].iflux_coar[0].bie;
      int jl = rbuf[n].iflux_coar[0].bjs;
      int ju = rbuf[n].iflux_coar[0].bje;
      int kl = rbuf[n].iflux_coar[0].bks;
      int ku = rbuf[n].iflux_coar[0].bke;
      const int ni = iu - il + 1;
      const int nj = ju - jl + 1;
      const int nk = ku - kl + 1;

      if (n < 8) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nk*nj),
        [&](const int idx) {
          int k = idx/nj;
          int j = idx - k*nj;
          k += kl;
          j += jl;
          const int q = j-jl + nj*(k-kl + nk*v);
          rbuf[n].flux_reg(m,q) -= stage_weight*flx.x1f(m,v,k,j,il);
        });
      } else if (n < 16) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nk*ni),
        [&](const int idx) {
          int k = idx/ni;
          int i = idx - k*ni;
          k += kl;
          i += il;
          const int q = i-il + ni*(k-kl + nk*v);
          rbuf[n].flux_reg(m,q) -= stage_weight*flx.x2f(m,v,k,jl,i);
        });
      } else {
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nj*ni),
        [&](const int idx) {
          int j = idx/ni;
          int i = idx - j*ni;
          j += jl;
          i += il;
          const int q = i-il + ni*(j-jl + nj*v);
          rbuf[n].flux_reg(m,q) -= stage_weight*flx.x3f(m,v,kl,j,i);
        });
      }

    // Fine contribution: area-restrict onto the matching destination coarse-side key.
    } else if (neighbor_level + 1 == source_level) {
      const bool remote_destination = (nghbr.d_view(m,n).rank != my_rank);
      const int dm = remote_destination ? m : nghbr.d_view(m,n).gid - mbgid.d_view(0);
      const int dn = nghbr.d_view(m,n).dest;
      int il = sbuf[n].iflux_coar[0].bis;
      int iu = sbuf[n].iflux_coar[0].bie;
      int jl = sbuf[n].iflux_coar[0].bjs;
      int ju = sbuf[n].iflux_coar[0].bje;
      int kl = sbuf[n].iflux_coar[0].bks;
      int ku = sbuf[n].iflux_coar[0].bke;
      const int ni = iu - il + 1;
      const int nj = ju - jl + 1;
      const int nk = ku - kl + 1;

      if (n < 8) {
        const int fi = 2*il - cis;
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nk*nj),
        [&](const int idx) {
          int k = idx/nj;
          int j = idx - k*nj;
          k += kl;
          j += jl;
          const int fj = 2*j - cjs;
          const int fk = 2*k - cks;
          Real restricted_flux;
          if (one_d) {
            restricted_flux = flx.x1f(m,v,0,0,fi);
          } else if (two_d) {
            restricted_flux = 0.5*(flx.x1f(m,v,0,fj,fi) +
                                   flx.x1f(m,v,0,fj+1,fi));
          } else {
            restricted_flux = 0.25*(flx.x1f(m,v,fk  ,fj  ,fi) +
                                    flx.x1f(m,v,fk  ,fj+1,fi) +
                                    flx.x1f(m,v,fk+1,fj  ,fi) +
                                    flx.x1f(m,v,fk+1,fj+1,fi));
          }
          const int q = j-jl + nj*(k-kl + nk*v);
          if (remote_destination) {
            sbuf[n].flux_reg(m,q) += stage_weight*restricted_flux;
          } else {
            rbuf[dn].flux_reg(dm,q) += stage_weight*restricted_flux;
          }
        });
      } else if (n < 16) {
        const int fj = 2*jl - cjs;
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nk*ni),
        [&](const int idx) {
          int k = idx/ni;
          int i = idx - k*ni;
          k += kl;
          i += il;
          const int fi = 2*i - cis;
          const int fk = 2*k - cks;
          Real restricted_flux;
          if (two_d) {
            restricted_flux = 0.5*(flx.x2f(m,v,0,fj,fi) +
                                   flx.x2f(m,v,0,fj,fi+1));
          } else {
            restricted_flux = 0.25*(flx.x2f(m,v,fk  ,fj,fi  ) +
                                    flx.x2f(m,v,fk  ,fj,fi+1) +
                                    flx.x2f(m,v,fk+1,fj,fi  ) +
                                    flx.x2f(m,v,fk+1,fj,fi+1));
          }
          const int q = i-il + ni*(k-kl + nk*v);
          if (remote_destination) {
            sbuf[n].flux_reg(m,q) += stage_weight*restricted_flux;
          } else {
            rbuf[dn].flux_reg(dm,q) += stage_weight*restricted_flux;
          }
        });
      } else {
        const int fk = 2*kl - cks;
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nj*ni),
        [&](const int idx) {
          int j = idx/ni;
          int i = idx - j*ni;
          j += jl;
          i += il;
          const int fi = 2*i - cis;
          const int fj = 2*j - cjs;
          const Real restricted_flux =
              0.25*(flx.x3f(m,v,fk,fj  ,fi  ) + flx.x3f(m,v,fk,fj  ,fi+1) +
                    flx.x3f(m,v,fk,fj+1,fi  ) + flx.x3f(m,v,fk,fj+1,fi+1));
          const int q = i-il + ni*(j-jl + nj*v);
          if (remote_destination) {
            sbuf[n].flux_reg(m,q) += stage_weight*restricted_flux;
          } else {
            rbuf[dn].flux_reg(dm,q) += stage_weight*restricted_flux;
          }
        });
      }
    }
    tmember.team_barrier();
  });

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MeshBoundaryValuesCC::ExchangeFluxRegistersCC()
//! \brief Transfer integrated fine-side registers to remote coarse owners.
//!
//! Every MPI rank must call this routine at the synchronization point for the same level
//! pair, including ranks with no local coarse blocks.  The transient receive flux buffer is
//! used only as staging so an MPI receive cannot overwrite the coarse contribution already
//! stored in the persistent receive-side register.  All requests complete before return,
//! allowing the normal flux buffers and tags to be reused by the next stage.

TaskStatus MeshBoundaryValuesCC::ExchangeFluxRegistersCC(int coarse_level) {
  if (flux_reg_nvar_ <= 0) {
    FluxRegisterCCError("ExchangeFluxRegistersCC called before initialization");
  }

#if MPI_PARALLEL_ENABLED
  const int nmb = pmy_pack->nmb_thispack;
  const int nnghbr = pmy_pack->pmb->nnghbr;
  const int my_rank = global_variable::my_rank;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mblev = pmy_pack->pmb->mb_lev;

  Kokkos::fence();
  bool no_errors = true;

  // Post receives first.  One receive-side key has exactly one reciprocal fine sender.
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      if (IsCCFluxFace(n) &&
          mblev.h_view(m) == coarse_level &&
          nghbr.h_view(m,n).gid >= 0 &&
          nghbr.h_view(m,n).lev == coarse_level + 1 &&
          nghbr.h_view(m,n).rank != my_rank) {
        if (recvbuf[n].flux_req[m] != MPI_REQUEST_NULL) {
          FluxRegisterCCError("receive request was still active at register exchange");
        }
        const int data_size = flux_reg_nvar_*recvbuf[n].iflxc_ndat;
        auto recv_ptr = Kokkos::subview(recvbuf[n].flux, m, Kokkos::ALL);
        const int tag = CreateBvals_MPI_Tag(m, n);
        const int ierr = MPI_Irecv(recv_ptr.data(), data_size, MPI_ATHENA_REAL,
                                   nghbr.h_view(m,n).rank, tag, comm_flux,
                                   &(recvbuf[n].flux_req[m]));
        if (ierr != MPI_SUCCESS) {no_errors = false;}
      }
    }
  }

  // Send the time-integrated restricted fine flux directly from its persistent register.
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      if (IsCCFluxFace(n) &&
          mblev.h_view(m) == coarse_level + 1 &&
          nghbr.h_view(m,n).gid >= 0 &&
          nghbr.h_view(m,n).lev == coarse_level &&
          nghbr.h_view(m,n).rank != my_rank) {
        if (sendbuf[n].flux_req[m] != MPI_REQUEST_NULL) {
          FluxRegisterCCError("send request was still active at register exchange");
        }
        const int drank = nghbr.h_view(m,n).rank;
        const int lid = nghbr.h_view(m,n).gid - pmy_pack->pmesh->gids_eachrank[drank];
        const int tag = CreateBvals_MPI_Tag(lid, nghbr.h_view(m,n).dest);
        const int data_size = flux_reg_nvar_*sendbuf[n].iflxc_ndat;
        auto send_ptr = Kokkos::subview(sendbuf[n].flux_reg, m, Kokkos::ALL);
        const int ierr = MPI_Isend(send_ptr.data(), data_size, MPI_ATHENA_REAL, drank, tag,
                                   comm_flux, &(sendbuf[n].flux_req[m]));
        if (ierr != MPI_SUCCESS) {no_errors = false;}
      }
    }
  }
  if (!no_errors) {
    FluxRegisterCCError("MPI error while posting register exchange");
  }

  // Complete all incoming transfers before adding staging data on the device.
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      if (IsCCFluxFace(n) &&
          mblev.h_view(m) == coarse_level &&
          nghbr.h_view(m,n).gid >= 0 &&
          nghbr.h_view(m,n).lev == coarse_level + 1 &&
          nghbr.h_view(m,n).rank != my_rank) {
        const int ierr = MPI_Wait(&(recvbuf[n].flux_req[m]), MPI_STATUS_IGNORE);
        if (ierr != MPI_SUCCESS) {no_errors = false;}
      }
    }
  }
  if (!no_errors) {
    FluxRegisterCCError("MPI error while receiving register exchange");
  }

  auto &rbuf = recvbuf;
  for (int n=0; n<nnghbr; ++n) {
    if (!IsCCFluxFace(n)) {continue;}
    const int ndata = flux_reg_nvar_*rbuf[n].iflxc_ndat;
    auto staging = rbuf[n].flux;
    auto reg = rbuf[n].flux_reg;
    Kokkos::TeamPolicy<> policy(DevExeSpace(), nmb, Kokkos::AUTO);
    Kokkos::parallel_for("add_remote_cc_flux_register", policy,
    KOKKOS_LAMBDA(TeamMember_t tmember) {
      const int m = tmember.league_rank();
      if ((mblev.d_view(m) == coarse_level) &&
          (nghbr.d_view(m,n).gid >= 0) &&
          (nghbr.d_view(m,n).lev == coarse_level + 1) &&
          (nghbr.d_view(m,n).rank != my_rank)) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, ndata),
        [&](const int q) {
          reg(m,q) += staging(m,q);
        });
      }
    });
  }

  // Outbound storage remains immutable until every nonblocking send has completed.
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      if (IsCCFluxFace(n) &&
          mblev.h_view(m) == coarse_level + 1 &&
          nghbr.h_view(m,n).gid >= 0 &&
          nghbr.h_view(m,n).lev == coarse_level &&
          nghbr.h_view(m,n).rank != my_rank) {
        const int ierr = MPI_Wait(&(sendbuf[n].flux_req[m]), MPI_STATUS_IGNORE);
        if (ierr != MPI_SUCCESS) {no_errors = false;}
      }
    }
  }
  if (!no_errors) {
    FluxRegisterCCError("MPI error while completing register exchange");
  }
  Kokkos::fence();
#endif

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MeshBoundaryValuesCC::ApplyFluxRegistersCC()
//! \brief Materialize the integrated mismatch in flux_scratch and reflux coarse cells.
//!
//! Building all three face fields first and applying one divergence per cell avoids races
//! where corrections from two differently oriented coarse/fine faces meet at an edge or
//! corner.  flux_scratch may be the module's normal flux array because it is no longer
//! needed after both RK2 flux evaluations have been accumulated.

TaskStatus MeshBoundaryValuesCC::ApplyFluxRegistersCC(
    DvceArray5D<Real> &cons, DvceFaceFld5D<Real> &flux_scratch, int coarse_level) {
  if (flux_reg_nvar_ <= 0) {
    FluxRegisterCCError("ApplyFluxRegistersCC called before initialization");
  }
  if (cons.extent_int(1) != flux_reg_nvar_ ||
      flux_scratch.x1f.extent_int(1) != flux_reg_nvar_) {
    FluxRegisterCCError("state/flux variable count does not match the allocated register");
  }

  const int nmb = pmy_pack->nmb_thispack;
  const int nnghbr = pmy_pack->pmb->nnghbr;
  const int nvar = flux_reg_nvar_;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mblev = pmy_pack->pmb->mb_lev;
  auto &rbuf = recvbuf;

  // The flux arrays are disposable scratch at a synchronization point.  Clearing all
  // blocks is simpler and safer than leaving stale non-interface faces on another level.
  Kokkos::deep_copy(flux_scratch.x1f, 0.0);
  Kokkos::deep_copy(flux_scratch.x2f, 0.0);
  Kokkos::deep_copy(flux_scratch.x3f, 0.0);

  Kokkos::TeamPolicy<> policy(DevExeSpace(), nmb*nnghbr*nvar, Kokkos::AUTO);
  Kokkos::parallel_for("materialize_cc_flux_register", policy,
  KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = tmember.league_rank()/(nnghbr*nvar);
    const int n = (tmember.league_rank() - m*(nnghbr*nvar))/nvar;
    const int v = tmember.league_rank() - m*(nnghbr*nvar) - n*nvar;

    if ((mblev.d_view(m) != coarse_level) ||
        (nghbr.d_view(m,n).gid < 0) ||
        (nghbr.d_view(m,n).lev != coarse_level + 1) ||
        !((n < 16) || ((n >= 24) && (n < 32)))) {
      return;
    }

    const int il = rbuf[n].iflux_coar[0].bis;
    const int iu = rbuf[n].iflux_coar[0].bie;
    const int jl = rbuf[n].iflux_coar[0].bjs;
    const int ju = rbuf[n].iflux_coar[0].bje;
    const int kl = rbuf[n].iflux_coar[0].bks;
    const int ku = rbuf[n].iflux_coar[0].bke;
    const int ni = iu - il + 1;
    const int nj = ju - jl + 1;
    const int nk = ku - kl + 1;

    if (n < 8) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nk*nj),
      [&](const int idx) {
        int k = idx/nj;
        int j = idx - k*nj;
        k += kl;
        j += jl;
        const int q = j-jl + nj*(k-kl + nk*v);
        flux_scratch.x1f(m,v,k,j,il) = rbuf[n].flux_reg(m,q);
      });
    } else if (n < 16) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nk*ni),
      [&](const int idx) {
        int k = idx/ni;
        int i = idx - k*ni;
        k += kl;
        i += il;
        const int q = i-il + ni*(k-kl + nk*v);
        flux_scratch.x2f(m,v,k,jl,i) = rbuf[n].flux_reg(m,q);
      });
    } else {
      Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nj*ni),
      [&](const int idx) {
        int j = idx/ni;
        int i = idx - j*ni;
        j += jl;
        i += il;
        const int q = i-il + ni*(j-jl + nj*v);
        flux_scratch.x3f(m,v,kl,j,i) = rbuf[n].flux_reg(m,q);
      });
    }
    tmember.team_barrier();
  });

  auto &indcs = pmy_pack->pmesh->mb_indcs;
  const int is = indcs.is;
  const int ie = indcs.ie;
  const int js = indcs.js;
  const int je = indcs.je;
  const int ks = indcs.ks;
  const int ke = indcs.ke;
  const bool multi_d = pmy_pack->pmesh->multi_d;
  const bool three_d = pmy_pack->pmesh->three_d;
  auto &mbsize = pmy_pack->pmb->mb_size;
  auto u = cons;
  auto df1 = flux_scratch.x1f;
  auto df2 = flux_scratch.x2f;
  auto df3 = flux_scratch.x3f;

  par_for("apply_cc_flux_register", DevExeSpace(), 0, nmb-1, 0, nvar-1,
          ks, ke, js, je, is, ie,
  KOKKOS_LAMBDA(const int m, const int v, const int k, const int j, const int i) {
    if (mblev.d_view(m) == coarse_level) {
      Real div_mismatch =
          (df1(m,v,k,j,i+1) - df1(m,v,k,j,i))/mbsize.d_view(m).dx1;
      if (multi_d) {
        div_mismatch +=
            (df2(m,v,k,j+1,i) - df2(m,v,k,j,i))/mbsize.d_view(m).dx2;
      }
      if (three_d) {
        div_mismatch +=
            (df3(m,v,k+1,j,i) - df3(m,v,k,j,i))/mbsize.d_view(m).dx3;
      }
      u(m,v,k,j,i) -= div_mismatch;
    }
  });

  // Prevent accidental double application and make this synchronization API host-complete.
  (void) ResetFluxRegistersCC(coarse_level);
  Kokkos::fence();
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBoundaryValuesCC::PackAndSendFlux()
//! \brief Pack restricted fluxes of cell-centered variables at fine/coarse boundaries
//! into boundary buffers and send to neighbors for flux-correction step.  These fluxes
//! (e.g. for the conserved hydro variables) live at cell faces.
//!
//! This routine packs ALL the buffers on ALL the faces simultaneously for ALL the
//! MeshBlocks. Buffer data are then sent (via MPI) or copied directly for periodic or
//! block boundaries.

TaskStatus MeshBoundaryValuesCC::PackAndSendFluxCC(DvceFaceFld5D<Real> &flx) {
  // create local references for variables in kernel
  int nmb = pmy_pack->nmb_thispack;
  int nnghbr = pmy_pack->pmb->nnghbr;
  int nvar = flx.x1f.extent_int(1);  // TODO(@user): 2nd idx from L of in arr must be NVAR

  auto &cis = pmy_pack->pmesh->mb_indcs.cis;
  auto &cjs = pmy_pack->pmesh->mb_indcs.cjs;
  auto &cks = pmy_pack->pmesh->mb_indcs.cks;

  int my_rank = global_variable::my_rank;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mbgid = pmy_pack->pmb->mb_gid;
  auto &mblev = pmy_pack->pmb->mb_lev;
  auto &sbuf = sendbuf;
  auto &rbuf = recvbuf;
  auto &one_d = pmy_pack->pmesh->one_d;
  auto &two_d = pmy_pack->pmesh->two_d;

  // Outer loop over (# of MeshBlocks)*(# of neighbors)*(# of variables)
  Kokkos::TeamPolicy<> policy(DevExeSpace(), (nmb*nnghbr*nvar), Kokkos::AUTO);
  Kokkos::parallel_for("RecvBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = (tmember.league_rank())/(nnghbr*nvar);
    const int n = (tmember.league_rank() - m*(nnghbr*nvar))/nvar;
    const int v = (tmember.league_rank() - m*(nnghbr*nvar) - n*nvar);

    // Note send buffer flux indices are for the coarse mesh
    int il = sbuf[n].iflux_coar[0].bis;
    int iu = sbuf[n].iflux_coar[0].bie;
    int jl = sbuf[n].iflux_coar[0].bjs;
    int ju = sbuf[n].iflux_coar[0].bje;
    int kl = sbuf[n].iflux_coar[0].bks;
    int ku = sbuf[n].iflux_coar[0].bke;
    const int ni = iu - il + 1;
    const int nj = ju - jl + 1;
    const int nk = ku - kl + 1;
    const int nji  = nj*ni;
    const int nkj  = nk*nj;
    const int nki  = nk*ni;

    // indices of recv'ing (destination) MB and buffer: MB IDs are stored sequentially
    // in MeshBlockPacks, so array index equals (target_id - first_id)
    int dm = nghbr.d_view(m,n).gid - mbgid.d_view(0);
    int dn = nghbr.d_view(m,n).dest;

    // only pack buffers when neighbor is at coarser level
    if ((nghbr.d_view(m,n).gid >=0) && (nghbr.d_view(m,n).lev < mblev.d_view(m))) {
      // x1faces
      if (n<8) {
        // i-index is fixed for flux correction on x1faces
        int fi = 2*il - cis;
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkj), [&](const int idx) {
          int k = idx / nj;
          int j = (idx - k * nj) + jl;
          k += kl;
          int fj = 2*j - cjs;
          int fk = 2*k - cks;
          Real rflx;
          if (one_d) {
            rflx = flx.x1f(m,v,0,0,fi);
          } else if (two_d) {
            rflx = 0.5*(flx.x1f(m,v,0,fj,fi) + flx.x1f(m,v,0,fj+1,fi));
          } else {
            rflx = 0.25*(flx.x1f(m,v,fk  ,fj,fi) + flx.x1f(m,v,fk  ,fj+1,fi) +
                         flx.x1f(m,v,fk+1,fj,fi) + flx.x1f(m,v,fk+1,fj+1,fi));
          }
          // copy directly into recv buffer if MeshBlocks on same rank
          if (nghbr.d_view(m,n).rank == my_rank) {
            rbuf[dn].flux(dm, (j-jl + nj*(k-kl + nk*v)) ) = rflx;
          // else copy into send buffer for MPI communication below
          } else {
            sbuf[n].flux(m, (j-jl + nj*(k-kl + nk*v)) ) = rflx;
          }
        });

      // x2faces
      } else if (n<16) {
        // j-index is fixed for flux correction on x2faces
        int fj = 2*jl - cjs;
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nki), [&](const int idx) {
          int k = idx / ni;
          int i = (idx - k * ni) + il;
          k += kl;
          int fi = 2*i - cis;
          int fk = 2*k - cks;
          Real rflx;
          if (two_d) {
            rflx = 0.5*(flx.x2f(m,v,0,fj,fi) + flx.x2f(m,v,0,fj,fi+1));
          } else {
            rflx = 0.25*(flx.x2f(m,v,fk  ,fj,fi) + flx.x2f(m,v,fk  ,fj,fi+1) +
                         flx.x2f(m,v,fk+1,fj,fi) + flx.x2f(m,v,fk+1,fj,fi+1));
          }
          // copy directly into recv buffer if MeshBlocks on same rank
          if (nghbr.d_view(m,n).rank == my_rank) {
            rbuf[dn].flux(dm, (i-il + ni*(k-kl + nk*v)) ) = rflx;
          // else copy into send buffer for MPI communication below
          } else {
            sbuf[n].flux(m, (i-il + ni*(k-kl + nk*v)) ) = rflx;
          }
        });

      // x3faces
      } else if ((n>=24) && (n<32)) {
        // k-index is fixed for flux correction on x3faces
        int fk = 2*kl - cks;
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nji), [&](const int idx) {
          int j = idx / ni;
          int i = (idx - j * ni) + il;
          j += jl;
          int fi = 2*i - cis;
          int fj = 2*j - cjs;
          Real rflx = 0.25*(flx.x3f(m,v,fk,fj  ,fi) + flx.x3f(m,v,fk,fj  ,fi+1) +
                            flx.x3f(m,v,fk,fj+1,fi) + flx.x3f(m,v,fk,fj+1,fi+1));
          // copy directly into recv buffer if MeshBlocks on same rank
          if (nghbr.d_view(m,n).rank == my_rank) {
            rbuf[dn].flux(dm, (i-il + ni*(j-jl + nj*v)) ) = rflx;
          // else copy into send buffer for MPI communication below
          } else {
            sbuf[n].flux(m, (i-il + ni*(j-jl + nj*v)) ) = rflx;
          }
        });
      }
    }  // end if-neighbor-exists block
    tmember.team_barrier();
  });  // end par_for_outer

#if MPI_PARALLEL_ENABLED
  // Send boundary buffer to neighboring MeshBlocks using MPI
  // Sends only occur to neighbors on FACES at a COARSER level
  Kokkos::fence();
  bool no_errors=true;
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      if ( (nghbr.h_view(m,n).gid >=0) &&
           (nghbr.h_view(m,n).lev < mblev.h_view(m)) &&
           ((n<16) || ((n>=24) && (n<32))) ) {
        // index and rank of destination Neighbor
        int dn = nghbr.h_view(m,n).dest;
        int drank = nghbr.h_view(m,n).rank;

        if (drank != my_rank) {
          // create tag using local ID and buffer index of *receiving* MeshBlock
          int lid = nghbr.h_view(m,n).gid - pmy_pack->pmesh->gids_eachrank[drank];
          int tag = CreateBvals_MPI_Tag(lid, dn);

          // get ptr to send buffer for fluxes
          int data_size = nvar*(sendbuf[n].iflxc_ndat);
          auto send_ptr = Kokkos::subview(sendbuf[n].flux, m, Kokkos::ALL);

          int ierr = MPI_Isend(send_ptr.data(), data_size, MPI_ATHENA_REAL, drank, tag,
                               comm_flux, &(sendbuf[n].flux_req[m]));
          if (ierr != MPI_SUCCESS) {no_errors=false;}
        }
      }
    }
  }
  // Quit if MPI error detected
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
       << std::endl << "MPI error in posting sends" << std::endl;
    std::exit(EXIT_FAILURE);
  }
#endif
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void RecvBuffers()
//! \brief Unpack boundary buffers for flux correction of CC variables.

TaskStatus MeshBoundaryValuesCC::RecvAndUnpackFluxCC(DvceFaceFld5D<Real> &flx) {
  // create local references for variables in kernel
  int nmb = pmy_pack->nmb_thispack;
  int nnghbr = pmy_pack->pmb->nnghbr;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mblev = pmy_pack->pmb->mb_lev;
  auto &rbuf = recvbuf;
#if MPI_PARALLEL_ENABLED
  //----- STEP 1: check that recv boundary buffer communications have all completed
  // receives only occur for neighbors on faces at a FINER level

  bool bflag = false;
  bool no_errors=true;
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      if ( (nghbr.h_view(m,n).gid >=0) &&
           (nghbr.h_view(m,n).lev > mblev.h_view(m)) &&
           ((n<16) || ((n>=24) && (n<32))) ) {
        if (nghbr.h_view(m,n).rank != global_variable::my_rank) {
          int test;
          int ierr = MPI_Test(&(rbuf[n].flux_req[m]), &test, MPI_STATUS_IGNORE);
          if (ierr != MPI_SUCCESS) {no_errors=false;}
          if (!(static_cast<bool>(test))) {
            bflag = true;
          }
        }
      }
    }
  }
  // Quit if MPI error detected
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "MPI error in testing non-blocking receives"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  // exit if recv boundary buffer communications have not completed
  if (bflag) {return TaskStatus::incomplete;}
#endif

  //----- STEP 2: buffers have all completed, so unpack

  int nvar = flx.x1f.extent_int(1); // TODO(@user): 2nd idx from L of in arr must be NVAR

  // Outer loop over (# of MeshBlocks)*(# of neighbors)*(# of variables)
  Kokkos::TeamPolicy<> policy(DevExeSpace(), (nmb*nnghbr*nvar), Kokkos::AUTO);
  Kokkos::parallel_for("RecvBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = (tmember.league_rank())/(nnghbr*nvar);
    const int n = (tmember.league_rank() - m*(nnghbr*nvar))/nvar;
    const int v = (tmember.league_rank() - m*(nnghbr*nvar) - n*nvar);

    // Recv buffer flux indices are for the regular mesh
    int il = rbuf[n].iflux_coar[0].bis;
    int iu = rbuf[n].iflux_coar[0].bie;
    int jl = rbuf[n].iflux_coar[0].bjs;
    int ju = rbuf[n].iflux_coar[0].bje;
    int kl = rbuf[n].iflux_coar[0].bks;
    int ku = rbuf[n].iflux_coar[0].bke;
    const int ni = iu - il + 1;
    const int nj = ju - jl + 1;
    const int nk = ku - kl + 1;
    const int nji  = nj*ni;
    const int nkj  = nk*nj;
    const int nki  = nk*ni;

    // only unpack buffers for faces when neighbor is at finer level
    if ((nghbr.d_view(m,n).gid >=0) && (nghbr.d_view(m,n).lev > mblev.d_view(m))) {
      //x1 faces
      if (n<8) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkj), [&](const int idx) {
          int k = idx / nj;
          int j = (idx - k * nj) + jl;
          k += kl;
          flx.x1f(m,v,k,j,il) = rbuf[n].flux(m,(j-jl + nj*(k-kl + nk*v)));
        });
      // x2faces
      } else if (n<16) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nki), [&](const int idx) {
          int k = idx / ni;
          int i = (idx - k * ni) + il;
          k += kl;
          flx.x2f(m,v,k,jl,i) = rbuf[n].flux(m,(i-il + ni*(k-kl + nk*v)));
        });
      // x3faces
      } else if ((n>=24) && (n<32)) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nji), [&](const int idx) {
          int j = idx / ni;
          int i = (idx - j * ni) + il;
          j += jl;
          flx.x3f(m,v,kl,j,i) = rbuf[n].flux(m,(i-il + ni*(j-jl + nj*v)));
        });
      }
    }  // end if-neighbor-exists block
    tmember.team_barrier();
  });  // end par_for_outer

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn  void BoundaryValuesCC::InitRecvFlux
//! \brief Posts non-blocking receives (with MPI) for boundary communication of fluxes of
//! cell-centered variables, which are communicated at FACES of MeshBlocks at the SAME
//! levels.  This is different than for fluxes of face-centered vars.

TaskStatus MeshBoundaryValuesCC::InitFluxRecv(const int nvars) {
  // Level subcycling integrates fine fluxes in persistent registers and exchanges them
  // once per level-pair synchronization.  Posting the legacy per-stage receive here would
  // leave it unmatched because MHD::SendFlux does not use PackAndSendFluxCC in this mode.
  if (!(pmy_pack->all_blocks_active)) {
    return TaskStatus::complete;
  }
#if MPI_PARALLEL_ENABLED
  int &nmb = pmy_pack->nmb_thispack;
  int &nnghbr = pmy_pack->pmb->nnghbr;
  auto &nghbr = pmy_pack->pmb->nghbr;

  // Initialize communications of fluxes
  bool no_errors=true;
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      // only post receives for neighbors on FACES at FINER level
      // this is the only thing different from BoundaryValuesFC::InitRecvFlux()
      if ( (nghbr.h_view(m,n).gid >=0) &&
           (nghbr.h_view(m,n).lev > pmy_pack->pmb->mb_lev.h_view(m)) &&
           ((n<16) || ((n>=24) && (n<32))) ) {
        // rank of destination buffer
        int drank = nghbr.h_view(m,n).rank;

        // post non-blocking receive if neighboring MeshBlock on a different rank
        if (drank != global_variable::my_rank) {
          // create tag using local ID and buffer index of *receiving* MeshBlock
          int tag = CreateBvals_MPI_Tag(m, n);

          // calculate amount of data to be passed, get pointer to variables
          int data_size = nvars*(recvbuf[n].iflxc_ndat);
          auto recv_ptr = Kokkos::subview(recvbuf[n].flux, m, Kokkos::ALL);

          // Post non-blocking receive for this buffer on this MeshBlock
          int ierr = MPI_Irecv(recv_ptr.data(), data_size, MPI_ATHENA_REAL, drank, tag,
                               comm_flux, &(recvbuf[n].flux_req[m]));
          if (ierr != MPI_SUCCESS) {no_errors=false;}
        }
      }
    }
  }
  // Quit if MPI error detected
  if (!(no_errors)) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
       << std::endl << "MPI error in posting non-blocking receives" << std::endl;
    std::exit(EXIT_FAILURE);
  }
#endif
  return TaskStatus::complete;
}
