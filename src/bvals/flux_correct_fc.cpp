//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file flux_correction_fc.cpp
//! \brief functions to pack/send and recv/unpack fluxes (emfs) for face-centered fields
//! (magnetic fields) at fine/coarse boundaries for the flux correction step.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "bvals.hpp"

namespace {

bool IsFCFluxBuffer(const int n) {
  return n >= 0 && n < 48;
}

[[noreturn]] void FluxRegisterFCError(const char *message) {
  if (global_variable::my_rank == 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << std::endl
              << "Face-centered level EMF register: " << message << std::endl;
  }
  std::exit(EXIT_FAILURE);
}

void ValidateFluxRegistersFC(MeshBoundaryBuffer *rbuf, const int nnghbr,
                             const int nmb) {
  for (int n=0; n<nnghbr && IsFCFluxBuffer(n); ++n) {
    const int ndata = 3*rbuf[n].iflxc_ndat;
    if (ndata > 0 &&
        (rbuf[n].flux_reg.extent_int(0) < nmb ||
         rbuf[n].flux_reg.extent_int(1) < ndata)) {
      FluxRegisterFCError("register operation called before initialization");
    }
  }
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void MeshBoundaryValuesFC::InitializeFluxRegistersFC()
//! \brief Explicitly allocate persistent receive-side EMF registers for subcycling.
//!
//! Registers are keyed by (coarse MeshBlock, coarse-side neighbor buffer).  They are
//! deliberately separate from the transient same-level EMF communication buffers, and
//! are allocated only when level subcycling opts in.

void MeshBoundaryValuesFC::InitializeFluxRegistersFC() {
  if (global_variable::nranks != 1 || pmy_pack->pmesh->nmb_packs_thisrank != 1) {
    FluxRegisterFCError("the initial implementation requires one rank and one MeshBlockPack");
  }
  if (!(pmy_pack->pmesh->multilevel) || pmy_pack->pmesh->adaptive) {
    FluxRegisterFCError("the initial implementation requires static mesh refinement");
  }

  const int nmb = std::max(pmy_pack->nmb_thispack,
                           pmy_pack->pmesh->nmb_maxperrank);
  const int nnghbr = pmy_pack->pmb->nnghbr;

  bool any_initialized = false;
  for (int n=0; n<nnghbr && IsFCFluxBuffer(n); ++n) {
    any_initialized = any_initialized || (recvbuf[n].flux_reg.extent_int(0) > 0);
  }
  if (any_initialized) {
    ValidateFluxRegistersFC(recvbuf, nnghbr, nmb);
    return;
  }

  for (int n=0; n<nnghbr && IsFCFluxBuffer(n); ++n) {
    const int ndata = 3*recvbuf[n].iflxc_ndat;
    if (ndata > 0) {
      Kokkos::realloc(recvbuf[n].flux_reg, nmb, ndata);
      Kokkos::deep_copy(recvbuf[n].flux_reg, 0.0);
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MeshBoundaryValuesFC::ResetFluxRegistersFC()
//! \brief Clear registers owned by coarse_level at interfaces with coarse_level+1.
//!
//! Registers owned by a coarser recursive caller remain intact while this level advances
//! through either of its child steps.

TaskStatus MeshBoundaryValuesFC::ResetFluxRegistersFC(int coarse_level) {
  const int nmb = pmy_pack->nmb_thispack;
  const int nnghbr = pmy_pack->pmb->nnghbr;
  ValidateFluxRegistersFC(recvbuf, nnghbr, nmb);

  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mblev = pmy_pack->pmb->mb_lev;
  auto &rbuf = recvbuf;

  for (int n=0; n<nnghbr && IsFCFluxBuffer(n); ++n) {
    const int ndata = 3*rbuf[n].iflxc_ndat;
    if (ndata <= 0) {
      continue;
    }
    auto reg = rbuf[n].flux_reg;
    Kokkos::TeamPolicy<> policy(DevExeSpace(), nmb, Kokkos::AUTO);
    Kokkos::parallel_for("reset_fc_emf_register", policy,
    KOKKOS_LAMBDA(TeamMember_t tmember) {
      const int m = tmember.league_rank();
      if ((mblev.d_view(m) == coarse_level) &&
          (nghbr.d_view(m,n).gid >= 0) &&
          (nghbr.d_view(m,n).lev == coarse_level + 1)) {
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, ndata),
        [&](const int q) {
          reg(m,q) = 0.0;
        });
      }
    });
  }
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MeshBoundaryValuesFC::AccumulateFluxRegistersFC()
//! \brief Add one RK2 stage to persistent coarse/fine EMF mismatch registers.
//!
//! stage_weight includes both the RK quadrature coefficient and the active-level dt; it
//! is dt/2 for each stage of Heun RK2.  Active coarse blocks add -weight*E_coarse to
//! their own coarse-side keys.  Active fine blocks restrict edge EMFs exactly as the
//! legacy PackAndSendFluxFC path and atomically add +weight*E_fine to the destination
//! coarse key.  The register therefore stores the time-integrated fine-minus-coarse EMF.

TaskStatus MeshBoundaryValuesFC::AccumulateFluxRegistersFC(
    DvceEdgeFld4D<Real> &flx, Real stage_weight) {
  const int nmb = pmy_pack->nmb_thispack;
  const int nnghbr = pmy_pack->pmb->nnghbr;
  ValidateFluxRegistersFC(recvbuf, nnghbr, nmb);
  if (!(stage_weight > 0.0) || !std::isfinite(stage_weight)) {
    FluxRegisterFCError("stage weight must be positive, finite, and include the RK time weight");
  }
  if (pmy_pack->all_blocks_active || pmy_pack->active_level < 0) {
    FluxRegisterFCError("an explicit active logical level must be selected before accumulation");
  }
  if (pmy_pack->nmb_active <= 0) {
    return TaskStatus::complete;
  }

  const int nactive = pmy_pack->nmb_active;
  const int active_offset = pmy_pack->active_offset;
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

  Kokkos::TeamPolicy<> policy(DevExeSpace(), 3*nactive*nnghbr, Kokkos::AUTO);
  Kokkos::parallel_for("accumulate_fc_emf_register", policy,
  KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int a = tmember.league_rank()/(3*nnghbr);
    const int n = (tmember.league_rank() - a*(3*nnghbr))/3;
    const int v = tmember.league_rank() - a*(3*nnghbr) - 3*n;
    const int m = active_lids(active_offset + a);

    if (n >= 48 || nghbr.d_view(m,n).gid < 0) {
      return;
    }
    const bool component_used =
        ((n < 8) && (v == 1 || v == 2)) ||
        ((n >= 8 && n < 16) && (v == 0 || v == 2)) ||
        ((n >= 16 && n < 24) && v == 2) ||
        ((n >= 24 && n < 32) && (v == 0 || v == 1)) ||
        ((n >= 32 && n < 40) && v == 1) ||
        ((n >= 40) && v == 0);
    if (!component_used) {
      return;
    }

    const int source_level = mblev.d_view(m);
    const int neighbor_level = nghbr.d_view(m,n).lev;
    const bool source_is_coarse = (neighbor_level == source_level + 1);
    const bool source_is_fine = (neighbor_level + 1 == source_level);
    if (!(source_is_coarse || source_is_fine)) {
      return;
    }
    if (source_is_fine && nghbr.d_view(m,n).rank != my_rank) {
      return;
    }

    int dm = m;
    int dn = n;
    int il, iu, jl, ju, kl, ku, ndat;
    if (source_is_coarse) {
      il = rbuf[n].iflux_coar[v].bis;
      iu = rbuf[n].iflux_coar[v].bie;
      jl = rbuf[n].iflux_coar[v].bjs;
      ju = rbuf[n].iflux_coar[v].bje;
      kl = rbuf[n].iflux_coar[v].bks;
      ku = rbuf[n].iflux_coar[v].bke;
      ndat = rbuf[n].iflxc_ndat;
    } else {
      dm = nghbr.d_view(m,n).gid - mbgid.d_view(0);
      dn = nghbr.d_view(m,n).dest;
      il = sbuf[n].iflux_coar[v].bis;
      iu = sbuf[n].iflux_coar[v].bie;
      jl = sbuf[n].iflux_coar[v].bjs;
      ju = sbuf[n].iflux_coar[v].bje;
      kl = sbuf[n].iflux_coar[v].bks;
      ku = sbuf[n].iflux_coar[v].bke;
      ndat = sbuf[n].iflxc_ndat;
    }
    const int ni = iu - il + 1;
    const int nj = ju - jl + 1;
    const int nk = ku - kl + 1;
    const int nji = nj*ni;
    const int nkj = nk*nj;
    const int nki = nk*ni;
    const Real coefficient = source_is_coarse ? -stage_weight : stage_weight;

    // x1faces: only x2e and x3e are tangential to the interface.
    if (n < 8) {
      const int fi = 2*il - cis;
      Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkj),
      [&](const int idx) {
        int k = idx/nj;
        int j = idx - k*nj;
        k += kl;
        j += jl;
        const int fj = 2*j - cjs;
        const int fk = 2*k - cks;
        Real emf;
        if (v == 1) {
          if (source_is_coarse) {
            emf = flx.x2e(m,k,j,il);
          } else if (one_d) {
            emf = flx.x2e(m,0,0,fi);
          } else if (two_d) {
            emf = 0.5*(flx.x2e(m,0,fj,fi) + flx.x2e(m,0,fj+1,fi));
          } else {
            emf = 0.5*(flx.x2e(m,fk,fj,fi) + flx.x2e(m,fk,fj+1,fi));
          }
        } else {
          if (source_is_coarse) {
            emf = flx.x3e(m,k,j,il);
          } else if (one_d) {
            emf = flx.x3e(m,0,0,fi);
          } else if (two_d) {
            emf = flx.x3e(m,0,fj,fi);
          } else {
            emf = 0.5*(flx.x3e(m,fk,fj,fi) + flx.x3e(m,fk+1,fj,fi));
          }
        }
        const int q = ndat*v + (j-jl + nj*(k-kl));
        Kokkos::atomic_add(&rbuf[dn].flux_reg(dm,q), coefficient*emf);
      });

    // x2faces: only x1e and x3e are tangential to the interface.
    } else if (n < 16) {
      const int fj = 2*jl - cjs;
      Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nki),
      [&](const int idx) {
        int k = idx/ni;
        int i = idx - k*ni;
        k += kl;
        i += il;
        const int fk = 2*k - cks;
        const int fi = 2*i - cis;
        Real emf;
        if (v == 0) {
          if (source_is_coarse) {
            emf = flx.x1e(m,k,jl,i);
          } else if (two_d) {
            emf = 0.5*(flx.x1e(m,0,fj,fi) + flx.x1e(m,0,fj,fi+1));
          } else {
            emf = 0.5*(flx.x1e(m,fk,fj,fi) + flx.x1e(m,fk,fj,fi+1));
          }
        } else {
          if (source_is_coarse) {
            emf = flx.x3e(m,k,jl,i);
          } else if (two_d) {
            emf = flx.x3e(m,0,fj,fi);
          } else {
            emf = 0.5*(flx.x3e(m,fk,fj,fi) + flx.x3e(m,fk+1,fj,fi));
          }
        }
        const int q = ndat*v + i-il + ni*(k-kl);
        Kokkos::atomic_add(&rbuf[dn].flux_reg(dm,q), coefficient*emf);
      });

    // x1x2 edges: only x3e lives on the edge.
    } else if (n < 24) {
      const int fi = 2*il - cis;
      const int fj = 2*jl - cjs;
      Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nk),
      [&](const int idx) {
        const int k = idx + kl;
        const int fk = 2*k - cks;
        Real emf;
        if (source_is_coarse) {
          emf = flx.x3e(m,k,jl,il);
        } else if (two_d) {
          emf = flx.x3e(m,0,fj,fi);
        } else {
          emf = 0.5*(flx.x3e(m,fk,fj,fi) + flx.x3e(m,fk+1,fj,fi));
        }
        const int q = ndat*v + (k-kl);
        Kokkos::atomic_add(&rbuf[dn].flux_reg(dm,q), coefficient*emf);
      });

    // x3faces: only x1e and x2e are tangential to the interface.
    } else if (n < 32) {
      const int fk = 2*kl - cks;
      Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nji),
      [&](const int idx) {
        int j = idx/ni;
        int i = idx - j*ni;
        j += jl;
        i += il;
        const int fi = 2*i - cis;
        const int fj = 2*j - cjs;
        Real emf;
        if (v == 0) {
          emf = source_is_coarse
              ? flx.x1e(m,kl,j,i)
              : 0.5*(flx.x1e(m,fk,fj,fi) + flx.x1e(m,fk,fj,fi+1));
        } else {
          emf = source_is_coarse
              ? flx.x2e(m,kl,j,i)
              : 0.5*(flx.x2e(m,fk,fj,fi) + flx.x2e(m,fk,fj+1,fi));
        }
        const int q = ndat*v + i-il + ni*(j-jl);
        Kokkos::atomic_add(&rbuf[dn].flux_reg(dm,q), coefficient*emf);
      });

    // x3x1 edges: only x2e lives on the edge.
    } else if (n < 40) {
      const int fi = 2*il - cis;
      const int fk = 2*kl - cks;
      Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nj),
      [&](const int idx) {
        const int j = idx + jl;
        const int fj = 2*j - cjs;
        const Real emf = source_is_coarse
            ? flx.x2e(m,kl,j,il)
            : 0.5*(flx.x2e(m,fk,fj,fi) + flx.x2e(m,fk,fj+1,fi));
        const int q = ndat*v + (j-jl);
        Kokkos::atomic_add(&rbuf[dn].flux_reg(dm,q), coefficient*emf);
      });

    // x2x3 edges: only x1e lives on the edge.
    } else {
      const int fj = 2*jl - cjs;
      const int fk = 2*kl - cks;
      Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, ni),
      [&](const int idx) {
        const int i = idx + il;
        const int fi = 2*i - cis;
        const Real emf = source_is_coarse
            ? flx.x1e(m,kl,jl,i)
            : 0.5*(flx.x1e(m,fk,fj,fi) + flx.x1e(m,fk,fj,fi+1));
        const int q = ndat*v + i-il;
        Kokkos::atomic_add(&rbuf[dn].flux_reg(dm,q), coefficient*emf);
      });
    }
    tmember.team_barrier();
  });

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn void MeshBoundaryValuesFC::PackAndSendFluxFC()
//! \brief Pack restricted fluxes of face-centered fields at fine/coarse boundaries
//! into boundary buffers and send to neighbors for flux-correction step. These fluxes
//! (e.g. EMFs) live at cell edges.
//!
//! This routine packs ALL the buffers on ALL the faces simultaneously for ALL the
//! MeshBlocks. Buffer data are then sent (via MPI) or copied directly for periodic or
//! block boundaries.

TaskStatus MeshBoundaryValuesFC::PackAndSendFluxFC(DvceEdgeFld4D<Real> &flx) {
  return PackAndSendFluxFC(flx, false);
}

//----------------------------------------------------------------------------------------
//! \brief Pack edge fields for the level-local path.
//!
//! With same_level_only=true this performs only the EMF synchronization required before
//! CT on the active level.  Coarse/fine EMFs are deliberately excluded: subcycling must
//! combine those through a time-integrated EMF register instead of replacing a coarse
//! RK-stage value with a fine RK-stage value at a different physical time.

TaskStatus MeshBoundaryValuesFC::PackAndSendFluxFC(DvceEdgeFld4D<Real> &flx,
                                                    bool same_level_only) {
  // create local references for variables in kernel
  int nmb = pmy_pack->nmb_thispack;
  int nnghbr = pmy_pack->pmb->nnghbr;

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
  const int active_level = pmy_pack->active_level;

  // Outer loop over (# of MeshBlocks)*(# of neighbors)*(3 field components)
  Kokkos::TeamPolicy<> policy(DevExeSpace(), (3*nmb*nnghbr), Kokkos::AUTO);
  Kokkos::parallel_for("RecvBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = (tmember.league_rank())/(3*nnghbr);
    const int n = (tmember.league_rank() - m*(3*nnghbr))/3;
    const int v = (tmember.league_rank() - m*(3*nnghbr) - 3*n);

    // Legacy mode loads same/coarser destinations.  The subcycling path loads only
    // same-level buffers sourced by the active level.
    const bool source_is_active = (mblev.d_view(m) == active_level);
    const bool eligible_neighbor = same_level_only
        ? (source_is_active && nghbr.d_view(m,n).lev == mblev.d_view(m))
        : (nghbr.d_view(m,n).lev <= mblev.d_view(m));
    if ((nghbr.d_view(m,n).gid >= 0) && eligible_neighbor) {
      // if neighbor is at coarser level, use cindices to pack buffer
      // Note indices can be different for each component of flux
      int il, iu, jl, ju, kl, ku, ndat;
      if (nghbr.d_view(m,n).lev < mblev.d_view(m)) {
        il = sbuf[n].iflux_coar[v].bis;
        iu = sbuf[n].iflux_coar[v].bie;
        jl = sbuf[n].iflux_coar[v].bjs;
        ju = sbuf[n].iflux_coar[v].bje;
        kl = sbuf[n].iflux_coar[v].bks;
        ku = sbuf[n].iflux_coar[v].bke;
        ndat = sbuf[n].iflxc_ndat;
      // if neighbor is at same level, use sindices to pack buffer
      } else if (nghbr.d_view(m,n).lev == mblev.d_view(m)) {
        il = sbuf[n].iflux_same[v].bis;
        iu = sbuf[n].iflux_same[v].bie;
        jl = sbuf[n].iflux_same[v].bjs;
        ju = sbuf[n].iflux_same[v].bje;
        kl = sbuf[n].iflux_same[v].bks;
        ku = sbuf[n].iflux_same[v].bke;
        ndat = sbuf[n].iflxs_ndat;
      }
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

      // x1faces (only load x2e and x3e)
      if (n<8) {
        // i-index is fixed for flux correction on x1faces
        const int fi = 2*il - cis;
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nkj), [&](const int idx) {
          int k = idx / nj;
          int j = (idx - k * nj) + jl;
          k += kl;
          int fj = 2*j - cjs;
          int fk = 2*k - cks;
          if (v==1) {
            Real rflx;
            // if neighbor is at same level, load x2e directly
            if (nghbr.d_view(m,n).lev == mblev.d_view(m)) {
              rflx = flx.x2e(m,k,j,il);
            // if neighbor is at coarser level, restrict x2e
            } else {
              if (one_d) {
                rflx = flx.x2e(m,0,0,fi);
              } else if (two_d) {
                rflx = 0.5*(flx.x2e(m,0,fj,fi) + flx.x2e(m,0,fj+1,fi));
              } else {
                rflx = 0.5*(flx.x2e(m,fk,fj,fi) + flx.x2e(m,fk,fj+1,fi));
              }
            }
            // copy directly into recv buffer if MeshBlocks on same rank
            if (nghbr.d_view(m,n).rank == my_rank) {
              rbuf[dn].flux(dm, ndat*v + (j-jl + nj*(k-kl))) = rflx;
            // else copy into send buffer for MPI communication below
            } else {
              sbuf[n].flux(m, ndat*v + (j-jl + nj*(k-kl))) = rflx;
            }
          } else if (v==2) {
            Real rflx;
            // if neighbor is at same level, load x3e directly
            if (nghbr.d_view(m,n).lev == mblev.d_view(m)) {
              rflx = flx.x3e(m,k,j,il);
            // if neighbor is at coarser level, restrict x3e
            } else {
              if (one_d) {
                rflx = flx.x3e(m,0,0,fi);
              } else if (two_d) {
                rflx = flx.x3e(m,0,fj,fi);
              } else {
                rflx = 0.5*(flx.x3e(m,fk,fj,fi) + flx.x3e(m,fk+1,fj,fi));
              }
            }
            if (nghbr.d_view(m,n).rank == my_rank) {
              rbuf[dn].flux(dm, ndat*v + (j-jl + nj*(k-kl))) = rflx;
            } else {
              sbuf[n].flux(m, ndat*v + (j-jl + nj*(k-kl))) = rflx;
            }
          }
        });

      // x2faces (only load x1e and x3e)
      } else if (n<16) {
        // j-index is fixed for flux correction on x2faces
        const int j = jl;
        const int fj = 2*jl - cjs;
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nki), [&](const int idx) {
          int k = idx / ni;
          int i = (idx - k * ni) + il;
          k += kl;
          int fk = 2*k - cks;
          int fi = 2*i - cis;
          if (v==0) {
            Real rflx;
            // if neighbor is at same level, load x1e directly
            if (nghbr.d_view(m,n).lev == mblev.d_view(m)) {
              rflx = flx.x1e(m,k,j,i);
            // if neighbor is at coarser level, restrict x1e
            } else {
              if (two_d) {
                rflx = 0.5*(flx.x1e(m,0,fj,fi) + flx.x1e(m,0,fj,fi+1));
              } else {
                rflx = 0.5*(flx.x1e(m,fk,fj,fi) + flx.x1e(m,fk,fj,fi+1));
              }
            }
            if (nghbr.d_view(m,n).rank == my_rank) {
              rbuf[dn].flux(dm, ndat*v + i-il + ni*(k-kl)) = rflx;
            } else {
              sbuf[n].flux(m, ndat*v + i-il + ni*(k-kl)) = rflx;
            }
          } else if (v==2) {
            Real rflx;
            // if neighbor is at same level, load x3e directly
            if (nghbr.d_view(m,n).lev == mblev.d_view(m)) {
              rflx = flx.x3e(m,k,j,i);
            // if neighbor is at coarser level, restrict x3e
            } else {
              if (two_d) {
                rflx = flx.x3e(m,0,fj,fi);
              } else {
                rflx = 0.5*(flx.x3e(m,fk,fj,fi) + flx.x3e(m,fk+1,fj,fi));
              }
            }
            if (nghbr.d_view(m,n).rank == my_rank) {
              rbuf[dn].flux(dm, ndat*v + i-il + ni*(k-kl)) = rflx;
            } else {
              sbuf[n].flux(m, ndat*v + i-il + ni*(k-kl)) = rflx;
            }
          }
        });

      // x1x2 edges (only load x3e)
      } else if (n<24) {
        // i/j-index is fixed for flux correction on x1x2 edges
        const int i = il;
        const int j = jl;
        const int fi = 2*il - cis;
        const int fj = 2*jl - cjs;
        if (v==2) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember,nk),[&](const int idx) {
            int k = idx + kl;
            int fk = 2*k - cks;
            Real rflx;
            // if neighbor is at same level, load x3e directly
            if (nghbr.d_view(m,n).lev == mblev.d_view(m)) {
              rflx = flx.x3e(m,k,j,i);
            // if neighbor is at coarser level, restrict x3e
            } else {
              if (two_d) {
                rflx = flx.x3e(m,0,fj,fi);
              } else {
                rflx = 0.5*(flx.x3e(m,fk,fj,fi) + flx.x3e(m,fk+1,fj,fi));
              }
            }
            if (nghbr.d_view(m,n).rank == my_rank) {
              rbuf[dn].flux(dm, ndat*v + (k-kl)) = rflx;
            } else {
              sbuf[n].flux(m, ndat*v + (k-kl)) = rflx;
            }
          });
        }

      // x3faces (only load x1e and x2e)
      } else if (n<32) {
        // k-index is fixed for flux correction on x3faces
        const int k = kl;
        const int fk = 2*kl - cks;
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nji), [&](const int idx) {
          int j = idx / ni;
          int i = (idx - j * ni) + il;
          j += jl;
          int fi = 2*i - cis;
          int fj = 2*j - cjs;
          if (v==0) {
            Real rflx;
            // if neighbor is at same level, load x1e directly
            if (nghbr.d_view(m,n).lev == mblev.d_view(m)) {
              rflx = flx.x1e(m,k,j,i);
            // if neighbor is at coarser level, restrict x1e
            } else {
              rflx = 0.5*(flx.x1e(m,fk,fj,fi) + flx.x1e(m,fk,fj,fi+1));
            }
            if (nghbr.d_view(m,n).rank == my_rank) {
              rbuf[dn].flux(dm, ndat*v + i-il + ni*(j-jl)) = rflx;
            } else {
              sbuf[n].flux(m, ndat*v + i-il + ni*(j-jl)) = rflx;
            }
          } else if (v==1) {
            Real rflx;
            // if neighbor is at same level, load x2e directly
            if (nghbr.d_view(m,n).lev == mblev.d_view(m)) {
              rflx = flx.x2e(m,k,j,i);
            // if neighbor is at coarser level, restrict x2e
            } else {
              rflx = 0.5*(flx.x2e(m,fk,fj,fi) + flx.x2e(m,fk,fj+1,fi));
            }
            if (nghbr.d_view(m,n).rank == my_rank) {
              rbuf[dn].flux(dm, ndat*v + i-il + ni*(j-jl)) = rflx;
            } else {
              sbuf[n].flux(m, ndat*v + i-il + ni*(j-jl)) = rflx;
            }
          }
        });

      // x3x1 edges (only load x2e)
      } else if (n<40) {
        const int i = il;
        const int k = kl;
        const int fi = 2*il - cis;
        const int fk = 2*kl - cks;
        if (v==1) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember,nj),[&](const int idx) {
            int j = idx + jl;
            int fj = 2*j - cjs;
            Real rflx;
            // if neighbor is at same level, load x2e directly
            if (nghbr.d_view(m,n).lev == mblev.d_view(m)) {
              rflx = flx.x2e(m,k,j,i);
            // if neighbor is at coarser level, restrict x2e
            } else {
              rflx = 0.5*(flx.x2e(m,fk,fj,fi) + flx.x2e(m,fk,fj+1,fi));
            }
            if (nghbr.d_view(m,n).rank == my_rank) {
              rbuf[dn].flux(dm, ndat*v + (j-jl)) = rflx;
            } else {
              sbuf[n].flux(m, ndat*v + (j-jl)) = rflx;
            }
          });
        }

      // x2x3 edges (only load x1e)
      } else if (n<48) {
        const int j = jl;
        const int k = kl;
        const int fj = 2*jl - cjs;
        const int fk = 2*kl - cks;
        if (v==0) {
          Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember,ni),[&](const int idx) {
            int i = idx + il;
            int fi = 2*i - cis;
            Real rflx;
            // if neighbor is at same level, load x1e directly
            if (nghbr.d_view(m,n).lev == mblev.d_view(m)) {
              rflx = flx.x1e(m,k,j,i);
            // if neighbor is at coarser level, restrict x1e
            } else {
              rflx = 0.5*(flx.x1e(m,fk,fj,fi) + flx.x1e(m,fk,fj,fi+1));
            }
            if (nghbr.d_view(m,n).rank == my_rank) {
              rbuf[dn].flux(dm, ndat*v + i-il) = rflx;
            } else {
              sbuf[n].flux(m, ndat*v + i-il) = rflx;
            }
          });
        }
      }
    }  // end if-neighbor-exists block
    tmember.team_barrier();
  });  // end par_for_outer

#if MPI_PARALLEL_ENABLED
  // Send boundary buffer to neighboring MeshBlocks using MPI
  // Sends only occur to neighbors on FACES and EDGES at COARSER or SAME level
  Kokkos::fence();
  bool no_errors=true;
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      const bool eligible_neighbor = same_level_only
          ? (mblev.h_view(m) == active_level &&
             nghbr.h_view(m,n).lev == mblev.h_view(m))
          : (nghbr.h_view(m,n).lev <= mblev.h_view(m));
      if ( (nghbr.h_view(m,n).gid >=0) && eligible_neighbor &&
           (n<48) ) {
        // index and rank of destination Neighbor
        int dn = nghbr.h_view(m,n).dest;
        int drank = nghbr.h_view(m,n).rank;

        if (drank != my_rank) {
          // create tag using local ID and buffer index of *receiving* MeshBlock
          int lid = nghbr.h_view(m,n).gid - pmy_pack->pmesh->gids_eachrank[drank];
          int tag = CreateBvals_MPI_Tag(lid, dn);

          // get ptr to send buffer for fluxes
          int data_size = 3;
          if ( nghbr.h_view(m,n).lev < pmy_pack->pmb->mb_lev.h_view(m) ) {
            data_size *= sendbuf[n].iflxc_ndat;
          } else if ( nghbr.h_view(m,n).lev == pmy_pack->pmb->mb_lev.h_view(m) ) {
            data_size *= sendbuf[n].iflxs_ndat;
          }
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
//! \fn void RecvAndUnpackFluxFC()
//! \brief Unpack boundary buffers for flux correction of FC variables.  This requires
//! averaging together fluxes from MeshBlocks at the same level, or replacing the fluxes
//! with the average from MeshBlocks at finer levels.

TaskStatus MeshBoundaryValuesFC::RecvAndUnpackFluxFC(DvceEdgeFld4D<Real> &flx) {
  return RecvAndUnpackFluxFC(flx, false);
}

//----------------------------------------------------------------------------------------
//! \brief Receive and average edge fields, optionally restricting work to active peers.

TaskStatus MeshBoundaryValuesFC::RecvAndUnpackFluxFC(DvceEdgeFld4D<Real> &flx,
                                                      bool same_level_only) {
  // create local references for variables in kernel
  int nmb = pmy_pack->nmb_thispack;
#if MPI_PARALLEL_ENABLED
  int nnghbr = pmy_pack->pmb->nnghbr;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &rbuf = recvbuf;
  auto &mblev = pmy_pack->pmb->mb_lev;
  //----- STEP 1: check that recv boundary buffer communications have all completed
  // receives only occur for neighbors on faces and edges at FINER or SAME level

  bool bflag = false;
  bool no_errors=true;
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      const bool target_is_active = (mblev.h_view(m) == pmy_pack->active_level);
      const bool eligible_neighbor = same_level_only
          ? (target_is_active && nghbr.h_view(m,n).lev == mblev.h_view(m))
          : (nghbr.h_view(m,n).lev >= mblev.h_view(m));
      if ( (nghbr.h_view(m,n).gid >=0) && eligible_neighbor &&
           (n<48) ) {
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

  //----- STEP 2: buffers have all completed, so unpack and perform appropriate averaging

  // 2D array to store number of fluxes summed into corner buffers
  DvceArray2D<int> nflx("nflx",nmb,48);
  par_for("init_nflx", DevExeSpace(), 0, (nmb-1), 0, 47,
  KOKKOS_LAMBDA(const int m, const int n) {
    nflx(m,n) = 1;
  });

  // Unpack and sum fluxes from the same level
  SumBoundaryFluxes(flx, true, nflx);

  // Zero EMFs at boundary that overlap with finer MeshBlocks (only use fine fluxes there)
  // Then unpack and sum fluxes from finer levels
  if (pmy_pack->pmesh->multilevel && !same_level_only) {
    ZeroFluxesAtBoundaryWithFiner(flx, nflx);
    SumBoundaryFluxes(flx, false, nflx);
  }

  // perform appropriate averaging depending on how many fluxes contributed to sums
  AverageBoundaryFluxes(flx, nflx, same_level_only);

  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn TaskStatus MeshBoundaryValuesFC::LoadFluxRegistersFC()
//! \brief Materialize the time-integrated EMF mismatch in an edge-field scratch array.
//!
//! The persistent registers are copied into the transient receive buffers so the legacy
//! fine-neighbor summation and T-junction averaging machinery can be reused exactly.
//! flx_scratch is cleared first and contains only the fine-minus-coarse integrated EMF on
//! coarse_level interfaces on return.  This routine does not apply a curl to B; its caller
//! owns the CT stencil.  Registers for this level are reset after the scratch is complete.

TaskStatus MeshBoundaryValuesFC::LoadFluxRegistersFC(
    DvceEdgeFld4D<Real> &flx_scratch, int coarse_level) {
  const int nmb = pmy_pack->nmb_thispack;
  const int nnghbr = pmy_pack->pmb->nnghbr;
  ValidateFluxRegistersFC(recvbuf, nnghbr, nmb);

  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &mblev = pmy_pack->pmb->mb_lev;
  auto &rbuf = recvbuf;

  // The edge fields and transient receive buffers are disposable scratch at a level
  // synchronization point.  Clearing every row prevents a non-active level from seeing
  // stale stage communication when the legacy summation kernels traverse the whole pack.
  Kokkos::deep_copy(flx_scratch.x1e, 0.0);
  Kokkos::deep_copy(flx_scratch.x2e, 0.0);
  Kokkos::deep_copy(flx_scratch.x3e, 0.0);
  for (int n=0; n<nnghbr && IsFCFluxBuffer(n); ++n) {
    Kokkos::deep_copy(rbuf[n].flux, 0.0);
  }

  Kokkos::TeamPolicy<> policy(DevExeSpace(), nmb*nnghbr, Kokkos::AUTO);
  Kokkos::parallel_for("load_fc_emf_register", policy,
  KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = tmember.league_rank()/nnghbr;
    const int n = tmember.league_rank() - m*nnghbr;
    if (n >= 48 ||
        mblev.d_view(m) != coarse_level ||
        nghbr.d_view(m,n).gid < 0 ||
        nghbr.d_view(m,n).lev != coarse_level + 1) {
      return;
    }
    const int ndata = 3*rbuf[n].iflxc_ndat;
    Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, ndata),
    [&](const int q) {
      rbuf[n].flux(m,q) = rbuf[n].flux_reg(m,q);
    });
    tmember.team_barrier();
  });
  Kokkos::fence();

  // Start with the same self-contribution count as the legacy path, then zero it exactly
  // where a finer neighbor replaces the coarse value.  Since each register key already
  // contains its own coarse subtraction, averaging overlapping fine keys yields
  // average(fine)-coarse at face seams, edges, and T-junctions.
  DvceArray2D<int> nflx("fc_register_nflx", nmb, 48);
  par_for("init_fc_register_nflx", DevExeSpace(), 0, nmb-1, 0, 47,
  KOKKOS_LAMBDA(const int m, const int n) {
    nflx(m,n) = 1;
  });
  ZeroFluxesAtBoundaryWithFiner(flx_scratch, nflx);
  SumBoundaryFluxes(flx_scratch, false, nflx);
  AverageBoundaryFluxes(flx_scratch, nflx, false);

  // Prevent accidental double loading, and make the synchronization API host-complete.
  (void) ResetFluxRegistersFC(coarse_level);
  Kokkos::fence();
  return TaskStatus::complete;
}

//----------------------------------------------------------------------------------------
//! \fn  void BoundaryValuesFC::SumBoundaryFluxes
//! \brief Sums boundary buffer fluxes from neighboring MeshBlocks at the same level into
//! flux (e.g. EMF) array if input argument 'same_level=true', or sums boundary buffer
//! fluxes from neighboring MeshBlocks at a finer level into flux array otherwise.

void MeshBoundaryValuesFC::SumBoundaryFluxes(DvceEdgeFld4D<Real> &flx,
                                          const bool same_level, DvceArray2D<int> &nflx) {
  // create local references for variables in kernel
  int nmb = pmy_pack->nmb_thispack;
  int nnghbr = pmy_pack->pmb->nnghbr;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &rbuf = recvbuf;
  auto &mblev = pmy_pack->pmb->mb_lev;
  auto &mbbcs = pmy_pack->pmb->mb_bcs;

  // Sum receive buffers into EMFs stored on MeshBlocks
  // Outer loop over (# of MeshBlocks)*(3 field components)
  Kokkos::TeamPolicy<> policy(DevExeSpace(), (3*nmb), Kokkos::AUTO);
  Kokkos::parallel_for("RecvBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = tmember.league_rank()/3;
    const int v = tmember.league_rank()%3;

    // scalar loop over neighbors (except corners) to prevent race condition in sums
    for (int n=0; n<nnghbr; ++n) {
      // only unpack buffers when neighbor exists AND
      // (neighbor at same level when same_level=true on input) OR
      // (neighbor at finer level when same_level=false on input)
      if (    (nghbr.d_view(m,n).gid >= 0) &&
           (( (same_level) && (nghbr.d_view(m,n).lev == mblev.d_view(m))) ||
            (!(same_level) && (nghbr.d_view(m,n).lev >  mblev.d_view(m)))) ) {
        bool inner_bufs=false, outer_bufs=false;
        if ((n==0)||(n==16)||(n==20)||(n==32)||(n==33)||(n==36)||(n==37)) {
          inner_bufs=true;
        }
        if ((n==4)||(n==18)||(n==22)||(n==34)||(n==35)||(n==38)||(n==39)) {
          outer_bufs=true;
        }

      // only unpack buffers when
      //   (both innerx1 AND outerx1 BCs are NOT shear_periodic) OR
      //   (only innerx1 BC is shear_periodic AND n!=0,16,20,32,33,36,37) OR
      //   (only outerx1 BC is shear_periodic AND n!=4,18,22,34,35,38,39) OR
      //   (both innerx1 AND outerx1 BC are shear_periodic AND n!=...)
      if ( ((mbbcs.d_view(m,0)!=BoundaryFlag::shear_periodic) &&
            (mbbcs.d_view(m,1)!=BoundaryFlag::shear_periodic)) ||
           ((mbbcs.d_view(m,0)==BoundaryFlag::shear_periodic) && !(inner_bufs) &&
            (mbbcs.d_view(m,1)!=BoundaryFlag::shear_periodic)) ||
           ((mbbcs.d_view(m,0)!=BoundaryFlag::shear_periodic) &&
            (mbbcs.d_view(m,1)==BoundaryFlag::shear_periodic) && !(outer_bufs)) ||
           ((mbbcs.d_view(m,0)==BoundaryFlag::shear_periodic) && !(inner_bufs) &&
            (mbbcs.d_view(m,1)==BoundaryFlag::shear_periodic) && !(outer_bufs)) ) {
        int il, iu, jl, ju, kl, ku, ndat;
        // if neighbor is at same level, use same indices to unpack buffer
        if (same_level) {
          il = rbuf[n].iflux_same[v].bis;
          iu = rbuf[n].iflux_same[v].bie;
          jl = rbuf[n].iflux_same[v].bjs;
          ju = rbuf[n].iflux_same[v].bje;
          kl = rbuf[n].iflux_same[v].bks;
          ku = rbuf[n].iflux_same[v].bke;
          ndat = rbuf[n].iflxs_ndat;
        // else neighbor is at finer level, use flux_coar indices to unpack buffer
        } else {
          il = rbuf[n].iflux_coar[v].bis;
          iu = rbuf[n].iflux_coar[v].bie;
          jl = rbuf[n].iflux_coar[v].bjs;
          ju = rbuf[n].iflux_coar[v].bje;
          kl = rbuf[n].iflux_coar[v].bks;
          ku = rbuf[n].iflux_coar[v].bke;
          ndat = rbuf[n].iflxc_ndat;
        }
        const int ni = iu - il + 1;
        const int nj = ju - jl + 1;
        const int nk = ku - kl + 1;
        const int nji  = nj*ni;
        const int nkj  = nk*nj;
        const int nki  = nk*ni;

        // x1faces
        if (n<8) {
          // always use v=0 thread index in sums to avoid race condition
          if (v==0) {
            if (n==0) {
              Kokkos::single(Kokkos::PerTeam(tmember), [&] () {
                nflx(m,16) += 1; nflx(m,20) += 1; nflx(m,32) += 1; nflx(m,36) += 1;
              });
            }
            if (n==4) {
              Kokkos::single(Kokkos::PerTeam(tmember), [&] () {
                nflx(m,18) += 1; nflx(m,22) += 1; nflx(m,34) += 1; nflx(m,38) += 1;
              });
            }
          } else {
            Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember,nkj),
            [&](const int idx){
              int k = idx / nj;
              int j = (idx - k * nj) + jl;
              k += kl;
              if (v==1) {
                flx.x2e(m,k,j,il) += rbuf[n].flux(m,ndat*v + (j-jl + nj*(k-kl)));
              } else if (v==2) {
                flx.x3e(m,k,j,il) += rbuf[n].flux(m,ndat*v + (j-jl + nj*(k-kl)));
              }
            });
          }

        // x2faces
        } else if (n<16) {
          // always use v=0 thread index in sums to avoid race condition
          if (v==0) {
            if (n==8) {
              Kokkos::single(Kokkos::PerTeam(tmember), [&] () {
                nflx(m,16) += 1; nflx(m,18) += 1; nflx(m,40) += 1; nflx(m,44) += 1;
              });
            }
            if (n==12) {
              Kokkos::single(Kokkos::PerTeam(tmember), [&] () {
                nflx(m,20) += 1; nflx(m,22) += 1; nflx(m,42) += 1; nflx(m,46) += 1;
              });
            }
          }
          Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember,nki),
          [&](const int idx){
            int k = idx/ni;
            int i = (idx - k * ni) + il;
            k += kl;
            if (v==0) {
              flx.x1e(m,k,jl,i) += rbuf[n].flux(m,ndat*v + i-il + ni*(k-kl));
            } else if (v==2) {
              flx.x3e(m,k,jl,i) += rbuf[n].flux(m,ndat*v + i-il + ni*(k-kl));
            }
          });

        // x1x2 edges
        } else if (n<24) {
          // always use v=0 thread index in sums to avoid race condition
          if (v==0) {
            Kokkos::single(Kokkos::PerTeam(tmember), [&] () {
              nflx(m,n) += 1;
            });
          } else if (v==2) {
            Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember,nk),[&](const int idx){
              int k = idx + kl;
              flx.x3e(m,k,jl,il) += rbuf[n].flux(m,ndat*v + (k-kl));
            });
          }

        // x3faces
        } else if (n<32)  {
          // always use v=0 thread index in sums to avoid race condition
          if (v==0) {
            if (n==24) {
              Kokkos::single(Kokkos::PerTeam(tmember), [&] () {
                nflx(m,32) += 1; nflx(m,34) += 1; nflx(m,40) += 1; nflx(m,42) += 1;
              });
            }
            if (n==28) {
              Kokkos::single(Kokkos::PerTeam(tmember), [&] () {
                nflx(m,36) += 1; nflx(m,38) += 1; nflx(m,44) += 1; nflx(m,46) += 1;
              });
            }
          }
          Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember,nji),
          [&](const int idx){
            int j = idx / ni;
            int i = (idx - j * ni) + il;
            j += jl;
            if (v==0) {
              flx.x1e(m,kl,j,i) += rbuf[n].flux(m,ndat*v + i-il + ni*(j-jl));
            } else if (v==1) {
              flx.x2e(m,kl,j,i) += rbuf[n].flux(m,ndat*v + i-il + ni*(j-jl));
            }
          });

        // x3x1 edges
        } else if (n<40) {
          // always use v=0 thread index in sums to avoid race condition
          if (v==0) {
            Kokkos::single(Kokkos::PerTeam(tmember), [&] () {
              nflx(m,n) += 1;
            });
          } else if (v==1) {
            Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember,nj),[&](const int idx){
              int j = idx + jl;
              flx.x2e(m,kl,j,il) += rbuf[n].flux(m,ndat*v + (j-jl));
            });
          }

        // x2x3 edges
        } else if (n<48) {
          if (v==0) {
            Kokkos::single(Kokkos::PerTeam(tmember), [&] () {
              nflx(m,n) += 1;
            });
            Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember,ni),[&](const int idx){
              int i = idx + il;
              flx.x1e(m,kl,jl,i) += rbuf[n].flux(m,ndat*v + i-il);
            });
          }
        }
      }  // end if-neighbor-exists block
      }
      tmember.team_barrier();
    }    // end for loop over n
  });    // end par_for_outer

  return;
}

//----------------------------------------------------------------------------------------
//! \fn  void MeshBoundaryValuesFC::ZeroFluxesAtBoundaryWithFiner
//! \brief Zeroes out fluxes of face-centered variables (e.g. EMFs) at boundaries with
//! MeshBlocks at a finer level, so that boundary buffer fluxes from finer level can be
//! summed (averaged) in place.

void MeshBoundaryValuesFC::ZeroFluxesAtBoundaryWithFiner(DvceEdgeFld4D<Real> &flx,
                                                         DvceArray2D<int> &nflx) {
  // create local references for variables in kernel
  int nmb = pmy_pack->nmb_thispack;
  int nnghbr = pmy_pack->pmb->nnghbr;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &rbuf = recvbuf;
  auto &mblev = pmy_pack->pmb->mb_lev;

  // Outer loop over (# of MeshBlocks)*(# of neighbors)*(3 field components)
  Kokkos::TeamPolicy<> policy(DevExeSpace(), (3*nmb*nnghbr), Kokkos::AUTO);
  Kokkos::parallel_for("RecvBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = (tmember.league_rank())/(3*nnghbr);
    const int n = (tmember.league_rank() - m*(3*nnghbr))/3;
    const int v = (tmember.league_rank() - m*(3*nnghbr) - 3*n);

    // only zero EMFs when neighbor exists and is at finer level
    if ((nghbr.d_view(m,n).gid >= 0) && (nghbr.d_view(m,n).lev > mblev.d_view(m))) {
      int il, iu, jl, ju, kl, ku;
      il = rbuf[n].iflux_coar[v].bis;
      iu = rbuf[n].iflux_coar[v].bie;
      jl = rbuf[n].iflux_coar[v].bjs;
      ju = rbuf[n].iflux_coar[v].bje;
      kl = rbuf[n].iflux_coar[v].bks;
      ku = rbuf[n].iflux_coar[v].bke;
      const int ni = iu - il + 1;
      const int nj = ju - jl + 1;
      const int nk = ku - kl + 1;
      const int nji  = nj*ni;
      const int nkj  = nk*nj;
      const int nki  = nk*ni;

      // x1faces
      if (n<8) {
        // use idle thread index to zero number of fluxes at corners of x1faces
        if (v==0) {
          if (n==0) {
            nflx(m,16) = 0; nflx(m,20) = 0; nflx(m,32) = 0; nflx(m,36) = 0;
          }
          if (n==4) {
            nflx(m,18) = 0; nflx(m,22) = 0; nflx(m,34) = 0; nflx(m,38) = 0;
          }
        // else zero fluxes
        } else {
          Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember,nkj),[&](const int idx){
            int k = idx / nj;
            int j = (idx - k * nj) + jl;
            k += kl;
            if (v==1) {
              flx.x2e(m,k,j,il) = 0.0;
            } else if (v==2) {
              flx.x3e(m,k,j,il) = 0.0;
            }
          });
        }

      // x2faces
      } else if (n<16) {
        // use idle thread index to zero number of fluxes at corners of x2faces
        if (v==1) {
          if (n==8) {
            nflx(m,16) = 0; nflx(m,18) = 0; nflx(m,40) = 0; nflx(m,44) = 0;
          }
          if (n==12) {
            nflx(m,20) = 0; nflx(m,22) = 0; nflx(m,42) = 0; nflx(m,46) = 0;
          }
        // else zero fluxes
        } else {
          Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember,nki),[&](const int idx){
            int k = idx/ni;
            int i = (idx - k * ni) + il;
            k += kl;
            if (v==0) {
              flx.x1e(m,k,jl,i) = 0.0;
            } else if (v==2) {
              flx.x3e(m,k,jl,i) = 0.0;
            }
          });
        }

      // x1x2 edges
      } else if (n<24) {
        if (v==2) {
          nflx(m,n) = 0;
          Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember,nk),[&](const int idx) {
            int k = idx + kl;
            flx.x3e(m,k,jl,il) = 0.0;
          });
        }

      // x3faces
      } else if (n<32)  {
        // use idle thread index to zero number of fluxes at corners of x2faces
        if (v==2) {
          if (n==24) {
            nflx(m,32) = 0; nflx(m,34) = 0; nflx(m,40) = 0; nflx(m,42) = 0;
          }
          if (n==28) {
            nflx(m,36) = 0; nflx(m,38) = 0; nflx(m,44) = 0; nflx(m,46) = 0;
          }
        // else zero fluxes
        } else {
          Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember,nji),[&](const int idx){
            int j = idx / ni;
            int i = (idx - j * ni) + il;
            j += jl;
            if (v==0) {
              flx.x1e(m,kl,j,i) = 0.0;
            } else if (v==1) {
              flx.x2e(m,kl,j,i) = 0.0;
            }
          });
        }

      // x3x1 edges
      } else if (n<40) {
        if (v==1) {
          nflx(m,n) = 0;
          Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember,nj),[&](const int idx){
            int j = idx + jl;
              flx.x2e(m,kl,j,il) = 0.0;
          });
        }

      // x2x3 edges
      } else if (n<48) {
        if (v==0) {
          nflx(m,n) = 0;
          Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember,ni),[&](const int idx){
            int i = idx + il;
              flx.x1e(m,kl,jl,i) = 0.0;
          });
        }
      }
    }  // end if-neighbor-exists block
    tmember.team_barrier();
  });  // end par_for_outer

  return;
}

//----------------------------------------------------------------------------------------
//! \fn  void MeshBoundaryValuesFC::AverageBoundaryFluxes
//! \brief Applies appropriate average to summed boundary fluxes, depending on number of
//! elements being averaged together.

void MeshBoundaryValuesFC::AverageBoundaryFluxes(DvceEdgeFld4D<Real> &flx,
                                                 DvceArray2D<int> &nflx,
                                                 const bool same_level_only) {
  // create local references for variables in kernel
  int nmb = pmy_pack->nmb_thispack;
  int nnghbr = pmy_pack->pmb->nnghbr;
  auto &nghbr = pmy_pack->pmb->nghbr;
  auto &rbuf = recvbuf;
  auto &mblev = pmy_pack->pmb->mb_lev;
  auto &mbbcs = pmy_pack->pmb->mb_bcs;
  bool &multi_d = pmy_pack->pmesh->multi_d;
  bool &three_d = pmy_pack->pmesh->three_d;

  // Outer loop over (# of MeshBlocks)*(# of neighbors)*(3 field components)
  Kokkos::TeamPolicy<> policy(DevExeSpace(), (3*nmb*nnghbr), Kokkos::AUTO);
  Kokkos::parallel_for("RecvBuff", policy, KOKKOS_LAMBDA(TeamMember_t tmember) {
    const int m = (tmember.league_rank())/(3*nnghbr);
    const int n = (tmember.league_rank() - m*(3*nnghbr))/3;
    const int v = (tmember.league_rank() - m*(3*nnghbr) - 3*n);
    // only average when
    //   (both innerx1 AND outerx1 BCs are NOT shear_periodic) OR
    //   (only innerx1 BC is shear_periodic AND n!=0) OR
    //   (only outerx1 BC is shear_periodic AND n!=4) OR
    //   (both innerx1 AND outerx1 BC are shear_periodic AND n!=0,4)
    if ( ((mbbcs.d_view(m,0)!=BoundaryFlag::shear_periodic) &&
          (mbbcs.d_view(m,1)!=BoundaryFlag::shear_periodic)) ||
         ((mbbcs.d_view(m,0)==BoundaryFlag::shear_periodic) && (n!=0) &&
          (mbbcs.d_view(m,1)!=BoundaryFlag::shear_periodic)) ||
         ((mbbcs.d_view(m,0)!=BoundaryFlag::shear_periodic) &&
          (mbbcs.d_view(m,1)==BoundaryFlag::shear_periodic) && (n!=4)) ||
         ((mbbcs.d_view(m,0)==BoundaryFlag::shear_periodic) && (n!=0) &&
          (mbbcs.d_view(m,1)==BoundaryFlag::shear_periodic) && (n!=4)) ) {
      int il, iu, jl, ju, kl, ku;
      il = rbuf[n].iflux_same[v].bis;
      iu = rbuf[n].iflux_same[v].bie;
      jl = rbuf[n].iflux_same[v].bjs;
      ju = rbuf[n].iflux_same[v].bje;
      kl = rbuf[n].iflux_same[v].bks;
      ku = rbuf[n].iflux_same[v].bke;

      // x1faces
      if (n==0 || n==4) {
        if (v==1) {
          int nj = ju - jl + 1;
          // same level; divide EMFs on face by 2, excluding edges
          if (nghbr.d_view(m,n).lev == mblev.d_view(m)) {
            if (three_d) {
              kl += 1; ku -= 1;
            }
            int nk = ku - kl + 1;
            Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nk*nj),
            [&](const int idx) {
              int k = idx / nj;
              int j = (idx - k * nj) + jl;
              k += kl;
              flx.x2e(m,k,j,il) *= 0.5;
            });
          // finer level; divide EMFs that overlap at edges of fine faces by 2
          } else if (!same_level_only &&
                     nghbr.d_view(m,n).lev >= mblev.d_view(m)) {
            if (three_d) {
              int k = kl + (ku - kl + 1)/2;
              Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nj),
              [&](const int idx) {
                int j = idx + jl;
                flx.x2e(m,k,j,il) *= 0.5;
              });
              tmember.team_barrier();
            }
          }
        } else if (v==2) {
          int nk = ku - kl + 1;
          // same level; divide EMFs on face by 2, excluding edges
          if (nghbr.d_view(m,n).lev == mblev.d_view(m)) {
            if (multi_d) {
              jl += 1; ju -= 1;
            }
            int nj = ju - jl + 1;
            Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nk*nj),
            [&](const int idx) {
              int k = idx / nj;
              int j = (idx - k * nj) + jl;
              k += kl;
              flx.x3e(m,k,j,il) *= 0.5;
            });
          // finer level; divide EMFs that overlap at edges of fine faces by 2
          } else if (!same_level_only &&
                     nghbr.d_view(m,n).lev >= mblev.d_view(m)) {
            if (multi_d) {
              int j = jl + (ju - jl + 1)/2;
              Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nk),
              [&](const int idx) {
              int k = idx + kl;
                  flx.x3e(m,k,j,il) *= 0.5;
              });
            }
          }
        }

      // x2faces
      } else if (multi_d && (n==8 || n==12)) {
        if (v==0) {
          int ni = iu - il + 1;
          // same level; divide EMFs on face by 2, excluding edges
          if (nghbr.d_view(m,n).lev == mblev.d_view(m)) {
            if (three_d) {
              kl += 1; ku -= 1;
            }
            int nk = ku - kl + 1;
            Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nk*ni),
            [&](const int idx) {
              int k = idx/ni;
              int i = (idx - k * ni) + il;
              k += kl;
              flx.x1e(m,k,jl,i) *= 0.5;
            });
          // finer level; divide EMFs that overlap at edges of fine faces by 2
          } else if (!same_level_only &&
                     nghbr.d_view(m,n).lev >= mblev.d_view(m)) {
            if (three_d) {
              int k = kl + (ku - kl + 1)/2;
              Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, ni),
              [&](const int idx) {
                int i = idx + il;
                flx.x1e(m,k,jl,i) *= 0.5;
              });
              tmember.team_barrier();
            }
          }
        } else if (v==2) {
          int nk = ku - kl + 1;
          // same level; divide EMFs on face by 2, excluding edges
          if (nghbr.d_view(m,n).lev == mblev.d_view(m)) {
            il += 1; iu -= 1;
            int ni = iu - il + 1;
            Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nk*ni),
            [&](const int idx) {
              int k = idx/ni;
                int i = (idx - k * ni) + il;
              k += kl;
              flx.x3e(m,k,jl,i) *= 0.5;
            });
            tmember.team_barrier();
          // finer level; divide EMFs that overlap at edges of fine faces by 2
          } else if (!same_level_only &&
                     nghbr.d_view(m,n).lev >= mblev.d_view(m)) {
            int i = il + (iu - il + 1)/2;
            Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nk),
            [&](const int idx) {
              int k = idx + kl;
              flx.x3e(m,k,jl,i) *= 0.5;
            });
            tmember.team_barrier();
          }
        }

      // x1x2 edges
      } else if (multi_d && (n==16 || n==18 || n==20 || n==22)) {
        if (v==2) {
          int nk = ku - kl + 1;
          Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember,nk),[&](const int idx) {
            int k = idx + kl;
            flx.x3e(m,k,jl,il) /= static_cast<Real>(nflx(m,n));
          });
        }

      // x3faces
      } else if (three_d && (n==24 || n==28))  {
        if (v==0) {
          int ni = iu - il + 1;
          // same level; divide EMFs on face by 2, excluding edges
          if (nghbr.d_view(m,n).lev == mblev.d_view(m)) {
            jl += 1; ju -= 1;
            int nj = ju - jl + 1;
            Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nj*ni),
            [&](const int idx) {
              int j = idx / ni;
              int i = (idx - j * ni) + il;
              j += jl;
              flx.x1e(m,kl,j,i) *= 0.5;
            });
            tmember.team_barrier();
          // finer level; divide EMFs that overlap at edges of fine faces by 2
          } else if (!same_level_only &&
                     nghbr.d_view(m,n).lev >= mblev.d_view(m)) {
            int j = jl + (ju - jl + 1)/2;
            Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember,ni),[&](const int idx){
              int i = idx + il;
              flx.x1e(m,kl,j,i) *= 0.5;
            });
            tmember.team_barrier();
          }
        } else if (v==1) {
          int nj = ju - jl + 1;
          // same level; divide EMFs on face by 2, excluding edges
            if (nghbr.d_view(m,n).lev == mblev.d_view(m)) {
            il += 1; iu -= 1;
            int ni = iu - il + 1;
            Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember, nj*ni),
              [&](const int idx) {
                int j = idx / ni;
              int i = (idx - j * ni) + il;
              j += jl;
              flx.x2e(m,kl,j,i) *= 0.5;
            });
            tmember.team_barrier();
          // finer level; divide EMFs that overlap at edges of fine faces by 2
          } else if (!same_level_only &&
                     nghbr.d_view(m,n).lev >= mblev.d_view(m)) {
            int i = il + (iu - il + 1)/2;
            Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember,nj),[&](const int idx){
              int j = idx + jl;
              flx.x2e(m,kl,j,i) *= 0.5;
            });
            tmember.team_barrier();
          }
        }

      // x3x1 edges
      } else if (three_d && (n==32 || n==34 || n==36 || n==38)) {
        if (v==1) {
          int nj = ju - jl + 1;
          Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember,nj),[&](const int idx) {
            int j = idx + jl;
            flx.x2e(m,kl,j,il) /= static_cast<Real>(nflx(m,n));
          });
        }

    // x2x3 edges
    } else if (three_d && (n==40 || n==42 || n==44 || n==46)) {
      if (v==0) {
        int ni = iu - il + 1;
        Kokkos::parallel_for(Kokkos::TeamThreadRange<>(tmember,ni),[&](const int idx) {
          int i = idx + il;
          flx.x1e(m,kl,jl,i) /= static_cast<Real>(nflx(m,n));
        });
      }
    }
    }
    tmember.team_barrier();
  });    // end par_for_outer

  return;
}

//----------------------------------------------------------------------------------------
//! \fn  void MeshBoundaryValuesFC::InitRecvFlux
//! \brief Posts non-blocking receives (with MPI) for boundary communication of fluxes of
//! face-centered variables, which are communicated at FACES and EDGES of MeshBlocks at
//! the SAME or FINER levels.  This is different than for fluxes of cell-centered vars.

TaskStatus MeshBoundaryValuesFC::InitFluxRecv(const int nvars) {
#if MPI_PARALLEL_ENABLED
  int &nmb = pmy_pack->nmb_thispack;
  int &nnghbr = pmy_pack->pmb->nnghbr;
  auto &nghbr = pmy_pack->pmb->nghbr;

  // Initialize communications of fluxes
  bool no_errors=true;
  for (int m=0; m<nmb; ++m) {
    for (int n=0; n<nnghbr; ++n) {
      // only post receives for neighbors on FACES and EDGES at FINER and SAME levels
      // this is the only thing different from BoundaryValuesCC::InitRecvFlux()
      if ( (nghbr.h_view(m,n).gid >=0) &&
           (nghbr.h_view(m,n).lev >= pmy_pack->pmb->mb_lev.h_view(m)) &&
           (n<48) ) {
        // rank of destination buffer
        int drank = nghbr.h_view(m,n).rank;

        // post non-blocking receive if neighboring MeshBlock on a different rank
        if (drank != global_variable::my_rank) {
          // create tag using local ID and buffer index of *receiving* MeshBlock
          int tag = CreateBvals_MPI_Tag(m, n);

          // calculate amount of data to be passed, get pointer to variables
          int data_size = nvars;
          if ( nghbr.h_view(m,n).lev > pmy_pack->pmb->mb_lev.h_view(m) ) {
            data_size *= recvbuf[n].iflxc_ndat;
          } else if ( nghbr.h_view(m,n).lev == pmy_pack->pmb->mb_lev.h_view(m) ) {
            data_size *= recvbuf[n].iflxs_ndat;
          }
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
