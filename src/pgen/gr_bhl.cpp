//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file gr_bhl.cpp
//! \brief Uniform wind initial conditions for GR Bondi-Hoyle-Lyttleton accretion in
//! Cartesian Kerr-Schild coordinates.

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "coordinates/adm.hpp"
#include "coordinates/cartesian_ks.hpp"
#include "coordinates/cell_locations.hpp"
#include "dyn_grmhd/dyn_grmhd.hpp"
#include "eos/eos.hpp"
#include "geodesic-grid/spherical_grid.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"

namespace {

struct bhl_pgen {
  Real rho_inf;
  Real cs_inf;
  Real v1_inf;
  Real v2_inf;
  Real v3_inf;
  Real u1_inf;
  Real u2_inf;
  Real u3_inf;
  Real pgas_inf;
  Real eint_inf;
  Real beta_inf;
  Real b1_inf;
  Real b2_inf;
  Real b3_inf;
  Real r_acc;
  Real mdot_hl;
  Real drag_exclude_radius;
  bool grav_drag;
};

bhl_pgen bhl;

KOKKOS_INLINE_FUNCTION
void SetWindPrimitives(const bhl_pgen pgen, const bool is_mhd,
                       const Real gm1, DvceArray5D<Real> w0,
                       DvceArray5D<Real> bcc, const int m, const int k,
                       const int j, const int i) {
  w0(m,IDN,k,j,i) = pgen.rho_inf;
  w0(m,IEN,k,j,i) = (is_mhd ? pgen.eint_inf : pgen.eint_inf);
  w0(m,IVX,k,j,i) = pgen.u1_inf;
  w0(m,IVY,k,j,i) = pgen.u2_inf;
  w0(m,IVZ,k,j,i) = pgen.u3_inf;
  if (is_mhd) {
    bcc(m,IBX,k,j,i) = pgen.b1_inf;
    bcc(m,IBY,k,j,i) = pgen.b2_inf;
    bcc(m,IBZ,k,j,i) = pgen.b3_inf;
  }
}

} // namespace

void BHLWindInnerX1(Mesh *pm);
void BHLFluxes(HistoryData *pdata, Mesh *pm);

KOKKOS_INLINE_FUNCTION
Real KerrSchildRadius(const Real x, const Real y, const Real z, const Real spin) {
  const Real rad2 = SQR(x) + SQR(y) + SQR(z);
  const Real a2 = SQR(spin);
  return std::sqrt(0.5*(rad2 - a2 + std::sqrt(SQR(rad2 - a2) + 4.0*a2*SQR(z))));
}

//----------------------------------------------------------------------------------------
//! \fn void ProblemGenerator::UserProblem()
//! \brief Set uniform BHL wind primitives and a uniform optional magnetic field.

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
  MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;

  if (!pmbp->pcoord->is_general_relativistic) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "gr_bhl requires general_rel=true in the <coord> block."
              << std::endl;
    exit(EXIT_FAILURE);
  }
  if (pmbp->phydro == nullptr && pmbp->pmhd == nullptr) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "gr_bhl requires either a <hydro> or <mhd> block." << std::endl;
    exit(EXIT_FAILURE);
  }

  const bool is_mhd = (pmbp->pmhd != nullptr);
  user_bcs_func = BHLWindInnerX1;

  Real gamma = is_mhd ? pmbp->pmhd->peos->eos_data.gamma
                      : pmbp->phydro->peos->eos_data.gamma;
  Real gm1 = gamma - 1.0;

  bhl.rho_inf = pin->GetOrAddReal("problem", "rho_inf", 1.0);
  bhl.v1_inf = pin->GetOrAddReal("problem", "v1_inf", 0.1);
  bhl.v2_inf = pin->GetOrAddReal("problem", "v2_inf", 0.0);
  bhl.v3_inf = pin->GetOrAddReal("problem", "v3_inf", 0.0);

  Real v2 = SQR(bhl.v1_inf) + SQR(bhl.v2_inf) + SQR(bhl.v3_inf);
  if (v2 <= 0.0 || v2 >= 1.0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "BHL wind speed must be finite and subluminal."
              << std::endl;
    exit(EXIT_FAILURE);
  }

  // Hoyle-Lyttleton convention used by Kaaz+2023 and 2409.12359: R_a = 2GM/v_inf^2.
  bhl.r_acc = 2.0/v2;
  if (pin->DoesParameterExist("problem", "r_acc")) {
    Real input_r_acc = pin->GetReal("problem", "r_acc");
    Real rel_err = std::abs(input_r_acc - bhl.r_acc)/bhl.r_acc;
    if (rel_err > 1.0e-12) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl
                << "r_acc is derived from v_inf and must not be set independently. "
                << "For the supplied velocity, r_acc = " << bhl.r_acc << "."
                << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  Real mach = pin->GetOrAddReal("problem", "mach_inf", -1.0);
  if (mach > 0.0) {
    bhl.cs_inf = std::sqrt(v2)/mach;
  } else {
    bhl.cs_inf = pin->GetOrAddReal("problem", "cs_inf", 0.05);
  }
  if (SQR(bhl.cs_inf) >= gm1) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl
              << "cs_inf is too large for the ideal-gas pressure formula."
              << std::endl;
    exit(EXIT_FAILURE);
  }

  Real lor = 1.0/std::sqrt(1.0 - v2);
  bhl.u1_inf = lor*bhl.v1_inf;
  bhl.u2_inf = lor*bhl.v2_inf;
  bhl.u3_inf = lor*bhl.v3_inf;

  bhl.eint_inf = SQR(bhl.cs_inf)*bhl.rho_inf/(gamma*(gm1 - SQR(bhl.cs_inf)));
  bhl.pgas_inf = gm1*bhl.eint_inf;
  bhl.mdot_hl = M_PI*SQR(bhl.r_acc)*bhl.rho_inf*std::sqrt(v2);

  bhl.beta_inf = pin->GetOrAddReal("problem", "beta_inf", -1.0);
  bhl.b1_inf = bhl.b2_inf = bhl.b3_inf = 0.0;
  if (is_mhd && bhl.beta_inf > 0.0) {
    Real theta_b = pin->GetOrAddReal("problem", "theta_b", 0.0) * (M_PI/180.0);
    Real b0 = std::sqrt(2.0*bhl.pgas_inf/bhl.beta_inf*(1.0 - v2));
    bhl.b1_inf = pin->GetOrAddReal("problem", "b1_dir", 0.0) * b0;
    bhl.b2_inf = pin->GetOrAddReal("problem", "b2_dir", std::sin(theta_b)) * b0;
    bhl.b3_inf = pin->GetOrAddReal("problem", "b3_dir", std::cos(theta_b)) * b0;
  }

  if (pin->GetOrAddBoolean("problem", "user_hist", false)) {
    auto &grids = spherical_grids;
    Real r_hist = pin->GetOrAddReal("problem", "r_hist", 3.0);
    int nangle = pin->GetOrAddInteger("problem", "hist_nangle", 5);
    bhl.grav_drag = pin->GetOrAddBoolean("problem", "grav_drag", true);
    bhl.drag_exclude_radius = pin->GetOrAddReal("problem", "drag_exclude_radius", -1.0);
    grids.push_back(std::make_unique<SphericalGrid>(pmbp, nangle, r_hist));
    user_hist_func = BHLFluxes;
  }

  if (restart) return;

  auto &indcs = pmy_mesh_->mb_indcs;
  auto &size = pmbp->pmb->mb_size;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb = pmbp->nmb_thispack;

  DvceArray5D<Real> w0, u0, bcc;
  if (is_mhd) {
    w0 = pmbp->pmhd->w0;
    u0 = pmbp->pmhd->u0;
    bcc = pmbp->pmhd->bcc0;
  } else {
    w0 = pmbp->phydro->w0;
    u0 = pmbp->phydro->u0;
  }

  auto bhl_ = bhl;
  par_for("pgen_bhl_wind", DevExeSpace(), 0,nmb-1,ks,ke,js,je,is,ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    SetWindPrimitives(bhl_, is_mhd, gm1, w0, bcc, m, k, j, i);
  });

  if (is_mhd) {
    auto &b0 = pmbp->pmhd->b0;
    par_for("pgen_bhl_bfield", DevExeSpace(), 0,nmb-1,ks,ke,js,je,is,ie,
    KOKKOS_LAMBDA(int m, int k, int j, int i) {
      b0.x1f(m,k,j,i) = bhl_.b1_inf;
      b0.x2f(m,k,j,i) = bhl_.b2_inf;
      b0.x3f(m,k,j,i) = bhl_.b3_inf;
      if (i == ie) b0.x1f(m,k,j,i+1) = bhl_.b1_inf;
      if (j == je) b0.x2f(m,k,j+1,i) = bhl_.b2_inf;
      if (k == ke) b0.x3f(m,k+1,j,i) = bhl_.b3_inf;
    });
  }

  if (pmbp->padm != nullptr) {
    pmbp->padm->SetADMVariables(pmbp);
    pmbp->pdyngr->PrimToConInit(is, ie, js, je, ks, ke);
  } else if (is_mhd) {
    pmbp->pmhd->peos->PrimToCons(w0, bcc, u0, is, ie, js, je, ks, ke);
  } else {
    pmbp->phydro->peos->PrimToCons(w0, u0, is, ie, js, je, ks, ke);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void BHLWindInnerX1()
//! \brief Hold any x1 user boundary fixed to the asymptotic wind state.

void BHLWindInnerX1(Mesh *pm) {
  auto &indcs = pm->mb_indcs;
  int ng = indcs.ng;
  int n1 = indcs.nx1 + 2*ng;
  int n2 = (indcs.nx2 > 1) ? indcs.nx2 + 2*ng : 1;
  int n3 = (indcs.nx3 > 1) ? indcs.nx3 + 2*ng : 1;
  int is = indcs.is, ie = indcs.ie;
  int js = indcs.js, je = indcs.je;
  int ks = indcs.ks, ke = indcs.ke;
  int nmb = pm->pmb_pack->nmb_thispack;
  auto &mb_bcs = pm->pmb_pack->pmb->mb_bcs;
  const bool is_mhd = (pm->pmb_pack->pmhd != nullptr);

  DvceArray5D<Real> w0, u0, bcc;
  if (is_mhd) {
    w0 = pm->pmb_pack->pmhd->w0;
    u0 = pm->pmb_pack->pmhd->u0;
    bcc = pm->pmb_pack->pmhd->bcc0;
  } else {
    w0 = pm->pmb_pack->phydro->w0;
    u0 = pm->pmb_pack->phydro->u0;
  }
  auto bhl_ = bhl;
  Real gm1 = is_mhd ? pm->pmb_pack->pmhd->peos->eos_data.gamma - 1.0
                    : pm->pmb_pack->phydro->peos->eos_data.gamma - 1.0;

  par_for("bhl_wind_x1", DevExeSpace(), 0,nmb-1,0,n3-1,0,n2-1,0,ng-1,
  KOKKOS_LAMBDA(int m, int k, int j, int n) {
    if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
      SetWindPrimitives(bhl_, is_mhd, gm1, w0, bcc, m, k, j, is-ng+n);
    }
    if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
      SetWindPrimitives(bhl_, is_mhd, gm1, w0, bcc, m, k, j, ie+1+n);
    }
  });

  if (is_mhd) {
    auto &b0 = pm->pmb_pack->pmhd->b0;
    par_for("bhl_bfield_x1", DevExeSpace(), 0,nmb-1,0,n3-1,0,n2-1,0,ng-1,
    KOKKOS_LAMBDA(int m, int k, int j, int n) {
      if (mb_bcs.d_view(m,BoundaryFace::inner_x1) == BoundaryFlag::user) {
        int i = is-ng+n;
        b0.x1f(m,k,j,i) = bhl_.b1_inf;
        b0.x2f(m,k,j,i) = bhl_.b2_inf;
        b0.x3f(m,k,j,i) = bhl_.b3_inf;
      }
      if (mb_bcs.d_view(m,BoundaryFace::outer_x1) == BoundaryFlag::user) {
        int i = ie+1+n;
        b0.x1f(m,k,j,i+1) = bhl_.b1_inf;
        b0.x2f(m,k,j,i) = bhl_.b2_inf;
        b0.x3f(m,k,j,i) = bhl_.b3_inf;
      }
    });
  }

  if (pm->pmb_pack->padm != nullptr) {
    pm->pmb_pack->pdyngr->PrimToConInit(is-ng, is-1, 0, n2-1, 0, n3-1);
    pm->pmb_pack->pdyngr->PrimToConInit(ie+1, ie+ng, 0, n2-1, 0, n3-1);
  } else if (is_mhd) {
    pm->pmb_pack->pmhd->peos->PrimToCons(w0, bcc, u0, is-ng, is-1, 0, n2-1, 0, n3-1);
    pm->pmb_pack->pmhd->peos->PrimToCons(w0, bcc, u0, ie+1, ie+ng, 0, n2-1, 0, n3-1);
  } else {
    pm->pmb_pack->phydro->peos->PrimToCons(w0, u0, is-ng, is-1, 0, n2-1, 0, n3-1);
    pm->pmb_pack->phydro->peos->PrimToCons(w0, u0, ie+1, ie+ng, 0, n2-1, 0, n3-1);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void BHLFluxes()
//! \brief Surface accretion fluxes and BHL drag diagnostics.

void BHLFluxes(HistoryData *pdata, Mesh *pm) {
  MeshBlockPack *pmbp = pm->pmb_pack;
  bool flat = pmbp->pcoord->coord_data.is_minkowski;
  Real spin = pmbp->pcoord->coord_data.bh_spin;
  bool is_mhd = (pmbp->pmhd != nullptr);

  Real gamma = is_mhd ? pmbp->pmhd->peos->eos_data.gamma
                      : pmbp->phydro->peos->eos_data.gamma;
  Real to_ien = (pmbp->pdyngr != nullptr) ? 1.0/(gamma - 1.0) : 1.0;

  int nvars = is_mhd ? pmbp->pmhd->nmhd + pmbp->pmhd->nscalars
                     : pmbp->phydro->nhydro + pmbp->phydro->nscalars;
  DvceArray5D<Real> w0 = is_mhd ? pmbp->pmhd->w0 : pmbp->phydro->w0;
  DvceArray5D<Real> bcc0;
  if (is_mhd) bcc0 = pmbp->pmhd->bcc0;

  auto &grids = pm->pgen->spherical_grids;
  const int phi_offset = 4;
  const int drag_offset = is_mhd ? 6 : 4;
  const int nflux = drag_offset + 9;
  pdata->nhist = grids.size()*nflux;
  if (pdata->nhist > NHISTORY_VARIABLES) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "BHL user history requested " << pdata->nhist
              << " columns, larger than NHISTORY_VARIABLES=" << NHISTORY_VARIABLES
              << ". Reduce the number of history radii or increase NHISTORY_VARIABLES."
              << std::endl;
    exit(EXIT_FAILURE);
  }
  for (int g=0; g<static_cast<int>(grids.size()); ++g) {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(1) << grids[g]->radius;
    std::string rad = stream.str();
    pdata->label[nflux*g+0] = "mdot_" + rad;
    pdata->label[nflux*g+1] = "mdotHL_" + rad;
    pdata->label[nflux*g+2] = "edot_" + rad;
    pdata->label[nflux*g+3] = "jdot_" + rad;
    if (is_mhd) {
      pdata->label[nflux*g+phi_offset+0] = "phi_" + rad;
      pdata->label[nflux*g+phi_offset+1] = "varphi_" + rad;
    }
    pdata->label[nflux*g+drag_offset+0] = "fm_x_" + rad;
    pdata->label[nflux*g+drag_offset+1] = "fm_y_" + rad;
    pdata->label[nflux*g+drag_offset+2] = "fm_z_" + rad;
    pdata->label[nflux*g+drag_offset+3] = "fg_x_" + rad;
    pdata->label[nflux*g+drag_offset+4] = "fg_y_" + rad;
    pdata->label[nflux*g+drag_offset+5] = "fg_z_" + rad;
    pdata->label[nflux*g+drag_offset+6] = "fd_x_" + rad;
    pdata->label[nflux*g+drag_offset+7] = "fd_y_" + rad;
    pdata->label[nflux*g+drag_offset+8] = "fd_z_" + rad;
  }

  DualArray2D<Real> interpolated_bcc;
  for (int g=0; g<static_cast<int>(grids.size()); ++g) {
    for (int n=0; n<nflux; ++n) {
      pdata->hdata[nflux*g+n] = 0.0;
    }

    if (is_mhd) {
      grids[g]->InterpolateToSphere(3, bcc0);
      Kokkos::realloc(interpolated_bcc, grids[g]->nangles, 3);
      Kokkos::deep_copy(interpolated_bcc, grids[g]->interp_vals);
      interpolated_bcc.template modify<DevExeSpace>();
      interpolated_bcc.template sync<HostMemSpace>();
    }
    grids[g]->InterpolateToSphere(nvars, w0);

    for (int n=0; n<grids[g]->nangles; ++n) {
      Real r = grids[g]->radius;
      Real theta = grids[g]->polar_pos.h_view(n,0);
      Real phi = grids[g]->polar_pos.h_view(n,1);
      Real x1 = grids[g]->interp_coord.h_view(n,0);
      Real x2 = grids[g]->interp_coord.h_view(n,1);
      Real x3 = grids[g]->interp_coord.h_view(n,2);
      Real glower[4][4], gupper[4][4];
      ComputeMetricAndInverse(x1, x2, x3, flat, spin, glower, gupper);

      Real rho = grids[g]->interp_vals.h_view(n,IDN);
      Real eint = grids[g]->interp_vals.h_view(n,IEN)*to_ien;
      Real vx = grids[g]->interp_vals.h_view(n,IVX);
      Real vy = grids[g]->interp_vals.h_view(n,IVY);
      Real vz = grids[g]->interp_vals.h_view(n,IVZ);
      Real bx = 0.0, by = 0.0, bz = 0.0;
      if (is_mhd) {
        bx = interpolated_bcc.h_view(n,IBX);
        by = interpolated_bcc.h_view(n,IBY);
        bz = interpolated_bcc.h_view(n,IBZ);
      }

      Real q = glower[1][1]*vx*vx + 2.0*glower[1][2]*vx*vy +
               2.0*glower[1][3]*vx*vz + glower[2][2]*vy*vy +
               2.0*glower[2][3]*vy*vz + glower[3][3]*vz*vz;
      Real alpha = std::sqrt(-1.0/gupper[0][0]);
      Real lor = std::sqrt(1.0 + q);
      Real u0 = lor/alpha;
      Real u1 = vx - alpha*lor*gupper[0][1];
      Real u2 = vy - alpha*lor*gupper[0][2];
      Real u3 = vz - alpha*lor*gupper[0][3];

      Real u_0 = glower[0][0]*u0 + glower[0][1]*u1 + glower[0][2]*u2 + glower[0][3]*u3;
      Real u_1 = glower[1][0]*u0 + glower[1][1]*u1 + glower[1][2]*u2 + glower[1][3]*u3;
      Real u_2 = glower[2][0]*u0 + glower[2][1]*u1 + glower[2][2]*u2 + glower[2][3]*u3;
      Real u_3 = glower[3][0]*u0 + glower[3][1]*u1 + glower[3][2]*u2 + glower[3][3]*u3;

      Real a2 = SQR(spin);
      Real rad2 = SQR(x1) + SQR(x2) + SQR(x3);
      Real r2 = SQR(r);
      Real sth = std::sin(theta);
      Real sph = std::sin(phi);
      Real cph = std::cos(phi);
      Real drdx = r*x1/(2.0*r2 - rad2 + a2);
      Real drdy = r*x2/(2.0*r2 - rad2 + a2);
      Real drdz = (r*x3 + a2*x3/r)/(2.0*r2 - rad2 + a2);
      Real ur = drdx*u1 + drdy*u2 + drdz*u3;
      Real sqrtmdet = r2 + SQR(spin*std::cos(theta));
      Real domega = grids[g]->solid_angles.h_view(n);

      Real b0 = u_1*bx + u_2*by + u_3*bz;
      Real b1 = (u0 != 0.0) ? (bx + b0*u1)/u0 : 0.0;
      Real b2 = (u0 != 0.0) ? (by + b0*u2)/u0 : 0.0;
      Real b3 = (u0 != 0.0) ? (bz + b0*u3)/u0 : 0.0;
      Real b_0 = glower[0][0]*b0 + glower[0][1]*b1 + glower[0][2]*b2 + glower[0][3]*b3;
      Real b_1 = glower[1][0]*b0 + glower[1][1]*b1 + glower[1][2]*b2 + glower[1][3]*b3;
      Real b_2 = glower[2][0]*b0 + glower[2][1]*b1 + glower[2][2]*b2 + glower[2][3]*b3;
      Real b_3 = glower[3][0]*b0 + glower[3][1]*b1 + glower[3][2]*b2 + glower[3][3]*b3;
      Real b_sq = b0*b_0 + b1*b_1 + b2*b_2 + b3*b_3;
      Real br = drdx*b1 + drdy*b2 + drdz*b3;
      Real u_ph = (-r*sph - spin*cph)*sth*u_1 + (r*cph - spin*sph)*sth*u_2;
      Real b_ph = (-r*sph - spin*cph)*sth*b_1 + (r*cph - spin*sph)*sth*b_2;
      Real wtot = rho + gamma*eint + b_sq;

      pdata->hdata[nflux*g+0] += -rho*ur*sqrtmdet*domega;
      pdata->hdata[nflux*g+2] += (wtot*ur*u_0 - br*b_0)*sqrtmdet*domega;
      pdata->hdata[nflux*g+3] += (wtot*ur*u_ph - br*b_ph)*sqrtmdet*domega;
      pdata->hdata[nflux*g+drag_offset+0] += (wtot*ur*u_1 - br*b_1)*sqrtmdet*domega;
      pdata->hdata[nflux*g+drag_offset+1] += (wtot*ur*u_2 - br*b_2)*sqrtmdet*domega;
      pdata->hdata[nflux*g+drag_offset+2] += (wtot*ur*u_3 - br*b_3)*sqrtmdet*domega;

      if (is_mhd) {
        pdata->hdata[nflux*g+phi_offset] +=
            0.5*std::abs(br*u0 - b0*ur)*sqrtmdet*domega;
      }
    }

    pdata->hdata[nflux*g+1] = pdata->hdata[nflux*g+0]/bhl.mdot_hl;
    if (is_mhd) {
      Real mdot = pdata->hdata[nflux*g+0];
      pdata->hdata[nflux*g+phi_offset+1] =
          (mdot > 0.0) ? pdata->hdata[nflux*g+phi_offset]/std::sqrt(mdot) : 0.0;
    }

    if (bhl.grav_drag) {
      auto &indcs = pmbp->pmesh->mb_indcs;
      auto &size = pmbp->pmb->mb_size;
      int is = indcs.is, js = indcs.js, ks = indcs.ks;
      int nx1 = indcs.nx1, nx2 = indcs.nx2, nx3 = indcs.nx3;
      int nmkji = pmbp->nmb_thispack*nx3*nx2*nx1;
      int nkji = nx3*nx2*nx1;
      int nji = nx2*nx1;
      Real r_exclude = (bhl.drag_exclude_radius > 0.0) ? bhl.drag_exclude_radius
                                                       : grids[g]->radius;
      array_sum::GlobalSum grav_sum;
      Kokkos::parallel_reduce("bhl_grav_drag", Kokkos::RangePolicy<>(DevExeSpace(), 0, nmkji),
      KOKKOS_LAMBDA(const int &idx, array_sum::GlobalSum &mb_sum) {
        int m = idx/nkji;
        int k = (idx - m*nkji)/nji;
        int j = (idx - m*nkji - k*nji)/nx1;
        int i = (idx - m*nkji - k*nji - j*nx1) + is;
        k += ks;
        j += js;

        Real x1min = size.d_view(m).x1min;
        Real x1max = size.d_view(m).x1max;
        Real x2min = size.d_view(m).x2min;
        Real x2max = size.d_view(m).x2max;
        Real x3min = size.d_view(m).x3min;
        Real x3max = size.d_view(m).x3max;
        Real x1v = CellCenterX(i-is, nx1, x1min, x1max);
        Real x2v = CellCenterX(j-js, nx2, x2min, x2max);
        Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);
        Real r_ks = KerrSchildRadius(x1v, x2v, x3v, spin);
        Real r2 = SQR(x1v) + SQR(x2v) + SQR(x3v);

        array_sum::GlobalSum hvars;
        if (r_ks >= r_exclude && r2 > 0.0) {
          Real vol = size.d_view(m).dx1*size.d_view(m).dx2*size.d_view(m).dx3;
          Real fac = w0(m,IDN,k,j,i)*vol/(r2*std::sqrt(r2));
          hvars.the_array[0] = fac*x1v;
          hvars.the_array[1] = fac*x2v;
          hvars.the_array[2] = fac*x3v;
        }
        for (int n=3; n<NREDUCTION_VARIABLES; ++n) {
          hvars.the_array[n] = 0.0;
        }
        mb_sum += hvars;
      }, Kokkos::Sum<array_sum::GlobalSum>(grav_sum));
      Kokkos::fence();

      pdata->hdata[nflux*g+drag_offset+3] = grav_sum.the_array[0];
      pdata->hdata[nflux*g+drag_offset+4] = grav_sum.the_array[1];
      pdata->hdata[nflux*g+drag_offset+5] = grav_sum.the_array[2];
    }

    pdata->hdata[nflux*g+drag_offset+6] =
        -pdata->hdata[nflux*g+drag_offset+0] + pdata->hdata[nflux*g+drag_offset+3];
    pdata->hdata[nflux*g+drag_offset+7] =
        -pdata->hdata[nflux*g+drag_offset+1] + pdata->hdata[nflux*g+drag_offset+4];
    pdata->hdata[nflux*g+drag_offset+8] =
        -pdata->hdata[nflux*g+drag_offset+2] + pdata->hdata[nflux*g+drag_offset+5];
  }

  for (int n=pdata->nhist; n<NHISTORY_VARIABLES; ++n) {
    pdata->hdata[n] = 0.0;
  }
}
