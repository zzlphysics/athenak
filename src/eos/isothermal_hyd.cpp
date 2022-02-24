//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file isothermal_hyd.cpp
//! \brief derived class that implements isothermal EOS for nonrelativistic hydro

#include "athena.hpp"
#include "parameter_input.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"
#include "eos/eos.hpp"

//----------------------------------------------------------------------------------------
// ctor: also calls EOS base class constructor

IsothermalHydro::IsothermalHydro(MeshBlockPack *pp, ParameterInput *pin) :
    EquationOfState("hydro", pp, pin) {
  eos_data.is_ideal = false;
  eos_data.iso_cs = pin->GetReal("hydro","iso_sound_speed");
  eos_data.gamma = 0.0;
  eos_data.use_e = false;
  eos_data.use_t = false;
}

//----------------------------------------------------------------------------------------
//! \fn void ConsToPrim()
//! \brief Converts conserved into primitive variables. Operates over range of cells given
//! in argument list.

void IsothermalHydro::ConsToPrim(DvceArray5D<Real> &cons, DvceArray5D<Real> &prim,
                                 const int il, const int iu, const int jl, const int ju,
                                 const int kl, const int ku) {
  int &nhyd  = pmy_pack->phydro->nhydro;
  int &nscal = pmy_pack->phydro->nscalars;
  int &nmb = pmy_pack->nmb_thispack;

  Real &dfloor_ = eos_data.dfloor;

  par_for("isohyd_con2prim", DevExeSpace(), 0, (nmb-1), kl, ku, jl, ju, il, iu,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real& u_d  = cons(m,IDN,k,j,i);
    const Real& u_m1 = cons(m,IM1,k,j,i);
    const Real& u_m2 = cons(m,IM2,k,j,i);
    const Real& u_m3 = cons(m,IM3,k,j,i);

    Real& w_d  = prim(m,IDN,k,j,i);
    Real& w_vx = prim(m,IVX,k,j,i);
    Real& w_vy = prim(m,IVY,k,j,i);
    Real& w_vz = prim(m,IVZ,k,j,i);

    // apply density floor, without changing momentum or energy
    u_d = (u_d > dfloor_) ?  u_d : dfloor_;
    w_d = u_d;

    Real di = 1.0/u_d;
    w_vx = u_m1*di;
    w_vy = u_m2*di;
    w_vz = u_m3*di;

    // convert scalars (if any)
    for (int n=nhyd; n<(nhyd+nscal); ++n) {
      prim(m,n,k,j,i) = cons(m,n,k,j,i)/u_d;
    }
  });

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void PrimToCons()
//! \brief Converts primitive into conserved variables. Operates over range of cells given
//! in argument list.

void IsothermalHydro::PrimToCons(const DvceArray5D<Real> &prim, DvceArray5D<Real> &cons,
                                 const int il, const int iu, const int jl, const int ju,
                                 const int kl, const int ku) {
  int &nhyd  = pmy_pack->phydro->nhydro;
  int &nscal = pmy_pack->phydro->nscalars;
  int &nmb = pmy_pack->nmb_thispack;

  par_for("isohyd_prim2con", DevExeSpace(), 0, (nmb-1), kl, ku, jl, ju, il, iu,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    Real& u_d  = cons(m,IDN,k,j,i);
    Real& u_m1 = cons(m,IM1,k,j,i);
    Real& u_m2 = cons(m,IM2,k,j,i);
    Real& u_m3 = cons(m,IM3,k,j,i);

    const Real& w_d  = prim(m,IDN,k,j,i);
    const Real& w_vx = prim(m,IVX,k,j,i);
    const Real& w_vy = prim(m,IVY,k,j,i);
    const Real& w_vz = prim(m,IVZ,k,j,i);

    u_d  = w_d;
    u_m1 = w_vx*w_d;
    u_m2 = w_vy*w_d;
    u_m3 = w_vz*w_d;

    // convert scalars (if any)
    for (int n=nhyd; n<(nhyd+nscal); ++n) {
      cons(m,n,k,j,i) = prim(m,n,k,j,i)*w_d;
    }
  });

  return;
}
