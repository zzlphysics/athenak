#ifndef PGEN_EMRI_SRMHD_CHARACTERISTICS_HPP_
#define PGEN_EMRI_SRMHD_CHARACTERISTICS_HPP_
//========================================================================================
// AthenaK astrophysical plasma code
// Copyright(C) 2026 AthenaK collaboration
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file emri_srmhd_characteristics.hpp
//! \brief Device-portable seven-wave primitive SRMHD boundary projection.

#include <cmath>

#include "athena.hpp"

namespace emri_srmhd {

constexpr int kModes = 7;

KOKKOS_INLINE_FUNCTION
Real QuadraticRoot(Real a1, Real a0, bool greater) {
  if (a1*a1 < 4.0*a0) return -0.5*a1;
  if (greater) {
    return (a1 >= 0.0) ? -2.0*a0/(a1 + sqrt(a1*a1 - 4.0*a0))
                       : 0.5*(-a1 + sqrt(a1*a1 - 4.0*a0));
  }
  return (a1 >= 0.0) ? 0.5*(-a1 - sqrt(a1*a1 - 4.0*a0))
                     : -2.0*a0/(a1 - sqrt(a1*a1 - 4.0*a0));
}

KOKKOS_INLINE_FUNCTION
Real CubicRootReal(Real a2, Real a1, Real a0) {
  const Real q = (a2*a2 - 3.0*a1)/9.0;
  const Real r = (2.0*a2*a2*a2 - 9.0*a1*a2 + 27.0*a0)/54.0;
  if (r*r < q*q*q && q > 0.0) {
    const Real argument = fmax(-1.0, fmin(1.0, r/sqrt(q*q*q)));
    return -2.0*sqrt(q)*cos(acos(argument)/3.0) - a2/3.0;
  }
  const Real radical = fmax(r*r - q*q*q, Real(0.0));
  const Real a = -copysign(1.0, r)*cbrt(fabs(r) + sqrt(radical));
  const Real b = (a != 0.0) ? q/a : 0.0;
  return a + b - a2/3.0;
}

KOKKOS_INLINE_FUNCTION
void QuarticRoots(Real a3, Real a2, Real a1, Real a0, Real &x1, Real &x2,
                  Real &x3, Real &x4) {
  const Real b2 = a2 - 3.0/8.0*a3*a3;
  const Real b1 = a1 - 0.5*a2*a3 + 1.0/8.0*a3*a3*a3;
  const Real b0 = a0 - 0.25*a1*a3 + 1.0/16.0*a2*a3*a3
                        - 3.0/256.0*a3*a3*a3*a3;
  const Real z0 = CubicRootReal(-b2, -4.0*b0, 4.0*b0*b2 - b1*b1);
  const Real d1 = (z0 > b2) ? sqrt(z0 - b2) : 0.0;
  const Real e1 = -d1;
  const Real inner = sqrt(fmax(0.25*z0*z0 - b0, Real(0.0)));
  const Real d0 = 0.5*z0 + ((b1 < 0.0) ? inner : -inner);
  const Real e0 = 0.5*z0 - ((b1 < 0.0) ? inner : -inner);
  const Real y1 = QuadraticRoot(d1, d0, false);
  const Real y2 = QuadraticRoot(d1, d0, true);
  const Real y3 = QuadraticRoot(e1, e0, false);
  const Real y4 = QuadraticRoot(e1, e0, true);
  x1 = fmin(y1, y3) - 0.25*a3;
  x2 = fmin(fmax(y1, y3), fmin(y2, y4)) - 0.25*a3;
  x3 = fmax(fmax(y1, y3), fmin(y2, y4)) - 0.25*a3;
  x4 = fmax(y2, y4) - 0.25*a3;
}

//! Construct one column of the primitive basis q=(rho,p,u_n,u_u,u_v,B_u,B_v).
KOKKOS_INLINE_FUNCTION
bool Eigenvector(const Real q[kModes], Real normal_field, Real gamma_adi, int mode,
                 Real &lambda, Real column[kModes]) {
  const Real rho = q[0];
  const Real pgas = q[1];
  const Real u0 = sqrt(1.0 + q[2]*q[2] + q[3]*q[3] + q[4]*q[4]);
  const Real velocity[3] = {q[2]/u0, q[3]/u0, q[4]/u0};
  const Real magnetic[3] = {normal_field, q[5], q[6]};
  Real u[4] = {u0, q[2], q[3], q[4]};
  Real b[4];
  b[0] = magnetic[0]*u[1] + magnetic[1]*u[2] + magnetic[2]*u[3];
  for (int component = 1; component < 4; ++component) {
    b[component] = (magnetic[component - 1] + b[0]*u[component])/u[0];
  }
  const Real wgas = rho + gamma_adi*pgas/(gamma_adi - 1.0);
  const Real cs_sq = gamma_adi*pgas/wgas;
  const Real b_sq = -b[0]*b[0] + b[1]*b[1] + b[2]*b[2] + b[3]*b[3];
  const Real wtot = wgas + b_sq;
  Real delta_rho = 0.0;
  Real delta_pgas = 0.0;
  Real delta_u[4] = {0.0, 0.0, 0.0, 0.0};
  Real delta_b[4] = {0.0, 0.0, 0.0, 0.0};

  if (mode == 3) {
    lambda = velocity[0];
    delta_rho = 1.0;
  } else if (mode == 1 || mode == 5) {
    const Real root_wtot = sqrt(wtot);
    const Real lambda_ap = (b[1] + root_wtot*u[1])/(b[0] + root_wtot*u[0]);
    const Real lambda_am = (b[1] - root_wtot*u[1])/(b[0] - root_wtot*u[0]);
    Real sign = 1.0;
    if ((lambda_ap > lambda_am && mode == 1) ||
        (lambda_ap <= lambda_am && mode == 5)) sign = -1.0;
    lambda = (sign > 0.0) ? lambda_ap : lambda_am;
    const Real alpha1[4] = {u[3], lambda*u[3], 0.0, u[0] - lambda*u[1]};
    const Real alpha2[4] = {-u[2], -lambda*u[2], lambda*u[1] - u[0], 0.0};
    const Real denominator = 1.0 - lambda*velocity[0];
    const Real g1 = (magnetic[1] + lambda*velocity[1]/denominator*magnetic[0])/u[0];
    const Real g2 = (magnetic[2] + lambda*velocity[2]/denominator*magnetic[0])/u[0];
    Real f1 = 1.0/sqrt(2.0);
    Real f2 = f1;
    if (g1 != 0.0 || g2 != 0.0) {
      const Real norm = sqrt(g1*g1 + g2*g2);
      f1 = g1/norm;
      f2 = g2/norm;
    }
    for (int mu = 0; mu < 4; ++mu) {
      delta_u[mu] = f1*alpha1[mu] + f2*alpha2[mu];
      delta_b[mu] = -sign*root_wtot*delta_u[mu];
    }
  } else {
    const Real factor_a = wgas*(1.0/cs_sq - 1.0);
    const Real factor_b = -(wgas + b_sq/cs_sq);
    const Real gamma2 = u[0]*u[0];
    const Real gamma4 = gamma2*gamma2;
    const Real coeff4 = factor_a*gamma4 - factor_b*gamma2 - b[0]*b[0];
    const Real coeff3 = -4.0*factor_a*gamma4*velocity[0]
                        + 2.0*factor_b*gamma2*velocity[0] + 2.0*b[0]*b[1];
    const Real coeff2 = 6.0*factor_a*gamma4*velocity[0]*velocity[0]
                        + factor_b*gamma2*(1.0 - velocity[0]*velocity[0])
                        + b[0]*b[0] - b[1]*b[1];
    const Real coeff1 = -4.0*factor_a*gamma4*velocity[0]*velocity[0]*velocity[0]
                        - 2.0*factor_b*gamma2*velocity[0] - 2.0*b[0]*b[1];
    const Real coeff0 = factor_a*gamma4*velocity[0]*velocity[0]
                            *velocity[0]*velocity[0]
                        + factor_b*gamma2*velocity[0]*velocity[0] + b[1]*b[1];
    if (fabs(coeff4) < 1.0e-30) return false;
    Real lambda_fl, lambda_sl, lambda_sr, lambda_fr;
    QuarticRoots(coeff3/coeff4, coeff2/coeff4, coeff1/coeff4, coeff0/coeff4,
                 lambda_fl, lambda_sl, lambda_sr, lambda_fr);
    Real lambda_other = 0.0;
    if (mode == 0) { lambda = lambda_fl; lambda_other = lambda_sl; }
    if (mode == 2) { lambda = lambda_sl; lambda_other = lambda_fl; }
    if (mode == 4) { lambda = lambda_sr; lambda_other = lambda_fr; }
    if (mode == 6) { lambda = lambda_fr; lambda_other = lambda_sr; }
    const Real root_wtot = sqrt(wtot);
    const Real lambda_ap = (b[1] + root_wtot*u[1])/(b[0] + root_wtot*u[0]);
    const Real lambda_am = (b[1] - root_wtot*u[1])/(b[0] - root_wtot*u[0]);
    Real lambda_a = lambda_ap;
    Real sign = 1.0;
    if (lambda_ap > lambda_am) {
      if (mode < 3) { lambda_a = lambda_am; sign = -1.0; }
    } else if (mode > 3) {
      lambda_a = lambda_am;
      sign = -1.0;
    }
    const Real a = u[0]*(velocity[0] - lambda);
    const Real g = 1.0 - lambda*lambda;
    const Real radicand = -factor_b - factor_a*a*a/g;
    if (!(g > 0.0) || radicand < 0.0) return false;
    const Real b_over_a = -sign*sqrt(radicand);
    const Real alpha1[4] = {u[3], lambda*u[3], 0.0, u[0] - lambda*u[1]};
    const Real alpha2[4] = {-u[2], -lambda*u[2], lambda*u[1] - u[0], 0.0};
    Real alpha11 = -alpha1[0]*alpha1[0];
    Real alpha12 = -alpha1[0]*alpha2[0];
    Real alpha22 = -alpha2[0]*alpha2[0];
    for (int mu = 1; mu < 4; ++mu) {
      alpha11 += alpha1[mu]*alpha1[mu];
      alpha12 += alpha1[mu]*alpha2[mu];
      alpha22 += alpha2[mu]*alpha2[mu];
    }
    const Real velocity_denom = 1.0 - lambda*velocity[0];
    const Real g1 = (magnetic[1] + lambda*velocity[1]/velocity_denom*magnetic[0])/u[0];
    const Real g2 = (magnetic[2] + lambda*velocity[2]/velocity_denom*magnetic[0])/u[0];
    const Real determinant = alpha11*alpha22 - alpha12*alpha12;
    if (fabs(determinant) < 1.0e-30) return false;
    const Real c1 = (g1*alpha12 + g2*alpha22)/determinant*u[0]*velocity_denom;
    const Real c2 = -(g1*alpha11 + g2*alpha12)/determinant*u[0]*velocity_denom;
    Real bt[4];
    for (int mu = 0; mu < 4; ++mu) bt[mu] = c1*alpha1[mu] + c2*alpha2[mu];
    Real f1 = 1.0/sqrt(2.0);
    Real f2 = f1;
    if (g1 != 0.0 || g2 != 0.0) {
      const Real norm = sqrt(g1*g1 + g2*g2);
      f1 = g1/norm;
      f2 = g2/norm;
    }
    Real phi[4];
    for (int mu = 0; mu < 4; ++mu) phi[mu] = a*u[mu];
    phi[0] += lambda;
    phi[1] += 1.0;
    if (fabs(lambda - lambda_a) <= fabs(lambda_other - lambda_a)) {
      const Real normalization_squared = determinant*
          (f1*f1*alpha11 + 2.0*f1*f2*alpha12 + f2*f2*alpha22);
      if (!(normalization_squared > 0.0)) return false;
      const Real normalization = sqrt(normalization_squared);
      Real bt_normalized[4];
      for (int mu = 0; mu < 4; ++mu) {
        bt_normalized[mu] = ((f1*alpha12 + f2*alpha22)*alpha1[mu]
                            - (f1*alpha11 + f2*alpha12)*alpha2[mu])/normalization;
      }
      Real bt_norm_squared = -bt[0]*bt[0];
      for (int mu = 1; mu < 4; ++mu) bt_norm_squared += bt[mu]*bt[mu];
      const Real denominator = a*a - (g + a*a)*cs_sq;
      if (bt_norm_squared < 0.0 || fabs(denominator) < 1.0e-30) return false;
      delta_pgas = -(g + a*a)*cs_sq/denominator*sqrt(bt_norm_squared);
      delta_rho = rho/(gamma_adi*pgas)*delta_pgas;
      for (int mu = 0; mu < 4; ++mu) {
        delta_u[mu] = -a*delta_pgas/(wgas*cs_sq*(g + a*a))*phi[mu]
                      - b_over_a/wgas*bt_normalized[mu];
        delta_b[mu] = -b_over_a*delta_pgas/wgas*u[mu]
                      - (1.0 + a*a/g)*bt_normalized[mu];
      }
    } else {
      delta_pgas = -1.0;
      delta_rho = rho/(gamma_adi*pgas)*delta_pgas;
      const Real denominator = wgas*a*a - b_sq*g;
      Real bt_reduced[4] = {0.0, 0.0, 0.0, 0.0};
      if (fabs(denominator) > 1.0e-30) {
        for (int mu = 0; mu < 4; ++mu) bt_reduced[mu] = bt[mu]/denominator;
      }
      for (int mu = 0; mu < 4; ++mu) {
        delta_u[mu] = a/(wgas*cs_sq*(g + a*a))*phi[mu]
                      - b_over_a*g/wgas*bt_reduced[mu];
        delta_b[mu] = b_over_a/wgas*u[mu]
                      - (1.0 + a*a/g)*g*bt_reduced[mu];
      }
    }
  }

  column[0] = delta_rho;
  column[1] = delta_pgas;
  column[2] = delta_u[1];
  column[3] = delta_u[2];
  column[4] = delta_u[3];
  column[5] = b[2]*delta_u[0] - b[0]*delta_u[2]
              + delta_b[2]*u[0] - delta_b[0]*u[2];
  column[6] = b[3]*delta_u[0] - b[0]*delta_u[3]
              + delta_b[3]*u[0] - delta_b[0]*u[3];
  Real norm_squared = 0.0;
  for (int variable = 0; variable < kModes; ++variable) {
    if (!isfinite(column[variable])) return false;
    norm_squared += column[variable]*column[variable];
  }
  if (!(norm_squared > 0.0) || !isfinite(lambda)) return false;
  const Real inverse_norm = 1.0/sqrt(norm_squared);
  for (int variable = 0; variable < kModes; ++variable) {
    column[variable] *= inverse_norm;
  }
  return true;
}

KOKKOS_INLINE_FUNCTION
bool Invert(const Real matrix[kModes][kModes], Real inverse[kModes][kModes]) {
  Real work[kModes][2*kModes];
  for (int row = 0; row < kModes; ++row) {
    for (int column = 0; column < kModes; ++column) {
      work[row][column] = matrix[row][column];
      work[row][column + kModes] = (row == column) ? 1.0 : 0.0;
    }
  }
  for (int column = 0; column < kModes; ++column) {
    int pivot = column;
    Real pivot_size = fabs(work[pivot][column]);
    for (int row = column + 1; row < kModes; ++row) {
      const Real size = fabs(work[row][column]);
      if (size > pivot_size) { pivot = row; pivot_size = size; }
    }
    if (!(pivot_size > 1.0e-12)) return false;
    if (pivot != column) {
      for (int entry = 0; entry < 2*kModes; ++entry) {
        const Real temporary = work[column][entry];
        work[column][entry] = work[pivot][entry];
        work[pivot][entry] = temporary;
      }
    }
    const Real inverse_pivot = 1.0/work[column][column];
    for (int entry = 0; entry < 2*kModes; ++entry) {
      work[column][entry] *= inverse_pivot;
    }
    for (int row = 0; row < kModes; ++row) {
      if (row == column) continue;
      const Real factor = work[row][column];
      for (int entry = 0; entry < 2*kModes; ++entry) {
        work[row][entry] -= factor*work[column][entry];
      }
    }
  }
  for (int row = 0; row < kModes; ++row) {
    for (int column = 0; column < kModes; ++column) {
      inverse[row][column] = work[row][column + kModes];
    }
  }
  return true;
}

//! Replace only modes entering the inner domain (negative outward-frame speed).
KOKKOS_INLINE_FUNCTION
bool ProjectIncoming(const Real interior[kModes], const Real exterior[kModes],
                     Real normal_field, Real gamma_adi, Real speed_tolerance,
                     Real outward_grid_speed, Real density_floor, Real pressure_floor,
                     Real result[kModes]) {
  if (!(interior[0] > 0.0 && interior[1] > 0.0 && exterior[0] > 0.0 &&
        exterior[1] > 0.0)) return false;
  Real reference[kModes];
  reference[0] = sqrt(interior[0]*exterior[0]);
  reference[1] = sqrt(interior[1]*exterior[1]);
  for (int variable = 2; variable < kModes; ++variable) {
    reference[variable] = 0.5*(interior[variable] + exterior[variable]);
  }
  Real right[kModes][kModes];
  Real speeds[kModes];
  for (int mode = 0; mode < kModes; ++mode) {
    Real column[kModes];
    if (!Eigenvector(reference, normal_field, gamma_adi, mode, speeds[mode], column)) {
      return false;
    }
    for (int variable = 0; variable < kModes; ++variable) {
      right[variable][mode] = column[variable];
    }
  }
  Real left[kModes][kModes];
  if (!Invert(right, left)) return false;
  Real amplitudes[kModes];
  for (int mode = 0; mode < kModes; ++mode) {
    amplitudes[mode] = 0.0;
    for (int variable = 0; variable < kModes; ++variable) {
      amplitudes[mode] += left[mode][variable]
                          *(exterior[variable] - interior[variable]);
    }
    if (speeds[mode] - outward_grid_speed >= -speed_tolerance) amplitudes[mode] = 0.0;
  }
  Real applied[kModes] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  for (int variable = 0; variable < kModes; ++variable) {
    for (int mode = 0; mode < kModes; ++mode) {
      applied[variable] += right[variable][mode]*amplitudes[mode];
    }
  }
  Real fraction = 1.0;
  const Real floors[2] = {density_floor, pressure_floor};
  for (int variable = 0; variable < 2; ++variable) {
    if (interior[variable] + applied[variable] < floors[variable] &&
        applied[variable] < 0.0) {
      fraction = fmin(fraction, 0.99*(interior[variable] - floors[variable])
                                /(-applied[variable]));
    }
  }
  fraction = fmax(0.0, fmin(1.0, fraction));
  for (int variable = 0; variable < kModes; ++variable) {
    result[variable] = interior[variable] + fraction*applied[variable];
    if (!isfinite(result[variable])) return false;
  }
  return result[0] > 0.0 && result[1] > 0.0;
}

}  // namespace emri_srmhd

#endif  // PGEN_EMRI_SRMHD_CHARACTERISTICS_HPP_
