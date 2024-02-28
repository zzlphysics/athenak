//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file radiation_femn_basis.cpp
//  \brief implementation of the radiation FEM/FPN basis functions and helpers

#include <iostream>
#include <gsl/gsl_sf_legendre.h>
#include "athena.hpp"
#include "radiation_femn/radiation_femn_geodesic_grid_matrices.hpp"

namespace radiationfemn {

/* Convert Barycentric coordinates to Cartesian coordinates given vertices of triangle
 *
 * Inputs:
 * (x1,y1,z1), (x2,y2,z2), (x3,y3,z3):  the three triangle vertices in cartesian coordinates
 * (xi1, xi2, xi3):the barycentric coordinates of a point inside the triangle
 *
 * Output:
 * (xval, yval, zval): the cartesian coordinates of the point
 */
inline void BarycentricToCartesian(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3,
                                   Real xi1, Real xi2, Real xi3, Real &xval, Real &yval, Real &zval) {

  xval = xi1 * x1 + xi2 * x2 + xi3 * x3;
  yval = xi1 * y1 + xi2 * y2 + xi3 * y3;
  zval = xi1 * z1 + xi2 * z2 + xi3 * z3;

}

/* Given index numbers of two vertices, finds if they share an edge and if so, return triangle info
 *
 * If a = b, this return all triangles that share the vertex
 *
 * Inputs:
 * a, b: index number of vertices
 * triangles: triangle information of the geodesic grid
 *
 * Outputs:
 * edge_triangles: the vertex information of shared edge(s)
 * is_edge: bool for if the vertices share an edge or not
 */
void FindTriangles(int a, int b, const HostArray2D<int> &triangles, HostArray2D<int> &edge_triangles, bool &is_edge) {

  is_edge = false;
  Kokkos::realloc(edge_triangles, 6, 3);
  Kokkos::deep_copy(edge_triangles, -42.);

  if (a == b) {
    size_t index{0};
    for (size_t i = 0; i < triangles.size() / 3; i++) {
      {
        if (triangles(i, 0) == a || triangles(i, 1) == a || triangles(i, 2) == a) {
          is_edge = true;
          edge_triangles(index, 0) = triangles(i, 0);
          edge_triangles(index, 1) = triangles(i, 1);
          edge_triangles(index, 2) = triangles(i, 2);
          index++;
        }
      }
    }
  } else if (a != b) {
    size_t index{0};
    for (size_t i = 0; i < triangles.size() / 3; i++) {
      if ((triangles(i, 0) == a && triangles(i, 1) == b) || (triangles(i, 0) == a && triangles(i, 2) == b) || (triangles(i, 0) == b && triangles(i, 1) == a)
          || (triangles(i, 0) == b && triangles(i, 2) == a) || (triangles(i, 1) == a && triangles(i, 2) == b) || (triangles(i, 1) == b && triangles(i, 2) == a)) {
        is_edge = true;
        edge_triangles(index, 0) = triangles(i, 0);
        edge_triangles(index, 1) = triangles(i, 1);
        edge_triangles(index, 2) = triangles(i, 2);
        index++;
      }
    }
  }
}


/* FEM basis functions (in barycentric coordinates)
 *
 * basis_choice: 1 => overlapping tent, 2 => small tent, 3 => overlapping honeycomb, 4 => small honeycomb
 * basis_index: 1 => basis peaked at xi1 = 1, 2 => basis peaked at xi2 = 1, 3 => basis peaked at xi3 = 1
 */
using fem_basis_type = Real (*)(Real, Real, Real);
const fem_basis_type fem_basis_fn[3][4] = {{fem_overtent_index1, fem_smalltent_index1, fem_overhoney_index1, fem_smallhoney_index1},
                                           {fem_overtent_index2, fem_smalltent_index2, fem_overhoney_index2, fem_smallhoney_index2},
                                           {fem_overtent_index3, fem_smalltent_index3, fem_overhoney_index3, fem_smallhoney_index3}};
inline Real fem_basis(Real xi1, Real xi2, Real xi3, int basis_index, int basis_choice) {
  return fem_basis_fn[basis_index - 1][basis_choice - 1](xi1, xi2, xi3);
}

// Overlapping tent basis
inline Real fem_overtent_index1(Real xi1, Real xi2, Real xi3) {
  return 2. * xi1 + xi2 + xi3 - 1.;
}
inline Real fem_overtent_index2(Real xi1, Real xi2, Real xi3) {
  return xi1 + 2. * xi2 + xi3 - 1.;
}
inline Real fem_overtent_index3(Real xi1, Real xi2, Real xi3) {
  return xi1 + xi2 + 2. * xi3 - 1.;
}

// Small tent basis
inline Real fem_smalltent_index1(Real xi1, Real xi2, Real xi3) {
  return (xi1 >= 0.5) * (xi1 - xi2 - xi3);
}
inline Real fem_smalltent_index2(Real xi1, Real xi2, Real xi3) {
  return (xi2 >= 0.5) * (xi2 - xi3 - xi1);
}
inline Real fem_smalltent_index3(Real xi1, Real xi2, Real xi3) {
  return (xi3 >= 0.5) * (xi3 - xi1 - xi2);
}

// Overlapping honeycomb
inline Real fem_overhoney_index1(Real xi1, Real xi2, Real xi3) {
  return 1.;
}
inline Real fem_overhoney_index2(Real xi1, Real xi2, Real xi3) {
  return 1.;
}
inline Real fem_overhoney_index3(Real xi1, Real xi2, Real xi3) {
  return 1.;
}

// Small honeycomb
inline Real fem_smallhoney_index1(Real xi1, Real xi2, Real xi3) {
  return (xi1 >= xi2) * (xi1 > xi3) * 1.;
}
inline Real fem_smallhoney_index2(Real xi1, Real xi2, Real xi3) {
  return (xi2 >= xi3) * (xi2 > xi1) * 1.;
}
inline Real fem_smallhoney_index3(Real xi1, Real xi2, Real xi3) {
  return (xi3 >= xi1) * (xi3 > xi2) * 1.;
}

// product of two FEM basis functions: psi_a psi_b
Real fem_basis_ab(int a, int b, int t1, int t2, int t3, Real xi1, Real xi2, Real xi3, int basis_choice) {
  int basis_index_a = (a == t1) * 1 + (a == t2) * 2 + (a == t3) * 3;
  int basis_index_b = (b == t1) * 1 + (b == t2) * 2 + (b == t3) * 3;
  return fem_basis(xi1, xi2, xi3, basis_index_a, basis_choice) * fem_basis(xi1, xi2, xi3, basis_index_b, basis_choice);
}

// a fem basis function: psi_a
Real fem_basis_a(int a, int t1, int t2, int t3, Real xi1, Real xi2, Real xi3, int basis_choice) {
  int basis_index_a = (a == t1) * 1 + (a == t2) * 2 + (a == t3) * 3;
  return fem_basis(xi1, xi2, xi3, basis_index_a, basis_choice);
}

/* Derivative of FEM basis wrt barycentric coordinates: dpsi/dxi
 *
 * basis_index: 1 => basis peaked at xi1 = 1, 2 => basis peaked at xi2 = 1, 3 => basis peaked at xi3 = 1
 * xi_index: 1 => xi1, 2 => xi2, 3 => xi3
 */
const double fem_basis_derivative[3][3] = {{2., 1., 1.},
                                           {1., 2., 1.},
                                           {1., 1., 2.}};
Real dfem_dxi(Real xi1, Real xi2, Real xi3, int basis_index, int xi_index) {
  return fem_basis_derivative[basis_index - 1][xi_index - 1];
}

/* FPN basis: real spherical harmonics
 *
 * Note: gsl_sf_legendre_sphPlm computes sqrt((2l+1)/(4 pi)) sqrt((l-m)!/(l+m)!) P_m^l(x)
 */
Real fpn_basis_lm(int l, int m, Real phi, Real theta) {
  Real result = 0.;
  if (m > 0) {
    result = sqrt(2.) * cos(m * phi) * gsl_sf_legendre_sphPlm(l, m, cos(theta));
  } else if (m == 0) {
    result = gsl_sf_legendre_sphPlm(l, 0, cos(theta));
  } else {
    result = sqrt(2.) * sin(abs(m) * phi) * gsl_sf_legendre_sphPlm(l, abs(m), cos(theta));
  }

  return result;
}

// FPN basis derivatives: dYlm/dphi = - m Yl-m
inline Real dfpn_dphi(int l, int m, Real phi, Real theta) {
  return -m * fpn_basis_lm(l, -m, phi, theta);
}

/* FPN basis derivatives: dYlm/dtheta
 * Note: Uses (1-x^2) dP^m_l/dx = (m-l-1) P^m_l+1 + (l+1)xP^m_l
 */
inline Real dfpn_dtheta(int l, int m, Real phi, Real theta) {
  Real result = 0.;

  Real der_legendre =
      -sin(theta) * ((m - l - 1.) * gsl_sf_legendre_sphPlm(l + 1, abs(m), cos(theta)) + (l + 1) * cos(theta) * gsl_sf_legendre_sphPlm(l, abs(m), cos(theta)))
          / (1. - cos(theta) * cos(theta));

  if (m > 0) {
    result = sqrt(2.) * cos(m * phi) * der_legendre;
  } else if (m == 0) {
    result = der_legendre;
  } else {
    result = sqrt(2.) * sin(abs(m) * phi) * der_legendre;
  }
  return result;
}

/* Derivative of FPN basis wrt angle
 *
 * var_index: 1 => phi, 2 => theta)
 */
using dfpn_dOmega_fn = Real (*)(int, int, Real, Real);
const dfpn_dOmega_fn dfpn_domega_fns[2] = {dfpn_dphi, dfpn_dtheta};
Real dfpn_dOmega(int l, int m, Real phi, Real theta, int var_index) {
  return dfpn_domega_fns[var_index - 1](l, m, phi, theta);
}

/* Jacobian P^Omega_ihat (energy = 1)
 *
 * Omega = (phi, theta), ihat = (x, y, z)
 */
Real PtildehatJac(Real phi, Real theta, int tilde_index, int hat_index) {
  if (tilde_index == 1 && hat_index == 1) {
    return -sin(phi) / sin(theta);
  } else if (tilde_index == 1 && hat_index == 2) {
    return cos(phi) / sin(theta);
  } else if (tilde_index == 1 && hat_index == 3) {
    return 0.;
  } else if (tilde_index == 2 && hat_index == 1) {
    return cos(phi) * cos(theta);
  } else if (tilde_index == 2 && hat_index == 2) {
    return sin(phi) * cos(theta);
  } else if (tilde_index == 2 && hat_index == 3) {
    return -sin(theta);
  } else {
    // std::cout << "Incorrect choice of index in radiation-femn block!" << std::endl;
    exit(EXIT_FAILURE);
  }
}

// -------------------------------------------------------------------------
// Cos Phi Sin Theta
//KOKKOS_INLINE_FUNCTION
Real CosPhiSinTheta(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3, Real xi1, Real xi2, Real xi3) {
  Real xval, yval, zval;
  BarycentricToCartesian(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3, xval, yval, zval);

  Real rval = sqrt(xval * xval + yval * yval + zval * zval);
  Real thetaval = acos(zval / rval);
  Real phival = atan2(yval, xval);

  return cos(phival) * sin(thetaval);
}

// ------------------------------------------------------------------------
// Sin Phi Sin Theta
//KOKKOS_INLINE_FUNCTION
Real SinPhiSinTheta(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3, Real xi1, Real xi2, Real xi3) {
  Real xval, yval, zval;
  BarycentricToCartesian(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3, xval, yval, zval);

  Real rval = sqrt(xval * xval + yval * yval + zval * zval);
  Real thetaval = acos(zval / rval);
  Real phival = atan2(yval, xval);

  return sin(phival) * sin(thetaval);
}

// ------------------------------------------------------------------------
// Cos Theta
//KOKKOS_INLINE_FUNCTION
Real CosTheta(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3, Real xi1, Real xi2, Real xi3) {
  Real xval, yval, zval;
  BarycentricToCartesian(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3, xval, yval, zval);

  Real rval = sqrt(xval * xval + yval * yval + zval * zval);
  Real thetaval = acos(zval / rval);

  return cos(thetaval);
}

// -------------------------------------------------------------------------
// sin Phi Cosec Theta
//KOKKOS_INLINE_FUNCTION
Real SinPhiCosecTheta(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3, Real xi1, Real xi2, Real xi3) {
  Real xval, yval, zval;
  BarycentricToCartesian(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3, xval, yval, zval);

  Real rval = sqrt(xval * xval + yval * yval + zval * zval);
  Real thetaval = acos(zval / rval);
  Real phival = atan2(yval, xval);

  return sin(phival) / sin(thetaval);
}

// -------------------------------------------------------------------------
// Cos Phi Cos Theta
//KOKKOS_INLINE_FUNCTION
Real CosPhiCosTheta(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3, Real xi1, Real xi2, Real xi3) {
  Real xval, yval, zval;
  BarycentricToCartesian(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3, xval, yval, zval);

  Real rval = sqrt(xval * xval + yval * yval + zval * zval);
  Real thetaval = acos(zval / rval);
  Real phival = atan2(yval, xval);

  return cos(phival) * cos(thetaval);
}

// ------------------------------------------------------------------------
// Cos Phi Cosec Theta
//KOKKOS_INLINE_FUNCTION
Real CosPhiCosecTheta(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3, Real xi1, Real xi2, Real xi3) {
  Real xval, yval, zval;
  BarycentricToCartesian(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3, xval, yval, zval);

  Real rval = sqrt(xval * xval + yval * yval + zval * zval);
  Real thetaval = acos(zval / rval);
  Real phival = atan2(yval, xval);

  return cos(phival) / sin(thetaval);
}

// ------------------------------------------------------------------------
// Sin Phi Cos Theta
//KOKKOS_INLINE_FUNCTION
Real SinPhiCosTheta(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3, Real xi1, Real xi2, Real xi3) {
  Real xval, yval, zval;
  BarycentricToCartesian(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3, xval, yval, zval);

  Real rval = sqrt(xval * xval + yval * yval + zval * zval);
  Real thetaval = acos(zval / rval);
  Real phival = atan2(yval, xval);

  return sin(phival) * cos(thetaval);
}

// ------------------------------------------------------------------------
// Sin Theta
//KOKKOS_INLINE_FUNCTION
Real SinTheta(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3, Real xi1, Real xi2, Real xi3) {
  Real xval, yval, zval;
  BarycentricToCartesian(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3, xval, yval, zval);

  Real rval = sqrt(xval * xval + yval * yval + zval * zval);
  Real thetaval = acos(zval / rval);

  return sin(thetaval);
}

// ------------------------------------------------------------
// Momentum contra-vector divided by energy (in comoving frame)
Real MomentumUnitEnergy(int mu, Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3, Real xi1, Real xi2, Real xi3) {
  Real result = 0.;
  if (mu == 0) {
    result = 1.;
  } else if (mu == 1) {
    result = CosPhiSinTheta(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3);
  } else if (mu == 2) {
    result = SinPhiSinTheta(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3);
  } else if (mu == 3) {
    result = CosTheta(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3);
  } else {
    // std::cout << "Incorrect choice of index for p^mu/e!" << std::endl;
    exit(EXIT_FAILURE);
  }

  return result;
}

Real MomentumUnitEnergy(int mu, Real phi, Real theta) {
  Real result = 0.;
  if (mu == 0) {
    result = 1.;
  } else if (mu == 1) {
    result = cos(phi) * sin(theta);
  } else if (mu == 2) {
    result = sin(phi) * sin(theta);
  } else if (mu == 3) {
    result = cos(theta);
  } else {
    // std::cout << "Incorrect choice of index for p^mu/e!" << std::endl;
    exit(EXIT_FAILURE);
  }

  return result;
}

// ------------------------------------------------------------------------
// partial xi1 / partial phi
Real pXi1pPhi(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3, Real xi1, Real xi2, Real xi3) {
  return (pow(x1 * xi1 + x2 * xi2 + x3 * xi3, 2) + pow(xi1 * y1 + xi2 * y2 + xi3 * y3, 2)) / (x2 * xi2 * y1 + x3 * xi3 * y1 - x1 * (xi2 * y2 + xi3 * y3));
}

// ------------------------------------------------------------------------
// partial xi2 / partial phi
Real pXi2pPhi(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3, Real xi1, Real xi2, Real xi3) {
  return (pow(x1 * xi1 + x2 * xi2 + x3 * xi3, 2) + pow(xi1 * y1 + xi2 * y2 + xi3 * y3, 2)) / ((x1 * xi1 + x3 * xi3) * y2 - x2 * (xi1 * y1 + xi3 * y3));
}

// ------------------------------------------------------------------------
// partial xi3 / partial phi
Real pXi3pPhi(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3, Real xi1, Real xi2, Real xi3) {
  return (pow(x1 * xi1 + x2 * xi2 + x3 * xi3, 2) + pow(xi1 * y1 + xi2 * y2 + xi3 * y3, 2)) / (-(x3 * (xi1 * y1 + xi2 * y2)) + (x1 * xi1 + x2 * xi2) * y3);
}

// ------------------------------------------------------------------------
// partial xi1 / partial theta
Real pXi1pTheta(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3, Real xi1, Real xi2, Real xi3) {
  return (-2 * pow(pow(x1 * xi1 + x2 * xi2 + x3 * xi3, 2) + pow(xi1 * y1 + xi2 * y2 + xi3 * y3, 2) + pow(xi1 * z1 + xi2 * z2 + xi3 * z3, 2), 1.5) *
      sqrt(1 - pow(xi1 * z1 + xi2 * z2 + xi3 * z3, 2)
          / (pow(x1 * xi1 + x2 * xi2 + x3 * xi3, 2) + pow(xi1 * y1 + xi2 * y2 + xi3 * y3, 2) + pow(xi1 * z1 + xi2 * z2 + xi3 * z3, 2)))) /
      (-2 * (xi1 * z1 + xi2 * z2 + xi3 * z3) * (x1 * (x1 * xi1 + x2 * xi2 + x3 * xi3) + y1 * (xi1 * y1 + xi2 * y2 + xi3 * y3) + z1 * (xi1 * z1 + xi2 * z2 + xi3 * z3))
          +
              2 * z1 * (pow(x1 * xi1 + x2 * xi2 + x3 * xi3, 2) + pow(xi1 * y1 + xi2 * y2 + xi3 * y3, 2) + pow(xi1 * z1 + xi2 * z2 + xi3 * z3, 2)));
}

// ------------------------------------------------------------------------
// partial xi2 / partial theta
Real pXi2pTheta(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3, Real xi1, Real xi2, Real xi3) {
  return (-2 * pow(pow(x1 * xi1 + x2 * xi2 + x3 * xi3, 2) + pow(xi1 * y1 + xi2 * y2 + xi3 * y3, 2) + pow(xi1 * z1 + xi2 * z2 + xi3 * z3, 2), 1.5) *
      sqrt(1 - pow(xi1 * z1 + xi2 * z2 + xi3 * z3, 2)
          / (pow(x1 * xi1 + x2 * xi2 + x3 * xi3, 2) + pow(xi1 * y1 + xi2 * y2 + xi3 * y3, 2) + pow(xi1 * z1 + xi2 * z2 + xi3 * z3, 2)))) /
      (-2 * (xi1 * z1 + xi2 * z2 + xi3 * z3) * (x2 * (x1 * xi1 + x2 * xi2 + x3 * xi3) + y2 * (xi1 * y1 + xi2 * y2 + xi3 * y3) + z2 * (xi1 * z1 + xi2 * z2 + xi3 * z3))
          +
              2 * z2 * (pow(x1 * xi1 + x2 * xi2 + x3 * xi3, 2) + pow(xi1 * y1 + xi2 * y2 + xi3 * y3, 2) + pow(xi1 * z1 + xi2 * z2 + xi3 * z3, 2)));
}

// ------------------------------------------------------------------------
// partial xi3 / partial theta
Real pXi3pTheta(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3, Real xi1, Real xi2, Real xi3) {
  return (-2 * pow(pow(x1 * xi1 + x2 * xi2 + x3 * xi3, 2) + pow(xi1 * y1 + xi2 * y2 + xi3 * y3, 2) + pow(xi1 * z1 + xi2 * z2 + xi3 * z3, 2), 1.5) *
      sqrt(1 - pow(xi1 * z1 + xi2 * z2 + xi3 * z3, 2)
          / (pow(x1 * xi1 + x2 * xi2 + x3 * xi3, 2) + pow(xi1 * y1 + xi2 * y2 + xi3 * y3, 2) + pow(xi1 * z1 + xi2 * z2 + xi3 * z3, 2)))) /
      (-2 * (xi1 * z1 + xi2 * z2 + xi3 * z3) * (x3 * (x1 * xi1 + x2 * xi2 + x3 * xi3) + y3 * (xi1 * y1 + xi2 * y2 + xi3 * y3) + z3 * (xi1 * z1 + xi2 * z2 + xi3 * z3))
          +
              2 * z3 * (pow(x1 * xi1 + x2 * xi2 + x3 * xi3, 2) + pow(xi1 * y1 + xi2 * y2 + xi3 * y3, 2) + pow(xi1 * z1 + xi2 * z2 + xi3 * z3, 2)));
}

inline Real dFEMBasisdphi(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3, Real xi1, Real xi2, Real xi3, int basis_index) {
  return dfem_dxi(xi1, xi2, xi3, basis_index, 1) * pXi1pPhi(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3)
      + dfem_dxi(xi1, xi2, xi3, basis_index, 2) * pXi2pPhi(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3)
      + dfem_dxi(xi1, xi2, xi3, basis_index, 3) * pXi3pPhi(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3);
}

inline Real dFEMBasisdtheta(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3, Real xi1, Real xi2, Real xi3, int basis_index) {
  return dfem_dxi(xi1, xi2, xi3, basis_index, 1) * pXi1pTheta(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3)
      + dfem_dxi(xi1, xi2, xi3, basis_index, 2) * pXi2pTheta(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3)
      + dfem_dxi(xi1, xi2, xi3, basis_index, 3) * pXi3pTheta(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3);
}

inline Real dJacxIxiJ(int i, int j, Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3) {

  Real result;
  if (i == 1 && j == 1) {
    result = x1;
  } else if (i == 1 && j == 2) {
    result = x2;
  } else if (i == 1 && j == 3) {
    result = x3;
  } else if (i == 2 && j == 1) {
    result = y1;
  } else if (i == 2 && j == 2) {
    result = y2;
  } else if (i == 2 && j == 3) {
    result = y3;
  } else if (i == 3 && j == 1) {
    result = z1;
  } else if (i == 3 && j == 2) {
    result = z2;
  } else {
    result = z3;
  }

  return result;
}

inline Real dpIdxJ(int i, int j, Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3, Real xi1, Real xi2, Real xi3) {

  Real x = xi1 * x1 + xi2 * x2 + xi3 * x3;
  Real y = xi1 * y1 + xi2 * y2 + xi3 * y3;
  Real z = xi1 * z1 + xi2 * z2 + xi3 * z3;

  Real result;
  if (i == 1 && j == 1) {
    result = (pow(y, 2) + pow(z, 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5);
  } else if ((i == 1 && j == 2) || (i == 2 && j == 1)) {
    result = -((x * y) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
  } else if ((i == 1 && j == 3) || (i == 3 && j == 1)) {
    result = -((x * z) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
  } else if (i == 2 && j == 2) {
    result = (pow(x, 2) + pow(z, 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5);
  } else if ((i == 2 && j == 3) || (i == 3 && j == 2)) {
    result = -((y * z) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5));
  } else {
    result = (pow(x, 2) + pow(y, 2)) / pow(pow(x, 2) + pow(y, 2) + pow(z, 2), 1.5);
  }

  return result;
}

inline Real dxiIdpJ(int i, int j, Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3, Real xi1, Real xi2, Real xi3) {

  Real result = dJacxIxiJ(i, 1, x1, y1, z1, x2, y2, z2, x3, y3, z3) * dpIdxJ(1, j, x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3)
      + dJacxIxiJ(i, 2, x1, y1, z1, x2, y2, z2, x3, y3, z3) * dpIdxJ(2, j, x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3)
      + dJacxIxiJ(i, 3, x1, y1, z1, x2, y2, z2, x3, y3, z3) * dpIdxJ(3, j, x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3);

  return 1 / result;
}

// Derivative of FEM basis with respect to cartesian coordinates on the unit sphere
inline Real dFEMBasisdxI(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3, Real xi1, Real xi2, Real xi3, int basis_index_a, int idx) {

  Real result = 0.;
  switch (idx) {
    case 1:
      result =
          dfem_dxi(xi1, xi2, xi3, basis_index_a, 1) * (y3 * z2 - y2 * z3)
              / (x3 * y2 * z1 - x2 * y3 * z1 - x3 * y1 * z2 + x1 * y3 * z2 + x2 * y1 * z3 - x1 * y2 * z3)
              + dfem_dxi(xi1, xi2, xi3, basis_index_a, 2) * (y3 * z1 - y1 * z3)
                  / (-(x3 * y2 * z1) + x2 * y3 * z1 + x3 * y1 * z2 - x1 * y3 * z2 - x2 * y1 * z3 + x1 * y2 * z3)
              + dfem_dxi(xi1, xi2, xi3, basis_index_a, 3) * (y2 * z1 - y1 * z2)
                  / (x3 * y2 * z1 - x2 * y3 * z1 - x3 * y1 * z2 + x1 * y3 * z2 + x2 * y1 * z3 - x1 * y2 * z3);
      break;
    case 2:
      result =
          dfem_dxi(xi1, xi2, xi3, basis_index_a, 1) * (x3 * z2 - x2 * z3)
              / (-(x3 * y2 * z1) + x2 * y3 * z1 + x3 * y1 * z2 - x1 * y3 * z2 - x2 * y1 * z3 + x1 * y2 * z3)
              + dfem_dxi(xi1, xi2, xi3, basis_index_a, 2) * (x3 * z1 - x1 * z3)
                  / (x3 * y2 * z1 - x2 * y3 * z1 - x3 * y1 * z2 + x1 * y3 * z2 + x2 * y1 * z3 - x1 * y2 * z3)
              + dfem_dxi(xi1, xi2, xi3, basis_index_a, 3) * (x2 * z1 - x1 * z2)
                  / (-(x3 * y2 * z1) + x2 * y3 * z1 + x3 * y1 * z2 - x1 * y3 * z2 - x2 * y1 * z3 + x1 * y2 * z3);
      break;
    case 3:
      result =
          dfem_dxi(xi1, xi2, xi3, basis_index_a, 1) * (x3 * y2 - x2 * y3)
              / (x3 * y2 * z1 - x2 * y3 * z1 - x3 * y1 * z2 + x1 * y3 * z2 + x2 * y1 * z3 - x1 * y2 * z3)
              + dfem_dxi(xi1, xi2, xi3, basis_index_a, 2) * (x3 * y1 - x1 * y3)
                  / (-(x3 * y2 * z1) + x2 * y3 * z1 + x3 * y1 * z2 - x1 * y3 * z2 - x2 * y1 * z3 + x1 * y2 * z3)
              + dfem_dxi(xi1, xi2, xi3, basis_index_a, 3) * (x2 * y1 - x1 * y2)
                  / (x3 * y2 * z1 - x2 * y3 * z1 - x3 * y1 * z2 + x1 * y3 * z2 + x2 * y1 * z3 - x1 * y2 * z3);
      break;
  }
  return result;
}

inline Real Rmatrix(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3, Real xi1, Real xi2, Real xi3, int i, int ihat) {

  Real xval, yval, zval;
  BarycentricToCartesian(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3, xval, yval, zval);

  Real rval = sqrt(xval * xval + yval * yval + zval * zval);
  Real thetaval = acos(zval / rval);
  Real phival = atan2(yval, xval);

  Real result = 0.;
  if (i == 1 && ihat == 1) {
    result = 1. - sin(thetaval) * sin(thetaval) * cos(thetaval) * cos(thetaval);
  } else if (i == 1 && ihat == 2) {
    result = -sin(thetaval) * sin(thetaval) * sin(phival) * cos(phival);
  } else if (i == 1 && ihat == 3) {
    result = -sin(thetaval);
  } else if (i == 2 && ihat == 1) {
    result = -sin(thetaval) * sin(thetaval) * sin(phival) * cos(phival);
  } else if (i == 2 && ihat == 2) {
    result = 1. - sin(thetaval) * sin(thetaval) * sin(phival) * sin(phival);
  } else if (i == 2 && ihat == 3) {
    result = -sin(thetaval) * cos(thetaval) * sin(phival);
  } else if (i == 3 && ihat == 1) {
    result = -sin(thetaval) * cos(thetaval) * cos(phival);
  } else if (i == 3 && ihat == 2) {
    result = -sin(thetaval) * cos(thetaval) * sin(phival);
  } else if (i == 3 && ihat == 3) {
    result = sin(thetaval) * sin(thetaval);
  } else {
    result = -42.;
  }

  return result;
}

Real dFEMBasisdp(int ihat, int a, int t1, int t2, int t3, Real x1, Real y1, Real z1, Real x2, Real y2, Real z2, Real x3, Real y3, Real z3,
                 Real xi1, Real xi2, Real xi3, int basis_choice) {

  int basis_index_a = (a == t1) * 1 + (a == t2) * 2 + (a == t3) * 3;

  Real
      result = dFEMBasisdxI(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3, basis_index_a, 1) * Rmatrix(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3, 1, ihat)
      + dFEMBasisdxI(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3, basis_index_a, 2) * Rmatrix(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3, 1, ihat)
      + dFEMBasisdxI(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3, basis_index_a, 3) * Rmatrix(x1, y1, z1, x2, y2, z2, x3, y3, z3, xi1, xi2, xi3, 1, ihat);
  return result;
}
} // namespace radiationfemn