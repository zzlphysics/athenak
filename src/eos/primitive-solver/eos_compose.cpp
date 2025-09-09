//========================================================================================
// PrimitiveSolver equation-of-state framework
// Copyright(C) 2023 Jacob M. Fields <jmf6719@psu.edu>
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file eos_compose.cpp
//  \brief Implementation of EOSCompose

#include <math.h>

#include <cassert>
#include <cstdio>
#include <limits>
#include <iostream>
#include <cstddef>
#include <string>
#include <vector>

#include "eos_compose.hpp"
#include "utils/tr_table.hpp"

using namespace Primitive; // NOLINT

void EOSCompOSE::ReadTableFromFile(std::string fname) {
  if (m_initialized==false) {
    TableReader::Table table;
    auto read_result = table.ReadTable(fname);
    if (read_result.error != TableReader::ReadResult::SUCCESS) {
      std::cout << "Table could not be read.\n";
      assert (false);
    }
    // Make sure table has correct dimentions
    assert(table.GetNDimensions()==3);
    // TODO(PH) check that required fields are present?

    // Read baryon (neutron) mass
    auto& table_scalars = table.GetScalars();
    mb = table_scalars.at("mn");

    // Get table dimesnions
    auto& point_info = table.GetPointInfo();
    m_nn = point_info[0].second;
    m_ny = point_info[1].second;
    m_nt = point_info[2].second;

    // (Re)Allocate device storage
    Kokkos::realloc(m_log_nb, m_nn);
    Kokkos::realloc(m_yq,     m_ny);
    Kokkos::realloc(m_log_t,  m_nt);
    Kokkos::realloc(m_table, ECNVARS, m_nn, m_ny, m_nt);

    // Create host storage to read into
    HostArray1D<Real>::HostMirror host_log_nb = create_mirror_view(m_log_nb);
    HostArray1D<Real>::HostMirror host_yq =     create_mirror_view(m_yq);
    HostArray1D<Real>::HostMirror host_log_t =  create_mirror_view(m_log_t);
    HostArray4D<Real>::HostMirror host_table =  create_mirror_view(m_table);

    { // read nb
      double * table_nb_double = table["nb"];
#if SINGLE_PRECISION_ENABLED
      // For single precision, convert directly during access
      for (size_t in=0; in<m_nn; ++in) {
        host_log_nb(in) = log(static_cast<Real>(table_nb_double[in]));
      }
      m_id_log_nb = Real(1.0)/(host_log_nb(1) - host_log_nb(0));
      min_n = static_cast<Real>(table_nb_double[0]);
      max_n = static_cast<Real>(table_nb_double[m_nn-1]);
#else
      // For double precision, use direct access
      Real * table_nb = reinterpret_cast<Real*>(table_nb_double);
      for (size_t in=0; in<m_nn; ++in) {
        host_log_nb(in) = log(table_nb[in]);
      }
      m_id_log_nb = 1.0/(host_log_nb(1) - host_log_nb(0));
      min_n = table_nb[0];
      max_n = table_nb[m_nn-1];
#endif
    }

    { // read yq
      double * table_yq_double = table["yq"];
#if SINGLE_PRECISION_ENABLED
      // For single precision, convert directly during access
      for (size_t iy=0; iy<m_ny; ++iy) {
        host_yq(iy) = static_cast<Real>(table_yq_double[iy]);
      }
      m_id_yq = Real(1.0)/(host_yq(1) - host_yq(0));
      min_Y[0] = static_cast<Real>(table_yq_double[0]);
      max_Y[0] = static_cast<Real>(table_yq_double[m_ny-1]);
#else
      // For double precision, use direct access
      Real * table_yq = reinterpret_cast<Real*>(table_yq_double);
      for (size_t iy=0; iy<m_ny; ++iy) {
        host_yq(iy) = table_yq[iy];
      }
      m_id_yq = 1.0/(host_yq(1) - host_yq(0));
      min_Y[0] = table_yq[0];
      max_Y[0] = table_yq[m_ny-1];
#endif
    }

    { // read T
      double * table_t_double = table["t"];
#if SINGLE_PRECISION_ENABLED
      // For single precision, convert directly during access
      for (size_t it=0; it<m_nt; ++it) {
        host_log_t(it) = log(static_cast<Real>(table_t_double[it]));
      }
      m_id_log_t = Real(1.0)/(host_log_t(1) - host_log_t(0));
      min_T = static_cast<Real>(table_t_double[1]);      // These are different
      max_T = static_cast<Real>(table_t_double[m_nt-2]); // on purpose
#else
      // For double precision, use direct access
      Real * table_t = reinterpret_cast<Real*>(table_t_double);
      for (size_t it=0; it<m_nt; ++it) {
        host_log_t(it) = log(table_t[it]);
      }
      m_id_log_t = 1.0/(host_log_t(1) - host_log_t(0));
      min_T = table_t[1];      // These are different
      max_T = table_t[m_nt-2]; // on purpose
#endif
    }

    { // Read Q1 -> log(P)
      double * table_Q1_double = table["Q1"];
#if SINGLE_PRECISION_ENABLED
      // For single precision, convert directly during access
      for (size_t in=0; in<m_nn; ++in) {
        for (size_t iy=0; iy<m_ny; ++iy) {
          for (size_t it=0; it<m_nt; ++it) {
            size_t iflat = it + m_nt*(iy + m_ny*in);
            host_table(ECLOGP,in,iy,it) = log(static_cast<Real>(table_Q1_double[iflat])) + host_log_nb(in);
          }
        }
      }
#else
      // For double precision, use direct access
      Real * table_Q1 = reinterpret_cast<Real*>(table_Q1_double);
      for (size_t in=0; in<m_nn; ++in) {
        for (size_t iy=0; iy<m_ny; ++iy) {
          for (size_t it=0; it<m_nt; ++it) {
            size_t iflat = it + m_nt*(iy + m_ny*in);
            host_table(ECLOGP,in,iy,it) = log(table_Q1[iflat]) + host_log_nb(in);
          }
        }
      }
#endif
    }

    { // Read Q2 -> S
      double * table_Q2_double = table["Q2"];
#if SINGLE_PRECISION_ENABLED
      // For single precision, convert directly during access
      for (size_t in=0; in<m_nn; ++in) {
        for (size_t iy=0; iy<m_ny; ++iy) {
          for (size_t it=0; it<m_nt; ++it) {
            size_t iflat = it + m_nt*(iy + m_ny*in);
            host_table(ECENT,in,iy,it) = static_cast<Real>(table_Q2_double[iflat]);
          }
        }
      }
#else
      // For double precision, use direct access
      Real * table_Q2 = reinterpret_cast<Real*>(table_Q2_double);
      for (size_t in=0; in<m_nn; ++in) {
        for (size_t iy=0; iy<m_ny; ++iy) {
          for (size_t it=0; it<m_nt; ++it) {
            size_t iflat = it + m_nt*(iy + m_ny*in);
            host_table(ECENT,in,iy,it) = table_Q2[iflat];
          }
        }
      }
#endif
    }

    { // Read Q3-> mu_b
      double * table_Q3_double = table["Q3"];
#if SINGLE_PRECISION_ENABLED
      // For single precision, convert directly during access
      for (size_t in=0; in<m_nn; ++in) {
        for (size_t iy=0; iy<m_ny; ++iy) {
          for (size_t it=0; it<m_nt; ++it) {
            size_t iflat = it + m_nt*(iy + m_ny*in);
            host_table(ECMUB,in,iy,it) = (static_cast<Real>(table_Q3_double[iflat])+Real(1))*mb;
          }
        }
      }
#else
      // For double precision, use direct access
      Real * table_Q3 = reinterpret_cast<Real*>(table_Q3_double);
      for (size_t in=0; in<m_nn; ++in) {
        for (size_t iy=0; iy<m_ny; ++iy) {
          for (size_t it=0; it<m_nt; ++it) {
            size_t iflat = it + m_nt*(iy + m_ny*in);
            host_table(ECMUB,in,iy,it) = (table_Q3[iflat]+1)*mb;
          }
        }
      }
#endif
    }

    { // Read Q4-> mu_q
      double * table_Q4_double = table["Q4"];
#if SINGLE_PRECISION_ENABLED
      // For single precision, convert directly during access
      for (size_t in=0; in<m_nn; ++in) {
        for (size_t iy=0; iy<m_ny; ++iy) {
          for (size_t it=0; it<m_nt; ++it) {
            size_t iflat = it + m_nt*(iy + m_ny*in);
            host_table(ECMUQ,in,iy,it) = static_cast<Real>(table_Q4_double[iflat])*mb;
          }
        }
      }
#else
      // For double precision, use direct access
      Real * table_Q4 = reinterpret_cast<Real*>(table_Q4_double);
      for (size_t in=0; in<m_nn; ++in) {
        for (size_t iy=0; iy<m_ny; ++iy) {
          for (size_t it=0; it<m_nt; ++it) {
            size_t iflat = it + m_nt*(iy + m_ny*in);
            host_table(ECMUQ,in,iy,it) = table_Q4[iflat]*mb;
          }
        }
      }
#endif
    }

    { // Read Q5-> mu_le
      double * table_Q5_double = table["Q5"];
#if SINGLE_PRECISION_ENABLED
      // For single precision, convert directly during access
      for (size_t in=0; in<m_nn; ++in) {
        for (size_t iy=0; iy<m_ny; ++iy) {
          for (size_t it=0; it<m_nt; ++it) {
            size_t iflat = it + m_nt*(iy + m_ny*in);
            host_table(ECMUL,in,iy,it) = static_cast<Real>(table_Q5_double[iflat])*mb;
          }
        }
      }
#else
      // For double precision, use direct access
      Real * table_Q5 = reinterpret_cast<Real*>(table_Q5_double);
      for (size_t in=0; in<m_nn; ++in) {
        for (size_t iy=0; iy<m_ny; ++iy) {
          for (size_t it=0; it<m_nt; ++it) {
            size_t iflat = it + m_nt*(iy + m_ny*in);
            host_table(ECMUL,in,iy,it) = table_Q5[iflat]*mb;
          }
        }
      }
#endif
    }

    { // Read Q7-> log(e)
      double * table_Q7_double = table["Q7"];
#if SINGLE_PRECISION_ENABLED
      // For single precision, convert directly during access
      for (size_t in=0; in<m_nn; ++in) {
        for (size_t iy=0; iy<m_ny; ++iy) {
          for (size_t it=0; it<m_nt; ++it) {
            size_t iflat = it + m_nt*(iy + m_ny*in);
            host_table(ECLOGE,in,iy,it) = log(mb*(static_cast<Real>(table_Q7_double[iflat]) + Real(1))) + host_log_nb(in);
          }
        }
      }
#else
      // For double precision, use direct access
      Real * table_Q7 = reinterpret_cast<Real*>(table_Q7_double);
      for (size_t in=0; in<m_nn; ++in) {
        for (size_t iy=0; iy<m_ny; ++iy) {
          for (size_t it=0; it<m_nt; ++it) {
            size_t iflat = it + m_nt*(iy + m_ny*in);
            host_table(ECLOGE,in,iy,it) = log(mb*(table_Q7[iflat] + 1)) + host_log_nb(in);
          }
        }
      }
#endif
    }

    { // Read cs2-> cs
      double * table_cs2_double = table["cs2"];
#if SINGLE_PRECISION_ENABLED
      // For single precision, convert directly during access
      for (size_t in=0; in<m_nn; ++in) {
        for (size_t iy=0; iy<m_ny; ++iy) {
          for (size_t it=0; it<m_nt; ++it) {
            size_t iflat = it + m_nt*(iy + m_ny*in);
            host_table(ECCS,in,iy,it) = sqrt(static_cast<Real>(table_cs2_double[iflat]));
          }
        }
      }
#else
      // For double precision, use direct access
      Real * table_cs2 = reinterpret_cast<Real*>(table_cs2_double);
      for (size_t in=0; in<m_nn; ++in) {
        for (size_t iy=0; iy<m_ny; ++iy) {
          for (size_t it=0; it<m_nt; ++it) {
            size_t iflat = it + m_nt*(iy + m_ny*in);
            host_table(ECCS,in,iy,it) = sqrt(table_cs2[iflat]);
          }
        }
      }
#endif
    }

    // Copy from host to device
    Kokkos::deep_copy(m_log_nb, host_log_nb);
    Kokkos::deep_copy(m_yq,     host_yq);
    Kokkos::deep_copy(m_log_t,  host_log_t);
    Kokkos::deep_copy(m_table,  host_table);

    m_initialized = true;

    m_min_h = std::numeric_limits<Real>::max();
    // Compute minimum enthalpy
    for (int in = 0; in < m_nn; ++in) {
      Real const nb = exp(host_log_nb(in));
      for (int it = 0; it < m_nt; ++it) {
        for (int iy = 0; iy < m_ny; ++iy) {
          // This would use GPU memory, and we are currently on the CPU, so Enthalpy is
          // hardcoded
          Real e = exp(host_table(ECLOGE,in,iy,it));
          Real p = exp(host_table(ECLOGP,in,iy,it));
          Real h = (e + p) / nb;
          m_min_h = fmin(m_min_h, h);
        }
      }
    }
  } // if (m_initialized==false)
}
