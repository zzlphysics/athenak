//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file eventlog.cpp
//! \brief writes diagnostic data collected by various event counters implemented
//! throughout the code to a log file.  Checks whether there is data to be written
//! every time step, but only writes data if one or more counters are non-zero

#include <algorithm>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "mhd/fofc_telemetry.hpp"
#include "mhd/mhd.hpp"
#include "outputs.hpp"

//----------------------------------------------------------------------------------------
// ctor: also calls BaseTypeOutput base class constructor

EventLogOutput::EventLogOutput(ParameterInput *pin, Mesh *pm, OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  header_written = false;
  write_zeros = pin->GetOrAddBoolean(op.block_name, "write_zeros", false);
  if (pm->pmb_pack->pmhd != nullptr) {
    fofc_spatial_telemetry = pm->pmb_pack->pmhd->fofc_spatial_telemetry;
  }
  if (fofc_spatial_telemetry) {
    // A fixed one-root-cycle interval is part of the telemetry contract.  In particular,
    // a time-based cadence could silently combine a variable number of root steps.
    if (op.dcycle != 1) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "FOFC spatial telemetry requires its event log to use "
                << "dcycle=1" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // A dense row is required even on a root step with no events.  Otherwise an absent
    // summary would be ambiguous between "zero corrections" and an interrupted output.
    if (!write_zeros) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "FOFC spatial telemetry requires its event log to use "
                << "write_zeros=true" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // This is the only nfofc contribution allowed to lack process-local spatial bins:
    // a prefix restored before Outputs construction (or produced by fresh IC setup).
    // Once the first dense event row consumes it, any later mismatch is a telemetry bug.
    fofc_telemetry_allowed_prefix = pm->ecounter.nfofc;
    fofc_telemetry_bins.assign(mhd::fofc_telemetry::kHistogramSize, 0);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void EventLogOutput::LoadOutputData()
//! \brief sums event counter data across MPI ranks

void EventLogOutput::LoadOutputData(Mesh *pm) {
  // FOFC records on the device without synchronizing every subcycled RK stage.  Drain
  // once at an output boundary before combining this rank's interval with other ranks.
  pm->DrainDeviceEventCounters();
  if (fofc_spatial_telemetry) {
    auto pending = pm->pmb_pack->pmhd->fofc_telemetry_pending;
    auto host_pending = Kokkos::create_mirror_view_and_copy(HostMemSpace(), pending);
    for (std::size_t n=0; n<fofc_telemetry_bins.size(); ++n) {
      fofc_telemetry_bins[n] = host_pending(n);
    }
    // Event output is the ownership boundary: after this copy, a future root step starts
    // a fresh interval even when AMR changed the MeshBlock topology in between.
    Kokkos::deep_copy(pending, std::uint64_t{0});
  }
#if MPI_PARALLEL_ENABLED
  // perform in-place sum or max over all MPI ranks, depending on counter
  std::uint64_t* sum_counters[] = {
    &(pm->ecounter.neos_dfloor), &(pm->ecounter.neos_efloor),
    &(pm->ecounter.neos_tfloor), &(pm->ecounter.neos_vceil),
    &(pm->ecounter.neos_fail), &(pm->ecounter.nfofc),
    &(pm->ecounter.ncons_adjust), &(pm->ecounter.nmag_adjust),
    &(pm->ecounter.nc2p_calls), &(pm->ecounter.nfofc_tests)
  };
  for (auto *counter : sum_counters) {
    MPI_Allreduce(MPI_IN_PLACE, counter, 1, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  }
  MPI_Allreduce(MPI_IN_PLACE, &(pm->ecounter.maxit_c2p), 1, MPI_INT, MPI_MAX,
                MPI_COMM_WORLD);
  if (fofc_spatial_telemetry) {
    MPI_Allreduce(MPI_IN_PLACE, fofc_telemetry_bins.data(),
                  static_cast<int>(fofc_telemetry_bins.size()), MPI_UINT64_T, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &fofc_telemetry_allowed_prefix, 1, MPI_UINT64_T,
                  MPI_SUM, MPI_COMM_WORLD);
  }
#endif

  if (fofc_spatial_telemetry) {
    fofc_telemetry_total = 0;
    for (const std::uint64_t count : fofc_telemetry_bins) {
      if (count > std::numeric_limits<std::uint64_t>::max()-fofc_telemetry_total) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "FOFC spatial telemetry total overflow" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      fofc_telemetry_total += count;
    }
    if (fofc_telemetry_total > pm->ecounter.nfofc) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "FOFC spatial histogram exceeds the authoritative "
                << "nfofc counter" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // Pending event counters are restart-persistent, while this diagnostic histogram is
    // intentionally not.  Attribute any carried prefix to an explicit overflow/unknown
    // bin so the histogram remains exactly conservative instead of pretending to know
    // its old level, stage, cause, or location.
    fofc_telemetry_unattributed = pm->ecounter.nfofc-fofc_telemetry_total;
    if (fofc_telemetry_unattributed != fofc_telemetry_allowed_prefix) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "FOFC spatial histogram mismatch is not the captured "
                << "pre-output prefix: missing=" << fofc_telemetry_unattributed
                << " allowed=" << fofc_telemetry_allowed_prefix << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (fofc_telemetry_unattributed > 0) {
      const std::size_t unknown_index = mhd::fofc_telemetry::HistogramIndex(
          mhd::fofc_telemetry::kLevelOverflow, mhd::fofc_telemetry::kStageOther,
          static_cast<int>(mhd::fofc_telemetry::Reason::unknown),
          mhd::fofc_telemetry::kRadiusBins-1, mhd::fofc_telemetry::kAbsZBins-1,
          mhd::fofc_telemetry::kLapseBins-1);
      if (fofc_telemetry_unattributed >
          std::numeric_limits<std::uint64_t>::max()-fofc_telemetry_bins[unknown_index]) {
        std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                  << std::endl << "FOFC unattributed telemetry overflow" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      fofc_telemetry_bins[unknown_index] += fofc_telemetry_unattributed;
      fofc_telemetry_total += fofc_telemetry_unattributed;
    }
    if (fofc_telemetry_total != pm->ecounter.nfofc) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "FOFC spatial histogram does not match nfofc"
                << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // check if there is any data to be written
  no_output = true;
  if (pm->ecounter.neos_dfloor > 0 ||
      pm->ecounter.neos_efloor > 0 ||
      pm->ecounter.neos_tfloor > 0 ||
      pm->ecounter.neos_vceil  > 0 ||
      pm->ecounter.neos_fail   > 0 ||
      pm->ecounter.nfofc > 0 ||
      pm->ecounter.ncons_adjust > 0 ||
      pm->ecounter.nmag_adjust > 0 ||
      pm->ecounter.nc2p_calls > 0 ||
      pm->ecounter.nfofc_tests > 0 ||
      pm->ecounter.maxit_c2p > 0) {
    no_output=false;
  }
}

//----------------------------------------------------------------------------------------
//! \fn void EventLogOutput::WriteOutputFile()
//! \brief writes event counter data to log file

void EventLogOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  if (header_written && no_output && !write_zeros) {
    // A sparse log still consumed this cadence point.  Persist it so Finalize and a
    // zero-step restart cannot sample the same interval a second time.
    UpdateOutputParameters(pm, pin, false);
    return;
  }

  // only the master rank writes the file
  if (global_variable::my_rank == 0) {
    // create filename: "file_basename" + ".log"
    // There is no file number or id in event log output filenames.
    std::string fname;
    fname.assign(out_params.file_basename);
    fname.append(".log");

    // open file for output
    FILE *pfile;
    if ((pfile = std::fopen(fname.c_str(),"a")) == nullptr) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "Output file '" << fname << "' could not be opened" <<std::endl;
      exit(EXIT_FAILURE);
    }

    // Write header, if it has not been written already
    if (!(header_written)) {
      std::fprintf(pfile,"# Athena event counter data\n");
      if (pm->pmb_pack->pdyngr != nullptr) {
        std::fprintf(pfile,
            "# DynGRMHD event totals count physical active zones only; ghost-zone "
            "cache work is excluded.\n");
      }
      std::fprintf(pfile,"#  cycle eos_dfloor eos_efloor eos_tfloor eos_vceil");
      std::fprintf(pfile," eos_fail c2p_it fofc cons_adjust mag_adjust");
      std::fprintf(pfile," c2p_calls fofc_tests");
      std::fprintf(pfile,"\n");  // terminate line
      if (fofc_spatial_telemetry) {
        std::fprintf(pfile,
            "# fofc_spatial_v1 kind=schema level_bins=0..31,overflow "
            "stage_bins=other,1,2,3 "
            "reason_bins=unknown,dmp_preflag,scalar,cons_density_floor,"
            "cons_energy_floor,prim_density_floor,prim_temperature_floor,"
            "rho_too_big,rho_too_small,nans_in_cons,mag_too_big,"
            "bracketing_failed,no_solution,invalid_geometry,other_c2p "
            "r_cyl_edges=2,4,8,16,32,64 abs_z_edges=0.5,1,2,4,8,16 "
            "lapse_edges=0.2,0.4,0.6,0.8,1 "
            "center1=%.17g center2=%.17g center3=%.17g\n",
            pm->pmb_pack->pmhd->fofc_telemetry_center[0],
            pm->pmb_pack->pmhd->fofc_telemetry_center[1],
            pm->pmb_pack->pmhd->fofc_telemetry_center[2]);
      }
    }

    // write event counters
    if (!(no_output) || write_zeros) {
      if (fofc_spatial_telemetry) {
        std::fprintf(pfile,
            "# fofc_spatial_v1 kind=summary cycle=%d count=%" PRIu64
            " nfofc=%" PRIu64 " unattributed=%" PRIu64 "\n",
            pm->ncycle, fofc_telemetry_total, pm->ecounter.nfofc,
            fofc_telemetry_unattributed);
        for (int level=0; level<mhd::fofc_telemetry::kLevelBins; ++level) {
          for (int stage=0; stage<mhd::fofc_telemetry::kStageBins; ++stage) {
            for (int reason=0; reason<mhd::fofc_telemetry::kReasonBins; ++reason) {
              for (int radius=0; radius<mhd::fofc_telemetry::kRadiusBins; ++radius) {
                for (int abs_z=0; abs_z<mhd::fofc_telemetry::kAbsZBins; ++abs_z) {
                  for (int lapse=0; lapse<mhd::fofc_telemetry::kLapseBins; ++lapse) {
                    const std::size_t index = mhd::fofc_telemetry::HistogramIndex(
                        level, stage, reason, radius, abs_z, lapse);
                    const std::uint64_t count = fofc_telemetry_bins[index];
                    if (count == 0) continue;
                    std::fprintf(pfile,
                        "# fofc_spatial_v1 kind=bin cycle=%d level_bin=%d "
                        "stage_bin=%d reason=%s r_cyl_bin=%d abs_z_bin=%d "
                        "lapse_bin=%d count=%" PRIu64 "\n",
                        pm->ncycle, level, stage,
                        mhd::fofc_telemetry::ReasonName(reason), radius, abs_z, lapse,
                        count);
                  }
                }
              }
            }
          }
        }
      }
      std::fprintf(pfile, "%8d", pm->ncycle);
      std::fprintf(pfile, " %8" PRIu64, pm->ecounter.neos_dfloor);
      std::fprintf(pfile, " %8" PRIu64, pm->ecounter.neos_efloor);
      std::fprintf(pfile, " %8" PRIu64, pm->ecounter.neos_tfloor);
      std::fprintf(pfile, " %8" PRIu64, pm->ecounter.neos_vceil);
      std::fprintf(pfile, " %8" PRIu64, pm->ecounter.neos_fail);
      std::fprintf(pfile, " %6d", pm->ecounter.maxit_c2p);
      std::fprintf(pfile, " %8" PRIu64, pm->ecounter.nfofc);
      std::fprintf(pfile, " %8" PRIu64, pm->ecounter.ncons_adjust);
      std::fprintf(pfile, " %8" PRIu64, pm->ecounter.nmag_adjust);
      std::fprintf(pfile, " %8" PRIu64, pm->ecounter.nc2p_calls);
      std::fprintf(pfile, " %8" PRIu64, pm->ecounter.nfofc_tests);
      std::fprintf(pfile,"\n"); // terminate line
    }
    std::fclose(pfile);
  }
  // Keep control flow and persisted output-cadence state identical on every rank.  Only
  // rank zero owns the file, but all ranks participated in this first header attempt.
  header_written = true;

  // reset counters
  pm->ecounter.neos_dfloor = 0;
  pm->ecounter.neos_efloor = 0;
  pm->ecounter.neos_tfloor = 0;
  pm->ecounter.neos_vceil = 0;
  pm->ecounter.neos_fail = 0;
  pm->ecounter.maxit_c2p = 0;
  pm->ecounter.nfofc = 0;
  pm->ecounter.ncons_adjust = 0;
  pm->ecounter.nmag_adjust = 0;
  pm->ecounter.nc2p_calls = 0;
  pm->ecounter.nfofc_tests = 0;
  if (fofc_spatial_telemetry) {
    std::fill(fofc_telemetry_bins.begin(), fofc_telemetry_bins.end(), 0);
    fofc_telemetry_total = 0;
    fofc_telemetry_unattributed = 0;
    fofc_telemetry_allowed_prefix = 0;
  }

  UpdateOutputParameters(pm, pin, false);
  return;
}
