//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file eventlog.cpp
//! \brief writes diagnostic data collected by various event counters implemented
//! throughout the code to a log file.  Checks whether there is data to be written
//! every time step, but only writes data if one or more counters are non-zero

#include <cstdio>
#include <cinttypes>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "athena.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "outputs.hpp"

//----------------------------------------------------------------------------------------
// ctor: also calls BaseTypeOutput base class constructor

EventLogOutput::EventLogOutput(ParameterInput *pin, Mesh *pm, OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  header_written = false;
  write_zeros = pin->GetOrAddBoolean(op.block_name, "write_zeros", false);
}

//----------------------------------------------------------------------------------------
//! \fn void EventLogOutput::LoadOutputData()
//! \brief sums event counter data across MPI ranks

void EventLogOutput::LoadOutputData(Mesh *pm) {
  // FOFC records on the device without synchronizing every subcycled RK stage.  Drain
  // once at an output boundary before combining this rank's interval with other ranks.
  pm->DrainDeviceEventCounters();
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
#endif

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
    }

    // write event counters
    if (!(no_output) || write_zeros) {
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

  UpdateOutputParameters(pm, pin, false);
  return;
}
