//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file restart.cpp
//! \brief writes restart files

#include <sys/stat.h>  // mkdir

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>      // fwrite(), fclose(), fopen(), fnprintf(), snprintf()
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <utility> // make_pair

#include "athena.hpp"
#include "coordinates/cell_locations.hpp"
#include "geodesic-grid/geodesic_grid.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"
#include "coordinates/adm.hpp"
#include "z4c/compact_object_tracker.hpp"
#include "z4c/z4c.hpp"
#include "radiation/radiation.hpp"
#include "srcterms/turb_driver.hpp"
//#include "outputs.hpp"

namespace {

constexpr std::uint64_t kAmrCycleCounterMagic = UINT64_C(0x41544b414d524331);
constexpr int kAmrCycleCounterVersion = 1;
constexpr std::uint64_t kEventCounterMagic = UINT64_C(0x41544b4556543031);
constexpr int kEventCounterVersion = 1;
constexpr int kEventSumCounterCount = 10;
constexpr std::uint64_t kPerRankPartitionMagic = UINT64_C(0x41544b5052543031);
constexpr int kPerRankPartitionVersion = 2;
constexpr int kCheckpointIdentityWords = 2;

std::uint64_t MixCheckpointIdentity(std::uint64_t value) {
  value ^= value >> 30;
  value *= UINT64_C(0xbf58476d1ce4e5b9);
  value ^= value >> 27;
  value *= UINT64_C(0x94d049bb133111eb);
  value ^= value >> 31;
  return value;
}

void GenerateCheckpointIdentity(int cycle, int file_number,
                                std::uint64_t identity[kCheckpointIdentityWords]) {
  static std::uint64_t sequence = 0;
  std::random_device entropy;
  std::uint64_t seed_low =
      (static_cast<std::uint64_t>(entropy()) << 32) ^
      static_cast<std::uint64_t>(entropy());
  std::uint64_t seed_high =
      (static_cast<std::uint64_t>(entropy()) << 32) ^
      static_cast<std::uint64_t>(entropy());
  const std::uint64_t high_clock = static_cast<std::uint64_t>(
      std::chrono::high_resolution_clock::now().time_since_epoch().count());
  const std::uint64_t steady_clock = static_cast<std::uint64_t>(
      std::chrono::steady_clock::now().time_since_epoch().count());
  const std::uint64_t current_sequence = ++sequence;
  seed_low ^= high_clock;
  seed_low ^= (static_cast<std::uint64_t>(cycle) << 32);
  seed_low ^= current_sequence * UINT64_C(0x9e3779b97f4a7c15);
  seed_high ^= steady_clock;
  seed_high ^= (static_cast<std::uint64_t>(file_number) << 32);
  seed_high ^= current_sequence * UINT64_C(0xd1b54a32d192ed03);
  identity[0] = MixCheckpointIdentity(seed_low);
  identity[1] = MixCheckpointIdentity(seed_high);
  if ((identity[0] | identity[1]) == 0) {
    identity[1] = UINT64_C(0x9e3779b97f4a7c15);
  }
}

}  // namespace

//----------------------------------------------------------------------------------------
// constructor: also calls BaseTypeOutput base class constructor

RestartOutput::RestartOutput(ParameterInput *pin, Mesh *pm, OutputParameters op) :
  BaseTypeOutput(pin, pm, op) {
  // create directories for outputs. Comments in binary.cpp constructor explain why
  mkdir("rst",0775);
  bool single_file_per_rank = op.single_file_per_rank;
  if (single_file_per_rank) {
    char rank_dir[20];
    std::snprintf(rank_dir, sizeof(rank_dir), "rst/rank_%08d/", global_variable::my_rank);
    mkdir(rank_dir, 0775);
  }
}

//----------------------------------------------------------------------------------------
// RestartOutput::LoadOutputData()
// overload of standard load data function specific to restarts.  Loads dependent
// variables, including ghost zones.

void RestartOutput::LoadOutputData(Mesh *pm) {
  // get spatial dimensions of arrays, including ghost zones
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int nout1 = indcs.nx1 + 2*(indcs.ng);
  int nout2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int nout3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  int nmb = pm->pmb_pack->nmb_thispack;

  // calculate total number of CC variables
  hydro::Hydro* phydro = pm->pmb_pack->phydro;
  mhd::MHD* pmhd = pm->pmb_pack->pmhd;
  adm::ADM* padm = pm->pmb_pack->padm;
  z4c::Z4c* pz4c = pm->pmb_pack->pz4c;
  radiation::Radiation* prad = pm->pmb_pack->prad;
  TurbulenceDriver* pturb=pm->pmb_pack->pturb;
  int nhydro=0, nmhd=0, nrad=0, nforce=3, nadm=0, nz4c=0;
  if (phydro != nullptr) {
    nhydro = phydro->nhydro + phydro->nscalars;
  }
  if (pmhd != nullptr) {
    nmhd = pmhd->nmhd + pmhd->nscalars;
  }
  if (pz4c != nullptr) {
    nz4c = pz4c->nz4c;
  } else if (padm != nullptr) {
    nadm = padm->nadm;
  }
  // if the spacetime is evolved, we do not need to checkpoint/recover the ADM variables
  if (prad != nullptr) {
    nrad = prad->prgeo->nangles;
  }

  // Note for restarts, outarrays are dimensioned (m,n,k,j,i)
  if (phydro != nullptr) {
    Kokkos::realloc(outarray_hyd, nmb, nhydro, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_hyd, Kokkos::subview(phydro->u0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (pmhd != nullptr) {
    Kokkos::realloc(outarray_mhd, nmb, nmhd, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_mhd, Kokkos::subview(pmhd->u0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    Kokkos::realloc(outfield.x1f, nmb, nout3, nout2, nout1+1);
    Kokkos::deep_copy(outfield.x1f, Kokkos::subview(pmhd->b0.x1f, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    Kokkos::realloc(outfield.x2f, nmb, nout3, nout2+1, nout1);
    Kokkos::deep_copy(outfield.x2f, Kokkos::subview(pmhd->b0.x2f, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
    Kokkos::realloc(outfield.x3f, nmb, nout3+1, nout2, nout1);
    Kokkos::deep_copy(outfield.x3f, Kokkos::subview(pmhd->b0.x3f, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (prad != nullptr) {
    Kokkos::realloc(outarray_rad, nmb, nrad, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_rad, Kokkos::subview(prad->i0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (pturb != nullptr) {
    Kokkos::realloc(outarray_force, nmb, nforce, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_force, Kokkos::subview(pturb->force, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }
  if (pz4c != nullptr) {
    Kokkos::realloc(outarray_z4c, nmb, nz4c, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_z4c, Kokkos::subview(pz4c->u0, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  } else if (padm != nullptr) {
    Kokkos::realloc(outarray_adm, nmb, nadm, nout3, nout2, nout1);
    Kokkos::deep_copy(outarray_adm, Kokkos::subview(padm->u_adm, std::make_pair(0,nmb),
                      Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL));
  }

  // calculate max/min number of MeshBlocks across all ranks
  noutmbs_max = pm->nmb_eachrank[0];
  noutmbs_min = pm->nmb_eachrank[0];
  for (int i=0; i<(global_variable::nranks); ++i) {
    noutmbs_max = std::max(noutmbs_max,pm->nmb_eachrank[i]);
    noutmbs_min = std::min(noutmbs_min,pm->nmb_eachrank[i]);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void RestartOutput:::WriteOutputFile(Mesh *pm)
//  \brief Cycles over all MeshBlocks and writes everything to a single restart file

void RestartOutput::WriteOutputFile(Mesh *pm, ParameterInput *pin) {
  // get spatial dimensions of arrays, including ghost zones
  auto &indcs = pm->pmb_pack->pmesh->mb_indcs;
  int nout1 = indcs.nx1 + 2*(indcs.ng);
  int nout2 = (indcs.nx2 > 1)? (indcs.nx2 + 2*(indcs.ng)) : 1;
  int nout3 = (indcs.nx3 > 1)? (indcs.nx3 + 2*(indcs.ng)) : 1;
  hydro::Hydro* phydro = pm->pmb_pack->phydro;
  mhd::MHD* pmhd = pm->pmb_pack->pmhd;
  radiation::Radiation* prad = pm->pmb_pack->prad;
  TurbulenceDriver* pturb=pm->pmb_pack->pturb;
  z4c::Z4c* pz4c = pm->pmb_pack->pz4c;
  adm::ADM* padm = pm->pmb_pack->padm;
  int nhydro=0, nmhd=0, nrad=0, nforce=3, nz4c=0, nadm=0, nco=0;
  if (phydro != nullptr) {
    nhydro = phydro->nhydro + phydro->nscalars;
  }
  if (pmhd != nullptr) {
    nmhd = pmhd->nmhd + pmhd->nscalars;
  }
  if (prad != nullptr) {
    nrad = prad->prgeo->nangles;
  }
  if (pz4c != nullptr) {
    nz4c = pz4c->nz4c;
    nco = pz4c->ptracker.size();
  } else if (padm != nullptr) {
    nadm = padm->nadm;
  }
  bool single_file_per_rank = out_params.single_file_per_rank;
  // FOFC counting is deferred on the device to avoid a fence at every subcycled RK
  // stage.  A checkpoint is a diagnostic boundary too: move pending increments into
  // the host interval before serializing it.  Draining clears only the device scalar;
  // the host accumulator remains live until an event-log output consumes it.
  pm->DrainDeviceEventCounters();
  // A shared restart must contain one global diagnostic accumulator.  Restore that
  // value on rank zero only (see BuildTreeFromRestart), so the next event-log
  // MPI_Allreduce adds it exactly once.  Per-rank restarts retain local accumulators;
  // their writer/current rank layout is already required to match exactly.
  const std::uint64_t local_event_sums[kEventSumCounterCount] = {
    pm->ecounter.neos_dfloor, pm->ecounter.neos_efloor,
    pm->ecounter.neos_tfloor, pm->ecounter.neos_vceil,
    pm->ecounter.neos_fail, pm->ecounter.nfofc,
    pm->ecounter.ncons_adjust, pm->ecounter.nmag_adjust,
    pm->ecounter.nc2p_calls, pm->ecounter.nfofc_tests
  };
  std::uint64_t restart_event_sums[kEventSumCounterCount] = {};
  int restart_event_maxit = 0;
  if (single_file_per_rank) {
    std::copy(local_event_sums, local_event_sums + kEventSumCounterCount,
              restart_event_sums);
    restart_event_maxit = pm->ecounter.maxit_c2p;
  } else {
#if MPI_PARALLEL_ENABLED
    MPI_Reduce(local_event_sums, restart_event_sums, kEventSumCounterCount,
               MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(pm->ecounter.maxit_c2p), &restart_event_maxit, 1, MPI_INT,
               MPI_MAX, 0, MPI_COMM_WORLD);
#else
    std::copy(local_event_sums, local_event_sums + kEventSumCounterCount,
              restart_event_sums);
    restart_event_maxit = pm->ecounter.maxit_c2p;
#endif
  }
  std::uint64_t checkpoint_identity[kCheckpointIdentityWords] = {0, 0};
  if (single_file_per_rank) {
    if (global_variable::my_rank == 0) {
      GenerateCheckpointIdentity(pm->ncycle, out_params.file_number,
                                 checkpoint_identity);
    }
#if MPI_PARALLEL_ENABLED
    MPI_Bcast(checkpoint_identity, sizeof(checkpoint_identity), MPI_BYTE, 0,
              MPI_COMM_WORLD);
#endif
  }
  std::string fname;
  if (single_file_per_rank) {
    // Generate a directory and filename for each rank
    // create filename: "rst/rank_YYYYYYY/file_basename" + "." + XXXXX + ".rst"
    // where YYYYYYY = 8-digit rank number
    // where XXXXX = 5-digit file_number
    char rank_dir[20];
    char number[7];
    std::snprintf(number, sizeof(number), ".%05d", out_params.file_number);
    std::snprintf(rank_dir, sizeof(rank_dir), "rank_%08d/", global_variable::my_rank);
    fname = std::string("rst/") + std::string(rank_dir) + out_params.file_basename
      + number + ".rst";

    // Debugging output to check directory and filename
    // std::cout << "Rank " << global_variable::my_rank << " generated filename: "
    //           << fname << std::endl;
  } else {
    // Existing behavior: single restart file
    // create filename: "rst/file_basename" + "." + XXXXX + ".rst"
    // where XXXXX = 5-digit file_number
    char number[7];
    std::snprintf(number, sizeof(number), ".%05d", out_params.file_number);
    fname = std::string("rst/") + out_params.file_basename + number + ".rst";
  }
  // Update before serializing ParameterInput so the checkpoint contains the next unique
  // file number and either the advanced or retained cadence selected by the Driver.
  UpdateOutputParameters(pm, pin, true);

  // Advertise the optional AMR metadata extension in the serialized ParameterInput.
  // Old restart files lack this key and retain their original byte layout.
  if (pm->adaptive) {
    pin->SetInteger("mesh_refinement", "restart_cycle_counters_version",
                    kAmrCycleCounterVersion);
  }

  // The binary header retains the actual last completed step for time-derived output.
  // Store the independent growth-limiter reference in ParameterInput so a timestep
  // clipped only to hit tlim does not force a roundoff-scale restart staircase.
  const Real restart_growth =
      (std::isfinite(pm->dt_restart_growth) && pm->dt_restart_growth > 0.0)
      ? pm->dt_restart_growth : std::numeric_limits<float>::max();
  pin->SetReal("time", "restart_dt_growth", restart_growth);

  // create string holding input parameters (copy of input file)
  std::stringstream ost;
  pin->ParameterDump(ost);
  std::string sbuf = ost.str();

  //--- STEP 1.  Root process writes header data (input file, critical variables)
  // Input file data is read by ParameterInput on restart, and the remaining header
  // variables are read in Mesh::BuildTreeFromRestart()

  // open file and  write the header; this part is serial
  IOWrapper resfile;
  resfile.Open(fname.c_str(), IOWrapper::FileMode::write, single_file_per_rank);
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    // The legacy header has one timestep slot.  Keep its unambiguous physical meaning:
    // the actual last completed step used by time-derived diagnostics.  The separate
    // restart_dt_growth ParameterInput value controls the post-restart growth limiter.
    // Cycle-zero checkpoints use the fresh-Mesh sentinel in this legacy slot.
    const Real restart_dt =
        (pm->ncycle > 0 && std::isfinite(pm->dt_last_completed) &&
         pm->dt_last_completed > 0.0)
        ? pm->dt_last_completed : std::numeric_limits<float>::max();

    // output the input parameters (input file)
    resfile.Write_any_type(sbuf.c_str(), sbuf.size(), "byte", single_file_per_rank);

    // output Mesh information
    resfile.Write_any_type(&(pm->nmb_total), (sizeof(int)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->root_level), (sizeof(int)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->mesh_size), (sizeof(RegionSize)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->mesh_indcs), (sizeof(RegionIndcs)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->mb_indcs), (sizeof(RegionIndcs)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&(pm->time), (sizeof(Real)), "byte",
                            single_file_per_rank);
    resfile.Write_any_type(&restart_dt, sizeof(Real), "byte", single_file_per_rank);
    resfile.Write_any_type(&(pm->ncycle), (sizeof(int)), "byte",
                            single_file_per_rank);
  }
  //--- STEP 2.  Root process writes list of logical locations and cost of MeshBlocks
  // This data read in Mesh::BuildTreeFromRestart()

  if (global_variable::my_rank == 0 || single_file_per_rank) {
    resfile.Write_any_type(&(pm->lloc_eachmb[0]),(pm->nmb_total)*sizeof(LogicalLocation),
                           "byte", single_file_per_rank);
    resfile.Write_any_type(&(pm->cost_eachmb[0]), (pm->nmb_total)*sizeof(float),
                           "byte", single_file_per_rank);
    // A per-rank file stores only the field data owned by its writer.  Persist that
    // exact partition so restart can reject a different MPI layout instead of silently
    // attaching the local data stream to the wrong global MeshBlock IDs.
    if (single_file_per_rank) {
      const int writer_nranks = global_variable::nranks;
      const int writer_rank = global_variable::my_rank;
      const int writer_gid_start = pm->gids_eachrank[writer_rank];
      const int writer_nmb = pm->nmb_eachrank[writer_rank];
      resfile.Write_any_type(&kPerRankPartitionMagic, sizeof(kPerRankPartitionMagic),
                             "byte", single_file_per_rank);
      resfile.Write_any_type(&kPerRankPartitionVersion,
                             sizeof(kPerRankPartitionVersion), "byte",
                             single_file_per_rank);
      resfile.Write_any_type(&writer_nranks, sizeof(writer_nranks), "byte",
                             single_file_per_rank);
      resfile.Write_any_type(&writer_rank, sizeof(writer_rank), "byte",
                             single_file_per_rank);
      resfile.Write_any_type(&writer_gid_start, sizeof(writer_gid_start), "byte",
                             single_file_per_rank);
      resfile.Write_any_type(&writer_nmb, sizeof(writer_nmb), "byte",
                             single_file_per_rank);
      resfile.Write_any_type(checkpoint_identity, sizeof(checkpoint_identity), "byte",
                             single_file_per_rank);
    }
    if (pm->adaptive) {
      const int counter_count = pm->nmb_total;
      resfile.Write_any_type(&kAmrCycleCounterMagic, sizeof(kAmrCycleCounterMagic),
                             "byte", single_file_per_rank);
      resfile.Write_any_type(&kAmrCycleCounterVersion, sizeof(kAmrCycleCounterVersion),
                             "byte", single_file_per_rank);
      resfile.Write_any_type(&counter_count, sizeof(counter_count), "byte",
                             single_file_per_rank);
      resfile.Write_any_type(pm->pmr->ncyc_since_ref.data(),
                             counter_count*sizeof(int), "byte", single_file_per_rank);
    }
    // Persist event counters even when every value is zero.  Presence is selected by
    // the serialized marker, so new readers remain compatible with legacy restarts and
    // parameter overrides cannot accidentally change the binary layout.
    resfile.Write_any_type(&kEventCounterMagic, sizeof(kEventCounterMagic),
                           "byte", single_file_per_rank);
    resfile.Write_any_type(&kEventCounterVersion, sizeof(kEventCounterVersion),
                           "byte", single_file_per_rank);
    resfile.Write_any_type(&kEventSumCounterCount, sizeof(kEventSumCounterCount),
                           "byte", single_file_per_rank);
    resfile.Write_any_type(restart_event_sums, sizeof(restart_event_sums),
                           "byte", single_file_per_rank);
    resfile.Write_any_type(&restart_event_maxit, sizeof(restart_event_maxit),
                           "byte", single_file_per_rank);
  }

  //--- STEP 3.  Root process writes internal state of objects that require it
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    // store z4c information
    if (pz4c != nullptr) {
      resfile.Write_any_type(&(pz4c->last_output_time), sizeof(Real), "byte",
                             single_file_per_rank);
    }
    // output puncture tracker data
    if (nco > 0) {
      for (auto & pt : pz4c->ptracker) {
        resfile.Write_any_type(pt->GetPos(), 3*sizeof(Real), "byte",
                               single_file_per_rank);
      }
    }
    // turbulence driver internal RNG
    if (pturb != nullptr) {
      resfile.Write_any_type(&(pturb->rstate), sizeof(RNG_State), "byte",
                             single_file_per_rank);
    }
  }

  //--- STEP 4.  All ranks write data over all MeshBlocks (5D arrays) in parallel
  // This data read in ProblemGenerator constructor for restarts

  // total size of all cell-centered variables and face-centered fields to be written by
  // this rank
  IOWrapperSizeT data_size = 0;
  if (phydro != nullptr) {
    data_size += nout1*nout2*nout3*nhydro*sizeof(Real); // hydro u0
  }
  if (pmhd != nullptr) {
    data_size += nout1*nout2*nout3*nmhd*sizeof(Real);   // mhd u0
    data_size += (nout1+1)*nout2*nout3*sizeof(Real);    // mhd b0.x1f
    data_size += nout1*(nout2+1)*nout3*sizeof(Real);    // mhd b0.x2f
    data_size += nout1*nout2*(nout3+1)*sizeof(Real);    // mhd b0.x3f
  }
  if (prad != nullptr) {
    data_size += nout1*nout2*nout3*nrad*sizeof(Real);   // radiation i0
  }
  if (pturb != nullptr) {
    data_size += nout1*nout2*nout3*nforce*sizeof(Real); // forcing
  }
  if (pz4c != nullptr) {
    data_size += nout1*nout2*nout3*nz4c*sizeof(Real);   // z4c u0
  } else if (padm != nullptr) {
    data_size += nout1*nout2*nout3*nadm*sizeof(Real);   // adm u_adm
  }
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    resfile.Write_any_type(&(data_size), sizeof(IOWrapperSizeT), "byte",
                            single_file_per_rank);
  }

  // calculate size of data written in Steps 1-2 above
  IOWrapperSizeT step1size = sbuf.size()*sizeof(char) + 3*sizeof(int) + 2*sizeof(Real) +
                             sizeof(RegionSize) + 2*sizeof(RegionIndcs);
  IOWrapperSizeT step2size = (pm->nmb_total)*(sizeof(LogicalLocation) + sizeof(float));
  if (single_file_per_rank) {
    step2size += sizeof(kPerRankPartitionMagic) + 5*sizeof(int)
               + kCheckpointIdentityWords*sizeof(std::uint64_t);
  }
  if (pm->adaptive) {
    step2size += sizeof(kAmrCycleCounterMagic) + 2*sizeof(int)
               + pm->nmb_total*sizeof(int);
  }
  step2size += sizeof(kEventCounterMagic) + 2*sizeof(int)
             + kEventSumCounterCount*sizeof(std::uint64_t) + sizeof(int);

  IOWrapperSizeT step3size = 3*nco*sizeof(Real);
  if (pz4c != nullptr) step3size += sizeof(Real);
  if (pturb != nullptr) step3size += sizeof(RNG_State);

  // write cell-centered variables in parallel
  IOWrapperSizeT offset_myrank = (step1size + step2size + step3size
                                  + sizeof(IOWrapperSizeT));

  if (!single_file_per_rank) {
    offset_myrank += data_size*(pm->gids_eachrank[global_variable::my_rank]);
  }

  IOWrapperSizeT myoffset = offset_myrank;

  // write cell-centered variables, one MeshBlock at a time (but parallelized over all
  // ranks). MeshBlocks are written seperately to reduce number of data elements per write
  // call, to avoid exceeding 2^31 limit for very large grids per MPI rank.
  if (phydro != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_hyd, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered hydro data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_hyd, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered hydro data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nhydro*sizeof(Real); // hydro u0
    myoffset = offset_myrank;
  }
  if (pmhd != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_mhd, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered mhd data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_mhd, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered mhd data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nmhd*sizeof(Real);   // mhd u0
    myoffset = offset_myrank;

    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to x1-face field
        auto x1fptr = Kokkos::subview(outfield.x1f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        int fldcnt = x1fptr.size();
        if (resfile.Write_any_type_at_all(x1fptr.data(),fldcnt,myoffset,"Real",
                                          single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x1f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        // get ptr to x2-face field
        auto x2fptr = Kokkos::subview(outfield.x2f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = x2fptr.size();
        if (resfile.Write_any_type_at_all(x2fptr.data(),fldcnt,myoffset,"Real",
                                          single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x2f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        // get ptr to x3-face field
        auto x3fptr = Kokkos::subview(outfield.x3f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = x3fptr.size();
        if (resfile.Write_any_type_at_all(x3fptr.data(),fldcnt,myoffset,"Real",
                                          single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x3f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        myoffset += data_size-(x1fptr.size()+x2fptr.size()+x3fptr.size())*sizeof(Real);

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to x1-face field
        auto x1fptr = Kokkos::subview(outfield.x1f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        int fldcnt = x1fptr.size();
        if (resfile.Write_any_type_at(x1fptr.data(),fldcnt,myoffset,"Real",
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x1f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        // get ptr to x2-face field
        auto x2fptr = Kokkos::subview(outfield.x2f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = x2fptr.size();
        if (resfile.Write_any_type_at(x2fptr.data(),fldcnt,myoffset,"Real",
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x2f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        // get ptr to x3-face field
        auto x3fptr = Kokkos::subview(outfield.x3f,m,Kokkos::ALL,Kokkos::ALL,Kokkos::ALL);
        fldcnt = x3fptr.size();
        if (resfile.Write_any_type_at(x3fptr.data(),fldcnt,myoffset,"Real",
                                      single_file_per_rank) != fldcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "b0.x3f data not written correctly to rst file, "
                    << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += fldcnt*sizeof(Real);

        myoffset += data_size-(x1fptr.size()+x2fptr.size()+x3fptr.size())*sizeof(Real);
      }
    }
    offset_myrank += (nout1+1)*nout2*nout3*sizeof(Real);    // mhd b0.x1f
    offset_myrank += nout1*(nout2+1)*nout3*sizeof(Real);    // mhd b0.x2f
    offset_myrank += nout1*nout2*(nout3+1)*sizeof(Real);    // mhd b0.x3f
    myoffset = offset_myrank;
  }

  if (prad != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_rad, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered rad data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_rad, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(),mbcnt,myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered rad data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nrad*sizeof(Real);   // radiation i0
    myoffset = offset_myrank;
  }

  if (pturb != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_force, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
          << std::endl << "cell-centered turb data not written correctly to rst file, "
          << "restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_force, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered turb data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nforce*sizeof(Real); // forcing
    myoffset = offset_myrank;
  }

  if (pz4c != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_z4c, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered z4c data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_z4c, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered z4c data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nz4c*sizeof(Real); // z4c u0
    myoffset = offset_myrank;
  } else if (padm != nullptr) {
    for (int m=0;  m<noutmbs_max; ++m) {
      // every rank has a MB to write, so write collectively
      if (m < noutmbs_min) {
        // get ptr to cell-centered MeshBlock data
        auto mbptr = Kokkos::subview(outarray_adm, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at_all(mbptr.data(),mbcnt,myoffset,"Real",
                                          single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered adm data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;

      // some ranks are finished writing, so use non-collective write
      } else if (m < pm->nmb_thisrank) {
        // get ptr to MeshBlock data
        auto mbptr = Kokkos::subview(outarray_adm, m, Kokkos::ALL, Kokkos::ALL,
                                     Kokkos::ALL, Kokkos::ALL);
        int mbcnt = mbptr.size();
        if (resfile.Write_any_type_at(mbptr.data(), mbcnt, myoffset,"Real",
                                      single_file_per_rank) != mbcnt) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                    << std::endl << "cell-centered adm data not written correctly"
                    << " to rst file, restart file is broken." << std::endl;
          exit(EXIT_FAILURE);
        }
        myoffset += data_size;
      }
    }
    offset_myrank += nout1*nout2*nout3*nadm*sizeof(Real); // adm u_adm
    myoffset = offset_myrank;
  }

  // close file, clean up
  resfile.Close(single_file_per_rank);

  return;
}
