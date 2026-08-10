//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file build_tree.cpp
//! \brief Functions to build MeshBlock, both for new runs and restarts

#include <algorithm>
#include <iostream>
#include <cinttypes>
#include <cmath>
#include <cstdlib>
#include <limits> // numeric_limits<>
#include <memory> // make_unique<>
#include <string>
#include <vector>

#include "athena.hpp"
#include "globals.hpp"
#include "parameter_input.hpp"
#include "mesh.hpp"
#include "coordinates/cell_locations.hpp"
#include "hydro/hydro.hpp"
#include "mhd/mhd.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

constexpr std::uint64_t kAmrCycleCounterMagic = UINT64_C(0x41544b414d524331);
constexpr int kAmrCycleCounterVersion = 1;
constexpr std::uint64_t kPerRankPartitionMagic = UINT64_C(0x41544b5052543031);
constexpr int kPerRankPartitionVersion = 2;
constexpr int kCheckpointIdentityWords = 2;

bool LevelSubcyclingRequested(ParameterInput *pin) {
  return pin->DoesParameterExist("time", "subcycling") &&
         pin->GetString("time", "subcycling") == "level";
}

float LevelSubcyclingBlockCost(const LogicalLocation &lloc, int root_level) {
  return std::ldexp(1.0f, lloc.level - root_level);
}

[[noreturn]] void BuildTreeError(const std::string &message) {
  if (global_variable::my_rank == 0) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << message << std::endl;
  }
#if MPI_PARALLEL_ENABLED
  if (global_variable::nranks > 1) {
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }
#endif
  std::exit(EXIT_FAILURE);
}

int ValidatedMaxBlocksPerRank(ParameterInput *pin, bool adaptive,
                              const int *nmb_eachrank) {
  int max_blocks = nmb_eachrank[global_variable::my_rank];
  if (adaptive) {
    if (!pin->DoesParameterExist("mesh_refinement", "max_nmb_per_rank")) {
      BuildTreeError("With AMR, <mesh_refinement>/max_nmb_per_rank must be specified");
    }
    max_blocks = pin->GetInteger("mesh_refinement", "max_nmb_per_rank");
    if (max_blocks < 1) {
      BuildTreeError("<mesh_refinement>/max_nmb_per_rank must be a positive integer");
    }
    for (int rank=0; rank<global_variable::nranks; ++rank) {
      if (nmb_eachrank[rank] > max_blocks) {
        BuildTreeError("Initial weighted partition assigns " +
                       std::to_string(nmb_eachrank[rank]) + " MeshBlocks to rank " +
                       std::to_string(rank) + ", exceeding max_nmb_per_rank=" +
                       std::to_string(max_blocks));
      }
    }
  }

#if MPI_PARALLEL_ENABLED
  int largest_partition = 0;
  for (int rank=0; rank<global_variable::nranks; ++rank) {
    largest_partition = std::max(largest_partition, nmb_eachrank[rank]);
  }
  if (std::max(max_blocks, largest_partition) > (1 << NUM_BITS_LID)) {
    BuildTreeError("Maximum number of MeshBlocks per rank exceeds the MPI tag limit");
  }
#endif
  return max_blocks;
}

}  // namespace

//----------------------------------------------------------------------------------------
//! \fn void Mesh::BuildTreeFromScratch():
//! Constructs MeshBlockTree, creates MeshBlockPack (containing the physics modules), and
//! divides grid into MeshBlock(s) for new runs (starting from scratch), using parameters
//! read from input file.  Also does initial load balance based on simple cost estimate.

void Mesh::BuildTreeFromScratch(ParameterInput *pin) {
  // calculate the number of MeshBlocks at root level in each dir
  nmb_rootx1 = mesh_indcs.nx1/mb_indcs.nx1;
  nmb_rootx2 = mesh_indcs.nx2/mb_indcs.nx2;
  nmb_rootx3 = mesh_indcs.nx3/mb_indcs.nx3;

  // find maximum number of MeshBlocks at root level in any dir
  int nmbmax = (nmb_rootx1 > nmb_rootx2) ? nmb_rootx1 : nmb_rootx2;
  nmbmax = (nmbmax > nmb_rootx3) ? nmbmax : nmb_rootx3;

  // Find smallest N such that 2^N > max number of MeshBlocks in any dimension (nmbmax)
  // Then N is logical level of root grid.  2^N implemented as left-shift (1<<root_level)
  for (root_level=0; ((1<<root_level) < nmbmax); root_level++) {}
  int current_level = root_level;

  // Construct tree and create root grid
  ptree = std::make_unique<MeshBlockTree>(this);
  ptree->CreateRootGrid();

  // Error check properties of input paraemters for SMR/AMR meshes.
  if (adaptive) {
    max_level = pin->GetOrAddInteger("mesh_refinement", "num_levels", 1) + root_level - 1;
    if (max_level > 31) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Number of refinement levels must be smaller than "
                << 31 - root_level + 1 << std::endl;
      std::exit(EXIT_FAILURE);
    }
  } else {
    max_level = 31;
  }

  // Read <refined_region> blocks and construct tree accordingly
  // These regions can be used with both SMR (in which case they will remain fixed) and
  // AMR (in which case they may be defined, unless the location refinement criteria used)
  if (multilevel) {
    // error check that number of cells in MeshBlock divisible by two
    if (mb_indcs.nx1 % 2 != 0 ||
       (mb_indcs.nx2 % 2 != 0 && multi_d) ||
       (mb_indcs.nx3 % 2 != 0 && three_d)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Number of cells in MeshBlock must be divisible by 2 "
                << "with SMR or AMR." << std::endl;
      std::exit(EXIT_FAILURE);
    }

    // cycle through ParameterInput list and find "refined_region" blocks, extract data
    // and expand MeshBlockTree
    for (auto it = pin->block.begin(); it != pin->block.end(); ++it) {
      if (it->block_name.compare(0, 14, "refined_region") == 0) {
        RegionSize ref_size;
        ref_size.x1min = pin->GetReal(it->block_name, "x1min");
        ref_size.x1max = pin->GetReal(it->block_name, "x1max");
        if (multi_d) {
          ref_size.x2min = pin->GetReal(it->block_name, "x2min");
          ref_size.x2max = pin->GetReal(it->block_name, "x2max");
        } else {
          ref_size.x2min = mesh_size.x2min;
          ref_size.x2max = mesh_size.x2max;
        }
        if (three_d) {
          ref_size.x3min = pin->GetReal(it->block_name, "x3min");
          ref_size.x3max = pin->GetReal(it->block_name, "x3max");
        } else {
          ref_size.x3min = mesh_size.x3min;
          ref_size.x3max = mesh_size.x3max;
        }
        int phy_ref_lev = pin->GetInteger(it->block_name, "level");
        int log_ref_lev = phy_ref_lev + root_level;
        if (log_ref_lev > current_level) current_level = log_ref_lev;

        // error check parameters in "refinement" blocks
        if (phy_ref_lev < 1) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl <<"<refined_region> level must be larger than 0 (root level=0)"
              << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (log_ref_lev > max_level) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<refined_region> level exceeds maximum allowed ("
              << max_level << ")" << std::endl << "Reduce/specify 'num_levels' in "
              << "<mesh_refinement> input block if using AMR" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (   ref_size.x1min > ref_size.x1max
            || ref_size.x2min > ref_size.x2max
            || ref_size.x3min > ref_size.x3max)  {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "Invalid <refined_region> (xmax < xmin in one direction)."
              << std::endl;
          std::exit(EXIT_FAILURE);
        }
        if (   ref_size.x1min < mesh_size.x1min || ref_size.x1max > mesh_size.x1max
            || ref_size.x2min < mesh_size.x2min || ref_size.x2max > mesh_size.x2max
            || ref_size.x3min < mesh_size.x3min || ref_size.x3max > mesh_size.x3max) {
          std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
              << std::endl << "<refined_region> must be fully contained within root mesh"
              << std::endl;
          std::exit(EXIT_FAILURE);
        }

        // note: if following is too slow, it could be replaced with bi-section search.
        // Suppose entire root domain is tiled with MeshBlocks at the desired refinement
        // level. Find range of x1-integer indices of such MeshBlocks that cover the
        // refinement region
        std::int32_t lx1min = 0, lx1max = 0;
        std::int32_t lx2min = 0, lx2max = 0;
        std::int32_t lx3min = 0, lx3max = 0;
        std::int32_t lxmax = nmb_rootx1*(1<<phy_ref_lev);
        for (lx1min=0; lx1min<lxmax; lx1min++) {
          if (LeftEdgeX(lx1min+1,lxmax,mesh_size.x1min,mesh_size.x1max) > ref_size.x1min)
            break;
        }
        for (lx1max=lx1min; lx1max<lxmax; lx1max++) {
          if (LeftEdgeX(lx1max+1,lxmax,mesh_size.x1min,mesh_size.x1max) >= ref_size.x1max)
            break;
        }
        if (lx1min % 2 == 1) lx1min--;
        if (lx1max % 2 == 0) lx1max++;

        // Find range of x2-indices of such MeshBlocks that cover the refinement region
        if (multi_d) { // 2D or 3D
          lxmax = nmb_rootx2*(1<<phy_ref_lev);
          for (lx2min=0; lx2min<lxmax; lx2min++) {
            if (LeftEdgeX(lx2min+1, lxmax, mesh_size.x2min, mesh_size.x2max) >
                ref_size.x2min)
            break;
          }
          for (lx2max=lx2min; lx2max<lxmax; lx2max++) {
            if (LeftEdgeX(lx2max+1, lxmax, mesh_size.x2min, mesh_size.x2max) >=
                ref_size.x2max)
            break;
          }
          if (lx2min % 2 == 1) lx2min--;
          if (lx2max % 2 == 0) lx2max++;
        }

        // Find range of x3-indices of such MeshBlocks that cover the refinement region
        if (three_d) { // 3D
          lxmax = nmb_rootx3*(1<<phy_ref_lev);
          for (lx3min=0; lx3min<lxmax; lx3min++) {
            if (LeftEdgeX(lx3min+1, lxmax, mesh_size.x3min, mesh_size.x3max) >
                ref_size.x3min)
            break;
          }
          for (lx3max=lx3min; lx3max<lxmax; lx3max++) {
            if (LeftEdgeX(lx3max+1, lxmax, mesh_size.x3min, mesh_size.x3max) >=
                ref_size.x3max)
            break;
          }
          if (lx3min % 2 == 1) lx3min--;
          if (lx3max % 2 == 0) lx3max++;
        }

        // Now add nodes to the MeshBlockTree corresponding to these MeshBlocks
        if (one_d) {  // 1D
          for (std::int32_t i=lx1min; i<lx1max; i+=2) {
            LogicalLocation nlloc;
            nlloc.level = log_ref_lev;
            nlloc.lx1 = i;
            nlloc.lx2 = 0;
            nlloc.lx3 = 0;
            int nnew;
            ptree->AddNode(nlloc, nnew);
          }
        }
        if (two_d) {  // 2D
          for (std::int32_t j=lx2min; j<lx2max; j+=2) {
            for (std::int32_t i=lx1min; i<lx1max; i+=2) {
              LogicalLocation nlloc;
              nlloc.level = log_ref_lev;
              nlloc.lx1 = i;
              nlloc.lx2 = j;
              nlloc.lx3 = 0;
              int nnew;
              ptree->AddNode(nlloc, nnew);
            }
          }
        }
        if (three_d) {  // 3D
          for (std::int32_t k=lx3min; k<lx3max; k+=2) {
            for (std::int32_t j=lx2min; j<lx2max; j+=2) {
              for (std::int32_t i=lx1min; i<lx1max; i+=2) {
                LogicalLocation nlloc;
                nlloc.level = log_ref_lev;
                nlloc.lx1 = i;
                nlloc.lx2 = j;
                nlloc.lx3 = k;
                int nnew;
                ptree->AddNode(nlloc, nnew);
              }
            }
          }
        }
      }
    }
  } // if (multilevel)

  if (!adaptive) max_level = current_level;

  // initial mesh hierarchy construction is completed here
  ptree->CountMeshBlocks(nmb_total);

  cost_eachmb = new float[nmb_total];
  rank_eachmb = new int[nmb_total];
  lloc_eachmb = new LogicalLocation[nmb_total];
  gids_eachrank = new int[global_variable::nranks];
  nmb_eachrank = new int[global_variable::nranks];

  // following returns LogicalLocation list sorted by Z-ordering, and total # of MBs
  ptree->CreateZOrderedLLList(lloc_eachmb, nullptr, nmb_total);

#if MPI_PARALLEL_ENABLED
  // check there is at least one MeshBlock per MPI rank
  if (nmb_total < global_variable::nranks) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
        << "Fewer MeshBlocks (nmb_total=" << nmb_total << ") than MPI ranks (nranks="
        << global_variable::nranks << ")" << std::endl;
    std::exit(EXIT_FAILURE);
  }
#endif

  // A level-l block advances 2^(l-root) times per synchronized root step.  Weighting
  // static partitions by that update frequency avoids assigning equal cost to blocks
  // whose strict-subcycling work differs by powers of two.
  const bool level_subcycling = LevelSubcyclingRequested(pin);
  for (int i=0; i<nmb_total; i++) {
    cost_eachmb[i] = level_subcycling
        ? LevelSubcyclingBlockCost(lloc_eachmb[i], root_level) : 1.0f;
  }
  LoadBalance(cost_eachmb, rank_eachmb, gids_eachrank, nmb_eachrank, nmb_total);
  nmb_maxperrank = ValidatedMaxBlocksPerRank(pin, adaptive, nmb_eachrank);

  // create MeshBlockPack for this rank
  int mbp_gids = gids_eachrank[global_variable::my_rank];
  int mbp_gide = mbp_gids + nmb_eachrank[global_variable::my_rank] - 1;
  nmb_thisrank = nmb_eachrank[global_variable::my_rank];

  pmb_pack = new MeshBlockPack(this, mbp_gids, mbp_gide);
  nmb_packs_thisrank = 1;
  pmb_pack->AddMeshBlocks(pin);
  pmb_pack->pmb->SetNeighbors(ptree, rank_eachmb);

  // Create new MeshRefinement object with either SMR or AMR (SMR needs Restrict fns)
  if (multilevel) {
    pmr = new MeshRefinement(this, pin);
  }

  // set initial time/cycle parameters, output diagnostics
  time = pin->GetOrAddReal("time", "start_time", 0.0);
  dt   = std::numeric_limits<float>::max();
  cfl_no = pin->GetReal("time", "cfl_number");
  ncycle = 0;
  if (global_variable::my_rank == 0) {PrintMeshDiagnostics();}

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void Mesh::BuildTreeFromRestart():
//! Constructs MeshBlockTree, creates MeshBlockPack (containing the physics modules), and
//! divides grid into MeshBlock(s) for restart runs, using parameters and data read from
//! restart file.

void Mesh::BuildTreeFromRestart(ParameterInput *pin, IOWrapper &resfile,
                                                     bool single_file_per_rank) {
  // At this point, the restartfile is already open and the ParameterInput (input file)
  // data has already been read in main(). Thus the file pointer is set to after <par_end>

  // following must be identical to calculation of headeroffset (excluding size of
  // ParameterInput data) in restart.cpp
  IOWrapperSizeT headersize = 3*sizeof(int) + 2*sizeof(Real)
    + sizeof(RegionSize) + 2*sizeof(RegionIndcs);
  char *headerdata = new char[headersize];

  // the master process reads the header data if single_file_per_rank is false
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    IOWrapperSizeT read_size = resfile.Read_bytes(headerdata, 1, headersize,
                                                  single_file_per_rank, true);
    if (read_size != headersize) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Header size read from restart file is incorrect, "
                << "expected " << headersize << ", got " << read_size << std::endl;
      BuildTreeError("Restart header is truncated");
    }
  }

#if MPI_PARALLEL_ENABLED
  // then broadcast the header data
  if (!single_file_per_rank) {
    int mpi_err = MPI_Bcast(headerdata, headersize, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (mpi_err != MPI_SUCCESS) {
      char error_string[1024];
      int length_of_error_string;
      MPI_Error_string(mpi_err, error_string, &length_of_error_string);
      std::cout << "MPI_Bcast failed with error: " << error_string << std::endl;
      BuildTreeError("Could not broadcast the restart header");
    }
  }
#endif

  // Now copy mesh data read from restart file into Mesh variables. Order of variables
  // set by Write()'s in restart.cpp
  // Note this overwrites size and indices initialized in Mesh constructor.
  IOWrapperSizeT hdos = 0;
  std::memcpy(&nmb_total, &(headerdata[hdos]), sizeof(int));
  hdos += sizeof(int);
  std::memcpy(&root_level, &(headerdata[hdos]), sizeof(int));
  hdos += sizeof(int);
  std::memcpy(&mesh_size, &(headerdata[hdos]), sizeof(RegionSize));
  hdos += sizeof(RegionSize);
  std::memcpy(&mesh_indcs, &(headerdata[hdos]), sizeof(RegionIndcs));
  hdos += sizeof(RegionIndcs);
  std::memcpy(&mb_indcs, &(headerdata[hdos]), sizeof(RegionIndcs));
  hdos += sizeof(RegionIndcs);
  std::memcpy(&time, &(headerdata[hdos]), sizeof(Real));
  hdos += sizeof(Real);
  std::memcpy(&dt, &(headerdata[hdos]), sizeof(Real));
  hdos += sizeof(Real);
  std::memcpy(&ncycle, &(headerdata[hdos]), sizeof(int));
  delete [] headerdata;
  const Real checkpoint_dt = dt;
  dt_last_completed =
      (ncycle > 0 && std::isfinite(checkpoint_dt) && checkpoint_dt > 0.0 &&
       checkpoint_dt != std::numeric_limits<float>::max()) ? checkpoint_dt : 0.0;

  // New checkpoints serialize the pre-tlim CFL/growth reference independently from
  // the actual last completed step above.  This keeps dt_last_completed physically
  // correct for time-derived diagnostics while preventing an endpoint roundoff tail
  // from throttling the resumed run.  For an old checkpoint, discard only a timestep
  // at floating-point roundoff scale; the freshly recomputed CFL limit still bounds the
  // next step, so this compatibility path cannot hide a genuinely restrictive CFL.
  const Real fresh_dt = std::numeric_limits<float>::max();
  dt_restart_growth = fresh_dt;
  if (ncycle > 0) {
    if (pin->DoesParameterExist("time", "restart_dt_growth")) {
      dt_restart_growth = pin->GetReal("time", "restart_dt_growth");
      if (!std::isfinite(dt_restart_growth) || !(dt_restart_growth > 0.0)) {
        BuildTreeError("Restart contains an invalid <time>/restart_dt_growth");
      }
    } else if (dt_last_completed > 0.0) {
      const Real scale = std::max(static_cast<Real>(1.0), std::abs(time));
      const Real roundoff_tail =
          static_cast<Real>(64.0)*std::numeric_limits<Real>::epsilon()*scale;
      if (dt_last_completed > roundoff_tail) {
        dt_restart_growth = dt_last_completed;
      }
    }
  }
  dt = dt_restart_growth;

  // calculate the number of MeshBlocks at root level in each dir
  nmb_rootx1 = mesh_indcs.nx1/mb_indcs.nx1;
  nmb_rootx2 = mesh_indcs.nx2/mb_indcs.nx2;
  nmb_rootx3 = mesh_indcs.nx3/mb_indcs.nx3;
  int current_level = root_level;
  std::vector<int> restart_cycle_counters;

  // Error check properties of input paraemters for SMR/AMR meshes.
  if (adaptive) {
    max_level = pin->GetOrAddInteger("mesh_refinement", "num_levels", 1) + root_level - 1;
    if (max_level > 31) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Number of refinement levels must be smaller than "
                << 31 - root_level + 1 << std::endl;
      BuildTreeError("Restart requests too many refinement levels");
    }
  } else {
    max_level = 31;
  }

  // allocate memory for lists read from restart
  cost_eachmb = new float[nmb_total];
  rank_eachmb = new int[nmb_total];
  lloc_eachmb = new LogicalLocation[nmb_total];
  gids_eachrank = new int[global_variable::nranks];
  nmb_eachrank = new int[global_variable::nranks];

  // allocate idlist buffer and read list of logical locations and cost
  IOWrapperSizeT listsize = sizeof(LogicalLocation) + sizeof(float);
  char *idlist = new char[listsize*nmb_total];
  // only the master process reads the ID list
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    if (resfile.Read_bytes(idlist, listsize, nmb_total, single_file_per_rank, true) !=
        static_cast<unsigned int>(nmb_total)) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
                << std::endl << "Incorrect number of MeshBlocks in restart file; "
                << "restart file is broken." << std::endl;
      BuildTreeError("Restart MeshBlock list is truncated");
    }
  }
#if MPI_PARALLEL_ENABLED
  // then broadcast the ID list
  if (!single_file_per_rank) {
    MPI_Bcast(idlist, listsize*nmb_total, MPI_CHAR, 0, MPI_COMM_WORLD);
  }
#endif

  // everyone sets the logical location and cost lists based on bradcasted data
  int os = 0;
  for (int i=0; i<nmb_total; i++) {
    std::memcpy(&(lloc_eachmb[i]), &(idlist[os]), sizeof(LogicalLocation));
    os += sizeof(LogicalLocation);
  }
  for (int i=0; i<nmb_total; i++) {
    std::memcpy(&(cost_eachmb[i]), &(idlist[os]), sizeof(float));
    os += sizeof(float);
    if (lloc_eachmb[i].level > current_level) current_level = lloc_eachmb[i].level;
  }
  delete [] idlist;

  // Each per-rank file contains only its writer's MeshBlock field stream.  Bind that
  // stream to the exact writer partition.  Legacy files without this metadata cannot
  // be restored safely because the original rank count and GID range are unknowable.
  int writer_gid_start = -1;
  int writer_nmb = -1;
  std::uint64_t checkpoint_identity[kCheckpointIdentityWords] = {0, 0};
  if (single_file_per_rank) {
    std::uint64_t partition_magic = 0;
    int partition_version = 0;
    int writer_nranks = 0;
    int writer_rank = -1;
    const bool partition_read_ok =
        resfile.Read_bytes(&partition_magic, 1, sizeof(partition_magic), true, true) ==
            sizeof(partition_magic) &&
        resfile.Read_bytes(&partition_version, 1, sizeof(partition_version), true,
                           true) ==
            sizeof(partition_version) &&
        resfile.Read_bytes(&writer_nranks, 1, sizeof(writer_nranks), true, true) ==
            sizeof(writer_nranks) &&
        resfile.Read_bytes(&writer_rank, 1, sizeof(writer_rank), true, true) ==
            sizeof(writer_rank) &&
        resfile.Read_bytes(&writer_gid_start, 1, sizeof(writer_gid_start), true, true) ==
            sizeof(writer_gid_start) &&
        resfile.Read_bytes(&writer_nmb, 1, sizeof(writer_nmb), true, true) ==
            sizeof(writer_nmb) &&
        resfile.Read_bytes(checkpoint_identity, 1, sizeof(checkpoint_identity), true,
                           true) == sizeof(checkpoint_identity);
    if (!partition_read_ok || partition_magic != kPerRankPartitionMagic ||
        partition_version != kPerRankPartitionVersion) {
      BuildTreeError("Per-rank restart lacks valid writer-partition metadata; "
                     "legacy per-rank files cannot be restored safely");
    }
    if (writer_nranks != global_variable::nranks ||
        writer_rank != global_variable::my_rank || writer_gid_start < 0 ||
        writer_nmb < 1 || writer_gid_start + writer_nmb > nmb_total) {
      BuildTreeError("Per-rank restart writer/current MPI rank layout is incompatible");
    }
#if MPI_PARALLEL_ENABLED
    std::uint64_t rank_zero_identity[kCheckpointIdentityWords] = {0, 0};
    if (global_variable::my_rank == 0) {
      rank_zero_identity[0] = checkpoint_identity[0];
      rank_zero_identity[1] = checkpoint_identity[1];
    }
    MPI_Bcast(rank_zero_identity, sizeof(rank_zero_identity), MPI_BYTE, 0,
              MPI_COMM_WORLD);
    int identity_matches =
        (checkpoint_identity[0] == rank_zero_identity[0] &&
         checkpoint_identity[1] == rank_zero_identity[1]) ? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &identity_matches, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    if (identity_matches == 0) {
      BuildTreeError("Per-rank checkpoint identity mismatch across rank files");
    }
#endif
  }

  // Versioned metadata is appended to the legacy location/cost list.  Probe the file
  // marker itself instead of a ParameterInput key: restart parameters can legitimately
  // be overridden from an input file or the command line.  If the marker is absent,
  // restore the old-format file position before reading the remaining object state.
  int counter_extension_present = 0;
  if (global_variable::my_rank == 0 || single_file_per_rank) {
    const IOWrapperSizeT extension_offset =
        resfile.GetPosition(single_file_per_rank);
    std::uint64_t counter_magic = 0;
    const std::size_t magic_bytes =
        resfile.Read_bytes(&counter_magic, 1, sizeof(counter_magic),
                           single_file_per_rank);
    if (magic_bytes == sizeof(counter_magic) &&
        counter_magic == kAmrCycleCounterMagic) {
      counter_extension_present = 1;
    } else if (resfile.Seek(extension_offset, single_file_per_rank) != 0) {
      BuildTreeError("Could not restore the legacy restart file position");
    }
  }
#if MPI_PARALLEL_ENABLED
  if (single_file_per_rank) {
    int extension_min = 0;
    int extension_max = 0;
    MPI_Allreduce(&counter_extension_present, &extension_min, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&counter_extension_present, &extension_max, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
    if (extension_min != extension_max) {
      BuildTreeError("Per-rank restart files disagree on AMR counter metadata");
    }
  } else {
    MPI_Bcast(&counter_extension_present, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
#endif

  if (counter_extension_present != 0) {
    int stored_version = 0;
    int counter_count = 0;
    int metadata_read_ok = 1;
    restart_cycle_counters.resize(nmb_total);
    if (global_variable::my_rank == 0 || single_file_per_rank) {
      metadata_read_ok =
          (resfile.Read_bytes(&stored_version, 1, sizeof(stored_version),
                              single_file_per_rank, true) == sizeof(stored_version)) &&
          (resfile.Read_bytes(&counter_count, 1, sizeof(counter_count),
                              single_file_per_rank, true) == sizeof(counter_count));
      if (metadata_read_ok && stored_version == kAmrCycleCounterVersion &&
          counter_count == nmb_total) {
        metadata_read_ok =
            resfile.Read_bytes(restart_cycle_counters.data(), sizeof(int), nmb_total,
                               single_file_per_rank, true) ==
            static_cast<std::size_t>(nmb_total);
      }
    }
#if MPI_PARALLEL_ENABLED
    if (!single_file_per_rank) {
      MPI_Bcast(&stored_version, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&counter_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&metadata_read_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
#endif
    if (!metadata_read_ok || stored_version != kAmrCycleCounterVersion ||
        counter_count != nmb_total) {
      BuildTreeError("AMR cycle-counter metadata in restart is invalid");
    }
#if MPI_PARALLEL_ENABLED
    if (!single_file_per_rank) {
      MPI_Bcast(restart_cycle_counters.data(), nmb_total, MPI_INT, 0, MPI_COMM_WORLD);
    }
#endif
  }
  if (!adaptive) max_level = current_level;

  // Checkpoints made before level-aware costs were introduced contain a uniform unit
  // cost list.  Upgrade only that recognizable legacy case; otherwise preserve any
  // measured or already-weighted costs stored in the restart.
  if (LevelSubcyclingRequested(pin)) {
    bool legacy_unit_costs = true;
    for (int i=0; i<nmb_total; ++i) {
      legacy_unit_costs = legacy_unit_costs && (cost_eachmb[i] == 1.0f);
    }
    if (legacy_unit_costs) {
      for (int i=0; i<nmb_total; ++i) {
        cost_eachmb[i] = LevelSubcyclingBlockCost(lloc_eachmb[i], root_level);
      }
    }
  }

  // rebuild the MeshBlockTree
  ptree = std::make_unique<MeshBlockTree>(this);
  ptree->CreateRootGrid();
  for (int i=0; i<nmb_total; i++) {ptree->AddNodeWithoutRefinement(lloc_eachmb[i]);}

  // check the tree structure by making sure total # of MBs counted in tree same as the
  // number read from the restart file.
  {
    int nnb;
    ptree->CreateZOrderedLLList(lloc_eachmb, nullptr, nnb);
    if (nnb != nmb_total) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__
        << std::endl << "Tree reconstruction failed. Total number of blocks in "
        << "reconstructed tree=" << nnb << ", number in file=" << nmb_total << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

#ifdef MPI_PARALLEL_ENABLED
  // check there is at least one MeshBlock per MPI rank
  if (!single_file_per_rank) {
    if (nmb_total < global_variable::nranks) {
      std::cout << "### FATAL ERROR in " << __FILE__ << " at line "
        << __LINE__ << std::endl
        << "Fewer MeshBlocks (nmb_total=" << nmb_total << ") than MPI ranks (nranks="
        << global_variable::nranks << ")" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
#endif

  LoadBalance(cost_eachmb, rank_eachmb, gids_eachrank, nmb_eachrank, nmb_total);
  if (single_file_per_rank) {
    int partition_matches =
        (gids_eachrank[global_variable::my_rank] == writer_gid_start &&
         nmb_eachrank[global_variable::my_rank] == writer_nmb) ? 1 : 0;
#if MPI_PARALLEL_ENABLED
    MPI_Allreduce(MPI_IN_PLACE, &partition_matches, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
#endif
    if (partition_matches == 0) {
      BuildTreeError("Per-rank restart writer GID partition differs from the current "
                     "load-balanced partition");
    }
  }
  nmb_maxperrank = ValidatedMaxBlocksPerRank(pin, adaptive, nmb_eachrank);

  // create MeshBlockPack for this rank
  int mbp_gids = gids_eachrank[global_variable::my_rank];
  int mbp_gide = mbp_gids + nmb_eachrank[global_variable::my_rank] - 1;
  nmb_thisrank = nmb_eachrank[global_variable::my_rank];

  pmb_pack = new MeshBlockPack(this, mbp_gids, mbp_gide);
  pmb_pack->AddMeshBlocks(pin);
  pmb_pack->pmb->SetNeighbors(ptree, rank_eachmb);

  // Create new MeshRefinement object with either SMR or AMR (SMR needs Restrict fns)
  if (multilevel) {
    pmr = new MeshRefinement(this, pin);
    if (!restart_cycle_counters.empty()) {
      for (int m=0; m<nmb_total; ++m) {
        pmr->ncyc_since_ref(m) = restart_cycle_counters[m];
      }
    }
  }

  // set remaining parameters, output diagnostics
  cfl_no = pin->GetReal("time", "cfl_number");
  if (global_variable::my_rank == 0) {PrintMeshDiagnostics();}
}
