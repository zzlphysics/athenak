//========================================================================================
// AthenaK astrophysical plasma code
// Copyright(C) 2026 AthenaK collaboration
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file emri_outer_worldtube.cpp
//! \brief Record exact CT face fluxes and RK-integrated edge EMFs on a fixed cube.

#include "pgen/emri_outer_worldtube.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "athena.hpp"
#include "driver/driver.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "mesh/meshblock.hpp"
#include "mesh/meshblock_pack.hpp"
#include "mhd/mhd.hpp"
#include "parameter_input.hpp"
#include "pgen/pgen.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

constexpr int kNumFaces = 6;
constexpr int kCellColumns = 7;
constexpr int kEdgeColumns = 7;
constexpr std::array<const char*, kNumFaces> kFaceNames = {
    "x1m", "x1p", "x2m", "x2p", "x3m", "x3p"};
constexpr std::array<int, kNumFaces> kNormalAxis = {0, 0, 1, 1, 2, 2};
constexpr std::array<int, kNumFaces> kNormalSign = {-1, 1, -1, 1, -1, 1};
constexpr std::array<int, kNumFaces> kUAxis = {1, 1, 2, 2, 0, 0};
constexpr std::array<int, kNumFaces> kUSign = {1, 1, 1, 1, 1, 1};
constexpr std::array<int, kNumFaces> kVAxis = {2, 2, 0, 0, 1, 1};
constexpr std::array<int, kNumFaces> kVSign = {-1, 1, -1, 1, -1, 1};

[[noreturn]] void WorldtubeFatal(const std::string &message) {
  std::cerr << "### FATAL ERROR in EMRI outer worldtube recorder on rank "
            << global_variable::my_rank << std::endl << message << std::endl;
#if MPI_PARALLEL_ENABLED
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
  std::exit(EXIT_FAILURE);
}

bool NearlyEqual(Real a, Real b, Real scale) {
  const Real tolerance = 256.0*std::numeric_limits<Real>::epsilon()*
      std::max({std::abs(a), std::abs(b), std::abs(scale), Real(1.0)});
  return std::abs(a - b) <= tolerance;
}

int AlignedInteger(Real value, Real scale, const std::string &description) {
  const Real rounded = std::round(value);
  if (!NearlyEqual(value, rounded, scale)) {
    std::ostringstream message;
    message << description << " is not grid aligned: index="
            << std::setprecision(std::numeric_limits<Real>::max_digits10) << value;
    WorldtubeFatal(message.str());
  }
  if (rounded < static_cast<Real>(std::numeric_limits<int>::min()) ||
      rounded > static_cast<Real>(std::numeric_limits<int>::max())) {
    WorldtubeFatal(description + " exceeds the supported integer index range");
  }
  return static_cast<int>(rounded);
}

std::string JsonEscape(const std::string &input) {
  std::ostringstream result;
  for (unsigned char character : input) {
    switch (character) {
      case '\"': result << "\\\""; break;
      case '\\': result << "\\\\"; break;
      case '\b': result << "\\b"; break;
      case '\f': result << "\\f"; break;
      case '\n': result << "\\n"; break;
      case '\r': result << "\\r"; break;
      case '\t': result << "\\t"; break;
      default:
        if (character < 0x20) {
          result << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                 << static_cast<int>(character) << std::dec << std::setfill(' ');
        } else {
          result << static_cast<char>(character);
        }
    }
  }
  return result.str();
}

std::string FileNameOnly(const std::string &path) {
  return std::filesystem::path(path).filename().string();
}

std::uint64_t ByteSwap64(std::uint64_t value) {
  value = ((value & 0x00000000ffffffffULL) << 32) |
          ((value & 0xffffffff00000000ULL) >> 32);
  value = ((value & 0x0000ffff0000ffffULL) << 16) |
          ((value & 0xffff0000ffff0000ULL) >> 16);
  value = ((value & 0x00ff00ff00ff00ffULL) << 8) |
          ((value & 0xff00ff00ff00ff00ULL) >> 8);
  return value;
}

bool HostIsLittleEndian() {
  const std::uint16_t value = 1;
  return *reinterpret_cast<const unsigned char*>(&value) == 1;
}

void WriteLittleEndianDouble(std::ofstream &stream, double value) {
  static_assert(sizeof(double) == sizeof(std::uint64_t),
                "worldtube binary format requires a 64-bit double");
  std::uint64_t bits = 0;
  std::memcpy(&bits, &value, sizeof(bits));
  if (!HostIsLittleEndian()) bits = ByteSwap64(bits);
  stream.write(reinterpret_cast<const char*>(&bits), sizeof(bits));
  if (!stream.good()) WorldtubeFatal("failed while writing a worldtube binary stream");
}

void WriteLittleEndianVector(std::ofstream &stream, const std::vector<double> &values,
                             std::size_t begin, std::size_t count, double divisor = 1.0) {
  if (!(std::isfinite(divisor) && divisor != 0.0) || begin + count > values.size()) {
    WorldtubeFatal("invalid range or divisor in worldtube binary write");
  }
  for (std::size_t index = begin; index < begin + count; ++index) {
    WriteLittleEndianDouble(stream, values[index]/divisor);
  }
}

void ReduceDoubleVector(const std::vector<double> &local, std::vector<double> &global) {
  global.resize(local.size());
#if MPI_PARALLEL_ENABLED
  std::size_t offset = 0;
  while (offset < local.size()) {
    const std::size_t remaining = local.size() - offset;
    const int count = static_cast<int>(std::min<std::size_t>(
        remaining, static_cast<std::size_t>(std::numeric_limits<int>::max())));
    MPI_Reduce(local.data() + offset, global.data() + offset, count, MPI_DOUBLE, MPI_SUM,
               0, MPI_COMM_WORLD);
    offset += static_cast<std::size_t>(count);
  }
#else
  global = local;
#endif
}

void AllreduceIntVector(std::vector<int> &values) {
#if MPI_PARALLEL_ENABLED
  std::size_t offset = 0;
  while (offset < values.size()) {
    const std::size_t remaining = values.size() - offset;
    const int count = static_cast<int>(std::min<std::size_t>(
        remaining, static_cast<std::size_t>(std::numeric_limits<int>::max())));
    MPI_Allreduce(MPI_IN_PLACE, values.data() + offset, count, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    offset += static_cast<std::size_t>(count);
  }
#endif
}

int CellIndexOnFace(Real lower, Real upper, Real cube_min, Real cube_max, Real dx,
                    int direction_sign, const std::string &description) {
  if (lower < cube_min && !NearlyEqual(lower, cube_min, dx)) return -1;
  if (upper > cube_max && !NearlyEqual(upper, cube_max, dx)) return -1;
  int result = 0;
  if (direction_sign > 0) {
    result = AlignedInteger((lower - cube_min)/dx, 1.0, description);
  } else {
    result = AlignedInteger((cube_max - upper)/dx, 1.0, description);
  }
  return result;
}

}  // namespace

EmriOuterWorldtubeWriter::EmriOuterWorldtubeWriter(ParameterInput *pin, Mesh *pm,
                                                   bool is_restart) :
    pmesh_(pm), is_restart_(is_restart) {
  if (pm->pmb_pack == nullptr || pm->pmb_pack->pmhd == nullptr) {
    WorldtubeFatal("<emri_worldtube> requires the MHD module");
  }
  if (!pm->three_d) {
    WorldtubeFatal("<emri_worldtube> requires a three-dimensional mesh");
  }
  if (pm->adaptive || pm->multilevel) {
    WorldtubeFatal(
        "the first outer writer requires a fixed, single-level mesh; AMR and SMR are "
        "rejected rather than interpolated across the extraction surface");
  }

  center_[0] = pin->GetOrAddReal("emri_worldtube", "center_x1", 0.0);
  center_[1] = pin->GetOrAddReal("emri_worldtube", "center_x2", 0.0);
  center_[2] = pin->GetOrAddReal("emri_worldtube", "center_x3", 0.0);
  half_width_ = pin->GetReal("emri_worldtube", "half_width");
  dcycle_ = pin->GetOrAddInteger("emri_worldtube", "dcycle", 1);
  overwrite_ = pin->GetOrAddBoolean("emri_worldtube", "overwrite", false);
  basename_ = pin->GetOrAddString(
      "emri_worldtube", "file_basename",
      pin->GetOrAddString("job", "basename", "athena") + ".outer_worldtube");
  if (!(std::isfinite(half_width_) && half_width_ > 0.0)) {
    WorldtubeFatal("<emri_worldtube>/half_width must be finite and positive");
  }
  if (dcycle_ < 1) {
    WorldtubeFatal("<emri_worldtube>/dcycle must be at least one");
  }

  const RegionSize &mesh = pm->mesh_size;
  dx_ = mesh.dx1;
  if (!(std::isfinite(dx_) && dx_ > 0.0) || !NearlyEqual(mesh.dx2, dx_, dx_) ||
      !NearlyEqual(mesh.dx3, dx_, dx_)) {
    WorldtubeFatal("<emri_worldtube> currently requires isotropic Cartesian cells");
  }
  cells_per_edge_ = AlignedInteger(2.0*half_width_/dx_, 1.0,
                                   "worldtube side length");
  if (cells_per_edge_ < 1) {
    WorldtubeFatal("worldtube cube must span at least one grid cell per edge");
  }
  for (int axis = 0; axis < 3; ++axis) {
    const Real minimum =
        (axis == 0) ? mesh.x1min : ((axis == 1) ? mesh.x2min : mesh.x3min);
    AlignedInteger((center_[axis] - half_width_ - minimum)/dx_, 1.0,
                   "worldtube lower face");
    AlignedInteger((center_[axis] + half_width_ - minimum)/dx_, 1.0,
                   "worldtube upper face");
  }

  mhd::MHD *pmhd = pm->pmb_pack->pmhd;
  nvar_ = pmhd->nmhd + pmhd->nscalars;
  state_names_ = {"rho", "u1", "u2", "u3"};
  if (pmhd->nmhd == 5) state_names_.push_back("pgas");
  for (int scalar = 0; scalar < pmhd->nscalars; ++scalar) {
    state_names_.push_back("scalar" + std::to_string(scalar));
  }
  if (static_cast<int>(state_names_.size()) != nvar_) {
    WorldtubeFatal("unexpected MHD primitive layout in outer worldtube writer");
  }

  BuildTopology(pm);
}

EmriOuterWorldtubeWriter::~EmriOuterWorldtubeWriter() {
  CloseFiles();
}

void EmriOuterWorldtubeWriter::BuildTopology(Mesh *pm) {
  const RegionIndcs &indcs = pm->mb_indcs;
  const int axis_cells[3] = {indcs.nx1, indcs.nx2, indcs.nx3};
  const int axis_start[3] = {indcs.is, indcs.js, indcs.ks};
  MeshBlockPack *pack = pm->pmb_pack;
  auto &sizes = pack->pmb->mb_size;
  const int cells = cells_per_edge_;
  const Real cube_min[3] = {center_[0] - half_width_, center_[1] - half_width_,
                            center_[2] - half_width_};
  const Real cube_max[3] = {center_[0] + half_width_, center_[1] + half_width_,
                            center_[2] + half_width_};
  std::vector<int> cell_coverage(static_cast<std::size_t>(kNumFaces)*cells*cells, 0);

  for (int m = 0; m < pack->nmb_thispack; ++m) {
    const RegionSize &size = sizes.h_view(m);
    const Real block_min[3] = {size.x1min, size.x2min, size.x3min};
    const Real block_max[3] = {size.x1max, size.x2max, size.x3max};
    const Real block_dx[3] = {size.dx1, size.dx2, size.dx3};
    for (int axis = 0; axis < 3; ++axis) {
      if (!NearlyEqual(block_dx[axis], dx_, dx_)) {
        WorldtubeFatal("a local MeshBlock does not share the worldtube grid spacing");
      }
    }

    for (int face = 0; face < kNumFaces; ++face) {
      const int normal_axis = kNormalAxis[face];
      const int normal_sign = kNormalSign[face];
      const int u_axis = kUAxis[face];
      const int v_axis = kVAxis[face];
      const Real plane = center_[normal_axis] + normal_sign*half_width_;
      if (plane < block_min[normal_axis] &&
          !NearlyEqual(plane, block_min[normal_axis], dx_)) continue;
      if (plane > block_max[normal_axis] &&
          !NearlyEqual(plane, block_max[normal_axis], dx_)) continue;
      const int plane_index = AlignedInteger(
          (plane - block_min[normal_axis])/dx_, 1.0, "worldtube block-face plane");
      const int normal_offset = (normal_sign < 0) ? plane_index : plane_index - 1;
      if (normal_offset < 0 || normal_offset >= axis_cells[normal_axis]) continue;

      for (int v_offset = 0; v_offset < axis_cells[v_axis]; ++v_offset) {
        const Real v_lower = block_min[v_axis] + v_offset*dx_;
        const Real v_upper = v_lower + dx_;
        const int v = CellIndexOnFace(v_lower, v_upper, cube_min[v_axis],
                                      cube_max[v_axis], dx_, kVSign[face],
                                      "worldtube v cell");
        if (v < 0 || v >= cells) continue;
        for (int u_offset = 0; u_offset < axis_cells[u_axis]; ++u_offset) {
          const Real u_lower = block_min[u_axis] + u_offset*dx_;
          const Real u_upper = u_lower + dx_;
          const int u = CellIndexOnFace(u_lower, u_upper, cube_min[u_axis],
                                        cube_max[u_axis], dx_, kUSign[face],
                                        "worldtube u cell");
          if (u < 0 || u >= cells) continue;
          int offsets[3] = {0, 0, 0};
          offsets[normal_axis] = normal_offset;
          offsets[u_axis] = u_offset;
          offsets[v_axis] = v_offset;
          CellRecord record{face, m,
                            axis_start[2] + offsets[2],
                            axis_start[1] + offsets[1],
                            axis_start[0] + offsets[0], v, u};
          host_cells_.push_back(record);
          const std::size_t output =
              (static_cast<std::size_t>(face)*cells + v)*cells + u;
          ++cell_coverage[output];
        }
      }
    }
  }

  AllreduceIntVector(cell_coverage);
  for (std::size_t index = 0; index < cell_coverage.size(); ++index) {
    if (cell_coverage[index] != 1) {
      std::ostringstream message;
      message << "worldtube face-cell topology has coverage " << cell_coverage[index]
              << " at flattened cell " << index << "; expected exactly one owner";
      WorldtubeFatal(message.str());
    }
  }

  const std::size_t edges_per_face = 2ULL*cells*(cells + 1ULL);
  const std::size_t total_edges = kNumFaces*edges_per_face;
  std::vector<int> edge_coverage(total_edges, 0);
  auto add_edge = [&](const CellRecord &cell, bool along_u, bool local_min,
                      int edge_v, int edge_u, std::size_t output) {
    const int face = cell.face;
    const int tangent_axis = along_u ? kUAxis[face] : kVAxis[face];
    const int tangent_sign = along_u ? kUSign[face] : kVSign[face];
    const int cross_axis = along_u ? kVAxis[face] : kUAxis[face];
    const int cross_sign = along_u ? kVSign[face] : kUSign[face];
    int indices[3] = {cell.i, cell.j, cell.k};
    indices[kNormalAxis[face]] += (kNormalSign[face] > 0) ? 1 : 0;
    const bool physical_lower = local_min ? (cross_sign > 0) : (cross_sign < 0);
    indices[cross_axis] += physical_lower ? 0 : 1;
    host_edges_.push_back(
        EdgeRecord{static_cast<int>(output), cell.m, indices[2], indices[1],
                   indices[0], tangent_axis, tangent_sign});
    ++edge_coverage[output];
    (void) edge_v;
    (void) edge_u;
  };

  for (const CellRecord &cell : host_cells_) {
    const std::size_t base = static_cast<std::size_t>(cell.face)*edges_per_face;
    const std::size_t u_base = base;
    const std::size_t v_base = base + static_cast<std::size_t>(cells + 1)*cells;
    add_edge(cell, true, true, cell.v, cell.u,
             u_base + static_cast<std::size_t>(cell.v)*cells + cell.u);
    add_edge(cell, false, true, cell.v, cell.u,
             v_base + static_cast<std::size_t>(cell.v)*(cells + 1) + cell.u);
    if (cell.v == cells - 1) {
      add_edge(cell, true, false, cells, cell.u,
               u_base + static_cast<std::size_t>(cells)*cells + cell.u);
    }
    if (cell.u == cells - 1) {
      add_edge(cell, false, false, cell.v, cells,
               v_base + static_cast<std::size_t>(cell.v)*(cells + 1) + cells);
    }
  }
  AllreduceIntVector(edge_coverage);
  for (std::size_t index = 0; index < edge_coverage.size(); ++index) {
    if (edge_coverage[index] != 1) {
      std::ostringstream message;
      message << "worldtube face-edge topology has coverage " << edge_coverage[index]
              << " at flattened edge " << index << "; expected exactly one owner";
      WorldtubeFatal(message.str());
    }
  }

  const int local_cells = static_cast<int>(host_cells_.size());
  const int local_edges = static_cast<int>(host_edges_.size());
  cell_records_ = DvceArray2D<int>("emri_worldtube_cells", std::max(local_cells, 1),
                                   kCellColumns);
  edge_records_ = DvceArray2D<int>("emri_worldtube_edges", std::max(local_edges, 1),
                                   kEdgeColumns);
  auto host_cell_view = Kokkos::create_mirror_view(cell_records_);
  for (int index = 0; index < local_cells; ++index) {
    const CellRecord &record = host_cells_[index];
    const int values[kCellColumns] = {record.face, record.m, record.k, record.j,
                                      record.i, record.v, record.u};
    for (int column = 0; column < kCellColumns; ++column) {
      host_cell_view(index, column) = values[column];
    }
  }
  Kokkos::deep_copy(cell_records_, host_cell_view);
  auto host_edge_view = Kokkos::create_mirror_view(edge_records_);
  for (int index = 0; index < local_edges; ++index) {
    const EdgeRecord &record = host_edges_[index];
    const int values[kEdgeColumns] = {record.output, record.m, record.k, record.j,
                                      record.i, record.component, record.sign};
    for (int column = 0; column < kEdgeColumns; ++column) {
      host_edge_view(index, column) = values[column];
    }
  }
  Kokkos::deep_copy(edge_records_, host_edge_view);

  const std::size_t total_face_cells =
      static_cast<std::size_t>(kNumFaces)*cells*cells;
  step_emf_integral_ = DvceArray1D<Real>("emri_worldtube_step_emf", total_edges);
  interval_emf_integral_ =
      DvceArray1D<Real>("emri_worldtube_interval_emf", total_edges);
  endpoint_state_ =
      DvceArray1D<Real>("emri_worldtube_state", total_face_cells*nvar_);
  endpoint_flux_ = DvceArray1D<Real>("emri_worldtube_flux", total_face_cells);
  Kokkos::deep_copy(step_emf_integral_, 0.0);
  Kokkos::deep_copy(interval_emf_integral_, 0.0);
}

void EmriOuterWorldtubeWriter::OpenFiles(Mesh *pm, Driver *pdrive) {
  if (initialized_) return;
  if (pdrive->LevelSubcyclingRequested() || pdrive->nimp_stages != 0 ||
      (pdrive->integrator != "rk1" && pdrive->integrator != "rk2" &&
       pdrive->integrator != "rk3")) {
    WorldtubeFatal(
        "the first outer writer supports explicit rk1/rk2/rk3 with "
        "<time>/subcycling=none; other recurrences need separate verification");
  }

  std::ostringstream suffix;
  suffix << ".cycle" << std::setw(8) << std::setfill('0') << pm->ncycle;
  stem_ = basename_ + suffix.str();
  manifest_path_ = stem_ + ".manifest.json";
  times_path_ = stem_ + ".times.bin";
  const std::filesystem::path parent = std::filesystem::path(stem_).parent_path();
  if (global_variable::my_rank == 0 && !parent.empty()) {
    std::error_code error;
    std::filesystem::create_directories(parent, error);
    if (error) WorldtubeFatal("could not create worldtube output directory: " +
                              error.message());
  }
#if MPI_PARALLEL_ENABLED
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if (global_variable::my_rank == 0) {
    std::vector<std::string> paths = {manifest_path_, times_path_};
    for (const char *face : kFaceNames) {
      paths.push_back(stem_ + "." + face + ".cell_state.bin");
      paths.push_back(stem_ + "." + face + ".normal_flux.bin");
      paths.push_back(stem_ + "." + face + ".emf_u.bin");
      paths.push_back(stem_ + "." + face + ".emf_v.bin");
    }
    if (!overwrite_) {
      for (const std::string &path : paths) {
        if (std::filesystem::exists(path)) {
          WorldtubeFatal("refusing to overwrite existing worldtube file " + path +
                         "; set <emri_worldtube>/overwrite=true explicitly");
        }
      }
    }
    times_stream_.open(times_path_, std::ios::binary | std::ios::trunc);
    if (!times_stream_.is_open()) WorldtubeFatal("could not open " + times_path_);
    for (int face = 0; face < kNumFaces; ++face) {
      const std::string prefix = stem_ + "." + kFaceNames[face];
      face_streams_[face].cell_state.open(prefix + ".cell_state.bin",
                                          std::ios::binary | std::ios::trunc);
      face_streams_[face].normal_flux.open(prefix + ".normal_flux.bin",
                                           std::ios::binary | std::ios::trunc);
      face_streams_[face].emf_u.open(prefix + ".emf_u.bin",
                                     std::ios::binary | std::ios::trunc);
      face_streams_[face].emf_v.open(prefix + ".emf_v.bin",
                                     std::ios::binary | std::ios::trunc);
      if (!face_streams_[face].cell_state.is_open() ||
          !face_streams_[face].normal_flux.is_open() ||
          !face_streams_[face].emf_u.is_open() ||
          !face_streams_[face].emf_v.is_open()) {
        WorldtubeFatal("could not open one or more worldtube face streams");
      }
    }
  }
  interval_start_time_ = pm->time;
  initialized_ = true;
  CaptureAndWriteEndpoint(pm, pm->time);
  WriteManifest(false);
}

void EmriOuterWorldtubeWriter::ObserveEField(Mesh *pm, Driver *pdrive, int stage) {
  if (finalized_) WorldtubeFatal("worldtube EMF observer called after finalization");
  if (!initialized_) {
    if (stage != 1) {
      WorldtubeFatal("first worldtube EMF observation was not RK stage one");
    }
    OpenFiles(pm, pdrive);
  }
  if (stage != last_stage_ + 1 || stage < 1 || stage > pdrive->nexp_stages) {
    WorldtubeFatal("worldtube recorder observed an invalid RK-stage sequence");
  }
  if (!(std::isfinite(pm->dt) && pm->dt > 0.0)) {
    WorldtubeFatal("worldtube recorder observed a non-positive timestep");
  }

  mhd::MHD *pmhd = pm->pmb_pack->pmhd;
  auto records = edge_records_;
  auto integral = step_emf_integral_;
  auto e1 = pmhd->efld.x1e;
  auto e2 = pmhd->efld.x2e;
  auto e3 = pmhd->efld.x3e;
  const int local_edges = static_cast<int>(host_edges_.size());
  const Real gam0 = pdrive->gam0[stage - 1];
  const Real beta_dt = pdrive->beta[stage - 1]*pm->dt;
  const Real dx = dx_;
  if (stage == 1) Kokkos::deep_copy(integral, 0.0);
  if (local_edges > 0) {
    par_for("emri_worldtube_accumulate_emf", DevExeSpace(), 0, local_edges - 1,
    KOKKOS_LAMBDA(int record_index) {
      const int output = records(record_index, 0);
      const int m = records(record_index, 1);
      const int k = records(record_index, 2);
      const int j = records(record_index, 3);
      const int i = records(record_index, 4);
      const int component = records(record_index, 5);
      const Real sign = static_cast<Real>(records(record_index, 6));
      Real electric = 0.0;
      if (component == 0) electric = e1(m, k, j, i);
      if (component == 1) electric = e2(m, k, j, i);
      if (component == 2) electric = e3(m, k, j, i);
      integral(output) = gam0*integral(output) + beta_dt*sign*electric*dx;
    });
  }
  last_stage_ = stage;
}

void EmriOuterWorldtubeWriter::CompleteStep(Mesh *pm, Driver *pdrive) {
  if (!initialized_ || finalized_) return;
  if (last_stage_ != pdrive->nexp_stages) {
    WorldtubeFatal("worldtube step callback did not follow the final RK-stage EMF");
  }
  auto interval = interval_emf_integral_;
  auto step = step_emf_integral_;
  const int total_edges = interval.extent_int(0);
  par_for("emri_worldtube_join_step", DevExeSpace(), 0, total_edges - 1,
  KOKKOS_LAMBDA(int index) {
    interval(index) += step(index);
  });
  ++interval_steps_;
  last_stage_ = 0;
  if (interval_steps_ >= dcycle_) {
    const Real interval_dt = pm->time - interval_start_time_;
    WriteInterval(pm, interval_dt);
    CaptureAndWriteEndpoint(pm, pm->time);
    interval_start_time_ = pm->time;
    interval_steps_ = 0;
    Kokkos::deep_copy(interval_emf_integral_, 0.0);
    WriteManifest(false);
  }
}

void EmriOuterWorldtubeWriter::CaptureAndWriteEndpoint(Mesh *pm, Real time) {
  mhd::MHD *pmhd = pm->pmb_pack->pmhd;
  auto records = cell_records_;
  auto state = endpoint_state_;
  auto flux = endpoint_flux_;
  auto w0 = pmhd->w0;
  auto b1 = pmhd->b0.x1f;
  auto b2 = pmhd->b0.x2f;
  auto b3 = pmhd->b0.x3f;
  const int local_cells = static_cast<int>(host_cells_.size());
  const int cells = cells_per_edge_;
  const int nvar = nvar_;
  const Real area = dx_*dx_;
  Kokkos::deep_copy(state, 0.0);
  Kokkos::deep_copy(flux, 0.0);
  if (local_cells > 0) {
    par_for("emri_worldtube_capture_endpoint", DevExeSpace(), 0, local_cells - 1,
    KOKKOS_LAMBDA(int record_index) {
      const int face = records(record_index, 0);
      const int m = records(record_index, 1);
      const int k = records(record_index, 2);
      const int j = records(record_index, 3);
      const int i = records(record_index, 4);
      const int v = records(record_index, 5);
      const int u = records(record_index, 6);
      const int face_cell = (face*cells + v)*cells + u;
      for (int variable = 0; variable < nvar; ++variable) {
        const int state_index = ((face*nvar + variable)*cells + v)*cells + u;
        state(state_index) = w0(m, variable, k, j, i);
      }
      const int normal_axis = face/2;
      const int normal_sign = (face % 2 == 0) ? -1 : 1;
      Real normal_field = 0.0;
      if (normal_axis == 0) normal_field = b1(m, k, j, i + (normal_sign > 0));
      if (normal_axis == 1) normal_field = b2(m, k, j + (normal_sign > 0), i);
      if (normal_axis == 2) normal_field = b3(m, k + (normal_sign > 0), j, i);
      flux(face_cell) = normal_sign*normal_field*area;
    });
  }
  auto host_state = Kokkos::create_mirror_view_and_copy(HostMemSpace(), state);
  auto host_flux = Kokkos::create_mirror_view_and_copy(HostMemSpace(), flux);
  std::vector<double> local_state(host_state.extent(0));
  std::vector<double> local_flux(host_flux.extent(0));
  for (std::size_t index = 0; index < local_state.size(); ++index) {
    local_state[index] = static_cast<double>(host_state(index));
  }
  for (std::size_t index = 0; index < local_flux.size(); ++index) {
    local_flux[index] = static_cast<double>(host_flux(index));
  }
  std::vector<double> global_state;
  std::vector<double> global_flux;
  ReduceDoubleVector(local_state, global_state);
  ReduceDoubleVector(local_flux, global_flux);
  if (global_variable::my_rank == 0) {
    WriteLittleEndianDouble(times_stream_, static_cast<double>(time));
    const std::size_t face_cells = static_cast<std::size_t>(cells)*cells;
    for (int face = 0; face < kNumFaces; ++face) {
      WriteLittleEndianVector(face_streams_[face].cell_state, global_state,
                              static_cast<std::size_t>(face)*nvar_*face_cells,
                              static_cast<std::size_t>(nvar_)*face_cells);
      WriteLittleEndianVector(face_streams_[face].normal_flux, global_flux,
                              static_cast<std::size_t>(face)*face_cells, face_cells);
    }
  }
  ++endpoints_written_;
}

void EmriOuterWorldtubeWriter::WriteInterval(Mesh *pm, Real interval_dt) {
  (void) pm;
  if (!(std::isfinite(interval_dt) && interval_dt > 0.0)) {
    WorldtubeFatal("worldtube output interval must be finite and positive");
  }
  auto host_integral = Kokkos::create_mirror_view_and_copy(
      HostMemSpace(), interval_emf_integral_);
  std::vector<double> local(host_integral.extent(0));
  for (std::size_t index = 0; index < local.size(); ++index) {
    local[index] = static_cast<double>(host_integral(index));
  }
  std::vector<double> global;
  ReduceDoubleVector(local, global);
  if (global_variable::my_rank == 0) {
    const std::size_t component_size =
        static_cast<std::size_t>(cells_per_edge_ + 1)*cells_per_edge_;
    const std::size_t face_size = 2*component_size;
    for (int face = 0; face < kNumFaces; ++face) {
      const std::size_t base = static_cast<std::size_t>(face)*face_size;
      WriteLittleEndianVector(face_streams_[face].emf_u, global, base, component_size,
                              static_cast<double>(interval_dt));
      WriteLittleEndianVector(face_streams_[face].emf_v, global,
                              base + component_size, component_size,
                              static_cast<double>(interval_dt));
    }
  }
  ++intervals_written_;
}

void EmriOuterWorldtubeWriter::WriteManifest(bool complete) const {
  if (global_variable::my_rank != 0) return;
  const std::string temporary = manifest_path_ + ".tmp";
  std::ofstream stream(temporary, std::ios::trunc);
  if (!stream.is_open()) WorldtubeFatal("could not open temporary worldtube manifest");
  stream << std::setprecision(std::numeric_limits<Real>::max_digits10);
  stream << "{\n"
         << "  \"classification\": \"athenak-emri-outer-worldtube-stream-v1\",\n"
         << "  \"target_classification\": "
         << "\"athenak-emri-cubical-flux-emf-worldtube-v1\",\n"
         << "  \"complete\": " << (complete ? "true" : "false") << ",\n"
         << "  \"binary_dtype\": \"<f8\",\n"
         << "  \"times_file\": \"" << JsonEscape(FileNameOnly(times_path_))
         << "\",\n"
         << "  \"nt\": " << endpoints_written_ << ",\n"
         << "  \"ninterval\": " << intervals_written_ << ",\n"
         << "  \"nvar\": " << nvar_ << ",\n"
         << "  \"cells_per_face_axis\": " << cells_per_edge_ << ",\n"
         << "  \"center\": [" << center_[0] << ", " << center_[1] << ", "
         << center_[2] << "],\n"
         << "  \"half_width\": " << half_width_ << ",\n"
         << "  \"grid_spacing\": " << dx_ << ",\n"
         << "  \"fixed_grid_aligned_cube\": true,\n"
         << "  \"restart_segment\": " << (is_restart_ ? "true" : "false") << ",\n"
         << "  \"state_sampling\": "
         << "\"interior-adjacent primitive cell average\",\n"
         << "  \"emf_sampling\": "
         << "\"synchronized CT line average with exact explicit-RK recurrence\",\n"
         << "  \"state_variables\": [";
  for (std::size_t index = 0; index < state_names_.size(); ++index) {
    if (index > 0) stream << ", ";
    stream << "\"" << JsonEscape(state_names_[index]) << "\"";
  }
  stream << "],\n  \"faces\": {\n";
  for (int face = 0; face < kNumFaces; ++face) {
    const std::string prefix = stem_ + "." + kFaceNames[face];
    stream << "    \"" << kFaceNames[face] << "\": {\n"
           << "      \"cell_state\": \""
           << JsonEscape(FileNameOnly(prefix + ".cell_state.bin")) << "\",\n"
           << "      \"normal_flux\": \""
           << JsonEscape(FileNameOnly(prefix + ".normal_flux.bin")) << "\",\n"
           << "      \"emf_u\": \""
           << JsonEscape(FileNameOnly(prefix + ".emf_u.bin")) << "\",\n"
           << "      \"emf_v\": \""
           << JsonEscape(FileNameOnly(prefix + ".emf_v.bin")) << "\"\n"
           << "    }" << ((face + 1 < kNumFaces) ? "," : "") << "\n";
  }
  stream << "  }\n}\n";
  stream.close();
  if (!stream.good()) WorldtubeFatal("failed to write worldtube manifest");
  std::error_code error;
  std::filesystem::remove(manifest_path_, error);
  error.clear();
  std::filesystem::rename(temporary, manifest_path_, error);
  if (error) WorldtubeFatal("could not install worldtube manifest: " + error.message());
}

void EmriOuterWorldtubeWriter::Finalize(Mesh *pm, Driver *pdrive) {
  (void) pdrive;
  if (finalized_) return;
  if (!initialized_) {
    finalized_ = true;
    return;
  }
  if (last_stage_ != 0) {
    WorldtubeFatal("worldtube finalization occurred inside an RK step");
  }
  if (interval_steps_ > 0) {
    const Real interval_dt = pm->time - interval_start_time_;
    WriteInterval(pm, interval_dt);
    CaptureAndWriteEndpoint(pm, pm->time);
    interval_start_time_ = pm->time;
    interval_steps_ = 0;
    Kokkos::deep_copy(interval_emf_integral_, 0.0);
  }
  WriteManifest(endpoints_written_ >= 2 && intervals_written_ + 1 == endpoints_written_);
  CloseFiles();
  finalized_ = true;
}

void EmriOuterWorldtubeWriter::CloseFiles() {
  if (times_stream_.is_open()) times_stream_.close();
  for (FaceStreams &streams : face_streams_) {
    if (streams.cell_state.is_open()) streams.cell_state.close();
    if (streams.normal_flux.is_open()) streams.normal_flux.close();
    if (streams.emf_u.is_open()) streams.emf_u.close();
    if (streams.emf_v.is_open()) streams.emf_v.close();
  }
}

void EmriOuterWorldtubeObserveEField(Mesh *pm, Driver *pdrive, int stage) {
  if (pm == nullptr || pm->pgen == nullptr ||
      pm->pgen->emri_outer_worldtube_ == nullptr) {
    WorldtubeFatal("invalid outer-worldtube EMF callback state");
  }
  pm->pgen->emri_outer_worldtube_->ObserveEField(pm, pdrive, stage);
}

void EmriOuterWorldtubeCompleteStep(Mesh *pm, Driver *pdrive) {
  if (pm == nullptr || pm->pgen == nullptr ||
      pm->pgen->emri_outer_worldtube_ == nullptr) {
    WorldtubeFatal("invalid outer-worldtube step callback state");
  }
  pm->pgen->emri_outer_worldtube_->CompleteStep(pm, pdrive);
}

void EmriOuterWorldtubeFinalize(Mesh *pm, Driver *pdrive) {
  if (pm == nullptr || pm->pgen == nullptr ||
      pm->pgen->emri_outer_worldtube_ == nullptr) {
    WorldtubeFatal("invalid outer-worldtube finalize callback state");
  }
  pm->pgen->emri_outer_worldtube_->Finalize(pm, pdrive);
}
