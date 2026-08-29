//========================================================================================
// AthenaK astrophysical plasma code
// Copyright(C) 2026 AthenaK collaboration
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file emri_inner_worldtube.cpp
//! \brief Replay validated face fluxes and interval edge EMFs on an inner cube.

#include "pgen/emri_inner_worldtube.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
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

constexpr std::uint64_t kHeaderBytes = 72;
constexpr std::uint32_t kBinaryVersion = 1;
constexpr int kFaceColumns = 7;
constexpr int kEdgeColumns = 7;
constexpr std::array<unsigned char, 16> kMagic = {
    'A', 'E', 'M', 'R', 'I', 'W', 'T', 'B', 'I', 'N', '0', '0', '0', '1', 0, 0};
constexpr std::array<int, 6> kNormalAxis = {0, 0, 1, 1, 2, 2};
constexpr std::array<int, 6> kNormalSign = {-1, 1, -1, 1, -1, 1};
constexpr std::array<int, 6> kUAxis = {1, 1, 2, 2, 0, 0};
constexpr std::array<int, 6> kUSign = {1, 1, 1, 1, 1, 1};
constexpr std::array<int, 6> kVAxis = {2, 2, 0, 0, 1, 1};
constexpr std::array<int, 6> kVSign = {-1, 1, -1, 1, -1, 1};

[[noreturn]] void InnerFatal(const std::string &message) {
  std::cerr << "### FATAL ERROR in EMRI inner worldtube replay on rank "
            << global_variable::my_rank << std::endl << message << std::endl;
#if MPI_PARALLEL_ENABLED
  MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
#endif
  std::exit(EXIT_FAILURE);
}

bool NearlyEqual(Real a, Real b, Real scale = 1.0) {
  const Real tolerance = 256.0*std::numeric_limits<Real>::epsilon()*
      std::max({std::abs(a), std::abs(b), std::abs(scale), Real(1.0)});
  return std::abs(a - b) <= tolerance;
}

std::uint32_t ReadLE32(const unsigned char *bytes) {
  return static_cast<std::uint32_t>(bytes[0]) |
      (static_cast<std::uint32_t>(bytes[1]) << 8) |
      (static_cast<std::uint32_t>(bytes[2]) << 16) |
      (static_cast<std::uint32_t>(bytes[3]) << 24);
}

std::uint64_t ReadLE64(const unsigned char *bytes) {
  std::uint64_t result = 0;
  for (int byte = 0; byte < 8; ++byte) {
    result |= static_cast<std::uint64_t>(bytes[byte]) << (8*byte);
  }
  return result;
}

double ReadLEDouble(const unsigned char *bytes) {
  const std::uint64_t bits = ReadLE64(bytes);
  double result = 0.0;
  std::memcpy(&result, &bits, sizeof(result));
  return result;
}

std::uint32_t UpdateCRC32(std::uint32_t crc, const unsigned char *bytes,
                          std::size_t count) {
  crc = ~crc;
  for (std::size_t index = 0; index < count; ++index) {
    crc ^= bytes[index];
    for (int bit = 0; bit < 8; ++bit) {
      const std::uint32_t mask = -(crc & 1U);
      crc = (crc >> 1) ^ (0xedb88320U & mask);
    }
  }
  return ~crc;
}

std::uint64_t CheckedProduct(std::initializer_list<std::uint64_t> factors,
                             const std::string &description) {
  std::uint64_t result = 1;
  for (std::uint64_t factor : factors) {
    if (factor != 0 && result > std::numeric_limits<std::uint64_t>::max()/factor) {
      InnerFatal(description + " overflows a 64-bit size");
    }
    result *= factor;
  }
  return result;
}

std::vector<Real> ReadDoublesAt(const std::string &path, std::uint64_t offset,
                                std::uint64_t count) {
  if (count > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max())) {
    InnerFatal("worldtube slab is too large for host address space");
  }
  std::ifstream stream(path, std::ios::binary);
  if (!stream.is_open()) InnerFatal("could not open inner replay file " + path);
  stream.seekg(static_cast<std::streamoff>(offset));
  if (!stream.good()) InnerFatal("could not seek in inner replay file");
  std::vector<unsigned char> bytes(static_cast<std::size_t>(count)*8);
  stream.read(reinterpret_cast<char*>(bytes.data()),
              static_cast<std::streamsize>(bytes.size()));
  if (stream.gcount() != static_cast<std::streamsize>(bytes.size())) {
    InnerFatal("inner replay file ended while reading a slab");
  }
  std::vector<Real> result(static_cast<std::size_t>(count));
  for (std::size_t index = 0; index < result.size(); ++index) {
    const double value = ReadLEDouble(bytes.data() + 8*index);
    if (!std::isfinite(value)) {
      InnerFatal("inner replay slab contains a non-finite value");
    }
    result[index] = static_cast<Real>(value);
  }
  return result;
}

}  // namespace

EmriInnerWorldtubeReplay::EmriInnerWorldtubeReplay(ParameterInput *pin, Mesh *pm,
                                                   bool is_restart) :
    pmesh_(pm), is_restart_(is_restart) {
  if (pm->pmb_pack == nullptr || pm->pmb_pack->pmhd == nullptr ||
      pm->pmb_pack->pdyngr == nullptr) {
    InnerFatal("inner worldtube replay currently requires dynamical GRMHD");
  }
  if (!pm->three_d || pm->adaptive || pm->multilevel) {
    InnerFatal("inner worldtube replay currently requires a fixed single-level 3D mesh");
  }
  path_ = pin->GetString("emri_worldtube", "file");
  flux_tolerance_ = pin->GetOrAddReal("emri_worldtube", "flux_tolerance", 1.0e-10);
  if (!(std::isfinite(flux_tolerance_) && flux_tolerance_ > 0.0)) {
    InnerFatal("<emri_worldtube>/flux_tolerance must be finite and positive");
  }
  ReadHeaderAndTimes(pin, pm);
  BuildBoundaryTopology(pm);
  LoadInterval(interval_);
  if (!is_restart_) {
    SetInitialNormalFlux(pm);
  } else {
    const Real residual = BoundaryFluxResidual(pm, flux_left_);
    if (residual > flux_tolerance_) {
      std::ostringstream message;
      message << "restart boundary flux differs from replay endpoint by " << residual
              << ", exceeding tolerance " << flux_tolerance_;
      InnerFatal(message.str());
    }
  }
}

void EmriInnerWorldtubeReplay::ReadHeaderAndTimes(ParameterInput *pin, Mesh *pm) {
  std::ifstream stream(path_, std::ios::binary | std::ios::ate);
  if (!stream.is_open()) InnerFatal("could not open inner replay file " + path_);
  const std::streamoff file_size_signed = stream.tellg();
  if (file_size_signed < static_cast<std::streamoff>(kHeaderBytes)) {
    InnerFatal("inner replay file is shorter than its header");
  }
  const std::uint64_t file_size = static_cast<std::uint64_t>(file_size_signed);
  stream.seekg(0);
  std::array<unsigned char, kHeaderBytes> header{};
  stream.read(reinterpret_cast<char*>(header.data()), header.size());
  if (stream.gcount() != static_cast<std::streamsize>(header.size())) {
    InnerFatal("failed to read inner replay header");
  }
  if (!std::equal(kMagic.begin(), kMagic.end(), header.begin())) {
    InnerFatal("inner replay magic is invalid");
  }
  const std::uint32_t version = ReadLE32(header.data() + 16);
  cells_per_edge_ = static_cast<int>(ReadLE32(header.data() + 20));
  nvar_ = static_cast<int>(ReadLE32(header.data() + 24));
  nt_ = static_cast<int>(ReadLE32(header.data() + 28));
  for (int axis = 0; axis < 3; ++axis) {
    data_center_[axis] = static_cast<Real>(ReadLEDouble(header.data() + 32 + 8*axis));
  }
  half_width_ = static_cast<Real>(ReadLEDouble(header.data() + 56));
  const std::uint64_t stored_crc = ReadLE64(header.data() + 64);
  if (version != kBinaryVersion || cells_per_edge_ < 1 || nvar_ < 1 || nt_ < 2 ||
      !(std::isfinite(half_width_) && half_width_ > 0.0) ||
      stored_crc > std::numeric_limits<std::uint32_t>::max()) {
    InnerFatal("inner replay header values are invalid or unsupported");
  }

  const std::uint64_t cells = static_cast<std::uint64_t>(cells_per_edge_);
  const std::uint64_t nt = static_cast<std::uint64_t>(nt_);
  const std::uint64_t state_count = CheckedProduct(
      {nt, static_cast<std::uint64_t>(nvar_), cells, cells}, "state array");
  const std::uint64_t flux_count = CheckedProduct({nt, cells, cells}, "flux array");
  const std::uint64_t emf_count = CheckedProduct(
      {nt - 1, cells + 1, cells}, "EMF array");
  std::uint64_t cursor = kHeaderBytes + 8*nt;
  for (int face = 0; face < 6; ++face) {
    state_offsets_[face] = cursor;
    cursor += 8*state_count;
    flux_offsets_[face] = cursor;
    cursor += 8*flux_count;
    emf_u_offsets_[face] = cursor;
    cursor += 8*emf_count;
    emf_v_offsets_[face] = cursor;
    cursor += 8*emf_count;
  }
  if (cursor != file_size) {
    InnerFatal("inner replay file size does not match its declared array dimensions");
  }

  stream.seekg(kHeaderBytes);
  std::array<unsigned char, 1 << 20> buffer{};
  std::uint32_t crc = 0;
  while (stream.good()) {
    stream.read(reinterpret_cast<char*>(buffer.data()), buffer.size());
    const std::streamsize count = stream.gcount();
    if (count > 0) {
      crc = UpdateCRC32(crc, buffer.data(), static_cast<std::size_t>(count));
    }
  }
  if (crc != static_cast<std::uint32_t>(stored_crc)) {
    InnerFatal("inner replay payload CRC32 mismatch");
  }

  const std::vector<Real> table_times = ReadDoublesAt(path_, kHeaderBytes, nt);
  times_ = table_times;
  for (int index = 0; index < nt_; ++index) {
    if (!std::isfinite(times_[index]) ||
        (index > 0 && !(times_[index] > times_[index - 1]))) {
      InnerFatal("inner replay times are not finite and strictly increasing");
    }
  }

  const std::string offset_text =
      pin->GetOrAddString("emri_worldtube", "time_offset", "auto");
  if (offset_text == "auto") {
    time_offset_ = times_.front() - pm->time;
  } else {
    std::size_t consumed = 0;
    try {
      time_offset_ = static_cast<Real>(std::stod(offset_text, &consumed));
    } catch (const std::exception&) {
      InnerFatal("<emri_worldtube>/time_offset is neither auto nor a number");
    }
    if (consumed != offset_text.size() || !std::isfinite(time_offset_)) {
      InnerFatal("<emri_worldtube>/time_offset is neither auto nor a finite number");
    }
  }
  const Real start_time = pm->time + time_offset_;
  int endpoint = -1;
  for (int index = 0; index < nt_; ++index) {
    if (NearlyEqual(start_time, times_[index], times_.back() - times_.front())) {
      endpoint = index;
      break;
    }
  }
  if (endpoint < 0 || endpoint >= nt_ - 1) {
    InnerFatal("inner run must start at a replay endpoint before the final table time");
  }
  interval_ = endpoint;

  mhd::MHD *pmhd = pm->pmb_pack->pmhd;
  const int fluid_nvar = pmhd->nmhd + pmhd->nscalars;
  has_cell_centered_magnetic_state_ = (nvar_ == fluid_nvar + 3);
  if (nvar_ != fluid_nvar && !has_cell_centered_magnetic_state_) {
    InnerFatal(
        "inner replay state must contain the compiled MHD primitives, with optional "
        "trailing bcc1,bcc2,bcc3");
  }
  const RegionSize &mesh = pm->mesh_size;
  dx_ = mesh.dx1;
  const Real lengths[3] = {mesh.x1max - mesh.x1min, mesh.x2max - mesh.x2min,
                           mesh.x3max - mesh.x3min};
  const Real mesh_center[3] = {0.5*(mesh.x1min + mesh.x1max),
                               0.5*(mesh.x2min + mesh.x2max),
                               0.5*(mesh.x3min + mesh.x3max)};
  for (int axis = 0; axis < 3; ++axis) {
    if (!NearlyEqual(data_center_[axis], mesh_center[axis], half_width_)) {
      InnerFatal(
          "replay and inner cube centers differ; prepare-inner repacks topology but "
          "does not perform a global-to-local coordinate/tetrad transformation");
    }
  }
  if (!NearlyEqual(mesh.dx2, dx_, dx_) || !NearlyEqual(mesh.dx3, dx_, dx_) ||
      !NearlyEqual(lengths[0], 2.0*half_width_, half_width_) ||
      !NearlyEqual(lengths[1], 2.0*half_width_, half_width_) ||
      !NearlyEqual(lengths[2], 2.0*half_width_, half_width_) ||
      !NearlyEqual(lengths[0]/dx_, cells_per_edge_, cells_per_edge_)) {
    InnerFatal("inner mesh must be the replay cube with matching isotropic resolution");
  }
}

void EmriInnerWorldtubeReplay::BuildBoundaryTopology(Mesh *pm) {
  const RegionIndcs &indcs = pm->mb_indcs;
  MeshBlockPack *pack = pm->pmb_pack;
  auto &sizes = pack->pmb->mb_size;
  const int cells = cells_per_edge_;
  std::vector<int> coverage(static_cast<std::size_t>(6)*cells*cells, 0);
  for (int m = 0; m < pack->nmb_thispack; ++m) {
    const RegionSize &size = sizes.h_view(m);
    const Real minimum[3] = {size.x1min, size.x2min, size.x3min};
    const Real maximum[3] = {size.x1max, size.x2max, size.x3max};
    const Real domain_minimum[3] = {pm->mesh_size.x1min, pm->mesh_size.x2min,
                                    pm->mesh_size.x3min};
    const Real domain_maximum[3] = {pm->mesh_size.x1max, pm->mesh_size.x2max,
                                    pm->mesh_size.x3max};
    const int axis_cells[3] = {indcs.nx1, indcs.nx2, indcs.nx3};
    const int axis_start[3] = {indcs.is, indcs.js, indcs.ks};
    for (int face = 0; face < 6; ++face) {
      const int normal_axis = kNormalAxis[face];
      const int normal_sign = kNormalSign[face];
      const Real boundary = (normal_sign < 0) ? domain_minimum[normal_axis]
                                               : domain_maximum[normal_axis];
      const Real block_boundary = (normal_sign < 0) ? minimum[normal_axis]
                                                     : maximum[normal_axis];
      if (!NearlyEqual(boundary, block_boundary, dx_)) continue;
      const int normal_offset = (normal_sign < 0) ? 0 : axis_cells[normal_axis] - 1;
      const int u_axis = kUAxis[face];
      const int v_axis = kVAxis[face];
      for (int v_offset = 0; v_offset < axis_cells[v_axis]; ++v_offset) {
        const Real v_lower = minimum[v_axis] + v_offset*dx_;
        int v = static_cast<int>(std::llround(
            ((kVSign[face] > 0) ? v_lower - domain_minimum[v_axis]
                                : domain_maximum[v_axis] - (v_lower + dx_))/dx_));
        for (int u_offset = 0; u_offset < axis_cells[u_axis]; ++u_offset) {
          const Real u_lower = minimum[u_axis] + u_offset*dx_;
          int u = static_cast<int>(std::llround(
              ((kUSign[face] > 0) ? u_lower - domain_minimum[u_axis]
                                  : domain_maximum[u_axis] - (u_lower + dx_))/dx_));
          if (u < 0 || u >= cells || v < 0 || v >= cells) {
            InnerFatal("inner MeshBlock produced an out-of-range face index");
          }
          int offsets[3] = {0, 0, 0};
          offsets[normal_axis] = normal_offset;
          offsets[u_axis] = u_offset;
          offsets[v_axis] = v_offset;
          host_faces_.push_back(FaceRecord{
              face, m, axis_start[2] + offsets[2], axis_start[1] + offsets[1],
              axis_start[0] + offsets[0], v, u});
          ++coverage[(static_cast<std::size_t>(face)*cells + v)*cells + u];
        }
      }
    }
  }
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, coverage.data(), static_cast<int>(coverage.size()), MPI_INT,
                MPI_SUM, MPI_COMM_WORLD);
#endif
  for (int count : coverage) {
    if (count != 1) InnerFatal("inner boundary face topology lacks a unique owner");
  }

  const std::size_t component_size = static_cast<std::size_t>(cells + 1)*cells;
  const std::size_t face_size = 2*component_size;
  using Target = std::tuple<int, int, int, int, int>;
  std::map<Target, EdgeRecord> unique_edges;
  auto add_edge = [&](const FaceRecord &cell, bool along_u, bool local_min,
                      int edge_v, int edge_u, std::size_t input) {
    const int face = cell.face;
    const int tangent_axis = along_u ? kUAxis[face] : kVAxis[face];
    const int tangent_sign = along_u ? kUSign[face] : kVSign[face];
    const int cross_axis = along_u ? kVAxis[face] : kUAxis[face];
    const int cross_sign = along_u ? kVSign[face] : kUSign[face];
    int indices[3] = {cell.i, cell.j, cell.k};
    indices[kNormalAxis[face]] += (kNormalSign[face] > 0) ? 1 : 0;
    const bool physical_lower = local_min ? (cross_sign > 0) : (cross_sign < 0);
    indices[cross_axis] += physical_lower ? 0 : 1;
    const Target target(cell.m, tangent_axis, indices[2], indices[1], indices[0]);
    const EdgeRecord record{static_cast<int>(input), cell.m, indices[2], indices[1],
                            indices[0], tangent_axis, tangent_sign};
    auto found = unique_edges.find(target);
    if (found == unique_edges.end() || record.input < found->second.input) {
      unique_edges[target] = record;
    }
    (void) edge_v;
    (void) edge_u;
  };
  for (const FaceRecord &cell : host_faces_) {
    const std::size_t base = static_cast<std::size_t>(cell.face)*face_size;
    const std::size_t vbase = base + component_size;
    add_edge(cell, true, true, cell.v, cell.u,
             base + static_cast<std::size_t>(cell.v)*cells + cell.u);
    add_edge(cell, true, false, cell.v + 1, cell.u,
             base + static_cast<std::size_t>(cell.v + 1)*cells + cell.u);
    add_edge(cell, false, true, cell.v, cell.u,
             vbase + static_cast<std::size_t>(cell.v)*(cells + 1) + cell.u);
    add_edge(cell, false, false, cell.v, cell.u + 1,
             vbase + static_cast<std::size_t>(cell.v)*(cells + 1) + cell.u + 1);
  }
  for (const auto &entry : unique_edges) host_edges_.push_back(entry.second);

  face_records_ = DvceArray2D<int>("emri_inner_faces",
                                   std::max<int>(host_faces_.size(), 1), kFaceColumns);
  edge_records_ = DvceArray2D<int>("emri_inner_edges",
                                   std::max<int>(host_edges_.size(), 1), kEdgeColumns);
  auto host_face_view = Kokkos::create_mirror_view(face_records_);
  for (std::size_t index = 0; index < host_faces_.size(); ++index) {
    const FaceRecord &record = host_faces_[index];
    const int values[kFaceColumns] = {record.face, record.m, record.k, record.j,
                                      record.i, record.v, record.u};
    for (int column = 0; column < kFaceColumns; ++column) {
      host_face_view(index, column) = values[column];
    }
  }
  Kokkos::deep_copy(face_records_, host_face_view);
  auto host_edge_view = Kokkos::create_mirror_view(edge_records_);
  for (std::size_t index = 0; index < host_edges_.size(); ++index) {
    const EdgeRecord &record = host_edges_[index];
    const int values[kEdgeColumns] = {record.input, record.m, record.k, record.j,
                                      record.i, record.component, record.sign};
    for (int column = 0; column < kEdgeColumns; ++column) {
      host_edge_view(index, column) = values[column];
    }
  }
  Kokkos::deep_copy(edge_records_, host_edge_view);

  const std::size_t face_cells = static_cast<std::size_t>(6)*cells*cells;
  const std::size_t face_edges = static_cast<std::size_t>(12)*cells*(cells + 1);
  state_left_ = DvceArray1D<Real>("emri_inner_state_left", face_cells*nvar_);
  state_right_ = DvceArray1D<Real>("emri_inner_state_right", face_cells*nvar_);
  flux_left_ = DvceArray1D<Real>("emri_inner_flux_left", face_cells);
  flux_right_ = DvceArray1D<Real>("emri_inner_flux_right", face_cells);
  interval_emf_ = DvceArray1D<Real>("emri_inner_interval_emf", face_edges);
}

void EmriInnerWorldtubeReplay::LoadInterval(int interval) {
  if (interval < 0 || interval >= nt_ - 1) {
    InnerFatal("requested inner replay interval is out of range");
  }
  const std::uint64_t cells = static_cast<std::uint64_t>(cells_per_edge_);
  const std::uint64_t state_slab = static_cast<std::uint64_t>(nvar_)*cells*cells;
  const std::uint64_t flux_slab = cells*cells;
  const std::uint64_t emf_slab = (cells + 1)*cells;
  auto host_state_left = Kokkos::create_mirror_view(state_left_);
  auto host_state_right = Kokkos::create_mirror_view(state_right_);
  auto host_flux_left = Kokkos::create_mirror_view(flux_left_);
  auto host_flux_right = Kokkos::create_mirror_view(flux_right_);
  auto host_emf = Kokkos::create_mirror_view(interval_emf_);
  for (int face = 0; face < 6; ++face) {
    const std::vector<Real> left_state = ReadDoublesAt(
        path_, state_offsets_[face] + 8*state_slab*interval, state_slab);
    const std::vector<Real> right_state = ReadDoublesAt(
        path_, state_offsets_[face] + 8*state_slab*(interval + 1), state_slab);
    const std::vector<Real> left_flux = ReadDoublesAt(
        path_, flux_offsets_[face] + 8*flux_slab*interval, flux_slab);
    const std::vector<Real> right_flux = ReadDoublesAt(
        path_, flux_offsets_[face] + 8*flux_slab*(interval + 1), flux_slab);
    const std::vector<Real> emf_u = ReadDoublesAt(
        path_, emf_u_offsets_[face] + 8*emf_slab*interval, emf_slab);
    const std::vector<Real> emf_v = ReadDoublesAt(
        path_, emf_v_offsets_[face] + 8*emf_slab*interval, emf_slab);
    const std::size_t state_base = static_cast<std::size_t>(face)*state_slab;
    const std::size_t flux_base = static_cast<std::size_t>(face)*flux_slab;
    const std::size_t emf_base = static_cast<std::size_t>(face)*2*emf_slab;
    for (std::size_t index = 0; index < state_slab; ++index) {
      host_state_left(state_base + index) = left_state[index];
      host_state_right(state_base + index) = right_state[index];
    }
    for (std::size_t index = 0; index < flux_slab; ++index) {
      host_flux_left(flux_base + index) = left_flux[index];
      host_flux_right(flux_base + index) = right_flux[index];
    }
    for (std::size_t index = 0; index < emf_slab; ++index) {
      host_emf(emf_base + index) = emf_u[index];
      host_emf(emf_base + emf_slab + index) = emf_v[index];
    }
  }
  Kokkos::deep_copy(state_left_, host_state_left);
  Kokkos::deep_copy(state_right_, host_state_right);
  Kokkos::deep_copy(flux_left_, host_flux_left);
  Kokkos::deep_copy(flux_right_, host_flux_right);
  Kokkos::deep_copy(interval_emf_, host_emf);
  interval_ = interval;
}

void EmriInnerWorldtubeReplay::SetInitialNormalFlux(Mesh *pm) {
  auto records = face_records_;
  auto expected = flux_left_;
  auto b1 = pm->pmb_pack->pmhd->b0.x1f;
  auto b2 = pm->pmb_pack->pmhd->b0.x2f;
  auto b3 = pm->pmb_pack->pmhd->b0.x3f;
  const int cells = cells_per_edge_;
  const Real inverse_area = 1.0/(dx_*dx_);
  const int count = static_cast<int>(host_faces_.size());
  if (count > 0) {
    par_for("emri_inner_set_flux", DevExeSpace(), 0, count - 1,
    KOKKOS_LAMBDA(int record_index) {
      const int face = records(record_index, 0);
      const int m = records(record_index, 1);
      const int k = records(record_index, 2);
      const int j = records(record_index, 3);
      const int i = records(record_index, 4);
      const int v = records(record_index, 5);
      const int u = records(record_index, 6);
      const int input = (face*cells + v)*cells + u;
      const int sign = (face % 2 == 0) ? -1 : 1;
      const Real field = sign*expected(input)*inverse_area;
      const int axis = face/2;
      if (axis == 0) b1(m, k, j, i + (sign > 0)) = field;
      if (axis == 1) b2(m, k, j + (sign > 0), i) = field;
      if (axis == 2) b3(m, k + (sign > 0), j, i) = field;
    });
  }
}

Real EmriInnerWorldtubeReplay::BoundaryFluxResidual(
    Mesh *pm, const DvceArray1D<Real> &expected) const {
  auto records = face_records_;
  auto b1 = pm->pmb_pack->pmhd->b0.x1f;
  auto b2 = pm->pmb_pack->pmhd->b0.x2f;
  auto b3 = pm->pmb_pack->pmhd->b0.x3f;
  const int cells = cells_per_edge_;
  const Real area = dx_*dx_;
  const int count = static_cast<int>(host_faces_.size());
  Real maximum = 0.0;
  if (count > 0) {
    Kokkos::parallel_reduce("emri_inner_flux_residual",
        Kokkos::RangePolicy<DevExeSpace>(0, count),
    KOKKOS_LAMBDA(int record_index, Real &thread_maximum) {
      const int face = records(record_index, 0);
      const int m = records(record_index, 1);
      const int k = records(record_index, 2);
      const int j = records(record_index, 3);
      const int i = records(record_index, 4);
      const int v = records(record_index, 5);
      const int u = records(record_index, 6);
      const int input = (face*cells + v)*cells + u;
      const int sign = (face % 2 == 0) ? -1 : 1;
      const int axis = face/2;
      Real field = 0.0;
      if (axis == 0) field = b1(m, k, j, i + (sign > 0));
      if (axis == 1) field = b2(m, k, j + (sign > 0), i);
      if (axis == 2) field = b3(m, k + (sign > 0), j, i);
      const Real difference = std::abs(sign*field*area - expected(input))/
                              (1.0 + std::abs(expected(input)));
      if (difference > thread_maximum) thread_maximum = difference;
    }, Kokkos::Max<Real>(maximum));
  }
#if MPI_PARALLEL_ENABLED
  MPI_Allreduce(MPI_IN_PLACE, &maximum, 1, MPI_ATHENA_REAL, MPI_MAX, MPI_COMM_WORLD);
#endif
  return maximum;
}

void EmriInnerWorldtubeReplay::InjectEField(Mesh *pm, Driver *pdrive, int stage) {
  if (exhausted_) InnerFatal("inner evolution advanced beyond the replay table");
  if (pdrive->LevelSubcyclingRequested() || pdrive->nimp_stages != 0 ||
      (pdrive->integrator != "rk1" && pdrive->integrator != "rk2" &&
       pdrive->integrator != "rk3")) {
    InnerFatal("inner replay supports explicit rk1/rk2/rk3 without level subcycling");
  }
  if (stage != last_stage_ + 1 || stage > pdrive->nexp_stages) {
    InnerFatal("inner replay observed an invalid RK-stage sequence");
  }
  auto records = edge_records_;
  auto input = interval_emf_;
  auto e1 = pm->pmb_pack->pmhd->efld.x1e;
  auto e2 = pm->pmb_pack->pmhd->efld.x2e;
  auto e3 = pm->pmb_pack->pmhd->efld.x3e;
  const Real inverse_length = 1.0/dx_;
  const int count = static_cast<int>(host_edges_.size());
  if (count > 0) {
    par_for("emri_inner_inject_emf", DevExeSpace(), 0, count - 1,
    KOKKOS_LAMBDA(int record_index) {
      const int source = records(record_index, 0);
      const int m = records(record_index, 1);
      const int k = records(record_index, 2);
      const int j = records(record_index, 3);
      const int i = records(record_index, 4);
      const int component = records(record_index, 5);
      const Real value = records(record_index, 6)*input(source)*inverse_length;
      if (component == 0) e1(m, k, j, i) = value;
      if (component == 1) e2(m, k, j, i) = value;
      if (component == 2) e3(m, k, j, i) = value;
    });
  }
  last_stage_ = stage;
}

void EmriInnerWorldtubeReplay::CompleteStep(Mesh *pm, Driver *pdrive) {
  if (last_stage_ != pdrive->nexp_stages) {
    InnerFatal("inner replay step ended before its final RK-stage injection");
  }
  last_stage_ = 0;
  const Real data_time = pm->time + time_offset_;
  const Real next_time = times_[interval_ + 1];
  if (data_time > next_time && !NearlyEqual(data_time, next_time, next_time)) {
    InnerFatal("inner timestep crossed a replay interval endpoint");
  }
  if (!NearlyEqual(data_time, next_time, next_time)) return;
  const Real residual = BoundaryFluxResidual(pm, flux_right_);
  if (global_variable::my_rank == 0) {
    std::cout << "EMRI inner worldtube: data_time=" << next_time
              << ", maximum normalized boundary-flux residual=" << residual
              << std::endl;
  }
  if (residual > flux_tolerance_) {
    std::ostringstream message;
    message << "replayed boundary flux residual " << residual
            << " exceeds tolerance " << flux_tolerance_
            << " at data time " << next_time;
    InnerFatal(message.str());
  }
  if (interval_ + 1 < nt_ - 1) {
    LoadInterval(interval_ + 1);
  } else {
    exhausted_ = true;
  }
}

void EmriInnerWorldtubeReplay::CapTimestep(Mesh *pm) {
  if (!(std::isfinite(pm->dt) && pm->dt > 0.0)) {
    InnerFatal("inner replay received a non-positive CFL timestep");
  }
  const Real data_time = pm->time + time_offset_;
  if (data_time >= times_.back() ||
      NearlyEqual(data_time, times_.back(), times_.back())) {
    return;
  }
  const Real next_time = times_[interval_ + 1];
  const Real remaining = next_time - data_time;
  if (!(remaining > 0.0)) InnerFatal("inner replay interval cursor is inconsistent");
  if (pm->dt > remaining && !NearlyEqual(pm->dt, remaining, remaining)) {
    pm->dt = remaining;
  }
}

void EmriInnerWorldtubeInjectEField(Mesh *pm, Driver *pdrive, int stage) {
  if (pm == nullptr || pm->pgen == nullptr ||
      pm->pgen->emri_inner_worldtube_ == nullptr) {
    InnerFatal("invalid inner-worldtube EMF callback state");
  }
  pm->pgen->emri_inner_worldtube_->InjectEField(pm, pdrive, stage);
}

void EmriInnerWorldtubeCompleteStep(Mesh *pm, Driver *pdrive) {
  if (pm == nullptr || pm->pgen == nullptr ||
      pm->pgen->emri_inner_worldtube_ == nullptr) {
    InnerFatal("invalid inner-worldtube step callback state");
  }
  pm->pgen->emri_inner_worldtube_->CompleteStep(pm, pdrive);
}

void EmriInnerWorldtubeCapTimestep(Mesh *pm) {
  if (pm == nullptr || pm->pgen == nullptr ||
      pm->pgen->emri_inner_worldtube_ == nullptr) {
    InnerFatal("invalid inner-worldtube timestep callback state");
  }
  pm->pgen->emri_inner_worldtube_->CapTimestep(pm);
}
