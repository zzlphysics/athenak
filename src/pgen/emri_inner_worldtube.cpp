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
#include "coordinates/adm.hpp"
#include "coordinates/cell_locations.hpp"
#include "driver/driver.hpp"
#include "dyn_grmhd/dyn_grmhd.hpp"
#include "globals.hpp"
#include "mesh/mesh.hpp"
#include "mesh/meshblock.hpp"
#include "mesh/meshblock_pack.hpp"
#include "eos/eos.hpp"
#include "mhd/mhd.hpp"
#include "parameter_input.hpp"
#include "pgen/emri_grmhd_tetrad.hpp"
#include "pgen/pgen.hpp"
#include "pgen/emri_srmhd_characteristics.hpp"

#if MPI_PARALLEL_ENABLED
#include <mpi.h>
#endif

namespace {

constexpr std::uint64_t kHeaderBytes = 72;
constexpr std::uint64_t kADMHeaderBytes = 96;
constexpr std::uint32_t kLegacyBinaryVersion = 1;
constexpr std::uint32_t kBinaryVersion = 2;
constexpr std::uint32_t kADMBinaryVersion = 1;
constexpr int kADMFields = 16;
constexpr int kFaceColumns = 7;
constexpr int kEdgeColumns = 7;
constexpr int kVolumeFaceColumns = 6;
constexpr std::array<unsigned char, 16> kMagic = {
    'A', 'E', 'M', 'R', 'I', 'W', 'T', 'B', 'I', 'N', '0', '0', '0', '1', 0, 0};
constexpr std::array<unsigned char, 16> kADMMagic = {
    'A', 'E', 'M', 'R', 'I', 'A', 'D', 'M', 'V', 'O', 'L', '0', '0', '1', 0, 0};
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

Real StageStateTime(Mesh *pm, Driver *pdrive, int stage) {
  Real start = pm->time;
  Real dt = pm->dt;
  if (pdrive->LevelSubcyclingRequested()) {
    start = pdrive->level_subcycling.time;
    dt = pdrive->level_subcycling.dt;
  }
  if (stage <= 0) return start;
  if (pdrive->integrator == "rk1" || pdrive->integrator == "rk2") {
    return start + dt;
  }
  if (pdrive->integrator == "rk3") {
    return start + ((stage == 2) ? 0.5 : 1.0)*dt;
  }
  InnerFatal("fluid worldtube replay encountered an unsupported RK integrator");
}

void ValidateSRCharacteristicReferenceBasis() {
  const Real lorentz = 1.0/std::sqrt(1.0 - (0.1*0.1 + 0.3*0.3 + 0.05*0.05));
  const Real primitive[emri_srmhd::kModes] = {
      4.0, 1.0, -0.1*lorentz, 0.3*lorentz, -0.05*lorentz, 1.8, -1.2};
  const Real expected_speeds[emri_srmhd::kModes] = {
      -0.7848755474024283, -0.6775953662109416, -0.3958024698080614,
      -0.1, 0.1681193457906640, 0.3894059410810621, 0.7189814956447750};
  const Real expected_right[emri_srmhd::kModes][emri_srmhd::kModes] = {
      {0.5782928355320754, 0.0, 0.9412713619603429, 1.0,
       0.9407335001623420, 0.0, 0.7452975389722657},
      {0.1927642785106916, 0.0, 0.3137571206534470, 0.0,
       0.3135778333874469, 0.0, 0.2484325129907554},
      {-0.1374857352758566, -0.0281591855092123, -0.0680845482723261, 0.0,
       0.0688822004851069, -0.0084968436978090, 0.1277011368582980},
      {0.0790192264949366, 0.1591482519396228, -0.0614738525924780, 0.0,
       0.0624394842443209, -0.1046879458929723, -0.1049287522806895},
      {-0.0910428544104706, 0.1800574080985748, 0.0397054032977907, 0.0,
       -0.0393066042700111, -0.1747336765432273, 0.0433738843454690},
      {0.5615454558826656, 0.6723271183537233, -0.0591435123941768, 0.0,
       -0.0680824551096638, 0.4377355357997133, 0.5194032697198726},
      {-0.5288202562124930, 0.6995959333008033, 0.0455596094936173, 0.0,
       0.0431942995194695, 0.8756848637561974, -0.2895413459974056}};
  const Real tolerance = (sizeof(Real) == sizeof(float)) ? 3.0e-4 : 3.0e-9;
  for (int mode = 0; mode < emri_srmhd::kModes; ++mode) {
    Real speed = 0.0;
    Real column[emri_srmhd::kModes];
    if (!emri_srmhd::Eigenvector(
            primitive, 2.5, 4.0/3.0, mode, speed, column) ||
        std::abs(speed - expected_speeds[mode]) > tolerance) {
      InnerFatal("SR characteristic analytic reference speed regression failed");
    }
    Real overlap = 0.0;
    for (int variable = 0; variable < emri_srmhd::kModes; ++variable) {
      overlap += column[variable]*expected_right[variable][mode];
    }
    if (std::abs(std::abs(overlap) - 1.0) > 10.0*tolerance) {
      InnerFatal("SR characteristic analytic eigenvector regression failed");
    }
  }
}

void ValidateGRFaceFrameReference() {
  const Real metric[NSPMETRIC] = {2.0, 0.3, -0.2, 1.5, 0.25, 1.2};
  const Real expected_basis[3][3] = {
      {-0.7281555866966863, 0.1718237643428152, -0.1571558820208675},
      {0.0, 0.8164965809277260, 0.0},
      {0.0, 0.1548574032706278, -0.9291444196237665}};
  const Real expected_dual[3][3] = {
      {-1.373332867686355, 0.0, 0.0},
      {0.2449489742783178, 1.224744871391589, 0.2041241452319315},
      {0.2322861049059417, 0.0, -1.076258952730863}};
  Real basis[3][3], dual[3][3], sqrt_determinant, sqrt_inverse_normal_metric;
  if (!emri_grmhd::BuildFaceFrame(metric, 0, -1, basis, dual, sqrt_determinant,
                                   sqrt_inverse_normal_metric)) {
    InnerFatal("GR characteristic face-frame reference construction failed");
  }
  const Real tolerance = (sizeof(Real) == sizeof(float)) ? 3.0e-5 : 3.0e-12;
  if (std::abs(sqrt_determinant - std::sqrt(3.277)) > tolerance ||
      std::abs(sqrt_inverse_normal_metric - 0.7281555866966863) > tolerance) {
    InnerFatal("GR characteristic face-frame metric regression failed");
  }
  for (int frame = 0; frame < 3; ++frame) {
    for (int component = 0; component < 3; ++component) {
      if (std::abs(basis[frame][component] - expected_basis[frame][component]) >
              tolerance ||
          std::abs(dual[frame][component] - expected_dual[frame][component]) >
              tolerance) {
        InnerFatal("GR characteristic face-frame basis regression failed");
      }
    }
  }
  const Real coordinate[3] = {0.4, -0.2, 0.7};
  Real local[3], recovered[3];
  emri_grmhd::CoordinateToFrame(dual, coordinate, local);
  emri_grmhd::FrameToCoordinate(basis, local, recovered);
  for (int component = 0; component < 3; ++component) {
    if (std::abs(recovered[component] - coordinate[component]) > 10.0*tolerance) {
      InnerFatal("GR characteristic face-frame round-trip regression failed");
    }
  }
  const Real grid_speed = emri_grmhd::OutwardGridSpeed(
      0.83, 0.15, -1, sqrt_inverse_normal_metric);
  if (std::abs(grid_speed + 0.2481926869312688) > tolerance) {
    InnerFatal("GR characteristic shift-speed regression failed");
  }
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

KOKKOS_INLINE_FUNCTION
Real SampleADMField(const DvceArray1D<Real> &slab, int field, int nx, int ny,
                    int nz, Real qx, Real qy, Real qz) {
  int i0 = static_cast<int>(floor(qx));
  int j0 = static_cast<int>(floor(qy));
  int k0 = static_cast<int>(floor(qz));
  i0 = (i0 < 0) ? 0 : ((i0 >= nx - 1) ? nx - 2 : i0);
  j0 = (j0 < 0) ? 0 : ((j0 >= ny - 1) ? ny - 2 : j0);
  k0 = (k0 < 0) ? 0 : ((k0 >= nz - 1) ? nz - 2 : k0);
  const Real fi = fmin(fmax(qx - i0, Real(0.0)), Real(1.0));
  const Real fj = fmin(fmax(qy - j0, Real(0.0)), Real(1.0));
  const Real fk = fmin(fmax(qz - k0, Real(0.0)), Real(1.0));
  const std::size_t field_offset =
      static_cast<std::size_t>(field)*nx*ny*nz;
  Real result = 0.0;
  for (int dk = 0; dk < 2; ++dk) {
    const Real wk = (dk == 0) ? 1.0 - fk : fk;
    for (int dj = 0; dj < 2; ++dj) {
      const Real wj = (dj == 0) ? 1.0 - fj : fj;
      for (int di = 0; di < 2; ++di) {
        const Real wi = (di == 0) ? 1.0 - fi : fi;
        const std::size_t index = field_offset
            + (static_cast<std::size_t>(k0 + dk)*ny + j0 + dj)*nx + i0 + di;
        result += wi*wj*wk*slab(index);
      }
    }
  }
  return result;
}

}  // namespace

EmriInnerWorldtubeReplay::EmriInnerWorldtubeReplay(ParameterInput *pin, Mesh *pm,
                                                   bool is_restart) :
    pmesh_(pm), is_restart_(is_restart) {
  if (pm->pmb_pack == nullptr || pm->pmb_pack->pmhd == nullptr) {
    InnerFatal("inner worldtube replay requires the MHD module");
  }
  if (!pm->three_d || pm->adaptive || pm->multilevel) {
    InnerFatal("inner worldtube replay currently requires a fixed single-level 3D mesh");
  }
  path_ = pin->GetString("emri_worldtube", "file");
  flux_tolerance_ = pin->GetOrAddReal("emri_worldtube", "flux_tolerance", 1.0e-10);
  if (!(std::isfinite(flux_tolerance_) && flux_tolerance_ > 0.0)) {
    InnerFatal("<emri_worldtube>/flux_tolerance must be finite and positive");
  }
  const std::string fluid_boundary =
      pin->GetOrAddString("emri_worldtube", "fluid_boundary", "off");
  if (fluid_boundary == "off") {
    fluid_boundary_enabled_ = false;
  } else if (fluid_boundary == "riemann") {
    fluid_boundary_enabled_ = true;
  } else if (fluid_boundary == "characteristic_sr") {
    fluid_boundary_enabled_ = true;
    characteristic_sr_boundary_ = true;
  } else if (fluid_boundary == "characteristic_gr") {
    fluid_boundary_enabled_ = true;
    characteristic_gr_boundary_ = true;
  } else {
    InnerFatal(
        "<emri_worldtube>/fluid_boundary must be off, riemann, characteristic_sr, "
        "or characteristic_gr");
  }
  characteristic_speed_tolerance_ = pin->GetOrAddReal(
      "emri_worldtube", "characteristic_speed_tolerance", 1.0e-10);
  if (!(std::isfinite(characteristic_speed_tolerance_) &&
        characteristic_speed_tolerance_ >= 0.0)) {
    InnerFatal("characteristic_speed_tolerance must be finite and nonnegative");
  }
  ReadHeaderAndTimes(pin, pm);
  if (pin->DoesBlockExist("emri_adm_replay") &&
      pin->GetOrAddBoolean("emri_adm_replay", "enabled", false)) {
    ReadADMVolume(pin, pm);
    pm->pmb_pack->padm->SetADMVariables = &EmriInnerWorldtubeSetADMVariables;
  }
  if (fluid_boundary_enabled_) {
    mhd::MHD *pmhd = pm->pmb_pack->pmhd;
    if (!has_cell_centered_magnetic_state_) {
      InnerFatal(
          "fluid_boundary=riemann requires trailing bcc1,bcc2,bcc3 state data");
    }
    if (pmhd->nmhd != 5 || !pmhd->peos->eos_data.is_ideal) {
      InnerFatal("fluid worldtube replay currently requires ideal-gas MHD");
    }
    if (characteristic_sr_boundary_ &&
        !pm->pmb_pack->pcoord->is_special_relativistic) {
      InnerFatal("fluid_boundary=characteristic_sr requires special-relativistic MHD");
    }
    if (characteristic_gr_boundary_ &&
        (pm->pmb_pack->pdyngr == nullptr || pm->pmb_pack->padm == nullptr ||
         pm->pmb_pack->pdyngr->eos_policy != DynGRMHD_EOS::eos_ideal)) {
      InnerFatal(
          "fluid_boundary=characteristic_gr requires ideal-gas DynGRMHD with ADM data");
    }
    if (characteristic_gr_boundary_ &&
        pin->GetOrAddString("emri_worldtube", "fluid_state_frame", "unverified") !=
            "inner_coordinate") {
      InnerFatal(
          "fluid_boundary=characteristic_gr requires fluid_state_frame=inner_coordinate; "
          "prepare-inner does not perform the global-to-inner coordinate transform");
    }
    if (characteristic_sr_boundary_ || characteristic_gr_boundary_) {
      ValidateSRCharacteristicReferenceBasis();
      ValidateGRFaceFrameReference();
    }
  }
  BuildBoundaryTopology(pm);
  LoadInterval(interval_);
  if (adm_volume_enabled_) {
    SetADMVariables(pm->pmb_pack);
    if (pm->pmb_pack->pcoord->coord_data.bh_excise) {
      pm->pmb_pack->pcoord->UpdateExcisionMasks();
    }
  }
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
  if ((version != kLegacyBinaryVersion && version != kBinaryVersion) ||
      cells_per_edge_ < 1 || nvar_ < 1 || nt_ < 2 ||
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
  has_initial_volume_flux_ = (version >= kBinaryVersion);
  if (has_initial_volume_flux_) {
    const std::uint64_t volume_face_count = CheckedProduct(
        {cells + 1, cells, cells}, "initial volume face-flux array");
    for (int component = 0; component < 3; ++component) {
      initial_volume_flux_offsets_[component] = cursor;
      cursor += 8*volume_face_count;
    }
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

void EmriInnerWorldtubeReplay::ReadADMVolume(ParameterInput *pin, Mesh *pm) {
  if (pm->pmb_pack->padm == nullptr || pm->pmb_pack->pdyngr == nullptr) {
    InnerFatal("ADM volume replay requires prescribed DynGRMHD");
  }
  if (!pm->pmb_pack->padm->is_dynamic) {
    InnerFatal("ADM volume replay requires <adm> dynamic=true for RK-stage updates");
  }
  if (pin->GetOrAddBoolean("problem", "user_hist", false)) {
    InnerFatal(
        "ADM volume replay does not yet support EMRI force history because that "
        "diagnostic still evaluates the analytic metric off-grid");
  }
  adm_path_ = pin->GetString("emri_adm_replay", "file");
  std::ifstream stream(adm_path_, std::ios::binary | std::ios::ate);
  if (!stream.is_open()) InnerFatal("could not open ADM volume replay " + adm_path_);
  const std::streamoff signed_size = stream.tellg();
  if (signed_size < static_cast<std::streamoff>(kADMHeaderBytes)) {
    InnerFatal("ADM volume replay is shorter than its header");
  }
  const std::uint64_t file_size = static_cast<std::uint64_t>(signed_size);
  stream.seekg(0);
  std::array<unsigned char, kADMHeaderBytes> header{};
  stream.read(reinterpret_cast<char*>(header.data()), header.size());
  if (stream.gcount() != static_cast<std::streamsize>(header.size()) ||
      !std::equal(kADMMagic.begin(), kADMMagic.end(), header.begin())) {
    InnerFatal("ADM volume replay header magic is invalid");
  }
  const std::uint32_t version = ReadLE32(header.data() + 16);
  const int nt = static_cast<int>(ReadLE32(header.data() + 20));
  adm_nvar_ = static_cast<int>(ReadLE32(header.data() + 24));
  adm_nx_ = static_cast<int>(ReadLE32(header.data() + 28));
  adm_ny_ = static_cast<int>(ReadLE32(header.data() + 32));
  adm_nz_ = static_cast<int>(ReadLE32(header.data() + 36));
  for (int axis = 0; axis < 3; ++axis) {
    adm_lower_[axis] = static_cast<Real>(ReadLEDouble(header.data() + 40 + 8*axis));
    adm_spacing_[axis] =
        static_cast<Real>(ReadLEDouble(header.data() + 64 + 8*axis));
  }
  const std::uint64_t stored_crc = ReadLE64(header.data() + 88);
  if (version != kADMBinaryVersion || nt != nt_ || adm_nvar_ != kADMFields ||
      adm_nx_ < 2 || adm_ny_ < 2 || adm_nz_ < 2 ||
      stored_crc > std::numeric_limits<std::uint32_t>::max()) {
    InnerFatal("ADM volume replay dimensions or version are unsupported");
  }
  for (int axis = 0; axis < 3; ++axis) {
    if (!(std::isfinite(adm_lower_[axis]) &&
          std::isfinite(adm_spacing_[axis]) && adm_spacing_[axis] > 0.0)) {
      InnerFatal("ADM volume replay grid geometry is invalid");
    }
  }
  adm_slab_values_ = CheckedProduct(
      {static_cast<std::uint64_t>(adm_nvar_),
       static_cast<std::uint64_t>(adm_nx_),
       static_cast<std::uint64_t>(adm_ny_),
       static_cast<std::uint64_t>(adm_nz_)},
      "ADM volume slab");
  const std::uint64_t payload_values = CheckedProduct(
      {static_cast<std::uint64_t>(nt), adm_slab_values_},
      "ADM volume payload");
  const std::uint64_t payload_bytes = CheckedProduct(
      {8, static_cast<std::uint64_t>(nt) + payload_values},
      "ADM volume payload bytes");
  if (payload_bytes > std::numeric_limits<std::uint64_t>::max() - kADMHeaderBytes ||
      kADMHeaderBytes + payload_bytes != file_size) {
    InnerFatal("ADM volume replay file size does not match its header");
  }
  stream.seekg(kADMHeaderBytes);
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
    InnerFatal("ADM volume replay payload CRC32 mismatch");
  }
  const std::vector<Real> adm_times = ReadDoublesAt(adm_path_, kADMHeaderBytes, nt);
  for (int index = 0; index < nt_; ++index) {
    if (!NearlyEqual(adm_times[index], times_[index], times_.back() - times_.front())) {
      InnerFatal("ADM volume and fluid worldtube time tables differ");
    }
  }
  adm_data_offset_ = kADMHeaderBytes + 8*static_cast<std::uint64_t>(nt);
  adm_left_ = DvceArray1D<Real>("emri_adm_left", adm_slab_values_);
  adm_right_ = DvceArray1D<Real>("emri_adm_right", adm_slab_values_);
  adm_interval_ = interval_;
  adm_volume_enabled_ = true;

  const RegionIndcs &indcs = pm->mb_indcs;
  auto &sizes = pm->pmb_pack->pmb->mb_size;
  const Real grid_upper[3] = {
      adm_lower_[0] + (adm_nx_ - 1)*adm_spacing_[0],
      adm_lower_[1] + (adm_ny_ - 1)*adm_spacing_[1],
      adm_lower_[2] + (adm_nz_ - 1)*adm_spacing_[2]};
  for (int m = 0; m < pm->pmb_pack->nmb_thispack; ++m) {
    const RegionSize &size = sizes.h_view(m);
    const Real minimum[3] = {
        size.x1min + (0.5 - indcs.ng)*size.dx1,
        size.x2min + (0.5 - indcs.ng)*size.dx2,
        size.x3min + (0.5 - indcs.ng)*size.dx3};
    const Real maximum[3] = {
        size.x1max + (indcs.ng - 0.5)*size.dx1,
        size.x2max + (indcs.ng - 0.5)*size.dx2,
        size.x3max + (indcs.ng - 0.5)*size.dx3};
    for (int axis = 0; axis < 3; ++axis) {
      const Real tolerance = 256.0*std::numeric_limits<Real>::epsilon()
          *std::max({std::abs(adm_lower_[axis]), std::abs(grid_upper[axis]),
                     std::abs(minimum[axis]), std::abs(maximum[axis]), Real(1.0)});
      if (minimum[axis] < adm_lower_[axis] - tolerance ||
          maximum[axis] > grid_upper[axis] + tolerance) {
        InnerFatal("ADM volume replay does not cover every MeshBlock ghost center");
      }
    }
  }
  LoadADMInterval(adm_interval_);
  if (global_variable::my_rank == 0) {
    std::cout << "EMRI numerical ADM volume replay: grid="
              << adm_nx_ << "x" << adm_ny_ << "x" << adm_nz_
              << ", samples=" << nt_ << std::endl;
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

  if (has_initial_volume_flux_) {
    const Real domain_minimum[3] = {pm->mesh_size.x1min, pm->mesh_size.x2min,
                                    pm->mesh_size.x3min};
    const std::size_t component_size =
        static_cast<std::size_t>(cells + 1)*cells*cells;
    if (3*component_size > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      InnerFatal("initial volume face-flux table exceeds 32-bit device indexing");
    }
    for (int m = 0; m < pack->nmb_thispack; ++m) {
      const RegionSize &size = sizes.h_view(m);
      const Real minimum[3] = {size.x1min, size.x2min, size.x3min};
      const int axis_cells[3] = {indcs.nx1, indcs.nx2, indcs.nx3};
      const int global_lower[3] = {
          static_cast<int>(std::llround((minimum[0] - domain_minimum[0])/dx_)),
          static_cast<int>(std::llround((minimum[1] - domain_minimum[1])/dx_)),
          static_cast<int>(std::llround((minimum[2] - domain_minimum[2])/dx_))};
      for (int axis = 0; axis < 3; ++axis) {
        if (global_lower[axis] < 0 || global_lower[axis] + axis_cells[axis] > cells) {
          InnerFatal("inner MeshBlock lies outside the initial volume-flux table");
        }
      }
      for (int k = 0; k < indcs.nx3; ++k) {
        const int global_k = global_lower[2] + k;
        for (int j = 0; j < indcs.nx2; ++j) {
          const int global_j = global_lower[1] + j;
          for (int i = 0; i <= indcs.nx1; ++i) {
            const int global_i = global_lower[0] + i;
            const std::size_t input =
                (static_cast<std::size_t>(global_k)*cells + global_j)
                *(cells + 1) + global_i;
            host_volume_faces_.push_back(VolumeFaceRecord{
                static_cast<int>(input), m, indcs.ks + k, indcs.js + j,
                indcs.is + i, 0});
          }
        }
      }
      for (int k = 0; k < indcs.nx3; ++k) {
        const int global_k = global_lower[2] + k;
        for (int j = 0; j <= indcs.nx2; ++j) {
          const int global_j = global_lower[1] + j;
          for (int i = 0; i < indcs.nx1; ++i) {
            const int global_i = global_lower[0] + i;
            const std::size_t input = component_size
                + (static_cast<std::size_t>(global_k)*(cells + 1) + global_j)
                *cells + global_i;
            host_volume_faces_.push_back(VolumeFaceRecord{
                static_cast<int>(input), m, indcs.ks + k, indcs.js + j,
                indcs.is + i, 1});
          }
        }
      }
      for (int k = 0; k <= indcs.nx3; ++k) {
        const int global_k = global_lower[2] + k;
        for (int j = 0; j < indcs.nx2; ++j) {
          const int global_j = global_lower[1] + j;
          for (int i = 0; i < indcs.nx1; ++i) {
            const int global_i = global_lower[0] + i;
            const std::size_t input = 2*component_size
                + (static_cast<std::size_t>(global_k)*cells + global_j)
                *cells + global_i;
            host_volume_faces_.push_back(VolumeFaceRecord{
                static_cast<int>(input), m, indcs.ks + k, indcs.js + j,
                indcs.is + i, 2});
          }
        }
      }
    }
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
  volume_face_records_ = DvceArray2D<int>(
      "emri_inner_volume_faces", std::max<int>(host_volume_faces_.size(), 1),
      kVolumeFaceColumns);
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
  auto host_volume_face_view = Kokkos::create_mirror_view(volume_face_records_);
  for (std::size_t index = 0; index < host_volume_faces_.size(); ++index) {
    const VolumeFaceRecord &record = host_volume_faces_[index];
    const int values[kVolumeFaceColumns] = {
        record.input, record.m, record.k, record.j, record.i, record.component};
    for (int column = 0; column < kVolumeFaceColumns; ++column) {
      host_volume_face_view(index, column) = values[column];
    }
  }
  Kokkos::deep_copy(volume_face_records_, host_volume_face_view);

  const std::size_t face_cells = static_cast<std::size_t>(6)*cells*cells;
  const std::size_t face_edges = static_cast<std::size_t>(12)*cells*(cells + 1);
  state_left_ = DvceArray1D<Real>("emri_inner_state_left", face_cells*nvar_);
  state_right_ = DvceArray1D<Real>("emri_inner_state_right", face_cells*nvar_);
  flux_left_ = DvceArray1D<Real>("emri_inner_flux_left", face_cells);
  flux_right_ = DvceArray1D<Real>("emri_inner_flux_right", face_cells);
  interval_emf_ = DvceArray1D<Real>("emri_inner_interval_emf", face_edges);
  const std::size_t volume_faces = static_cast<std::size_t>(3)*(cells + 1)*cells*cells;
  initial_volume_flux_ = DvceArray1D<Real>(
      "emri_inner_initial_volume_flux", std::max<std::size_t>(volume_faces, 1));
  if (has_initial_volume_flux_) {
    auto host_volume_flux = Kokkos::create_mirror_view(initial_volume_flux_);
    const std::uint64_t component_count =
        static_cast<std::uint64_t>(cells + 1)*cells*cells;
    for (int component = 0; component < 3; ++component) {
      const std::vector<Real> values = ReadDoublesAt(
          path_, initial_volume_flux_offsets_[component], component_count);
      const std::size_t base = static_cast<std::size_t>(component)*component_count;
      for (std::size_t index = 0; index < values.size(); ++index) {
        host_volume_flux(base + index) = values[index];
      }
    }
    Kokkos::deep_copy(initial_volume_flux_, host_volume_flux);
  }
  characteristic_diagnostics_ = DvceArray1D<int>(
      "emri_inner_characteristic_diagnostics", 2);
  Kokkos::deep_copy(characteristic_diagnostics_, 0);
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

void EmriInnerWorldtubeReplay::LoadADMInterval(int interval) {
  if (!adm_volume_enabled_ || interval < 0 || interval >= nt_ - 1) {
    InnerFatal("requested ADM replay interval is out of range");
  }
  const std::uint64_t left_offset = adm_data_offset_
      + 8*adm_slab_values_*static_cast<std::uint64_t>(interval);
  const std::uint64_t right_offset = left_offset + 8*adm_slab_values_;
  const std::vector<Real> left = ReadDoublesAt(
      adm_path_, left_offset, adm_slab_values_);
  const std::vector<Real> right = ReadDoublesAt(
      adm_path_, right_offset, adm_slab_values_);
  auto host_left = Kokkos::create_mirror_view(adm_left_);
  auto host_right = Kokkos::create_mirror_view(adm_right_);
  for (std::size_t index = 0; index < left.size(); ++index) {
    host_left(index) = left[index];
    host_right(index) = right[index];
  }
  Kokkos::deep_copy(adm_left_, host_left);
  Kokkos::deep_copy(adm_right_, host_right);
  adm_interval_ = interval;
}

void EmriInnerWorldtubeReplay::SetADMVariables(MeshBlockPack *pack) {
  if (!adm_volume_enabled_ || pack == nullptr || pack->padm == nullptr) {
    InnerFatal("invalid numerical ADM replay callback state");
  }
  if (adm_interval_ != interval_) {
    InnerFatal("fluid and ADM replay interval cursors differ");
  }
  const Real table_time = pack->pmesh->time + time_offset_;
  const Real left_time = times_[adm_interval_];
  const Real right_time = times_[adm_interval_ + 1];
  Real fraction = (table_time - left_time)/(right_time - left_time);
  const Real tolerance = 512.0*std::numeric_limits<Real>::epsilon()
      *std::max({std::abs(table_time), std::abs(left_time),
                 std::abs(right_time), Real(1.0)});
  if (fraction < 0.0 && table_time >= left_time - tolerance) fraction = 0.0;
  if (fraction > 1.0 && table_time <= right_time + tolerance) fraction = 1.0;
  if (!(fraction >= 0.0 && fraction <= 1.0)) {
    InnerFatal("ADM metric time lies outside the loaded replay interval");
  }

  auto left = adm_left_;
  auto right = adm_right_;
  auto &adm_vars = pack->padm->adm;
  auto &sizes = pack->pmb->mb_size;
  const RegionIndcs &indcs = pack->pmesh->mb_indcs;
  const int is = indcs.is;
  const int js = indcs.js;
  const int ks = indcs.ks;
  const int n1 = indcs.nx1 + 2*indcs.ng;
  const int n2 = indcs.nx2 + 2*indcs.ng;
  const int n3 = indcs.nx3 + 2*indcs.ng;
  const int nx = adm_nx_;
  const int ny = adm_ny_;
  const int nz = adm_nz_;
  const Real lower_x = adm_lower_[0];
  const Real lower_y = adm_lower_[1];
  const Real lower_z = adm_lower_[2];
  const Real inverse_dx = 1.0/adm_spacing_[0];
  const Real inverse_dy = 1.0/adm_spacing_[1];
  const Real inverse_dz = 1.0/adm_spacing_[2];
  const int nmb = pack->nmb_thispack;
  par_for("emri_numerical_adm_replay", DevExeSpace(), 0, nmb - 1,
          0, n3 - 1, 0, n2 - 1, 0, n1 - 1,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    const Real x = CellCenterX(
        i - is, indcs.nx1, sizes.d_view(m).x1min, sizes.d_view(m).x1max);
    const Real y = CellCenterX(
        j - js, indcs.nx2, sizes.d_view(m).x2min, sizes.d_view(m).x2max);
    const Real z = CellCenterX(
        k - ks, indcs.nx3, sizes.d_view(m).x3min, sizes.d_view(m).x3max);
    const Real qx = (x - lower_x)*inverse_dx;
    const Real qy = (y - lower_y)*inverse_dy;
    const Real qz = (z - lower_z)*inverse_dz;
    Real values[kADMFields];
    for (int field = 0; field < kADMFields; ++field) {
      const Real value_left = SampleADMField(left, field, nx, ny, nz, qx, qy, qz);
      const Real value_right = SampleADMField(
          right, field, nx, ny, nz, qx, qy, qz);
      values[field] = (1.0 - fraction)*value_left + fraction*value_right;
    }
    const Real gxx = values[4];
    const Real gxy = values[5];
    const Real gxz = values[6];
    const Real gyy = values[7];
    const Real gyz = values[8];
    const Real gzz = values[9];
    const Real determinant = adm::SpatialDet(gxx, gxy, gxz, gyy, gyz, gzz);
    const Real minor2 = gxx*gyy - gxy*gxy;
    if (!(gxx > 0.0 && minor2 > 0.0 && determinant > 0.0)) {
      Kokkos::abort("interpolated numerical ADM spatial metric is not positive");
    }
    Real inverse[3][3];
    adm::SpatialInv(1.0/determinant, gxx, gxy, gxz, gyy, gyz, gzz,
        &inverse[0][0], &inverse[0][1], &inverse[0][2],
        &inverse[1][1], &inverse[1][2], &inverse[2][2]);
    inverse[1][0] = inverse[0][1];
    inverse[2][0] = inverse[0][2];
    inverse[2][1] = inverse[1][2];
    const Real beta_lower[3] = {values[1], values[2], values[3]};
    Real beta[3] = {0.0, 0.0, 0.0};
    Real beta_squared = 0.0;
    for (int a = 0; a < 3; ++a) {
      for (int b = 0; b < 3; ++b) beta[a] += inverse[a][b]*beta_lower[b];
      beta_squared += beta[a]*beta_lower[a];
    }
    const Real lapse_squared = beta_squared - values[0];
    if (!(lapse_squared > 0.0)) {
      Kokkos::abort("interpolated numerical ADM lapse is not positive");
    }
    const Real gamma[6] = {gxx, gxy, gxz, gyy, gyz, gzz};
    const int row[6] = {0, 0, 0, 1, 1, 2};
    const int column[6] = {0, 1, 2, 1, 2, 2};
    for (int a = 0; a < 3; ++a) adm_vars.beta_u(m, a, k, j, i) = beta[a];
    for (int component = 0; component < 6; ++component) {
      adm_vars.g_dd(m, row[component], column[component], k, j, i) =
          gamma[component];
      adm_vars.vK_dd(m, row[component], column[component], k, j, i) =
          values[10 + component];
    }
    adm_vars.alpha(m, k, j, i) = sqrt(lapse_squared);
    adm_vars.psi4(m, k, j, i) = cbrt(determinant);
  });
}

void EmriInnerWorldtubeReplay::SetInitialNormalFlux(Mesh *pm) {
  auto b1 = pm->pmb_pack->pmhd->b0.x1f;
  auto b2 = pm->pmb_pack->pmhd->b0.x2f;
  auto b3 = pm->pmb_pack->pmhd->b0.x3f;
  const Real inverse_area = 1.0/(dx_*dx_);
  if (has_initial_volume_flux_) {
    auto records = volume_face_records_;
    auto expected = initial_volume_flux_;
    const int count = static_cast<int>(host_volume_faces_.size());
    if (count > 0) {
      par_for("emri_inner_set_volume_flux", DevExeSpace(), 0, count - 1,
      KOKKOS_LAMBDA(int record_index) {
        const int input = records(record_index, 0);
        const int m = records(record_index, 1);
        const int k = records(record_index, 2);
        const int j = records(record_index, 3);
        const int i = records(record_index, 4);
        const int component = records(record_index, 5);
        const Real field = expected(input)*inverse_area;
        if (component == 0) b1(m, k, j, i) = field;
        if (component == 1) b2(m, k, j, i) = field;
        if (component == 2) b3(m, k, j, i) = field;
      });
    }
  } else {
    auto records = face_records_;
    auto expected = flux_left_;
    const int cells = cells_per_edge_;
    const int count = static_cast<int>(host_faces_.size());
    if (count > 0) {
      par_for("emri_inner_set_boundary_flux", DevExeSpace(), 0, count - 1,
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

  MeshBlockPack *pack = pm->pmb_pack;
  auto bcc = pack->pmhd->bcc0;
  const RegionIndcs &indcs = pm->mb_indcs;
  const int nmb = pack->nmb_thispack;
  par_for("emri_inner_set_cell_field", DevExeSpace(), 0, nmb - 1,
          indcs.ks, indcs.ke, indcs.js, indcs.je, indcs.is, indcs.ie,
  KOKKOS_LAMBDA(int m, int k, int j, int i) {
    bcc(m, IBX, k, j, i) = 0.5*(b1(m, k, j, i) + b1(m, k, j, i + 1));
    bcc(m, IBY, k, j, i) = 0.5*(b2(m, k, j, i) + b2(m, k, j + 1, i));
    bcc(m, IBZ, k, j, i) = 0.5*(b3(m, k, j, i) + b3(m, k + 1, j, i));
  });
  if (pack->pdyngr == nullptr) {
    InnerFatal("inner volume-flux initialization currently requires DynGRMHD");
  }
  pack->pdyngr->PrimToConInit(
      indcs.is, indcs.ie, indcs.js, indcs.je, indcs.ks, indcs.ke);
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

void EmriInnerWorldtubeReplay::ApplyPrimitiveBoundary(Mesh *pm, Driver *pdrive,
                                                       int stage) {
  if (!fluid_boundary_enabled_ || exhausted_) return;
  const Real data_time = StageStateTime(pm, pdrive, stage) + time_offset_;
  const Real left_time = times_[interval_];
  const Real right_time = times_[interval_ + 1];
  const Real time_scale = times_.back() - times_.front();
  if (data_time < left_time && !NearlyEqual(data_time, left_time, time_scale)) {
    InnerFatal("fluid replay stage time precedes the loaded worldtube interval");
  }
  if (data_time > right_time && !NearlyEqual(data_time, right_time, time_scale)) {
    InnerFatal("fluid replay stage time exceeds the loaded worldtube interval");
  }
  const Real theta = std::clamp(
      (data_time - left_time)/(right_time - left_time), Real(0.0), Real(1.0));
  mhd::MHD *pmhd = pm->pmb_pack->pmhd;
  auto records = face_records_;
  auto left = state_left_;
  auto right = state_right_;
  auto left_flux = flux_left_;
  auto right_flux = flux_right_;
  auto w0 = pmhd->w0;
  auto bcc0 = pmhd->bcc0;
  const EOS_Data eos = pmhd->peos->eos_data;
  auto characteristic_diagnostics = characteristic_diagnostics_;
  const int local_faces = static_cast<int>(host_faces_.size());
  const int cells = cells_per_edge_;
  const int nvar = nvar_;
  const int fluid_nvar = pmhd->nmhd + pmhd->nscalars;
  const int nghost = pm->mb_indcs.ng;
  const Real inverse_area = 1.0/(dx_*dx_);
  const bool characteristic_sr = characteristic_sr_boundary_;
  const bool characteristic_gr = characteristic_gr_boundary_;
  const bool characteristic = characteristic_sr || characteristic_gr;
  const Real speed_tolerance = characteristic_speed_tolerance_;
  adm::ADM::ADM_vars adm_fields;
  if (characteristic_gr) adm_fields = pm->pmb_pack->padm->adm;
  if (local_faces == 0) return;
  par_for("emri_inner_fluid_boundary", DevExeSpace(), 0, local_faces - 1,
  KOKKOS_LAMBDA(int record_index) {
    const int face = records(record_index, 0);
    const int m = records(record_index, 1);
    const int k = records(record_index, 2);
    const int j = records(record_index, 3);
    const int i = records(record_index, 4);
    const int v = records(record_index, 5);
    const int u = records(record_index, 6);
    const int normal_axis = face/2;
    const int normal_sign = (face % 2 == 0) ? -1 : 1;
    const int face_cell = (face*cells + v)*cells + u;
    Real exterior_b[3];
    for (int component = 0; component < 3; ++component) {
      const int variable = fluid_nvar + component;
      const int input = ((face*nvar + variable)*cells + v)*cells + u;
      exterior_b[component] = (1.0 - theta)*left(input) + theta*right(input);
    }
    const Real outward_flux =
        (1.0 - theta)*left_flux(face_cell) + theta*right_flux(face_cell);
    exterior_b[normal_axis] = normal_sign*outward_flux*inverse_area;
    Real interior_q[emri_srmhd::kModes] = {0.0};
    Real boundary_q[emri_srmhd::kModes] = {0.0};
    Real basis[3][3] = {{0.0}};
    Real normal_field = 0.0;
    Real sqrt_determinant = 1.0;
    bool projected = false;
    if (characteristic) {
      Real metric[NSPMETRIC] = {1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
      Real shift[3] = {0.0, 0.0, 0.0};
      Real lapse = 1.0;
      if (characteristic_gr) {
        if (normal_axis == 0) {
          adm::Face1Metric(m, k, j, i + (normal_sign > 0), adm_fields.g_dd,
                           adm_fields.beta_u, adm_fields.alpha, metric, shift, lapse);
        } else if (normal_axis == 1) {
          adm::Face2Metric(m, k, j + (normal_sign > 0), i, adm_fields.g_dd,
                           adm_fields.beta_u, adm_fields.alpha, metric, shift, lapse);
        } else {
          adm::Face3Metric(m, k + (normal_sign > 0), j, i, adm_fields.g_dd,
                           adm_fields.beta_u, adm_fields.alpha, metric, shift, lapse);
        }
      }
      Real dual[3][3];
      Real sqrt_inverse_normal_metric = 1.0;
      bool frame_valid = emri_grmhd::BuildFaceFrame(
          metric, normal_axis, normal_sign, basis, dual, sqrt_determinant,
          sqrt_inverse_normal_metric);
      Real outward_grid_speed = 0.0;
      if (characteristic_gr && frame_valid) {
        frame_valid = isfinite(lapse) && lapse > 0.0;
        if (frame_valid) {
          outward_grid_speed = emri_grmhd::OutwardGridSpeed(
              lapse, shift[normal_axis], normal_sign, sqrt_inverse_normal_metric);
          frame_valid = isfinite(outward_grid_speed) && fabs(outward_grid_speed) < 1.0;
        }
      }
      if (frame_valid) {
        Real exterior_q[emri_srmhd::kModes];
        Real exterior_u[3], exterior_b_undensitized[3];
        for (int component = 0; component < 3; ++component) {
          const int input = ((face*nvar + IVX + component)*cells + v)*cells + u;
          exterior_u[component] = (1.0 - theta)*left(input) + theta*right(input);
          exterior_b_undensitized[component] = exterior_b[component]/sqrt_determinant;
        }
        Real exterior_u_frame[3], exterior_b_frame[3];
        emri_grmhd::CoordinateToFrame(dual, exterior_u, exterior_u_frame);
        emri_grmhd::CoordinateToFrame(
            dual, exterior_b_undensitized, exterior_b_frame);
        const int density_input = ((face*nvar + IDN)*cells + v)*cells + u;
        const int energy_input = ((face*nvar + IEN)*cells + v)*cells + u;
        exterior_q[0] = (1.0 - theta)*left(density_input) + theta*right(density_input);
        const Real exterior_thermodynamic =
            (1.0 - theta)*left(energy_input) + theta*right(energy_input);
        exterior_q[1] = characteristic_gr ? exterior_thermodynamic
                                          : eos.IdealGasPressure(exterior_thermodynamic);
        exterior_q[2] = exterior_u_frame[0];
        exterior_q[3] = exterior_u_frame[1];
        exterior_q[4] = exterior_u_frame[2];
        exterior_q[5] = exterior_b_frame[1];
        exterior_q[6] = exterior_b_frame[2];
        normal_field = exterior_b_frame[0];

        Real interior_u[3], interior_b_undensitized[3];
        for (int component = 0; component < 3; ++component) {
          interior_u[component] = w0(m, IVX + component, k, j, i);
          interior_b_undensitized[component] =
              bcc0(m, component, k, j, i)/sqrt_determinant;
        }
        Real interior_u_frame[3], interior_b_frame[3];
        emri_grmhd::CoordinateToFrame(dual, interior_u, interior_u_frame);
        emri_grmhd::CoordinateToFrame(
            dual, interior_b_undensitized, interior_b_frame);
        interior_q[0] = w0(m, IDN, k, j, i);
        interior_q[1] = characteristic_gr ? w0(m, IEN, k, j, i)
                                          : eos.IdealGasPressure(w0(m, IEN, k, j, i));
        interior_q[2] = interior_u_frame[0];
        interior_q[3] = interior_u_frame[1];
        interior_q[4] = interior_u_frame[2];
        interior_q[5] = interior_b_frame[1];
        interior_q[6] = interior_b_frame[2];

        int inward[3] = {i, j, k};
        inward[normal_axis] -= normal_sign;
        Real inward_q[emri_srmhd::kModes];
        Real inward_u[3], inward_b_undensitized[3];
        for (int component = 0; component < 3; ++component) {
          inward_u[component] =
              w0(m, IVX + component, inward[2], inward[1], inward[0]);
          inward_b_undensitized[component] =
              bcc0(m, component, inward[2], inward[1], inward[0])/sqrt_determinant;
        }
        Real inward_u_frame[3], inward_b_frame[3];
        emri_grmhd::CoordinateToFrame(dual, inward_u, inward_u_frame);
        emri_grmhd::CoordinateToFrame(dual, inward_b_undensitized, inward_b_frame);
        inward_q[0] = w0(m, IDN, inward[2], inward[1], inward[0]);
        inward_q[1] = characteristic_gr
            ? w0(m, IEN, inward[2], inward[1], inward[0])
            : eos.IdealGasPressure(w0(m, IEN, inward[2], inward[1], inward[0]));
        inward_q[2] = inward_u_frame[0];
        inward_q[3] = inward_u_frame[1];
        inward_q[4] = inward_u_frame[2];
        inward_q[5] = inward_b_frame[1];
        inward_q[6] = inward_b_frame[2];
        Real predicted_q[emri_srmhd::kModes];
        for (int variable = 0; variable < emri_srmhd::kModes; ++variable) {
          // Retain the outgoing one-sided physical-frame gradient.  Using the active
          // boundary cell itself would flatten every outgoing characteristic mode.
          predicted_q[variable] = 2.0*interior_q[variable] - inward_q[variable];
        }
        projected = emri_srmhd::ProjectIncoming(
            predicted_q, exterior_q, normal_field, eos.gamma, speed_tolerance,
            outward_grid_speed, fmax(eos.dfloor, Real(1.0e-14)),
            fmax(eos.pfloor, Real(1.0e-14)), boundary_q);
      }
      Kokkos::atomic_increment(&characteristic_diagnostics(projected ? 0 : 1));
    }
    for (int layer = 1; layer <= nghost; ++layer) {
      int ghost[3] = {i, j, k};
      ghost[normal_axis] += normal_sign*layer;
      for (int variable = 0; variable < fluid_nvar; ++variable) {
        const int input = ((face*nvar + variable)*cells + v)*cells + u;
        w0(m, variable, ghost[2], ghost[1], ghost[0]) =
            (1.0 - theta)*left(input) + theta*right(input);
      }
      for (int component = 0; component < 3; ++component) {
        bcc0(m, component, ghost[2], ghost[1], ghost[0]) = exterior_b[component];
      }
      if (projected) {
        Real layer_q[emri_srmhd::kModes];
        for (int variable = 0; variable < emri_srmhd::kModes; ++variable) {
          layer_q[variable] = interior_q[variable]
                              + layer*(boundary_q[variable] - interior_q[variable]);
        }
        layer_q[0] = fmax(layer_q[0], fmax(eos.dfloor, Real(1.0e-14)));
        layer_q[1] = fmax(layer_q[1], fmax(eos.pfloor, Real(1.0e-14)));
        w0(m, IDN, ghost[2], ghost[1], ghost[0]) = layer_q[0];
        w0(m, IEN, ghost[2], ghost[1], ghost[0]) = characteristic_gr
            ? layer_q[1] : layer_q[1]/(eos.gamma - 1.0);
        const Real velocity_frame[3] = {layer_q[2], layer_q[3], layer_q[4]};
        Real velocity_coordinate[3];
        emri_grmhd::FrameToCoordinate(basis, velocity_frame, velocity_coordinate);
        const Real magnetic_frame[3] = {normal_field, layer_q[5], layer_q[6]};
        Real magnetic_coordinate[3];
        emri_grmhd::FrameToCoordinate(basis, magnetic_frame, magnetic_coordinate);
        for (int component = 0; component < 3; ++component) {
          w0(m, IVX + component, ghost[2], ghost[1], ghost[0]) =
              velocity_coordinate[component];
          bcc0(m, component, ghost[2], ghost[1], ghost[0]) =
              sqrt_determinant*magnetic_coordinate[component];
        }
        bcc0(m, normal_axis, ghost[2], ghost[1], ghost[0]) =
            normal_sign*outward_flux*inverse_area;
      }
    }
  });
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
    if (adm_volume_enabled_) LoadADMInterval(interval_);
  } else {
    exhausted_ = true;
    if ((characteristic_sr_boundary_ || characteristic_gr_boundary_) &&
        global_variable::my_rank == 0) {
      auto diagnostics = Kokkos::create_mirror_view_and_copy(
          HostMemSpace(), characteristic_diagnostics_);
      std::cout << "EMRI inner characteristic "
                << (characteristic_gr_boundary_ ? "GR" : "SR")
                << " boundary: projections="
                << diagnostics(0) << ", fallbacks=" << diagnostics(1) << std::endl;
    }
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

void EmriInnerWorldtubeApplyPrimitiveBoundary(Mesh *pm, Driver *pdrive, int stage) {
  if (pm == nullptr || pm->pgen == nullptr ||
      pm->pgen->emri_inner_worldtube_ == nullptr) {
    InnerFatal("invalid inner-worldtube primitive-boundary callback state");
  }
  pm->pgen->emri_inner_worldtube_->ApplyPrimitiveBoundary(pm, pdrive, stage);
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

void EmriInnerWorldtubeSetADMVariables(MeshBlockPack *pack) {
  if (pack == nullptr || pack->pmesh == nullptr || pack->pmesh->pgen == nullptr ||
      pack->pmesh->pgen->emri_inner_worldtube_ == nullptr ||
      !pack->pmesh->pgen->emri_inner_worldtube_->ADMVolumeEnabled()) {
    InnerFatal("invalid inner-worldtube numerical ADM callback state");
  }
  pack->pmesh->pgen->emri_inner_worldtube_->SetADMVariables(pack);
}
