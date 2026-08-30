#ifndef PGEN_EMRI_INNER_WORLDTUBE_HPP_
#define PGEN_EMRI_INNER_WORLDTUBE_HPP_
//========================================================================================
// AthenaK astrophysical plasma code
// Copyright(C) 2026 AthenaK collaboration
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file emri_inner_worldtube.hpp
//! \brief Bounded-slab magnetic replay for a cubical EMRI inner worldtube.

#include <array>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

#include "athena.hpp"

class Driver;
class Mesh;
class MeshBlockPack;
class ParameterInput;

class EmriInnerWorldtubeReplay {
 public:
  EmriInnerWorldtubeReplay(ParameterInput *pin, Mesh *pm, bool is_restart);
  ~EmriInnerWorldtubeReplay() = default;

  void InjectEField(Mesh *pm, Driver *pdrive, int stage);
  void ApplyPrimitiveBoundary(Mesh *pm, Driver *pdrive, int stage);
  void CompleteStep(Mesh *pm, Driver *pdrive);
  void CapTimestep(Mesh *pm);
  void SetADMVariables(MeshBlockPack *pack);
  bool FluidBoundaryEnabled() const { return fluid_boundary_enabled_; }
  bool ADMVolumeEnabled() const { return adm_volume_enabled_; }
  // CUDA extended lambdas require their enclosing member to be publicly accessible.
  void SetInitialNormalFlux(Mesh *pm);
  Real BoundaryFluxResidual(Mesh *pm, const DvceArray1D<Real> &expected) const;

 private:
  struct FaceRecord {
    int face;
    int m;
    int k;
    int j;
    int i;
    int v;
    int u;
  };

  struct EdgeRecord {
    int input;
    int m;
    int k;
    int j;
    int i;
    int component;
    int sign;
  };

  struct VolumeFaceRecord {
    int input;
    int m;
    int k;
    int j;
    int i;
    int component;
  };

  void ReadHeaderAndTimes(ParameterInput *pin, Mesh *pm);
  void ReadADMVolume(ParameterInput *pin, Mesh *pm);
  void BuildBoundaryTopology(Mesh *pm);
  void LoadInterval(int interval);
  void LoadADMInterval(int interval);
  Mesh *pmesh_ = nullptr;
  bool is_restart_ = false;
  bool exhausted_ = false;
  bool fluid_boundary_enabled_ = false;
  bool characteristic_sr_boundary_ = false;
  bool characteristic_gr_boundary_ = false;
  bool adm_volume_enabled_ = false;
  int cells_per_edge_ = 0;
  int nvar_ = 0;
  bool has_cell_centered_magnetic_state_ = false;
  bool has_initial_volume_flux_ = false;
  int nt_ = 0;
  int interval_ = 0;
  int last_stage_ = 0;
  int adm_nx_ = 0;
  int adm_ny_ = 0;
  int adm_nz_ = 0;
  int adm_nvar_ = 0;
  int adm_interval_ = 0;
  Real data_center_[3] = {0.0, 0.0, 0.0};
  Real half_width_ = 0.0;
  Real dx_ = 0.0;
  Real time_offset_ = 0.0;
  Real flux_tolerance_ = 1.0e-10;
  Real characteristic_speed_tolerance_ = 1.0e-10;
  std::string path_;
  std::string adm_path_;
  std::uint64_t adm_data_offset_ = 0;
  std::uint64_t adm_slab_values_ = 0;
  Real adm_lower_[3] = {0.0, 0.0, 0.0};
  Real adm_spacing_[3] = {0.0, 0.0, 0.0};
  std::vector<Real> times_;
  std::array<std::uint64_t, 6> state_offsets_{};
  std::array<std::uint64_t, 6> flux_offsets_{};
  std::array<std::uint64_t, 6> emf_u_offsets_{};
  std::array<std::uint64_t, 6> emf_v_offsets_{};
  std::array<std::uint64_t, 3> initial_volume_flux_offsets_{};
  std::vector<FaceRecord> host_faces_;
  std::vector<EdgeRecord> host_edges_;
  std::vector<VolumeFaceRecord> host_volume_faces_;
  DvceArray2D<int> face_records_;
  DvceArray2D<int> edge_records_;
  DvceArray2D<int> volume_face_records_;
  DvceArray1D<Real> state_left_;
  DvceArray1D<Real> state_right_;
  DvceArray1D<Real> flux_left_;
  DvceArray1D<Real> flux_right_;
  DvceArray1D<Real> interval_emf_;
  DvceArray1D<Real> initial_volume_flux_;
  DvceArray1D<Real> adm_left_;
  DvceArray1D<Real> adm_right_;
  DvceArray1D<int> characteristic_diagnostics_;
};

void EmriInnerWorldtubeInjectEField(Mesh *pm, Driver *pdrive, int stage);
void EmriInnerWorldtubeApplyPrimitiveBoundary(Mesh *pm, Driver *pdrive, int stage);
void EmriInnerWorldtubeCompleteStep(Mesh *pm, Driver *pdrive);
void EmriInnerWorldtubeCapTimestep(Mesh *pm);
void EmriInnerWorldtubeSetADMVariables(MeshBlockPack *pack);

#endif  // PGEN_EMRI_INNER_WORLDTUBE_HPP_
