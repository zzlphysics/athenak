#ifndef PGEN_EMRI_OUTER_WORLDTUBE_HPP_
#define PGEN_EMRI_OUTER_WORLDTUBE_HPP_
//========================================================================================
// AthenaK astrophysical plasma code
// Copyright(C) 2026 AthenaK collaboration
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file emri_outer_worldtube.hpp
//! \brief Fixed, grid-aligned outer-worldtube recorder for EMRI multiscale matching.

#include <array>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

#include "athena.hpp"

class Driver;
class Mesh;
class ParameterInput;

class EmriOuterWorldtubeWriter {
 public:
  EmriOuterWorldtubeWriter(ParameterInput *pin, Mesh *pm, bool is_restart);
  ~EmriOuterWorldtubeWriter();

  void ObserveEField(Mesh *pm, Driver *pdrive, int stage);
  void CompleteStep(Mesh *pm, Driver *pdrive);
  void Finalize(Mesh *pm, Driver *pdrive);

 private:
  struct CellRecord {
    int face;
    int m;
    int k;
    int j;
    int i;
    int v;
    int u;
  };

  struct EdgeRecord {
    int output;
    int m;
    int k;
    int j;
    int i;
    int component;
    int sign;
  };

  struct FaceStreams {
    std::ofstream cell_state;
    std::ofstream normal_flux;
    std::ofstream emf_u;
    std::ofstream emf_v;
  };

  void BuildTopology(Mesh *pm);
  void OpenFiles(Mesh *pm, Driver *pdrive);
  void CaptureAndWriteEndpoint(Mesh *pm, Real time);
  void WriteInterval(Mesh *pm, Real interval_dt);
  void WriteManifest(bool complete) const;
  void CloseFiles();

  Mesh *pmesh_ = nullptr;
  bool is_restart_ = false;
  bool initialized_ = false;
  bool finalized_ = false;
  bool overwrite_ = false;
  int cells_per_edge_ = 0;
  int nvar_ = 0;
  int dcycle_ = 1;
  int last_stage_ = 0;
  int interval_steps_ = 0;
  std::int64_t endpoints_written_ = 0;
  std::int64_t intervals_written_ = 0;
  Real center_[3] = {0.0, 0.0, 0.0};
  Real half_width_ = 0.0;
  Real dx_ = 0.0;
  Real interval_start_time_ = 0.0;
  std::string basename_;
  std::string stem_;
  std::string manifest_path_;
  std::string times_path_;
  std::vector<std::string> state_names_;
  std::vector<CellRecord> host_cells_;
  std::vector<EdgeRecord> host_edges_;
  DvceArray2D<int> cell_records_;
  DvceArray2D<int> edge_records_;
  DvceArray1D<Real> step_emf_integral_;
  DvceArray1D<Real> interval_emf_integral_;
  DvceArray1D<Real> endpoint_state_;
  DvceArray1D<Real> endpoint_flux_;
  std::ofstream times_stream_;
  std::array<FaceStreams, 6> face_streams_;
};

void EmriOuterWorldtubeObserveEField(Mesh *pm, Driver *pdrive, int stage);
void EmriOuterWorldtubeCompleteStep(Mesh *pm, Driver *pdrive);
void EmriOuterWorldtubeFinalize(Mesh *pm, Driver *pdrive);

#endif  // PGEN_EMRI_OUTER_WORLDTUBE_HPP_
