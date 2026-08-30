#ifndef PGEN_PGEN_HPP_
#define PGEN_PGEN_HPP_
//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file pgen.hpp
//  \brief definitions for ProblemGenerator class

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "geodesic-grid/spherical_grid.hpp"
#include "parameter_input.hpp"

class Driver;
class EmriInnerWorldtubeReplay;
struct EmriADMReplayGeometry;
class EmriOuterWorldtubeWriter;

void EmriOuterWorldtubeObserveEField(Mesh *pm, Driver *pdrive, int stage);
void EmriOuterWorldtubeCompleteStep(Mesh *pm, Driver *pdrive);
void EmriOuterWorldtubeFinalize(Mesh *pm, Driver *pdrive);
void EmriInnerWorldtubeInjectEField(Mesh *pm, Driver *pdrive, int stage);
void EmriInnerWorldtubeApplyPrimitiveBoundary(Mesh *pm, Driver *pdrive, int stage);
void EmriInnerWorldtubeCompleteStep(Mesh *pm, Driver *pdrive);
void EmriInnerWorldtubeCapTimestep(Mesh *pm);
void EmriInnerWorldtubeSetADMVariables(MeshBlockPack *pmbp);
bool EmriInnerWorldtubeADMReplayGeometry(Mesh *pm, EmriADMReplayGeometry &geometry);

using ProblemFinalizeFnPtr = void (*)(ParameterInput *pin, Mesh *pm);
using UserBoundaryFnPtr = void (*)(Mesh* pm);
using UserPrimitiveBoundaryFnPtr = void (*)(Mesh* pm, Driver* pdrive, int stage);
using UserEFieldFnPtr = void (*)(Mesh* pm, Driver* pdrive, int stage);
using UserStepFnPtr = void (*)(Mesh* pm, Driver* pdrive);
using UserTimestepFnPtr = void (*)(Mesh* pm);
using UserSrctermFnPtr = void (*)(Mesh* pm, const Real bdt);
using UserRefinementFnPtr = void (*)(MeshBlockPack* pmbp);
using UserHistoryFnPtr = void (*)(HistoryData *pdata, Mesh *pm);
using UserOutputRegionFnPtr = bool (*)(const std::string &name, Real time,
                                       Real center[3]);

//----------------------------------------------------------------------------------------
//! \class ProblemGenerator

class ProblemGenerator {
 public:
  // constructor for new problems
  ProblemGenerator(ParameterInput *pin, Mesh *pmesh);
  // constructor for restarts
  ProblemGenerator(ParameterInput *pin, Mesh *pmesh, IOWrapper resfile,
                   bool single_file_per_rank=false);
  ~ProblemGenerator();

  // true if user BCs are specified on any face
  bool user_bcs;

  // true if user srcterms are specified
  bool user_srcs;

  // true if user history outputs are specified
  bool user_hist;

  // vector of SphericalGrid objects for analysis
  std::vector<std::unique_ptr<SphericalGrid>> spherical_grids;

  // function pointer for final work after main loop (e.g. compute errors).  Called by
  // Driver::Finalize()
  ProblemFinalizeFnPtr pgen_final_func=nullptr;
  // function pointer for user-enrolled BCs.  Called in ApplyPhysicalBCs in task list
  UserBoundaryFnPtr user_bcs_func=nullptr;
  // Optional primitive-cache boundary hook.  It runs after C2P, so a boundary Riemann
  // problem can supply an exterior primitive state without rewriting outgoing conserved
  // modes before the ordinary solver sees the interface.
  UserPrimitiveBoundaryFnPtr user_primitive_bcs_func=nullptr;
  // Optional edge-electric-field hook, called after built-in E-field source terms and
  // before EMF communication/CT at every explicit RK stage.
  UserEFieldFnPtr user_efld_func=nullptr;
  // Read the synchronized edge field after EMF communication and immediately before CT.
  UserEFieldFnPtr user_efld_observer_func=nullptr;
  // Observe a fully completed root step after Mesh::time and Mesh::ncycle are advanced.
  UserStepFnPtr user_step_observer_func=nullptr;
  // Flush problem-owned streaming diagnostics before final problem analysis runs.
  UserStepFnPtr user_finalize_observer_func=nullptr;
  // Apply a problem-specific upper bound after the ordinary CFL timestep is selected.
  UserTimestepFnPtr user_timestep_func=nullptr;
  UserSrctermFnPtr user_srcs_func=nullptr;
  UserRefinementFnPtr user_ref_func=nullptr;
  UserHistoryFnPtr user_hist_func=nullptr;
  // Resolve a problem-defined, time-dependent center used by filtered mesh output.
  // Returning false rejects an unknown region name before any partial file is written.
  UserOutputRegionFnPtr user_output_region_func=nullptr;

  // predefined problem generator functions (default test suite)
  void CallProblemGenerator(ParameterInput *pin, bool is_restart);
  void Advection(ParameterInput *pin, const bool restart);
  void AlfvenWave(ParameterInput *pin, const bool restart);
  void BondiAccretion(ParameterInput *pin, const bool restart);
  void CShock(ParameterInput *pin, const bool restart);
  void Diffusion(ParameterInput *pin, const bool restart);
  void LinearWave(ParameterInput *pin, const bool restart);
  void LWImplode(ParameterInput *pin, const bool restart);
  void Monopole(ParameterInput *pin, const bool restart);
  void MRI3d(ParameterInput *pin, const bool restart);
  void OrszagTang(ParameterInput *pin, const bool restart);
  void ShockTube(ParameterInput *pin, const bool restart);
  void Shwave(ParameterInput *pin, const bool restart);
  void SphericalCollapse(ParameterInput *pin, const bool restart);
  void RadiationLinearWave(ParameterInput *pin, const bool restart);
  void RadiationBeam(ParameterInput *pin, const bool restart);
  void Z4cBoostedPuncture(ParameterInput *pin, const bool restart);
  void Z4cLinearWave(ParameterInput *pin, const bool restart);

  // Generic error output function (using difference u0-u1)
  void OutputErrors(ParameterInput *pin, Mesh *pm);

  // template for user-specified problem generator
  void UserProblem(ParameterInput *pin, const bool restart);

 private:
  friend void EmriOuterWorldtubeObserveEField(Mesh*, Driver*, int);
  friend void EmriOuterWorldtubeCompleteStep(Mesh*, Driver*);
  friend void EmriOuterWorldtubeFinalize(Mesh*, Driver*);
  friend void EmriInnerWorldtubeInjectEField(Mesh*, Driver*, int);
  friend void EmriInnerWorldtubeApplyPrimitiveBoundary(Mesh*, Driver*, int);
  friend void EmriInnerWorldtubeCompleteStep(Mesh*, Driver*);
  friend void EmriInnerWorldtubeCapTimestep(Mesh*);
  friend void EmriInnerWorldtubeSetADMVariables(MeshBlockPack*);
  friend bool EmriInnerWorldtubeADMReplayGeometry(Mesh*, EmriADMReplayGeometry&);

  void ConfigureEmriOuterWorldtube(ParameterInput *pin, bool is_restart);

  std::unique_ptr<EmriInnerWorldtubeReplay> emri_inner_worldtube_;
  std::unique_ptr<EmriOuterWorldtubeWriter> emri_outer_worldtube_;
  bool single_file_per_rank; // for restart file naming
  Mesh* pmy_mesh_;
};

#endif // PGEN_PGEN_HPP_
