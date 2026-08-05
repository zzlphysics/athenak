# Effective binary-black-hole GRMHD

Build AthenaK with `-DPROBLEM=dynbbh` and run
`inputs/dyngr/effective_bbh.athinput`.  The problem generator supplies the boosted,
superposed Kerr–Schild metric of [Combi & Ressler](https://arxiv.org/abs/2403.13308)
through AthenaK's prescribed dynamic-ADM interface.  The default initializes a uniform
gas and an optional uniform magnetic field.  The optional `reference_fm_torus` initial
data maps a Fishbone--Moncrief reference torus into the instantaneous binary ADM metric;
it is useful circumbinary initial data, but is not an exact equilibrium of the
time-dependent binary spacetime.

The metric kernel is a clean implementation of the equations in the paper, not a copy of
the generated C in the public archive.  The reference trajectory workflow is the authors'
[AnalyticalBBH repository](https://gitlab.com/combi.luciano/analyticalbbh), inspected at
commit [`2ce68e3d`](https://gitlab.com/combi.luciano/analyticalbbh/-/commit/2ce68e3d49e8758b32efc8841d239354d8d619d6).
The older [Zenodo v1 archive](https://zenodo.org/records/10841021) should not be used as a
golden arbitrary-spin implementation; its generated source predates corrections present
in the current repository.

## Trajectories

`trajectory_mode=circular` is an analytic Newtonian circular orbit intended for smoke
tests.  An inspiral, merger transition, remnant spin, and kick compatible with the
reference workflow require `trajectory_mode=table`.  Generate the HDF5 trajectory with
the authors' current AnalyticalBBH/CBwaves workflow, then convert it:

```bash
python3 scripts/bbh_trajectory_h5_to_ascii.py trajectory.h5 trajectory.dat
```

The converter normalizes all dimensional columns by the initial total mass unless
`--mass-unit` is supplied.  Legacy files with missing spin datasets are rejected; use
`--assume-missing-spin-zero` only after confirming that the omitted components really are
zero.  Existing outputs are not overwritten unless `--force` is passed.

The ASCII schema is one row per time and 21 columns:

```text
t x1 y1 z1 x2 y2 z2 vx1 vy1 vz1 vx2 vy2 vz2
a1x a1y a1z a2x a2y a2z m1_full m2_full
```

Here `aAi=SAi/MA` has dimensions of length.  Its components are expressed in the
instantaneous moving-hole spatial axes used by the paper; this implementation chooses
those axes parallel to the global axes before each pure boost.  The table must cover the
full simulation interval plus `metric_fd_step` on both ends.  AthenaK uses cubic Hermite
interpolation for each position/velocity pair and linear interpolation for masses and
spin components once per RK substage; it never extrapolates.
`trajectory_time_offset` maps simulation time zero into the table as
`t_table=t_simulation+trajectory_time_offset`.  Individual component masses may reach
zero after merger, but their sum must remain positive.  The converter checks dimensions,
finiteness, mass signs, subluminal velocities, and consistency between reported velocity
and finite-difference `dx/dt`, but cannot certify the physical consistency of a merger
fit.  Inspect the mass split, remnant velocity, and both spin vectors, especially for
unequal-mass runs.  On input, AthenaK additionally uses the quadratic Bézier control
vectors of every Hermite interval as a conservative, rigorous subluminal-speed
certificate; a badly undersampled position/velocity pair is rejected before evolution.

The direct two-term SKS superposition is not guaranteed to retain Lorentz signature when
two distinct terms overlap strongly.  AthenaK therefore checks representative points at
every table interval used by the run, and retains exact per-cell ADM checks.  A table must
start its smooth inspiral-to-remnant transition before this check fails; a pre-merger-only
table cannot simply be evolved through close approach.  The circular mode is deliberately
rejected at unsafe small separations.  This dense host-side scan is a one-time startup
cost; expect large trajectory tables with tens of thousands of rows to take tens of
seconds to validate on one CPU rank.

The generated reference C hides a negative lapse squared with an absolute value.  This
implementation does not do that, because it would no longer be the ADM decomposition of
the supplied four-metric.  Instead, each subextremal term is smoothly tapered only inside
its Kerr `r_+` surface and geometrically excised; the exact paper metric is unchanged in
the causal exterior.  `spin_buffer1/2` and `singularity_floor` provide an auxiliary ring
cutoff for superextremal half-mass terms during the canonical remnant interpolation.
`singularity_floor` is dimensionless and relative to each component mass (the legacy
`cutoff_floor` alias has the same revised meaning), and each
`spin_buffer+singularity_floor` sum must not exceed one.

## Required integration settings

- Use `<adm> dynamic=true` and an `<mhd>` block.  AthenaK does not support `<adm>+<hydro>`.
- Excision is required; use `<coord> excise=true` and `excision_scheme=lapse`.
- Only RK1, RK2, and RK3 are enabled for this time-dependent prescribed metric.
- With AMR, use `nghost=4` (or another supported even value) and
  `<mesh_refinement> prolong_primitives=false`.  The example keeps `nghost=3` because AMR
  is disabled by default.
- Radiation and Kerr-specific derived outputs such as `mhd_jcon` are not supported by
  this setup because those paths still assume a single stationary Kerr metric.

## Reference circumbinary torus initial data

Set `initial_data=reference_fm_torus` to construct a constant-angular-momentum
Fishbone--Moncrief torus using a stationary reference potential and then normalize its
four-velocity in the actual effective-BBH ADM fields at simulation time zero.  The
construction is therefore a controlled, reproducible initial condition, not a claim of
hydrostatic equilibrium in a nonstationary, nonaxisymmetric spacetime.  Its initial
transient and dependence on the reference potential must be included in a scientific
error budget.

The compact reference setup uses:

```text
initial_data           = reference_fm_torus
torus_reference_mass   = 1
torus_reference_center1 = 0
torus_reference_center2 = 0
torus_reference_center3 = 0
torus_reference_velocity1 = 0
torus_reference_velocity2 = 0
torus_reference_velocity3 = 0
torus_r_edge           = 18
torus_r_peak           = 29
torus_rho_max          = 1
torus_rho_min          = 1e-5
torus_rho_pow          = -1.5
torus_pgas_min         = 3.333333333333e-8
torus_pgas_pow         = -2.5
torus_rho_atm          = 1e-10
torus_temp_atm         = 3.333333333333e-13
torus_pert_amp         = 0.02
torus_seed             = 1
torus_magnetic_field   = single_loop
torus_potential_cutoff = 0.2
torus_mag_norm         = density_weighted_beta
torus_mag_target       = 100
torus_min_grid_peak_fraction = 0.9
torus_min_magnetic_cells = 64
torus_require_full_domain = true
```

The binary, `r_edge=18M`, `r_peak=29M`, single-loop field, and density-weighted
`beta=100` follow the compact disk experiment in
[Combi & Ressler (2024), section V](https://arxiv.org/pdf/2403.13308#page=11).
That paper does not specify the atmosphere or perturbation for this experiment.  The
power-law background, vector-potential cutoff, and two-percent pressure perturbation used
here instead follow the independently documented
[AthenaK magnetized-disk test](https://arxiv.org/html/2409.10384#S4.SS6); this mixed
provenance must be stated rather than described as an exact reproduction.  For the
Valencia solver, set both `torus_temp_atm` and
`mhd/tfloor=3.333333333333e-13`; these are temperature floors, not the HARM-like pressure
floor.  The atmosphere pressure is
`max(torus_pgas_min*(r/M_ref)^torus_pgas_pow, rho_background*torus_temp_atm)`.

`torus_rho_min*(r/M_ref)^torus_rho_pow` and
`torus_pgas_min*(r/M_ref)^torus_pgas_pow` define the radial atmosphere profile.
`torus_rho_atm` is the absolute density lower bound and `torus_temp_atm` is the absolute
temperature lower bound.  The
EOS floors must not exceed them.  Only `mhd/eos=ideal` and `mhd/dyn_eos=ideal` are
currently accepted, so the Gamma-law thermodynamics and internal-energy normalization
cannot be silently applied to another EOS.  `torus_seed` makes the pressure perturbation
reproducible across restarts and MPI layouts.

`torus_reference_mass`, `torus_reference_center{1,2,3}`, and
`torus_reference_velocity{1,2,3}` define the stationary Schwarzschild reference potential
and its initial coordinate frame.  If omitted, they are inferred from the mass-scaled,
mass-weighted trajectory state at the run's initial time.  The translational velocity is
added before normalizing the four-velocity in the actual BBH ADM metric; this is a useful
general initial frame, not an exact Lorentz-boosted equilibrium.  `torus_r_edge` and
`torus_r_peak` are absolute code-coordinate lengths.  If omitted they default to
`18*M_ref` and `29*M_ref`; explicitly supplied values must be scaled by the user when
changing the reference mass.  Keep the approximate reference translation small; input is
rejected unless its Euclidean coordinate speed is below one, and every torus cell must
still pass the actual-metric timelike check.

`torus_magnetic_field=none` gives a hydrodynamic reference case.  `single_loop` constructs
a cell-edge vector potential from density above `torus_potential_cutoff`, then takes its
discrete constrained-transport curl.  This is the supported path for preserving
divergence at static- and adaptive-refinement interfaces; do not initialize cell-centered
magnetic components independently.

The magnetic normalization choices are:

- `density_weighted_beta`: ratio
  `sum(rho*P*sqrt(det(gamma_ij))*dV)/sum(rho*b^2/2*sqrt(det(gamma_ij))*dV)` over
  the analytic torus;
- `peak_beta`: ratio based on the separately measured pressure and magnetic-pressure
  extrema;
- `integrated_pressure`: volume-integrated gas pressure divided by volume-integrated
  magnetic pressure;
- `integrated_internal_energy`: volume-integrated internal energy divided by
  volume-integrated magnetic pressure.

The three sum-based choices form compensated double-precision sums within each
MeshBlock, gather the block contributions in global-GID order, and compute and broadcast
one field scale from rank zero.  This keeps fresh initial data invariant when the same
mesh and homogeneous platform are repartitioned across MPI ranks.  It is not a promise of
bitwise identity across different CPU/GPU architectures or compiler math modes.

The paper says density-averaged beta but does not publish its exact averaging formula.
`density_weighted_beta` above is therefore AthenaK's explicit reproducible convention,
not a claim that the private paper implementation is algebraically identical.  The
meaning of `torus_mag_target` follows the selected normalization.  Values from different
choices are not interchangeable, so record the choice, unnormalized ratio, target, and
applied field scale printed at initialization in run metadata.

Initialization fails if the resolved density peak is below
`torus_min_grid_peak_fraction*torus_rho_max`, or if fewer than
`torus_min_magnetic_cells` active cells support the magnetic loop.  With
`torus_require_full_domain=true`, the analytic outer equatorial radius must fit inside
every domain direction; setting it false is an explicit opt-out for deliberately
truncated software tests, not a production recommendation.

Two inputs deliberately separate software coverage from production planning:

- `effective_bbh_torus_smoke.athinput` uses a `160M`-wide root domain with one static
  refinement box.  Its `dx=1.25M` inner spacing only makes the `r_edge=18M`, `r_peak=29M`
  torus visible; it is not adequate for MRI, horizon, flux, turbulence, or accretion
  measurements.
- `effective_bbh_torus_production_template.athinput` records the equal-mass, zero-spin,
  separation-`20M`, Gamma-`4/3`, RK2/CFL-`0.3`, HLLE/WENO-Z reference setup on a
  `[-1000M,1000M]^3` domain.  It intentionally contains a missing trajectory, a `32^3`
  unrefined resource guard, and a two-cycle limit, so it is not directly runnable as a
  production calculation.

The production resolution goal is `dx <= M/64` near both horizons and throughout every
region used for science diagnostics.  On the exact `[-1000M,1000M]` domain, a planned
`128^3` root plus ten dyadic refinement levels gives
`dx=2000/(128*2^10) M=0.0152587890625M`, slightly finer than `M/64`.  The template does
not enable these boxes: their extents must follow the verified trajectory while retaining
the cavity and disk science region, and the resulting MeshBlock count, GPU memory,
checkpoint size, output cadence, and NAS transfer rate must be measured first.  Merely
resolving the horizons does not resolve the MRI or justify a magnetic-structure claim.

Before a paper run, replace the placeholder with a provenance-recorded PN-to-remnant
trajectory, test `metric_fd_step`, add validated moving-hole accretion and magnetic-flux
diagnostics, address horizon-following fluid/magnetic forcing, and demonstrate spatial,
temporal, floor/excision, trajectory, and initial-transient convergence.

### Strict 2:1 level subcycling

Set `<time> subcycling=level` to advance every physical refinement level with exactly
twice as many steps as its parent.  The recursive driver time-interpolates coarse boundary
data, accumulates conserved-variable flux registers, and applies a constrained-transport
EMF reflux-curl whenever a child catches its parent.  Regridding, diagnostics, output, and
restart are permitted only at synchronized root-level endpoints.

The validated path supports multilevel static SMR and synchronized moving AMR, including
weighted load balancing over multiple MPI ranks.  Shared restart files may be resumed with
a different rank count.  Per-rank restart files record their writer partition, AMR
derefinement counters, and a 128-bit checkpoint identity, and therefore require the same
rank layout.  `inputs/dyngr/effective_bbh_smr3_smoke.athinput` exercises three physical
levels; `inputs/dyngr/effective_bbh_amr_subcycle_smoke.athinput` creates, moves, and removes
three-level refinement regions around the holes.

The regression suite also moves and fully de-refines a 3D magnetic octet on cells with a
`1:2:2` aspect ratio.  It checks divergence- and curl-preserving Toth-Roe face-flux
prolongation, same-rank and MPI de-refinement, and a three-to-two-rank restart at every
synchronized root cycle.

The implementation remains deliberately fail-fast outside its validated physics envelope:
one MeshBlockPack per MPI rank, RK2, single-fluid MHD+CT with a prescribed dynamic ADM
metric, conserved-variable prolongation, and no radiation, particles, Z4c, diffusion,
built-in source terms, turbulence, shearing-box, orbital-advection, or full-mesh user
boundary/source modules.  `mhd_jcon` output is also unavailable because it is not
time-consistent on intermediate levels and assumes a single Kerr metric.  Use
`subcycling=none` for the historical global-timestep path.

Reducing the number of MeshBlock updates does not guarantee a GPU wall-clock speedup on a
small mesh: level-local kernel launches and synchronization can dominate when individual
levels underfill the device.  Benchmark with the intended block count and MPI layout.  The
smoke inputs are far too coarse for scientific horizon or magnetic-structure measurements.

The implementation uses second-order centered spacetime derivatives to construct
`K_ij`.  Test `metric_fd_step` at `h/2`, `h`, and `2h`; making it arbitrarily small can
increase floating-point and table-interpolation error.  The metric is an approximation
suitable for low-self-gravity gas dynamics, not a replacement for numerical relativity or
an accurate gravitational-wave model.

The example's `48^3` root grid has `dx=1` and does not resolve either horizon for
scientific work.  It is a CPU integration test.  Enabling the included tracking criterion
requires changing `refinement=adaptive`, switching to an even ghost-zone count, and then
performing resolution convergence.

The excision path combines AthenaK's lapse mask with a conservative geometric mask around
the compact metric-regularization regions and uses atmosphere recovery.  It does not yet
implement the horizon-following fluid forcing used to suppress magnetic leakage in long,
strongly magnetized production runs; see the time-dependent GRMHD treatment of
[Ressler et al.](https://arxiv.org/abs/2404.02193).  Such runs need resolution and
excision-convergence tests before scientific use.

## Related applications

The published applications motivating this setup are
[Ressler et al. 2025](https://arxiv.org/abs/2410.10944) (inspiral-only Athena++),
[Fedrigo et al. 2026](https://arxiv.org/abs/2603.00224) (GIZMO GRHD, without magnetic
fields), and [Jiang et al. 2026](https://arxiv.org/abs/2603.10792) (AthenaK with a separate
PN setup).  They establish useful physical and numerical targets, but do not provide a
drop-in AthenaK production input and implementation for this repository.  Reproducing
their scientific runs therefore also requires their initial conditions, refinement
strategy, diagnostics, and a large-resource convergence campaign.
