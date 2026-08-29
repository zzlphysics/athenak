# Local EMRI GRMHD wind tunnel

Build with

```sh
cmake -S . -B build_emri -DPROBLEM=emri_windtunnel
cmake --build build_emri -j
```

and run the small CPU smoke problem with

```sh
./build_emri/src/athena -i inputs/emri/emri_windtunnel_smoke.athinput
```

## Spacetime model

The numerical domain is centered on the secondary black hole and does not contain the
primary.  The prescribed external metric nevertheless retains the primary Kerr field.
At `t=0` the primary is at global Cartesian Kerr-Schild position `(0,0,0)` and the
secondary orbit is anchored at `(R0,0,0)`, where

```text
R0 = sqrt(orbital_radius^2 + (primary_chi*primary_mass)^2).
```

Local axes rotate at `orbital_omega` and translate with the secondary.  The exact metric
pullback supplies the rotating-frame shift.  Circular equatorial motion and spins aligned
with the orbital z axis make the local metric stationary.  `omega_mode=kerr_geodesic`
uses the test-particle equatorial Kerr frequency; `omega_mode=custom` reads
`orbital_omega` directly.

The default `secondary_embedding=tangent_tetrad` constructs the small-hole Kerr-Schild
perturbation in an orthonormal tetrad of the external metric at the orbital anchor, then
maps that covariant perturbation into the numerical chart.  This is the appropriate
local construction: the constant part of the primary field changes the tetrad, while
only its gradients and curvature across the box can produce locally measurable
corrections.  `secondary_embedding=global_boost` retains the earlier global
Minkowski-boost prescription as a legacy control; it should not be used for EMRI science
because it can turn an order-`M/r` coordinate potential into a spurious order-unity
change of the secondary near field.

Use `<adm> dynamic=false` for ordinary global-timestep evolution: the stationary metric
is cached and is regenerated automatically after AMR topology changes.  The repository's
strict `time/subcycling=level` scheduler currently requires `dynamic=true`; that mode
re-evaluates the same stationary metric at RK stages and is therefore more expensive.

This is an effective prescribed spacetime, not a matched-asymptotic or constraint-solved
EMRI metric.  It keeps the full primary Kerr field in the external sector and an exact
small-hole Kerr perturbation in the anchor tetrad, but neglects nonlinear cross terms,
inspiral, radiation reaction, primary recoil, and disk self-gravity.  Those assumptions
are intentional for a one-way local wind tunnel.

## Controlled background modes

`background_mode` defines three runs that keep the secondary, wind, magnetic field,
domain, and diagnostic radii fixed:

| mode | primary Kerr term | orbital pullback | role |
|---|---:|---:|---|
| `full` | yes | yes | effective EMRI wind tunnel |
| `frame_only` | no | yes | fixed-chart counterfactual for adding the primary metric |
| `isolated` | no | no | inertial isolated-Kerr BHL control |

The useful contrasts are `full-frame_only` and `frame_only-isolated`.  The first changes
only the primary Kerr contribution while holding the raw orbital chart fixed.  The second
shows the effect of the orbital non-inertial chart.  With the default source-tetrad wind,
force components, and physical extraction spheres, the three runs specify the same
locally measured upstream state and diagnostics.  In `frame_only`, centrifugal
acceleration is intentionally not balanced by primary gravity, so its nonzero
`geo_resid` is expected.  The contrasts remain attribution experiments between
prescribed effective spacetimes, not globally gauge-invariant binary observables.

Run the same input in three separate directories, changing only the mode and basename:

```sh
mkdir -p runs/emri/{full,frame_only,isolated}
./build_emri/src/athena -d runs/emri/full \
  -i inputs/emri/emri_windtunnel_smoke.athinput \
  problem/background_mode=full job/basename=full
./build_emri/src/athena -d runs/emri/frame_only \
  -i inputs/emri/emri_windtunnel_smoke.athinput \
  problem/background_mode=frame_only job/basename=frame_only
./build_emri/src/athena -d runs/emri/isolated \
  -i inputs/emri/emri_windtunnel_smoke.athinput \
  problem/background_mode=isolated job/basename=isolated
```

The selected mode, secondary embedding, and wind/force frame choices are part of the
restart contract, so a checkpoint cannot silently be continued with another model.
Prototype checkpoints with older contracts are deliberately rejected.

## Units and scale checks

All masses and lengths use `G=c=1`.  The natural local choice is
`secondary_mass=1`, so `primary_mass=1/q`.  The startup diagnostic reports

```text
q, orbital_radius/primary_mass, orbital_omega*primary_mass,
r_H/secondary_mass, box_size/orbital_radius.
```

Interpretation as a local model requires the box to be much smaller than the orbital
radius and not to intersect the primary horizon.  Extreme mass ratios require a
double-precision build.  Metric derivatives use separate finite-difference scales for
the secondary near field and the slowly varying primary/rotating-frame contribution;
these can be overridden with `metric_fd_step` and `external_metric_fd_step`.  Before
production, compare the accretion scale

```text
r_a ~ 2 m / (v_rel^2 + c_s^2 + v_A^2)
```

with the Hill radius

```text
r_H ~ orbital_radius * (q/3)^(1/3).
```

The present uniform-wind stage is most informative when runs are organized by `r_a/r_H`,
not by a brute-force scan in `q`.

## Capture-scale and cost planner

`plan_bhl_hierarchy.py` turns an upstream source-tetrad state into a reproducible decision
between a direct relativistic calculation, a scale-separated outer--inner calculation,
and a disk-scale outer model.  Run it before choosing a box:

```bash
python3 inputs/emri/plan_bhl_hierarchy.py \
  --secondary-mass 1 \
  --primary-mass 100000 \
  --orbital-radius 1000000 \
  --rho 1 --pgas 1e-4 --adiabatic-index 1.3333333333333333 \
  --four-velocity 0.05 0 0 \
  --magnetic-field 0.01 0 0 \
  --coherence-scale disk_H=100000 \
  --coherence-scale L_rho=200000 \
  --output-prefix plans/slow-wind
```

All inputs must already use one common `G=c=1` unit system.  With
`wind_frame=source_tetrad`, the velocity and magnetic arguments have exactly the same
meaning as `u1..u3` and `b1..b3` in the local problem.  For an ideal relativistic gas the
planner computes

```text
h = 1 + gamma p/[(gamma-1) rho],       cs^2 = gamma p/(rho h),
b^2 = [B^2 + (B dot u)^2]/W^2,        vA^2 = b^2/(rho h+b^2),
cf,proxy^2 = cs^2 + vA^2 - cs^2 vA^2.
```

It reports both `m/(v^2+cf,proxy^2)` and the conservative factor-two cost radius
`2m/(v^2+cf,proxy^2)`.  The latter is used to size the default domain.  This is a cost
proxy, not the exact directional GRMHD fast characteristic, which requires the full
dispersion relation.  The default upstream extent of four capture radii follows the
minimum scale commonly used in relativistic BHL calculations; the downstream side is
longer to retain the wake.  See the
[relativistic BHL domain discussion](https://academic.oup.com/mnras/article/471/3/3127/3979475)
and the magnetized calculations in [arXiv:2201.11753](https://arxiv.org/abs/2201.11753)
and [arXiv:2409.12359](https://arxiv.org/abs/2409.12359).

The output consists of strict JSON plus a readable Markdown report.  It gives a uniform
horizon-resolving cell count, an optimistic nested-box AMR count, the number of finest
steps, and global-timestep zone updates.  The last quantity prevents a misleading result
in which a modest resident mesh is evolved for millions of secondary-scale timesteps.
All budgets are command-line controls and should be replaced by measured AthenaK
throughput on the intended machine.

The decision tree is:

1. If the cost radius is smaller than the supplied disk/gradient scales and Hill radius,
   and the direct cell, level, step, and zone-update budgets pass, use direct GRMHD.
2. If direct evolution fails but
   `r_inner << r_match=sqrt(r_inner*r_a) << r_a`, evolve a nonrelativistic or SRMHD outer
   BHL problem with a finite sink, and replay a worldtube into the relativistic inner
   problem.  Vary `r_match`; the geometric mean is only the first trial.
3. If `r_a` approaches `H`, a fitted gradient scale, or `r_H`, the uniform BHL premise has
   already failed.  Use a global-disk or shearing-patch outer calculation and retain the
   relativistic inner calculation only if a smaller overlap still exists.  Enlarging a
   uniform box cannot repair the missing disk geometry or tidal shear.
4. For a weakly magnetized, sub-fast, nearly spherical state, a Bondi/Michel outer
   solution is a useful reduced control.  It is not a generic replacement for a
   magnetized disk wind.

For analytic or extracted gradients, provide the smallest physically relevant coherence
scale.  Useful estimates are `L_rho=1/|grad ln rho|`,
`L_p=1/|grad ln p|`, `L_u=|u|/||grad u||`, and
`L_B=|B|/||grad B||`, together with the disk height.  Run the auditor on every state in a
time-dependent table and design for the largest capture radius and smallest coherence
scale, not for the temporal mean.

### Outer--inner matching contract

The planner writes an explicit spatial-worldtube contract.  Its host-side flux/EMF
format, validator, conservative regridder, and first fixed-grid AthenaK outer writer are
implemented below.  The inner CT magnetic driver is also implemented; source-tetrad
transformation and incoming-characteristic fluid replay remain.  The outer calculation has a
secondary-centered cubical extraction worldtube and a sink strictly inside it.  The
inner calculation replaces the sink region.  A trustworthy implementation must transfer:

- incoming fluid characteristic amplitudes while leaving outgoing modes unconstrained;
- face-integrated normal magnetic flux, conservatively regridded so fine-face sums equal
  their parent-face flux;
- one edge-integrated, time-centered EMF shared by every incident CT face;
- the source tetrad, coordinate maps, and time-centering metadata.

Pointwise interpolation of cell-centered `B` is insufficient: it neither preserves the
discrete divergence constraint across different grids nor supplies Faraday-consistent
time evolution.  For the initial inner volume, reconstruct a field from transferred face
fluxes or a compatible vector potential; then use the imported edge EMFs at the six
worldtube faces.  A closed six-face worldtube also removes the current assumption that a
fixed Cartesian face remains upstream when the wind direction changes.

The outer sink must be converged by reducing its radius or shown to lie downstream of a
fast-magnetosonic critical surface; otherwise it can communicate with and alter the
matching data.  This causal issue is demonstrated explicitly in magnetized BHL sink
studies such as [Lee et al. 2017](https://academic.oup.com/mnras/article/468/1/717/3043736).

For a mean force, first settle the cheap outer problem for several `r_a/v_eff`, then run
the inner problem on a time-averaged boundary or on several separated stationary
windows, each several `r_match/v_eff` long.  The second option measures nonlinear
snapshot-to-snapshot scatter at much lower cost than replaying the complete outer
settling interval.  It does not retain low-frequency outer--inner correlations; spectra,
coherent variability, and feedback require continuous replay or a response model.

Do not add the force on the outer sink to the inner accreted-momentum force.  A
partition-of-unity overlap gives the intended bookkeeping,

```text
F_total = -Fmom_inner
          + integral w_inner f_inner,rel dV
          + integral w_outer f_outer,grav dV,
w_inner + w_outer = 1.
```

Move the matching surface and overlap width to expose double counting or a gap.  If an
inner jet or magnetic eruption reaches the matching surface and changes the outer wake,
the one-way hierarchy is physically inconsistent; iterate outer and inner solutions or
perform a two-way/global calculation.  The strong magnetic feedback seen in
[arXiv:2201.11753](https://arxiv.org/abs/2201.11753) and
[arXiv:2409.12359](https://arxiv.org/abs/2409.12359) makes this a production gate rather
than a formal caveat.

## Wind and boundaries

`rho0`, `pgas0`, `u1..u3`, and `b1..b3` set the magnetized wind at the source
anchor.  With all `d*_dxh*` parameters zero this is the original uniform-wind
configuration.
With the default `wind_frame=source_tetrad`, `u1..u3` are the Eulerian orthonormal spatial
four-velocity components and `b1..b3` are magnetic-field components in the
secondary-comoving external-background tetrad at the anchor.  The code reconstructs the
same tangent-chart GRMHD state point by point and maps it into each numerical slicing;
the Lorentz factor therefore remains timelike near the small hole without confusing a
coordinate boost with a physical change of wind.  `wind_frame=normal_frame` is the
legacy coordinate-primitive interpretation.

Source-frame wind uses a `user` flag on every face that supplies upstream data, because
the coordinate primitive varies over the boundary even when the tangent-frame state is
uniform.  The sample has positive `u1` and therefore uses `ix1_bc=user`; downstream faces
remain `outflow`.  For oblique or sub-fast flows, more than one user face may be required.
The zero-gradient face field has constant densitized coordinate components.  A future
replay boundary can replace this state with tetrad-projected data from a global
single-SMBH GRMHD run.

### Controlled gradients and shear

The analytic profile uses source-tangent coordinates
`xh1,xh2,xh3 = radial, prograde tangential, vertical`.  Density and pressure are
positive exponential profiles,

```text
rho  = rho0   * exp(sum_i dlnrho_dxhi  * xhi)
pgas = pgas0 * exp(sum_i dlnpgas_dxhi * xhi).
```

`dui_dxhj` is the full gradient tensor of the source-frame Eulerian spatial
four-velocity.  For a slow Newtonian Keplerian control, the leading shear is
approximately `du2_dxh1 = -q_sh*Omega`, with `q_sh=3/2`; strong-field production
values should instead be obtained from the chosen relativistic disk model.  These
analytic profiles are not automatically in pressure/tidal equilibrium, so startup
transients are part of the controlled experiment unless a balanced profile is supplied.
`max_log_contrast` rejects an accidentally excessive density or pressure contrast at a
domain corner.

`dbi_dxhj` specifies a linear source-tangent magnetic-flux-density gradient.  It must be
trace-free.  If `db3_dxh3` is omitted, its default is
`-db1_dxh1-db2_dxh2`; if all nine entries are present, the code checks the trace rather
than silently projecting it.  The source gradient is transformed as a contravariant
vector density into the numerical chart.  The constant anchor field retains the exact
legacy tetrad-to-slicing transformation.

The coordinate densitized field is initialized from

```text
A = 0.5 * (B0 cross x) - (1/3) * x cross (G x),
Bdens = curl(A) = B0 + G x,
```

using edge-centered vector potential values and the discrete CT curl.  Coarse edges
touching finer neighbors use fine-edge averages, so shared coarse/fine face fluxes are
consistent.  User-boundary ghost faces sample the same trace-free analytic field; active
face fluxes remain CT evolved.  A spatially varying profile therefore requires `user`
rather than the built-in constant `inflow` boundary on every supplying face.

## Static global-GRMHD Taylor extraction

`extract_static_taylor_worldtube.py` provides the first one-way global-to-local path.  A
global single-SMBH run must write co-temporal three-dimensional `mhd_w_bcc` and `adm`
binary dumps.  For example,

```bash
python3 inputs/emri/extract_static_taylor_worldtube.py \
  --state global.mhd_w_bcc.00100.bin \
  --adm global.adm.00100.bin \
  --output-prefix profiles/orbit-r10-t1000 \
  --anchor 10 0 0 \
  --primary-center 0 0 0 \
  --disk-normal 0 0 1 \
  --orbital-omega 0.03162277660168379 \
  --fit-radius 0.5
```

The extractor defaults to identical global/local code units.  When the global disk uses
`M_primary=1` but the local tunnel uses `m_secondary=1`, pass
`--global-length-in-local-units M_primary/m_secondary` (normally `1/q`).  If `L` is this
ratio and `D` is `--density-renormalization`, conversion gives
`t_local=L t_global`, `(rho,p)_local=D (rho,p)_global/L^2`,
`B_local=sqrt(D) B_global/L`, scalar and velocity gradients divided by `L`, and magnetic
gradients multiplied by `sqrt(D)/L^2`.  Choosing `D=L^2` keeps the density normalization
near its global numerical value while preserving pressure-to-density ratio and
magnetization.  Both factors are recorded in the manifest; relying on an implicit unit
convention is not supported.

`--source-velocity vx vy vz` can replace `--orbital-omega`; both are global coordinate
velocities `dx^i/dt`, not physical three-velocities.  The extractor:

1. volume-weight fits the ADM metric at the worldline anchor;
2. constructs an orthonormal source tetrad with radial, prograde, and vertical axes;
3. matches state and ADM leaf MeshBlocks by logical location, so differing file order is
   harmless;
4. projects AthenaK's `W v^i` and densitized `B^i` through four-vectors into the source
   tetrad;
5. Gaussian- and cell-volume-weight fits the scalar and velocity Taylor coefficients;
6. jointly fits all magnetic components with `tr(dB_i/dxh_j)=0` imposed exactly.

The output prefix produces an `.athinput` fragment and a `.json` manifest containing
hashes, time/cycle, tetrad/coframe, fitted parameters, residuals, cell-volume range, and
ready-to-use command-line overrides.  Copy the fragment values into the local
`<problem>` block.  If the residuals are large or the fit spans a wide range of AMR cell
volumes, the first-order model is not an adequate representation; use a smaller fitting
radius or proceed to a spatial worldtube replay rather than treating interpolation as
additional physical resolution.

This extractor is intentionally static and one-way.  It holds the anchor tetrad fixed
across the fitting sphere and replaces cell-centered magnetic samples with their best
trace-free linear approximation.  The profile-series stage below removes the
time-independence; the one-way and first-order spatial limitations remain explicit in
every manifest.

## Time-dependent Taylor worldtube replay

`build_taylor_worldtube_series.py` repeats the constrained extraction along a supplied
circular-equatorial worldline.  Its input is a JSON manifest; dump paths may be absolute
or relative to the manifest:

```json
{
  "classification": "athenak-emri-worldline-v1",
  "samples": [
    {
      "state": "dumps/global.mhd_w_bcc.00100.bin",
      "adm": "dumps/global.adm.00100.bin",
      "anchor": [10.0, 0.0, 0.0],
      "source_velocity": [0.0, 0.0316, 0.0]
    },
    {
      "state": "dumps/global.mhd_w_bcc.00101.bin",
      "adm": "dumps/global.adm.00101.bin",
      "anchor": [9.999, 0.141, 0.0],
      "source_velocity": [-0.00045, 0.0316, 0.0]
    }
  ]
}
```

Build the replay table and its provenance manifest with

```bash
python3 inputs/emri/build_taylor_worldtube_series.py \
  --manifest profiles/worldline.json \
  --output-prefix profiles/orbit-r10 \
  --primary-center 0 0 0 \
  --disk-normal 0 0 1 \
  --fit-radius 0.5 \
  --global-length-in-local-units 100000 \
  --density-renormalization 1.0e10
```

The builder rejects non-increasing dump times and a worldline whose radius, height,
angular frequency, radial speed, or vertical speed violates `--orbit-tolerance`.  This is
not merely a convenience check: the current local metric is circular and equatorial, so
feeding it an eccentric or inclined trajectory would combine incompatible frames.  The
converted coordinate radius and angular frequency are embedded in the table; the C++
loader rejects a table that does not match the configured local orbit within the same
tolerance.

Enable online replay with

```text
<adm>
dynamic = true

<problem>
profile_file          = profiles/orbit-r10.dat
profile_time_offset   = 1000.0
profile_extrapolation = error
```

The table time is `simulation_time + profile_time_offset`; use `auto` (the template
default) to align local startup with the first table row.  The resolved numeric offset is
written into restarts.  `profile_extrapolation=error` is the production-safe default,
while `hold` clamps to the first or last state for controlled tests.  Dynamic ADM refresh
supplies the RK-stage time, so the replay is synchronized with RK1/2/3.
Level subcycling is rejected for now because a single process-global boundary profile
cannot yet represent two simultaneously cached level times.

The restart contract stores a compact versioned digest of all metric/profile controls
and an FNV-1a digest of the table bytes.  A changed table is therefore rejected before
restart evolution, while the short digest keeps AthenaK's parameter header below its
fixed read limit.

Density and pressure are interpolated logarithmically.  Velocity, magnetic field, and
all gradients are interpolated linearly; interpolation of two trace-free magnetic
gradients remains trace-free.  Initial magnetic faces still come from the discrete curl
of the analytic vector potential, and replay only prescribes user ghost faces, leaving
active faces under CT evolution.  Consequently the active-grid divergence constraint is
preserved even when the upstream field changes in time.

This divergence statement is numerical, not a claim that arbitrary snapshot-to-snapshot
changes satisfy Faraday's law at the boundary.  The coefficient replay does not import a
global edge EMF.  Use a dump cadence that resolves the shortest relevant eddy/advection
time, check convergence after halving that cadence, and monitor boundary Poynting flux.
Rapid magnetic variability requires the later face-flux-plus-EMF worldtube path.

Every face marked `user` is prescribed by the analytic profile.  Mark only faces that
remain upstream over the complete table and verify that the fitted normal velocity never
reverses there.  If turbulence changes which face supplies the domain, the current fixed
face contract is insufficient; making every face prescribed would also corrupt the
downstream wake.  A characteristic inflow/outflow switch should be implemented before
using such a dataset.

This coefficient replay deliberately avoids pointwise interpolation between the global
and local meshes.  Each global AMR leaf cell contributes with its physical volume to a
local source-tetrad fit, and the small simulation evaluates that continuous fit on its
own mesh.  If the local box spans enough global cells that first-order residuals are no
longer small, the next extension should replay a spatial worldtube (preferably vector
potential or face flux plus EMF), not silently increase the polynomial order.

## Cubical face-flux/edge-EMF worldtube

`worldtube_flux_emf.py` implements the topology and host-side transfer layer for the
spatial replay.  Its strict NPZ container is classified as
`athenak-emri-cubical-flux-emf-worldtube-v1` and stores, on each of the six outward
oriented faces:

```text
cell_state[t, variable, v, u]
normal_flux[t, v, u]
emf_u[time_interval, v_edge, u_segment]
emf_v[time_interval, v_segment, u_edge]
```

`normal_flux` is a face-cell integral, not a point value of `B`.  Each EMF is a line
integral oriented along the face-local direction, and `u cross v` is the outward normal.
The current outer writer declares and records `cell_state` in the order
`rho,u1,u2,u3[,pgas][,scalars...],bcc1,bcc2,bcc3`.  These final three cell-centered
samples supply the two tangential fields needed by an MHD characteristic projection;
the CT-normal field is still taken from `normal_flux`, not from `bcc`.
The EMF over one stored interval is the time average actually used by the outer CT
update.  Therefore every face cell must satisfy

```text
Phi[n+1] = Phi[n] - dt * (Eu_bottom + Ev_right - Eu_top - Ev_left).
```

The validator checks this relation on every cell and interval, pairs the two copies of
all twelve cube-edge EMFs with the correct orientation, and checks zero net outward flux
over the closed cube.  It rejects a file before any inner evolution if one of these
topological identities fails:

```bash
python3 inputs/emri/worldtube_flux_emf.py validate outer_worldtube.npz
```

Different outer and inner face resolutions are handled without pointwise `B`
interpolation:

```bash
python3 inputs/emri/worldtube_flux_emf.py resample \
  outer_worldtube.npz inner_worldtube.npz --cells-per-face 128
```

Normal fluxes are transferred by exact area overlap; line EMFs use conservative overlap
along their edge and nodal interpolation transverse to it.  These tensor-product
operators commute with the discrete face curl, including non-integer changes in
resolution.  Fluid face-cell values use area-average transfer.  Regridding preserves
both total magnetic flux and the stored Faraday update to roundoff, but it assumes the
same physical cube; changing the matching radius is a new outer extraction, not a grid
interpolation.

AthenaK now exposes two optional problem hooks at the required RK-stage locations.
`user_efld_func` runs after built-in corner-EMF sources but before EMF communication; an
inner injector writes there so block and AMR synchronization still see its data.
`user_efld_observer_func` runs after `RecvE` and immediately before CT; an outer writer
reads there so it records the synchronized EMF actually used by the update.  The writer
must accumulate those `mhd::MHD::efld` values with the integrator recurrence.  Neither
`mhd_w_bcc` snapshots nor `-v cross B` reconstructed afterward are equivalent.  In
particular, constructing some `E=-dA/dt` from snapshot changes can satisfy a continuous
Faraday equation yet produce a large first-cell artifact when stitched to the inner
Riemann EMF.  That shortcut is deliberately not provided.

The current container uses one interval-average EMF, which permits conservative
piecewise-constant temporal replay at different local timesteps.  Higher-order temporal
replay will require outer RK-stage samples or temporal moments in a later schema.  The
implemented inner driver loads bounded slabs and injects stored line EMFs in
`user_efld_func` at every local RK stage.  It also has validated flat-spacetime and
ADM-face GRMHD seven-wave incoming-mode boundaries.  The offline global-snapshot path
below now performs the upstream state/two-form transformation.  An online moving
cut-surface writer remains preferable when exact source RK EMFs are available.

### First AthenaK outer writer

The pgen-independent outer writer is enabled by an optional input block, so it can be
attached to a global single-black-hole GRMHD problem without modifying that problem's
initial-data code:

```text
<emri_worldtube>
enabled       = true
center_x1     = 20.0
center_x2     = 0.0
center_x3     = 0.0
half_width    = 4.0
dcycle        = 8
overwrite     = false
file_basename = disk_patch.outer_worldtube
```

At RK stage one it records the initial face flux and exterior-adjacent primitive state.
After EMF synchronization at every stage it advances a device-resident line-integral
register with

```text
I_stage = gamma0_stage * I_previous + beta_stage * dt * E_stage * edge_length.
```

At each `dcycle` endpoint it sums the completed-step integrals, writes their interval
average, and writes the new state and normal flux.  A final off-cadence interval is
flushed by the driver finalize hook.  MPI ranks contribute disjoint surface cells/edges
and reduce only the small surface arrays to rank zero; the full volume is never copied to
the host.  Output is an incrementally updated JSON manifest plus little-endian float64
streams.  Convert and validate it with

```bash
python3 inputs/emri/worldtube_flux_emf.py pack-outer \
  disk_patch.outer_worldtube.cycle00000000.manifest.json disk_patch.worldtube.npz
```

This first writer deliberately enforces a fixed, grid-aligned cube on an isotropic,
single-level Cartesian mesh and supports explicit RK1/RK2/RK3 without level subcycling.
It rejects AMR/SMR, a moving cube, non-aligned faces, or an existing output segment
unless `overwrite=true`.  These are correctness gates, not fundamental limitations.
A continuously moving worldtube requires the motional edge term and a discrete geometric
conservation law; merely changing `center_x*` each output would break the same Faraday
identity the format is designed to preserve.  Static refinement can be added next by
extracting on one uniform leaf level and using the existing mimetic regridder.

### Moving/source-tetrad frame contract

`worldtube_frame.py` supplies the geometry and discrete-constraint layer used by the
offline cut-surface sampler below.  For the affine local map

```text
x^mu(T,X) = z^mu(T) + e^mu_a(T) X^a,
```

it constructs the full Jacobian, including the position-dependent time leg
`dz^mu/dT + (de^mu_a/dT) X^a`, and pulls back the Faraday two-form.  When the global and
local time coordinates agree, the electric edge one-form reduces to

```text
E_local,a = e^i_a [E_i + (w cross B)_i].
```

Thus both the translation of the secondary and rotation of its spatial tetrad enter the
EMF.  The same module transforms contravariant four-vectors with the inverse spacetime
Jacobian.  Magnetic components come from the spatial part of the pulled-back two-form,
so no separate, ambiguous vector-component rotation is used.

An instantaneous frame contract can be audited before a run:

```json
{
  "classification": "athenak-emri-affine-worldtube-frame-v1",
  "worldline_tangent": [1.0, 0.0, 0.2, 0.0],
  "spatial_legs": [[0,0,0], [1,0,0], [0,1,0], [0,0,1]],
  "spatial_leg_derivative": [[0,0,0], [0,-0.01,0], [0.01,0,0], [0,0,0]]
}
```

```bash
python3 inputs/emri/worldtube_frame.py audit frame.json
```

Only a stationary signed permutation of the Cartesian axes is an exact relabeling of
the current fixed cube.  A translating axis-aligned cube needs an ALE surface sampler;
a general tetrad rotation needs cut-face/cut-edge geometry as well.  The audit reports
these cases separately and never treats ordinary array interpolation as a tensor
transformation.

Even a correct moving-surface quadrature will generally leave small inconsistencies
between interpolated endpoint fluxes and interpolated edge EMFs.  Given a raw NPZ with
the usual arrays, apply the unique-edge surface-cochain projection with

```bash
python3 inputs/emri/worldtube_frame.py project-moving \
  sampled_moving_worldtube.npz projected_worldtube.npz
```

The projection first merges the two samples of each of the twelve physical cube-edge
segments in a common orientation.  On every interval it then finds the minimum-norm
edge correction whose discrete curl equals `(Phi[n]-Phi[n+1])/dt` on all six face
grids.  Endpoint fluxes and fluid states are not changed.  A changing net flux through
the closed surface is rejected because no edge field can repair it.  The output records
the raw Faraday error, shared-edge mismatch, iteration residual, and maximum/RMS edge
correction.  These corrections must decrease under outer spatial and temporal
refinement; the projection restores topology but cannot rescue an under-resolved
physical sample.

For an aligned-spin circular equatorial Kerr orbit, generate the full Hermite table
directly on the binary-dump time axis:

```bash
python3 inputs/emri/build_kerr_circular_frame.py \
  --snapshot-manifest global_worldtube.json \
  --primary-mass 1.0 --primary-chi 0.7 --orbital-radius 8.0 \
  --orbit-direction 1 --output orbit.frame.json
```

Set the extraction manifest `frame` field to `orbit.frame.json`.  The radius is the
Kerr-Schild/Boyer-Lindquist radial coordinate; the Cartesian equatorial coordinate
radius is `sqrt(r^2+a^2)`.  The generator uses the same test-particle Kerr frequency
and ISCO convention as `emri_windtunnel`, rejects a superluminal or inside-ISCO orbit,
and constructs an orthonormal radial, prograde-tangential, vertical source tetrad.  It
uses a co-rotating tetrad convention, not a parallel-transported or Fermi-Walker frame.
When a snapshot manifest is supplied, it reads `global_length_in_local_units`, converts
dump times with `T_local=L t_global`, and scales every Jacobian column by `1/L`; this is
required for an extreme mass ratio and is not merely provenance metadata.
The output reports Hermite midpoint worldline error, interpolated tetrad Gram error,
timelike margin, spatial-leg time-row norm, and Jacobian determinant/condition number.
An orthonormal moving tetrad generally has nonzero `e^0_a`, so source dumps must extend
beyond the requested local sample-time endpoints by the resulting worldtube time tilt;
using every source time as an output time is not automatically valid.  The preflight
below calculates this exactly.  Hermite interpolation errors must be reduced by adding
frame knots when the orbital phase advance per dump is large.
Generic inclined, eccentric, or inspiralling worldlines still require an external
ephemeris and transport prescription; they must not be approximated by this circular
generator.

For such a numerical worldline, start from `kerr_ephemeris.example.json`.  Its
`global_times`, Cartesian Kerr-Schild `positions`, and `coordinate_velocities` define a
cubic-Hermite trajectory.  Build the transported frame with

```bash
python3 inputs/emri/build_kerr_ephemeris_frame.py \
  --ephemeris inputs/emri/kerr_ephemeris.example.json \
  --output orbit.ephemeris.frame.json
```

The default `transport_mode=fermi_walker` computes the aligned-spin Cartesian Kerr-
Schild connection by centered metric differences, reconstructs the four-velocity and
four-acceleration, and integrates the spatial legs with fixed-substep RK4.  It is the
appropriate nonrotating convention for an accelerated prescribed worldline.  The
spatial legs are Gram-projected after each RK substep; the maximum projection correction
is recorded and must remain small.  The generator also repeats the transport with half
as many substeps and records the maximum coarse/fine leg difference.  Repeat production
runs with a smaller `metric_fd_step_global_units` as an independent connection-
differencing convergence check.

`transport_mode=parallel` is reserved for a resolved geodesic ephemeris.  It hard-fails
when the dimensionless proper acceleration `M|a|` exceeds
`geodesic_acceleration_tolerance_Ma`; relaxing that threshold does not turn an
accelerated inspiral into a geodesic.  Hermite acceleration can jump at ephemeris knots,
so both its absolute and relative jump are reported.  Refine the ephemeris cadence until
the acceleration jump, interpolated tetrad Gram error, RK coarse/fine difference, and
science observables all converge.

If `initial_spatial_basis` is omitted, the initial legs use cylindrical radial,
prograde-tangential, and disk-normal coordinate directions before metric Gram-Schmidt.
For an inclined or otherwise geometry-specific orbit, provide the desired three basis
columns explicitly; Fermi-Walker transport then fixes their subsequent nonrotating
evolution.  As in the circular generator, all output times and Jacobian columns include
the declared extreme-mass-ratio unit conversion.

This generator assumes a stationary Kerr background whose signed spin is along the
Cartesian Kerr-Schild coordinate z axis; `disk_normal` controls only the default initial
spatial basis.  It is not valid when the global ADM metric contains material time
dependence or a non-Kerr contribution large enough to affect tetrad transport; that case
requires connection data extracted from the actual ADM time series.

For that numerical-geometry case, use the same snapshot manifest that will later drive
the fluid worldtube together with `adm_ephemeris.example.json`:

```bash
python3 inputs/emri/build_adm_ephemeris_frame.py \
  --manifest global_worldtube.json \
  --ephemeris inputs/emri/adm_ephemeris.example.json \
  --output orbit.adm.frame.json
```

The builder reads only the ADM variables from each binary pair; the state dump is opened
for its time and MeshBlock topology but its fluid arrays are not loaded.  It applies the
same fixed-leaf-level and trilinear spatial interpolation rules as the fluid extractor.
At every trajectory event it samples the metric at the center and at three centered
offset pairs.  The default offset is one quarter of the smallest selected-level source
cell; set `adm_metric_fd_step_global_units` to override it.  A second connection estimate
at half that offset is reported as an independent finite-difference diagnostic.

ADM fields are linear between source snapshots.  Their metric time derivative is
computed analytically rather than by subtracting two nearly coincident times.  The
transport grid is the union of ephemeris knots and all enclosed source-snapshot times,
so an RK4 step never crosses a jump in the piecewise temporal derivative.  The output
reports the connection jump across every such source knot.  Refine the dump cadence
until that jump, the interpolated tetrad Gram error, and the resulting force/accretion
observables converge.

`fermi_walker` remains the default for a prescribed accelerated inspiral.  Numerical ADM
data can make an analytically geodesic ephemeris weakly nongeodesic, so `parallel` still
hard-fails against `geodesic_acceleration_tolerance_scaled`; the multiplier used to make
`|a|` dimensionless is `acceleration_scale_global_units`, normally the primary mass in
global units.

Afterward, replace the extraction manifest's `frame` object with the path to
`orbit.adm.frame.json`, then run the existing metadata-only preflight and worldtube
extractor.  The frame, fluid state, Faraday tensor, face flux, and motional edge EMF are
then all derived from the same ADM/GRMHD series.  The usual restriction remains: every
orbit finite-difference stencil and every worldtube quadrature stencil must stay on the
declared fixed leaf level; neither tool silently interpolates across an AMR interface.

The implementation was exercised on the existing six-snapshot `16^3` DynGRMHD closure
series.  With a two-snapshot cache, chronological fine/coarse transport plus diagnostics
used 12 cache misses (two bounded passes through the series).  The generated tilted
frame passed the exact worldtube preflight, and a subsequent real fluid extraction left
the closed-flux and Faraday residuals at roundoff; its largest relative edge projection
was `1.58e-3`.  A separate level-one AMR series loaded only its selected ADM leaf blocks
and completed with exactly two cache misses for two snapshots.

For a future exact online path, the remaining gap is narrow and explicit: sample `F` and
fluid four-vectors on the moving cut surface during the global RK recurrence, then feed
those raw integrals through this pullback/projection layer.  The fixed-grid observer is
still the only sampler that records the source CT EMF recurrence exactly.

### Offline global-snapshot cut-surface extraction

`extract_global_worldtube.py` provides the complete one-way path for an existing global
GRMHD time series.  Every sample pairs co-temporal `mhd_w_bcc` and `adm` binary dumps.
A level-zero uniform tiling is assembled into a dense Cartesian source grid.  For an
AMR/SMR output, set the manifest `source_level` to a leaf level that covers every sampled
locus, including the exterior-adjacent state centers, plus its one-cell interpolation
halo.  The extractor retains only that level's leaf MeshBlocks, interpolates across
neighboring same-level block boundaries, and rejects any stencil that touches a
coarse-fine interface.  It deliberately does not invent a cell-centered magnetic
prolongation rule from `mhd_w_bcc` dumps.  A globally uniformly refined leaf tiling is
recognized automatically even when its sole physical level is nonzero.

The target cube may have a different spacing, center, translation, or rotation and need
not align with the source cells.  The source and ADM leaf topology may change between
snapshots, provided the selected level continues to cover every sampled stencil.  Each
output records the selected level, all available leaf levels, the number of retained
blocks, and whether dense or sparse source storage was used.  Production should place
the matching surface at least one selected-level source cell away from an AMR interface;
additional margin is advisable when the source refinement region moves between dumps.

The frame table is cubic Hermite data for `z^mu(T)` and `e^mu_a(T)`: it contains their
values and first derivatives at every knot.  At each sampled event the extractor
reconstructs

```text
W = sqrt(1 + gamma_ij u^i u^j),
U^0 = W/alpha,
U^i = u^i - beta^i U^0,
v_transport^i = U^i/U^0,
E_i = -(v_transport cross bcc)_i,
F_0i = -E_i,  F_ij = epsilon_ijk bcc^k,
```

where `bcc^i=sqrt(gamma) B^i`.  It then transforms the spacetime metric, four-velocity,
and Faraday two-form rather than rotating three-component arrays independently.  If one
global length unit contains `L` local units and the optional density renormalization is
`D`, the numerical-unit conversion is

```text
g_local = L^2 J^T g_global J,
U_local = J^-1 U_global/L,
F_local = L sqrt(D) J^T F_global J,
(rho,p)_local = D (rho,p)_global/L^2.
```

For the common map `x_global=X_local/L`, this leaves a Minkowski metric unchanged and
gives `B_local=sqrt(D) B_global/L`.  This separation is important: treating the
mass-ratio unit conversion as only a coordinate pullback would introduce an extra
erroneous power of `L`.

The fluid state is sampled at the exterior-adjacent local cell center, exactly matching
the operational outer-writer convention.  Magnetic flux is Gauss-integrated from the
spatial pulled-back two-form on the geometrical face.  The time-edge part of the same
two-form is Gauss-integrated along each edge and over each stored time interval, so
translation and rotation automatically contribute the motional EMF.  Finally, the tool
performs two small topological projections:

1. subtract the equal-area minimum-norm mean face flux at each endpoint so the closed
   cubical surface has zero magnetic monopole;
2. minimally correct the unique edge cochain so its curl equals the endpoint flux
   change on every face cell.

Both relative corrections are written to the output metadata.  They are error
estimators, not free accuracy: repeat the extraction with finer source dumps, target
faces, temporal cadence, and quadrature, and require the corrections and observables to
converge.  The complete frame contract, resolved source paths, sizes, and manifest hash
are retained as provenance.  Set `hash_source_files=true` for a final production product
to add SHA-256 for every large dump; it defaults to false so exploratory extraction does
not perform a second full read of the source series.

Long series are scanned once in metadata-only mode and loaded through a bounded LRU
cache.  `snapshot_cache_size` defaults to two, the minimum needed for linear time
interpolation.  Sampling is ordered by time interval across all six faces, so a
time-orthogonal frame whose output times follow the dump times normally loads every
source snapshot exactly once.  A frame with temporal tilt `e^0_a X^a` can span several
source intervals on one face; increase the cache if the recorded miss and eviction
counts show repeated loads.  The output metadata records capacity, hits, misses,
evictions, and peak resident snapshots.

Consequently, resident source data scale as
`O(snapshot_cache_size * selected-level snapshot size)` rather than with the number of
dumps.  This does not yet make a single dump block-streaming: the binary reader loads
the requested state and ADM variables for one complete snapshot pair before the
extractor keeps only the selected leaf-level blocks.  Production memory estimates must
therefore include that one-pair transient peak in addition to the cache and output
worldtube.

Start from `global_worldtube_manifest.example.json`, replace its dump paths and affine
frame with the physical worldline/tetrad map.  Before loading any fluid arrays, traverse
the exact state-center, face-quadrature, edge-quadrature, and temporal-quadrature loci:

```bash
python3 inputs/emri/preflight_global_worldtube.py \
  --manifest global_worldtube.json \
  --output disk_patch.preflight.json
```

The preflight reads only headers and MeshBlock topology.  It checks both temporal
interpolation endpoints for every mapped event, including temporal tetrad tilt, and
hard-fails on the source cell-center envelope, a coarse-fine stencil, or a singular
affine Jacobian.  `minimum_additional_stencil_halo_cells=0` is accepted but warned: the
requested stencil is exactly covered with no spare selected-level cell.  Production
should normally retain at least one additional cell, and more when the AMR region moves
between dumps.  The reported envelope margin is also measured in selected-level source
cells.

After a passing preflight, run

```bash
python3 inputs/emri/extract_global_worldtube.py \
  --manifest global_worldtube.json \
  --output disk_patch.worldtube.npz \
  --diagnostics disk_patch.extraction.json
python3 inputs/emri/worldtube_flux_emf.py validate disk_patch.worldtube.npz
python3 inputs/emri/worldtube_flux_emf.py prepare-inner \
  disk_patch.worldtube.npz disk_patch.inner.bin
```

All mapped event times must lie inside the source snapshot range, including any temporal
offset `e^0_a X^a`, and all mapped spatial points must lie inside the source cell-center
envelope.  The latter is stricter than merely lying inside the domain because binary
dumps do not contain source ghost cells.  The output declares
`fluid_state_frame=inner_coordinate`, so it can be used with
`fluid_boundary=characteristic_gr`; the inner total metric should nevertheless approach
the transformed global metric at the replay boundary, otherwise the one-way matching
surface is too close to the secondary.

On the standard `16^3` DynGRMHD closure snapshot pair, a static `8x8` face extraction
reconstructed the exterior state to `2.55e-8` relative L2 (the dump is float32), face
flux to `1.61e-5`, and the exact-RK writer EMF to about `8e-3`.  Its closed-flux
correction was `2.22e-17` relative and the maximum Faraday edge correction was
`2.90e-5` in line-integral units.  The state/flux agreement validates coordinates and
GR variable reconstruction; the larger EMF difference quantifies why sparse snapshots
are an approximate fallback rather than a replacement for an online RK observer.
Replaying this offline product through `characteristic_gr` completed 1,536 mode
projections with zero fallback and a `3.42e-18` CT endpoint residual.  Relative to the
coincident global subvolume, its density, pressure, velocity-vector, and magnetic-vector
L2 errors were `1.84e-5`, `6.78e-5`, `3.94e-5`, and `3.35e-5`.  The first three retain
the online closure accuracy; the larger magnetic error isolates the snapshot-EMF
approximation that must be converged with dump cadence.

A real two-level DynGRMHD smoke output provides a separate AMR addressing test.  The
source has 56 level-zero and 64 level-one leaf blocks; a `12^3`, half-width-three
worldtube and its halo lie inside the centered fine region.  Relative to a uniform
`32^3` source with the same local `dx=0.5`, the sparse level-one extraction differs by
`4.69e-8` in state L2, `7.40e-9` in face-flux L2, and `9.57e-5` in EMF L2 after one
`1e-4` step.  The maximum final edge-projection residual is `1.20e-18`.  Inner replay
performs 3,456 characteristic projections with zero fallback and an `8.65e-19` boundary
CT residual; density, pressure, velocity-vector, and magnetic-vector closure errors are
`1.03e-7`, `3.11e-7`, `2.35e-7`, and `1.71e-7`.  Moving a stencil across the fine-region
boundary is a hard error rather than an implicit coarse-fine interpolation.
The metadata-only preflight visits 16,128 local quadrature events.  It measures nine
cells of source-domain envelope margin but zero additional level-one halo cells for the
deliberately tight AMR regression region, and therefore emits the expected production
margin warning; the equivalent uniform source has eight additional cells.

Reproduce the complete AMR source, equivalent uniform source, extraction comparison,
and inner replay with

```bash
python3 inputs/emri/run_global_worldtube_amr_check.py \
  --athena build_emri/src/athena \
  --workdir runs/emri/global_worldtube_amr_check \
  --fail-on-gate
```

The default check uses the optional centered `refined_region1` in
`emri_windtunnel_smoke.athinput`; ordinary smoke runs retain `refinement=none`.

### Offline extraction convergence campaign

`analyze_global_worldtube_convergence.py` compares an offline NPZ against the online
outer-writer stream.  Candidate endpoint times may be a subset of the online times.  In
that case the reference edge EMF is not subsampled: each group of fine intervals is
combined with its exact `dt`-weighted temporal mean.  Consequently the coarsened
reference retains the integrated discrete Faraday update and shared-edge orientation.

The complete CPU campaign can be reproduced with

```bash
python3 inputs/emri/run_global_worldtube_convergence.py \
  --athena build_emri/src/athena \
  --workdir runs/emri/offline_convergence \
  --outer-cells 16 24 32 \
  --cadence-strides 1 2 4 \
  --quadrature-orders 1 2 3 \
  --replay-quadrature 2 \
  --tlim 0.08 --cfl 0.02 --fail-on-gate
```

For every source resolution, the driver evolves one DynGRMHD reference, writes
co-temporal `mhd_w_bcc` and ADM dumps at every step, and records the exact online RK
worldtube.  It then extracts all cadence/quadrature combinations, performs the two CT
projections, compares against the conservatively coarsened online stream, replays the
second-order products through `characteristic_gr`, and writes both per-case JSON and a
global `summary.json`.  The input is a smooth short-duration regression problem; the
reported gates are software acceptance criteria, not universal physics tolerances.

The standard 27-case CPU result gives the following finest-cadence errors and extraction
costs (wall time depends on the host):

| source cells | Gauss order | flux relative L2 | EMF relative L2 | CPU seconds |
|---:|---:|---:|---:|---:|
| 16 | 2 | 5.61e-5 | 7.97e-3 | 1.93 |
| 16 | 3 | 4.75e-5 | 5.32e-3 | 4.04 |
| 24 | 2 | 2.76e-5 | 3.44e-3 | 6.48 |
| 24 | 3 | 2.37e-5 | 2.29e-3 | 13.79 |
| 32 | 2 | 1.59e-5 | 1.91e-3 | 14.02 |
| 32 | 3 | 1.36e-5 | 1.26e-3 | 29.91 |

At third order the measured adjacent spatial orders are `1.72, 1.93` for face flux and
`2.08, 2.06` for EMF.  Cadence strides through four steps keep the EMF-error change
below `1.5e-3` relative in this test; its largest snapshot interval is only about
`0.066 dx` in light-crossing units, so this result must not be extrapolated to sparse
disk dumps.  Across all nine second-order replays, 43 of 226,037 characteristic
projections used the safe exterior-state fallback (`1.90e-4`), while the maximum CT
boundary residual was `1.37e-17`.  The worst density, pressure, velocity-vector, and
magnetic-vector closure errors were `7.22e-5`, `2.65e-4`, `1.56e-4`, and `1.34e-4`.

For a global-disk production series, start with the following policy:

1. Keep the matching cube on one fixed source level and use target spacing no finer than
   the source spacing.  Repeat with both spacings reduced together and demand convergence
   of the force, accretion rate, magnetic flux, and variability statistics used in the
   paper.
2. Choose dumps so the largest local characteristic displacement between snapshots is
   initially below about `0.1 dx`; then halve the interval.  The halved-cadence change in
   science observables, rather than the smoke-test field error, is the acceptance test.
3. Use second-order Gauss integration for the broad survey and repeat representative
   epochs at third order.  Here third order reduces the finest EMF error by one third but
   costs about twice as much as second order.
4. Require closed-flux and final Faraday residuals near roundoff.  Treat the true
   relative edge correction (correction divided by the largest raw sampled edge EMF) as
   an interpolation-error estimator: it should decrease under refinement and remain
   below the tolerance dictated by magnetic-force observables.  The smooth campaign
   remains below one percent.
5. Track the characteristic fallback fraction and its spatial distribution.  A small
   global fraction is insufficient if fallback cells form a coherent patch near the
   downstream wake or a current sheet.

### Inner CT magnetic replay

Prepare the validated NPZ for bounded-memory C++ loading:

```bash
python3 inputs/emri/worldtube_flux_emf.py prepare-inner \
  disk_patch.worldtube.npz disk_patch.inner.bin
python3 inputs/emri/worldtube_flux_emf.py inspect-inner disk_patch.inner.bin
```

The binary has a fixed little-endian header, exact declared dimensions, and a CRC32 over
the entire payload.  The C++ reader verifies all of these, the strictly increasing time
array, MHD variable count, coincident cube center, and exact cube/grid match before
changing a field.  `prepare-inner` only repacks already transformed data; it does not
translate or rotate a global worldtube.  The driver keeps
only the two endpoint state slabs, two endpoint flux slabs, and one interval EMF slab on
the device; the full time series remains on disk.  Enable it with

```text
<emri_worldtube>
enabled        = true
mode           = inner
file           = /absolute/path/disk_patch.inner.bin
time_offset    = auto
flux_tolerance = 1.0e-10
fluid_boundary = off
characteristic_speed_tolerance = 1.0e-10
fluid_state_frame = unverified
```

On a fresh run the replay initializes the six active boundary-normal face fields from
the first flux slab.  At each local RK stage it replaces every physical-boundary edge
EMF before communication, including every MeshBlock copy of an edge.  The ordinary CFL
timestep is clipped so no step crosses an outer-data endpoint.  An interval-average EMF
can therefore drive any number of smaller inner timesteps while preserving the exact
integrated CT update.  At every outer endpoint the code compares the evolved boundary
flux to the stored next slab and aborts if the normalized discrepancy exceeds
`flux_tolerance`.  Single-block, eight-block, equal-timestep, and five-inner-step per
outer-interval tests close this check.

The state slabs are loaded but remain inactive by default.  Use the offline extractor
above or a future equivalent online transformer before enabling a fluid boundary;
copying global coordinate primitives into ghost cells uses the wrong frame, and a moving
cube also needs the motional EMF.  Set `fluid_boundary=riemann` only after `cell_state`
has the inner-coordinate primitive convention and trailing `bcc1,bcc2,bcc3` fields.
This enables the first operational fluid replay:
after every C2P cache refresh it linearly interpolates the exterior state in worldtube
time, fills the face-normal ghost columns, replaces the ghost normal `bcc` by the CT
flux value, and lets the ordinary GRMHD HLLE solver pose the boundary Riemann problem.
It never overwrites active cells or the evolved CT face field.

This HLL boundary is characteristic-compatible in the two-wave sense.  If the entire
HLL fan leaves the box, the solver uses the interior flux exactly; if the entire fan
enters, it uses the exterior flux; a mixed fan gets the standard HLL combination.  It
is substantially safer than overwriting active primitives, but it is more diffusive
than the seven-wave projector below and does not separately preserve every outgoing
Alfvén/contact/slow amplitude.  `fluid_boundary=off` remains the default until reflection
and frame-transformation tests are passed.

`fluid_boundary=characteristic_sr` enables the flat-spacetime device implementation;
`fluid_boundary=characteristic_gr` applies the same analytic seven-wave projector in an
Eulerian orthonormal face frame constructed from the current ADM metric.  Both first
linearly extrapolate the active and next-inward cells to preserve every outgoing spatial
gradient, then replace only incoming amplitudes with exterior worldtube data.  Deeper
ghost layers continue the projected first-ghost slope.  Density and pressure are
positivity limited, and every failed/degenerate basis or non-timelike fixed-coordinate
face falls back to the full exterior state and is counted in the final log.  Startup
regressions compare the analytic RMHD basis with the independent numerical-Jacobian
reference and the GR frame with a non-diagonal-metric reference.

The GR option is deliberately opt-in twice:

```text
fluid_boundary    = characteristic_gr
fluid_state_frame = inner_coordinate
```

The second line asserts that the stored `rho,u^i,pgas,sqrt(gamma)B^i` values have already
been transformed into the inner simulation coordinates.  `prepare-inner` validates and
repacks topology but does not perform that tensor transformation, so leaving the value
as `unverified` aborts before evolution.

### Incoming-characteristic reference projector

`worldtube_characteristics.py` fixes the mathematical convention for the remaining
fluid boundary.  At each face point the global/source-frame state must first be mapped
to a local orthonormal frame whose first spatial axis is the outward normal.  In that
frame the independent ideal-RMHD primitive vector is

```text
q = (rho, pgas, u_normal, u_tangent_u, u_tangent_v,
     B_tangent_u, B_tangent_v),
```

while `B_normal` is supplied by the CT face flux and is not an independently prescribed
wave variable.  The module differentiates the SRMHD conserved variables and normal
fluxes, constructs the seven-wave primitive Jacobian, and projects
`q_external-q_interior` onto its eigenvectors.  With an outward normal, only eigenvalues
`lambda < 0` enter the domain in a zero-shift flat frame.  The boundary state is
therefore

```text
q_predict = 2 q_boundary-cell - q_next-inward,
q_ghost1 = q_predict + R P_incoming L (q_external - q_predict).
```

This is the local-characteristic strategy for GRMHD: spacetime curvature enters through
the boundary tetrad and coordinate wave-speed conversion, while the local wave fan is
special relativistic.  The reference implementation uses complex-step Jacobians and a
dense eigensolve so it is transparent and independently testable.  It recovers the exact
relativistic hydrodynamic sound speeds when `B=0`, returns the complete exterior state
for super-fast inflow, retains the interior state for super-fast outflow, and replaces
only the incoming amplitudes in a mixed fan.

The operational SR and GR paths use the same closed-form RMHD eigenbasis on CPU/GPU.
For coordinate face `x^n=constant`, the GR path constructs

```text
e_normal^i = normal_sign gamma^{in}/sqrt(gamma^{nn}),
B^i = bcc^i/sqrt(gamma),
lambda_relative = lambda_Eulerian
                  - normal_sign beta^n/(alpha sqrt(gamma^{nn})).
```

It injects a mode only when `lambda_relative` is negative, transforms the projected
velocity and field back to coordinate components, re-densitizes the field, and then
sets its normal component exactly from the CT flux.  The two tangent legs remain in the
coordinate face and are Gram-Schmidt orthogonalized, so all six faces share a
right-handed convention even for a non-diagonal spatial metric.  The metric is frozen at
the Riemann face across the ghost reconstruction stencil; its truncation error is
therefore controlled by `dx` relative to the curvature scale.  The host reference
aborts on an ill-conditioned eigenbasis, while the device path records a conservative
exterior-state fallback.  These choices follow the local-Riemann construction of Antón
et al. (2006) and prevent the replay from overconstraining outgoing modes.

Before deciding whether that more expensive projector is necessary for a particular
disk patch, audit the actual packed worldtube rather than relying only on the upstream
Mach number:

```bash
python3 inputs/emri/audit_worldtube_boundary.py disk_patch.inner.bin \
  --gamma 1.3333333333333333 --output disk_patch.boundary-audit.json
```

The audit reports the incoming-mode histogram separately on all six faces, the fraction
of samples with a mixed HLL fan, the numerical eigenbasis condition number, and
`max(|lambda| dt_worldtube/dx_face)` for the fluid sampling cadence.  It also reports a
linear per-mode HLL flux-gain error.  The extremal fast modes and every genuinely
super-fast face have zero error in this local test; slower modes in a mixed fan need not.
This gain error is a screening diagnostic, not a wave-packet reflection coefficient.
Inspect its median, 95th percentile, and maximum together: a very large maximum can be
caused by an almost stationary mode.  Repeat the audit with a physically justified
`--speed-tolerance` to expose that sensitivity, but do not discard a slow mode merely
to improve the reported number.
Any production choice of `fluid_boundary=riemann` must still pass a propagated
outgoing-wave reflection test at the intended resolution and worldtube cadence.

An end-to-end domain-decomposition closure test is now available before that more
selective wave test.  Evolve a larger reference domain while recording a cubical
worldtube, replay its coincident inner cube at identical `dx`, output `mhd_w_bcc` at the
same time, then run

```bash
python3 inputs/emri/compare_worldtube_closure.py \
  outer.mhd_w_bcc.00001.bin inner.mhd_w_bcc.00001.bin \
  --output closure.json
```

For a mode-resolved boundary-reflection gate, build the built-in problem generators and
run the seven native special-relativistic MHD wave families through the same outer to
inner path:

```bash
python3 inputs/emri/run_worldtube_reflection_campaign.py \
  --athena build_worldtube_wave/src/athena \
  --workdir /tmp/emri-worldtube-reflection \
  --modes 0 1 2 3 4 5 6 --inner-cells 8 --cfl 0.3 --dcycle 1 \
  --fluid-boundary characteristic_sr
```

The campaign compares the replayed inner cube with the coincident part of the periodic
outer run, projects their primitive-state difference onto the seven local RMHD
characteristics, and reports the modes whose x1 speed is opposite to the injected wave.
Its `reflected_amplitude_coefficient` is an amplitude-like closure norm in a
Euclidean-normalized primitive eigenbasis.  It is deliberately not labeled as an
energy reflection coefficient; that interpretation requires a physical symmetrizer.
The binary plot dumps are float32 by format design, so the default perturbation amplitude
is `1e-3`, well inside the linear regime but safely above output quantization.

For the standard test background at eight inner cells per axis, exterior-adjacent
sampling plus `riemann` gives reflected-amplitude coefficients of 3.04%, 1.91%, 0.276%,
machine zero, 1.04%, 3.53%, and 7.21% for modes 0 through 6.  The seven-wave option gives
0.741%, 0.363%, 0.167%, machine zero, 0.264%, 0.608%, and 2.36%.  The worst mode 6 then
falls to 1.10% and 0.619% at 12 and 16 cells per axis, corresponding to empirical orders
1.89 and 1.99.  These are primitive-eigenmode amplitude diagnostics, not energy
reflection coefficients.

The tool assembles arbitrary fixed-level MeshBlock tilings, extracts the geometrically
coincident reference subvolume, and reports componentwise plus vector-group relative
norms.  The CPU regression exercises a one-block `16^3` outer run against an eight-block
`8^3` inner replay for five RK3 steps.  This measures the total decomposition error; a
mode-resolved reflected wave amplitude remains a separate validation gate.

The curved-spacetime closure campaign automates the corresponding DynGRMHD A/B test:

```bash
python3 inputs/emri/run_worldtube_gr_closure.py \
  --athena build_emri/src/athena \
  --workdir /tmp/emri-worldtube-gr-closure
```

On the standard full effective-EMRI metric, the `16^3` reference versus eight-block
`8^3` replay gives relative L2 errors

| boundary | density | pressure | velocity vector | magnetic vector |
|---|---:|---:|---:|---:|
| `riemann` | 3.48e-5 | 1.03e-4 | 5.75e-5 | 8.90e-6 |
| `characteristic_gr` | 1.84e-5 | 6.77e-5 | 3.94e-5 | 8.89e-6 |

The characteristic run performs 1,536 projections with zero fallback, while the CT
endpoint residual remains approximately `3e-18`.  This is a curved-background domain
decomposition closure result, not yet a physical reflected-energy coefficient.

For that mode-resolved gate, `worldtube_linear_wave.athinput` uses AthenaK's native
seven-wave SRMHD problem in a flat local orthonormal frame.  The worldtube transport is
problem-independent, so the same outer writer and inner CT/fluid replay can be tested
without black-hole curvature mixing the manufactured modes.  Inner replay supports
SRMHD, fixed-background GRMHD, and DynGRMHD; production EMRI runs still use the latter.

## Validation gates before science production

1. Recover the isolated Kerr BHL result as `primary_mass/orbital_radius -> 0` at fixed
   secondary parameters.
2. Check convergence with `metric_fd_step`, box size, resolution, and excision threshold.
3. Verify that the secondary-centered metric is stationary to roundoff and that its
   primary-only curvature agrees with Kerr at the orbital point.
4. Separate horizon momentum flux from near-wake and unresolved far-wake forces; a local
   domain alone does not contain the full dynamical-friction Coulomb logarithm.
5. Compare uniform wind, controlled density/shear gradients, and global-GRMHD boundary
   replay before mapping the measured four-force to changes in `E`, `Lz`, and `Q`.

## Force diagnostics

Set `problem/user_hist=true` and add a history output with
`user_hist_only=true`.  The online diagnostics follow the decomposition used in
arXiv:2201.11753 and arXiv:2409.12359, but retain all three force components and add a
relativistic source-force estimator.  The default `force_frame=source_tetrad` defines
the extraction and cutoff radii by

```text
rhat^2 = sum_a (theta^a_i x^i)^2,
```

where `theta^a_i` is the spatial source co-basis at the orbital anchor.  A source-frame
sphere is generally an ellipsoid in the raw rotating coordinates.  The geodesic-grid
sampling points and its covariant surface-area vector are both transformed; changing
only the points is not sufficient.  The accretion rate is reported per source proper
time, and the momentum flux uses the complete covariant four-momentum flux before
projection:

```text
mdot_hat   = -(dt/dtau) integral rho u^i dSigma_i,
Fmom_hat_a =  (dt/dtau) e_a^mu integral T^i_mu dSigma_i.
```

The `mu=0` energy-flux term in the second expression is essential: omitting it creates a
large spurious tangential force in a rotating chart.  Set `force_frame=coordinate` only
for legacy coordinate-sphere diagnostics.

The paper-compatible far-field estimator is

```text
Fnewt_hat_a = m integral (rho-rho0) xhat_a/rhat^3 d^3xhat,
```

over `r_s < r < force_outer_radius_3`.  Subtracting the uniform background removes a
finite-box cancellation error and is the default.  Set
`force_subtract_background=false` to use `rho` literally as in the cited work.  The
factor of secondary mass `m`, implicit when `G=m=c=1`, is retained here.

The relativistic estimator treats the secondary position as a source parameter:

```text
Frel_i(R) = -1/2 integral T^{mu nu} d_i hsec_{mu nu} sqrt(-g) d^3x,
```

where `d_i` differentiates the secondary Kerr-Schild perturbation with respect to
field-point displacement while holding the rotating-coordinate Jacobian fixed.  Thus
`-d_i hsec` is the derivative with respect to the secondary position.  This construction
isolates the force associated with the small hole: derivatives of the primary metric and
of the rotating basis are not included.  It also includes gas pressure and magnetic
stress, whereas `Fnewt_i` uses density alone.  Its weak-field limit is
`m integral rho x_i/r^3 dV`.

For source-frame output, the generalized-force covector is converted to force per source
proper time and projected on the source tetrad.  `Frel1H_i`, `Frel2H_i`, and `Frel3H_i`
are accumulated to the three configured physical outer radii.  The corresponding total
force estimates are formed in post-processing:

```text
FtotalkH_i = -FmomH_i + FrelkH_i.
```

The axes at the reported orbital phase are `x` radial from the primary through the
secondary, `y` prograde tangential, and `z` normal to the orbital plane.  Consequently,
`FtotalH_y` is the leading orbital-energy/angular-momentum drag channel, while `x` and
`z` can drive eccentricity and inclination.  For a differently oriented imposed wind,
use the full vector rather than identifying a component by name.

Use `force_surface_radius >= 3m` as a starting point and repeat the calculation at other
extraction radii: arXiv:2409.12359 showed that horizon-adjacent momentum fluxes can be
contaminated by density floors.  The extraction sphere must enclose the secondary
horizon in the selected frame, and its coordinate ellipsoid must fit inside the box.
Dependence on `force_outer_radius_*` measures the missing far wake and is a physical
systematic of a local wind tunnel, not merely a numerical error.  The force remains a
quasi-local diagnostic in a prescribed effective metric; it is not fed back into the
orbit.  Also converge `force_surface_nlevel`, which controls the geodesic-grid angular
quadrature independently of the Cartesian fluid resolution.

The history columns are

```text
mass_ratio orbit_r_M omega_M mdot_hat
FmomH_x..z FnewtH_x..z
Frel1H_x..z Frel2H_x..z Frel3H_x..z geo_resid
```

The shorter `H` spelling is used because AthenaK history labels are limited to ten
characters.  Coordinate-frame output retains the un-hatted legacy column names.

## Force averaging and controlled contrasts

`analyze_force_history.py` forms all three
`Ftotal_k=-Fmom+Frel_k` estimates, the incremental outer-wake contributions
`Frel2-Frel1` and `Frel3-Frel2`, trapezoidal time averages, and equal-duration block
standard errors.  With all three conventional labels it also forms the three contrasts
automatically:

```sh
python3 inputs/emri/analyze_force_history.py \
  full=runs/emri/full/full.user.hst \
  frame_only=runs/emri/frame_only/frame_only.user.hst \
  isolated=runs/emri/isolated/isolated.user.hst \
  --tmin 100 --blocks 8
```

Use `--quantity Ftotal3H_y` repeatedly for source-frame histories (or `Ftotal3_y` for
legacy coordinate histories), and `--format csv` or `--format json` for downstream
analysis.  The script detects either schema and rejects a mixed-schema comparison.  The
default averaging interval is the common time overlap.  Duplicate times from appended
restart histories retain their last value.  A block error is meaningful only after the
flow is statistically stationary and each block is longer than the force autocorrelation
time; the script caps the number of blocks at the number of sampled time intervals and
reports `nan` when only one block is possible.

## Decision gate for the next science runs

At a local scale of order the secondary mass, true primary-curvature effects are
parametrically small:

```text
Omega*m ~ q/(r/M)^(3/2),
|Riemann_primary|*m^2 ~ q^2/(r/M)^3.
```

Large differences seen in a `q=10^-3` smoke test therefore must not be extrapolated to a
LISA EMRI; metric components can differ at order `M/r` even when the locally measurable
tidal field is tiny.  Treat a primary-metric signal as resolved only if

```text
|mean(full-frame_only)| > max(2*block_error,
                              resolution_change,
                              outer-radius_change).
```

The one-step CPU consistency check at `q=10^-5`, `r/M=10` illustrates the required
resolution test.  With the default source-tetrad embedding, state, and diagnostics,
increasing the uniform grid from `16^3` to `64^3` reduced the initial
`full-frame_only` contrast in `Ftotal3H_x` from about `8.7e-1` to `6.9e-3`, while the
individual force was about `18.4`.  The large coarse-grid contrast was therefore
coordinate-grid sampling of differently shaped physical ellipsoids, not a detected
primary-curvature force.  These startup values are validation numbers, not stationary
BHL results, but they show that realistic EMRI tidal corrections can lie far below the
truncation floor of a cheap local run.

The economical production sequence is:

1. Evolve one physically motivated, mostly tangential disk-wind setup in all three modes
   until the wake has crossed the largest force radius, then establish a stationary
   averaging window.
2. Repeat `full` and `isolated` at one higher resolution and with a larger box/outer force
   radius.  Converge `force_surface_nlevel` separately because it controls only the inner
   momentum-flux quadrature.
3. Only if the primary-metric contrast clears that numerical/systematic floor, repeat it
   at a genuine EMRI ratio such as `q=10^-5` or `10^-6`.  If it does not, move effort to
   tetrad-matched shear/gradient or global-GRMHD replay boundaries, where the disk
   environment supplies information absent from a uniform wind.

This gate avoids an expensive broad scan whose leading result could be a coordinate or
finite-box effect.  A spin-sign pair for the secondary and a magnetic-field orientation
pair are the next targeted tests for transverse/Magnus-like forces after the baseline is
converged.
