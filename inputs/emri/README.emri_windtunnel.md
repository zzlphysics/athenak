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
`user_efld_func` at every local RK stage.  Its remaining physics step is to inject only
incoming fluid modes after the required global-to-source-tetrad transformation.

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

At RK stage one it records the initial face flux and interior-adjacent primitive state.
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

`worldtube_frame.py` now supplies the geometry and discrete-constraint layer needed by a
future cut-surface sampler.  For the affine local map

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

The remaining implementation gap is consequently narrow and explicit: sample `F` and
fluid four-vectors on the moving cut surface during the global run, then feed those raw
integrals through this pullback/projection layer.  The existing fixed-grid observer is
still the only production outer sampler.

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

The state slabs are loaded but are intentionally not imposed yet.  A production global
disk to local-EMRI coupling must first transform density/pressure, four-velocity, magnetic
flux 2-forms, and EMF 1-forms into the same comoving source-tetrad chart, then prescribe
only incoming GRMHD characteristic amplitudes.  Copying global coordinate primitives
into all ghost cells would both use the wrong frame and overconstrain outgoing modes.
Likewise, a moving extraction cube needs the motional EMF in the discrete map.  The
current magnetic replay is consequently a complete CT/topology transport layer and an
end-to-end test bed, not yet the final physical fluid boundary condition.

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
