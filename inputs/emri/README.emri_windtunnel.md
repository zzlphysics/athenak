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
trace-free linear approximation.  Those limitations are recorded in every manifest and
will be lifted by the time-dependent worldtube stage.

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
