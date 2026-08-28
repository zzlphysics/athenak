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
primary.  The prescribed metric nevertheless contains both Kerr-Schild terms.  At `t=0`
the primary is at global Cartesian Kerr-Schild position `(0,0,0)` and the secondary is at
`(R0,0,0)`, where

```text
R0 = sqrt(orbital_radius^2 + (primary_chi*primary_mass)^2).
```

Local axes rotate at `orbital_omega` and translate with the secondary.  The exact metric
pullback supplies the rotating-frame shift.  Circular equatorial motion and spins aligned
with the orbital z axis make the local metric stationary.  `omega_mode=kerr_geodesic`
uses the test-particle equatorial Kerr frequency; `omega_mode=custom` reads
`orbital_omega` directly.

Use `<adm> dynamic=false` for ordinary global-timestep evolution: the stationary metric
is cached and is regenerated automatically after AMR topology changes.  The repository's
strict `time/subcycling=level` scheduler currently requires `dynamic=true`; that mode
re-evaluates the same stationary metric at RK stages and is therefore more expensive.

This is an effective prescribed spacetime, not a constraint-solved EMRI metric.  It keeps
the full primary Kerr field and all orders in local position present in the superposed
Kerr-Schild approximation, but it neglects inspiral, radiation reaction, primary recoil,
and the disk's gravitational backreaction.  Those assumptions are intentional for a
one-way local wind tunnel.

## Controlled background modes

`background_mode` defines three runs that keep the secondary, wind, magnetic field,
domain, and diagnostic radii fixed:

| mode | primary Kerr term | orbital pullback and secondary boost | role |
|---|---:|---:|---|
| `full` | yes | yes | effective EMRI wind tunnel |
| `frame_only` | no | yes | fixed-chart counterfactual for adding the primary metric |
| `isolated` | no | no | inertial isolated-Kerr BHL control |

The useful contrasts are `full-frame_only` and `frame_only-isolated`.  The first changes
only the primary Kerr contribution while holding the coordinate map fixed.  The second
shows the effect of the orbital non-inertial chart and boost.  In `frame_only`, centrifugal
acceleration is intentionally not balanced by primary gravity, so its nonzero
`geo_resid` is expected.  These are attribution experiments between prescribed
spacetimes, not gauge-invariant observables that may be subtracted without qualification.
In particular, the reported forces are coordinate components of a quasi-local generalized
force.

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

The selected mode is part of the restart contract, so a checkpoint cannot silently be
continued with another background.  Version-1 prototype checkpoints predate this field
and are deliberately rejected.

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

`rho0`, `pgas0`, `u1..u3`, and `b1..b3` initialize a uniform magnetized wind.
The velocity inputs are AthenaK's normal-frame spatial four-velocity components, not
coordinate three-velocities.  Constant face-centered magnetic fields are divergence-free.

Use the built-in `inflow` flag on every upstream face and `outflow` downstream.  The
sample has positive `u1` and therefore uses `ix1_bc=inflow`.  For oblique or sub-fast
flows, more than one inflow face may be physically required.  A future replay boundary
will replace these constants with tetrad-projected data from a global single-SMBH GRMHD
run while preserving constrained-transport face fluxes.

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
relativistic source-force estimator.  On the extraction sphere `r=r_s`, the code computes

```text
mdot    = - integral rho u^j n_j sqrt(-g) r_s^2 dOmega,
Fmom_i  =   integral T^j_i n_j sqrt(-g) r_s^2 dOmega.
```

The paper-compatible far-field estimator is

```text
Fnewt_i = m integral (rho-rho0) x_i/r^3 sqrt(gamma) d^3x,
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

`Frel1_i`, `Frel2_i`, and `Frel3_i` are accumulated to the three configured outer radii.
The corresponding total force estimates are formed in post-processing:

```text
Ftotal_i(Rk) = -Fmom_i + Frelk_i.
```

The axes at the reported orbital phase are `x` radial from the primary through the
secondary, `y` prograde tangential, and `z` normal to the orbital plane.  Consequently,
`Ftotal_y` is the leading orbital-energy/angular-momentum drag channel, while `x` and `z`
can drive eccentricity and inclination.  For a differently oriented imposed wind, use
the full vector rather than identifying a component by name.

Use `force_surface_radius >= 3m` as a starting point and repeat the calculation at other
extraction radii: arXiv:2409.12359 showed that horizon-adjacent momentum fluxes can be
contaminated by density floors.  The extraction sphere must enclose the boosted secondary
horizon, and all outer radii must fit inside the largest origin-centered sphere in the
box.  Dependence on `force_outer_radius_*` measures the missing far wake and is a physical
systematic of a local wind tunnel, not merely a numerical error.  The force remains a
quasi-local diagnostic in a prescribed effective metric; it is not fed back into the
orbit.  Also converge `force_surface_nlevel`, which controls the geodesic-grid angular
quadrature independently of the Cartesian fluid resolution.

The history columns are

```text
mass_ratio orbit_r_M omega_M mdot
Fmom_x..z Fnewt_x..z
Frel1_x..z Frel2_x..z Frel3_x..z geo_resid
```

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

Use `--quantity Ftotal3_y` repeatedly to select columns and `--format csv` or
`--format json` for downstream analysis.  The default averaging interval is the common
time overlap.  Duplicate times from appended restart histories retain their last value.
A block error is meaningful only after the flow is statistically stationary and each
block is longer than the force autocorrelation time; the script caps the number of blocks
at the number of sampled time intervals and reports `nan` when only one block is possible.

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
