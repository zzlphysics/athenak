# NAS `r=56M1`, `t=20000M1` scale audit

## Source

- Run: `20251102_mad_a09375_192-7_gamma16`
- Dump: `bin/torus.mhd_w_bcc.02000.bin`
- SHA-256: `92d31c1e319f97bb3f315f659e57d3fb71e208d24e6fa487d26c6fa21f2aef90`
- Dump time/cycle: `20000M1` / `800000`
- Global method: RK2, WENO-Z, HLLE, FOFC, `gamma=5/3`
- Orbit audit: prograde circular orbit, `a/M1=0.9375`, `r_BL=56M1`,
  `Omega M1=0.002380934652520433`, `T_orb=2638.96M1`

The source root mesh is `192^3` on `[-1024,1024]^3`.  The orbit lies in static
refinement level 4 (`|x_i|<64M1`), so

```text
dx_source = (2048/192)/2^4 = 0.6666667 M1 = 66666.7 m2  (q=1e-5).
```

## Resolved global-disk scale

A density-weighted cylindrical-annulus measurement over `52 <= R < 60M1` gives

```text
<R>_rho = 56.0865 M1
H_rms   = 36.5142 M1
H/R     = 0.6510
H/dx    = 54.8 cells
midplane density coefficient of variation = 0.490
```

The exact value of `H` depends on the vertical window and includes the thick
MAD atmosphere.  It is nevertheless sufficient to show that the global dump resolves
large-scale, non-axisymmetric disk structure; it does not establish MRI convergence,
which still requires directional quality factors and a global-resolution comparison.

## Eight-phase Taylor audit

The frozen-profile builder used fitting radii `[2,3,4]M1`.  None of the eight phases
passed all strict local Taylor gates.  Every phase failed fitting-scale sensitivity;
several also failed pressure, density, velocity, or magnetic residual gates.

| phase | max scale sensitivity | other failed gates | `r_a/m2` | nominal duration `m2` | levels |
|---:|---:|---|---:|---:|---:|
| 0 deg | 0.559 | pressure | 48.78 | 1204.6 | 6 |
| 45 deg | 1.314 | density, pressure, velocity, magnetic | 51.79 | 1317.7 | 6 |
| 90 deg | 0.979 | density, pressure, velocity, magnetic | 80.32 | 2545.0 | 7 |
| 135 deg | 0.462 | density, pressure, magnetic | 141.68 | 5962.3 | 8 |
| 180 deg | 0.306 | pressure, magnetic | 138.14 | 5740.1 | 8 |
| 225 deg | 0.719 | density, pressure, velocity, magnetic | 131.70 | 5343.4 | 8 |
| 270 deg | 0.435 | pressure | 90.73 | 3055.5 | 7 |
| 315 deg | 0.802 | pressure, magnetic | 125.77 | 4986.9 | 7 |

This failure means that first derivatives inferred from only a few `0.667M1` cells are
not robust to the fitting radius.  It does not mean that the global disk is uniform.

## Multiscale interpretation

The earlier `t=5000M1` local production box was

```text
4735.39 x 4490.20 x 2638.26 m2
= 0.04735 x 0.04490 x 0.02638 M1.
```

Thus its longest side is only `0.071` of one source cell.  Its two-crossing duration,
`8440.93m2`, is `0.08441M1`, only `3.20e-5` of an orbit.  The phase-zero late-state
cost model is smaller still: approximately
`428 x 620 x 430m2 = 0.00428 x 0.00620 x 0.00430M1`, with a nominal duration
`1204.6m2 = 0.01205M1`.

Consequently a frozen local run cannot contain source-resolved time-dependent disk
eddies.  At genuine EMRI mass ratio it is a relativistic magnetized BHL response problem
conditioned on a global-disk state.  The scientifically defensible disk information is
the distribution across independently selected phase/time snapshots and robustly
smoothed large-scale gradients, not cell-by-cell temporal replay inside one local run.

The next data step is therefore an ensemble audit over several late times and phases,
with a uniform-state baseline and separately gated smoothed gradients.  A claim about
live turbulent forcing requires either a much more highly refined global patch, a
matched worldtube simulation, or controlled divergence-free synthetic fluctuations.

## Independent `t=30000M1` check

Repeating the same eight-phase audit on `torus.mhd_w_bcc.03000.bin` found one strict
Taylor pass, at phase 135 degrees.  Its maximum fitting-scale sensitivity is `0.1393`
against the `0.15` gate; its density, pressure, velocity, and magnetic residuals all
pass.  The dump SHA-256 is
`359f978b34f18c9f8a50e36e6d37328447082672b9b8f35c0d16a73c30dce9d4`.
This state has `r_a=189.4m2`, an eight-level direct-AMR estimate, and a nominal
five-crossing duration `9217m2`.  It is not marked execution-ready because the current
cost ceiling rejects the direct calculation while no clean two-sided matching overlap
exists.  The other seven phases again fail fitting-scale sensitivity.

The second epoch therefore strengthens the interpretation: robust first gradients are
occasionally available, not generic.  A first PPM4 cloud rerun should use a separately
defined uniform/intercept baseline selected from the late ensemble; a gradient contrast
should be added only for a state that passes the derivative gates and the direct-cost
gate, rather than weakening either test after seeing the result.
