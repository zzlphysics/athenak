# Effective binary-black-hole GRMHD

Build AthenaK with `-DPROBLEM=dynbbh` and run
`inputs/dyngr/effective_bbh.athinput`.  The problem generator supplies the boosted,
superposed Kerr–Schild metric of [Combi & Ressler](https://arxiv.org/abs/2403.13308)
through AthenaK's prescribed dynamic-ADM interface.  It initializes a uniform gas and an
optional uniform magnetic field; it is a metric/GRMHD integration setup, not an
equilibrium circumbinary torus model.

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
