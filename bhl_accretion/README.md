# GR BHL Accretion Runs

This directory contains notes and helper files for reproducing the two `a=0`
Bondi-Hoyle-Lyttleton runs from Kaaz et al. 2023 with AthenaK.

## Files

- `2201.11753_paper.pdf` / `2201.11753_paper_mineru.md`: Kaaz et al. 2023.
- `2409.12359_paper.pdf` / `2409.12359_paper_mineru.md`: AthenaK GRMHD BHL reference.
- `gr_bhl_kaaz_b10nr53_debug.athinput`: magnetized `B10NR53`, `a=0`, `beta_inf=10`, `gamma=5/3`.
- `gr_bhl_kaaz_nbnr53_debug.athinput`: unmagnetized `NBNR53`, `a=0`, `gamma=5/3`.
- `plot_bhl.py`: quick-look plotting wrapper for history and binary dumps.
- `gr_fm_torus_mad_192_7_1024.athinput`: previous torus run used as the SMR/domain template.

The problem generator is `src/pgen/gr_bhl.cpp` and is compiled as a custom pgen with
`-D PROBLEM=gr_bhl`.

## Physical Setup

The simulations use code units `G = M = c = rho_inf = 1`.

The Hoyle-Lyttleton accretion radius is not an independent parameter in the current
problem generator. It is derived from the wind speed:

```text
R_a = 2GM / v_inf^2
```

In code units:

```text
r_acc = 2 / (v1_inf^2 + v2_inf^2 + v3_inf^2)
```

For the Kaaz+2023 setup:

```text
v_inf = 0.1
r_acc = 200
mach_inf = 2.45
cs_inf = v_inf / mach_inf ~= 0.0408
tau_a = r_acc / v_inf = 2000
```

If `r_acc` is explicitly present in an input file, `gr_bhl.cpp` checks that it matches
`2/v_inf^2` and aborts if it does not. Prefer setting `v1_inf`, `v2_inf`, `v3_inf`,
and `mach_inf`.

## Build

Use the configured Spack environment:

```bash
spack env activate grmhd
cmake -S . -B build-gr_bhl -D PROBLEM=gr_bhl -D CMAKE_BUILD_TYPE=Release
cmake --build build-gr_bhl -j 8
```

The executable is:

```bash
build-gr_bhl/src/athena
```

## Smoke Tests

These commands reduce the mesh and stop after one cycle. They only validate parsing,
initialization, boundary conditions, and the evolution task path.

```bash
spack env activate grmhd

build-gr_bhl/src/athena \
  -i bhl_accretion/gr_bhl_kaaz_b10nr53_debug.athinput \
  mesh/nx1=16 mesh/nx2=16 mesh/nx3=16 \
  meshblock/nx1=8 meshblock/nx2=8 meshblock/nx3=8 \
  time/nlim=1 \
  output1/dt=1.0e99 output2/dt=1.0e99 output3/dt=1.0e99 \
  job/basename=bhl_b10nr53_smoke

build-gr_bhl/src/athena \
  -i bhl_accretion/gr_bhl_kaaz_nbnr53_debug.athinput \
  mesh/nx1=16 mesh/nx2=16 mesh/nx3=16 \
  meshblock/nx1=8 meshblock/nx2=8 meshblock/nx3=8 \
  time/nlim=1 \
  output1/dt=1.0e99 output2/dt=1.0e99 output3/dt=1.0e99 \
  job/basename=bhl_nbnr53_smoke
```

## Debug Runs

Run the magnetized Schwarzschild model:

```bash
spack env activate grmhd
build-gr_bhl/src/athena -i bhl_accretion/gr_bhl_kaaz_b10nr53_debug.athinput
```

Run the unmagnetized Schwarzschild model:

```bash
spack env activate grmhd
build-gr_bhl/src/athena -i bhl_accretion/gr_bhl_kaaz_nbnr53_debug.athinput
```

The current inputs use a 6-level SMR debug grid based on the torus template:

```text
domain: [-1024, 1024]^3
root resolution: 192^3
meshblock: 48^3
finest refined box: [-16, 16]^3
```

The magnetized debug input currently stops at `tlim = 2000 = 1 tau_a`. This is useful
for early debugging, not for comparison to the paper. The unmagnetized input is set to
`tlim = 40000 = 20 tau_a`, which is the start of the quasi-steady averaging window
used in Kaaz+2023 Appendix A. A closer paper comparison needs `t >= 20 tau_a`, and
preferably longer.

## Outputs

History output:

```text
<basename>.user.hst
```

The BHL pgen writes:

```text
mdot_<r_hist>
mdotHL_<r_hist>
phi_<r_hist>      # MHD only
```

Binary dumps:

```text
bin/<basename>.<variable>.00000.bin
```

For the provided inputs:

```text
B10NR53: variable = mhd_w_bcc
NBNR53:  variable = hydro_w
```

## Visualization

Quick-look latest dump and history:

```bash
spack env activate grmhd
bhl_accretion/plot_bhl.py --run-dir . --basename bhl_b10nr53 --latest-only
```

Plot every fifth dump:

```bash
bhl_accretion/plot_bhl.py --run-dir . --basename bhl_b10nr53 --stride 5 --r-max 300
```

For the hydro run:

```bash
bhl_accretion/plot_bhl.py \
  --run-dir . \
  --basename bhl_nbnr53 \
  --variables dens,derived:pgas,derived:vx \
  --latest-only
```

Plots are written to:

```text
bhl_plots/
```

The script creates:

- `history_bhl.png`: `mdot/mdot_HL` and, for MHD, `phi`.
- `*.xz.*.png`: x-z slices.
- `*.xy.*.png`: x-y slices.

The slice plotting is delegated to `vis/python/plot_slice.py`, so any variable or
`derived:*` quantity supported there can be passed through `--variables`.

## Current Caveats

- The upstream boundary is fixed wind on `inner_x1`; the other outer faces are outflow.
- The current implementation initializes a uniform magnetic field in Cartesian
coordinates. For Kaaz+2023 aligned-field runs, use `theta_b = 0` for `B || +z`.
- `phi` is the unsigned magnetic flux through `r_hist`, not yet normalized by
`sqrt(mdot)` into the dimensionless `varphi` convention used in MAD papers.
- Drag-force diagnostics are not implemented yet.

