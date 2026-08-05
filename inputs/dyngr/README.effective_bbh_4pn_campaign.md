# 4PN + multilevel-AMR convergence campaign

This directory contains a fail-closed generator and a machine-readable L/M/H matrix for
the equal-mass, initially nonspinning, `d=20M` effective-BBH torus baseline.  It is a
candidate qualification campaign, not by itself a claim that a run is publication ready.
The 4PN trajectory remains an external data product: commit its manifest and SHA-256, not
the roughly 40 MB ASCII table.

The three spatial tiers use the same `[-1024M,1024M]^3` domain, `128^3` root grid,
`16^3` MeshBlocks, RK2/CFL 0.3, and strict per-level 2:1 subcycling.  A persistent central
floor keeps `r<=80M` at least at physical L4 (`dx=M`); moving-hole shell radii are listed
explicitly at every level.  Only the inner levels change between tiers.

| tier | finest level | finest `dx` | MeshBlock gate | generated `max_nmb/rank` | A100-40G aggregate-memory lower bound / qualified launch | streaming scratch | undrained scratch for `10000M` |
|---|---:|---:|---:|---:|---:|---:|---:|
| L | L9 | `M/32` | 8,690 | 2,176 | 4 / 4 | 128 GiB | 940 GiB |
| M | L10 | `M/64` | 12,470 | 1,600 | 6 / 8 | 256 GiB | 1,350 GiB |
| H | L11 | `M/128` | 55,380 | 2,048 | 25 / 32 | 512 GiB | 6,000 GiB |

The gate is 1.25 times the larger of (a) the actual initial hierarchy reported by
AthenaK `-m`, (b) an alignment probe, and (c) the post-merger `2.5M` horizon-guard topology
model.  For the frozen v3 table at its reference offset, initial preseeding is
6,952/9,976/37,024 MeshBlocks for L/M/H; a nearby H alignment probe reaches 40,440.  These
are artifact-and-phase-qualified measurements, not universal counts for arbitrary q=1
tables.  The central L4 cube is intentionally larger than the persistent spherical floor
so the torus is never initialized on the root grid.  The post-merger model still controls
H (44,304 MeshBlocks before margin), and every newly generated input must be rechecked
with AthenaK `-m` before allocation.

The 4/6/25 values are arithmetic aggregate-memory lower bounds only; they are **not**
approved launch layouts.  This frozen matrix accepts exactly L/4, M/8, or H/32 MPI
ranks/GPUs.  Any other count fails closed until a real GID-ordered topology partition and
peak-memory run at that count are archived and the matrix is revised.

`max_nmb_per_rank` is not computed from `total/ranks` alone.  The generated value is the
larger of the rounded global campaign gate per rank and an artifact-qualified contiguous-
GID (Z-order) partition floor.  At 4/8 ranks, the frozen L/M initial artifacts measured
1,738--1,738 blocks at cost 231,008--231,008 and 1,247--1,247 blocks at cost
530,224--530,224 per rank.  The H 32-rank alignment artifact measured 985--1,893 blocks
and weighted cost 1,968,128--1,970,584 per rank.  Adding 8% to each measured block-count
peak and rounding to 64 gives hard floors 1,920/1,408/2,048 for L/M/H.  Therefore H uses
2,048 even though `55,380/32` rounded to 64 is only 1,792.  The matrix pins each
`mesh_structure.dat` SHA-256; file-order replay is invalid because that file is grouped by
physical level and must first be sorted by its `#MeshBlock GID`.

The aggregate GPU lower bounds are not performance or launch recommendations.  H is a
multi-node-size job on ordinary 8-GPU hosts; do not launch it over slow Ethernet merely
because aggregate memory is sufficient.  Benchmark the actual MPI fabric and complete M
first; unqualified rank counts are rejected by the generator.

## The horizon lower bound matters

At physical level `l`, the tracker uses

```text
r_effective(l) = max(r_explicit(l), 1.25 r_horizon,enclosing).
```

The horizon term is a physical lower bound at each level and is **not** multiplied by
`2^(Lmax-l)`.  For the normalized q=1 nonspinning inspiral it is about `1.25M` around each
component; after stitching to a roughly `0.95M` spinning remnant it is commonly about
`2.2M`.  Thus a nominal finest radius of `1M` does not imply a `1M` finest region.  The
matrix budgets a `2.5M` post-merger bound and then adds 25% topology margin.  The generator
rejects a trajectory whose certified guard exceeds `2.5M`; update and re-audit the matrix
instead of bypassing that failure.  Certification is interval-aware: it uses the same
Hermite middle velocity control as AthenaK, the convex-hull speed bound, linear mass/spin
bounds, and the analytic Kerr enclosing-radius monotonicity.  For the frozen table the
largest row value is `2.356581M`, while the conservative all-times certificate is
`2.398361M`.  Endpoint-only scans are not accepted because a Hermite interval can hide a
large Lorentz factor between otherwise harmless rows.

## Generate an input

Generate on the compute host so the embedded trajectory path is valid there.  The example
uses streaming output; omit `--streaming-drain --drain-mib-s 8` only when the declared
scratch disk can hold the full projected campaign.

```bash
python3 inputs/dyngr/generate_effective_bbh_4pn_campaign.py L \
  --trajectory /data/trajectories/q1_4pn_to_remnant.dat \
  --trajectory-provenance /data/trajectories/q1_4pn_to_remnant.dat.provenance.json \
  --source-artifact /data/trajectories/circular_orbit_PN_sep20.h5 \
  --expected-source-revision 2ce68e3d49e8758b32efc8841d239354d8d619d6 \
  --gpus 4 --gpu-memory-gib 40 --scratch-gib 128 \
  --streaming-drain --drain-mib-s 8 \
  --output /data/athenak/run-L/effective_bbh_4pn_L.athinput
```

Use `--gpus 8 --scratch-gib 256` for M and `--gpus 32 --scratch-gib 512` for H.
`--validate-only` performs every trajectory/resource/storage check without writing an
input.  `tlim` defaults to the table endpoint minus the metric finite-difference padding;
use `--tlim` for a short qualification segment.  Unless explicitly supplied,
`trajectory_time_offset` is derived as the first table time plus `metric_fd_step` (rounded
one representable value upward), so the initial centered stencil is covered; explicit
offsets are checked against both ends of the full `[offset-h, offset+tlim+h]` interval.
The default
`metric_fd_step=5e-5M` and the measured 14.33 MiB/MeshBlock memory model both correspond
to the validated double-precision A100+MPI build (`Athena_SINGLE_PRECISION=OFF`); the
matrix records `2.5e-5`, `5e-5`, and `1e-4M` as the required finite-difference sensitivity
triplet.  The restart storage estimate is also double precision (ordinary `bin` outputs
remain float32 by file-format design).
Use `--cfl-number 0.15` to generate the half-CFL temporal-convergence variant; the
generator permits decreasing, but never silently increasing, the audited baseline CFL.

The generator verifies syntax-independent campaign invariants, the complete 21-column
table, monotonic times, finite and subluminal states, time coverage, trajectory SHA-256,
the q=1/nonspinning/center-of-mass/20M baseline, the horizon bound, GPU capacity, and
scratch capacity.  `--trajectory-provenance` is not a label: it must be the stitcher's v2
JSON sidecar.  Its `output.sha256`, row/time metadata, 21-column schema,
`frozen-CBwaves`/`local-instantaneous-4PN` model fields, canonical-single-term remnant
representation, mass, spin, and kick are cross-checked against the hashed table.  The
sidecar must also contain `declaration_source=CLI`, the exact audited 40-hex AnalyticalBBH
commit, and a `source_provenance` object with path and 64-hex SHA-256.  The campaign matrix
currently pins revision `2ce68e3d49e8758b32efc8841d239354d8d619d6`.  The generator re-hashes the source file;
use `--source-artifact` when its local copy moved since stitching.  Hard-coded labels or a
declared-but-unavailable source artifact do not pass.  Archive the sidecar itself by the
SHA-256 copied into every generated input header.

This matrix is intentionally frozen to trajectory SHA-256
`42575c6b6a07f4a4fad22b6fecedbf4719aada11093b91556c46475379a9cac0`, sidecar SHA-256
`64247d4d8b018d48d3d8ba9b4c7048d86418e36ba4923f8fecc48d10b489bf1e`, and source
artifact SHA-256 `94f9f7cf848ba51f9186fa27506f526ffd3466fefb0d4a8af5337cf4e8864751`.
A new trajectory or stitch is a new campaign-matrix revision with fresh topology and
resource measurements, even when its human-readable labels are unchanged.

Publication-certificate runs require double precision.  Single precision is permitted
only as an explicitly acknowledged short systematic (`--real-precision single
--allow-single-precision-systematic`) with an explicit `--metric-fd-step`.  The generator
rejects it if `t-h`, `t`, and `t+h` are not distinct float32 values at either endpoint;
at late times the `h/2,h,2h` triplet can otherwise collapse onto the same float ULP and
produce a false finite-difference convergence result.  Always confirm the executable with
`athena -c` before launching.

## Qualification order

1. Convert and validate the trajectory with
   `scripts/bbh_trajectory_h5_to_ascii.py --validation-profile q1-nonspinning`, stitch a
   single-remnant tail when needed, and archive both manifests.
2. Generate all three inputs from the same table and matrix.  Keep their generated headers
   and SHA-256 values with the run metadata.
3. Run a short L segment through AMR creation, restart/resume, and synchronized output.
   Require finite conserved/primitive variables, bounded atmosphere resets, stable
   MeshBlock counts, restart equivalence, and constrained-transport `divB` at roundoff
   scale before extending it.
4. Complete L, then repeat the short and full gates for M.  Start H only after L/M science
   diagnostics show a usable convergence trend and a measured H wall-time/storage
   projection has been approved.
5. Repeat at least one resolved interval with half CFL, all three `metric_fd_step` values,
   and varied floor/excision parameters.  Compare time-aligned accretion rate, magnetic
   flux, shell-integrated stress/luminosity proxies, disk mass, and topology-independent
   volume integrals; report differences alongside L/M/H spatial convergence.

Segment cloud runs into 2--4 hour restart-safe units and use
`scripts/mark_output_ready.py` plus `scripts/pull_ready_outputs.py`.  The intended local
destination is `/home/zhangzelin/UGreenNAS/Projects/GRMHD_AthenaK`.  Verify its write access
before starting.  Follow `scripts/README.cloud-output.md`: warn at 65% remote-disk use,
start no new segment at 75%, stop after a synchronized checkpoint at 80%, emergency-stop
before another write at 85%, and always reserve the larger of 50 GiB or two restart files.

## Scientific limitation of this matrix

L/M/H converge the moving-horizon hierarchy while holding the central disk floor fixed at
L4.  That protects torus initialization and the large-scale field from a root-grid start,
but it does not establish MRI or small-scale magnetic-turbulence convergence throughout
the disk.  A paper making such claims still needs a geometry- or physics-aware persistent
disk criterion and its own disk-resolution sequence.  Moving-hole accretion/flux
diagnostics, horizon-following magnetic treatment, initial-transient tests, and
floor/excision sensitivity also remain explicit publication gates.
