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

| tier | finest level | finest `dx` | MeshBlock gate | generated `max_nmb/rank` | A100-40G aggregate lower bound | accepted launch layout | streaming scratch | undrained scratch for `10000M` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| L | L9 | `M/32` | 18,088 | 4,544 | 8 | 4 x A100-80G | 192 GiB | 10,700 GiB |
| M | L10 | `M/64` | 28,255 | 3,584 | 13 | candidate 8 x A100-80G | 320 GiB | 16,600 GiB |
| H | L11 | `M/128` | 55,380 | 2,048 | 25 | topology-only 32 x A100-40G | 640 GiB | 32,500 GiB |

The gate is 1.25 times the larger of (a) the actual initial hierarchy reported by
AthenaK `-m`, (b) an alignment probe, (c) the post-merger `2.5M` horizon-guard topology
model, and (d) an archived root-step-capped runtime or capacity-sweep observation when
one exists.  For the frozen v3 table at its reference offset, initial preseeding is
6,952/9,976/37,024 MeshBlocks for L/M/H; a nearby H alignment probe reaches 40,440.  These
are artifact-and-phase-qualified measurements, not universal counts for arbitrary q=1
tables.  The central L4 cube is intentionally larger than the persistent spherical floor
so the torus is never initialized on the root grid.  The post-merger model still controls
H (44,304 MeshBlocks before margin), and every newly generated input must be rechecked
with AthenaK `-m` before allocation.

The L rootcap run reached 14,372 MeshBlocks at cycle 2; its trajectory-wide capacity
sweep reached 14,470 total blocks and a 3,899-block rank peak.  The M four-rank probe
reached 22,604 MeshBlocks at cycle 2.  Both exceeded the old static gates, so schema 4
archives the evidence hashes and uses 18,088/28,255 blocks after the 25% margin.  The M
observation is a capacity input, not an eight-rank qualification result.

The 8/13/25 values are arithmetic A100-40G aggregate-memory lower bounds only; they are
**not** approved launch layouts.  The per-rank allocation rejects 40G for L and M.  This
frozen matrix accepts exactly L/4 on 80G, M/8 on 80G, or H/32 on 40G MPI ranks/GPUs.  M/8
remains a candidate layout until the real eight-rank runtime and partition evidence are
archived; H remains topology-only.  Any other count fails closed until a real GID-ordered
topology partition and peak-memory run at that count are archived and the matrix is
revised.

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

For L/M the larger gate-derived allocations, 4,544 and 3,584 blocks/rank, now dominate
those initial-partition floors.  L's separately archived full-sweep rank peak of 3,899 is
also below the new allocation.  M/8 must still verify the actual weighted partition; a
global total divided by eight is not treated as evidence of rank balance.

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
  --gpus 4 --gpu-memory-gib 80 --scratch-gib 192 \
  --streaming-drain --drain-mib-s 8 \
  --output /data/athenak/run-L/effective_bbh_4pn_L.athinput
```

Use `--gpus 8 --gpu-memory-gib 80 --scratch-gib 320` for the M qualification candidate
and `--gpus 32 --gpu-memory-gib 40 --scratch-gib 640` for H topology qualification.
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
The default closed-file layout is hierarchical.  Generic and BBH-specific history are
written at every synchronized root cycle.  Global three-dimensional primitive state and
native metric-aware diagnostics are written every `50M`, global `divB` every `25M`, and
restart every `250M`.  Two additional `80M`-wide equatorial streams follow the
mass-weighted `bbh_com`: `bbh_local_w` stores primitive/MHD state and `bbh_local_gr`
stores covariant GRMHD diagnostics every root cycle.  With the baseline root-step cap the
local-frame and history spacing is at most `4.8M`, not `1M`: synchronized AthenaK output
cannot observe intermediate fine-level substeps.  The event log also checks and resets
EOS floors, velocity ceilings, C2P failures/iteration maxima, and FOFC events every root
cycle.

Local output selects complete AMR MeshBlocks intersecting the moving cube and then writes
only one cell plane through its instantaneous center.  Thus it preserves MPI-IO record
shape and AMR provenance while reducing each selected `16^3` block to `16^2` cells.  The
storage gate pessimistically assumes every gated MeshBlock intersects the local window;
the measured campaign analyzer uses actual closed-file sizes instead.
The streaming scratch gate budgets a conservative `100M` cloud segment, two retained
restart generations, forced end-of-segment outputs, and 25% headroom.  The corresponding
undrained `10000M` archive projections are 5.94/9.27/18.17 TiB for L/M/H, including
the restart forced at every `100M` segment boundary.  Full campaigns must therefore
stream checksum-verified closed files instead of accumulating them on an instance disk.
Every generated input explicitly sets `time/root_dt_max=cfl_number*root_dx` (`4.8M` at
the baseline CFL).  This bounds moving-BBH AMR lookahead and makes the cap available for
fail-closed command-line overrides; AthenaK cannot override a parameter that is absent
from the input file.  The half-CFL variant therefore receives a `2.4M` root-step cap.
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

## Two-A100 infrastructure qualification (2026-08-05)

A compact infrastructure run, deliberately smaller than tier L, exercised four physical
AMR levels, strict 2:1 level subcycling, the frozen v3 4PN table, reference FM-torus
initialization, adaptive refinement, shared restart output, and MPI restart on two
A100-SXM4-40GB GPUs.  The full fresh/restart run used commit `85d552d3`; commit
`d03ac167` then fixed and requalified a restart-segment-only load-balance diagnostic.
This does **not** qualify the L/M/H spatial convergence tiers or the H post-merger
partition.

The host OpenMPI 4.1.2 build was rejected after it reported no CUDA support and failed on
the first device-buffer send.  An isolated OpenMPI 5.0.9 configured with CUDA reported
`mpi_built_with_cuda_support=true`; an independent two-rank device-pointer `Sendrecv`
passed before AthenaK was launched.  AthenaK was linked to that installation rather than
the system MPI.

The fresh hierarchy began with 1,296 MeshBlocks, split exactly 648/648 at weighted cost
2,416/2,416.  After two root cycles it had created 1,232 blocks (2,528 total), with
1,264 blocks and weighted cost 7,636 on each rank.  Fresh and cycle-2-to-4 restart both
exited zero; the table content fingerprint was identical across segments.  The stable
post-refinement hierarchy requires 15,272 MeshBlock updates per root cycle under level
subcycling, versus 40,448 if every block were advanced at the L4 cadence: a 62.2% update
reduction (2.65x fewer updates), not a claim of equal wall-clock speedup.  Across the two
fresh topologies the measured count was 20,104 rather than 61,184 updates (67.1%
reduction).

Peak allocated memory was 27,429 MiB per GPU with `max_nmb_per_rank=2000`, consistent
with and slightly below the matrix's conservative 14.33 MiB per configured MeshBlock
slot.  A 2,528-block double-precision restart was 7,024,456,205 bytes, validating the
restart-byte model to file-header accuracy.  The qualification directory occupied 38
GiB; only 1.1 MiB of logs, inputs, build/MPI/GPU provenance, histories, and hashes was
downloaded.  Large checkpoints remain on the retained cloud disk.

The complete SHA-verified metadata is archived under
`/home/zhangzelin/UGreenNAS/Projects/GRMHD_AthenaK/4pn_amr_campaign/cloud_636963/` in
`20260805-a100x2-4level-subcycling-85d552d3-metadata-v2` and
`20260805-a100x2-restart-d03ac167-diagnostic`.  The latter records why an initial wrapper
marker was a decimal-format false negative and retains the subsequent strict verification.

## Four-A100 controlled subcycling comparison (2026-08-08)

A same-binary, same-restart, same-hardware comparison advanced the L hierarchy from
`t=48M` to `52.8M` on four A100-80G GPUs.  Strict level subcycling took `513.286s`; the
historical uniform-finest-step path took `697.485s`, a measured `1.3589x` wall-clock
speedup and `26.4%` time reduction.  At the common starting topology, the analytical
MeshBlock-update counts are 1,796,672 versus 4,448,256, a `59.6%` reduction (`2.4758x`
fewer updates).  The observed uniform run regridded after every fine step and deleted
2,772 blocks, so its realized update reduction was only `41.5%`.  Subcycling achieved
`79.5%` of the uniform path's updates/s because recursive synchronization and smaller
level-local batches reduce GPU occupancy; therefore the update-count ratio must not be
reported as the wall-clock speedup.

The 68 KiB evidence bundle is archived below `cloud_638770` as
`uniform-cycle10-to52p8`.  Its `comparison-summary.json` SHA-256 is
`09d69c0a0a48bdb1d62a019476df6c6cc204ffe3fe47f5f7eb7333676c171f82`.

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
