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
| L | L9 | `M/32` | 18,088 | 4,544 | 8 | 4 x A100-80G | 192 GiB | 11,200 GiB |
| M | L10 | `M/64` | 28,255 | 3,584 | 13 | candidate 8 x A100-80G | 320 GiB | 17,500 GiB |
| H | L11 | `M/128` | 55,380 | 2,048 | 25 | topology-only 32 x A100-40G | 640 GiB | 34,300 GiB |

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
reached 22,604 MeshBlocks at cycle 2.  Both exceeded the old static gates, so schema 5
archives the evidence hashes, the `10M` global-output policy, and uses 18,088/28,255
blocks after the 25% margin.  The M
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

## Non-publication 8xV100 long-run qualification

`effective_bbh_4pn_v100_qualification.athinput` is a deliberately reduced, standalone
software-qualification input for periods when the accepted A100/A800 L/M/H layouts are
unavailable.  It is not generated from, and does not modify, the frozen L/M/H matrix.
It must not be used as another spatial convergence tier or as evidence for resolved disk
or magnetic turbulence.

The input retains the full `[-1024M,1024M]^3` domain, the frozen v3 4PN-to-remnant table,
reference FM torus and magnetic loop, dynamic ADM spacetime, ten physical AMR levels,
moving-hole tracking, and strict level 2:1 subcycling.  It reduces the root grid from
`128^3` to `64^3`, giving `dx0=32M` and `dx9=M/16`; the persistent disk floor is `64M`
through L4 and the moving shells are `12,12,12,12,12,8,4,2,1M`.  The `4.8M` root-step
cap is intentionally retained even though the coarser root grid has a `9.6M` CFL limit,
so root synchronization, output scheduling, and AMR lookahead use the production cadence.

The cold-start hierarchy deliberately pre-applies the two reachable tracker windows
`(0,2.4,4.8)M` and `(4.8,7.2,9.6)M`, including the tracker's `1.2M` light-speed padding.
The old finest-only seed began with 1,184 MeshBlocks and then created 2,744 blocks at the
first regrid; it did not provide the complete moving-shell lookahead before the first
`4.8M` step.  The replacement encodes only the directly selected L7 and L8 parents as
disjoint cuboids and lets `MeshBlockTree` supply the strict 2:1 propagation.  A single
bounding box per shell is forbidden: that alternative was measured at 4,992 blocks and
828,544 weighted updates.

With the pinned trajectory at its reference offset and the current hierarchy builder, a
double-precision `athena -m` now reports exactly 4,320 initial MeshBlocks, with physical-
level counts `4,420,420,420,420,420,396,588,592,640`.  Their strict-subcycling cost is
605,884 MeshBlock updates per root step: 34.4% below the L initial topology and 72.6%
below advancing every block at the finest cadence.  A GID-sorted replay of the production
load balancer at eight ranks assigns `351--724` blocks and weighted cost `75,683--75,792`
per rank, leaving at least 300 slots below the configured 1,024-block cap.  Confirm that
replay with an actual eight-rank `-m` before evolution.  These are initial-topology facts,
not a trajectory-wide capacity certificate: abort rather than increasing the cap if a
later partition exceeds 1,024 blocks/rank.  Do not infer the MPI partition by reading
`mesh_structure.dat` in file order: it is grouped by physical level and must be sorted by
the `#MeshBlock GID` before any independent replay.

The deployed table path is
`/data/athenak-l4/assets/bbh-dense-single-r7-d12-v3.dat`.  Verify its SHA-256 against the
input header before launch.  A local parse may override only the path.  A serial `-m`
check additionally needs a topology-only capacity override because its one process must
temporarily hold all 4,320 initial blocks; this override must never be used for evolution:

```bash
build_bbh_cpu/src/athena -n \
  -i inputs/dyngr/effective_bbh_4pn_v100_qualification.athinput \
  problem/trajectory_file=/absolute/path/to/bbh-dense-single-r7-d12-v3.dat

mkdir -p /tmp/athenak-v100q-mesh
cd /tmp/athenak-v100q-mesh
/absolute/path/to/athenak/build_bbh_cpu/src/athena -m \
  -i /absolute/path/to/athenak/inputs/dyngr/effective_bbh_4pn_v100_qualification.athinput \
  problem/trajectory_file=/absolute/path/to/bbh-dense-single-r7-d12-v3.dat \
  mesh_refinement/max_nmb_per_rank=4320
```

The production eight-rank command must not contain that serial-only override: it uses the
input's fail-closed `max_nmb_per_rank=1024` on every GPU.

The authoritative partition check uses eight MPI processes and no capacity override:

```bash
mpirun --bind-to none -np 8 /absolute/path/to/athenak/build_bbh_mpi/src/athena -m \
  -i /absolute/path/to/athenak/inputs/dyngr/effective_bbh_4pn_v100_qualification.athinput \
  problem/trajectory_file=/absolute/path/to/bbh-dense-single-r7-d12-v3.dat
```

Before GPU evolution, require a double-precision Volta70 build, ECC enabled, one MPI rank
per V100, a CUDA-aware MPI device-pointer `Sendrecv`, the exact 4,320-block topology, and
an observed partition no larger than 1,024 blocks/rank.  The first runtime gates are
`0--9.6M`, split-versus-continuous restart equivalence through `19.2M`, and then the
`96--115.2M` checkpoint-boundary test of moving-AMR lookahead.  Only after those pass
should restart-safe segments advance toward `604.8M`; a completed orbit is established
from an unwrapped BBH phase of at least `2*pi`, not from final time alone.

The first 8xV100-32G qualification allocated only 14.07--14.14 GiB per GPU over
3,928--4,068 live MeshBlocks, leaving about 17.9 GiB per device.  A deliberately loose
`max_bsq=1e6` branch remained finite through `105.6M`, held
`|divB|_max <= 5.0e-14`, and showed no GPU-memory growth, but its outside-excision
`sigma_max` jumped after `86.4M` and reached `1.05e4`; C2P floor and FOFC rates rose at
the same time.  A fair one-root-step restart gate at `96M` found that `max_bsq=100`
reduced outside-excision `sigma_max` from 6,589 to 4,403, maximum C2P iterations from 17
to 14, and FOFC corrections by 14%, with no C2P failure and only a `1.1e-12` relative
change in baryon mass.  The standalone V100 input therefore uses 100 for long software
soaks.

This is not a science choice: production runs must include at least a `25/100/1000`
floor-sensitivity group and report injected-floor diagnostics before interpreting the
funnel or jet.

The first cap-100 soak restarted from the `96M` checkpoint and reached `134.4M` in eight
root steps with exit status zero.  It used 26m15.8s wall time (about `87.9M/hour`),
14.07--14.14 GiB/GPU, and no volatile corrected or uncorrected ECC event.  The final
checkpoint held 4,012 MeshBlocks, at most 673 blocks/rank, 351 spare slots on the most
populated rank, and exact level-subcycling costs.  At the final synchronized endpoint all
16,433,152 cells in both the native-GR and `divB` audits were finite,
`|divB|_max=1.46e-13`, there were no C2P failures or velocity-ceiling corrections, and
the maximum outside coordinate distance `r>1M` from either hole was `sigma=44.96` (none
above 100).  Baryon mass changed by -0.342% over the `96--134.4M` segment.  The raw
unmasked `sigma=1.43e9` maximum was at distance `0.974M` from a hole and is explicitly
inside the evolution's excision-floor region; this is the case that motivated the exact
mask output below.  The endpoint's weighted level work is 582,904 block-updates versus
2,054,144 for finest-cadence stepping of every live block, a 71.6% update reduction
(3.52x ideal work ratio).  This is a software-soak result, not a resolution-convergence
or physical-floor validation.

The current V100 qualification output policy is history and event log plus both
COM-following local slices every root endpoint, global primitive and native-GR state
every `48M`, the strict full-domain `divB` topology reference every root endpoint
(`4.8M`), and the source restart stream every `48M` (not `19.2M`).  Before the next
long segment, the strict planner may make the one-time audited `48M -> 19.2M` restart
transition with `--target-restart-dt 19.2`; it preserves the serialized `last_time`
phase, records the exact `output4/dt=19.2` runtime token, and serializes `19.2M` into the
new checkpoints.  The same plan may make the paired `48M -> 10M` full-domain 3-D
transition with `--target-global-dt 10`; this atomically emits
`output2/dt=10.0` for `mhd_w_bcc` and `output5/dt=10.0` for
`mhd_gr_diagnostics`, preserves their common serialized phase, and binds both old and
new schedules and endpoint counters.  Their optional `id` fields must equal their own
variables, all numbered output path templates must remain disjoint, and every planned
write plus the endpoint's next counter must fit AthenaK's five-digit filename range
(`<100000`).  A mismatch, alias, collision, or exhausted counter is fatal before
launch.  Later segments omit both one-time options and inherit the serialized tighter
cadences without another override.  The event log writes a row even when
all fault counters are zero and records 64-bit totals for conserved/primitive floors,
C2P failures and actual iteration maxima, conserved and magnetization adjustments,
non-excision FOFC corrections, normal C2P calls, and FOFC trial solves.  Detailed C2P
failure dumps are independently capped at eight per rank, so aggregate failure counts
remain exact without allowing a first-failure GPU log storm.  Transfer only
checksum-verified closed files and retain at least the latest three restart generations
until NAS acknowledgement.
These totals count the physical active zones of every leaf MeshBlock; ghost-zone C2P and
FOFC cache work still executes but is excluded, so rates are invariant to ghost width,
MPI decomposition, and restart-only cache hydration.  Excision-floor cells bypass the
normal C2P solver and are not part of the `c2p_calls` denominator.  Restart event-counter
format v2 stores this active-zone definition.  Legacy v1 pending totals included ghosts
and cannot be converted: an explicit one-time
`allow_legacy_ghost_event_counters=true` qualification discards them and writes v2.
Native GR diagnostic files produced after this qualification append an exact
`gr_excision_mask`; the dashboard masks those cells automatically while remaining able
to read earlier four-field files.  This prevents horizon-interior regularization values
(`sigma` reached `2.3e8` in the raw unmasked field) from being confused with the
outside-excision physical maximum.

### Raw event-ratio policy v2

The original campaign used per-row hard maxima of 1% for `fofc/fofc_tests` and 0.5% for
each C2P adjustment divided by `c2p_calls`.  Those values were conservative operational
heuristics, not bounds derived from a convergence theorem, a physical-volume fraction,
or a literature-standard publication criterion.  They are especially misleading with
strict level subcycling: repeated Runge--Kutta and fine-level opportunities weight the
smallest near-hole cells far more heavily than their coordinate volume or mass.

Policy v2 therefore uses 5% FOFC and 10% conserved/magnetization adjustment as per-row
**emergency guards**.  It emits a nonfatal yellow diagnostic when FOFC exceeds 1%, or
either adjustment exceeds 2%, for three consecutive root rows.  Exact zero C2P failure,
zero velocity-ceiling failure, finite fields, `c2p_it<25`, divB, baryon history, restart
integrity, GPU ECC, and capacity checks remain hard requirements.  The raw-rate change
does not waive those structural tests and does not turn a yellow interval into a
publication-quality physical result.

The change is motivated by the independently replayed cycle-602 peak: FOFC was 1.669%,
the two C2P adjustments were 3.239% and 3.045%, EOS failure and velocity-ceiling counts
were zero, and the maximum C2P iteration count was 12.  Spatial telemetry placed 93.1%
of FOFC events and at least 97.9% of both adjustment classes on logical L11, predominantly
at cylindrical radius 8--16M.  Full-state analysis showed that the dominant L11 region
occupied only about `3.55e-8` of the non-excised coordinate volume and contained roughly
`4.11e-5` of its coordinate-volume-weighted conserved density; the global baryon history
remained smooth across the peak.  This establishes why a raw event count cannot be read
as “1.669% of the simulated volume is bad.”  It does not establish that the interventions
are physically negligible.

Publication acceptance consequently remains incomplete until the campaign binds (1)
spatial intervention telemetry, (2) signed and absolute conservative correction budgets
for mass and energy, and (3) L/M/H resolution and floor/magnetization sensitivity.  An
`event_policy_v2_requalification_v1` report may qualify an immutable old endpoint only as
a continuation source.  It preserves the old predeclared failure in the report, claims
no publication acceptance, and must never replace or overwrite the original evidence.
Large closed-output inventories may be checked with `--audit-workers N` (`1<=N<=32`).
Workers audit independent files in parallel, while results are consumed in immutable
plan order and every binding is rechecked after all workers finish.  The report records
the worker count; concurrency changes throughput, never thresholds or accepted bytes.

### Opt-in FOFC spatial telemetry

For a short diagnostic replay, set `mhd/fofc_spatial_telemetry=true` in the input and use
exactly one active event-log output with `dcycle=1` and `write_zeros=true`.  The dense row
makes every root step explicit, including intervals with zero corrections.  Telemetry is
deliberately disabled by default: the disabled path allocates neither its per-cell reason
byte nor its fixed histogram, adds no atomics, and retains the existing numeric event-row
schema.  Enabling it adds one byte per allocated FOFC cell, a 582,120-bin (`4.44MiB`)
`uint64` histogram per rank, and one histogram MPI reduction per root cycle; it is not
intended for an entire production campaign.

Each physical active cell counted by `fofc` is added once to a joint histogram of:

- absolute logical AMR level (0--31 plus overflow);
- RK stage (1--3 plus `other`);
- C2P/floor cause, DMP preflag, scalar-only correction, or `unknown`;
- cylindrical radius from `fofc_telemetry_center1/2` (default 0,0), absolute height from
  `fofc_telemetry_center3` (default 0), and local lapse.

The fixed edges are `R_cyl={2,4,8,16,32,64}M`,
`|z|={0.5,1,2,4,8,16}M`, and
`alpha={0.2,0.4,0.6,0.8,1}`.  Radius and height distinguish disk/funnel location; lapse
acts as a trajectory-independent moving-hole proximity proxy.  Rank-local device
histograms survive AMR topology changes and are summed over MPI at every synchronized
root cycle.  Sparse records are written immediately before the unchanged numeric row as
comments beginning with `# fofc_spatial_v1`; the campaign checker and prefix recovery
therefore ignore them as event data.  Every `kind=summary` record asserts
`count=nfofc`, and the following `kind=bin` counts sum exactly to that value.
Non-finite coordinates and non-finite or negative lapse are tagged with the dedicated
`invalid_geometry` reason instead of being silently mixed into a physical overflow bin.
These are raw contributions to global `nfofc`, not per-bin failure probabilities.  A
fine level has more active cells and, under strict 2:1 subcycling, approximately twice as
many cell-stage exposures per added level during each root cycle.  Cross-level or spatial
rates therefore require the matching active cell-stage exposure as denominator; the
global `fofc_tests` value is not a joint per-bin denominator.  This diagnostic deliberately
does not add exposure atomics to the billions of tested cells per production step.

The histogram is diagnostic state, not evolution state: its pending bins are reset only
after event output and are not serialized in restart files.  Restart-persistent `nfofc`
totals that predate the current process are conservatively placed in the joint
overflow/`other`/`unknown` bin and reported as `unattributed`; no old location or cause is
invented.  The per-cell reason scratch is also not migrated by AMR, because FOFC consumes
and clears it inside each RK stage before root-boundary regridding.  These semantics make
one-step restart replays useful for localization while keeping the authoritative event
counter and all numerical updates unchanged.

For an older checkpoint that predates this parameter (for example the retained `c322`
checkpoint), a command-line `mhd/fofc_spatial_telemetry=true` override is intentionally
rejected because the key is absent.  Keep the checkpoint read-only and launch the reviewed
binary with both `-r` and a minimal read-only `-i` overlay containing only
`<mhd>/fofc_spatial_telemetry=true`; `LoadFromFile` adds that key without modifying the
source.  Bind the exact checkpoint, binary, trajectory, and overlay SHA-256 values in the
diagnostic manifest, use a separate output directory, and retain `dcycle=1` plus
`write_zeros=true`.  First stop after one root step and require `count=nfofc`.  Any pending
prefix counter restored from `c322` is explicitly reported as `unattributed`; if clean
attribution of just the new step is required, first account for that prefix in a separate
default-off zero-step child, but start the diagnostic itself from the unchanged source and
record the subtraction rather than promoting either child into the production chain.  A
default-off restart never inserts this telemetry key into an old checkpoint header.

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
native metric-aware diagnostics are written every `10M`, global `divB` every `25M`, and
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
undrained `10000M` archive projections are 8.74/13.65/26.75 TiB for L/M/H, including
the restart forced at every `100M` segment boundary.  Full campaigns must therefore
stream checksum-verified closed files instead of accumulating them on an instance disk.
Every generated input explicitly sets `time/root_dt_max=cfl_number*root_dx` (`4.8M` at
the baseline CFL).  This bounds moving-BBH AMR lookahead and makes the cap available for
fail-closed command-line overrides; AthenaK cannot override a parameter that is absent
from the input file.  The half-CFL variant therefore receives a `2.4M` root-step cap.
Use `--cfl-number 0.15` to generate the half-CFL temporal-convergence variant; the
generator permits decreasing, but never silently increasing, the audited baseline CFL.

At a segment endpoint, the driver treats a residual within eight adjacent representable
times of `tlim` as reached without changing the stored evolution time.  This avoids recursively advancing
every refined level through a physically meaningless roundoff-sized root step while keeping
prescribed dynamic metrics restart-equivalent to an uninterrupted run: a run-specific
`tlim` must not rewrite the arithmetic root endpoint.  Restart files preserve the actual
last completed timestep for
time-derived diagnostics and separately serialize `time/restart_dt_growth`, the
CFL/growth-limited value before any genuine final `tlim` clip.  Old checkpoints without
this internal parameter treat only a machine-roundoff-scale header timestep as an endpoint
artifact; the newly evaluated CFL and `root_dt_max` still bound the resumed step.

AMR topology is intentionally immutable inside one root step.  All finer levels complete
their recursive 2:1 substeps and reflux to the root endpoint before regridding, output,
or restart is legal.  The BBH spacetime itself is not frozen: ADM variables are evaluated
at each local RK stage time.  To keep the moving holes resolved between root-sync
regrids, the tracker samples each trajectory at the beginning, midpoint, and end of the
reachable next root interval and pads those shells by a subluminal light-speed bound.
With `ncycle_check=1`, `refinement_interval=1`, and the `4.8M` root cap, every root
endpoint therefore installs a path-covering hierarchy for the next at-most-`4.8M` step;
regridding independently inside fine substeps is neither required nor supported.

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
