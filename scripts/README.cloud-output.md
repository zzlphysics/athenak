# Cloud output transfer

AthenaK production output must be segmented so a full remote disk cannot corrupt a long
run.  The intended protocol is:

1. Run AthenaK for a 2--4 hour segment and stop at a synchronized restart point.
2. On the compute host, publish every closed file in that segment with
   `mark_output_ready.py`.  The script refuses files younger than 120 seconds, scans
   `/proc/*/fd` before and after hashing to prove that no process has the file open,
   verifies size/mtime stability while computing SHA256, and atomically creates
   an immutable `SEGMENT.manifest.ready`.  A segment identifier can never replace an
   existing ready manifest, so receivers cannot ACK one payload and later observe another.
3. On the workstation, run `pull_ready_outputs.py`.  It resumes partial `rsync` transfers,
   verifies size and SHA256 in a hidden incoming directory, atomically installs each
   verified file, and writes matching local and remote ACK records.
4. Delete cloud files only after their ACK exists.  Keep at least the two newest restart
   files on the compute host.  Remote deletion is intentionally not implemented by the
   puller, so a transfer failure cannot destroy the only copy.

Never publish the newest in-progress output merely because its size has stopped changing.
For mid-segment streaming, publish a binary dump only after the following numbered dump
exists and the closure checks pass.  History files are append-only and remain unpublished
until the AthenaK segment exits.  Restart files are published only after synchronized
shutdown; retain two newer verified restart generations before considering an ACK-gated
cleanup.

## Strict long-run segment gate

For a chained BBH run, terminate each production segment by root-cycle count.  Do not use
the decimal time limit as the primary endpoint: repeated binary64 additions can put the
stored time several ULPs away from the printed decimal value.  Create an immutable plan
from the closed source restart before launch:

```bash
python3 scripts/plan_athenak_segment.py \
  --repo /data/athenak/repo-clean \
  --source-restart /scratch/run/previous/state/rst/effective_bbh.00010.rst \
  --source-history /scratch/run/previous/state/effective_bbh.user.hst \
  --anchor \
  --binary /data/athenak/build/src/athena \
  --trajectory /data/athenak/trajectory/trajectory.dat \
  --root-steps 100 --root-dt 4.8 --tlim-guard-steps 2 --ranks 8 \
  --launcher /opt/openmpi/bin/mpirun \
  --state-dir /scratch/run/segment-0002/state \
  --staging-dir /scratch/run/segment-0002/staging \
  --evidence-dir /scratch/run/segment-0002/evidence \
  --wall-time-seconds 28800 \
  --planned-peak-output-gib 110 \
  --required-unnumbered effective_bbh.mhd.hst \
  --required-unnumbered effective_bbh.user.hst \
  --required-unnumbered effective_bbh.log \
  --output /scratch/run/segment-0002/evidence/segment.plan.json
```

Use `--anchor` only for the first segment imported from the legacy campaign.  Every
subsequent segment must instead supply the preceding immutable
`--parent-segment-pass .../segment.pass.ready` report.  The planner binds the clean Git
commit, executable and tool hashes, source history/restart, trajectory, directories,
MPI launcher, the exact `nvidia-smi` executable, GPU/rank count, and wall limit.  It also
binds a minimal explicit child environment and fixed directory descriptors, so inherited
`LD_*`, `OMPI_MCA_*`, `PATH`, and other ambient variables cannot alter the run.  The plan
filename is fixed to the direct evidence child `segment.plan.json`.  The planner also
proves that the sole runtime
cadence override `output3/dt=4.8` produces one full-domain, ghost-free `divB` topology
file on every root cycle; a serialized `output3/dcycle` is rejected.

`--planned-peak-output-gib` is a conservative peak for files newly created under the
segment's state directory.  It does not include the staged source-restart and trajectory
copies or the reserved recovery space; the launcher adds those separately.  The default
is 200 GiB;
pass an explicit measured upper bound for production (for example 30--40 GiB for a
10-root-step qualification segment or 110 GiB for a 50-root-step segment).  Before any
staging copy and again through fixed directory descriptors immediately before spawning
MPI, the launcher rejects every involved filesystem when its exact used-space fraction
reaches 75% (independent of `df` display rounding), or when it has less than
the larger of the computed budget and 180 GiB free.  State and staging contributions are
combined once when they share a filesystem and budgeted independently when they do not.

Use the plan's exact `time/nlim` value as the primary stop and its `time/tlim` value only
as a two-root-step guard.  Launch exactly that immutable plan; the launcher holds private
read-only descriptors for the staged restart and trajectory and fixed descriptors for
the state/evidence directories for the complete MPI lifetime, proves the rank-to-GPU
mapping, and writes immutable launch and completion records.  A failed launch proof
terminates the whole new MPI process group.  The launcher waits for this one segment but
never monitors or chains another segment:

```bash
python3 scripts/launch_athenak_segment.py \
  --plan /scratch/run/segment-0002/evidence/segment.plan.json \
  --state-dir /scratch/run/segment-0002/state \
  --launch-record /scratch/run/segment-0002/evidence/segment.launch.ready \
  --completion-record /scratch/run/segment-0002/evidence/segment.completion.ready \
  --run-log /scratch/run/segment-0002/evidence/run.log \
  --exit-status /scratch/run/segment-0002/evidence/exit.status \
  --gpu-before /scratch/run/segment-0002/evidence/gpu-before.csv \
  --gpu-after /scratch/run/segment-0002/evidence/gpu-after.csv
```

After a clean process exit and a 120-second settling interval, qualify the entire closed
segment:

```bash
python3 scripts/check_athenak_segment.py \
  --plan /scratch/run/segment-0002/evidence/segment.plan.json \
  --launch-record /scratch/run/segment-0002/evidence/segment.launch.ready \
  --completion-record /scratch/run/segment-0002/evidence/segment.completion.ready \
  --endpoint-restart /scratch/run/segment-0002/state/rst/effective_bbh.00020.rst \
  --state-dir /scratch/run/segment-0002/state \
  --run-log /scratch/run/segment-0002/evidence/run.log \
  --event-log /scratch/run/segment-0002/state/effective_bbh.log \
  --exit-status /scratch/run/segment-0002/evidence/exit.status \
  --gpu-before /scratch/run/segment-0002/evidence/gpu-before.csv \
  --gpu-after /scratch/run/segment-0002/evidence/gpu-after.csv \
  --output /scratch/run/segment-0002/evidence/segment.pass.ready
```

Exit status 75 means that at least one file is still inside the settling interval; wait
and retry manually.  Any other nonzero status is a failed gate.  Only a newly created,
immutable `segment.pass.ready` permits the endpoint restart to become the source for the
next segment.  The gate rechecks the source hash, exact cycle/time endpoint, strict-level
subcycling restart contract, capacity headroom, event counters, GPU UUID/ECC state, every
stored endpoint-restart Real, every planned binary field, `divB`, histories, and per-step
baryon-mass loss.  Publish and transfer the closed segment only after this qualification;
do not let a launcher silently continue past a missing or failed report.

Each immutable pass also contains `scientific_advisories`.  These advisories do not
weaken or replace the hard gate: they mark `yellow` only after three consecutive root
cycles above 0.001 for FOFC/tests, conservative-adjust/C2P, or magnetic-adjust/C2P;
after a baryon-mass loss above 0.0025 in one root step; after a loss above 0.02 in any
sliding ten-root-step window (48M when `root_dt=4.8`); or when `divB` exceeds its yellow
level.  The report names the concrete cycles, intervals, and observed rates.  Density,
energy, and temperature floor counts are always summarized per C2P call; chained
segments additionally record their normalized change from the bound parent pass as a
trend, without silently converting that trend into a pass/fail threshold.

Example on the compute host after a segment has closed:

```bash
python3 scripts/mark_output_ready.py \
  --root /data/athenak/run01 \
  --manifest-dir /data/athenak/run01/manifests \
  --segment segment-0001 \
  bin/effective_bbh.mhd_w.00001.bin \
  effective_bbh.mhd.hst effective_bbh.user.hst effective_bbh.log \
  rst/effective_bbh.00001.rst
```

The generic and user history files and the event log are append-only and must be
published only after the segment process exits.  Keep all three in the same immutable
segment manifest so post-processing can align physics diagnostics with numerical events.

Example on the workstation:

```bash
python3 scripts/pull_ready_outputs.py \
  --host root@GPU_HOST \
  --remote-root /data/athenak/run01 \
  --remote-manifest-dir /data/athenak/run01/manifests \
  --destination ~/UGreenNAS/Projects/GRMHD_AthenaK/run01 \
  --expected-mount-source 192.168.99.198:/volume1/Projects \
  --expected-mount-fstype nfs \
  --bwlimit-kib 9000 --poll-seconds 60
```

Before creating the destination, the puller uses `findmnt` to require the exact declared
mount source and an `rw` mount; it then probes with a real create/delete operation.  This
prevents a missing NAS mount from silently redirecting a large transfer onto the local
filesystem.  A managed sandbox may expose the same target as `ro`; perform the pull and
its write probe in the host mount namespace.  On 2026-08-10 the host view of this NAS was
verified `rw` by create, `fsync`, readback, SHA-256, and cleanup.

Zhixing Cloud instances default to 32 Mbps, which caps an SSH transfer at about
3.8 MiB/s even when the NAS is faster.  Before a large pull, query `MaxBandwidth` from
the instance list, call
[`change_bandwidth_query_price`](https://s.apifox.cn/b0fc397f-c455-4c9a-9d82-875fc48ae106/api-319535385),
and then call
[`change_bandwidth`](https://s.apifox.cn/b0fc397f-c455-4c9a-9d82-875fc48ae106/api-242525904)
on the running instance.  A 100 Mbps target is appropriate for a NAS path measured at
about 10 MB/s.  On 2026-08-08, the quoted incremental price was 0.204 CNY/hour (68 paid
Mbps above the free 32 Mbps), but always query the current price rather than hard-coding
it.  Stop the transfer instance promptly after both ACKs are verified so GPU, disk, and
bandwidth billing all end together.

Use the following remote disk policy for the segment launcher:

- below 65%: normal operation;
- 65%: warning and disable optional/high-cadence dumps;
- 75%: do not start another segment;
- 80%: finish the current synchronized checkpoint and stop;
- 85%: emergency stop before another output write.

Always reserve the larger of 50 GiB or twice the largest restart.  At a measured 10 MB/s,
the theoretical transfer ceiling is 36 GB/hour (0.864 TB/day); scheduling to 25 GB/hour
leaves useful retry margin.  A 1 TB backlog needs at least about 28 hours to drain.
