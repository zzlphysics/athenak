# Cloud output transfer

AthenaK production output must be segmented so a full remote disk cannot corrupt a long
run.  The intended protocol is:

1. Run AthenaK for a 2--4 hour segment and stop at a synchronized restart point.
2. On the compute host, publish every closed file in that segment with
   `mark_output_ready.py`.  The script refuses files younger than 120 seconds, scans
   `/proc/*/fd` and writable mappings before and after hashing, fails closed on a process
   visibility gap, and hashes through one `O_NOFOLLOW` descriptor while binding the
   file's device, inode, link count, size, mtime, and ctime.  It atomically creates a
   read-only `SEGMENT.manifest.ready`; a segment identifier can never replace an existing
   ready manifest.
3. On a Zhixing workstation/transfer host, run `zhixing_pull_manifest.py`.  It requires
   an independently pinned manifest SHA256 and SSH host-key SHA256, resumes partial
   transfers, verifies every final NAS path, fsyncs the files and directories, rereads the
   remote manifest, and create-no-replace publishes matching local and remote ACKs.  The
   generic `pull_ready_outputs.py` remains useful for non-destructive legacy pulls, but
   its ACK likewise never authorizes cloud deletion.
4. Treat every local or remote transfer ACK as a receipt only.  ACK v2 carries
   `authorizes_remote_deletion: false`; old ACKs are equally non-authorizing.  No cloud
   file may be deleted from an ACK.  The current multi-client read-write NFS topology has
   no verifiable server-side read-only snapshot/write barrier, so
   `authorize_zhixing_cleanup.py` fails closed and emits no cleanup authorization.
   Remote deletion is intentionally not implemented by either transfer program.

Never publish the newest in-progress output merely because its size has stopped changing.
For mid-segment streaming, publish a binary dump only after the following numbered dump
exists and the closure checks pass.  History files are append-only and remain unpublished
until the AthenaK segment exits.  Restart files are published only after synchronized
shutdown.  A future cleanup design must retain at least the newest three independently
verified restart generations, but cleanup remains disabled until both a verifiable
storage snapshot/write barrier and a transactional per-file authorization consumer are
implemented.

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
  --target-max-nmb-per-rank 1280 \
  --launcher /data/athenak-l4/opt/strict-prrte-5.0.9/prterun \
  --mca-prefix /data/athenak-l4/opt/openmpi-5.0.9-cuda \
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
also binds the explicit, canonical `--mca-prefix` directory plus existence or absence
of the six Open MPI 5 default MCA files under
`$HOME/.{openmpi,prte,pmix}/mca-params.conf` and the explicit installation prefix's
`etc/{openmpi,prte,pmix}-mca-params.conf`.  Present files are bound by canonical path,
device/inode, owner/mode, timestamps, size, and SHA-256.  Do not create, remove, chmod,
or edit any of those paths between planning, launch qualification, and checking.  The
filename is fixed to the direct evidence child `segment.plan.json`.  The planner also
proves that the sole runtime
cadence override `output3/dt=4.8` produces one full-domain, ghost-free `divB` topology
file on every root cycle; a serialized `output3/dcycle` is rejected.

For this deployment, bind the absolute, regular copied PRRTE executable shown above
directly and bind its actual Open MPI installation/configuration tree separately with
`--mca-prefix /data/athenak-l4/opt/openmpi-5.0.9-cuda`.  Never derive the prefix from
the copied launcher's parent directory.  Do not pass an `mpirun` shell/personality
wrapper that self-execs another program: that can change the live executable or argv
after planning and invalidates the launch proof.  The strict launch environment uses
`PRTE_MCA_schizo_proxy=ompi`.

`--target-max-nmb-per-rank` is the only authorized Mesh runtime transition.  Omit it
to preserve the source restart value.  When present it must be a strict increase, may
not exceed 16384, and produces exactly one
`mesh_refinement/max_nmb_per_rank=<target>` token; equal, decreasing, duplicate, or
additional Mesh overrides fail closed.  The plan binds a conservative fixed model of
14.33 MiB per MeshBlock slot and requires the target capacity to consume no more than
80% of every rank's reported `memory.total`; it also requires the observed
`memory.used` plus the modeled requirement to fit within `memory.total`.  The launcher
obtains both total and used GPU memory from the plan-bound `nvidia-smi` and applies both
gates before spawning MPI, then repeats them immediately before process creation;
the checker rebinds the source restart to the old capacity and every normal or recovered
endpoint restart to the target capacity.  For an unchanged segment, omit the option
rather than repeating the current value.

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

Launch qualification rechecks the MCA snapshot before staging completes, immediately
before process creation, and after all ranks are live.  It reads every rank's complete
environment from `/proc` and requires the exact 39-key Open MPI/PRRTE/PMIx 5.0.9 V100
profile: fixed values are plan-bound, rank/world/PID/hostname/state/argv values are
independently derived, and all five PMIx URI aliases must name the same loopback server
and agree across ranks.  Any additional key, including an unreviewed `OMPI_*`, `PMIX_*`,
or `PRTE_*` variable, fails the proof.  The checker repeats these derivations and
rehashes the live MCA paths independently.

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

### Interrupted-segment scheduled-prefix recovery

Do not run the complete-segment checker against a partial run.  If a launched segment
ends nonzero, or if its launcher is lost while the same Linux boot can still prove every
recorded MPI/rank identity has ended, the separately plan-bound recovery tool may
qualify a shorter prefix:

```bash
install -d -m 0700 /scratch/run/recovery-0002/state \
  /scratch/run/recovery-0002/evidence
python3 scripts/recover_athenak_segment_prefix.py \
  --plan /scratch/run/segment-0002/evidence/segment.plan.json \
  --state-dir /scratch/run/segment-0002/state \
  --launch-record /scratch/run/segment-0002/evidence/segment.launch.ready \
  --completion-record /scratch/run/segment-0002/evidence/segment.completion.ready \
  --recovery-state-dir /scratch/run/recovery-0002/state \
  --recovery-evidence-dir /scratch/run/recovery-0002/evidence \
  --output /scratch/run/recovery-0002/evidence/segment.prefix.pass.ready
```

The optional `--completion-record` must be the original plan path and must contain a
nonzero exit.  If it is absent, recovery is allowed only on the launch host during the
same boot, after a bounded quiescence audit proves every recorded launcher/holder/rank
identity gone, the entire managed MPI process group gone, the exact launch GPU identity
set unchanged with zero volatile ECC errors, and no remaining compute contexts.  A boot
change is a permanent recovery failure; use an earlier already-qualified segment
instead.

An interrupted run need not contain the `tlim=... nlim=...` line printed only at the end
of Driver finalization.  For a canonically audited nonzero completion or the same-boot
closed-process case above, prefix recovery may replace that absent line only with the
immutable plan's `expected.final_cycle`/`expected.tlim` and the actual launch record's
single, byte-exact `time/nlim=...` and `time/tlim=...` argv tokens.  The recovery record
preserves that alternate binding as `run_log_prefix_audit.original_limit_evidence`.
If any Driver limit-state line is present, it remains authoritative: a contradictory
line or more than one line is fatal and cannot fall back to argv evidence.  This exception
is confined to scheduled-prefix recovery; the complete-segment checker still requires
its normal unique limit-state and cycle-limit termination evidence.

Recovery considers only restart writes whose original immutable plan labels
`scheduled`.  It chooses the highest complete candidate.  An obviously truncated later
write may be retained as suffix evidence while an older complete candidate is used, but
a later structurally complete restart that fails cadence, restart-contract, or scientific
checks forbids fallback.  Cadence is replayed through `Outputs::Execute` only—no
`Finalize` write is invented.  Every numbered prefix output must exist, unknown files or
directories fail the recovery, and planned post-prefix artifacts are SHA-256 inventoried.

History and event outputs are never edited in place.  Their exact newline-terminated
byte prefixes are copied into the new owner-only recovery state tree; discarded suffix
bytes remain in the original read-only tree and are recorded by size and SHA-256.  The
tool neither deletes nor modifies the original segment.  It first publishes an
immutable `.recovery.ready` record with `status=prepared`; that record alone is not a
qualification and is deliberately non-consumable.  The commit point is the later
immutable `.pass.ready` in the same bound evidence directory, containing an
`athenak_segment_pass` with the distinct
`qualification_mode=scheduled_prefix_recovery_v1`.  The next planner and checker
require both artifacts and explicitly revalidate this provenance.  A normal successful
segment instead carries
`qualification_mode=complete_segment_v1`; recovery does not weaken that checker.

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

Receipt-only Zhixing transfer workflow on the workstation (the manifest SHA and SSH
host-key fingerprint must come from independent trusted channels, not from this SSH
session):

```bash
python3 scripts/zhixing_pull_manifest.py \
  --ssh-state /tmp/zhixing_l4_ssh.private.json \
  --host-key-sha256 SHA256:PINNED_GATEWAY_FINGERPRINT \
  --remote-root /data/athenak/run01 \
  --remote-manifest /data/athenak/run01/manifests/segment-0001.manifest.ready \
  --expected-manifest-sha256 PINNED_MANIFEST_SHA256 \
  --segment segment-0001 \
  --destination ~/UGreenNAS/Projects/GRMHD_AthenaK/run01/segments/segment-0001 \
  --expected-mount-source 192.168.99.198:/volume1/Projects \
  --expected-mount-fstype nfs
```

The receiver requires the ready manifest to be non-writable and creates a dedicated
segment lock.  Re-running the same transaction is safe: an existing manifest or ACK is
accepted only when its immutable bytes match exactly.  If only one ACK was committed
before a network interruption, its original `verified_unix` and bytes are reused to
complete the other side.  The ACK is a non-authorizing transfer receipt even if a later
payload identity check fails and leaves one side of that receipt in place.  Neither
transfer program deletes cloud files.

Payload streams use at most 256 MiB per SSH transport, below Paramiko 2.12's 512 MiB
rekey threshold.  Every planned transport rotation and every bounded reconnect repeats
the pinned host-key check, byte-for-byte remote-manifest verification, remote-source
identity check, held destination-directory identity check, exact `findmnt` snapshot,
and NAS write probe.  Sequential partials and each parallel range are fsynced and resume
from their pinned local inode size.  Transport failures use at most four exponential-
backoff reconnects per window; during payload transfer, a changed manifest/source/mount
or a local I/O failure stops immediately before the ACK phase.  The source and already
committed local payloads are never removed by recovery.

`authorize_zhixing_cleanup.py` documents the pinned inputs and exact per-file delete-set
shape required by a future cleanup gate, while deliberately having no successful path
today.  On the deployed multi-client read-write NFS mount it exits nonzero before making
an SSH connection or creating any authorization artifact.  Local mode bits, local
`/proc` scans, advisory locks, and a client-side read-only remount cannot revoke a writer
already open on another NFS client.  There is no `--force` escape hatch.

A future implementation may enable authorization only after code verifies an immutable
server-side snapshot/generation, rehashes the complete remote source and NAS replica,
binds the exact pinned manifest plus byte-identical local and remote transfer ACKs,
enforces the newest-three-restart retention policy, and emits a short-lived immutable
authorization for explicit individual files.  A fixed transactional consumer must then
re-establish the barrier and repeat identity/SHA/retention checks immediately before
quarantining and unlinking each file.  An authorization must never be interpreted as
permission for manual, recursive, prefix, or glob deletion.  Until that backend and
consumer exist, all cloud source files remain retained.

Legacy/non-destructive example on the workstation:

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
