# Cloud output transfer

AthenaK production output must be segmented so a full remote disk cannot corrupt a long
run.  The intended protocol is:

1. Run AthenaK for a 2--4 hour segment and stop at a synchronized restart point.
2. On the compute host, publish every closed file in that segment with
   `mark_output_ready.py`.  The script checks that files are stable, computes SHA256, and
   atomically creates `SEGMENT.manifest.ready`.
3. On the workstation, run `pull_ready_outputs.py`.  It resumes partial `rsync` transfers,
   verifies size and SHA256 in a hidden incoming directory, atomically installs each
   verified file, and writes matching local and remote ACK records.
4. Delete cloud files only after their ACK exists.  Keep at least the two newest restart
   files on the compute host.  Remote deletion is intentionally not implemented by the
   puller, so a transfer failure cannot destroy the only copy.

Example on the compute host after a segment has closed:

```bash
python3 scripts/mark_output_ready.py \
  --root /data/athenak/run01 \
  --manifest-dir /data/athenak/run01/manifests \
  --segment segment-0001 \
  bin/effective_bbh.mhd_w.00001.bin rst/effective_bbh.00001.rst
```

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
filesystem.  The current NAS must be remounted read-write before use.

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

Always reserve the larger of 50 GB or twice the largest restart.  At a measured 10 MB/s,
the theoretical transfer ceiling is 36 GB/hour (0.864 TB/day); scheduling to 25 GB/hour
leaves useful retry margin.  A 1 TB backlog needs at least about 28 hours to drain.
