# BBH GRMHD post-processing

`plot_bbh_grmhd.py` creates AMR-aware multi-panel figures from AthenaK binary
outputs.  It uses the slice-selection and block-rendering approach from
`plot_slice.py`, but reads all fields for a dashboard in one file scan.  This is
substantially cheaper than scanning a gigabyte-scale dump once for every panel.

## L-tier example

```bash
python3 vis/python/plot_bbh_grmhd.py \
  /path/to/state/bin/effective_bbh_4pn_L_q8.mhd_w_bcc.00005.bin \
  --output-dir output/postprocess/L-cycle25 \
  --plane z --location 0 --extent 80 --grid \
  --trajectory /path/to/bbh-dense-single-r7-d12-v3.dat \
  --history /path/to/state/effective_bbh_4pn_L_q8.mhd.hst
```

The default dashboard contains density, gas pressure, temperature, coordinate
magnetic-field magnitude, an inverse-beta proxy, and AMR level.  Available panels
can be listed with `--help`; for example:

```bash
python3 vis/python/plot_bbh_grmhd.py "state/bin/*.bin" \
  --output-dir frames \
  --panels dens,press,bmag,beta_inv_proxy,sigma_proxy,velmag_proxy
```

Every PNG has a neighboring JSON summary, and the output directory receives a
`postprocess-manifest.json` describing all processed frames.

## Closed-file campaign workflow

For cloud campaign data, use `analyze_bbh_grmhd_campaign.py` instead of passing an
unverified transfer directory directly to the plotter.  Every segment must have a local
`.acks/*.ack` written by `scripts/pull_ready_outputs.py`.  By default the workflow
rechecks both size and SHA256 before reading any `.bin` or `.hst` file, merges overlapping
restart histories in command-line order, inventories binary times/cycles/MeshBlock
counts, measures output cadence, projects storage to the requested target time, and
renders every verified frame plus the final frame:

```bash
python3 vis/python/analyze_bbh_grmhd_campaign.py \
  /nas/campaign/L-segment-000 /nas/campaign/L-segment-001 \
  --output-dir output/L-full --target-time 3500 \
  --trajectory /path/to/q1_4pn_to_remnant.dat
```

Use `--render-every 5` for a lower-cadence preview, or `--render-every 0` to produce only
the verified inventory, merged history, cadence plot, and JSON report.  `--no-sha256`
exists only for a quick local debugging pass; its report is explicitly marked unverified
and it must not be used for a paper artifact.

The workflow never discovers arbitrary files with a glob.  A science file absent from a
local ACK is ignored, and a missing, size-changed, or hash-changed ACK record aborts the
whole analysis.  This is the downstream half of the closed-file protocol documented in
`scripts/README.cloud-output.md`.

After separate L/M/H workflows exist, build the time-aligned history comparison with:

```bash
python3 vis/python/compare_bbh_grmhd_convergence.py \
  L=output/L-full/merged-history.csv \
  M=output/M-full/merged-history.csv \
  H=output/H-full/merged-history.csv \
  --output-dir output/LMH-convergence
```

It writes aligned CSV, a comparison figure, pairwise L2 differences, and an empirical
2:1 tier-difference order.  This order only probes diagnostics affected by the moving-hole
inner hierarchy: the campaign deliberately keeps the bulk disk at physical L4, so it is
not a claim of MRI or small-scale magnetic-turbulence convergence.

Proxy panels mask cells below `1e-8` of the maximum slice density by default so
that atmosphere floors do not dominate their color scale.  Change this with
`--rho-mask-fraction`, or pass `--rho-mask-fraction 0` to disable the mask.

The separate divergence output can be rendered with a symmetric-log scale:

```bash
python3 vis/python/plot_bbh_grmhd.py state/bin/run.mhd_divb.00005.bin \
  --output-dir output/divb --panels divb,level --plane z --extent 80
```

## Dynamic-spacetime limitation

The current `mhd_w_bcc` output stores primitive fluid variables and coordinate
magnetic components, but not the dynamical ADM metric.  Therefore density,
`press`, `temperature`, `vel*`, and `bcc*` are exact output values, while `bmag`,
`beta_inv_proxy`, `sigma_proxy`, and `velmag_proxy` are explicitly labelled
coordinate-component proxies.  They are useful for morphology and debugging but
must not be quoted as covariant GR diagnostics.

Publication analysis of relativistic magnetization, Lorentz factor, Bernoulli
parameter, and horizon flux should add synchronized lapse, shift, spatial-metric,
and preferably determinant data to the science output.  The post-processor can
then be extended to perform metric contractions without reconstructing the
effective BBH metric independently.
