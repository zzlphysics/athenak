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

The default dashboard contains density, gas pressure, temperature, the magnitude of
the stored densitized `bcc` components, an inverse-beta proxy, and AMR level.  Available
panels can be listed with `--help`; for example:

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
rechecks both size and SHA256 before reading any `.bin`, `.hst`, or event-log file, merges
overlapping restart histories in command-line order, keeps generic `.mhd.hst` and BBH
`.user.hst` columns in separate products, inventories binary times/cycles/MeshBlock
counts by physical AMR level, measures the strict-subcycling MeshBlock-update reduction,
measures output cadence, projects storage to the requested target time, and renders every
verified frame plus the final frame:

```bash
python3 vis/python/analyze_bbh_grmhd_campaign.py \
  /nas/campaign/L-segment-000 /nas/campaign/L-segment-001 \
  --output-dir output/L-full --target-time 3500 \
  --segment-span 100 --root-step-seconds 522 --drain-mib-s 8 \
  --trajectory /path/to/q1_4pn_to_remnant.dat
```

`--segment-span` includes the synchronized checkpoint forced at every cloud-segment end,
even when the configured restart cadence is longer.  With a measured root-step wall time,
the report converts primitive/divB/checkpoint production into average MiB/s and compares
it with the conservative sustained NAS drain rate.  It separately reports total NAS
archive growth and the remote working set obtained by retaining only the two newest
restart generations.  Both a machine-readable `campaign-analysis.json` and a concise
human-readable `campaign-analysis.md` are written.

The analyzer distinguishes low-cadence global `mhd_w_bcc`/
`mhd_gr_diagnostics` files from high-cadence `bbh_local_w`/`bbh_local_gr` files even
though they contain the same variable groups.  It reads already-sliced variable-size
MeshBlock records directly, derives root-cycle cadence from `dcycle` when only one frame
is present, and plots a local stream using its stored moving-window extent rather than a
fixed origin-centered crop.

The generic history remains `merged-history.{csv,hst,png}`.  Native BBH diagnostics are
written as `merged-user-history.{csv,hst,png}` and include the two 3D positions,
separation, orbital angular frequency, trajectory mass terms, outside-excision baryon
mass, proper-volume gas/magnetic integrals, mass-weighted Lorentz factor and
magnetization, coordinate `J_z`, inner-region mass, and outside-excision density and
magnetization maxima.  The analyzer also merges the root-cycle Athena event log and
reports total floor/ceiling/C2P/FOFC events plus the maximum C2P iteration count.

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
  L=output/L-full/merged-user-history.csv \
  M=output/M-full/merged-user-history.csv \
  H=output/H-full/merged-user-history.csv \
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

Fresh DynGRMHD builds also provide a compact native diagnostic stream.  It evaluates
the contractions against the synchronized ADM spatial metric inside AthenaK, including
the required `bcc/sqrt(det(gamma))` conversion.  New files append
`gr_excision_mask` (one inside the evolution's excision-floor region, zero outside),
and the dashboard automatically masks those cells from the four physical GR panels.
Legacy four-field files remain readable; use the optional `excision_mask` panel to audit
new output explicitly:

```bash
python3 vis/python/plot_bbh_grmhd.py \
  state/bin/run.mhd_gr_diagnostics.00005.bin \
  --output-dir output/gr-diagnostics \
  --panels gr_bsq,gr_lorentz,gr_sigma,gr_beta_inv,excision_mask,level \
  --plane z --extent 80
```

## Fixed-color movies and equatorial interpolation

Use `make_bbh_grmhd_movie.py` for a reproducible time sequence.  It selects only
ACK-bound campaign files, samples the complete selected sequence before rendering,
chooses one robust color range per panel, and applies those limits to every frame.
The resulting movie therefore does not flicker because an individual frame changed
its colorbar range.  For example:

```bash
python3 vis/python/make_bbh_grmhd_movie.py /nas/campaign \
  --path-prefix segment-000/ --path-prefix segment-001/ \
  --output-dir output/movie --stream mhd_gr_diagnostics \
  --panels gr_bsq,gr_sigma,gr_beta_inv \
  --plane z --location 0 --extent 40 --interpolate-plane \
  --trajectory /path/to/q1_4pn_to_remnant.dat --fps 12 --workers 4
```

`--interpolate-plane` is intended for an exact cell face such as `z=0` in an
even-cell domain.  It reads the cell-center planes immediately below and above the
requested face, requires identical projected AMR topology, and averages corresponding
fields with equal weights.  This avoids choosing one neighboring cell at a refinement
boundary.  It is deliberately rejected for `bbh_local_*` files that were already
reduced to one stored cell along the slice axis: the discarded second plane cannot be
reconstructed after the simulation.  Use full-3D `mhd_w_bcc` or
`mhd_gr_diagnostics` output when two-sided interpolation is required.

For moving-window local output, set `--extent` to the configured
`region_half_width`.  MeshBlocks intersecting the requested window may extend beyond
it; plotting their complete extents can otherwise show asymmetric blank corners that
are outside the stored region.  Cropping to the configured half-width does not discard
data from inside the requested window.

## Dynamic-spacetime limitation

The `mhd_w_bcc` output stores primitive fluid variables and densitized magnetic
components, `bcc=sqrt(det(gamma_ij)) B^i`, but not the dynamical ADM metric.  Therefore
density, `press`, `temperature`, `vel*`, and `bcc*` are exact stored values, while
`bmag`, `beta_inv_proxy`, `sigma_proxy`, and `velmag_proxy` are explicitly labelled
stored-component proxies.  They are useful for morphology and debugging but must not
be quoted as covariant GR diagnostics.

Use `mhd_gr_diagnostics` for publication measurements of `b^2`, Lorentz factor,
`b^2/rho`, and `b^2/(2p)`.  Bernoulli parameters, MRI quality factors, and moving-horizon
fluxes still require additional native diagnostics; they must not be reconstructed from
the proxy panels.

The generic history plot is likewise a global coordinate-volume diagnostic, not a
closed-system conservation proof.  The BBH user history improves the definitions and
excludes atmosphere values inside the excision-floor mask, but its integration domain
still moves with that mask.  The prescribed time-dependent metric exchanges coordinate
energy and momentum with the gas; atmosphere recovery, accretion, and boundary flux can
change domain mass.  Quote moving-surface fluxes and controlled systematics rather than
interpreting `tot-E` or `baryon_m` drift alone as a numerical error norm.
