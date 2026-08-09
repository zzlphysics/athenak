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
