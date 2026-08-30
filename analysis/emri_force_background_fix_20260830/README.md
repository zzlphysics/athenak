# EMRI force-background and resolution audit (2026-08-30)

This audit accompanies the cellwise relativistic-force background subtraction and the
fixed-AMR-tree production hand-off.

The local Taylor-gradient zero test used the same initial GRMHD state twice.  Contracting
the full stress tensor gave `||Frel||=21.7290`; contracting
`delta T = T - T_upstream(x,t)` gave `||Frel||=4.32e-15`, with a largest component of
`3.16e-15`.  The remaining background-subtracted momentum flux is a surface
interpolation/quadrature residual, so the production campaign measures it by resetting
the analytic wind on a copy of the exact warmup restart tree and taking one negligible
probe step.

For `r56_p000_t5000`, the actual anisotropic root spacings and nine binary refinements
give finest coordinate spacings `(0.07226, 0.07830, 0.06441) m2`.  With a nonspinning
secondary, `r_+=2m2`, giving `(27.7, 25.5, 31.1)` cells per horizon radius, or
`(55.4, 51.1, 62.1)` per diameter.  The AthenaK SANE-disk paper used `a=0.9375`,
so its `0.125M`, `0.0833M`, and `0.0625M` runs had about `10.8`, `16.2`, and `21.6`
cells per Kerr horizon radius.  The current coordinate-axis count is therefore above
even that paper's highest near-horizon count, but that does not replace a convergence test:
the present setup uses PLM rather than the paper's PPM4, and magnetic flux, force, and
wake statistics can converge differently from density or mass accretion.

The cloud I/O qualification did create a final `mhd_w_bcc` file of 410,039,066 bytes and
recorded its SHA-256.  It was not copied into the small local evidence bundle, and the
provider disk was subsequently released.  The scheduled science cadence is one field
dump every quarter capture-crossing (`dt=1057.1925 m2`); the completed `55.3767 m2`
qualification was only `0.0131` crossing.  Forced terminal output, rather than the
science cadence, is why the I/O qualification nevertheless had one final primitive file.

The local upstream state is not an artificial uniform torus atmosphere: it is the
first-order source-tetrad Taylor fit at `r_BL=56M`, phase zero, to the evolved global
MAD snapshot `torus.mhd_w_bcc.00250.bin` at `t=5000.01M`.  It carries the instantaneous
density, pressure, three velocity components, magnetic field, and their fitted spatial
gradients.  It does **not** carry the full turbulent eddy spectrum or short-time
variability: the local volume and boundaries are smooth analytic Taylor fields and are
frozen for this campaign.  Thus it samples a real turbulent MAD environment only through
one local instantaneous value and gradient.  Frozen snapshots at other phases/times form
an ensemble; true time-dependent turbulence requires sufficiently cadenced global
worldtube data or a separately controlled stochastic/turbulent boundary model.
