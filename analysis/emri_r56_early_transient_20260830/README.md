# EMRI r56 production-grid early-transient analysis

This directory analyzes the completed `r56_p000_t5000` production-grid qualification.
It contains 11 de-duplicated history samples through `t=55.3767 m2`.  The run was built
to qualify AMR, force diagnostics, restart continuity, CT, memory, and I/O.  It was not
long enough to produce a stationary Bondi-Hoyle-Lyttleton solution.

**Force-definition note.**  The histories analyzed here predate the cellwise
`delta T = T-T_upstream` correction.  Their raw `Frel` values contain the analytic
upstream stress and must not be reused as physical drag.  The replacement zero test and
fixed-tree control are documented in
`analysis/emri_force_background_fix_20260830/README.md`; a production-grid rerun is
required for corrected force time series.

## Main physical conclusion

The run is numerically healthy but physically still in startup.  Its capture-crossing
time is `4220.47 m2`, so the saved interval covers only `0.01312` crossing.  Using the
planner's causal-speed proxy, the disturbance can have propagated only `4.32 m2` by the
last sample.  This barely exceeds the `3 m2` momentum-flux surface and is much smaller
than the three volume-force radii, `164.5`, `329.0`, and `658.1 m2`.

The upstream state is weakly magnetized and sub-fast rather than a classic supersonic
BHL wind:

- `beta=3388`, `sigma=1.70e-6`;
- three-speed `v=0.04365`;
- sound-speed squared `c_s^2=0.004172`;
- fast-Mach proxy `M_fast=0.676`.

Pressure therefore matters more than magnetic stress in the upstream capture scale, and
a strong, narrow BHL shock cone is not expected for this case.  Its eventual morphology
should be closer to a pressure-supported Bondi-Hoyle transition with an oblique drift.

This parameter point is not directly comparable to the published magnetized BHL runs.
[Kaaz et al. (2023)](https://arxiv.org/abs/2201.11753) used a Mach-2.45 wind and mostly
`beta=1--200`, evolved as long as 50 accretion times, and noted that the flow commonly
requires roughly 10--20 accretion times to reach steady state.  [Kim & Most
(2025)](https://arxiv.org/abs/2409.12359) used Mach 2 and `beta=10` in their fiducial
model; their tangential drag settled only around 25 accretion times.  Our `M_fast=0.676`,
`beta=3388` case is both more pressure dominated and much less magnetized.  The planned
two crossings are therefore an initial physics target with a stationarity gate, not a
literature-supported guarantee of relaxation.

## Accretion rate

The recorded proper-time accretion diagnostic rises from `44.13` to `488.93`, a raw
factor of `11.1`.  This factor is not a physical growth rate because the AMR hierarchy is
also being constructed: the calculation begins with 287 MeshBlocks and one physical
refinement level and reaches 2779 MeshBlocks and nine physical refinement levels (ten
levels including the root) only near the end.

For scale only, the planner's classical proxy is

```text
Mdot_BHL,proxy = 4 pi rho m2^2 / (v_eff^2)^(3/2) = 1.9814e4.
```

The final value is only `2.47%` of this proxy.  Even between the last two saved samples,
after the full topology has been reached, `Mdot` changes by `+8.10%`.  The accretion rate
therefore has not saturated.  The classical proxy is a normalization, not a precision
prediction for a relativistic, magnetized, sub-fast flow.

## Force history

The largest-shell raw total-force estimator at the final time is

```text
Ftotal(2 r_a) = (125.77, 95.50, -0.215)
```

in source-tetrad code units.  The small final vertical fraction, `1.36e-3` of the vector
magnitude, is consistent with recovery toward approximate equatorial symmetry.  It is
not yet a stationary inclination-force measurement.

The outer increment remains large:

```text
|Frel(2 r_a)-Frel(r_a)| / |Ftotal(2 r_a)| = 0.576,
```

well above the provisional `0.15` closure gate.  Conversely, that outer increment
changes by only `1.82%` relative to its cycle-0 value.  This is exactly what the causal
scale predicts: the outer shells have not yet responded, so their raw values mainly
measure the initialized background and discretized volume integral rather than a formed
wake.

Subtracting the cycle-0 value makes the three shell estimates look similar at the final
time,

```text
Delta Ftotal(0.5 r_a) = (-24.96, -5.47, 3.62)
Delta Ftotal(1.0 r_a) = (-23.93, -9.19, 2.85)
Delta Ftotal(2.0 r_a) = (-24.99, -7.93, 2.67),
```

which is consistent with a common near-zone response.  This subtraction is diagnostic
only: it also includes the change in AMR quadrature between the 287-block and 2779-block
meshes, so it cannot be reported as a physical drag force.

There is a second baseline systematic that should be resolved before interpreting a long
production run.  `Fnewt` subtracts `rho0` in the current implementation, while `Frel`
contracts the full stress-energy tensor with the secondary-metric derivative.  At cycle
zero, `Fnewt` is close to zero but `Frel(2 r_a)` is already about
`(152.1, 103.2, -2.87)`.  The raw relativistic estimator therefore contains a large
unperturbed-flow/background term.  This is not evidence for a large GR drag correction.
For a wake-induced-force measurement, retain the raw estimator but add either:

1. a relativistic perturbation estimator using the analytic upstream Taylor state,
   `delta T = T - T_upstream(x)`, with the same metric, tetrad, and volume quadrature; or
2. a companion fully refined background calculation on exactly the same AMR topology,
   followed by pointwise/history subtraction.

The new estimator should pass uniform-wind zero-baseline, weak-field `Frel-Fnewt`, outer-
radius, and resolution tests before it is used for an EMRI dephasing estimate.

## Numerical result

The numerical qualification itself is strong:

- restarted and continuous cycle-11 endpoints agree exactly in MHD cells, CT faces, ADM
  fields, and the de-duplicated 40-column history;
- `max|divB|=1.02e-13`;
- the restart-cache proposed change is only `2.44e-4` of its tolerance;
- the geodesic residual stays at `4.50e-12`;
- full-topology root cycles settle to roughly `109--120 s` on one A100;
- peak A100 memory is `47675 MiB`.

Thus the solver, metric cache, CT, restart, force-history continuity, memory envelope, and
runtime model are ready.  What is not ready is a stationary physical inference.

## Required time coverage

At the current causal-speed proxy, the earliest possible response times are approximately

```text
momentum surface (3 m2):       38.5 m2
0.5 r_a shell (164.5 m2):   2110  m2 = 0.5 crossing
r_a shell (329.0 m2):       4220  m2 = 1.0 crossing
2 r_a shell (658.1 m2):     8441  m2 = 2.0 crossings
```

Consequently, near-hole accretion can begin responding in the present run, but an outer-
wake force-closure statement requires essentially the full planned two crossings and
may require a substantial stationarity extension.  Science averaging should start only
after the wake has reached the relevant shell and the adjacent-window drift and
autocorrelation gates pass.  Given the 10--25-crossing relaxation behavior in the
supersonic literature, the runtime planner should permit multi-crossing continuation
rather than budgeting only one quarter-crossing extension.

## Files

- `early_transient_physics.png/pdf`: accretion, raw force, baseline difference, and shell
  closure.
- `force_decomposition.png/pdf`: momentum, Newtonian, relativistic, and total-force terms
  projected on the upstream velocity.
- `scale_and_numerical_health.png/pdf`: causal/spatial scales and exact-restart health.
- `runtime_amr_ramp.png/pdf`: AMR construction cost and full-topology cycle time.
- `derived_history.csv`: derived force and closure time series.
- `analysis_summary.json`: machine-readable numerical values.
- `plot_early_transient.py`: reproducible plotting and analysis script.

No primitive or magnetic field dump was retained locally from the cloud qualification,
so this analysis cannot show density, Mach, magnetization, or field-line slices.  Those
plots require a subsequent field-output checkpoint.
