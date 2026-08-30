# Partial-subcycling design note

## Semantics

Let physical AMR levels be `L0...Lmax`.  A useful control is a synchronization anchor
`S`, not merely the number of finest levels:

- levels `L0...LS` advance as one synchronous AMR group with timestep `dt_S`;
- levels above `S` retain strict recursive 2:1 subcycling;
- ordinary outputs, full-volume force integrals, regridding, and restart occur only at
  group synchronization points.

An optional user spelling `subcycling_finest_levels=N` maps to `S=Lmax-N`.  For five
levels `L0...L4`, `N=3` gives `S=1`: `L0,L1` step together, while `L2,L3,L4` subcycle.

For the current ten-level hierarchy `L0...L9`, however, `N=3` gives `S=6`.  With the
measured strict-subcycling root interval `5.03425m2`, this would synchronize every

```text
5.03425 / 2^6 = 0.07866 m2,
```

far more often than the requested `1m2`.  A target near `1m2` needs `S=2` or `S=3`
(`1.2586m2` or `0.6293m2`), meaning that the finest seven or six level transitions
remain subcycled.

## Work model

Over one old root interval, strict subcycling advances level `l` `2^l` times.  The
hybrid anchor advances each level `l<=S` `2^S` times and retains `2^l` advances above
`S`.  At `S=3`, the additional factors on `L0,L1,L2` are `8,4,2`; `L3...L9` retain
their existing CFL work.  This can be substantially cheaper than setting
`root_dt_max=1m2`, which over-resolves every level, including the finest level, by about
a factor five.

## Why this is not an output-only switch

The current recursive driver activates one exact level at a time.  A hybrid anchor needs
a new multi-level active range, synchronous flux correction inside the coarse group,
and the existing temporal interpolation/reflux contract only across the `S/S+1`
interface.  Writing a full-field dump or the present EMRI force history during an
ordinary fine substep would be invalid: coarse states and the outer force shells are at
different RK times.

Required qualification before production:

1. conservative linear-wave and shock tests across interfaces inside and above the
   anchor group;
2. anisotropic 3-D CT and `divB` tests;
3. bitwise continuous/restart equivalence at anchor synchronization;
4. MPI-decomposition invariance and load-balance weighting;
5. adaptive-tree creation/freeze and FOFC event continuity;
6. prescribed dynamic-ADM and active user-boundary tests.

Until that implementation is qualified, full force history is limited to the measured
root synchronization.  A separate near-hole history sampled at a fine-level sync is a
possible lower-cost feature for `mdot` and horizon-scale fluxes, but it must not be mixed
with the volume force shells extending to `O(r_a)`.
