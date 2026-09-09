# V6-6-1 Batch17 Surface Interpolation Optimization

## 1. Background

The CPU optimization work identified `AnalysisStuff` as a major part of the
evolution cost.  The relevant path was:

```text
Evolve
  -> AnalysisStuff
    -> surf_MassPAng
      -> Patch::Interp_Points
        -> f_global_interp
          -> global_interp
            -> decide3d
            -> polin3
```

V6-4-1 had already replaced the dense surface-interpolation collective with a
specialized sparse owner-computes path:

```text
consumer's local point range
  -> determine owner
  -> Alltoallv(point index) to owner
  -> interpolate on owner
  -> Alltoallv(17 results) to consumer
  -> integrate on consumer
```

This removed dense `NN * 17` result exchange, but profiling still attributed
most of `surf_MassPAng` to interpolation on the owner ranks.

## 2. Profiling discovery

For every surface point, `surf_MassPAng` interpolates the same 17 fields in
this order:

1. `Sfx_rhs`
2. `Sfy_rhs`
3. `Sfz_rhs`
4. `chi`
5. `trK`
6. `gxx`
7. `gxy`
8. `gxz`
9. `gyy`
10. `gyz`
11. `gzz`
12. `Axx`
13. `Axy`
14. `Axz`
15. `Ayy`
16. `Ayz`
17. `Azz`

All 17 fields use the same point, owner Block, Cell-centered geometry,
6-by-6-by-6 stencil geometry, normalized interpolation coordinates, source
index mapping, reflection geometry, and six-point interpolation weights.

## 3. Root cause

The V6-4-1 owner loop nevertheless called the generic interpolation path once
per field:

```text
for each point:
  for each of 17 fields:
    global_interp
      -> recompute geometry and stencil mapping
      -> repeat reflection decisions
      -> decide3d
      -> polin3
      -> 43 polint calls
```

The necessary field traffic was therefore accompanied by 17 repetitions of
point-invariant preparation and generic routine overhead.

## 4. Failed V6-5 attempt

V6-5-1 tested precomputed indices and coefficients through
`f_global_interpind`.  It regressed significantly and was reverted.
`global_interpind` is not just a 216-element weighted dot product: it still
constructs `ya(6,6,6)`, and boundary/symmetry cases can construct and
initialize a full-Block extension `fh` through `symmetry_tbd`.  For some
surface points this produces memory traffic far beyond the local stencil.

This result motivated V6-6-1 to retain the local 6-by-6-by-6 symmetry model of
`decide3d` and fuse only the 17 fields belonging to one point.

## 5. V6-6-1 design

V6-6-1 introduces a Cell-only `global_interp_batch17` routine.  Vertex builds
do not receive a dummy or unsupported implementation.

For each point, the routine computes once:

- Cell-centered `cxI`, `cxB`, `cxT`, clamps, and normalized `cx`;
- `src_i(6,6,6)`, `src_j(6,6,6)`, and `src_k(6,6,6)`;
- `reflect_mask(6,6,6)`, with bit 0/1/2 representing x/y/z reflection;
- `wx(6)`, `wy(6)`, and `wz(6)` Lagrange weights.

It then evaluates each field independently:

```text
for each of 17 fields:
  gather 216 values using the shared source mapping
  apply reflection signs using this field's own SoA(1:3)
  reduce z -> y -> x with the shared weights
  store one result
```

Final signs are not shared between fields.  Only reflection geometry is
shared.  The routine does not call `global_interpind`, generic `polin3`, or
generic `polint`, and does not construct a full-Block `fh` array.

Each field still requires 216 double loads.  The optimization removes repeated
geometry, index arithmetic, reflection tests, weight construction, and generic
call overhead; it does not claim to remove necessary field data access.

## 6. Implementation

The implementation consists of:

- the Cell-only Fortran `global_interp_batch17` kernel and C++ ABI declaration;
- one batch call for the 17 fields in the sparse-owner owner loop;
- fallback to the original 17 `f_global_interp` calls unless the build is Cell,
  `num_var == 17`, `ORDN == 6`, and the batch routine returns success;
- unchanged sparse request/result `Alltoallv`, consumer ordering, and Weight
  behavior;
- unchanged generic `Patch::Interp_Points`, `global_interp`, and `polin3` APIs.

## 7. Profiling result

The following results are profiling builds and must not be confused with the
final clean measurements:

| Measurement (max rank) | V6-4-1 | V6-6-1 | Improvement |
|---|---:|---:|---:|
| Total Evolve Time | 924.599 s | 607.589 s | about 1.52x overall |
| `surf_masspang` | 302.606 s | 22.926 s | about 13.2x |
| `interp_points` | 301.226 s | 21.630 s | about 13.9x |
| `interp_global_interp` | 300.324 s | 20.443 s | strongly reduced |

Level-6 and Level-7 surface-path samples are empty in V6-6-1 because the
optimized `surf_MassPAng` path no longer enters the generic
`f_global_interp -> global_interp -> polin3` chain.

## 8. Clean performance setup

Level-1 through Level-7 profiling remains in the source but is controlled by:

```c
#ifndef ENABLE_ABE_PROFILING
#define ENABLE_ABE_PROFILING 0
#endif
```

The default clean build also has `VERIFY_BATCH17=0`.  Consequently it performs
no profiling `MPI_Wtime` calls, counter bookkeeping, profiling reductions or
gathers, profile printing, or dual-path correctness evaluation.  A clean
`make clean && make -j2 ABE` succeeded, and the resulting binary contained none
of the Level-1 through Level-7 report strings.

## 9. Baseline environment change

An older baseline measurement reported 846.867 seconds of evolution, or
42.34335 seconds per step.  A later clean baseline measurement in the current
device/WSL/runtime environment produced 1004.78 seconds instead.

The cause of this change has not been established.  Plausible but unproven
factors include WSL runtime state, host load, CPU frequency/power/thermal
state, Windows or hybrid-CPU scheduling, MPI placement, background tasks, and
memory/cache state.

Therefore the final comparison does **not** use the historical 846.867-second
baseline.  It uses the clean 1004.78-second baseline remeasured in the same
current environment as the final V6-6-1 run, making the final A/B comparison
more comparable without claiming that the environments are perfectly
controlled.

## 10. Final clean benchmark

Both final values below are clean, 20-step measurements:

| Version | Evolve | Total running | Evolve per step |
|---|---:|---:|---:|
| Current-environment clean baseline | 1004.78 s | 1008.06 s | 50.239 s |
| V6-6-1 clean run 2 | 586.906 s | 589.825 s | 29.3453 s |

The first V6-6-1 clean run, 681.167 seconds, is retained only as evidence of
runtime variability and is not used as the final result.

Formal comparison:

```text
Evolve speedup = 1004.78 / 586.906 = approximately 1.712x
Evolve reduction = (1004.78 - 586.906) / 1004.78 = approximately 41.59%
Total-running speedup = 1008.06 / 589.825 = approximately 1.709x
```

Thus the reported result is 50.239 to 29.345 seconds per step, a 1.712x
evolution speedup and 41.59% reduction in single-step evolution time.

## 11. Correctness verification

The Cell batch kernel was compared directly with 17 calls to the original
production `f_global_interp` path on real owner Blocks and surface points.

The current `Symmetry=1` run covered reflection masks `[0, 4]`: no reflection
and z reflection.  It did not cover x, y, or multi-axis reflection, and no such
claim is made.

Across 23 sampled points and 391 field samples:

| Statistic | Value |
|---|---:|
| maximum absolute difference | 6.93889390390722838e-18 |
| mean absolute difference | 1.33154764489875235e-19 |
| maximum relative difference | 4.43827200454890002e-15 |
| mean relative difference | 1.87151296325769851e-16 |

These differences are at double-precision rounding scale for the tested
configuration.  Further x/y/multi-axis reflection validation would require an
appropriate `Symmetry=2` run.

## 12. Final conclusion

V6-6-1 removes the dominant repeated interpolation preparation in
`surf_MassPAng` while preserving the V6-4 sparse owner-computes data flow and
per-field symmetry semantics.  In the current clean A/B environment it reduces
evolution time from 1004.78 to 586.906 seconds over 20 steps, or from 50.239 to
29.345 seconds per step.  The measured speedup is 1.712x, with a 41.59% time
reduction.  Correctness testing for the active `Symmetry=1` configuration found
only double-precision-scale differences.

## 13. Remaining gap to 2x

A 2x target relative to the current clean baseline is:

```text
50.239 / 2 = 25.1195 seconds per step
```

The current result is 29.3453 seconds per step, so it has not reached 2x.  The
remaining factor required from the current version is:

```text
29.3453 / 25.1195 = approximately 1.168x
```
