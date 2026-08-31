# V5 Lopsided SIMD Optimization

## 1. Motivation

Profiling attributes about 8.04% of runtime to `lopsided_`. The active fourth-order advection kernel executes 21 times per BSSN RHS and previously failed to vectorize.

## 2. Original Kernel

`lopsided(ex,X,Y,Z,f,f_rhs,Sfx,Sfy,Sfz,Symmetry,SoA)` adds three directional derivatives for one grid function to `f_rhs`. `Sfx/Sfy/Sfz` are the shift components. For positive shift in direction `q`, the regular stencil is `s/(12*dq) * (-3*f[-1] - 10*f[0] + 18*f[+1] - 6*f[+2] + f[+3])`. For negative shift it is `-s/(12*dq) * (-3*f[+1] - 10*f[0] + 18*f[-1] - 6*f[-2] + f[-3])`. Exactly zero shift performs no update.

The active configuration is `ghost_width=3`, `GAUGE=0`, and cell centering. Points without a complete upwind stencil retain the original centered fourth-order fallback and one-sided opposite-stencil fallback. `symmetry_bd(3,...)` supplies parity-dependent negative-index ghost values; symmetry 1 extends z and symmetry 2 also extends x/y.

## 3. Vectorization Failure

GCC reported the original inner loop as `not vectorized: control flow in loop` and vectorized zero loops in the function. The blocker was the per-point shift-sign branch together with directional boundary/fallback branches, not dependency, aliasing, unknown stride, or a function call inside the hot loop.

## 4. Candidate Designs

Computing both complete positive and negative stencils would nearly double stencil arithmetic. `MERGE` on real values remained scalar control flow under GCC. Sign-region segmentation would require an extra scan and is poorly suited to arbitrary sign patterns. The selected design uses integer bit masks and `MERGE_BITS` to select the required neighbor values, then evaluates one stencil.

## 5. Selected SIMD Design

The common regular 3-D interior is split from boundary/fallback points. Within it, sign masks select the positive or negative neighbor ordering, `abs(shift)` supplies the original magnitude, and a nonzero mask preserves exact zero-shift behavior. x, y, and z are accumulated in their original order inside one vector loop to preserve floating-point rounding. Every non-common-interior point executes the unchanged original directional logic.

## 6. Numerical Equivalence

The complete production candidate was compared with a renamed object built from the original source. Tests covered 11x10x9 and 17x15x13 grids; symmetry 0/1/2; all seven production parity triples; and all-positive, all-negative, all-zero, alternating, random, tiny-positive, and tiny-negative shifts. All 294 cases were bitwise equal over complete arrays: `max_abs_error=0`, `max_rel_error=0`, NaN=0, Inf=0.

## 7. Compiler Vectorization

With the unchanged release-compatible flags `-O3 -march=native -finline-functions`, GCC vectorized the regular inner loop using 32-byte vectors plus a 16-byte remainder. The exact boundary loop remains scalar because it intentionally retains control flow.

## 8. Assembly Analysis

`objdump -d -Mintel` shows YMM `vcmpltpd`/`vcmpneqpd`, `vpblendvb`, vector multiply/add/FMA, and vector stores in the hot loop. There is no per-element branch in that loop. Register pressure causes two YMM stack slots for broadcast scale constants; no full stencil temporary is spilled or materialized.

## 9. Microbenchmark

An independent 96x96x48, random-shift, 40-call benchmark ran seven interleaved trials. Original times were 0.268--0.286 s (mean 0.275000 s); candidate times were 0.102--0.107 s (mean 0.104429 s). Mean kernel speedup was 2.633x, with per-trial speedups 2.576--2.750x.

## 10. 20-Step Benchmark

All runs used the same executable inputs, eight MPI ranks, and the temporary launcher setting `HWLOC_COMPONENTS=-gl`. `Total Evolve Time` is the primary metric; shell `real` is reported separately and is not mixed with it.

| Version/run | Total Evolve (s) | Total Running (s) | real (s) | user (s) | sys (s) |
|---|---:|---:|---:|---:|---:|
| V4 run 1 | 852.284 | 855.118 | 875.78 | 6913.07 | 85.91 |
| V5 run 1 | 871.981 | 875.893 | 897.23 | 7074.24 | 94.59 |
| V5 run 2 retry | 935.481 | 941.534 | 979.14 | 7694.73 | 119.07 |
| V4 run 2 | 885.629 | 889.002 | 930.96 | 7328.93 | 102.23 |

The intended order was V4 run 1, V5 run 1, V5 run 2, V4 run 2. The first attempt at V5 run 2 was externally interrupted after timestep 15 and is excluded; it was restarted from a clean `run2_retry` directory. The incomplete directory was retained rather than mixed with the valid result.

| Total Evolve statistic | V4 (s) | V5 (s) |
|---|---:|---:|
| mean | 868.9565 | 903.7310 |
| median | 868.9565 | 903.7310 |
| min | 852.284 | 871.981 |
| max | 885.629 | 935.481 |
| range | 33.345 | 63.500 |
| population CV | 1.919% | 3.513% |

The paired differences are unfavorable in both directions: V5 run 1 minus V4 run 1 is +19.697 s (+2.311%), and V5 run 2 minus V4 run 2 is +49.852 s (+5.629%). The ratio of the aggregate means is `V4/V5 = 0.9615x`, equivalently V5 takes 4.002% more Evolve time on average.

Per-timestep distributions also favor V4. V4 run 1 has mean/median 43.5185/44.1152 s and V5 run 1 44.5197/44.4973 s. Under the later, slower machine period, V5 run 2 has mean/median 48.4209/47.6215 s and the immediately following V4 run 2 46.1906/46.1492 s. The corresponding timestep population CVs are 4.435%, 0.993%, 4.255%, and 1.344%. Thus the result is not an artifact of comparing only the two-version means: both paired trend and timestep medians reject an end-to-end gain.

## 11. Correctness

For both complete V4/V5 pairs, `bssn_BH.dat` and `bssn_constraint.dat` are textually identical after excluding the file-creation timestamp. Consequently the numerical-column comparisons give BH trajectory `max_abs_error=0`, RMS=0 and constraint `max_abs_error=0`, RMS=0. A case-insensitive scan of all eight comparison files found NaN=0 and Inf=0.

Both versions finish at the same puncture positions:

- BH 0: `(-0.595187, 4.40883, 1.77552e-05)`
- BH 1: `(0.779982, -5.45474, 2.96487e-05)`

## 12. Decision

**REVERT.** The candidate is bitwise correct, genuinely vectorizes, and accelerates the isolated random-shift kernel by 2.633x. It nevertheless slows both contemporary 20-step comparisons, by 2.311% and 5.629%, with a 4.002% mean Evolve-time increase. The isolated SIMD kernel speedup therefore did not translate into an application speedup.

Possible explanations include changed cache behavior, additional simultaneous load streams, memory-bandwidth pressure, vector mask/blend overhead, and interaction with the surrounding BSSN computation. These are hypotheses only; the measurements do not establish any one of them as the root cause.

The experimental `lopsidediff.f90` implementation is preserved in Git history and then reverted, leaving the production source at its pre-experiment state. No push was performed.
