# V5 Metric Derivative Temporary Elimination

## 1. Motivation

`compute_rhs_bssn` previously called `fdderivs` six times for conformal metric fields. Each call wrote six full 3-D second-derivative arrays, and the caller immediately contracted those arrays into one Ricci component. The next call overwrote all six arrays. The change removes this write/read round trip without combining different metric fields into a mega-kernel.

## 2. Original Data Flow

The original calls in `AMSS_NCKU_source/bssn_rhs.f90` were:

| Input field | Parity | Temporary outputs | Immediate consumer | Final value |
|---|---|---|---|---|
| `dxx` | SSS | `fxx,fxy,fxz,fyy,fyz,fzz` | old lines 379--380 | `Rxx` |
| `dyy` | SSS | same six arrays | old lines 383--384 | `Ryy` |
| `dzz` | SSS | same six arrays | old lines 387--388 | `Rzz` |
| `gxy` | AAS | same six arrays | old lines 391--392 | `Rxy` |
| `gxz` | ASA | same six arrays | old lines 395--396 | `Rxz` |
| `gyz` | SAA | same six arrays | old lines 399--400 | `Ryz` |

For each call, `fdderivs` writes all six arrays. Their first and last reads are both in the immediately following contraction. None is reused by another RHS term; the next metric call overwrites all six. The later chi and lapse `fdderivs` calls begin new lifetimes and remain unchanged.

The exact retained discrete expression is:

`contract = gupxx*Dxx(field) + gupyy*Dyy(field) + gupzz*Dzz(field) + (gupxy*Dxy(field) + gupxz*Dxz(field) + gupyz*Dyz(field))*TWO`

No continuum identity or alternative Ricci formula was substituted.

## 3. Temporary Memory Traffic

Old path per metric field and grid point: six derivative-array stores, followed by six derivative-array loads for the contraction, plus one Ricci store. New path: stencil values remain scalar/vector temporaries and only the contracted Ricci array is stored. Across six fields this removes 36 full-array derivative outputs and their immediate reads per RHS invocation. The extended symmetry buffer `fh` remains because it carries the existing boundary/parity semantics.

## 4. Prototype Component

The prototype was `dxx -> Rxx`, chosen because its SSS parity and immediate single contraction make the lifetime unambiguous. A reusable narrow helper, `fdderivs_metric_contract`, was added to `AMSS_NCKU_source/diff_new.f90`; it processes one field per call. The call in `bssn_rhs.f90` supplies the six inverse-metric coefficient arrays and receives one contracted array.

The helper retains the current fourth-order centered pure and mixed stencils. Its regular interior is `i,j,k = 3 .. ex(axis)-2`. A second loop covers `1 .. ex(axis)-1`, skips that regular interior, and uses the exact old fourth-order/second-order fallback conditions. `contract` is initialized to zero over the full array, preserving the old zero value at uncomputed points including the final index planes. `symmetry_bd(2,...)`, parity arguments, coordinate spacing, and the symmetry-dependent `imin/jmin/kmin` values are unchanged.

## 5. Numerical Equivalence

An offline Fortran harness in `/tmp/v5_metric_equiv.f90` linked both the production `fdderivs` and new helper from the same compiled `diff_new.f90`. It used a fixed seed, grids 11x10x9 and 17x15x13, symmetry values 0/1/2, and the real metric parity set SSS/SSS/SSS/AAS/ASA/SAA.

| Metric component | max abs error | max rel error | NaN | Inf | Result |
|---|---:|---:|---:|---:|---|
| Rxx | 0 | 0 | 0 | 0 | PASS |
| Ryy | 0 | 0 | 0 | 0 | PASS |
| Rzz | 0 | 0 | 0 | 0 | PASS |
| Rxy | 0 | 0 | 0 | 0 | PASS |
| Rxz | 0 | 0 | 0 | 0 | PASS |
| Ryz | 0 | 0 | 0 | 0 | PASS |

All 36 shape/symmetry/component cases were bitwise equal over the complete arrays, including fallback and zero-output regions.

## 6. Vectorization

GNU Fortran 11.4.0 was invoked with `-O3 -march=native -finline-functions` and `-fopt-info-vec-all`; LTO was omitted only from this diagnostic compile because GCC defers its report under LTO. The production build retains LTO.

- New regular-interior `i` loop: vectorized with 32-byte AVX vectors and a 16-byte remainder version.
- Boundary/fallback loop: not vectorized because of control flow, as expected.
- Original monolithic active `fdderivs` loop: not vectorized because of control flow.
- No dependency or alias diagnostic blocked the new interior loop.
- An initial implementation kept six derivative vectors live and spilled four vector temporaries. Shortening lifetimes through diagonal and mixed partial sums reduced this to one observed YMM stack temporary while retaining bitwise equivalence.

## 7. Assembly Comparison

`objdump -d -Mintel` confirms AVX/YMM arithmetic in the new regular-interior loop. The old source path necessarily stores six output arrays and reloads them in `bssn_rhs`; the new helper has one result-array store per point. The final vector loop uses YMM0--YMM15 except YMM14 and has one vector stack temporary; this is mild register pressure, not the four-spill pattern of the first prototype. No repeated six-array store/reload stream remains.

## 8. Full Metric Conversion

After the prototype passed numerical, boundary, SIMD, and assembly checks, the same one-field helper replaced the other five metric calls. This is six independent narrow calls, not one all-metric kernel. Shift, chi, lapse, AMR, MPI, physical parameters, compiler flags, and floating-point flags were not changed.

The normal release build completed successfully with `-O3 -march=native -flto -finline-functions`. No `-pg` or `-ffast-math` is present.

## MPI Benchmark Issue

The launcher hang was independent of ABE and the V5 numerical changes. On this WSL host, even `mpirun -np 1 hostname` hung before spawning a child. A system-call trace stopped in the external hwloc OpenGL component while it probed X11/NVIDIA displays: after successful local X0/NV-CONTROL requests, it blocked connecting to `127.0.0.1:6001`. Unsetting `DISPLAY` did not help because the hwloc GL component actively scans displays.

The minimal per-command workaround is `HWLOC_COMPONENTS=-gl`; no shell or Open MPI configuration was changed. With that environment variable, `mpirun -np 1 hostname`, `-np 2 hostname`, `-np 8 hostname`, and `-np 8 true` all returned normally, with the hostname commands producing exactly 1, 2, and 8 workers. `--mca btl self,vader,tcp` also worked when combined with the hwloc workaround, but produced no behavioral difference and is not required. The benchmark executables link the same system Open MPI 4.1.2 libraries used by `/usr/bin/mpirun`, and `ldd` reported no missing libraries.

An ABE `-np 1` smoke test with `HWLOC_COMPONENTS=-gl` completed AMSS initialization, read all 65,000 Ansorg records, printed `Before Evolve`, and was then stopped at the intended smoke-test boundary. Both formal 8-rank executables subsequently initialized and completed 20 steps normally.

## 20-Step Benchmark

A clean V4 worktree at commit `8433b47` was built with the same GCC/Open MPI toolchain and release flags as V5: `-O3 -march=native -flto -finline-functions`. Each run used eight ranks, identical `input.par` and `Ansorg.psid` hashes, the same host, and the same launcher workaround: `HWLOC_COMPONENTS=-gl mpirun -np 8 ./ABE`. The order was V4 then V5 for round 1 and V5 then V4 for round 2.

| Version | Round | Total Evolve (s) | Total Running (s) | real (s) | user (s) | sys (s) | Result |
|---|---:|---:|---:|---:|---:|---:|---|
| V4 | 1 | 886.979 | 890.072 | 906.96 | 7156.71 | 92.74 | successful |
| V5 | 1 | 906.133 | 909.515 | 928.16 | 7317.92 | 99.22 | successful |
| V5 | 2 | 909.730 | 914.337 | 933.46 | 7348.94 | 110.96 | successful |
| V4 | 2 | 976.959 | 981.120 | 1002.82 | 7905.13 | 104.74 | successful |

For Total Evolve, V4 mean/min/max are 931.969/886.979/976.959 seconds and range/mean variation is 9.655%. V5 mean/min/max are 907.932/906.133/909.730 seconds and variation is 0.396%. Round 1 gives V5 time reduction -2.159% and speedup 0.9789x; round 2 gives +6.881% and 1.0739x. Comparing the two-version means gives +2.579% and 1.0265x, but this apparent gain is smaller than the V4 temporal variation and the paired rounds disagree in direction. Per-timestep times likewise moved between roughly 41 and 52 seconds. The data therefore do not demonstrate stable performance improvement.

## 10. Correctness

Localized equivalence remains bitwise exact for all required grids, symmetries, parities, regular interior, fallback, and zero-output regions. End-to-end comparison used the complete numeric rows generated by each 20-step run. For each round, BH trajectory compared 20 rows/140 values and constraint compared 180 rows/1,440 values. Both datasets have `max_abs_error = 0` and `RMS = 0`; no NaN or Inf was found. The final printed timestep-20 puncture positions agree exactly: BH0 `(-0.595187, 4.40883, 1.77552e-05)` and BH1 `(0.779982, -5.45474, 2.96487e-05)`.

## 11. Decision

The implementation passes the numerical acceptance gates but does not pass the stable-performance gate: V5 was slower in round 1 and faster in round 2, while host variation dominated the mean difference. The recommendation from these contemporary samples is **REVERT**, unless a more tightly controlled benchmark environment is available for another decision run. The V5 source changes remain uncommitted.
