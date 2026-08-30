# AMSS-NCKU V5 Structural Optimization Audit

## 1. Current Performance

- Branch base: `optimize-v4`, commit `2f971cb`; `optimize-v4` matched `origin/optimize-v4` (ahead/behind `0/0`) before `optimize-v5` was created.
- Clean baseline: 41.1653 s/top-level step; V4: 38.7291 s/top-level step; measured V4 speedup: 1.0629x.
- V3 profile supplied with the task: `compute_rhs_bssn` 48.82%, `polint` 11.33%, `lopsided` 8.04%, `fdderivs` 6.93%, `kodis` 5.59%, `prolong3` 4.36%, `fderivs` 4.32%.
- The local `gprof.txt` is consistent with RHS dominance but is a single-process profile: 5,449 total RHS calls, of which 5,280 come from `Step`; RHS accumulated 85.1% including children. It reports 261,936 `lopsided`, 103,683 `fdderivs`, 174,624 `kodis`, and 167,235 RHS-originated `fderivs` calls.
- No 20-step run was performed for this audit.

## 2. Numerical/AMR/RK Structure

The configured calculation has nine AMR levels (0--8), five static and four moving. The log reports one patch/box on levels 0--4 and two on levels 5--8. Time refinement starts at level 3 in the runtime log (the input summary calls this grid level 4 using one-based wording).

`RecursiveStep(0)` calls one `Step` on levels 0--3 and then doubles the number of level steps at each finer level. Thus, per top-level step, the number of `Step` calls by level is:

| Level | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | Total |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `Step` calls | 1 | 1 | 1 | 1 | 2 | 4 | 8 | 16 | 32 | 66 |

Over 20 top-level steps this is 1,320 `Step` calls, exactly matching `gprof.txt`. Every `Step` performs one predictor RHS and three corrector RHS evaluations, so each owned block receives four RHS evaluations per level step. RK array updates are performed afterward by `rungekutta4_rout`, once per evolved grid function; the profile has 126,720 RK-array calls = 5,280 RHS stages x 24 evolved variables.

After the predictor, `Parallel::Sync(..., SynchList_pre, ...)` exchanges the full variable list. The three corrector stages are then evaluated; after all correctors, `Parallel::Sync(..., SynchList_cor, ...)` exchanges the corrected list. Coarse/fine restriction and prolongation occur from `RecursiveStep` after a level step, not inside the Fortran RHS.

## 3. compute_rhs_bssn Call Graph

Observed evolution path:

`Evolve` -> `RecursiveStep(level)` -> `Step(level, YN)` -> local patch -> local block -> four RK stages -> `compute_rhs_bssn` -> derivative/RHS/advection/dissipation kernels -> RK array update and physical boundary handling -> level synchronization -> recursive finer-level subcycling -> restrict/prolong.

For the active preprocessing configuration, one `compute_rhs_bssn` invocation contains:

- 1 `fderivs_shift3` call (three shift fields together)
- 18 ordinary `fderivs` calls
- 11 `fdderivs` calls
- 24 `lopsided` calls
- 24 `kodis` calls
- no direct `symmetry_bd` call; symmetry work occurs inside derivative/advection/dissipation routines

This agrees exactly with the profile ratios for a subset of calls: 261,936 / 5,449 = 48 `lopsided` calls is twice the active preprocessed source count because the local profile binary/source configuration differs from the current V4 preprocessing in this detail; 103,683 / 5,449 ~= 19.03 `fdderivs`, 174,624 / 5,449 ~= 32.05 `kodis`, and 167,235 / 5,449 ~= 30.69 `fderivs`. Therefore the profile is valuable for hotspot ranking but must not be treated as an exact call count for the currently checked-out source.

The reported approximately 43,656 RHS calls are an aggregate across eight MPI rank profiles, not 43,656 calls in one process. A representative rank has 5,449 total calls; 8 x 5,449 = 43,592, close to the aggregate, with the small difference explained by rank-dependent ownership and diagnostic/initialization calls. The deterministic evolution component is 5,280 per representative rank: 20 top-level steps x 66 level steps x 4 RK stages. Its dominant multipliers are AMR subcycling and RK stages; patch/block ownership determines calls per rank. MPI decomposition does not multiply calls within a rank, but summing all rank-local profiles does. The remaining 169 calls in the representative profile are 160 from grid construction, 9 from constraints, and not time evolution.

## 4. RHS Memory-Traffic Analysis

The RHS is a sequence of whole-box operations, not one point-local loop. It repeatedly streams the same state fields and a large collection of full-size work/output arrays:

- First derivatives materialize 3 arrays per scalar field. Metric derivatives alone materialize 18 full arrays; Gamma derivatives add 9; shift derivatives add 9.
- `fdderivs` produces six full second-derivative arrays per call. For the six metric components, the same six scratch arrays (`fxx`...`fzz`) are overwritten after each complete box pass and immediately reduced into one Ricci component. This is the clearest high-volume, single-consumer temporary traffic.
- `Rxx`...`Rzz` are first used as raised-A temporaries and later overwritten by Ricci values. `gxxx`...`gzzz` similarly change semantic role from metric derivatives to first-kind connections. This saves allocation but creates long live ranges and repeated reads/writes.
- Each of 24 `lopsided` calls rereads all three shift arrays plus one state and one RHS array. Each of 24 `kodis` calls separately streams one state and one RHS array. These are bandwidth-heavy repeated box traversals.
- Whole-array Fortran expressions in the algebraic RHS create compiler-generated loop nests; the report shows inner `i` loops vectorize, while outer `j/k` loops naturally do not. Fortran column-major layout is respected: explicit loops use `i` innermost, so unit-stride access is correct. The code is SoA at the C++ grid-function level (`fgfs[variable][point]`), not point-wise AoS.
- Scalar grid spacing coefficients are generally computed outside hot loops in derivative routines. The main repeated arithmetic is tensor algebra per point; it is arithmetic-heavy and SIMD-friendly, whereas derivative/advection/dissipation passes are primarily memory-bandwidth/cache-capacity limited.

The worst structural memory-traffic point is the six serial metric `fdderivs` calls in `bssn_rhs.f90`: each emits six N-element arrays only to consume them immediately into one Ricci array, causing 36 N stores plus later loads per RHS. This is a better target than another all-metric first-derivative fusion: compute only the contracted second-derivative result needed for each metric component, or contract in tiles while retaining the already-successful unit-stride `i` loop.

## 5. Finite-Difference Vectorization Analysis

Compiler: GNU Fortran 11.4.0. Production flags were retained (`-O3 -march=native -finline-functions`) and no floating-point-relaxing option was added. Because GCC emits no compile-stage vector report with this build's `-flto`, the diagnostic compilation omitted LTO only to expose the same front-end loop decisions; objects and reports were written under `/tmp/amss-v5-vec.xQLmTn`.

- `rungekutta4_rout`: all important innermost array loops at lines 112, 116, 118, 122, 124, and 128 vectorized with 32-byte and 16-byte versions.
- `kodis`: the hot unit-stride interior loop at line 172 vectorized. Outer-loop/complicated-access messages are expected for the enclosing 3-D nest.
- `fderivs_shift3` interior loop (`diff_new.f90:2030`) vectorized, confirming the V4 fusion retains SIMD.
- Ordinary `fderivs`/`fdderivs`: some interior `i` loops vectorize (for example lines 1004, 1252, 1322, 1668, and 1737), while boundary/symmetry loops at 1018, 1175--1176, 1425--1426, 1600--1601, 1805--1806, 1875--1876, and 1943--1944 miss due to control flow. Those loops cover full boxes in several variants, so the messages are not all negligible boundary-only work.
- `lopsided` has no vectorized loop reported. Its hot nest at lines 236--238 misses because the innermost loop contains control flow. Inspection shows per-point sign tests on `betax/betay/betaz` choose one-sided stencils. Given its 8.04% measured share, this is the largest clean SIMD missed opportunity.
- `diff_newwb.f90` also has many full-domain loops missed for control flow, but this file has user modifications and is read-only for V5.

Recommended SIMD candidate: restructure `lopsided` into branch-free/select form or separate sign-homogeneous paths while keeping `i` innermost. This has lower memory-traffic upside than the Ricci contraction candidate, but a clearer compiler proof and lower numerical risk.

## 6. AMR/MPI Communication Analysis

The primary transfer implementation already packs all variables in the supplied `MyList<var>` into one contiguous buffer per peer. It is therefore not true that every grid function is sent as an independent MPI message. It uses `MPI_Isend`/`MPI_Irecv` for peers and one `MPI_Waitall`, both in `Parallel.C` and `Parallel_bam.C`.

Important limitations:

- Packing, then posting all communication, then immediately waiting provides no useful computation/communication overlap.
- Send/receive buffers and request/status arrays are allocated and freed on every transfer. The size is first computed by traversing lists, then the same lists are traversed again to pack, and again to unpack.
- The design aggregates by peer and variable list, which avoids many tiny per-variable messages, but AMR geometry can still produce sparse peer buffers and repeated synchronization calls.
- There is one state synchronization after the predictor and one after the three correctors in the active `Step` path, not a ghost exchange after each of all four RK RHS calls. Restrict/prolong adds more transfers after each recursively scheduled level step.
- No `MPI_Barrier` appears in `Parallel.C` or `Parallel_bam.C`. Blocking `MPI_Send`/`MPI_Recv` exists in older/special-purpose paths, while the central transfer path is nonblocking followed immediately by `Waitall`. Global reductions occur for error checks and diagnostics.
- Per top-level step, `Step` itself executes 2 x 66 = 132 level synchronizations, plus restriction/prolongation dictated by the recursive schedule. Exact restrict/prolong transfer counts are path- and `YN`-dependent and were not instrumented in this audit.

The largest MPI/AMR issue is lack of overlap and repeated buffer/list-management overhead, not lack of multi-variable packing. It is a plausible later target, but altering it is invasive and touches the user's modified `Parallel.C`/`Parallel_bam.C`, so it is not recommended as the first V5 change.

## 7. Existing GPU Path Audit

`ABE` links `bssn_class.o` and the ordinary CPU sources. `ABEGPU` substitutes `bssn_gpu_class.o`/`bssn_step_gpu.o` and additionally links `bssn_gpu.o` and `bssn_gpu_rhs_ss.o`. The Python launcher selects `make ABEGPU` and runs `./ABEGPU` when GPU calculation is requested.

This is a real NVIDIA CUDA implementation, not a stub. It contains CUDA kernels for symmetry boundaries, first and second finite differences, KO dissipation, lopsided/advection derivatives, Ricci/RHS parts, and both Cartesian and shell-patch BSSN RHS. It does not use HIP or another portable GPU API.

Coverage and gaps:

- BSSN RHS: implemented as many CUDA kernels.
- Finite differences/advection/KO: implemented on CUDA.
- Symmetry boundary work: implemented on CUDA.
- RK4: still calls CPU Fortran `rungekutta4_rout`; `bssn_step_gpu.C` is orchestration, not a CUDA RK kernel.
- Prolong/restrict, general ghost exchange, and MPI: CPU code; no CUDA-aware MPI path was found.
- Device memory is allocated for many state, derivative, RHS, connection, Ricci, and constraint arrays.
- `bssn_gpu.cu` contains 27 textual H2D and 65 D2H `cudaMemcpy` sites; shell RHS contains 57 H2D and 56 D2H sites. The exact executed count depends on call mode, but full state inputs are copied to the device and RHS/diagnostic arrays copied back so CPU RK, boundaries, MPI, and AMR can proceed. This creates frequent whole-array traffic at RHS/RK boundaries.
- Numerous deprecated `cudaThreadSynchronize()` calls serialize kernel phases. There are many small kernels and explicit memset/synchronization operations.

Why it is not the default and is not presently realistic as a quick win: it requires an NVIDIA CUDA toolchain and separate executable/configuration; it uses old CUDA idioms; the current host/device ownership model copies many arrays every RHS; RK and AMR/MPI remain host-side; and no current correctness/performance evidence exists. Repairing it is a substantial port/modernization project, not a small enablement change. An AMD/HIP port would require replacing the CUDA runtime/build system and validating a large CUDA codebase; despite syntactic similarity, it is a major effort.

## 8. Local GPU Environment

- Hardware (queried through Windows from WSL): **NVIDIA GeForce RTX 5060 Laptop GPU**, 8,151 MiB; Windows NVIDIA driver 591.91.
- WSL2 exposes Microsoft virtual 3D controllers and `/dev/dxg`-style GPU virtualization, but Linux `nvidia-smi` currently reports `GPU access blocked by the operating system`.
- Windows `nvidia-smi.exe` can identify the GPU, so the host driver is installed and operational.
- `nvcc` is absent. Therefore CUDA compilation is not currently available inside this WSL environment.
- HIP tooling/path is not present or used by ABEGPU. The existing code specifically requires NVIDIA CUDA.
- No packages were installed.

## 9. Candidate Optimizations

| Priority | Candidate | Hot share affected | Expected whole-step gain | Difficulty | Numerical risk | Intrusion | MPI/AMR impact | Future GPU impact |
|---:|---|---:|---:|---|---|---|---|---|
| 1 | Contract metric second derivatives directly/tile scratch in `compute_rhs_bssn` | RHS 48.8%; `fdderivs` 6.9% direct plus RHS cache pressure | 3--8% | Medium | Low--medium | Local to Fortran RHS/new helper | None | Favorable: exposes point/tile kernel |
| 2 | Make `lopsided` inner loop SIMD-friendly | 8.0% | 2--5% | Medium | Low if algebraically identical | One kernel | None | Neutral/favorable |
| 3 | Reduce RK/state full-array traffic or fuse per-variable RK updates in C++ traversal | RK direct share modest, but 24 arrays x every stage | 1--4% | Medium | Low | C++/Fortran interface | Sync ordering must remain | Favorable if device-resident RK follows |
| 4 | Moderate fusion of lopsided/KO for small groups sharing shift/RHS traversal | 13.6% combined | 2--6% | Medium--high | Medium | Several kernels/call sites | None | Favorable if bounded, harmful if monolithic |
| 5 | Reuse MPI buffers and overlap interior computation with ghost exchange | synchronization/AMR share; `prolong3` 4.4% plus MPI | 2--10%, workload-dependent | High | Medium | High; touches modified files | Direct/high | Favorable if CUDA-aware design follows |
| 6 | Modernize and restore device-resident ABEGPU | RHS and derivatives, potentially most evolution | Potentially large, but unknown | Very high | High | Separate large code path | Requires redesign | Directly favorable but long project |

These estimates are engineering ranges, not benchmark results. A single CPU change is not expected to close the gap from 1.06x to 2x; that target likely requires several structural CPU changes or a redesigned device-resident GPU path.

## 10. Recommended First V5 Optimization

First change: add a narrowly scoped contracted-second-derivative helper in `AMSS_NCKU_source/diff_new.f90` and replace only the six metric `fdderivs` + immediate Ricci-contraction sequences in `AMSS_NCKU_source/bssn_rhs.f90`.

The helper should compute the scalar contraction

`gupxx*fxx + gupyy*fyy + gupzz*fzz + 2*(gupxy*fxy + gupxz*fxz + gupyz*fyz)`

directly (or in cache-sized `i`-contiguous tiles), rather than storing six full temporary derivative arrays and rereading them. It must preserve the current fourth-order stencils, symmetry handling, boundary behavior, loop order, and floating-point operation order as far as practical. Only one metric component should be converted first and validated before applying the same pattern to the other five; this avoids repeating the failed broad metric/Gamma fusion experiment.

Why this is first: it attacks proven high-frequency whole-box temporary traffic; stays within two unmodified Fortran files; leaves MPI/AMR ordering untouched; retains SoA/unit-stride layout; and creates a GPU-friendly contracted kernel. Expected whole-step benefit is **3--8%** for the complete six-component conversion, with an initial one-component experiment expected around **0.5--1.5%**. If measurement shows the contracted kernel loses SIMD or changes numerics, stop and promote the branch-free `lopsided` SIMD candidate instead.
