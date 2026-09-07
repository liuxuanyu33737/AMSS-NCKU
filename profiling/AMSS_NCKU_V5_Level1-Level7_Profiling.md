# AMSS-NCKU V5 CPU Profiling：Level-1 ～ Level-7

> **归档说明**
>
> - 归档分支：`optimize-v5`
> - 归档文件：`profiling/AMSS_NCKU_V5_Level1-Level7_Profiling.md`
> - 归档日期：2026-09-07
> - 内容性质：Level-1 ～ Level-7 逐层 CPU profiling 的过程记录与最终定位结论，作为后续 V5 针对性优化的依据。
> - 重要前提：以下所有 Level-1～Level-7 数据均来自**带 profiling 插桩的运行**，这些运行时间不能作为最终性能 benchmark。

## 1. Profiling 目的

本轮 profiling 的目的：

- 不是直接优化 `compute_rhs_bssn`；
- 而是按照“逐层定位”的方式，寻找 ABE 中真正执行时间异常/耗时较高的函数；
- 从 ABE Evolve 逐层向下定位；
- 最终找到可以进行针对性优化的具体函数。

明确说明：

- 所有 Level-1～Level-7 运行均带有 profiling 插桩，因此这些运行时间不能作为最终性能 benchmark。

## 2. Profiling 环境

记录本次 profiling 运行时的确认环境：

| 项目 | 值 |
| --- | --- |
| CPU | Intel Core i7-14650HX |
| CPU(s) | 24 |
| cores/socket | 12 |
| threads/core | 2 |
| MPI | OpenMPI 4.1.2 |
| GCC/GFortran | 11.4.0 |
| MPI ranks | 8 |
| GPU_Calculation | no |
| Evolution_Step_Number | 20 |
| OS | WSL Ubuntu 22.04.5 LTS |

编译优化中使用的关键选项（与 `AMSS_NCKU_source/makefile.inc` 中实际配置一致）：

- `-O3`
- `-march=native`
- `-flto`
- `-finline-functions`

（C++ 侧关键选项为 `-O3 -march=native -flto`；Fortran 侧为 `-O3 -march=native -flto -finline-functions`。）

## 3. Level-1：ABE 粗粒度 profiling

### 3.1 结果

| 阶段 | count | total(s) | avg(s) | maxrank(s) |
| --- | ---: | ---: | ---: | ---: |
| compute_rhs_bssn | 2112 | 129.446321 | 6.129087e-02 | 16.430768 |
| rk4 | 50688 | 3.315873 | 6.541732e-05 | 0.429344 |
| boundary | 52800 | 0.390180 | 9.438133e-06 | 0.057718 |
| enforce_ga | 2112 | 4.234896 | 2.005159e-03 | 0.539117 |
| sync | 2112 | 3.889536 | 1.841637e-03 | 0.529202 |
| amr | 808 | 21.528627 | 2.664434e-02 | 2.704379 |
| other | 2640 | 78.725026 | 2.982009e-02 | 9.851318 |

### 3.2 解释

- `compute_rhs_bssn` 仍然是主要阶段之一；
- 但 `other` 出现明显 long-tail；
- `other` 的 maxrank 接近 9.85 s；
- 因此不能只盯着 `compute_rhs_bssn`；
- 下一层继续拆分 `other`。

## 4. Level-2：拆分 other

### 4.1 结果

| 阶段 | count | total(s) | avg(s) | maxrank(s) |
| --- | ---: | ---: | ---: | ---: |
| bh_porg_rhs | 256 | 0.225053 | — | 0.030065 |
| bh_euler_pred | 256 | 0.000030 | — | 0.000005 |
| analysis | 8 | 78.153983 | 9.769248 | 9.769877 |
| swap_pre_cor | 1056 | 0.001640 | — | 0.000279 |
| swap_state_old | 528 | 0.001440 | — | 0.000207 |
| err_allreduce | 2112 | 19.848933 | 9.398169e-03 | 2.869079 |

（表中 “—” 表示该行插桩未单独记录该项。）

### 4.2 重点结论

**`analysis` 是 `other` 的主要异常来源。**

- `analysis` 的 total（≈78.15 s）与 maxrank（≈9.77 s）都显著高于 `other` 内其它子阶段；
- `err_allreduce` 虽也有一定 total（≈19.85 s），但其 maxrank 约 2.87 s，且属于 Allreduce 同步类开销，与 `analysis` 性质不同；
- 因此继续对 `analysis`（即 AnalysisStuff）进行下钻。

## 5. Level-3：拆分 AnalysisStuff

### 5.1 结果

| 阶段 | count | total(s) | avg(s) | maxrank(s) |
| --- | ---: | ---: | ---: | ---: |
| compute_psi4 | 8 | 0.417166 | 0.05214572 | 0.053106 |
| surf_wave | 96 | 11.105772 | 0.1156851 | 1.390047 |
| surf_masspang | 96 | 72.530521 | 0.7555263 | 9.066649 |
| psi4_write | 96 | 0.001579 | — | 0.001543 |
| map_write | 96 | 0.000360 | — | 0.000350 |
| bh_write | 8 | 0.000646 | — | 0.000645 |

### 5.2 重点结论

**`surf_MassPAng` 是 AnalysisStuff 的主要热点。**

说明：

- `surf_MassPAng` maxrank ≈ 9.07 s；
- AnalysisStuff maxrank ≈ 10.51 s；
- 因此 `surf_MassPAng` 约占 AnalysisStuff maxrank 的 86%（9.066649 / 10.51 ≈ 0.863）。

## 6. Level-4：拆分 surf_MassPAng

### 6.1 结果

| 阶段 | count | total(s) | avg(s) | maxrank(s) |
| --- | ---: | ---: | ---: | ---: |
| surf_masspang_total | 96 | 64.952504 | 6.765886e-01 | 8.119232 |
| volume_f_admmass | 96 | 0.488794 | 5.091602e-03 | 0.071750 |
| prep_arrays | 96 | 0.016583 | — | 0.002591 |
| interp_points | 96 | 64.423222 | 6.710752e-01 | 8.063309 |
| sphere_local_loop | 96 | 0.008597 | — | 0.001277 |
| mpi_allreduce_mass | 96 | 0.012032 | — | 0.002240 |
| mpi_allreduce_sx | — | 0.000435 | — | 0.000093 |
| mpi_allreduce_sy | — | 0.000184 | — | 0.000026 |
| mpi_allreduce_sz | — | 0.000121 | — | 0.000024 |
| mpi_allreduce_px | — | 0.000184 | — | 0.000025 |
| mpi_allreduce_py | — | 0.000098 | — | 0.000012 |
| mpi_allreduce_pz | — | 0.000100 | — | 0.000013 |
| post_process | — | 0.000013 | — | 0.000002 |
| cleanup | — | 0.001938 | — | 0.000394 |

（表中 “—” 表示该行插桩未单独记录该项。）

### 6.2 重点结论

**`Interp_Points` maxrank = 8.063309 s**

- `surf_MassPAng`（本层计时口径）maxrank = 8.119232 s；
- `Interp_Points` ≈ 8.063309 / 8.119232 ≈ 99.3%；
- 因此继续向 `Patch::Interp_Points()` 下钻。

## 7. Level-5：Patch::Interp_Points

### 7.1 运行概要

- nprocs = 8
- marked_interp_calls(sum) = 96
- points_fed(sum) = 3,538,944
- points_interp_local(sum) = 442,368

### 7.2 各阶段结果

| 阶段 | total(s) | avg(s) | maxrank(s) |
| --- | ---: | ---: | ---: |
| interp_setup | 0.101892 | 1.061371e-03 | 0.016060 |
| interp_point_loop | 29.613511 | 0.3084741 | 8.735395 |
| interp_block_search | 0.071251 | 7.422029e-04 | 0.012226 |
| interp_global_interp | 29.542260 | 0.3077319 | 8.723168 |
| interp_ar_shellf | 40.704748 | 0.4240078 | 8.790552 |
| interp_ar_weight | 0.036985 | 3.852645e-04 | 0.006306 |
| interp_post | 0.003421 | — | 0.000486 |
| interp_cleanup | 0.048506 | — | 0.006322 |
| interp_total | 70.509526 | 0.7344742 | 8.820580 |

（表中 “—” 表示该行插桩未单独记录该项；`interp_total` 为本次标记调用路径的整体计时口径。）

### 7.3 重点结论

- `interp_global_interp` 是主要计算热点；
- `interp_ar_shellf` 是主要通信热点；
- `interp_block_search` 几乎可以忽略；
- `shellf` Allreduce 是大数组 Allreduce，不能与 Level-4 的几个标量 Allreduce 混为一谈；
- 这一层同时发现计算热点和通信热点。

## 8. Level-6：f_global_interp

### 8.1 结果

`f_global_interp` 调用总次数 count = 7,520,256。

| 阶段 | total(s) | avg(s) | maxrank(s) |
| --- | ---: | ---: | ---: |
| gi_index | 0.348275 | 4.631162e-08 | 0.103672 |
| gi_decide3d | 1.484201 | 1.973605e-07 | 0.447901 |
| gi_polin3 | 29.998879 | 3.989077e-06 | 8.951208 |
| gi_total | 31.831356 | 4.232749e-06 | 9.502782 |

### 8.2 重点结论

**`polin3` 是 `global_interp` 的主要计算热点。**

- `gi_polin3` maxrank = 8.951208 s；
- `gi_polin3` 的 total（≈30.00 s）与 maxrank（≈8.95 s）均远高于同层其它阶段；
- 因此继续向 `polin3` 内部下钻。

## 9. Level-7：polin3 内部

### 9.1 运行概要

- polint_calls(sum) = 323,371,008
- polint_calls/polin3(avg) = 43.0000

说明（order 6）：36 + 6 + 1 = 43。

### 9.2 各阶段结果

| 阶段 | count | total(s) | avg(s) | maxrank(s) |
| --- | ---: | ---: | ---: | ---: |
| p3_prep | 7520256 | 0.301934 | 4.014944e-08 | 0.087499 |
| p3_z_x3a | 7520256 | 39.373214 | 5.235621e-06 | 11.465719 |
| p3_y_x2a | 7520256 | 8.155290 | 1.084443e-06 | 2.375968 |
| p3_x_x1a | 7520256 | 1.338965 | 1.780479e-07 | 0.389693 |
| p3_total | 7520256 | 53.702625 | 7.141063e-06 | 15.640531 |

### 9.3 重点结论

- `z` phase 是 polin3 内绝对主要阶段；
- `z` phase 进行 36 次 polint；
- `y` phase 6 次；
- `x` phase 1 次；
- 共 43 次；
- 36/43 ≈ 83.7%；
- 不是简单的 z 方向异常访存；
- 当前 Cell 分支下 `ya` 是 6×6×6 的小数组；
- 其大小约 1.7 KB，能够驻留 L1。

配套的 Level-8 源码分析（源码走读确认，非新一轮运行插桩 profiling）确认：

- polint 三个阶段的单次调用体基本相同；
- z/y/x 的主要区别是调用次数；
- z 方向主要由 polin3 的张量积循环结构导致；
- 当前 `surf_MassPAng → global_interp → polin3` 路径没有使用 `polint6_fast`；
- 当前实际调用的是通用 `polint`；
- `polint6_fast` 只在 prolong/restrict 路径使用。

## 10. 最终定位链

用一个清晰的树状结构总结：

```
ABE
└── AnalysisStuff
    └── surf_MassPAng
        └── Patch::Interp_Points
            └── global_interp
                └── polin3
                    └── 36 × polint (z)
                       + 6 × polint (y)
                       + 1 × polint (x)
```

最终结论：

**当前 V5 CPU profiling 已经把 AnalysisStuff 的主要计算热点逐层定位到了 `polin3()`，进一步定位到 z-phase 中重复执行的 36 次 order-6 `polint`。**

## 11. 下一步优化候选

> 以下只记录候选，不实施。

首选候选：

针对固定 order=6、相同 x3a/x3 的情况，为 z/y/x 三个阶段设计专用 order-6 插值路径：

- 共享插值权重；
- 减少重复最近节点搜索；
- 减少重复构造 c/d/ho/den；
- 中间阶段不计算无用 dy；
- 对 36 条 z 线复用同一组权重；
- 保留原 polint fallback。

明确说明：

- 这是候选优化方向，**尚未实施，也尚未 benchmark 验证**。

同时记录：

- Level-5 还发现 `interp_ar_shellf` 是一个独立的大数组 MPI_Allreduce 通信热点，未来可以作为第二条优化路线研究。

## 12. Profiling 结论

本轮不是简单地“看到 compute_rhs_bssn 占比高就优化 compute_rhs_bssn”。

而是通过：

Level-1 → Level-2 → Level-3 → Level-4 → Level-5 → Level-6 → Level-7

逐层缩小范围，最终定位到具体函数和具体计算阶段。

这份报告用于记录 V5 的 profiling 过程和后续优化依据。




