# AMSS-NCKU V1-V4 优化总结与正式性能复测

## 1. 优化背景与测试说明

本项目针对 AMSS-NCKU 的 CPU 演化路径进行性能优化。V1-V4 覆盖编译配置、插值、有限差分、BSSN RHS 导数计算和 KO 数值耗散；每项修改均以数学等价和端到端正确运行为前提，实测无稳定收益的实验予以回退。

正式复测统一使用以下条件：

- CPU，`GPU Calculation = no`
- MPI processes = 8
- Evolution Step Number = 20
- 同一机器、同一 WSL 环境
- 官方 Python 入口 `AMSS_NCKU_Program.py`
- Clean Baseline 使用原始 `-O3` release build，不包含 `-pg`
- V4 使用 `optimize-v4` 的正式 release build

本文严格区分四种计时：

1. `Total Evolve Time`：主要反映演化阶段耗时。
2. `Average Step Time = Total Evolve Time / 20`：作为暑期训练“单步耗时”的主要指标。
3. `This Program Cost`：官方 Python driver 的内部计时，作为 ASC 风格整体性能指标。
4. shell `real`：完整执行 `python3 AMSS_NCKU_Program.py` 的墙钟时间，包含更多 Python 前后处理、编译和绘图等过程，只作为辅助参考。

上述指标覆盖的时间区间不同，不能混合比较。

## 2. V1 优化

### 2.1 编译参数优化

原始 baseline 的 C++ 和 Fortran 编译均以 `-O3` 为基础。提交 `d4cf697` 在正式优化版本中加入：

- C++：`-march=native`
- Fortran：`-march=native -finline-functions`

这些参数使编译器可针对当前 CPU 指令集生成代码，并更积极地处理 Fortran 内联。若早期 profiling build 使用过 `-pg`，删除 `-pg` 只是从带插桩构建恢复到正常 release benchmark 环境，不属于算法或源码优化收益。

### 2.2 固定六点 `polint` 优化

早期 profiling 显示 `polint` 是主要热点之一。实际运行中，`ghost_width=3` 对应大量固定六点插值；V1 在通用 Neville 插值中为六点情形展开计算，减少循环控制和数组操作，其余阶数仍走原有通用路径。修改保持 Neville 插值方法、选点逻辑和数学结果等价，目的是降低超高频小函数的成本。

## 3. V2 优化

### 3.1 `polint6_fast`

提交 `7765e9c` 将固定六点标量 Neville 计算整理为专用 kernel `polint6_fast`，并在 cell-centered 与 vertex-centered prolongation/restriction 的适用位置调用 fast path；不满足固定六点条件时仍保留 generic fallback。

### 3.2 LTO

提交 `119ab9e` 在 C++、Fortran 编译参数和链接阶段加入 `-flto`，目标是让编译器跨 Fortran 源文件内联小函数并进行整体优化。最终二进制中未再看到独立的 `polint6_fast` 符号，可作为其可能已被跨文件内联的侧面证据，但不能单独视为绝对证明。

### 3.3 `fderivs` interior/boundary 分离

提交 `b706ed5` 调整 `diff_new.f90` 中的 `fderivs`：规则内部区域采用固定四阶中心差分 fast path，边界区域继续使用原有合法性检查和降阶处理。数值格式没有改变，结构调整的目的是减少 interior 热循环中的逐点条件判断。

## 4. V3 优化探索

V4 实际从提交 `7316333` 之后分叉。以下 metric、Gamma 和 lopsided 实验保存在 `optimize-v3` 分支，用于记录 V3 阶段的性能探索；由于最终测试无稳定收益，相关源码修改被回退，因此没有进入 `optimize-v4` 的实际祖先链。`328d889` 是该独立分支上的 V3 总结提交。

### 4.1 最终进入 V4：shift derivative fusion

提交 `7316333`（`Fuse shift derivative passes in BSSN RHS`）融合 `betax`、`betay`、`betaz` 三个 shift 分量的一阶导数计算，减少重复网格遍历和循环调度。该修改通过随机数组、不同 symmetry 条件的离线等价测试以及完整 20-step end-to-end 运行，最终保留，并成为 V4 的祖先提交。

### 4.2 metric derivative fusion

- 实验提交：`871ab07`
- 回退提交：`3018a69`
- 分支位置：`optimize-v3`，不在 `optimize-v4` 祖先链

该实验数学结果正确，但实际 benchmark 变慢，因此回退。较大的循环体、寄存器压力、cache 行为或编译器生成代码均可能影响结果，但这些只是可能因素，尚无严格归因。

### 4.3 Gamma derivative fusion

- 实验提交：`e7bc6c8`
- 回退提交：`b9e0edb`
- 分支位置：`optimize-v3`，不在 `optimize-v4` 祖先链

该实验通过数值检查，但 benchmark 没有性能收益，因此回退。

### 4.4 lopsided interior fast path

- 实验提交：`7c2a357`
- 回退提交：`24f03a4`
- 分支位置：`optimize-v3`，不在 `optimize-v4` 祖先链

离线正确性测试通过，但完整 20-step 测试中的 `Total Running Time` 没有稳定收益。考虑到 WSL 环境下 shell `real` 存在波动，最终依据更稳定的内部演化指标决定不保留该修改。

V3 的原则是：正确但没有性能收益的修改同样回退，并保留 Git 负实验记录。

## 5. V4 优化

### 5.1 `fdderivs` center cache 尝试

- 实验提交：`448e828`（`Reduce redundant work in fdderivs`）
- 回退提交：`91d56ef`

该实验尝试缓存 `fh(i,j,k)`，避免纯二阶导数 `fxx/fyy/fzz` 表面上的重复中心点读取。检查 GCC `-O3` 生成的汇编后，修改前后的关键代码基本等价，说明编译器已经完成相应的公共子表达式消除。由于没有形成机器码层面的实际优化价值，未再耗费一次长期 20-step benchmark，直接回退。这是一项 V4 历史实验，最终没有净源码改动。

### 5.2 KO dissipation kernel

提交 `8433b47`（`Optimize KO dissipation kernel`）是 V4 当前保留的主要源码优化。原 `kodis` 主循环在每个网格点重复判断 stencil 是否落在有效范围；修改后预先计算 `ilo/ihi`、`jlo/jhi`、`klo/khi`，只遍历 radius-3 stencil 合法区域，将边界合法性判断移出 hot loop。

KO stencil 保持为：

```text
1, -6, 15, -20, 15, -6, 1
```

修改不改变数学 stencil、系数、边界语义，也不主动改变除法运算顺序。汇编分析显示原内层循环中的多次边界条件分支被移出主循环，GCC 因而更容易生成 AVX vector loop。随机数组和 symmetry/parity 离线测试得到 `max_abs_error = 0`，完整 20-step 运行成功，因此保留。

## 6. Profiling 演变

| Routine | Early baseline | V3 |
| --- | ---: | ---: |
| `compute_rhs_bssn_` | 39.66% | 48.82% |
| `polint_` | 23.67% | 11.33% |
| `lopsided_` | 6.89% | 8.04% |
| `fdderivs_` | 5.81% | 6.93% |
| `kodis_` | 4.52% | 5.59% |
| `prolong3_` | 4.95% | 4.36% |
| `fderivs_` | 4.24% | 4.32% |

`polint_` 的相对占比由约 23.67% 降至约 11.33%，说明前两阶段针对插值路径的优化明显改变了热点结构；与此同时，`compute_rhs_bssn_` 的占比达到约 48.82%，后续最大优化目标已转向 BSSN RHS。

gprof 百分比是单次 profiling 中的相对占比。不同 profiling run 的 self seconds 和百分比不能代替同口径的正式 benchmark，也不能仅凭占比变化推导绝对耗时变化。

## 7. 前期性能记录错误与正式修正

项目早期曾记录 baseline `Total Evolve Time` 约 3852.32 s、某优化版本约 800 s，并一度据此形成约 4-5x speedup 的说法。该结论现已作废，不能作为正式性能结果，原因包括：

1. 不同阶段曾混用 profiling build 与 release build。
2. 部分旧数据可能包含 `-pg` profiling overhead。
3. `Total Evolve Time`、`Total Running Time`、shell `real` 和 `This Program Cost` 的计时范围不同。
4. 早期记录曾把不同计时口径放在一起比较。
5. WSL + MPI 环境下 shell `real` 存在明显波动。
6. 不同日期和时段的 benchmark 不足以可靠判断几个百分点以内的小优化。
7. 删除 `-pg` 属于恢复正常 benchmark 环境，不能计入算法或源码优化收益。

从 V4 收尾开始，项目重新建立 clean baseline。正式 Clean Baseline 使用提交 `2b54b0e`：原始源码、`-O3`、无 `-pg`、8 MPI、CPU、20 steps、官方 Python driver，并与 V4 在同一机器和 WSL 环境运行。V4 对应源码提交 `8433b47`，使用相同 MPI 数、设备、步数、driver 和环境。今后的正式性能结果均以这组 clean benchmark 为准。

## 8. Clean Baseline 与 V4 正式性能对比

### Clean Baseline（`2b54b0e`）

- `Total Evolve Time`：823.306 s
- `Average Step Time`：823.306 / 20 = 41.1653 s/step
- `Total Running Time`：826.269 s
- `This run used`：826.343 s
- `This Program Cost`：3744.7120990753174 s
- shell time：`real 64m33.304s`，`user 124m53.245s`，`sys 7m14.091s`

### V4（`8433b47`）

- `Total Evolve Time`：774.582 s
- `Average Step Time`：774.582 / 20 = 38.7291 s/step
- `Total Running Time`：777.368 s
- `This run used`：777.455 s
- `This Program Cost`：3675.8069472312927 s
- shell time：`real 63m19.342s`，`user 116m58.127s`，`sys 7m35.982s`

| Metric | Clean Baseline | V4 | Time Reduction |
| --- | ---: | ---: | ---: |
| Total Evolve Time | 823.306 s | 774.582 s | 5.92% |
| Average Step Time | 41.1653 s | 38.7291 s | 5.92% |
| Total Running Time | 826.269 s | 777.368 s | 5.92% |
| This Program Cost | 3744.712 s | 3675.807 s | 1.84% |
| Full Python real | 64m33.304s | 63m19.342s | 1.91% |

对应 speedup factor 为：

- Average Step：`41.1653 / 38.7291 ≈ 1.0629x`
- Total Evolve：`823.306 / 774.582 ≈ 1.0629x`
- Total Running：`826.269 / 777.368 ≈ 1.0629x`
- This Program Cost：`3744.712099 / 3675.806947 ≈ 1.0187x`
- Full Python real：`3873.304 / 3799.342 ≈ 1.0195x`

其中 `64m33.304s = 3873.304 s`，`63m19.342s = 3799.342 s`。时间降低比例与 speedup factor 是不同概念；5.92% 的时间降低不能表述为“提升 5.92 倍”。

## 9. 当前正式结论

clean benchmark 下，平均单步时间由 41.1653 s 降至 38.7291 s，时间降低约 5.92%，speedup 约为 1.0629x。目前尚未达到训练任务要求的 2x 单步加速。

ASC 风格的 `This Program Cost` 由 3744.712 s 降至 3675.807 s，时间降低约 1.84%，speedup 约为 1.0187x。

因此，当前可复现、可解释的正式结果是真实约 1.063x 单步加速，而不是前期因口径混用形成的 4-5x 记录。

## 10. V1-V4 涉及文件和 Git 历史

文件清单依据各提交的实际 diff 核对。`git diff --name-only 2b54b0e..8433b47` 还会包含 baseline 记录提交 `e51f6c6` 引入的 `BENCHMARK.md`、`macrodef.*`、`microdef.fh` 等文件；它们不应仅因出现在区间 diff 中就被误记为 V1-V4 的优化源码。本节按实际优化提交归类。

### 10.1 最终进入 `optimize-v4` 的修改

| 阶段 | 提交 | 实际修改的 tracked 文件 |
| --- | --- | --- |
| V1 | `d4cf697` | `AMSS_NCKU_source/makefile.inc`、`AMSS_NCKU_source/fmisc.f90`、`AMSS_NCKU_source/.gitignore`、`reports/OPTIMIZATION_V1_REPORT.md` |
| V2 | `7765e9c` | `AMSS_NCKU_source/fmisc.f90`、`AMSS_NCKU_source/prolongrestrict_cell.f90`、`AMSS_NCKU_source/prolongrestrict_vertex.f90` |
| V2 | `119ab9e` | `AMSS_NCKU_source/makefile.inc` |
| V2 | `b706ed5` | `AMSS_NCKU_source/diff_new.f90` |
| V2 | `58e16b7` | `OPTIMIZATION_V2_REPORT.md` |
| V3 | `7316333` | `AMSS_NCKU_source/bssn_rhs.f90`、`AMSS_NCKU_source/diff_new.f90` |
| V4 | `8433b47` | `AMSS_NCKU_source/kodiss.f90` |

最终进入 V4 的净优化源码文件为：

```text
AMSS_NCKU_source/bssn_rhs.f90
AMSS_NCKU_source/diff_new.f90
AMSS_NCKU_source/fmisc.f90
AMSS_NCKU_source/kodiss.f90
AMSS_NCKU_source/makefile.inc
AMSS_NCKU_source/prolongrestrict_cell.f90
AMSS_NCKU_source/prolongrestrict_vertex.f90
```

### 10.2 历史实验与回退

| 阶段 | 实验与回退 | 实际修改的 tracked 文件 | 状态 |
| --- | --- | --- | --- |
| V3 | `871ab07` / `3018a69` | `AMSS_NCKU_source/bssn_rhs.f90`、`AMSS_NCKU_source/diff_new.f90` | metric fusion 历史实验，最终回退，且未进入 V4 祖先链 |
| V3 | `e7bc6c8` / `b9e0edb` | `AMSS_NCKU_source/bssn_rhs.f90`、`AMSS_NCKU_source/diff_new.f90` | Gamma fusion 历史实验，最终回退，且未进入 V4 祖先链 |
| V3 | `7c2a357` / `24f03a4` | `AMSS_NCKU_source/lopsidediff.f90` | lopsided fast path 历史实验，最终回退，且未进入 V4 祖先链 |
| V4 | `448e828` / `91d56ef` | `AMSS_NCKU_source/diff_new.f90` | center cache 历史实验，最终回退 |

当前工作区未提交的 `Parallel.C`、`Parallel_bam.C`、`diff_newwb.f90` 以及测试目录、日志、profiling 输出和备份文件均不是上述正式优化提交的成果。

## 11. 后续优化方向

V3 profiling 中 `compute_rhs_bssn_` 已占约 48.82%，是下一阶段最值得分析的热点。后续工作建议包括：

- 对 `compute_rhs_bssn` 内部进一步分段计时。
- 检查重复网格遍历和临时数组成本。
- 分析 memory bandwidth、cache 和向量化情况。
- 寻找能够安全融合、又不会重现 metric/Gamma 大循环性能回退的计算。
- 分析 MPI / AMR synchronization 的占比与等待时间。

这些内容仅作为 V5 方向记录；本次不修改任何 V5 源码。所有后续实验仍应以数学等价、离线正确性测试和完整端到端验证为前提，并继续沿用本文定义的 clean benchmark 口径。
