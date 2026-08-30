# AMSS-NCKU V1–V4 优化总结与正式基准结果

## 1. 文档目的

本文总结 `optimize-v1` 至 `optimize-v4` 的主要优化、实验回退和正式性能结果。

当前 V4 以提交 `8433b47` 为关键版本。本文只把经过验证并保留在当前 `optimize-v4` 分支中的改动视为最终优化；曾经提交但因无收益或负收益而回退的实验，会单独标明，不计入最终 V4。

## 2. V1–V4 优化内容

### 2.1 V1：编译参数与 `polint` 六点专用优化

主要提交：

- `d4cf697`：`Optimization V1: compiler flags and polint optimization`
- `a73064c`：`Move and translate V1 optimization report`（V1 分支上的报告整理提交，不在当前 V4 的祖先链中）

主要工作：

- 调整 `AMSS_NCKU_source/makefile.inc` 中的编译优化参数，使正式运行使用面向性能的编译配置。
- 针对程序中高频出现的固定六点插值路径优化 `polint`，减少通用实现带来的循环、分支和临时工作。
- 保持原有数值算法、插值阶数和调用语义不变。

涉及的主要文件：

- `AMSS_NCKU_source/makefile.inc`
- `AMSS_NCKU_source/fmisc.f90`
- `AMSS_NCKU_source/.gitignore`
- `reports/OPTIMIZATION_V1_REPORT.md`

### 2.2 V2：`polint6_fast`、LTO 与 `fderivs` 常规内部区域拆分

主要提交：

- `7765e9c`：`Optimize fixed-order interpolation with scalar Neville kernel`
- `119ab9e`：`Enable LTO for cross-file Fortran inlining`
- `b706ed5`：`Split regular fderivs interior from boundary handling`
- `58e16b7`：`Document and finalize optimization V2`

主要工作：

- 新增并使用六点定阶标量 Neville 插值核心 `polint6_fast`，将已知阶数的热点调用从通用数组式 `polint` 路径中拆出。
- 在 cell-centered 和 vertex-centered 的 prolongation/restriction 调用处接入六点快速路径。
- 启用 LTO（链接时优化），帮助编译器跨 Fortran 源文件内联和消除额外调用开销。
- 将 `fderivs` 的常规内部网格与边界处理拆分：内部区域使用无逐点边界判断的紧凑热循环，边界区域继续保留原有合法性检查和降阶逻辑，从而减少热点循环中的条件判断并改善向量化机会。

涉及的主要文件：

- `AMSS_NCKU_source/fmisc.f90`
- `AMSS_NCKU_source/prolongrestrict_cell.f90`
- `AMSS_NCKU_source/prolongrestrict_vertex.f90`
- `AMSS_NCKU_source/makefile.inc`
- `AMSS_NCKU_source/diff_new.f90`
- `AMSS_NCKU_source/macrodef.fh`
- `AMSS_NCKU_source/macrodef.h`
- `AMSS_NCKU_source/microdef.fh`
- `OPTIMIZATION_V2_REPORT.md`

### 2.3 V3：保留 shift 三分量导数融合，回退负优化

最终保留提交：

- `7316333`：`Fuse shift derivative passes in BSSN RHS`

保留内容：

- 在 `compute_rhs_bssn` 路径中，将 `betax`、`betay`、`betaz` 三个 shift 分量原先彼此独立的导数调用融合为 `fderivs_shift3`。
- 在同一组网格循环中生成三个分量的九个一阶导数，减少重复的网格遍历、子程序调度与边界准备工作，提高数据局部性。
- 常规内部区域与边界区域仍分别处理，以保持原有边界和对称性语义。

涉及的主要文件：

- `AMSS_NCKU_source/bssn_rhs.f90`
- `AMSS_NCKU_source/diff_new.f90`

V3 还验证了以下结构性尝试，但实测没有带来稳定收益或出现负收益，因此均已回退，不属于最终 V4：

- metric 导数融合：`871ab07`，由 `3018a69` 回退。
- Gamma 导数融合：`e7bc6c8`，由 `b9e0edb` 回退。
- `lopsided` 内部有限差分优化：`7c2a357`，由 `24f03a4` 回退。
- `328d889` 为 V3 分支的整理报告提交，不在当前 V4 的祖先链中。

这些回退说明：减少循环次数并不必然缩短总时间；寄存器压力、内存带宽、编译器向量化结果和缓存行为都可能抵消理论收益。最终版本只保留经过正式测试有稳定正收益的 shift 三分量融合。

### 2.4 V4：回退 `fdderivs/fcenter` 等价尝试，保留 KO 热循环优化

实验并回退：

- `448e828`：`Reduce redundant work in fdderivs`
- `91d56ef`：`Revert "Reduce redundant work in fdderivs"`

该尝试希望减少 `fdderivs`/`fcenter` 路径中的表面重复计算，但检查 GCC 生成的汇编后发现关键代码与原实现等价，没有形成可验证的有效优化，因此回退，避免保留无收益改动。

最终保留提交：

- `8433b47`：`Optimize KO dissipation kernel`

保留内容：

- 优化 `kodis`/Kreiss–Oliger 数值耗散热点循环。
- 预先计算满足六阶 stencil 的有效 `i/j/k` 上下界，将原来每个网格点上的六个边界条件判断外提到循环外。
- 热循环只遍历合法内部区域，减少分支开销，并为编译器生成更好的向量化代码创造条件。
- 数值 stencil、耗散系数及边界语义保持不变。

涉及文件：

- `AMSS_NCKU_source/kodiss.f90`

## 3. 前期性能记录口径错误说明

前期曾记录或讨论过类似：

```text
3852 s -> 800 s
约 4–5x 加速
```

这组数字不能作为 V1–V4 的正式加速结论，主要原因如下：

1. 混用了旧 baseline、profiling 构建和不同阶段的运行数据，比较对象并非严格一致。
2. 混淆了程序内部的 `Total Evolve Time`、`Total Running Time`，shell 的 `time real`，以及官方流程输出的 `This Program Cost`。这些指标覆盖的代码范围不同，不能直接交叉相除。
3. 一些结果来自不同运行环境、不同代码时间点或不同系统负载，运行噪声和环境差异没有被排除。

因此，旧的 `3852 s -> 800 s` 和 `4–5x` 只能作为排查过程中的历史记录，不能写成正式性能结论，也不能用来判断是否达到训练任务要求。

## 4. 正式基准口径

正式比较必须满足以下条件：

- Clean Baseline 使用原始源码。
- 编译参数为 `-O3`，不使用 `-pg` profiling 插桩。
- Baseline 与 V4 均使用 8 个 MPI 进程。
- 均运行 20 个演化步。
- 使用同一个 Python 入口和相同完整流程。
- V4 的全部指标取自同一次带 shell `time` 的正式运行，不能拼接不同运行的数字。

建议从同一次日志中核对内部指标：

```bash
grep -Ei "Total Evolve|Total Running|This run used|This Program Cost" RUN_LOG
```

完整 Python 墙钟时间由同一次命令外层的 shell `time` 记录：

```bash
/usr/bin/time -p python3 YOUR_OFFICIAL_ENTRY.py 2>&1 | tee RUN_LOG
```

其中入口脚本和参数应按正式任务实际使用的命令替换，Baseline 与 V4 必须完全一致。

## 5. Clean Baseline 与 V4 正式结果

| 指标 | Clean Baseline | V4 | 变化 |
|---|---:|---:|---:|
| `Total Evolve Time` | 823.306 s | 774.582 s | 降低约 5.92% |
| 平均单步时间（20 steps） | 41.1653 s/step | 38.7291 s/step | 降低约 5.92%，speedup 约 1.0629x |
| `Total Running Time` | 826.269 s | 777.368 s | 降低约 5.92% |
| `This run used` | 826.343 s | 777.455 s | 降低约 5.92% |
| `This Program Cost` | 3744.712099 s | 3675.806947 s | 降低约 1.84%，speedup 约 1.0187x |
| 完整 Python `real` | 64m33.304s | 63m19.342s | 减少 73.962 s，降低约 1.91% |

计算方法：

```text
Baseline 平均单步 = 823.306 / 20 = 41.1653 s/step
V4 平均单步       = 774.582 / 20 = 38.7291 s/step

单步 speedup      = 823.306 / 774.582 = 1.0629x
演化时间降低比例  = (823.306 - 774.582) / 823.306 ≈ 5.92%

Program Cost speedup = 3744.712099 / 3675.806947 ≈ 1.0187x
Program Cost 降低比例 = (3744.712099 - 3675.806947) / 3744.712099 ≈ 1.84%

完整 Python real：
64m33.304s = 3873.304 s
63m19.342s = 3799.342 s
(3873.304 - 3799.342) / 3873.304 ≈ 1.91%
```

## 6. 当前结论与后续方向

在严格一致的 8 MPI、20 steps、`-O3`、无 `-pg`、相同 Python 完整流程下，V4 将平均单步时间从 `41.1653 s` 降至 `38.7291 s`，获得约 `1.0629x` 单步加速，演化阶段耗时降低约 `5.92%`。

该结果说明 V1–V4 的保留优化具有稳定正收益，但当前正式结果并未达到训练要求的 `2x` 加速。完整官方流程的 `This Program Cost` 仅降低约 `1.84%`，也表明演化内核之外仍存在大量耗时。

后续不宜继续只做零散小改，应重新进行同口径 profiling，优先针对 `compute_rhs_bssn` 等主要热点开展结构级优化，例如：

- 细分 `compute_rhs_bssn` 内部各类导数、张量运算和临时数组的耗时。
- 检查重复网格遍历、重复表达式与不必要的中间数组。
- 分析内存访问、缓存命中、向量化报告和寄存器压力。
- 区分计算热点与 MPI/边界同步开销，并分别验证。

所有后续优化都应继续使用本文定义的 clean baseline 和同流程测试方法，分别报告 `Total Evolve Time`、平均单步时间、`Total Running Time`、`This Program Cost` 与完整 Python `real`，避免再次混淆统计口径。

## 7. 当前 `optimize-v4` 已提交的主要优化文件

相对原始主分支，当前 `optimize-v4` 已提交的优化源码和报告主要包括：

```text
AMSS_NCKU_source/.gitignore
AMSS_NCKU_source/bssn_rhs.f90
AMSS_NCKU_source/diff_new.f90
AMSS_NCKU_source/fmisc.f90
AMSS_NCKU_source/kodiss.f90
AMSS_NCKU_source/macrodef.fh
AMSS_NCKU_source/macrodef.h
AMSS_NCKU_source/makefile.inc
AMSS_NCKU_source/microdef.fh
AMSS_NCKU_source/prolongrestrict_cell.f90
AMSS_NCKU_source/prolongrestrict_vertex.f90
OPTIMIZATION_V2_REPORT.md
reports/OPTIMIZATION_V1_REPORT.md
```

这些文件的优化改动已经存在于分支历史中，无需在本次文档提交中重复暂存。工作区中的 `Parallel.C`、`Parallel_bam.C`、`diff_newwb.f90`，以及测试目录、日志、profiling 输出和备份文件，不属于本次总结提交。
