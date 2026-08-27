# AMSS-NCKU 第二阶段优化报告

## 1. 优化目标

第二阶段优化（V2）在 V1 的基础上继续优化高频插值和有限差分路径，而不是从 baseline 重新开始。目标是在保持物理配置、数值算法、插值阶数和差分阶数不变的前提下，减少固定六点插值与规则内部网格差分的执行开销，并通过链接时优化扩大跨源文件优化范围。

## 2. 测试环境

- 操作系统：Ubuntu 22.04 WSL
- MPI：OpenMPI
- 编译器：GCC / GFortran 11.4
- MPI 进程数：8
- Evolution steps：20
- 网格与物理配置：与 baseline 相同

## 3. V1 基础

V2 完整继承 V1 已有优化：

- 使用 `-O3`。
- 使用 `-march=native` 启用本机 CPU 指令优化。
- 使用 `-finline-functions` 加强 Fortran 函数内联。
- 保留 V1 的 `polint` 固定六点路径优化以及其他阶数的通用 fallback。

## 4. V2 优化内容

### 4.1 `polint6_fast` 专用六点插值

当前生产配置为 `ghost_width=3`，因此高频插值路径的阶数固定为 `2*ghost_width=6`。V2 新增固定长度 6 的 `polint6_fast` 内核，并将 `prolongrestrict_cell.f90` 和 `prolongrestrict_vertex.f90` 中对应的固定六点调用切换到该内核。

其他阶数仍使用通用 `polint` fallback，`polin2` 和 `polin3` 的通用路径也继续保留。该改动不改变 Neville 插值算法、插值阶数或物理计算逻辑。

### 4.2 标量 Neville 展开

`polint6_fast` 进一步将固定六点路径中的 `c`、`d`、`ho` 等中间量展开为标量，并展开固定次数的 Neville 递推。这样可减少数组切片、动态索引、内部循环和函数内部管理开销，同时保持原有数学计算过程不变。

### 4.3 LTO

在 C++、Fortran 编译及最终链接阶段统一启用 `-flto`，允许编译器进行跨源文件优化和内联。构建结果中已观察到 `polint6_fast` 被内联，最终二进制中不再保留独立符号。

### 4.4 `fderivs` 内部区域拆分

V2 修改实际参与构建的 `diff_new.f90`。在 `ghost_width=3` 对应的四阶中心差分配置下，将规则内部区域与边界区域分开处理：内部区域直接执行四阶中心差分，避免每个内部网格点重复进行边界条件判断；边界区域仍保留原有条件、对称边界处理和 fallback。

该优化不改变差分阶数、数学公式或边界计算逻辑。

## 5. 正确性保证

- 不改变物理算法、插值算法和差分公式。
- 不改变插值阶数与差分阶数。
- 不改变 MPI 进程数（8）与 Evolution steps（20）。
- 正式生效的配置已与 baseline 对照确认：`ABEtype=0`（BSSN vacuum）、`GAUGE=0`、`ghost_width=3`、`RPB=0`、`PSTR=0`，未启用错误的 Z4c 配置。
- `microdef.fh` 为正常文件，不是指向 `macrodef.fh` 的错误软链接；当前正式配置由 `macrodef.h` 引用 `macrodef.fh`。
- 完整 20 steps 运行成功，puncture position 按正确运行轨迹持续输出。
- 运行最终输出：`Simulation is successfully done!!`

## 6. 性能结果

| 版本 | Total Evolve | Total Running | real | 说明 |
| --- | ---: | ---: | ---: | --- |
| Baseline | 3852.320 s | 4009.230 s | 66m49.233s | 沿用 V1 正式报告数据 |
| V1 | 821.391 s | 824.376 s | 57m35.947s | V1 正式报告数据 |
| V2 Run #1 | 892.440 s | 895.866 s | 55m28.511s | 带 MPI profiling |
| V2 Run #2 | 861.726 s | 864.966 s | 56m07.120s | 带 MPI profiling |
| V2 average | 877.083 s | 880.416 s | 55m47.816s | 两次带 MPI profiling 运行的平均值 |

Run #1 的 `This run used` 为 895.928 s，Run #2 为 865.031 s。一次表现较好的 V2 内部计时为 `Total Running Time = 864.966 s`，但由于运行环境波动及 profiling 插桩会影响内部统计，不据此与 V1 的 824.376 s 得出 V2 计算更慢的绝对结论。

**V2 当前性能结果为带 MPI profiling 插桩的测试结果，不作为最终无插桩极限成绩。**

## 7. 当前加速效果

以两次 V2 profiling 的平均 real time 3347.816 s 计算：

- 当前带 profiling 插桩情况下，相比 Baseline：`4009.233 / 3347.816 = 1.198x`，约为 **1.20x** 加速，墙钟时间下降约 **16.5%**。
- 当前带 profiling 插桩情况下，相比 V1：`3455.947 / 3347.816 = 1.032x`，约为 **1.03x** 加速，墙钟时间下降约 **3.1%**。

这些数据是包含 profiling 额外开销的保守结果，不代表移除插桩后的 V2 最终极限性能。

## 8. Profiling 影响

测试时 `Parallel.C` 和 `Parallel_bam.C` 中仍包含 MPI PROFILE 插桩，涉及：

- `MPI_Wtime` 计时；
- `MPI_Gather` 汇总各 rank 数据；
- `MPI_Reduce` 统计 max / min / average；
- rank 0 输出 `MPI PROFILE REPORT`、`L2Norm` 及 rank 0 到 rank 7 的统计文本。

这些操作会引入额外通信、同步和文本输出开销，并干扰 real time。因此，本次归档不把 profiling instrumentation 作为正式性能优化代码提交；后续应移除插桩后重新进行干净 benchmark，作为 V2 的最终极限成绩。

## 9. 后续优化方向

1. 移除 MPI profiling 插桩，进行无插桩的干净 benchmark。
2. 继续分析 `compute_rhs_bssn_` 热点。
3. 探索多场 derivative kernel 融合，减少重复遍历。
4. 减少三维临时数组及其带来的内存带宽压力。
5. 探索通信与计算重叠，并复用 MPI 缓冲区。
