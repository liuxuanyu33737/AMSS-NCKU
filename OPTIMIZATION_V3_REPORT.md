# AMSS-NCKU CPU 优化 V3 报告

## 1. V3 优化目标

V3 建立在 V1/V2 的正式优化基础上，重点继续分析 BSSN RHS 中高频有限差分路径。目标是在不改变物理公式、MPI 配置、边界条件、差分阶数、`GAUGE`、`ABEtype` 和 `ghost_width` 的前提下，通过减少重复网格遍历降低计算开销，并使用数值对比和 20-step A/B benchmark 判断优化是否值得保留。

V3 不以“减少 loop 数量”作为成功标准。每项实验必须同时满足：

- 离线数值结果与原实现一致；
- 完整 20-step 模拟成功；
- 性能数据能够支持保留该改动。

## 2. V3 开始前状态

V3 从提交 `58e16b7 Document and finalize optimization V2` 继续。生产配置保持：

- Ubuntu 22.04 WSL；
- OpenMPI，8 ranks；
- GCC / GFortran 11.4；
- `ABEtype=0`；
- `GAUGE=0`；
- `ghost_width=3`；
- 20 evolution steps；
- 167 个 grid functions。

V1/V2 已经完成编译参数、固定六点 `polint`、LTO 和 `fderivs` 规则内部区域等优化。V3 主要研究 BSSN RHS 中进一步融合 derivative pass 的可行性，并在新 profiling 结果上检查 `lopsided`。

## 3. 第一刀：Shift 导数融合

### 3.1 热点背景

`compute_rhs_bssn` 原先分别对 `betax`、`betay`、`betaz` 调用三次 `fderivs`。三个场采用相同网格和差分阶数，因此可在保持各自 symmetry 处理的前提下，合并主要导数遍历。

该优化对应提交：

```text
7316333 Fuse shift derivative passes in BSSN RHS
```

### 3.2 原实现

三个 shift 分量的 symmetry parity 分别为：

| 场 | x | y | z |
| --- | --- | --- | --- |
| `betax` | ANTI | SYM | SYM |
| `betay` | SYM | ANTI | SYM |
| `betaz` | SYM | SYM | ANTI |

原实现分别调用三次 `fderivs`。规则内部区域有 3 个 loop nest，边界区域也有 3 个 loop nest。

### 3.3 优化方式

新增专用 kernel `fderivs_shift3`，一次处理三个 shift 场，输出：

```text
betaxx betaxy betaxz
betayx betayy betayz
betazx betazy betazz
```

新实现保持：

- 三次独立 `symmetry_bd`；
- 独立扩展数组 `fhx`、`fhy`、`fhz`；
- 原四阶中心差分 stencil；
- 原边界判断和二阶 fallback；
- 原数值计算路径和物理公式。

循环结构由内部 `3 -> 1`、边界 `3 -> 1`。该改动减少了重复 loop control 和网格遍历，但没有合并不同场的 symmetry 操作。

### 3.4 数值正确性

在 `11 x 10 x 9` 随机三维网格上，将原三次独立 `fderivs` 与 `fderivs_shift3` 的 9 个输出逐项比较：

| Symmetry | max_abs_error |
| ---: | ---: |
| 0 | 0 |
| 1 | 0 |
| 2 | 0 |

### 3.5 20-step 性能结果

第一次正式测试：

| 指标 | 结果 |
| --- | ---: |
| Total Evolve Time | 756.876 s |
| Total Running Time | 759.967 s |
| This run used | 760.010 s |
| real | 55m48.388s |
| user | 105m48.673s |
| sys | 1m14.477s |

最终 puncture position：

```text
no.0 (-0.595187 4.40883 1.77552e-05)
no.1 (0.779982 -5.45474 2.96487e-05)
Simulation is successfully done!!
```

在第四刀评估阶段，又从 `7316333` 建立独立 worktree，并在同期环境重新 clean build 和运行：

| 指标 | 同期 shift3 复测 |
| --- | ---: |
| Total Evolve Time | 803.028 s |
| Total Running Time | 805.984 s |
| This run used | 806.034 s |
| real | 54m39.262s |
| user | 107m57.314s |
| sys | 1m16.036s |

### 3.6 结论

Shift derivative fusion 数值完全一致，并作为 V3 唯一新增的正式有效优化保留。报告同时列出两次 shift3 结果，以呈现运行环境波动，不只选择其中最快的一次。

## 4. 第二刀：Metric6 融合实验

### 4.1 思路

提交 `871ab07 Fuse metric derivative passes in BSSN RHS` 融合以下六个 metric 场的一阶导数：

```text
gxx gxy gxz gyy gyz gzz
```

共生成 18 个输出，将内部 loop nest 从 `6 -> 1`，边界 loop nest 从 `6 -> 1`，同时保持每个 metric 分量原有的 symmetry 参数。

### 4.2 数值验证

`11 x 10 x 9` 随机网格的验证结果为：

| Symmetry | max_abs_error |
| ---: | ---: |
| 0 | 0 |
| 1 | 0 |
| 2 | 0 |

数值正确性通过。

### 4.3 性能结果

| 指标 | shift3 | metric6 | 变化 |
| --- | ---: | ---: | ---: |
| Total Evolve | 756.876 s | 786.662 s | +29.786 s |
| Total Running | 759.967 s | 789.869 s | +29.902 s，约 +3.9% |
| This run | 760.010 s | 789.907 s | +29.897 s |
| real | 55m48.388s | 56m59.996s | +71.608 s，约 +2.1% |
| user | 105m48.673s | 110m25.017s | +4m36.344s |
| sys | 1m14.477s | 1m23.001s | +8.524 s |

### 4.4 回滚原因

该实验虽然减少了 loop nest 且数值正确，但内部计时和 wall-clock 均退化，因此通过提交 `3018a69 Revert "Fuse metric derivative passes in BSSN RHS"` 回滚。

可能原因包括寄存器压力增加、单个 kernel 过大、数组读写或写回压力增加，以及 cache locality 收益不足以抵消这些成本。这些是基于实现结构和性能现象的可能解释，尚未由硬件性能计数器证明。

## 5. 第三刀：Gamma3 融合实验

### 5.1 思路

提交 `e7bc6c8 Fuse Gamma derivative passes in BSSN RHS` 融合 `Gamx`、`Gamy`、`Gamz` 三个 Gamma 场。

其 symmetry parity 为：

| 场 | x | y | z |
| --- | --- | --- | --- |
| `Gamx` | ANTI | SYM | SYM |
| `Gamy` | SYM | ANTI | SYM |
| `Gamz` | SYM | SYM | ANTI |

新增 `fderivs_gamma3`，输出：

```text
Gamxx Gamxy Gamxz
Gamyx Gamyy Gamyz
Gamzx Gamzy Gamzz
```

内部和边界 loop nest 均由 `3 -> 1`，三次 `symmetry_bd` 保持不变。

### 5.2 数值验证

| Symmetry | max_abs_error |
| ---: | ---: |
| 0 | 0 |
| 1 | 0 |
| 2 | 0 |

数值正确性通过。

### 5.3 性能结果

| 指标 | shift3 | Gamma3 | 变化 |
| --- | ---: | ---: | ---: |
| Total Evolve | 756.876 s | 781.571 s | +24.695 s |
| Total Running | 759.967 s | 784.410 s | +24.443 s，约 +3.2% |
| This run | 760.010 s | 784.466 s | +24.456 s |
| real | 55m48.388s | 56m54.659s | +66.271 s，约 +2.0% |
| user | 105m48.673s | 109m32.599s | +3m43.926s |
| sys | 1m14.477s | 1m15.135s | +0.658 s |

### 5.4 回滚原因

Gamma3 同样在数值正确的情况下出现性能退化，因此通过提交 `b9e0edb Revert "Fuse Gamma derivative passes in BSSN RHS"` 回滚。该结果进一步说明，小于 metric6 的三场融合也不能仅凭 loop 数量减少推断会加速；寄存器压力、kernel 大小、cache 和 memory traffic 都可能抵消融合收益。

## 6. V3 重新 Profiling

### 6.1 gprof 方法

在 V3 最优的 shift3 状态重新进行 gprof：

- 编译参数为 `-O3 -march=native -flto -finline-functions -pg`；
- 使用 8 个 MPI ranks；
- 使用 `GMON_OUT_PREFIX` 分别收集各 rank 数据；
- 将 8 个 rank 的 gmon 数据合并为 `gmon.sum`。

gprof 的 self seconds 是多个 MPI rank 合并后的采样结果。因此不能直接将合并 profile 中约 3838 秒的累计时间与程序输出的约 760 秒 `Total Running` 作绝对时间比较。本节主要比较百分比、热点排序和调用次数。

### 6.2 新 Top Hotspots

| 排名 | 函数 | 占比 | self seconds | calls |
| ---: | --- | ---: | ---: | ---: |
| 1 | `compute_rhs_bssn_` | 48.82% | 1873.86 | 43,656 |
| 2 | `polint_` | 11.33% | 434.91 | 2,939,811,944 |
| 3 | `lopsided_` | 8.04% | 308.75 | 1,047,744 |
| 4 | `fdderivs_` | 6.93% | 266.12 | 480,216 |
| 5 | `kodis_` | 5.59% | 214.61 | 1,047,744 |
| 6 | `prolong3_` | 4.36% | 167.23 | 1,892,582 |
| 7 | `fderivs_` | 4.32% | 165.79 | 598,928 |
| 8 | `enforce_ga_` | 2.09% | 80.28 | - |
| 9 | `symmetry_bd_` | 1.75% | 67.32 | 6,256,902 |
| 10 | `rungekutta4_rout_` | 1.74% | 66.75 | - |
| 11 | `global_interp_` | 0.96% | 36.70 | - |
| 12 | `restrict3_` | 0.76% | 28.99 | - |

### 6.3 和 baseline 对比

| 函数 | baseline 占比 | V3 占比 | 观察 |
| --- | ---: | ---: | --- |
| `compute_rhs_bssn_` | 39.66% | 48.82% | 相对占比提高 |
| `polint_` | 23.67% | 11.33% | 占比显著下降 |
| `lopsided_` | 6.89% | 8.04% | 成为第三大热点 |
| `fdderivs_` | 5.81% | 6.93% | 相对占比提高 |
| `prolong3_` | 4.95% | 4.36% | 相对占比略降 |
| `kodis_` | 4.52% | 5.59% | 相对占比提高 |
| `fderivs_` | 4.24% | 4.32% | 占比接近 |

`polint_` 从 23.67% 降到 11.33%，表明 V1/V2 的 fixed-order interpolation 优化有效。两次 profiling 的配置和 MPI 数据合并方式可能不同，因此不简单比较 self seconds，也不从 self seconds 推导绝对加速。

### 6.4 热点转移分析

`compute_rhs_bssn_` 从 39.66% 上升到 48.82% 不表示它已经被证明绝对变慢。更合理的解释是：外围热点下降后，BSSN RHS 在总 CPU 时间中的相对占比提高，成为更集中的核心瓶颈。

新的主要外围热点为 `lopsided_`、`fdderivs_` 和 `kodis_`，占比分别为 8.04%、6.93% 和 5.59%，合计 20.56%。此外，`polint_` 仍有约 29.4 亿次调用，即使单次成本下降，累计开销仍值得继续调查。

## 7. 第四刀：Lopsided Interior Fast Path

### 7.1 思路

提交 `7c2a357 Optimize lopsided interior finite differences` 针对 `lopsided_` 的 `ghost_width=3` 路径，将规则内部区域与通用边界区域分开。

正向 stencil：

```text
S/(12h) * [-3f(i-1) - 10f(i) + 18f(i+1) - 6f(i+2) + f(i+3)]
```

负向 stencil：

```text
-S/(12h) * [-3f(i+1) - 10f(i) + 18f(i-1) - 6f(i-2) + f(i-3)]
```

规则内部区域为：

```text
i = 4 : ex(1)-3
j = 4 : ex(2)-3
k = 4 : ex(3)-3
```

内部区域先走无 stencil 边界判断的 fast loop；原通用 boundary loop、`symmetry_bd` 和全部 fallback 保留，内部点在通用 loop 中 `cycle`。

### 7.2 数值正确性

离线测试使用 `11 x 10 x 9` 随机数组，覆盖：

- `Symmetry=0/1/2`；
- 7 种实际 parity 组合；
- 三个方向全部 `-1/0/+1` 模式；
- 共 567 组对比。

结果：

| Symmetry | max_abs_error |
| ---: | ---: |
| 0 | 0 |
| 1 | 0 |
| 2 | 0 |

### 7.3 两次 benchmark

| 指标 | lopsided Run #1 | lopsided Run #2 | 平均 |
| --- | ---: | ---: | ---: |
| Total Evolve | 829.889 s | 823.473 s | 826.681 s |
| Total Running | 833.032 s | 826.760 s | 829.896 s |
| This run | 833.091 s | 826.834 s | 829.963 s |
| real | 54m26.619s | 54m19.326s | 约 54m22.973s |
| user | 109m46.998s | 108m51.191s | 约 109m19.095s |
| sys | 1m17.918s | 1m22.763s | 约 1m20.341s |

### 7.4 同期 shift3 A/B

为避免直接使用不同阶段的旧数据，测试从 `7316333` 建立独立 worktree，在同期环境重新执行 `make clean`、`make -j8 ABE` 和相同 20-step workload。

| 指标 | shift3 同期复测 | lopsided 两次平均 | lopsided 相对变化 |
| --- | ---: | ---: | ---: |
| Total Running | 805.984 s | 829.896 s | +23.912 s，约 +2.97% |
| real | 54m39.262s | 约 54m22.973s | -16.289 s，约 -0.50% |

### 7.5 为什么最终回滚

第四刀的 wall-clock 平均值只快约 0.50%，该差异较小，可能落在当前 WSL/MPI 环境波动范围内；与此同时，程序内部 `Total Running` 稳定退化约 2.97%。因此没有足够证据将它认定为有效优化。

该实验通过提交 `24f03a4 Revert "Optimize lopsided interior finite differences"` 正式回滚。第四刀只作为通过正确性验证但未通过性能验收的实验记录，不属于 V3 最终有效优化。

## 8. V3 最终版本

V3 保留完整实验历史，不 reset、不 squash。最终源码状态为：

- 保留 `7316333` 引入的 `fderivs_shift3`；
- metric6 已由 `3018a69` 回滚，正式源码不存在 `fderivs_metric6`；
- Gamma3 已由 `b9e0edb` 回滚，正式源码不存在 `fderivs_gamma3`；
- lopsided interior fast path 已由 `24f03a4` 回滚，`lopsidediff.f90` 恢复到第四刀之前的实现。

最终 V3 唯一新增并保留的核心优化是 shift derivative fusion。

## 9. Baseline / V1 / V2 / V3 性能总览

| 版本 | Total Evolve | Total Running | real | 说明 |
| --- | ---: | ---: | ---: | --- |
| 官方 Baseline | 3852.320 s | 4009.230 s | 66m49.233s | 原始正式数据 |
| V1 | 821.391 s | 824.376 s | 57m35.947s | V1 正式数据 |
| V2 profiling Run #1 | 892.440 s | 895.866 s | 55m28.511s | 带 MPI profiling instrumentation |
| V2 profiling Run #2 | 861.726 s | 864.966 s | 56m07.120s | 带 MPI profiling instrumentation |
| V2 profiling average | 877.083 s | 880.416 s | 55m47.816s | 两次 profiling 平均 |
| V3 shift3 Run #1 | 756.876 s | 759.967 s | 55m48.388s | clean V3 正式测试 |
| V3 shift3 同期复测 | 803.028 s | 805.984 s | 54m39.262s | 第四刀同期 A/B 基准 |

V2 数据包含 MPI profiling instrumentation，而 clean V3 已移除这些插桩，因此两者的 CPU 内部计时不能作为严格 apples-to-apples 对比。两个 V3 shift3 run 均列入报告，以避免单次运行波动造成选择性结论。

从 wall-clock 看，V3 两次 shift3 结果均显著低于官方 baseline，但相对 V1/V2 的小幅差异仍需结合测试条件谨慎解释。

### Wall-clock 波动说明

当前 WSL + OpenMPI 环境中，程序内部 `Total Running Time` 与 shell `time` 的 `real` 存在明显差异。例如内部计时约为 760 至 830 秒，而 wall-clock 约为 54 至 57 分钟。

现有数据不足以确定差异的确切来源。wall-clock 可能受到程序内部计时未覆盖阶段、MPI runtime、WSL 调度和系统负载等因素影响，需要进一步的专门 instrumentation 才能区分各部分贡献。不能将差异简单断言为全部来自 MPI 等待、I/O 或 WSL 调度。

因此代码微优化需要同时参考：

- `Total Evolve`；
- `Total Running`；
- `user`；
- `real`。

比赛最终仍关注 `real`，但小于约 1% 的单次 wall-clock 差异不能直接认定为有效优化。

## 10. V3 的经验与下一阶段方向

V3 的价值不仅是成功保留了 shift derivative fusion，也通过三个负收益实验确定了多项实现边界：

1. loop nest 数量减少不等于整体运行一定更快。
2. 扩大 derivative fusion 时，寄存器压力、cache 行为、memory traffic 和 kernel size 可能抵消循环融合收益。
3. 数值正确性是必要条件，但仍必须由同期 A/B benchmark 决定性能改动是否保留。
4. wall-clock 小幅变化需要重复测试，并与程序内部计时共同判断。
5. 失败实验应保留 commit 和 revert 历史，避免后续重复尝试同一路径。

下一阶段优先研究：

1. `compute_rhs_bssn` 内部计算结构；
2. `fdderivs`；
3. `kodis`；
4. `polint` 仍高达约 29 亿次调用，继续调查固定阶二维/三维插值链。

这些方向仅作为后续工作建议，V3 收尾不再修改相应源码。
