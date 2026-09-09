# AMSS-NCKU V6-6-1 优化记录

## 1. 优化背景

逐层 MPI profiling 将 CPU 路径中的异常耗时定位到：

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

`surf_MassPAng` 是 `AnalysisStuff` 中最突出的热点，其中主要成本又集中在曲面采样点插值。早期分析将这部分拆成两类开销：一是 `global_interp/polin3` 的实际插值计算，二是汇总完整 `shellf` 的 MPI 通信。V6-6-1 针对前者中重复执行的工作进行优化。

## 2. 前期优化基础

V6-4-1 已经把 `surf_MassPAng` 从 dense collective 改成专用的 sparse owner-computes 数据流：

```text
consumer 的本地点段
  -> 判断每个点的 owner
  -> 通过 Alltoallv 把点索引发给 owner
  -> owner 完成插值
  -> 通过 Alltoallv 返回 17 个变量的结果
  -> consumer 按原顺序积分
```

这避免了对完整 `NN * 17` 结果数组进行 dense collective，同时保留了 generic `Patch::Interp_Points`，因此不会改变 `Prolongint`、`surf_Wave` 等其他调用路径。不过，V6-4-1 的 owner 端仍然要对每个曲面点调用 17 次通用插值函数，插值计算本身依旧很重。

## 3. 问题分析

同一个 surface point 需要插值以下 17 个变量：

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

对于同一个点，这 17 个变量共享相同的：

- `xyz` 坐标；
- owner Block；
- Cell-centered 几何信息；
- `6 * 6 * 6` 插值 stencil；
- source index 映射；
- reflection 几何；
- 三个方向的六点插值权重。

V6-4-1 中仍按变量逐次执行：

```text
for each point:
  for each of 17 variables:
    f_global_interp
      -> global_interp
        -> 重算几何与 stencil 映射
        -> 重做 reflection 判断
        -> decide3d
        -> polin3
          -> 43 次 polint
```

也就是说，必要的 17 组变量数据读取之外，同一点不随变量变化的准备工作也重复了 17 次。这正是 V6-6-1 要消除的部分。

## 4. V6-5-1 失败尝试

V6-5-1 曾尝试缓存插值索引与系数，再通过 `f_global_interpind` 计算：

```text
cached inds + coef
  -> f_global_interpind
```

性能测试出现明显退化，因此该方案已经回退。源码分析表明，`global_interpind` 并非简单执行 216 项加权求和：它仍会构造 `ya(6,6,6)`；在 boundary/symmetry 情况下，还可能通过 `symmetry_tbd` 生成并初始化完整 Block 大小的 `fh`，带来远超局部 stencil 的内存访问和复制成本。

因此 V6-6-1 没有继续做 persistent cache，也没有沿用 `global_interpind`，而是保留 `decide3d` 的局部 `6 * 6 * 6` symmetry 语义，只合并同一点的 17 个变量。

## 5. V6-6-1 核心优化

V6-6-1 新增 Cell-only 的 `global_interp_batch17` kernel。原来的结构是：

```text
for each point:
  for each of 17 variables:
    geometry
    decide3d
    polin3
```

新的结构是：

```text
for each point:
  geometry              * 1
  source mapping        * 1
  reflection mapping    * 1
  wx / wy / wz          * 1

  for each of 17 variables:
    按该变量自己的 SoA(1:3) 处理 reflection sign
    gather 216 values
    z -> y -> x tensor reduction
```

每个点只计算一次 Cell-centered 坐标、clamp、归一化坐标、`src_i/src_j/src_k`、reflection mask 以及 `wx/wy/wz`。之后 17 个变量分别读取自己的 216 个 double，并使用各自的 symmetry parity；不同变量的最终符号没有被错误共享，共享的只有 reflection 几何。

该实现的边界很明确：

- 不调用 `global_interpind`；
- 不构造完整 Block 的 `fh`；
- 不调用 generic `polin3/polint`；
- 不维护跨点或跨 timestep 的 persistent cache；
- 不减少每个变量必须读取的 216 个数值，只消除同一点重复的几何、索引、reflection 和权重准备；
- generic `Patch::Interp_Points`、`global_interp` 和 `polin3` API 保持不变。

专用 sparse-owner 路径仅在 Cell build、`num_var == 17`、`ORDN == 6` 且 batch routine 成功时使用 batch17，否则回退到原来的 17 次 `f_global_interp`。

## 6. Profiling 结果

以下数据来自开启 profiling 的构建，只用于定位和比较热点，不与最终 clean performance 数据混用：

| 指标（max rank） | V6-4-1 profiling | V6-6-1 profiling | 变化 |
|---|---:|---:|---:|
| Total Evolve Time | 924.599 s | 607.589 s | 整体约 1.52 倍 |
| `surf_masspang` | 302.606 s | 22.926 s | 约 13.2 倍 |
| `interp_points` | 301.226 s | 21.630 s | 约 13.9 倍 |
| `interp_global_interp` | 300.324 s | 20.443 s | 显著下降 |

其中：

```text
surf_MassPAng: 302.606 s -> 22.926 s，约 13.2x
Interp_Points: 301.226 s -> 21.630 s，约 13.9x
```

V6-6-1 的 surface 优化路径不再进入 generic `f_global_interp -> global_interp -> polin3` 调用链，因此该路径的 Level-6、Level-7 样本为空。这是调用路径被绕开后的预期现象，不代表计时数据丢失。

## 7. Profiling 清理

Level-1 到 Level-7 profiling 源码仍然保留，以便后续重新启用，但统一受编译期开关控制：

```c
#ifndef ENABLE_ABE_PROFILING
#define ENABLE_ABE_PROFILING 0
#endif
```

默认 clean performance 构建同时设置：

```text
ENABLE_ABE_PROFILING=0
VERIFY_BATCH17=0
```

因此默认构建不会执行：

- profiling 用的 `MPI_Wtime`；
- profiling counters 与 min/max/count bookkeeping；
- profiling 汇总所需的 `MPI_Reduce/MPI_Gather`；
- Level-1 到 Level-7 输出；
- batch17 与 generic 路径的 correctness 双路径计算。

这些开关只裁剪诊断代码，不改变 V6-6-1 的数学实现。此前 clean `make clean && make -j2 ABE` 已成功，生成的默认 ABE 中也未检出 Level-1 到 Level-7 的 profile 输出字符串。

## 8. Baseline 环境变化说明

历史 baseline 曾测得：

```text
Total Evolve Time = 846.867 s
单步 = 42.34335 s/step
```

之后在当前设备和运行环境中重新 clean 测量 baseline，结果变为：

```text
Total Evolve Time = 1004.78 s
Total Running Time = 1008.06 s
```

目前无法明确判断这一变化的具体原因。可能相关的因素包括 WSL 状态、Windows host 调度、CPU 频率/功耗/温度、hybrid CPU scheduling、MPI placement、后台负载，以及 memory/cache/system state；这些都只是可能性，没有经过严格验证。

因此最终对比不使用历史的 846.867 s，而采用与 V6-6-1 最终测试处于当前同一环境下重新测得的 clean baseline 1004.78 s。这样能减少环境差异造成的不公平比较，但不声称两次运行环境得到了完全控制。

## 9. 最终 Clean 性能结果

最终比较的两组数据均为 clean 20-step 测量：

| 指标 | Baseline | V6-6-1 |
|---|---:|---:|
| Total Evolve Time | 1004.78 s | 586.906 s |
| Total Running Time | 1008.06 s | 589.825 s |
| steps | 20 | 20 |
| Evolve 单步耗时 | 50.239 s/step | 29.3453 s/step |

计算如下：

```text
Baseline 单步 = 1004.78 / 20 = 50.239 s/step
V6-6-1 单步 = 586.906 / 20 = 29.3453 s/step

Evolve 加速比 = 1004.78 / 586.906 = 约 1.712x
Evolve 耗时下降 = (1004.78 - 586.906) / 1004.78 = 约 41.59%
Total Running Time 加速比 = 1008.06 / 589.825 = 约 1.709x
```

V6-6-1 clean 第一次运行的 Evolve 时间为 681.167 s，说明当前环境仍存在运行波动；该数据只作为波动证据保留，不作为最终结果。

## 10. Correctness 验证

correctness 检查在真实生产 Block 和真实 surface point 上直接比较：

```text
global_interp_batch17
vs
原生产 f_global_interp * 17
```

当前 `Symmetry=1` 测试覆盖 23 个点、391 个 variable samples，reflection masks 为 `[0, 4]`，即：

- 无 reflection；
- z reflection。

结果如下：

| 指标 | 数值 |
|---|---:|
| max_abs | 6.93889390390722838e-18 |
| mean_abs | 1.33154764489875235e-19 |
| max_rel | 4.43827200454890002e-15 |
| mean_rel | 1.87151296325769851e-16 |

这些差异属于双精度浮点运算的舍入误差量级。本次验证没有覆盖 x reflection、y reflection 或多轴 reflection，因此不对这些情况作 correctness 声明；如需验证，必须使用能够实际覆盖相应路径的 `Symmetry=2` 配置。

## 11. 当前结果与 2x 目标

以当前 clean baseline 的 50.239 s/step 为基准，2 倍加速目标为：

```text
50.239 / 2 = 25.1195 s/step
```

V6-6-1 当前为 29.3453 s/step，尚未达到 2x。相对当前版本还需要：

```text
29.3453 / 25.1195 = 约 1.168x
```

的额外加速。

## 12. 总结

本轮工作通过逐层 profiling 从 `AnalysisStuff` 定位到 `surf_MassPAng`，最终确认同一个曲面点对 17 个变量重复执行几何、stencil、reflection 和权重准备是主要问题。Cell-only `global_interp_batch17` 将这些不随变量变化的工作从 17 次降为 1 次，同时保留每个变量独立的数据读取和 symmetry sign 语义。

结果上，profiling 构建中的 `surf_MassPAng` 从 302.606 s 降到 22.926 s，约加速 13.2 倍；当前环境下的 clean 整体 Evolve 从 1004.78 s 降到 586.906 s，约加速 1.712 倍，耗时下降约 41.59%。当前 `Symmetry=1` 覆盖范围内，batch17 与原生产路径的差异处于双精度舍入误差量级。

目前整体性能尚未达到 2 倍目标。由于 `surf_MassPAng` 热点已经显著缩小，后续更合理的方向是重新 profiling 并转向 `compute_rhs_bssn`、`surf_Wave` 等剩余热点，而不是继续集中优化已经降到较低占比的 `surf_MassPAng`。
