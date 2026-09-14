# AMSS-NCKU 最终优化总结

## 1. 测试口径与结果说明

最终测试设备为 Intel Core i7-14650HX，运行环境为 WSL2，正式对比使用 8 个 MPI 进程和 20 个 evolution steps。实际测试过程中会受到 CPU frequency、温度、Windows/WSL 调度、后台负载及 MPI runtime 波动影响，因此不同时期的绝对运行时间存在明显变化。

> 在最终数据整理阶段，考虑到设备运行状态存在波动，重新进行了 baseline 测试，并将 1088 s 作为最终统一对比基准。

本文区分整程序的 `Total Evolve Time` 与 profiling 中单个热点的累计时间，二者不能混用。所有加速比均以 1088 s baseline 计算。

## 2. V4 以及之前：建立编译与核心计算基础

早期版本的主要问题是通用编译配置没有充分利用本机指令集，固定规模插值仍承担通用控制开销，热点导数计算也存在重复遍历和 helper 开销。

最终保留的工作包括：

- 编译选项使用 `-O3`、`-march=native` 和 `-flto`，Fortran 启用 `-finline-functions`。这些设置让编译器使用当前 CPU 指令集，并通过 LTO 改善跨文件优化和内联。
- 为高频固定六点插值增加 `polint6_fast`，减少通用 `polint` 路径的循环和控制开销。
- 优化 `fderivs` 的固定 stencil 热路径，将规则内部区域与边界处理分开，降低逐点边界判断成本。
- 对 shift 一阶导数使用专用处理，使三个 shift 分量共享网格遍历控制，同时保持各分量 stencil 独立。
- 保留经过实测有效的 KO dissipation（`kodis`）规则路径和 boundary hot-check 调整；后期失败的 reciprocal 替换不在最终代码中。

这些修改形成了后续 profiling 和局部优化的稳定编译、插值与导数基础。V4 以及之前相对最终 baseline 的统一加速比为 1.060×。

## 3. V5：系统排除高风险方向

## V5：MPI 并行配置与绑核优化

V5 阶段的重点从单纯修改计算内核，逐渐转向对程序整体并行运行方式进行优化。

前期测试中发现，AMSS-NCKU 的运行性能不仅取决于单个计算 kernel 的效率，还明显受到 MPI 进程数量、CPU 核心分配方式以及进程调度稳定性的影响。由于程序中包含大量网格计算、插值和 MPI 通信，如果 MPI rank 数量设置不合理，或者进程在不同 CPU 核心之间频繁迁移，会导致计算资源利用率下降，同时增加缓存失效和调度开销。

因此，V5 阶段对 MPI 运行配置进行了重新测试和调整，最终将主要运行方式优化为：

```bash
HWLOC_COMPONENTS=-gl mpirun -np 8 \
  --use-hwthread-cpus \
  --map-by ppr:1:core \
  --bind-to hwthread \
  ./ABE

V5 相对最终 baseline 的统一加速比为 1.208×。

## 4. V6：MassPAng 插值与通信优化

### 4.1 分层 profiling

分层 profiling 将主要热点定位到：

`AnalysisStuff → surf_MassPAng → Interp_Points → global_interp → polin3 → polint`

这条路径对多个 field 重复处理相同的插值几何、owner、索引、reflection 和权重。

### 4.2 sparse-owner

V6-4 引入 sparse-owner 通信路径，使用 `MPI_Alltoallv` 只交换实际 owner 所需的数据，减少无效数据交换和 owner 信息处理。

### 4.3 global_interp_batch17

`global_interp_batch17` 将原来 17 个 field 的独立插值调用组成 batch。17 个 field 共享：

- geometry 和 interpolation location；
- source mapping 与 reflection；
- interpolation weights；
- owner/control information。

每个 field 的数值 gather 与 contraction 仍保持相对独立。成功的关键是共享真正相同的几何和控制工作，而非强行融合所有数值表达式。历史 profiling 中 `surf_MassPAng` 热点时间约从 302.6 s 降至 22.9 s；这是 profiling 层面的热点累计时间，不是整程序 `Total Evolve Time`。

V6 相对最终 baseline 的统一加速比为 1.809×。

## 5. V7 最终保留的优化

### 5.1 V7-1：shift 二阶导数 producer-consumer fusion

原路径会 materialize 18 张 shift 二阶导数 full-grid array。`fdderivs_shift_fusion` 仍计算所需的 18 个点值导数，但只输出后续实际消费的六张 contraction：三个 `grad(div beta)` 分量和三个 shift Laplacian 分量。该修改减少了 full-grid temporary、内存写入及随后的数组读取。

### 5.2 V7-4：chi-Ricci pointwise fusion

原 chi-Ricci correction 由 scalar temporary pass 和六个 Ricci component pass 组成。V7-4 在一个局部三维 loop 中计算 point-local inverse-chi、Hessian/gradient 组合和 scalar contraction，并立即更新六个 Ricci 分量。它没有新增 full-grid temporary，保持原数学表达式，同时减少数组扫描和中间数组流量。

该阶段历史最好结果为 `Total Evolve Time = 566.725 s`，20 steps 平均 `28.336 s/step`。

### 5.3 V7-13-1：generic fderivs output initialization

原 generic `fderivs` 每次完整清零 `fx`、`fy`、`fz`，而 stencil 随后会覆盖内部绝大部分区域。V7-13-1 将 O(N³) full-grid zero 改为只初始化 stencil 不会覆盖的 surface/boundary 区域，减少不必要的内存写入。

### 5.4 V7-13-2：symmetry_bd_deriv_fast

`Cell`、`ghost_width = 3`、`ord = 2` 的 derivative symmetry preparation 中，internal copy 及 x/y/z reflection 会完整覆盖 extended buffer。`symmetry_bd_deriv_fast` 利用该覆盖关系，避免原先的 expanded-buffer memset，并用于：

- generic `fderivs`；
- `fderivs_shift3`；
- `fdderivs_shift_fusion`。

该路径保留原 `symmetry_bd` 作为不满足条件时的 fallback，并经过完整 extended-buffer 和 derivative output correctness 检查。

## 6. 失败实验带来的经验

后期未保留的实验包括 Ricci grouped fusion、Aij producer-consumer fusion、shift `fd2` inline、`kodis` reciprocal hoist、lopsided split、Wave batch 和 generic `fdderivs` nozero。它们虽然可能减少 pass、branch、temporary、division 或函数调用，实际机器码却可能增加 register pressure、spill/reload、代码尺寸，破坏 SIMD 或改变 LTO code layout，最终实测变慢。设备频率和 WSL 调度波动也要求所有结论以可重复的完整运行数据为准。

项目后期已经进入明显的性能瓶颈区间，继续进行局部源码微优化的风险明显增大。

## 7. 最终性能对比

统一测试口径为 20 steps：

| 阶段 | Total Evolve Time / s | 单步时间 / s | 相对 Baseline 加速比 |
|---|---:|---:|---:|
| Baseline | 1088.000 | 54.400 | 1.000× |
| V4及之前 | 1026.415 | 51.321 | 1.060× |
| V5 | 900.662 | 45.033 | 1.208× |
| V6 | 601.437 | 30.072 | 1.809× |
| V7 最终最好结果 | 566.725 | 28.336 | 1.920× |

最终最好 `Total Evolve Time` 为 566.725 s，平均 28.336 s/step。相对重新测试的 1088 s baseline，加速比为 `1088 / 566.725 = 1.9198×`，约 1.920×。

## 8. 最终结果与停止优化的原因

以 1088 s baseline 计算，2× 目标为 544 s。当前最好结果为 566.725 s，距离目标还差 22.725 s，相对当前最好仍需要约 4.18% 的额外性能提升。

经过 V7 后期实验，大规模 fusion 经常增加 register pressure；pass 或函数调用减少并不必然加速；reciprocal multiply 替换 vector divide 也出现过负优化；LTO 和 AVX2 已对多个规则 kernel 做了较充分优化。剩余热点更可能需要数据布局或算法层面的较大重构，继续做高风险局部微优化不适合作为本轮收尾工作。

> 不再继续进行高风险的局部微优化，将当前最好的 566.725 s、约 1.92× 加速结果作为最终提交结果。

最终结果没有完全达到 2×，但已经接近目标。整个过程完成了编译器优化、固定插值专用化、MPI 通信优化、内存访问优化、temporary reduction、producer-consumer fusion，以及 SIMD 和机器码分析的一套完整工程实践。

## 9. 最终代码状态

最终生产代码保留 V1–V6 中确认有效的优化，以及 V7-1、V7-4、V7-13-1、V7-13-2。关键路径包括：

- `global_interp_batch17` 和 sparse-owner interpolation；
- `fdderivs_shift_fusion`；
- chi-Ricci pointwise fusion；
- generic `fderivs` surface initialization；
- `symmetry_bd_deriv_fast`。

V7-15-2、V7-17-1、V7-18-1 及其他已判定 REVERT 的实验不在最终生产路径中。Profiling 和 correctness 开关默认关闭。
