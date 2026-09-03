# AMSS-NCKU V5 优化记录

## 1. V4 状态

V4 是进入当前 V5 阶段前的稳定源码优化版本。已有主要源码优化包括：

- 编译器优化：`-O3`、`-march=native`、`-flto` 以及已有的 Fortran inline 配置；
- `polint` 六点插值专用快速路径；
- `fderivs` 内部区域直接 stencil 优化；
- `shift3` 中 `betax`、`betay`、`betaz` 导数融合；
- `kodis` 边界范围优化，去除热点循环中的逐点边界判断并改善向量化。

此前 V4 的历史测试曾得到：

- Total Evolve Time ≈ 774.582 s；
- 单步 ≈ 38.729 s/step。

但后续在同一台机器和当前 WSL 环境重新测试时，性能数字发生了明显波动。因此，该 V4 数字只作为历史记录；后续 V5 对比优先使用近期同环境测试数据。

## 2. Baseline 重新测试与异常变化说明

以前保存的 baseline 历史测试为：

- Total Evolve Time = 823.306 s；
- 单步时间 = 41.1653 s/step。

近期重新运行保存的 clean baseline，使用 4 MPI rank、20 steps，得到：

- Total Evolve Time = 886.447 s；
- Total Running Time = 889.867 s；
- shell real = 953.19 s；
- 单步时间 = 44.32235 s/step。

baseline 相比早期记录明显变慢。目前无法确认造成 baseline 变化的唯一原因。

**该变化目前原因未知，可能与 WSL 调度、Windows 主机负载、CPU 功耗/频率、后台进程、温度、MPI/hwloc 运行环境等因素有关，但目前没有证据证明具体是哪一个原因，因此不能把其中任何一种可能性写成已确认原因。**

后续当前阶段统一使用近期重新测试得到的 **44.322 s/step** 作为当前机器上的 baseline 参考值。

## 3. MPI rank 数量测试

以下均为近期同环境测试结果。

### 当前源码，4 MPI rank

- Total Evolve Time = 785.343 s；
- Total Running Time = 788.746 s；
- shell real = 834.87 s；
- 单步：785.343 / 20 = 39.26715 s/step；
- 相对当前 4-rank baseline：44.32235 / 39.26715 ≈ 1.129x，即单步时间约减少 11.4%。

### 当前源码，8 MPI rank，不显式绑核

- Total Evolve Time = 757.344 s；
- Total Running Time = 760.150 s；
- shell real = 800.07 s；
- 单步时间 = 37.8672 s/step。

### 12 MPI rank

两次测试的 Total Evolve Time 分别为 986.879 s 和 1046.220 s，对应约 49.344 s/step 和 52.311 s/step，明显比 8 rank 慢。这说明增加 MPI rank 并不会持续提高性能。

### 10 MPI rank + 当前 hwthread 绑定

- Total Evolve Time = 974.349 s；
- Total Running Time = 977.103 s；
- 单步时间 = 48.71745 s/step。

该结果明显劣于 8 rank。因此，在当前测试范围内，**8 MPI rank 是最佳选择**。

## 4. Process Pinning / MPI 绑核优化

这是 V5 阶段根据 Zero Point ASC26 优化思路进行的系统级优化。

### 8 rank 默认方式

- Total Evolve Time = 757.344 s；
- 单步时间 = 37.8672 s/step。

### 8 rank，绑定整个 core

命令：

```bash
HWLOC_COMPONENTS=-gl time mpirun -np 8 \
  --bind-to core \
  --map-by core \
  ./ABE
```

结果：

- Total Evolve Time = 739.600 s；
- Total Running Time = 743.048 s；
- shell real = 789.21 s；
- 单步时间 = 36.9800 s/step；
- 相比未绑核约提升 2.34%。

### 8 rank，每 rank 固定一个不同 core 的单个 hardware thread

命令：

```bash
HWLOC_COMPONENTS=-gl time mpirun -np 8 \
  --use-hwthread-cpus \
  --map-by ppr:1:core \
  --bind-to hwthread \
  ./ABE
```

OpenMPI 映射确认：rank 0 到 rank 7 分别映射至 core 0 到 core 7，并且各自使用 hwt 0。

20-step 结果：

- Total Evolve Time = 733.762 s；
- Total Running Time = 736.838 s；
- shell real = 799.48 s；
- 单步：733.762 / 20 = 36.6881 s/step。

这是目前 V5 阶段的最好成绩。相对当前 baseline：44.32235 / 36.6881 ≈ 1.208x，单步时间下降约 17.2%。因此，当前推荐运行配置就是上述单 hwthread pinning 命令。

## 5. WSL / hwloc 环境问题

在当前 WSL 环境下，OpenMPI/hwloc 默认可能因为 OpenGL/NVIDIA display probing 在程序启动前卡住。目前使用临时环境变量：

```bash
HWLOC_COMPONENTS=-gl
```

绕过该问题。不修改系统配置，也不将其写入 shell profile。

`lscpu` 在 WSL 中显示 24 logical CPUs、12 cores、2 threads/core；`/sys/devices/system/cpu/cpu*/cpu_capacity` 全部为 1024。因此，WSL 没有可靠暴露 Intel i7-14650HX 的 P-core/E-core 差异，当前不能可靠判断具体哪些 Linux CPU ID 对应真实 P-core/E-core。当前绑核策略只根据 WSL 暴露的 core/hardware thread 拓扑设计。

## 6. Thread Control 检查

检查结果如下：

```text
OMP_NUM_THREADS=
OPENBLAS_NUM_THREADS=
MKL_NUM_THREADS=
NUMEXPR_NUM_THREADS=
```

执行以下动态库检查没有输出：

```bash
ldd ./ABE | egrep 'gomp|omp|openblas|mkl|blas'
```

因此，当前 ABE 没有发现明显的 OpenMP/OpenBLAS/MKL/BLAS 隐式多线程问题，这一方向暂时没有继续优化。

## 7. 当前 V5 阶段结论

| 配置 | 单步时间 | 相对当前 baseline 加速比 |
| --- | ---: | ---: |
| Baseline，4 rank | 44.322 s/step | 1.000x |
| 当前源码，4 rank | 39.267 s/step | 约 1.129x |
| 当前源码，8 rank | 37.867 s/step | 约 1.170x |
| 当前源码，8 rank + bind-to-core | 36.980 s/step | 约 1.199x |
| 当前源码，8 rank + 单 hwthread pinning | 36.688 s/step | 约 1.208x |

训练目标约为 baseline 的 2x，对应目标单步约为：44.322 / 2 = 22.161 s/step。因此当前仍未达到 2x。

当前 V5 的系统级 MPI 调优已经基本确定：

- 8 rank 优于 4 / 10 / 12 rank；
- process pinning 有小幅正收益；
- 每 rank 固定到不同 core 的单一 hardware thread 是目前最好配置；
- 仅靠 MPI 参数优化不足以达到 2x；
- 后续需要继续进行源码级、通信级或 GPU 架构优化。


## 8. V4 → 当前 V5 源码状态审查

以 `optimize-v4` 与当前 `optimize-v5` 的已提交历史进行比较，V5 阶段进行过若干结构性源码实验，包括 metric temporary、lopsided SIMD、RK4 traversal、prolong3 paired-x reuse、MPI Sync batching 和 conformal metric RHS traversal。上述实验均已在提交历史中回滚；当前从 V4 到 V5 的已提交净差异只有相关审计和回滚说明文档，没有保留下来的净源码差异。

工作区仍存在输入参数、OpenMP、GPU、构建与诊断性质的未提交改动，以及测试输出、备份、二进制文件和日志。由于这些内容没有被当前稳定 V5 提交历史确认，部分 GPU 实验也尚未完成运行正确性验证，因此不纳入本次 CPU V5 稳定提交。
