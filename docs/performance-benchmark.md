# 性能基准测试

性能基准程序用于测量 v1.0 稀疏求解链在不同模型规模下的组装、静力和 Newmark 瞬态计算表现。

Benchmark 使用轴向 Beam3D 链式模型，主要用于观察：

- 稀疏矩阵非零元增长
- 稀疏存储优势
- 组装耗时
- 静力求解耗时
- Newmark 求解耗时
- 模型规模扩展趋势

它不是针对某台机器给出固定实时性能承诺。

## 构建与运行

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel

./build/solver_benchmark \
  --sizes 100,500,2000 \
  --steps 50 \
  --repeats 5 \
  --csv benchmark.csv
```

## CSV 输出字段

```text
elements
global_dofs
free_dofs
k_nnz
m_nnz
k_fill_ratio
sparse_k_bytes
dense_k_bytes
dense_to_sparse_k_memory_ratio
assembly_ms
static_total_ms
newmark_total_ms
newmark_steps
static_tip_relative_error
```

含义如下：

- `elements`：梁单元数量
- `global_dofs`：全局自由度数量
- `free_dofs`：自由自由度数量
- `k_nnz`：总体刚度矩阵 K 的非零元数量
- `m_nnz`：总体质量矩阵 M 的非零元数量
- `k_fill_ratio`：K 的非零元占完整稠密矩阵元素总数的比例
- `sparse_k_bytes`：K 的近似稀疏存储量
- `dense_k_bytes`：同尺寸稠密 K 所需存储量
- `dense_to_sparse_k_memory_ratio`：稠密/稀疏存储量比值
- `assembly_ms`：组装耗时
- `static_total_ms`：静力全过程耗时
- `newmark_total_ms`：Newmark 全过程耗时
- `newmark_steps`：Newmark 时间步数
- `static_tip_relative_error`：静力端部位移相对解析误差

## 时间定义

### assembly_ms

`Model::assemble()` 的中位数运行时间。

包括总体：

- K
- M
- F

的组装。

### static_total_ms

`LinearStaticSolver::solve()` 的端到端中位数时间，包括：

- 总体组装
- 自由 DOF 提取
- 稀疏约化
- 稀疏矩阵分解
- 方程求解
- 支座反力恢复

### newmark_total_ms

`NewmarkBetaSolver::solve()` 的端到端中位数时间，包括：

- 总体组装
- 稀疏约化 K/M
- Rayleigh 阻尼
- 质量矩阵分解
- 初始加速度
- 有效刚度构造
- `Keff` 一次性稀疏分解
- 全部时间步求解
- 位移/速度/加速度历史存储

## 为什么采用多次运行中位数

CI 运行环境会受到：

- CPU 调度
- 同机其他任务
- 虚拟化
- 缓存状态

影响。

因此 benchmark 不使用单次运行时间，而使用多次运行结果的中位数降低随机噪声。

## 稀疏存储估算

`sparse_k_bytes` 根据压缩稀疏矩阵中的：

- 数值
- 索引
- 外层指针

进行近似估算。

`dense_k_bytes` 按：

```text
global_dofs * global_dofs * sizeof(double)
```

计算。

这个比值表示 K 本身的理论存储差异，不代表整个进程的实际常驻内存。

## 当前 CI 基准规模

CI 自动运行：

```text
100 个梁单元
500 个梁单元
2000 个梁单元

50 个 Newmark 时间步
每个规模重复 5 次
```

并上传：

```text
benchmark.csv
```

作为 `solver-benchmark` artifact。

## 正确性保护

Benchmark 不只计时，还会将轴向梁链的静力端部位移与解析解比较：

```text
u_tip = P L_total / EA
```

如果相对误差超过基准测试允许范围，benchmark 会直接失败。

因此性能测试必须建立在结果正确的前提下。

## 如何理解 CI 时间

GitHub Hosted Runner 的绝对毫秒数不应视为固定性能承诺。

合理的比较方式是：

- 使用相同 Build Type
- 使用相同 runner 类型
- 使用相同模型
- 使用相同时间步
- 使用多次运行
- 观察长期趋势

## 模态分析说明

生产 `ModalSolver` 已采用稀疏 `Kff/Mff` 与 Shift-Invert Lanczos 求解低阶广义特征值。

当前 `solver_benchmark` 可执行程序仍主要测量：

- 稀疏总体组装
- 稀疏静力求解
- 稀疏 Newmark 瞬态求解

模态正确性和稀疏求解路径由 `v1_verification` 与 `sparse_modal` 测试覆盖；其中 `sparse_modal` 会将低阶 Lanczos 结果与小规模 dense 参考解交叉验证。
