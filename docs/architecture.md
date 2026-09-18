# 求解器架构

本项目采用分层结构，将模型数据、单元力学、自由度编号、总体组装、数值求解和结果后处理分离。

设计思想参考了 OpenSees 的 Domain/Node/Element 分层、deal.II 的独立 DOF 管理，以及 MFEM、CalculiX 中“单元局部计算 + 总体组装”的有限元组织方式，但未复制上述项目的实现代码。

## 总体计算流程

```text
Node + Material + Section
          │
          ▼
       Beam3D
   局部 Ke / Me
          │
          ▼
      坐标变换
          │
          ▼
   全局单元 Ke / Me
          │
          ▼
    Model::assemble
          │
          ▼
      DofManager
          │
          ▼
 Sparse Global K / M
          │
          ├────────────────────┐
          ▼                    ▼
LinearStaticSolver         ModalSolver
          │                    │
          │               Sparse Kff/Mff
          │                    │
          │             Shift-Invert Lanczos
          │
          ▼
NewmarkBetaSolver
          │
          ▼
   U / V / A / Reaction
          │
          ▼
Element / Section Response
          │
          ▼
       Recorder
          │
          ▼
       Envelope
```

## 核心模块

### 1. Model

`Model` 负责保存：

- 节点
- 单元
- 节点荷载
- 单元荷载
- 时变单元荷载

同时负责总体：

- 刚度矩阵 K
- 质量矩阵 M
- 荷载向量 F

的组装。

K/M 采用 `Eigen::SparseMatrix<double>`。

## 2. DofManager

`DofManager` 负责：

```text
(NodeId, Dof) -> global equation number
```

节点 ID 不要求连续，因此不能简单采用：

```text
equation = node_id * 6 + offset
```

独立 DOF 映射可以保证 Model、Solver、Recorder 使用同一套全局编号。

## 3. Element / Beam3D

`Element` 负责定义单元级接口。

`Beam3D` 是当前主要单元类型：

- 两节点
- 每节点 6 DOF
- 12×12 局部刚度矩阵
- 12×12 一致质量矩阵
- 轴向
- Saint-Venant 扭转
- 双向 Euler-Bernoulli 弯曲

局部坐标 x 轴由：

```text
node i -> node j
```

确定。

刚度矩阵通过：

```text
K_global = T^T K_local T
```

转换到全局坐标系。

## 4. ElementLoad

单元荷载独立于单元本体，用于描述：

- 均布荷载
- 线性分布荷载
- 局部分布荷载
- 任意位置集中力/力矩
- 时变单元荷载

连续荷载通过形函数转换为一致等效节点荷载。

时变荷载采用空间荷载 + 标量时间函数的组合：

```text
f(x,t) = scale(t) * f_spatial(x)
```

## 5. 稀疏总体组装

单元矩阵仍使用小型稠密矩阵。

总体组装采用：

```text
Element Ke / Me
     ↓
DOF mapping
     ↓
Eigen::Triplet
     ↓
setFromTriplets
     ↓
Sparse Global K / M
```

共享节点产生的重复矩阵项自动累加。

## 6. LinearStaticSolver

静力方程：

```text
K U = F
```

约束后求解：

```text
Kff Uf = Ff
```

`SparseDofReducer` 从总体稀疏矩阵中直接提取稀疏约化矩阵。

求解器使用：

```text
Eigen::SimplicialLDLT
```

支座反力通过：

```text
R = K U - F
```

恢复。

## 7. ModalSolver

模态分析求解：

```text
K phi = omega^2 M phi
```

总体和约化矩阵均保持稀疏：

```text
Sparse Global K/M
      ↓
SparseDofReducer
      ↓
Sparse Kff/Mff
      ↓
Spectra Shift-Invert
      ↓
Lanczos
```

求解器通过 Spectra 的对称广义特征值 Shift-Invert Lanczos 算法提取最低阶结构模态。

对于自由-自由模型，为避免 `K` 因刚体模态奇异而无法在零点做 shift-invert，程序采用一个很小的负 shift。由于非负特征值到该负 shift 的距离仍按特征值从小到大排列，因此首先得到刚体零频模态，然后得到最低阶正结构模态；程序会过滤零/非正特征值。

返回的振型按：

```text
phi^T M phi = 1
```

进行质量归一化，同时保留 Lanczos 迭代次数和矩阵算子调用次数用于诊断。

## 8. NewmarkBetaSolver

动力平衡方程：

```text
M a + C v + K u = F(t)
```

Rayleigh 阻尼：

```text
C = alpha M + beta K
```

Newmark 有效刚度：

```text
Keff = K + a0 M + a1 C
```

对于当前线性时不变系统，`Keff` 在整个时程中保持不变，因此：

1. 在时间循环前构造 `Keff`
2. 使用 `SimplicialLDLT` 分解一次
3. 每个时间步只更新有效荷载
4. 重复使用同一分解结果求解

## 9. 响应恢复

梁端力恢复采用：

```text
u_local = T u_global

f_local = K_local u_local - f_eq_local
```

对于梁任意位置的截面内力，不使用端力简单插值，而是根据真实荷载分布做截面平衡，得到：

```text
N(x)
Vy(x)
Vz(x)
T(x)
My(x)
Mz(x)
```

集中荷载位置使用 Left/Right 两侧极限保存内力跳变。

## 10. 截面应力

当前支持：

```text
sigma_x = N/A - Mz*y/Iz + My*z/Iy
```

并根据具体截面类型提供部分剪切/扭转应力模型。

## 11. Recorder

结果访问层包括：

- `ElementRecorder`
- `NodeRecorder`
- `SectionRecorder`
- `EnvelopeRecorder`

Solver 只负责求解全局状态，Recorder 负责将全局结果转换为工程上关心的节点、截面、应力和包络结果。

## 验证体系

项目使用三层验证：

### 单元测试

验证局部数学实现，例如：

- Beam3D 刚度/质量
- 荷载等效
- 坐标转换
- SparseDofReducer
- 应力公式

### 回归测试

防止已修复问题再次出现，例如：

- 集中荷载左右极限
- 单元荷载固定端力
- 非支持应力分量显式报错
- 稀疏求解奇异系统检测

### v1.0 Verification Suite

使用解析解检查完整计算链，包括：

- 静力
- 均布荷载
- 模态
- Newmark
- Rayleigh 阻尼
- 截面应力
- Recorder / Envelope

## 参考项目

- OpenSees: https://github.com/OpenSees/OpenSees
- deal.II: https://github.com/dealii/dealii
- MFEM: https://github.com/mfem/mfem
- CalculiX: http://www.calculix.de/
