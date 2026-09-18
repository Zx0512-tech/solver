# solver

一个基于 **C++17 + Eigen + Spectra** 实现的三维梁结构有限元求解器，重点关注清晰的单元公式、可扩展的分析架构、完整的结果恢复，以及可验证的数值计算流程。

当前 v1.0 已形成从建模、组装、求解到后处理和验证的完整链路，支持线性静力、模态分析和 Newmark-β 瞬态动力分析。

## 已实现功能

### 结构模型

- 节点、边界约束和节点荷载
- 线弹性材料与梁截面
- 确定性的全局自由度编号
- 使用 `Eigen::SparseMatrix<double>` 组装总体刚度矩阵和质量矩阵
- 总体荷载向量保持为稠密向量
- 单元荷载抽象，单元负责自身等效节点荷载

### Beam3D 三维梁单元

- 两节点三维 Euler-Bernoulli 梁/框架单元
- 每个节点 6 个自由度：
  `ux, uy, uz, rx, ry, rz`
- 轴向、扭转和双向弯曲刚度
- 12×12 一致质量矩阵
- 包含扭转转动惯量
- 稳健的局部坐标系构造
- 刚度、质量、荷载和内力的局部/全局坐标转换
- 局部变形和梁端内力恢复

### 荷载与响应

当前支持：

- 节点力与节点力矩
- `BeamUniformLoad3D`：均布线荷载
- `BeamLinearLoad3D`：三角形/梯形线荷载
- `BeamPointLoad3D`：梁任意位置集中力和集中力矩
- `BeamPartialLinearLoad3D`：作用于指定梁段的局部均布/三角形/梯形荷载
- `TimeDependentElementLoad`：将任意空间单元荷载包装为时变荷载，支持 `P(t)` 和 `q(x,t)`
- 梁荷载可在局部坐标系或全局坐标系中定义
- 基于梁形函数计算一致等效节点荷载
- 梁端内力恢复保留单元荷载固定端效应
- 静力单元响应记录
- 动力单元内力时程记录
- 梁任意位置截面内力恢复：

```text
N(x), Vy(x), Vz(x), T(x), My(x), Mz(x)
```

- `BeamResponseSampler`：均匀采样并自动插入荷载分界点
- 集中荷载位置保留左右极限，避免剪力/内力跳变被平滑掉
- CSV 输出：

```text
x, side, N, Vy, Vz, T, My, Mz
```

- 无额外绘图库依赖的 SVG 内力图：
  - 轴力图
  - 剪力图
  - 弯矩图
  - 扭矩图

Beam3D 的局部梁端内力顺序为：

```text
[N_i, Vy_i, Vz_i, T_i, My_i, Mz_i,
 N_j, Vy_j, Vz_j, T_j, My_j, Mz_j]
```

### 梁荷载示例

```cpp
// 局部坐标系均布荷载
model.addElementLoad<fem::BeamUniformLoad3D>(
    1, Eigen::Vector3d(0.0, -2000.0, 0.0));

// 任意空间方向梁上的全局竖向荷载
model.addElementLoad<fem::BeamUniformLoad3D>(
    1,
    Eigen::Vector3d(0.0, 0.0, -5000.0),
    fem::BeamLoadCoordinateSystem::Global);

// 从 qi 线性变化到 qj 的梯形荷载
model.addElementLoad<fem::BeamLinearLoad3D>(
    1,
    Eigen::Vector3d(0.0, -1000.0, 0.0),
    Eigen::Vector3d(0.0, -3000.0, 0.0));

// 距节点 i 为 1.5 m 的集中力
model.addElementLoad<fem::BeamPointLoad3D>(
    1, 1.5, Eigen::Vector3d(0.0, -10000.0, 0.0));

// 仅作用于 x=1.0 m 到 x=4.0 m 的局部梯形荷载
model.addElementLoad<fem::BeamPartialLinearLoad3D>(
    1,
    1.0,
    4.0,
    Eigen::Vector3d(0.0, -1000.0, 0.0),
    Eigen::Vector3d(0.0, -3000.0, 0.0));

// 时变集中荷载 P(t) = P0 sin(omega t)
model.addTimeDependentElementLoad<fem::BeamPointLoad3D>(
    [omega](double t) { return std::sin(omega * t); },
    1,
    1.5,
    Eigen::Vector3d(0.0, -10000.0, 0.0));

// 时变局部分布荷载 q(x,t) = q0(x) * scale(t)
model.addTimeDependentElementLoad<fem::BeamPartialLinearLoad3D>(
    [](double t) { return 1.0 + 0.5 * std::sin(20.0 * t); },
    1,
    1.0,
    4.0,
    Eigen::Vector3d(0.0, -1000.0, 0.0),
    Eigen::Vector3d(0.0, -3000.0, 0.0));
```

## 梁内力图与 CSV 输出

```cpp
const auto result = fem::LinearStaticSolver{}.solve(model);

const auto samples =
    fem::BeamResponseSampler{}.sample(
        model,
        beam_id,
        result.displacement,
        31);

fem::BeamForceCsvWriter{}.writeFile(
    "results/beam_1_forces.csv",
    samples);

fem::BeamForceDiagramSvgWriter{}.writeSet(
    "results",
    "beam_1",
    samples);
```

输出文件：

```text
beam_1_axial.svg
beam_1_shear.svg
beam_1_bending.svg
beam_1_torsion.svg
```

对于三维梁：

- 剪力图同时包含 `Vy` 和 `Vz`
- 弯矩图同时包含 `My` 和 `Mz`
- SVG 文件可直接使用浏览器打开

## 截面体系与截面应力恢复

截面层保留原有的属性输入方式，同时提供几何截面类型：

```cpp
fem::GeneralSection general(A, Iy, Iz, J);
fem::RectangleSection rectangle(size_y, size_z);
fem::CircularSection circle(radius);
```

### 截面属性

`RectangleSection` 根据截面尺寸自动计算：

- `A`
- `Iy`
- `Iz`
- Saint-Venant 扭转常数工程近似值

`CircularSection` 自动计算实心圆截面的：

- `A`
- `Iy = Iz`
- 极惯性矩/扭转常数 `J`

### 正应力恢复

任意截面均支持轴力和双向弯曲正应力：

```text
sigma_x = N/A - Mz*y/Iz + My*z/Iy
```

对于具有明确几何信息的截面，程序还会检查请求的 `(y,z)` 点是否位于截面内部。

v1.0 当前支持：

- `RectangleSection`
  - 轴力 + 双向弯曲正应力
  - 经典矩形 `Vy/Vz` 横向剪应力
- `CircularSection`
  - 轴力 + 双向弯曲正应力
  - 实心圆 Saint-Venant 扭转剪应力
- `GeneralSection`
  - 轴力 + 双向弯曲正应力

对于当前没有可靠几何公式支持的组合，程序会明确抛出异常，而不是静默忽略某个应力分量。

示例：

```cpp
const double sigma_x =
    model.beamNormalStressAt(
        beam_id,
        x,
        y,
        z,
        result.displacement);

const auto stress =
    model.beamStressAt(
        beam_id,
        x,
        y,
        z,
        result.displacement);

const auto extrema =
    model.beamNormalStressExtrema(
        beam_id,
        x,
        result.displacement);
```

`beamNormalStressAt` 始终可以恢复轴力和弯曲产生的正应力。

`beamNormalStressExtrema` 可返回矩形截面和实心圆截面的：

- `sigma_min`
- `sigma_max`
- 对应截面坐标

## 稀疏矩阵组装与求解

总体结构矩阵直接组装为压缩稀疏矩阵：

```text
单元 Ke / Me
      ↓
全局 DOF 映射
      ↓
Eigen::Triplet<double>
      ↓
setFromTriplets()
      ↓
SparseMatrix K / M
```

共享节点产生的重复 Triplet 会在组装阶段自动累加，因此不需要先创建稠密总体矩阵。

### 静力稀疏求解

```text
Sparse Global K
      ↓
SparseDofReducer
      ↓
Sparse Kff
      ↓
Eigen::SimplicialLDLT
      ↓
Uf
```

静力平衡方程为：

```text
Kff * Uf = Ff
```

求得位移后，通过：

```text
R = K U - F
```

恢复支座反力。

### Newmark-β 稀疏瞬态求解

动力方程：

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

当前实现中：

- `Kff` 为稀疏矩阵
- `Mff` 为稀疏矩阵
- `Cff` 为稀疏矩阵
- `Keff` 为稀疏矩阵
- `Keff` 在时间循环外只分解一次
- 所有时间步复用同一个稀疏分解结果
- 质量矩阵同样只在计算初始加速度时分解一次

### 模态分析

模态方程：

```text
K phi = omega^2 M phi
```

当前模态求解链保持稀疏：

```text
Sparse Global K / M
      ↓
SparseDofReducer
      ↓
Sparse Kff / Mff
      ↓
Shift-Invert
      ↓
Lanczos 子空间迭代
      ↓
最低阶结构模态
```

求解器使用 Spectra 1.2.0 的对称广义特征值 `SymGEigsShiftSolver`，对 `Kff/Mff` 采用 Shift-Invert Lanczos 方法提取低阶特征值。

为兼容自由-自由结构的刚体零频模态，使用一个相对于系统特征值尺度很小的负 shift，使 `K - sigma M` 可稳定分解；求得候选特征值后过滤零/非正特征值并返回最低阶正结构模态。

返回的振型按质量矩阵进行归一化，并记录 Lanczos 的迭代次数和算子调用次数。单自由度问题作为解析特例直接计算 `lambda = K/M`。

## 分析能力

- 线性静力分析
- 支座反力恢复
- 稀疏 Shift-Invert Lanczos 广义特征值模态分析
- Rayleigh 阻尼
- Newmark-β 瞬态动力分析
- 任意节点时变荷载
- 任意单元时变荷载
- 初始位移和初始速度
- 时间相关梁端内力恢复
- 节点、截面力和应力时程记录
- 动力响应包络提取

## 结果记录与动力包络

结果层提供：

- `ElementRecorder`
- `NodeRecorder`
- `SectionRecorder`
- `EnvelopeRecorder`

示例：

```cpp
const auto node_history =
    fem::NodeRecorder{}.record(node_id, newmark_result);

const auto section_history =
    fem::SectionRecorder{}.record(
        model, beam_id, x, newmark_result);

const auto stress_history =
    fem::SectionRecorder{}.recordNormalStress(
        model, beam_id, x, y, z, newmark_result);
```

动力包络：

```cpp
const auto uy_envelope =
    fem::EnvelopeRecorder{}.record(
        node_history,
        fem::NodeResponseQuantity::Displacement,
        fem::Dof::UY);

const auto mz_envelope =
    fem::EnvelopeRecorder{}.record(
        section_history,
        fem::BeamSectionForceComponent::Mz);

const auto sigma_envelope =
    fem::EnvelopeRecorder{}.record(stress_history);
```

每个包络包含：

- 最小值、发生时间和步号
- 最大值、发生时间和步号
- 最大绝对值
- 最大绝对值对应的带符号原值
- 发生时间和步号

## 响应架构

```text
ElementLoad
   │
   ├──> 等效节点荷载 ──> 全局 F
   │
Solver ──> 全局位移 U
                     │
                     ▼
                  Element
                     │
                     ▼
              ElementResponse
           局部变形 / 梁端内力
                     │
                     ▼
                 Recorder
             静力 / 动力时程
```

基本设计原则是：

> 单元负责自身力学行为，Solver 负责总体方程求解，Recorder 只负责结果访问和整理，不重复实现单元力学。

## 构建与测试

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

项目使用 Eigen 3.4 进行线性代数计算，并使用 Spectra 1.2.0 完成稀疏 Lanczos 模态求解。

### 示例程序

梁内力图：

```bash
./build/beam_force_diagrams_example
```

截面应力恢复：

```bash
./build/beam_stress_recovery_example
```

动力包络：

```bash
./build/result_envelope_example
```

## v1.0 验证体系

发布级端到端验证：

```bash
./build/test_v1_verification
```

Verification Suite 使用结构力学解析解检查：

- 三维悬臂梁组合静力响应
- 支座反力平衡
- 均布荷载挠度和转角
- 梁任意位置剪力和弯矩
- Euler-Bernoulli 悬臂梁一阶频率
- Newmark 自由振动
- Rayleigh 阻尼目标值
- 截面正应力
- Recorder / Envelope 完整结果链路

当前完整 CI 包含 14 组 CTest，并额外运行 v1.0 Verification Suite 和性能基准测试。

详细说明：

- `docs/v1-verification.md`
- `docs/performance-benchmark.md`

## 性能基准测试

```bash
./build/solver_benchmark \
  --sizes 100,500,2000 \
  --steps 50 \
  --repeats 5 \
  --csv benchmark.csv
```

Benchmark 输出：

- 单元数量
- 全局 DOF 数
- 自由 DOF 数
- K/M 非零元数量
- 稀疏率
- 稀疏 K 近似存储量
- 等效稠密 K 存储量
- 稀疏/稠密存储比
- 组装耗时
- 静力总耗时
- Newmark 总耗时
- 静力解析解相对误差

CI 会自动上传 `benchmark.csv` 作为构建产物。

## 文档

- `docs/architecture.md`：求解器架构
- `docs/v1-verification.md`：v1.0 端到端验证
- `docs/performance-benchmark.md`：性能基准测试
