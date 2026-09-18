# v1.0 验证套件

v1.0 Verification Suite 是求解器的发布级端到端验收层。

它不是对已有单元测试的简单重复，而是使用结构力学解析解验证从模型输入、有限元组装、求解到工程结果恢复的完整计算链。

## 运行方式

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
./build/test_v1_verification
```

该套件同时注册为 CTest：

```text
v1_verification
```

## 验收算例

### 1. 三维悬臂梁组合静力响应

对单根 Beam3D 同时施加：

- 轴向力
- 局部 y 方向横向力
- 局部 z 方向横向力
- 扭矩

验证：

```text
ux = Fx L / EA
uy = Fy L^3 / (3 E Iz)
uz = Fz L^3 / (3 E Iy)
rx = Mx L / GJ
```

并检查相应支座反力和支座弯矩是否满足整体平衡。

### 2. 均布荷载与截面内力平衡

对悬臂梁施加全跨局部 y 向均布荷载。

验证：

```text
u_tip = q L^4 / (8 E Iz)
theta_tip = q L^3 / (6 E Iz)

Vy(x) = q (L - x)
Mz(x) = q (L - x)^2 / 2
```

该算例同时验证：

```text
单元荷载
  ↓
一致等效节点荷载
  ↓
总体稀疏求解
  ↓
位移
  ↓
截面内力恢复
```

### 3. Euler-Bernoulli 悬臂梁一阶弯曲模态

使用 8 个 Beam3D 单元离散均匀悬臂梁。

理论一阶圆频率：

```text
omega_1 = beta_1^2 sqrt(E Iz / (rho A L^4))

beta_1 = 1.875104068711961
```

由于有限元模态振型是对连续解析振型的空间离散逼近，因此允许存在合理的离散误差。

### 4. Newmark-β 自由振动

构造一个轴向单自由度梁体系：

- 初始位移非零
- 初始速度为零
- 无外荷载
- 无阻尼

积分一个解析振动周期，验证：

```text
u(T/2) = -u0
u(T)   =  u0
```

用于检查：

- 一致质量
- 初始加速度
- Newmark 时间积分
- 稀疏有效刚度求解

### 5. Recorder 与 Envelope 完整链路

同一动力分析结果进一步经过：

```text
NewmarkResult
  ↓
NodeRecorder
  ↓
EnvelopeRecorder
```

以及：

```text
NewmarkResult
  ↓
SectionRecorder
  ↓
EnvelopeRecorder
```

检查节点位移包络和截面轴力包络是否与理论振幅一致。

### 6. Rayleigh 阻尼目标值

使用：

```text
RayleighDamping::fromModalTargets
```

根据两个目标频率和阻尼比拟合 Rayleigh 参数。

随后在两个目标频率处重新计算阻尼比，检查是否回到给定值。

### 7. 矩形截面正应力

验证：

```text
静力求解
  ↓
截面 N / M
  ↓
BeamSectionStressRecovery
  ↓
sigma_x
```

并与：

```text
sigma_x = N/A - Mz*y/Iz + My*z/Iy
```

解析公式比较。

## 与普通回归测试的区别

Verification Suite 不替代已有测试。

现有测试仍负责检查局部边界问题，包括：

- Beam3D 局部/全局坐标转换
- 一致质量矩阵
- 均布、线性、局部和时变荷载
- 集中荷载位置内力跳变
- CSV / SVG 内力输出
- 截面几何和应力恢复
- Recorder
- 稀疏矩阵组装
- 稀疏约化系统与静力/Newmark 求解

Verification Suite 负责更高一层的：

> 完整分析流程是否最终回到正确的结构力学答案。

## 当前验收状态

当前 CI 中：

```text
14 / 14 CTest 通过
22 / 22 v1.0 端到端检查通过
```

当前最大相对误差来自有限元离散后的悬臂梁一阶弯曲模态，而不是静力或动力线性方程求解误差。

## 稀疏 Lanczos 回归验证

除发布级解析模态验证外，新增 `sparse_modal` 回归测试：

- 20 单元悬臂梁提取前 4 阶模态
- 与独立 dense generalized eigensolve 参考结果交叉比较
- 检查自由 DOF 子空间残差

```text
Kff phi - lambda Mff phi
```

- 验证返回后端为 `SparseLanczosShiftInvert`
- 检查 Lanczos 迭代次数和算子调用次数
- 自由-自由 Beam3D 模型包含刚体零频模态时，验证程序能够跳过零模态并恢复最低阶正结构模态

dense 求解仅存在于测试代码中作为小规模参考答案，不再属于生产 `ModalSolver` 的求解路径。
