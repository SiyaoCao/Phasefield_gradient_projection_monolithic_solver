# `main.cc` 按主函数运行顺序的详细讲解（含行号与公式对应）

> 说明：
> 1. 按 `main()` 的真实执行顺序讲解。
> 2. 每段先给**代码行号**与代码含义，再给对应公式。
> 3. 对“涉及计算”的代码给出公式；简单 I/O、日志、容器操作给逻辑说明。
> 4. 文末给出**所有类方法 / 辅助函数 / 工具方法**清单（含未必每次运行都会走到的分支方法）。

---

## 1. 主入口：`main()`

### 1.1 参数检查与维度分流
- 代码行号：`6468-6494`
- 关键逻辑：
  - `argc != 2` 时报错。
  - `dim = stoi(argv[1])`。
  - `dim==2`：创建 `PhaseFieldMonolithicSolve<2>` 并 `run()`。
  - `dim==3`：创建 `PhaseFieldMonolithicSolve<3>` 并 `run()`。
  - 否则抛错。

对应公式（参数映射）：
[
\text{dim} = \operatorname{int}(\text{argv}[1]),\quad \text{dim}\in\{2,3\}
]

[
\text{ProgramFlow}=
\begin{cases}
\text{Solve}_{2D}, & \text{dim}=2 \\
\text{Solve}_{3D}, & \text{dim}=3
\end{cases}
]

---

## 2. 构造阶段（由 `main()` 触发）

### 2.1 `PhaseFieldMonolithicSolve` 构造函数
- 代码行号：`2035-2054`
- 关键逻辑：
  - 读取参数对象 `m_parameters(input_file)`。
  - 初始化网格、时间、日志、计时器、有限元系统、积分规则与块结构。
  - `m_fe = [\mathbf{u}, d]` 两个场（位移向量 + 相场标量）。

对应公式（离散未知向量）：
[
\mathbf{x} = [\mathbf{u}, d]^T
]

[
\mathbf{u}\in \mathbb{R}^{n_u},\quad d\in\mathbb{R}^{n_d},\quad \mathbf{x}\in\mathbb{R}^{n_u+n_d}
]

---

## 3. 运行主流程：`run()`

### 3.1 初始化与首帧输出
- 代码行号：`6356-6370`
- 调用顺序：
  1. `print_parameter_information()`
  2. `read_material_data(...)`
  3. `read_time_data(...)`
  4. `make_grid()`
  5. `setup_system()`
  6. `output_results()`（初始输出）

#### 3.1.1 材料读取
- 代码行号：`1517-1571`
- 从文件读取每个材料区域：`\lambda, \mu, l, g_c, \eta, k`。

对应公式（泊松比校验）：
[
\nu = \frac{\lambda}{2(\lambda+\mu)}\in[-1,0.5]
]

#### 3.1.2 时间表读取
- 代码行号：`1574-1617`
- 每行读取 `[t_0,t_1,\Delta t,\text{magnitude}]`，并校验最终 `t_1=end_time`。

对应公式（时间段）：
[
t\in[t_0,t_1],\quad \Delta t>0
]

#### 3.1.3 时间推进器增量规则
- 代码行号：`684-703`（`Time::increment`）

对应公式：
[
t_{n+1}=t_n+\Delta t_n,
\quad n\leftarrow n+1
]

### 3.2 网格生成 `make_grid()`
- 代码行号：`2057-2097`
- 根据 `scenario` 路由到 `make_grid_case_1...11`。
- 输出初始网格并计算参考体积。

对应公式（体积）：
[
V_0 = \int_{\Omega_0} 1\,\mathrm{d}\Omega
]

### 3.3 系统建立 `setup_system()`
- 代码行号：`3288-3354`
- 分配 DoF、块结构、稀疏模式、矩阵/向量；调用 `setup_qph()` 初始化积分点历史变量。

对应公式（块矩阵结构）：
[
\mathbf{K}=
\begin{bmatrix}
\mathbf{K}_{uu} & \mathbf{K}_{ud}\\
\mathbf{K}_{du} & \mathbf{K}_{dd}
\end{bmatrix},\quad
\mathbf{r}=
\begin{bmatrix}
\mathbf{r}_u\\
\mathbf{r}_d
\end{bmatrix}
]

### 3.4 积分点物理场初始化 `setup_qph()` + 材料更新
- 代码行号：`1620-1667`, `908-939`, `820-894`
- 在每个单元每个积分点创建 `PointHistory`，并设置材料。
- `update_material_data()` 完成谱分解、拉压分裂、退化应力与切线、能量更新。

对应公式（退化函数）：
[
g(d)=(1-d)^2,
\quad g'(d)=2(d-1),
\quad g''(d)=2
]

对应公式（应力分裂）：
[
\boldsymbol{\sigma}=\big(g(d)+k\big)\,\boldsymbol{\sigma}^+ + \boldsymbol{\sigma}^-
]

对应公式（总应变能密度）：
[
\psi = \big(g(d)+k\big)\psi^+ + \psi^-
]

对应公式（裂纹耗散+黏性正则）：
[
\gamma(d,\nabla d)=\frac{1}{2l}d^2+\frac{l}{2}|\nabla d|^2
]

[
\psi_\text{crack}=g_c\,\gamma(d,\nabla d) + \frac{\eta}{2\Delta t}(d-d_{n})^2
]

---

## 4. 时间步循环（`run()` 的 while 主体）

- 代码行号：`6371-6464`
- 每步流程：
  1. 非线性求解（LBFGS 或 LBFGSB）
  2. 可选自适应局部加密与解转移
  3. 输出结果
  4. 计算能量、反力、写历史
  5. 时间推进

对应判据：
[
t < t_{\text{end}} + 10^{-6}\Delta t
]

---

## 5. 非线性求解核心（按运行顺序）

## 5.1 约束施加 `make_constraints()`
- 代码行号：`3357-3638`
- 根据 `scenario` 施加位移边界、支撑、旋转位移等；首迭代含非齐次约束，后续转为齐次增量形式。

对应公式（增量 Dirichlet）：
[
\Delta \mathbf{u}_{\Gamma_D}=\hat{\mathbf{u}}(t_{n+1})-\hat{\mathbf{u}}(t_n)
]

## 5.2 残量装配 `assemble_system_rhs_BFGS_parallel()`
- 代码行号：`3681-3719`（并行入口）
- 单元核函数：`3722-3857`

对应公式（位移方程弱式残量）：
[
r_u^i = \int_\Omega \nabla^s \mathbf{N}_i : \boldsymbol{\sigma}\,\mathrm{d}\Omega
      - \int_\Omega \mathbf{N}_i\cdot\mathbf{b}\,\mathrm{d}\Omega
      - \int_{\Gamma_t} \mathbf{N}_i\cdot\bar{\mathbf{t}}\,\mathrm{d}\Gamma
]

对应公式（相场方程弱式残量）：
[
r_d^i = \int_\Omega
\Big[g_c l\,\nabla N_i\!\cdot\!\nabla d
+\big(\frac{g_c}{l}d + \frac{\eta}{\Delta t}(d-d_n) + g'(d)\psi^+\big)N_i\Big] \,\mathrm{d}\Omega
]

## 5.3 切线装配 `assemble_system_B0()`
- 代码行号：`3641-3678`（并行入口）
- 单元核函数：`3860-3959`

对应公式（位移切线块）：
[
K_{uu}^{ij}=\int_\Omega (\nabla^s\mathbf{N}_i):\mathbb{C}:(\nabla^s\mathbf{N}_j)\,\mathrm{d}\Omega
]

对应公式（相场切线块）：
[
K_{dd}^{ij}=\int_\Omega
\Big[\big(\frac{g_c}{l}+\frac{\eta}{\Delta t}+g''(d)\psi^+\big)N_iN_j
+ g_c l\,\nabla N_i\cdot\nabla N_j\Big]\,\mathrm{d}\Omega
]

## 5.4 标准 L-BFGS 分支 `solve_nonlinear_timestep_LBFGS()`
- 代码行号：`4471-4771`

对应公式（两环递推）：
[
\rho_k=\frac{1}{\mathbf{y}_k^T\mathbf{s}_k},\quad
\mathbf{q}=\nabla f_k,
\quad
\mathbf{p}_k=-H_k\nabla f_k
]

对应公式（更新）：
[
\mathbf{x}_{k+1}=\mathbf{x}_k+\alpha_k\mathbf{p}_k
]

## 5.5 受限 L-BFGS-B 分支 `solve_nonlinear_timestep_LBFGS_B()`
- 代码行号：`4953-5805`
- 包含：断点、广义 Cauchy 点、自由变量子空间解、可行性修正、线搜索。

### 5.5.1 点投影与断点
- 代码行号：`1379-1395`, `1401-1443`

对应公式（箱约束投影）：
[
d^{n}\le d \le 1
]

[
\Pi_{[d^n,1]}(x)=\min\big(1,\max(d^n,x)\big)
]

### 5.5.2 广义 Cauchy 点
- 代码行号：`4775-4949`

对应公式（沿投影梯度分段线搜索）：
[
\mathbf{x}(t)=\Pi_{[l,u]}\big(\mathbf{x}_k - t\nabla f_k\big)
]

对应公式（分段二次模型导数）：
[
\phi'(t)=\nabla m(\mathbf{x}(t))^T\mathbf{x}'(t),\quad
\phi''(t)=\mathbf{x}'(t)^T\nabla^2m\,\mathbf{x}'(t)
]

### 5.5.3 子空间最小化
- 代码行号：`5113-5623`（矩阵构造与 CG/Direct 子问题）

对应公式（紧致表示）：
[
B_k = B_k^0 - W_k M_k W_k^T
]

对应公式（自由子空间系统）：
[
Z^T B_k Z\,p = -Z^T\nabla f(\mathbf{x}^c)
]

### 5.5.4 可行性修正
- 代码行号：`5634-5663`

对应公式：
[
d_{k+1}=\Pi_{[d^n,1]}(d_k+\Delta d)
]

---

## 6. 线搜索模块

### 6.1 梯度型步长
- 代码行号：`4129-4193`

对应公式（割线思想）：
[
\Delta \alpha = -\Delta\alpha_{old}\frac{g(\alpha)^T p}{(g(\alpha)-g(\alpha_{old}))^T p}
]

### 6.2 Strong Wolfe
- 代码行号：`4196-4256`, `4260-4311`, `4315-4345`, `4349-4368`

对应公式（Wolfe 条件）：
[
\phi(\alpha) \le \phi(0)+c_1\alpha\phi'(0)
]

[
|\phi'(\alpha)|\le c_2|\phi'(0)|
]

对应公式（三次插值）：
[
\alpha = \alpha_1 - (\alpha_1-\alpha_0)
\frac{\phi_1'+d_2-d_1}{\phi_1'-\phi_0'+2d_2}
]

---

## 7. 结果输出与后处理

### 7.1 VTK 输出 `output_results()`
- 代码行号：`5808-5901`
- 输出位移、相场、活跃集状态、材料 ID、L2 投影应力。

对应公式（L2 投影）：
[
\int_\Omega N_i\,\sigma_h\,\mathrm{d}\Omega
=
\int_\Omega N_i\,\sigma\,\mathrm{d}\Omega
]

### 7.2 反力计算 `calculate_reaction_force()`
- 代码行号：`5904-6058`

对应公式（边界反力分量求和）：
[
R_j = \sum_{i\in\Gamma_{ID}} r_{u,i}^{(j)}
]

### 7.3 历史写出 `write_history_data()`
- 代码行号：`6061-6114`
- 写 `Reaction_force.hist` 与 `Energy.hist`。

### 7.4 能量积分
- 代码行号：`6117-6140`, `6144-6168`

对应公式：
[
\Pi = \int_\Omega \big(\psi_{\text{strain}} + \psi_{\text{crack}}\big)\,\mathrm{d}\Omega
]

[
E_{\text{strain}}=\int_\Omega \psi_{\text{strain}}\,\mathrm{d}\Omega,
\quad
E_{\text{crack}}=\int_\Omega \psi_{\text{crack}}\,\mathrm{d}\Omega
]

### 7.5 自适应网格与解转移
- 代码行号：`6171-6281`

对应公式（相场阈值触发）：
[
\text{if } d_h > d_{\text{th}} \text{ and } h > l\cdot r_{h/l},\ \text{then refine}
]

---

## 8. 辅助函数与工具函数（运行中被间接使用）

### 8.1 顶点自由度映射
- 代码行号：`103-119`（`get_vertex_dofs`）
- 功能：给一个顶点提取该顶点上的局部 DoF 全局编号。

### 8.2 体力场构造
- 代码行号：`218-242`（`right_hand_side`）

对应公式：
[
\mathbf{b}(\mathbf{x})=(f_x,f_y[,f_z])
]

### 8.3 预条件器包装
- 代码行号：
  - Jacobi：`138-149`
  - SparseLU：`168-176`
  - SparseILU：`198-214`

对应公式（预条件）：
[
P^{-1}A\mathbf{x}=P^{-1}\mathbf{b}
]

---

## 9. 参数系统与小工具类

### 9.1 参数声明/解析
- `Scenario`：`289-425`
- `FESystem`：`438-463`
- `BodyForce`：`477-508`
- `NonlinearSolver`：`525-575`
- `TimeInfo`：`587-609`
- `AllParameters`：`624-648`

### 9.2 时间类 `Time`
- 代码行号：`664-703`（访问器 + `increment`）

### 9.3 误差归一化
- 代码行号：`1084-1099`（`Errors::reset/normalize`）

对应公式（相对化）：
[
\varepsilon_{rel}=\frac{\varepsilon}{\varepsilon_0}
]

---

## 10. 所有类方法 / 辅助函数 / 工具方法总清单（`main.cc`）

> 下面按“功能簇”汇总，满足“包括所有类方法、辅助函数和工具方法”的要求。

1. **辅助函数**：`get_vertex_dofs(103)`, `right_hand_side(218)`, `degradation_function(244)`, `degradation_function_derivative(249)`, `degradation_function_2nd_order_derivative(254)`。
2. **预条件器类方法**：`usr_Jacobi_preconditioner::{ctor,vmult}`, `usr_sparseLU_preconditioner::{ctor,vmult}`, `usr_sparseILU_preconditioner::{ctor,vmult}`（`138-214`）。
3. **参数类方法**：`Scenario/FESystem/BodyForce/NonlinearSolver/TimeInfo/AllParameters` 的 `declare_parameters` 与 `parse_parameters`（`289-648`）。
4. **Time 类方法**：`current,end,get_delta_t,get_magnitude,get_timestep,increment`（`664-703`）。
5. **材料类方法**（`LinearIsotropicElasticityAdditiveSplit`）：访问器 + `update_material_data`（`747-894`）。
6. **积分点历史类**（`PointHistory`）：`setup_lqp`, `update_field_values` 与所有 getter（`908-994`）。
7. **主求解类核心方法**（`PhaseFieldMonolithicSolve`）：
   - 网格：`make_grid`, `make_grid_case_1...11`（`2057-3285`）
   - 系统：`setup_system`, `make_constraints`, `assemble_system_B0`, `assemble_system_rhs_BFGS_parallel`, `assemble_system_B0_one_cell`, `assemble_system_rhs_BFGS_one_cell`, `assemble_system_rhs_BFGS`（`3288-4124`）
   - 线搜索：`line_search_stepsize_gradient_based`, `line_search_stepsize_strong_wolfe`, `line_search_zoom_strong_wolfe`, `line_search_interpolation_cubic`, `calculate_phi_and_phi_prime`（`4129-4368`）
   - 优化器：`LBFGS_B0`, `solve_nonlinear_timestep_LBFGS`, `calculate_cauchy_point`, `solve_nonlinear_timestep_LBFGS_B`（`4371-5805`）
   - 后处理：`output_results`, `calculate_reaction_force`, `write_history_data`, `calculate_energy_functional`, `calculate_total_strain_energy_and_crack_energy_dissipation`, `local_refine_and_solution_transfer`, `print_parameter_information`, `run`（`5808-6465`）
   - L-BFGS-B工具子程序：`point_projection`, `calculate_break_points`, `zT_B0_z`, `z_x_vector`, `zT_x_vector`, `ebT_x_B0_x_v`（`1260-1443`）
8. **并行装配/更新用数据结构方法**：`PerTaskData_*::reset`, `ScratchData_*::reset`, `copy_local_to_global_UQPH`（`1722-2030`, `1201-1202`）。
9. **主函数**：`main`（`6468-6494`）。

---

## 11. 运行顺序总览（简图）

- `main(6468)`
  - `PhaseFieldMonolithicSolve::<ctor>(2035)`
  - `run(6356)`
    - `print_parameter_information(6284)`
    - `read_material_data(1517)`
    - `read_time_data(1574)`
    - `make_grid(2057)` -> `make_grid_case_i`
    - `setup_system(3288)` -> `setup_qph(1620)`
    - 时间循环：
      - `make_constraints(3357)`
      - `assemble_system_rhs_BFGS_parallel(3681)`
      - `assemble_system_B0(3641)`
      - `solve_nonlinear_timestep_LBFGS(4471)` 或 `solve_nonlinear_timestep_LBFGS_B(4953)`
      - `local_refine_and_solution_transfer(6171)`（自适应策略）
      - `output_results(5808)`
      - `calculate_energy_functional(6117)`
      - `calculate_total_strain_energy_and_crack_energy_dissipation(6144)`
      - `calculate_reaction_force(5904)`（若启用）
      - `write_history_data(6061)`

---

## 12. 备注

- `PhaseFieldMonolithicSolve` 类内有 `determine_component_extractors` 的声明（`1129`）但在 `main.cc` 中无对应实现与调用，不参与当前运行路径。
- `assemble_system_rhs_BFGS`（串行版）与 `assemble_system_rhs_BFGS_parallel`（并行版）数学形式一致，默认流程使用并行版。



---

## 13. 按“方法”逐项解释（覆盖全部类方法/辅助函数/工具方法）

> 说明：以下按“主流程先、分支后、工具最后”的顺序列出，确保 `main.cc` 中定义的方法均被覆盖。

### 13.1 入口与总控
- `6468-6494 main`：解析维度参数并分派 2D/3D 求解。
- `2035-2054 PhaseFieldMonolithicSolve::PhaseFieldMonolithicSolve`：构造离散系统对象、积分规则、日志和时间对象。
- `6356-6465 PhaseFieldMonolithicSolve::run`：驱动完整求解生命周期（读参→建网格→装配/求解→输出→时间推进）。

### 13.2 参数子系统（全部方法）
- `289-397 Scenario::declare_parameters`：声明场景/求解器/网格策略参数。
- `399-425 Scenario::parse_parameters`：读取场景参数到成员变量。
- `438-453 FESystem::declare_parameters`；`455-463 FESystem::parse_parameters`：有限元阶次与积分阶。
- `477-497 BodyForce::declare_parameters`；`499-508 BodyForce::parse_parameters`：体力参数。
- `525-560 NonlinearSolver::declare_parameters`；`562-575 NonlinearSolver::parse_parameters`：非线性容限与最大迭代。
- `587-599 TimeInfo::declare_parameters`；`601-609 TimeInfo::parse_parameters`：终止时间与时间表文件。
- `624-630 AllParameters::AllParameters`：统一读取 `parameters.prm`。
- `632-639 AllParameters::declare_parameters`；`641-648 AllParameters::parse_parameters`：组合调用各子类参数接口。

### 13.3 时间与误差小类（全部方法）
- `664-667 current`，`668-671 end`，`672-675 get_delta_t`，`676-679 get_magnitude`，`680-683 get_timestep`：时间访问器。
- `684-703 increment`：按时间表更新 [ t_{n+1}=t_n+\Delta t_n ]。
- `1084-1089 Errors::reset`：重置误差量。
- `1091-1099 Errors::normalize`：归一化 [ \varepsilon_{rel}=\varepsilon/\varepsilon_0 ]。

### 13.4 本构与积分点历史（全部方法）
- `747-790`（`LinearIsotropicElasticityAdditiveSplit` 全部 getter）：返回应力/切线/能量/相场状态。
- `820-894 update_material_data`：谱分解、拉压分裂、应力切线和能量更新；核心关系：
  [ \boldsymbol{\sigma}=(g(d)+k)\boldsymbol{\sigma}^+ + \boldsymbol{\sigma}^- ]。
- `908-929 PointHistory::setup_lqp`：创建积分点材料对象。
- `931-939 PointHistory::update_field_values`：把当前增量场喂给本构更新。
- `941-994 PointHistory` 全部 getter：读取积分点能量、应力、材料常数等。

### 13.5 网格与系统（全部方法）
- `2057-2097 make_grid`：按场景分发网格构造并计算参考体积 [ V_0=\int_{\Omega_0}1\,\mathrm{d}\Omega ]。
- `2100-3285 make_grid_case_1 ... make_grid_case_11`：11个算例的几何、边界编号、预加密/自适应初始加密。
- `3288-3354 setup_system`：DoF 分配、稀疏结构、矩阵向量初始化、`setup_qph()`。
- `3357-3638 make_constraints`：各场景位移边界与增量约束装配。

### 13.6 装配（全部方法）
- `3641-3678 assemble_system_B0`：并行装配切线矩阵入口。
- `3860-3959 assemble_system_B0_one_cell`：单元切线贡献，公式见第5.3节。
- `3681-3719 assemble_system_rhs_BFGS_parallel`：并行装配残量入口。
- `3722-3857 assemble_system_rhs_BFGS_one_cell`：单元残量贡献，公式见第5.2节。
- `3962-4124 assemble_system_rhs_BFGS`：串行残量装配版本（同数学模型）。

### 13.7 L-BFGS / L-BFGS-B 工具方法（全部方法）
- `1260-1298 zT_B0_z`：对活跃自由度做子空间矩阵压缩。
- `1301-1331 z_x_vector`、`1334-1360 zT_x_vector`：自由变量与全量变量之间的映射。
- `1363-1376 ebT_x_B0_x_v`：取第 `b` 行矩阵-向量乘积。
- `1379-1395 point_projection`：相场投影到可行域 [ d\in[d^n,1] ]。
- `1401-1443 calculate_break_points`：沿投影梯度计算断点时刻。

### 13.8 误差、输入、QPH更新（全部方法）
- `1446-1457 get_error_residual`：无约束残量范数统计。
- `1460-1500 get_error_residual_LBFGSB`：考虑箱约束的投影残量范数。
- `1503-1514 get_error_update`：增量范数统计。
- `1517-1571 read_material_data`：材料表读入与校验。
- `1574-1617 read_time_data`：时间表读入与校验。
- `1620-1667 setup_qph`：初始化所有单元积分点历史。
- `1670-1676 get_total_solution`：返回 [ x_{tot}=x+\Delta x ]。
- `1680-1717 update_qph_incremental`：并行更新 QPH。
- `1785-1817 update_qph_incremental_one_cell`：单元级 QPH 更新。
- `1201-1202 copy_local_to_global_UQPH`：空拷贝器（该流程不需要本地到全局聚合）。

### 13.9 并行任务数据结构方法（全部）
- `1722-1723 PerTaskData_UQPH::reset`
- `1771-1781 ScratchData_UQPH::reset`
- `1832-1836 PerTaskData_ASM::reset`
- `1850-1853 PerTaskData_ASM_RHS_BFGS::reset`
- `1910-1941 ScratchData_ASM::reset`
- `1999-2030 ScratchData_ASM_RHS_BFGS::reset`

### 13.10 线搜索与核心优化（全部方法）
- `4129-4193 line_search_stepsize_gradient_based`
- `4196-4256 line_search_stepsize_strong_wolfe`
- `4260-4311 line_search_zoom_strong_wolfe`
- `4315-4345 line_search_interpolation_cubic`
- `4349-4368 calculate_phi_and_phi_prime`
- `4371-4426 LBFGS_B0`
- `4471-4771 solve_nonlinear_timestep_LBFGS`
- `4775-4949 calculate_cauchy_point`
- `4953-5805 solve_nonlinear_timestep_LBFGS_B`

### 13.11 输出与后处理（全部方法）
- `5808-5901 output_results`
- `5904-6058 calculate_reaction_force`
- `6061-6114 write_history_data`
- `6117-6140 calculate_energy_functional`
- `6144-6168 calculate_total_strain_energy_and_crack_energy_dissipation`
- `6171-6281 local_refine_and_solution_transfer`
- `6284-6353 print_parameter_information`

### 13.12 其他辅助函数（全部）
- `103-119 get_vertex_dofs`：顶点自由度索引提取。
- `218-242 right_hand_side`：生成体力向量场 [ \mathbf{b}=(f_x,f_y[,f_z]) ]。
- `244-258 degradation_function*`：退化函数及其导数。
- `138-214` 三个预条件器类的构造与 `vmult`：CG 预条件封装。

