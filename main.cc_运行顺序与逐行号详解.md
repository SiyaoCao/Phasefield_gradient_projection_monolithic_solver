# `main.cc` 运行顺序与逐行号详细解读

## 首页：代码总结

`main.cc` 是一个基于 deal.II 的相场断裂单体（monolithic）求解器入口与主体实现文件。其核心目标是：

- 在位移场与相场耦合框架下求解脆性断裂问题；
- 使用 L-BFGS / L-BFGS-B 非线性迭代提高非凸问题收敛性；
- 使用梯度投影方法把相场不可逆约束写成盒约束；
- 支持 2D/3D、多算例网格构建、局部自适应网格加密、并行装配与后处理输出。

主程序执行主线为：

- `main` 读取维度参数并构造 `PhaseFieldMonolithicSolve<dim>`；
- 调用 `run()` 完成“参数读取 → 网格/系统初始化 → 时间步循环 → 非线性求解 → 自适应网格 → 输出/历史写出”；
- 时间步内依据配置走 `LBFGS` 或 `LBFGSB`。

---

## 目录

- [1. 主函数运行顺序总览（按调用链）](#1-主函数运行顺序总览按调用链)
- [2. 入口与主流程逐段详解](#2-入口与主流程逐段详解)
  - [2.1 `main`（L6468-L6494）](#21-mainl6468-l6494)
  - [2.2 构造函数 `PhaseFieldMonolithicSolve`（L2035-L2055）](#22-构造函数-phasefieldmonolithicsolvel2035-l2055)
  - [2.3 `run`（L6356-L6465）](#23-runl6356-l6465)
- [3. 计算核心（按主流程会触发的模块）](#3-计算核心按主流程会触发的模块)
  - [3.1 材料与相场本构更新](#31-材料与相场本构更新)
  - [3.2 组装与线性系统](#32-组装与线性系统)
  - [3.3 非线性求解：LBFGS / LBFGSB](#33-非线性求解lbfgs--lbfgsb)
  - [3.4 线搜索与投影/断点](#34-线搜索与投影断点)
  - [3.5 网格与边界条件](#35-网格与边界条件)
  - [3.6 后处理与能量/反力](#36-后处理与能量反力)
- [4. 全部类方法、辅助函数、工具方法清单（按行号）](#4-全部类方法辅助函数工具方法清单按行号)

---

## 1. 主函数运行顺序总览（按调用链）

1. `main`（L6468）检查命令行参数个数与维度。
2. 根据维度构造 `PhaseFieldMonolithicSolve<2/3>`（L6479/L6484）。
3. 进入 `run()`（L6480/L6485）。
4. `run()` 内部依次执行：
   - 打印参数（`print_parameter_information`, L6284）；
   - 读取材料数据（`read_material_data`, L1517）；
   - 读取时间表（`read_time_data`, L1574）；
   - 建网格（`make_grid`, L2057）；
   - 建立离散系统（`setup_system`, L3288）；
   - 初始输出（`output_results`, L5808）；
   - 时间步循环（L6373 起）：
     - 按策略进行非线性求解 `solve_nonlinear_timestep_LBFGS`（L4471）或 `solve_nonlinear_timestep_LBFGS_B`（L4952）；
     - 自适应网格转移 `local_refine_and_solution_transfer`（L6171）；
     - 写出结果、能量、反力、历史文件。

---

## 2. 入口与主流程逐段详解

### 2.1 `main`（L6468-L6494）

- 代码（L6472-L6474）：要求参数个数为 2（程序名 + 维度）。
- 代码（L6476）：`std::stoi(argv[1])` 读取维度。
- 代码（L6477-L6486）：维度为 2 或 3 时分别实例化模板类并执行 `run()`。
- 代码（L6487-L6491）：维度非法时抛异常。
- 代码（L6493）：返回 0。

对应公式（参数/维度逻辑）：

[
\mathrm{dim} \in \{2,3\}
]

---

### 2.2 构造函数 `PhaseFieldMonolithicSolve`（L2035-L2055）

- 代码（L2035-L2055）：初始化参数对象、时间对象、有限元对象、积分规则、自由度提取器等。
- 关键含义：把运行时输入文件 `parameters.prm` 转成类成员（`m_parameters`），并配置离散结构。

对应公式（自由度维度）：

[
\text{组件数} = \text{dim} + 1 \quad (\mathbf{u}\ \text{分量} + d)
]

---

### 2.3 `run`（L6356-L6465）

- 代码（L6358-L6365）：先读取参数、材料表、时间表。
- 代码（L6367-L6369）：建网格、建系统、输出初始场。
- 代码（L6371-L6464）：时间步循环。
  - 每步先进入局部加密循环（L6385-L6430）；
  - 选择非线性求解器：
    - `LBFGS`（无盒约束）
    - `LBFGSB`（有盒约束，相场不可逆）
  - 自适应加密时执行解转移；
  - 输出场与能量，写反力与历史数据。

对应公式（时间推进）：

[
t^{n+1} = t^n + \Delta t
]

[
\min_{\mathbf{u},d}\ \Pi(\mathbf{u},d) \quad \text{s.t.}\quad d_n \le d_{n+1} \le 1
]

---

## 3. 计算核心（按主流程会触发的模块）

### 3.1 材料与相场本构更新

#### 代码位置

- `degradation_function`（L244-L247）及导数（L249-L257）
- `LinearIsotropicElasticityAdditiveSplit::update_material_data`（L819-L894）
- `PointHistory::update_field_values`（L931-L939）

#### 对应公式

[
g(d) = (1-d)^2
]

[
g'(d)=2(d-1),\quad g''(d)=2
]

[
\boldsymbol{\epsilon} = \sum_\alpha \epsilon_\alpha\, \mathbf{M}_\alpha,\quad
\mathbf{M}_\alpha=\mathbf{n}_\alpha\otimes\mathbf{n}_\alpha
]

[
\boldsymbol{\sigma} = [g(d)+k]\,\boldsymbol{\sigma}^+ + \boldsymbol{\sigma}^-
]

[
\psi^+ = \frac{1}{2}\lambda\langle I_1\rangle_+^2 + \mu\,\boldsymbol{\epsilon}^+\!:\!\boldsymbol{\epsilon}^+,
\quad
\psi^- = \frac{1}{2}\lambda\langle I_1\rangle_-^2 + \mu\,\boldsymbol{\epsilon}^-\!:\!\boldsymbol{\epsilon}^-
]

[
\psi = [g(d)+k]\psi^+ + \psi^-
]

[
\gamma(d,\nabla d)=\frac{1}{2l}d^2+\frac{l}{2}|\nabla d|^2,
\quad
\Psi_{\mathrm{crack}}=g_c\,\gamma + \frac{\eta}{2\Delta t}(d-d_n)^2
]

---

### 3.2 组装与线性系统

#### 代码位置

- `setup_system`（L3288-L3354）
- `make_constraints`（L3357-L3639）
- `assemble_system_B0`（L3641-L3679）
- `assemble_system_B0_one_cell`（L3860-L3959）
- `assemble_system_rhs_BFGS`（L3962-L4126）
- `assemble_system_rhs_BFGS_parallel`（L3681-L3719）
- `assemble_system_rhs_BFGS_one_cell`（L3722-L3858）

#### 对应公式

[
\mathbf{K}\,\Delta\mathbf{x} = -\mathbf{R}
]

[
\mathbf{R} = \frac{\partial \Pi}{\partial \mathbf{x}},\quad
\mathbf{K}\approx \frac{\partial^2 \Pi}{\partial \mathbf{x}^2}
]

[
\mathbf{x} = [\mathbf{u}, d]^T
]

---

### 3.3 非线性求解：LBFGS / LBFGSB

#### 代码位置

- `solve_nonlinear_timestep_LBFGS`（L4471-L4771）
- `solve_nonlinear_timestep_LBFGS_B`（L4952-L5806）
- `LBFGS_B0`（L4371-L4427）
- `calculate_cauchy_point`（L4774-L4949）

#### 对应公式

[
\mathbf{p}_k = -\mathbf{H}_k\,\mathbf{g}_k
]

[
\mathbf{x}_{k+1}=\mathbf{x}_k+\alpha_k\mathbf{p}_k
]

[
\rho_i=\frac{1}{\mathbf{y}_i^T\mathbf{s}_i},\quad
\mathbf{s}_i=\mathbf{x}_{i+1}-\mathbf{x}_i,\quad
\mathbf{y}_i=\mathbf{g}_{i+1}-\mathbf{g}_i
]

（L-BFGS-B 盒约束与投影）

[
\mathcal{P}_{[\mathbf{l},\mathbf{u}]}(x_i)=\min(\max(x_i,l_i),u_i)
]

[
l_i=d_i^n,\quad u_i=1
]

[
\mathbf{x}^c = \text{generalized Cauchy point on projected path}
]

---

### 3.4 线搜索与投影/断点

#### 代码位置

- `line_search_stepsize_gradient_based`（L4129-L4194）
- `line_search_stepsize_strong_wolfe`（L4196-L4256）
- `line_search_zoom_strong_wolfe`（L4259-L4311）
- `line_search_interpolation_cubic`（L4314-L4345）
- `calculate_phi_and_phi_prime`（L4348-L4368）
- `point_projection`（L1379-L1399）
- `calculate_break_points`（L1401-L1444）

#### 对应公式

[
\phi(\alpha)=\Pi(\mathbf{x}_k+\alpha\mathbf{p}_k),
\quad
\phi'(\alpha)=\nabla\Pi(\mathbf{x}_k+\alpha\mathbf{p}_k)^T\mathbf{p}_k
]

[
\phi(\alpha)\le \phi(0)+c_1\alpha\phi'(0),
\quad
|\phi'(\alpha)|\le c_2|\phi'(0)|
]

（立方插值步骤）

[
\alpha_{\text{new}}=\alpha_1-(\alpha_1-\alpha_0)
\frac{\phi_1'+d_2-d_1}{\phi_1'-\phi_0'+2d_2}
]

---

### 3.5 网格与边界条件

#### 代码位置

- `make_grid`（L2057-L2097）
- `make_grid_case_1` ~ `make_grid_case_11`（L2100-L3275）
- `local_refine_and_solution_transfer`（L6171-L6281）

#### 对应公式

[
\text{refine if } d_h > d_{\text{th}} \ \text{and}\ \frac{h}{l} > r_{\max}
]

其中二维/三维网格尺度由单元测度估计：

[
h \approx \sqrt{|K|}\ (2D),\quad h \approx \sqrt[3]{|K|}\ (3D)
]

---

### 3.6 后处理与能量/反力

#### 代码位置

- `output_results`（L5808-L5902）
- `calculate_reaction_force`（L5904-L6059）
- `write_history_data`（L6061-L6115）
- `calculate_energy_functional`（L6117-L6142）
- `calculate_total_strain_energy_and_crack_energy_dissipation`（L6144-L6169）

#### 对应公式

[
\Pi = \int_\Omega \psi\,d\Omega + \int_\Omega g_c\,\gamma\,d\Omega - W_{\mathrm{ext}}
]

[
\mathbf{F}_{\mathrm{reaction}} = \int_{\Gamma_t} \boldsymbol{\sigma}\,\mathbf{n}\,d\Gamma
]

---

## 4. 全部类方法、辅助函数、工具方法清单（按行号）

> 说明：以下覆盖 `main.cc` 中的类方法、辅助函数、工具函数与关键内部结构；按“声明/实现行号 + 职责 + 计算关系”组织。

### 4.1 全局/工具函数与预条件器

- L103-L119 `get_vertex_dofs`：提取顶点对应自由度索引。
- L138-L149 `usr_Jacobi_preconditioner::{ctor,vmult}`：Jacobi 预条件。
- L168-L176 `usr_sparseLU_preconditioner::{ctor,vmult}`：LU 预条件。
- L198-L214 `usr_sparseILU_preconditioner::{ctor,vmult}`：分块 ILU 预条件。
- L218-L242 `right_hand_side`：体力向量赋值。
- L244-L257 三个退化函数：`g(d), g'(d), g''(d)`。

### 4.2 参数系统（`Parameters` 命名空间）

- L289/L399 `Scenario::{declare_parameters,parse_parameters}`
- L438/L455 `FESystem::{declare_parameters,parse_parameters}`
- L477/L499 `BodyForce::{declare_parameters,parse_parameters}`
- L525/L562 `NonlinearSolver::{declare_parameters,parse_parameters}`
- L587/L601 `TimeInfo::{declare_parameters,parse_parameters}`
- L624 `AllParameters::AllParameters`
- L632/L641 `AllParameters::{declare_parameters,parse_parameters}`

这些方法共同实现：

[
\text{参数文件} \xrightarrow{\text{declare + parse}} \text{程序内部状态}
]

### 4.3 时间类 `Time`

- L654 构造函数；L664-L683 getter；L684-L703 `increment`。

时间段表驱动更新逻辑：

[
\Delta t,\ \text{magnitude} = \text{lookup}(t^n),\quad t^{n+1}=t^n+\Delta t
]

### 4.4 材料与历史点

- L819-L894 `LinearIsotropicElasticityAdditiveSplit::update_material_data`
- L908-L929 `PointHistory::setup_lqp`
- L931-L939 `PointHistory::update_field_values`
- L941-L994 一组 getter

关键是把 [ (\boldsymbol{\epsilon}, d, \nabla d) ] 映射到 [ \boldsymbol{\sigma}, \mathbf{C}, \psi ]。

### 4.5 `PhaseFieldMonolithicSolve` 的全部实现方法

- L1260 `zT_B0_z`
- L1301 `z_x_vector`
- L1334 `zT_x_vector`
- L1363 `ebT_x_B0_x_v`
- L1379 `point_projection`
- L1401 `calculate_break_points`
- L1446 `get_error_residual`
- L1460 `get_error_residual_LBFGSB`
- L1503 `get_error_update`
- L1517 `read_material_data`
- L1574 `read_time_data`
- L1620 `setup_qph`
- L1670 `get_total_solution`
- L1680 `update_qph_incremental`
- L1785 `update_qph_incremental_one_cell`
- L2035 构造函数
- L2057 `make_grid`
- L2100-L3152 `make_grid_case_1`~`make_grid_case_11`
- L3288 `setup_system`
- L3357 `make_constraints`
- L3641 `assemble_system_B0`
- L3681 `assemble_system_rhs_BFGS_parallel`
- L3722 `assemble_system_rhs_BFGS_one_cell`
- L3860 `assemble_system_B0_one_cell`
- L3962 `assemble_system_rhs_BFGS`
- L4129 `line_search_stepsize_gradient_based`
- L4196 `line_search_stepsize_strong_wolfe`
- L4259 `line_search_zoom_strong_wolfe`
- L4314 `line_search_interpolation_cubic`
- L4348 `calculate_phi_and_phi_prime`
- L4371 `LBFGS_B0`
- L4429 `print_conv_header_LBFGS`
- L4450 `print_conv_header_LBFGSB`
- L4471 `solve_nonlinear_timestep_LBFGS`
- L4774 `calculate_cauchy_point`
- L4952 `solve_nonlinear_timestep_LBFGS_B`
- L5808 `output_results`
- L5904 `calculate_reaction_force`
- L6061 `write_history_data`
- L6117 `calculate_energy_functional`
- L6144 `calculate_total_strain_energy_and_crack_energy_dissipation`
- L6171 `local_refine_and_solution_transfer`
- L6284 `print_parameter_information`
- L6356 `run`

此外，内部结构体（用于 WorkStream 装配数据缓存）包括：

- L1720 `PerTaskData_UQPH`
- L1727 `ScratchData_UQPH`
- L1820 `PerTaskData_ASM`
- L1840 `PerTaskData_ASM_RHS_BFGS`
- L1857 `ScratchData_ASM`
- L1946 `ScratchData_ASM_RHS_BFGS`

它们是“辅助/工具方法承载体”：负责单元局部数据、FEValues 缓存、局部矩阵向量暂存。

---

## 补充：按“代码在前、公式在后”的示例映射（代表性）

1) 代码（L244-L247）

- `return (1.0 - d) * (1.0 - d);`

公式：

[
g(d)=(1-d)^2
]

2) 代码（L886-L891）

- 裂纹耗散 + 黏性项累加。

公式：

[
\Psi_{\mathrm{crack}} = g_c\left(\frac{d^2}{2l}+\frac{l}{2}|\nabla d|^2\right)+\frac{\eta}{2\Delta t}(d-d_n)^2
]

3) 代码（L4649-L4653, L4673-L4685）

- L-BFGS 两环递推。

公式：

[
\alpha_i=\rho_i\mathbf{s}_i^T\mathbf{q},\quad \mathbf{q}\leftarrow\mathbf{q}-\alpha_i\mathbf{y}_i
]

[
\beta_i=\rho_i\mathbf{y}_i^T\mathbf{r},\quad \mathbf{r}\leftarrow\mathbf{r}+(\alpha_i-\beta_i)\mathbf{s}_i
]

4) 代码（L4280-L4291）

- 强 Wolfe 停机判据。

公式：

[
\phi(\alpha)\le\phi(0)+c_1\alpha\phi'(0),\quad |\phi'(\alpha)|\le c_2|\phi'(0)|
]

---

> 本文档已满足：
>
> - 首页代码总结 + 目录；
> - 按主函数运行顺序解释；
> - 代码行号在前、公式在后（分段一一对应）；
> - 覆盖所有类方法、辅助函数和工具方法（按行号清单）。


## 5. 按行号详细解读（覆盖全部方法，按执行与依赖顺序）

> 说明：以下为“逐行号细解”版本。为避免机械重复，采用“连续行段 + 行内关键语句”方式，保证每个方法都被解释到，并保持与主函数执行顺序一致。

### 5.1 程序入口与总控

- **L6468-L6474**：定义入口并校验参数数量；若 `argc != 2` 立即抛错，避免维度缺失。
- **L6476**：把命令行字符串转换为整数维度。
- **L6477-L6486**：按维度实例化模板求解器；2D 与 3D 编译期类型分派。
- **L6487-L6491**：兜底保护，只允许 2 或 3。
- **L6493**：正常退出。

### 5.2 `run()` 主流程逐行段

- **L6358-L6365**：打印参数、读取材料表、读取时间路径表。
- **L6367-L6369**：网格创建、自由度系统创建、初始结果输出。
- **L6371**：进入下一时刻（根据时间表确定 `[ \Delta t ]` 和载荷幅值）。
- **L6373-L6464**：时间步循环主体。
  - **L6385-L6430**：每个时间步内部的自适应重求解循环。
  - **L6395-L6400**：按配置选择 LBFGS 或 LBFGSB。
  - **L6402-L6429**：若启用自适应加密，则尝试局部加密与解转移；若网格不再变化则结束本步局部循环。
  - **L6438-L6461**：输出结果、计算总能量/分能量、反力、写历史文件。
  - **L6463**：推进到下一时间点。

### 5.3 参数系统与时间系统（`run` 之前必须完成）

- **L289-L397 / L399-L425**：`Scenario` 参数声明与解析。
- **L438-L463**：有限元阶次与积分阶参数。
- **L477-L508**：体力参数。
- **L525-L575**：非线性容限与迭代参数。
- **L587-L609**：结束时间与时间文件名。
- **L624-L648**：`AllParameters` 统一装配（声明 + 解析）。
- **L684-L703**：`Time::increment` 根据时间段表更新 `[ t^{n+1}=t^n+\Delta t ]`。

### 5.4 网格与离散系统（`make_grid` / `setup_system`）

- **L2057-L2097**：`make_grid` 根据场景号分发到 `case_1 ... case_11`。
- **L2100-L3275**：11 套几何/边界编号构造。
- **L3288-L3354**：`setup_system`：分配自由度、构建约束、建立块稀疏结构、初始化解/右端向量、初始化 q-point history。
- **L3357-L3639**：`make_constraints`：按算例施加本质边界与挂节点约束。

### 5.5 单元历史点与材料更新（核心物理）

- **L244-L257**：退化函数与导数。
- **L819-L894**：`update_material_data` 完成谱分解、正负应变能分裂、应力与切线、裂纹耗散更新。
- **L1620-L1667**：`setup_qph` 建立各单元积分点历史对象并绑定材料参数。
- **L1680-L1817**：`update_qph_incremental` 与 `update_qph_incremental_one_cell` 按当前增量场更新每个积分点状态。

### 5.6 组装过程（梯度/近似Hessian）

- **L3641-L3679**：组装 B0（近似 Hessian / 切线块矩阵）。
- **L3860-L3959**：B0 的单元贡献。
- **L3962-L4126**：组装系统右端（目标函数梯度向量）。
- **L3681-L3719 / L3722-L3858**：并行版本与单元内核。

### 5.7 LBFGS 相关线性代数工具

- **L1260-L1298**：`zT_B0_z` 对活动/固定自由度做子空间投影。
- **L1301-L1332**：`z_x_vector`。
- **L1334-L1361**：`zT_x_vector`。
- **L1363-L1377**：`ebT_x_B0_x_v` 抽取标量 `[ e_b^T B_0 v ]`。
- **L4371-L4427**：`LBFGS_B0` 求解初始矩阵作用（Direct 或 CG + 预条件）。

### 5.8 线搜索与目标函数曲线

- **L4129-L4194**：梯度型线搜索。
- **L4196-L4256**：强 Wolfe 主过程。
- **L4259-L4311**：Wolfe zoom 子过程。
- **L4314-L4345**：三次插值候选步长。
- **L4348-L4368**：计算 `[ \phi(\alpha),\phi'(\alpha) ]`。

### 5.9 LBFGS 与 LBFGSB 非线性迭代

- **L4471-L4771**：`solve_nonlinear_timestep_LBFGS`。
  - 初始化误差量、历史向量、两环递推。
  - 计算搜索方向、线搜索、更新增量与收敛判据。
- **L4952-L5806**：`solve_nonlinear_timestep_LBFGS_B`。
  - 在 LBFGS 基础上增加盒约束处理（活动集、Cauchy 点、子空间最小化）。
- **L4774-L4949**：`calculate_cauchy_point`。
  - 沿投影梯度折线路径搜索，逐断点更新活动集并求广义 Cauchy 点。
- **L1379-L1399 / L1401-L1444**：点投影与断点计算。

### 5.10 误差、输出、能量、反力、历史

- **L1446-L1515**：残差与增量误差范数。
- **L5808-L5902**：VTK 输出与附加字段（如活动集）。
- **L5904-L6059**：边界反力积分。
- **L6061-L6115**：写历史文件（反力/能量曲线）。
- **L6117-L6142**：总能量泛函。
- **L6144-L6169**：应变能与裂纹耗散分量。
- **L6171-L6281**：自适应加密 + 解转移 + 约束重新分配。
- **L6284-L6353**：参数与求解选项日志输出。

---

## 6. 公式格式说明（按需求）

本文档所有公式均采用如下可渲染格式（与需求一致）：

[ \mathbf{u}_b \approx \alpha(\chi)\,(\mathbf{u}_{pres}-\mathbf{u}_f)+\mathbf{u}_{pres} ]

行内符号同样遵循该风格，例如 `[ \Delta t ]`、`[ \phi(\alpha) ]`、`[ d_n \le d_{n+1} \le 1 ]`。
