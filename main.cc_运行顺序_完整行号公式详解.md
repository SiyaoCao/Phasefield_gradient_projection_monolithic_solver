# main.cc 运行顺序、完整代码与公式对应详解

## 首页：代码总结

本文档针对 `/home/runner/work/Phasefield_gradient_projection_monolithic_solver/Phasefield_gradient_projection_monolithic_solver/main.cc` 做**按主函数运行顺序**的完整说明，满足以下目标：

1. 给出主流程（`main()` → `PhaseFieldMonolithicSolve::run()`）的分段解释，并明确对应**行号**。
2. 对涉及计算的关键代码段给出数学表达，并采用统一格式：`[ ... ]`。
3. 以一个 Markdown 文档记录全部说明内容，包含：
   - 类方法（class methods）
   - 辅助函数（helper functions）
   - 工具函数（utility functions）
4. 提供 `main.cc` 的**完整代码展示（带行号）**，确保可逐行定位。

---

## 目录

- [1. 主函数运行顺序总览（含行号）](#1-主函数运行顺序总览含行号)
- [2. 按运行顺序的分段代码与解释（代码在前、公式在后）](#2-按运行顺序的分段代码与解释代码在前公式在后)
  - [2.1 程序入口 main()](#21-程序入口-main)
  - [2.2 run() 顶层时间步驱动](#22-run-顶层时间步驱动)
  - [2.3 参数/时间/材料读取与网格系统初始化](#23-参数时间材料读取与网格系统初始化)
  - [2.4 非线性求解主干：LBFGS 与 LBFGS-B](#24-非线性求解主干lbfgs-与-lbfgs-b)
  - [2.5 L-BFGS-B 约束处理：投影、断点、广义 Cauchy 点与子空间求解](#25-l-bfgs-b-约束处理投影断点广义-cauchy-点与子空间求解)
  - [2.6 线搜索、能量与后处理输出](#26-线搜索能量与后处理输出)
  - [2.7 本构更新与谱分解相关计算](#27-本构更新与谱分解相关计算)
- [3. 全部类方法/辅助函数/工具函数覆盖表（含行号）](#3-全部类方法辅助函数工具函数覆盖表含行号)
- [4. 完整 main.cc 代码（逐行行号）](#4-完整-maincc-代码逐行行号)

---

## 1. 主函数运行顺序总览（含行号）

- 程序入口：`main()`，行 **6468-6494**。
- 根据命令行维度参数创建对象并调用 `run()`：行 **6476-6485**。
- `run()` 主驱动：行 **6356-6465**。
- 在每个时间步中：
  1. 非线性求解（`solve_nonlinear_timestep_LBFGS` 或 `solve_nonlinear_timestep_LBFGS_B`）
  2. 可能触发局部自适应加密与解转移（`local_refine_and_solution_transfer`）
  3. 输出结果与能量统计（`output_results`, `calculate_energy_functional`, `calculate_total_strain_energy_and_crack_energy_dissipation`）
  4. 反力计算与历史写出（`calculate_reaction_force`, `write_history_data`）

---

## 2. 按运行顺序的分段代码与解释（代码在前、公式在后）

> 说明：以下按主流程执行顺序解释。每段先给代码（含原始行号），后给公式与说明。根据要求采用“分段对应”，不强制逐行一一对应。

### 2.1 程序入口 main()

**代码段（行 6468-6494）：**

```cpp
6468 | int main(int argc, char* argv[])
6469 | {
6470 |   using namespace dealii;
6471 | 
6472 |   if (argc != 2)
6473 |     AssertThrow(false,
6474 |                 ExcMessage("The number of arguments provided to the program has to be 2!"));
6475 | 
6476 |   const unsigned int dim = std::stoi(argv[1]);
6477 |   if (dim == 2 )
6478 |     {
6479 |       PhaseField::PhaseFieldMonolithicSolve<2> problem_2D("parameters.prm");
6480 |       problem_2D.run();
6481 |     }
6482 |   else if (dim == 3)
6483 |     {
6484 |       PhaseField::PhaseFieldMonolithicSolve<3> problem_3D("parameters.prm");
6485 |       problem_3D.run();
6486 |     }
6487 |   else
6488 |     {
6489 |       AssertThrow(false,
6490 |                   ExcMessage("Dimension has to be either 2 or 3"));
6491 |     }
6492 | 
6493 |   return 0;
6494 | }
```

**对应公式与说明：**

- 维度选择本质是离散空间维度分支：
  [ \dim(\Omega) \in \{2,3\} ]
- `run()` 触发整个时间增量求解流程：
  [ (\mathbf{u}_{n+1}, d_{n+1}) = \arg\min\limits_{\mathbf{u},d} \Pi(\mathbf{u}, d) \ \text{s.t.}\ d_n \le d_{n+1} \le 1 ]

---

### 2.2 run() 顶层时间步驱动

**代码段（行 6356-6465）：**

```cpp
6356 | void PhaseFieldMonolithicSolve<dim>::run()
6357 | {
6358 |   print_parameter_information();
6359 |   read_material_data(...);
6360 |   read_time_data(...);
6361 |   make_grid();
6362 |   setup_system();
6363 |   output_results();
6364 | 
6365 |   m_time.increment(time_table);
6366 |   while(m_time.current() < m_time.end() + m_time.get_delta_t()*1.0e-6)
6367 |   {
6368 |     ... 非线性求解 + 自适应加密 + 输出 + 能量 + 反力 + 历史写出 ...
6369 |     m_time.increment(time_table);
6370 |   }
6371 | }
```

**对应公式与说明：**

- 时间步推进：
  [ t_{n+1} = t_n + \Delta t_n ]
- 每步最优化问题：
  [ (\mathbf{u}_{n+1}, d_{n+1}) = \arg\min \Pi(\mathbf{u}, d) ]
- 不可逆约束（盒约束）：
  [ d_n \le d_{n+1} \le 1 ]

---

### 2.3 参数/时间/材料读取与网格系统初始化

**代码段（行 624-648, 651-711, 1517-1618, 2057-3356，节选）：**

```cpp
 624 | AllParameters::AllParameters(const std::string &input_file)
 ...
 641 | void AllParameters::parse_parameters(ParameterHandler &prm)
 ...
 651 | class Time { ... void increment(std::vector<std::array<double, 4>> time_table) ... }
 ...
1517 | void PhaseFieldMonolithicSolve<dim>::read_material_data(...)
1574 | void PhaseFieldMonolithicSolve<dim>::read_time_data(...)
2057 | void PhaseFieldMonolithicSolve<dim>::make_grid()
3288 | void PhaseFieldMonolithicSolve<dim>::setup_system()
3357 | void PhaseFieldMonolithicSolve<dim>::make_constraints(const unsigned int it_nr)
```

**对应公式与说明：**

- 材料弹性常数与断裂参数被读入并映射到材料区域：
  [ (\lambda,\mu,l,g_c,\eta,k_{res})_m \leftarrow \text{material file} ]
- 载荷幅值与时间表驱动：
  [ \Delta t = \Delta t(t),\ \text{magnitude}=\mathcal{M}(t) ]
- 离散后自由度构建：
  [ \mathbf{x} = [\mathbf{u}, d]^T \in \mathbb{R}^N ]

---

### 2.4 非线性求解主干：LBFGS 与 LBFGS-B

**代码段（行 4471-4771 与 4952-5805，节选）：**

```cpp
4471 | void PhaseFieldMonolithicSolve<dim>::solve_nonlinear_timestep_LBFGS(...)
 ...
4640 | LBFGS_q_vector = m_system_rhs;
 ...
4670 | LBFGS_B0(LBFGS_r_vector, LBFGS_q_vector);
 ...
4687 | LBFGS_r_vector *= -1.0; // search direction
 ...
4714 | LBFGS_r_vector *= line_search_parameter;
 ...
4774 | void PhaseFieldMonolithicSolve<dim>::calculate_cauchy_point(...)
 ...
4952 | void PhaseFieldMonolithicSolve<dim>::solve_nonlinear_timestep_LBFGS_B(...)
 ...
```

**对应公式与说明：**

- L-BFGS 两循环核心：
  [ \mathbf{q} \leftarrow \mathbf{g}_k,\ \alpha_i = \rho_i \mathbf{s}_i^T\mathbf{q},\ \mathbf{q}\leftarrow \mathbf{q}-\alpha_i\mathbf{y}_i ]
  [ \mathbf{r}\leftarrow H_k^0\mathbf{q},\ \beta_i=\rho_i\mathbf{y}_i^T\mathbf{r},\ \mathbf{r}\leftarrow \mathbf{r}+(\alpha_i-\beta_i)\mathbf{s}_i ]
  [ \mathbf{p}_k = -\mathbf{r} ]
- 更新步：
  [ \mathbf{x}_{k+1} = \mathbf{x}_k + \alpha_k \mathbf{p}_k ]
- BFGS 曲率条件（代码中 `s_dot_y` 判据）：
  [ \mathbf{s}_k^T\mathbf{y}_k > \varepsilon\, \|\mathbf{y}_k\|^2 ]

---

### 2.5 L-BFGS-B 约束处理：投影、断点、广义 Cauchy 点与子空间求解

**代码段（行 1379-1444, 1260-1378, 4774-4950, 5113-5628，节选）：**

```cpp
1379 | void PhaseFieldMonolithicSolve<dim>::point_projection(BlockVector<double> & solution_delta)
 ...
1260 | void PhaseFieldMonolithicSolve<dim>::zT_B0_z(...)
1301 | void PhaseFieldMonolithicSolve<dim>::z_x_vector(...)
1334 | void PhaseFieldMonolithicSolve<dim>::zT_x_vector(...)
1363 | double PhaseFieldMonolithicSolve<dim>::ebT_x_B0_x_v(...)
 ...
4774 | void PhaseFieldMonolithicSolve<dim>::calculate_cauchy_point(...)
 ...
5280 | if (m_parameters.m_type_linear_solver == "CG") ...
5481 | else if (m_parameters.m_type_linear_solver == "Direct") ...
```

**对应公式与说明：**

- 点投影（盒约束）：
  [ d_i \leftarrow \min\{1,\max\{d_i^n, d_i\}\} ]
- 文中要求的符号风格示例（同格式）：
  [ \mathbf{u}b \approx \alpha(\chi),(\mathbf{u}{pres}-\mathbf{u}f)+\mathbf{u}{pres} ]
- 广义 Cauchy 点沿分段线性可行路径：
  [ \mathbf{x}^c = \mathbf{x} + t^* \mathbf{d},\ t^* = \arg\min_{t\ge0} m_k(\mathbf{x}+t\mathbf{d}) ]
- 子空间方程（自由变量）：
  [ Z^T B_k Z\,\mathbf{p} = -Z^T\hat{g} ]
- 紧凑 L-BFGS-B 形式：
  [ B_k = B_k^0 - W M W^T ]

---

### 2.6 线搜索、能量与后处理输出

**代码段（行 4129-4370, 5808-6168, 5904-6114，节选）：**

```cpp
4129 | double line_search_stepsize_gradient_based(...)
4196 | double line_search_stepsize_strong_wolfe(...)
4259 | double line_search_zoom_strong_wolfe(...)
4314 | double line_search_interpolation_cubic(...)
4348 | std::pair<double, double> calculate_phi_and_phi_prime(...)
 ...
5808 | void output_results() const
5904 | void calculate_reaction_force(unsigned int face_ID)
6061 | void write_history_data()
6117 | double calculate_energy_functional() const
6143 | std::pair<double, double> calculate_total_strain_energy_and_crack_energy_dissipation() const
```

**对应公式与说明：**

- 强 Wolfe 条件：
  [ \phi(\alpha) \le \phi(0) + c_1\alpha\phi'(0) ]
  [ |\phi'(\alpha)| \le c_2|\phi'(0)| ]
- 三次插值（`line_search_interpolation_cubic`）使用两点函数值与导数构造：
  [ \alpha_{new} = \text{cubic}(\alpha_0,\phi_0,\phi'_0,\alpha_1,\phi_1,\phi'_1) ]
- 总能量泛函离散积分：
  [ \Pi_h = \sum_{K\in\mathcal{T}_h} \sum_q (\psi_{tot} + \gamma_{crack}) J_xW ]
- 反力汇总：
  [ \mathbf{R} = \sum_{i\in\Gamma_{ID}} r_i\,\mathbf{e}_{(i\bmod dim)} ]

---

### 2.7 本构更新与谱分解相关计算

**代码段（行 244-258, 792-894, 877-884，节选）：**

```cpp
244 | double degradation_function(const double d) { return (1.0 - d) * (1.0 - d); }
249 | double degradation_function_derivative(const double d) { return 2.0 * (d - 1.0); }
254 | double degradation_function_2nd_order_derivative(const double d) { return 2.0; }
 ...
820 | void LinearIsotropicElasticityAdditiveSplit<dim>::update_material_data(...)
857 | stress_positive = ...
864 | m_stress = degradation * stress_positive + stress_negative;
876 | m_strain_energy_positive = ...
884 | m_strain_energy_total = degradation * m_strain_energy_positive + m_strain_energy_negative;
886 | m_crack_energy_dissipation = ... + viscosity_term;
```

**对应公式与说明：**

- 退化函数与导数：
  [ g(d)=(1-d)^2 ]
  [ g'(d)=2(d-1) ]
  [ g''(d)=2 ]
- 加性分裂应力：
  [ \boldsymbol{\sigma} = (g(d)+k)\boldsymbol{\sigma}^+ + \boldsymbol{\sigma}^- ]
- 正负能量：
  [ \psi^+ = \frac{\lambda}{2}\langle \operatorname{tr}(\varepsilon)\rangle_+^2 + \mu\,\varepsilon^+:\varepsilon^+ ]
  [ \psi^- = \frac{\lambda}{2}\langle \operatorname{tr}(\varepsilon)\rangle_-^2 + \mu\,\varepsilon^-:\varepsilon^- ]
- 裂纹耗散+黏性项：
  [ \gamma = g_c\left(\frac{d^2}{2l}+\frac{l}{2}|\nabla d|^2\right) + \frac{\eta}{2\Delta t}(d-d_n)^2 ]

---

## 3. 全部类方法/辅助函数/工具函数覆盖表（含行号）

> 下表覆盖 main.cc 中程序运行涉及的类、方法、辅助函数和工具函数；其中“主流程调用”表示是否在 `main()->run()` 路径中直接或间接使用。

| 类别 | 名称 | 行号（定义起点） | 主流程调用 | 作用摘要 |
|---|---:|---:|---|---|
| 工具函数 | `get_vertex_dofs` | 103 | 是 | 获取顶点相关自由度索引 |
| 预条件器 | `usr_Jacobi_preconditioner::vmult` | 142 | 是 | Jacobi 预条件作用 |
| 预条件器 | `usr_sparseLU_preconditioner::vmult` | 172 | 是 | 稀疏 LU 预条件作用 |
| 预条件器 | `usr_sparseILU_preconditioner::vmult` | 204 | 是 | 分块 ILU 预条件作用 |
| 辅助函数 | `right_hand_side` | 218 | 是 | 体力向量赋值 |
| 工具函数 | `degradation_function` | 244 | 是 | 相场退化函数 |
| 工具函数 | `degradation_function_derivative` | 249 | 是 | 退化函数一阶导 |
| 工具函数 | `degradation_function_2nd_order_derivative` | 254 | 是 | 退化函数二阶导 |
| 参数类 | `Scenario/FESystem/BodyForce/NonlinearSolver/TimeInfo` 解析函数 | 289-607 | 是 | 参数声明与解析 |
| 参数类 | `AllParameters` 构造/解析 | 624-648 | 是 | 总参数对象初始化 |
| 时间类 | `Time::increment` | 684 | 是 | 时间步和幅值更新 |
| 本构类 | `LinearIsotropicElasticityAdditiveSplit::update_material_data` | 820 | 是 | 谱分解+应力应变能更新 |
| 历史点类 | `PointHistory` 一组 getter/setter | 908-994 | 是 | 积分点历史量读写 |
| 求解器类 | `zT_B0_z/z_x_vector/zT_x_vector/ebT_x_B0_x_v` | 1260/1301/1334/1363 | 是 | L-BFGS-B 子空间线代工具 |
| 求解器类 | `point_projection` | 1379 | 是 | 盒约束投影 |
| 求解器类 | `get_error_residual/get_error_residual_LBFGSB/get_error_update` | 1446/1460/1503 | 是 | 收敛指标计算 |
| IO | `read_material_data/read_time_data` | 1517/1574 | 是 | 材料表与时间表读取 |
| QPH | `setup_qph/get_total_solution/update_qph_incremental_one_cell` | 1620/1670/1785 | 是 | 积分点场更新 |
| 网格 | `make_grid` + `make_grid_case_1..11` | 2057-3287 | 是 | 场景网格与边界设置 |
| 系统 | `setup_system/make_constraints` | 3288/3357 | 是 | DoF、约束、矩阵向量初始化 |
| 装配 | `assemble_system_B0/_one_cell` | 3641/3860 | 是 | 近似 Hessian 初始矩阵装配 |
| 装配 | `assemble_system_rhs_BFGS_parallel/_one_cell` | 3681/3722 | 是 | 残量（梯度）装配 |
| 装配 | `assemble_system_rhs_BFGS` | 3962 | 是 | 串行残量装配路径 |
| 线搜索 | `line_search_stepsize_gradient_based` | 4129 | 是 | 梯度型线搜索 |
| 线搜索 | `line_search_stepsize_strong_wolfe` | 4196 | 是 | 强 Wolfe 线搜索 |
| 线搜索 | `line_search_zoom_strong_wolfe` | 4259 | 是 | Wolfe zoom 子过程 |
| 线搜索 | `line_search_interpolation_cubic` | 4314 | 是 | 三次插值 |
| 线搜索 | `calculate_phi_and_phi_prime` | 4348 | 是 | 计算 [ \phi(\alpha),\phi'(\alpha) ] |
| 线代 | `LBFGS_B0` | 4371 | 是 | 初始逆 Hessian 作用 |
| 输出 | `print_conv_header_LBFGS/LBFGSB` | 4429/4450 | 是 | 收敛日志表头 |
| 主算法 | `solve_nonlinear_timestep_LBFGS` | 4471 | 是 | 无约束 L-BFGS |
| 主算法 | `calculate_cauchy_point` | 4774 | 是 | 计算广义 Cauchy 点 |
| 主算法 | `solve_nonlinear_timestep_LBFGS_B` | 4952 | 是 | 有界约束 L-BFGS-B |
| 后处理 | `output_results` | 5808 | 是 | VTU 输出 |
| 后处理 | `calculate_reaction_force` | 5904 | 是 | 反力计算 |
| 后处理 | `write_history_data` | 6061 | 是 | 历史文件写出 |
| 能量 | `calculate_energy_functional` | 6117 | 是 | 总能量积分 |
| 能量 | `calculate_total_strain_energy_and_crack_energy_dissipation` | 6143 | 是 | 分项能量积分 |
| 自适应 | `local_refine_and_solution_transfer` | 6171 | 是 | 局部加密与解转移 |
| 输出 | `print_parameter_information` | 6284 | 是 | 运行参数打印 |
| 主流程 | `run` | 6356 | 是 | 时间步主循环 |
| 入口 | `main` | 6468 | 是 | 程序入口 |

---

## 4. 完整 main.cc 代码（逐行行号）

> 以下为 `main.cc` 完整代码，已逐行标注行号，便于与上文解释一一定位。

```cpp
   1 | /* ---------------------------------------------------------------------
   2 |  *
   3 |  * Copyright (C) 2006 - 2020 by the deal.II authors
   4 |  *
   5 |  * This file is part of the deal.II library.
   6 |  *
   7 |  * The deal.II library is free software; you can use it, redistribute
   8 |  * it, and/or modify it under the terms of the GNU Lesser General
   9 |  * Public License as published by the Free Software Foundation; either
  10 |  * version 2.1 of the License, or (at your option) any later version.
  11 |  * The full text of the license can be found in the file LICENSE.md at
  12 |  * the top level directory of deal.II.
  13 |  *
  14 |  * ---------------------------------------------------------------------
  15 | 
  16 |  *
  17 |  * Author: Tao Jin, PhD
  18 |  *         University of Ottawa, Ottawa, Ontario, Canada
  19 |  *         July 2024
  20 |  */
  21 | 
  22 | /* A monolithic scheme based on the L-BFGS method and the gradient projection method
  23 |  *  to solve the phase-field crack problem
  24 |  * 1. The phase-field method treats the phasefield irreversibility using the gradient
  25 |  *    projection method. During a load step [t_n, t_n+1], let d_n represent the phasefield
  26 |  *    at the beginning of the load step (known), then the inequality constraints
  27 |  *    d_n <= d_n+1 <= 1.0 are treated as box constraints.
  28 |  * 2. Various direct and iterative linear solvers are designed to improve the wall-clock
  29 |  *    run time.
  30 |  * 3. Using TBB for stiffness assembly and Gauss point calculation.
  31 |  * 4. Using adaptive mesh refinement.
  32 |  */
  33 | 
  34 | #include <deal.II/grid/tria.h>
  35 | #include <deal.II/grid/grid_generator.h>
  36 | #include <deal.II/grid/grid_refinement.h>
  37 | #include <deal.II/grid/grid_out.h>
  38 | #include <deal.II/grid/grid_in.h>
  39 | #include <deal.II/grid/manifold_lib.h>
  40 | 
  41 | #include <deal.II/dofs/dof_handler.h>
  42 | #include <deal.II/dofs/dof_tools.h>
  43 | #include <deal.II/dofs/dof_renumbering.h>
  44 | 
  45 | #include <deal.II/fe/fe_values.h>
  46 | #include <deal.II/fe/fe_system.h>
  47 | #include <deal.II/fe/fe_q.h>
  48 | #include <deal.II/fe/fe_dgp_monomial.h>
  49 | #include <deal.II/fe/mapping_q_eulerian.h>
  50 | 
  51 | #include <deal.II/base/timer.h>
  52 | #include <deal.II/base/quadrature_point_data.h>
  53 | #include <deal.II/base/parameter_handler.h>
  54 | 
  55 | #include <deal.II/lac/affine_constraints.h>
  56 | #include <deal.II/lac/vector.h>
  57 | #include <deal.II/lac/full_matrix.h>
  58 | #include <deal.II/lac/sparse_matrix.h>
  59 | #include <deal.II/lac/dynamic_sparsity_pattern.h>
  60 | #include <deal.II/lac/block_sparse_matrix.h>
  61 | #include <deal.II/lac/block_vector.h>
  62 | 
  63 | #include <deal.II/numerics/vector_tools.h>
  64 | #include <deal.II/numerics/matrix_tools.h>
  65 | #include <deal.II/numerics/data_out.h>
  66 | 
  67 | #include <deal.II/lac/solver_cg.h>
  68 | #include <deal.II/lac/precondition.h>
  69 | #include <deal.II/lac/linear_operator.h>
  70 | #include <deal.II/lac/packaged_operation.h>
  71 | #include <deal.II/lac/precondition_selector.h>
  72 | #include <deal.II/lac/solver_selector.h>
  73 | #include <deal.II/lac/sparse_direct.h>
  74 | 
  75 | #include <deal.II/numerics/error_estimator.h>
  76 | 
  77 | #include <deal.II/physics/elasticity/standard_tensors.h>
  78 | 
  79 | #include <deal.II/base/quadrature_point_data.h>
  80 | 
  81 | #include <deal.II/grid/grid_tools.h>
  82 | 
  83 | #include <deal.II/base/work_stream.h>
  84 | 
  85 | #include <deal.II/numerics/solution_transfer.h>
  86 | 
  87 | #include <deal.II/lac/linear_operator_tools.h>
  88 | #include <deal.II/lac/sparse_ilu.h>
  89 | 
  90 | #include <fstream>
  91 | #include <iostream>
  92 | 
  93 | #include <deal.II/base/logstream.h>
  94 | 
  95 | #include "SpectrumDecomposition.h"
  96 | #include "Utilities.h"
  97 | 
  98 | namespace PhaseField
  99 | {
 100 |   using namespace dealii;
 101 | 
 102 |   template <int dim>
 103 |   std::vector<types::global_dof_index> get_vertex_dofs(
 104 |     const typename Triangulation<dim>::active_vertex_iterator &vertex,
 105 |     const DoFHandler<dim> &dof_handler)
 106 |   {
 107 |     DoFAccessor<0, dim, dim, false> vertex_dofs(
 108 |         &(dof_handler.get_triangulation()),
 109 |         vertex->level(),
 110 |         vertex->index(),
 111 |         &dof_handler);
 112 |     const unsigned int n_dofs = dof_handler.get_fe().dofs_per_vertex;
 113 |     std::vector<types::global_dof_index> dofs(n_dofs);
 114 |     for (unsigned int i = 0; i < n_dofs; ++i)
 115 |     {
 116 |       dofs[i] = vertex_dofs.vertex_dof_index(0, i);
 117 |     }
 118 |     return dofs;
 119 |   }
 120 | 
 121 |   // Jacobi preconditioner
 122 |   class usr_Jacobi_preconditioner : public Subscriptor
 123 |   {
 124 |   public:
 125 |     usr_Jacobi_preconditioner(const BlockSparseMatrix<double> & S);
 126 | 
 127 |     void vmult(BlockVector<double> & dst,
 128 | 	       const BlockVector<double> & src) const;
 129 | 
 130 |   private:
 131 | #  if DEAL_II_VERSION_GTE(9, 7, 0)
 132 |     const ObserverPointer<const BlockSparseMatrix<double> > m_system_matrix;
 133 | #  else
 134 |     const SmartPointer<const BlockSparseMatrix<double> > m_system_matrix;
 135 | #  endif
 136 |   };
 137 | 
 138 |   usr_Jacobi_preconditioner::usr_Jacobi_preconditioner(const BlockSparseMatrix<double> & S)
 139 |   : m_system_matrix(&S)
 140 |   {}
 141 | 
 142 |   void usr_Jacobi_preconditioner::vmult(BlockVector<double> & dst,
 143 | 					const BlockVector<double> & src) const
 144 |   {
 145 |     PreconditionJacobi<BlockSparseMatrix<double>> preconditioner;
 146 |     preconditioner.initialize(*m_system_matrix, 1.0);
 147 | 
 148 |     preconditioner.vmult(dst, src);
 149 |   }
 150 | 
 151 |   // LU preconditioner
 152 |   class usr_sparseLU_preconditioner : public Subscriptor
 153 |   {
 154 |   public:
 155 |     usr_sparseLU_preconditioner(const SparseDirectUMFPACK & matrix_factorization);
 156 | 
 157 |     void vmult(BlockVector<double> & dst,
 158 | 	       const BlockVector<double> & src) const;
 159 | 
 160 |   private:
 161 | #  if DEAL_II_VERSION_GTE(9, 7, 0)
 162 |     const ObserverPointer<const SparseDirectUMFPACK > m_matrix_LU;
 163 | #  else
 164 |     const SmartPointer<const SparseDirectUMFPACK > m_matrix_LU;
 165 | #  endif
 166 |   };
 167 | 
 168 |   usr_sparseLU_preconditioner::usr_sparseLU_preconditioner(const SparseDirectUMFPACK & matrix_factorization)
 169 |   : m_matrix_LU(&matrix_factorization)
 170 |   {}
 171 | 
 172 |   void usr_sparseLU_preconditioner::vmult(BlockVector<double> & dst,
 173 | 					       const BlockVector<double> & src) const
 174 |   {
 175 |     (*m_matrix_LU).vmult(dst, src);
 176 |   }
 177 | 
 178 |   // Incomplete LU preconditioner
 179 |   class usr_sparseILU_preconditioner : public Subscriptor
 180 |   {
 181 |   public:
 182 |     usr_sparseILU_preconditioner(const SparseILU<double> & ILU_factorization_disp,
 183 | 				 const SparseILU<double> & ILU_factorization_phasefield);
 184 | 
 185 |     void vmult(BlockVector<double> & dst,
 186 | 	     const BlockVector<double> & src) const;
 187 | 
 188 |   private:
 189 | #  if DEAL_II_VERSION_GTE(9, 7, 0)
 190 |     const ObserverPointer<const SparseILU<double> > m_ILU_factorization_disp;
 191 |     const ObserverPointer<const SparseILU<double> > m_ILU_factorization_phasefield;
 192 | #  else
 193 |     const SmartPointer<const SparseILU<double> > m_ILU_factorization_disp;
 194 |     const SmartPointer<const SparseILU<double> > m_ILU_factorization_phasefield;
 195 | #  endif
 196 |   };
 197 | 
 198 |   usr_sparseILU_preconditioner::usr_sparseILU_preconditioner(const SparseILU<double> & ILU_factorization_disp,
 199 | 							     const SparseILU<double> & ILU_factorization_phasefield)
 200 |   : m_ILU_factorization_disp(& ILU_factorization_disp)
 201 |   , m_ILU_factorization_phasefield(& ILU_factorization_phasefield)
 202 |   {}
 203 | 
 204 |   void usr_sparseILU_preconditioner::vmult(BlockVector<double> & dst,
 205 | 					   const BlockVector<double> & src) const
 206 |   {
 207 |     std::vector<types::global_dof_index> sizes_per_block(src.n_blocks());
 208 |     for (unsigned int i = 0; i < src.n_blocks(); ++i)
 209 |       sizes_per_block[i] = src.block(i).size();
 210 |     dst.reinit(sizes_per_block);
 211 | 
 212 |     (*m_ILU_factorization_disp).vmult(dst.block(0), src.block(0));
 213 |     (*m_ILU_factorization_phasefield).vmult(dst.block(1), src.block(1));
 214 |   }
 215 | 
 216 |   // body force
 217 |   template <int dim>
 218 |   void right_hand_side(const std::vector<Point<dim>> &points,
 219 | 		       std::vector<Tensor<1, dim>> &  values,
 220 | 		       const double fx,
 221 | 		       const double fy,
 222 | 		       const double fz)
 223 |   {
 224 |     Assert(values.size() == points.size(),
 225 |            ExcDimensionMismatch(values.size(), points.size()));
 226 |     Assert(dim >= 2, ExcNotImplemented());
 227 | 
 228 |     for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
 229 |       {
 230 | 	if (dim == 2)
 231 | 	  {
 232 | 	    values[point_n][0] = fx;
 233 | 	    values[point_n][1] = fy;
 234 | 	  }
 235 | 	else
 236 | 	  {
 237 | 	    values[point_n][0] = fx;
 238 | 	    values[point_n][1] = fy;
 239 | 	    values[point_n][2] = fz;
 240 | 	  }
 241 |       }
 242 |   }
 243 | 
 244 |   double degradation_function(const double d)
 245 |   {
 246 |     return (1.0 - d) * (1.0 - d);
 247 |   }
 248 | 
 249 |   double degradation_function_derivative(const double d)
 250 |   {
 251 |     return 2.0 * (d - 1.0);
 252 |   }
 253 | 
 254 |   double degradation_function_2nd_order_derivative(const double d)
 255 |   {
 256 |     (void) d;
 257 |     return 2.0;
 258 |   }
 259 | 
 260 |   namespace Parameters
 261 |   {
 262 |     struct Scenario
 263 |     {
 264 |       unsigned int m_scenario;
 265 |       std::string m_logfile_name;
 266 |       bool m_output_iteration_history;
 267 |       bool m_plane_stress;
 268 |       std::string m_type_nonlinear_solver;
 269 |       std::string m_type_line_search;
 270 |       std::string m_type_linear_solver;
 271 |       std::string m_type_preconditioner;
 272 |       double m_CG_tolerace;
 273 |       std::string m_refinement_strategy;
 274 |       unsigned int m_LBFGS_m;
 275 |       unsigned int m_global_refine_times;
 276 |       unsigned int m_local_prerefine_times;
 277 |       unsigned int m_max_adaptive_refine_times;
 278 |       int m_max_allowed_refinement_level;
 279 |       double m_phasefield_refine_threshold;
 280 |       double m_allowed_max_h_l_ratio;
 281 |       unsigned int m_total_material_regions;
 282 |       std::string m_material_file_name;
 283 |       int m_reaction_force_face_id;
 284 | 
 285 |       static void declare_parameters(ParameterHandler &prm);
 286 |       void parse_parameters(ParameterHandler &prm);
 287 |     };
 288 | 
 289 |     void Scenario::declare_parameters(ParameterHandler &prm)
 290 |     {
 291 |       prm.enter_subsection("Scenario");
 292 |       {
 293 |         prm.declare_entry("Scenario number",
 294 |                           "1",
 295 |                           Patterns::Integer(0),
 296 |                           "Geometry, loading and boundary conditions scenario");
 297 | 
 298 |         prm.declare_entry("Log file name",
 299 | 			  "Output.log",
 300 |                           Patterns::FileName(Patterns::FileName::input),
 301 | 			  "Name of the file for log");
 302 | 
 303 |         prm.declare_entry("Output iteration history",
 304 | 			  "yes",
 305 |                           Patterns::Selection("yes|no"),
 306 | 			  "Shall we write iteration history to the log file?");
 307 | 
 308 |         prm.declare_entry("Plane stress",
 309 | 			  "no",
 310 | 			  Patterns::Selection("yes|no"),
 311 | 			  "If it is 2D, is it plane-stress?");
 312 | 
 313 |         prm.declare_entry("Nonlinear solver type",
 314 |                           "LBFGSB",
 315 |                           Patterns::Selection("LBFGS|LBFGSB"),
 316 |                           "Type of solver used to solve the nonlinear system");
 317 | 
 318 |         prm.declare_entry("Line search type",
 319 |                           "GradientBased",
 320 |                           Patterns::Selection("GradientBased|StrongWolfe"),
 321 |                           "Type of line search method, the gradient-based method "
 322 |                           "should be preferred since it is generally faster");
 323 | 
 324 |         prm.declare_entry("Linear solver type",
 325 |                           "CG",
 326 |                           Patterns::Selection("Direct|CG"),
 327 |                           "Type of solver used to solve the linear system");
 328 | 
 329 |         prm.declare_entry("Preconditioner type for CG",
 330 |                           "ILU",
 331 |                           Patterns::Selection("None|Jacobi|LU|ILU"),
 332 |                           "Type of preconditioner used to solve the linear system");
 333 | 
 334 |         prm.declare_entry("CG tolerance",
 335 |                           "1.0e-6",
 336 |                           Patterns::Double(0.0),
 337 |                           "Convergence tolerance of CG iterations");
 338 | 
 339 |         prm.declare_entry("Mesh refinement strategy",
 340 |                           "adaptive-refine",
 341 |                           Patterns::Selection("pre-refine|adaptive-refine"),
 342 |                           "Mesh refinement strategy: pre-refine or adaptive-refine");
 343 | 
 344 |         prm.declare_entry("LBFGS m",
 345 |                           "40",
 346 |                           Patterns::Integer(0),
 347 |                           "Number of vectors used for LBFGS");
 348 | 
 349 |         prm.declare_entry("Global refinement times",
 350 |                           "0",
 351 |                           Patterns::Integer(0),
 352 |                           "Global refinement times (across the entire domain)");
 353 | 
 354 |         prm.declare_entry("Local prerefinement times",
 355 |                           "0",
 356 |                           Patterns::Integer(0),
 357 |                           "Local pre-refinement times (assume crack path is known a priori), "
 358 |                           "only refine along the crack path.");
 359 | 
 360 |         prm.declare_entry("Max adaptive refinement times",
 361 |                           "100",
 362 |                           Patterns::Integer(0),
 363 |                           "Maximum number of adaptive refinement times allowed in each step");
 364 | 
 365 |         prm.declare_entry("Max allowed refinement level",
 366 |                           "100",
 367 |                           Patterns::Integer(0),
 368 |                           "Maximum allowed cell refinement level");
 369 | 
 370 |         prm.declare_entry("Phasefield refine threshold",
 371 | 			  "0.8",
 372 | 			  Patterns::Double(),
 373 | 			  "Phasefield-based refinement threshold value");
 374 | 
 375 |         prm.declare_entry("Allowed max hl ratio",
 376 | 			  "0.25",
 377 | 			  Patterns::Double(),
 378 | 			  "Allowed maximum ratio between mesh size h and length scale l");
 379 | 
 380 |         prm.declare_entry("Material regions",
 381 |                           "1",
 382 |                           Patterns::Integer(0),
 383 |                           "Number of material regions");
 384 | 
 385 |         prm.declare_entry("Material data file",
 386 |                           "1",
 387 |                           Patterns::FileName(Patterns::FileName::input),
 388 |                           "Material data file");
 389 | 
 390 |         prm.declare_entry("Reaction force face ID",
 391 |                           "1",
 392 |                           Patterns::Integer(),
 393 |                           "Face id where reaction forces should be calculated "
 394 |                           "(negative integer means not to calculate reaction force)");
 395 |       }
 396 |       prm.leave_subsection();
 397 |     }
 398 | 
 399 |     void Scenario::parse_parameters(ParameterHandler &prm)
 400 |     {
 401 |       prm.enter_subsection("Scenario");
 402 |       {
 403 |         m_scenario = prm.get_integer("Scenario number");
 404 |         m_logfile_name = prm.get("Log file name");
 405 |         m_output_iteration_history = prm.get_bool("Output iteration history");
 406 |         m_plane_stress = prm.get_bool("Plane stress");
 407 |         m_type_nonlinear_solver = prm.get("Nonlinear solver type");
 408 |         m_type_line_search = prm.get("Line search type");
 409 |         m_type_linear_solver = prm.get("Linear solver type");
 410 |         m_type_preconditioner = prm.get("Preconditioner type for CG");
 411 |         m_CG_tolerace = prm.get_double("CG tolerance");
 412 |         m_refinement_strategy = prm.get("Mesh refinement strategy");
 413 |         m_LBFGS_m = prm.get_integer("LBFGS m");
 414 |         m_global_refine_times = prm.get_integer("Global refinement times");
 415 |         m_local_prerefine_times = prm.get_integer("Local prerefinement times");
 416 |         m_max_adaptive_refine_times = prm.get_integer("Max adaptive refinement times");
 417 |         m_max_allowed_refinement_level = prm.get_integer("Max allowed refinement level");
 418 |         m_phasefield_refine_threshold = prm.get_double("Phasefield refine threshold");
 419 |         m_allowed_max_h_l_ratio = prm.get_double("Allowed max hl ratio");
 420 |         m_total_material_regions = prm.get_integer("Material regions");
 421 |         m_material_file_name = prm.get("Material data file");
 422 |         m_reaction_force_face_id = prm.get_integer("Reaction force face ID");
 423 |       }
 424 |       prm.leave_subsection();
 425 |     }
 426 | 
 427 |     struct FESystem
 428 |     {
 429 |       unsigned int m_poly_degree;
 430 |       unsigned int m_quad_order;
 431 | 
 432 |       static void declare_parameters(ParameterHandler &prm);
 433 | 
 434 |       void parse_parameters(ParameterHandler &prm);
 435 |     };
 436 | 
 437 | 
 438 |     void FESystem::declare_parameters(ParameterHandler &prm)
 439 |     {
 440 |       prm.enter_subsection("Finite element system");
 441 |       {
 442 |         prm.declare_entry("Polynomial degree",
 443 |                           "1",
 444 |                           Patterns::Integer(0),
 445 |                           "Phase field polynomial order");
 446 | 
 447 |         prm.declare_entry("Quadrature order",
 448 |                           "2",
 449 |                           Patterns::Integer(0),
 450 |                           "Gauss quadrature order");
 451 |       }
 452 |       prm.leave_subsection();
 453 |     }
 454 | 
 455 |     void FESystem::parse_parameters(ParameterHandler &prm)
 456 |     {
 457 |       prm.enter_subsection("Finite element system");
 458 |       {
 459 |         m_poly_degree = prm.get_integer("Polynomial degree");
 460 |         m_quad_order  = prm.get_integer("Quadrature order");
 461 |       }
 462 |       prm.leave_subsection();
 463 |     }
 464 | 
 465 |     // body force (N/m^3)
 466 |     struct BodyForce
 467 |     {
 468 |       double m_x_component;
 469 |       double m_y_component;
 470 |       double m_z_component;
 471 | 
 472 |       static void declare_parameters(ParameterHandler &prm);
 473 | 
 474 |       void parse_parameters(ParameterHandler &prm);
 475 |     };
 476 | 
 477 |     void BodyForce::declare_parameters(ParameterHandler &prm)
 478 |     {
 479 |       prm.enter_subsection("Body force");
 480 |       {
 481 |         prm.declare_entry("Body force x component",
 482 | 			  "0.0",
 483 | 			  Patterns::Double(),
 484 | 			  "Body force x-component (N/m^3)");
 485 | 
 486 |         prm.declare_entry("Body force y component",
 487 | 			  "0.0",
 488 | 			  Patterns::Double(),
 489 | 			  "Body force y-component (N/m^3)");
 490 | 
 491 |         prm.declare_entry("Body force z component",
 492 | 			  "0.0",
 493 | 			  Patterns::Double(),
 494 | 			  "Body force z-component (N/m^3)");
 495 |       }
 496 |       prm.leave_subsection();
 497 |     }
 498 | 
 499 |     void BodyForce::parse_parameters(ParameterHandler &prm)
 500 |     {
 501 |       prm.enter_subsection("Body force");
 502 |       {
 503 |         m_x_component = prm.get_double("Body force x component");
 504 |         m_y_component = prm.get_double("Body force y component");
 505 |         m_z_component = prm.get_double("Body force z component");
 506 |       }
 507 |       prm.leave_subsection();
 508 |     }
 509 | 
 510 |     struct NonlinearSolver
 511 |     {
 512 |       unsigned int m_max_iterations_BFGS;
 513 |       bool m_relative_residual;
 514 | 
 515 |       double       m_tol_u_residual;
 516 |       double       m_tol_d_residual;
 517 |       double       m_tol_u_incr;
 518 |       double       m_tol_d_incr;
 519 | 
 520 |       static void declare_parameters(ParameterHandler &prm);
 521 | 
 522 |       void parse_parameters(ParameterHandler &prm);
 523 |     };
 524 | 
 525 |     void NonlinearSolver::declare_parameters(ParameterHandler &prm)
 526 |     {
 527 |       prm.enter_subsection("Nonlinear solver");
 528 |       {
 529 |         prm.declare_entry("Max iterations BFGS",
 530 |                           "20",
 531 |                           Patterns::Integer(0),
 532 |                           "Number of BFGS iterations allowed");
 533 | 
 534 |         prm.declare_entry("Relative residual",
 535 | 			  "yes",
 536 |                           Patterns::Selection("yes|no"),
 537 | 			  "Shall we use relative residual for convergence?");
 538 | 
 539 |         prm.declare_entry("Tolerance displacement residual",
 540 |                           "1.0e-9",
 541 |                           Patterns::Double(0.0),
 542 |                           "Displacement residual tolerance");
 543 | 
 544 |         prm.declare_entry("Tolerance phasefield residual",
 545 |                           "1.0e-9",
 546 |                           Patterns::Double(0.0),
 547 |                           "Phasefield residual tolerance");
 548 | 
 549 |         prm.declare_entry("Tolerance displacement increment",
 550 |                           "1.0e-9",
 551 |                           Patterns::Double(0.0),
 552 |                           "Displacement increment tolerance");
 553 | 
 554 |         prm.declare_entry("Tolerance phasefield increment",
 555 |                           "1.0e-9",
 556 |                           Patterns::Double(0.0),
 557 |                           "Phasefield increment tolerance");
 558 |       }
 559 |       prm.leave_subsection();
 560 |     }
 561 | 
 562 |     void NonlinearSolver::parse_parameters(ParameterHandler &prm)
 563 |     {
 564 |       prm.enter_subsection("Nonlinear solver");
 565 |       {
 566 |         m_max_iterations_BFGS = prm.get_integer("Max iterations BFGS");
 567 |         m_relative_residual = prm.get_bool("Relative residual");
 568 | 
 569 |         m_tol_u_residual           = prm.get_double("Tolerance displacement residual");
 570 |         m_tol_d_residual           = prm.get_double("Tolerance phasefield residual");
 571 |         m_tol_u_incr               = prm.get_double("Tolerance displacement increment");
 572 |         m_tol_d_incr               = prm.get_double("Tolerance phasefield increment");
 573 |       }
 574 |       prm.leave_subsection();
 575 |     }
 576 | 
 577 |     struct TimeInfo
 578 |     {
 579 |       double m_end_time;
 580 |       std::string m_time_file_name;
 581 | 
 582 |       static void declare_parameters(ParameterHandler &prm);
 583 | 
 584 |       void parse_parameters(ParameterHandler &prm);
 585 |     };
 586 | 
 587 |     void TimeInfo::declare_parameters(ParameterHandler &prm)
 588 |     {
 589 |       prm.enter_subsection("Time");
 590 |       {
 591 |         prm.declare_entry("End time", "1", Patterns::Double(), "End time");
 592 | 
 593 |         prm.declare_entry("Time data file",
 594 |                           "1",
 595 |                           Patterns::FileName(Patterns::FileName::input),
 596 |                           "Time data file");
 597 |       }
 598 |       prm.leave_subsection();
 599 |     }
 600 | 
 601 |     void TimeInfo::parse_parameters(ParameterHandler &prm)
 602 |     {
 603 |       prm.enter_subsection("Time");
 604 |       {
 605 |         m_end_time = prm.get_double("End time");
 606 |         m_time_file_name = prm.get("Time data file");
 607 |       }
 608 |       prm.leave_subsection();
 609 |     }
 610 | 
 611 |     struct AllParameters : public Scenario,
 612 | 	                   public FESystem,
 613 | 	                   public BodyForce,
 614 | 			   public NonlinearSolver,
 615 | 			   public TimeInfo
 616 |     {
 617 |       AllParameters(const std::string &input_file);
 618 | 
 619 |       static void declare_parameters(ParameterHandler &prm);
 620 | 
 621 |       void parse_parameters(ParameterHandler &prm);
 622 |     };
 623 | 
 624 |     AllParameters::AllParameters(const std::string &input_file)
 625 |     {
 626 |       ParameterHandler prm;
 627 |       declare_parameters(prm);
 628 |       prm.parse_input(input_file);
 629 |       parse_parameters(prm);
 630 |     }
 631 | 
 632 |     void AllParameters::declare_parameters(ParameterHandler &prm)
 633 |     {
 634 |       Scenario::declare_parameters(prm);
 635 |       FESystem::declare_parameters(prm);
 636 |       BodyForce::declare_parameters(prm);
 637 |       NonlinearSolver::declare_parameters(prm);
 638 |       TimeInfo::declare_parameters(prm);
 639 |     }
 640 | 
 641 |     void AllParameters::parse_parameters(ParameterHandler &prm)
 642 |     {
 643 |       Scenario::parse_parameters(prm);
 644 |       FESystem::parse_parameters(prm);
 645 |       BodyForce::parse_parameters(prm);
 646 |       NonlinearSolver::parse_parameters(prm);
 647 |       TimeInfo::parse_parameters(prm);
 648 |     }
 649 |   } // namespace Parameters
 650 | 
 651 |   class Time
 652 |   {
 653 |   public:
 654 |     Time(const double time_end)
 655 |       : m_timestep(0)
 656 |       , m_time_current(0.0)
 657 |       , m_time_end(time_end)
 658 |       , m_delta_t(0.0)
 659 |       , m_magnitude(1.0)
 660 |     {}
 661 | 
 662 |     virtual ~Time() = default;
 663 | 
 664 |     double current() const
 665 |     {
 666 |       return m_time_current;
 667 |     }
 668 |     double end() const
 669 |     {
 670 |       return m_time_end;
 671 |     }
 672 |     double get_delta_t() const
 673 |     {
 674 |       return m_delta_t;
 675 |     }
 676 |     double get_magnitude() const
 677 |     {
 678 |       return m_magnitude;
 679 |     }
 680 |     unsigned int get_timestep() const
 681 |     {
 682 |       return m_timestep;
 683 |     }
 684 |     void increment(std::vector<std::array<double, 4>> time_table)
 685 |     {
 686 |       double t_1, t_delta, t_magnitude;
 687 |       for (auto & time_group : time_table)
 688 |         {
 689 | 	  t_1 = time_group[1];
 690 | 	  t_delta = time_group[2];
 691 | 	  t_magnitude = time_group[3];
 692 | 
 693 | 	  if (m_time_current < t_1 - 1.0e-6*t_delta)
 694 | 	    {
 695 | 	      m_delta_t = t_delta;
 696 | 	      m_magnitude = t_magnitude;
 697 | 	      break;
 698 | 	    }
 699 |         }
 700 | 
 701 |       m_time_current += m_delta_t;
 702 |       ++m_timestep;
 703 |     }
 704 | 
 705 |   private:
 706 |     unsigned int m_timestep;
 707 |     double       m_time_current;
 708 |     const double m_time_end;
 709 |     double m_delta_t;
 710 |     double m_magnitude;
 711 |   };
 712 | 
 713 |   template <int dim>
 714 |   class LinearIsotropicElasticityAdditiveSplit
 715 |   {
 716 |   public:
 717 |     LinearIsotropicElasticityAdditiveSplit(const double lame_lambda,
 718 | 			                   const double lame_mu,
 719 | 				           const double residual_k,
 720 | 					   const double length_scale,
 721 | 					   const double viscosity,
 722 | 					   const double gc,
 723 | 					   const bool   plane_stress_flag)
 724 |       : m_lame_lambda(lame_lambda)
 725 |       , m_lame_mu(lame_mu)
 726 |       , m_residual_k(residual_k)
 727 |       , m_length_scale(length_scale)
 728 |       , m_eta(viscosity)
 729 |       , m_gc(gc)
 730 |       , m_plane_stress(plane_stress_flag)
 731 |       , m_phase_field_value(0.0)
 732 |       , m_grad_phasefield(Tensor<1, dim>())
 733 |       , m_strain(SymmetricTensor<2, dim>())
 734 |       , m_stress(SymmetricTensor<2, dim>())
 735 |       , m_stress_positive(SymmetricTensor<2, dim>())
 736 |       , m_mechanical_C(SymmetricTensor<4, dim>())
 737 |       , m_strain_energy_positive(0.0)
 738 |       , m_strain_energy_negative(0.0)
 739 |       , m_strain_energy_total(0.0)
 740 |       , m_crack_energy_dissipation(0.0)
 741 |     {
 742 |       Assert(  ( lame_lambda / (2*(lame_lambda + lame_mu)) <= 0.5)
 743 | 	     & ( lame_lambda / (2*(lame_lambda + lame_mu)) >=-1.0),
 744 | 	     ExcInternalError() );
 745 |     }
 746 | 
 747 |     const SymmetricTensor<4, dim> & get_mechanical_C() const
 748 |     {
 749 |       return m_mechanical_C;
 750 |     }
 751 | 
 752 |     const SymmetricTensor<2, dim> & get_cauchy_stress() const
 753 |     {
 754 |       return m_stress;
 755 |     }
 756 | 
 757 |     const SymmetricTensor<2, dim> & get_cauchy_stress_positive() const
 758 |     {
 759 |       return m_stress_positive;
 760 |     }
 761 | 
 762 |     double get_positive_strain_energy() const
 763 |     {
 764 |       return m_strain_energy_positive;
 765 |     }
 766 | 
 767 |     double get_negative_strain_energy() const
 768 |     {
 769 |       return m_strain_energy_negative;
 770 |     }
 771 | 
 772 |     double get_total_strain_energy() const
 773 |     {
 774 |       return m_strain_energy_total;
 775 |     }
 776 | 
 777 |     double get_crack_energy_dissipation() const
 778 |     {
 779 |       return m_crack_energy_dissipation;
 780 |     }
 781 | 
 782 |     double get_phase_field_value() const
 783 |     {
 784 |       return m_phase_field_value;
 785 |     }
 786 | 
 787 |     const Tensor<1, dim> get_phase_field_gradient() const
 788 |     {
 789 |       return m_grad_phasefield;
 790 |     }
 791 | 
 792 |     void update_material_data(const SymmetricTensor<2, dim> & strain,
 793 | 			      const double phase_field_value,
 794 | 			      const Tensor<1, dim> & grad_phasefield,
 795 | 			      const double phase_field_value_previous_step,
 796 | 			      const double delta_time);
 797 | 
 798 |   private:
 799 |     const double m_lame_lambda;
 800 |     const double m_lame_mu;
 801 |     const double m_residual_k;
 802 |     const double m_length_scale;
 803 |     const double m_eta;
 804 |     const double m_gc;
 805 |     const bool   m_plane_stress;
 806 |     double m_phase_field_value;
 807 |     Tensor<1, dim> m_grad_phasefield;
 808 |     SymmetricTensor<2, dim> m_strain;
 809 |     SymmetricTensor<2, dim> m_stress;
 810 |     SymmetricTensor<2, dim> m_stress_positive;
 811 |     SymmetricTensor<4, dim> m_mechanical_C;
 812 |     double m_strain_energy_positive;
 813 |     double m_strain_energy_negative;
 814 |     double m_strain_energy_total;
 815 |     double m_crack_energy_dissipation;
 816 |   };
 817 | 
 818 |   template <int dim>
 819 |   void LinearIsotropicElasticityAdditiveSplit<dim>::
 820 |    update_material_data(const SymmetricTensor<2, dim> & strain,
 821 | 			const double phase_field_value,
 822 | 			const Tensor<1, dim> & grad_phasefield,
 823 | 			const double phase_field_value_previous_step,
 824 | 			const double delta_time)
 825 |   {
 826 |     m_strain = strain;
 827 |     m_phase_field_value = phase_field_value;
 828 |     m_grad_phasefield = grad_phasefield;
 829 |     Vector<double>              eigenvalues(dim);
 830 |     std::vector<Tensor<1, dim>> eigenvectors(dim);
 831 |     usr_spectrum_decomposition::spectrum_decomposition<dim>(m_strain,
 832 |   							      eigenvalues,
 833 |   							      eigenvectors);
 834 | 
 835 |     SymmetricTensor<2, dim> strain_positive, strain_negative;
 836 |     strain_positive = usr_spectrum_decomposition::positive_tensor(eigenvalues, eigenvectors);
 837 |     strain_negative = usr_spectrum_decomposition::negative_tensor(eigenvalues, eigenvectors);
 838 | 
 839 |     SymmetricTensor<4, dim> projector_positive, projector_negative;
 840 |     usr_spectrum_decomposition::positive_negative_projectors(eigenvalues,
 841 |   							       eigenvectors,
 842 | 							       projector_positive,
 843 | 							       projector_negative);
 844 | 
 845 |     SymmetricTensor<2, dim> stress_positive, stress_negative;
 846 |     const double degradation = degradation_function(m_phase_field_value) + m_residual_k;
 847 |     const double I_1 = trace(m_strain);
 848 | 
 849 |     // 2D plane strain and 3D cases
 850 |     double my_lambda = m_lame_lambda;
 851 | 
 852 |     // 2D plane stress case
 853 |     if (    dim == 2
 854 | 	   && m_plane_stress)
 855 |       my_lambda = 2 * m_lame_mu * m_lame_lambda / (m_lame_lambda + 2 * m_lame_mu);
 856 | 
 857 |     stress_positive = my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)
 858 |                                     * Physics::Elasticity::StandardTensors<dim>::I
 859 |                     + 2 * m_lame_mu * strain_positive;
 860 |     stress_negative = my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)
 861 |                                     * Physics::Elasticity::StandardTensors<dim>::I
 862 |     		      + 2 * m_lame_mu * strain_negative;
 863 | 
 864 |     m_stress = degradation * stress_positive + stress_negative;
 865 |     m_stress_positive = stress_positive;
 866 | 
 867 |     SymmetricTensor<4, dim> C_positive, C_negative;
 868 |     C_positive = my_lambda * usr_spectrum_decomposition::heaviside_function(I_1)
 869 |                                * Physics::Elasticity::StandardTensors<dim>::IxI
 870 | 		 + 2 * m_lame_mu * projector_positive;
 871 |     C_negative = my_lambda * usr_spectrum_decomposition::heaviside_function(-I_1)
 872 |                                * Physics::Elasticity::StandardTensors<dim>::IxI
 873 |     		 + 2 * m_lame_mu * projector_negative;
 874 |     m_mechanical_C = degradation * C_positive + C_negative;
 875 | 
 876 |     m_strain_energy_positive = 0.5 * my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)
 877 |                                                    * usr_spectrum_decomposition::positive_ramp_function(I_1)
 878 |                              + m_lame_mu * strain_positive * strain_positive;
 879 | 
 880 |     m_strain_energy_negative = 0.5 * my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)
 881 |                                                    * usr_spectrum_decomposition::negative_ramp_function(I_1)
 882 |                              + m_lame_mu * strain_negative * strain_negative;
 883 | 
 884 |     m_strain_energy_total = degradation * m_strain_energy_positive + m_strain_energy_negative;
 885 | 
 886 |     m_crack_energy_dissipation = m_gc * (  0.5 / m_length_scale * m_phase_field_value * m_phase_field_value
 887 | 	                                   + 0.5 * m_length_scale * m_grad_phasefield * m_grad_phasefield)
 888 | 	                                   // the term due to viscosity regularization
 889 | 	                                   + (m_phase_field_value - phase_field_value_previous_step)
 890 | 					   * (m_phase_field_value - phase_field_value_previous_step)
 891 | 				           * 0.5 * m_eta / delta_time;
 892 |     //(void)delta_time;
 893 |     //(void)phase_field_value_previous_step;
 894 |   }
 895 | 
 896 |   template <int dim>
 897 |   class PointHistory
 898 |   {
 899 |   public:
 900 |     PointHistory()
 901 |       : m_length_scale(0.0)
 902 |       , m_gc(0.0)
 903 |       , m_viscosity(0.0)
 904 |     {}
 905 | 
 906 |     virtual ~PointHistory() = default;
 907 | 
 908 |     void setup_lqp(const double lame_lambda,
 909 | 		   const double lame_mu,
 910 | 		   const double length_scale,
 911 | 		   const double gc,
 912 | 		   const double viscosity,
 913 | 		   const double residual_k,
 914 | 		   const bool   plane_stress_flag)
 915 |     {
 916 |       m_material =
 917 |               std::make_shared<LinearIsotropicElasticityAdditiveSplit<dim>>(lame_lambda,
 918 |         	                                                            lame_mu,
 919 | 								            residual_k,
 920 | 									    length_scale,
 921 | 									    viscosity,
 922 | 									    gc,
 923 | 									    plane_stress_flag);
 924 |       m_length_scale = length_scale;
 925 |       m_gc = gc;
 926 |       m_viscosity = viscosity;
 927 | 
 928 |       update_field_values(SymmetricTensor<2, dim>(), 0.0, Tensor<1, dim>(), 0.0, 1.0);
 929 |     }
 930 | 
 931 |     void update_field_values(const SymmetricTensor<2, dim> & strain,
 932 | 		             const double phase_field_value,
 933 | 			     const Tensor<1, dim> & grad_phasefield,
 934 | 			     const double phase_field_value_previous_step,
 935 | 			     const double delta_time)
 936 |     {
 937 |       m_material->update_material_data(strain, phase_field_value, grad_phasefield,
 938 | 				       phase_field_value_previous_step, delta_time);
 939 |     }
 940 | 
 941 |     double get_current_positive_strain_energy() const
 942 |     {
 943 |       return m_material->get_positive_strain_energy();
 944 |     }
 945 | 
 946 |     const SymmetricTensor<4, dim> & get_mechanical_C() const
 947 |     {
 948 |       return m_material->get_mechanical_C();
 949 |     }
 950 | 
 951 |     const SymmetricTensor<2, dim> & get_cauchy_stress() const
 952 |     {
 953 |       return m_material->get_cauchy_stress();
 954 |     }
 955 | 
 956 |     const SymmetricTensor<2, dim> & get_cauchy_stress_positive() const
 957 |     {
 958 |       return m_material->get_cauchy_stress_positive();
 959 |     }
 960 | 
 961 |     double get_total_strain_energy() const
 962 |     {
 963 |       return m_material->get_total_strain_energy();
 964 |     }
 965 | 
 966 |     double get_crack_energy_dissipation() const
 967 |     {
 968 |       return m_material->get_crack_energy_dissipation();
 969 |     }
 970 | 
 971 |     double get_phase_field_value() const
 972 |     {
 973 |       return m_material->get_phase_field_value();
 974 |     }
 975 | 
 976 |     const Tensor<1, dim> get_phase_field_gradient() const
 977 |     {
 978 |       return m_material->get_phase_field_gradient();
 979 |     }
 980 | 
 981 |     double get_length_scale() const
 982 |     {
 983 |       return m_length_scale;
 984 |     }
 985 | 
 986 |     double get_critical_energy_release_rate() const
 987 |     {
 988 |       return m_gc;
 989 |     }
 990 | 
 991 |     double get_viscosity() const
 992 |     {
 993 |       return m_viscosity;
 994 |     }
 995 |   private:
 996 |     std::shared_ptr<LinearIsotropicElasticityAdditiveSplit<dim>> m_material;
 997 |     double m_length_scale;
 998 |     double m_gc;
 999 |     double m_viscosity;
1000 |   };
1001 | 
1002 |   template <int dim>
1003 |   class PhaseFieldMonolithicSolve
1004 |   {
1005 |   public:
1006 |     PhaseFieldMonolithicSolve(const std::string &input_file);
1007 | 
1008 |     virtual ~PhaseFieldMonolithicSolve() = default;
1009 |     void run();
1010 | 
1011 |   private:
1012 |     struct PerTaskData_ASM;
1013 |     struct ScratchData_ASM;
1014 | 
1015 |     struct PerTaskData_ASM_RHS_BFGS;
1016 |     struct ScratchData_ASM_RHS_BFGS;
1017 | 
1018 |     struct PerTaskData_UQPH;
1019 |     struct ScratchData_UQPH;
1020 | 
1021 |     Parameters::AllParameters m_parameters;
1022 |     Triangulation<dim> m_triangulation;
1023 | 
1024 |     CellDataStorage<typename Triangulation<dim>::cell_iterator,
1025 |                     PointHistory<dim>>
1026 |       m_quadrature_point_history;
1027 | 
1028 |     Time                m_time;
1029 |     std::ofstream m_logfile;
1030 |     mutable TimerOutput m_timer;
1031 | 
1032 |     DoFHandler<dim>                  m_dof_handler;
1033 |     FESystem<dim>                    m_fe;
1034 |     const unsigned int               m_dofs_per_cell;
1035 |     const FEValuesExtractors::Vector m_u_fe;
1036 |     const FEValuesExtractors::Scalar m_d_fe;
1037 | 
1038 |     static const unsigned int m_n_blocks          = 2;
1039 |     static const unsigned int m_n_components      = dim + 1;
1040 |     static const unsigned int m_first_u_component = 0;
1041 |     static const unsigned int m_d_component       = dim;
1042 | 
1043 |     enum
1044 |     {
1045 |       m_u_dof = 0,
1046 |       m_d_dof = 1
1047 |     };
1048 | 
1049 |     std::vector<types::global_dof_index> m_dofs_per_block;
1050 | 
1051 |     const QGauss<dim>     m_qf_cell;
1052 |     const QGauss<dim - 1> m_qf_face;
1053 |     const unsigned int    m_n_q_points;
1054 | 
1055 |     double m_vol_reference;
1056 | 
1057 |     AffineConstraints<double> m_constraints;
1058 |     BlockSparsityPattern      m_sparsity_pattern;
1059 |     BlockSparseMatrix<double> m_tangent_matrix;
1060 |     BlockVector<double>       m_system_rhs;
1061 |     BlockVector<double>       m_solution;
1062 |     SparseDirectUMFPACK       m_A_direct;
1063 |     // m_active_set_phasefield has 0 (inactive constraint)
1064 |     //                          or 1 (active constraint lower bound)
1065 |     //                          or 2 (active constraint upper bound)
1066 |     // In order to add active set into the VTK output, we have to declare
1067 |     // it as double, not int or unsigned int
1068 |     Vector<double> m_active_set_phasefield;
1069 | 
1070 |     std::map<unsigned int, std::vector<double>> m_material_data;
1071 | 
1072 |     std::vector<std::pair<double, std::vector<double>>> m_history_reaction_force;
1073 |     std::vector<std::pair<double, std::array<double, 3>>> m_history_energy;
1074 | 
1075 | 
1076 |     struct Errors
1077 |     {
1078 |       Errors()
1079 |         : m_norm(1.0)
1080 |         , m_u(1.0)
1081 |         , m_d(1.0)
1082 |       {}
1083 | 
1084 |       void reset()
1085 |       {
1086 |         m_norm = 1.0;
1087 |         m_u    = 1.0;
1088 |         m_d    = 1.0;
1089 |       }
1090 | 
1091 |       void normalize(const Errors &rhs)
1092 |       {
1093 |         if (rhs.m_norm != 0.0)
1094 |           m_norm /= rhs.m_norm;
1095 |         if (rhs.m_u != 0.0)
1096 |           m_u /= rhs.m_u;
1097 |         if (rhs.m_d != 0.0)
1098 |           m_d /= rhs.m_d;
1099 |       }
1100 | 
1101 |       double m_norm, m_u, m_d;
1102 |     };
1103 | 
1104 |     Errors m_error_residual, m_error_residual_0, m_error_residual_norm, m_error_update,
1105 |       m_error_update_0, m_error_update_norm;
1106 | 
1107 |     void get_error_residual(Errors &error_residual);
1108 |     void get_error_residual_LBFGSB(Errors &error_residual,
1109 | 				   const BlockVector<double> & solution_delta);
1110 | 
1111 |     void get_error_update(const BlockVector<double> &newton_update,
1112 |                           Errors & error_update);
1113 | 
1114 |     void make_grid();
1115 |     void make_grid_case_1();
1116 |     void make_grid_case_2();
1117 |     void make_grid_case_3();
1118 |     void make_grid_case_4();
1119 |     void make_grid_case_5();
1120 |     void make_grid_case_6();
1121 |     void make_grid_case_7();
1122 |     void make_grid_case_8();
1123 |     void make_grid_case_9();
1124 |     void make_grid_case_10();
1125 |     void make_grid_case_11();
1126 | 
1127 |     void setup_system();
1128 | 
1129 |     void determine_component_extractors();
1130 | 
1131 |     void make_constraints(const unsigned int it_nr);
1132 | 
1133 |     void assemble_system_B0(const BlockVector<double> & solution_old);
1134 | 
1135 |     void assemble_system_B0_one_cell(
1136 |       const typename DoFHandler<dim>::active_cell_iterator &cell,
1137 |       ScratchData_ASM &                                     scratch,
1138 |       PerTaskData_ASM &                                     data) const;
1139 | 
1140 |     void assemble_system_rhs_BFGS_one_cell(
1141 |       const typename DoFHandler<dim>::active_cell_iterator &cell,
1142 |       ScratchData_ASM_RHS_BFGS &                           scratch,
1143 |       PerTaskData_ASM_RHS_BFGS &                           data) const;
1144 | 
1145 |     void assemble_system_rhs_BFGS(const BlockVector<double> & solution_old,
1146 | 				  BlockVector<double> & system_rhs);
1147 | 
1148 |     void assemble_system_rhs_BFGS_parallel(const BlockVector<double> & solution_old,
1149 |     				           BlockVector<double> & system_rhs);
1150 | 
1151 |     void solve_nonlinear_timestep_LBFGS(BlockVector<double> &solution_delta,
1152 | 				        BlockVector<double> & LBFGS_update_refine);
1153 | 
1154 |     void solve_nonlinear_timestep_LBFGS_B(BlockVector<double> &solution_delta,
1155 |     				          BlockVector<double> & LBFGS_update_refine);
1156 | 
1157 |     void calculate_cauchy_point(const BlockSparseMatrix<double> & B0_matrix,
1158 | 	                        const std::list<BlockVector<double>> & y_vector_list,
1159 | 				const std::list<BlockVector<double>> & b0xs_vector_list,
1160 | 				const FullMatrix<double> & M_matrix,
1161 | 				const BlockVector<double> & gradient_g,
1162 | 				const BlockVector<double> & solution_delta,
1163 | 				BlockVector<double> & solution_delta_cauchy_point);
1164 | 
1165 |     double line_search_stepsize_gradient_based(const BlockVector<double> & BFGS_p_vector,
1166 | 					       const BlockVector<double> & solution_delta);
1167 | 
1168 |     double line_search_stepsize_strong_wolfe(const double phi_0,
1169 | 				             const double phi_0_prime,
1170 | 				             const BlockVector<double> & BFGS_p_vector,
1171 | 				             const BlockVector<double> & solution_delta);
1172 | 
1173 |     double line_search_zoom_strong_wolfe(double phi_low, double phi_low_prime, double alpha_low,
1174 | 					 double phi_high, double phi_high_prime, double alpha_high,
1175 | 					 double phi_0, double phi_0_prime, const BlockVector<double> & BFGS_p_vector,
1176 | 					 double c1, double c2, unsigned int max_iter,
1177 | 					 const BlockVector<double> & solution_delta);
1178 | 
1179 |     double line_search_interpolation_cubic(const double alpha_0, const double phi_0, const double phi_0_prime,
1180 | 					   const double alpha_1, const double phi_1, const double phi_1_prime);
1181 | 
1182 |     std::pair<double, double> calculate_phi_and_phi_prime(const double alpha,
1183 | 							  const BlockVector<double> & BFGS_p_vector,
1184 | 							  const BlockVector<double> & solution_delta);
1185 | 
1186 |     void LBFGS_B0(BlockVector<double> & LBFGS_r_vector,
1187 | 		  BlockVector<double> & LBFGS_q_vector);
1188 | 
1189 |     void output_results() const;
1190 | 
1191 |     void setup_qph();
1192 | 
1193 |     void update_qph_incremental(const BlockVector<double> &solution_delta,
1194 | 				const BlockVector<double> &solution_old);
1195 | 
1196 |     void update_qph_incremental_one_cell(
1197 |       const typename DoFHandler<dim>::active_cell_iterator &cell,
1198 |       ScratchData_UQPH &                                    scratch,
1199 |       PerTaskData_UQPH &                                    data);
1200 | 
1201 |     void copy_local_to_global_UQPH(const PerTaskData_UQPH & /*data*/)
1202 |     {}
1203 | 
1204 |     BlockVector<double>
1205 |     get_total_solution(const BlockVector<double> &solution_delta) const;
1206 | 
1207 |     // Should not make this function const
1208 |     void read_material_data(const std::string &data_file,
1209 | 			    const unsigned int total_material_regions);
1210 | 
1211 |     void read_time_data(const std::string &data_file,
1212 |     		        std::vector<std::array<double, 4>> & time_table);
1213 | 
1214 |     void print_conv_header_LBFGS();
1215 | 
1216 |     void print_conv_header_LBFGSB();
1217 | 
1218 |     void print_parameter_information();
1219 | 
1220 |     void calculate_reaction_force(unsigned int face_ID);
1221 | 
1222 |     void write_history_data();
1223 | 
1224 |     double calculate_energy_functional() const;
1225 | 
1226 |     std::pair<double, double> calculate_total_strain_energy_and_crack_energy_dissipation() const;
1227 | 
1228 |     bool local_refine_and_solution_transfer(BlockVector<double> & solution_delta,
1229 | 					    BlockVector<double> & LBFGS_update_refine);
1230 | 
1231 |     // L-BFGS-B subroutines
1232 |     void point_projection(BlockVector<double> & solution_delta);
1233 | 
1234 |     std::priority_queue< std::pair<double, unsigned int>,
1235 |                          std::vector<std::pair<double, unsigned int>>,
1236 |     		         std::greater<std::pair<double, unsigned int>> >
1237 |       calculate_break_points(const BlockVector<double> & solution_delta,
1238 | 			     const BlockVector<double> & gradient_g,
1239 | 			     BlockVector<double> & gradient_d);
1240 | 
1241 |     double ebT_x_B0_x_v(const unsigned int b,
1242 | 			const BlockSparseMatrix<double> & B0_matrix,
1243 | 			const BlockVector<double> & v);
1244 | 
1245 |     void zT_x_vector(const BlockVector<double> & z,
1246 | 		     const BlockVector<double> & src_vector,
1247 | 		     BlockVector<double> & target_vector);
1248 | 
1249 |     void z_x_vector(const BlockVector<double> & z,
1250 | 		    const BlockVector<double> & src_vector,
1251 | 		    BlockVector<double> & target_vector);
1252 | 
1253 |     void zT_B0_z(const BlockVector<double> & z,
1254 | 		     BlockSparseMatrix<double> & B0_matrix);
1255 | 
1256 | 
1257 |   }; // class PhaseFieldSplitSolve
1258 | 
1259 |   template <int dim>
1260 |   void PhaseFieldMonolithicSolve<dim>::zT_B0_z(const BlockVector<double> & z,
1261 | 					       BlockSparseMatrix<double> & B0_matrix)
1262 |   {
1263 |     // block 1: displacement
1264 |     for (unsigned int i = 0; i < z.block(m_u_dof).size(); ++i)
1265 |       {
1266 | 	if (z.block(m_u_dof)[i] < 0)
1267 | 	  {
1268 | 	    for (auto itr = B0_matrix.block(m_u_dof, m_u_dof).begin(i);
1269 | 		      itr != B0_matrix.block(m_u_dof, m_u_dof).end(i);
1270 | 		      ++itr)
1271 | 	      {
1272 | 		if (itr->column() != itr->row())
1273 | 		  {
1274 | 		    B0_matrix.block(m_u_dof, m_u_dof).set(itr->row(),    itr->column(), 0.0);
1275 | 		    B0_matrix.block(m_u_dof, m_u_dof).set(itr->column(), itr->row(),    0.0);
1276 | 		  }
1277 | 	      }
1278 | 	  }
1279 |       } // for (unsigned int i = 0; i < z.block(m_u_dof).size(); ++i)
1280 | 
1281 |     // block 2: phasefield
1282 |     for (unsigned int i = 0; i < z.block(m_d_dof).size(); ++i)
1283 |       {
1284 | 	if (z.block(m_d_dof)[i] < 0)
1285 | 	  {
1286 | 	    for (auto itr = B0_matrix.block(m_d_dof, m_d_dof).begin(i);
1287 | 		      itr != B0_matrix.block(m_d_dof, m_d_dof).end(i);
1288 | 		      ++itr)
1289 | 	      {
1290 | 		if (itr->column() != itr->row())
1291 | 		  {
1292 | 		    B0_matrix.block(m_d_dof, m_d_dof).set(itr->row(),    itr->column(), 0.0);
1293 | 		    B0_matrix.block(m_d_dof, m_d_dof).set(itr->column(), itr->row(),    0.0);
1294 | 		  }
1295 | 	      }
1296 | 	  }
1297 |       } // for (unsigned int i = 0; i < z.block(m_d_dof).size(); ++i)
1298 |   }
1299 | 
1300 |   template <int dim>
1301 |   void PhaseFieldMonolithicSolve<dim>::z_x_vector(const BlockVector<double> & z,
1302 | 						  const BlockVector<double> & src_vector,
1303 | 						  BlockVector<double> & target_vector)
1304 |   {
1305 |     //We assume that the dimensions of all the block matrices are correct
1306 |     //block 1: displacement
1307 |     unsigned int target_vector_index = 0;
1308 |     for (unsigned int i = 0; i < z.block(m_u_dof).size(); ++i)
1309 |       {
1310 | 	if (z.block(m_u_dof)[i] > 0)
1311 | 	  {
1312 | 	    target_vector.block(m_u_dof)[i] = src_vector.block(m_u_dof)[target_vector_index];
1313 | 	    ++target_vector_index;
1314 | 	  }
1315 | 	else
1316 | 	  target_vector.block(m_u_dof)[i] = 0;
1317 |       }
1318 | 
1319 |     //block 2: phasefield
1320 |     target_vector_index = 0;
1321 |     for (unsigned int i = 0; i < z.block(m_d_dof).size(); ++i)
1322 |       {
1323 | 	if (z.block(m_d_dof)[i] > 0)
1324 | 	  {
1325 | 	    target_vector.block(m_d_dof)[i] = src_vector.block(m_d_dof)[target_vector_index];
1326 | 	    ++target_vector_index;
1327 | 	  }
1328 | 	else
1329 | 	  target_vector.block(m_d_dof)[i] = 0;
1330 |       }
1331 |   }
1332 | 
1333 |   template <int dim>
1334 |   void PhaseFieldMonolithicSolve<dim>::zT_x_vector(const BlockVector<double> & z,
1335 | 						   const BlockVector<double> & src_vector,
1336 | 						   BlockVector<double> & target_vector)
1337 |   {
1338 |     //We assume that the dimensions of all the block matrices are correct
1339 |     //block 1: displacement
1340 |     unsigned int target_vector_index = 0;
1341 |     for (unsigned int i = 0; i < z.block(m_u_dof).size(); ++i)
1342 |       {
1343 | 	if (z.block(m_u_dof)[i] > 0)
1344 | 	  {
1345 | 	    target_vector.block(m_u_dof)[target_vector_index] = src_vector.block(m_u_dof)[i];
1346 | 	    ++target_vector_index;
1347 | 	  }
1348 |       }
1349 | 
1350 |     //block 2: phasefield
1351 |     target_vector_index = 0;
1352 |     for (unsigned int i = 0; i < z.block(m_d_dof).size(); ++i)
1353 |       {
1354 |         if (z.block(m_d_dof)[i] > 0)
1355 |         {
1356 | 	  target_vector.block(m_d_dof)[target_vector_index] = src_vector.block(m_d_dof)[i];
1357 | 	  ++target_vector_index;
1358 |         }
1359 |       }
1360 |   }
1361 | 
1362 |   template <int dim>
1363 |   double PhaseFieldMonolithicSolve<dim>::ebT_x_B0_x_v(const unsigned int b,
1364 | 						      const BlockSparseMatrix<double> & B0_matrix,
1365 | 						      const BlockVector<double> & v)
1366 |   {
1367 |     double row_sum = 0.0;
1368 |     for (auto itr = B0_matrix.block(m_d_dof, m_d_dof).begin(b);
1369 |               itr != B0_matrix.block(m_d_dof, m_d_dof).end(b);
1370 |               ++itr)
1371 |       {
1372 |         row_sum += itr->value() * v.block(m_d_dof)(itr->column());
1373 |       }
1374 | 
1375 |     return row_sum;
1376 |   }
1377 | 
1378 |   template <int dim>
1379 |   void PhaseFieldMonolithicSolve<dim>::point_projection(BlockVector<double> & solution_delta)
1380 |   {
1381 |     // Phase-field value cannot exceed 1.0
1382 |     const double upper_limit = 1.0;
1383 | 
1384 |     Vector<double> solution_phasefield_total(m_solution.block(m_d_dof));
1385 |     solution_phasefield_total += solution_delta.block(m_d_dof);
1386 | 
1387 |     for (unsigned int i = 0; i < solution_phasefield_total.size(); ++i)
1388 |       {
1389 | 	if (solution_delta.block(m_d_dof)[i] < 0.0)
1390 | 	  solution_delta.block(m_d_dof)[i] = 0.0;
1391 | 
1392 | 	if (solution_phasefield_total[i] > upper_limit)
1393 | 	  solution_delta.block(m_d_dof)[i] = upper_limit - m_solution.block(m_d_dof)[i];
1394 |       }
1395 |   }
1396 | 
1397 |   template <int dim>
1398 |   std::priority_queue< std::pair<double, unsigned int>,
1399 |                        std::vector<std::pair<double, unsigned int>>,
1400 |   		       std::greater<std::pair<double, unsigned int>> >
1401 |     PhaseFieldMonolithicSolve<dim>::calculate_break_points(const BlockVector<double> & solution_delta,
1402 | 		       				           const BlockVector<double> & gradient_g,
1403 | 							   BlockVector<double> & gradient_d)
1404 |   {
1405 |     // Creates a min heap of break points
1406 |     std::priority_queue< std::pair<double, unsigned int>,
1407 |                          std::vector<std::pair<double, unsigned int>>,
1408 |     		         std::greater<std::pair<double, unsigned int>> >
1409 |     break_points_sorted;
1410 | 
1411 |     double t = 0.0;
1412 | 
1413 |     Vector<double> solution_phasefield_total(m_solution.block(m_d_dof));
1414 |     solution_phasefield_total += solution_delta.block(m_d_dof);
1415 | 
1416 |     // upper bound is 1.0, lower bound is the solution at the previous step.
1417 |     for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
1418 |       {
1419 | 	if (gradient_g.block(m_d_dof)[i] < 0)
1420 | 	  t = (solution_phasefield_total[i] - 1.0 ) / gradient_g.block(m_d_dof)[i];
1421 | 	else if (gradient_g.block(m_d_dof)[i] > 0)
1422 | 	  t = solution_delta.block(m_d_dof)[i] / gradient_g.block(m_d_dof)[i];
1423 | 	else
1424 | 	  t = std::numeric_limits<double>::max();
1425 | 
1426 |         //AssertThrow(t >= 0, ExcMessage("Break point has to be a non-negative t value"));
1427 | 
1428 |         if (t > 0)
1429 |           {
1430 | 	    break_points_sorted.push(std::make_pair(t, i));
1431 |           }
1432 |         else // if t == 0, i is in the active set
1433 |           {
1434 |             gradient_d.block(m_d_dof)[i] = 0;
1435 |             if (gradient_g.block(m_d_dof)[i] > 0)
1436 |               m_active_set_phasefield(i) = 1; //lower bound
1437 |             else
1438 |               m_active_set_phasefield(i) = 2; //upper bound
1439 |           }
1440 |       }
1441 | 
1442 |     return break_points_sorted;
1443 |   }
1444 | 
1445 |   template <int dim>
1446 |   void PhaseFieldMonolithicSolve<dim>::get_error_residual(Errors &error_residual)
1447 |   {
1448 |     BlockVector<double> error_res(m_dofs_per_block);
1449 | 
1450 |     for (unsigned int i = 0; i < m_dof_handler.n_dofs(); ++i)
1451 |       if (!m_constraints.is_constrained(i))
1452 |         error_res(i) = m_system_rhs(i);
1453 | 
1454 |     error_residual.m_norm = error_res.l2_norm();
1455 |     error_residual.m_u    = error_res.block(m_u_dof).l2_norm();
1456 |     error_residual.m_d    = error_res.block(m_d_dof).l2_norm();
1457 |   }
1458 | 
1459 |   template <int dim>
1460 |   void PhaseFieldMonolithicSolve<dim>::get_error_residual_LBFGSB(Errors &error_residual,
1461 | 								 const BlockVector<double> & solution_delta)
1462 |   {
1463 |     // We use L_2 norm
1464 |     BlockVector<double> error_res(m_dofs_per_block);
1465 | 
1466 |     // For displacement DOFs, except essential boundary conditions
1467 |     // and hanging-node constraints, there are no box constraints
1468 |     for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
1469 |       {
1470 |         if (!m_constraints.is_constrained(i))
1471 | 	  error_res.block(m_u_dof)[i] = m_system_rhs.block(m_u_dof)[i];
1472 |       }
1473 | 
1474 |     // For phasefield DOFs, there are points with active box constraints
1475 |     // and points with inactive box constraints
1476 |     const double upper_limit = 1.0;
1477 |     Vector<double> solution_phasefield_total(m_solution.block(m_d_dof));
1478 |     solution_phasefield_total += solution_delta.block(m_d_dof);
1479 | 
1480 |     double trial_solution = 0.0;
1481 |     for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
1482 |       {
1483 | 	// phasefield DOFs can still be constrained due to hanging-nodes
1484 |         if (!m_constraints.is_constrained(i + m_dofs_per_block[m_u_dof]))
1485 |           {
1486 |             trial_solution = solution_phasefield_total(i) - m_system_rhs.block(m_d_dof)[i];
1487 | 
1488 |             if (trial_solution < m_solution.block(m_d_dof)[i])
1489 |               error_res.block(m_d_dof)[i] = m_solution.block(m_d_dof)[i] - solution_phasefield_total(i);
1490 |             else if (trial_solution > upper_limit)
1491 |               error_res.block(m_d_dof)[i] = upper_limit - solution_phasefield_total(i);
1492 |             else
1493 |               error_res.block(m_d_dof)[i] = (-m_system_rhs.block(m_d_dof)[i]);
1494 |           }
1495 |       }
1496 | 
1497 |     error_residual.m_norm = error_res.l2_norm();
1498 |     error_residual.m_u    = error_res.block(m_u_dof).l2_norm();
1499 |     error_residual.m_d    = error_res.block(m_d_dof).l2_norm();
1500 |   }
1501 | 
1502 |   template <int dim>
1503 |   void PhaseFieldMonolithicSolve<dim>::get_error_update(const BlockVector<double> &newton_update,
1504 | 							Errors & error_update)
1505 |   {
1506 |     BlockVector<double> error_ud(m_dofs_per_block);
1507 |     for (unsigned int i = 0; i < m_dof_handler.n_dofs(); ++i)
1508 |       if (!m_constraints.is_constrained(i))
1509 | 	error_ud(i) = newton_update(i);
1510 | 
1511 |     error_update.m_norm = error_ud.l2_norm();
1512 |     error_update.m_u    = error_ud.block(m_u_dof).l2_norm();
1513 |     error_update.m_d    = error_ud.block(m_d_dof).l2_norm();
1514 |   }
1515 | 
1516 |   template <int dim>
1517 |   void PhaseFieldMonolithicSolve<dim>::read_material_data(const std::string &data_file,
1518 | 				                     const unsigned int total_material_regions)
1519 |   {
1520 |     std::ifstream myfile (data_file);
1521 | 
1522 |     double lame_lambda, lame_mu, length_scale, gc, viscosity, residual_k;
1523 |     int material_region;
1524 |     double poisson_ratio;
1525 |     if (myfile.is_open())
1526 |       {
1527 |         m_logfile << "Reading material data file ..." << std::endl;
1528 | 
1529 |         while ( myfile >> material_region
1530 |                        >> lame_lambda
1531 | 		       >> lame_mu
1532 | 		       >> length_scale
1533 | 		       >> gc
1534 | 		       >> viscosity
1535 | 		       >> residual_k)
1536 |           {
1537 |             m_material_data[material_region] = {lame_lambda,
1538 |         	                                lame_mu,
1539 | 						length_scale,
1540 | 						gc,
1541 | 						viscosity,
1542 |                                                 residual_k};
1543 |             poisson_ratio = lame_lambda / (2*(lame_lambda + lame_mu));
1544 |             Assert( (poisson_ratio <= 0.5)&(poisson_ratio >=-1.0) , ExcInternalError());
1545 | 
1546 |             m_logfile << "\tRegion " << material_region << " : " << std::endl;
1547 |             m_logfile << "\t\tLame lambda = " << lame_lambda << std::endl;
1548 |             m_logfile << "\t\tLame mu = "  << lame_mu << std::endl;
1549 |             m_logfile << "\t\tPoisson ratio = "  << poisson_ratio << std::endl;
1550 |             m_logfile << "\t\tPhase field length scale (l) = " << length_scale << std::endl;
1551 |             m_logfile << "\t\tCritical energy release rate (gc) = "  << gc << std::endl;
1552 |             m_logfile << "\t\tViscosity for regularization (eta) = "  << viscosity << std::endl;
1553 |             m_logfile << "\t\tResidual_k (k) = "  << residual_k << std::endl;
1554 |           }
1555 | 
1556 |         if (m_material_data.size() != total_material_regions)
1557 |           {
1558 |             m_logfile << "Material data file has " << m_material_data.size() << " rows. However, "
1559 |         	      << "the mesh has " << total_material_regions << " material regions."
1560 | 		      << std::endl;
1561 |             Assert(m_material_data.size() == total_material_regions,
1562 |                        ExcDimensionMismatch(m_material_data.size(), total_material_regions));
1563 |           }
1564 |         myfile.close();
1565 |       }
1566 |     else
1567 |       {
1568 | 	m_logfile << "Material data file : " << data_file << " not exist!" << std::endl;
1569 | 	Assert(false, ExcMessage("Failed to read material data file"));
1570 |       }
1571 |   }
1572 | 
1573 |   template <int dim>
1574 |   void PhaseFieldMonolithicSolve<dim>::read_time_data(const std::string &data_file,
1575 | 				                 std::vector<std::array<double, 4>> & time_table)
1576 |   {
1577 |     std::ifstream myfile (data_file);
1578 | 
1579 |     double t_0, t_1, delta_t, t_magnitude;
1580 | 
1581 |     if (myfile.is_open())
1582 |       {
1583 | 	m_logfile << "Reading time data file ..." << std::endl;
1584 | 
1585 | 	while ( myfile >> t_0
1586 | 		       >> t_1
1587 | 		       >> delta_t
1588 | 		       >> t_magnitude)
1589 | 	  {
1590 | 	    Assert( t_0 < t_1,
1591 | 		    ExcMessage("For each time pair, "
1592 | 			       "the start time should be smaller than the end time"));
1593 | 	    time_table.push_back({{t_0, t_1, delta_t, t_magnitude}});
1594 | 	  }
1595 | 
1596 | 	Assert(std::fabs(t_1 - m_parameters.m_end_time) < 1.0e-9,
1597 | 	       ExcMessage("End time in time table is inconsistent with input data in parameters.prm"));
1598 | 
1599 | 	Assert(time_table.size() > 0,
1600 | 	       ExcMessage("Time data file is empty."));
1601 | 	myfile.close();
1602 |       }
1603 |     else
1604 |       {
1605 |         m_logfile << "Time data file : " << data_file << " not exist!" << std::endl;
1606 |         Assert(false, ExcMessage("Failed to read time data file"));
1607 |       }
1608 | 
1609 |     for (auto & time_group : time_table)
1610 |       {
1611 | 	m_logfile << "\t\t"
1612 | 	          << time_group[0] << ",\t"
1613 | 	          << time_group[1] << ",\t"
1614 | 		  << time_group[2] << ",\t"
1615 | 		  << time_group[3] << std::endl;
1616 |       }
1617 |   }
1618 | 
1619 |   template <int dim>
1620 |   void PhaseFieldMonolithicSolve<dim>::setup_qph()
1621 |   {
1622 |     m_logfile << "\t\tSetting up quadrature point data ("
1623 | 	      << m_n_q_points
1624 | 	      << " points per cell)" << std::endl;
1625 | 
1626 |     m_quadrature_point_history.clear();
1627 |     for (auto const & cell : m_triangulation.active_cell_iterators())
1628 |       {
1629 | 	m_quadrature_point_history.initialize(cell, m_n_q_points);
1630 |       }
1631 | 
1632 |     unsigned int material_id;
1633 |     double lame_lambda = 0.0;
1634 |     double lame_mu = 0.0;
1635 |     double length_scale = 0.0;
1636 |     double gc = 0.0;
1637 |     double viscosity = 0.0;
1638 |     double residual_k = 0.0;
1639 | 
1640 |     for (const auto &cell : m_triangulation.active_cell_iterators())
1641 |       {
1642 |         material_id = cell->material_id();
1643 |         if (m_material_data.find(material_id) != m_material_data.end())
1644 |           {
1645 |             lame_lambda                = m_material_data[material_id][0];
1646 |             lame_mu                    = m_material_data[material_id][1];
1647 |             length_scale               = m_material_data[material_id][2];
1648 |             gc                         = m_material_data[material_id][3];
1649 |             viscosity                  = m_material_data[material_id][4];
1650 |             residual_k                 = m_material_data[material_id][5];
1651 | 	  }
1652 |         else
1653 |           {
1654 |             m_logfile << "Could not find material data for material id: " << material_id << std::endl;
1655 |             AssertThrow(false, ExcMessage("Could not find material data for material id."));
1656 |           }
1657 | 
1658 |         const std::vector<std::shared_ptr<PointHistory<dim>>> lqph =
1659 |           m_quadrature_point_history.get_data(cell);
1660 |         Assert(lqph.size() == m_n_q_points, ExcInternalError());
1661 | 
1662 |         for (unsigned int q_point = 0; q_point < m_n_q_points; ++q_point)
1663 |           lqph[q_point]->setup_lqp(lame_lambda, lame_mu, length_scale,
1664 | 				   gc, viscosity, residual_k,
1665 | 				   m_parameters.m_plane_stress);
1666 |       }
1667 |   }
1668 | 
1669 |   template <int dim>
1670 |   BlockVector<double> PhaseFieldMonolithicSolve<dim>::get_total_solution(
1671 |     const BlockVector<double> &solution_delta) const
1672 |   {
1673 |     BlockVector<double> solution_total(m_solution);
1674 |     solution_total += solution_delta;
1675 |     return solution_total;
1676 |   }
1677 | 
1678 |   template <int dim>
1679 |   void
1680 |   PhaseFieldMonolithicSolve<dim>::update_qph_incremental(const BlockVector<double> &solution_delta,
1681 | 							 const BlockVector<double> &solution_old)
1682 |   {
1683 |     m_timer.enter_subsection("Update QPH data");
1684 | 
1685 |     const BlockVector<double> solution_total(get_total_solution(solution_delta));
1686 | 
1687 |     const UpdateFlags uf_UQPH(update_values | update_gradients);
1688 |     PerTaskData_UQPH  per_task_data_UQPH;
1689 |     ScratchData_UQPH  scratch_data_UQPH(m_fe,
1690 | 					m_qf_cell,
1691 | 					uf_UQPH,
1692 | 					solution_total,
1693 | 					solution_old,
1694 | 					m_time.get_delta_t());
1695 | 
1696 |     auto worker = [this](const typename DoFHandler<dim>::active_cell_iterator &cell,
1697 | 	                 ScratchData_UQPH & scratch,
1698 | 	                 PerTaskData_UQPH & data)
1699 |       {
1700 |         this->update_qph_incremental_one_cell(cell, scratch, data);
1701 |       };
1702 | 
1703 |     auto copier = [this](const PerTaskData_UQPH &data)
1704 |       {
1705 |         this->copy_local_to_global_UQPH(data);
1706 |       };
1707 | 
1708 |     WorkStream::run(
1709 | 	m_dof_handler.begin_active(),
1710 | 	m_dof_handler.end(),
1711 | 	worker,
1712 | 	copier,
1713 | 	scratch_data_UQPH,
1714 | 	per_task_data_UQPH);
1715 | 
1716 |     m_timer.leave_subsection();
1717 |   }
1718 | 
1719 |   template <int dim>
1720 |   struct PhaseFieldMonolithicSolve<dim>::PerTaskData_UQPH
1721 |   {
1722 |     void reset()
1723 |     {}
1724 |   };
1725 | 
1726 |   template <int dim>
1727 |   struct PhaseFieldMonolithicSolve<dim>::ScratchData_UQPH
1728 |   {
1729 |     const BlockVector<double> & m_solution_UQPH;
1730 | 
1731 |     std::vector<SymmetricTensor<2, dim>> m_solution_symm_grads_u_cell;
1732 |     std::vector<double>         m_solution_values_phasefield_cell;
1733 |     std::vector<Tensor<1, dim>> m_solution_grad_phasefield_cell;
1734 | 
1735 |     FEValues<dim> m_fe_values;
1736 | 
1737 |     const BlockVector<double>&       m_solution_previous_step;
1738 |     std::vector<double>              m_phasefield_previous_step_cell;
1739 | 
1740 |     const double                     m_delta_time;
1741 | 
1742 |     ScratchData_UQPH(const FiniteElement<dim> & fe_cell,
1743 |                      const QGauss<dim> &        qf_cell,
1744 |                      const UpdateFlags          uf_cell,
1745 |                      const BlockVector<double> &solution_total,
1746 | 		     const BlockVector<double> &solution_old,
1747 | 		     const double delta_time)
1748 |       : m_solution_UQPH(solution_total)
1749 |       , m_solution_symm_grads_u_cell(qf_cell.size())
1750 |       , m_solution_values_phasefield_cell(qf_cell.size())
1751 |       , m_solution_grad_phasefield_cell(qf_cell.size())
1752 |       , m_fe_values(fe_cell, qf_cell, uf_cell)
1753 |       , m_solution_previous_step(solution_old)
1754 |       , m_phasefield_previous_step_cell(qf_cell.size())
1755 |       , m_delta_time(delta_time)
1756 |     {}
1757 | 
1758 |     ScratchData_UQPH(const ScratchData_UQPH &rhs)
1759 |       : m_solution_UQPH(rhs.m_solution_UQPH)
1760 |       , m_solution_symm_grads_u_cell(rhs.m_solution_symm_grads_u_cell)
1761 |       , m_solution_values_phasefield_cell(rhs.m_solution_values_phasefield_cell)
1762 |       , m_solution_grad_phasefield_cell(rhs.m_solution_grad_phasefield_cell)
1763 |       , m_fe_values(rhs.m_fe_values.get_fe(),
1764 |                     rhs.m_fe_values.get_quadrature(),
1765 |                     rhs.m_fe_values.get_update_flags())
1766 |       , m_solution_previous_step(rhs.m_solution_previous_step)
1767 |       , m_phasefield_previous_step_cell(rhs.m_phasefield_previous_step_cell)
1768 |       , m_delta_time(rhs.m_delta_time)
1769 |     {}
1770 | 
1771 |     void reset()
1772 |     {
1773 |       const unsigned int n_q_points = m_solution_symm_grads_u_cell.size();
1774 |       for (unsigned int q = 0; q < n_q_points; ++q)
1775 |         {
1776 |           m_solution_symm_grads_u_cell[q]  = 0.0;
1777 |           m_solution_values_phasefield_cell[q] = 0.0;
1778 |           m_solution_grad_phasefield_cell[q] = 0.0;
1779 |           m_phasefield_previous_step_cell[q] = 0.0;
1780 |         }
1781 |     }
1782 |   };
1783 | 
1784 |   template <int dim>
1785 |   void PhaseFieldMonolithicSolve<dim>::update_qph_incremental_one_cell(
1786 |     const typename DoFHandler<dim>::active_cell_iterator &cell,
1787 |     ScratchData_UQPH & scratch,
1788 |     PerTaskData_UQPH & /*data*/)
1789 |   {
1790 |     scratch.reset();
1791 | 
1792 |     scratch.m_fe_values.reinit(cell);
1793 | 
1794 |     const std::vector<std::shared_ptr<PointHistory<dim>>> lqph =
1795 |       m_quadrature_point_history.get_data(cell);
1796 |     Assert(lqph.size() == m_n_q_points, ExcInternalError());
1797 | 
1798 |     const FEValuesExtractors::Vector displacement(0);
1799 | 
1800 |     scratch.m_fe_values[m_u_fe].get_function_symmetric_gradients(
1801 |       scratch.m_solution_UQPH, scratch.m_solution_symm_grads_u_cell);
1802 |     scratch.m_fe_values[m_d_fe].get_function_values(
1803 |       scratch.m_solution_UQPH, scratch.m_solution_values_phasefield_cell);
1804 |     scratch.m_fe_values[m_d_fe].get_function_gradients(
1805 |       scratch.m_solution_UQPH, scratch.m_solution_grad_phasefield_cell);
1806 | 
1807 |     scratch.m_fe_values[m_d_fe].get_function_values(
1808 |       scratch.m_solution_previous_step, scratch.m_phasefield_previous_step_cell);
1809 | 
1810 |     for (const unsigned int q_point :
1811 |          scratch.m_fe_values.quadrature_point_indices())
1812 |       lqph[q_point]->update_field_values(scratch.m_solution_symm_grads_u_cell[q_point],
1813 |                                          scratch.m_solution_values_phasefield_cell[q_point],
1814 | 					 scratch.m_solution_grad_phasefield_cell[q_point],
1815 | 					 scratch.m_phasefield_previous_step_cell[q_point],
1816 | 					 scratch.m_delta_time);
1817 |   }
1818 | 
1819 |   template <int dim>
1820 |   struct PhaseFieldMonolithicSolve<dim>::PerTaskData_ASM
1821 |   {
1822 |     FullMatrix<double>                   m_cell_matrix;
1823 |     Vector<double>                       m_cell_rhs;
1824 |     std::vector<types::global_dof_index> m_local_dof_indices;
1825 | 
1826 |     PerTaskData_ASM(const unsigned int dofs_per_cell)
1827 |       : m_cell_matrix(dofs_per_cell, dofs_per_cell)
1828 |       , m_cell_rhs(dofs_per_cell)
1829 |       , m_local_dof_indices(dofs_per_cell)
1830 |     {}
1831 | 
1832 |     void reset()
1833 |     {
1834 |       m_cell_matrix = 0.0;
1835 |       m_cell_rhs    = 0.0;
1836 |     }
1837 |   };
1838 | 
1839 |   template <int dim>
1840 |   struct PhaseFieldMonolithicSolve<dim>::PerTaskData_ASM_RHS_BFGS
1841 |   {
1842 |     Vector<double>                       m_cell_rhs;
1843 |     std::vector<types::global_dof_index> m_local_dof_indices;
1844 | 
1845 |     PerTaskData_ASM_RHS_BFGS(const unsigned int dofs_per_cell)
1846 |       : m_cell_rhs(dofs_per_cell)
1847 |       , m_local_dof_indices(dofs_per_cell)
1848 |     {}
1849 | 
1850 |     void reset()
1851 |     {
1852 |       m_cell_rhs    = 0.0;
1853 |     }
1854 |   };
1855 | 
1856 |   template <int dim>
1857 |   struct PhaseFieldMonolithicSolve<dim>::ScratchData_ASM
1858 |   {
1859 |     FEValues<dim>     m_fe_values;
1860 |     FEFaceValues<dim> m_fe_face_values;
1861 | 
1862 |     std::vector<std::vector<double>>                  m_Nx_phasefield;      // shape function values for phase-field
1863 |     std::vector<std::vector<Tensor<1, dim>>>          m_grad_Nx_phasefield; // gradient of shape function values for phase field
1864 | 
1865 |     std::vector<std::vector<Tensor<1, dim>>>          m_Nx_disp;       // shape function values for displacement
1866 |     std::vector<std::vector<Tensor<2, dim>>>          m_grad_Nx_disp;  // gradient of shape function values for displacement
1867 |     std::vector<std::vector<SymmetricTensor<2, dim>>> m_symm_grad_Nx_disp;  // symmetric gradient of shape function values for displacement
1868 | 
1869 |     const BlockVector<double>&       m_solution_previous_step;
1870 |     std::vector<double>              m_phasefield_previous_step_cell;
1871 | 
1872 |     ScratchData_ASM(const FiniteElement<dim> & fe_cell,
1873 |                     const QGauss<dim> &        qf_cell,
1874 |                     const UpdateFlags          uf_cell,
1875 | 		    const QGauss<dim - 1> &    qf_face,
1876 | 		    const UpdateFlags          uf_face,
1877 | 		    const BlockVector<double>& solution_old)
1878 |       : m_fe_values(fe_cell, qf_cell, uf_cell)
1879 |       , m_fe_face_values(fe_cell, qf_face, uf_face)
1880 |       , m_Nx_phasefield(qf_cell.size(),
1881 | 	                std::vector<double>(fe_cell.n_dofs_per_cell()))
1882 |       , m_grad_Nx_phasefield(qf_cell.size(),
1883 | 		             std::vector<Tensor<1, dim>>(fe_cell.n_dofs_per_cell()))
1884 |       , m_Nx_disp(qf_cell.size(),
1885 | 		  std::vector<Tensor<1, dim>>(fe_cell.n_dofs_per_cell()))
1886 |       , m_grad_Nx_disp(qf_cell.size(),
1887 |                        std::vector<Tensor<2, dim>>(fe_cell.n_dofs_per_cell()))
1888 |       , m_symm_grad_Nx_disp(qf_cell.size(),
1889 |                             std::vector<SymmetricTensor<2, dim>>(fe_cell.n_dofs_per_cell()))
1890 |       , m_solution_previous_step(solution_old)
1891 |       , m_phasefield_previous_step_cell(qf_cell.size())
1892 |     {}
1893 | 
1894 |     ScratchData_ASM(const ScratchData_ASM &rhs)
1895 |       : m_fe_values(rhs.m_fe_values.get_fe(),
1896 |                     rhs.m_fe_values.get_quadrature(),
1897 |                     rhs.m_fe_values.get_update_flags())
1898 |       , m_fe_face_values(rhs.m_fe_face_values.get_fe(),
1899 | 	                 rhs.m_fe_face_values.get_quadrature(),
1900 | 	                 rhs.m_fe_face_values.get_update_flags())
1901 |       , m_Nx_phasefield(rhs.m_Nx_phasefield)
1902 |       , m_grad_Nx_phasefield(rhs.m_grad_Nx_phasefield)
1903 |       , m_Nx_disp(rhs.m_Nx_disp)
1904 |       , m_grad_Nx_disp(rhs.m_grad_Nx_disp)
1905 |       , m_symm_grad_Nx_disp(rhs.m_symm_grad_Nx_disp)
1906 |       , m_solution_previous_step(rhs.m_solution_previous_step)
1907 |       , m_phasefield_previous_step_cell(rhs.m_phasefield_previous_step_cell)
1908 |     {}
1909 | 
1910 |     void reset()
1911 |     {
1912 |       const unsigned int n_q_points      = m_Nx_phasefield.size();
1913 |       const unsigned int n_dofs_per_cell = m_Nx_phasefield[0].size();
1914 |       for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
1915 |         {
1916 |           Assert(m_Nx_phasefield[q_point].size() == n_dofs_per_cell,
1917 | 		 ExcInternalError());
1918 | 
1919 |           Assert(m_grad_Nx_phasefield[q_point].size() == n_dofs_per_cell,
1920 |                  ExcInternalError());
1921 | 
1922 |           Assert(m_Nx_disp[q_point].size() == n_dofs_per_cell,
1923 | 		 ExcInternalError());
1924 | 
1925 |           Assert(m_grad_Nx_disp[q_point].size() == n_dofs_per_cell,
1926 |                  ExcInternalError());
1927 | 
1928 |           Assert(m_symm_grad_Nx_disp[q_point].size() == n_dofs_per_cell,
1929 |                  ExcInternalError());
1930 | 
1931 |           m_phasefield_previous_step_cell[q_point] = 0.0;
1932 |           for (unsigned int k = 0; k < n_dofs_per_cell; ++k)
1933 |             {
1934 |               m_Nx_phasefield[q_point][k]           = 0.0;
1935 |               m_grad_Nx_phasefield[q_point][k]      = 0.0;
1936 |               m_Nx_disp[q_point][k]                 = 0.0;
1937 |               m_grad_Nx_disp[q_point][k]            = 0.0;
1938 |               m_symm_grad_Nx_disp[q_point][k]       = 0.0;
1939 |             }
1940 |         }
1941 |     }
1942 |   };
1943 | 
1944 | 
1945 |   template <int dim>
1946 |   struct PhaseFieldMonolithicSolve<dim>::ScratchData_ASM_RHS_BFGS
1947 |   {
1948 |     FEValues<dim>     m_fe_values;
1949 |     FEFaceValues<dim> m_fe_face_values;
1950 | 
1951 |     std::vector<std::vector<double>>                  m_Nx_phasefield;      // shape function values for phase-field
1952 |     std::vector<std::vector<Tensor<1, dim>>>          m_grad_Nx_phasefield; // gradient of shape function values for phase field
1953 | 
1954 |     std::vector<std::vector<Tensor<1, dim>>>          m_Nx_disp;       // shape function values for displacement
1955 |     std::vector<std::vector<Tensor<2, dim>>>          m_grad_Nx_disp;  // gradient of shape function values for displacement
1956 |     std::vector<std::vector<SymmetricTensor<2, dim>>> m_symm_grad_Nx_disp;  // symmetric gradient of shape function values for displacement
1957 | 
1958 |     const BlockVector<double>&       m_solution_previous_step;
1959 |     std::vector<double>              m_phasefield_previous_step_cell;
1960 | 
1961 |     ScratchData_ASM_RHS_BFGS(const FiniteElement<dim> & fe_cell,
1962 |                              const QGauss<dim> &        qf_cell,
1963 |                              const UpdateFlags          uf_cell,
1964 | 		             const QGauss<dim - 1> &    qf_face,
1965 | 		             const UpdateFlags          uf_face,
1966 | 		             const BlockVector<double>& solution_old)
1967 |       : m_fe_values(fe_cell, qf_cell, uf_cell)
1968 |       , m_fe_face_values(fe_cell, qf_face, uf_face)
1969 |       , m_Nx_phasefield(qf_cell.size(),
1970 | 	                std::vector<double>(fe_cell.n_dofs_per_cell()))
1971 |       , m_grad_Nx_phasefield(qf_cell.size(),
1972 | 		             std::vector<Tensor<1, dim>>(fe_cell.n_dofs_per_cell()))
1973 |       , m_Nx_disp(qf_cell.size(),
1974 | 		  std::vector<Tensor<1, dim>>(fe_cell.n_dofs_per_cell()))
1975 |       , m_grad_Nx_disp(qf_cell.size(),
1976 |                        std::vector<Tensor<2, dim>>(fe_cell.n_dofs_per_cell()))
1977 |       , m_symm_grad_Nx_disp(qf_cell.size(),
1978 |                             std::vector<SymmetricTensor<2, dim>>(fe_cell.n_dofs_per_cell()))
1979 |       , m_solution_previous_step(solution_old)
1980 |       , m_phasefield_previous_step_cell(qf_cell.size())
1981 |     {}
1982 | 
1983 |     ScratchData_ASM_RHS_BFGS(const ScratchData_ASM_RHS_BFGS &rhs)
1984 |       : m_fe_values(rhs.m_fe_values.get_fe(),
1985 |                     rhs.m_fe_values.get_quadrature(),
1986 |                     rhs.m_fe_values.get_update_flags())
1987 |       , m_fe_face_values(rhs.m_fe_face_values.get_fe(),
1988 | 	                 rhs.m_fe_face_values.get_quadrature(),
1989 | 	                 rhs.m_fe_face_values.get_update_flags())
1990 |       , m_Nx_phasefield(rhs.m_Nx_phasefield)
1991 |       , m_grad_Nx_phasefield(rhs.m_grad_Nx_phasefield)
1992 |       , m_Nx_disp(rhs.m_Nx_disp)
1993 |       , m_grad_Nx_disp(rhs.m_grad_Nx_disp)
1994 |       , m_symm_grad_Nx_disp(rhs.m_symm_grad_Nx_disp)
1995 |       , m_solution_previous_step(rhs.m_solution_previous_step)
1996 |       , m_phasefield_previous_step_cell(rhs.m_phasefield_previous_step_cell)
1997 |     {}
1998 | 
1999 |     void reset()
2000 |     {
2001 |       const unsigned int n_q_points      = m_Nx_phasefield.size();
2002 |       const unsigned int n_dofs_per_cell = m_Nx_phasefield[0].size();
2003 |       for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
2004 |         {
2005 |           Assert(m_Nx_phasefield[q_point].size() == n_dofs_per_cell,
2006 | 		 ExcInternalError());
2007 | 
2008 |           Assert(m_grad_Nx_phasefield[q_point].size() == n_dofs_per_cell,
2009 |                  ExcInternalError());
2010 | 
2011 |           Assert(m_Nx_disp[q_point].size() == n_dofs_per_cell,
2012 | 		 ExcInternalError());
2013 | 
2014 |           Assert(m_grad_Nx_disp[q_point].size() == n_dofs_per_cell,
2015 |                  ExcInternalError());
2016 | 
2017 |           Assert(m_symm_grad_Nx_disp[q_point].size() == n_dofs_per_cell,
2018 |                  ExcInternalError());
2019 | 
2020 |           m_phasefield_previous_step_cell[q_point] = 0.0;
2021 |           for (unsigned int k = 0; k < n_dofs_per_cell; ++k)
2022 |             {
2023 |               m_Nx_phasefield[q_point][k]           = 0.0;
2024 |               m_grad_Nx_phasefield[q_point][k]      = 0.0;
2025 |               m_Nx_disp[q_point][k]                 = 0.0;
2026 |               m_grad_Nx_disp[q_point][k]            = 0.0;
2027 |               m_symm_grad_Nx_disp[q_point][k]       = 0.0;
2028 |             }
2029 |         }
2030 |     }
2031 |   };
2032 | 
2033 |   // constructor has no return type
2034 |   template <int dim>
2035 |   PhaseFieldMonolithicSolve<dim>::PhaseFieldMonolithicSolve(const std::string &input_file)
2036 |     : m_parameters(input_file)
2037 |     , m_triangulation(Triangulation<dim>::maximum_smoothing)
2038 |     , m_time(m_parameters.m_end_time)
2039 |     , m_logfile(m_parameters.m_logfile_name)
2040 |     , m_timer(m_logfile, TimerOutput::summary, TimerOutput::wall_times)
2041 |     , m_dof_handler(m_triangulation)
2042 |     , m_fe(FE_Q<dim>(m_parameters.m_poly_degree),
2043 | 	   dim, // displacement
2044 | 	   FE_Q<dim>(m_parameters.m_poly_degree),
2045 | 	   1)   // phasefield
2046 |     , m_dofs_per_cell(m_fe.n_dofs_per_cell())
2047 |     , m_u_fe(m_first_u_component)
2048 |     , m_d_fe(m_d_component)
2049 |     , m_dofs_per_block(m_n_blocks)
2050 |     , m_qf_cell(m_parameters.m_quad_order)
2051 |     , m_qf_face(m_parameters.m_quad_order)
2052 |     , m_n_q_points(m_qf_cell.size())
2053 |     , m_vol_reference(0.0)
2054 |   {}
2055 | 
2056 |   template <int dim>
2057 |   void PhaseFieldMonolithicSolve<dim>::make_grid()
2058 |   {
2059 |     if (m_parameters.m_scenario == 1)
2060 |       make_grid_case_1();
2061 |     else if (m_parameters.m_scenario == 2)
2062 |       make_grid_case_2();
2063 |     else if (m_parameters.m_scenario == 3)
2064 |       make_grid_case_3();
2065 |     else if (m_parameters.m_scenario == 4)
2066 |       make_grid_case_4();
2067 |     else if (m_parameters.m_scenario == 5)
2068 |       make_grid_case_5();
2069 |     else if (m_parameters.m_scenario == 6)
2070 |       make_grid_case_6();
2071 |     else if (m_parameters.m_scenario == 7)
2072 |       make_grid_case_7();
2073 |     else if (m_parameters.m_scenario == 8)
2074 |       make_grid_case_8();
2075 |     else if (m_parameters.m_scenario == 9)
2076 |       make_grid_case_9();
2077 |     else if (m_parameters.m_scenario == 10)
2078 |       make_grid_case_10();
2079 |     else if (m_parameters.m_scenario == 11)
2080 |       make_grid_case_11();
2081 |     else
2082 |       Assert(false, ExcMessage("The scenario has not been implemented!"));
2083 | 
2084 |     m_logfile << "\t\tTriangulation:"
2085 |               << "\n\t\t\tNumber of active cells: "
2086 |               << m_triangulation.n_active_cells()
2087 |               << "\n\t\t\tNumber of used vertices: "
2088 |               << m_triangulation.n_used_vertices()
2089 | 	      << std::endl;
2090 | 
2091 |     std::ofstream out("original_mesh.vtu");
2092 |     GridOut       grid_out;
2093 |     grid_out.write_vtu(m_triangulation, out);
2094 | 
2095 |     m_vol_reference = GridTools::volume(m_triangulation);
2096 |     m_logfile << "\t\tGrid:\n\t\t\tReference volume: " << m_vol_reference << std::endl;
2097 |   }
2098 | 
2099 |   template <int dim>
2100 |   void PhaseFieldMonolithicSolve<dim>::make_grid_case_1()
2101 |   {
2102 |     for (unsigned int i = 0; i < 80; ++i)
2103 |       m_logfile << "*";
2104 |     m_logfile << std::endl;
2105 |     m_logfile << "\t\t\tSquare tension (unstructured)" << std::endl;
2106 |     for (unsigned int i = 0; i < 80; ++i)
2107 |       m_logfile << "*";
2108 |     m_logfile << std::endl;
2109 | 
2110 |     AssertThrow(dim==2, ExcMessage("The dimension has to be 2D!"));
2111 | 
2112 |     GridIn<dim> gridin;
2113 |     gridin.attach_triangulation(m_triangulation);
2114 |     std::ifstream f("square_tension_unstructured.msh");
2115 |     gridin.read_msh(f);
2116 | 
2117 |     for (const auto &cell : m_triangulation.active_cell_iterators())
2118 |       for (const auto &face : cell->face_iterators())
2119 | 	{
2120 | 	  if (face->at_boundary() == true)
2121 | 	    {
2122 | 	      if (std::fabs(face->center()[1] + 0.5 ) < 1.0e-9 )
2123 | 		face->set_boundary_id(0);
2124 | 	      else if (std::fabs(face->center()[1] - 0.5 ) < 1.0e-9)
2125 | 	        face->set_boundary_id(1);
2126 | 	      else
2127 | 	        face->set_boundary_id(2);
2128 | 	    }
2129 | 	}
2130 | 
2131 |     m_triangulation.refine_global(m_parameters.m_global_refine_times);
2132 | 
2133 |     if (m_parameters.m_refinement_strategy == "pre-refine")
2134 |       {
2135 | 	unsigned int material_id;
2136 | 	double length_scale;
2137 | 	for (unsigned int i = 0; i < m_parameters.m_local_prerefine_times; i++)
2138 | 	  {
2139 | 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2140 | 	      {
2141 | 		if (   std::fabs(cell->center()[1]) < 0.01
2142 | 		    && cell->center()[0] > 0.495)
2143 | 		  {
2144 | 		    material_id = cell->material_id();
2145 | 		    length_scale = m_material_data[material_id][2];
2146 | 		    if (  std::sqrt(cell->measure())
2147 | 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2148 | 		      cell->set_refine_flag();
2149 | 		  }
2150 | 	      }
2151 | 	    m_triangulation.execute_coarsening_and_refinement();
2152 | 	  }
2153 |       }
2154 |     else if (m_parameters.m_refinement_strategy == "adaptive-refine")
2155 |       {
2156 | 	unsigned int material_id;
2157 | 	double length_scale;
2158 | 	bool initiation_point_refine_unfinished = true;
2159 | 	while (initiation_point_refine_unfinished)
2160 | 	  {
2161 | 	    initiation_point_refine_unfinished = false;
2162 | 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2163 | 	      {
2164 | 		if (   std::fabs(cell->center()[1] - 0.0) < 0.05
2165 | 		    && std::fabs(cell->center()[0] - 0.5) < 0.05)
2166 | 		  {
2167 | 		    material_id = cell->material_id();
2168 | 		    length_scale = m_material_data[material_id][2];
2169 | 		    if (  std::sqrt(cell->measure())
2170 | 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2171 | 		      {
2172 | 		        cell->set_refine_flag();
2173 | 		        initiation_point_refine_unfinished = true;
2174 | 		      }
2175 | 		  }
2176 | 	      }
2177 | 	    m_triangulation.execute_coarsening_and_refinement();
2178 | 	  }
2179 |       }
2180 |     else
2181 |       {
2182 | 	AssertThrow(false,
2183 | 	            ExcMessage("Selected mesh refinement strategy not implemented!"));
2184 |       }
2185 |   }
2186 | 
2187 | 
2188 |   template <int dim>
2189 |   void PhaseFieldMonolithicSolve<dim>::make_grid_case_2()
2190 |   {
2191 |     for (unsigned int i = 0; i < 80; ++i)
2192 |       m_logfile << "*";
2193 |     m_logfile << std::endl;
2194 |     m_logfile << "\t\t\t\tSquare shear (unstructured)" << std::endl;
2195 |     for (unsigned int i = 0; i < 80; ++i)
2196 |       m_logfile << "*";
2197 |     m_logfile << std::endl;
2198 | 
2199 |     AssertThrow(dim==2, ExcMessage("The dimension has to be 2D!"));
2200 | 
2201 |     GridIn<dim> gridin;
2202 |     gridin.attach_triangulation(m_triangulation);
2203 |     std::ifstream f("square_shear_unstructured.msh");
2204 |     gridin.read_msh(f);
2205 | 
2206 |     for (const auto &cell : m_triangulation.active_cell_iterators())
2207 |       for (const auto &face : cell->face_iterators())
2208 | 	{
2209 | 	  if (face->at_boundary() == true)
2210 | 	    {
2211 | 	      if (std::fabs(face->center()[1] + 0.5 ) < 1.0e-9 )
2212 | 		face->set_boundary_id(0);
2213 | 	      else if (std::fabs(face->center()[1] - 0.5 ) < 1.0e-9)
2214 | 	        face->set_boundary_id(1);
2215 | 	      else if (   (std::fabs(face->center()[0] - 0.0 ) < 1.0e-9)
2216 | 		       || (std::fabs(face->center()[0] - 1.0 ) < 1.0e-9))
2217 | 	        face->set_boundary_id(2);
2218 | 	      else
2219 | 	        face->set_boundary_id(3);
2220 | 	    }
2221 | 	}
2222 | 
2223 |     m_triangulation.refine_global(m_parameters.m_global_refine_times);
2224 | 
2225 |     if (m_parameters.m_refinement_strategy == "pre-refine")
2226 |       {
2227 | 	unsigned int material_id;
2228 | 	double length_scale;
2229 | 	for (unsigned int i = 0; i < m_parameters.m_local_prerefine_times; i++)
2230 | 	  {
2231 | 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2232 | 	      {
2233 | 		if (    (cell->center()[0] > 0.45)
2234 | 		     && (cell->center()[1] < 0.05) )
2235 | 		  {
2236 | 		    material_id = cell->material_id();
2237 | 		    length_scale = m_material_data[material_id][2];
2238 | 		    if (  std::sqrt(cell->measure())
2239 | 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2240 | 		      cell->set_refine_flag();
2241 | 		  }
2242 | 	      }
2243 | 	    m_triangulation.execute_coarsening_and_refinement();
2244 | 	  }
2245 |       }
2246 |     else if (m_parameters.m_refinement_strategy == "adaptive-refine")
2247 |       {
2248 | 	unsigned int material_id;
2249 | 	double length_scale;
2250 | 	bool initiation_point_refine_unfinished = true;
2251 | 	while (initiation_point_refine_unfinished)
2252 | 	  {
2253 | 	    initiation_point_refine_unfinished = false;
2254 | 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2255 | 	      {
2256 | 		if (    std::fabs(cell->center()[0] - 0.5) < 0.025
2257 | 		     && cell->center()[1] < 0.0 && cell->center()[1] > -0.025)
2258 | 		  {
2259 | 		    material_id = cell->material_id();
2260 | 		    length_scale = m_material_data[material_id][2];
2261 | 		    if (  std::sqrt(cell->measure())
2262 | 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2263 | 		      {
2264 | 		        cell->set_refine_flag();
2265 | 		        initiation_point_refine_unfinished = true;
2266 | 		      }
2267 | 		  }
2268 | 	      }
2269 | 	    m_triangulation.execute_coarsening_and_refinement();
2270 | 	  }
2271 |       }
2272 |     else
2273 |       {
2274 | 	AssertThrow(false,
2275 | 	            ExcMessage("Selected mesh refinement strategy not implemented!"));
2276 |       }
2277 |   }
2278 | 
2279 |   template <int dim>
2280 |   void PhaseFieldMonolithicSolve<dim>::make_grid_case_3()
2281 |   {
2282 |     for (unsigned int i = 0; i < 80; ++i)
2283 |       m_logfile << "*";
2284 |     m_logfile << std::endl;
2285 |     m_logfile << "\t\t\tSquare tension (structured)" << std::endl;
2286 |     for (unsigned int i = 0; i < 80; ++i)
2287 |       m_logfile << "*";
2288 |     m_logfile << std::endl;
2289 | 
2290 |     AssertThrow(dim==2, ExcMessage("The dimension has to be 2D!"));
2291 | 
2292 |     GridIn<dim> gridin;
2293 |     gridin.attach_triangulation(m_triangulation);
2294 |     std::ifstream f("square_tension_structured.msh");
2295 |     gridin.read_msh(f);
2296 | 
2297 |     for (const auto &cell : m_triangulation.active_cell_iterators())
2298 |       for (const auto &face : cell->face_iterators())
2299 | 	{
2300 | 	  if (face->at_boundary() == true)
2301 | 	    {
2302 | 	      if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9 )
2303 | 		face->set_boundary_id(0);
2304 | 	      else if (std::fabs(face->center()[1] - 1.0 ) < 1.0e-9)
2305 | 	        face->set_boundary_id(1);
2306 | 	      else
2307 | 	        face->set_boundary_id(2);
2308 | 	    }
2309 | 	}
2310 | 
2311 |     m_triangulation.refine_global(m_parameters.m_global_refine_times);
2312 | 
2313 |     if (m_parameters.m_refinement_strategy == "pre-refine")
2314 |       {
2315 | 	unsigned int material_id;
2316 | 	double length_scale;
2317 | 	for (unsigned int i = 0; i < m_parameters.m_local_prerefine_times; i++)
2318 | 	  {
2319 | 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2320 | 	      {
2321 | 		if (    (std::fabs(cell->center()[1] - 0.5) < 0.025)
2322 | 		     && (cell->center()[0] > 0.475) )
2323 | 		  {
2324 | 		    material_id = cell->material_id();
2325 | 		    length_scale = m_material_data[material_id][2];
2326 | 		    if (  std::sqrt(cell->measure())
2327 | 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2328 | 		      cell->set_refine_flag();
2329 | 		  }
2330 | 	      }
2331 | 	    m_triangulation.execute_coarsening_and_refinement();
2332 | 	  }
2333 |       }
2334 |     else if (m_parameters.m_refinement_strategy == "adaptive-refine")
2335 |       {
2336 | 	unsigned int material_id;
2337 | 	double length_scale;
2338 | 	bool initiation_point_refine_unfinished = true;
2339 | 	while (initiation_point_refine_unfinished)
2340 | 	  {
2341 | 	    initiation_point_refine_unfinished = false;
2342 | 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2343 | 	      {
2344 | 		if (    std::fabs(cell->center()[0] - 0.5) < 0.025
2345 | 		     && std::fabs(cell->center()[1] - 0.5) < 0.025 )
2346 | 		  {
2347 | 		    material_id = cell->material_id();
2348 | 		    length_scale = m_material_data[material_id][2];
2349 | 		    if (  std::sqrt(cell->measure())
2350 | 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2351 | 		      {
2352 | 		        cell->set_refine_flag();
2353 | 		        initiation_point_refine_unfinished = true;
2354 | 		      }
2355 | 		  }
2356 | 	      }
2357 | 	    m_triangulation.execute_coarsening_and_refinement();
2358 | 	  }
2359 |       }
2360 |     else
2361 |       {
2362 | 	AssertThrow(false,
2363 | 	            ExcMessage("Selected mesh refinement strategy not implemented!"));
2364 |       }
2365 |   }
2366 | 
2367 |   template <int dim>
2368 |   void PhaseFieldMonolithicSolve<dim>::make_grid_case_4()
2369 |   {
2370 |     for (unsigned int i = 0; i < 80; ++i)
2371 |       m_logfile << "*";
2372 |     m_logfile << std::endl;
2373 |     m_logfile << "\t\t\t\tSquare shear (structured)" << std::endl;
2374 |     for (unsigned int i = 0; i < 80; ++i)
2375 |       m_logfile << "*";
2376 |     m_logfile << std::endl;
2377 | 
2378 |     AssertThrow(dim==2, ExcMessage("The dimension has to be 2D!"));
2379 | 
2380 |     GridIn<dim> gridin;
2381 |     gridin.attach_triangulation(m_triangulation);
2382 |     std::ifstream f("square_shear_structured.msh");
2383 |     gridin.read_msh(f);
2384 | 
2385 |     for (const auto &cell : m_triangulation.active_cell_iterators())
2386 |       for (const auto &face : cell->face_iterators())
2387 | 	{
2388 | 	  if (face->at_boundary() == true)
2389 | 	    {
2390 | 	      if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9 )
2391 | 		face->set_boundary_id(0);
2392 | 	      else if (std::fabs(face->center()[1] - 1.0 ) < 1.0e-9)
2393 | 	        face->set_boundary_id(1);
2394 | 	      else if (   (std::fabs(face->center()[0] - 0.0 ) < 1.0e-9)
2395 | 		       || (std::fabs(face->center()[0] - 1.0 ) < 1.0e-9))
2396 | 	        face->set_boundary_id(2);
2397 | 	      else
2398 | 	        face->set_boundary_id(3);
2399 | 	    }
2400 | 	}
2401 | 
2402 |     m_triangulation.refine_global(m_parameters.m_global_refine_times);
2403 | 
2404 |     if (m_parameters.m_refinement_strategy == "pre-refine")
2405 |       {
2406 | 	unsigned int material_id;
2407 | 	double length_scale;
2408 | 	for (unsigned int i = 0; i < m_parameters.m_local_prerefine_times; i++)
2409 | 	  {
2410 | 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2411 | 	      {
2412 | 		if (   (   (cell->center()[0] > 0.475)
2413 | 		        && (cell->center()[1] < 0.525) )
2414 | 		    //|| (    cell->center()[1] > 0.975)
2415 | 		   )
2416 | 		  {
2417 | 		    material_id = cell->material_id();
2418 | 		    length_scale = m_material_data[material_id][2];
2419 | 		    if (  std::sqrt(cell->measure())
2420 | 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2421 | 		      cell->set_refine_flag();
2422 | 		  }
2423 | 	      }
2424 | 	    m_triangulation.execute_coarsening_and_refinement();
2425 | 	  }
2426 |       }
2427 |     else if (m_parameters.m_refinement_strategy == "adaptive-refine")
2428 |       {
2429 | 	unsigned int material_id;
2430 | 	double length_scale;
2431 | 	bool initiation_point_refine_unfinished = true;
2432 | 	while (initiation_point_refine_unfinished)
2433 | 	  {
2434 | 	    initiation_point_refine_unfinished = false;
2435 | 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2436 | 	      {
2437 | 		if (   (   std::fabs(cell->center()[0] - 0.5) < 0.025
2438 | 		        && cell->center()[1] < 0.5 && cell->center()[1] > 0.475 )
2439 | 		    // we also need to refine the top edge, since there might be a conflict between
2440 | 		    // inhomogeneous boundary conditions and the hanging-node constraints at the
2441 | 		    // top edge
2442 | 		    || (   cell->center()[1] > 0.975 ) )
2443 | 		  {
2444 | 		    material_id = cell->material_id();
2445 | 		    length_scale = m_material_data[material_id][2];
2446 | 		    if (  std::sqrt(cell->measure())
2447 | 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2448 | 		      {
2449 | 		        cell->set_refine_flag();
2450 | 		        initiation_point_refine_unfinished = true;
2451 | 		      }
2452 | 		  }
2453 | 	      }
2454 | 	    m_triangulation.execute_coarsening_and_refinement();
2455 | 	  }
2456 |       }
2457 |     else
2458 |       {
2459 | 	AssertThrow(false,
2460 | 	            ExcMessage("Selected mesh refinement strategy not implemented!"));
2461 |       }
2462 |   }
2463 | 
2464 |   template <int dim>
2465 |   void PhaseFieldMonolithicSolve<dim>::make_grid_case_5()
2466 |   {
2467 |     for (unsigned int i = 0; i < 80; ++i)
2468 |       m_logfile << "*";
2469 |     m_logfile << std::endl;
2470 |     m_logfile << "\t\t\t\tThree-point bending (structured)" << std::endl;
2471 |     for (unsigned int i = 0; i < 80; ++i)
2472 |       m_logfile << "*";
2473 |     m_logfile << std::endl;
2474 | 
2475 |     AssertThrow(dim==2, ExcMessage("The dimension has to be 2D!"));
2476 | 
2477 |     GridIn<dim> gridin;
2478 |     gridin.attach_triangulation(m_triangulation);
2479 |     std::ifstream f("three_point_bending_structured.msh");
2480 |     gridin.read_msh(f);
2481 | 
2482 |     for (const auto &cell : m_triangulation.active_cell_iterators())
2483 |       for (const auto &face : cell->face_iterators())
2484 | 	{
2485 | 	  if (face->at_boundary() == true)
2486 | 	    {
2487 | 	      if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9 )
2488 | 		face->set_boundary_id(0);
2489 | 	      else if (std::fabs(face->center()[1] - 2.0 ) < 1.0e-9)
2490 | 	        face->set_boundary_id(1);
2491 | 	      else
2492 | 	        face->set_boundary_id(2);
2493 | 	    }
2494 | 	}
2495 | 
2496 |     m_triangulation.refine_global(m_parameters.m_global_refine_times);
2497 | 
2498 |     if (m_parameters.m_refinement_strategy == "pre-refine")
2499 |       {
2500 | 	for (const auto &cell : m_triangulation.active_cell_iterators())
2501 | 	  {
2502 | 	    if (    std::fabs(cell->center()[0] - 4.0) < 0.075
2503 | 		 && cell->center()[1] < 1.6)
2504 | 	      {
2505 | 		cell->set_refine_flag();
2506 | 	      }
2507 | 	  }
2508 | 	m_triangulation.execute_coarsening_and_refinement();
2509 | 
2510 | 	unsigned int material_id;
2511 | 	double length_scale;
2512 | 	for (unsigned int i = 0; i < m_parameters.m_local_prerefine_times; i++)
2513 | 	  {
2514 | 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2515 | 	      {
2516 | 		if (    std::fabs(cell->center()[0] - 4.0) < 0.05
2517 | 		     && cell->center()[1] < 1.6)
2518 | 		  {
2519 | 		    material_id = cell->material_id();
2520 | 		    length_scale = m_material_data[material_id][2];
2521 | 		    if (  std::sqrt(cell->measure())
2522 | 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2523 | 		      cell->set_refine_flag();
2524 | 		  }
2525 | 	      }
2526 | 	    m_triangulation.execute_coarsening_and_refinement();
2527 | 	  }
2528 |       }
2529 |     else if (m_parameters.m_refinement_strategy == "adaptive-refine")
2530 |       {
2531 | 	unsigned int material_id;
2532 | 	double length_scale;
2533 | 	bool initiation_point_refine_unfinished = true;
2534 | 	while (initiation_point_refine_unfinished)
2535 | 	  {
2536 | 	    initiation_point_refine_unfinished = false;
2537 | 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2538 | 	      {
2539 | 		if (    std::fabs(cell->center()[0] - 4.0) < 0.075
2540 | 		     && std::fabs(cell->center()[1] - 0.4) < 0.075 )
2541 | 		  {
2542 | 		    material_id = cell->material_id();
2543 | 		    length_scale = m_material_data[material_id][2];
2544 | 		    if (  std::sqrt(cell->measure())
2545 | 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2546 | 		      {
2547 | 		        cell->set_refine_flag();
2548 | 		        initiation_point_refine_unfinished = true;
2549 | 		      }
2550 | 		  }
2551 | 	      }
2552 | 	    m_triangulation.execute_coarsening_and_refinement();
2553 | 	  }
2554 |       }
2555 |     else
2556 |       {
2557 | 	AssertThrow(false,
2558 | 	            ExcMessage("Selected mesh refinement strategy not implemented!"));
2559 |       }
2560 |   }
2561 | 
2562 |   template <int dim>
2563 |   void PhaseFieldMonolithicSolve<dim>::make_grid_case_6()
2564 |   {
2565 |     AssertThrow(dim==3, ExcMessage("The dimension has to be 3D!"));
2566 | 
2567 |     for (unsigned int i = 0; i < 80; ++i)
2568 |       m_logfile << "*";
2569 |     m_logfile << std::endl;
2570 |     m_logfile << "\t\t\t\tSphere inclusion (3D structured)" << std::endl;
2571 |     for (unsigned int i = 0; i < 80; ++i)
2572 |       m_logfile << "*";
2573 |     m_logfile << std::endl;
2574 | 
2575 |     Triangulation<dim> tria_inner;
2576 |     GridGenerator::hyper_ball(tria_inner, Point<dim>(), 0.5);
2577 | 
2578 |     Triangulation<dim> tria_outer;
2579 |     GridGenerator::hyper_shell(
2580 |       tria_outer, Point<dim>(), 0.5, std::sqrt(dim), 2 * dim);
2581 | 
2582 |     Triangulation<dim> tmp_triangulation;
2583 | 
2584 |     GridGenerator::merge_triangulations(tria_inner, tria_outer, tmp_triangulation);
2585 | 
2586 |     tmp_triangulation.reset_all_manifolds();
2587 |     tmp_triangulation.set_all_manifold_ids(0);
2588 | 
2589 |     for (const auto &cell : tmp_triangulation.cell_iterators())
2590 |       {
2591 |         for (const auto &face : cell->face_iterators())
2592 |           {
2593 |             bool face_at_sphere_boundary = true;
2594 |             for (const auto v : face->vertex_indices())
2595 |               {
2596 |                 if (std::abs(face->vertex(v).norm_square() - 0.25) > 1e-12)
2597 |                   face_at_sphere_boundary = false;
2598 |               }
2599 |             if (face_at_sphere_boundary)
2600 |               face->set_all_manifold_ids(1);
2601 |           }
2602 |         if (cell->center().norm_square() < 0.25)
2603 |           cell->set_material_id(1);
2604 |         else
2605 |           cell->set_material_id(0);
2606 |       }
2607 | 
2608 |     tmp_triangulation.set_manifold(1, SphericalManifold<dim>());
2609 | 
2610 |     TransfiniteInterpolationManifold<dim> transfinite_manifold;
2611 |     transfinite_manifold.initialize(tmp_triangulation);
2612 |     tmp_triangulation.set_manifold(0, transfinite_manifold);
2613 | 
2614 |     tmp_triangulation.refine_global(m_parameters.m_global_refine_times);
2615 | 
2616 |     std::set<typename Triangulation< dim >::active_cell_iterator >
2617 |       cells_to_remove;
2618 | 
2619 |     for (const auto &cell : tmp_triangulation.active_cell_iterators())
2620 |       {
2621 | 	if (   cell->center()[0] < 0.0
2622 | 	    || cell->center()[1] < 0.0
2623 | 	    || cell->center()[2] < 0.0)
2624 | 	  {
2625 | 	    cells_to_remove.insert(cell);
2626 | 	  }
2627 |       }
2628 | 
2629 |     GridGenerator::create_triangulation_with_removed_cells(tmp_triangulation,
2630 | 							   cells_to_remove,
2631 | 							   m_triangulation);
2632 | 
2633 |     for (const auto &cell : m_triangulation.active_cell_iterators())
2634 |       for (const auto &face : cell->face_iterators())
2635 | 	{
2636 | 	  if (face->at_boundary() == true)
2637 | 	    {
2638 | 	      if (std::fabs(face->center()[0] - 0.0 ) < 1.0e-9 )
2639 | 		face->set_boundary_id(0);
2640 | 	      else if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9)
2641 | 		face->set_boundary_id(1);
2642 | 	      else if (std::fabs(face->center()[2] - 0.0 ) < 1.0e-9)
2643 | 		face->set_boundary_id(2);
2644 | 	      else if (std::fabs(face->center()[2] - 1.0 ) < 1.0e-9)
2645 | 		face->set_boundary_id(3);
2646 | 	      else
2647 | 		face->set_boundary_id(4);
2648 | 	    }
2649 | 	}
2650 | 
2651 |     if (m_parameters.m_refinement_strategy == "adaptive-refine")
2652 |       {
2653 | 	unsigned int material_id;
2654 | 	double length_scale;
2655 | 	bool initiation_point_refine_unfinished = true;
2656 | 	while (initiation_point_refine_unfinished)
2657 | 	  {
2658 | 	    initiation_point_refine_unfinished = false;
2659 | 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2660 | 	      {
2661 | 		if (    cell->center()[2] > 0.525
2662 | 		     && cell->center()[2] < 0.575
2663 | 		     && cell->center()[0] < 0.05
2664 | 		     && cell->center()[1] < 0.05 )
2665 | 		  {
2666 | 		    material_id = cell->material_id();
2667 | 		    length_scale = m_material_data[material_id][2];
2668 | 		    if (  std::cbrt(cell->measure())
2669 | 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2670 | 		      {
2671 | 			cell->set_refine_flag();
2672 | 			initiation_point_refine_unfinished = true;
2673 | 		      }
2674 | 		  }
2675 | 	      }
2676 | 	    m_triangulation.execute_coarsening_and_refinement();
2677 | 	  }
2678 |       }
2679 |     else
2680 |       {
2681 | 	AssertThrow(false,
2682 | 		    ExcMessage("Selected mesh refinement strategy not implemented!"));
2683 |       }
2684 |   }
2685 | 
2686 |   template <int dim>
2687 |   void PhaseFieldMonolithicSolve<dim>::make_grid_case_7()
2688 |   {
2689 |     AssertThrow(dim==3, ExcMessage("The dimension has to be 3D!"));
2690 | 
2691 |     for (unsigned int i = 0; i < 80; ++i)
2692 |       m_logfile << "*";
2693 |     m_logfile << std::endl;
2694 |     m_logfile << "\t\t\t\tSphere inclusion (3D structured version 2)" << std::endl;
2695 |     for (unsigned int i = 0; i < 80; ++i)
2696 |       m_logfile << "*";
2697 |     m_logfile << std::endl;
2698 | 
2699 |     Triangulation<dim> tria_inner;
2700 |     GridGenerator::hyper_ball(tria_inner, Point<dim>(), 0.49);
2701 | 
2702 |     Triangulation<dim> tria_outer;
2703 |     GridGenerator::hyper_shell(
2704 |       tria_outer, Point<dim>(), 0.49, std::sqrt(dim)*0.5, 2 * dim);
2705 | 
2706 |     Triangulation<dim> cube1;
2707 |     GridGenerator::hyper_rectangle(cube1, Point<dim>(0, 0, 0.5), Point<dim>(1, 1, 1.5));
2708 |     Triangulation<dim> cube2;
2709 |     GridGenerator::hyper_rectangle(cube2, Point<dim>(0, 0.5, -0.5), Point<dim>(1, 1.5, 0.5));
2710 |     Triangulation<dim> cube3;
2711 |     GridGenerator::hyper_rectangle(cube3, Point<dim>(0.5, -0.5, -0.5), Point<dim>(1.5, 0.5, 0.5));
2712 | 
2713 |     Triangulation<dim> tmp_triangulation;
2714 |     GridGenerator::merge_triangulations({&tria_inner, &tria_outer,
2715 |                                          &cube1, &cube2, &cube3}, tmp_triangulation);
2716 | 
2717 |     tmp_triangulation.reset_all_manifolds();
2718 |     tmp_triangulation.set_all_manifold_ids(0);
2719 | 
2720 |     for (const auto &cell : tmp_triangulation.cell_iterators())
2721 |       {
2722 |         for (const auto &face : cell->face_iterators())
2723 |           {
2724 |             bool face_at_sphere_boundary = true;
2725 |             for (const auto v : face->vertex_indices())
2726 |               {
2727 |                 if (std::abs(face->vertex(v).norm_square() - 0.49 * 0.49) > 1e-12)
2728 |                   face_at_sphere_boundary = false;
2729 |               }
2730 |             if (face_at_sphere_boundary)
2731 |               face->set_all_manifold_ids(1);
2732 |           }
2733 |         if (cell->center().norm_square() < 0.1)
2734 |           cell->set_material_id(1);
2735 |         else
2736 |           cell->set_material_id(0);
2737 |       }
2738 | 
2739 |     tmp_triangulation.set_manifold(1, SphericalManifold<dim>());
2740 | 
2741 |     TransfiniteInterpolationManifold<dim> transfinite_manifold;
2742 |     transfinite_manifold.initialize(tmp_triangulation);
2743 |     tmp_triangulation.set_manifold(0, transfinite_manifold);
2744 | 
2745 |     tmp_triangulation.refine_global(m_parameters.m_global_refine_times);
2746 | 
2747 |     std::set<typename Triangulation< dim >::active_cell_iterator >
2748 |       cells_to_remove;
2749 | 
2750 |     for (const auto &cell : tmp_triangulation.active_cell_iterators())
2751 |       {
2752 | 	if (   cell->center()[0] < 0.0
2753 | 	    || cell->center()[1] < 0.0
2754 | 	    || cell->center()[2] < 0.0
2755 | 	    || cell->center()[0] > 1.0
2756 | 	    || cell->center()[1] > 1.0
2757 | 	    || cell->center()[2] > 1.0)
2758 | 	  {
2759 | 	    cells_to_remove.insert(cell);
2760 | 	  }
2761 |       }
2762 | 
2763 |     GridGenerator::create_triangulation_with_removed_cells(tmp_triangulation,
2764 | 							   cells_to_remove,
2765 | 							   m_triangulation);
2766 | 
2767 |     for (const auto &cell : m_triangulation.active_cell_iterators())
2768 |       for (const auto &face : cell->face_iterators())
2769 | 	{
2770 | 	  if (face->at_boundary() == true)
2771 | 	    {
2772 | 	      if (std::fabs(face->center()[0] - 0.0 ) < 1.0e-9 )
2773 | 		face->set_boundary_id(0);
2774 | 	      else if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9)
2775 | 		face->set_boundary_id(1);
2776 | 	      else if (std::fabs(face->center()[2] - 0.0 ) < 1.0e-9)
2777 | 		face->set_boundary_id(2);
2778 | 	      else if (std::fabs(face->center()[2] - 1.0 ) < 1.0e-9)
2779 | 		face->set_boundary_id(3);
2780 | 	      else
2781 | 		face->set_boundary_id(4);
2782 | 	    }
2783 | 	}
2784 | 
2785 |     if (m_parameters.m_refinement_strategy == "adaptive-refine")
2786 |       {
2787 | 	unsigned int material_id;
2788 | 	double length_scale;
2789 | 	bool initiation_point_refine_unfinished = true;
2790 | 	while (initiation_point_refine_unfinished)
2791 | 	  {
2792 | 	    initiation_point_refine_unfinished = false;
2793 | 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2794 | 	      {
2795 | 		if (    cell->center()[2] > 0.505
2796 | 		     && cell->center()[2] < 0.575
2797 | 		     && cell->center()[0] < 0.05
2798 | 		     && cell->center()[1] < 0.05 )
2799 | 		  {
2800 | 		    material_id = cell->material_id();
2801 | 		    length_scale = m_material_data[material_id][2];
2802 | 		    if (  std::cbrt(cell->measure())
2803 | 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2804 | 		      {
2805 | 			cell->set_refine_flag();
2806 | 			initiation_point_refine_unfinished = true;
2807 | 		      }
2808 | 		  }
2809 | 	      }
2810 | 	    m_triangulation.execute_coarsening_and_refinement();
2811 | 	  }
2812 |       }
2813 |     else
2814 |       {
2815 | 	AssertThrow(false,
2816 | 		    ExcMessage("Selected mesh refinement strategy not implemented!"));
2817 |       }
2818 |   }
2819 | 
2820 | 
2821 |   template <int dim>
2822 |   void PhaseFieldMonolithicSolve<dim>::make_grid_case_8()
2823 |   {
2824 |     AssertThrow(dim==3, ExcMessage("The dimension has to be 3D!"));
2825 | 
2826 |     for (unsigned int i = 0; i < 80; ++i)
2827 |       m_logfile << "*";
2828 |     m_logfile << std::endl;
2829 |     m_logfile << "\t\t\t\tSphere inclusion (3D structured version 2 with barriers)" << std::endl;
2830 |     for (unsigned int i = 0; i < 80; ++i)
2831 |       m_logfile << "*";
2832 |     m_logfile << std::endl;
2833 | 
2834 |     Triangulation<dim> tria_inner;
2835 |     GridGenerator::hyper_ball(tria_inner, Point<dim>(), 0.49);
2836 | 
2837 |     Triangulation<dim> tria_outer;
2838 |     GridGenerator::hyper_shell(
2839 |       tria_outer, Point<dim>(), 0.49, std::sqrt(dim)*0.5, 2 * dim);
2840 | 
2841 |     Triangulation<dim> cube1;
2842 |     GridGenerator::hyper_rectangle(cube1, Point<dim>(0, 0, 0.5), Point<dim>(1, 1, 1.5));
2843 |     Triangulation<dim> cube2;
2844 |     GridGenerator::hyper_rectangle(cube2, Point<dim>(0, 0.5, -0.5), Point<dim>(1, 1.5, 0.5));
2845 |     Triangulation<dim> cube3;
2846 |     GridGenerator::hyper_rectangle(cube3, Point<dim>(0.5, -0.5, -0.5), Point<dim>(1.5, 0.5, 0.5));
2847 | 
2848 |     Triangulation<dim> tmp_triangulation;
2849 |     GridGenerator::merge_triangulations({&tria_inner, &tria_outer,
2850 |                                          &cube1, &cube2, &cube3}, tmp_triangulation);
2851 | 
2852 |     tmp_triangulation.reset_all_manifolds();
2853 |     tmp_triangulation.set_all_manifold_ids(0);
2854 | 
2855 |     for (const auto &cell : tmp_triangulation.cell_iterators())
2856 |       {
2857 |         for (const auto &face : cell->face_iterators())
2858 |           {
2859 |             bool face_at_sphere_boundary = true;
2860 |             for (const auto v : face->vertex_indices())
2861 |               {
2862 |                 if (std::abs(face->vertex(v).norm_square() - 0.49 * 0.49) > 1e-12)
2863 |                   face_at_sphere_boundary = false;
2864 |               }
2865 |             if (face_at_sphere_boundary)
2866 |               face->set_all_manifold_ids(1);
2867 |           }
2868 |         if (cell->center().norm_square() < 0.1)
2869 |           cell->set_material_id(1);
2870 |         else
2871 |           cell->set_material_id(0);
2872 |       }
2873 | 
2874 |     tmp_triangulation.set_manifold(1, SphericalManifold<dim>());
2875 | 
2876 |     TransfiniteInterpolationManifold<dim> transfinite_manifold;
2877 |     transfinite_manifold.initialize(tmp_triangulation);
2878 |     tmp_triangulation.set_manifold(0, transfinite_manifold);
2879 | 
2880 |     tmp_triangulation.refine_global(m_parameters.m_global_refine_times);
2881 | 
2882 |     // some extra barriers
2883 |     for (const auto &cell : tmp_triangulation.cell_iterators())
2884 |       {
2885 |         if (    std::fabs(cell->center()[1] - 0.75) < 0.05
2886 |              && std::fabs(cell->center()[2] - 0.5625) < 0.05
2887 |              && std::fabs(cell->center()[0] - 0.0) < 0.2)
2888 |           cell->set_material_id(1);
2889 | 
2890 |         if (    std::fabs(cell->center()[1] - 0.0) < 0.2
2891 |              && std::fabs(cell->center()[2] - 0.5) < 0.1
2892 |              && std::fabs(cell->center()[0] - 0.75) < 0.05)
2893 |           cell->set_material_id(1);
2894 |       }
2895 | 
2896 |     std::set<typename Triangulation< dim >::active_cell_iterator >
2897 |       cells_to_remove;
2898 | 
2899 |     for (const auto &cell : tmp_triangulation.active_cell_iterators())
2900 |       {
2901 | 	if (   cell->center()[0] < 0.0
2902 | 	    || cell->center()[1] < 0.0
2903 | 	    || cell->center()[2] < 0.0
2904 | 	    || cell->center()[0] > 1.0
2905 | 	    || cell->center()[1] > 1.0
2906 | 	    || cell->center()[2] > 1.0)
2907 | 	  {
2908 | 	    cells_to_remove.insert(cell);
2909 | 	  }
2910 |       }
2911 | 
2912 |     GridGenerator::create_triangulation_with_removed_cells(tmp_triangulation,
2913 | 							   cells_to_remove,
2914 | 							   m_triangulation);
2915 | 
2916 |     for (const auto &cell : m_triangulation.active_cell_iterators())
2917 |       for (const auto &face : cell->face_iterators())
2918 | 	{
2919 | 	  if (face->at_boundary() == true)
2920 | 	    {
2921 | 	      if (std::fabs(face->center()[0] - 0.0 ) < 1.0e-9 )
2922 | 		face->set_boundary_id(0);
2923 | 	      else if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9)
2924 | 		face->set_boundary_id(1);
2925 | 	      else if (std::fabs(face->center()[2] - 0.0 ) < 1.0e-9)
2926 | 		face->set_boundary_id(2);
2927 | 	      else if (std::fabs(face->center()[2] - 1.0 ) < 1.0e-9)
2928 | 		face->set_boundary_id(3);
2929 | 	      else
2930 | 		face->set_boundary_id(4);
2931 | 	    }
2932 | 	}
2933 | 
2934 |     if (m_parameters.m_refinement_strategy == "adaptive-refine")
2935 |       {
2936 | 	unsigned int material_id;
2937 | 	double length_scale;
2938 | 	bool initiation_point_refine_unfinished = true;
2939 | 	while (initiation_point_refine_unfinished)
2940 | 	  {
2941 | 	    initiation_point_refine_unfinished = false;
2942 | 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2943 | 	      {
2944 | 		if (    cell->center()[2] > 0.505
2945 | 		     && cell->center()[2] < 0.575
2946 | 		     && cell->center()[0] < 0.05
2947 | 		     && cell->center()[1] < 0.05 )
2948 | 		  {
2949 | 		    material_id = cell->material_id();
2950 | 		    length_scale = m_material_data[material_id][2];
2951 | 		    if (  std::cbrt(cell->measure())
2952 | 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2953 | 		      {
2954 | 			cell->set_refine_flag();
2955 | 			initiation_point_refine_unfinished = true;
2956 | 		      }
2957 | 		  }
2958 | 	      }
2959 | 	    m_triangulation.execute_coarsening_and_refinement();
2960 | 	  }
2961 |       }
2962 |     else
2963 |       {
2964 | 	AssertThrow(false,
2965 | 		    ExcMessage("Selected mesh refinement strategy not implemented!"));
2966 |       }
2967 |   }
2968 | 
2969 |   template <int dim>
2970 |   void PhaseFieldMonolithicSolve<dim>::make_grid_case_9()
2971 |   {
2972 |     AssertThrow(dim==2, ExcMessage("The dimension has to be 2D!"));
2973 | 
2974 |     for (unsigned int i = 0; i < 80; ++i)
2975 |       m_logfile << "*";
2976 |     m_logfile << std::endl;
2977 |     m_logfile << "\t\t\t\tL-shape bending (2D structured)" << std::endl;
2978 |     for (unsigned int i = 0; i < 80; ++i)
2979 |       m_logfile << "*";
2980 |     m_logfile << std::endl;
2981 | 
2982 |     GridIn<dim> gridin;
2983 |     gridin.attach_triangulation(m_triangulation);
2984 |     std::ifstream f("L-Shape.msh");
2985 |     gridin.read_msh(f);
2986 | 
2987 |     for (const auto &cell : m_triangulation.active_cell_iterators())
2988 |       for (const auto &face : cell->face_iterators())
2989 | 	{
2990 | 	  if (face->at_boundary() == true)
2991 | 	    {
2992 | 	      if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9 )
2993 | 		face->set_boundary_id(0);
2994 | 	      else
2995 | 	        face->set_boundary_id(1);
2996 | 	    }
2997 | 	}
2998 | 
2999 |     m_triangulation.refine_global(m_parameters.m_global_refine_times);
3000 | 
3001 |     if (m_parameters.m_refinement_strategy == "pre-refine")
3002 |       {
3003 | 	unsigned int material_id;
3004 | 	double length_scale;
3005 | 	for (unsigned int i = 0; i < m_parameters.m_local_prerefine_times; i++)
3006 | 	  {
3007 | 	    for (const auto &cell : m_triangulation.active_cell_iterators())
3008 | 	      {
3009 | 		if (    (cell->center()[1] > 242.0)
3010 | 		     && (cell->center()[1] < 312.5)
3011 | 		     && (cell->center()[0] < 258.0) )
3012 | 		  {
3013 | 		    material_id = cell->material_id();
3014 | 		    length_scale = m_material_data[material_id][2];
3015 | 		    if (  std::sqrt(cell->measure())
3016 | 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
3017 | 		      cell->set_refine_flag();
3018 | 		  }
3019 | 	      }
3020 | 	    m_triangulation.execute_coarsening_and_refinement();
3021 | 	  }
3022 |       }
3023 |     else if (m_parameters.m_refinement_strategy == "adaptive-refine")
3024 |       {
3025 | 	unsigned int material_id;
3026 | 	double length_scale;
3027 | 	bool initiation_point_refine_unfinished = true;
3028 | 	while (initiation_point_refine_unfinished)
3029 | 	  {
3030 | 	    initiation_point_refine_unfinished = false;
3031 | 	    for (const auto &cell : m_triangulation.active_cell_iterators())
3032 | 	      {
3033 | 		if (             (cell->center()[0] - 250) < 0.0
3034 | 		     &&          (cell->center()[0] - 240) > 0.0
3035 | 		     && std::fabs(cell->center()[1] - 250) < 10.0 )
3036 | 		  {
3037 | 		    material_id = cell->material_id();
3038 | 		    length_scale = m_material_data[material_id][2];
3039 | 		    if (  std::sqrt(cell->measure())
3040 | 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
3041 | 		      {
3042 | 		        cell->set_refine_flag();
3043 | 		        initiation_point_refine_unfinished = true;
3044 | 		      }
3045 | 		  }
3046 | 	      }
3047 | 	    m_triangulation.execute_coarsening_and_refinement();
3048 | 	  }
3049 |       }
3050 |     else
3051 |       {
3052 | 	AssertThrow(false,
3053 | 	            ExcMessage("Selected mesh refinement strategy not implemented!"));
3054 |       }
3055 |   }
3056 | 
3057 |   template <int dim>
3058 |   void PhaseFieldMonolithicSolve<dim>::make_grid_case_10()
3059 |   {
3060 |     AssertThrow(dim==3, ExcMessage("The dimension has to be 3D!"));
3061 | 
3062 |     for (unsigned int i = 0; i < 80; ++i)
3063 |       m_logfile << "*";
3064 |     m_logfile << std::endl;
3065 |     m_logfile << "\t\t\t\tL-shape bending (3D structured)" << std::endl;
3066 |     for (unsigned int i = 0; i < 80; ++i)
3067 |       m_logfile << "*";
3068 |     m_logfile << std::endl;
3069 | 
3070 |     Triangulation<2> triangulation_2d;
3071 | 
3072 |     GridIn<2> gridin;
3073 |     gridin.attach_triangulation(triangulation_2d);
3074 |     std::ifstream f("L-Shape.msh");
3075 |     gridin.read_msh(f);
3076 | 
3077 |     const double thickness = 150.0;
3078 |     const unsigned int n_layer = 11;
3079 |     GridGenerator::extrude_triangulation(triangulation_2d, n_layer, thickness, m_triangulation);
3080 | 
3081 |     for (const auto &cell : m_triangulation.active_cell_iterators())
3082 |       for (const auto &face : cell->face_iterators())
3083 | 	{
3084 | 	  if (face->at_boundary() == true)
3085 | 	    {
3086 | 	      if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9 )
3087 | 		face->set_boundary_id(0);
3088 | 	      else
3089 | 	        face->set_boundary_id(1);
3090 | 	    }
3091 | 	}
3092 | 
3093 |     m_triangulation.refine_global(m_parameters.m_global_refine_times);
3094 | 
3095 |     if (m_parameters.m_refinement_strategy == "pre-refine")
3096 |       {
3097 | 	unsigned int material_id;
3098 | 	double length_scale;
3099 | 	for (unsigned int i = 0; i < m_parameters.m_local_prerefine_times; i++)
3100 | 	  {
3101 | 	    for (const auto &cell : m_triangulation.active_cell_iterators())
3102 | 	      {
3103 | 		if (    (std::fabs(cell->center()[1] - 250.0) < 10.0)
3104 | 		     && (cell->center()[0] < 250.0) )
3105 | 		  {
3106 | 		    material_id = cell->material_id();
3107 | 		    length_scale = m_material_data[material_id][2];
3108 | 		    if (  std::cbrt(cell->measure())
3109 | 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
3110 | 		      cell->set_refine_flag();
3111 | 		  }
3112 | 	      }
3113 | 	    m_triangulation.execute_coarsening_and_refinement();
3114 | 	  }
3115 |       }
3116 |     else if (m_parameters.m_refinement_strategy == "adaptive-refine")
3117 |       {
3118 | 	unsigned int material_id;
3119 | 	double length_scale;
3120 | 	bool initiation_point_refine_unfinished = true;
3121 | 	while (initiation_point_refine_unfinished)
3122 | 	  {
3123 | 	    initiation_point_refine_unfinished = false;
3124 | 	    for (const auto &cell : m_triangulation.active_cell_iterators())
3125 | 	      {
3126 | 		if (             (cell->center()[0] - 250) < 0.0
3127 | 		     &&          (cell->center()[0] - 240) > 0.0
3128 | 		     && std::fabs(cell->center()[1] - 250) < 10.0 )
3129 | 		  {
3130 | 		    material_id = cell->material_id();
3131 | 		    length_scale = m_material_data[material_id][2];
3132 | 		    if (  std::cbrt(cell->measure())
3133 | 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
3134 | 		      {
3135 | 		        cell->set_refine_flag();
3136 | 		        initiation_point_refine_unfinished = true;
3137 | 		      }
3138 | 		  }
3139 | 	      }
3140 | 	    m_triangulation.execute_coarsening_and_refinement();
3141 | 	  }
3142 |       }
3143 |     else
3144 |       {
3145 | 	AssertThrow(false,
3146 | 	            ExcMessage("Selected mesh refinement strategy not implemented!"));
3147 |       }
3148 |   }
3149 | 
3150 | 
3151 |   template <int dim>
3152 |   void PhaseFieldMonolithicSolve<dim>::make_grid_case_11()
3153 |   {
3154 |     AssertThrow(dim==3, ExcMessage("The dimension has to be 3D!"));
3155 | 
3156 |     for (unsigned int i = 0; i < 80; ++i)
3157 |       m_logfile << "*";
3158 |     m_logfile << std::endl;
3159 |     m_logfile << "\t\t\t\tBrokenshire torsion (3D structured)" << std::endl;
3160 |     for (unsigned int i = 0; i < 80; ++i)
3161 |       m_logfile << "*";
3162 |     m_logfile << std::endl;
3163 | 
3164 |     Triangulation<2> triangulation_2d;
3165 | 
3166 |     double const length = 200.0;
3167 |     double const width = 50.0;
3168 |     double const height = 50.0;
3169 |     double const delta_L = 25.0;
3170 |     double const tan_theta = delta_L / (0.5*width);
3171 | 
3172 |     std::vector<unsigned int> repetitions(2, 1);
3173 |     repetitions[0] = 20;
3174 |     repetitions[1] = 5;
3175 | 
3176 |     Point<2> point1(0.0, 0.0);
3177 |     Point<2> point2(length, width);
3178 | 
3179 |     GridGenerator::subdivided_hyper_rectangle(triangulation_2d,
3180 | 					      repetitions,
3181 | 					      point1,
3182 | 					      point2 );
3183 | 
3184 |     typename Triangulation<2>::vertex_iterator vertex_ptr;
3185 |     vertex_ptr = triangulation_2d.begin_active_vertex();
3186 |     while (vertex_ptr != triangulation_2d.end_vertex())
3187 |       {
3188 | 	Point<2> & vertex_point = vertex_ptr->vertex();
3189 | 
3190 | 	const double delta_x = (vertex_point(1) - 0.5*width) * tan_theta;
3191 | 
3192 | 	if (std::fabs(vertex_point(0) - 0.5*length) < 1.0e-6)
3193 | 	  {
3194 | 	    vertex_point(0) += delta_x;
3195 | 	  }
3196 | 	else if (std::fabs(vertex_point(0) + length/repetitions[0] - 0.5*length) < 1.0e-6)
3197 | 	  {
3198 | 	    vertex_point(0) += (delta_x + length/repetitions[0]*0.5);
3199 | 	  }
3200 | 	else if (std::fabs(vertex_point(0) - length/repetitions[0] - 0.5*length) < 1.0e-6)
3201 | 	  {
3202 | 	    vertex_point(0) += (delta_x - length/repetitions[0]*0.5);
3203 | 	  }
3204 | 	else if (vertex_point(0) < 0.5*length - length/repetitions[0] - 1.0e-6)
3205 | 	  {
3206 | 	    vertex_point(0) += (delta_x + length/repetitions[0]*0.5) * vertex_point(0)/(0.5*length - length/repetitions[0]);
3207 | 	  }
3208 | 	else if (vertex_point(0) > 0.5*length + length/repetitions[0] + 1.0e-6)
3209 | 	  {
3210 | 	    vertex_point(0) += (delta_x - length/repetitions[0]*0.5) * (length - vertex_point(0))/(0.5*length - length/repetitions[0]);
3211 | 	  }
3212 | 
3213 | 	++vertex_ptr;
3214 |       }
3215 | 
3216 |     Triangulation<dim> tmp_triangulation;
3217 |     const unsigned int n_layer = repetitions[1] + 1;
3218 |     GridGenerator::extrude_triangulation(triangulation_2d, n_layer, height, tmp_triangulation);
3219 | 
3220 |     tmp_triangulation.refine_global(m_parameters.m_global_refine_times);
3221 | 
3222 |     std::set<typename Triangulation< dim >::active_cell_iterator >
3223 |       cells_to_remove;
3224 | 
3225 |     for (const auto &cell : tmp_triangulation.active_cell_iterators())
3226 |       {
3227 | 	if (    (std::fabs(cell->center()[0] - (cell->center()[1] - 0.5*width)*tan_theta - 0.5*length) < 2.5)
3228 | 	     && cell->center()[2] > 0.5* height  )
3229 | 	  {
3230 | 	    cells_to_remove.insert(cell);
3231 | 	  }
3232 |       }
3233 | 
3234 |     GridGenerator::create_triangulation_with_removed_cells(tmp_triangulation,
3235 | 							   cells_to_remove,
3236 | 							   m_triangulation);
3237 | 
3238 |     if (m_parameters.m_refinement_strategy == "adaptive-refine")
3239 |       {
3240 | 	unsigned int material_id;
3241 | 	double length_scale;
3242 | 	bool initiation_point_refine_unfinished = true;
3243 | 	while (initiation_point_refine_unfinished)
3244 | 	  {
3245 | 	    initiation_point_refine_unfinished = false;
3246 | 	    for (const auto &cell : m_triangulation.active_cell_iterators())
3247 | 	      {
3248 | 		if (  (std::fabs(cell->center()[0] - (cell->center()[1] - 0.5*width)*tan_theta - 0.5*length) < 5.0)
3249 | 		    && cell->center()[2] <= 0.5*height
3250 | 		    && cell->center()[2] > 0.5*height - 5.0 )
3251 | 		  {
3252 | 		    material_id = cell->material_id();
3253 | 		    length_scale = m_material_data[material_id][2];
3254 | 		    if (  std::cbrt(cell->measure())
3255 | 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
3256 | 		      {
3257 | 			cell->set_refine_flag();
3258 | 			initiation_point_refine_unfinished = true;
3259 | 		      }
3260 | 		  }
3261 | 	      }
3262 | 	    m_triangulation.execute_coarsening_and_refinement();
3263 | 	  }
3264 |       }
3265 |     else
3266 |       {
3267 | 	AssertThrow(false,
3268 | 		    ExcMessage("Selected mesh refinement strategy not implemented!"));
3269 |       }
3270 | 
3271 | 
3272 |     for (const auto &cell : m_triangulation.active_cell_iterators())
3273 |       for (const auto &face : cell->face_iterators())
3274 | 	{
3275 | 	  if (face->at_boundary() == true)
3276 | 	    {
3277 | 	      if (std::fabs(face->center()[0] - length) < 1.0e-6 )
3278 | 		face->set_boundary_id(0);
3279 | 	      else if (std::fabs(face->center()[0] - 0.0) < 1.0e-6 )
3280 | 		face->set_boundary_id(1);
3281 | 	      else
3282 | 		face->set_boundary_id(2);
3283 | 	    }
3284 | 	}
3285 |   }
3286 | 
3287 |   template <int dim>
3288 |   void PhaseFieldMonolithicSolve<dim>::setup_system()
3289 |   {
3290 |     m_timer.enter_subsection("Setup system");
3291 | 
3292 |     std::vector<unsigned int> block_component(m_n_components,
3293 |                                               m_u_dof); // displacement
3294 |     block_component[m_d_component] = m_d_dof;           // phasefield
3295 | 
3296 |     m_dof_handler.distribute_dofs(m_fe);
3297 |     DoFRenumbering::Cuthill_McKee(m_dof_handler);
3298 |     DoFRenumbering::component_wise(m_dof_handler, block_component);
3299 | 
3300 |     m_constraints.clear();
3301 |     DoFTools::make_hanging_node_constraints(m_dof_handler, m_constraints);
3302 |     m_constraints.close();
3303 | 
3304 |     m_dofs_per_block =
3305 |       DoFTools::count_dofs_per_fe_block(m_dof_handler, block_component);
3306 | 
3307 |     m_logfile << "\t\tTriangulation:"
3308 |               << "\n\t\t\t Number of active cells: "
3309 |               << m_triangulation.n_active_cells()
3310 |               << "\n\t\t\t Number of used vertices: "
3311 |               << m_triangulation.n_used_vertices()
3312 |               << "\n\t\t\t Number of active edges: "
3313 |               << m_triangulation.n_active_lines()
3314 |               << "\n\t\t\t Number of active faces: "
3315 |               << m_triangulation.n_active_faces()
3316 |               << "\n\t\t\t Number of degrees of freedom (total): "
3317 | 	      << m_dof_handler.n_dofs()
3318 | 	      << "\n\t\t\t Number of degrees of freedom (disp): "
3319 | 	      << m_dofs_per_block[m_u_dof]
3320 | 	      << "\n\t\t\t Number of degrees of freedom (phasefield): "
3321 | 	      << m_dofs_per_block[m_d_dof]
3322 |               << std::endl;
3323 | 
3324 |     m_tangent_matrix.clear();
3325 |     {
3326 |       BlockDynamicSparsityPattern dsp(m_dofs_per_block, m_dofs_per_block);
3327 | 
3328 |       Table<2, DoFTools::Coupling> coupling(m_n_components, m_n_components);
3329 |       for (unsigned int ii = 0; ii < m_n_components; ++ii)
3330 |         for (unsigned int jj = 0; jj < m_n_components; ++jj)
3331 |           {
3332 |             if (   ((ii < m_d_component) && (jj == m_d_component))
3333 |                 || ((ii == m_d_component) && (jj < m_d_component)) )
3334 |               coupling[ii][jj] = DoFTools::none;
3335 |             else
3336 |               coupling[ii][jj] = DoFTools::always;
3337 |           }
3338 | 
3339 |       DoFTools::make_sparsity_pattern(
3340 |         m_dof_handler, coupling, dsp, m_constraints, false);
3341 |       m_sparsity_pattern.copy_from(dsp);
3342 |     }
3343 | 
3344 |     m_tangent_matrix.reinit(m_sparsity_pattern);
3345 | 
3346 |     m_system_rhs.reinit(m_dofs_per_block);
3347 |     m_solution.reinit(m_dofs_per_block);
3348 | 
3349 |     m_active_set_phasefield.reinit(m_dofs_per_block[m_d_dof]);
3350 | 
3351 |     setup_qph();
3352 | 
3353 |     m_timer.leave_subsection();
3354 |   }
3355 | 
3356 |   template <int dim>
3357 |   void PhaseFieldMonolithicSolve<dim>::make_constraints(const unsigned int it_nr)
3358 |   {
3359 |     const bool apply_dirichlet_bc = (it_nr == 0);
3360 | 
3361 |     if (it_nr > 1)
3362 |       {
3363 |         return;
3364 |       }
3365 | 
3366 |     if (apply_dirichlet_bc)
3367 |       {
3368 | 	m_constraints.clear();
3369 | 	DoFTools::make_hanging_node_constraints(m_dof_handler,
3370 | 						m_constraints);
3371 | 
3372 | 	const FEValuesExtractors::Scalar x_displacement(0);
3373 | 	const FEValuesExtractors::Scalar y_displacement(1);
3374 | 	const FEValuesExtractors::Scalar z_displacement(2);
3375 | 
3376 | 	const FEValuesExtractors::Vector displacements(0);
3377 | 
3378 | 	if (   m_parameters.m_scenario == 1
3379 | 	    || m_parameters.m_scenario == 3)
3380 | 	  {
3381 | 	    // Dirichlet B,C. bottom surface
3382 | 	    const int boundary_id_bottom_surface = 0;
3383 | 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3384 | 						     boundary_id_bottom_surface,
3385 | 						     Functions::ZeroFunction<dim>(m_n_components),
3386 | 						     m_constraints,
3387 | 						     m_fe.component_mask(y_displacement));
3388 | 
3389 | 	    typename Triangulation<dim>::active_vertex_iterator vertex_itr;
3390 | 	    vertex_itr = m_triangulation.begin_active_vertex();
3391 | 	    std::vector<types::global_dof_index> node_xy(m_fe.dofs_per_vertex);
3392 | 
3393 | 	    for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)
3394 | 	      {
3395 | 		if (   (std::fabs(vertex_itr->vertex()[0] - 0.0) < 1.0e-9)
3396 | 		    && (std::fabs(vertex_itr->vertex()[1] - 0.0) < 1.0e-9) )
3397 | 		  {
3398 | 		    node_xy = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);
3399 | 		  }
3400 | 	      }
3401 | 	    m_constraints.add_line(node_xy[0]);
3402 | 	    m_constraints.set_inhomogeneity(node_xy[0], 0.0);
3403 | 
3404 | 	    m_constraints.add_line(node_xy[1]);
3405 | 	    m_constraints.set_inhomogeneity(node_xy[1], 0.0);
3406 | 
3407 | 	    const int boundary_id_top_surface = 1;
3408 | 	    /*
3409 | 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3410 | 						     boundary_id_top_surface,
3411 | 						     Functions::ZeroFunction<dim>(m_n_components),
3412 | 						     m_constraints,
3413 | 						     m_fe.component_mask(x_displacement));
3414 | 	    */
3415 |             const double time_inc = m_time.get_delta_t();
3416 |             double disp_magnitude = m_time.get_magnitude();
3417 | 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3418 | 						     boundary_id_top_surface,
3419 | 						     Functions::ConstantFunction<dim>(
3420 | 						       disp_magnitude*time_inc, m_n_components),
3421 | 						     m_constraints,
3422 | 						     m_fe.component_mask(y_displacement));
3423 | 	  }
3424 | 	else if (   m_parameters.m_scenario == 2
3425 | 	         || m_parameters.m_scenario == 4)
3426 | 	  {
3427 | 	    // Dirichlet B,C. bottom surface
3428 | 	    const int boundary_id_bottom_surface = 0;
3429 | 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3430 | 						     boundary_id_bottom_surface,
3431 | 						     Functions::ZeroFunction<dim>(m_n_components),
3432 | 						     m_constraints,
3433 | 						     m_fe.component_mask(displacements));
3434 | 
3435 | 	    const int boundary_id_top_surface = 1;
3436 | 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3437 | 						     boundary_id_top_surface,
3438 | 						     Functions::ZeroFunction<dim>(m_n_components),
3439 | 						     m_constraints,
3440 | 						     m_fe.component_mask(y_displacement));
3441 | 
3442 | 	    const double time_inc = m_time.get_delta_t();
3443 | 	    double disp_magnitude = m_time.get_magnitude();
3444 | 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3445 | 						     boundary_id_top_surface,
3446 | 						     Functions::ConstantFunction<dim>(
3447 | 						       disp_magnitude*time_inc, m_n_components),
3448 | 						     m_constraints,
3449 | 						     m_fe.component_mask(x_displacement));
3450 | 
3451 | 	    const int boundary_id_side_surfaces = 2;
3452 | 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3453 | 						     boundary_id_side_surfaces,
3454 | 						     Functions::ZeroFunction<dim>(m_n_components),
3455 | 						     m_constraints,
3456 | 						     m_fe.component_mask(y_displacement));
3457 | 	  }
3458 | 	else if (m_parameters.m_scenario == 5)
3459 | 	  {
3460 | 	    typename Triangulation<dim>::active_vertex_iterator vertex_itr;
3461 | 	    vertex_itr = m_triangulation.begin_active_vertex();
3462 | 	    std::vector<types::global_dof_index> node_bottomleft(m_fe.dofs_per_vertex);
3463 | 	    std::vector<types::global_dof_index> node_bottomright(m_fe.dofs_per_vertex);
3464 | 	    std::vector<types::global_dof_index> node_topcenter(m_fe.dofs_per_vertex);
3465 | 
3466 | 	    for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)
3467 | 	      {
3468 | 		if (   (std::fabs(vertex_itr->vertex()[0] - 0.0) < 1.0e-9)
3469 | 		    && (std::fabs(vertex_itr->vertex()[1] - 0.0) < 1.0e-9) )
3470 | 		  {
3471 | 		    node_bottomleft = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);
3472 | 		  }
3473 | 		if (   (std::fabs(vertex_itr->vertex()[0] - 8.0) < 1.0e-9)
3474 | 		    && (std::fabs(vertex_itr->vertex()[1] - 0.0) < 1.0e-9) )
3475 | 		  {
3476 | 		    node_bottomright = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);
3477 | 		  }
3478 | 		if (   (std::fabs(vertex_itr->vertex()[0] - 4.0) < 1.0e-9)
3479 | 		    && (std::fabs(vertex_itr->vertex()[1] - 2.0) < 1.0e-9) )
3480 | 		  {
3481 | 		    node_topcenter = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);
3482 | 		  }
3483 | 	      }
3484 | 	    // bottom-left node fixed in both x- and y-directions
3485 | 	    m_constraints.add_line(node_bottomleft[0]);
3486 | 	    m_constraints.set_inhomogeneity(node_bottomleft[0], 0.0);
3487 | 
3488 | 	    m_constraints.add_line(node_bottomleft[1]);
3489 | 	    m_constraints.set_inhomogeneity(node_bottomleft[1], 0.0);
3490 | 
3491 | 	    // bottom-right node only fixed in y-direction
3492 | 	    m_constraints.add_line(node_bottomright[1]);
3493 | 	    m_constraints.set_inhomogeneity(node_bottomright[1], 0.0);
3494 | 
3495 | 	    // top-center node applied with y-displacement
3496 | 	    const double time_inc = m_time.get_delta_t();
3497 | 	    double disp_magnitude = m_time.get_magnitude();
3498 | 
3499 | 	    m_constraints.add_line(node_topcenter[1]);
3500 | 	    m_constraints.set_inhomogeneity(node_topcenter[1], disp_magnitude*time_inc);
3501 | 	  }
3502 | 	else if (   m_parameters.m_scenario == 6
3503 | 	         || m_parameters.m_scenario == 7
3504 | 		 || m_parameters.m_scenario == 8)
3505 | 	  {
3506 | 	    const int x0_surface = 0;
3507 | 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3508 | 						     x0_surface,
3509 | 						     Functions::ZeroFunction<dim>(m_n_components),
3510 | 						     m_constraints,
3511 | 						     m_fe.component_mask(x_displacement));
3512 | 	    const int y0_surface = 1;
3513 | 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3514 | 						     y0_surface,
3515 | 						     Functions::ZeroFunction<dim>(m_n_components),
3516 | 						     m_constraints,
3517 | 						     m_fe.component_mask(y_displacement));
3518 | 	    const int z0_surface = 2;
3519 | 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3520 | 						     z0_surface,
3521 | 						     Functions::ZeroFunction<dim>(m_n_components),
3522 | 						     m_constraints,
3523 | 						     m_fe.component_mask(z_displacement));
3524 | 
3525 | 	    const int z1_surface = 3;
3526 | 	    const double time_inc = m_time.get_delta_t();
3527 | 	    double disp_magnitude = 1.0;
3528 | 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3529 | 						     z1_surface,
3530 | 						     Functions::ConstantFunction<dim>(
3531 | 						       disp_magnitude*time_inc, m_n_components),
3532 | 						     m_constraints,
3533 | 						     m_fe.component_mask(z_displacement));
3534 | 	  }
3535 | 	else if (   m_parameters.m_scenario == 9
3536 | 	         || m_parameters.m_scenario == 10)
3537 | 	  {
3538 | 	    // Dirichlet B,C. bottom surface
3539 | 	    const int boundary_id_bottom_surface = 0;
3540 | 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3541 | 						     boundary_id_bottom_surface,
3542 | 						     Functions::ZeroFunction<dim>(m_n_components),
3543 | 						     m_constraints,
3544 | 						     m_fe.component_mask(displacements));
3545 | 
3546 | 	    typename Triangulation<dim>::active_vertex_iterator vertex_itr;
3547 | 	    vertex_itr = m_triangulation.begin_active_vertex();
3548 | 	    std::vector<types::global_dof_index> node_disp_control(m_fe.dofs_per_vertex);
3549 | 
3550 | 	    for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)
3551 | 	      {
3552 | 		if (   (std::fabs(vertex_itr->vertex()[0] - 470.0) < 1.0e-9)
3553 | 		    && (std::fabs(vertex_itr->vertex()[1] - 250.0) < 1.0e-9) )
3554 | 		  {
3555 | 		    node_disp_control = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);
3556 | 	            // node applied with y-displacement
3557 | 		    const double time_inc = m_time.get_delta_t();
3558 | 		    double disp_magnitude = m_time.get_magnitude();
3559 | 
3560 | 		    m_constraints.add_line(node_disp_control[1]);
3561 | 		    m_constraints.set_inhomogeneity(node_disp_control[1], disp_magnitude*time_inc);
3562 | 		  }
3563 | 	      }
3564 | 	  }
3565 | 	else if (m_parameters.m_scenario == 11)
3566 | 	  {
3567 | 	    // Dirichlet B,C. right surface
3568 | 	    const int boundary_id_right_surface = 0;
3569 | 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3570 | 						     boundary_id_right_surface,
3571 | 						     Functions::ZeroFunction<dim>(m_n_components),
3572 | 						     m_constraints,
3573 | 						     m_fe.component_mask(displacements));
3574 | 
3575 | 	    // Dirichlet B,C. left surface
3576 | 	    const int boundary_id_left_surface = 1;
3577 | 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3578 | 						     boundary_id_left_surface,
3579 | 						     Functions::ZeroFunction<dim>(m_n_components),
3580 | 						     m_constraints,
3581 | 						     m_fe.component_mask(x_displacement));
3582 | 
3583 | 	    typename Triangulation<dim>::active_vertex_iterator vertex_itr;
3584 | 	    vertex_itr = m_triangulation.begin_active_vertex();
3585 | 	    std::vector<types::global_dof_index> node_rotate(m_fe.dofs_per_vertex);
3586 | 	    double node_dist = 0.0;
3587 | 	    double disp_mag = 0.0;
3588 | 	    double angle_theta = 0.0;
3589 | 	    double disp_y = 0;
3590 | 	    double disp_z = 0;
3591 | 
3592 | 	    for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)
3593 | 	      {
3594 | 		if (std::fabs(vertex_itr->vertex()[0] - 0.0) < 1.0e-9)
3595 | 		  {
3596 | 		    node_rotate = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);
3597 | 		    node_dist = std::sqrt(  vertex_itr->vertex()[1] * vertex_itr->vertex()[1]
3598 | 			                  + vertex_itr->vertex()[2] * vertex_itr->vertex()[2]);
3599 | 
3600 | 		    angle_theta = m_time.get_delta_t() * m_time.get_magnitude();
3601 | 		    disp_mag = node_dist * std::tan(angle_theta);
3602 | 
3603 | 		    if (node_dist > 0)
3604 | 		      {
3605 | 		        disp_y = vertex_itr->vertex()[2]/node_dist * disp_mag;
3606 | 		        disp_z = -vertex_itr->vertex()[1]/node_dist * disp_mag;
3607 | 		      }
3608 | 		    else
3609 | 		      {
3610 | 			disp_y = 0.0;
3611 | 			disp_z = 0.0;
3612 | 		      }
3613 | 
3614 | 		    m_constraints.add_line(node_rotate[1]);
3615 | 		    m_constraints.set_inhomogeneity(node_rotate[1], disp_y);
3616 | 
3617 | 		    m_constraints.add_line(node_rotate[2]);
3618 | 		    m_constraints.set_inhomogeneity(node_rotate[2], disp_z);
3619 | 		  }
3620 | 	      }
3621 | 	  }
3622 | 	else
3623 | 	  Assert(false, ExcMessage("The scenario has not been implemented!"));
3624 |       }
3625 |     else  // inhomogeneous constraints
3626 |       {
3627 |         if (m_constraints.has_inhomogeneities())
3628 |           {
3629 |             AffineConstraints<double> homogeneous_constraints(m_constraints);
3630 |             for (unsigned int dof = 0; dof != m_dof_handler.n_dofs(); ++dof)
3631 |               if (homogeneous_constraints.is_inhomogeneously_constrained(dof))
3632 |                 homogeneous_constraints.set_inhomogeneity(dof, 0.0);
3633 |             m_constraints.clear();
3634 |             m_constraints.copy_from(homogeneous_constraints);
3635 |           }
3636 |       }
3637 |     m_constraints.close();
3638 |   }
3639 | 
3640 |   template <int dim>
3641 |   void PhaseFieldMonolithicSolve<dim>::assemble_system_B0(const BlockVector<double> & solution_old)
3642 |   {
3643 |     m_timer.enter_subsection("Assemble B0");
3644 | 
3645 |     m_tangent_matrix = 0.0;
3646 | 
3647 |     const UpdateFlags uf_cell(update_values | update_gradients |
3648 | 			      update_quadrature_points | update_JxW_values);
3649 |     const UpdateFlags uf_face(update_values | update_normal_vectors |
3650 |                               update_JxW_values);
3651 | 
3652 |     PerTaskData_ASM per_task_data(m_fe.n_dofs_per_cell());
3653 |     ScratchData_ASM scratch_data(m_fe, m_qf_cell, uf_cell, m_qf_face, uf_face, solution_old);
3654 | 
3655 |     auto worker =
3656 |       [this](const typename DoFHandler<dim>::active_cell_iterator &cell,
3657 | 	     ScratchData_ASM & scratch,
3658 | 	     PerTaskData_ASM & data)
3659 |       {
3660 |         this->assemble_system_B0_one_cell(cell, scratch, data);
3661 |       };
3662 | 
3663 |     auto copier = [this](const PerTaskData_ASM &data)
3664 |       {
3665 |         this->m_constraints.distribute_local_to_global(data.m_cell_matrix,
3666 |                                                        data.m_local_dof_indices,
3667 | 						       m_tangent_matrix);
3668 |       };
3669 | 
3670 |     WorkStream::run(
3671 |       m_dof_handler.active_cell_iterators(),
3672 |       worker,
3673 |       copier,
3674 |       scratch_data,
3675 |       per_task_data);
3676 | 
3677 |     m_timer.leave_subsection();
3678 |   }
3679 | 
3680 |   template <int dim>
3681 |   void PhaseFieldMonolithicSolve<dim>::assemble_system_rhs_BFGS_parallel(const BlockVector<double> & solution_old,
3682 | 								         BlockVector<double> & system_rhs)
3683 |   {
3684 |     m_timer.enter_subsection("Assemble RHS");
3685 | 
3686 |     system_rhs = 0.0;
3687 | 
3688 |     const UpdateFlags uf_cell(update_values | update_gradients |
3689 | 			      update_quadrature_points | update_JxW_values);
3690 |     const UpdateFlags uf_face(update_values | update_normal_vectors |
3691 | 			      update_JxW_values);
3692 | 
3693 |     PerTaskData_ASM_RHS_BFGS per_task_data(m_fe.n_dofs_per_cell());
3694 |     ScratchData_ASM_RHS_BFGS scratch_data(m_fe, m_qf_cell, uf_cell, m_qf_face, uf_face, solution_old);
3695 | 
3696 |     auto worker =
3697 |       [this](const typename DoFHandler<dim>::active_cell_iterator &cell,
3698 | 	     ScratchData_ASM_RHS_BFGS & scratch,
3699 | 	     PerTaskData_ASM_RHS_BFGS & data)
3700 |       {
3701 |         this->assemble_system_rhs_BFGS_one_cell(cell, scratch, data);
3702 |       };
3703 | 
3704 |     auto copier = [this, &system_rhs](const PerTaskData_ASM_RHS_BFGS &data)
3705 |       {
3706 |         this->m_constraints.distribute_local_to_global(data.m_cell_rhs,
3707 |                                                        data.m_local_dof_indices,
3708 | 						       system_rhs);
3709 |       };
3710 | 
3711 |     WorkStream::run(
3712 |       m_dof_handler.active_cell_iterators(),
3713 |       worker,
3714 |       copier,
3715 |       scratch_data,
3716 |       per_task_data);
3717 | 
3718 |     m_timer.leave_subsection();
3719 |   }
3720 | 
3721 |   template <int dim>
3722 |   void PhaseFieldMonolithicSolve<dim>::assemble_system_rhs_BFGS_one_cell(
3723 |       const typename DoFHandler<dim>::active_cell_iterator &cell,
3724 |       ScratchData_ASM_RHS_BFGS & scratch,
3725 |       PerTaskData_ASM_RHS_BFGS & data) const
3726 |   {
3727 |     data.reset();
3728 |     scratch.reset();
3729 |     scratch.m_fe_values.reinit(cell);
3730 |     cell->get_dof_indices(data.m_local_dof_indices);
3731 | 
3732 |     scratch.m_fe_values[m_d_fe].get_function_values(
3733 |       scratch.m_solution_previous_step, scratch.m_phasefield_previous_step_cell);
3734 | 
3735 |     const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
3736 |       m_quadrature_point_history.get_data(cell);
3737 |     Assert(lqph.size() == m_n_q_points, ExcInternalError());
3738 | 
3739 |     const double time_ramp = (m_time.current() / m_time.end());
3740 |     std::vector<Tensor<1, dim>> rhs_values(m_n_q_points);
3741 | 
3742 |     right_hand_side(scratch.m_fe_values.get_quadrature_points(),
3743 | 		    rhs_values,
3744 | 		    m_parameters.m_x_component*1.0,
3745 | 		    m_parameters.m_y_component*1.0,
3746 | 		    m_parameters.m_z_component*1.0);
3747 | 
3748 |     const double delta_time = m_time.get_delta_t();
3749 | 
3750 |     for (const unsigned int q_point : scratch.m_fe_values.quadrature_point_indices())
3751 |       {
3752 |         for (const unsigned int k : scratch.m_fe_values.dof_indices())
3753 |           {
3754 |             const unsigned int k_group = m_fe.system_to_base_index(k).first.first;
3755 | 
3756 |             if (k_group == m_u_dof)
3757 |               {
3758 |                 scratch.m_Nx_disp[q_point][k] =
3759 |                   scratch.m_fe_values[m_u_fe].value(k, q_point);
3760 |                 scratch.m_grad_Nx_disp[q_point][k] =
3761 |                   scratch.m_fe_values[m_u_fe].gradient(k, q_point);
3762 |                 scratch.m_symm_grad_Nx_disp[q_point][k] =
3763 |                   symmetrize(scratch.m_grad_Nx_disp[q_point][k]);
3764 |               }
3765 |             else if (k_group == m_d_dof)
3766 |               {
3767 | 		scratch.m_Nx_phasefield[q_point][k] =
3768 | 		  scratch.m_fe_values[m_d_fe].value(k, q_point);
3769 | 		scratch.m_grad_Nx_phasefield[q_point][k] =
3770 | 		  scratch.m_fe_values[m_d_fe].gradient(k, q_point);
3771 |               }
3772 |             else
3773 |               Assert(k_group <= m_d_dof, ExcInternalError());
3774 |           }
3775 |       }
3776 | 
3777 |     for (const unsigned int q_point : scratch.m_fe_values.quadrature_point_indices())
3778 |       {
3779 | 	const double length_scale            = lqph[q_point]->get_length_scale();
3780 | 	const double gc                      = lqph[q_point]->get_critical_energy_release_rate();
3781 | 	const double eta                     = lqph[q_point]->get_viscosity();
3782 | 	const double current_positive_strain_energy = lqph[q_point]->get_current_positive_strain_energy();
3783 | 
3784 | 	const double phasefield_value        = lqph[q_point]->get_phase_field_value();
3785 | 	const Tensor<1, dim> phasefield_grad = lqph[q_point]->get_phase_field_gradient();
3786 | 
3787 |         const std::vector<double>         &      N_phasefield = scratch.m_Nx_phasefield[q_point];
3788 |         const std::vector<Tensor<1, dim>> & grad_N_phasefield = scratch.m_grad_Nx_phasefield[q_point];
3789 |         const double                old_phasefield = scratch.m_phasefield_previous_step_cell[q_point];
3790 | 
3791 |         const SymmetricTensor<2, dim> & cauchy_stress = lqph[q_point]->get_cauchy_stress();
3792 | 
3793 |         const std::vector<Tensor<1,dim>> & N_disp = scratch.m_Nx_disp[q_point];
3794 |         const std::vector<SymmetricTensor<2, dim>> & symm_grad_N_disp =
3795 |           scratch.m_symm_grad_Nx_disp[q_point];
3796 |         const double JxW = scratch.m_fe_values.JxW(q_point);
3797 | 
3798 |         SymmetricTensor<2, dim> symm_grad_Nx_i_x_C;
3799 | 
3800 |         for (const unsigned int i : scratch.m_fe_values.dof_indices())
3801 |           {
3802 |             const unsigned int i_group = m_fe.system_to_base_index(i).first.first;
3803 | 
3804 |             if (i_group == m_u_dof)
3805 |               {
3806 |                 data.m_cell_rhs(i) += (symm_grad_N_disp[i] * cauchy_stress) * JxW;
3807 | 
3808 | 		// contributions from the body force to right-hand side
3809 | 		data.m_cell_rhs(i) -= N_disp[i] * rhs_values[q_point] * JxW;
3810 |               }
3811 |             else if (i_group == m_d_dof)
3812 |               {
3813 |     	        data.m_cell_rhs(i) += (    gc * length_scale * grad_N_phasefield[i] * phasefield_grad
3814 |     	                                +  (   gc / length_scale * phasefield_value
3815 | 					     + eta / delta_time  * (phasefield_value - old_phasefield)
3816 | 					     + degradation_function_derivative(phasefield_value)
3817 | 					     * current_positive_strain_energy )
3818 | 					  * N_phasefield[i]
3819 | 				      ) * JxW;
3820 |               }
3821 |             else
3822 |               Assert(i_group <= m_d_dof, ExcInternalError());
3823 |           }  // i
3824 |       }  // q_point
3825 | 
3826 |     // if there is surface pressure, this surface pressure always applied to the
3827 |     // reference configuration
3828 |     const unsigned int face_pressure_id = 100;
3829 |     const double p0 = 0.0;
3830 | 
3831 |     for (const auto &face : cell->face_iterators())
3832 |       if (face->at_boundary() && face->boundary_id() == face_pressure_id)
3833 |         {
3834 |           scratch.m_fe_face_values.reinit(cell, face);
3835 | 
3836 |           for (const unsigned int f_q_point : scratch.m_fe_face_values.quadrature_point_indices())
3837 |             {
3838 |               const Tensor<1, dim> &N = scratch.m_fe_face_values.normal_vector(f_q_point);
3839 | 
3840 |               const double         pressure  = p0 * time_ramp;
3841 |               const Tensor<1, dim> traction  = pressure * N;
3842 | 
3843 |               for (const unsigned int i : scratch.m_fe_values.dof_indices())
3844 |                 {
3845 |                   const unsigned int i_group = m_fe.system_to_base_index(i).first.first;
3846 | 
3847 |                   if (i_group == m_u_dof)
3848 |                     {
3849 |     		      const unsigned int component_i = m_fe.system_to_component_index(i).first;
3850 |     		      const double Ni = scratch.m_fe_face_values.shape_value(i, f_q_point);
3851 |     		      const double JxW = scratch.m_fe_face_values.JxW(f_q_point);
3852 |     		      data.m_cell_rhs(i) -= (Ni * traction[component_i]) * JxW;
3853 |                     }
3854 |                 }
3855 |             }
3856 |         }
3857 |   }
3858 | 
3859 |   template <int dim>
3860 |   void PhaseFieldMonolithicSolve<dim>::assemble_system_B0_one_cell(
3861 |       const typename DoFHandler<dim>::active_cell_iterator &cell,
3862 |       ScratchData_ASM & scratch,
3863 |       PerTaskData_ASM & data) const
3864 |   {
3865 |     data.reset();
3866 |     scratch.reset();
3867 |     scratch.m_fe_values.reinit(cell);
3868 |     cell->get_dof_indices(data.m_local_dof_indices);
3869 | 
3870 |     scratch.m_fe_values[m_d_fe].get_function_values(
3871 |       scratch.m_solution_previous_step, scratch.m_phasefield_previous_step_cell);
3872 | 
3873 |     const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
3874 |       m_quadrature_point_history.get_data(cell);
3875 |     Assert(lqph.size() == m_n_q_points, ExcInternalError());
3876 | 
3877 |     const double delta_time = m_time.get_delta_t();
3878 | 
3879 |     for (const unsigned int q_point : scratch.m_fe_values.quadrature_point_indices())
3880 |       {
3881 |         for (const unsigned int k : scratch.m_fe_values.dof_indices())
3882 |           {
3883 |             const unsigned int k_group = m_fe.system_to_base_index(k).first.first;
3884 | 
3885 |             if (k_group == m_u_dof)
3886 |               {
3887 |                 scratch.m_Nx_disp[q_point][k] =
3888 |                   scratch.m_fe_values[m_u_fe].value(k, q_point);
3889 |                 scratch.m_grad_Nx_disp[q_point][k] =
3890 |                   scratch.m_fe_values[m_u_fe].gradient(k, q_point);
3891 |                 scratch.m_symm_grad_Nx_disp[q_point][k] =
3892 |                   symmetrize(scratch.m_grad_Nx_disp[q_point][k]);
3893 |               }
3894 |             else if (k_group == m_d_dof)
3895 |               {
3896 | 		scratch.m_Nx_phasefield[q_point][k] =
3897 | 		  scratch.m_fe_values[m_d_fe].value(k, q_point);
3898 | 		scratch.m_grad_Nx_phasefield[q_point][k] =
3899 | 		  scratch.m_fe_values[m_d_fe].gradient(k, q_point);
3900 |               }
3901 |             else
3902 |               Assert(k_group <= m_d_dof, ExcInternalError());
3903 |           }
3904 |       }
3905 | 
3906 |     for (const unsigned int q_point : scratch.m_fe_values.quadrature_point_indices())
3907 |       {
3908 | 	const double length_scale            = lqph[q_point]->get_length_scale();
3909 | 	const double gc                      = lqph[q_point]->get_critical_energy_release_rate();
3910 | 	const double eta                     = lqph[q_point]->get_viscosity();
3911 | 	const double current_positive_strain_energy = lqph[q_point]->get_current_positive_strain_energy();
3912 | 
3913 | 	const double phasefield_value        = lqph[q_point]->get_phase_field_value();
3914 | 
3915 |         const std::vector<double>         &      N_phasefield = scratch.m_Nx_phasefield[q_point];
3916 |         const std::vector<Tensor<1, dim>> & grad_N_phasefield = scratch.m_grad_Nx_phasefield[q_point];
3917 | 
3918 |         //const SymmetricTensor<2, dim> & cauchy_stress_positive = lqph[q_point]->get_cauchy_stress_positive();
3919 |         const SymmetricTensor<4, dim> & mechanical_C  = lqph[q_point]->get_mechanical_C();
3920 | 
3921 |         const std::vector<SymmetricTensor<2, dim>> & symm_grad_N_disp =
3922 |           scratch.m_symm_grad_Nx_disp[q_point];
3923 |         const double JxW = scratch.m_fe_values.JxW(q_point);
3924 | 
3925 |         SymmetricTensor<2, dim> symm_grad_Nx_i_x_C;
3926 | 
3927 |         for (const unsigned int i : scratch.m_fe_values.dof_indices())
3928 |           {
3929 |             const unsigned int i_group = m_fe.system_to_base_index(i).first.first;
3930 | 
3931 |             if (i_group == m_u_dof)
3932 |               {
3933 |                 symm_grad_Nx_i_x_C = symm_grad_N_disp[i] * mechanical_C;
3934 |               }
3935 | 
3936 |             for (const unsigned int j : scratch.m_fe_values.dof_indices())
3937 |               {
3938 |                 const unsigned int j_group = m_fe.system_to_base_index(j).first.first;
3939 | 
3940 |                 if ((i_group == j_group) && (i_group == m_u_dof))
3941 |                   {
3942 |                     data.m_cell_matrix(i, j) += symm_grad_Nx_i_x_C * symm_grad_N_disp[j] * JxW;
3943 |                   }
3944 |                 else if ((i_group == j_group) && (i_group == m_d_dof))
3945 |                   {
3946 |                     data.m_cell_matrix(i, j) += (  (   gc/length_scale + eta/delta_time
3947 |                 	                             + degradation_function_2nd_order_derivative(phasefield_value)
3948 | 						     * current_positive_strain_energy  )
3949 |                 	                          * N_phasefield[i] * N_phasefield[j]
3950 | 					          + gc * length_scale * grad_N_phasefield[i] * grad_N_phasefield[j]
3951 | 					        ) * JxW;
3952 |                   }
3953 |                 else
3954 |                   Assert((i_group <= m_d_dof) && (j_group <= m_d_dof),
3955 |                          ExcInternalError());
3956 |               } // j
3957 |           }  // i
3958 |       }  // q_point
3959 |   }
3960 | 
3961 |   template <int dim>
3962 |   void PhaseFieldMonolithicSolve<dim>::assemble_system_rhs_BFGS(const BlockVector<double> & solution_old,
3963 | 								BlockVector<double> & system_rhs)
3964 |   {
3965 |     m_timer.enter_subsection("Assemble RHS");
3966 | 
3967 |     system_rhs = 0.0;
3968 | 
3969 |     Vector<double> cell_rhs(m_dofs_per_cell);
3970 |     std::vector<types::global_dof_index> local_dof_indices(m_dofs_per_cell);
3971 | 
3972 |     const double time_ramp = (m_time.current() / m_time.end());
3973 |     const double delta_time = m_time.get_delta_t();
3974 | 
3975 |     std::vector<Tensor<1, dim>> rhs_values(m_n_q_points);
3976 |     const UpdateFlags uf_cell(update_values | update_gradients |
3977 | 			      update_quadrature_points | update_JxW_values);
3978 |     const UpdateFlags uf_face(update_values | update_normal_vectors |
3979 | 			      update_JxW_values);
3980 | 
3981 |     FEValues<dim> fe_values(m_fe, m_qf_cell, uf_cell);
3982 |     FEFaceValues<dim> fe_face_values(m_fe, m_qf_face, uf_face);
3983 | 
3984 |     // shape function values for displacement field
3985 |     std::vector<std::vector<Tensor<1, dim>>>
3986 |       Nx_disp(m_qf_cell.size(), std::vector<Tensor<1, dim>>(m_dofs_per_cell));
3987 |     std::vector<std::vector<Tensor<2, dim>>>
3988 |       grad_Nx_disp(m_qf_cell.size(), std::vector<Tensor<2, dim>>(m_dofs_per_cell));
3989 |     std::vector<std::vector<SymmetricTensor<2, dim>>>
3990 |       symm_grad_Nx_disp(m_qf_cell.size(), std::vector<SymmetricTensor<2, dim>>(m_dofs_per_cell));
3991 | 
3992 |     // shape function values for phase field
3993 |     std::vector<std::vector<double>>
3994 |       Nx_phasefield(m_qf_cell.size(), std::vector<double>(m_dofs_per_cell));
3995 |     std::vector<std::vector<Tensor<1, dim>>>
3996 |       grad_Nx_phasefield(m_qf_cell.size(), std::vector<Tensor<1, dim>>(m_dofs_per_cell));
3997 | 
3998 |     std::vector<double> phasefield_previous_step_cell(m_qf_cell.size());
3999 | 
4000 |     for (const auto &cell : m_dof_handler.active_cell_iterators())
4001 |       {
4002 | 	const std::vector<std::shared_ptr< PointHistory<dim>>> lqph =
4003 | 	  m_quadrature_point_history.get_data(cell);
4004 | 	Assert(lqph.size() == m_n_q_points, ExcInternalError());
4005 | 
4006 | 	cell_rhs = 0.0;
4007 | 	fe_values.reinit(cell);
4008 | 	right_hand_side(fe_values.get_quadrature_points(),
4009 | 			rhs_values,
4010 | 			m_parameters.m_x_component*time_ramp,
4011 | 			m_parameters.m_y_component*time_ramp,
4012 | 			m_parameters.m_z_component*time_ramp);
4013 | 
4014 | 	fe_values[m_d_fe].get_function_values(
4015 | 	    solution_old, phasefield_previous_step_cell);
4016 | 
4017 | 	for (const unsigned int q_point : fe_values.quadrature_point_indices())
4018 | 	  {
4019 | 	    for (const unsigned int k : fe_values.dof_indices())
4020 | 	      {
4021 | 		const unsigned int k_group = m_fe.system_to_base_index(k).first.first;
4022 | 
4023 | 		if (k_group == m_u_dof)
4024 | 		  {
4025 | 		    Nx_disp[q_point][k] = fe_values[m_u_fe].value(k, q_point);
4026 | 		    grad_Nx_disp[q_point][k] = fe_values[m_u_fe].gradient(k, q_point);
4027 | 		    symm_grad_Nx_disp[q_point][k] = symmetrize(grad_Nx_disp[q_point][k]);
4028 | 		  }
4029 | 		else if (k_group == m_d_dof)
4030 | 		  {
4031 | 		    Nx_phasefield[q_point][k] = fe_values[m_d_fe].value(k, q_point);
4032 | 		    grad_Nx_phasefield[q_point][k] = fe_values[m_d_fe].gradient(k, q_point);
4033 | 		  }
4034 | 		else
4035 | 		  Assert(k_group <= m_d_dof, ExcInternalError());
4036 | 	      }
4037 | 	  }
4038 | 
4039 | 	for (const unsigned int q_point : fe_values.quadrature_point_indices())
4040 | 	  {
4041 | 	    const double length_scale            = lqph[q_point]->get_length_scale();
4042 | 	    const double gc                      = lqph[q_point]->get_critical_energy_release_rate();
4043 | 	    const double eta                     = lqph[q_point]->get_viscosity();
4044 | 	    const double current_positive_strain_energy = lqph[q_point]->get_current_positive_strain_energy();
4045 | 
4046 | 	    const double phasefield_value        = lqph[q_point]->get_phase_field_value();
4047 | 	    const Tensor<1, dim> phasefield_grad = lqph[q_point]->get_phase_field_gradient();
4048 | 
4049 | 	    const std::vector<double>         &      N_phasefield = Nx_phasefield[q_point];
4050 | 	    const std::vector<Tensor<1, dim>> & grad_N_phasefield = grad_Nx_phasefield[q_point];
4051 | 	    const double                old_phasefield = phasefield_previous_step_cell[q_point];
4052 | 
4053 | 	    const SymmetricTensor<2, dim> & cauchy_stress = lqph[q_point]->get_cauchy_stress();
4054 | 
4055 | 	    const std::vector<Tensor<1,dim>> & N = Nx_disp[q_point];
4056 | 	    const std::vector<SymmetricTensor<2, dim>> & symm_grad_N = symm_grad_Nx_disp[q_point];
4057 | 	    const double JxW = fe_values.JxW(q_point);
4058 | 
4059 | 	    for (const unsigned int i : fe_values.dof_indices())
4060 | 	      {
4061 | 		const unsigned int i_group = m_fe.system_to_base_index(i).first.first;
4062 | 
4063 | 		if (i_group == m_u_dof)
4064 | 		  {
4065 | 		    cell_rhs(i) += (symm_grad_N[i] * cauchy_stress) * JxW;
4066 | 		    // contributions from the body force to right-hand side
4067 | 		    cell_rhs(i) -= N[i] * rhs_values[q_point] * JxW;
4068 | 		  }
4069 | 		else if (i_group == m_d_dof)
4070 | 		  {
4071 | 		    cell_rhs(i) += (    gc * length_scale * grad_N_phasefield[i] * phasefield_grad
4072 | 	    	                     +  (   gc / length_scale * phasefield_value
4073 | 			                  + eta / delta_time  * (phasefield_value - old_phasefield)
4074 | 				          + degradation_function_derivative(phasefield_value)
4075 | 					  * current_positive_strain_energy )
4076 | 				     * N_phasefield[i]
4077 | 				   ) * JxW;
4078 | 		  }
4079 | 		else
4080 | 		  Assert(i_group <= m_d_dof, ExcInternalError());
4081 | 	      }
4082 | 	  }
4083 | 
4084 | 	// if there is surface pressure, this surface pressure always applied to the
4085 | 	// reference configuration
4086 | 	const unsigned int face_pressure_id = 100;
4087 | 	const double p0 = 0.0;
4088 | 
4089 | 	for (const auto &face : cell->face_iterators())
4090 | 	  {
4091 | 	    if (face->at_boundary() && face->boundary_id() == face_pressure_id)
4092 | 	      {
4093 | 		fe_face_values.reinit(cell, face);
4094 | 
4095 | 		for (const unsigned int f_q_point : fe_face_values.quadrature_point_indices())
4096 | 		  {
4097 | 		    const Tensor<1, dim> &N = fe_face_values.normal_vector(f_q_point);
4098 | 
4099 | 		    const double         pressure  = p0 * time_ramp;
4100 | 		    const Tensor<1, dim> traction  = pressure * N;
4101 | 
4102 | 		    for (const unsigned int i : fe_values.dof_indices())
4103 | 		      {
4104 | 			const unsigned int i_group = m_fe.system_to_base_index(i).first.first;
4105 | 
4106 | 			if (i_group == m_u_dof)
4107 | 			  {
4108 | 			    const unsigned int component_i = m_fe.system_to_component_index(i).first;
4109 | 			    const double Ni = fe_face_values.shape_value(i, f_q_point);
4110 | 			    const double JxW = fe_face_values.JxW(f_q_point);
4111 | 			    cell_rhs(i) -= (Ni * traction[component_i]) * JxW;
4112 | 			  }
4113 | 		      }
4114 | 		  }
4115 | 	      }
4116 | 	  }
4117 | 
4118 | 	cell->get_dof_indices(local_dof_indices);
4119 | 	for (const unsigned int i : fe_values.dof_indices())
4120 | 	  system_rhs(local_dof_indices[i]) += cell_rhs(i);
4121 |       } // for (const auto &cell : m_dof_handler.active_cell_iterators())
4122 | 
4123 |     m_timer.leave_subsection();
4124 |   }
4125 | 
4126 | 
4127 | 
4128 |   template <int dim>
4129 |   double PhaseFieldMonolithicSolve<dim>::line_search_stepsize_gradient_based(const BlockVector<double> & BFGS_p_vector,
4130 | 				                                             const BlockVector<double> & solution_delta)
4131 |   {
4132 |     BlockVector<double> g_old(m_system_rhs);
4133 | 
4134 |     // BFGS_p_vector is the search direction
4135 |     BlockVector<double> solution_delta_trial(solution_delta);
4136 |     // take a full step size 1.0
4137 |     solution_delta_trial.add(1.0, BFGS_p_vector);
4138 | 
4139 |     update_qph_incremental(solution_delta_trial, m_solution);
4140 | 
4141 |     BlockVector<double> g_new(m_dofs_per_block);
4142 |     assemble_system_rhs_BFGS_parallel(m_solution, g_new);
4143 | 
4144 |     BlockVector<double> y_old(m_dofs_per_block);
4145 | 
4146 |     y_old = g_new - g_old;
4147 | 
4148 |     double alpha = 1.0;
4149 | 
4150 |     double alpha_old = 0.0;
4151 | 
4152 |     double delta_alpha_old = alpha - alpha_old;
4153 | 
4154 |     double delta_alpha_new;
4155 | 
4156 |     unsigned int ls_max = 10;
4157 | 
4158 |     unsigned int i = 1;
4159 | 
4160 |     for (; i <= ls_max; ++i)
4161 |       {
4162 | 	delta_alpha_new = -delta_alpha_old
4163 | 	                * (g_new * BFGS_p_vector)/(y_old * BFGS_p_vector);
4164 | 	alpha += delta_alpha_new;
4165 | 
4166 | 	if (std::fabs(delta_alpha_new) < 1.0e-5)
4167 | 	  break;
4168 | 
4169 |         if (i == ls_max)
4170 |           {
4171 |             alpha = 1.0;
4172 |             break;
4173 |           }
4174 | 
4175 |         g_old = g_new;
4176 | 
4177 |         // BFGS_p_vector is the search direction
4178 |         solution_delta_trial = solution_delta;
4179 |         solution_delta_trial.add(alpha, BFGS_p_vector);
4180 |         update_qph_incremental(solution_delta_trial, m_solution);
4181 |         assemble_system_rhs_BFGS_parallel(m_solution, g_new);
4182 | 
4183 |         y_old = g_new - g_old;
4184 | 
4185 |         delta_alpha_old = delta_alpha_new;
4186 |       }
4187 | 
4188 |     if (alpha < 1.0e-3)
4189 |       alpha = 1.0;
4190 | 
4191 |     //num_ls = i;
4192 |     return alpha;
4193 |   }
4194 | 
4195 |   template <int dim>
4196 |   double PhaseFieldMonolithicSolve<dim>::line_search_stepsize_strong_wolfe(const double phi_0,
4197 | 				                                           const double phi_0_prime,
4198 | 				                                           const BlockVector<double> & BFGS_p_vector,
4199 | 				                                           const BlockVector<double> & solution_delta)
4200 |   {
4201 |     //AssertThrow(phi_0_prime < 0,
4202 |     //            ExcMessage("The derivative of phi at alpha = 0 should be negative!"));
4203 | 
4204 |     // Some line search parameters
4205 |     const double c1 = 0.0001;
4206 |     const double c2 = 0.9;
4207 |     const double alpha_max = 1.0;
4208 |     const unsigned int max_iter = 20;
4209 |     double alpha = 1.0;
4210 | 
4211 |     double phi_old = phi_0;
4212 |     double phi_prime_old = phi_0_prime;
4213 |     double alpha_old = 0.0;
4214 | 
4215 |     double phi, phi_prime;
4216 | 
4217 |     std::pair<double, double> current_phi_phi_prime;
4218 | 
4219 |     unsigned int i = 0;
4220 |     for (; i < max_iter; ++i)
4221 |       {
4222 | 	current_phi_phi_prime = calculate_phi_and_phi_prime(alpha, BFGS_p_vector, solution_delta);
4223 | 	phi = current_phi_phi_prime.first;
4224 | 	phi_prime = current_phi_phi_prime.second;
4225 | 
4226 | 	if (   ( phi > (phi_0 + c1 * alpha * phi_0_prime) )
4227 | 	    || ( i > 0 && phi > phi_old ) )
4228 | 	  {
4229 | 	    return line_search_zoom_strong_wolfe(phi_old, phi_prime_old, alpha_old,
4230 | 						 phi,     phi_prime,     alpha,
4231 | 						 phi_0,   phi_0_prime,   BFGS_p_vector,
4232 | 						 c1,      c2,            max_iter, solution_delta);
4233 | 	  }
4234 | 
4235 | 	if (std::fabs(phi_prime) <= c2 * std::fabs(phi_0_prime))
4236 | 	  {
4237 | 	    return alpha;
4238 | 	  }
4239 | 
4240 | 	if (phi_prime >= 0)
4241 | 	  {
4242 | 	    return line_search_zoom_strong_wolfe(phi,     phi_prime,     alpha,
4243 | 						 phi_old, phi_prime_old, alpha_old,
4244 | 						 phi_0,   phi_0_prime,   BFGS_p_vector,
4245 | 						 c1,      c2,            max_iter, solution_delta);
4246 | 	  }
4247 | 
4248 | 	phi_old = phi;
4249 | 	phi_prime_old = phi_prime;
4250 | 	alpha_old = alpha;
4251 | 
4252 | 	alpha = std::min(0.6*alpha, alpha_max);
4253 |       }
4254 | 
4255 |     return alpha;
4256 |   }
4257 | 
4258 |   template <int dim>
4259 |   double PhaseFieldMonolithicSolve<dim>::
4260 |     line_search_zoom_strong_wolfe(double phi_low, double phi_low_prime, double alpha_low,
4261 | 				  double phi_high, double phi_high_prime, double alpha_high,
4262 | 				  double phi_0, double phi_0_prime, const BlockVector<double> & BFGS_p_vector,
4263 | 				  double c1, double c2, unsigned int max_iter, const BlockVector<double> & solution_delta)
4264 |   {
4265 |     double alpha = 0;
4266 |     std::pair<double, double> current_phi_phi_prime;
4267 |     double phi, phi_prime;
4268 | 
4269 |     unsigned int i = 0;
4270 |     for (; i < max_iter; ++i)
4271 |       {
4272 | 	// a simple bisection is faster than cubic interpolation
4273 | 	alpha = 0.5 * (alpha_low + alpha_high);
4274 | 	//alpha = line_search_interpolation_cubic(alpha_low, phi_low, phi_low_prime,
4275 | 	//					alpha_high, phi_high, phi_high_prime);
4276 | 	current_phi_phi_prime = calculate_phi_and_phi_prime(alpha, BFGS_p_vector, solution_delta);
4277 | 	phi = current_phi_phi_prime.first;
4278 | 	phi_prime = current_phi_phi_prime.second;
4279 | 
4280 | 	if (   (phi > phi_0 + c1 * alpha * phi_0_prime)
4281 | 	    || (phi > phi_low) )
4282 | 	  {
4283 | 	    alpha_high = alpha;
4284 | 	    phi_high = phi;
4285 | 	    phi_high_prime = phi_prime;
4286 | 	  }
4287 | 	else
4288 | 	  {
4289 | 	    if (std::fabs(phi_prime) <= c2 * std::fabs(phi_0_prime))
4290 | 	      {
4291 | 		return alpha;
4292 | 	      }
4293 | 
4294 | 	    if (phi_prime * (alpha_high - alpha_low) >= 0.0)
4295 | 	      {
4296 | 		alpha_high = alpha_low;
4297 | 		phi_high_prime = phi_low_prime;
4298 | 		phi_high = phi_low;
4299 | 	      }
4300 | 
4301 | 	    alpha_low = alpha;
4302 | 	    phi_low_prime = phi_prime;
4303 | 	    phi_low = phi;
4304 | 	  }
4305 |       }
4306 | 
4307 |     // avoid unused variable warnings from compiler
4308 |     (void)phi_high;
4309 |     (void)phi_high_prime;
4310 |     return alpha;
4311 |   }
4312 | 
4313 |   template <int dim>
4314 |   double PhaseFieldMonolithicSolve<dim>::
4315 |     line_search_interpolation_cubic(const double alpha_0, const double phi_0, const double phi_0_prime,
4316 |   			            const double alpha_1, const double phi_1, const double phi_1_prime)
4317 |   {
4318 |     const double d1 = phi_0_prime + phi_1_prime - 3.0 * (phi_0 - phi_1) / (alpha_0 - alpha_1);
4319 | 
4320 |     const double temp = d1 * d1 - phi_0_prime * phi_1_prime;
4321 | 
4322 |     if (temp < 0.0)
4323 |       return 0.5 * (alpha_0 + alpha_1);
4324 | 
4325 |     int sign;
4326 |     if (alpha_1 > alpha_0)
4327 |       sign = 1;
4328 |     else
4329 |       sign = -1;
4330 | 
4331 |     const double d2 = sign * std::sqrt(temp);
4332 | 
4333 |     const double alpha = alpha_1 - (alpha_1 - alpha_0)
4334 | 	               * (phi_1_prime + d2 - d1) / (phi_1_prime - phi_0_prime + 2*d2);
4335 | 
4336 |     if (    (alpha_1 > alpha_0)
4337 | 	 && (alpha > alpha_1 || alpha < alpha_0))
4338 |       return 0.5 * (alpha_0 + alpha_1);
4339 | 
4340 |     if (    (alpha_0 > alpha_1)
4341 | 	 && (alpha > alpha_0 || alpha < alpha_1))
4342 |       return 0.5 * (alpha_0 + alpha_1);
4343 | 
4344 |     return alpha;
4345 |   }
4346 | 
4347 |   template <int dim>
4348 |   std::pair<double, double> PhaseFieldMonolithicSolve<dim>::
4349 |     calculate_phi_and_phi_prime(const double alpha,
4350 | 				const BlockVector<double> & BFGS_p_vector,
4351 | 				const BlockVector<double> & solution_delta)
4352 |   {
4353 |     // the first component is phi(alpha), the second component is phi_prime(alpha),
4354 |     std::pair<double, double> phi_values;
4355 | 
4356 |     BlockVector<double> solution_delta_trial(solution_delta);
4357 |     solution_delta_trial.add(alpha, BFGS_p_vector);
4358 | 
4359 |     update_qph_incremental(solution_delta_trial, m_solution);
4360 | 
4361 |     BlockVector<double> system_rhs(m_dofs_per_block);
4362 |     assemble_system_rhs_BFGS_parallel(m_solution, system_rhs);
4363 |     //m_constraints.condense(system_rhs);
4364 | 
4365 |     phi_values.first = calculate_energy_functional();
4366 |     phi_values.second = system_rhs * BFGS_p_vector;
4367 |     return phi_values;
4368 |   }
4369 | 
4370 |   template <int dim>
4371 |   void PhaseFieldMonolithicSolve<dim>::LBFGS_B0(BlockVector<double> & LBFGS_r_vector,
4372 | 						BlockVector<double> & LBFGS_q_vector)
4373 |   {
4374 |     m_timer.enter_subsection("Solve B0");
4375 | 
4376 |     assemble_system_B0(m_solution);
4377 | 
4378 |     if (m_parameters.m_type_linear_solver == "Direct")
4379 |       {
4380 | 	SparseDirectUMFPACK A_direct;
4381 | 	A_direct.initialize(m_tangent_matrix);
4382 | 	A_direct.vmult(LBFGS_r_vector,
4383 | 		       LBFGS_q_vector);
4384 |       }
4385 |     else if (m_parameters.m_type_linear_solver == "CG")
4386 |       {
4387 | /*
4388 | 	SolverControl            solver_control(1e6, 1e-9);
4389 | 	SolverCG<BlockVector<double>> cg(solver_control);
4390 | 
4391 | 	PreconditionJacobi<BlockSparseMatrix<double>> preconditioner;
4392 | 	preconditioner.initialize(m_tangent_matrix, 1.0);
4393 | 
4394 | 	cg.solve(m_tangent_matrix,
4395 | 		 LBFGS_r_vector,
4396 | 		 LBFGS_q_vector,
4397 | 		 preconditioner);
4398 | */
4399 | 	SolverControl            solver_control_uu(1e6, 1e-9);
4400 | 	SolverCG<Vector<double>> cg_uu(solver_control_uu);
4401 | 
4402 | 	PreconditionJacobi<SparseMatrix<double>> preconditioner_uu;
4403 | 	preconditioner_uu.initialize(m_tangent_matrix.block(m_u_dof, m_u_dof), 1.0);
4404 | 	cg_uu.solve(m_tangent_matrix.block(m_u_dof, m_u_dof),
4405 | 	            LBFGS_r_vector.block(m_u_dof),
4406 | 	            LBFGS_q_vector.block(m_u_dof),
4407 | 	            preconditioner_uu);
4408 | 
4409 | 	SolverControl            solver_control_dd(1e6, 1e-15);
4410 | 	SolverCG<Vector<double>> cg_dd(solver_control_dd);
4411 | 
4412 | 	PreconditionJacobi<SparseMatrix<double>> preconditioner_dd;
4413 | 	preconditioner_dd.initialize(m_tangent_matrix.block(m_d_dof, m_d_dof), 1.0);
4414 | 	cg_dd.solve(m_tangent_matrix.block(m_d_dof, m_d_dof),
4415 | 	            LBFGS_r_vector.block(m_d_dof),
4416 | 	            LBFGS_q_vector.block(m_d_dof),
4417 | 	            preconditioner_dd);
4418 |       }
4419 |     else
4420 |       {
4421 | 	AssertThrow(false,
4422 | 	            ExcMessage("Selected linear solver not implemented!"));
4423 |       }
4424 | 
4425 |     m_timer.leave_subsection();
4426 |   }
4427 | 
4428 |   template <int dim>
4429 |   void PhaseFieldMonolithicSolve<dim>::print_conv_header_LBFGS()
4430 |   {
4431 |     m_logfile << "\t\t" << "L-BFGS (warning: without phasefield irreversibility)" << std::endl;;
4432 |     static const unsigned int l_width = 100;
4433 |     m_logfile << '\t' << '\t';
4434 |     for (unsigned int i = 0; i < l_width; ++i)
4435 |       m_logfile << '_';
4436 |     m_logfile << std::endl;
4437 | 
4438 |     m_logfile << "\t\t itr "
4439 |               << " |  LS-alpha     Energy      Res_Norm    "
4440 |               << " Res_u      Res_d    Inc_Norm   "
4441 |               << " Inc_u      Inc_d" << std::endl;
4442 | 
4443 |     m_logfile << '\t' << '\t';
4444 |     for (unsigned int i = 0; i < l_width; ++i)
4445 |       m_logfile << '_';
4446 |     m_logfile << std::endl;
4447 |   }
4448 | 
4449 |   template <int dim>
4450 |   void PhaseFieldMonolithicSolve<dim>::print_conv_header_LBFGSB()
4451 |   {
4452 |     m_logfile << '\t' << "L-BFGS-B" << std::endl;;
4453 |     static const unsigned int l_width = 130;
4454 |     m_logfile << '\t';
4455 |     for (unsigned int i = 0; i < l_width; ++i)
4456 |       m_logfile << '_';
4457 |     m_logfile << std::endl;
4458 | 
4459 |     m_logfile << "\t itr  "
4460 |               << " |    LB    UB    LUB    CG-itr    LS-alpha     Energy      Res_Norm    "
4461 |               << " Res_u      Res_d    Inc_Norm   "
4462 |               << " Inc_u      Inc_d" << std::endl;
4463 | 
4464 |     m_logfile << '\t';
4465 |     for (unsigned int i = 0; i < l_width; ++i)
4466 |       m_logfile << '_';
4467 |     m_logfile << std::endl;
4468 |   }
4469 | 
4470 |   template <int dim>
4471 |   void PhaseFieldMonolithicSolve<dim>::
4472 |   solve_nonlinear_timestep_LBFGS(BlockVector<double> & solution_delta,
4473 | 				 BlockVector<double> & LBFGS_update_refine)
4474 |   {
4475 |     BlockVector<double> LBFGS_update(m_dofs_per_block);
4476 | 
4477 |     LBFGS_update = 0.0;
4478 | 
4479 |     m_error_residual.reset();
4480 |     m_error_residual_0.reset();
4481 |     m_error_residual_norm.reset();
4482 |     m_error_update.reset();
4483 |     m_error_update_0.reset();
4484 |     m_error_update_norm.reset();
4485 | 
4486 |     if (m_parameters.m_output_iteration_history)
4487 |       print_conv_header_LBFGS();
4488 | 
4489 |     unsigned int LBFGS_iteration = 0;
4490 | 
4491 |     BlockVector<double> LBFGS_r_vector(m_dofs_per_block);
4492 |     BlockVector<double> LBFGS_y_vector(m_dofs_per_block);
4493 |     BlockVector<double> LBFGS_q_vector(m_dofs_per_block);
4494 |     BlockVector<double> LBFGS_s_vector(m_dofs_per_block);
4495 |     std::list<std::pair< std::pair<BlockVector<double>,
4496 |                                    BlockVector<double>>,
4497 |                          double>> LBFGS_vector_list;
4498 | 
4499 |     const unsigned int LBFGS_m = m_parameters.m_LBFGS_m;
4500 |     std::list<double> LBFGS_alpha_list;
4501 | 
4502 |     double line_search_parameter = 0.0;
4503 |     double LBFGS_beta = 0.0;
4504 |     double rho = 0.0;
4505 | 
4506 |     for (; LBFGS_iteration < m_parameters.m_max_iterations_BFGS; ++LBFGS_iteration)
4507 |       {
4508 | 	if (m_parameters.m_output_iteration_history)
4509 | 	  m_logfile << '\t' << '\t' << std::setw(4) << LBFGS_iteration << ' '
4510 |                     << std::flush;
4511 | 
4512 |         make_constraints(LBFGS_iteration);
4513 | 
4514 |         // At the first step, we simply distribute the inhomogeneous part of
4515 |         // the constraints
4516 |         if (LBFGS_iteration == 0)
4517 |           {
4518 |             // use the solution from the previous solve on the
4519 |             // refined mesh as initial guess
4520 |             LBFGS_update = LBFGS_update_refine;
4521 | 
4522 |             m_constraints.distribute(LBFGS_update);
4523 |             solution_delta += LBFGS_update;
4524 | 
4525 |             update_qph_incremental(solution_delta, m_solution);
4526 |             if (m_parameters.m_output_iteration_history)
4527 |               {
4528 |                 m_logfile << " | " << std::flush;
4529 |                 m_logfile << std::endl;
4530 |               }
4531 |             continue;
4532 |           }
4533 |         else if (LBFGS_iteration == 1)
4534 |           {
4535 | 	    // Calculate the residual vector r. NOTICE that in the context of
4536 | 	    // BFGS, this r is the gradient of the energy functional (objective function),
4537 | 	    // NOT the negative gradient of the energy functional
4538 | 	    assemble_system_rhs_BFGS_parallel(m_solution, m_system_rhs);
4539 | 
4540 | 	    // We cannot simply zero out the dofs that are constrained, since we might
4541 | 	    // have hanging node constraints. In this case, we need to modify the RHS
4542 | 	    // as C^T * b, which C contains entries of 0.5 (x_3 = 0.5*x_1 + 0.5*x_2)
4543 | 	    //for (unsigned int i = 0; i < m_dof_handler.n_dofs(); ++i)
4544 | 	      //if (m_constraints.is_constrained(i))
4545 | 		//m_system_rhs(i) = 0.0;
4546 | 
4547 | 	    // if m_constraints has inhomogeneity, we cannot call m_constraints.condense(m_system_rhs),
4548 | 	    // since the m_system_matrix needs to be provided to modify the RHS properly. However, this
4549 | 	    // error will not be detected in the release mode and only will be detected on the debug mode
4550 | 	    // if we use assemble_system_rhs_BFGS_parallel, then condense() is not necessary
4551 | 	    //m_constraints.condense(m_system_rhs);
4552 |           }
4553 | 
4554 |         get_error_residual(m_error_residual);
4555 |         if (LBFGS_iteration == 1)
4556 |           m_error_residual_0 = m_error_residual;
4557 | 
4558 |         m_error_residual_norm = m_error_residual;
4559 |         // For three-point bending problem and 3D problem, we use absolute residual
4560 |         // for convergence test
4561 |         if (m_parameters.m_relative_residual)
4562 |           m_error_residual_norm.normalize(m_error_residual_0);
4563 | 
4564 |         if (LBFGS_iteration > 1 && m_error_update_norm.m_u <= m_parameters.m_tol_u_incr
4565 |                                 && m_error_residual_norm.m_u <= m_parameters.m_tol_u_residual
4566 | 			        && m_error_update_norm.m_d <= m_parameters.m_tol_d_incr
4567 | 			        && m_error_residual_norm.m_d <= m_parameters.m_tol_d_residual
4568 | 				)
4569 |           {
4570 |             if (m_parameters.m_output_iteration_history)
4571 |               {
4572 | 		m_logfile << " | ";
4573 | 		m_logfile << " CONVERGED! " << std::fixed << std::setprecision(3) << std::setw(7)
4574 | 			  << std::scientific
4575 | 		      << "    ----    "
4576 | 		      << "  " << m_error_residual_norm.m_norm
4577 | 		      << "  " << m_error_residual_norm.m_u
4578 | 		      << "  " << m_error_residual_norm.m_d
4579 | 		      << "  " << m_error_update_norm.m_norm
4580 | 		      << "  " << m_error_update_norm.m_u
4581 | 		      << "  " << m_error_update_norm.m_d
4582 | 		      << "  " << std::endl;
4583 | 
4584 | 		m_logfile << '\t' << '\t';
4585 | 		for (unsigned int i = 0; i < 100; ++i)
4586 | 		  m_logfile << '_';
4587 | 		m_logfile << std::endl;
4588 |               }
4589 | 
4590 |             m_logfile << "\t\tConvergence is reached after "
4591 |         	      << LBFGS_iteration << " L-BFGS iterations."<< std::endl;
4592 | 
4593 |             m_logfile << "\t\tResidual information of convergence:" << std::endl;
4594 | 
4595 |             if (m_parameters.m_relative_residual)
4596 |               {
4597 | 		m_logfile << "\t\t\tRelative residual of disp. equation: "
4598 | 			  << m_error_residual_norm.m_u << std::endl;
4599 | 
4600 | 		m_logfile << "\t\t\tAbsolute residual of disp. equation: "
4601 | 			  << m_error_residual_norm.m_u * m_error_residual_0.m_u << std::endl;
4602 | 
4603 | 		m_logfile << "\t\t\tRelative residual of phasefield equation: "
4604 | 			  << m_error_residual_norm.m_d << std::endl;
4605 | 
4606 | 		m_logfile << "\t\t\tAbsolute residual of phasefield equation: "
4607 | 			  << m_error_residual_norm.m_d * m_error_residual_0.m_d << std::endl;
4608 | 
4609 | 		m_logfile << "\t\t\tRelative increment of disp.: "
4610 | 			  << m_error_update_norm.m_u << std::endl;
4611 | 
4612 | 		m_logfile << "\t\t\tAbsolute increment of disp.: "
4613 | 			  << m_error_update_norm.m_u * m_error_update_0.m_u << std::endl;
4614 | 
4615 | 		m_logfile << "\t\t\tRelative increment of phasefield: "
4616 | 			  << m_error_update_norm.m_d << std::endl;
4617 | 
4618 | 		m_logfile << "\t\t\tAbsolute increment of phasefield: "
4619 | 			  << m_error_update_norm.m_d * m_error_update_0.m_d << std::endl;
4620 |               }
4621 |             else
4622 |               {
4623 | 		m_logfile << "\t\t\tAbsolute residual of disp. equation: "
4624 | 			  << m_error_residual_norm.m_u << std::endl;
4625 | 
4626 | 		m_logfile << "\t\t\tAbsolute residual of phasefield equation: "
4627 | 			  << m_error_residual_norm.m_d << std::endl;
4628 | 
4629 | 		m_logfile << "\t\t\tAbsolute increment of disp.: "
4630 | 			  << m_error_update_norm.m_u << std::endl;
4631 | 
4632 | 		m_logfile << "\t\t\tAbsolute increment of phasefield: "
4633 | 			  << m_error_update_norm.m_d << std::endl;
4634 |               }
4635 | 
4636 |             break;
4637 |           }
4638 | 
4639 |         // LBFGS algorithm
4640 |         LBFGS_q_vector = m_system_rhs;
4641 | 
4642 |         LBFGS_alpha_list.clear();
4643 |         for (auto itr = LBFGS_vector_list.begin(); itr != LBFGS_vector_list.end(); ++itr)
4644 |           {
4645 |             LBFGS_s_vector = (itr->first).first;
4646 |             LBFGS_y_vector = (itr->first).second;
4647 |             rho = itr->second;
4648 | 
4649 |             const double alpha = rho * (LBFGS_s_vector * LBFGS_q_vector);
4650 |             LBFGS_alpha_list.push_back(alpha);
4651 | 
4652 |             LBFGS_q_vector.add(-alpha, LBFGS_y_vector);
4653 |           }
4654 | /*
4655 |         double scale_gamma = 0.0;
4656 |         if (LBFGS_iteration == 1)
4657 |           {
4658 |             scale_gamma = 1.0;
4659 |           }
4660 |         else
4661 |           {
4662 |             LBFGS_s_vector = LBFGS_vector_list.front().first.first;
4663 |             LBFGS_y_vector = LBFGS_vector_list.front().first.second;
4664 |             scale_gamma = (LBFGS_s_vector * LBFGS_y_vector)/(LBFGS_y_vector * LBFGS_y_vector);
4665 |           }
4666 | 
4667 |         LBFGS_q_vector *= scale_gamma;
4668 |         LBFGS_r_vector = LBFGS_q_vector;
4669 | */
4670 |         LBFGS_B0(LBFGS_r_vector,
4671 | 		 LBFGS_q_vector);
4672 | 
4673 |         for (auto itr = LBFGS_vector_list.rbegin(); itr != LBFGS_vector_list.rend(); ++itr)
4674 |           {
4675 |             LBFGS_s_vector = (itr->first).first;
4676 |             LBFGS_y_vector = (itr->first).second;
4677 |             rho = itr->second;
4678 | 
4679 |             LBFGS_beta = rho * (LBFGS_y_vector * LBFGS_r_vector);
4680 | 
4681 |             const double alpha = LBFGS_alpha_list.back();
4682 |             LBFGS_alpha_list.pop_back();
4683 | 
4684 |             LBFGS_r_vector.add(alpha - LBFGS_beta, LBFGS_s_vector);
4685 |           }
4686 | 
4687 |         LBFGS_r_vector *= -1.0; // this is the p_vector (search direction)
4688 | 
4689 |         m_constraints.distribute(LBFGS_r_vector);
4690 | 
4691 |         // We need a line search algorithm to decide line_search_parameter
4692 | 
4693 |         if(m_parameters.m_type_line_search == "StrongWolfe")
4694 |           {
4695 |             const double phi_0 = calculate_energy_functional();
4696 |             const double phi_0_prime = m_system_rhs * LBFGS_r_vector;
4697 | 
4698 |             line_search_parameter = line_search_stepsize_strong_wolfe(phi_0,
4699 | 		    				                      phi_0_prime,
4700 | 								      LBFGS_r_vector,
4701 | 						                      solution_delta);
4702 |           }
4703 |         else if(m_parameters.m_type_line_search == "GradientBased")
4704 |           {
4705 | 	    // LBFGS_r_vector is the search direction
4706 | 	    line_search_parameter = line_search_stepsize_gradient_based(LBFGS_r_vector,
4707 | 									solution_delta);
4708 |           }
4709 |         else
4710 |           {
4711 |             Assert(false, ExcMessage("An unknown line search method is called!"));
4712 |           }
4713 | 
4714 |         LBFGS_r_vector *= line_search_parameter;
4715 |         LBFGS_update = LBFGS_r_vector;
4716 | 
4717 |         get_error_update(LBFGS_update, m_error_update);
4718 |         if (LBFGS_iteration == 1)
4719 |           m_error_update_0 = m_error_update;
4720 | 
4721 |         m_error_update_norm = m_error_update;
4722 |         // For three-point bending problem and the sphere inclusion problem,
4723 |         // we use absolute residual for convergence test
4724 |         if (m_parameters.m_relative_residual)
4725 |           m_error_update_norm.normalize(m_error_update_0);
4726 | 
4727 |         solution_delta += LBFGS_update;
4728 |         update_qph_incremental(solution_delta, m_solution);
4729 | 
4730 |         LBFGS_y_vector = m_system_rhs;
4731 |         LBFGS_y_vector *= -1.0;
4732 |         assemble_system_rhs_BFGS_parallel(m_solution, m_system_rhs);
4733 |         // if we use assemble_system_rhs_BFGS_parallel, then condense() is not necessary
4734 |         //m_constraints.condense(m_system_rhs);
4735 |         LBFGS_y_vector += m_system_rhs;
4736 | 
4737 |         LBFGS_s_vector = LBFGS_update;
4738 | 
4739 |         if (LBFGS_iteration > LBFGS_m)
4740 |           LBFGS_vector_list.pop_back();
4741 | 
4742 |         rho = 1.0 / (LBFGS_y_vector * LBFGS_s_vector);
4743 | 
4744 |         LBFGS_vector_list.push_front(std::make_pair(std::make_pair(LBFGS_s_vector,
4745 | 								   LBFGS_y_vector),
4746 | 						    rho));
4747 |         if (m_parameters.m_output_iteration_history)
4748 |           {
4749 | 	    const double energy_functional = calculate_energy_functional();
4750 | 
4751 | 	    m_logfile << " | " << std::fixed << std::setprecision(3) << std::setw(1)
4752 | 		      << std::scientific
4753 | 		      << "" << line_search_parameter
4754 | 		      << std::fixed << std::setprecision(6) << std::setw(1)
4755 | 					<< std::scientific
4756 | 		      << "  " << energy_functional
4757 | 		      << std::fixed << std::setprecision(3) << std::setw(1)
4758 | 					<< std::scientific
4759 | 		      << "  " << m_error_residual_norm.m_norm
4760 | 		      << "  " << m_error_residual_norm.m_u
4761 | 		      << "  " << m_error_residual_norm.m_d
4762 | 		      << "  " << m_error_update_norm.m_norm
4763 | 		      << "  " << m_error_update_norm.m_u
4764 | 		      << "  " << m_error_update_norm.m_d
4765 | 		      << "  " << std::endl;
4766 |           }
4767 |       }
4768 | 
4769 |     AssertThrow(LBFGS_iteration < m_parameters.m_max_iterations_BFGS,
4770 |                 ExcMessage("No convergence in L-BFGS nonlinear solver!"));
4771 |   }
4772 | 
4773 |   template <int dim>
4774 |   void PhaseFieldMonolithicSolve<dim>::
4775 |   calculate_cauchy_point(const BlockSparseMatrix<double> & B0_matrix,
4776 | 	                 const std::list<BlockVector<double>> & y_vector_list,
4777 | 		         const std::list<BlockVector<double>> & b0xs_vector_list,
4778 | 			 const FullMatrix<double> & M_matrix,
4779 | 			 const BlockVector<double> & gradient_g,
4780 | 			 const BlockVector<double> & solution_delta,
4781 | 			 BlockVector<double> & solution_delta_cauchy_point)
4782 |   {
4783 |     m_timer.enter_subsection("Calculate Cauchy point");
4784 | 
4785 |     solution_delta_cauchy_point = 0.0;
4786 |     BlockVector<double> gradient_d(gradient_g);
4787 |     gradient_d *= -1;
4788 | 
4789 |     const unsigned int list_size = y_vector_list.size();
4790 |     const auto itr_y_begin    = y_vector_list.begin();
4791 |     const auto itr_b0xs_begin = b0xs_vector_list.begin();
4792 | 
4793 |     // t_series only contains t > 0
4794 |     std::priority_queue< std::pair<double, unsigned int>,
4795 |                          std::vector<std::pair<double, unsigned int>>,
4796 |         		 std::greater<std::pair<double, unsigned int>> >
4797 |     t_series = calculate_break_points(solution_delta,
4798 |     			              gradient_g,
4799 | 				      gradient_d);
4800 | 
4801 |     // m_active_set_phasefield contains 1 or 2 for active set and 0 for inactive set
4802 |     for (unsigned int i = 0; i < m_active_set_phasefield.size(); ++i)
4803 |       {
4804 | 	if (m_active_set_phasefield(i) > 0.5)
4805 | 	  solution_delta_cauchy_point.block(m_d_dof)[i] = solution_delta.block(m_d_dof)[i];
4806 |       }
4807 | 
4808 |     // p = W^T * d
4809 |     Vector<double> p(2 * list_size);
4810 |     for (unsigned int i = 0; i < list_size; ++i)
4811 |       {
4812 |         p(i)             = (*std::next(itr_y_begin,    i)) * gradient_d;
4813 |         p(i + list_size) = (*std::next(itr_b0xs_begin, i)) * gradient_d;
4814 |       }
4815 | 
4816 |     Vector<double> c(2 * list_size);
4817 |     c = 0.0;
4818 | 
4819 |     double f_prime = -(gradient_d * gradient_d);
4820 | 
4821 |     // M * p
4822 |     Vector<double> Mp(2 * list_size);
4823 |     if (list_size > 0)
4824 |       M_matrix.vmult(Mp, p);
4825 | 
4826 |     // B_0 * d
4827 |     BlockVector<double> B0_grandient_d(m_dofs_per_block);
4828 |     B0_matrix.vmult(B0_grandient_d, gradient_d);
4829 | 
4830 |     double f_prime_prime = gradient_d * B0_grandient_d;
4831 |     if (list_size > 0)
4832 |       f_prime_prime -= (p * Mp);
4833 | 
4834 |     double delta_t_min = -f_prime / f_prime_prime;
4835 | 
4836 |     double t_old = 0.0;
4837 | 
4838 |     std::pair<double, unsigned int> top_pair = t_series.top();
4839 |     double t = top_pair.first;
4840 |     unsigned int b = top_pair.second;
4841 | 
4842 |     double delta_t = t - t_old;
4843 | 
4844 |     BlockVector<double> z(m_dofs_per_block);
4845 |     z = 0.0;
4846 | 
4847 |     // w_b = W^T * e_b
4848 |     Vector<double> w_b(2 * list_size);
4849 | 
4850 |     // w_b_T_x_M = w_b^T * M
4851 |     Vector<double> w_b_T_x_M(2 * list_size);
4852 | 
4853 |     Vector<double> temp_vector(m_dofs_per_block[m_d_dof]);
4854 | 
4855 |     while (delta_t_min >= delta_t)
4856 |       {
4857 | 	t_series.pop();
4858 | 
4859 | 	if (gradient_d.block(m_d_dof)[b] > 0)
4860 | 	  solution_delta_cauchy_point.block(m_d_dof)[b] = 1.0 - m_solution.block(m_d_dof)[b];
4861 | 	else if (gradient_d.block(m_d_dof)[b] < 0)
4862 | 	  solution_delta_cauchy_point.block(m_d_dof)[b] = 0.0;
4863 | 	else
4864 | 	  AssertThrow(false,
4865 | 	              ExcMessage("gradient_d(b) cannot be zero!"));
4866 | 
4867 | 	if (gradient_d.block(m_d_dof)[b] < 0)
4868 | 	  m_active_set_phasefield[b] = 1; //lower bound
4869 | 	else
4870 | 	  m_active_set_phasefield[b] = 2; //upper bound
4871 | 
4872 |         // z_{j} = z_{j-1} + delta_t_{j-1} * gradient_d_{j-1};
4873 | 	z.sadd(1.0, delta_t, gradient_d);
4874 | 
4875 | 	// c_{j} = c_{j-1} + delta_t_{j-1} * p_{j-1};
4876 | 	if (list_size > 0)
4877 | 	  c.sadd(1.0, delta_t, p);
4878 | 
4879 |         double temp_scalar = ebT_x_B0_x_v(b, B0_matrix, z);
4880 | 
4881 |         // w_b = W^T * e_b
4882 |         for (unsigned int i = 0; i < list_size; ++i)
4883 |           {
4884 |             w_b(i)             = (*std::next(itr_y_begin,    i)).block(m_d_dof)[b];
4885 |             w_b(i + list_size) = (*std::next(itr_b0xs_begin, i)).block(m_d_dof)[b];
4886 |           }
4887 | 
4888 |         if (list_size > 0)
4889 |           M_matrix.vmult(w_b_T_x_M, w_b);
4890 | 
4891 | 	f_prime += delta_t * f_prime_prime
4892 | 	         + gradient_g.block(m_d_dof)[b] * gradient_g.block(m_d_dof)[b]
4893 | 	         + temp_scalar * gradient_g.block(m_d_dof)[b];
4894 | 
4895 | 	if (list_size > 0)
4896 | 	  f_prime -= (w_b_T_x_M * c) * gradient_g.block(m_d_dof)[b];
4897 | 
4898 | 	temp_scalar = ebT_x_B0_x_v(b, B0_matrix, gradient_d);
4899 | 
4900 | 	f_prime_prime += 2.0 * gradient_g.block(m_d_dof)[b] * temp_scalar
4901 | 	               + gradient_g.block(m_d_dof)[b] * gradient_g.block(m_d_dof)[b]
4902 | 		             * B0_matrix.block(m_d_dof, m_d_dof)(b, b);
4903 | 
4904 | 	if (list_size > 0)
4905 | 	  {
4906 | 	    f_prime_prime -= 2.0 * gradient_g.block(m_d_dof)[b] * (w_b_T_x_M * p);
4907 | 	    f_prime_prime -= gradient_g.block(m_d_dof)[b] * gradient_g.block(m_d_dof)[b] * (w_b_T_x_M * w_b);
4908 | 	  }
4909 | 
4910 | 	// p_{j} = p_{j-1} + g_b * w_b;
4911 | 	if (list_size > 0)
4912 | 	  p.sadd(1.0, gradient_g.block(m_d_dof)[b], w_b);
4913 | 
4914 | 	gradient_d.block(m_d_dof)[b] = 0.0;
4915 | 
4916 | 	delta_t_min = -f_prime / f_prime_prime;
4917 | 
4918 | 	t_old = t;
4919 | 
4920 | 	top_pair = t_series.top();
4921 | 	t = top_pair.first;
4922 | 	b = top_pair.second;
4923 | 
4924 | 	delta_t = t - t_old;
4925 |       }
4926 | 
4927 |     if (delta_t_min < 0)
4928 |       delta_t_min = 0;
4929 | 
4930 |     t_old += delta_t_min;
4931 | 
4932 |     for(unsigned int i = 0; i < m_active_set_phasefield.size(); ++i)
4933 |       {
4934 | 	// inactive phasefield dof
4935 | 	if (m_active_set_phasefield(i) < 0.5)
4936 | 	  solution_delta_cauchy_point.block(m_d_dof)[i] = solution_delta.block(m_d_dof)[i]
4937 | 						+ t_old * gradient_d.block(m_d_dof)[i];
4938 |       }
4939 | 
4940 |     // There are no active constraints in the displacement field
4941 |     solution_delta_cauchy_point.block(m_u_dof) = solution_delta.block(m_u_dof);
4942 |     (solution_delta_cauchy_point.block(m_u_dof)).add(t_old, gradient_d.block(m_u_dof));
4943 | 
4944 |     // We need to make sure the solution_delta_cauchy_point satisfies the essential
4945 |     // boundary conditions and the hanging-node constraints
4946 |     m_constraints.distribute(solution_delta_cauchy_point);
4947 | 
4948 |     m_timer.leave_subsection();
4949 |   }
4950 | 
4951 |   template <int dim>
4952 |   void PhaseFieldMonolithicSolve<dim>::
4953 |   solve_nonlinear_timestep_LBFGS_B(BlockVector<double> & solution_delta,
4954 | 				   BlockVector<double> & LBFGS_update_refine)
4955 |   {
4956 |     BlockVector<double> LBFGS_update(m_dofs_per_block);
4957 |     BlockVector<double> solution_delta_cauchy_point(m_dofs_per_block);
4958 |     LBFGS_update = 0.0;
4959 | 
4960 |     const unsigned int LBFGS_m = m_parameters.m_LBFGS_m;
4961 | 
4962 |     unsigned int LBFGS_iteration = 0;
4963 | 
4964 |     m_error_residual.reset();
4965 |     m_error_residual_0.reset();
4966 |     m_error_residual_norm.reset();
4967 |     m_error_update.reset();
4968 |     m_error_update_0.reset();
4969 |     m_error_update_norm.reset();
4970 | 
4971 |     if (m_parameters.m_output_iteration_history)
4972 |       print_conv_header_LBFGSB();
4973 | 
4974 |     BlockVector<double> LBFGS_s_vector(m_dofs_per_block);
4975 |     BlockVector<double> LBFGS_y_vector(m_dofs_per_block);
4976 |     BlockVector<double> free_dofs(m_dofs_per_block);
4977 |     BlockVector<double> b0xs_vector(m_dofs_per_block);
4978 | 
4979 |     // all the list goes from k-m to k-1
4980 |     // the front is the oldest quantity,and the end is
4981 |     // newest quantity
4982 |     std::list<BlockVector<double>> s_vector_list;
4983 |     std::list<BlockVector<double>> y_vector_list;
4984 |     std::list<double> s_dot_y_list;
4985 |     std::list<BlockVector<double>> b0xs_vector_list;
4986 | 
4987 |     double line_search_parameter = 0.0;
4988 | 
4989 |     unsigned int lower_bound_number_old = 0;
4990 |     unsigned int upper_bound_number_old = 0;
4991 |     unsigned int lowerupper_bound_number_old = 0;
4992 | 
4993 |     unsigned int lower_bound_number_new = 0;
4994 |     unsigned int upper_bound_number_new = 0;
4995 |     unsigned int lowerupper_bound_number_new = 0;
4996 | 
4997 |     for (; LBFGS_iteration < m_parameters.m_max_iterations_BFGS; ++LBFGS_iteration)
4998 |       {
4999 | 	if (m_parameters.m_output_iteration_history)
5000 | 	  m_logfile << '\t' << std::setw(4) << LBFGS_iteration << ' '
5001 |                     << std::flush;
5002 | 
5003 |         make_constraints(LBFGS_iteration);
5004 | 
5005 |         // At the first step, we simply distribute the inhomogeneous part of
5006 |         // the constraints
5007 |         if (LBFGS_iteration == 0)
5008 |           {
5009 |             // use the solution from the previous solve on the
5010 |             // refined mesh as initial guess
5011 |             LBFGS_update = LBFGS_update_refine;
5012 | 
5013 |             m_constraints.distribute(LBFGS_update);
5014 |             solution_delta += LBFGS_update;
5015 |             update_qph_incremental(solution_delta, m_solution);
5016 |             assemble_system_rhs_BFGS_parallel(m_solution, m_system_rhs);
5017 |             m_logfile << "  | " << std::endl;
5018 | 
5019 |             continue;
5020 |           }
5021 | 
5022 |         get_error_residual_LBFGSB(m_error_residual,
5023 | 				  solution_delta);
5024 | 
5025 |         if (LBFGS_iteration == 1)
5026 |           m_error_residual_0 = m_error_residual;
5027 | 
5028 |         m_error_residual_norm = m_error_residual;
5029 | 
5030 |         if (m_parameters.m_relative_residual)
5031 |           m_error_residual_norm.normalize(m_error_residual_0);
5032 | 
5033 |         if (LBFGS_iteration > 1 && m_error_update_norm.m_u <= m_parameters.m_tol_u_incr
5034 |                                 && m_error_residual_norm.m_u <= m_parameters.m_tol_u_residual
5035 | 			        && m_error_update_norm.m_d <= m_parameters.m_tol_d_incr
5036 | 			        && m_error_residual_norm.m_d <= m_parameters.m_tol_d_residual
5037 | 				&& lower_bound_number_new == lower_bound_number_old
5038 | 				&& upper_bound_number_new == upper_bound_number_old
5039 | 				&& lowerupper_bound_number_new == lowerupper_bound_number_old)
5040 |           {
5041 |             if (m_parameters.m_output_iteration_history)
5042 |               {
5043 | 		m_logfile << "  | ";
5044 | 		m_logfile << " CONVERGED! " << std::fixed << std::setprecision(3) << std::setw(7)
5045 | 			  << std::scientific
5046 | 		      << "           ---      "
5047 | 		      << "\t\t\t"
5048 | 		      << "  " << m_error_residual_norm.m_norm
5049 | 		      << "  " << m_error_residual_norm.m_u
5050 | 		      << "  " << m_error_residual_norm.m_d
5051 | 		      << "  " << m_error_update_norm.m_norm
5052 | 		      << "  " << m_error_update_norm.m_u
5053 | 		      << "  " << m_error_update_norm.m_d
5054 | 		      << "  " << std::endl;
5055 | 
5056 | 		m_logfile << '\t';
5057 | 		for (unsigned int i = 0; i < 130; ++i)
5058 | 		  m_logfile << '_';
5059 | 		m_logfile << std::endl;
5060 |               }
5061 | 
5062 |             m_logfile << "\t\tThe current L-BFGS-B step converges in "
5063 |         	      << LBFGS_iteration
5064 |         	      << " iterations." << std::endl;
5065 |             m_logfile << "\t\tNumber of active lower bounds not changed." << std::endl;
5066 |     	    m_logfile << "\t\tNumber of active upper bounds not changed." << std::endl;
5067 |     	    m_logfile << "\t\tNumber of active lower-upper bounds not changed." << std::endl;
5068 | 
5069 |             if (m_parameters.m_relative_residual)
5070 |               {
5071 | 		m_logfile << "\t\tProjected gradient of disp. (relative): "
5072 | 			  << m_error_residual_norm.m_u << std::endl;
5073 | 
5074 | 		m_logfile << "\t\tProjected gradient of disp. (absolute): "
5075 | 			  << m_error_residual_norm.m_u * m_error_residual_0.m_u << std::endl;
5076 | 
5077 | 		m_logfile << "\t\tProjected gradient of phasefield (relative): "
5078 | 			  << m_error_residual_norm.m_d << std::endl;
5079 | 
5080 | 		m_logfile << "\t\tProjected gradient of phasefield (absolute): "
5081 | 			  << m_error_residual_norm.m_d * m_error_residual_0.m_d << std::endl;
5082 | 
5083 | 		m_logfile << "\t\tRelative increment of disp.: "
5084 | 			  << m_error_update_norm.m_u << std::endl;
5085 | 
5086 | 		m_logfile << "\t\tAbsolute increment of disp.: "
5087 | 			  << m_error_update_norm.m_u * m_error_update_0.m_u << std::endl;
5088 | 
5089 | 		m_logfile << "\t\tRelative increment of phasefield: "
5090 | 			  << m_error_update_norm.m_d << std::endl;
5091 | 
5092 | 		m_logfile << "\t\tAbsolute increment of phasefield: "
5093 | 			  << m_error_update_norm.m_d * m_error_update_0.m_d << std::endl;
5094 |               }
5095 |             else
5096 |               {
5097 | 		m_logfile << "\t\tProjected gradient of disp. (absolute): "
5098 | 			  << m_error_residual_norm.m_u << std::endl;
5099 | 
5100 | 		m_logfile << "\t\tProjected gradient of phasefield (absolute): "
5101 | 			  << m_error_residual_norm.m_d << std::endl;
5102 | 
5103 | 		m_logfile << "\t\tAbsolute increment of disp.: "
5104 | 			  << m_error_update_norm.m_u << std::endl;
5105 | 
5106 | 		m_logfile << "\t\tAbsolute increment of phasefield: "
5107 | 			  << m_error_update_norm.m_d << std::endl;
5108 |               }
5109 | 
5110 |             break;
5111 |           }
5112 | 
5113 | 	// assemble the initial B_0 matrix at the k-th L-BFGS iteration
5114 | 	// m_solution is the old solution from the previous converged step
5115 | 	// it is needed only for the viscosity term
5116 | 	// the output is m_tangent_matrix (B^0_k)
5117 | 	assemble_system_B0(m_solution);
5118 | 
5119 | 	// B^0_k * s_vector has to be completely recalculated from scratch
5120 | 	// at each L-BFGS iteration, since B^0_k is different
5121 | 	b0xs_vector_list.clear();
5122 | 	for (auto itr = s_vector_list.begin(); itr != s_vector_list.end(); ++itr)
5123 | 	  {
5124 | 	    m_tangent_matrix.vmult(b0xs_vector, *itr);
5125 | 	    b0xs_vector_list.push_back(b0xs_vector);
5126 | 	  }
5127 | 
5128 | 	// In the iteration LBFGS_iteration = 0, only the essential boundary conditions
5129 | 	// are applied.
5130 | 	// WHen LBFGS_iteration = 1, it is the first step of LBFGS update, and the
5131 | 	// s_vector_list and y_vector_list are empty.
5132 | 	// Since the pair of s and y will only be added to the list if s dot y > tol,
5133 | 	// it is safer to decide the matrix dimension by the size of the list.
5134 | 	const unsigned int list_size = s_vector_list.size();
5135 | 	const auto itr_s_begin    = s_vector_list.begin();
5136 | 	const auto itr_y_begin    = y_vector_list.begin();
5137 | 	const auto itr_b0xs_begin = b0xs_vector_list.begin();
5138 | 	const auto itr_s_dot_y_begin = s_dot_y_list.begin();
5139 | 
5140 | 	FullMatrix<double> sTxBxs_matrix(list_size);
5141 | 	sTxBxs_matrix = 0;
5142 | 	for (unsigned int i = 0; i < list_size; ++i)
5143 | 	  for (unsigned int j = 0; j <= i; ++j)
5144 | 	    {
5145 | 	      sTxBxs_matrix(i, j) = (*std::next(itr_s_begin,    i))
5146 | 		                  * (*std::next(itr_b0xs_begin, j));
5147 | 	    }
5148 | 	for (unsigned int i = 0; i < list_size; ++i)
5149 | 	  for (unsigned int j = i + 1; j < list_size; ++j)
5150 | 	    {
5151 | 	      sTxBxs_matrix(i, j) = sTxBxs_matrix(j, i);
5152 | 	    }
5153 | 
5154 | 	FullMatrix<double> D_matrix(list_size);
5155 | 	D_matrix = 0;
5156 | 	for (unsigned int i = 0; i < list_size; ++i)
5157 | 	  D_matrix(i, i) = (*std::next(itr_s_dot_y_begin, i));
5158 | 
5159 | 	FullMatrix<double> L_matrix(list_size);
5160 | 	L_matrix = 0;
5161 | 	for (unsigned int i = 0; i < list_size; ++i)
5162 | 	  for (unsigned int j = 0; j < i; ++j)
5163 | 	    L_matrix(i, j) = (*std::next(itr_s_begin, i))
5164 |                            * (*std::next(itr_y_begin, j));
5165 | 
5166 | 	FullMatrix<double> M_matrix_inv(2 * list_size);
5167 | 	FullMatrix<double> M_matrix(2 * list_size);
5168 | 
5169 | 	M_matrix_inv = 0;
5170 | 	for (unsigned int i = 0; i < list_size; ++i)
5171 | 	  M_matrix_inv(i, i) = -D_matrix(i, i);
5172 | 
5173 | 	for (unsigned int i = 0; i < list_size; ++i)
5174 |           for (unsigned int j = 0; j < list_size; ++j)
5175 |             {
5176 |               M_matrix_inv(i + list_size, j + list_size) = sTxBxs_matrix(i, j);
5177 |               M_matrix_inv(i + list_size, j            ) = L_matrix(i, j);
5178 |               M_matrix_inv(i            , j + list_size) = L_matrix(j, i);
5179 |             }
5180 | 
5181 | 	if (!M_matrix_inv.empty())
5182 | 	  M_matrix.invert(M_matrix_inv);
5183 | 
5184 | 	m_active_set_phasefield = 0;
5185 | 	calculate_cauchy_point(m_tangent_matrix,
5186 | 			       y_vector_list,
5187 | 			       b0xs_vector_list,
5188 | 			       M_matrix,
5189 | 			       m_system_rhs,
5190 | 			       solution_delta,
5191 | 			       solution_delta_cauchy_point);
5192 | 
5193 | 	// We need to find out which DOFs are free:
5194 | 	// no essential boundary conditions, no hanging node constraints
5195 | 	// no active box constraints
5196 | 	unsigned int free_disp_number = 0;
5197 | 	for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
5198 | 	  {
5199 | 	    if (m_constraints.is_constrained(i))
5200 | 	      free_dofs.block(m_u_dof)[i] = -1;
5201 | 	    else
5202 | 	      {
5203 | 	        free_dofs.block(m_u_dof)[i] = 1;
5204 | 	        ++free_disp_number;
5205 | 	      }
5206 | 	  }
5207 | 
5208 | 	unsigned int free_phasefield_number = 0;
5209 | 	for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
5210 | 	  {
5211 | 	    if (   m_constraints.is_constrained(i + m_dofs_per_block[m_u_dof])
5212 | 		|| m_active_set_phasefield(i) > 0.5)
5213 | 	      free_dofs.block(m_d_dof)[i] = -1;
5214 | 	    else
5215 | 	      {
5216 | 	        free_dofs.block(m_d_dof)[i] = 1;
5217 | 	        ++free_phasefield_number;
5218 | 	      }
5219 | 	  }
5220 | 
5221 | 	// temp_vector_1 = x^c - x_k
5222 | 	BlockVector<double> temp_vector_1(solution_delta_cauchy_point);
5223 | 	temp_vector_1 -= solution_delta;
5224 | 
5225 | 	// temp_vector_2 = B_0 * (x^c - x_k)
5226 | 	BlockVector<double> temp_vector_2(m_dofs_per_block);
5227 | 	m_tangent_matrix.vmult(temp_vector_2, temp_vector_1);
5228 | 
5229 | 	// temp_vector_3 = W^T * (x^c - x_k)
5230 | 	Vector<double> temp_vector_3(2 * list_size);
5231 | 	for (unsigned int i = 0; i < list_size; ++i)
5232 | 	  {
5233 | 	    temp_vector_3(i)             = (*std::next(itr_y_begin,    i)) * temp_vector_1;
5234 | 	    temp_vector_3(i + list_size) = (*std::next(itr_b0xs_begin, i)) * temp_vector_1;
5235 | 	  }
5236 | 
5237 | 	// temp_vector_4 = M * W^T * (x^c - x_k)
5238 | 	Vector<double> temp_vector_4(2 * list_size);
5239 | 	if (list_size > 0)
5240 | 	  M_matrix.vmult(temp_vector_4, temp_vector_3);
5241 | 
5242 | 	// temp_vector_5 = W * M * W^T * (x^c - x_k)
5243 | 	BlockVector<double> temp_vector_5(m_dofs_per_block);
5244 | 	for (unsigned int i = 0; i < list_size; ++i)
5245 | 	  {
5246 | 	    temp_vector_5.add(temp_vector_4(i),             (*std::next(itr_y_begin,    i)));
5247 | 	    temp_vector_5.add(temp_vector_4(i + list_size), (*std::next(itr_b0xs_begin, i)));
5248 | 	  }
5249 | 
5250 | 	// temp_vector_2 = B_0 * (x^c - x_k) - W * M * W^T * (x^c - x_k)
5251 | 	if (list_size > 0)
5252 | 	  temp_vector_2 -= temp_vector_5;
5253 | 
5254 | 	// temp_vector_2 = g + B_0 * (x^c - x_k) - W * M * W^T * (x^c - x_k)
5255 | 	temp_vector_2 += m_system_rhs;
5256 | 
5257 | 	// temp_vector_2 = Z^T * [g + B_0 * (x^c - x_k) - W * M * W^T * (x^c - x_k)]
5258 | 	for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
5259 | 	  {
5260 | 	    if (free_dofs.block(m_u_dof)[i] < 0)
5261 | 	      temp_vector_2.block(m_u_dof)[i] = 0;
5262 | 	  }
5263 | 
5264 | 	for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
5265 | 	  {
5266 | 	    if (free_dofs.block(m_d_dof)[i] < 0)
5267 | 	      temp_vector_2.block(m_d_dof)[i] = 0;
5268 | 	  }
5269 | 
5270 | 	BlockVector<double> rhs_vector(temp_vector_2);
5271 | 	rhs_vector *= -1;
5272 | 
5273 | 	BlockVector<double> search_direction(m_dofs_per_block);
5274 | 	search_direction = 0.0;
5275 | 
5276 | 	unsigned int cg_iterations = 0;
5277 | 
5278 | 	double alpha_backtrack = 1.0;
5279 | 
5280 | 	if (m_parameters.m_type_linear_solver == "CG")
5281 | 	  {
5282 | 	    m_timer.enter_subsection("Subspace CG solve (inverse operator)");
5283 | 
5284 | 	    //const double rc_hat_norm = rhs_vector.l2_norm();
5285 | 	    const double cg_tol = m_parameters.m_CG_tolerace; //std::min( 0.1, std::sqrt(rc_hat_norm) ) * rc_hat_norm;
5286 | 
5287 | 	    zT_B0_z(free_dofs, m_tangent_matrix);
5288 | 
5289 | 	    const auto op_zT_B0_z = block_operator(m_tangent_matrix);
5290 | 
5291 | 	    if (list_size > 0)
5292 | 	      {
5293 | 		std::list<BlockVector<double>> zT_y_list;
5294 | 		BlockVector<double> zT_y_vector(m_dofs_per_block);
5295 | 		for (auto itr = y_vector_list.begin(); itr != y_vector_list.end(); ++itr)
5296 | 		  {
5297 | 		    zT_y_vector = (*itr);
5298 | 		    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
5299 | 		      {
5300 | 			if (free_dofs.block(m_u_dof)[i] < 0)
5301 | 			  zT_y_vector.block(m_u_dof)[i] = 0;
5302 | 		      }
5303 | 
5304 | 		    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
5305 | 		      {
5306 | 			if (free_dofs.block(m_d_dof)[i] < 0)
5307 | 			  zT_y_vector.block(m_d_dof)[i] = 0;
5308 | 		      }
5309 | 
5310 | 		    zT_y_list.push_back(zT_y_vector);
5311 | 		  }
5312 | 
5313 | 		std::list<BlockVector<double>> zT_b0xs_list;
5314 | 		BlockVector<double> zT_b0xs_vector(m_dofs_per_block);
5315 | 		for (auto itr = b0xs_vector_list.begin(); itr != b0xs_vector_list.end(); ++itr)
5316 | 		  {
5317 | 		    zT_b0xs_vector = (*itr);
5318 | 		    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
5319 | 		      {
5320 | 			if (free_dofs.block(m_u_dof)[i] < 0)
5321 | 			  zT_b0xs_vector.block(m_u_dof)[i] = 0;
5322 | 		      }
5323 | 
5324 | 		    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
5325 | 		      {
5326 | 			if (free_dofs.block(m_d_dof)[i] < 0)
5327 | 			  zT_b0xs_vector.block(m_d_dof)[i] = 0;
5328 | 		      }
5329 | 
5330 | 		    zT_b0xs_list.push_back(zT_b0xs_vector);
5331 | 		  }
5332 | 
5333 | 		const auto op_M_matrix = linear_operator(M_matrix);
5334 | 
5335 | 		FullMatrix<double> zT_W_matrix_u(m_dofs_per_block[m_u_dof], 2*list_size);
5336 | 		unsigned int j = 0;
5337 | 		for (auto itr = zT_y_list.begin(); itr != zT_y_list.end(); ++itr)
5338 | 		  {
5339 | 		    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
5340 | 		      zT_W_matrix_u(i, j) = (*itr).block(m_u_dof)[i];
5341 | 		    ++j;
5342 | 		  }
5343 | 		j = 0;
5344 | 		for (auto itr = zT_b0xs_list.begin(); itr != zT_b0xs_list.end(); ++itr)
5345 | 		  {
5346 | 		    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
5347 | 		      zT_W_matrix_u(i, j + list_size) = (*itr).block(m_u_dof)[i];
5348 | 		    ++j;
5349 | 		  }
5350 | 
5351 | 		FullMatrix<double> zT_W_matrix_d(m_dofs_per_block[m_d_dof], 2*list_size);
5352 | 		j = 0;
5353 | 		for (auto itr = zT_y_list.begin(); itr != zT_y_list.end(); ++itr)
5354 | 		  {
5355 | 		    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
5356 | 		      zT_W_matrix_d(i, j) = (*itr).block(m_d_dof)[i];
5357 | 		    ++j;
5358 | 		  }
5359 | 		j = 0;
5360 | 		for (auto itr = zT_b0xs_list.begin(); itr != zT_b0xs_list.end(); ++itr)
5361 | 		  {
5362 | 		    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
5363 | 		      zT_W_matrix_d(i, j + list_size) = (*itr).block(m_d_dof)[i];
5364 | 		    ++j;
5365 | 		  }
5366 | 
5367 | 		const auto op_zT_W_matrix_u = linear_operator(zT_W_matrix_u);
5368 | 		const auto op_zT_W_matrix_d = linear_operator(zT_W_matrix_d);
5369 | 
5370 | 		const auto op_uMuT = op_zT_W_matrix_u * op_M_matrix * transpose_operator(op_zT_W_matrix_u);
5371 | 
5372 | 		const auto op_uMdT = op_zT_W_matrix_u * op_M_matrix * transpose_operator(op_zT_W_matrix_d);
5373 | 
5374 | 		const auto op_dMuT = op_zT_W_matrix_d * op_M_matrix * transpose_operator(op_zT_W_matrix_u);
5375 | 
5376 | 		const auto op_dMdT = op_zT_W_matrix_d * op_M_matrix * transpose_operator(op_zT_W_matrix_d);
5377 | 
5378 | 		const auto op_zT_wMwT_z = block_operator<2, 2, BlockVector<double>>({op_uMuT, op_uMdT,
5379 | 										     op_dMuT, op_dMdT});
5380 | 
5381 | 		const auto op_total = op_zT_B0_z - op_zT_wMwT_z;
5382 | 
5383 | 		SolverControl            solver_control(1e5, cg_tol);
5384 | 		SolverCG<BlockVector<double>> cg(solver_control);
5385 | 
5386 | 		if (m_parameters.m_type_preconditioner == "None")
5387 | 		  {
5388 | 		    // somehow op_total_inv has to be made const, or the
5389 | 		    // program will have compliation error
5390 | 		    const auto op_total_inv = inverse_operator(op_total, cg);
5391 | 		    op_total_inv.vmult(search_direction, rhs_vector);
5392 | 		  }
5393 | 		else if (m_parameters.m_type_preconditioner == "Jacobi")
5394 | 		  {
5395 | 		    usr_Jacobi_preconditioner preconditioner(m_tangent_matrix);
5396 | 
5397 | 		    const auto op_total_inv = inverse_operator(op_total, cg, preconditioner);
5398 | 		    op_total_inv.vmult(search_direction, rhs_vector);
5399 | 		  }
5400 | 		else if (m_parameters.m_type_preconditioner == "LU")
5401 | 		  {
5402 | 		    SparseDirectUMFPACK matrix_factorization;
5403 | 		    matrix_factorization.initialize(m_tangent_matrix);
5404 | 
5405 | 		    usr_sparseLU_preconditioner preconditioner(matrix_factorization);
5406 | 		    const auto op_total_inv = inverse_operator(op_total, cg, preconditioner);
5407 | 		    op_total_inv.vmult(search_direction, rhs_vector);
5408 | 		  }
5409 | 		else if (m_parameters.m_type_preconditioner == "ILU")
5410 | 		  {
5411 | 		    SparseILU<double> SparseILU_disp;
5412 | 		    SparseILU_disp.initialize(m_tangent_matrix.block(0, 0));
5413 | 		    SparseILU<double> SparseILU_phasefield;
5414 | 		    SparseILU_phasefield.initialize(m_tangent_matrix.block(1, 1));
5415 | 
5416 | 		    usr_sparseILU_preconditioner preconditioner(SparseILU_disp,
5417 | 								SparseILU_phasefield);
5418 | 
5419 | 		    const auto op_total_inv = inverse_operator(op_total, cg, preconditioner);
5420 | 		    op_total_inv.vmult(search_direction, rhs_vector);
5421 | 		  }
5422 | 		else
5423 | 		  {
5424 | 		    AssertThrow(false, ExcMessage("Preconditioner type not implemented"));
5425 | 		  }
5426 | 
5427 | 		cg_iterations = solver_control.last_step();
5428 | 	      } // if (list_size > 0)
5429 | 	    else
5430 | 	      {
5431 | 		const auto op_total = op_zT_B0_z;
5432 | 		SolverControl            solver_control(1e5, cg_tol);
5433 | 		SolverCG<BlockVector<double>> cg(solver_control);
5434 | 
5435 | 		if (m_parameters.m_type_preconditioner == "None")
5436 | 		  {
5437 | 		    // somehow op_total_inv has to be made const, or the
5438 | 		    // program will have compliation error
5439 | 		    const auto op_total_inv = inverse_operator(op_total, cg);
5440 | 		    op_total_inv.vmult(search_direction, rhs_vector);
5441 | 		  }
5442 | 		else if (m_parameters.m_type_preconditioner == "Jacobi")
5443 | 		  {
5444 | 		    usr_Jacobi_preconditioner preconditioner(m_tangent_matrix);
5445 | 
5446 | 		    const auto op_total_inv = inverse_operator(op_total, cg, preconditioner);
5447 | 		    op_total_inv.vmult(search_direction, rhs_vector);
5448 | 		  }
5449 | 		else if (m_parameters.m_type_preconditioner == "LU")
5450 | 		  {
5451 | 		    SparseDirectUMFPACK matrix_factorization;
5452 | 		    matrix_factorization.initialize(m_tangent_matrix);
5453 | 
5454 | 		    usr_sparseLU_preconditioner preconditioner(matrix_factorization);
5455 | 		    const auto op_total_inv = inverse_operator(op_total, cg, preconditioner);
5456 | 		    op_total_inv.vmult(search_direction, rhs_vector);
5457 | 		  }
5458 | 		else if (m_parameters.m_type_preconditioner == "ILU")
5459 | 		  {
5460 | 		    SparseILU<double> SparseILU_disp;
5461 | 		    SparseILU_disp.initialize(m_tangent_matrix.block(0, 0));
5462 | 		    SparseILU<double> SparseILU_phasefield;
5463 | 		    SparseILU_phasefield.initialize(m_tangent_matrix.block(1, 1));
5464 | 
5465 | 		    usr_sparseILU_preconditioner preconditioner(SparseILU_disp,
5466 | 								SparseILU_phasefield);
5467 | 
5468 | 		    const auto op_total_inv = inverse_operator(op_total, cg, preconditioner);
5469 | 		    op_total_inv.vmult(search_direction, rhs_vector);
5470 | 		  }
5471 | 		else
5472 | 		  {
5473 | 		    AssertThrow(false, ExcMessage("Preconditioner type not implemented"));
5474 | 		  }
5475 | 
5476 | 		cg_iterations = solver_control.last_step();
5477 | 	      } // // if (list_size == 0)
5478 | 
5479 |             m_timer.leave_subsection();
5480 | 	  } // if (m_parameters.m_type_linear_solver == "CG")
5481 | 	else if (m_parameters.m_type_linear_solver == "Direct")
5482 | 	  {
5483 | 	    m_timer.enter_subsection("Subspace direct solve (LU factorization)");
5484 | 
5485 | 	    zT_B0_z(free_dofs, m_tangent_matrix);
5486 | 
5487 | 	    SparseDirectUMFPACK zT_B0_z_inv;
5488 | 	    zT_B0_z_inv.initialize(m_tangent_matrix);
5489 | 
5490 | 	    //SparseDirectUMFPACK zT_B0_z_inv_disp;
5491 | 	    //zT_B0_z_inv_disp.initialize(m_tangent_matrix.block(m_u_dof, m_u_dof));
5492 | 
5493 | 	    //SparseDirectUMFPACK zT_B0_z_inv_phasefield;
5494 | 	    //zT_B0_z_inv_phasefield.initialize(m_tangent_matrix.block(m_d_dof, m_d_dof));
5495 | 
5496 | 	    m_timer.leave_subsection();
5497 | 
5498 | 	    m_timer.enter_subsection("Subspace direct solve (LU matrix-vector multiplication)");
5499 | 
5500 | 	    zT_B0_z_inv.vmult(search_direction, rhs_vector);
5501 | 	    //zT_B0_z_inv_disp.vmult(search_direction.block(m_u_dof), rhs_vector.block(m_u_dof));
5502 | 	    //zT_B0_z_inv_phasefield.vmult(search_direction.block(m_d_dof), rhs_vector.block(m_d_dof));
5503 | 
5504 | 	    BlockVector<double> update_vector(m_dofs_per_block);
5505 | 	    update_vector = 0;
5506 | 	    if (list_size > 0)
5507 | 	      {
5508 | 		std::list<BlockVector<double>> zT_B0_z_inv_zT_y_list;
5509 | 		std::list<BlockVector<double>> zT_y_list;
5510 | 		BlockVector<double> zT_y_vector(m_dofs_per_block);
5511 | 		BlockVector<double> zT_B0_z_inv_zT_y_vector(m_dofs_per_block);
5512 | 		for (auto itr = y_vector_list.begin(); itr != y_vector_list.end(); ++itr)
5513 | 		  {
5514 | 		    zT_y_vector = (*itr);
5515 | 		    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
5516 | 		      {
5517 | 			if (free_dofs.block(m_u_dof)[i] < 0)
5518 | 			  zT_y_vector.block(m_u_dof)[i] = 0;
5519 | 		      }
5520 | 
5521 | 		    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
5522 | 		      {
5523 | 			if (free_dofs.block(m_d_dof)[i] < 0)
5524 | 			  zT_y_vector.block(m_d_dof)[i] = 0;
5525 | 		      }
5526 | 
5527 | 		    zT_y_list.push_back(zT_y_vector);
5528 | 
5529 | 		    //zT_B0_z_inv_disp.vmult(zT_B0_z_inv_zT_y_vector.block(m_u_dof), zT_y_vector.block(m_u_dof));
5530 | 		    //zT_B0_z_inv_phasefield.vmult(zT_B0_z_inv_zT_y_vector.block(m_d_dof), zT_y_vector.block(m_d_dof));
5531 | 		    zT_B0_z_inv.vmult(zT_B0_z_inv_zT_y_vector, zT_y_vector);
5532 | 
5533 | 		    zT_B0_z_inv_zT_y_list.push_back(zT_B0_z_inv_zT_y_vector);
5534 | 		  }
5535 | 
5536 | 		std::list<BlockVector<double>> zT_B0_z_inv_zT_b0xs_list;
5537 | 		std::list<BlockVector<double>> zT_b0xs_list;
5538 | 		BlockVector<double> zT_b0xs_vector(m_dofs_per_block);
5539 | 		BlockVector<double> zT_B0_z_inv_zT_b0xs_vector(m_dofs_per_block);
5540 | 		for (auto itr = b0xs_vector_list.begin(); itr != b0xs_vector_list.end(); ++itr)
5541 | 		  {
5542 | 		    zT_b0xs_vector = (*itr);
5543 | 		    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
5544 | 		      {
5545 | 			if (free_dofs.block(m_u_dof)[i] < 0)
5546 | 			  zT_b0xs_vector.block(m_u_dof)[i] = 0;
5547 | 		      }
5548 | 
5549 | 		    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
5550 | 		      {
5551 | 			if (free_dofs.block(m_d_dof)[i] < 0)
5552 | 			  zT_b0xs_vector.block(m_d_dof)[i] = 0;
5553 | 		      }
5554 | 
5555 | 		    zT_b0xs_list.push_back(zT_b0xs_vector);
5556 | 
5557 | 		    //zT_B0_z_inv_disp.vmult(zT_B0_z_inv_zT_b0xs_vector.block(m_u_dof), zT_b0xs_vector.block(m_u_dof));
5558 | 		    //zT_B0_z_inv_phasefield.vmult(zT_B0_z_inv_zT_b0xs_vector.block(m_d_dof), zT_b0xs_vector.block(m_d_dof));
5559 | 		    zT_B0_z_inv.vmult(zT_B0_z_inv_zT_b0xs_vector, zT_b0xs_vector);
5560 | 
5561 | 		    zT_B0_z_inv_zT_b0xs_list.push_back(zT_B0_z_inv_zT_b0xs_vector);
5562 | 		  }
5563 | 
5564 | 		FullMatrix<double> wT_z_zT_B0_z_inv_zT_w(2 * list_size);
5565 | 		const auto itr_zT_y_list_begin = zT_y_list.begin();
5566 | 		const auto itr_zT_b0xs_list_begin = zT_b0xs_list.begin();
5567 | 		const auto itr_zT_B0_z_inv_zT_y_list_begin = zT_B0_z_inv_zT_y_list.begin();
5568 | 		const auto itr_zT_B0_z_inv_zT_b0xs_list_begin = zT_B0_z_inv_zT_b0xs_list.begin();
5569 | 		for (unsigned int i = 0; i < list_size; ++i)
5570 | 		  for (unsigned int j = 0; j < list_size; ++j)
5571 | 		    {
5572 | 		      wT_z_zT_B0_z_inv_zT_w(i          , j          ) = (*std::next(itr_zT_y_list_begin            , i))
5573 | 								      * (*std::next(itr_zT_B0_z_inv_zT_y_list_begin, j));
5574 | 
5575 | 		      wT_z_zT_B0_z_inv_zT_w(i          , j+list_size) = (*std::next(itr_zT_y_list_begin               , i))
5576 | 								      * (*std::next(itr_zT_B0_z_inv_zT_b0xs_list_begin, j));
5577 | 
5578 | 		      wT_z_zT_B0_z_inv_zT_w(i+list_size, j          ) = (*std::next(itr_zT_b0xs_list_begin         , i))
5579 | 								      * (*std::next(itr_zT_B0_z_inv_zT_y_list_begin, j));
5580 | 
5581 | 		      wT_z_zT_B0_z_inv_zT_w(i+list_size, j+list_size) = (*std::next(itr_zT_b0xs_list_begin            , i))
5582 | 								      * (*std::next(itr_zT_B0_z_inv_zT_b0xs_list_begin, j));
5583 | 		    }
5584 | 
5585 | 		FullMatrix<double> temp_matrix(2 * list_size);
5586 | 		M_matrix.mmult(temp_matrix, wT_z_zT_B0_z_inv_zT_w);
5587 | 
5588 | 		FullMatrix<double> middle_matrix(IdentityMatrix(2*list_size));
5589 | 		middle_matrix.add(-1.0, temp_matrix);
5590 | 
5591 | 		FullMatrix<double> middle_matrix_inv(2 * list_size);
5592 | 		middle_matrix_inv.invert(middle_matrix);
5593 | 
5594 | 		middle_matrix_inv.mmult(middle_matrix, M_matrix);
5595 | 
5596 | 		Vector<double> wT_z_zT_B0_z_inv_rhs(2 * list_size);
5597 | 		for (unsigned int i = 0; i < list_size; ++i)
5598 | 		  {
5599 | 		    wT_z_zT_B0_z_inv_rhs(i            ) = (*std::next(itr_zT_B0_z_inv_zT_y_list_begin   , i)) * rhs_vector;
5600 | 		    wT_z_zT_B0_z_inv_rhs(i + list_size) = (*std::next(itr_zT_B0_z_inv_zT_b0xs_list_begin, i)) * rhs_vector;
5601 | 		  }
5602 | 
5603 | 		Vector<double> middle_matrix_wT_z_zT_B0_z_inv_rhs(2 * list_size);
5604 | 		middle_matrix.vmult(middle_matrix_wT_z_zT_B0_z_inv_rhs,
5605 | 				    wT_z_zT_B0_z_inv_rhs);
5606 | 
5607 | 		unsigned int index = 0;
5608 | 		for (auto itr = zT_B0_z_inv_zT_y_list.begin(); itr != zT_B0_z_inv_zT_y_list.end(); ++itr)
5609 | 		  {
5610 | 		    update_vector.add(middle_matrix_wT_z_zT_B0_z_inv_rhs(index), *itr);
5611 | 		    ++index;
5612 | 		  }
5613 | 		for (auto itr = zT_B0_z_inv_zT_b0xs_list.begin(); itr != zT_B0_z_inv_zT_b0xs_list.end(); ++itr)
5614 | 		  {
5615 | 		    update_vector.add(middle_matrix_wT_z_zT_B0_z_inv_rhs(index), *itr);
5616 | 		    ++index;
5617 | 		  }
5618 | 	      } //	if (list_size > 0)
5619 | 
5620 | 	    search_direction += update_vector;
5621 | 
5622 | 	    m_timer.leave_subsection();
5623 | 	  } // else if (m_parameters.m_type_linear_solver == "Direct")
5624 | 	else
5625 | 	  {
5626 | 	    AssertThrow(false, ExcMessage("Linear solver type not implemented"));
5627 | 	  }
5628 | 
5629 | 	// We don't do backtrack yet. We will make sure phasefield
5630 | 	// remains feasible later
5631 | 	alpha_backtrack = 1.0;
5632 | 	search_direction *= alpha_backtrack;
5633 | 
5634 | 	LBFGS_update.block(m_u_dof) = solution_delta_cauchy_point.block(m_u_dof);
5635 | 	LBFGS_update.block(m_u_dof) += search_direction.block(m_u_dof);
5636 | 	LBFGS_update.block(m_u_dof) -= solution_delta.block(m_u_dof);
5637 | 
5638 | 	for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
5639 | 	  {
5640 | 	    // phasefield active constraints
5641 | 	    if (m_active_set_phasefield(i) > 0.5)
5642 | 	      LBFGS_update.block(m_d_dof)[i] = solution_delta_cauchy_point.block(m_d_dof)[i]
5643 | 					     - solution_delta.block(m_d_dof)[i];
5644 | 	    else
5645 | 	      LBFGS_update.block(m_d_dof)[i] = solution_delta_cauchy_point.block(m_d_dof)[i]
5646 | 					     + search_direction.block(m_d_dof)[i]
5647 | 					     - solution_delta.block(m_d_dof)[i];
5648 | 	  }
5649 | 
5650 | 	// make sure the phasefield solutions are feasible
5651 | 	for(unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
5652 | 	  {
5653 | 	    if (solution_delta.block(m_d_dof)[i] + LBFGS_update.block(m_d_dof)[i] < 0.0)
5654 | 	      LBFGS_update.block(m_d_dof)[i] = -solution_delta.block(m_d_dof)[i];
5655 | 
5656 | 	    if (  solution_delta.block(m_d_dof)[i]
5657 | 		+ m_solution.block(m_d_dof)[i]
5658 | 		+ LBFGS_update.block(m_d_dof)[i] > 1.0)
5659 | 	      LBFGS_update.block(m_d_dof)[i] = 1.0 - m_solution.block(m_d_dof)[i]
5660 | 						   - solution_delta.block(m_d_dof)[i];
5661 | 	  }
5662 | 
5663 | 	m_constraints.distribute(LBFGS_update);
5664 | 
5665 | 	// We need a line search algorithm to decide line_search_parameter
5666 | 
5667 |         if(m_parameters.m_type_line_search == "StrongWolfe")
5668 |           {
5669 | 	    const double phi_0 = calculate_energy_functional();
5670 | 	    const double phi_0_prime = m_system_rhs * LBFGS_update;
5671 | 
5672 | 	    line_search_parameter = line_search_stepsize_strong_wolfe(phi_0,
5673 | 								      phi_0_prime,
5674 | 								      LBFGS_update,
5675 | 								      solution_delta);
5676 |           }
5677 |         else if(m_parameters.m_type_line_search == "GradientBased")
5678 |           {
5679 | 	    // LBFGS_r_vector is the search direction
5680 | 	    line_search_parameter = line_search_stepsize_gradient_based(LBFGS_update,
5681 | 									solution_delta);
5682 |           }
5683 |         else
5684 |           {
5685 |             Assert(false, ExcMessage("An unknown line search method is called!"));
5686 |           }
5687 | 
5688 | 	LBFGS_update *= line_search_parameter;
5689 | 
5690 |         get_error_update(LBFGS_update, m_error_update);
5691 |         if (LBFGS_iteration == 1)
5692 |           m_error_update_0 = m_error_update;
5693 | 
5694 |         m_error_update_norm = m_error_update;
5695 |         // For three-point bending problem and the sphere inclusion problem,
5696 |         // we use absolute residual for convergence test
5697 |         if (m_parameters.m_relative_residual)
5698 |           m_error_update_norm.normalize(m_error_update_0);
5699 | 
5700 | 	solution_delta += LBFGS_update;
5701 | 
5702 | 	update_qph_incremental(solution_delta, m_solution);
5703 | 
5704 |         LBFGS_y_vector = m_system_rhs;
5705 |         LBFGS_y_vector *= -1.0;
5706 |         assemble_system_rhs_BFGS_parallel(m_solution, m_system_rhs);
5707 |         // if we use assemble_system_rhs_BFGS_parallel, then condense() is not necessary
5708 |         //m_constraints.condense(m_system_rhs);
5709 |         LBFGS_y_vector += m_system_rhs;
5710 | 
5711 |         LBFGS_s_vector = LBFGS_update;
5712 | 
5713 | 	// s_vector_list, y_vector_list, s_dot_y_list only need to discard
5714 | 	// the front (oldest) item and add the newest item to the end at
5715 | 	// each L-BFGS iteration
5716 | 	double s_dot_y = LBFGS_s_vector * LBFGS_y_vector;
5717 | 	if (s_dot_y > 1.0e-16 * LBFGS_y_vector.norm_sqr())
5718 | 	  {
5719 | 	    if (list_size >= LBFGS_m)
5720 | 	      {
5721 | 		s_vector_list.pop_front();
5722 | 		y_vector_list.pop_front();
5723 | 		s_dot_y_list.pop_front();
5724 | 	      }
5725 | 
5726 | 	    s_vector_list.push_back(LBFGS_s_vector);
5727 | 	    y_vector_list.push_back(LBFGS_y_vector);
5728 | 	    s_dot_y_list.push_back(s_dot_y);
5729 | 	  }
5730 | 
5731 | 	Vector<double> solution_phasefield_total(m_solution.block(m_d_dof));
5732 | 	solution_phasefield_total += solution_delta.block(m_d_dof);
5733 | 
5734 | 	// Since line search parameter might be less than one, we need update
5735 | 	// the phasefield active set status
5736 | 	// upper bound is 1.0, lower bound is the solution at the previous step.
5737 | 	unsigned int number_active_constraint_lower_bound = 0;
5738 | 	unsigned int number_active_constraint_upper_bound = 0;
5739 | 	unsigned int number_active_constraint_lowerupper_bound = 0;
5740 | 
5741 | 	for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
5742 | 	  {
5743 | 	    if (   solution_delta.block(m_d_dof)[i] == 0.0
5744 | 		&& solution_phasefield_total[i] == 1.0)
5745 | 	      {
5746 | 		m_active_set_phasefield(i) = 3; //lower bound overlaps with upper bound
5747 | 		++number_active_constraint_lowerupper_bound;
5748 | 	      }
5749 | 	    else if (   solution_delta.block(m_d_dof)[i] == 0.0
5750 | 		     && solution_phasefield_total[i] != 1.0)
5751 | 	      {
5752 | 		m_active_set_phasefield(i) = 1; //lower bound
5753 | 		++number_active_constraint_lower_bound;
5754 | 	      }
5755 | 	    else if (   solution_phasefield_total[i] == 1.0
5756 | 		     && solution_delta.block(m_d_dof)[i] != 0.0)
5757 | 	      {
5758 | 	        m_active_set_phasefield(i) = 2; //upper bound
5759 | 	        ++number_active_constraint_upper_bound;
5760 | 	      }
5761 | 	    else
5762 | 	      {
5763 | 	        m_active_set_phasefield(i) = 0;
5764 | 	      }
5765 | 	  }
5766 | 
5767 | 	lower_bound_number_old = lower_bound_number_new;
5768 | 	upper_bound_number_old = upper_bound_number_new;
5769 | 	lowerupper_bound_number_old = lowerupper_bound_number_new;
5770 | 
5771 | 	lower_bound_number_new = number_active_constraint_lower_bound;
5772 | 	upper_bound_number_new = number_active_constraint_upper_bound;
5773 | 	lowerupper_bound_number_new = number_active_constraint_lowerupper_bound;
5774 | 
5775 | 	if (m_parameters.m_output_iteration_history)
5776 |           {
5777 | 	    const double energy_functional = calculate_energy_functional();
5778 | 
5779 | 	    m_logfile << "  | "
5780 | 		      << std::setw(6) << number_active_constraint_lower_bound
5781 | 		      << std::setw(6) << number_active_constraint_upper_bound
5782 | 		      << std::setw(6) << number_active_constraint_lowerupper_bound;
5783 |             if (m_parameters.m_type_linear_solver == "CG")
5784 |               m_logfile << std::setw(8) << cg_iterations;
5785 |             else
5786 |               m_logfile << std::setw(8) << "---";
5787 |             m_logfile << "      "
5788 |         	      << std::fixed << std::setprecision(3) << std::scientific << line_search_parameter
5789 | 		      << std::fixed << std::setprecision(6) << std::setw(1) << std::scientific
5790 | 		      << "  " << energy_functional
5791 | 		      << std::fixed << std::setprecision(3) << std::setw(1)
5792 | 					<< std::scientific
5793 | 		      << "  " << m_error_residual_norm.m_norm
5794 | 		      << "  " << m_error_residual_norm.m_u
5795 | 		      << "  " << m_error_residual_norm.m_d
5796 | 		      << "  " << m_error_update_norm.m_norm
5797 | 		      << "  " << m_error_update_norm.m_u
5798 | 		      << "  " << m_error_update_norm.m_d
5799 | 		      << "  " << std::endl;
5800 |           }
5801 |       } // for (; LBFGS_iteration < m_parameters.m_max_iterations_BFGS; ++LBFGS_iteration)
5802 | 
5803 |     AssertThrow(LBFGS_iteration < m_parameters.m_max_iterations_BFGS,
5804 |                 ExcMessage("No convergence in L-BFGS-B nonlinear solver!"));
5805 |   }
5806 | 
5807 |   template <int dim>
5808 |   void PhaseFieldMonolithicSolve<dim>::output_results() const
5809 |   {
5810 |     m_timer.enter_subsection("Output results");
5811 | 
5812 |     DataOut<dim> data_out;
5813 | 
5814 |     std::vector<DataComponentInterpretation::DataComponentInterpretation>
5815 |       data_component_interpretation(
5816 |         dim, DataComponentInterpretation::component_is_part_of_vector);
5817 | 
5818 |     data_component_interpretation.push_back(
5819 |       DataComponentInterpretation::component_is_scalar);
5820 | 
5821 |     std::vector<std::string> solution_name(dim, "displacement");
5822 |     solution_name.emplace_back("phasefield");
5823 | 
5824 |     data_out.attach_dof_handler(m_dof_handler);
5825 |     data_out.add_data_vector(m_solution,
5826 |                              solution_name,
5827 |                              DataOut<dim>::type_dof_data,
5828 |                              data_component_interpretation);
5829 | 
5830 |     // output phasefield active set status
5831 |     BlockVector<double> active_set_status(m_dofs_per_block);
5832 |     active_set_status.block(m_d_dof) = m_active_set_phasefield;
5833 |     std::vector<DataComponentInterpretation::DataComponentInterpretation>
5834 |       data_component_interpretation_active_set(
5835 |         dim+1, DataComponentInterpretation::component_is_scalar);
5836 |     std::vector<std::string> solution_name_active_set;
5837 |     solution_name_active_set.emplace_back("disp_x_active_set");
5838 |     solution_name_active_set.emplace_back("disp_y_active_set");
5839 |     if (dim ==3)
5840 |       solution_name_active_set.emplace_back("disp_z_active_set");
5841 |     solution_name_active_set.emplace_back("phasefield_active_set");
5842 |     data_out.add_data_vector(active_set_status,
5843 | 			     solution_name_active_set,
5844 | 			     DataOut<dim>::type_dof_data,
5845 | 			     data_component_interpretation_active_set);
5846 | 
5847 |     Vector<double> cell_material_id(m_triangulation.n_active_cells());
5848 |     // output material ID for each cell
5849 |     for (const auto &cell : m_triangulation.active_cell_iterators())
5850 |       {
5851 | 	cell_material_id(cell->active_cell_index()) = cell->material_id();
5852 |       }
5853 |     data_out.add_data_vector(cell_material_id, "materialID");
5854 | 
5855 |     // Stress L2 projection
5856 |     DoFHandler<dim> stresses_dof_handler_L2(m_triangulation);
5857 |     FE_Q<dim>     stresses_fe_L2(m_parameters.m_poly_degree); //FE_Q element is continuous
5858 |     stresses_dof_handler_L2.distribute_dofs(stresses_fe_L2);
5859 |     AffineConstraints<double> constraints;
5860 |     constraints.clear();
5861 |     DoFTools::make_hanging_node_constraints(stresses_dof_handler_L2, constraints);
5862 |     constraints.close();
5863 |     std::vector<DataComponentInterpretation::DataComponentInterpretation>
5864 | 	  data_component_interpretation_stress(1,
5865 | 					       DataComponentInterpretation::component_is_scalar);
5866 | 
5867 |     for (unsigned int i = 0; i < dim; ++i)
5868 |       for (unsigned int j = i; j < dim; ++j)
5869 | 	{
5870 | 	  Vector<double> stress_field_L2;
5871 | 	  stress_field_L2.reinit(stresses_dof_handler_L2.n_dofs());
5872 | 
5873 | 	  MappingQ<dim> mapping(m_parameters.m_poly_degree + 1);
5874 | 	  VectorTools::project(mapping,
5875 | 			       stresses_dof_handler_L2,
5876 | 			       constraints,
5877 | 			       m_qf_cell,
5878 | 			       [&] (const typename DoFHandler<dim>::active_cell_iterator & cell,
5879 | 				    const unsigned int q) -> double
5880 | 			       {
5881 | 				 return m_quadrature_point_history.get_data(cell)[q]->get_cauchy_stress()[i][j];
5882 | 			       },
5883 | 			       stress_field_L2);
5884 | 
5885 | 	  std::string stress_name = "Cauchy_stress_" + std::to_string(i+1) + std::to_string(j+1)
5886 | 				  + "_L2";
5887 | 
5888 | 	  data_out.add_data_vector(stresses_dof_handler_L2,
5889 | 				   stress_field_L2,
5890 | 				   stress_name,
5891 | 				   data_component_interpretation_stress);
5892 | 	}
5893 | 
5894 |     data_out.build_patches(m_parameters.m_poly_degree);
5895 | 
5896 |     std::ofstream output("Solution-" + std::to_string(dim) + "d-" +
5897 | 			 Utilities::int_to_string(m_time.get_timestep(),4) + ".vtu");
5898 | 
5899 |     data_out.write_vtu(output);
5900 |     m_timer.leave_subsection();
5901 |   }
5902 | 
5903 |   template <int dim>
5904 |   void PhaseFieldMonolithicSolve<dim>::calculate_reaction_force(unsigned int face_ID)
5905 |   {
5906 |     m_timer.enter_subsection("Calculate reaction force");
5907 | 
5908 |     BlockVector<double>       system_rhs;
5909 |     system_rhs.reinit(m_dofs_per_block);
5910 | 
5911 |     Vector<double> cell_rhs(m_dofs_per_cell);
5912 |     std::vector<types::global_dof_index> local_dof_indices(m_dofs_per_cell);
5913 | 
5914 |     const double time_ramp = (m_time.current() / m_time.end());
5915 |     std::vector<Tensor<1, dim>> rhs_values(m_n_q_points);
5916 |     const UpdateFlags uf_cell(update_values | update_gradients |
5917 | 			      update_quadrature_points | update_JxW_values);
5918 |     const UpdateFlags uf_face(update_values | update_normal_vectors |
5919 |                               update_JxW_values);
5920 | 
5921 |     FEValues<dim> fe_values(m_fe, m_qf_cell, uf_cell);
5922 |     FEFaceValues<dim> fe_face_values(m_fe, m_qf_face, uf_face);
5923 | 
5924 |     // shape function values for displacement field
5925 |     std::vector<std::vector<Tensor<1, dim>>>
5926 |       Nx(m_qf_cell.size(), std::vector<Tensor<1, dim>>(m_dofs_per_cell));
5927 |     std::vector<std::vector<Tensor<2, dim>>>
5928 |       grad_Nx(m_qf_cell.size(), std::vector<Tensor<2, dim>>(m_dofs_per_cell));
5929 |     std::vector<std::vector<SymmetricTensor<2, dim>>>
5930 |       symm_grad_Nx(m_qf_cell.size(), std::vector<SymmetricTensor<2, dim>>(m_dofs_per_cell));
5931 | 
5932 |     for (const auto &cell : m_dof_handler.active_cell_iterators())
5933 |       {
5934 | 	// if calculate_reaction_force() is defined as const, then
5935 | 	// we also need to put a const in std::shared_ptr,
5936 | 	// that is, std::shared_ptr<const PointHistory<dim>>
5937 | 	const std::vector<std::shared_ptr< PointHistory<dim>>> lqph =
5938 | 	  m_quadrature_point_history.get_data(cell);
5939 | 	Assert(lqph.size() == m_n_q_points, ExcInternalError());
5940 |         cell_rhs = 0.0;
5941 |         fe_values.reinit(cell);
5942 |         right_hand_side(fe_values.get_quadrature_points(),
5943 |     		        rhs_values,
5944 |     		        m_parameters.m_x_component*time_ramp,
5945 |     		        m_parameters.m_y_component*time_ramp,
5946 |     		        m_parameters.m_z_component*time_ramp);
5947 | 
5948 |         for (const unsigned int q_point : fe_values.quadrature_point_indices())
5949 |           {
5950 |             for (const unsigned int k : fe_values.dof_indices())
5951 |               {
5952 |                 const unsigned int k_group = m_fe.system_to_base_index(k).first.first;
5953 | 
5954 |                 if (k_group == m_u_dof)
5955 |                   {
5956 |     		    Nx[q_point][k] = fe_values[m_u_fe].value(k, q_point);
5957 |     		    grad_Nx[q_point][k] = fe_values[m_u_fe].gradient(k, q_point);
5958 |     		    symm_grad_Nx[q_point][k] = symmetrize(grad_Nx[q_point][k]);
5959 |                   }
5960 |               }
5961 |           }
5962 | 
5963 |         for (const unsigned int q_point : fe_values.quadrature_point_indices())
5964 |           {
5965 |             const SymmetricTensor<2, dim> & cauchy_stress = lqph[q_point]->get_cauchy_stress();
5966 | 
5967 |             const std::vector<Tensor<1,dim>> & N = Nx[q_point];
5968 |             const std::vector<SymmetricTensor<2, dim>> & symm_grad_N = symm_grad_Nx[q_point];
5969 |             const double JxW = fe_values.JxW(q_point);
5970 | 
5971 |             for (const unsigned int i : fe_values.dof_indices())
5972 |               {
5973 |                 const unsigned int i_group = m_fe.system_to_base_index(i).first.first;
5974 | 
5975 |                 if (i_group == m_u_dof)
5976 |                   {
5977 |                     cell_rhs(i) -= (symm_grad_N[i] * cauchy_stress) * JxW;
5978 |     		    // contributions from the body force to right-hand side
5979 |     		    cell_rhs(i) += N[i] * rhs_values[q_point] * JxW;
5980 |                   }
5981 |               }
5982 |           }
5983 | 
5984 |         // if there is surface pressure, this surface pressure always applied to the
5985 |         // reference configuration
5986 |         const unsigned int face_pressure_id = 100;
5987 |         const double p0 = 0.0;
5988 | 
5989 |         for (const auto &face : cell->face_iterators())
5990 |           {
5991 | 	    if (face->at_boundary() && face->boundary_id() == face_pressure_id)
5992 | 	      {
5993 | 		fe_face_values.reinit(cell, face);
5994 | 
5995 | 		for (const unsigned int f_q_point : fe_face_values.quadrature_point_indices())
5996 | 		  {
5997 | 		    const Tensor<1, dim> &N = fe_face_values.normal_vector(f_q_point);
5998 | 
5999 | 		    const double         pressure  = p0 * time_ramp;
6000 | 		    const Tensor<1, dim> traction  = pressure * N;
6001 | 
6002 | 		    for (const unsigned int i : fe_values.dof_indices())
6003 | 		      {
6004 | 			const unsigned int i_group = m_fe.system_to_base_index(i).first.first;
6005 | 
6006 | 			if (i_group == m_u_dof)
6007 | 			  {
6008 | 			    const unsigned int component_i = m_fe.system_to_component_index(i).first;
6009 | 			    const double Ni = fe_face_values.shape_value(i, f_q_point);
6010 | 			    const double JxW = fe_face_values.JxW(f_q_point);
6011 | 			    cell_rhs(i) += (Ni * traction[component_i]) * JxW;
6012 | 			  }
6013 | 		      }
6014 | 		  }
6015 | 	      }
6016 |           }
6017 | 
6018 |         cell->get_dof_indices(local_dof_indices);
6019 |         for (const unsigned int i : fe_values.dof_indices())
6020 |           system_rhs(local_dof_indices[i]) += cell_rhs(i);
6021 |       } // for (const auto &cell : m_dof_handler.active_cell_iterators())
6022 | 
6023 |     // The difference between the above assembled system_rhs and m_system_rhs
6024 |     // is that m_system_rhs is condensed by the m_constraints, which zero out
6025 |     // the rhs values associated with the constrained DOFs and modify the rhs
6026 |     // values associated with the unconstrained DOFs.
6027 | 
6028 |     std::vector< types::global_dof_index > mapping;
6029 |     std::set<types::boundary_id> boundary_ids;
6030 |     boundary_ids.insert(face_ID);
6031 |     DoFTools::map_dof_to_boundary_indices(m_dof_handler,
6032 | 					  boundary_ids,
6033 | 					  mapping);
6034 | 
6035 |     std::vector<double> reaction_force(dim, 0.0);
6036 | 
6037 |     for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
6038 |       {
6039 | 	if (mapping[i] != numbers::invalid_dof_index)
6040 | 	  {
6041 | 	    reaction_force[i % dim] += system_rhs.block(m_u_dof)(i);
6042 | 	  }
6043 |       }
6044 | 
6045 |     for (unsigned int i = 0; i < dim; i++)
6046 |       m_logfile << "\t\tReaction force in direction " << i << " on boundary ID " << face_ID
6047 |                 << " = "
6048 | 		<< std::fixed << std::setprecision(3) << std::setw(1)
6049 |                 << std::scientific
6050 | 		<< reaction_force[i] << std::endl;
6051 | 
6052 |     std::pair<double, std::vector<double>> time_force;
6053 |     time_force.first = m_time.current();
6054 |     time_force.second = reaction_force;
6055 |     m_history_reaction_force.push_back(time_force);
6056 | 
6057 |     m_timer.leave_subsection();
6058 |   }
6059 | 
6060 |   template <int dim>
6061 |   void PhaseFieldMonolithicSolve<dim>::write_history_data()
6062 |   {
6063 |     m_logfile << "\t\tWrite history data ... \n"<<std::endl;
6064 | 
6065 |     std::ofstream myfile_reaction_force ("Reaction_force.hist");
6066 |     if (myfile_reaction_force.is_open())
6067 |     {
6068 |       myfile_reaction_force << 0.0 << "\t";
6069 |       if (dim == 2)
6070 | 	myfile_reaction_force << 0.0 << "\t"
6071 | 	       << 0.0 << std::endl;
6072 |       if (dim == 3)
6073 | 	myfile_reaction_force << 0.0 << "\t"
6074 | 	       << 0.0 << "\t"
6075 | 	       << 0.0 << std::endl;
6076 | 
6077 |       for (auto const & time_force : m_history_reaction_force)
6078 | 	{
6079 | 	  myfile_reaction_force << time_force.first << "\t";
6080 | 	  if (dim == 2)
6081 | 	    myfile_reaction_force << time_force.second[0] << "\t"
6082 | 	           << time_force.second[1] << std::endl;
6083 | 	  if (dim == 3)
6084 | 	    myfile_reaction_force << time_force.second[0] << "\t"
6085 | 	           << time_force.second[1] << "\t"
6086 | 		   << time_force.second[2] << std::endl;
6087 | 	}
6088 |       myfile_reaction_force.close();
6089 |     }
6090 |     else
6091 |       m_logfile << "Unable to open file";
6092 | 
6093 |     std::ofstream myfile_energy ("Energy.hist");
6094 |     if (myfile_energy.is_open())
6095 |     {
6096 |       myfile_energy << std::fixed << std::setprecision(10) << std::scientific
6097 |                     << 0.0 << "\t"
6098 |                     << 0.0 << "\t"
6099 | 	            << 0.0 << "\t"
6100 | 	            << 0.0 << std::endl;
6101 | 
6102 |       for (auto const & time_energy : m_history_energy)
6103 | 	{
6104 | 	  myfile_energy << std::fixed << std::setprecision(10) << std::scientific
6105 | 	                << time_energy.first     << "\t"
6106 |                         << time_energy.second[0] << "\t"
6107 | 	                << time_energy.second[1] << "\t"
6108 | 		        << time_energy.second[2] << std::endl;
6109 | 	}
6110 |       myfile_energy.close();
6111 |     }
6112 |     else
6113 |       m_logfile << "Unable to open file";
6114 |   }
6115 | 
6116 |   template <int dim>
6117 |   double PhaseFieldMonolithicSolve<dim>::calculate_energy_functional() const
6118 |   {
6119 |     double energy_functional = 0.0;
6120 | 
6121 |     FEValues<dim> fe_values(m_fe, m_qf_cell, update_JxW_values);
6122 | 
6123 |     for (const auto &cell : m_dof_handler.active_cell_iterators())
6124 |       {
6125 |         fe_values.reinit(cell);
6126 | 
6127 |         const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
6128 |           m_quadrature_point_history.get_data(cell);
6129 |         Assert(lqph.size() == m_n_q_points, ExcInternalError());
6130 | 
6131 |         for (unsigned int q_point = 0; q_point < m_n_q_points; ++q_point)
6132 |           {
6133 |             const double JxW = fe_values.JxW(q_point);
6134 |             energy_functional += lqph[q_point]->get_total_strain_energy() * JxW;
6135 |             energy_functional += lqph[q_point]->get_crack_energy_dissipation() * JxW;
6136 |           }
6137 |       }
6138 | 
6139 |     return energy_functional;
6140 |   }
6141 | 
6142 |   template <int dim>
6143 |   std::pair<double, double>
6144 |     PhaseFieldMonolithicSolve<dim>::calculate_total_strain_energy_and_crack_energy_dissipation() const
6145 |   {
6146 |     double total_strain_energy = 0.0;
6147 |     double crack_energy_dissipation = 0.0;
6148 | 
6149 |     FEValues<dim> fe_values(m_fe, m_qf_cell, update_JxW_values);
6150 | 
6151 |     for (const auto &cell : m_dof_handler.active_cell_iterators())
6152 |       {
6153 |         fe_values.reinit(cell);
6154 | 
6155 |         const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
6156 |           m_quadrature_point_history.get_data(cell);
6157 |         Assert(lqph.size() == m_n_q_points, ExcInternalError());
6158 | 
6159 |         for (unsigned int q_point = 0; q_point < m_n_q_points; ++q_point)
6160 |           {
6161 |             const double JxW = fe_values.JxW(q_point);
6162 |             total_strain_energy += lqph[q_point]->get_total_strain_energy() * JxW;
6163 |             crack_energy_dissipation += lqph[q_point]->get_crack_energy_dissipation() * JxW;
6164 |           }
6165 |       }
6166 | 
6167 |     return std::make_pair(total_strain_energy, crack_energy_dissipation);
6168 |   }
6169 | 
6170 |   template <int dim>
6171 |   bool PhaseFieldMonolithicSolve<dim>::local_refine_and_solution_transfer(BlockVector<double> & solution_delta,
6172 | 									  BlockVector<double> & LBFGS_update_refine)
6173 |   {
6174 |     // This is the solution at (n+1) obtained from the old (coarse) mesh
6175 |     BlockVector<double> solution_next_step(m_dofs_per_block);
6176 |     solution_next_step = m_solution + solution_delta;
6177 |     bool mesh_is_same = true;
6178 |     bool cell_refine_flag = true;
6179 | 
6180 |     unsigned int material_id;
6181 |     double length_scale;
6182 |     double cell_length;
6183 |     while(cell_refine_flag)
6184 |       {
6185 | 	cell_refine_flag = false;
6186 | 
6187 | 	std::vector<types::global_dof_index> local_dof_indices(m_fe.dofs_per_cell);
6188 | 	for (const auto &cell : m_dof_handler.active_cell_iterators())
6189 | 	  {
6190 | 	    cell->get_dof_indices(local_dof_indices);
6191 | 
6192 | 	    for (unsigned int i = 0; i< m_fe.dofs_per_cell; ++i)
6193 | 	      {
6194 | 		const unsigned int comp_i = m_fe.system_to_component_index(i).first;
6195 | 		if (comp_i == m_d_component) //phasefield component
6196 | 		  {
6197 | 		    if (  solution_next_step(local_dof_indices[i])
6198 | 			> m_parameters.m_phasefield_refine_threshold )
6199 | 		      {
6200 | 			material_id = cell->material_id();
6201 | 	                length_scale = m_material_data[material_id][2];
6202 | 	                if (dim == 2)
6203 | 	                  cell_length = std::sqrt(cell->measure());
6204 | 	                else
6205 | 	                  cell_length = std::cbrt(cell->measure());
6206 | 			if (  cell_length
6207 | 			    > length_scale * m_parameters.m_allowed_max_h_l_ratio )
6208 | 			  {
6209 | 			    if (cell->level() < m_parameters.m_max_allowed_refinement_level)
6210 | 			      {
6211 | 			        cell->set_refine_flag();
6212 | 			        break;
6213 | 			      }
6214 | 			  }
6215 | 		      }
6216 | 		  }
6217 | 	      }
6218 | 	  }
6219 | 
6220 | 	for (const auto &cell : m_dof_handler.active_cell_iterators())
6221 | 	  {
6222 | 	    if (cell->refine_flag_set())
6223 | 	      {
6224 | 		cell_refine_flag = true;
6225 | 		break;
6226 | 	      }
6227 | 	  }
6228 | 
6229 | 	// if any cell is refined, we need to project the solution
6230 | 	// to the newly refined mesh
6231 | 	if (cell_refine_flag)
6232 | 	  {
6233 | 	    mesh_is_same = false;
6234 | 
6235 | 	    std::vector<BlockVector<double> > old_solutions(2);
6236 | 	    old_solutions[0] = solution_next_step;
6237 | 	    old_solutions[1] = m_solution;
6238 | 
6239 | 	    m_triangulation.prepare_coarsening_and_refinement();
6240 | 	    SolutionTransfer<dim, BlockVector<double>> solution_transfer(m_dof_handler);
6241 | 	    solution_transfer.prepare_for_coarsening_and_refinement(old_solutions);
6242 | 	    m_triangulation.execute_coarsening_and_refinement();
6243 | 
6244 | 	    setup_system();
6245 | 
6246 | 	    std::vector<BlockVector<double>> tmp_solutions(2);
6247 | 	    tmp_solutions[0].reinit(m_dofs_per_block);
6248 | 	    tmp_solutions[1].reinit(m_dofs_per_block);
6249 | 
6250 |             #  if DEAL_II_VERSION_GTE(9, 7, 0)
6251 | 	    solution_transfer.interpolate(tmp_solutions);
6252 | 	    #  else
6253 | 	    // If an older version of dealII is used, for example, 9.4.0, interpolate()
6254 |             // needs to use the following interface.
6255 |             solution_transfer.interpolate(old_solutions, tmp_solutions);
6256 |             #  endif
6257 | 	    solution_next_step = tmp_solutions[0];
6258 | 	    m_solution = tmp_solutions[1];
6259 | 
6260 | 	    // make sure the projected solutions still satisfy
6261 | 	    // hanging node constraints
6262 | 	    m_constraints.distribute(solution_next_step);
6263 | 	    m_constraints.distribute(m_solution);
6264 | 	  } // if (cell_refine_flag)
6265 |       } // while(cell_refine_flag)
6266 | 
6267 |     // calculate field variables for newly refined cells
6268 |     if (!mesh_is_same)
6269 |       {
6270 | 	BlockVector<double> temp_solution_delta(m_dofs_per_block);
6271 | 	BlockVector<double> temp_previous_solution(m_dofs_per_block);
6272 | 	temp_solution_delta = 0.0;
6273 | 	temp_previous_solution = 0.0;
6274 | 	update_qph_incremental(temp_solution_delta, temp_previous_solution);
6275 | 
6276 | 	// initial guess for the resolve on the refined mesh
6277 | 	LBFGS_update_refine = solution_next_step - m_solution;
6278 |       }
6279 | 
6280 |     return mesh_is_same;
6281 |   }
6282 | 
6283 |   template <int dim>
6284 |   void PhaseFieldMonolithicSolve<dim>::print_parameter_information()
6285 |   {
6286 |     if (m_parameters.m_type_nonlinear_solver == "LBFGS")
6287 |       {
6288 | 	m_logfile << "WARNING: this version of LBFGS does not enforce"
6289 | 	    " phase-field irreversibility." << std::endl;
6290 | 	m_logfile << "It should only be used for demonstrating the importance of"
6291 | 	    " inequality constraints." << std::endl;
6292 | 	m_logfile << "The obtained result is meaningless!" << std::endl;
6293 |       }
6294 | 
6295 |     m_logfile << "Scenario number = " << m_parameters.m_scenario << std::endl;
6296 |     m_logfile << "Log file = " << m_parameters.m_logfile_name << std::endl;
6297 |     m_logfile << "Write iteration history to log file? = " << std::boolalpha
6298 | 	      << m_parameters.m_output_iteration_history << std::endl;
6299 | 
6300 |     if (dim == 2)
6301 |       {
6302 | 	if (m_parameters.m_plane_stress)
6303 | 	  m_logfile << "2D plane-stress case" << std::endl;
6304 | 	else
6305 | 	  m_logfile << "2D plane-strain case" << std::endl;
6306 |       }
6307 | 
6308 |     m_logfile << "Nonlinear solver type = " << m_parameters.m_type_nonlinear_solver << std::endl;
6309 |     m_logfile << "Line search type = " << m_parameters.m_type_line_search << std::endl;
6310 |     m_logfile << "Linear solver type = " << m_parameters.m_type_linear_solver << std::endl;
6311 | 
6312 |     if (m_parameters.m_type_linear_solver == "CG")
6313 |       {
6314 |         m_logfile << "Preconditioner type for CG = " << m_parameters.m_type_preconditioner << std::endl;
6315 |         m_logfile << "Convergence tolerance for CG iterations = " << m_parameters.m_CG_tolerace << std::endl;
6316 |       }
6317 | 
6318 |     m_logfile << "Mesh refinement strategy = " << m_parameters.m_refinement_strategy << std::endl;
6319 | 
6320 |     if (m_parameters.m_refinement_strategy == "adaptive-refine")
6321 |       {
6322 | 	m_logfile << "\tMaximum adaptive refinement times allowed in each step = "
6323 | 		  << m_parameters.m_max_adaptive_refine_times << std::endl;
6324 | 	m_logfile << "\tMaximum allowed cell refinement level = "
6325 | 		  << m_parameters.m_max_allowed_refinement_level << std::endl;
6326 | 	m_logfile << "\tPhasefield-based refinement threshold value = "
6327 | 		  << m_parameters.m_phasefield_refine_threshold << std::endl;
6328 |       }
6329 | 
6330 |     m_logfile << "L-BFGS_m = " << m_parameters.m_LBFGS_m << std::endl;
6331 |     m_logfile << "Global refinement times = " << m_parameters.m_global_refine_times << std::endl;
6332 |     m_logfile << "Local prerefinement times = " <<m_parameters. m_local_prerefine_times << std::endl;
6333 |     m_logfile << "Allowed maximum h/l ratio = " << m_parameters.m_allowed_max_h_l_ratio << std::endl;
6334 |     m_logfile << "total number of material types = " << m_parameters.m_total_material_regions << std::endl;
6335 |     m_logfile << "material data file name = " << m_parameters.m_material_file_name << std::endl;
6336 |     if (m_parameters.m_reaction_force_face_id >= 0)
6337 |       m_logfile << "Calculate reaction forces on Face ID = " << m_parameters.m_reaction_force_face_id << std::endl;
6338 |     else
6339 |       m_logfile << "No need to calculate reaction forces." << std::endl;
6340 | 
6341 |     if (m_parameters.m_relative_residual)
6342 |       m_logfile << "Relative residual for convergence." << std::endl;
6343 |     else
6344 |       m_logfile << "Absolute residual for convergence." << std::endl;
6345 | 
6346 |     m_logfile << "Body force = (" << m_parameters.m_x_component << ", "
6347 |                                   << m_parameters.m_y_component << ", "
6348 | 	                          << m_parameters.m_z_component << ") (N/m^3)"
6349 | 				  << std::endl;
6350 | 
6351 |     m_logfile << "End time = " << m_parameters.m_end_time << std::endl;
6352 |     m_logfile << "Time data file name = " << m_parameters.m_time_file_name << std::endl;
6353 |   }
6354 | 
6355 |   template <int dim>
6356 |   void PhaseFieldMonolithicSolve<dim>::run()
6357 |   {
6358 |     print_parameter_information();
6359 | 
6360 |     read_material_data(m_parameters.m_material_file_name,
6361 |     		       m_parameters.m_total_material_regions);
6362 | 
6363 |     std::vector<std::array<double, 4>> time_table;
6364 | 
6365 |     read_time_data(m_parameters.m_time_file_name, time_table);
6366 | 
6367 |     make_grid();
6368 |     setup_system();
6369 |     output_results();
6370 | 
6371 |     m_time.increment(time_table);
6372 | 
6373 |     while(m_time.current() < m_time.end() + m_time.get_delta_t()*1.0e-6)
6374 |       {
6375 | 	m_logfile << std::endl
6376 | 		  << "Timestep " << m_time.get_timestep() << " @ " << m_time.current()
6377 | 		  << 's' << std::endl;
6378 | 
6379 |         bool mesh_is_same = false;
6380 | 
6381 |         // initial guess for the resolve on the refined mesh
6382 | 	BlockVector<double> LBFGS_update_refine(m_dofs_per_block);
6383 | 	LBFGS_update_refine = 0.0;
6384 | 
6385 |         // local adaptive mesh refinement loop
6386 | 	unsigned int adp_refine_iteration = 0;
6387 |         for (; adp_refine_iteration < m_parameters.m_max_adaptive_refine_times + 1; ++adp_refine_iteration)
6388 |           {
6389 | 	    if (m_parameters.m_refinement_strategy == "adaptive-refine")
6390 | 	      m_logfile << "\tAdaptive refinement-"<< adp_refine_iteration << ": " << std::endl;
6391 | 
6392 | 	    BlockVector<double> solution_delta(m_dofs_per_block);
6393 | 	    solution_delta = 0.0;
6394 | 
6395 |             if (m_parameters.m_type_nonlinear_solver == "LBFGS")
6396 | 	      solve_nonlinear_timestep_LBFGS(solution_delta, LBFGS_update_refine);
6397 |             else if (m_parameters.m_type_nonlinear_solver == "LBFGSB")
6398 |               solve_nonlinear_timestep_LBFGS_B(solution_delta, LBFGS_update_refine);
6399 | 	    else
6400 | 	      AssertThrow(false, ExcMessage("Nonlinear solver type not implemented"));
6401 | 
6402 | 	    if (m_parameters.m_refinement_strategy == "adaptive-refine")
6403 | 	      {
6404 | 
6405 | 		if (adp_refine_iteration == m_parameters.m_max_adaptive_refine_times)
6406 | 		  {
6407 | 		    m_solution += solution_delta;
6408 | 		    break;
6409 | 		  }
6410 | 
6411 | 		mesh_is_same = local_refine_and_solution_transfer(solution_delta,
6412 | 								  LBFGS_update_refine);
6413 | 
6414 | 		if (mesh_is_same)
6415 | 		  {
6416 | 		    m_solution += solution_delta;
6417 | 		    break;
6418 | 		  }
6419 | 	      }
6420 | 	    else if (m_parameters.m_refinement_strategy == "pre-refine")
6421 | 	      {
6422 | 		m_solution += solution_delta;
6423 | 	        break;
6424 | 	      }
6425 | 	    else
6426 | 	      {
6427 | 		AssertThrow(false,
6428 | 		            ExcMessage("Selected mesh refinement strategy not implemented!"));
6429 | 	      }
6430 |           } // for (; adp_refine_iteration < m_parameters.m_max_adaptive_refine_times; ++adp_refine_iteration)
6431 | 
6432 |         //AssertThrow(adp_refine_iteration < m_parameters.m_max_adaptive_refine_times,
6433 |         //            ExcMessage("Number of local adaptive mesh refinement exceeds allowed maximum times!"));
6434 | 
6435 | 	// output vtk files every 10 steps if there are too
6436 | 	// many time steps
6437 | 	//if (m_time.get_timestep() % 10 == 0)
6438 |         output_results();
6439 | 
6440 | 	double energy_functional_current = calculate_energy_functional();
6441 | 	m_logfile << "\t\tEnergy functional (J) = " << std::fixed << std::setprecision(10) << std::scientific
6442 | 	          << energy_functional_current << std::endl;
6443 | 
6444 | 	std::pair<double, double> energy_pair = calculate_total_strain_energy_and_crack_energy_dissipation();
6445 | 	m_logfile << "\t\tTotal strain energy (J) = " << std::fixed << std::setprecision(10) << std::scientific
6446 | 		  << energy_pair.first << std::endl;
6447 | 	m_logfile << "\t\tCrack energy dissipation (J) = " << std::fixed << std::setprecision(10) << std::scientific
6448 | 		  << energy_pair.second << std::endl;
6449 | 
6450 | 	std::pair<double, std::array<double, 3>> time_energy;
6451 | 	time_energy.first = m_time.current();
6452 | 	time_energy.second[0] = energy_pair.first;
6453 | 	time_energy.second[1] = energy_pair.second;
6454 | 	time_energy.second[2] = energy_pair.first + energy_pair.second;
6455 | 	m_history_energy.push_back(time_energy);
6456 | 
6457 | 	int face_ID = m_parameters.m_reaction_force_face_id;
6458 | 	if (face_ID >= 0)
6459 | 	  calculate_reaction_force(face_ID);
6460 | 
6461 |         write_history_data();
6462 | 
6463 | 	m_time.increment(time_table);
6464 |       } // while(m_time.current() < m_time.end() + m_time.get_delta_t()*1.0e-6)
6465 |   }
6466 | } // namespace PhaseField
6467 | 
6468 | int main(int argc, char* argv[])
6469 | {
6470 |   using namespace dealii;
6471 | 
6472 |   if (argc != 2)
6473 |     AssertThrow(false,
6474 |     		ExcMessage("The number of arguments provided to the program has to be 2!"));
6475 | 
6476 |   const unsigned int dim = std::stoi(argv[1]);
6477 |   if (dim == 2 )
6478 |     {
6479 |       PhaseField::PhaseFieldMonolithicSolve<2> problem_2D("parameters.prm");
6480 |       problem_2D.run();
6481 |     }
6482 |   else if (dim == 3)
6483 |     {
6484 |       PhaseField::PhaseFieldMonolithicSolve<3> problem_3D("parameters.prm");
6485 |       problem_3D.run();
6486 |     }
6487 |   else
6488 |     {
6489 |       AssertThrow(false,
6490 |                   ExcMessage("Dimension has to be either 2 or 3"));
6491 |     }
6492 | 
6493 |   return 0;
6494 | }
```
