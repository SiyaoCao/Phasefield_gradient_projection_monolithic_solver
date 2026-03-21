# main.cc 运行顺序、逐行号与公式对应详解

> 文档目标：按 **主函数运行顺序** 解释 `main.cc`，并覆盖所有类方法、辅助函数、工具函数。

## 首页：代码总结

`main.cc` 实现了一个基于 [L-BFGS / L-BFGS-B + 梯度投影] 的相场断裂单体（monolithic）求解器。核心流程为：
1. 读取参数、材料与时间表；
2. 按场景构建 2D/3D 网格并建立有限元离散；
3. 在每个时间步执行非线性迭代（含约束投影、线搜索、BFGS 更新）；
4. 进行自适应网格加密与解传递；
5. 输出场变量、能量与反力历史。

**文件规模**：
- `main.cc` 总行数：`6494`
- 文档覆盖：`1~6494` 全部代码行

## 目录

- [1. 主函数运行顺序总览（含代码行号）](#1-主函数运行顺序总览含代码行号)
- [2. 运行期关键计算与公式对应](#2-运行期关键计算与公式对应)
- [3. 全量代码（完整显示，带行号）](#3-全量代码完整显示带行号)
- [4. 按行号详细解读（覆盖全文件，逐行）](#4-按行号详细解读覆盖全文件逐行)
- [5. 全部方法/函数索引（起止行号）](#5-全部方法函数索引起止行号)

## 1. 主函数运行顺序总览（含代码行号）

- **[6468~6494] 程序入口 main(argc, argv)**：读取维度参数，构造 2D/3D 求解器并调用 run()。
- **[2035~2055] 构造函数 PhaseFieldMonolithicSolve(...)**：初始化参数、网格、自由度对象、时间对象、日志、有限元系统和积分规则。
- **[6356~6465] run() 总控流程**：按时间步推进：读材料与时间表、建模、求解、后处理、写历史数据。
- **[1517~1568] read_material_data(...)**：读取材料参数表并存入 m_material_data。
- **[1574~1615] read_time_data(...)**：读取分段时间步长和载荷幅值表。
- **[2057~2097] make_grid() 场景分发**：根据 Scenario number 调用 make_grid_case_1~11。
- **[2100~3279] make_grid_case_1~11**：生成/读取几何网格并设置边界、局部加密等。
- **[3288~3353] setup_system()**：分配 DoF、块结构矩阵/向量、约束与稀疏模式，初始化 QPH。
- **[1620~1667] setup_qph()**：为每个积分点创建 PointHistory 与材料状态。
- **[4471~5806] 非线性求解核心（LBFGS/LBFGSB 与 Cauchy 点）**：构造 BFGS 近似、投影、线搜索、子空间优化并更新解。
- **[3641~4118] 离散系统组装（B0、RHS）**：单元积分并装配切线块矩阵及残量向量。
- **[4129~4368] 线搜索与步长控制**：梯度型线搜索/强 Wolfe 线搜索及 zoom 插值。
- **[5904~6170] 后处理（反力、能量）**：计算反力、总能、裂纹耗散并记录历史。
- **[6171~6281] 自适应网格与解转移**：局部加密、解传递并决定是否继续自适应循环。
- **[5808~5903] output_results()**：输出位移/相场/活动集等 VTK 数据。

## 2. 运行期关键计算与公式对应

> 按“代码在前、公式在后”给出。为避免逐行重复，采用分段一一对应。

### 2.1 [244~257] 退化函数及导数（相场损伤）

对应代码行：

```cpp
// main.cc:244~257
```

对应公式：

[ g(d) = (1-d)^2,\quad g'(d)=2(d-1),\quad g''(d)=2 ]

### 2.2 [684~703] 时间推进（分段时间表）

对应代码行：

```cpp
// main.cc:684~703
```

对应公式：

[ t_{n+1}=t_n+\Delta t_n,\quad m_{n+1}=m(t_n) ]

### 2.3 [1379~1398] 相场点投影（箱约束）

对应代码行：

```cpp
// main.cc:1379~1398
```

对应公式：

[ d_i^{k+1}=\min\left(1,\max\left(d_i^{(n)},\,d_i^{k}+\Delta d_iight)ight) ]

### 2.4 [1401~1444] 断点（break points）计算

对应代码行：

```cpp
// main.cc:1401~1444
```

对应公式：

[ t_i = egin{cases}\dfrac{d_i-1}{g_i}, & g_i<0\[2pt]\dfrac{d_i-d_i^{(n)}}{g_i}, & g_i>0\[2pt]+\infty,& g_i=0\end{cases} ]

### 2.5 [3641~4118] 有限元残量与切线组装

对应代码行：

```cpp
// main.cc:3641~4118
```

对应公式：

[ \mathbf{R}(\mathbf{u},d)=\int_{\Omega} \mathbf{B}^Toldsymbol\sigma\,\mathrm{d}\Omega + \int_{\Omega}\left( g'(d)\psi^+ + \dfrac{g_c}{l}d ight)N_d\,\mathrm{d}\Omega + \int_{\Omega} g_c l\,
abla N_d\cdot
abla d\,\mathrm{d}\Omega ]

### 2.6 [4129~4368] 线搜索（Armijo/Wolfe）

对应代码行：

```cpp
// main.cc:4129~4368
```

对应公式：

[ \phi(lpha)=\Pi(\mathbf{x}+lpha\mathbf{p}),\quad \phi'(lpha)=
abla\Pi(\mathbf{x}+lpha\mathbf{p})^T\mathbf{p} ]

### 2.7 [4196~4311] 强 Wolfe 条件

对应代码行：

```cpp
// main.cc:4196~4311
```

对应公式：

[ \phi(lpha)\le\phi(0)+c_1lpha\phi'(0),\quad |\phi'(lpha)|\le c_2|\phi'(0)| ]

### 2.8 [6117~6169] 总能量与裂纹耗散后处理

对应代码行：

```cpp
// main.cc:6117~6169
```

对应公式：

[ \Pi=\int_{\Omega}\psi(arepsilon,d)\,\mathrm{d}\Omega+\int_{\Omega}g_c\gamma(d,
abla d)\,\mathrm{d}\Omega-\mathcal{W}_{\mathrm{ext}} ]

### 2.9 [5904~6060] 反力积分（边界）

对应代码行：

```cpp
// main.cc:5904~6060
```

对应公式：

[ \mathbf{F}_{\Gamma}=\int_{\Gamma_t}oldsymbol\sigma\mathbf{n}\,\mathrm{d}\Gamma ]

## 3. 全量代码（完整显示，带行号）

> 以下完整粘贴 `main.cc`，并为每行添加行号。

```cpp
   1: /* ---------------------------------------------------------------------
   2:  *
   3:  * Copyright (C) 2006 - 2020 by the deal.II authors
   4:  *
   5:  * This file is part of the deal.II library.
   6:  *
   7:  * The deal.II library is free software; you can use it, redistribute
   8:  * it, and/or modify it under the terms of the GNU Lesser General
   9:  * Public License as published by the Free Software Foundation; either
  10:  * version 2.1 of the License, or (at your option) any later version.
  11:  * The full text of the license can be found in the file LICENSE.md at
  12:  * the top level directory of deal.II.
  13:  *
  14:  * ---------------------------------------------------------------------
  15: 
  16:  *
  17:  * Author: Tao Jin, PhD
  18:  *         University of Ottawa, Ottawa, Ontario, Canada
  19:  *         July 2024
  20:  */
  21: 
  22: /* A monolithic scheme based on the L-BFGS method and the gradient projection method
  23:  *  to solve the phase-field crack problem
  24:  * 1. The phase-field method treats the phasefield irreversibility using the gradient
  25:  *    projection method. During a load step [t_n, t_n+1], let d_n represent the phasefield
  26:  *    at the beginning of the load step (known), then the inequality constraints
  27:  *    d_n <= d_n+1 <= 1.0 are treated as box constraints.
  28:  * 2. Various direct and iterative linear solvers are designed to improve the wall-clock
  29:  *    run time.
  30:  * 3. Using TBB for stiffness assembly and Gauss point calculation.
  31:  * 4. Using adaptive mesh refinement.
  32:  */
  33: 
  34: #include <deal.II/grid/tria.h>
  35: #include <deal.II/grid/grid_generator.h>
  36: #include <deal.II/grid/grid_refinement.h>
  37: #include <deal.II/grid/grid_out.h>
  38: #include <deal.II/grid/grid_in.h>
  39: #include <deal.II/grid/manifold_lib.h>
  40: 
  41: #include <deal.II/dofs/dof_handler.h>
  42: #include <deal.II/dofs/dof_tools.h>
  43: #include <deal.II/dofs/dof_renumbering.h>
  44: 
  45: #include <deal.II/fe/fe_values.h>
  46: #include <deal.II/fe/fe_system.h>
  47: #include <deal.II/fe/fe_q.h>
  48: #include <deal.II/fe/fe_dgp_monomial.h>
  49: #include <deal.II/fe/mapping_q_eulerian.h>
  50: 
  51: #include <deal.II/base/timer.h>
  52: #include <deal.II/base/quadrature_point_data.h>
  53: #include <deal.II/base/parameter_handler.h>
  54: 
  55: #include <deal.II/lac/affine_constraints.h>
  56: #include <deal.II/lac/vector.h>
  57: #include <deal.II/lac/full_matrix.h>
  58: #include <deal.II/lac/sparse_matrix.h>
  59: #include <deal.II/lac/dynamic_sparsity_pattern.h>
  60: #include <deal.II/lac/block_sparse_matrix.h>
  61: #include <deal.II/lac/block_vector.h>
  62: 
  63: #include <deal.II/numerics/vector_tools.h>
  64: #include <deal.II/numerics/matrix_tools.h>
  65: #include <deal.II/numerics/data_out.h>
  66: 
  67: #include <deal.II/lac/solver_cg.h>
  68: #include <deal.II/lac/precondition.h>
  69: #include <deal.II/lac/linear_operator.h>
  70: #include <deal.II/lac/packaged_operation.h>
  71: #include <deal.II/lac/precondition_selector.h>
  72: #include <deal.II/lac/solver_selector.h>
  73: #include <deal.II/lac/sparse_direct.h>
  74: 
  75: #include <deal.II/numerics/error_estimator.h>
  76: 
  77: #include <deal.II/physics/elasticity/standard_tensors.h>
  78: 
  79: #include <deal.II/base/quadrature_point_data.h>
  80: 
  81: #include <deal.II/grid/grid_tools.h>
  82: 
  83: #include <deal.II/base/work_stream.h>
  84: 
  85: #include <deal.II/numerics/solution_transfer.h>
  86: 
  87: #include <deal.II/lac/linear_operator_tools.h>
  88: #include <deal.II/lac/sparse_ilu.h>
  89: 
  90: #include <fstream>
  91: #include <iostream>
  92: 
  93: #include <deal.II/base/logstream.h>
  94: 
  95: #include "SpectrumDecomposition.h"
  96: #include "Utilities.h"
  97: 
  98: namespace PhaseField
  99: {
 100:   using namespace dealii;
 101: 
 102:   template <int dim>
 103:   std::vector<types::global_dof_index> get_vertex_dofs(
 104:     const typename Triangulation<dim>::active_vertex_iterator &vertex,
 105:     const DoFHandler<dim> &dof_handler)
 106:   {
 107:     DoFAccessor<0, dim, dim, false> vertex_dofs(
 108:         &(dof_handler.get_triangulation()),
 109:         vertex->level(),
 110:         vertex->index(),
 111:         &dof_handler);
 112:     const unsigned int n_dofs = dof_handler.get_fe().dofs_per_vertex;
 113:     std::vector<types::global_dof_index> dofs(n_dofs);
 114:     for (unsigned int i = 0; i < n_dofs; ++i)
 115:     {
 116:       dofs[i] = vertex_dofs.vertex_dof_index(0, i);
 117:     }
 118:     return dofs;
 119:   }
 120: 
 121:   // Jacobi preconditioner
 122:   class usr_Jacobi_preconditioner : public Subscriptor
 123:   {
 124:   public:
 125:     usr_Jacobi_preconditioner(const BlockSparseMatrix<double> & S);
 126: 
 127:     void vmult(BlockVector<double> & dst,
 128: 	       const BlockVector<double> & src) const;
 129: 
 130:   private:
 131: #  if DEAL_II_VERSION_GTE(9, 7, 0)
 132:     const ObserverPointer<const BlockSparseMatrix<double> > m_system_matrix;
 133: #  else
 134:     const SmartPointer<const BlockSparseMatrix<double> > m_system_matrix;
 135: #  endif
 136:   };
 137: 
 138:   usr_Jacobi_preconditioner::usr_Jacobi_preconditioner(const BlockSparseMatrix<double> & S)
 139:   : m_system_matrix(&S)
 140:   {}
 141: 
 142:   void usr_Jacobi_preconditioner::vmult(BlockVector<double> & dst,
 143: 					const BlockVector<double> & src) const
 144:   {
 145:     PreconditionJacobi<BlockSparseMatrix<double>> preconditioner;
 146:     preconditioner.initialize(*m_system_matrix, 1.0);
 147: 
 148:     preconditioner.vmult(dst, src);
 149:   }
 150: 
 151:   // LU preconditioner
 152:   class usr_sparseLU_preconditioner : public Subscriptor
 153:   {
 154:   public:
 155:     usr_sparseLU_preconditioner(const SparseDirectUMFPACK & matrix_factorization);
 156: 
 157:     void vmult(BlockVector<double> & dst,
 158: 	       const BlockVector<double> & src) const;
 159: 
 160:   private:
 161: #  if DEAL_II_VERSION_GTE(9, 7, 0)
 162:     const ObserverPointer<const SparseDirectUMFPACK > m_matrix_LU;
 163: #  else
 164:     const SmartPointer<const SparseDirectUMFPACK > m_matrix_LU;
 165: #  endif
 166:   };
 167: 
 168:   usr_sparseLU_preconditioner::usr_sparseLU_preconditioner(const SparseDirectUMFPACK & matrix_factorization)
 169:   : m_matrix_LU(&matrix_factorization)
 170:   {}
 171: 
 172:   void usr_sparseLU_preconditioner::vmult(BlockVector<double> & dst,
 173: 					       const BlockVector<double> & src) const
 174:   {
 175:     (*m_matrix_LU).vmult(dst, src);
 176:   }
 177: 
 178:   // Incomplete LU preconditioner
 179:   class usr_sparseILU_preconditioner : public Subscriptor
 180:   {
 181:   public:
 182:     usr_sparseILU_preconditioner(const SparseILU<double> & ILU_factorization_disp,
 183: 				 const SparseILU<double> & ILU_factorization_phasefield);
 184: 
 185:     void vmult(BlockVector<double> & dst,
 186: 	     const BlockVector<double> & src) const;
 187: 
 188:   private:
 189: #  if DEAL_II_VERSION_GTE(9, 7, 0)
 190:     const ObserverPointer<const SparseILU<double> > m_ILU_factorization_disp;
 191:     const ObserverPointer<const SparseILU<double> > m_ILU_factorization_phasefield;
 192: #  else
 193:     const SmartPointer<const SparseILU<double> > m_ILU_factorization_disp;
 194:     const SmartPointer<const SparseILU<double> > m_ILU_factorization_phasefield;
 195: #  endif
 196:   };
 197: 
 198:   usr_sparseILU_preconditioner::usr_sparseILU_preconditioner(const SparseILU<double> & ILU_factorization_disp,
 199: 							     const SparseILU<double> & ILU_factorization_phasefield)
 200:   : m_ILU_factorization_disp(& ILU_factorization_disp)
 201:   , m_ILU_factorization_phasefield(& ILU_factorization_phasefield)
 202:   {}
 203: 
 204:   void usr_sparseILU_preconditioner::vmult(BlockVector<double> & dst,
 205: 					   const BlockVector<double> & src) const
 206:   {
 207:     std::vector<types::global_dof_index> sizes_per_block(src.n_blocks());
 208:     for (unsigned int i = 0; i < src.n_blocks(); ++i)
 209:       sizes_per_block[i] = src.block(i).size();
 210:     dst.reinit(sizes_per_block);
 211: 
 212:     (*m_ILU_factorization_disp).vmult(dst.block(0), src.block(0));
 213:     (*m_ILU_factorization_phasefield).vmult(dst.block(1), src.block(1));
 214:   }
 215: 
 216:   // body force
 217:   template <int dim>
 218:   void right_hand_side(const std::vector<Point<dim>> &points,
 219: 		       std::vector<Tensor<1, dim>> &  values,
 220: 		       const double fx,
 221: 		       const double fy,
 222: 		       const double fz)
 223:   {
 224:     Assert(values.size() == points.size(),
 225:            ExcDimensionMismatch(values.size(), points.size()));
 226:     Assert(dim >= 2, ExcNotImplemented());
 227: 
 228:     for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
 229:       {
 230: 	if (dim == 2)
 231: 	  {
 232: 	    values[point_n][0] = fx;
 233: 	    values[point_n][1] = fy;
 234: 	  }
 235: 	else
 236: 	  {
 237: 	    values[point_n][0] = fx;
 238: 	    values[point_n][1] = fy;
 239: 	    values[point_n][2] = fz;
 240: 	  }
 241:       }
 242:   }
 243: 
 244:   double degradation_function(const double d)
 245:   {
 246:     return (1.0 - d) * (1.0 - d);
 247:   }
 248: 
 249:   double degradation_function_derivative(const double d)
 250:   {
 251:     return 2.0 * (d - 1.0);
 252:   }
 253: 
 254:   double degradation_function_2nd_order_derivative(const double d)
 255:   {
 256:     (void) d;
 257:     return 2.0;
 258:   }
 259: 
 260:   namespace Parameters
 261:   {
 262:     struct Scenario
 263:     {
 264:       unsigned int m_scenario;
 265:       std::string m_logfile_name;
 266:       bool m_output_iteration_history;
 267:       bool m_plane_stress;
 268:       std::string m_type_nonlinear_solver;
 269:       std::string m_type_line_search;
 270:       std::string m_type_linear_solver;
 271:       std::string m_type_preconditioner;
 272:       double m_CG_tolerace;
 273:       std::string m_refinement_strategy;
 274:       unsigned int m_LBFGS_m;
 275:       unsigned int m_global_refine_times;
 276:       unsigned int m_local_prerefine_times;
 277:       unsigned int m_max_adaptive_refine_times;
 278:       int m_max_allowed_refinement_level;
 279:       double m_phasefield_refine_threshold;
 280:       double m_allowed_max_h_l_ratio;
 281:       unsigned int m_total_material_regions;
 282:       std::string m_material_file_name;
 283:       int m_reaction_force_face_id;
 284: 
 285:       static void declare_parameters(ParameterHandler &prm);
 286:       void parse_parameters(ParameterHandler &prm);
 287:     };
 288: 
 289:     void Scenario::declare_parameters(ParameterHandler &prm)
 290:     {
 291:       prm.enter_subsection("Scenario");
 292:       {
 293:         prm.declare_entry("Scenario number",
 294:                           "1",
 295:                           Patterns::Integer(0),
 296:                           "Geometry, loading and boundary conditions scenario");
 297: 
 298:         prm.declare_entry("Log file name",
 299: 			  "Output.log",
 300:                           Patterns::FileName(Patterns::FileName::input),
 301: 			  "Name of the file for log");
 302: 
 303:         prm.declare_entry("Output iteration history",
 304: 			  "yes",
 305:                           Patterns::Selection("yes|no"),
 306: 			  "Shall we write iteration history to the log file?");
 307: 
 308:         prm.declare_entry("Plane stress",
 309: 			  "no",
 310: 			  Patterns::Selection("yes|no"),
 311: 			  "If it is 2D, is it plane-stress?");
 312: 
 313:         prm.declare_entry("Nonlinear solver type",
 314:                           "LBFGSB",
 315:                           Patterns::Selection("LBFGS|LBFGSB"),
 316:                           "Type of solver used to solve the nonlinear system");
 317: 
 318:         prm.declare_entry("Line search type",
 319:                           "GradientBased",
 320:                           Patterns::Selection("GradientBased|StrongWolfe"),
 321:                           "Type of line search method, the gradient-based method "
 322:                           "should be preferred since it is generally faster");
 323: 
 324:         prm.declare_entry("Linear solver type",
 325:                           "CG",
 326:                           Patterns::Selection("Direct|CG"),
 327:                           "Type of solver used to solve the linear system");
 328: 
 329:         prm.declare_entry("Preconditioner type for CG",
 330:                           "ILU",
 331:                           Patterns::Selection("None|Jacobi|LU|ILU"),
 332:                           "Type of preconditioner used to solve the linear system");
 333: 
 334:         prm.declare_entry("CG tolerance",
 335:                           "1.0e-6",
 336:                           Patterns::Double(0.0),
 337:                           "Convergence tolerance of CG iterations");
 338: 
 339:         prm.declare_entry("Mesh refinement strategy",
 340:                           "adaptive-refine",
 341:                           Patterns::Selection("pre-refine|adaptive-refine"),
 342:                           "Mesh refinement strategy: pre-refine or adaptive-refine");
 343: 
 344:         prm.declare_entry("LBFGS m",
 345:                           "40",
 346:                           Patterns::Integer(0),
 347:                           "Number of vectors used for LBFGS");
 348: 
 349:         prm.declare_entry("Global refinement times",
 350:                           "0",
 351:                           Patterns::Integer(0),
 352:                           "Global refinement times (across the entire domain)");
 353: 
 354:         prm.declare_entry("Local prerefinement times",
 355:                           "0",
 356:                           Patterns::Integer(0),
 357:                           "Local pre-refinement times (assume crack path is known a priori), "
 358:                           "only refine along the crack path.");
 359: 
 360:         prm.declare_entry("Max adaptive refinement times",
 361:                           "100",
 362:                           Patterns::Integer(0),
 363:                           "Maximum number of adaptive refinement times allowed in each step");
 364: 
 365:         prm.declare_entry("Max allowed refinement level",
 366:                           "100",
 367:                           Patterns::Integer(0),
 368:                           "Maximum allowed cell refinement level");
 369: 
 370:         prm.declare_entry("Phasefield refine threshold",
 371: 			  "0.8",
 372: 			  Patterns::Double(),
 373: 			  "Phasefield-based refinement threshold value");
 374: 
 375:         prm.declare_entry("Allowed max hl ratio",
 376: 			  "0.25",
 377: 			  Patterns::Double(),
 378: 			  "Allowed maximum ratio between mesh size h and length scale l");
 379: 
 380:         prm.declare_entry("Material regions",
 381:                           "1",
 382:                           Patterns::Integer(0),
 383:                           "Number of material regions");
 384: 
 385:         prm.declare_entry("Material data file",
 386:                           "1",
 387:                           Patterns::FileName(Patterns::FileName::input),
 388:                           "Material data file");
 389: 
 390:         prm.declare_entry("Reaction force face ID",
 391:                           "1",
 392:                           Patterns::Integer(),
 393:                           "Face id where reaction forces should be calculated "
 394:                           "(negative integer means not to calculate reaction force)");
 395:       }
 396:       prm.leave_subsection();
 397:     }
 398: 
 399:     void Scenario::parse_parameters(ParameterHandler &prm)
 400:     {
 401:       prm.enter_subsection("Scenario");
 402:       {
 403:         m_scenario = prm.get_integer("Scenario number");
 404:         m_logfile_name = prm.get("Log file name");
 405:         m_output_iteration_history = prm.get_bool("Output iteration history");
 406:         m_plane_stress = prm.get_bool("Plane stress");
 407:         m_type_nonlinear_solver = prm.get("Nonlinear solver type");
 408:         m_type_line_search = prm.get("Line search type");
 409:         m_type_linear_solver = prm.get("Linear solver type");
 410:         m_type_preconditioner = prm.get("Preconditioner type for CG");
 411:         m_CG_tolerace = prm.get_double("CG tolerance");
 412:         m_refinement_strategy = prm.get("Mesh refinement strategy");
 413:         m_LBFGS_m = prm.get_integer("LBFGS m");
 414:         m_global_refine_times = prm.get_integer("Global refinement times");
 415:         m_local_prerefine_times = prm.get_integer("Local prerefinement times");
 416:         m_max_adaptive_refine_times = prm.get_integer("Max adaptive refinement times");
 417:         m_max_allowed_refinement_level = prm.get_integer("Max allowed refinement level");
 418:         m_phasefield_refine_threshold = prm.get_double("Phasefield refine threshold");
 419:         m_allowed_max_h_l_ratio = prm.get_double("Allowed max hl ratio");
 420:         m_total_material_regions = prm.get_integer("Material regions");
 421:         m_material_file_name = prm.get("Material data file");
 422:         m_reaction_force_face_id = prm.get_integer("Reaction force face ID");
 423:       }
 424:       prm.leave_subsection();
 425:     }
 426: 
 427:     struct FESystem
 428:     {
 429:       unsigned int m_poly_degree;
 430:       unsigned int m_quad_order;
 431: 
 432:       static void declare_parameters(ParameterHandler &prm);
 433: 
 434:       void parse_parameters(ParameterHandler &prm);
 435:     };
 436: 
 437: 
 438:     void FESystem::declare_parameters(ParameterHandler &prm)
 439:     {
 440:       prm.enter_subsection("Finite element system");
 441:       {
 442:         prm.declare_entry("Polynomial degree",
 443:                           "1",
 444:                           Patterns::Integer(0),
 445:                           "Phase field polynomial order");
 446: 
 447:         prm.declare_entry("Quadrature order",
 448:                           "2",
 449:                           Patterns::Integer(0),
 450:                           "Gauss quadrature order");
 451:       }
 452:       prm.leave_subsection();
 453:     }
 454: 
 455:     void FESystem::parse_parameters(ParameterHandler &prm)
 456:     {
 457:       prm.enter_subsection("Finite element system");
 458:       {
 459:         m_poly_degree = prm.get_integer("Polynomial degree");
 460:         m_quad_order  = prm.get_integer("Quadrature order");
 461:       }
 462:       prm.leave_subsection();
 463:     }
 464: 
 465:     // body force (N/m^3)
 466:     struct BodyForce
 467:     {
 468:       double m_x_component;
 469:       double m_y_component;
 470:       double m_z_component;
 471: 
 472:       static void declare_parameters(ParameterHandler &prm);
 473: 
 474:       void parse_parameters(ParameterHandler &prm);
 475:     };
 476: 
 477:     void BodyForce::declare_parameters(ParameterHandler &prm)
 478:     {
 479:       prm.enter_subsection("Body force");
 480:       {
 481:         prm.declare_entry("Body force x component",
 482: 			  "0.0",
 483: 			  Patterns::Double(),
 484: 			  "Body force x-component (N/m^3)");
 485: 
 486:         prm.declare_entry("Body force y component",
 487: 			  "0.0",
 488: 			  Patterns::Double(),
 489: 			  "Body force y-component (N/m^3)");
 490: 
 491:         prm.declare_entry("Body force z component",
 492: 			  "0.0",
 493: 			  Patterns::Double(),
 494: 			  "Body force z-component (N/m^3)");
 495:       }
 496:       prm.leave_subsection();
 497:     }
 498: 
 499:     void BodyForce::parse_parameters(ParameterHandler &prm)
 500:     {
 501:       prm.enter_subsection("Body force");
 502:       {
 503:         m_x_component = prm.get_double("Body force x component");
 504:         m_y_component = prm.get_double("Body force y component");
 505:         m_z_component = prm.get_double("Body force z component");
 506:       }
 507:       prm.leave_subsection();
 508:     }
 509: 
 510:     struct NonlinearSolver
 511:     {
 512:       unsigned int m_max_iterations_BFGS;
 513:       bool m_relative_residual;
 514: 
 515:       double       m_tol_u_residual;
 516:       double       m_tol_d_residual;
 517:       double       m_tol_u_incr;
 518:       double       m_tol_d_incr;
 519: 
 520:       static void declare_parameters(ParameterHandler &prm);
 521: 
 522:       void parse_parameters(ParameterHandler &prm);
 523:     };
 524: 
 525:     void NonlinearSolver::declare_parameters(ParameterHandler &prm)
 526:     {
 527:       prm.enter_subsection("Nonlinear solver");
 528:       {
 529:         prm.declare_entry("Max iterations BFGS",
 530:                           "20",
 531:                           Patterns::Integer(0),
 532:                           "Number of BFGS iterations allowed");
 533: 
 534:         prm.declare_entry("Relative residual",
 535: 			  "yes",
 536:                           Patterns::Selection("yes|no"),
 537: 			  "Shall we use relative residual for convergence?");
 538: 
 539:         prm.declare_entry("Tolerance displacement residual",
 540:                           "1.0e-9",
 541:                           Patterns::Double(0.0),
 542:                           "Displacement residual tolerance");
 543: 
 544:         prm.declare_entry("Tolerance phasefield residual",
 545:                           "1.0e-9",
 546:                           Patterns::Double(0.0),
 547:                           "Phasefield residual tolerance");
 548: 
 549:         prm.declare_entry("Tolerance displacement increment",
 550:                           "1.0e-9",
 551:                           Patterns::Double(0.0),
 552:                           "Displacement increment tolerance");
 553: 
 554:         prm.declare_entry("Tolerance phasefield increment",
 555:                           "1.0e-9",
 556:                           Patterns::Double(0.0),
 557:                           "Phasefield increment tolerance");
 558:       }
 559:       prm.leave_subsection();
 560:     }
 561: 
 562:     void NonlinearSolver::parse_parameters(ParameterHandler &prm)
 563:     {
 564:       prm.enter_subsection("Nonlinear solver");
 565:       {
 566:         m_max_iterations_BFGS = prm.get_integer("Max iterations BFGS");
 567:         m_relative_residual = prm.get_bool("Relative residual");
 568: 
 569:         m_tol_u_residual           = prm.get_double("Tolerance displacement residual");
 570:         m_tol_d_residual           = prm.get_double("Tolerance phasefield residual");
 571:         m_tol_u_incr               = prm.get_double("Tolerance displacement increment");
 572:         m_tol_d_incr               = prm.get_double("Tolerance phasefield increment");
 573:       }
 574:       prm.leave_subsection();
 575:     }
 576: 
 577:     struct TimeInfo
 578:     {
 579:       double m_end_time;
 580:       std::string m_time_file_name;
 581: 
 582:       static void declare_parameters(ParameterHandler &prm);
 583: 
 584:       void parse_parameters(ParameterHandler &prm);
 585:     };
 586: 
 587:     void TimeInfo::declare_parameters(ParameterHandler &prm)
 588:     {
 589:       prm.enter_subsection("Time");
 590:       {
 591:         prm.declare_entry("End time", "1", Patterns::Double(), "End time");
 592: 
 593:         prm.declare_entry("Time data file",
 594:                           "1",
 595:                           Patterns::FileName(Patterns::FileName::input),
 596:                           "Time data file");
 597:       }
 598:       prm.leave_subsection();
 599:     }
 600: 
 601:     void TimeInfo::parse_parameters(ParameterHandler &prm)
 602:     {
 603:       prm.enter_subsection("Time");
 604:       {
 605:         m_end_time = prm.get_double("End time");
 606:         m_time_file_name = prm.get("Time data file");
 607:       }
 608:       prm.leave_subsection();
 609:     }
 610: 
 611:     struct AllParameters : public Scenario,
 612: 	                   public FESystem,
 613: 	                   public BodyForce,
 614: 			   public NonlinearSolver,
 615: 			   public TimeInfo
 616:     {
 617:       AllParameters(const std::string &input_file);
 618: 
 619:       static void declare_parameters(ParameterHandler &prm);
 620: 
 621:       void parse_parameters(ParameterHandler &prm);
 622:     };
 623: 
 624:     AllParameters::AllParameters(const std::string &input_file)
 625:     {
 626:       ParameterHandler prm;
 627:       declare_parameters(prm);
 628:       prm.parse_input(input_file);
 629:       parse_parameters(prm);
 630:     }
 631: 
 632:     void AllParameters::declare_parameters(ParameterHandler &prm)
 633:     {
 634:       Scenario::declare_parameters(prm);
 635:       FESystem::declare_parameters(prm);
 636:       BodyForce::declare_parameters(prm);
 637:       NonlinearSolver::declare_parameters(prm);
 638:       TimeInfo::declare_parameters(prm);
 639:     }
 640: 
 641:     void AllParameters::parse_parameters(ParameterHandler &prm)
 642:     {
 643:       Scenario::parse_parameters(prm);
 644:       FESystem::parse_parameters(prm);
 645:       BodyForce::parse_parameters(prm);
 646:       NonlinearSolver::parse_parameters(prm);
 647:       TimeInfo::parse_parameters(prm);
 648:     }
 649:   } // namespace Parameters
 650: 
 651:   class Time
 652:   {
 653:   public:
 654:     Time(const double time_end)
 655:       : m_timestep(0)
 656:       , m_time_current(0.0)
 657:       , m_time_end(time_end)
 658:       , m_delta_t(0.0)
 659:       , m_magnitude(1.0)
 660:     {}
 661: 
 662:     virtual ~Time() = default;
 663: 
 664:     double current() const
 665:     {
 666:       return m_time_current;
 667:     }
 668:     double end() const
 669:     {
 670:       return m_time_end;
 671:     }
 672:     double get_delta_t() const
 673:     {
 674:       return m_delta_t;
 675:     }
 676:     double get_magnitude() const
 677:     {
 678:       return m_magnitude;
 679:     }
 680:     unsigned int get_timestep() const
 681:     {
 682:       return m_timestep;
 683:     }
 684:     void increment(std::vector<std::array<double, 4>> time_table)
 685:     {
 686:       double t_1, t_delta, t_magnitude;
 687:       for (auto & time_group : time_table)
 688:         {
 689: 	  t_1 = time_group[1];
 690: 	  t_delta = time_group[2];
 691: 	  t_magnitude = time_group[3];
 692: 
 693: 	  if (m_time_current < t_1 - 1.0e-6*t_delta)
 694: 	    {
 695: 	      m_delta_t = t_delta;
 696: 	      m_magnitude = t_magnitude;
 697: 	      break;
 698: 	    }
 699:         }
 700: 
 701:       m_time_current += m_delta_t;
 702:       ++m_timestep;
 703:     }
 704: 
 705:   private:
 706:     unsigned int m_timestep;
 707:     double       m_time_current;
 708:     const double m_time_end;
 709:     double m_delta_t;
 710:     double m_magnitude;
 711:   };
 712: 
 713:   template <int dim>
 714:   class LinearIsotropicElasticityAdditiveSplit
 715:   {
 716:   public:
 717:     LinearIsotropicElasticityAdditiveSplit(const double lame_lambda,
 718: 			                   const double lame_mu,
 719: 				           const double residual_k,
 720: 					   const double length_scale,
 721: 					   const double viscosity,
 722: 					   const double gc,
 723: 					   const bool   plane_stress_flag)
 724:       : m_lame_lambda(lame_lambda)
 725:       , m_lame_mu(lame_mu)
 726:       , m_residual_k(residual_k)
 727:       , m_length_scale(length_scale)
 728:       , m_eta(viscosity)
 729:       , m_gc(gc)
 730:       , m_plane_stress(plane_stress_flag)
 731:       , m_phase_field_value(0.0)
 732:       , m_grad_phasefield(Tensor<1, dim>())
 733:       , m_strain(SymmetricTensor<2, dim>())
 734:       , m_stress(SymmetricTensor<2, dim>())
 735:       , m_stress_positive(SymmetricTensor<2, dim>())
 736:       , m_mechanical_C(SymmetricTensor<4, dim>())
 737:       , m_strain_energy_positive(0.0)
 738:       , m_strain_energy_negative(0.0)
 739:       , m_strain_energy_total(0.0)
 740:       , m_crack_energy_dissipation(0.0)
 741:     {
 742:       Assert(  ( lame_lambda / (2*(lame_lambda + lame_mu)) <= 0.5)
 743: 	     & ( lame_lambda / (2*(lame_lambda + lame_mu)) >=-1.0),
 744: 	     ExcInternalError() );
 745:     }
 746: 
 747:     const SymmetricTensor<4, dim> & get_mechanical_C() const
 748:     {
 749:       return m_mechanical_C;
 750:     }
 751: 
 752:     const SymmetricTensor<2, dim> & get_cauchy_stress() const
 753:     {
 754:       return m_stress;
 755:     }
 756: 
 757:     const SymmetricTensor<2, dim> & get_cauchy_stress_positive() const
 758:     {
 759:       return m_stress_positive;
 760:     }
 761: 
 762:     double get_positive_strain_energy() const
 763:     {
 764:       return m_strain_energy_positive;
 765:     }
 766: 
 767:     double get_negative_strain_energy() const
 768:     {
 769:       return m_strain_energy_negative;
 770:     }
 771: 
 772:     double get_total_strain_energy() const
 773:     {
 774:       return m_strain_energy_total;
 775:     }
 776: 
 777:     double get_crack_energy_dissipation() const
 778:     {
 779:       return m_crack_energy_dissipation;
 780:     }
 781: 
 782:     double get_phase_field_value() const
 783:     {
 784:       return m_phase_field_value;
 785:     }
 786: 
 787:     const Tensor<1, dim> get_phase_field_gradient() const
 788:     {
 789:       return m_grad_phasefield;
 790:     }
 791: 
 792:     void update_material_data(const SymmetricTensor<2, dim> & strain,
 793: 			      const double phase_field_value,
 794: 			      const Tensor<1, dim> & grad_phasefield,
 795: 			      const double phase_field_value_previous_step,
 796: 			      const double delta_time);
 797: 
 798:   private:
 799:     const double m_lame_lambda;
 800:     const double m_lame_mu;
 801:     const double m_residual_k;
 802:     const double m_length_scale;
 803:     const double m_eta;
 804:     const double m_gc;
 805:     const bool   m_plane_stress;
 806:     double m_phase_field_value;
 807:     Tensor<1, dim> m_grad_phasefield;
 808:     SymmetricTensor<2, dim> m_strain;
 809:     SymmetricTensor<2, dim> m_stress;
 810:     SymmetricTensor<2, dim> m_stress_positive;
 811:     SymmetricTensor<4, dim> m_mechanical_C;
 812:     double m_strain_energy_positive;
 813:     double m_strain_energy_negative;
 814:     double m_strain_energy_total;
 815:     double m_crack_energy_dissipation;
 816:   };
 817: 
 818:   template <int dim>
 819:   void LinearIsotropicElasticityAdditiveSplit<dim>::
 820:    update_material_data(const SymmetricTensor<2, dim> & strain,
 821: 			const double phase_field_value,
 822: 			const Tensor<1, dim> & grad_phasefield,
 823: 			const double phase_field_value_previous_step,
 824: 			const double delta_time)
 825:   {
 826:     m_strain = strain;
 827:     m_phase_field_value = phase_field_value;
 828:     m_grad_phasefield = grad_phasefield;
 829:     Vector<double>              eigenvalues(dim);
 830:     std::vector<Tensor<1, dim>> eigenvectors(dim);
 831:     usr_spectrum_decomposition::spectrum_decomposition<dim>(m_strain,
 832:   							      eigenvalues,
 833:   							      eigenvectors);
 834: 
 835:     SymmetricTensor<2, dim> strain_positive, strain_negative;
 836:     strain_positive = usr_spectrum_decomposition::positive_tensor(eigenvalues, eigenvectors);
 837:     strain_negative = usr_spectrum_decomposition::negative_tensor(eigenvalues, eigenvectors);
 838: 
 839:     SymmetricTensor<4, dim> projector_positive, projector_negative;
 840:     usr_spectrum_decomposition::positive_negative_projectors(eigenvalues,
 841:   							       eigenvectors,
 842: 							       projector_positive,
 843: 							       projector_negative);
 844: 
 845:     SymmetricTensor<2, dim> stress_positive, stress_negative;
 846:     const double degradation = degradation_function(m_phase_field_value) + m_residual_k;
 847:     const double I_1 = trace(m_strain);
 848: 
 849:     // 2D plane strain and 3D cases
 850:     double my_lambda = m_lame_lambda;
 851: 
 852:     // 2D plane stress case
 853:     if (    dim == 2
 854: 	   && m_plane_stress)
 855:       my_lambda = 2 * m_lame_mu * m_lame_lambda / (m_lame_lambda + 2 * m_lame_mu);
 856: 
 857:     stress_positive = my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)
 858:                                     * Physics::Elasticity::StandardTensors<dim>::I
 859:                     + 2 * m_lame_mu * strain_positive;
 860:     stress_negative = my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)
 861:                                     * Physics::Elasticity::StandardTensors<dim>::I
 862:     		      + 2 * m_lame_mu * strain_negative;
 863: 
 864:     m_stress = degradation * stress_positive + stress_negative;
 865:     m_stress_positive = stress_positive;
 866: 
 867:     SymmetricTensor<4, dim> C_positive, C_negative;
 868:     C_positive = my_lambda * usr_spectrum_decomposition::heaviside_function(I_1)
 869:                                * Physics::Elasticity::StandardTensors<dim>::IxI
 870: 		 + 2 * m_lame_mu * projector_positive;
 871:     C_negative = my_lambda * usr_spectrum_decomposition::heaviside_function(-I_1)
 872:                                * Physics::Elasticity::StandardTensors<dim>::IxI
 873:     		 + 2 * m_lame_mu * projector_negative;
 874:     m_mechanical_C = degradation * C_positive + C_negative;
 875: 
 876:     m_strain_energy_positive = 0.5 * my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)
 877:                                                    * usr_spectrum_decomposition::positive_ramp_function(I_1)
 878:                              + m_lame_mu * strain_positive * strain_positive;
 879: 
 880:     m_strain_energy_negative = 0.5 * my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)
 881:                                                    * usr_spectrum_decomposition::negative_ramp_function(I_1)
 882:                              + m_lame_mu * strain_negative * strain_negative;
 883: 
 884:     m_strain_energy_total = degradation * m_strain_energy_positive + m_strain_energy_negative;
 885: 
 886:     m_crack_energy_dissipation = m_gc * (  0.5 / m_length_scale * m_phase_field_value * m_phase_field_value
 887: 	                                   + 0.5 * m_length_scale * m_grad_phasefield * m_grad_phasefield)
 888: 	                                   // the term due to viscosity regularization
 889: 	                                   + (m_phase_field_value - phase_field_value_previous_step)
 890: 					   * (m_phase_field_value - phase_field_value_previous_step)
 891: 				           * 0.5 * m_eta / delta_time;
 892:     //(void)delta_time;
 893:     //(void)phase_field_value_previous_step;
 894:   }
 895: 
 896:   template <int dim>
 897:   class PointHistory
 898:   {
 899:   public:
 900:     PointHistory()
 901:       : m_length_scale(0.0)
 902:       , m_gc(0.0)
 903:       , m_viscosity(0.0)
 904:     {}
 905: 
 906:     virtual ~PointHistory() = default;
 907: 
 908:     void setup_lqp(const double lame_lambda,
 909: 		   const double lame_mu,
 910: 		   const double length_scale,
 911: 		   const double gc,
 912: 		   const double viscosity,
 913: 		   const double residual_k,
 914: 		   const bool   plane_stress_flag)
 915:     {
 916:       m_material =
 917:               std::make_shared<LinearIsotropicElasticityAdditiveSplit<dim>>(lame_lambda,
 918:         	                                                            lame_mu,
 919: 								            residual_k,
 920: 									    length_scale,
 921: 									    viscosity,
 922: 									    gc,
 923: 									    plane_stress_flag);
 924:       m_length_scale = length_scale;
 925:       m_gc = gc;
 926:       m_viscosity = viscosity;
 927: 
 928:       update_field_values(SymmetricTensor<2, dim>(), 0.0, Tensor<1, dim>(), 0.0, 1.0);
 929:     }
 930: 
 931:     void update_field_values(const SymmetricTensor<2, dim> & strain,
 932: 		             const double phase_field_value,
 933: 			     const Tensor<1, dim> & grad_phasefield,
 934: 			     const double phase_field_value_previous_step,
 935: 			     const double delta_time)
 936:     {
 937:       m_material->update_material_data(strain, phase_field_value, grad_phasefield,
 938: 				       phase_field_value_previous_step, delta_time);
 939:     }
 940: 
 941:     double get_current_positive_strain_energy() const
 942:     {
 943:       return m_material->get_positive_strain_energy();
 944:     }
 945: 
 946:     const SymmetricTensor<4, dim> & get_mechanical_C() const
 947:     {
 948:       return m_material->get_mechanical_C();
 949:     }
 950: 
 951:     const SymmetricTensor<2, dim> & get_cauchy_stress() const
 952:     {
 953:       return m_material->get_cauchy_stress();
 954:     }
 955: 
 956:     const SymmetricTensor<2, dim> & get_cauchy_stress_positive() const
 957:     {
 958:       return m_material->get_cauchy_stress_positive();
 959:     }
 960: 
 961:     double get_total_strain_energy() const
 962:     {
 963:       return m_material->get_total_strain_energy();
 964:     }
 965: 
 966:     double get_crack_energy_dissipation() const
 967:     {
 968:       return m_material->get_crack_energy_dissipation();
 969:     }
 970: 
 971:     double get_phase_field_value() const
 972:     {
 973:       return m_material->get_phase_field_value();
 974:     }
 975: 
 976:     const Tensor<1, dim> get_phase_field_gradient() const
 977:     {
 978:       return m_material->get_phase_field_gradient();
 979:     }
 980: 
 981:     double get_length_scale() const
 982:     {
 983:       return m_length_scale;
 984:     }
 985: 
 986:     double get_critical_energy_release_rate() const
 987:     {
 988:       return m_gc;
 989:     }
 990: 
 991:     double get_viscosity() const
 992:     {
 993:       return m_viscosity;
 994:     }
 995:   private:
 996:     std::shared_ptr<LinearIsotropicElasticityAdditiveSplit<dim>> m_material;
 997:     double m_length_scale;
 998:     double m_gc;
 999:     double m_viscosity;
1000:   };
1001: 
1002:   template <int dim>
1003:   class PhaseFieldMonolithicSolve
1004:   {
1005:   public:
1006:     PhaseFieldMonolithicSolve(const std::string &input_file);
1007: 
1008:     virtual ~PhaseFieldMonolithicSolve() = default;
1009:     void run();
1010: 
1011:   private:
1012:     struct PerTaskData_ASM;
1013:     struct ScratchData_ASM;
1014: 
1015:     struct PerTaskData_ASM_RHS_BFGS;
1016:     struct ScratchData_ASM_RHS_BFGS;
1017: 
1018:     struct PerTaskData_UQPH;
1019:     struct ScratchData_UQPH;
1020: 
1021:     Parameters::AllParameters m_parameters;
1022:     Triangulation<dim> m_triangulation;
1023: 
1024:     CellDataStorage<typename Triangulation<dim>::cell_iterator,
1025:                     PointHistory<dim>>
1026:       m_quadrature_point_history;
1027: 
1028:     Time                m_time;
1029:     std::ofstream m_logfile;
1030:     mutable TimerOutput m_timer;
1031: 
1032:     DoFHandler<dim>                  m_dof_handler;
1033:     FESystem<dim>                    m_fe;
1034:     const unsigned int               m_dofs_per_cell;
1035:     const FEValuesExtractors::Vector m_u_fe;
1036:     const FEValuesExtractors::Scalar m_d_fe;
1037: 
1038:     static const unsigned int m_n_blocks          = 2;
1039:     static const unsigned int m_n_components      = dim + 1;
1040:     static const unsigned int m_first_u_component = 0;
1041:     static const unsigned int m_d_component       = dim;
1042: 
1043:     enum
1044:     {
1045:       m_u_dof = 0,
1046:       m_d_dof = 1
1047:     };
1048: 
1049:     std::vector<types::global_dof_index> m_dofs_per_block;
1050: 
1051:     const QGauss<dim>     m_qf_cell;
1052:     const QGauss<dim - 1> m_qf_face;
1053:     const unsigned int    m_n_q_points;
1054: 
1055:     double m_vol_reference;
1056: 
1057:     AffineConstraints<double> m_constraints;
1058:     BlockSparsityPattern      m_sparsity_pattern;
1059:     BlockSparseMatrix<double> m_tangent_matrix;
1060:     BlockVector<double>       m_system_rhs;
1061:     BlockVector<double>       m_solution;
1062:     SparseDirectUMFPACK       m_A_direct;
1063:     // m_active_set_phasefield has 0 (inactive constraint)
1064:     //                          or 1 (active constraint lower bound)
1065:     //                          or 2 (active constraint upper bound)
1066:     // In order to add active set into the VTK output, we have to declare
1067:     // it as double, not int or unsigned int
1068:     Vector<double> m_active_set_phasefield;
1069: 
1070:     std::map<unsigned int, std::vector<double>> m_material_data;
1071: 
1072:     std::vector<std::pair<double, std::vector<double>>> m_history_reaction_force;
1073:     std::vector<std::pair<double, std::array<double, 3>>> m_history_energy;
1074: 
1075: 
1076:     struct Errors
1077:     {
1078:       Errors()
1079:         : m_norm(1.0)
1080:         , m_u(1.0)
1081:         , m_d(1.0)
1082:       {}
1083: 
1084:       void reset()
1085:       {
1086:         m_norm = 1.0;
1087:         m_u    = 1.0;
1088:         m_d    = 1.0;
1089:       }
1090: 
1091:       void normalize(const Errors &rhs)
1092:       {
1093:         if (rhs.m_norm != 0.0)
1094:           m_norm /= rhs.m_norm;
1095:         if (rhs.m_u != 0.0)
1096:           m_u /= rhs.m_u;
1097:         if (rhs.m_d != 0.0)
1098:           m_d /= rhs.m_d;
1099:       }
1100: 
1101:       double m_norm, m_u, m_d;
1102:     };
1103: 
1104:     Errors m_error_residual, m_error_residual_0, m_error_residual_norm, m_error_update,
1105:       m_error_update_0, m_error_update_norm;
1106: 
1107:     void get_error_residual(Errors &error_residual);
1108:     void get_error_residual_LBFGSB(Errors &error_residual,
1109: 				   const BlockVector<double> & solution_delta);
1110: 
1111:     void get_error_update(const BlockVector<double> &newton_update,
1112:                           Errors & error_update);
1113: 
1114:     void make_grid();
1115:     void make_grid_case_1();
1116:     void make_grid_case_2();
1117:     void make_grid_case_3();
1118:     void make_grid_case_4();
1119:     void make_grid_case_5();
1120:     void make_grid_case_6();
1121:     void make_grid_case_7();
1122:     void make_grid_case_8();
1123:     void make_grid_case_9();
1124:     void make_grid_case_10();
1125:     void make_grid_case_11();
1126: 
1127:     void setup_system();
1128: 
1129:     void determine_component_extractors();
1130: 
1131:     void make_constraints(const unsigned int it_nr);
1132: 
1133:     void assemble_system_B0(const BlockVector<double> & solution_old);
1134: 
1135:     void assemble_system_B0_one_cell(
1136:       const typename DoFHandler<dim>::active_cell_iterator &cell,
1137:       ScratchData_ASM &                                     scratch,
1138:       PerTaskData_ASM &                                     data) const;
1139: 
1140:     void assemble_system_rhs_BFGS_one_cell(
1141:       const typename DoFHandler<dim>::active_cell_iterator &cell,
1142:       ScratchData_ASM_RHS_BFGS &                           scratch,
1143:       PerTaskData_ASM_RHS_BFGS &                           data) const;
1144: 
1145:     void assemble_system_rhs_BFGS(const BlockVector<double> & solution_old,
1146: 				  BlockVector<double> & system_rhs);
1147: 
1148:     void assemble_system_rhs_BFGS_parallel(const BlockVector<double> & solution_old,
1149:     				           BlockVector<double> & system_rhs);
1150: 
1151:     void solve_nonlinear_timestep_LBFGS(BlockVector<double> &solution_delta,
1152: 				        BlockVector<double> & LBFGS_update_refine);
1153: 
1154:     void solve_nonlinear_timestep_LBFGS_B(BlockVector<double> &solution_delta,
1155:     				          BlockVector<double> & LBFGS_update_refine);
1156: 
1157:     void calculate_cauchy_point(const BlockSparseMatrix<double> & B0_matrix,
1158: 	                        const std::list<BlockVector<double>> & y_vector_list,
1159: 				const std::list<BlockVector<double>> & b0xs_vector_list,
1160: 				const FullMatrix<double> & M_matrix,
1161: 				const BlockVector<double> & gradient_g,
1162: 				const BlockVector<double> & solution_delta,
1163: 				BlockVector<double> & solution_delta_cauchy_point);
1164: 
1165:     double line_search_stepsize_gradient_based(const BlockVector<double> & BFGS_p_vector,
1166: 					       const BlockVector<double> & solution_delta);
1167: 
1168:     double line_search_stepsize_strong_wolfe(const double phi_0,
1169: 				             const double phi_0_prime,
1170: 				             const BlockVector<double> & BFGS_p_vector,
1171: 				             const BlockVector<double> & solution_delta);
1172: 
1173:     double line_search_zoom_strong_wolfe(double phi_low, double phi_low_prime, double alpha_low,
1174: 					 double phi_high, double phi_high_prime, double alpha_high,
1175: 					 double phi_0, double phi_0_prime, const BlockVector<double> & BFGS_p_vector,
1176: 					 double c1, double c2, unsigned int max_iter,
1177: 					 const BlockVector<double> & solution_delta);
1178: 
1179:     double line_search_interpolation_cubic(const double alpha_0, const double phi_0, const double phi_0_prime,
1180: 					   const double alpha_1, const double phi_1, const double phi_1_prime);
1181: 
1182:     std::pair<double, double> calculate_phi_and_phi_prime(const double alpha,
1183: 							  const BlockVector<double> & BFGS_p_vector,
1184: 							  const BlockVector<double> & solution_delta);
1185: 
1186:     void LBFGS_B0(BlockVector<double> & LBFGS_r_vector,
1187: 		  BlockVector<double> & LBFGS_q_vector);
1188: 
1189:     void output_results() const;
1190: 
1191:     void setup_qph();
1192: 
1193:     void update_qph_incremental(const BlockVector<double> &solution_delta,
1194: 				const BlockVector<double> &solution_old);
1195: 
1196:     void update_qph_incremental_one_cell(
1197:       const typename DoFHandler<dim>::active_cell_iterator &cell,
1198:       ScratchData_UQPH &                                    scratch,
1199:       PerTaskData_UQPH &                                    data);
1200: 
1201:     void copy_local_to_global_UQPH(const PerTaskData_UQPH & /*data*/)
1202:     {}
1203: 
1204:     BlockVector<double>
1205:     get_total_solution(const BlockVector<double> &solution_delta) const;
1206: 
1207:     // Should not make this function const
1208:     void read_material_data(const std::string &data_file,
1209: 			    const unsigned int total_material_regions);
1210: 
1211:     void read_time_data(const std::string &data_file,
1212:     		        std::vector<std::array<double, 4>> & time_table);
1213: 
1214:     void print_conv_header_LBFGS();
1215: 
1216:     void print_conv_header_LBFGSB();
1217: 
1218:     void print_parameter_information();
1219: 
1220:     void calculate_reaction_force(unsigned int face_ID);
1221: 
1222:     void write_history_data();
1223: 
1224:     double calculate_energy_functional() const;
1225: 
1226:     std::pair<double, double> calculate_total_strain_energy_and_crack_energy_dissipation() const;
1227: 
1228:     bool local_refine_and_solution_transfer(BlockVector<double> & solution_delta,
1229: 					    BlockVector<double> & LBFGS_update_refine);
1230: 
1231:     // L-BFGS-B subroutines
1232:     void point_projection(BlockVector<double> & solution_delta);
1233: 
1234:     std::priority_queue< std::pair<double, unsigned int>,
1235:                          std::vector<std::pair<double, unsigned int>>,
1236:     		         std::greater<std::pair<double, unsigned int>> >
1237:       calculate_break_points(const BlockVector<double> & solution_delta,
1238: 			     const BlockVector<double> & gradient_g,
1239: 			     BlockVector<double> & gradient_d);
1240: 
1241:     double ebT_x_B0_x_v(const unsigned int b,
1242: 			const BlockSparseMatrix<double> & B0_matrix,
1243: 			const BlockVector<double> & v);
1244: 
1245:     void zT_x_vector(const BlockVector<double> & z,
1246: 		     const BlockVector<double> & src_vector,
1247: 		     BlockVector<double> & target_vector);
1248: 
1249:     void z_x_vector(const BlockVector<double> & z,
1250: 		    const BlockVector<double> & src_vector,
1251: 		    BlockVector<double> & target_vector);
1252: 
1253:     void zT_B0_z(const BlockVector<double> & z,
1254: 		     BlockSparseMatrix<double> & B0_matrix);
1255: 
1256: 
1257:   }; // class PhaseFieldSplitSolve
1258: 
1259:   template <int dim>
1260:   void PhaseFieldMonolithicSolve<dim>::zT_B0_z(const BlockVector<double> & z,
1261: 					       BlockSparseMatrix<double> & B0_matrix)
1262:   {
1263:     // block 1: displacement
1264:     for (unsigned int i = 0; i < z.block(m_u_dof).size(); ++i)
1265:       {
1266: 	if (z.block(m_u_dof)[i] < 0)
1267: 	  {
1268: 	    for (auto itr = B0_matrix.block(m_u_dof, m_u_dof).begin(i);
1269: 		      itr != B0_matrix.block(m_u_dof, m_u_dof).end(i);
1270: 		      ++itr)
1271: 	      {
1272: 		if (itr->column() != itr->row())
1273: 		  {
1274: 		    B0_matrix.block(m_u_dof, m_u_dof).set(itr->row(),    itr->column(), 0.0);
1275: 		    B0_matrix.block(m_u_dof, m_u_dof).set(itr->column(), itr->row(),    0.0);
1276: 		  }
1277: 	      }
1278: 	  }
1279:       } // for (unsigned int i = 0; i < z.block(m_u_dof).size(); ++i)
1280: 
1281:     // block 2: phasefield
1282:     for (unsigned int i = 0; i < z.block(m_d_dof).size(); ++i)
1283:       {
1284: 	if (z.block(m_d_dof)[i] < 0)
1285: 	  {
1286: 	    for (auto itr = B0_matrix.block(m_d_dof, m_d_dof).begin(i);
1287: 		      itr != B0_matrix.block(m_d_dof, m_d_dof).end(i);
1288: 		      ++itr)
1289: 	      {
1290: 		if (itr->column() != itr->row())
1291: 		  {
1292: 		    B0_matrix.block(m_d_dof, m_d_dof).set(itr->row(),    itr->column(), 0.0);
1293: 		    B0_matrix.block(m_d_dof, m_d_dof).set(itr->column(), itr->row(),    0.0);
1294: 		  }
1295: 	      }
1296: 	  }
1297:       } // for (unsigned int i = 0; i < z.block(m_d_dof).size(); ++i)
1298:   }
1299: 
1300:   template <int dim>
1301:   void PhaseFieldMonolithicSolve<dim>::z_x_vector(const BlockVector<double> & z,
1302: 						  const BlockVector<double> & src_vector,
1303: 						  BlockVector<double> & target_vector)
1304:   {
1305:     //We assume that the dimensions of all the block matrices are correct
1306:     //block 1: displacement
1307:     unsigned int target_vector_index = 0;
1308:     for (unsigned int i = 0; i < z.block(m_u_dof).size(); ++i)
1309:       {
1310: 	if (z.block(m_u_dof)[i] > 0)
1311: 	  {
1312: 	    target_vector.block(m_u_dof)[i] = src_vector.block(m_u_dof)[target_vector_index];
1313: 	    ++target_vector_index;
1314: 	  }
1315: 	else
1316: 	  target_vector.block(m_u_dof)[i] = 0;
1317:       }
1318: 
1319:     //block 2: phasefield
1320:     target_vector_index = 0;
1321:     for (unsigned int i = 0; i < z.block(m_d_dof).size(); ++i)
1322:       {
1323: 	if (z.block(m_d_dof)[i] > 0)
1324: 	  {
1325: 	    target_vector.block(m_d_dof)[i] = src_vector.block(m_d_dof)[target_vector_index];
1326: 	    ++target_vector_index;
1327: 	  }
1328: 	else
1329: 	  target_vector.block(m_d_dof)[i] = 0;
1330:       }
1331:   }
1332: 
1333:   template <int dim>
1334:   void PhaseFieldMonolithicSolve<dim>::zT_x_vector(const BlockVector<double> & z,
1335: 						   const BlockVector<double> & src_vector,
1336: 						   BlockVector<double> & target_vector)
1337:   {
1338:     //We assume that the dimensions of all the block matrices are correct
1339:     //block 1: displacement
1340:     unsigned int target_vector_index = 0;
1341:     for (unsigned int i = 0; i < z.block(m_u_dof).size(); ++i)
1342:       {
1343: 	if (z.block(m_u_dof)[i] > 0)
1344: 	  {
1345: 	    target_vector.block(m_u_dof)[target_vector_index] = src_vector.block(m_u_dof)[i];
1346: 	    ++target_vector_index;
1347: 	  }
1348:       }
1349: 
1350:     //block 2: phasefield
1351:     target_vector_index = 0;
1352:     for (unsigned int i = 0; i < z.block(m_d_dof).size(); ++i)
1353:       {
1354:         if (z.block(m_d_dof)[i] > 0)
1355:         {
1356: 	  target_vector.block(m_d_dof)[target_vector_index] = src_vector.block(m_d_dof)[i];
1357: 	  ++target_vector_index;
1358:         }
1359:       }
1360:   }
1361: 
1362:   template <int dim>
1363:   double PhaseFieldMonolithicSolve<dim>::ebT_x_B0_x_v(const unsigned int b,
1364: 						      const BlockSparseMatrix<double> & B0_matrix,
1365: 						      const BlockVector<double> & v)
1366:   {
1367:     double row_sum = 0.0;
1368:     for (auto itr = B0_matrix.block(m_d_dof, m_d_dof).begin(b);
1369:               itr != B0_matrix.block(m_d_dof, m_d_dof).end(b);
1370:               ++itr)
1371:       {
1372:         row_sum += itr->value() * v.block(m_d_dof)(itr->column());
1373:       }
1374: 
1375:     return row_sum;
1376:   }
1377: 
1378:   template <int dim>
1379:   void PhaseFieldMonolithicSolve<dim>::point_projection(BlockVector<double> & solution_delta)
1380:   {
1381:     // Phase-field value cannot exceed 1.0
1382:     const double upper_limit = 1.0;
1383: 
1384:     Vector<double> solution_phasefield_total(m_solution.block(m_d_dof));
1385:     solution_phasefield_total += solution_delta.block(m_d_dof);
1386: 
1387:     for (unsigned int i = 0; i < solution_phasefield_total.size(); ++i)
1388:       {
1389: 	if (solution_delta.block(m_d_dof)[i] < 0.0)
1390: 	  solution_delta.block(m_d_dof)[i] = 0.0;
1391: 
1392: 	if (solution_phasefield_total[i] > upper_limit)
1393: 	  solution_delta.block(m_d_dof)[i] = upper_limit - m_solution.block(m_d_dof)[i];
1394:       }
1395:   }
1396: 
1397:   template <int dim>
1398:   std::priority_queue< std::pair<double, unsigned int>,
1399:                        std::vector<std::pair<double, unsigned int>>,
1400:   		       std::greater<std::pair<double, unsigned int>> >
1401:     PhaseFieldMonolithicSolve<dim>::calculate_break_points(const BlockVector<double> & solution_delta,
1402: 		       				           const BlockVector<double> & gradient_g,
1403: 							   BlockVector<double> & gradient_d)
1404:   {
1405:     // Creates a min heap of break points
1406:     std::priority_queue< std::pair<double, unsigned int>,
1407:                          std::vector<std::pair<double, unsigned int>>,
1408:     		         std::greater<std::pair<double, unsigned int>> >
1409:     break_points_sorted;
1410: 
1411:     double t = 0.0;
1412: 
1413:     Vector<double> solution_phasefield_total(m_solution.block(m_d_dof));
1414:     solution_phasefield_total += solution_delta.block(m_d_dof);
1415: 
1416:     // upper bound is 1.0, lower bound is the solution at the previous step.
1417:     for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
1418:       {
1419: 	if (gradient_g.block(m_d_dof)[i] < 0)
1420: 	  t = (solution_phasefield_total[i] - 1.0 ) / gradient_g.block(m_d_dof)[i];
1421: 	else if (gradient_g.block(m_d_dof)[i] > 0)
1422: 	  t = solution_delta.block(m_d_dof)[i] / gradient_g.block(m_d_dof)[i];
1423: 	else
1424: 	  t = std::numeric_limits<double>::max();
1425: 
1426:         //AssertThrow(t >= 0, ExcMessage("Break point has to be a non-negative t value"));
1427: 
1428:         if (t > 0)
1429:           {
1430: 	    break_points_sorted.push(std::make_pair(t, i));
1431:           }
1432:         else // if t == 0, i is in the active set
1433:           {
1434:             gradient_d.block(m_d_dof)[i] = 0;
1435:             if (gradient_g.block(m_d_dof)[i] > 0)
1436:               m_active_set_phasefield(i) = 1; //lower bound
1437:             else
1438:               m_active_set_phasefield(i) = 2; //upper bound
1439:           }
1440:       }
1441: 
1442:     return break_points_sorted;
1443:   }
1444: 
1445:   template <int dim>
1446:   void PhaseFieldMonolithicSolve<dim>::get_error_residual(Errors &error_residual)
1447:   {
1448:     BlockVector<double> error_res(m_dofs_per_block);
1449: 
1450:     for (unsigned int i = 0; i < m_dof_handler.n_dofs(); ++i)
1451:       if (!m_constraints.is_constrained(i))
1452:         error_res(i) = m_system_rhs(i);
1453: 
1454:     error_residual.m_norm = error_res.l2_norm();
1455:     error_residual.m_u    = error_res.block(m_u_dof).l2_norm();
1456:     error_residual.m_d    = error_res.block(m_d_dof).l2_norm();
1457:   }
1458: 
1459:   template <int dim>
1460:   void PhaseFieldMonolithicSolve<dim>::get_error_residual_LBFGSB(Errors &error_residual,
1461: 								 const BlockVector<double> & solution_delta)
1462:   {
1463:     // We use L_2 norm
1464:     BlockVector<double> error_res(m_dofs_per_block);
1465: 
1466:     // For displacement DOFs, except essential boundary conditions
1467:     // and hanging-node constraints, there are no box constraints
1468:     for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
1469:       {
1470:         if (!m_constraints.is_constrained(i))
1471: 	  error_res.block(m_u_dof)[i] = m_system_rhs.block(m_u_dof)[i];
1472:       }
1473: 
1474:     // For phasefield DOFs, there are points with active box constraints
1475:     // and points with inactive box constraints
1476:     const double upper_limit = 1.0;
1477:     Vector<double> solution_phasefield_total(m_solution.block(m_d_dof));
1478:     solution_phasefield_total += solution_delta.block(m_d_dof);
1479: 
1480:     double trial_solution = 0.0;
1481:     for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
1482:       {
1483: 	// phasefield DOFs can still be constrained due to hanging-nodes
1484:         if (!m_constraints.is_constrained(i + m_dofs_per_block[m_u_dof]))
1485:           {
1486:             trial_solution = solution_phasefield_total(i) - m_system_rhs.block(m_d_dof)[i];
1487: 
1488:             if (trial_solution < m_solution.block(m_d_dof)[i])
1489:               error_res.block(m_d_dof)[i] = m_solution.block(m_d_dof)[i] - solution_phasefield_total(i);
1490:             else if (trial_solution > upper_limit)
1491:               error_res.block(m_d_dof)[i] = upper_limit - solution_phasefield_total(i);
1492:             else
1493:               error_res.block(m_d_dof)[i] = (-m_system_rhs.block(m_d_dof)[i]);
1494:           }
1495:       }
1496: 
1497:     error_residual.m_norm = error_res.l2_norm();
1498:     error_residual.m_u    = error_res.block(m_u_dof).l2_norm();
1499:     error_residual.m_d    = error_res.block(m_d_dof).l2_norm();
1500:   }
1501: 
1502:   template <int dim>
1503:   void PhaseFieldMonolithicSolve<dim>::get_error_update(const BlockVector<double> &newton_update,
1504: 							Errors & error_update)
1505:   {
1506:     BlockVector<double> error_ud(m_dofs_per_block);
1507:     for (unsigned int i = 0; i < m_dof_handler.n_dofs(); ++i)
1508:       if (!m_constraints.is_constrained(i))
1509: 	error_ud(i) = newton_update(i);
1510: 
1511:     error_update.m_norm = error_ud.l2_norm();
1512:     error_update.m_u    = error_ud.block(m_u_dof).l2_norm();
1513:     error_update.m_d    = error_ud.block(m_d_dof).l2_norm();
1514:   }
1515: 
1516:   template <int dim>
1517:   void PhaseFieldMonolithicSolve<dim>::read_material_data(const std::string &data_file,
1518: 				                     const unsigned int total_material_regions)
1519:   {
1520:     std::ifstream myfile (data_file);
1521: 
1522:     double lame_lambda, lame_mu, length_scale, gc, viscosity, residual_k;
1523:     int material_region;
1524:     double poisson_ratio;
1525:     if (myfile.is_open())
1526:       {
1527:         m_logfile << "Reading material data file ..." << std::endl;
1528: 
1529:         while ( myfile >> material_region
1530:                        >> lame_lambda
1531: 		       >> lame_mu
1532: 		       >> length_scale
1533: 		       >> gc
1534: 		       >> viscosity
1535: 		       >> residual_k)
1536:           {
1537:             m_material_data[material_region] = {lame_lambda,
1538:         	                                lame_mu,
1539: 						length_scale,
1540: 						gc,
1541: 						viscosity,
1542:                                                 residual_k};
1543:             poisson_ratio = lame_lambda / (2*(lame_lambda + lame_mu));
1544:             Assert( (poisson_ratio <= 0.5)&(poisson_ratio >=-1.0) , ExcInternalError());
1545: 
1546:             m_logfile << "\tRegion " << material_region << " : " << std::endl;
1547:             m_logfile << "\t\tLame lambda = " << lame_lambda << std::endl;
1548:             m_logfile << "\t\tLame mu = "  << lame_mu << std::endl;
1549:             m_logfile << "\t\tPoisson ratio = "  << poisson_ratio << std::endl;
1550:             m_logfile << "\t\tPhase field length scale (l) = " << length_scale << std::endl;
1551:             m_logfile << "\t\tCritical energy release rate (gc) = "  << gc << std::endl;
1552:             m_logfile << "\t\tViscosity for regularization (eta) = "  << viscosity << std::endl;
1553:             m_logfile << "\t\tResidual_k (k) = "  << residual_k << std::endl;
1554:           }
1555: 
1556:         if (m_material_data.size() != total_material_regions)
1557:           {
1558:             m_logfile << "Material data file has " << m_material_data.size() << " rows. However, "
1559:         	      << "the mesh has " << total_material_regions << " material regions."
1560: 		      << std::endl;
1561:             Assert(m_material_data.size() == total_material_regions,
1562:                        ExcDimensionMismatch(m_material_data.size(), total_material_regions));
1563:           }
1564:         myfile.close();
1565:       }
1566:     else
1567:       {
1568: 	m_logfile << "Material data file : " << data_file << " not exist!" << std::endl;
1569: 	Assert(false, ExcMessage("Failed to read material data file"));
1570:       }
1571:   }
1572: 
1573:   template <int dim>
1574:   void PhaseFieldMonolithicSolve<dim>::read_time_data(const std::string &data_file,
1575: 				                 std::vector<std::array<double, 4>> & time_table)
1576:   {
1577:     std::ifstream myfile (data_file);
1578: 
1579:     double t_0, t_1, delta_t, t_magnitude;
1580: 
1581:     if (myfile.is_open())
1582:       {
1583: 	m_logfile << "Reading time data file ..." << std::endl;
1584: 
1585: 	while ( myfile >> t_0
1586: 		       >> t_1
1587: 		       >> delta_t
1588: 		       >> t_magnitude)
1589: 	  {
1590: 	    Assert( t_0 < t_1,
1591: 		    ExcMessage("For each time pair, "
1592: 			       "the start time should be smaller than the end time"));
1593: 	    time_table.push_back({{t_0, t_1, delta_t, t_magnitude}});
1594: 	  }
1595: 
1596: 	Assert(std::fabs(t_1 - m_parameters.m_end_time) < 1.0e-9,
1597: 	       ExcMessage("End time in time table is inconsistent with input data in parameters.prm"));
1598: 
1599: 	Assert(time_table.size() > 0,
1600: 	       ExcMessage("Time data file is empty."));
1601: 	myfile.close();
1602:       }
1603:     else
1604:       {
1605:         m_logfile << "Time data file : " << data_file << " not exist!" << std::endl;
1606:         Assert(false, ExcMessage("Failed to read time data file"));
1607:       }
1608: 
1609:     for (auto & time_group : time_table)
1610:       {
1611: 	m_logfile << "\t\t"
1612: 	          << time_group[0] << ",\t"
1613: 	          << time_group[1] << ",\t"
1614: 		  << time_group[2] << ",\t"
1615: 		  << time_group[3] << std::endl;
1616:       }
1617:   }
1618: 
1619:   template <int dim>
1620:   void PhaseFieldMonolithicSolve<dim>::setup_qph()
1621:   {
1622:     m_logfile << "\t\tSetting up quadrature point data ("
1623: 	      << m_n_q_points
1624: 	      << " points per cell)" << std::endl;
1625: 
1626:     m_quadrature_point_history.clear();
1627:     for (auto const & cell : m_triangulation.active_cell_iterators())
1628:       {
1629: 	m_quadrature_point_history.initialize(cell, m_n_q_points);
1630:       }
1631: 
1632:     unsigned int material_id;
1633:     double lame_lambda = 0.0;
1634:     double lame_mu = 0.0;
1635:     double length_scale = 0.0;
1636:     double gc = 0.0;
1637:     double viscosity = 0.0;
1638:     double residual_k = 0.0;
1639: 
1640:     for (const auto &cell : m_triangulation.active_cell_iterators())
1641:       {
1642:         material_id = cell->material_id();
1643:         if (m_material_data.find(material_id) != m_material_data.end())
1644:           {
1645:             lame_lambda                = m_material_data[material_id][0];
1646:             lame_mu                    = m_material_data[material_id][1];
1647:             length_scale               = m_material_data[material_id][2];
1648:             gc                         = m_material_data[material_id][3];
1649:             viscosity                  = m_material_data[material_id][4];
1650:             residual_k                 = m_material_data[material_id][5];
1651: 	  }
1652:         else
1653:           {
1654:             m_logfile << "Could not find material data for material id: " << material_id << std::endl;
1655:             AssertThrow(false, ExcMessage("Could not find material data for material id."));
1656:           }
1657: 
1658:         const std::vector<std::shared_ptr<PointHistory<dim>>> lqph =
1659:           m_quadrature_point_history.get_data(cell);
1660:         Assert(lqph.size() == m_n_q_points, ExcInternalError());
1661: 
1662:         for (unsigned int q_point = 0; q_point < m_n_q_points; ++q_point)
1663:           lqph[q_point]->setup_lqp(lame_lambda, lame_mu, length_scale,
1664: 				   gc, viscosity, residual_k,
1665: 				   m_parameters.m_plane_stress);
1666:       }
1667:   }
1668: 
1669:   template <int dim>
1670:   BlockVector<double> PhaseFieldMonolithicSolve<dim>::get_total_solution(
1671:     const BlockVector<double> &solution_delta) const
1672:   {
1673:     BlockVector<double> solution_total(m_solution);
1674:     solution_total += solution_delta;
1675:     return solution_total;
1676:   }
1677: 
1678:   template <int dim>
1679:   void
1680:   PhaseFieldMonolithicSolve<dim>::update_qph_incremental(const BlockVector<double> &solution_delta,
1681: 							 const BlockVector<double> &solution_old)
1682:   {
1683:     m_timer.enter_subsection("Update QPH data");
1684: 
1685:     const BlockVector<double> solution_total(get_total_solution(solution_delta));
1686: 
1687:     const UpdateFlags uf_UQPH(update_values | update_gradients);
1688:     PerTaskData_UQPH  per_task_data_UQPH;
1689:     ScratchData_UQPH  scratch_data_UQPH(m_fe,
1690: 					m_qf_cell,
1691: 					uf_UQPH,
1692: 					solution_total,
1693: 					solution_old,
1694: 					m_time.get_delta_t());
1695: 
1696:     auto worker = [this](const typename DoFHandler<dim>::active_cell_iterator &cell,
1697: 	                 ScratchData_UQPH & scratch,
1698: 	                 PerTaskData_UQPH & data)
1699:       {
1700:         this->update_qph_incremental_one_cell(cell, scratch, data);
1701:       };
1702: 
1703:     auto copier = [this](const PerTaskData_UQPH &data)
1704:       {
1705:         this->copy_local_to_global_UQPH(data);
1706:       };
1707: 
1708:     WorkStream::run(
1709: 	m_dof_handler.begin_active(),
1710: 	m_dof_handler.end(),
1711: 	worker,
1712: 	copier,
1713: 	scratch_data_UQPH,
1714: 	per_task_data_UQPH);
1715: 
1716:     m_timer.leave_subsection();
1717:   }
1718: 
1719:   template <int dim>
1720:   struct PhaseFieldMonolithicSolve<dim>::PerTaskData_UQPH
1721:   {
1722:     void reset()
1723:     {}
1724:   };
1725: 
1726:   template <int dim>
1727:   struct PhaseFieldMonolithicSolve<dim>::ScratchData_UQPH
1728:   {
1729:     const BlockVector<double> & m_solution_UQPH;
1730: 
1731:     std::vector<SymmetricTensor<2, dim>> m_solution_symm_grads_u_cell;
1732:     std::vector<double>         m_solution_values_phasefield_cell;
1733:     std::vector<Tensor<1, dim>> m_solution_grad_phasefield_cell;
1734: 
1735:     FEValues<dim> m_fe_values;
1736: 
1737:     const BlockVector<double>&       m_solution_previous_step;
1738:     std::vector<double>              m_phasefield_previous_step_cell;
1739: 
1740:     const double                     m_delta_time;
1741: 
1742:     ScratchData_UQPH(const FiniteElement<dim> & fe_cell,
1743:                      const QGauss<dim> &        qf_cell,
1744:                      const UpdateFlags          uf_cell,
1745:                      const BlockVector<double> &solution_total,
1746: 		     const BlockVector<double> &solution_old,
1747: 		     const double delta_time)
1748:       : m_solution_UQPH(solution_total)
1749:       , m_solution_symm_grads_u_cell(qf_cell.size())
1750:       , m_solution_values_phasefield_cell(qf_cell.size())
1751:       , m_solution_grad_phasefield_cell(qf_cell.size())
1752:       , m_fe_values(fe_cell, qf_cell, uf_cell)
1753:       , m_solution_previous_step(solution_old)
1754:       , m_phasefield_previous_step_cell(qf_cell.size())
1755:       , m_delta_time(delta_time)
1756:     {}
1757: 
1758:     ScratchData_UQPH(const ScratchData_UQPH &rhs)
1759:       : m_solution_UQPH(rhs.m_solution_UQPH)
1760:       , m_solution_symm_grads_u_cell(rhs.m_solution_symm_grads_u_cell)
1761:       , m_solution_values_phasefield_cell(rhs.m_solution_values_phasefield_cell)
1762:       , m_solution_grad_phasefield_cell(rhs.m_solution_grad_phasefield_cell)
1763:       , m_fe_values(rhs.m_fe_values.get_fe(),
1764:                     rhs.m_fe_values.get_quadrature(),
1765:                     rhs.m_fe_values.get_update_flags())
1766:       , m_solution_previous_step(rhs.m_solution_previous_step)
1767:       , m_phasefield_previous_step_cell(rhs.m_phasefield_previous_step_cell)
1768:       , m_delta_time(rhs.m_delta_time)
1769:     {}
1770: 
1771:     void reset()
1772:     {
1773:       const unsigned int n_q_points = m_solution_symm_grads_u_cell.size();
1774:       for (unsigned int q = 0; q < n_q_points; ++q)
1775:         {
1776:           m_solution_symm_grads_u_cell[q]  = 0.0;
1777:           m_solution_values_phasefield_cell[q] = 0.0;
1778:           m_solution_grad_phasefield_cell[q] = 0.0;
1779:           m_phasefield_previous_step_cell[q] = 0.0;
1780:         }
1781:     }
1782:   };
1783: 
1784:   template <int dim>
1785:   void PhaseFieldMonolithicSolve<dim>::update_qph_incremental_one_cell(
1786:     const typename DoFHandler<dim>::active_cell_iterator &cell,
1787:     ScratchData_UQPH & scratch,
1788:     PerTaskData_UQPH & /*data*/)
1789:   {
1790:     scratch.reset();
1791: 
1792:     scratch.m_fe_values.reinit(cell);
1793: 
1794:     const std::vector<std::shared_ptr<PointHistory<dim>>> lqph =
1795:       m_quadrature_point_history.get_data(cell);
1796:     Assert(lqph.size() == m_n_q_points, ExcInternalError());
1797: 
1798:     const FEValuesExtractors::Vector displacement(0);
1799: 
1800:     scratch.m_fe_values[m_u_fe].get_function_symmetric_gradients(
1801:       scratch.m_solution_UQPH, scratch.m_solution_symm_grads_u_cell);
1802:     scratch.m_fe_values[m_d_fe].get_function_values(
1803:       scratch.m_solution_UQPH, scratch.m_solution_values_phasefield_cell);
1804:     scratch.m_fe_values[m_d_fe].get_function_gradients(
1805:       scratch.m_solution_UQPH, scratch.m_solution_grad_phasefield_cell);
1806: 
1807:     scratch.m_fe_values[m_d_fe].get_function_values(
1808:       scratch.m_solution_previous_step, scratch.m_phasefield_previous_step_cell);
1809: 
1810:     for (const unsigned int q_point :
1811:          scratch.m_fe_values.quadrature_point_indices())
1812:       lqph[q_point]->update_field_values(scratch.m_solution_symm_grads_u_cell[q_point],
1813:                                          scratch.m_solution_values_phasefield_cell[q_point],
1814: 					 scratch.m_solution_grad_phasefield_cell[q_point],
1815: 					 scratch.m_phasefield_previous_step_cell[q_point],
1816: 					 scratch.m_delta_time);
1817:   }
1818: 
1819:   template <int dim>
1820:   struct PhaseFieldMonolithicSolve<dim>::PerTaskData_ASM
1821:   {
1822:     FullMatrix<double>                   m_cell_matrix;
1823:     Vector<double>                       m_cell_rhs;
1824:     std::vector<types::global_dof_index> m_local_dof_indices;
1825: 
1826:     PerTaskData_ASM(const unsigned int dofs_per_cell)
1827:       : m_cell_matrix(dofs_per_cell, dofs_per_cell)
1828:       , m_cell_rhs(dofs_per_cell)
1829:       , m_local_dof_indices(dofs_per_cell)
1830:     {}
1831: 
1832:     void reset()
1833:     {
1834:       m_cell_matrix = 0.0;
1835:       m_cell_rhs    = 0.0;
1836:     }
1837:   };
1838: 
1839:   template <int dim>
1840:   struct PhaseFieldMonolithicSolve<dim>::PerTaskData_ASM_RHS_BFGS
1841:   {
1842:     Vector<double>                       m_cell_rhs;
1843:     std::vector<types::global_dof_index> m_local_dof_indices;
1844: 
1845:     PerTaskData_ASM_RHS_BFGS(const unsigned int dofs_per_cell)
1846:       : m_cell_rhs(dofs_per_cell)
1847:       , m_local_dof_indices(dofs_per_cell)
1848:     {}
1849: 
1850:     void reset()
1851:     {
1852:       m_cell_rhs    = 0.0;
1853:     }
1854:   };
1855: 
1856:   template <int dim>
1857:   struct PhaseFieldMonolithicSolve<dim>::ScratchData_ASM
1858:   {
1859:     FEValues<dim>     m_fe_values;
1860:     FEFaceValues<dim> m_fe_face_values;
1861: 
1862:     std::vector<std::vector<double>>                  m_Nx_phasefield;      // shape function values for phase-field
1863:     std::vector<std::vector<Tensor<1, dim>>>          m_grad_Nx_phasefield; // gradient of shape function values for phase field
1864: 
1865:     std::vector<std::vector<Tensor<1, dim>>>          m_Nx_disp;       // shape function values for displacement
1866:     std::vector<std::vector<Tensor<2, dim>>>          m_grad_Nx_disp;  // gradient of shape function values for displacement
1867:     std::vector<std::vector<SymmetricTensor<2, dim>>> m_symm_grad_Nx_disp;  // symmetric gradient of shape function values for displacement
1868: 
1869:     const BlockVector<double>&       m_solution_previous_step;
1870:     std::vector<double>              m_phasefield_previous_step_cell;
1871: 
1872:     ScratchData_ASM(const FiniteElement<dim> & fe_cell,
1873:                     const QGauss<dim> &        qf_cell,
1874:                     const UpdateFlags          uf_cell,
1875: 		    const QGauss<dim - 1> &    qf_face,
1876: 		    const UpdateFlags          uf_face,
1877: 		    const BlockVector<double>& solution_old)
1878:       : m_fe_values(fe_cell, qf_cell, uf_cell)
1879:       , m_fe_face_values(fe_cell, qf_face, uf_face)
1880:       , m_Nx_phasefield(qf_cell.size(),
1881: 	                std::vector<double>(fe_cell.n_dofs_per_cell()))
1882:       , m_grad_Nx_phasefield(qf_cell.size(),
1883: 		             std::vector<Tensor<1, dim>>(fe_cell.n_dofs_per_cell()))
1884:       , m_Nx_disp(qf_cell.size(),
1885: 		  std::vector<Tensor<1, dim>>(fe_cell.n_dofs_per_cell()))
1886:       , m_grad_Nx_disp(qf_cell.size(),
1887:                        std::vector<Tensor<2, dim>>(fe_cell.n_dofs_per_cell()))
1888:       , m_symm_grad_Nx_disp(qf_cell.size(),
1889:                             std::vector<SymmetricTensor<2, dim>>(fe_cell.n_dofs_per_cell()))
1890:       , m_solution_previous_step(solution_old)
1891:       , m_phasefield_previous_step_cell(qf_cell.size())
1892:     {}
1893: 
1894:     ScratchData_ASM(const ScratchData_ASM &rhs)
1895:       : m_fe_values(rhs.m_fe_values.get_fe(),
1896:                     rhs.m_fe_values.get_quadrature(),
1897:                     rhs.m_fe_values.get_update_flags())
1898:       , m_fe_face_values(rhs.m_fe_face_values.get_fe(),
1899: 	                 rhs.m_fe_face_values.get_quadrature(),
1900: 	                 rhs.m_fe_face_values.get_update_flags())
1901:       , m_Nx_phasefield(rhs.m_Nx_phasefield)
1902:       , m_grad_Nx_phasefield(rhs.m_grad_Nx_phasefield)
1903:       , m_Nx_disp(rhs.m_Nx_disp)
1904:       , m_grad_Nx_disp(rhs.m_grad_Nx_disp)
1905:       , m_symm_grad_Nx_disp(rhs.m_symm_grad_Nx_disp)
1906:       , m_solution_previous_step(rhs.m_solution_previous_step)
1907:       , m_phasefield_previous_step_cell(rhs.m_phasefield_previous_step_cell)
1908:     {}
1909: 
1910:     void reset()
1911:     {
1912:       const unsigned int n_q_points      = m_Nx_phasefield.size();
1913:       const unsigned int n_dofs_per_cell = m_Nx_phasefield[0].size();
1914:       for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
1915:         {
1916:           Assert(m_Nx_phasefield[q_point].size() == n_dofs_per_cell,
1917: 		 ExcInternalError());
1918: 
1919:           Assert(m_grad_Nx_phasefield[q_point].size() == n_dofs_per_cell,
1920:                  ExcInternalError());
1921: 
1922:           Assert(m_Nx_disp[q_point].size() == n_dofs_per_cell,
1923: 		 ExcInternalError());
1924: 
1925:           Assert(m_grad_Nx_disp[q_point].size() == n_dofs_per_cell,
1926:                  ExcInternalError());
1927: 
1928:           Assert(m_symm_grad_Nx_disp[q_point].size() == n_dofs_per_cell,
1929:                  ExcInternalError());
1930: 
1931:           m_phasefield_previous_step_cell[q_point] = 0.0;
1932:           for (unsigned int k = 0; k < n_dofs_per_cell; ++k)
1933:             {
1934:               m_Nx_phasefield[q_point][k]           = 0.0;
1935:               m_grad_Nx_phasefield[q_point][k]      = 0.0;
1936:               m_Nx_disp[q_point][k]                 = 0.0;
1937:               m_grad_Nx_disp[q_point][k]            = 0.0;
1938:               m_symm_grad_Nx_disp[q_point][k]       = 0.0;
1939:             }
1940:         }
1941:     }
1942:   };
1943: 
1944: 
1945:   template <int dim>
1946:   struct PhaseFieldMonolithicSolve<dim>::ScratchData_ASM_RHS_BFGS
1947:   {
1948:     FEValues<dim>     m_fe_values;
1949:     FEFaceValues<dim> m_fe_face_values;
1950: 
1951:     std::vector<std::vector<double>>                  m_Nx_phasefield;      // shape function values for phase-field
1952:     std::vector<std::vector<Tensor<1, dim>>>          m_grad_Nx_phasefield; // gradient of shape function values for phase field
1953: 
1954:     std::vector<std::vector<Tensor<1, dim>>>          m_Nx_disp;       // shape function values for displacement
1955:     std::vector<std::vector<Tensor<2, dim>>>          m_grad_Nx_disp;  // gradient of shape function values for displacement
1956:     std::vector<std::vector<SymmetricTensor<2, dim>>> m_symm_grad_Nx_disp;  // symmetric gradient of shape function values for displacement
1957: 
1958:     const BlockVector<double>&       m_solution_previous_step;
1959:     std::vector<double>              m_phasefield_previous_step_cell;
1960: 
1961:     ScratchData_ASM_RHS_BFGS(const FiniteElement<dim> & fe_cell,
1962:                              const QGauss<dim> &        qf_cell,
1963:                              const UpdateFlags          uf_cell,
1964: 		             const QGauss<dim - 1> &    qf_face,
1965: 		             const UpdateFlags          uf_face,
1966: 		             const BlockVector<double>& solution_old)
1967:       : m_fe_values(fe_cell, qf_cell, uf_cell)
1968:       , m_fe_face_values(fe_cell, qf_face, uf_face)
1969:       , m_Nx_phasefield(qf_cell.size(),
1970: 	                std::vector<double>(fe_cell.n_dofs_per_cell()))
1971:       , m_grad_Nx_phasefield(qf_cell.size(),
1972: 		             std::vector<Tensor<1, dim>>(fe_cell.n_dofs_per_cell()))
1973:       , m_Nx_disp(qf_cell.size(),
1974: 		  std::vector<Tensor<1, dim>>(fe_cell.n_dofs_per_cell()))
1975:       , m_grad_Nx_disp(qf_cell.size(),
1976:                        std::vector<Tensor<2, dim>>(fe_cell.n_dofs_per_cell()))
1977:       , m_symm_grad_Nx_disp(qf_cell.size(),
1978:                             std::vector<SymmetricTensor<2, dim>>(fe_cell.n_dofs_per_cell()))
1979:       , m_solution_previous_step(solution_old)
1980:       , m_phasefield_previous_step_cell(qf_cell.size())
1981:     {}
1982: 
1983:     ScratchData_ASM_RHS_BFGS(const ScratchData_ASM_RHS_BFGS &rhs)
1984:       : m_fe_values(rhs.m_fe_values.get_fe(),
1985:                     rhs.m_fe_values.get_quadrature(),
1986:                     rhs.m_fe_values.get_update_flags())
1987:       , m_fe_face_values(rhs.m_fe_face_values.get_fe(),
1988: 	                 rhs.m_fe_face_values.get_quadrature(),
1989: 	                 rhs.m_fe_face_values.get_update_flags())
1990:       , m_Nx_phasefield(rhs.m_Nx_phasefield)
1991:       , m_grad_Nx_phasefield(rhs.m_grad_Nx_phasefield)
1992:       , m_Nx_disp(rhs.m_Nx_disp)
1993:       , m_grad_Nx_disp(rhs.m_grad_Nx_disp)
1994:       , m_symm_grad_Nx_disp(rhs.m_symm_grad_Nx_disp)
1995:       , m_solution_previous_step(rhs.m_solution_previous_step)
1996:       , m_phasefield_previous_step_cell(rhs.m_phasefield_previous_step_cell)
1997:     {}
1998: 
1999:     void reset()
2000:     {
2001:       const unsigned int n_q_points      = m_Nx_phasefield.size();
2002:       const unsigned int n_dofs_per_cell = m_Nx_phasefield[0].size();
2003:       for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
2004:         {
2005:           Assert(m_Nx_phasefield[q_point].size() == n_dofs_per_cell,
2006: 		 ExcInternalError());
2007: 
2008:           Assert(m_grad_Nx_phasefield[q_point].size() == n_dofs_per_cell,
2009:                  ExcInternalError());
2010: 
2011:           Assert(m_Nx_disp[q_point].size() == n_dofs_per_cell,
2012: 		 ExcInternalError());
2013: 
2014:           Assert(m_grad_Nx_disp[q_point].size() == n_dofs_per_cell,
2015:                  ExcInternalError());
2016: 
2017:           Assert(m_symm_grad_Nx_disp[q_point].size() == n_dofs_per_cell,
2018:                  ExcInternalError());
2019: 
2020:           m_phasefield_previous_step_cell[q_point] = 0.0;
2021:           for (unsigned int k = 0; k < n_dofs_per_cell; ++k)
2022:             {
2023:               m_Nx_phasefield[q_point][k]           = 0.0;
2024:               m_grad_Nx_phasefield[q_point][k]      = 0.0;
2025:               m_Nx_disp[q_point][k]                 = 0.0;
2026:               m_grad_Nx_disp[q_point][k]            = 0.0;
2027:               m_symm_grad_Nx_disp[q_point][k]       = 0.0;
2028:             }
2029:         }
2030:     }
2031:   };
2032: 
2033:   // constructor has no return type
2034:   template <int dim>
2035:   PhaseFieldMonolithicSolve<dim>::PhaseFieldMonolithicSolve(const std::string &input_file)
2036:     : m_parameters(input_file)
2037:     , m_triangulation(Triangulation<dim>::maximum_smoothing)
2038:     , m_time(m_parameters.m_end_time)
2039:     , m_logfile(m_parameters.m_logfile_name)
2040:     , m_timer(m_logfile, TimerOutput::summary, TimerOutput::wall_times)
2041:     , m_dof_handler(m_triangulation)
2042:     , m_fe(FE_Q<dim>(m_parameters.m_poly_degree),
2043: 	   dim, // displacement
2044: 	   FE_Q<dim>(m_parameters.m_poly_degree),
2045: 	   1)   // phasefield
2046:     , m_dofs_per_cell(m_fe.n_dofs_per_cell())
2047:     , m_u_fe(m_first_u_component)
2048:     , m_d_fe(m_d_component)
2049:     , m_dofs_per_block(m_n_blocks)
2050:     , m_qf_cell(m_parameters.m_quad_order)
2051:     , m_qf_face(m_parameters.m_quad_order)
2052:     , m_n_q_points(m_qf_cell.size())
2053:     , m_vol_reference(0.0)
2054:   {}
2055: 
2056:   template <int dim>
2057:   void PhaseFieldMonolithicSolve<dim>::make_grid()
2058:   {
2059:     if (m_parameters.m_scenario == 1)
2060:       make_grid_case_1();
2061:     else if (m_parameters.m_scenario == 2)
2062:       make_grid_case_2();
2063:     else if (m_parameters.m_scenario == 3)
2064:       make_grid_case_3();
2065:     else if (m_parameters.m_scenario == 4)
2066:       make_grid_case_4();
2067:     else if (m_parameters.m_scenario == 5)
2068:       make_grid_case_5();
2069:     else if (m_parameters.m_scenario == 6)
2070:       make_grid_case_6();
2071:     else if (m_parameters.m_scenario == 7)
2072:       make_grid_case_7();
2073:     else if (m_parameters.m_scenario == 8)
2074:       make_grid_case_8();
2075:     else if (m_parameters.m_scenario == 9)
2076:       make_grid_case_9();
2077:     else if (m_parameters.m_scenario == 10)
2078:       make_grid_case_10();
2079:     else if (m_parameters.m_scenario == 11)
2080:       make_grid_case_11();
2081:     else
2082:       Assert(false, ExcMessage("The scenario has not been implemented!"));
2083: 
2084:     m_logfile << "\t\tTriangulation:"
2085:               << "\n\t\t\tNumber of active cells: "
2086:               << m_triangulation.n_active_cells()
2087:               << "\n\t\t\tNumber of used vertices: "
2088:               << m_triangulation.n_used_vertices()
2089: 	      << std::endl;
2090: 
2091:     std::ofstream out("original_mesh.vtu");
2092:     GridOut       grid_out;
2093:     grid_out.write_vtu(m_triangulation, out);
2094: 
2095:     m_vol_reference = GridTools::volume(m_triangulation);
2096:     m_logfile << "\t\tGrid:\n\t\t\tReference volume: " << m_vol_reference << std::endl;
2097:   }
2098: 
2099:   template <int dim>
2100:   void PhaseFieldMonolithicSolve<dim>::make_grid_case_1()
2101:   {
2102:     for (unsigned int i = 0; i < 80; ++i)
2103:       m_logfile << "*";
2104:     m_logfile << std::endl;
2105:     m_logfile << "\t\t\tSquare tension (unstructured)" << std::endl;
2106:     for (unsigned int i = 0; i < 80; ++i)
2107:       m_logfile << "*";
2108:     m_logfile << std::endl;
2109: 
2110:     AssertThrow(dim==2, ExcMessage("The dimension has to be 2D!"));
2111: 
2112:     GridIn<dim> gridin;
2113:     gridin.attach_triangulation(m_triangulation);
2114:     std::ifstream f("square_tension_unstructured.msh");
2115:     gridin.read_msh(f);
2116: 
2117:     for (const auto &cell : m_triangulation.active_cell_iterators())
2118:       for (const auto &face : cell->face_iterators())
2119: 	{
2120: 	  if (face->at_boundary() == true)
2121: 	    {
2122: 	      if (std::fabs(face->center()[1] + 0.5 ) < 1.0e-9 )
2123: 		face->set_boundary_id(0);
2124: 	      else if (std::fabs(face->center()[1] - 0.5 ) < 1.0e-9)
2125: 	        face->set_boundary_id(1);
2126: 	      else
2127: 	        face->set_boundary_id(2);
2128: 	    }
2129: 	}
2130: 
2131:     m_triangulation.refine_global(m_parameters.m_global_refine_times);
2132: 
2133:     if (m_parameters.m_refinement_strategy == "pre-refine")
2134:       {
2135: 	unsigned int material_id;
2136: 	double length_scale;
2137: 	for (unsigned int i = 0; i < m_parameters.m_local_prerefine_times; i++)
2138: 	  {
2139: 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2140: 	      {
2141: 		if (   std::fabs(cell->center()[1]) < 0.01
2142: 		    && cell->center()[0] > 0.495)
2143: 		  {
2144: 		    material_id = cell->material_id();
2145: 		    length_scale = m_material_data[material_id][2];
2146: 		    if (  std::sqrt(cell->measure())
2147: 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2148: 		      cell->set_refine_flag();
2149: 		  }
2150: 	      }
2151: 	    m_triangulation.execute_coarsening_and_refinement();
2152: 	  }
2153:       }
2154:     else if (m_parameters.m_refinement_strategy == "adaptive-refine")
2155:       {
2156: 	unsigned int material_id;
2157: 	double length_scale;
2158: 	bool initiation_point_refine_unfinished = true;
2159: 	while (initiation_point_refine_unfinished)
2160: 	  {
2161: 	    initiation_point_refine_unfinished = false;
2162: 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2163: 	      {
2164: 		if (   std::fabs(cell->center()[1] - 0.0) < 0.05
2165: 		    && std::fabs(cell->center()[0] - 0.5) < 0.05)
2166: 		  {
2167: 		    material_id = cell->material_id();
2168: 		    length_scale = m_material_data[material_id][2];
2169: 		    if (  std::sqrt(cell->measure())
2170: 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2171: 		      {
2172: 		        cell->set_refine_flag();
2173: 		        initiation_point_refine_unfinished = true;
2174: 		      }
2175: 		  }
2176: 	      }
2177: 	    m_triangulation.execute_coarsening_and_refinement();
2178: 	  }
2179:       }
2180:     else
2181:       {
2182: 	AssertThrow(false,
2183: 	            ExcMessage("Selected mesh refinement strategy not implemented!"));
2184:       }
2185:   }
2186: 
2187: 
2188:   template <int dim>
2189:   void PhaseFieldMonolithicSolve<dim>::make_grid_case_2()
2190:   {
2191:     for (unsigned int i = 0; i < 80; ++i)
2192:       m_logfile << "*";
2193:     m_logfile << std::endl;
2194:     m_logfile << "\t\t\t\tSquare shear (unstructured)" << std::endl;
2195:     for (unsigned int i = 0; i < 80; ++i)
2196:       m_logfile << "*";
2197:     m_logfile << std::endl;
2198: 
2199:     AssertThrow(dim==2, ExcMessage("The dimension has to be 2D!"));
2200: 
2201:     GridIn<dim> gridin;
2202:     gridin.attach_triangulation(m_triangulation);
2203:     std::ifstream f("square_shear_unstructured.msh");
2204:     gridin.read_msh(f);
2205: 
2206:     for (const auto &cell : m_triangulation.active_cell_iterators())
2207:       for (const auto &face : cell->face_iterators())
2208: 	{
2209: 	  if (face->at_boundary() == true)
2210: 	    {
2211: 	      if (std::fabs(face->center()[1] + 0.5 ) < 1.0e-9 )
2212: 		face->set_boundary_id(0);
2213: 	      else if (std::fabs(face->center()[1] - 0.5 ) < 1.0e-9)
2214: 	        face->set_boundary_id(1);
2215: 	      else if (   (std::fabs(face->center()[0] - 0.0 ) < 1.0e-9)
2216: 		       || (std::fabs(face->center()[0] - 1.0 ) < 1.0e-9))
2217: 	        face->set_boundary_id(2);
2218: 	      else
2219: 	        face->set_boundary_id(3);
2220: 	    }
2221: 	}
2222: 
2223:     m_triangulation.refine_global(m_parameters.m_global_refine_times);
2224: 
2225:     if (m_parameters.m_refinement_strategy == "pre-refine")
2226:       {
2227: 	unsigned int material_id;
2228: 	double length_scale;
2229: 	for (unsigned int i = 0; i < m_parameters.m_local_prerefine_times; i++)
2230: 	  {
2231: 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2232: 	      {
2233: 		if (    (cell->center()[0] > 0.45)
2234: 		     && (cell->center()[1] < 0.05) )
2235: 		  {
2236: 		    material_id = cell->material_id();
2237: 		    length_scale = m_material_data[material_id][2];
2238: 		    if (  std::sqrt(cell->measure())
2239: 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2240: 		      cell->set_refine_flag();
2241: 		  }
2242: 	      }
2243: 	    m_triangulation.execute_coarsening_and_refinement();
2244: 	  }
2245:       }
2246:     else if (m_parameters.m_refinement_strategy == "adaptive-refine")
2247:       {
2248: 	unsigned int material_id;
2249: 	double length_scale;
2250: 	bool initiation_point_refine_unfinished = true;
2251: 	while (initiation_point_refine_unfinished)
2252: 	  {
2253: 	    initiation_point_refine_unfinished = false;
2254: 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2255: 	      {
2256: 		if (    std::fabs(cell->center()[0] - 0.5) < 0.025
2257: 		     && cell->center()[1] < 0.0 && cell->center()[1] > -0.025)
2258: 		  {
2259: 		    material_id = cell->material_id();
2260: 		    length_scale = m_material_data[material_id][2];
2261: 		    if (  std::sqrt(cell->measure())
2262: 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2263: 		      {
2264: 		        cell->set_refine_flag();
2265: 		        initiation_point_refine_unfinished = true;
2266: 		      }
2267: 		  }
2268: 	      }
2269: 	    m_triangulation.execute_coarsening_and_refinement();
2270: 	  }
2271:       }
2272:     else
2273:       {
2274: 	AssertThrow(false,
2275: 	            ExcMessage("Selected mesh refinement strategy not implemented!"));
2276:       }
2277:   }
2278: 
2279:   template <int dim>
2280:   void PhaseFieldMonolithicSolve<dim>::make_grid_case_3()
2281:   {
2282:     for (unsigned int i = 0; i < 80; ++i)
2283:       m_logfile << "*";
2284:     m_logfile << std::endl;
2285:     m_logfile << "\t\t\tSquare tension (structured)" << std::endl;
2286:     for (unsigned int i = 0; i < 80; ++i)
2287:       m_logfile << "*";
2288:     m_logfile << std::endl;
2289: 
2290:     AssertThrow(dim==2, ExcMessage("The dimension has to be 2D!"));
2291: 
2292:     GridIn<dim> gridin;
2293:     gridin.attach_triangulation(m_triangulation);
2294:     std::ifstream f("square_tension_structured.msh");
2295:     gridin.read_msh(f);
2296: 
2297:     for (const auto &cell : m_triangulation.active_cell_iterators())
2298:       for (const auto &face : cell->face_iterators())
2299: 	{
2300: 	  if (face->at_boundary() == true)
2301: 	    {
2302: 	      if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9 )
2303: 		face->set_boundary_id(0);
2304: 	      else if (std::fabs(face->center()[1] - 1.0 ) < 1.0e-9)
2305: 	        face->set_boundary_id(1);
2306: 	      else
2307: 	        face->set_boundary_id(2);
2308: 	    }
2309: 	}
2310: 
2311:     m_triangulation.refine_global(m_parameters.m_global_refine_times);
2312: 
2313:     if (m_parameters.m_refinement_strategy == "pre-refine")
2314:       {
2315: 	unsigned int material_id;
2316: 	double length_scale;
2317: 	for (unsigned int i = 0; i < m_parameters.m_local_prerefine_times; i++)
2318: 	  {
2319: 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2320: 	      {
2321: 		if (    (std::fabs(cell->center()[1] - 0.5) < 0.025)
2322: 		     && (cell->center()[0] > 0.475) )
2323: 		  {
2324: 		    material_id = cell->material_id();
2325: 		    length_scale = m_material_data[material_id][2];
2326: 		    if (  std::sqrt(cell->measure())
2327: 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2328: 		      cell->set_refine_flag();
2329: 		  }
2330: 	      }
2331: 	    m_triangulation.execute_coarsening_and_refinement();
2332: 	  }
2333:       }
2334:     else if (m_parameters.m_refinement_strategy == "adaptive-refine")
2335:       {
2336: 	unsigned int material_id;
2337: 	double length_scale;
2338: 	bool initiation_point_refine_unfinished = true;
2339: 	while (initiation_point_refine_unfinished)
2340: 	  {
2341: 	    initiation_point_refine_unfinished = false;
2342: 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2343: 	      {
2344: 		if (    std::fabs(cell->center()[0] - 0.5) < 0.025
2345: 		     && std::fabs(cell->center()[1] - 0.5) < 0.025 )
2346: 		  {
2347: 		    material_id = cell->material_id();
2348: 		    length_scale = m_material_data[material_id][2];
2349: 		    if (  std::sqrt(cell->measure())
2350: 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2351: 		      {
2352: 		        cell->set_refine_flag();
2353: 		        initiation_point_refine_unfinished = true;
2354: 		      }
2355: 		  }
2356: 	      }
2357: 	    m_triangulation.execute_coarsening_and_refinement();
2358: 	  }
2359:       }
2360:     else
2361:       {
2362: 	AssertThrow(false,
2363: 	            ExcMessage("Selected mesh refinement strategy not implemented!"));
2364:       }
2365:   }
2366: 
2367:   template <int dim>
2368:   void PhaseFieldMonolithicSolve<dim>::make_grid_case_4()
2369:   {
2370:     for (unsigned int i = 0; i < 80; ++i)
2371:       m_logfile << "*";
2372:     m_logfile << std::endl;
2373:     m_logfile << "\t\t\t\tSquare shear (structured)" << std::endl;
2374:     for (unsigned int i = 0; i < 80; ++i)
2375:       m_logfile << "*";
2376:     m_logfile << std::endl;
2377: 
2378:     AssertThrow(dim==2, ExcMessage("The dimension has to be 2D!"));
2379: 
2380:     GridIn<dim> gridin;
2381:     gridin.attach_triangulation(m_triangulation);
2382:     std::ifstream f("square_shear_structured.msh");
2383:     gridin.read_msh(f);
2384: 
2385:     for (const auto &cell : m_triangulation.active_cell_iterators())
2386:       for (const auto &face : cell->face_iterators())
2387: 	{
2388: 	  if (face->at_boundary() == true)
2389: 	    {
2390: 	      if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9 )
2391: 		face->set_boundary_id(0);
2392: 	      else if (std::fabs(face->center()[1] - 1.0 ) < 1.0e-9)
2393: 	        face->set_boundary_id(1);
2394: 	      else if (   (std::fabs(face->center()[0] - 0.0 ) < 1.0e-9)
2395: 		       || (std::fabs(face->center()[0] - 1.0 ) < 1.0e-9))
2396: 	        face->set_boundary_id(2);
2397: 	      else
2398: 	        face->set_boundary_id(3);
2399: 	    }
2400: 	}
2401: 
2402:     m_triangulation.refine_global(m_parameters.m_global_refine_times);
2403: 
2404:     if (m_parameters.m_refinement_strategy == "pre-refine")
2405:       {
2406: 	unsigned int material_id;
2407: 	double length_scale;
2408: 	for (unsigned int i = 0; i < m_parameters.m_local_prerefine_times; i++)
2409: 	  {
2410: 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2411: 	      {
2412: 		if (   (   (cell->center()[0] > 0.475)
2413: 		        && (cell->center()[1] < 0.525) )
2414: 		    //|| (    cell->center()[1] > 0.975)
2415: 		   )
2416: 		  {
2417: 		    material_id = cell->material_id();
2418: 		    length_scale = m_material_data[material_id][2];
2419: 		    if (  std::sqrt(cell->measure())
2420: 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2421: 		      cell->set_refine_flag();
2422: 		  }
2423: 	      }
2424: 	    m_triangulation.execute_coarsening_and_refinement();
2425: 	  }
2426:       }
2427:     else if (m_parameters.m_refinement_strategy == "adaptive-refine")
2428:       {
2429: 	unsigned int material_id;
2430: 	double length_scale;
2431: 	bool initiation_point_refine_unfinished = true;
2432: 	while (initiation_point_refine_unfinished)
2433: 	  {
2434: 	    initiation_point_refine_unfinished = false;
2435: 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2436: 	      {
2437: 		if (   (   std::fabs(cell->center()[0] - 0.5) < 0.025
2438: 		        && cell->center()[1] < 0.5 && cell->center()[1] > 0.475 )
2439: 		    // we also need to refine the top edge, since there might be a conflict between
2440: 		    // inhomogeneous boundary conditions and the hanging-node constraints at the
2441: 		    // top edge
2442: 		    || (   cell->center()[1] > 0.975 ) )
2443: 		  {
2444: 		    material_id = cell->material_id();
2445: 		    length_scale = m_material_data[material_id][2];
2446: 		    if (  std::sqrt(cell->measure())
2447: 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2448: 		      {
2449: 		        cell->set_refine_flag();
2450: 		        initiation_point_refine_unfinished = true;
2451: 		      }
2452: 		  }
2453: 	      }
2454: 	    m_triangulation.execute_coarsening_and_refinement();
2455: 	  }
2456:       }
2457:     else
2458:       {
2459: 	AssertThrow(false,
2460: 	            ExcMessage("Selected mesh refinement strategy not implemented!"));
2461:       }
2462:   }
2463: 
2464:   template <int dim>
2465:   void PhaseFieldMonolithicSolve<dim>::make_grid_case_5()
2466:   {
2467:     for (unsigned int i = 0; i < 80; ++i)
2468:       m_logfile << "*";
2469:     m_logfile << std::endl;
2470:     m_logfile << "\t\t\t\tThree-point bending (structured)" << std::endl;
2471:     for (unsigned int i = 0; i < 80; ++i)
2472:       m_logfile << "*";
2473:     m_logfile << std::endl;
2474: 
2475:     AssertThrow(dim==2, ExcMessage("The dimension has to be 2D!"));
2476: 
2477:     GridIn<dim> gridin;
2478:     gridin.attach_triangulation(m_triangulation);
2479:     std::ifstream f("three_point_bending_structured.msh");
2480:     gridin.read_msh(f);
2481: 
2482:     for (const auto &cell : m_triangulation.active_cell_iterators())
2483:       for (const auto &face : cell->face_iterators())
2484: 	{
2485: 	  if (face->at_boundary() == true)
2486: 	    {
2487: 	      if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9 )
2488: 		face->set_boundary_id(0);
2489: 	      else if (std::fabs(face->center()[1] - 2.0 ) < 1.0e-9)
2490: 	        face->set_boundary_id(1);
2491: 	      else
2492: 	        face->set_boundary_id(2);
2493: 	    }
2494: 	}
2495: 
2496:     m_triangulation.refine_global(m_parameters.m_global_refine_times);
2497: 
2498:     if (m_parameters.m_refinement_strategy == "pre-refine")
2499:       {
2500: 	for (const auto &cell : m_triangulation.active_cell_iterators())
2501: 	  {
2502: 	    if (    std::fabs(cell->center()[0] - 4.0) < 0.075
2503: 		 && cell->center()[1] < 1.6)
2504: 	      {
2505: 		cell->set_refine_flag();
2506: 	      }
2507: 	  }
2508: 	m_triangulation.execute_coarsening_and_refinement();
2509: 
2510: 	unsigned int material_id;
2511: 	double length_scale;
2512: 	for (unsigned int i = 0; i < m_parameters.m_local_prerefine_times; i++)
2513: 	  {
2514: 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2515: 	      {
2516: 		if (    std::fabs(cell->center()[0] - 4.0) < 0.05
2517: 		     && cell->center()[1] < 1.6)
2518: 		  {
2519: 		    material_id = cell->material_id();
2520: 		    length_scale = m_material_data[material_id][2];
2521: 		    if (  std::sqrt(cell->measure())
2522: 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2523: 		      cell->set_refine_flag();
2524: 		  }
2525: 	      }
2526: 	    m_triangulation.execute_coarsening_and_refinement();
2527: 	  }
2528:       }
2529:     else if (m_parameters.m_refinement_strategy == "adaptive-refine")
2530:       {
2531: 	unsigned int material_id;
2532: 	double length_scale;
2533: 	bool initiation_point_refine_unfinished = true;
2534: 	while (initiation_point_refine_unfinished)
2535: 	  {
2536: 	    initiation_point_refine_unfinished = false;
2537: 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2538: 	      {
2539: 		if (    std::fabs(cell->center()[0] - 4.0) < 0.075
2540: 		     && std::fabs(cell->center()[1] - 0.4) < 0.075 )
2541: 		  {
2542: 		    material_id = cell->material_id();
2543: 		    length_scale = m_material_data[material_id][2];
2544: 		    if (  std::sqrt(cell->measure())
2545: 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2546: 		      {
2547: 		        cell->set_refine_flag();
2548: 		        initiation_point_refine_unfinished = true;
2549: 		      }
2550: 		  }
2551: 	      }
2552: 	    m_triangulation.execute_coarsening_and_refinement();
2553: 	  }
2554:       }
2555:     else
2556:       {
2557: 	AssertThrow(false,
2558: 	            ExcMessage("Selected mesh refinement strategy not implemented!"));
2559:       }
2560:   }
2561: 
2562:   template <int dim>
2563:   void PhaseFieldMonolithicSolve<dim>::make_grid_case_6()
2564:   {
2565:     AssertThrow(dim==3, ExcMessage("The dimension has to be 3D!"));
2566: 
2567:     for (unsigned int i = 0; i < 80; ++i)
2568:       m_logfile << "*";
2569:     m_logfile << std::endl;
2570:     m_logfile << "\t\t\t\tSphere inclusion (3D structured)" << std::endl;
2571:     for (unsigned int i = 0; i < 80; ++i)
2572:       m_logfile << "*";
2573:     m_logfile << std::endl;
2574: 
2575:     Triangulation<dim> tria_inner;
2576:     GridGenerator::hyper_ball(tria_inner, Point<dim>(), 0.5);
2577: 
2578:     Triangulation<dim> tria_outer;
2579:     GridGenerator::hyper_shell(
2580:       tria_outer, Point<dim>(), 0.5, std::sqrt(dim), 2 * dim);
2581: 
2582:     Triangulation<dim> tmp_triangulation;
2583: 
2584:     GridGenerator::merge_triangulations(tria_inner, tria_outer, tmp_triangulation);
2585: 
2586:     tmp_triangulation.reset_all_manifolds();
2587:     tmp_triangulation.set_all_manifold_ids(0);
2588: 
2589:     for (const auto &cell : tmp_triangulation.cell_iterators())
2590:       {
2591:         for (const auto &face : cell->face_iterators())
2592:           {
2593:             bool face_at_sphere_boundary = true;
2594:             for (const auto v : face->vertex_indices())
2595:               {
2596:                 if (std::abs(face->vertex(v).norm_square() - 0.25) > 1e-12)
2597:                   face_at_sphere_boundary = false;
2598:               }
2599:             if (face_at_sphere_boundary)
2600:               face->set_all_manifold_ids(1);
2601:           }
2602:         if (cell->center().norm_square() < 0.25)
2603:           cell->set_material_id(1);
2604:         else
2605:           cell->set_material_id(0);
2606:       }
2607: 
2608:     tmp_triangulation.set_manifold(1, SphericalManifold<dim>());
2609: 
2610:     TransfiniteInterpolationManifold<dim> transfinite_manifold;
2611:     transfinite_manifold.initialize(tmp_triangulation);
2612:     tmp_triangulation.set_manifold(0, transfinite_manifold);
2613: 
2614:     tmp_triangulation.refine_global(m_parameters.m_global_refine_times);
2615: 
2616:     std::set<typename Triangulation< dim >::active_cell_iterator >
2617:       cells_to_remove;
2618: 
2619:     for (const auto &cell : tmp_triangulation.active_cell_iterators())
2620:       {
2621: 	if (   cell->center()[0] < 0.0
2622: 	    || cell->center()[1] < 0.0
2623: 	    || cell->center()[2] < 0.0)
2624: 	  {
2625: 	    cells_to_remove.insert(cell);
2626: 	  }
2627:       }
2628: 
2629:     GridGenerator::create_triangulation_with_removed_cells(tmp_triangulation,
2630: 							   cells_to_remove,
2631: 							   m_triangulation);
2632: 
2633:     for (const auto &cell : m_triangulation.active_cell_iterators())
2634:       for (const auto &face : cell->face_iterators())
2635: 	{
2636: 	  if (face->at_boundary() == true)
2637: 	    {
2638: 	      if (std::fabs(face->center()[0] - 0.0 ) < 1.0e-9 )
2639: 		face->set_boundary_id(0);
2640: 	      else if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9)
2641: 		face->set_boundary_id(1);
2642: 	      else if (std::fabs(face->center()[2] - 0.0 ) < 1.0e-9)
2643: 		face->set_boundary_id(2);
2644: 	      else if (std::fabs(face->center()[2] - 1.0 ) < 1.0e-9)
2645: 		face->set_boundary_id(3);
2646: 	      else
2647: 		face->set_boundary_id(4);
2648: 	    }
2649: 	}
2650: 
2651:     if (m_parameters.m_refinement_strategy == "adaptive-refine")
2652:       {
2653: 	unsigned int material_id;
2654: 	double length_scale;
2655: 	bool initiation_point_refine_unfinished = true;
2656: 	while (initiation_point_refine_unfinished)
2657: 	  {
2658: 	    initiation_point_refine_unfinished = false;
2659: 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2660: 	      {
2661: 		if (    cell->center()[2] > 0.525
2662: 		     && cell->center()[2] < 0.575
2663: 		     && cell->center()[0] < 0.05
2664: 		     && cell->center()[1] < 0.05 )
2665: 		  {
2666: 		    material_id = cell->material_id();
2667: 		    length_scale = m_material_data[material_id][2];
2668: 		    if (  std::cbrt(cell->measure())
2669: 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2670: 		      {
2671: 			cell->set_refine_flag();
2672: 			initiation_point_refine_unfinished = true;
2673: 		      }
2674: 		  }
2675: 	      }
2676: 	    m_triangulation.execute_coarsening_and_refinement();
2677: 	  }
2678:       }
2679:     else
2680:       {
2681: 	AssertThrow(false,
2682: 		    ExcMessage("Selected mesh refinement strategy not implemented!"));
2683:       }
2684:   }
2685: 
2686:   template <int dim>
2687:   void PhaseFieldMonolithicSolve<dim>::make_grid_case_7()
2688:   {
2689:     AssertThrow(dim==3, ExcMessage("The dimension has to be 3D!"));
2690: 
2691:     for (unsigned int i = 0; i < 80; ++i)
2692:       m_logfile << "*";
2693:     m_logfile << std::endl;
2694:     m_logfile << "\t\t\t\tSphere inclusion (3D structured version 2)" << std::endl;
2695:     for (unsigned int i = 0; i < 80; ++i)
2696:       m_logfile << "*";
2697:     m_logfile << std::endl;
2698: 
2699:     Triangulation<dim> tria_inner;
2700:     GridGenerator::hyper_ball(tria_inner, Point<dim>(), 0.49);
2701: 
2702:     Triangulation<dim> tria_outer;
2703:     GridGenerator::hyper_shell(
2704:       tria_outer, Point<dim>(), 0.49, std::sqrt(dim)*0.5, 2 * dim);
2705: 
2706:     Triangulation<dim> cube1;
2707:     GridGenerator::hyper_rectangle(cube1, Point<dim>(0, 0, 0.5), Point<dim>(1, 1, 1.5));
2708:     Triangulation<dim> cube2;
2709:     GridGenerator::hyper_rectangle(cube2, Point<dim>(0, 0.5, -0.5), Point<dim>(1, 1.5, 0.5));
2710:     Triangulation<dim> cube3;
2711:     GridGenerator::hyper_rectangle(cube3, Point<dim>(0.5, -0.5, -0.5), Point<dim>(1.5, 0.5, 0.5));
2712: 
2713:     Triangulation<dim> tmp_triangulation;
2714:     GridGenerator::merge_triangulations({&tria_inner, &tria_outer,
2715:                                          &cube1, &cube2, &cube3}, tmp_triangulation);
2716: 
2717:     tmp_triangulation.reset_all_manifolds();
2718:     tmp_triangulation.set_all_manifold_ids(0);
2719: 
2720:     for (const auto &cell : tmp_triangulation.cell_iterators())
2721:       {
2722:         for (const auto &face : cell->face_iterators())
2723:           {
2724:             bool face_at_sphere_boundary = true;
2725:             for (const auto v : face->vertex_indices())
2726:               {
2727:                 if (std::abs(face->vertex(v).norm_square() - 0.49 * 0.49) > 1e-12)
2728:                   face_at_sphere_boundary = false;
2729:               }
2730:             if (face_at_sphere_boundary)
2731:               face->set_all_manifold_ids(1);
2732:           }
2733:         if (cell->center().norm_square() < 0.1)
2734:           cell->set_material_id(1);
2735:         else
2736:           cell->set_material_id(0);
2737:       }
2738: 
2739:     tmp_triangulation.set_manifold(1, SphericalManifold<dim>());
2740: 
2741:     TransfiniteInterpolationManifold<dim> transfinite_manifold;
2742:     transfinite_manifold.initialize(tmp_triangulation);
2743:     tmp_triangulation.set_manifold(0, transfinite_manifold);
2744: 
2745:     tmp_triangulation.refine_global(m_parameters.m_global_refine_times);
2746: 
2747:     std::set<typename Triangulation< dim >::active_cell_iterator >
2748:       cells_to_remove;
2749: 
2750:     for (const auto &cell : tmp_triangulation.active_cell_iterators())
2751:       {
2752: 	if (   cell->center()[0] < 0.0
2753: 	    || cell->center()[1] < 0.0
2754: 	    || cell->center()[2] < 0.0
2755: 	    || cell->center()[0] > 1.0
2756: 	    || cell->center()[1] > 1.0
2757: 	    || cell->center()[2] > 1.0)
2758: 	  {
2759: 	    cells_to_remove.insert(cell);
2760: 	  }
2761:       }
2762: 
2763:     GridGenerator::create_triangulation_with_removed_cells(tmp_triangulation,
2764: 							   cells_to_remove,
2765: 							   m_triangulation);
2766: 
2767:     for (const auto &cell : m_triangulation.active_cell_iterators())
2768:       for (const auto &face : cell->face_iterators())
2769: 	{
2770: 	  if (face->at_boundary() == true)
2771: 	    {
2772: 	      if (std::fabs(face->center()[0] - 0.0 ) < 1.0e-9 )
2773: 		face->set_boundary_id(0);
2774: 	      else if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9)
2775: 		face->set_boundary_id(1);
2776: 	      else if (std::fabs(face->center()[2] - 0.0 ) < 1.0e-9)
2777: 		face->set_boundary_id(2);
2778: 	      else if (std::fabs(face->center()[2] - 1.0 ) < 1.0e-9)
2779: 		face->set_boundary_id(3);
2780: 	      else
2781: 		face->set_boundary_id(4);
2782: 	    }
2783: 	}
2784: 
2785:     if (m_parameters.m_refinement_strategy == "adaptive-refine")
2786:       {
2787: 	unsigned int material_id;
2788: 	double length_scale;
2789: 	bool initiation_point_refine_unfinished = true;
2790: 	while (initiation_point_refine_unfinished)
2791: 	  {
2792: 	    initiation_point_refine_unfinished = false;
2793: 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2794: 	      {
2795: 		if (    cell->center()[2] > 0.505
2796: 		     && cell->center()[2] < 0.575
2797: 		     && cell->center()[0] < 0.05
2798: 		     && cell->center()[1] < 0.05 )
2799: 		  {
2800: 		    material_id = cell->material_id();
2801: 		    length_scale = m_material_data[material_id][2];
2802: 		    if (  std::cbrt(cell->measure())
2803: 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2804: 		      {
2805: 			cell->set_refine_flag();
2806: 			initiation_point_refine_unfinished = true;
2807: 		      }
2808: 		  }
2809: 	      }
2810: 	    m_triangulation.execute_coarsening_and_refinement();
2811: 	  }
2812:       }
2813:     else
2814:       {
2815: 	AssertThrow(false,
2816: 		    ExcMessage("Selected mesh refinement strategy not implemented!"));
2817:       }
2818:   }
2819: 
2820: 
2821:   template <int dim>
2822:   void PhaseFieldMonolithicSolve<dim>::make_grid_case_8()
2823:   {
2824:     AssertThrow(dim==3, ExcMessage("The dimension has to be 3D!"));
2825: 
2826:     for (unsigned int i = 0; i < 80; ++i)
2827:       m_logfile << "*";
2828:     m_logfile << std::endl;
2829:     m_logfile << "\t\t\t\tSphere inclusion (3D structured version 2 with barriers)" << std::endl;
2830:     for (unsigned int i = 0; i < 80; ++i)
2831:       m_logfile << "*";
2832:     m_logfile << std::endl;
2833: 
2834:     Triangulation<dim> tria_inner;
2835:     GridGenerator::hyper_ball(tria_inner, Point<dim>(), 0.49);
2836: 
2837:     Triangulation<dim> tria_outer;
2838:     GridGenerator::hyper_shell(
2839:       tria_outer, Point<dim>(), 0.49, std::sqrt(dim)*0.5, 2 * dim);
2840: 
2841:     Triangulation<dim> cube1;
2842:     GridGenerator::hyper_rectangle(cube1, Point<dim>(0, 0, 0.5), Point<dim>(1, 1, 1.5));
2843:     Triangulation<dim> cube2;
2844:     GridGenerator::hyper_rectangle(cube2, Point<dim>(0, 0.5, -0.5), Point<dim>(1, 1.5, 0.5));
2845:     Triangulation<dim> cube3;
2846:     GridGenerator::hyper_rectangle(cube3, Point<dim>(0.5, -0.5, -0.5), Point<dim>(1.5, 0.5, 0.5));
2847: 
2848:     Triangulation<dim> tmp_triangulation;
2849:     GridGenerator::merge_triangulations({&tria_inner, &tria_outer,
2850:                                          &cube1, &cube2, &cube3}, tmp_triangulation);
2851: 
2852:     tmp_triangulation.reset_all_manifolds();
2853:     tmp_triangulation.set_all_manifold_ids(0);
2854: 
2855:     for (const auto &cell : tmp_triangulation.cell_iterators())
2856:       {
2857:         for (const auto &face : cell->face_iterators())
2858:           {
2859:             bool face_at_sphere_boundary = true;
2860:             for (const auto v : face->vertex_indices())
2861:               {
2862:                 if (std::abs(face->vertex(v).norm_square() - 0.49 * 0.49) > 1e-12)
2863:                   face_at_sphere_boundary = false;
2864:               }
2865:             if (face_at_sphere_boundary)
2866:               face->set_all_manifold_ids(1);
2867:           }
2868:         if (cell->center().norm_square() < 0.1)
2869:           cell->set_material_id(1);
2870:         else
2871:           cell->set_material_id(0);
2872:       }
2873: 
2874:     tmp_triangulation.set_manifold(1, SphericalManifold<dim>());
2875: 
2876:     TransfiniteInterpolationManifold<dim> transfinite_manifold;
2877:     transfinite_manifold.initialize(tmp_triangulation);
2878:     tmp_triangulation.set_manifold(0, transfinite_manifold);
2879: 
2880:     tmp_triangulation.refine_global(m_parameters.m_global_refine_times);
2881: 
2882:     // some extra barriers
2883:     for (const auto &cell : tmp_triangulation.cell_iterators())
2884:       {
2885:         if (    std::fabs(cell->center()[1] - 0.75) < 0.05
2886:              && std::fabs(cell->center()[2] - 0.5625) < 0.05
2887:              && std::fabs(cell->center()[0] - 0.0) < 0.2)
2888:           cell->set_material_id(1);
2889: 
2890:         if (    std::fabs(cell->center()[1] - 0.0) < 0.2
2891:              && std::fabs(cell->center()[2] - 0.5) < 0.1
2892:              && std::fabs(cell->center()[0] - 0.75) < 0.05)
2893:           cell->set_material_id(1);
2894:       }
2895: 
2896:     std::set<typename Triangulation< dim >::active_cell_iterator >
2897:       cells_to_remove;
2898: 
2899:     for (const auto &cell : tmp_triangulation.active_cell_iterators())
2900:       {
2901: 	if (   cell->center()[0] < 0.0
2902: 	    || cell->center()[1] < 0.0
2903: 	    || cell->center()[2] < 0.0
2904: 	    || cell->center()[0] > 1.0
2905: 	    || cell->center()[1] > 1.0
2906: 	    || cell->center()[2] > 1.0)
2907: 	  {
2908: 	    cells_to_remove.insert(cell);
2909: 	  }
2910:       }
2911: 
2912:     GridGenerator::create_triangulation_with_removed_cells(tmp_triangulation,
2913: 							   cells_to_remove,
2914: 							   m_triangulation);
2915: 
2916:     for (const auto &cell : m_triangulation.active_cell_iterators())
2917:       for (const auto &face : cell->face_iterators())
2918: 	{
2919: 	  if (face->at_boundary() == true)
2920: 	    {
2921: 	      if (std::fabs(face->center()[0] - 0.0 ) < 1.0e-9 )
2922: 		face->set_boundary_id(0);
2923: 	      else if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9)
2924: 		face->set_boundary_id(1);
2925: 	      else if (std::fabs(face->center()[2] - 0.0 ) < 1.0e-9)
2926: 		face->set_boundary_id(2);
2927: 	      else if (std::fabs(face->center()[2] - 1.0 ) < 1.0e-9)
2928: 		face->set_boundary_id(3);
2929: 	      else
2930: 		face->set_boundary_id(4);
2931: 	    }
2932: 	}
2933: 
2934:     if (m_parameters.m_refinement_strategy == "adaptive-refine")
2935:       {
2936: 	unsigned int material_id;
2937: 	double length_scale;
2938: 	bool initiation_point_refine_unfinished = true;
2939: 	while (initiation_point_refine_unfinished)
2940: 	  {
2941: 	    initiation_point_refine_unfinished = false;
2942: 	    for (const auto &cell : m_triangulation.active_cell_iterators())
2943: 	      {
2944: 		if (    cell->center()[2] > 0.505
2945: 		     && cell->center()[2] < 0.575
2946: 		     && cell->center()[0] < 0.05
2947: 		     && cell->center()[1] < 0.05 )
2948: 		  {
2949: 		    material_id = cell->material_id();
2950: 		    length_scale = m_material_data[material_id][2];
2951: 		    if (  std::cbrt(cell->measure())
2952: 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
2953: 		      {
2954: 			cell->set_refine_flag();
2955: 			initiation_point_refine_unfinished = true;
2956: 		      }
2957: 		  }
2958: 	      }
2959: 	    m_triangulation.execute_coarsening_and_refinement();
2960: 	  }
2961:       }
2962:     else
2963:       {
2964: 	AssertThrow(false,
2965: 		    ExcMessage("Selected mesh refinement strategy not implemented!"));
2966:       }
2967:   }
2968: 
2969:   template <int dim>
2970:   void PhaseFieldMonolithicSolve<dim>::make_grid_case_9()
2971:   {
2972:     AssertThrow(dim==2, ExcMessage("The dimension has to be 2D!"));
2973: 
2974:     for (unsigned int i = 0; i < 80; ++i)
2975:       m_logfile << "*";
2976:     m_logfile << std::endl;
2977:     m_logfile << "\t\t\t\tL-shape bending (2D structured)" << std::endl;
2978:     for (unsigned int i = 0; i < 80; ++i)
2979:       m_logfile << "*";
2980:     m_logfile << std::endl;
2981: 
2982:     GridIn<dim> gridin;
2983:     gridin.attach_triangulation(m_triangulation);
2984:     std::ifstream f("L-Shape.msh");
2985:     gridin.read_msh(f);
2986: 
2987:     for (const auto &cell : m_triangulation.active_cell_iterators())
2988:       for (const auto &face : cell->face_iterators())
2989: 	{
2990: 	  if (face->at_boundary() == true)
2991: 	    {
2992: 	      if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9 )
2993: 		face->set_boundary_id(0);
2994: 	      else
2995: 	        face->set_boundary_id(1);
2996: 	    }
2997: 	}
2998: 
2999:     m_triangulation.refine_global(m_parameters.m_global_refine_times);
3000: 
3001:     if (m_parameters.m_refinement_strategy == "pre-refine")
3002:       {
3003: 	unsigned int material_id;
3004: 	double length_scale;
3005: 	for (unsigned int i = 0; i < m_parameters.m_local_prerefine_times; i++)
3006: 	  {
3007: 	    for (const auto &cell : m_triangulation.active_cell_iterators())
3008: 	      {
3009: 		if (    (cell->center()[1] > 242.0)
3010: 		     && (cell->center()[1] < 312.5)
3011: 		     && (cell->center()[0] < 258.0) )
3012: 		  {
3013: 		    material_id = cell->material_id();
3014: 		    length_scale = m_material_data[material_id][2];
3015: 		    if (  std::sqrt(cell->measure())
3016: 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
3017: 		      cell->set_refine_flag();
3018: 		  }
3019: 	      }
3020: 	    m_triangulation.execute_coarsening_and_refinement();
3021: 	  }
3022:       }
3023:     else if (m_parameters.m_refinement_strategy == "adaptive-refine")
3024:       {
3025: 	unsigned int material_id;
3026: 	double length_scale;
3027: 	bool initiation_point_refine_unfinished = true;
3028: 	while (initiation_point_refine_unfinished)
3029: 	  {
3030: 	    initiation_point_refine_unfinished = false;
3031: 	    for (const auto &cell : m_triangulation.active_cell_iterators())
3032: 	      {
3033: 		if (             (cell->center()[0] - 250) < 0.0
3034: 		     &&          (cell->center()[0] - 240) > 0.0
3035: 		     && std::fabs(cell->center()[1] - 250) < 10.0 )
3036: 		  {
3037: 		    material_id = cell->material_id();
3038: 		    length_scale = m_material_data[material_id][2];
3039: 		    if (  std::sqrt(cell->measure())
3040: 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
3041: 		      {
3042: 		        cell->set_refine_flag();
3043: 		        initiation_point_refine_unfinished = true;
3044: 		      }
3045: 		  }
3046: 	      }
3047: 	    m_triangulation.execute_coarsening_and_refinement();
3048: 	  }
3049:       }
3050:     else
3051:       {
3052: 	AssertThrow(false,
3053: 	            ExcMessage("Selected mesh refinement strategy not implemented!"));
3054:       }
3055:   }
3056: 
3057:   template <int dim>
3058:   void PhaseFieldMonolithicSolve<dim>::make_grid_case_10()
3059:   {
3060:     AssertThrow(dim==3, ExcMessage("The dimension has to be 3D!"));
3061: 
3062:     for (unsigned int i = 0; i < 80; ++i)
3063:       m_logfile << "*";
3064:     m_logfile << std::endl;
3065:     m_logfile << "\t\t\t\tL-shape bending (3D structured)" << std::endl;
3066:     for (unsigned int i = 0; i < 80; ++i)
3067:       m_logfile << "*";
3068:     m_logfile << std::endl;
3069: 
3070:     Triangulation<2> triangulation_2d;
3071: 
3072:     GridIn<2> gridin;
3073:     gridin.attach_triangulation(triangulation_2d);
3074:     std::ifstream f("L-Shape.msh");
3075:     gridin.read_msh(f);
3076: 
3077:     const double thickness = 150.0;
3078:     const unsigned int n_layer = 11;
3079:     GridGenerator::extrude_triangulation(triangulation_2d, n_layer, thickness, m_triangulation);
3080: 
3081:     for (const auto &cell : m_triangulation.active_cell_iterators())
3082:       for (const auto &face : cell->face_iterators())
3083: 	{
3084: 	  if (face->at_boundary() == true)
3085: 	    {
3086: 	      if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9 )
3087: 		face->set_boundary_id(0);
3088: 	      else
3089: 	        face->set_boundary_id(1);
3090: 	    }
3091: 	}
3092: 
3093:     m_triangulation.refine_global(m_parameters.m_global_refine_times);
3094: 
3095:     if (m_parameters.m_refinement_strategy == "pre-refine")
3096:       {
3097: 	unsigned int material_id;
3098: 	double length_scale;
3099: 	for (unsigned int i = 0; i < m_parameters.m_local_prerefine_times; i++)
3100: 	  {
3101: 	    for (const auto &cell : m_triangulation.active_cell_iterators())
3102: 	      {
3103: 		if (    (std::fabs(cell->center()[1] - 250.0) < 10.0)
3104: 		     && (cell->center()[0] < 250.0) )
3105: 		  {
3106: 		    material_id = cell->material_id();
3107: 		    length_scale = m_material_data[material_id][2];
3108: 		    if (  std::cbrt(cell->measure())
3109: 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
3110: 		      cell->set_refine_flag();
3111: 		  }
3112: 	      }
3113: 	    m_triangulation.execute_coarsening_and_refinement();
3114: 	  }
3115:       }
3116:     else if (m_parameters.m_refinement_strategy == "adaptive-refine")
3117:       {
3118: 	unsigned int material_id;
3119: 	double length_scale;
3120: 	bool initiation_point_refine_unfinished = true;
3121: 	while (initiation_point_refine_unfinished)
3122: 	  {
3123: 	    initiation_point_refine_unfinished = false;
3124: 	    for (const auto &cell : m_triangulation.active_cell_iterators())
3125: 	      {
3126: 		if (             (cell->center()[0] - 250) < 0.0
3127: 		     &&          (cell->center()[0] - 240) > 0.0
3128: 		     && std::fabs(cell->center()[1] - 250) < 10.0 )
3129: 		  {
3130: 		    material_id = cell->material_id();
3131: 		    length_scale = m_material_data[material_id][2];
3132: 		    if (  std::cbrt(cell->measure())
3133: 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
3134: 		      {
3135: 		        cell->set_refine_flag();
3136: 		        initiation_point_refine_unfinished = true;
3137: 		      }
3138: 		  }
3139: 	      }
3140: 	    m_triangulation.execute_coarsening_and_refinement();
3141: 	  }
3142:       }
3143:     else
3144:       {
3145: 	AssertThrow(false,
3146: 	            ExcMessage("Selected mesh refinement strategy not implemented!"));
3147:       }
3148:   }
3149: 
3150: 
3151:   template <int dim>
3152:   void PhaseFieldMonolithicSolve<dim>::make_grid_case_11()
3153:   {
3154:     AssertThrow(dim==3, ExcMessage("The dimension has to be 3D!"));
3155: 
3156:     for (unsigned int i = 0; i < 80; ++i)
3157:       m_logfile << "*";
3158:     m_logfile << std::endl;
3159:     m_logfile << "\t\t\t\tBrokenshire torsion (3D structured)" << std::endl;
3160:     for (unsigned int i = 0; i < 80; ++i)
3161:       m_logfile << "*";
3162:     m_logfile << std::endl;
3163: 
3164:     Triangulation<2> triangulation_2d;
3165: 
3166:     double const length = 200.0;
3167:     double const width = 50.0;
3168:     double const height = 50.0;
3169:     double const delta_L = 25.0;
3170:     double const tan_theta = delta_L / (0.5*width);
3171: 
3172:     std::vector<unsigned int> repetitions(2, 1);
3173:     repetitions[0] = 20;
3174:     repetitions[1] = 5;
3175: 
3176:     Point<2> point1(0.0, 0.0);
3177:     Point<2> point2(length, width);
3178: 
3179:     GridGenerator::subdivided_hyper_rectangle(triangulation_2d,
3180: 					      repetitions,
3181: 					      point1,
3182: 					      point2 );
3183: 
3184:     typename Triangulation<2>::vertex_iterator vertex_ptr;
3185:     vertex_ptr = triangulation_2d.begin_active_vertex();
3186:     while (vertex_ptr != triangulation_2d.end_vertex())
3187:       {
3188: 	Point<2> & vertex_point = vertex_ptr->vertex();
3189: 
3190: 	const double delta_x = (vertex_point(1) - 0.5*width) * tan_theta;
3191: 
3192: 	if (std::fabs(vertex_point(0) - 0.5*length) < 1.0e-6)
3193: 	  {
3194: 	    vertex_point(0) += delta_x;
3195: 	  }
3196: 	else if (std::fabs(vertex_point(0) + length/repetitions[0] - 0.5*length) < 1.0e-6)
3197: 	  {
3198: 	    vertex_point(0) += (delta_x + length/repetitions[0]*0.5);
3199: 	  }
3200: 	else if (std::fabs(vertex_point(0) - length/repetitions[0] - 0.5*length) < 1.0e-6)
3201: 	  {
3202: 	    vertex_point(0) += (delta_x - length/repetitions[0]*0.5);
3203: 	  }
3204: 	else if (vertex_point(0) < 0.5*length - length/repetitions[0] - 1.0e-6)
3205: 	  {
3206: 	    vertex_point(0) += (delta_x + length/repetitions[0]*0.5) * vertex_point(0)/(0.5*length - length/repetitions[0]);
3207: 	  }
3208: 	else if (vertex_point(0) > 0.5*length + length/repetitions[0] + 1.0e-6)
3209: 	  {
3210: 	    vertex_point(0) += (delta_x - length/repetitions[0]*0.5) * (length - vertex_point(0))/(0.5*length - length/repetitions[0]);
3211: 	  }
3212: 
3213: 	++vertex_ptr;
3214:       }
3215: 
3216:     Triangulation<dim> tmp_triangulation;
3217:     const unsigned int n_layer = repetitions[1] + 1;
3218:     GridGenerator::extrude_triangulation(triangulation_2d, n_layer, height, tmp_triangulation);
3219: 
3220:     tmp_triangulation.refine_global(m_parameters.m_global_refine_times);
3221: 
3222:     std::set<typename Triangulation< dim >::active_cell_iterator >
3223:       cells_to_remove;
3224: 
3225:     for (const auto &cell : tmp_triangulation.active_cell_iterators())
3226:       {
3227: 	if (    (std::fabs(cell->center()[0] - (cell->center()[1] - 0.5*width)*tan_theta - 0.5*length) < 2.5)
3228: 	     && cell->center()[2] > 0.5* height  )
3229: 	  {
3230: 	    cells_to_remove.insert(cell);
3231: 	  }
3232:       }
3233: 
3234:     GridGenerator::create_triangulation_with_removed_cells(tmp_triangulation,
3235: 							   cells_to_remove,
3236: 							   m_triangulation);
3237: 
3238:     if (m_parameters.m_refinement_strategy == "adaptive-refine")
3239:       {
3240: 	unsigned int material_id;
3241: 	double length_scale;
3242: 	bool initiation_point_refine_unfinished = true;
3243: 	while (initiation_point_refine_unfinished)
3244: 	  {
3245: 	    initiation_point_refine_unfinished = false;
3246: 	    for (const auto &cell : m_triangulation.active_cell_iterators())
3247: 	      {
3248: 		if (  (std::fabs(cell->center()[0] - (cell->center()[1] - 0.5*width)*tan_theta - 0.5*length) < 5.0)
3249: 		    && cell->center()[2] <= 0.5*height
3250: 		    && cell->center()[2] > 0.5*height - 5.0 )
3251: 		  {
3252: 		    material_id = cell->material_id();
3253: 		    length_scale = m_material_data[material_id][2];
3254: 		    if (  std::cbrt(cell->measure())
3255: 			> length_scale * m_parameters.m_allowed_max_h_l_ratio )
3256: 		      {
3257: 			cell->set_refine_flag();
3258: 			initiation_point_refine_unfinished = true;
3259: 		      }
3260: 		  }
3261: 	      }
3262: 	    m_triangulation.execute_coarsening_and_refinement();
3263: 	  }
3264:       }
3265:     else
3266:       {
3267: 	AssertThrow(false,
3268: 		    ExcMessage("Selected mesh refinement strategy not implemented!"));
3269:       }
3270: 
3271: 
3272:     for (const auto &cell : m_triangulation.active_cell_iterators())
3273:       for (const auto &face : cell->face_iterators())
3274: 	{
3275: 	  if (face->at_boundary() == true)
3276: 	    {
3277: 	      if (std::fabs(face->center()[0] - length) < 1.0e-6 )
3278: 		face->set_boundary_id(0);
3279: 	      else if (std::fabs(face->center()[0] - 0.0) < 1.0e-6 )
3280: 		face->set_boundary_id(1);
3281: 	      else
3282: 		face->set_boundary_id(2);
3283: 	    }
3284: 	}
3285:   }
3286: 
3287:   template <int dim>
3288:   void PhaseFieldMonolithicSolve<dim>::setup_system()
3289:   {
3290:     m_timer.enter_subsection("Setup system");
3291: 
3292:     std::vector<unsigned int> block_component(m_n_components,
3293:                                               m_u_dof); // displacement
3294:     block_component[m_d_component] = m_d_dof;           // phasefield
3295: 
3296:     m_dof_handler.distribute_dofs(m_fe);
3297:     DoFRenumbering::Cuthill_McKee(m_dof_handler);
3298:     DoFRenumbering::component_wise(m_dof_handler, block_component);
3299: 
3300:     m_constraints.clear();
3301:     DoFTools::make_hanging_node_constraints(m_dof_handler, m_constraints);
3302:     m_constraints.close();
3303: 
3304:     m_dofs_per_block =
3305:       DoFTools::count_dofs_per_fe_block(m_dof_handler, block_component);
3306: 
3307:     m_logfile << "\t\tTriangulation:"
3308:               << "\n\t\t\t Number of active cells: "
3309:               << m_triangulation.n_active_cells()
3310:               << "\n\t\t\t Number of used vertices: "
3311:               << m_triangulation.n_used_vertices()
3312:               << "\n\t\t\t Number of active edges: "
3313:               << m_triangulation.n_active_lines()
3314:               << "\n\t\t\t Number of active faces: "
3315:               << m_triangulation.n_active_faces()
3316:               << "\n\t\t\t Number of degrees of freedom (total): "
3317: 	      << m_dof_handler.n_dofs()
3318: 	      << "\n\t\t\t Number of degrees of freedom (disp): "
3319: 	      << m_dofs_per_block[m_u_dof]
3320: 	      << "\n\t\t\t Number of degrees of freedom (phasefield): "
3321: 	      << m_dofs_per_block[m_d_dof]
3322:               << std::endl;
3323: 
3324:     m_tangent_matrix.clear();
3325:     {
3326:       BlockDynamicSparsityPattern dsp(m_dofs_per_block, m_dofs_per_block);
3327: 
3328:       Table<2, DoFTools::Coupling> coupling(m_n_components, m_n_components);
3329:       for (unsigned int ii = 0; ii < m_n_components; ++ii)
3330:         for (unsigned int jj = 0; jj < m_n_components; ++jj)
3331:           {
3332:             if (   ((ii < m_d_component) && (jj == m_d_component))
3333:                 || ((ii == m_d_component) && (jj < m_d_component)) )
3334:               coupling[ii][jj] = DoFTools::none;
3335:             else
3336:               coupling[ii][jj] = DoFTools::always;
3337:           }
3338: 
3339:       DoFTools::make_sparsity_pattern(
3340:         m_dof_handler, coupling, dsp, m_constraints, false);
3341:       m_sparsity_pattern.copy_from(dsp);
3342:     }
3343: 
3344:     m_tangent_matrix.reinit(m_sparsity_pattern);
3345: 
3346:     m_system_rhs.reinit(m_dofs_per_block);
3347:     m_solution.reinit(m_dofs_per_block);
3348: 
3349:     m_active_set_phasefield.reinit(m_dofs_per_block[m_d_dof]);
3350: 
3351:     setup_qph();
3352: 
3353:     m_timer.leave_subsection();
3354:   }
3355: 
3356:   template <int dim>
3357:   void PhaseFieldMonolithicSolve<dim>::make_constraints(const unsigned int it_nr)
3358:   {
3359:     const bool apply_dirichlet_bc = (it_nr == 0);
3360: 
3361:     if (it_nr > 1)
3362:       {
3363:         return;
3364:       }
3365: 
3366:     if (apply_dirichlet_bc)
3367:       {
3368: 	m_constraints.clear();
3369: 	DoFTools::make_hanging_node_constraints(m_dof_handler,
3370: 						m_constraints);
3371: 
3372: 	const FEValuesExtractors::Scalar x_displacement(0);
3373: 	const FEValuesExtractors::Scalar y_displacement(1);
3374: 	const FEValuesExtractors::Scalar z_displacement(2);
3375: 
3376: 	const FEValuesExtractors::Vector displacements(0);
3377: 
3378: 	if (   m_parameters.m_scenario == 1
3379: 	    || m_parameters.m_scenario == 3)
3380: 	  {
3381: 	    // Dirichlet B,C. bottom surface
3382: 	    const int boundary_id_bottom_surface = 0;
3383: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3384: 						     boundary_id_bottom_surface,
3385: 						     Functions::ZeroFunction<dim>(m_n_components),
3386: 						     m_constraints,
3387: 						     m_fe.component_mask(y_displacement));
3388: 
3389: 	    typename Triangulation<dim>::active_vertex_iterator vertex_itr;
3390: 	    vertex_itr = m_triangulation.begin_active_vertex();
3391: 	    std::vector<types::global_dof_index> node_xy(m_fe.dofs_per_vertex);
3392: 
3393: 	    for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)
3394: 	      {
3395: 		if (   (std::fabs(vertex_itr->vertex()[0] - 0.0) < 1.0e-9)
3396: 		    && (std::fabs(vertex_itr->vertex()[1] - 0.0) < 1.0e-9) )
3397: 		  {
3398: 		    node_xy = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);
3399: 		  }
3400: 	      }
3401: 	    m_constraints.add_line(node_xy[0]);
3402: 	    m_constraints.set_inhomogeneity(node_xy[0], 0.0);
3403: 
3404: 	    m_constraints.add_line(node_xy[1]);
3405: 	    m_constraints.set_inhomogeneity(node_xy[1], 0.0);
3406: 
3407: 	    const int boundary_id_top_surface = 1;
3408: 	    /*
3409: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3410: 						     boundary_id_top_surface,
3411: 						     Functions::ZeroFunction<dim>(m_n_components),
3412: 						     m_constraints,
3413: 						     m_fe.component_mask(x_displacement));
3414: 	    */
3415:             const double time_inc = m_time.get_delta_t();
3416:             double disp_magnitude = m_time.get_magnitude();
3417: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3418: 						     boundary_id_top_surface,
3419: 						     Functions::ConstantFunction<dim>(
3420: 						       disp_magnitude*time_inc, m_n_components),
3421: 						     m_constraints,
3422: 						     m_fe.component_mask(y_displacement));
3423: 	  }
3424: 	else if (   m_parameters.m_scenario == 2
3425: 	         || m_parameters.m_scenario == 4)
3426: 	  {
3427: 	    // Dirichlet B,C. bottom surface
3428: 	    const int boundary_id_bottom_surface = 0;
3429: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3430: 						     boundary_id_bottom_surface,
3431: 						     Functions::ZeroFunction<dim>(m_n_components),
3432: 						     m_constraints,
3433: 						     m_fe.component_mask(displacements));
3434: 
3435: 	    const int boundary_id_top_surface = 1;
3436: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3437: 						     boundary_id_top_surface,
3438: 						     Functions::ZeroFunction<dim>(m_n_components),
3439: 						     m_constraints,
3440: 						     m_fe.component_mask(y_displacement));
3441: 
3442: 	    const double time_inc = m_time.get_delta_t();
3443: 	    double disp_magnitude = m_time.get_magnitude();
3444: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3445: 						     boundary_id_top_surface,
3446: 						     Functions::ConstantFunction<dim>(
3447: 						       disp_magnitude*time_inc, m_n_components),
3448: 						     m_constraints,
3449: 						     m_fe.component_mask(x_displacement));
3450: 
3451: 	    const int boundary_id_side_surfaces = 2;
3452: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3453: 						     boundary_id_side_surfaces,
3454: 						     Functions::ZeroFunction<dim>(m_n_components),
3455: 						     m_constraints,
3456: 						     m_fe.component_mask(y_displacement));
3457: 	  }
3458: 	else if (m_parameters.m_scenario == 5)
3459: 	  {
3460: 	    typename Triangulation<dim>::active_vertex_iterator vertex_itr;
3461: 	    vertex_itr = m_triangulation.begin_active_vertex();
3462: 	    std::vector<types::global_dof_index> node_bottomleft(m_fe.dofs_per_vertex);
3463: 	    std::vector<types::global_dof_index> node_bottomright(m_fe.dofs_per_vertex);
3464: 	    std::vector<types::global_dof_index> node_topcenter(m_fe.dofs_per_vertex);
3465: 
3466: 	    for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)
3467: 	      {
3468: 		if (   (std::fabs(vertex_itr->vertex()[0] - 0.0) < 1.0e-9)
3469: 		    && (std::fabs(vertex_itr->vertex()[1] - 0.0) < 1.0e-9) )
3470: 		  {
3471: 		    node_bottomleft = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);
3472: 		  }
3473: 		if (   (std::fabs(vertex_itr->vertex()[0] - 8.0) < 1.0e-9)
3474: 		    && (std::fabs(vertex_itr->vertex()[1] - 0.0) < 1.0e-9) )
3475: 		  {
3476: 		    node_bottomright = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);
3477: 		  }
3478: 		if (   (std::fabs(vertex_itr->vertex()[0] - 4.0) < 1.0e-9)
3479: 		    && (std::fabs(vertex_itr->vertex()[1] - 2.0) < 1.0e-9) )
3480: 		  {
3481: 		    node_topcenter = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);
3482: 		  }
3483: 	      }
3484: 	    // bottom-left node fixed in both x- and y-directions
3485: 	    m_constraints.add_line(node_bottomleft[0]);
3486: 	    m_constraints.set_inhomogeneity(node_bottomleft[0], 0.0);
3487: 
3488: 	    m_constraints.add_line(node_bottomleft[1]);
3489: 	    m_constraints.set_inhomogeneity(node_bottomleft[1], 0.0);
3490: 
3491: 	    // bottom-right node only fixed in y-direction
3492: 	    m_constraints.add_line(node_bottomright[1]);
3493: 	    m_constraints.set_inhomogeneity(node_bottomright[1], 0.0);
3494: 
3495: 	    // top-center node applied with y-displacement
3496: 	    const double time_inc = m_time.get_delta_t();
3497: 	    double disp_magnitude = m_time.get_magnitude();
3498: 
3499: 	    m_constraints.add_line(node_topcenter[1]);
3500: 	    m_constraints.set_inhomogeneity(node_topcenter[1], disp_magnitude*time_inc);
3501: 	  }
3502: 	else if (   m_parameters.m_scenario == 6
3503: 	         || m_parameters.m_scenario == 7
3504: 		 || m_parameters.m_scenario == 8)
3505: 	  {
3506: 	    const int x0_surface = 0;
3507: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3508: 						     x0_surface,
3509: 						     Functions::ZeroFunction<dim>(m_n_components),
3510: 						     m_constraints,
3511: 						     m_fe.component_mask(x_displacement));
3512: 	    const int y0_surface = 1;
3513: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3514: 						     y0_surface,
3515: 						     Functions::ZeroFunction<dim>(m_n_components),
3516: 						     m_constraints,
3517: 						     m_fe.component_mask(y_displacement));
3518: 	    const int z0_surface = 2;
3519: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3520: 						     z0_surface,
3521: 						     Functions::ZeroFunction<dim>(m_n_components),
3522: 						     m_constraints,
3523: 						     m_fe.component_mask(z_displacement));
3524: 
3525: 	    const int z1_surface = 3;
3526: 	    const double time_inc = m_time.get_delta_t();
3527: 	    double disp_magnitude = 1.0;
3528: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3529: 						     z1_surface,
3530: 						     Functions::ConstantFunction<dim>(
3531: 						       disp_magnitude*time_inc, m_n_components),
3532: 						     m_constraints,
3533: 						     m_fe.component_mask(z_displacement));
3534: 	  }
3535: 	else if (   m_parameters.m_scenario == 9
3536: 	         || m_parameters.m_scenario == 10)
3537: 	  {
3538: 	    // Dirichlet B,C. bottom surface
3539: 	    const int boundary_id_bottom_surface = 0;
3540: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3541: 						     boundary_id_bottom_surface,
3542: 						     Functions::ZeroFunction<dim>(m_n_components),
3543: 						     m_constraints,
3544: 						     m_fe.component_mask(displacements));
3545: 
3546: 	    typename Triangulation<dim>::active_vertex_iterator vertex_itr;
3547: 	    vertex_itr = m_triangulation.begin_active_vertex();
3548: 	    std::vector<types::global_dof_index> node_disp_control(m_fe.dofs_per_vertex);
3549: 
3550: 	    for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)
3551: 	      {
3552: 		if (   (std::fabs(vertex_itr->vertex()[0] - 470.0) < 1.0e-9)
3553: 		    && (std::fabs(vertex_itr->vertex()[1] - 250.0) < 1.0e-9) )
3554: 		  {
3555: 		    node_disp_control = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);
3556: 	            // node applied with y-displacement
3557: 		    const double time_inc = m_time.get_delta_t();
3558: 		    double disp_magnitude = m_time.get_magnitude();
3559: 
3560: 		    m_constraints.add_line(node_disp_control[1]);
3561: 		    m_constraints.set_inhomogeneity(node_disp_control[1], disp_magnitude*time_inc);
3562: 		  }
3563: 	      }
3564: 	  }
3565: 	else if (m_parameters.m_scenario == 11)
3566: 	  {
3567: 	    // Dirichlet B,C. right surface
3568: 	    const int boundary_id_right_surface = 0;
3569: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3570: 						     boundary_id_right_surface,
3571: 						     Functions::ZeroFunction<dim>(m_n_components),
3572: 						     m_constraints,
3573: 						     m_fe.component_mask(displacements));
3574: 
3575: 	    // Dirichlet B,C. left surface
3576: 	    const int boundary_id_left_surface = 1;
3577: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3578: 						     boundary_id_left_surface,
3579: 						     Functions::ZeroFunction<dim>(m_n_components),
3580: 						     m_constraints,
3581: 						     m_fe.component_mask(x_displacement));
3582: 
3583: 	    typename Triangulation<dim>::active_vertex_iterator vertex_itr;
3584: 	    vertex_itr = m_triangulation.begin_active_vertex();
3585: 	    std::vector<types::global_dof_index> node_rotate(m_fe.dofs_per_vertex);
3586: 	    double node_dist = 0.0;
3587: 	    double disp_mag = 0.0;
3588: 	    double angle_theta = 0.0;
3589: 	    double disp_y = 0;
3590: 	    double disp_z = 0;
3591: 
3592: 	    for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)
3593: 	      {
3594: 		if (std::fabs(vertex_itr->vertex()[0] - 0.0) < 1.0e-9)
3595: 		  {
3596: 		    node_rotate = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);
3597: 		    node_dist = std::sqrt(  vertex_itr->vertex()[1] * vertex_itr->vertex()[1]
3598: 			                  + vertex_itr->vertex()[2] * vertex_itr->vertex()[2]);
3599: 
3600: 		    angle_theta = m_time.get_delta_t() * m_time.get_magnitude();
3601: 		    disp_mag = node_dist * std::tan(angle_theta);
3602: 
3603: 		    if (node_dist > 0)
3604: 		      {
3605: 		        disp_y = vertex_itr->vertex()[2]/node_dist * disp_mag;
3606: 		        disp_z = -vertex_itr->vertex()[1]/node_dist * disp_mag;
3607: 		      }
3608: 		    else
3609: 		      {
3610: 			disp_y = 0.0;
3611: 			disp_z = 0.0;
3612: 		      }
3613: 
3614: 		    m_constraints.add_line(node_rotate[1]);
3615: 		    m_constraints.set_inhomogeneity(node_rotate[1], disp_y);
3616: 
3617: 		    m_constraints.add_line(node_rotate[2]);
3618: 		    m_constraints.set_inhomogeneity(node_rotate[2], disp_z);
3619: 		  }
3620: 	      }
3621: 	  }
3622: 	else
3623: 	  Assert(false, ExcMessage("The scenario has not been implemented!"));
3624:       }
3625:     else  // inhomogeneous constraints
3626:       {
3627:         if (m_constraints.has_inhomogeneities())
3628:           {
3629:             AffineConstraints<double> homogeneous_constraints(m_constraints);
3630:             for (unsigned int dof = 0; dof != m_dof_handler.n_dofs(); ++dof)
3631:               if (homogeneous_constraints.is_inhomogeneously_constrained(dof))
3632:                 homogeneous_constraints.set_inhomogeneity(dof, 0.0);
3633:             m_constraints.clear();
3634:             m_constraints.copy_from(homogeneous_constraints);
3635:           }
3636:       }
3637:     m_constraints.close();
3638:   }
3639: 
3640:   template <int dim>
3641:   void PhaseFieldMonolithicSolve<dim>::assemble_system_B0(const BlockVector<double> & solution_old)
3642:   {
3643:     m_timer.enter_subsection("Assemble B0");
3644: 
3645:     m_tangent_matrix = 0.0;
3646: 
3647:     const UpdateFlags uf_cell(update_values | update_gradients |
3648: 			      update_quadrature_points | update_JxW_values);
3649:     const UpdateFlags uf_face(update_values | update_normal_vectors |
3650:                               update_JxW_values);
3651: 
3652:     PerTaskData_ASM per_task_data(m_fe.n_dofs_per_cell());
3653:     ScratchData_ASM scratch_data(m_fe, m_qf_cell, uf_cell, m_qf_face, uf_face, solution_old);
3654: 
3655:     auto worker =
3656:       [this](const typename DoFHandler<dim>::active_cell_iterator &cell,
3657: 	     ScratchData_ASM & scratch,
3658: 	     PerTaskData_ASM & data)
3659:       {
3660:         this->assemble_system_B0_one_cell(cell, scratch, data);
3661:       };
3662: 
3663:     auto copier = [this](const PerTaskData_ASM &data)
3664:       {
3665:         this->m_constraints.distribute_local_to_global(data.m_cell_matrix,
3666:                                                        data.m_local_dof_indices,
3667: 						       m_tangent_matrix);
3668:       };
3669: 
3670:     WorkStream::run(
3671:       m_dof_handler.active_cell_iterators(),
3672:       worker,
3673:       copier,
3674:       scratch_data,
3675:       per_task_data);
3676: 
3677:     m_timer.leave_subsection();
3678:   }
3679: 
3680:   template <int dim>
3681:   void PhaseFieldMonolithicSolve<dim>::assemble_system_rhs_BFGS_parallel(const BlockVector<double> & solution_old,
3682: 								         BlockVector<double> & system_rhs)
3683:   {
3684:     m_timer.enter_subsection("Assemble RHS");
3685: 
3686:     system_rhs = 0.0;
3687: 
3688:     const UpdateFlags uf_cell(update_values | update_gradients |
3689: 			      update_quadrature_points | update_JxW_values);
3690:     const UpdateFlags uf_face(update_values | update_normal_vectors |
3691: 			      update_JxW_values);
3692: 
3693:     PerTaskData_ASM_RHS_BFGS per_task_data(m_fe.n_dofs_per_cell());
3694:     ScratchData_ASM_RHS_BFGS scratch_data(m_fe, m_qf_cell, uf_cell, m_qf_face, uf_face, solution_old);
3695: 
3696:     auto worker =
3697:       [this](const typename DoFHandler<dim>::active_cell_iterator &cell,
3698: 	     ScratchData_ASM_RHS_BFGS & scratch,
3699: 	     PerTaskData_ASM_RHS_BFGS & data)
3700:       {
3701:         this->assemble_system_rhs_BFGS_one_cell(cell, scratch, data);
3702:       };
3703: 
3704:     auto copier = [this, &system_rhs](const PerTaskData_ASM_RHS_BFGS &data)
3705:       {
3706:         this->m_constraints.distribute_local_to_global(data.m_cell_rhs,
3707:                                                        data.m_local_dof_indices,
3708: 						       system_rhs);
3709:       };
3710: 
3711:     WorkStream::run(
3712:       m_dof_handler.active_cell_iterators(),
3713:       worker,
3714:       copier,
3715:       scratch_data,
3716:       per_task_data);
3717: 
3718:     m_timer.leave_subsection();
3719:   }
3720: 
3721:   template <int dim>
3722:   void PhaseFieldMonolithicSolve<dim>::assemble_system_rhs_BFGS_one_cell(
3723:       const typename DoFHandler<dim>::active_cell_iterator &cell,
3724:       ScratchData_ASM_RHS_BFGS & scratch,
3725:       PerTaskData_ASM_RHS_BFGS & data) const
3726:   {
3727:     data.reset();
3728:     scratch.reset();
3729:     scratch.m_fe_values.reinit(cell);
3730:     cell->get_dof_indices(data.m_local_dof_indices);
3731: 
3732:     scratch.m_fe_values[m_d_fe].get_function_values(
3733:       scratch.m_solution_previous_step, scratch.m_phasefield_previous_step_cell);
3734: 
3735:     const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
3736:       m_quadrature_point_history.get_data(cell);
3737:     Assert(lqph.size() == m_n_q_points, ExcInternalError());
3738: 
3739:     const double time_ramp = (m_time.current() / m_time.end());
3740:     std::vector<Tensor<1, dim>> rhs_values(m_n_q_points);
3741: 
3742:     right_hand_side(scratch.m_fe_values.get_quadrature_points(),
3743: 		    rhs_values,
3744: 		    m_parameters.m_x_component*1.0,
3745: 		    m_parameters.m_y_component*1.0,
3746: 		    m_parameters.m_z_component*1.0);
3747: 
3748:     const double delta_time = m_time.get_delta_t();
3749: 
3750:     for (const unsigned int q_point : scratch.m_fe_values.quadrature_point_indices())
3751:       {
3752:         for (const unsigned int k : scratch.m_fe_values.dof_indices())
3753:           {
3754:             const unsigned int k_group = m_fe.system_to_base_index(k).first.first;
3755: 
3756:             if (k_group == m_u_dof)
3757:               {
3758:                 scratch.m_Nx_disp[q_point][k] =
3759:                   scratch.m_fe_values[m_u_fe].value(k, q_point);
3760:                 scratch.m_grad_Nx_disp[q_point][k] =
3761:                   scratch.m_fe_values[m_u_fe].gradient(k, q_point);
3762:                 scratch.m_symm_grad_Nx_disp[q_point][k] =
3763:                   symmetrize(scratch.m_grad_Nx_disp[q_point][k]);
3764:               }
3765:             else if (k_group == m_d_dof)
3766:               {
3767: 		scratch.m_Nx_phasefield[q_point][k] =
3768: 		  scratch.m_fe_values[m_d_fe].value(k, q_point);
3769: 		scratch.m_grad_Nx_phasefield[q_point][k] =
3770: 		  scratch.m_fe_values[m_d_fe].gradient(k, q_point);
3771:               }
3772:             else
3773:               Assert(k_group <= m_d_dof, ExcInternalError());
3774:           }
3775:       }
3776: 
3777:     for (const unsigned int q_point : scratch.m_fe_values.quadrature_point_indices())
3778:       {
3779: 	const double length_scale            = lqph[q_point]->get_length_scale();
3780: 	const double gc                      = lqph[q_point]->get_critical_energy_release_rate();
3781: 	const double eta                     = lqph[q_point]->get_viscosity();
3782: 	const double current_positive_strain_energy = lqph[q_point]->get_current_positive_strain_energy();
3783: 
3784: 	const double phasefield_value        = lqph[q_point]->get_phase_field_value();
3785: 	const Tensor<1, dim> phasefield_grad = lqph[q_point]->get_phase_field_gradient();
3786: 
3787:         const std::vector<double>         &      N_phasefield = scratch.m_Nx_phasefield[q_point];
3788:         const std::vector<Tensor<1, dim>> & grad_N_phasefield = scratch.m_grad_Nx_phasefield[q_point];
3789:         const double                old_phasefield = scratch.m_phasefield_previous_step_cell[q_point];
3790: 
3791:         const SymmetricTensor<2, dim> & cauchy_stress = lqph[q_point]->get_cauchy_stress();
3792: 
3793:         const std::vector<Tensor<1,dim>> & N_disp = scratch.m_Nx_disp[q_point];
3794:         const std::vector<SymmetricTensor<2, dim>> & symm_grad_N_disp =
3795:           scratch.m_symm_grad_Nx_disp[q_point];
3796:         const double JxW = scratch.m_fe_values.JxW(q_point);
3797: 
3798:         SymmetricTensor<2, dim> symm_grad_Nx_i_x_C;
3799: 
3800:         for (const unsigned int i : scratch.m_fe_values.dof_indices())
3801:           {
3802:             const unsigned int i_group = m_fe.system_to_base_index(i).first.first;
3803: 
3804:             if (i_group == m_u_dof)
3805:               {
3806:                 data.m_cell_rhs(i) += (symm_grad_N_disp[i] * cauchy_stress) * JxW;
3807: 
3808: 		// contributions from the body force to right-hand side
3809: 		data.m_cell_rhs(i) -= N_disp[i] * rhs_values[q_point] * JxW;
3810:               }
3811:             else if (i_group == m_d_dof)
3812:               {
3813:     	        data.m_cell_rhs(i) += (    gc * length_scale * grad_N_phasefield[i] * phasefield_grad
3814:     	                                +  (   gc / length_scale * phasefield_value
3815: 					     + eta / delta_time  * (phasefield_value - old_phasefield)
3816: 					     + degradation_function_derivative(phasefield_value)
3817: 					     * current_positive_strain_energy )
3818: 					  * N_phasefield[i]
3819: 				      ) * JxW;
3820:               }
3821:             else
3822:               Assert(i_group <= m_d_dof, ExcInternalError());
3823:           }  // i
3824:       }  // q_point
3825: 
3826:     // if there is surface pressure, this surface pressure always applied to the
3827:     // reference configuration
3828:     const unsigned int face_pressure_id = 100;
3829:     const double p0 = 0.0;
3830: 
3831:     for (const auto &face : cell->face_iterators())
3832:       if (face->at_boundary() && face->boundary_id() == face_pressure_id)
3833:         {
3834:           scratch.m_fe_face_values.reinit(cell, face);
3835: 
3836:           for (const unsigned int f_q_point : scratch.m_fe_face_values.quadrature_point_indices())
3837:             {
3838:               const Tensor<1, dim> &N = scratch.m_fe_face_values.normal_vector(f_q_point);
3839: 
3840:               const double         pressure  = p0 * time_ramp;
3841:               const Tensor<1, dim> traction  = pressure * N;
3842: 
3843:               for (const unsigned int i : scratch.m_fe_values.dof_indices())
3844:                 {
3845:                   const unsigned int i_group = m_fe.system_to_base_index(i).first.first;
3846: 
3847:                   if (i_group == m_u_dof)
3848:                     {
3849:     		      const unsigned int component_i = m_fe.system_to_component_index(i).first;
3850:     		      const double Ni = scratch.m_fe_face_values.shape_value(i, f_q_point);
3851:     		      const double JxW = scratch.m_fe_face_values.JxW(f_q_point);
3852:     		      data.m_cell_rhs(i) -= (Ni * traction[component_i]) * JxW;
3853:                     }
3854:                 }
3855:             }
3856:         }
3857:   }
3858: 
3859:   template <int dim>
3860:   void PhaseFieldMonolithicSolve<dim>::assemble_system_B0_one_cell(
3861:       const typename DoFHandler<dim>::active_cell_iterator &cell,
3862:       ScratchData_ASM & scratch,
3863:       PerTaskData_ASM & data) const
3864:   {
3865:     data.reset();
3866:     scratch.reset();
3867:     scratch.m_fe_values.reinit(cell);
3868:     cell->get_dof_indices(data.m_local_dof_indices);
3869: 
3870:     scratch.m_fe_values[m_d_fe].get_function_values(
3871:       scratch.m_solution_previous_step, scratch.m_phasefield_previous_step_cell);
3872: 
3873:     const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
3874:       m_quadrature_point_history.get_data(cell);
3875:     Assert(lqph.size() == m_n_q_points, ExcInternalError());
3876: 
3877:     const double delta_time = m_time.get_delta_t();
3878: 
3879:     for (const unsigned int q_point : scratch.m_fe_values.quadrature_point_indices())
3880:       {
3881:         for (const unsigned int k : scratch.m_fe_values.dof_indices())
3882:           {
3883:             const unsigned int k_group = m_fe.system_to_base_index(k).first.first;
3884: 
3885:             if (k_group == m_u_dof)
3886:               {
3887:                 scratch.m_Nx_disp[q_point][k] =
3888:                   scratch.m_fe_values[m_u_fe].value(k, q_point);
3889:                 scratch.m_grad_Nx_disp[q_point][k] =
3890:                   scratch.m_fe_values[m_u_fe].gradient(k, q_point);
3891:                 scratch.m_symm_grad_Nx_disp[q_point][k] =
3892:                   symmetrize(scratch.m_grad_Nx_disp[q_point][k]);
3893:               }
3894:             else if (k_group == m_d_dof)
3895:               {
3896: 		scratch.m_Nx_phasefield[q_point][k] =
3897: 		  scratch.m_fe_values[m_d_fe].value(k, q_point);
3898: 		scratch.m_grad_Nx_phasefield[q_point][k] =
3899: 		  scratch.m_fe_values[m_d_fe].gradient(k, q_point);
3900:               }
3901:             else
3902:               Assert(k_group <= m_d_dof, ExcInternalError());
3903:           }
3904:       }
3905: 
3906:     for (const unsigned int q_point : scratch.m_fe_values.quadrature_point_indices())
3907:       {
3908: 	const double length_scale            = lqph[q_point]->get_length_scale();
3909: 	const double gc                      = lqph[q_point]->get_critical_energy_release_rate();
3910: 	const double eta                     = lqph[q_point]->get_viscosity();
3911: 	const double current_positive_strain_energy = lqph[q_point]->get_current_positive_strain_energy();
3912: 
3913: 	const double phasefield_value        = lqph[q_point]->get_phase_field_value();
3914: 
3915:         const std::vector<double>         &      N_phasefield = scratch.m_Nx_phasefield[q_point];
3916:         const std::vector<Tensor<1, dim>> & grad_N_phasefield = scratch.m_grad_Nx_phasefield[q_point];
3917: 
3918:         //const SymmetricTensor<2, dim> & cauchy_stress_positive = lqph[q_point]->get_cauchy_stress_positive();
3919:         const SymmetricTensor<4, dim> & mechanical_C  = lqph[q_point]->get_mechanical_C();
3920: 
3921:         const std::vector<SymmetricTensor<2, dim>> & symm_grad_N_disp =
3922:           scratch.m_symm_grad_Nx_disp[q_point];
3923:         const double JxW = scratch.m_fe_values.JxW(q_point);
3924: 
3925:         SymmetricTensor<2, dim> symm_grad_Nx_i_x_C;
3926: 
3927:         for (const unsigned int i : scratch.m_fe_values.dof_indices())
3928:           {
3929:             const unsigned int i_group = m_fe.system_to_base_index(i).first.first;
3930: 
3931:             if (i_group == m_u_dof)
3932:               {
3933:                 symm_grad_Nx_i_x_C = symm_grad_N_disp[i] * mechanical_C;
3934:               }
3935: 
3936:             for (const unsigned int j : scratch.m_fe_values.dof_indices())
3937:               {
3938:                 const unsigned int j_group = m_fe.system_to_base_index(j).first.first;
3939: 
3940:                 if ((i_group == j_group) && (i_group == m_u_dof))
3941:                   {
3942:                     data.m_cell_matrix(i, j) += symm_grad_Nx_i_x_C * symm_grad_N_disp[j] * JxW;
3943:                   }
3944:                 else if ((i_group == j_group) && (i_group == m_d_dof))
3945:                   {
3946:                     data.m_cell_matrix(i, j) += (  (   gc/length_scale + eta/delta_time
3947:                 	                             + degradation_function_2nd_order_derivative(phasefield_value)
3948: 						     * current_positive_strain_energy  )
3949:                 	                          * N_phasefield[i] * N_phasefield[j]
3950: 					          + gc * length_scale * grad_N_phasefield[i] * grad_N_phasefield[j]
3951: 					        ) * JxW;
3952:                   }
3953:                 else
3954:                   Assert((i_group <= m_d_dof) && (j_group <= m_d_dof),
3955:                          ExcInternalError());
3956:               } // j
3957:           }  // i
3958:       }  // q_point
3959:   }
3960: 
3961:   template <int dim>
3962:   void PhaseFieldMonolithicSolve<dim>::assemble_system_rhs_BFGS(const BlockVector<double> & solution_old,
3963: 								BlockVector<double> & system_rhs)
3964:   {
3965:     m_timer.enter_subsection("Assemble RHS");
3966: 
3967:     system_rhs = 0.0;
3968: 
3969:     Vector<double> cell_rhs(m_dofs_per_cell);
3970:     std::vector<types::global_dof_index> local_dof_indices(m_dofs_per_cell);
3971: 
3972:     const double time_ramp = (m_time.current() / m_time.end());
3973:     const double delta_time = m_time.get_delta_t();
3974: 
3975:     std::vector<Tensor<1, dim>> rhs_values(m_n_q_points);
3976:     const UpdateFlags uf_cell(update_values | update_gradients |
3977: 			      update_quadrature_points | update_JxW_values);
3978:     const UpdateFlags uf_face(update_values | update_normal_vectors |
3979: 			      update_JxW_values);
3980: 
3981:     FEValues<dim> fe_values(m_fe, m_qf_cell, uf_cell);
3982:     FEFaceValues<dim> fe_face_values(m_fe, m_qf_face, uf_face);
3983: 
3984:     // shape function values for displacement field
3985:     std::vector<std::vector<Tensor<1, dim>>>
3986:       Nx_disp(m_qf_cell.size(), std::vector<Tensor<1, dim>>(m_dofs_per_cell));
3987:     std::vector<std::vector<Tensor<2, dim>>>
3988:       grad_Nx_disp(m_qf_cell.size(), std::vector<Tensor<2, dim>>(m_dofs_per_cell));
3989:     std::vector<std::vector<SymmetricTensor<2, dim>>>
3990:       symm_grad_Nx_disp(m_qf_cell.size(), std::vector<SymmetricTensor<2, dim>>(m_dofs_per_cell));
3991: 
3992:     // shape function values for phase field
3993:     std::vector<std::vector<double>>
3994:       Nx_phasefield(m_qf_cell.size(), std::vector<double>(m_dofs_per_cell));
3995:     std::vector<std::vector<Tensor<1, dim>>>
3996:       grad_Nx_phasefield(m_qf_cell.size(), std::vector<Tensor<1, dim>>(m_dofs_per_cell));
3997: 
3998:     std::vector<double> phasefield_previous_step_cell(m_qf_cell.size());
3999: 
4000:     for (const auto &cell : m_dof_handler.active_cell_iterators())
4001:       {
4002: 	const std::vector<std::shared_ptr< PointHistory<dim>>> lqph =
4003: 	  m_quadrature_point_history.get_data(cell);
4004: 	Assert(lqph.size() == m_n_q_points, ExcInternalError());
4005: 
4006: 	cell_rhs = 0.0;
4007: 	fe_values.reinit(cell);
4008: 	right_hand_side(fe_values.get_quadrature_points(),
4009: 			rhs_values,
4010: 			m_parameters.m_x_component*time_ramp,
4011: 			m_parameters.m_y_component*time_ramp,
4012: 			m_parameters.m_z_component*time_ramp);
4013: 
4014: 	fe_values[m_d_fe].get_function_values(
4015: 	    solution_old, phasefield_previous_step_cell);
4016: 
4017: 	for (const unsigned int q_point : fe_values.quadrature_point_indices())
4018: 	  {
4019: 	    for (const unsigned int k : fe_values.dof_indices())
4020: 	      {
4021: 		const unsigned int k_group = m_fe.system_to_base_index(k).first.first;
4022: 
4023: 		if (k_group == m_u_dof)
4024: 		  {
4025: 		    Nx_disp[q_point][k] = fe_values[m_u_fe].value(k, q_point);
4026: 		    grad_Nx_disp[q_point][k] = fe_values[m_u_fe].gradient(k, q_point);
4027: 		    symm_grad_Nx_disp[q_point][k] = symmetrize(grad_Nx_disp[q_point][k]);
4028: 		  }
4029: 		else if (k_group == m_d_dof)
4030: 		  {
4031: 		    Nx_phasefield[q_point][k] = fe_values[m_d_fe].value(k, q_point);
4032: 		    grad_Nx_phasefield[q_point][k] = fe_values[m_d_fe].gradient(k, q_point);
4033: 		  }
4034: 		else
4035: 		  Assert(k_group <= m_d_dof, ExcInternalError());
4036: 	      }
4037: 	  }
4038: 
4039: 	for (const unsigned int q_point : fe_values.quadrature_point_indices())
4040: 	  {
4041: 	    const double length_scale            = lqph[q_point]->get_length_scale();
4042: 	    const double gc                      = lqph[q_point]->get_critical_energy_release_rate();
4043: 	    const double eta                     = lqph[q_point]->get_viscosity();
4044: 	    const double current_positive_strain_energy = lqph[q_point]->get_current_positive_strain_energy();
4045: 
4046: 	    const double phasefield_value        = lqph[q_point]->get_phase_field_value();
4047: 	    const Tensor<1, dim> phasefield_grad = lqph[q_point]->get_phase_field_gradient();
4048: 
4049: 	    const std::vector<double>         &      N_phasefield = Nx_phasefield[q_point];
4050: 	    const std::vector<Tensor<1, dim>> & grad_N_phasefield = grad_Nx_phasefield[q_point];
4051: 	    const double                old_phasefield = phasefield_previous_step_cell[q_point];
4052: 
4053: 	    const SymmetricTensor<2, dim> & cauchy_stress = lqph[q_point]->get_cauchy_stress();
4054: 
4055: 	    const std::vector<Tensor<1,dim>> & N = Nx_disp[q_point];
4056: 	    const std::vector<SymmetricTensor<2, dim>> & symm_grad_N = symm_grad_Nx_disp[q_point];
4057: 	    const double JxW = fe_values.JxW(q_point);
4058: 
4059: 	    for (const unsigned int i : fe_values.dof_indices())
4060: 	      {
4061: 		const unsigned int i_group = m_fe.system_to_base_index(i).first.first;
4062: 
4063: 		if (i_group == m_u_dof)
4064: 		  {
4065: 		    cell_rhs(i) += (symm_grad_N[i] * cauchy_stress) * JxW;
4066: 		    // contributions from the body force to right-hand side
4067: 		    cell_rhs(i) -= N[i] * rhs_values[q_point] * JxW;
4068: 		  }
4069: 		else if (i_group == m_d_dof)
4070: 		  {
4071: 		    cell_rhs(i) += (    gc * length_scale * grad_N_phasefield[i] * phasefield_grad
4072: 	    	                     +  (   gc / length_scale * phasefield_value
4073: 			                  + eta / delta_time  * (phasefield_value - old_phasefield)
4074: 				          + degradation_function_derivative(phasefield_value)
4075: 					  * current_positive_strain_energy )
4076: 				     * N_phasefield[i]
4077: 				   ) * JxW;
4078: 		  }
4079: 		else
4080: 		  Assert(i_group <= m_d_dof, ExcInternalError());
4081: 	      }
4082: 	  }
4083: 
4084: 	// if there is surface pressure, this surface pressure always applied to the
4085: 	// reference configuration
4086: 	const unsigned int face_pressure_id = 100;
4087: 	const double p0 = 0.0;
4088: 
4089: 	for (const auto &face : cell->face_iterators())
4090: 	  {
4091: 	    if (face->at_boundary() && face->boundary_id() == face_pressure_id)
4092: 	      {
4093: 		fe_face_values.reinit(cell, face);
4094: 
4095: 		for (const unsigned int f_q_point : fe_face_values.quadrature_point_indices())
4096: 		  {
4097: 		    const Tensor<1, dim> &N = fe_face_values.normal_vector(f_q_point);
4098: 
4099: 		    const double         pressure  = p0 * time_ramp;
4100: 		    const Tensor<1, dim> traction  = pressure * N;
4101: 
4102: 		    for (const unsigned int i : fe_values.dof_indices())
4103: 		      {
4104: 			const unsigned int i_group = m_fe.system_to_base_index(i).first.first;
4105: 
4106: 			if (i_group == m_u_dof)
4107: 			  {
4108: 			    const unsigned int component_i = m_fe.system_to_component_index(i).first;
4109: 			    const double Ni = fe_face_values.shape_value(i, f_q_point);
4110: 			    const double JxW = fe_face_values.JxW(f_q_point);
4111: 			    cell_rhs(i) -= (Ni * traction[component_i]) * JxW;
4112: 			  }
4113: 		      }
4114: 		  }
4115: 	      }
4116: 	  }
4117: 
4118: 	cell->get_dof_indices(local_dof_indices);
4119: 	for (const unsigned int i : fe_values.dof_indices())
4120: 	  system_rhs(local_dof_indices[i]) += cell_rhs(i);
4121:       } // for (const auto &cell : m_dof_handler.active_cell_iterators())
4122: 
4123:     m_timer.leave_subsection();
4124:   }
4125: 
4126: 
4127: 
4128:   template <int dim>
4129:   double PhaseFieldMonolithicSolve<dim>::line_search_stepsize_gradient_based(const BlockVector<double> & BFGS_p_vector,
4130: 				                                             const BlockVector<double> & solution_delta)
4131:   {
4132:     BlockVector<double> g_old(m_system_rhs);
4133: 
4134:     // BFGS_p_vector is the search direction
4135:     BlockVector<double> solution_delta_trial(solution_delta);
4136:     // take a full step size 1.0
4137:     solution_delta_trial.add(1.0, BFGS_p_vector);
4138: 
4139:     update_qph_incremental(solution_delta_trial, m_solution);
4140: 
4141:     BlockVector<double> g_new(m_dofs_per_block);
4142:     assemble_system_rhs_BFGS_parallel(m_solution, g_new);
4143: 
4144:     BlockVector<double> y_old(m_dofs_per_block);
4145: 
4146:     y_old = g_new - g_old;
4147: 
4148:     double alpha = 1.0;
4149: 
4150:     double alpha_old = 0.0;
4151: 
4152:     double delta_alpha_old = alpha - alpha_old;
4153: 
4154:     double delta_alpha_new;
4155: 
4156:     unsigned int ls_max = 10;
4157: 
4158:     unsigned int i = 1;
4159: 
4160:     for (; i <= ls_max; ++i)
4161:       {
4162: 	delta_alpha_new = -delta_alpha_old
4163: 	                * (g_new * BFGS_p_vector)/(y_old * BFGS_p_vector);
4164: 	alpha += delta_alpha_new;
4165: 
4166: 	if (std::fabs(delta_alpha_new) < 1.0e-5)
4167: 	  break;
4168: 
4169:         if (i == ls_max)
4170:           {
4171:             alpha = 1.0;
4172:             break;
4173:           }
4174: 
4175:         g_old = g_new;
4176: 
4177:         // BFGS_p_vector is the search direction
4178:         solution_delta_trial = solution_delta;
4179:         solution_delta_trial.add(alpha, BFGS_p_vector);
4180:         update_qph_incremental(solution_delta_trial, m_solution);
4181:         assemble_system_rhs_BFGS_parallel(m_solution, g_new);
4182: 
4183:         y_old = g_new - g_old;
4184: 
4185:         delta_alpha_old = delta_alpha_new;
4186:       }
4187: 
4188:     if (alpha < 1.0e-3)
4189:       alpha = 1.0;
4190: 
4191:     //num_ls = i;
4192:     return alpha;
4193:   }
4194: 
4195:   template <int dim>
4196:   double PhaseFieldMonolithicSolve<dim>::line_search_stepsize_strong_wolfe(const double phi_0,
4197: 				                                           const double phi_0_prime,
4198: 				                                           const BlockVector<double> & BFGS_p_vector,
4199: 				                                           const BlockVector<double> & solution_delta)
4200:   {
4201:     //AssertThrow(phi_0_prime < 0,
4202:     //            ExcMessage("The derivative of phi at alpha = 0 should be negative!"));
4203: 
4204:     // Some line search parameters
4205:     const double c1 = 0.0001;
4206:     const double c2 = 0.9;
4207:     const double alpha_max = 1.0;
4208:     const unsigned int max_iter = 20;
4209:     double alpha = 1.0;
4210: 
4211:     double phi_old = phi_0;
4212:     double phi_prime_old = phi_0_prime;
4213:     double alpha_old = 0.0;
4214: 
4215:     double phi, phi_prime;
4216: 
4217:     std::pair<double, double> current_phi_phi_prime;
4218: 
4219:     unsigned int i = 0;
4220:     for (; i < max_iter; ++i)
4221:       {
4222: 	current_phi_phi_prime = calculate_phi_and_phi_prime(alpha, BFGS_p_vector, solution_delta);
4223: 	phi = current_phi_phi_prime.first;
4224: 	phi_prime = current_phi_phi_prime.second;
4225: 
4226: 	if (   ( phi > (phi_0 + c1 * alpha * phi_0_prime) )
4227: 	    || ( i > 0 && phi > phi_old ) )
4228: 	  {
4229: 	    return line_search_zoom_strong_wolfe(phi_old, phi_prime_old, alpha_old,
4230: 						 phi,     phi_prime,     alpha,
4231: 						 phi_0,   phi_0_prime,   BFGS_p_vector,
4232: 						 c1,      c2,            max_iter, solution_delta);
4233: 	  }
4234: 
4235: 	if (std::fabs(phi_prime) <= c2 * std::fabs(phi_0_prime))
4236: 	  {
4237: 	    return alpha;
4238: 	  }
4239: 
4240: 	if (phi_prime >= 0)
4241: 	  {
4242: 	    return line_search_zoom_strong_wolfe(phi,     phi_prime,     alpha,
4243: 						 phi_old, phi_prime_old, alpha_old,
4244: 						 phi_0,   phi_0_prime,   BFGS_p_vector,
4245: 						 c1,      c2,            max_iter, solution_delta);
4246: 	  }
4247: 
4248: 	phi_old = phi;
4249: 	phi_prime_old = phi_prime;
4250: 	alpha_old = alpha;
4251: 
4252: 	alpha = std::min(0.6*alpha, alpha_max);
4253:       }
4254: 
4255:     return alpha;
4256:   }
4257: 
4258:   template <int dim>
4259:   double PhaseFieldMonolithicSolve<dim>::
4260:     line_search_zoom_strong_wolfe(double phi_low, double phi_low_prime, double alpha_low,
4261: 				  double phi_high, double phi_high_prime, double alpha_high,
4262: 				  double phi_0, double phi_0_prime, const BlockVector<double> & BFGS_p_vector,
4263: 				  double c1, double c2, unsigned int max_iter, const BlockVector<double> & solution_delta)
4264:   {
4265:     double alpha = 0;
4266:     std::pair<double, double> current_phi_phi_prime;
4267:     double phi, phi_prime;
4268: 
4269:     unsigned int i = 0;
4270:     for (; i < max_iter; ++i)
4271:       {
4272: 	// a simple bisection is faster than cubic interpolation
4273: 	alpha = 0.5 * (alpha_low + alpha_high);
4274: 	//alpha = line_search_interpolation_cubic(alpha_low, phi_low, phi_low_prime,
4275: 	//					alpha_high, phi_high, phi_high_prime);
4276: 	current_phi_phi_prime = calculate_phi_and_phi_prime(alpha, BFGS_p_vector, solution_delta);
4277: 	phi = current_phi_phi_prime.first;
4278: 	phi_prime = current_phi_phi_prime.second;
4279: 
4280: 	if (   (phi > phi_0 + c1 * alpha * phi_0_prime)
4281: 	    || (phi > phi_low) )
4282: 	  {
4283: 	    alpha_high = alpha;
4284: 	    phi_high = phi;
4285: 	    phi_high_prime = phi_prime;
4286: 	  }
4287: 	else
4288: 	  {
4289: 	    if (std::fabs(phi_prime) <= c2 * std::fabs(phi_0_prime))
4290: 	      {
4291: 		return alpha;
4292: 	      }
4293: 
4294: 	    if (phi_prime * (alpha_high - alpha_low) >= 0.0)
4295: 	      {
4296: 		alpha_high = alpha_low;
4297: 		phi_high_prime = phi_low_prime;
4298: 		phi_high = phi_low;
4299: 	      }
4300: 
4301: 	    alpha_low = alpha;
4302: 	    phi_low_prime = phi_prime;
4303: 	    phi_low = phi;
4304: 	  }
4305:       }
4306: 
4307:     // avoid unused variable warnings from compiler
4308:     (void)phi_high;
4309:     (void)phi_high_prime;
4310:     return alpha;
4311:   }
4312: 
4313:   template <int dim>
4314:   double PhaseFieldMonolithicSolve<dim>::
4315:     line_search_interpolation_cubic(const double alpha_0, const double phi_0, const double phi_0_prime,
4316:   			            const double alpha_1, const double phi_1, const double phi_1_prime)
4317:   {
4318:     const double d1 = phi_0_prime + phi_1_prime - 3.0 * (phi_0 - phi_1) / (alpha_0 - alpha_1);
4319: 
4320:     const double temp = d1 * d1 - phi_0_prime * phi_1_prime;
4321: 
4322:     if (temp < 0.0)
4323:       return 0.5 * (alpha_0 + alpha_1);
4324: 
4325:     int sign;
4326:     if (alpha_1 > alpha_0)
4327:       sign = 1;
4328:     else
4329:       sign = -1;
4330: 
4331:     const double d2 = sign * std::sqrt(temp);
4332: 
4333:     const double alpha = alpha_1 - (alpha_1 - alpha_0)
4334: 	               * (phi_1_prime + d2 - d1) / (phi_1_prime - phi_0_prime + 2*d2);
4335: 
4336:     if (    (alpha_1 > alpha_0)
4337: 	 && (alpha > alpha_1 || alpha < alpha_0))
4338:       return 0.5 * (alpha_0 + alpha_1);
4339: 
4340:     if (    (alpha_0 > alpha_1)
4341: 	 && (alpha > alpha_0 || alpha < alpha_1))
4342:       return 0.5 * (alpha_0 + alpha_1);
4343: 
4344:     return alpha;
4345:   }
4346: 
4347:   template <int dim>
4348:   std::pair<double, double> PhaseFieldMonolithicSolve<dim>::
4349:     calculate_phi_and_phi_prime(const double alpha,
4350: 				const BlockVector<double> & BFGS_p_vector,
4351: 				const BlockVector<double> & solution_delta)
4352:   {
4353:     // the first component is phi(alpha), the second component is phi_prime(alpha),
4354:     std::pair<double, double> phi_values;
4355: 
4356:     BlockVector<double> solution_delta_trial(solution_delta);
4357:     solution_delta_trial.add(alpha, BFGS_p_vector);
4358: 
4359:     update_qph_incremental(solution_delta_trial, m_solution);
4360: 
4361:     BlockVector<double> system_rhs(m_dofs_per_block);
4362:     assemble_system_rhs_BFGS_parallel(m_solution, system_rhs);
4363:     //m_constraints.condense(system_rhs);
4364: 
4365:     phi_values.first = calculate_energy_functional();
4366:     phi_values.second = system_rhs * BFGS_p_vector;
4367:     return phi_values;
4368:   }
4369: 
4370:   template <int dim>
4371:   void PhaseFieldMonolithicSolve<dim>::LBFGS_B0(BlockVector<double> & LBFGS_r_vector,
4372: 						BlockVector<double> & LBFGS_q_vector)
4373:   {
4374:     m_timer.enter_subsection("Solve B0");
4375: 
4376:     assemble_system_B0(m_solution);
4377: 
4378:     if (m_parameters.m_type_linear_solver == "Direct")
4379:       {
4380: 	SparseDirectUMFPACK A_direct;
4381: 	A_direct.initialize(m_tangent_matrix);
4382: 	A_direct.vmult(LBFGS_r_vector,
4383: 		       LBFGS_q_vector);
4384:       }
4385:     else if (m_parameters.m_type_linear_solver == "CG")
4386:       {
4387: /*
4388: 	SolverControl            solver_control(1e6, 1e-9);
4389: 	SolverCG<BlockVector<double>> cg(solver_control);
4390: 
4391: 	PreconditionJacobi<BlockSparseMatrix<double>> preconditioner;
4392: 	preconditioner.initialize(m_tangent_matrix, 1.0);
4393: 
4394: 	cg.solve(m_tangent_matrix,
4395: 		 LBFGS_r_vector,
4396: 		 LBFGS_q_vector,
4397: 		 preconditioner);
4398: */
4399: 	SolverControl            solver_control_uu(1e6, 1e-9);
4400: 	SolverCG<Vector<double>> cg_uu(solver_control_uu);
4401: 
4402: 	PreconditionJacobi<SparseMatrix<double>> preconditioner_uu;
4403: 	preconditioner_uu.initialize(m_tangent_matrix.block(m_u_dof, m_u_dof), 1.0);
4404: 	cg_uu.solve(m_tangent_matrix.block(m_u_dof, m_u_dof),
4405: 	            LBFGS_r_vector.block(m_u_dof),
4406: 	            LBFGS_q_vector.block(m_u_dof),
4407: 	            preconditioner_uu);
4408: 
4409: 	SolverControl            solver_control_dd(1e6, 1e-15);
4410: 	SolverCG<Vector<double>> cg_dd(solver_control_dd);
4411: 
4412: 	PreconditionJacobi<SparseMatrix<double>> preconditioner_dd;
4413: 	preconditioner_dd.initialize(m_tangent_matrix.block(m_d_dof, m_d_dof), 1.0);
4414: 	cg_dd.solve(m_tangent_matrix.block(m_d_dof, m_d_dof),
4415: 	            LBFGS_r_vector.block(m_d_dof),
4416: 	            LBFGS_q_vector.block(m_d_dof),
4417: 	            preconditioner_dd);
4418:       }
4419:     else
4420:       {
4421: 	AssertThrow(false,
4422: 	            ExcMessage("Selected linear solver not implemented!"));
4423:       }
4424: 
4425:     m_timer.leave_subsection();
4426:   }
4427: 
4428:   template <int dim>
4429:   void PhaseFieldMonolithicSolve<dim>::print_conv_header_LBFGS()
4430:   {
4431:     m_logfile << "\t\t" << "L-BFGS (warning: without phasefield irreversibility)" << std::endl;;
4432:     static const unsigned int l_width = 100;
4433:     m_logfile << '\t' << '\t';
4434:     for (unsigned int i = 0; i < l_width; ++i)
4435:       m_logfile << '_';
4436:     m_logfile << std::endl;
4437: 
4438:     m_logfile << "\t\t itr "
4439:               << " |  LS-alpha     Energy      Res_Norm    "
4440:               << " Res_u      Res_d    Inc_Norm   "
4441:               << " Inc_u      Inc_d" << std::endl;
4442: 
4443:     m_logfile << '\t' << '\t';
4444:     for (unsigned int i = 0; i < l_width; ++i)
4445:       m_logfile << '_';
4446:     m_logfile << std::endl;
4447:   }
4448: 
4449:   template <int dim>
4450:   void PhaseFieldMonolithicSolve<dim>::print_conv_header_LBFGSB()
4451:   {
4452:     m_logfile << '\t' << "L-BFGS-B" << std::endl;;
4453:     static const unsigned int l_width = 130;
4454:     m_logfile << '\t';
4455:     for (unsigned int i = 0; i < l_width; ++i)
4456:       m_logfile << '_';
4457:     m_logfile << std::endl;
4458: 
4459:     m_logfile << "\t itr  "
4460:               << " |    LB    UB    LUB    CG-itr    LS-alpha     Energy      Res_Norm    "
4461:               << " Res_u      Res_d    Inc_Norm   "
4462:               << " Inc_u      Inc_d" << std::endl;
4463: 
4464:     m_logfile << '\t';
4465:     for (unsigned int i = 0; i < l_width; ++i)
4466:       m_logfile << '_';
4467:     m_logfile << std::endl;
4468:   }
4469: 
4470:   template <int dim>
4471:   void PhaseFieldMonolithicSolve<dim>::
4472:   solve_nonlinear_timestep_LBFGS(BlockVector<double> & solution_delta,
4473: 				 BlockVector<double> & LBFGS_update_refine)
4474:   {
4475:     BlockVector<double> LBFGS_update(m_dofs_per_block);
4476: 
4477:     LBFGS_update = 0.0;
4478: 
4479:     m_error_residual.reset();
4480:     m_error_residual_0.reset();
4481:     m_error_residual_norm.reset();
4482:     m_error_update.reset();
4483:     m_error_update_0.reset();
4484:     m_error_update_norm.reset();
4485: 
4486:     if (m_parameters.m_output_iteration_history)
4487:       print_conv_header_LBFGS();
4488: 
4489:     unsigned int LBFGS_iteration = 0;
4490: 
4491:     BlockVector<double> LBFGS_r_vector(m_dofs_per_block);
4492:     BlockVector<double> LBFGS_y_vector(m_dofs_per_block);
4493:     BlockVector<double> LBFGS_q_vector(m_dofs_per_block);
4494:     BlockVector<double> LBFGS_s_vector(m_dofs_per_block);
4495:     std::list<std::pair< std::pair<BlockVector<double>,
4496:                                    BlockVector<double>>,
4497:                          double>> LBFGS_vector_list;
4498: 
4499:     const unsigned int LBFGS_m = m_parameters.m_LBFGS_m;
4500:     std::list<double> LBFGS_alpha_list;
4501: 
4502:     double line_search_parameter = 0.0;
4503:     double LBFGS_beta = 0.0;
4504:     double rho = 0.0;
4505: 
4506:     for (; LBFGS_iteration < m_parameters.m_max_iterations_BFGS; ++LBFGS_iteration)
4507:       {
4508: 	if (m_parameters.m_output_iteration_history)
4509: 	  m_logfile << '\t' << '\t' << std::setw(4) << LBFGS_iteration << ' '
4510:                     << std::flush;
4511: 
4512:         make_constraints(LBFGS_iteration);
4513: 
4514:         // At the first step, we simply distribute the inhomogeneous part of
4515:         // the constraints
4516:         if (LBFGS_iteration == 0)
4517:           {
4518:             // use the solution from the previous solve on the
4519:             // refined mesh as initial guess
4520:             LBFGS_update = LBFGS_update_refine;
4521: 
4522:             m_constraints.distribute(LBFGS_update);
4523:             solution_delta += LBFGS_update;
4524: 
4525:             update_qph_incremental(solution_delta, m_solution);
4526:             if (m_parameters.m_output_iteration_history)
4527:               {
4528:                 m_logfile << " | " << std::flush;
4529:                 m_logfile << std::endl;
4530:               }
4531:             continue;
4532:           }
4533:         else if (LBFGS_iteration == 1)
4534:           {
4535: 	    // Calculate the residual vector r. NOTICE that in the context of
4536: 	    // BFGS, this r is the gradient of the energy functional (objective function),
4537: 	    // NOT the negative gradient of the energy functional
4538: 	    assemble_system_rhs_BFGS_parallel(m_solution, m_system_rhs);
4539: 
4540: 	    // We cannot simply zero out the dofs that are constrained, since we might
4541: 	    // have hanging node constraints. In this case, we need to modify the RHS
4542: 	    // as C^T * b, which C contains entries of 0.5 (x_3 = 0.5*x_1 + 0.5*x_2)
4543: 	    //for (unsigned int i = 0; i < m_dof_handler.n_dofs(); ++i)
4544: 	      //if (m_constraints.is_constrained(i))
4545: 		//m_system_rhs(i) = 0.0;
4546: 
4547: 	    // if m_constraints has inhomogeneity, we cannot call m_constraints.condense(m_system_rhs),
4548: 	    // since the m_system_matrix needs to be provided to modify the RHS properly. However, this
4549: 	    // error will not be detected in the release mode and only will be detected on the debug mode
4550: 	    // if we use assemble_system_rhs_BFGS_parallel, then condense() is not necessary
4551: 	    //m_constraints.condense(m_system_rhs);
4552:           }
4553: 
4554:         get_error_residual(m_error_residual);
4555:         if (LBFGS_iteration == 1)
4556:           m_error_residual_0 = m_error_residual;
4557: 
4558:         m_error_residual_norm = m_error_residual;
4559:         // For three-point bending problem and 3D problem, we use absolute residual
4560:         // for convergence test
4561:         if (m_parameters.m_relative_residual)
4562:           m_error_residual_norm.normalize(m_error_residual_0);
4563: 
4564:         if (LBFGS_iteration > 1 && m_error_update_norm.m_u <= m_parameters.m_tol_u_incr
4565:                                 && m_error_residual_norm.m_u <= m_parameters.m_tol_u_residual
4566: 			        && m_error_update_norm.m_d <= m_parameters.m_tol_d_incr
4567: 			        && m_error_residual_norm.m_d <= m_parameters.m_tol_d_residual
4568: 				)
4569:           {
4570:             if (m_parameters.m_output_iteration_history)
4571:               {
4572: 		m_logfile << " | ";
4573: 		m_logfile << " CONVERGED! " << std::fixed << std::setprecision(3) << std::setw(7)
4574: 			  << std::scientific
4575: 		      << "    ----    "
4576: 		      << "  " << m_error_residual_norm.m_norm
4577: 		      << "  " << m_error_residual_norm.m_u
4578: 		      << "  " << m_error_residual_norm.m_d
4579: 		      << "  " << m_error_update_norm.m_norm
4580: 		      << "  " << m_error_update_norm.m_u
4581: 		      << "  " << m_error_update_norm.m_d
4582: 		      << "  " << std::endl;
4583: 
4584: 		m_logfile << '\t' << '\t';
4585: 		for (unsigned int i = 0; i < 100; ++i)
4586: 		  m_logfile << '_';
4587: 		m_logfile << std::endl;
4588:               }
4589: 
4590:             m_logfile << "\t\tConvergence is reached after "
4591:         	      << LBFGS_iteration << " L-BFGS iterations."<< std::endl;
4592: 
4593:             m_logfile << "\t\tResidual information of convergence:" << std::endl;
4594: 
4595:             if (m_parameters.m_relative_residual)
4596:               {
4597: 		m_logfile << "\t\t\tRelative residual of disp. equation: "
4598: 			  << m_error_residual_norm.m_u << std::endl;
4599: 
4600: 		m_logfile << "\t\t\tAbsolute residual of disp. equation: "
4601: 			  << m_error_residual_norm.m_u * m_error_residual_0.m_u << std::endl;
4602: 
4603: 		m_logfile << "\t\t\tRelative residual of phasefield equation: "
4604: 			  << m_error_residual_norm.m_d << std::endl;
4605: 
4606: 		m_logfile << "\t\t\tAbsolute residual of phasefield equation: "
4607: 			  << m_error_residual_norm.m_d * m_error_residual_0.m_d << std::endl;
4608: 
4609: 		m_logfile << "\t\t\tRelative increment of disp.: "
4610: 			  << m_error_update_norm.m_u << std::endl;
4611: 
4612: 		m_logfile << "\t\t\tAbsolute increment of disp.: "
4613: 			  << m_error_update_norm.m_u * m_error_update_0.m_u << std::endl;
4614: 
4615: 		m_logfile << "\t\t\tRelative increment of phasefield: "
4616: 			  << m_error_update_norm.m_d << std::endl;
4617: 
4618: 		m_logfile << "\t\t\tAbsolute increment of phasefield: "
4619: 			  << m_error_update_norm.m_d * m_error_update_0.m_d << std::endl;
4620:               }
4621:             else
4622:               {
4623: 		m_logfile << "\t\t\tAbsolute residual of disp. equation: "
4624: 			  << m_error_residual_norm.m_u << std::endl;
4625: 
4626: 		m_logfile << "\t\t\tAbsolute residual of phasefield equation: "
4627: 			  << m_error_residual_norm.m_d << std::endl;
4628: 
4629: 		m_logfile << "\t\t\tAbsolute increment of disp.: "
4630: 			  << m_error_update_norm.m_u << std::endl;
4631: 
4632: 		m_logfile << "\t\t\tAbsolute increment of phasefield: "
4633: 			  << m_error_update_norm.m_d << std::endl;
4634:               }
4635: 
4636:             break;
4637:           }
4638: 
4639:         // LBFGS algorithm
4640:         LBFGS_q_vector = m_system_rhs;
4641: 
4642:         LBFGS_alpha_list.clear();
4643:         for (auto itr = LBFGS_vector_list.begin(); itr != LBFGS_vector_list.end(); ++itr)
4644:           {
4645:             LBFGS_s_vector = (itr->first).first;
4646:             LBFGS_y_vector = (itr->first).second;
4647:             rho = itr->second;
4648: 
4649:             const double alpha = rho * (LBFGS_s_vector * LBFGS_q_vector);
4650:             LBFGS_alpha_list.push_back(alpha);
4651: 
4652:             LBFGS_q_vector.add(-alpha, LBFGS_y_vector);
4653:           }
4654: /*
4655:         double scale_gamma = 0.0;
4656:         if (LBFGS_iteration == 1)
4657:           {
4658:             scale_gamma = 1.0;
4659:           }
4660:         else
4661:           {
4662:             LBFGS_s_vector = LBFGS_vector_list.front().first.first;
4663:             LBFGS_y_vector = LBFGS_vector_list.front().first.second;
4664:             scale_gamma = (LBFGS_s_vector * LBFGS_y_vector)/(LBFGS_y_vector * LBFGS_y_vector);
4665:           }
4666: 
4667:         LBFGS_q_vector *= scale_gamma;
4668:         LBFGS_r_vector = LBFGS_q_vector;
4669: */
4670:         LBFGS_B0(LBFGS_r_vector,
4671: 		 LBFGS_q_vector);
4672: 
4673:         for (auto itr = LBFGS_vector_list.rbegin(); itr != LBFGS_vector_list.rend(); ++itr)
4674:           {
4675:             LBFGS_s_vector = (itr->first).first;
4676:             LBFGS_y_vector = (itr->first).second;
4677:             rho = itr->second;
4678: 
4679:             LBFGS_beta = rho * (LBFGS_y_vector * LBFGS_r_vector);
4680: 
4681:             const double alpha = LBFGS_alpha_list.back();
4682:             LBFGS_alpha_list.pop_back();
4683: 
4684:             LBFGS_r_vector.add(alpha - LBFGS_beta, LBFGS_s_vector);
4685:           }
4686: 
4687:         LBFGS_r_vector *= -1.0; // this is the p_vector (search direction)
4688: 
4689:         m_constraints.distribute(LBFGS_r_vector);
4690: 
4691:         // We need a line search algorithm to decide line_search_parameter
4692: 
4693:         if(m_parameters.m_type_line_search == "StrongWolfe")
4694:           {
4695:             const double phi_0 = calculate_energy_functional();
4696:             const double phi_0_prime = m_system_rhs * LBFGS_r_vector;
4697: 
4698:             line_search_parameter = line_search_stepsize_strong_wolfe(phi_0,
4699: 		    				                      phi_0_prime,
4700: 								      LBFGS_r_vector,
4701: 						                      solution_delta);
4702:           }
4703:         else if(m_parameters.m_type_line_search == "GradientBased")
4704:           {
4705: 	    // LBFGS_r_vector is the search direction
4706: 	    line_search_parameter = line_search_stepsize_gradient_based(LBFGS_r_vector,
4707: 									solution_delta);
4708:           }
4709:         else
4710:           {
4711:             Assert(false, ExcMessage("An unknown line search method is called!"));
4712:           }
4713: 
4714:         LBFGS_r_vector *= line_search_parameter;
4715:         LBFGS_update = LBFGS_r_vector;
4716: 
4717:         get_error_update(LBFGS_update, m_error_update);
4718:         if (LBFGS_iteration == 1)
4719:           m_error_update_0 = m_error_update;
4720: 
4721:         m_error_update_norm = m_error_update;
4722:         // For three-point bending problem and the sphere inclusion problem,
4723:         // we use absolute residual for convergence test
4724:         if (m_parameters.m_relative_residual)
4725:           m_error_update_norm.normalize(m_error_update_0);
4726: 
4727:         solution_delta += LBFGS_update;
4728:         update_qph_incremental(solution_delta, m_solution);
4729: 
4730:         LBFGS_y_vector = m_system_rhs;
4731:         LBFGS_y_vector *= -1.0;
4732:         assemble_system_rhs_BFGS_parallel(m_solution, m_system_rhs);
4733:         // if we use assemble_system_rhs_BFGS_parallel, then condense() is not necessary
4734:         //m_constraints.condense(m_system_rhs);
4735:         LBFGS_y_vector += m_system_rhs;
4736: 
4737:         LBFGS_s_vector = LBFGS_update;
4738: 
4739:         if (LBFGS_iteration > LBFGS_m)
4740:           LBFGS_vector_list.pop_back();
4741: 
4742:         rho = 1.0 / (LBFGS_y_vector * LBFGS_s_vector);
4743: 
4744:         LBFGS_vector_list.push_front(std::make_pair(std::make_pair(LBFGS_s_vector,
4745: 								   LBFGS_y_vector),
4746: 						    rho));
4747:         if (m_parameters.m_output_iteration_history)
4748:           {
4749: 	    const double energy_functional = calculate_energy_functional();
4750: 
4751: 	    m_logfile << " | " << std::fixed << std::setprecision(3) << std::setw(1)
4752: 		      << std::scientific
4753: 		      << "" << line_search_parameter
4754: 		      << std::fixed << std::setprecision(6) << std::setw(1)
4755: 					<< std::scientific
4756: 		      << "  " << energy_functional
4757: 		      << std::fixed << std::setprecision(3) << std::setw(1)
4758: 					<< std::scientific
4759: 		      << "  " << m_error_residual_norm.m_norm
4760: 		      << "  " << m_error_residual_norm.m_u
4761: 		      << "  " << m_error_residual_norm.m_d
4762: 		      << "  " << m_error_update_norm.m_norm
4763: 		      << "  " << m_error_update_norm.m_u
4764: 		      << "  " << m_error_update_norm.m_d
4765: 		      << "  " << std::endl;
4766:           }
4767:       }
4768: 
4769:     AssertThrow(LBFGS_iteration < m_parameters.m_max_iterations_BFGS,
4770:                 ExcMessage("No convergence in L-BFGS nonlinear solver!"));
4771:   }
4772: 
4773:   template <int dim>
4774:   void PhaseFieldMonolithicSolve<dim>::
4775:   calculate_cauchy_point(const BlockSparseMatrix<double> & B0_matrix,
4776: 	                 const std::list<BlockVector<double>> & y_vector_list,
4777: 		         const std::list<BlockVector<double>> & b0xs_vector_list,
4778: 			 const FullMatrix<double> & M_matrix,
4779: 			 const BlockVector<double> & gradient_g,
4780: 			 const BlockVector<double> & solution_delta,
4781: 			 BlockVector<double> & solution_delta_cauchy_point)
4782:   {
4783:     m_timer.enter_subsection("Calculate Cauchy point");
4784: 
4785:     solution_delta_cauchy_point = 0.0;
4786:     BlockVector<double> gradient_d(gradient_g);
4787:     gradient_d *= -1;
4788: 
4789:     const unsigned int list_size = y_vector_list.size();
4790:     const auto itr_y_begin    = y_vector_list.begin();
4791:     const auto itr_b0xs_begin = b0xs_vector_list.begin();
4792: 
4793:     // t_series only contains t > 0
4794:     std::priority_queue< std::pair<double, unsigned int>,
4795:                          std::vector<std::pair<double, unsigned int>>,
4796:         		 std::greater<std::pair<double, unsigned int>> >
4797:     t_series = calculate_break_points(solution_delta,
4798:     			              gradient_g,
4799: 				      gradient_d);
4800: 
4801:     // m_active_set_phasefield contains 1 or 2 for active set and 0 for inactive set
4802:     for (unsigned int i = 0; i < m_active_set_phasefield.size(); ++i)
4803:       {
4804: 	if (m_active_set_phasefield(i) > 0.5)
4805: 	  solution_delta_cauchy_point.block(m_d_dof)[i] = solution_delta.block(m_d_dof)[i];
4806:       }
4807: 
4808:     // p = W^T * d
4809:     Vector<double> p(2 * list_size);
4810:     for (unsigned int i = 0; i < list_size; ++i)
4811:       {
4812:         p(i)             = (*std::next(itr_y_begin,    i)) * gradient_d;
4813:         p(i + list_size) = (*std::next(itr_b0xs_begin, i)) * gradient_d;
4814:       }
4815: 
4816:     Vector<double> c(2 * list_size);
4817:     c = 0.0;
4818: 
4819:     double f_prime = -(gradient_d * gradient_d);
4820: 
4821:     // M * p
4822:     Vector<double> Mp(2 * list_size);
4823:     if (list_size > 0)
4824:       M_matrix.vmult(Mp, p);
4825: 
4826:     // B_0 * d
4827:     BlockVector<double> B0_grandient_d(m_dofs_per_block);
4828:     B0_matrix.vmult(B0_grandient_d, gradient_d);
4829: 
4830:     double f_prime_prime = gradient_d * B0_grandient_d;
4831:     if (list_size > 0)
4832:       f_prime_prime -= (p * Mp);
4833: 
4834:     double delta_t_min = -f_prime / f_prime_prime;
4835: 
4836:     double t_old = 0.0;
4837: 
4838:     std::pair<double, unsigned int> top_pair = t_series.top();
4839:     double t = top_pair.first;
4840:     unsigned int b = top_pair.second;
4841: 
4842:     double delta_t = t - t_old;
4843: 
4844:     BlockVector<double> z(m_dofs_per_block);
4845:     z = 0.0;
4846: 
4847:     // w_b = W^T * e_b
4848:     Vector<double> w_b(2 * list_size);
4849: 
4850:     // w_b_T_x_M = w_b^T * M
4851:     Vector<double> w_b_T_x_M(2 * list_size);
4852: 
4853:     Vector<double> temp_vector(m_dofs_per_block[m_d_dof]);
4854: 
4855:     while (delta_t_min >= delta_t)
4856:       {
4857: 	t_series.pop();
4858: 
4859: 	if (gradient_d.block(m_d_dof)[b] > 0)
4860: 	  solution_delta_cauchy_point.block(m_d_dof)[b] = 1.0 - m_solution.block(m_d_dof)[b];
4861: 	else if (gradient_d.block(m_d_dof)[b] < 0)
4862: 	  solution_delta_cauchy_point.block(m_d_dof)[b] = 0.0;
4863: 	else
4864: 	  AssertThrow(false,
4865: 	              ExcMessage("gradient_d(b) cannot be zero!"));
4866: 
4867: 	if (gradient_d.block(m_d_dof)[b] < 0)
4868: 	  m_active_set_phasefield[b] = 1; //lower bound
4869: 	else
4870: 	  m_active_set_phasefield[b] = 2; //upper bound
4871: 
4872:         // z_{j} = z_{j-1} + delta_t_{j-1} * gradient_d_{j-1};
4873: 	z.sadd(1.0, delta_t, gradient_d);
4874: 
4875: 	// c_{j} = c_{j-1} + delta_t_{j-1} * p_{j-1};
4876: 	if (list_size > 0)
4877: 	  c.sadd(1.0, delta_t, p);
4878: 
4879:         double temp_scalar = ebT_x_B0_x_v(b, B0_matrix, z);
4880: 
4881:         // w_b = W^T * e_b
4882:         for (unsigned int i = 0; i < list_size; ++i)
4883:           {
4884:             w_b(i)             = (*std::next(itr_y_begin,    i)).block(m_d_dof)[b];
4885:             w_b(i + list_size) = (*std::next(itr_b0xs_begin, i)).block(m_d_dof)[b];
4886:           }
4887: 
4888:         if (list_size > 0)
4889:           M_matrix.vmult(w_b_T_x_M, w_b);
4890: 
4891: 	f_prime += delta_t * f_prime_prime
4892: 	         + gradient_g.block(m_d_dof)[b] * gradient_g.block(m_d_dof)[b]
4893: 	         + temp_scalar * gradient_g.block(m_d_dof)[b];
4894: 
4895: 	if (list_size > 0)
4896: 	  f_prime -= (w_b_T_x_M * c) * gradient_g.block(m_d_dof)[b];
4897: 
4898: 	temp_scalar = ebT_x_B0_x_v(b, B0_matrix, gradient_d);
4899: 
4900: 	f_prime_prime += 2.0 * gradient_g.block(m_d_dof)[b] * temp_scalar
4901: 	               + gradient_g.block(m_d_dof)[b] * gradient_g.block(m_d_dof)[b]
4902: 		             * B0_matrix.block(m_d_dof, m_d_dof)(b, b);
4903: 
4904: 	if (list_size > 0)
4905: 	  {
4906: 	    f_prime_prime -= 2.0 * gradient_g.block(m_d_dof)[b] * (w_b_T_x_M * p);
4907: 	    f_prime_prime -= gradient_g.block(m_d_dof)[b] * gradient_g.block(m_d_dof)[b] * (w_b_T_x_M * w_b);
4908: 	  }
4909: 
4910: 	// p_{j} = p_{j-1} + g_b * w_b;
4911: 	if (list_size > 0)
4912: 	  p.sadd(1.0, gradient_g.block(m_d_dof)[b], w_b);
4913: 
4914: 	gradient_d.block(m_d_dof)[b] = 0.0;
4915: 
4916: 	delta_t_min = -f_prime / f_prime_prime;
4917: 
4918: 	t_old = t;
4919: 
4920: 	top_pair = t_series.top();
4921: 	t = top_pair.first;
4922: 	b = top_pair.second;
4923: 
4924: 	delta_t = t - t_old;
4925:       }
4926: 
4927:     if (delta_t_min < 0)
4928:       delta_t_min = 0;
4929: 
4930:     t_old += delta_t_min;
4931: 
4932:     for(unsigned int i = 0; i < m_active_set_phasefield.size(); ++i)
4933:       {
4934: 	// inactive phasefield dof
4935: 	if (m_active_set_phasefield(i) < 0.5)
4936: 	  solution_delta_cauchy_point.block(m_d_dof)[i] = solution_delta.block(m_d_dof)[i]
4937: 						+ t_old * gradient_d.block(m_d_dof)[i];
4938:       }
4939: 
4940:     // There are no active constraints in the displacement field
4941:     solution_delta_cauchy_point.block(m_u_dof) = solution_delta.block(m_u_dof);
4942:     (solution_delta_cauchy_point.block(m_u_dof)).add(t_old, gradient_d.block(m_u_dof));
4943: 
4944:     // We need to make sure the solution_delta_cauchy_point satisfies the essential
4945:     // boundary conditions and the hanging-node constraints
4946:     m_constraints.distribute(solution_delta_cauchy_point);
4947: 
4948:     m_timer.leave_subsection();
4949:   }
4950: 
4951:   template <int dim>
4952:   void PhaseFieldMonolithicSolve<dim>::
4953:   solve_nonlinear_timestep_LBFGS_B(BlockVector<double> & solution_delta,
4954: 				   BlockVector<double> & LBFGS_update_refine)
4955:   {
4956:     BlockVector<double> LBFGS_update(m_dofs_per_block);
4957:     BlockVector<double> solution_delta_cauchy_point(m_dofs_per_block);
4958:     LBFGS_update = 0.0;
4959: 
4960:     const unsigned int LBFGS_m = m_parameters.m_LBFGS_m;
4961: 
4962:     unsigned int LBFGS_iteration = 0;
4963: 
4964:     m_error_residual.reset();
4965:     m_error_residual_0.reset();
4966:     m_error_residual_norm.reset();
4967:     m_error_update.reset();
4968:     m_error_update_0.reset();
4969:     m_error_update_norm.reset();
4970: 
4971:     if (m_parameters.m_output_iteration_history)
4972:       print_conv_header_LBFGSB();
4973: 
4974:     BlockVector<double> LBFGS_s_vector(m_dofs_per_block);
4975:     BlockVector<double> LBFGS_y_vector(m_dofs_per_block);
4976:     BlockVector<double> free_dofs(m_dofs_per_block);
4977:     BlockVector<double> b0xs_vector(m_dofs_per_block);
4978: 
4979:     // all the list goes from k-m to k-1
4980:     // the front is the oldest quantity,and the end is
4981:     // newest quantity
4982:     std::list<BlockVector<double>> s_vector_list;
4983:     std::list<BlockVector<double>> y_vector_list;
4984:     std::list<double> s_dot_y_list;
4985:     std::list<BlockVector<double>> b0xs_vector_list;
4986: 
4987:     double line_search_parameter = 0.0;
4988: 
4989:     unsigned int lower_bound_number_old = 0;
4990:     unsigned int upper_bound_number_old = 0;
4991:     unsigned int lowerupper_bound_number_old = 0;
4992: 
4993:     unsigned int lower_bound_number_new = 0;
4994:     unsigned int upper_bound_number_new = 0;
4995:     unsigned int lowerupper_bound_number_new = 0;
4996: 
4997:     for (; LBFGS_iteration < m_parameters.m_max_iterations_BFGS; ++LBFGS_iteration)
4998:       {
4999: 	if (m_parameters.m_output_iteration_history)
5000: 	  m_logfile << '\t' << std::setw(4) << LBFGS_iteration << ' '
5001:                     << std::flush;
5002: 
5003:         make_constraints(LBFGS_iteration);
5004: 
5005:         // At the first step, we simply distribute the inhomogeneous part of
5006:         // the constraints
5007:         if (LBFGS_iteration == 0)
5008:           {
5009:             // use the solution from the previous solve on the
5010:             // refined mesh as initial guess
5011:             LBFGS_update = LBFGS_update_refine;
5012: 
5013:             m_constraints.distribute(LBFGS_update);
5014:             solution_delta += LBFGS_update;
5015:             update_qph_incremental(solution_delta, m_solution);
5016:             assemble_system_rhs_BFGS_parallel(m_solution, m_system_rhs);
5017:             m_logfile << "  | " << std::endl;
5018: 
5019:             continue;
5020:           }
5021: 
5022:         get_error_residual_LBFGSB(m_error_residual,
5023: 				  solution_delta);
5024: 
5025:         if (LBFGS_iteration == 1)
5026:           m_error_residual_0 = m_error_residual;
5027: 
5028:         m_error_residual_norm = m_error_residual;
5029: 
5030:         if (m_parameters.m_relative_residual)
5031:           m_error_residual_norm.normalize(m_error_residual_0);
5032: 
5033:         if (LBFGS_iteration > 1 && m_error_update_norm.m_u <= m_parameters.m_tol_u_incr
5034:                                 && m_error_residual_norm.m_u <= m_parameters.m_tol_u_residual
5035: 			        && m_error_update_norm.m_d <= m_parameters.m_tol_d_incr
5036: 			        && m_error_residual_norm.m_d <= m_parameters.m_tol_d_residual
5037: 				&& lower_bound_number_new == lower_bound_number_old
5038: 				&& upper_bound_number_new == upper_bound_number_old
5039: 				&& lowerupper_bound_number_new == lowerupper_bound_number_old)
5040:           {
5041:             if (m_parameters.m_output_iteration_history)
5042:               {
5043: 		m_logfile << "  | ";
5044: 		m_logfile << " CONVERGED! " << std::fixed << std::setprecision(3) << std::setw(7)
5045: 			  << std::scientific
5046: 		      << "           ---      "
5047: 		      << "\t\t\t"
5048: 		      << "  " << m_error_residual_norm.m_norm
5049: 		      << "  " << m_error_residual_norm.m_u
5050: 		      << "  " << m_error_residual_norm.m_d
5051: 		      << "  " << m_error_update_norm.m_norm
5052: 		      << "  " << m_error_update_norm.m_u
5053: 		      << "  " << m_error_update_norm.m_d
5054: 		      << "  " << std::endl;
5055: 
5056: 		m_logfile << '\t';
5057: 		for (unsigned int i = 0; i < 130; ++i)
5058: 		  m_logfile << '_';
5059: 		m_logfile << std::endl;
5060:               }
5061: 
5062:             m_logfile << "\t\tThe current L-BFGS-B step converges in "
5063:         	      << LBFGS_iteration
5064:         	      << " iterations." << std::endl;
5065:             m_logfile << "\t\tNumber of active lower bounds not changed." << std::endl;
5066:     	    m_logfile << "\t\tNumber of active upper bounds not changed." << std::endl;
5067:     	    m_logfile << "\t\tNumber of active lower-upper bounds not changed." << std::endl;
5068: 
5069:             if (m_parameters.m_relative_residual)
5070:               {
5071: 		m_logfile << "\t\tProjected gradient of disp. (relative): "
5072: 			  << m_error_residual_norm.m_u << std::endl;
5073: 
5074: 		m_logfile << "\t\tProjected gradient of disp. (absolute): "
5075: 			  << m_error_residual_norm.m_u * m_error_residual_0.m_u << std::endl;
5076: 
5077: 		m_logfile << "\t\tProjected gradient of phasefield (relative): "
5078: 			  << m_error_residual_norm.m_d << std::endl;
5079: 
5080: 		m_logfile << "\t\tProjected gradient of phasefield (absolute): "
5081: 			  << m_error_residual_norm.m_d * m_error_residual_0.m_d << std::endl;
5082: 
5083: 		m_logfile << "\t\tRelative increment of disp.: "
5084: 			  << m_error_update_norm.m_u << std::endl;
5085: 
5086: 		m_logfile << "\t\tAbsolute increment of disp.: "
5087: 			  << m_error_update_norm.m_u * m_error_update_0.m_u << std::endl;
5088: 
5089: 		m_logfile << "\t\tRelative increment of phasefield: "
5090: 			  << m_error_update_norm.m_d << std::endl;
5091: 
5092: 		m_logfile << "\t\tAbsolute increment of phasefield: "
5093: 			  << m_error_update_norm.m_d * m_error_update_0.m_d << std::endl;
5094:               }
5095:             else
5096:               {
5097: 		m_logfile << "\t\tProjected gradient of disp. (absolute): "
5098: 			  << m_error_residual_norm.m_u << std::endl;
5099: 
5100: 		m_logfile << "\t\tProjected gradient of phasefield (absolute): "
5101: 			  << m_error_residual_norm.m_d << std::endl;
5102: 
5103: 		m_logfile << "\t\tAbsolute increment of disp.: "
5104: 			  << m_error_update_norm.m_u << std::endl;
5105: 
5106: 		m_logfile << "\t\tAbsolute increment of phasefield: "
5107: 			  << m_error_update_norm.m_d << std::endl;
5108:               }
5109: 
5110:             break;
5111:           }
5112: 
5113: 	// assemble the initial B_0 matrix at the k-th L-BFGS iteration
5114: 	// m_solution is the old solution from the previous converged step
5115: 	// it is needed only for the viscosity term
5116: 	// the output is m_tangent_matrix (B^0_k)
5117: 	assemble_system_B0(m_solution);
5118: 
5119: 	// B^0_k * s_vector has to be completely recalculated from scratch
5120: 	// at each L-BFGS iteration, since B^0_k is different
5121: 	b0xs_vector_list.clear();
5122: 	for (auto itr = s_vector_list.begin(); itr != s_vector_list.end(); ++itr)
5123: 	  {
5124: 	    m_tangent_matrix.vmult(b0xs_vector, *itr);
5125: 	    b0xs_vector_list.push_back(b0xs_vector);
5126: 	  }
5127: 
5128: 	// In the iteration LBFGS_iteration = 0, only the essential boundary conditions
5129: 	// are applied.
5130: 	// WHen LBFGS_iteration = 1, it is the first step of LBFGS update, and the
5131: 	// s_vector_list and y_vector_list are empty.
5132: 	// Since the pair of s and y will only be added to the list if s dot y > tol,
5133: 	// it is safer to decide the matrix dimension by the size of the list.
5134: 	const unsigned int list_size = s_vector_list.size();
5135: 	const auto itr_s_begin    = s_vector_list.begin();
5136: 	const auto itr_y_begin    = y_vector_list.begin();
5137: 	const auto itr_b0xs_begin = b0xs_vector_list.begin();
5138: 	const auto itr_s_dot_y_begin = s_dot_y_list.begin();
5139: 
5140: 	FullMatrix<double> sTxBxs_matrix(list_size);
5141: 	sTxBxs_matrix = 0;
5142: 	for (unsigned int i = 0; i < list_size; ++i)
5143: 	  for (unsigned int j = 0; j <= i; ++j)
5144: 	    {
5145: 	      sTxBxs_matrix(i, j) = (*std::next(itr_s_begin,    i))
5146: 		                  * (*std::next(itr_b0xs_begin, j));
5147: 	    }
5148: 	for (unsigned int i = 0; i < list_size; ++i)
5149: 	  for (unsigned int j = i + 1; j < list_size; ++j)
5150: 	    {
5151: 	      sTxBxs_matrix(i, j) = sTxBxs_matrix(j, i);
5152: 	    }
5153: 
5154: 	FullMatrix<double> D_matrix(list_size);
5155: 	D_matrix = 0;
5156: 	for (unsigned int i = 0; i < list_size; ++i)
5157: 	  D_matrix(i, i) = (*std::next(itr_s_dot_y_begin, i));
5158: 
5159: 	FullMatrix<double> L_matrix(list_size);
5160: 	L_matrix = 0;
5161: 	for (unsigned int i = 0; i < list_size; ++i)
5162: 	  for (unsigned int j = 0; j < i; ++j)
5163: 	    L_matrix(i, j) = (*std::next(itr_s_begin, i))
5164:                            * (*std::next(itr_y_begin, j));
5165: 
5166: 	FullMatrix<double> M_matrix_inv(2 * list_size);
5167: 	FullMatrix<double> M_matrix(2 * list_size);
5168: 
5169: 	M_matrix_inv = 0;
5170: 	for (unsigned int i = 0; i < list_size; ++i)
5171: 	  M_matrix_inv(i, i) = -D_matrix(i, i);
5172: 
5173: 	for (unsigned int i = 0; i < list_size; ++i)
5174:           for (unsigned int j = 0; j < list_size; ++j)
5175:             {
5176:               M_matrix_inv(i + list_size, j + list_size) = sTxBxs_matrix(i, j);
5177:               M_matrix_inv(i + list_size, j            ) = L_matrix(i, j);
5178:               M_matrix_inv(i            , j + list_size) = L_matrix(j, i);
5179:             }
5180: 
5181: 	if (!M_matrix_inv.empty())
5182: 	  M_matrix.invert(M_matrix_inv);
5183: 
5184: 	m_active_set_phasefield = 0;
5185: 	calculate_cauchy_point(m_tangent_matrix,
5186: 			       y_vector_list,
5187: 			       b0xs_vector_list,
5188: 			       M_matrix,
5189: 			       m_system_rhs,
5190: 			       solution_delta,
5191: 			       solution_delta_cauchy_point);
5192: 
5193: 	// We need to find out which DOFs are free:
5194: 	// no essential boundary conditions, no hanging node constraints
5195: 	// no active box constraints
5196: 	unsigned int free_disp_number = 0;
5197: 	for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
5198: 	  {
5199: 	    if (m_constraints.is_constrained(i))
5200: 	      free_dofs.block(m_u_dof)[i] = -1;
5201: 	    else
5202: 	      {
5203: 	        free_dofs.block(m_u_dof)[i] = 1;
5204: 	        ++free_disp_number;
5205: 	      }
5206: 	  }
5207: 
5208: 	unsigned int free_phasefield_number = 0;
5209: 	for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
5210: 	  {
5211: 	    if (   m_constraints.is_constrained(i + m_dofs_per_block[m_u_dof])
5212: 		|| m_active_set_phasefield(i) > 0.5)
5213: 	      free_dofs.block(m_d_dof)[i] = -1;
5214: 	    else
5215: 	      {
5216: 	        free_dofs.block(m_d_dof)[i] = 1;
5217: 	        ++free_phasefield_number;
5218: 	      }
5219: 	  }
5220: 
5221: 	// temp_vector_1 = x^c - x_k
5222: 	BlockVector<double> temp_vector_1(solution_delta_cauchy_point);
5223: 	temp_vector_1 -= solution_delta;
5224: 
5225: 	// temp_vector_2 = B_0 * (x^c - x_k)
5226: 	BlockVector<double> temp_vector_2(m_dofs_per_block);
5227: 	m_tangent_matrix.vmult(temp_vector_2, temp_vector_1);
5228: 
5229: 	// temp_vector_3 = W^T * (x^c - x_k)
5230: 	Vector<double> temp_vector_3(2 * list_size);
5231: 	for (unsigned int i = 0; i < list_size; ++i)
5232: 	  {
5233: 	    temp_vector_3(i)             = (*std::next(itr_y_begin,    i)) * temp_vector_1;
5234: 	    temp_vector_3(i + list_size) = (*std::next(itr_b0xs_begin, i)) * temp_vector_1;
5235: 	  }
5236: 
5237: 	// temp_vector_4 = M * W^T * (x^c - x_k)
5238: 	Vector<double> temp_vector_4(2 * list_size);
5239: 	if (list_size > 0)
5240: 	  M_matrix.vmult(temp_vector_4, temp_vector_3);
5241: 
5242: 	// temp_vector_5 = W * M * W^T * (x^c - x_k)
5243: 	BlockVector<double> temp_vector_5(m_dofs_per_block);
5244: 	for (unsigned int i = 0; i < list_size; ++i)
5245: 	  {
5246: 	    temp_vector_5.add(temp_vector_4(i),             (*std::next(itr_y_begin,    i)));
5247: 	    temp_vector_5.add(temp_vector_4(i + list_size), (*std::next(itr_b0xs_begin, i)));
5248: 	  }
5249: 
5250: 	// temp_vector_2 = B_0 * (x^c - x_k) - W * M * W^T * (x^c - x_k)
5251: 	if (list_size > 0)
5252: 	  temp_vector_2 -= temp_vector_5;
5253: 
5254: 	// temp_vector_2 = g + B_0 * (x^c - x_k) - W * M * W^T * (x^c - x_k)
5255: 	temp_vector_2 += m_system_rhs;
5256: 
5257: 	// temp_vector_2 = Z^T * [g + B_0 * (x^c - x_k) - W * M * W^T * (x^c - x_k)]
5258: 	for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
5259: 	  {
5260: 	    if (free_dofs.block(m_u_dof)[i] < 0)
5261: 	      temp_vector_2.block(m_u_dof)[i] = 0;
5262: 	  }
5263: 
5264: 	for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
5265: 	  {
5266: 	    if (free_dofs.block(m_d_dof)[i] < 0)
5267: 	      temp_vector_2.block(m_d_dof)[i] = 0;
5268: 	  }
5269: 
5270: 	BlockVector<double> rhs_vector(temp_vector_2);
5271: 	rhs_vector *= -1;
5272: 
5273: 	BlockVector<double> search_direction(m_dofs_per_block);
5274: 	search_direction = 0.0;
5275: 
5276: 	unsigned int cg_iterations = 0;
5277: 
5278: 	double alpha_backtrack = 1.0;
5279: 
5280: 	if (m_parameters.m_type_linear_solver == "CG")
5281: 	  {
5282: 	    m_timer.enter_subsection("Subspace CG solve (inverse operator)");
5283: 
5284: 	    //const double rc_hat_norm = rhs_vector.l2_norm();
5285: 	    const double cg_tol = m_parameters.m_CG_tolerace; //std::min( 0.1, std::sqrt(rc_hat_norm) ) * rc_hat_norm;
5286: 
5287: 	    zT_B0_z(free_dofs, m_tangent_matrix);
5288: 
5289: 	    const auto op_zT_B0_z = block_operator(m_tangent_matrix);
5290: 
5291: 	    if (list_size > 0)
5292: 	      {
5293: 		std::list<BlockVector<double>> zT_y_list;
5294: 		BlockVector<double> zT_y_vector(m_dofs_per_block);
5295: 		for (auto itr = y_vector_list.begin(); itr != y_vector_list.end(); ++itr)
5296: 		  {
5297: 		    zT_y_vector = (*itr);
5298: 		    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
5299: 		      {
5300: 			if (free_dofs.block(m_u_dof)[i] < 0)
5301: 			  zT_y_vector.block(m_u_dof)[i] = 0;
5302: 		      }
5303: 
5304: 		    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
5305: 		      {
5306: 			if (free_dofs.block(m_d_dof)[i] < 0)
5307: 			  zT_y_vector.block(m_d_dof)[i] = 0;
5308: 		      }
5309: 
5310: 		    zT_y_list.push_back(zT_y_vector);
5311: 		  }
5312: 
5313: 		std::list<BlockVector<double>> zT_b0xs_list;
5314: 		BlockVector<double> zT_b0xs_vector(m_dofs_per_block);
5315: 		for (auto itr = b0xs_vector_list.begin(); itr != b0xs_vector_list.end(); ++itr)
5316: 		  {
5317: 		    zT_b0xs_vector = (*itr);
5318: 		    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
5319: 		      {
5320: 			if (free_dofs.block(m_u_dof)[i] < 0)
5321: 			  zT_b0xs_vector.block(m_u_dof)[i] = 0;
5322: 		      }
5323: 
5324: 		    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
5325: 		      {
5326: 			if (free_dofs.block(m_d_dof)[i] < 0)
5327: 			  zT_b0xs_vector.block(m_d_dof)[i] = 0;
5328: 		      }
5329: 
5330: 		    zT_b0xs_list.push_back(zT_b0xs_vector);
5331: 		  }
5332: 
5333: 		const auto op_M_matrix = linear_operator(M_matrix);
5334: 
5335: 		FullMatrix<double> zT_W_matrix_u(m_dofs_per_block[m_u_dof], 2*list_size);
5336: 		unsigned int j = 0;
5337: 		for (auto itr = zT_y_list.begin(); itr != zT_y_list.end(); ++itr)
5338: 		  {
5339: 		    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
5340: 		      zT_W_matrix_u(i, j) = (*itr).block(m_u_dof)[i];
5341: 		    ++j;
5342: 		  }
5343: 		j = 0;
5344: 		for (auto itr = zT_b0xs_list.begin(); itr != zT_b0xs_list.end(); ++itr)
5345: 		  {
5346: 		    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
5347: 		      zT_W_matrix_u(i, j + list_size) = (*itr).block(m_u_dof)[i];
5348: 		    ++j;
5349: 		  }
5350: 
5351: 		FullMatrix<double> zT_W_matrix_d(m_dofs_per_block[m_d_dof], 2*list_size);
5352: 		j = 0;
5353: 		for (auto itr = zT_y_list.begin(); itr != zT_y_list.end(); ++itr)
5354: 		  {
5355: 		    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
5356: 		      zT_W_matrix_d(i, j) = (*itr).block(m_d_dof)[i];
5357: 		    ++j;
5358: 		  }
5359: 		j = 0;
5360: 		for (auto itr = zT_b0xs_list.begin(); itr != zT_b0xs_list.end(); ++itr)
5361: 		  {
5362: 		    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
5363: 		      zT_W_matrix_d(i, j + list_size) = (*itr).block(m_d_dof)[i];
5364: 		    ++j;
5365: 		  }
5366: 
5367: 		const auto op_zT_W_matrix_u = linear_operator(zT_W_matrix_u);
5368: 		const auto op_zT_W_matrix_d = linear_operator(zT_W_matrix_d);
5369: 
5370: 		const auto op_uMuT = op_zT_W_matrix_u * op_M_matrix * transpose_operator(op_zT_W_matrix_u);
5371: 
5372: 		const auto op_uMdT = op_zT_W_matrix_u * op_M_matrix * transpose_operator(op_zT_W_matrix_d);
5373: 
5374: 		const auto op_dMuT = op_zT_W_matrix_d * op_M_matrix * transpose_operator(op_zT_W_matrix_u);
5375: 
5376: 		const auto op_dMdT = op_zT_W_matrix_d * op_M_matrix * transpose_operator(op_zT_W_matrix_d);
5377: 
5378: 		const auto op_zT_wMwT_z = block_operator<2, 2, BlockVector<double>>({op_uMuT, op_uMdT,
5379: 										     op_dMuT, op_dMdT});
5380: 
5381: 		const auto op_total = op_zT_B0_z - op_zT_wMwT_z;
5382: 
5383: 		SolverControl            solver_control(1e5, cg_tol);
5384: 		SolverCG<BlockVector<double>> cg(solver_control);
5385: 
5386: 		if (m_parameters.m_type_preconditioner == "None")
5387: 		  {
5388: 		    // somehow op_total_inv has to be made const, or the
5389: 		    // program will have compliation error
5390: 		    const auto op_total_inv = inverse_operator(op_total, cg);
5391: 		    op_total_inv.vmult(search_direction, rhs_vector);
5392: 		  }
5393: 		else if (m_parameters.m_type_preconditioner == "Jacobi")
5394: 		  {
5395: 		    usr_Jacobi_preconditioner preconditioner(m_tangent_matrix);
5396: 
5397: 		    const auto op_total_inv = inverse_operator(op_total, cg, preconditioner);
5398: 		    op_total_inv.vmult(search_direction, rhs_vector);
5399: 		  }
5400: 		else if (m_parameters.m_type_preconditioner == "LU")
5401: 		  {
5402: 		    SparseDirectUMFPACK matrix_factorization;
5403: 		    matrix_factorization.initialize(m_tangent_matrix);
5404: 
5405: 		    usr_sparseLU_preconditioner preconditioner(matrix_factorization);
5406: 		    const auto op_total_inv = inverse_operator(op_total, cg, preconditioner);
5407: 		    op_total_inv.vmult(search_direction, rhs_vector);
5408: 		  }
5409: 		else if (m_parameters.m_type_preconditioner == "ILU")
5410: 		  {
5411: 		    SparseILU<double> SparseILU_disp;
5412: 		    SparseILU_disp.initialize(m_tangent_matrix.block(0, 0));
5413: 		    SparseILU<double> SparseILU_phasefield;
5414: 		    SparseILU_phasefield.initialize(m_tangent_matrix.block(1, 1));
5415: 
5416: 		    usr_sparseILU_preconditioner preconditioner(SparseILU_disp,
5417: 								SparseILU_phasefield);
5418: 
5419: 		    const auto op_total_inv = inverse_operator(op_total, cg, preconditioner);
5420: 		    op_total_inv.vmult(search_direction, rhs_vector);
5421: 		  }
5422: 		else
5423: 		  {
5424: 		    AssertThrow(false, ExcMessage("Preconditioner type not implemented"));
5425: 		  }
5426: 
5427: 		cg_iterations = solver_control.last_step();
5428: 	      } // if (list_size > 0)
5429: 	    else
5430: 	      {
5431: 		const auto op_total = op_zT_B0_z;
5432: 		SolverControl            solver_control(1e5, cg_tol);
5433: 		SolverCG<BlockVector<double>> cg(solver_control);
5434: 
5435: 		if (m_parameters.m_type_preconditioner == "None")
5436: 		  {
5437: 		    // somehow op_total_inv has to be made const, or the
5438: 		    // program will have compliation error
5439: 		    const auto op_total_inv = inverse_operator(op_total, cg);
5440: 		    op_total_inv.vmult(search_direction, rhs_vector);
5441: 		  }
5442: 		else if (m_parameters.m_type_preconditioner == "Jacobi")
5443: 		  {
5444: 		    usr_Jacobi_preconditioner preconditioner(m_tangent_matrix);
5445: 
5446: 		    const auto op_total_inv = inverse_operator(op_total, cg, preconditioner);
5447: 		    op_total_inv.vmult(search_direction, rhs_vector);
5448: 		  }
5449: 		else if (m_parameters.m_type_preconditioner == "LU")
5450: 		  {
5451: 		    SparseDirectUMFPACK matrix_factorization;
5452: 		    matrix_factorization.initialize(m_tangent_matrix);
5453: 
5454: 		    usr_sparseLU_preconditioner preconditioner(matrix_factorization);
5455: 		    const auto op_total_inv = inverse_operator(op_total, cg, preconditioner);
5456: 		    op_total_inv.vmult(search_direction, rhs_vector);
5457: 		  }
5458: 		else if (m_parameters.m_type_preconditioner == "ILU")
5459: 		  {
5460: 		    SparseILU<double> SparseILU_disp;
5461: 		    SparseILU_disp.initialize(m_tangent_matrix.block(0, 0));
5462: 		    SparseILU<double> SparseILU_phasefield;
5463: 		    SparseILU_phasefield.initialize(m_tangent_matrix.block(1, 1));
5464: 
5465: 		    usr_sparseILU_preconditioner preconditioner(SparseILU_disp,
5466: 								SparseILU_phasefield);
5467: 
5468: 		    const auto op_total_inv = inverse_operator(op_total, cg, preconditioner);
5469: 		    op_total_inv.vmult(search_direction, rhs_vector);
5470: 		  }
5471: 		else
5472: 		  {
5473: 		    AssertThrow(false, ExcMessage("Preconditioner type not implemented"));
5474: 		  }
5475: 
5476: 		cg_iterations = solver_control.last_step();
5477: 	      } // // if (list_size == 0)
5478: 
5479:             m_timer.leave_subsection();
5480: 	  } // if (m_parameters.m_type_linear_solver == "CG")
5481: 	else if (m_parameters.m_type_linear_solver == "Direct")
5482: 	  {
5483: 	    m_timer.enter_subsection("Subspace direct solve (LU factorization)");
5484: 
5485: 	    zT_B0_z(free_dofs, m_tangent_matrix);
5486: 
5487: 	    SparseDirectUMFPACK zT_B0_z_inv;
5488: 	    zT_B0_z_inv.initialize(m_tangent_matrix);
5489: 
5490: 	    //SparseDirectUMFPACK zT_B0_z_inv_disp;
5491: 	    //zT_B0_z_inv_disp.initialize(m_tangent_matrix.block(m_u_dof, m_u_dof));
5492: 
5493: 	    //SparseDirectUMFPACK zT_B0_z_inv_phasefield;
5494: 	    //zT_B0_z_inv_phasefield.initialize(m_tangent_matrix.block(m_d_dof, m_d_dof));
5495: 
5496: 	    m_timer.leave_subsection();
5497: 
5498: 	    m_timer.enter_subsection("Subspace direct solve (LU matrix-vector multiplication)");
5499: 
5500: 	    zT_B0_z_inv.vmult(search_direction, rhs_vector);
5501: 	    //zT_B0_z_inv_disp.vmult(search_direction.block(m_u_dof), rhs_vector.block(m_u_dof));
5502: 	    //zT_B0_z_inv_phasefield.vmult(search_direction.block(m_d_dof), rhs_vector.block(m_d_dof));
5503: 
5504: 	    BlockVector<double> update_vector(m_dofs_per_block);
5505: 	    update_vector = 0;
5506: 	    if (list_size > 0)
5507: 	      {
5508: 		std::list<BlockVector<double>> zT_B0_z_inv_zT_y_list;
5509: 		std::list<BlockVector<double>> zT_y_list;
5510: 		BlockVector<double> zT_y_vector(m_dofs_per_block);
5511: 		BlockVector<double> zT_B0_z_inv_zT_y_vector(m_dofs_per_block);
5512: 		for (auto itr = y_vector_list.begin(); itr != y_vector_list.end(); ++itr)
5513: 		  {
5514: 		    zT_y_vector = (*itr);
5515: 		    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
5516: 		      {
5517: 			if (free_dofs.block(m_u_dof)[i] < 0)
5518: 			  zT_y_vector.block(m_u_dof)[i] = 0;
5519: 		      }
5520: 
5521: 		    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
5522: 		      {
5523: 			if (free_dofs.block(m_d_dof)[i] < 0)
5524: 			  zT_y_vector.block(m_d_dof)[i] = 0;
5525: 		      }
5526: 
5527: 		    zT_y_list.push_back(zT_y_vector);
5528: 
5529: 		    //zT_B0_z_inv_disp.vmult(zT_B0_z_inv_zT_y_vector.block(m_u_dof), zT_y_vector.block(m_u_dof));
5530: 		    //zT_B0_z_inv_phasefield.vmult(zT_B0_z_inv_zT_y_vector.block(m_d_dof), zT_y_vector.block(m_d_dof));
5531: 		    zT_B0_z_inv.vmult(zT_B0_z_inv_zT_y_vector, zT_y_vector);
5532: 
5533: 		    zT_B0_z_inv_zT_y_list.push_back(zT_B0_z_inv_zT_y_vector);
5534: 		  }
5535: 
5536: 		std::list<BlockVector<double>> zT_B0_z_inv_zT_b0xs_list;
5537: 		std::list<BlockVector<double>> zT_b0xs_list;
5538: 		BlockVector<double> zT_b0xs_vector(m_dofs_per_block);
5539: 		BlockVector<double> zT_B0_z_inv_zT_b0xs_vector(m_dofs_per_block);
5540: 		for (auto itr = b0xs_vector_list.begin(); itr != b0xs_vector_list.end(); ++itr)
5541: 		  {
5542: 		    zT_b0xs_vector = (*itr);
5543: 		    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
5544: 		      {
5545: 			if (free_dofs.block(m_u_dof)[i] < 0)
5546: 			  zT_b0xs_vector.block(m_u_dof)[i] = 0;
5547: 		      }
5548: 
5549: 		    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
5550: 		      {
5551: 			if (free_dofs.block(m_d_dof)[i] < 0)
5552: 			  zT_b0xs_vector.block(m_d_dof)[i] = 0;
5553: 		      }
5554: 
5555: 		    zT_b0xs_list.push_back(zT_b0xs_vector);
5556: 
5557: 		    //zT_B0_z_inv_disp.vmult(zT_B0_z_inv_zT_b0xs_vector.block(m_u_dof), zT_b0xs_vector.block(m_u_dof));
5558: 		    //zT_B0_z_inv_phasefield.vmult(zT_B0_z_inv_zT_b0xs_vector.block(m_d_dof), zT_b0xs_vector.block(m_d_dof));
5559: 		    zT_B0_z_inv.vmult(zT_B0_z_inv_zT_b0xs_vector, zT_b0xs_vector);
5560: 
5561: 		    zT_B0_z_inv_zT_b0xs_list.push_back(zT_B0_z_inv_zT_b0xs_vector);
5562: 		  }
5563: 
5564: 		FullMatrix<double> wT_z_zT_B0_z_inv_zT_w(2 * list_size);
5565: 		const auto itr_zT_y_list_begin = zT_y_list.begin();
5566: 		const auto itr_zT_b0xs_list_begin = zT_b0xs_list.begin();
5567: 		const auto itr_zT_B0_z_inv_zT_y_list_begin = zT_B0_z_inv_zT_y_list.begin();
5568: 		const auto itr_zT_B0_z_inv_zT_b0xs_list_begin = zT_B0_z_inv_zT_b0xs_list.begin();
5569: 		for (unsigned int i = 0; i < list_size; ++i)
5570: 		  for (unsigned int j = 0; j < list_size; ++j)
5571: 		    {
5572: 		      wT_z_zT_B0_z_inv_zT_w(i          , j          ) = (*std::next(itr_zT_y_list_begin            , i))
5573: 								      * (*std::next(itr_zT_B0_z_inv_zT_y_list_begin, j));
5574: 
5575: 		      wT_z_zT_B0_z_inv_zT_w(i          , j+list_size) = (*std::next(itr_zT_y_list_begin               , i))
5576: 								      * (*std::next(itr_zT_B0_z_inv_zT_b0xs_list_begin, j));
5577: 
5578: 		      wT_z_zT_B0_z_inv_zT_w(i+list_size, j          ) = (*std::next(itr_zT_b0xs_list_begin         , i))
5579: 								      * (*std::next(itr_zT_B0_z_inv_zT_y_list_begin, j));
5580: 
5581: 		      wT_z_zT_B0_z_inv_zT_w(i+list_size, j+list_size) = (*std::next(itr_zT_b0xs_list_begin            , i))
5582: 								      * (*std::next(itr_zT_B0_z_inv_zT_b0xs_list_begin, j));
5583: 		    }
5584: 
5585: 		FullMatrix<double> temp_matrix(2 * list_size);
5586: 		M_matrix.mmult(temp_matrix, wT_z_zT_B0_z_inv_zT_w);
5587: 
5588: 		FullMatrix<double> middle_matrix(IdentityMatrix(2*list_size));
5589: 		middle_matrix.add(-1.0, temp_matrix);
5590: 
5591: 		FullMatrix<double> middle_matrix_inv(2 * list_size);
5592: 		middle_matrix_inv.invert(middle_matrix);
5593: 
5594: 		middle_matrix_inv.mmult(middle_matrix, M_matrix);
5595: 
5596: 		Vector<double> wT_z_zT_B0_z_inv_rhs(2 * list_size);
5597: 		for (unsigned int i = 0; i < list_size; ++i)
5598: 		  {
5599: 		    wT_z_zT_B0_z_inv_rhs(i            ) = (*std::next(itr_zT_B0_z_inv_zT_y_list_begin   , i)) * rhs_vector;
5600: 		    wT_z_zT_B0_z_inv_rhs(i + list_size) = (*std::next(itr_zT_B0_z_inv_zT_b0xs_list_begin, i)) * rhs_vector;
5601: 		  }
5602: 
5603: 		Vector<double> middle_matrix_wT_z_zT_B0_z_inv_rhs(2 * list_size);
5604: 		middle_matrix.vmult(middle_matrix_wT_z_zT_B0_z_inv_rhs,
5605: 				    wT_z_zT_B0_z_inv_rhs);
5606: 
5607: 		unsigned int index = 0;
5608: 		for (auto itr = zT_B0_z_inv_zT_y_list.begin(); itr != zT_B0_z_inv_zT_y_list.end(); ++itr)
5609: 		  {
5610: 		    update_vector.add(middle_matrix_wT_z_zT_B0_z_inv_rhs(index), *itr);
5611: 		    ++index;
5612: 		  }
5613: 		for (auto itr = zT_B0_z_inv_zT_b0xs_list.begin(); itr != zT_B0_z_inv_zT_b0xs_list.end(); ++itr)
5614: 		  {
5615: 		    update_vector.add(middle_matrix_wT_z_zT_B0_z_inv_rhs(index), *itr);
5616: 		    ++index;
5617: 		  }
5618: 	      } //	if (list_size > 0)
5619: 
5620: 	    search_direction += update_vector;
5621: 
5622: 	    m_timer.leave_subsection();
5623: 	  } // else if (m_parameters.m_type_linear_solver == "Direct")
5624: 	else
5625: 	  {
5626: 	    AssertThrow(false, ExcMessage("Linear solver type not implemented"));
5627: 	  }
5628: 
5629: 	// We don't do backtrack yet. We will make sure phasefield
5630: 	// remains feasible later
5631: 	alpha_backtrack = 1.0;
5632: 	search_direction *= alpha_backtrack;
5633: 
5634: 	LBFGS_update.block(m_u_dof) = solution_delta_cauchy_point.block(m_u_dof);
5635: 	LBFGS_update.block(m_u_dof) += search_direction.block(m_u_dof);
5636: 	LBFGS_update.block(m_u_dof) -= solution_delta.block(m_u_dof);
5637: 
5638: 	for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
5639: 	  {
5640: 	    // phasefield active constraints
5641: 	    if (m_active_set_phasefield(i) > 0.5)
5642: 	      LBFGS_update.block(m_d_dof)[i] = solution_delta_cauchy_point.block(m_d_dof)[i]
5643: 					     - solution_delta.block(m_d_dof)[i];
5644: 	    else
5645: 	      LBFGS_update.block(m_d_dof)[i] = solution_delta_cauchy_point.block(m_d_dof)[i]
5646: 					     + search_direction.block(m_d_dof)[i]
5647: 					     - solution_delta.block(m_d_dof)[i];
5648: 	  }
5649: 
5650: 	// make sure the phasefield solutions are feasible
5651: 	for(unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
5652: 	  {
5653: 	    if (solution_delta.block(m_d_dof)[i] + LBFGS_update.block(m_d_dof)[i] < 0.0)
5654: 	      LBFGS_update.block(m_d_dof)[i] = -solution_delta.block(m_d_dof)[i];
5655: 
5656: 	    if (  solution_delta.block(m_d_dof)[i]
5657: 		+ m_solution.block(m_d_dof)[i]
5658: 		+ LBFGS_update.block(m_d_dof)[i] > 1.0)
5659: 	      LBFGS_update.block(m_d_dof)[i] = 1.0 - m_solution.block(m_d_dof)[i]
5660: 						   - solution_delta.block(m_d_dof)[i];
5661: 	  }
5662: 
5663: 	m_constraints.distribute(LBFGS_update);
5664: 
5665: 	// We need a line search algorithm to decide line_search_parameter
5666: 
5667:         if(m_parameters.m_type_line_search == "StrongWolfe")
5668:           {
5669: 	    const double phi_0 = calculate_energy_functional();
5670: 	    const double phi_0_prime = m_system_rhs * LBFGS_update;
5671: 
5672: 	    line_search_parameter = line_search_stepsize_strong_wolfe(phi_0,
5673: 								      phi_0_prime,
5674: 								      LBFGS_update,
5675: 								      solution_delta);
5676:           }
5677:         else if(m_parameters.m_type_line_search == "GradientBased")
5678:           {
5679: 	    // LBFGS_r_vector is the search direction
5680: 	    line_search_parameter = line_search_stepsize_gradient_based(LBFGS_update,
5681: 									solution_delta);
5682:           }
5683:         else
5684:           {
5685:             Assert(false, ExcMessage("An unknown line search method is called!"));
5686:           }
5687: 
5688: 	LBFGS_update *= line_search_parameter;
5689: 
5690:         get_error_update(LBFGS_update, m_error_update);
5691:         if (LBFGS_iteration == 1)
5692:           m_error_update_0 = m_error_update;
5693: 
5694:         m_error_update_norm = m_error_update;
5695:         // For three-point bending problem and the sphere inclusion problem,
5696:         // we use absolute residual for convergence test
5697:         if (m_parameters.m_relative_residual)
5698:           m_error_update_norm.normalize(m_error_update_0);
5699: 
5700: 	solution_delta += LBFGS_update;
5701: 
5702: 	update_qph_incremental(solution_delta, m_solution);
5703: 
5704:         LBFGS_y_vector = m_system_rhs;
5705:         LBFGS_y_vector *= -1.0;
5706:         assemble_system_rhs_BFGS_parallel(m_solution, m_system_rhs);
5707:         // if we use assemble_system_rhs_BFGS_parallel, then condense() is not necessary
5708:         //m_constraints.condense(m_system_rhs);
5709:         LBFGS_y_vector += m_system_rhs;
5710: 
5711:         LBFGS_s_vector = LBFGS_update;
5712: 
5713: 	// s_vector_list, y_vector_list, s_dot_y_list only need to discard
5714: 	// the front (oldest) item and add the newest item to the end at
5715: 	// each L-BFGS iteration
5716: 	double s_dot_y = LBFGS_s_vector * LBFGS_y_vector;
5717: 	if (s_dot_y > 1.0e-16 * LBFGS_y_vector.norm_sqr())
5718: 	  {
5719: 	    if (list_size >= LBFGS_m)
5720: 	      {
5721: 		s_vector_list.pop_front();
5722: 		y_vector_list.pop_front();
5723: 		s_dot_y_list.pop_front();
5724: 	      }
5725: 
5726: 	    s_vector_list.push_back(LBFGS_s_vector);
5727: 	    y_vector_list.push_back(LBFGS_y_vector);
5728: 	    s_dot_y_list.push_back(s_dot_y);
5729: 	  }
5730: 
5731: 	Vector<double> solution_phasefield_total(m_solution.block(m_d_dof));
5732: 	solution_phasefield_total += solution_delta.block(m_d_dof);
5733: 
5734: 	// Since line search parameter might be less than one, we need update
5735: 	// the phasefield active set status
5736: 	// upper bound is 1.0, lower bound is the solution at the previous step.
5737: 	unsigned int number_active_constraint_lower_bound = 0;
5738: 	unsigned int number_active_constraint_upper_bound = 0;
5739: 	unsigned int number_active_constraint_lowerupper_bound = 0;
5740: 
5741: 	for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
5742: 	  {
5743: 	    if (   solution_delta.block(m_d_dof)[i] == 0.0
5744: 		&& solution_phasefield_total[i] == 1.0)
5745: 	      {
5746: 		m_active_set_phasefield(i) = 3; //lower bound overlaps with upper bound
5747: 		++number_active_constraint_lowerupper_bound;
5748: 	      }
5749: 	    else if (   solution_delta.block(m_d_dof)[i] == 0.0
5750: 		     && solution_phasefield_total[i] != 1.0)
5751: 	      {
5752: 		m_active_set_phasefield(i) = 1; //lower bound
5753: 		++number_active_constraint_lower_bound;
5754: 	      }
5755: 	    else if (   solution_phasefield_total[i] == 1.0
5756: 		     && solution_delta.block(m_d_dof)[i] != 0.0)
5757: 	      {
5758: 	        m_active_set_phasefield(i) = 2; //upper bound
5759: 	        ++number_active_constraint_upper_bound;
5760: 	      }
5761: 	    else
5762: 	      {
5763: 	        m_active_set_phasefield(i) = 0;
5764: 	      }
5765: 	  }
5766: 
5767: 	lower_bound_number_old = lower_bound_number_new;
5768: 	upper_bound_number_old = upper_bound_number_new;
5769: 	lowerupper_bound_number_old = lowerupper_bound_number_new;
5770: 
5771: 	lower_bound_number_new = number_active_constraint_lower_bound;
5772: 	upper_bound_number_new = number_active_constraint_upper_bound;
5773: 	lowerupper_bound_number_new = number_active_constraint_lowerupper_bound;
5774: 
5775: 	if (m_parameters.m_output_iteration_history)
5776:           {
5777: 	    const double energy_functional = calculate_energy_functional();
5778: 
5779: 	    m_logfile << "  | "
5780: 		      << std::setw(6) << number_active_constraint_lower_bound
5781: 		      << std::setw(6) << number_active_constraint_upper_bound
5782: 		      << std::setw(6) << number_active_constraint_lowerupper_bound;
5783:             if (m_parameters.m_type_linear_solver == "CG")
5784:               m_logfile << std::setw(8) << cg_iterations;
5785:             else
5786:               m_logfile << std::setw(8) << "---";
5787:             m_logfile << "      "
5788:         	      << std::fixed << std::setprecision(3) << std::scientific << line_search_parameter
5789: 		      << std::fixed << std::setprecision(6) << std::setw(1) << std::scientific
5790: 		      << "  " << energy_functional
5791: 		      << std::fixed << std::setprecision(3) << std::setw(1)
5792: 					<< std::scientific
5793: 		      << "  " << m_error_residual_norm.m_norm
5794: 		      << "  " << m_error_residual_norm.m_u
5795: 		      << "  " << m_error_residual_norm.m_d
5796: 		      << "  " << m_error_update_norm.m_norm
5797: 		      << "  " << m_error_update_norm.m_u
5798: 		      << "  " << m_error_update_norm.m_d
5799: 		      << "  " << std::endl;
5800:           }
5801:       } // for (; LBFGS_iteration < m_parameters.m_max_iterations_BFGS; ++LBFGS_iteration)
5802: 
5803:     AssertThrow(LBFGS_iteration < m_parameters.m_max_iterations_BFGS,
5804:                 ExcMessage("No convergence in L-BFGS-B nonlinear solver!"));
5805:   }
5806: 
5807:   template <int dim>
5808:   void PhaseFieldMonolithicSolve<dim>::output_results() const
5809:   {
5810:     m_timer.enter_subsection("Output results");
5811: 
5812:     DataOut<dim> data_out;
5813: 
5814:     std::vector<DataComponentInterpretation::DataComponentInterpretation>
5815:       data_component_interpretation(
5816:         dim, DataComponentInterpretation::component_is_part_of_vector);
5817: 
5818:     data_component_interpretation.push_back(
5819:       DataComponentInterpretation::component_is_scalar);
5820: 
5821:     std::vector<std::string> solution_name(dim, "displacement");
5822:     solution_name.emplace_back("phasefield");
5823: 
5824:     data_out.attach_dof_handler(m_dof_handler);
5825:     data_out.add_data_vector(m_solution,
5826:                              solution_name,
5827:                              DataOut<dim>::type_dof_data,
5828:                              data_component_interpretation);
5829: 
5830:     // output phasefield active set status
5831:     BlockVector<double> active_set_status(m_dofs_per_block);
5832:     active_set_status.block(m_d_dof) = m_active_set_phasefield;
5833:     std::vector<DataComponentInterpretation::DataComponentInterpretation>
5834:       data_component_interpretation_active_set(
5835:         dim+1, DataComponentInterpretation::component_is_scalar);
5836:     std::vector<std::string> solution_name_active_set;
5837:     solution_name_active_set.emplace_back("disp_x_active_set");
5838:     solution_name_active_set.emplace_back("disp_y_active_set");
5839:     if (dim ==3)
5840:       solution_name_active_set.emplace_back("disp_z_active_set");
5841:     solution_name_active_set.emplace_back("phasefield_active_set");
5842:     data_out.add_data_vector(active_set_status,
5843: 			     solution_name_active_set,
5844: 			     DataOut<dim>::type_dof_data,
5845: 			     data_component_interpretation_active_set);
5846: 
5847:     Vector<double> cell_material_id(m_triangulation.n_active_cells());
5848:     // output material ID for each cell
5849:     for (const auto &cell : m_triangulation.active_cell_iterators())
5850:       {
5851: 	cell_material_id(cell->active_cell_index()) = cell->material_id();
5852:       }
5853:     data_out.add_data_vector(cell_material_id, "materialID");
5854: 
5855:     // Stress L2 projection
5856:     DoFHandler<dim> stresses_dof_handler_L2(m_triangulation);
5857:     FE_Q<dim>     stresses_fe_L2(m_parameters.m_poly_degree); //FE_Q element is continuous
5858:     stresses_dof_handler_L2.distribute_dofs(stresses_fe_L2);
5859:     AffineConstraints<double> constraints;
5860:     constraints.clear();
5861:     DoFTools::make_hanging_node_constraints(stresses_dof_handler_L2, constraints);
5862:     constraints.close();
5863:     std::vector<DataComponentInterpretation::DataComponentInterpretation>
5864: 	  data_component_interpretation_stress(1,
5865: 					       DataComponentInterpretation::component_is_scalar);
5866: 
5867:     for (unsigned int i = 0; i < dim; ++i)
5868:       for (unsigned int j = i; j < dim; ++j)
5869: 	{
5870: 	  Vector<double> stress_field_L2;
5871: 	  stress_field_L2.reinit(stresses_dof_handler_L2.n_dofs());
5872: 
5873: 	  MappingQ<dim> mapping(m_parameters.m_poly_degree + 1);
5874: 	  VectorTools::project(mapping,
5875: 			       stresses_dof_handler_L2,
5876: 			       constraints,
5877: 			       m_qf_cell,
5878: 			       [&] (const typename DoFHandler<dim>::active_cell_iterator & cell,
5879: 				    const unsigned int q) -> double
5880: 			       {
5881: 				 return m_quadrature_point_history.get_data(cell)[q]->get_cauchy_stress()[i][j];
5882: 			       },
5883: 			       stress_field_L2);
5884: 
5885: 	  std::string stress_name = "Cauchy_stress_" + std::to_string(i+1) + std::to_string(j+1)
5886: 				  + "_L2";
5887: 
5888: 	  data_out.add_data_vector(stresses_dof_handler_L2,
5889: 				   stress_field_L2,
5890: 				   stress_name,
5891: 				   data_component_interpretation_stress);
5892: 	}
5893: 
5894:     data_out.build_patches(m_parameters.m_poly_degree);
5895: 
5896:     std::ofstream output("Solution-" + std::to_string(dim) + "d-" +
5897: 			 Utilities::int_to_string(m_time.get_timestep(),4) + ".vtu");
5898: 
5899:     data_out.write_vtu(output);
5900:     m_timer.leave_subsection();
5901:   }
5902: 
5903:   template <int dim>
5904:   void PhaseFieldMonolithicSolve<dim>::calculate_reaction_force(unsigned int face_ID)
5905:   {
5906:     m_timer.enter_subsection("Calculate reaction force");
5907: 
5908:     BlockVector<double>       system_rhs;
5909:     system_rhs.reinit(m_dofs_per_block);
5910: 
5911:     Vector<double> cell_rhs(m_dofs_per_cell);
5912:     std::vector<types::global_dof_index> local_dof_indices(m_dofs_per_cell);
5913: 
5914:     const double time_ramp = (m_time.current() / m_time.end());
5915:     std::vector<Tensor<1, dim>> rhs_values(m_n_q_points);
5916:     const UpdateFlags uf_cell(update_values | update_gradients |
5917: 			      update_quadrature_points | update_JxW_values);
5918:     const UpdateFlags uf_face(update_values | update_normal_vectors |
5919:                               update_JxW_values);
5920: 
5921:     FEValues<dim> fe_values(m_fe, m_qf_cell, uf_cell);
5922:     FEFaceValues<dim> fe_face_values(m_fe, m_qf_face, uf_face);
5923: 
5924:     // shape function values for displacement field
5925:     std::vector<std::vector<Tensor<1, dim>>>
5926:       Nx(m_qf_cell.size(), std::vector<Tensor<1, dim>>(m_dofs_per_cell));
5927:     std::vector<std::vector<Tensor<2, dim>>>
5928:       grad_Nx(m_qf_cell.size(), std::vector<Tensor<2, dim>>(m_dofs_per_cell));
5929:     std::vector<std::vector<SymmetricTensor<2, dim>>>
5930:       symm_grad_Nx(m_qf_cell.size(), std::vector<SymmetricTensor<2, dim>>(m_dofs_per_cell));
5931: 
5932:     for (const auto &cell : m_dof_handler.active_cell_iterators())
5933:       {
5934: 	// if calculate_reaction_force() is defined as const, then
5935: 	// we also need to put a const in std::shared_ptr,
5936: 	// that is, std::shared_ptr<const PointHistory<dim>>
5937: 	const std::vector<std::shared_ptr< PointHistory<dim>>> lqph =
5938: 	  m_quadrature_point_history.get_data(cell);
5939: 	Assert(lqph.size() == m_n_q_points, ExcInternalError());
5940:         cell_rhs = 0.0;
5941:         fe_values.reinit(cell);
5942:         right_hand_side(fe_values.get_quadrature_points(),
5943:     		        rhs_values,
5944:     		        m_parameters.m_x_component*time_ramp,
5945:     		        m_parameters.m_y_component*time_ramp,
5946:     		        m_parameters.m_z_component*time_ramp);
5947: 
5948:         for (const unsigned int q_point : fe_values.quadrature_point_indices())
5949:           {
5950:             for (const unsigned int k : fe_values.dof_indices())
5951:               {
5952:                 const unsigned int k_group = m_fe.system_to_base_index(k).first.first;
5953: 
5954:                 if (k_group == m_u_dof)
5955:                   {
5956:     		    Nx[q_point][k] = fe_values[m_u_fe].value(k, q_point);
5957:     		    grad_Nx[q_point][k] = fe_values[m_u_fe].gradient(k, q_point);
5958:     		    symm_grad_Nx[q_point][k] = symmetrize(grad_Nx[q_point][k]);
5959:                   }
5960:               }
5961:           }
5962: 
5963:         for (const unsigned int q_point : fe_values.quadrature_point_indices())
5964:           {
5965:             const SymmetricTensor<2, dim> & cauchy_stress = lqph[q_point]->get_cauchy_stress();
5966: 
5967:             const std::vector<Tensor<1,dim>> & N = Nx[q_point];
5968:             const std::vector<SymmetricTensor<2, dim>> & symm_grad_N = symm_grad_Nx[q_point];
5969:             const double JxW = fe_values.JxW(q_point);
5970: 
5971:             for (const unsigned int i : fe_values.dof_indices())
5972:               {
5973:                 const unsigned int i_group = m_fe.system_to_base_index(i).first.first;
5974: 
5975:                 if (i_group == m_u_dof)
5976:                   {
5977:                     cell_rhs(i) -= (symm_grad_N[i] * cauchy_stress) * JxW;
5978:     		    // contributions from the body force to right-hand side
5979:     		    cell_rhs(i) += N[i] * rhs_values[q_point] * JxW;
5980:                   }
5981:               }
5982:           }
5983: 
5984:         // if there is surface pressure, this surface pressure always applied to the
5985:         // reference configuration
5986:         const unsigned int face_pressure_id = 100;
5987:         const double p0 = 0.0;
5988: 
5989:         for (const auto &face : cell->face_iterators())
5990:           {
5991: 	    if (face->at_boundary() && face->boundary_id() == face_pressure_id)
5992: 	      {
5993: 		fe_face_values.reinit(cell, face);
5994: 
5995: 		for (const unsigned int f_q_point : fe_face_values.quadrature_point_indices())
5996: 		  {
5997: 		    const Tensor<1, dim> &N = fe_face_values.normal_vector(f_q_point);
5998: 
5999: 		    const double         pressure  = p0 * time_ramp;
6000: 		    const Tensor<1, dim> traction  = pressure * N;
6001: 
6002: 		    for (const unsigned int i : fe_values.dof_indices())
6003: 		      {
6004: 			const unsigned int i_group = m_fe.system_to_base_index(i).first.first;
6005: 
6006: 			if (i_group == m_u_dof)
6007: 			  {
6008: 			    const unsigned int component_i = m_fe.system_to_component_index(i).first;
6009: 			    const double Ni = fe_face_values.shape_value(i, f_q_point);
6010: 			    const double JxW = fe_face_values.JxW(f_q_point);
6011: 			    cell_rhs(i) += (Ni * traction[component_i]) * JxW;
6012: 			  }
6013: 		      }
6014: 		  }
6015: 	      }
6016:           }
6017: 
6018:         cell->get_dof_indices(local_dof_indices);
6019:         for (const unsigned int i : fe_values.dof_indices())
6020:           system_rhs(local_dof_indices[i]) += cell_rhs(i);
6021:       } // for (const auto &cell : m_dof_handler.active_cell_iterators())
6022: 
6023:     // The difference between the above assembled system_rhs and m_system_rhs
6024:     // is that m_system_rhs is condensed by the m_constraints, which zero out
6025:     // the rhs values associated with the constrained DOFs and modify the rhs
6026:     // values associated with the unconstrained DOFs.
6027: 
6028:     std::vector< types::global_dof_index > mapping;
6029:     std::set<types::boundary_id> boundary_ids;
6030:     boundary_ids.insert(face_ID);
6031:     DoFTools::map_dof_to_boundary_indices(m_dof_handler,
6032: 					  boundary_ids,
6033: 					  mapping);
6034: 
6035:     std::vector<double> reaction_force(dim, 0.0);
6036: 
6037:     for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
6038:       {
6039: 	if (mapping[i] != numbers::invalid_dof_index)
6040: 	  {
6041: 	    reaction_force[i % dim] += system_rhs.block(m_u_dof)(i);
6042: 	  }
6043:       }
6044: 
6045:     for (unsigned int i = 0; i < dim; i++)
6046:       m_logfile << "\t\tReaction force in direction " << i << " on boundary ID " << face_ID
6047:                 << " = "
6048: 		<< std::fixed << std::setprecision(3) << std::setw(1)
6049:                 << std::scientific
6050: 		<< reaction_force[i] << std::endl;
6051: 
6052:     std::pair<double, std::vector<double>> time_force;
6053:     time_force.first = m_time.current();
6054:     time_force.second = reaction_force;
6055:     m_history_reaction_force.push_back(time_force);
6056: 
6057:     m_timer.leave_subsection();
6058:   }
6059: 
6060:   template <int dim>
6061:   void PhaseFieldMonolithicSolve<dim>::write_history_data()
6062:   {
6063:     m_logfile << "\t\tWrite history data ... \n"<<std::endl;
6064: 
6065:     std::ofstream myfile_reaction_force ("Reaction_force.hist");
6066:     if (myfile_reaction_force.is_open())
6067:     {
6068:       myfile_reaction_force << 0.0 << "\t";
6069:       if (dim == 2)
6070: 	myfile_reaction_force << 0.0 << "\t"
6071: 	       << 0.0 << std::endl;
6072:       if (dim == 3)
6073: 	myfile_reaction_force << 0.0 << "\t"
6074: 	       << 0.0 << "\t"
6075: 	       << 0.0 << std::endl;
6076: 
6077:       for (auto const & time_force : m_history_reaction_force)
6078: 	{
6079: 	  myfile_reaction_force << time_force.first << "\t";
6080: 	  if (dim == 2)
6081: 	    myfile_reaction_force << time_force.second[0] << "\t"
6082: 	           << time_force.second[1] << std::endl;
6083: 	  if (dim == 3)
6084: 	    myfile_reaction_force << time_force.second[0] << "\t"
6085: 	           << time_force.second[1] << "\t"
6086: 		   << time_force.second[2] << std::endl;
6087: 	}
6088:       myfile_reaction_force.close();
6089:     }
6090:     else
6091:       m_logfile << "Unable to open file";
6092: 
6093:     std::ofstream myfile_energy ("Energy.hist");
6094:     if (myfile_energy.is_open())
6095:     {
6096:       myfile_energy << std::fixed << std::setprecision(10) << std::scientific
6097:                     << 0.0 << "\t"
6098:                     << 0.0 << "\t"
6099: 	            << 0.0 << "\t"
6100: 	            << 0.0 << std::endl;
6101: 
6102:       for (auto const & time_energy : m_history_energy)
6103: 	{
6104: 	  myfile_energy << std::fixed << std::setprecision(10) << std::scientific
6105: 	                << time_energy.first     << "\t"
6106:                         << time_energy.second[0] << "\t"
6107: 	                << time_energy.second[1] << "\t"
6108: 		        << time_energy.second[2] << std::endl;
6109: 	}
6110:       myfile_energy.close();
6111:     }
6112:     else
6113:       m_logfile << "Unable to open file";
6114:   }
6115: 
6116:   template <int dim>
6117:   double PhaseFieldMonolithicSolve<dim>::calculate_energy_functional() const
6118:   {
6119:     double energy_functional = 0.0;
6120: 
6121:     FEValues<dim> fe_values(m_fe, m_qf_cell, update_JxW_values);
6122: 
6123:     for (const auto &cell : m_dof_handler.active_cell_iterators())
6124:       {
6125:         fe_values.reinit(cell);
6126: 
6127:         const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
6128:           m_quadrature_point_history.get_data(cell);
6129:         Assert(lqph.size() == m_n_q_points, ExcInternalError());
6130: 
6131:         for (unsigned int q_point = 0; q_point < m_n_q_points; ++q_point)
6132:           {
6133:             const double JxW = fe_values.JxW(q_point);
6134:             energy_functional += lqph[q_point]->get_total_strain_energy() * JxW;
6135:             energy_functional += lqph[q_point]->get_crack_energy_dissipation() * JxW;
6136:           }
6137:       }
6138: 
6139:     return energy_functional;
6140:   }
6141: 
6142:   template <int dim>
6143:   std::pair<double, double>
6144:     PhaseFieldMonolithicSolve<dim>::calculate_total_strain_energy_and_crack_energy_dissipation() const
6145:   {
6146:     double total_strain_energy = 0.0;
6147:     double crack_energy_dissipation = 0.0;
6148: 
6149:     FEValues<dim> fe_values(m_fe, m_qf_cell, update_JxW_values);
6150: 
6151:     for (const auto &cell : m_dof_handler.active_cell_iterators())
6152:       {
6153:         fe_values.reinit(cell);
6154: 
6155:         const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
6156:           m_quadrature_point_history.get_data(cell);
6157:         Assert(lqph.size() == m_n_q_points, ExcInternalError());
6158: 
6159:         for (unsigned int q_point = 0; q_point < m_n_q_points; ++q_point)
6160:           {
6161:             const double JxW = fe_values.JxW(q_point);
6162:             total_strain_energy += lqph[q_point]->get_total_strain_energy() * JxW;
6163:             crack_energy_dissipation += lqph[q_point]->get_crack_energy_dissipation() * JxW;
6164:           }
6165:       }
6166: 
6167:     return std::make_pair(total_strain_energy, crack_energy_dissipation);
6168:   }
6169: 
6170:   template <int dim>
6171:   bool PhaseFieldMonolithicSolve<dim>::local_refine_and_solution_transfer(BlockVector<double> & solution_delta,
6172: 									  BlockVector<double> & LBFGS_update_refine)
6173:   {
6174:     // This is the solution at (n+1) obtained from the old (coarse) mesh
6175:     BlockVector<double> solution_next_step(m_dofs_per_block);
6176:     solution_next_step = m_solution + solution_delta;
6177:     bool mesh_is_same = true;
6178:     bool cell_refine_flag = true;
6179: 
6180:     unsigned int material_id;
6181:     double length_scale;
6182:     double cell_length;
6183:     while(cell_refine_flag)
6184:       {
6185: 	cell_refine_flag = false;
6186: 
6187: 	std::vector<types::global_dof_index> local_dof_indices(m_fe.dofs_per_cell);
6188: 	for (const auto &cell : m_dof_handler.active_cell_iterators())
6189: 	  {
6190: 	    cell->get_dof_indices(local_dof_indices);
6191: 
6192: 	    for (unsigned int i = 0; i< m_fe.dofs_per_cell; ++i)
6193: 	      {
6194: 		const unsigned int comp_i = m_fe.system_to_component_index(i).first;
6195: 		if (comp_i == m_d_component) //phasefield component
6196: 		  {
6197: 		    if (  solution_next_step(local_dof_indices[i])
6198: 			> m_parameters.m_phasefield_refine_threshold )
6199: 		      {
6200: 			material_id = cell->material_id();
6201: 	                length_scale = m_material_data[material_id][2];
6202: 	                if (dim == 2)
6203: 	                  cell_length = std::sqrt(cell->measure());
6204: 	                else
6205: 	                  cell_length = std::cbrt(cell->measure());
6206: 			if (  cell_length
6207: 			    > length_scale * m_parameters.m_allowed_max_h_l_ratio )
6208: 			  {
6209: 			    if (cell->level() < m_parameters.m_max_allowed_refinement_level)
6210: 			      {
6211: 			        cell->set_refine_flag();
6212: 			        break;
6213: 			      }
6214: 			  }
6215: 		      }
6216: 		  }
6217: 	      }
6218: 	  }
6219: 
6220: 	for (const auto &cell : m_dof_handler.active_cell_iterators())
6221: 	  {
6222: 	    if (cell->refine_flag_set())
6223: 	      {
6224: 		cell_refine_flag = true;
6225: 		break;
6226: 	      }
6227: 	  }
6228: 
6229: 	// if any cell is refined, we need to project the solution
6230: 	// to the newly refined mesh
6231: 	if (cell_refine_flag)
6232: 	  {
6233: 	    mesh_is_same = false;
6234: 
6235: 	    std::vector<BlockVector<double> > old_solutions(2);
6236: 	    old_solutions[0] = solution_next_step;
6237: 	    old_solutions[1] = m_solution;
6238: 
6239: 	    m_triangulation.prepare_coarsening_and_refinement();
6240: 	    SolutionTransfer<dim, BlockVector<double>> solution_transfer(m_dof_handler);
6241: 	    solution_transfer.prepare_for_coarsening_and_refinement(old_solutions);
6242: 	    m_triangulation.execute_coarsening_and_refinement();
6243: 
6244: 	    setup_system();
6245: 
6246: 	    std::vector<BlockVector<double>> tmp_solutions(2);
6247: 	    tmp_solutions[0].reinit(m_dofs_per_block);
6248: 	    tmp_solutions[1].reinit(m_dofs_per_block);
6249: 
6250:             #  if DEAL_II_VERSION_GTE(9, 7, 0)
6251: 	    solution_transfer.interpolate(tmp_solutions);
6252: 	    #  else
6253: 	    // If an older version of dealII is used, for example, 9.4.0, interpolate()
6254:             // needs to use the following interface.
6255:             solution_transfer.interpolate(old_solutions, tmp_solutions);
6256:             #  endif
6257: 	    solution_next_step = tmp_solutions[0];
6258: 	    m_solution = tmp_solutions[1];
6259: 
6260: 	    // make sure the projected solutions still satisfy
6261: 	    // hanging node constraints
6262: 	    m_constraints.distribute(solution_next_step);
6263: 	    m_constraints.distribute(m_solution);
6264: 	  } // if (cell_refine_flag)
6265:       } // while(cell_refine_flag)
6266: 
6267:     // calculate field variables for newly refined cells
6268:     if (!mesh_is_same)
6269:       {
6270: 	BlockVector<double> temp_solution_delta(m_dofs_per_block);
6271: 	BlockVector<double> temp_previous_solution(m_dofs_per_block);
6272: 	temp_solution_delta = 0.0;
6273: 	temp_previous_solution = 0.0;
6274: 	update_qph_incremental(temp_solution_delta, temp_previous_solution);
6275: 
6276: 	// initial guess for the resolve on the refined mesh
6277: 	LBFGS_update_refine = solution_next_step - m_solution;
6278:       }
6279: 
6280:     return mesh_is_same;
6281:   }
6282: 
6283:   template <int dim>
6284:   void PhaseFieldMonolithicSolve<dim>::print_parameter_information()
6285:   {
6286:     if (m_parameters.m_type_nonlinear_solver == "LBFGS")
6287:       {
6288: 	m_logfile << "WARNING: this version of LBFGS does not enforce"
6289: 	    " phase-field irreversibility." << std::endl;
6290: 	m_logfile << "It should only be used for demonstrating the importance of"
6291: 	    " inequality constraints." << std::endl;
6292: 	m_logfile << "The obtained result is meaningless!" << std::endl;
6293:       }
6294: 
6295:     m_logfile << "Scenario number = " << m_parameters.m_scenario << std::endl;
6296:     m_logfile << "Log file = " << m_parameters.m_logfile_name << std::endl;
6297:     m_logfile << "Write iteration history to log file? = " << std::boolalpha
6298: 	      << m_parameters.m_output_iteration_history << std::endl;
6299: 
6300:     if (dim == 2)
6301:       {
6302: 	if (m_parameters.m_plane_stress)
6303: 	  m_logfile << "2D plane-stress case" << std::endl;
6304: 	else
6305: 	  m_logfile << "2D plane-strain case" << std::endl;
6306:       }
6307: 
6308:     m_logfile << "Nonlinear solver type = " << m_parameters.m_type_nonlinear_solver << std::endl;
6309:     m_logfile << "Line search type = " << m_parameters.m_type_line_search << std::endl;
6310:     m_logfile << "Linear solver type = " << m_parameters.m_type_linear_solver << std::endl;
6311: 
6312:     if (m_parameters.m_type_linear_solver == "CG")
6313:       {
6314:         m_logfile << "Preconditioner type for CG = " << m_parameters.m_type_preconditioner << std::endl;
6315:         m_logfile << "Convergence tolerance for CG iterations = " << m_parameters.m_CG_tolerace << std::endl;
6316:       }
6317: 
6318:     m_logfile << "Mesh refinement strategy = " << m_parameters.m_refinement_strategy << std::endl;
6319: 
6320:     if (m_parameters.m_refinement_strategy == "adaptive-refine")
6321:       {
6322: 	m_logfile << "\tMaximum adaptive refinement times allowed in each step = "
6323: 		  << m_parameters.m_max_adaptive_refine_times << std::endl;
6324: 	m_logfile << "\tMaximum allowed cell refinement level = "
6325: 		  << m_parameters.m_max_allowed_refinement_level << std::endl;
6326: 	m_logfile << "\tPhasefield-based refinement threshold value = "
6327: 		  << m_parameters.m_phasefield_refine_threshold << std::endl;
6328:       }
6329: 
6330:     m_logfile << "L-BFGS_m = " << m_parameters.m_LBFGS_m << std::endl;
6331:     m_logfile << "Global refinement times = " << m_parameters.m_global_refine_times << std::endl;
6332:     m_logfile << "Local prerefinement times = " <<m_parameters. m_local_prerefine_times << std::endl;
6333:     m_logfile << "Allowed maximum h/l ratio = " << m_parameters.m_allowed_max_h_l_ratio << std::endl;
6334:     m_logfile << "total number of material types = " << m_parameters.m_total_material_regions << std::endl;
6335:     m_logfile << "material data file name = " << m_parameters.m_material_file_name << std::endl;
6336:     if (m_parameters.m_reaction_force_face_id >= 0)
6337:       m_logfile << "Calculate reaction forces on Face ID = " << m_parameters.m_reaction_force_face_id << std::endl;
6338:     else
6339:       m_logfile << "No need to calculate reaction forces." << std::endl;
6340: 
6341:     if (m_parameters.m_relative_residual)
6342:       m_logfile << "Relative residual for convergence." << std::endl;
6343:     else
6344:       m_logfile << "Absolute residual for convergence." << std::endl;
6345: 
6346:     m_logfile << "Body force = (" << m_parameters.m_x_component << ", "
6347:                                   << m_parameters.m_y_component << ", "
6348: 	                          << m_parameters.m_z_component << ") (N/m^3)"
6349: 				  << std::endl;
6350: 
6351:     m_logfile << "End time = " << m_parameters.m_end_time << std::endl;
6352:     m_logfile << "Time data file name = " << m_parameters.m_time_file_name << std::endl;
6353:   }
6354: 
6355:   template <int dim>
6356:   void PhaseFieldMonolithicSolve<dim>::run()
6357:   {
6358:     print_parameter_information();
6359: 
6360:     read_material_data(m_parameters.m_material_file_name,
6361:     		       m_parameters.m_total_material_regions);
6362: 
6363:     std::vector<std::array<double, 4>> time_table;
6364: 
6365:     read_time_data(m_parameters.m_time_file_name, time_table);
6366: 
6367:     make_grid();
6368:     setup_system();
6369:     output_results();
6370: 
6371:     m_time.increment(time_table);
6372: 
6373:     while(m_time.current() < m_time.end() + m_time.get_delta_t()*1.0e-6)
6374:       {
6375: 	m_logfile << std::endl
6376: 		  << "Timestep " << m_time.get_timestep() << " @ " << m_time.current()
6377: 		  << 's' << std::endl;
6378: 
6379:         bool mesh_is_same = false;
6380: 
6381:         // initial guess for the resolve on the refined mesh
6382: 	BlockVector<double> LBFGS_update_refine(m_dofs_per_block);
6383: 	LBFGS_update_refine = 0.0;
6384: 
6385:         // local adaptive mesh refinement loop
6386: 	unsigned int adp_refine_iteration = 0;
6387:         for (; adp_refine_iteration < m_parameters.m_max_adaptive_refine_times + 1; ++adp_refine_iteration)
6388:           {
6389: 	    if (m_parameters.m_refinement_strategy == "adaptive-refine")
6390: 	      m_logfile << "\tAdaptive refinement-"<< adp_refine_iteration << ": " << std::endl;
6391: 
6392: 	    BlockVector<double> solution_delta(m_dofs_per_block);
6393: 	    solution_delta = 0.0;
6394: 
6395:             if (m_parameters.m_type_nonlinear_solver == "LBFGS")
6396: 	      solve_nonlinear_timestep_LBFGS(solution_delta, LBFGS_update_refine);
6397:             else if (m_parameters.m_type_nonlinear_solver == "LBFGSB")
6398:               solve_nonlinear_timestep_LBFGS_B(solution_delta, LBFGS_update_refine);
6399: 	    else
6400: 	      AssertThrow(false, ExcMessage("Nonlinear solver type not implemented"));
6401: 
6402: 	    if (m_parameters.m_refinement_strategy == "adaptive-refine")
6403: 	      {
6404: 
6405: 		if (adp_refine_iteration == m_parameters.m_max_adaptive_refine_times)
6406: 		  {
6407: 		    m_solution += solution_delta;
6408: 		    break;
6409: 		  }
6410: 
6411: 		mesh_is_same = local_refine_and_solution_transfer(solution_delta,
6412: 								  LBFGS_update_refine);
6413: 
6414: 		if (mesh_is_same)
6415: 		  {
6416: 		    m_solution += solution_delta;
6417: 		    break;
6418: 		  }
6419: 	      }
6420: 	    else if (m_parameters.m_refinement_strategy == "pre-refine")
6421: 	      {
6422: 		m_solution += solution_delta;
6423: 	        break;
6424: 	      }
6425: 	    else
6426: 	      {
6427: 		AssertThrow(false,
6428: 		            ExcMessage("Selected mesh refinement strategy not implemented!"));
6429: 	      }
6430:           } // for (; adp_refine_iteration < m_parameters.m_max_adaptive_refine_times; ++adp_refine_iteration)
6431: 
6432:         //AssertThrow(adp_refine_iteration < m_parameters.m_max_adaptive_refine_times,
6433:         //            ExcMessage("Number of local adaptive mesh refinement exceeds allowed maximum times!"));
6434: 
6435: 	// output vtk files every 10 steps if there are too
6436: 	// many time steps
6437: 	//if (m_time.get_timestep() % 10 == 0)
6438:         output_results();
6439: 
6440: 	double energy_functional_current = calculate_energy_functional();
6441: 	m_logfile << "\t\tEnergy functional (J) = " << std::fixed << std::setprecision(10) << std::scientific
6442: 	          << energy_functional_current << std::endl;
6443: 
6444: 	std::pair<double, double> energy_pair = calculate_total_strain_energy_and_crack_energy_dissipation();
6445: 	m_logfile << "\t\tTotal strain energy (J) = " << std::fixed << std::setprecision(10) << std::scientific
6446: 		  << energy_pair.first << std::endl;
6447: 	m_logfile << "\t\tCrack energy dissipation (J) = " << std::fixed << std::setprecision(10) << std::scientific
6448: 		  << energy_pair.second << std::endl;
6449: 
6450: 	std::pair<double, std::array<double, 3>> time_energy;
6451: 	time_energy.first = m_time.current();
6452: 	time_energy.second[0] = energy_pair.first;
6453: 	time_energy.second[1] = energy_pair.second;
6454: 	time_energy.second[2] = energy_pair.first + energy_pair.second;
6455: 	m_history_energy.push_back(time_energy);
6456: 
6457: 	int face_ID = m_parameters.m_reaction_force_face_id;
6458: 	if (face_ID >= 0)
6459: 	  calculate_reaction_force(face_ID);
6460: 
6461:         write_history_data();
6462: 
6463: 	m_time.increment(time_table);
6464:       } // while(m_time.current() < m_time.end() + m_time.get_delta_t()*1.0e-6)
6465:   }
6466: } // namespace PhaseField
6467: 
6468: int main(int argc, char* argv[])
6469: {
6470:   using namespace dealii;
6471: 
6472:   if (argc != 2)
6473:     AssertThrow(false,
6474:     		ExcMessage("The number of arguments provided to the program has to be 2!"));
6475: 
6476:   const unsigned int dim = std::stoi(argv[1]);
6477:   if (dim == 2 )
6478:     {
6479:       PhaseField::PhaseFieldMonolithicSolve<2> problem_2D("parameters.prm");
6480:       problem_2D.run();
6481:     }
6482:   else if (dim == 3)
6483:     {
6484:       PhaseField::PhaseFieldMonolithicSolve<3> problem_3D("parameters.prm");
6485:       problem_3D.run();
6486:     }
6487:   else
6488:     {
6489:       AssertThrow(false,
6490:                   ExcMessage("Dimension has to be either 2 or 3"));
6491:     }
6492: 
6493:   return 0;
6494: }
```

## 4. 按行号详细解读（覆盖全文件，逐行）

> 本节逐行给出说明；涉及核心计算的行可结合第 2 节公式。

- **[行 1]** `/* ---------------------------------------------------------------------`：注释行，用于说明算法背景、假设或实现意图。
- **[行 2]** ` *`：注释行，用于说明算法背景、假设或实现意图。
- **[行 3]** ` * Copyright (C) 2006 - 2020 by the deal.II authors`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4]** ` *`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5]** ` * This file is part of the deal.II library.`：注释行，用于说明算法背景、假设或实现意图。
- **[行 6]** ` *`：注释行，用于说明算法背景、假设或实现意图。
- **[行 7]** ` * The deal.II library is free software; you can use it, redistribute`：注释行，用于说明算法背景、假设或实现意图。
- **[行 8]** ` * it, and/or modify it under the terms of the GNU Lesser General`：注释行，用于说明算法背景、假设或实现意图。
- **[行 9]** ` * Public License as published by the Free Software Foundation; either`：注释行，用于说明算法背景、假设或实现意图。
- **[行 10]** ` * version 2.1 of the License, or (at your option) any later version.`：注释行，用于说明算法背景、假设或实现意图。
- **[行 11]** ` * The full text of the license can be found in the file LICENSE.md at`：注释行，用于说明算法背景、假设或实现意图。
- **[行 12]** ` * the top level directory of deal.II.`：注释行，用于说明算法背景、假设或实现意图。
- **[行 13]** ` *`：注释行，用于说明算法背景、假设或实现意图。
- **[行 14]** ` * ---------------------------------------------------------------------`：注释行，用于说明算法背景、假设或实现意图。
- **[行 15]** ``：空行，用于分隔逻辑块。
- **[行 16]** ` *`：注释行，用于说明算法背景、假设或实现意图。
- **[行 17]** ` * Author: Tao Jin, PhD`：注释行，用于说明算法背景、假设或实现意图。
- **[行 18]** ` *         University of Ottawa, Ottawa, Ontario, Canada`：注释行，用于说明算法背景、假设或实现意图。
- **[行 19]** ` *         July 2024`：注释行，用于说明算法背景、假设或实现意图。
- **[行 20]** ` */`：注释行，用于说明算法背景、假设或实现意图。
- **[行 21]** ``：空行，用于分隔逻辑块。
- **[行 22]** `/* A monolithic scheme based on the L-BFGS method and the gradient projection method`：注释行，用于说明算法背景、假设或实现意图。
- **[行 23]** ` *  to solve the phase-field crack problem`：注释行，用于说明算法背景、假设或实现意图。
- **[行 24]** ` * 1. The phase-field method treats the phasefield irreversibility using the gradient`：注释行，用于说明算法背景、假设或实现意图。
- **[行 25]** ` *    projection method. During a load step [t_n, t_n+1], let d_n represent the phasefield`：注释行，用于说明算法背景、假设或实现意图。
- **[行 26]** ` *    at the beginning of the load step (known), then the inequality constraints`：注释行，用于说明算法背景、假设或实现意图。
- **[行 27]** ` *    d_n <= d_n+1 <= 1.0 are treated as box constraints.`：注释行，用于说明算法背景、假设或实现意图。
- **[行 28]** ` * 2. Various direct and iterative linear solvers are designed to improve the wall-clock`：注释行，用于说明算法背景、假设或实现意图。
- **[行 29]** ` *    run time.`：注释行，用于说明算法背景、假设或实现意图。
- **[行 30]** ` * 3. Using TBB for stiffness assembly and Gauss point calculation.`：注释行，用于说明算法背景、假设或实现意图。
- **[行 31]** ` * 4. Using adaptive mesh refinement.`：注释行，用于说明算法背景、假设或实现意图。
- **[行 32]** ` */`：注释行，用于说明算法背景、假设或实现意图。
- **[行 33]** ``：空行，用于分隔逻辑块。
- **[行 34]** `#include <deal.II/grid/tria.h>`：引入头文件 <deal.II/grid/tria.h>，提供后续类型/函数定义。
- **[行 35]** `#include <deal.II/grid/grid_generator.h>`：引入头文件 <deal.II/grid/grid_generator.h>，提供后续类型/函数定义。
- **[行 36]** `#include <deal.II/grid/grid_refinement.h>`：引入头文件 <deal.II/grid/grid_refinement.h>，提供后续类型/函数定义。
- **[行 37]** `#include <deal.II/grid/grid_out.h>`：引入头文件 <deal.II/grid/grid_out.h>，提供后续类型/函数定义。
- **[行 38]** `#include <deal.II/grid/grid_in.h>`：引入头文件 <deal.II/grid/grid_in.h>，提供后续类型/函数定义。
- **[行 39]** `#include <deal.II/grid/manifold_lib.h>`：引入头文件 <deal.II/grid/manifold_lib.h>，提供后续类型/函数定义。
- **[行 40]** ``：空行，用于分隔逻辑块。
- **[行 41]** `#include <deal.II/dofs/dof_handler.h>`：引入头文件 <deal.II/dofs/dof_handler.h>，提供后续类型/函数定义。
- **[行 42]** `#include <deal.II/dofs/dof_tools.h>`：引入头文件 <deal.II/dofs/dof_tools.h>，提供后续类型/函数定义。
- **[行 43]** `#include <deal.II/dofs/dof_renumbering.h>`：引入头文件 <deal.II/dofs/dof_renumbering.h>，提供后续类型/函数定义。
- **[行 44]** ``：空行，用于分隔逻辑块。
- **[行 45]** `#include <deal.II/fe/fe_values.h>`：引入头文件 <deal.II/fe/fe_values.h>，提供后续类型/函数定义。
- **[行 46]** `#include <deal.II/fe/fe_system.h>`：引入头文件 <deal.II/fe/fe_system.h>，提供后续类型/函数定义。
- **[行 47]** `#include <deal.II/fe/fe_q.h>`：引入头文件 <deal.II/fe/fe_q.h>，提供后续类型/函数定义。
- **[行 48]** `#include <deal.II/fe/fe_dgp_monomial.h>`：引入头文件 <deal.II/fe/fe_dgp_monomial.h>，提供后续类型/函数定义。
- **[行 49]** `#include <deal.II/fe/mapping_q_eulerian.h>`：引入头文件 <deal.II/fe/mapping_q_eulerian.h>，提供后续类型/函数定义。
- **[行 50]** ``：空行，用于分隔逻辑块。
- **[行 51]** `#include <deal.II/base/timer.h>`：引入头文件 <deal.II/base/timer.h>，提供后续类型/函数定义。
- **[行 52]** `#include <deal.II/base/quadrature_point_data.h>`：引入头文件 <deal.II/base/quadrature_point_data.h>，提供后续类型/函数定义。
- **[行 53]** `#include <deal.II/base/parameter_handler.h>`：引入头文件 <deal.II/base/parameter_handler.h>，提供后续类型/函数定义。
- **[行 54]** ``：空行，用于分隔逻辑块。
- **[行 55]** `#include <deal.II/lac/affine_constraints.h>`：引入头文件 <deal.II/lac/affine_constraints.h>，提供后续类型/函数定义。
- **[行 56]** `#include <deal.II/lac/vector.h>`：引入头文件 <deal.II/lac/vector.h>，提供后续类型/函数定义。
- **[行 57]** `#include <deal.II/lac/full_matrix.h>`：引入头文件 <deal.II/lac/full_matrix.h>，提供后续类型/函数定义。
- **[行 58]** `#include <deal.II/lac/sparse_matrix.h>`：引入头文件 <deal.II/lac/sparse_matrix.h>，提供后续类型/函数定义。
- **[行 59]** `#include <deal.II/lac/dynamic_sparsity_pattern.h>`：引入头文件 <deal.II/lac/dynamic_sparsity_pattern.h>，提供后续类型/函数定义。
- **[行 60]** `#include <deal.II/lac/block_sparse_matrix.h>`：引入头文件 <deal.II/lac/block_sparse_matrix.h>，提供后续类型/函数定义。
- **[行 61]** `#include <deal.II/lac/block_vector.h>`：引入头文件 <deal.II/lac/block_vector.h>，提供后续类型/函数定义。
- **[行 62]** ``：空行，用于分隔逻辑块。
- **[行 63]** `#include <deal.II/numerics/vector_tools.h>`：引入头文件 <deal.II/numerics/vector_tools.h>，提供后续类型/函数定义。
- **[行 64]** `#include <deal.II/numerics/matrix_tools.h>`：引入头文件 <deal.II/numerics/matrix_tools.h>，提供后续类型/函数定义。
- **[行 65]** `#include <deal.II/numerics/data_out.h>`：引入头文件 <deal.II/numerics/data_out.h>，提供后续类型/函数定义。
- **[行 66]** ``：空行，用于分隔逻辑块。
- **[行 67]** `#include <deal.II/lac/solver_cg.h>`：引入头文件 <deal.II/lac/solver_cg.h>，提供后续类型/函数定义。
- **[行 68]** `#include <deal.II/lac/precondition.h>`：引入头文件 <deal.II/lac/precondition.h>，提供后续类型/函数定义。
- **[行 69]** `#include <deal.II/lac/linear_operator.h>`：引入头文件 <deal.II/lac/linear_operator.h>，提供后续类型/函数定义。
- **[行 70]** `#include <deal.II/lac/packaged_operation.h>`：引入头文件 <deal.II/lac/packaged_operation.h>，提供后续类型/函数定义。
- **[行 71]** `#include <deal.II/lac/precondition_selector.h>`：引入头文件 <deal.II/lac/precondition_selector.h>，提供后续类型/函数定义。
- **[行 72]** `#include <deal.II/lac/solver_selector.h>`：引入头文件 <deal.II/lac/solver_selector.h>，提供后续类型/函数定义。
- **[行 73]** `#include <deal.II/lac/sparse_direct.h>`：引入头文件 <deal.II/lac/sparse_direct.h>，提供后续类型/函数定义。
- **[行 74]** ``：空行，用于分隔逻辑块。
- **[行 75]** `#include <deal.II/numerics/error_estimator.h>`：引入头文件 <deal.II/numerics/error_estimator.h>，提供后续类型/函数定义。
- **[行 76]** ``：空行，用于分隔逻辑块。
- **[行 77]** `#include <deal.II/physics/elasticity/standard_tensors.h>`：引入头文件 <deal.II/physics/elasticity/standard_tensors.h>，提供后续类型/函数定义。
- **[行 78]** ``：空行，用于分隔逻辑块。
- **[行 79]** `#include <deal.II/base/quadrature_point_data.h>`：引入头文件 <deal.II/base/quadrature_point_data.h>，提供后续类型/函数定义。
- **[行 80]** ``：空行，用于分隔逻辑块。
- **[行 81]** `#include <deal.II/grid/grid_tools.h>`：引入头文件 <deal.II/grid/grid_tools.h>，提供后续类型/函数定义。
- **[行 82]** ``：空行，用于分隔逻辑块。
- **[行 83]** `#include <deal.II/base/work_stream.h>`：引入头文件 <deal.II/base/work_stream.h>，提供后续类型/函数定义。
- **[行 84]** ``：空行，用于分隔逻辑块。
- **[行 85]** `#include <deal.II/numerics/solution_transfer.h>`：引入头文件 <deal.II/numerics/solution_transfer.h>，提供后续类型/函数定义。
- **[行 86]** ``：空行，用于分隔逻辑块。
- **[行 87]** `#include <deal.II/lac/linear_operator_tools.h>`：引入头文件 <deal.II/lac/linear_operator_tools.h>，提供后续类型/函数定义。
- **[行 88]** `#include <deal.II/lac/sparse_ilu.h>`：引入头文件 <deal.II/lac/sparse_ilu.h>，提供后续类型/函数定义。
- **[行 89]** ``：空行，用于分隔逻辑块。
- **[行 90]** `#include <fstream>`：引入头文件 <fstream>，提供后续类型/函数定义。
- **[行 91]** `#include <iostream>`：引入头文件 <iostream>，提供后续类型/函数定义。
- **[行 92]** ``：空行，用于分隔逻辑块。
- **[行 93]** `#include <deal.II/base/logstream.h>`：引入头文件 <deal.II/base/logstream.h>，提供后续类型/函数定义。
- **[行 94]** ``：空行，用于分隔逻辑块。
- **[行 95]** `#include "SpectrumDecomposition.h"`：引入头文件 "SpectrumDecomposition.h"，提供后续类型/函数定义。
- **[行 96]** `#include "Utilities.h"`：引入头文件 "Utilities.h"，提供后续类型/函数定义。
- **[行 97]** ``：空行，用于分隔逻辑块。
- **[行 98]** `namespace PhaseField`：命名空间声明，组织符号并避免命名冲突。
- **[行 99]** `{`：作用域边界（代码块开始/结束）。
- **[行 100]** `  using namespace dealii;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 101]** ``：空行，用于分隔逻辑块。
- **[行 102]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 103]** `  std::vector<types::global_dof_index> get_vertex_dofs(`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 104]** `    const typename Triangulation<dim>::active_vertex_iterator &vertex,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 105]** `    const DoFHandler<dim> &dof_handler)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 106]** `  {`：作用域边界（代码块开始/结束）。
- **[行 107]** `    DoFAccessor<0, dim, dim, false> vertex_dofs(`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 108]** `        &(dof_handler.get_triangulation()),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 109]** `        vertex->level(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 110]** `        vertex->index(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 111]** `        &dof_handler);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 112]** `    const unsigned int n_dofs = dof_handler.get_fe().dofs_per_vertex;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 113]** `    std::vector<types::global_dof_index> dofs(n_dofs);`：调用 C++ 标准库工具函数/容器接口。
- **[行 114]** `    for (unsigned int i = 0; i < n_dofs; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 115]** `    {`：作用域边界（代码块开始/结束）。
- **[行 116]** `      dofs[i] = vertex_dofs.vertex_dof_index(0, i);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 117]** `    }`：作用域边界（代码块开始/结束）。
- **[行 118]** `    return dofs;`：返回当前函数结果。
- **[行 119]** `  }`：作用域边界（代码块开始/结束）。
- **[行 120]** ``：空行，用于分隔逻辑块。
- **[行 121]** `  // Jacobi preconditioner`：注释行，用于说明算法背景、假设或实现意图。
- **[行 122]** `  class usr_Jacobi_preconditioner : public Subscriptor`：类型声明（类/结构体），封装数据与行为。
- **[行 123]** `  {`：作用域边界（代码块开始/结束）。
- **[行 124]** `  public:`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 125]** `    usr_Jacobi_preconditioner(const BlockSparseMatrix<double> & S);`：函数调用语句，触发对应计算或操作。
- **[行 126]** ``：空行，用于分隔逻辑块。
- **[行 127]** `    void vmult(BlockVector<double> & dst,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 128]** `	       const BlockVector<double> & src) const;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 129]** ``：空行，用于分隔逻辑块。
- **[行 130]** `  private:`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 131]** `#  if DEAL_II_VERSION_GTE(9, 7, 0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 132]** `    const ObserverPointer<const BlockSparseMatrix<double> > m_system_matrix;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 133]** `#  else`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 134]** `    const SmartPointer<const BlockSparseMatrix<double> > m_system_matrix;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 135]** `#  endif`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 136]** `  };`：作用域边界（代码块开始/结束）。
- **[行 137]** ``：空行，用于分隔逻辑块。
- **[行 138]** `  usr_Jacobi_preconditioner::usr_Jacobi_preconditioner(const BlockSparseMatrix<double> & S)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 139]** `  : m_system_matrix(&S)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 140]** `  {}`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 141]** ``：空行，用于分隔逻辑块。
- **[行 142]** `  void usr_Jacobi_preconditioner::vmult(BlockVector<double> & dst,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 143]** `					const BlockVector<double> & src) const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 144]** `  {`：作用域边界（代码块开始/结束）。
- **[行 145]** `    PreconditionJacobi<BlockSparseMatrix<double>> preconditioner;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 146]** `    preconditioner.initialize(*m_system_matrix, 1.0);`：函数调用语句，触发对应计算或操作。
- **[行 147]** ``：空行，用于分隔逻辑块。
- **[行 148]** `    preconditioner.vmult(dst, src);`：函数调用语句，触发对应计算或操作。
- **[行 149]** `  }`：作用域边界（代码块开始/结束）。
- **[行 150]** ``：空行，用于分隔逻辑块。
- **[行 151]** `  // LU preconditioner`：注释行，用于说明算法背景、假设或实现意图。
- **[行 152]** `  class usr_sparseLU_preconditioner : public Subscriptor`：类型声明（类/结构体），封装数据与行为。
- **[行 153]** `  {`：作用域边界（代码块开始/结束）。
- **[行 154]** `  public:`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 155]** `    usr_sparseLU_preconditioner(const SparseDirectUMFPACK & matrix_factorization);`：函数调用语句，触发对应计算或操作。
- **[行 156]** ``：空行，用于分隔逻辑块。
- **[行 157]** `    void vmult(BlockVector<double> & dst,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 158]** `	       const BlockVector<double> & src) const;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 159]** ``：空行，用于分隔逻辑块。
- **[行 160]** `  private:`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 161]** `#  if DEAL_II_VERSION_GTE(9, 7, 0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 162]** `    const ObserverPointer<const SparseDirectUMFPACK > m_matrix_LU;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 163]** `#  else`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 164]** `    const SmartPointer<const SparseDirectUMFPACK > m_matrix_LU;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 165]** `#  endif`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 166]** `  };`：作用域边界（代码块开始/结束）。
- **[行 167]** ``：空行，用于分隔逻辑块。
- **[行 168]** `  usr_sparseLU_preconditioner::usr_sparseLU_preconditioner(const SparseDirectUMFPACK & matrix_factorization)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 169]** `  : m_matrix_LU(&matrix_factorization)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 170]** `  {}`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 171]** ``：空行，用于分隔逻辑块。
- **[行 172]** `  void usr_sparseLU_preconditioner::vmult(BlockVector<double> & dst,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 173]** `					       const BlockVector<double> & src) const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 174]** `  {`：作用域边界（代码块开始/结束）。
- **[行 175]** `    (*m_matrix_LU).vmult(dst, src);`：函数调用语句，触发对应计算或操作。
- **[行 176]** `  }`：作用域边界（代码块开始/结束）。
- **[行 177]** ``：空行，用于分隔逻辑块。
- **[行 178]** `  // Incomplete LU preconditioner`：注释行，用于说明算法背景、假设或实现意图。
- **[行 179]** `  class usr_sparseILU_preconditioner : public Subscriptor`：类型声明（类/结构体），封装数据与行为。
- **[行 180]** `  {`：作用域边界（代码块开始/结束）。
- **[行 181]** `  public:`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 182]** `    usr_sparseILU_preconditioner(const SparseILU<double> & ILU_factorization_disp,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 183]** `				 const SparseILU<double> & ILU_factorization_phasefield);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 184]** ``：空行，用于分隔逻辑块。
- **[行 185]** `    void vmult(BlockVector<double> & dst,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 186]** `	     const BlockVector<double> & src) const;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 187]** ``：空行，用于分隔逻辑块。
- **[行 188]** `  private:`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 189]** `#  if DEAL_II_VERSION_GTE(9, 7, 0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 190]** `    const ObserverPointer<const SparseILU<double> > m_ILU_factorization_disp;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 191]** `    const ObserverPointer<const SparseILU<double> > m_ILU_factorization_phasefield;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 192]** `#  else`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 193]** `    const SmartPointer<const SparseILU<double> > m_ILU_factorization_disp;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 194]** `    const SmartPointer<const SparseILU<double> > m_ILU_factorization_phasefield;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 195]** `#  endif`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 196]** `  };`：作用域边界（代码块开始/结束）。
- **[行 197]** ``：空行，用于分隔逻辑块。
- **[行 198]** `  usr_sparseILU_preconditioner::usr_sparseILU_preconditioner(const SparseILU<double> & ILU_factorization_disp,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 199]** `							     const SparseILU<double> & ILU_factorization_phasefield)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 200]** `  : m_ILU_factorization_disp(& ILU_factorization_disp)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 201]** `  , m_ILU_factorization_phasefield(& ILU_factorization_phasefield)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 202]** `  {}`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 203]** ``：空行，用于分隔逻辑块。
- **[行 204]** `  void usr_sparseILU_preconditioner::vmult(BlockVector<double> & dst,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 205]** `					   const BlockVector<double> & src) const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 206]** `  {`：作用域边界（代码块开始/结束）。
- **[行 207]** `    std::vector<types::global_dof_index> sizes_per_block(src.n_blocks());`：调用 C++ 标准库工具函数/容器接口。
- **[行 208]** `    for (unsigned int i = 0; i < src.n_blocks(); ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 209]** `      sizes_per_block[i] = src.block(i).size();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 210]** `    dst.reinit(sizes_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 211]** ``：空行，用于分隔逻辑块。
- **[行 212]** `    (*m_ILU_factorization_disp).vmult(dst.block(0), src.block(0));`：函数调用语句，触发对应计算或操作。
- **[行 213]** `    (*m_ILU_factorization_phasefield).vmult(dst.block(1), src.block(1));`：函数调用语句，触发对应计算或操作。
- **[行 214]** `  }`：作用域边界（代码块开始/结束）。
- **[行 215]** ``：空行，用于分隔逻辑块。
- **[行 216]** `  // body force`：注释行，用于说明算法背景、假设或实现意图。
- **[行 217]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 218]** `  void right_hand_side(const std::vector<Point<dim>> &points,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 219]** `		       std::vector<Tensor<1, dim>> &  values,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 220]** `		       const double fx,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 221]** `		       const double fy,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 222]** `		       const double fz)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 223]** `  {`：作用域边界（代码块开始/结束）。
- **[行 224]** `    Assert(values.size() == points.size(),`：运行期断言/检查，验证输入与状态合法性。
- **[行 225]** `           ExcDimensionMismatch(values.size(), points.size()));`：函数调用语句，触发对应计算或操作。
- **[行 226]** `    Assert(dim >= 2, ExcNotImplemented());`：运行期断言/检查，验证输入与状态合法性。
- **[行 227]** ``：空行，用于分隔逻辑块。
- **[行 228]** `    for (unsigned int point_n = 0; point_n < points.size(); ++point_n)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 229]** `      {`：作用域边界（代码块开始/结束）。
- **[行 230]** `	if (dim == 2)`：条件分支：根据当前状态选择执行路径。
- **[行 231]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 232]** `	    values[point_n][0] = fx;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 233]** `	    values[point_n][1] = fy;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 234]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 235]** `	else`：条件分支的兜底路径。
- **[行 236]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 237]** `	    values[point_n][0] = fx;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 238]** `	    values[point_n][1] = fy;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 239]** `	    values[point_n][2] = fz;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 240]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 241]** `      }`：作用域边界（代码块开始/结束）。
- **[行 242]** `  }`：作用域边界（代码块开始/结束）。
- **[行 243]** ``：空行，用于分隔逻辑块。
- **[行 244]** `  double degradation_function(const double d)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 245]** `  {`：作用域边界（代码块开始/结束）。
- **[行 246]** `    return (1.0 - d) * (1.0 - d);`：返回当前函数结果。
- **[行 247]** `  }`：作用域边界（代码块开始/结束）。
- **[行 248]** ``：空行，用于分隔逻辑块。
- **[行 249]** `  double degradation_function_derivative(const double d)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 250]** `  {`：作用域边界（代码块开始/结束）。
- **[行 251]** `    return 2.0 * (d - 1.0);`：返回当前函数结果。
- **[行 252]** `  }`：作用域边界（代码块开始/结束）。
- **[行 253]** ``：空行，用于分隔逻辑块。
- **[行 254]** `  double degradation_function_2nd_order_derivative(const double d)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 255]** `  {`：作用域边界（代码块开始/结束）。
- **[行 256]** `    (void) d;`：函数调用语句，触发对应计算或操作。
- **[行 257]** `    return 2.0;`：返回当前函数结果。
- **[行 258]** `  }`：作用域边界（代码块开始/结束）。
- **[行 259]** ``：空行，用于分隔逻辑块。
- **[行 260]** `  namespace Parameters`：命名空间声明，组织符号并避免命名冲突。
- **[行 261]** `  {`：作用域边界（代码块开始/结束）。
- **[行 262]** `    struct Scenario`：类型声明（类/结构体），封装数据与行为。
- **[行 263]** `    {`：作用域边界（代码块开始/结束）。
- **[行 264]** `      unsigned int m_scenario;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 265]** `      std::string m_logfile_name;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 266]** `      bool m_output_iteration_history;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 267]** `      bool m_plane_stress;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 268]** `      std::string m_type_nonlinear_solver;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 269]** `      std::string m_type_line_search;`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 270]** `      std::string m_type_linear_solver;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 271]** `      std::string m_type_preconditioner;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 272]** `      double m_CG_tolerace;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 273]** `      std::string m_refinement_strategy;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 274]** `      unsigned int m_LBFGS_m;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 275]** `      unsigned int m_global_refine_times;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 276]** `      unsigned int m_local_prerefine_times;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 277]** `      unsigned int m_max_adaptive_refine_times;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 278]** `      int m_max_allowed_refinement_level;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 279]** `      double m_phasefield_refine_threshold;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 280]** `      double m_allowed_max_h_l_ratio;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 281]** `      unsigned int m_total_material_regions;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 282]** `      std::string m_material_file_name;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 283]** `      int m_reaction_force_face_id;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 284]** ``：空行，用于分隔逻辑块。
- **[行 285]** `      static void declare_parameters(ParameterHandler &prm);`：函数调用语句，触发对应计算或操作。
- **[行 286]** `      void parse_parameters(ParameterHandler &prm);`：函数调用语句，触发对应计算或操作。
- **[行 287]** `    };`：作用域边界（代码块开始/结束）。
- **[行 288]** ``：空行，用于分隔逻辑块。
- **[行 289]** `    void Scenario::declare_parameters(ParameterHandler &prm)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 290]** `    {`：作用域边界（代码块开始/结束）。
- **[行 291]** `      prm.enter_subsection("Scenario");`：函数调用语句，触发对应计算或操作。
- **[行 292]** `      {`：作用域边界（代码块开始/结束）。
- **[行 293]** `        prm.declare_entry("Scenario number",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 294]** `                          "1",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 295]** `                          Patterns::Integer(0),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 296]** `                          "Geometry, loading and boundary conditions scenario");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 297]** ``：空行，用于分隔逻辑块。
- **[行 298]** `        prm.declare_entry("Log file name",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 299]** `			  "Output.log",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 300]** `                          Patterns::FileName(Patterns::FileName::input),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 301]** `			  "Name of the file for log");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 302]** ``：空行，用于分隔逻辑块。
- **[行 303]** `        prm.declare_entry("Output iteration history",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 304]** `			  "yes",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 305]** `                          Patterns::Selection("yes|no"),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 306]** `			  "Shall we write iteration history to the log file?");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 307]** ``：空行，用于分隔逻辑块。
- **[行 308]** `        prm.declare_entry("Plane stress",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 309]** `			  "no",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 310]** `			  Patterns::Selection("yes|no"),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 311]** `			  "If it is 2D, is it plane-stress?");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 312]** ``：空行，用于分隔逻辑块。
- **[行 313]** `        prm.declare_entry("Nonlinear solver type",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 314]** `                          "LBFGSB",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 315]** `                          Patterns::Selection("LBFGS|LBFGSB"),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 316]** `                          "Type of solver used to solve the nonlinear system");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 317]** ``：空行，用于分隔逻辑块。
- **[行 318]** `        prm.declare_entry("Line search type",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 319]** `                          "GradientBased",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 320]** `                          Patterns::Selection("GradientBased|StrongWolfe"),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 321]** `                          "Type of line search method, the gradient-based method "`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 322]** `                          "should be preferred since it is generally faster");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 323]** ``：空行，用于分隔逻辑块。
- **[行 324]** `        prm.declare_entry("Linear solver type",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 325]** `                          "CG",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 326]** `                          Patterns::Selection("Direct|CG"),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 327]** `                          "Type of solver used to solve the linear system");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 328]** ``：空行，用于分隔逻辑块。
- **[行 329]** `        prm.declare_entry("Preconditioner type for CG",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 330]** `                          "ILU",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 331]** `                          Patterns::Selection("None|Jacobi|LU|ILU"),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 332]** `                          "Type of preconditioner used to solve the linear system");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 333]** ``：空行，用于分隔逻辑块。
- **[行 334]** `        prm.declare_entry("CG tolerance",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 335]** `                          "1.0e-6",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 336]** `                          Patterns::Double(0.0),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 337]** `                          "Convergence tolerance of CG iterations");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 338]** ``：空行，用于分隔逻辑块。
- **[行 339]** `        prm.declare_entry("Mesh refinement strategy",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 340]** `                          "adaptive-refine",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 341]** `                          Patterns::Selection("pre-refine|adaptive-refine"),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 342]** `                          "Mesh refinement strategy: pre-refine or adaptive-refine");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 343]** ``：空行，用于分隔逻辑块。
- **[行 344]** `        prm.declare_entry("LBFGS m",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 345]** `                          "40",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 346]** `                          Patterns::Integer(0),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 347]** `                          "Number of vectors used for LBFGS");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 348]** ``：空行，用于分隔逻辑块。
- **[行 349]** `        prm.declare_entry("Global refinement times",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 350]** `                          "0",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 351]** `                          Patterns::Integer(0),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 352]** `                          "Global refinement times (across the entire domain)");`：函数调用语句，触发对应计算或操作。
- **[行 353]** ``：空行，用于分隔逻辑块。
- **[行 354]** `        prm.declare_entry("Local prerefinement times",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 355]** `                          "0",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 356]** `                          Patterns::Integer(0),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 357]** `                          "Local pre-refinement times (assume crack path is known a priori), "`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 358]** `                          "only refine along the crack path.");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 359]** ``：空行，用于分隔逻辑块。
- **[行 360]** `        prm.declare_entry("Max adaptive refinement times",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 361]** `                          "100",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 362]** `                          Patterns::Integer(0),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 363]** `                          "Maximum number of adaptive refinement times allowed in each step");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 364]** ``：空行，用于分隔逻辑块。
- **[行 365]** `        prm.declare_entry("Max allowed refinement level",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 366]** `                          "100",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 367]** `                          Patterns::Integer(0),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 368]** `                          "Maximum allowed cell refinement level");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 369]** ``：空行，用于分隔逻辑块。
- **[行 370]** `        prm.declare_entry("Phasefield refine threshold",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 371]** `			  "0.8",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 372]** `			  Patterns::Double(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 373]** `			  "Phasefield-based refinement threshold value");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 374]** ``：空行，用于分隔逻辑块。
- **[行 375]** `        prm.declare_entry("Allowed max hl ratio",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 376]** `			  "0.25",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 377]** `			  Patterns::Double(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 378]** `			  "Allowed maximum ratio between mesh size h and length scale l");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 379]** ``：空行，用于分隔逻辑块。
- **[行 380]** `        prm.declare_entry("Material regions",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 381]** `                          "1",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 382]** `                          Patterns::Integer(0),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 383]** `                          "Number of material regions");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 384]** ``：空行，用于分隔逻辑块。
- **[行 385]** `        prm.declare_entry("Material data file",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 386]** `                          "1",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 387]** `                          Patterns::FileName(Patterns::FileName::input),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 388]** `                          "Material data file");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 389]** ``：空行，用于分隔逻辑块。
- **[行 390]** `        prm.declare_entry("Reaction force face ID",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 391]** `                          "1",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 392]** `                          Patterns::Integer(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 393]** `                          "Face id where reaction forces should be calculated "`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 394]** `                          "(negative integer means not to calculate reaction force)");`：函数调用语句，触发对应计算或操作。
- **[行 395]** `      }`：作用域边界（代码块开始/结束）。
- **[行 396]** `      prm.leave_subsection();`：函数调用语句，触发对应计算或操作。
- **[行 397]** `    }`：作用域边界（代码块开始/结束）。
- **[行 398]** ``：空行，用于分隔逻辑块。
- **[行 399]** `    void Scenario::parse_parameters(ParameterHandler &prm)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 400]** `    {`：作用域边界（代码块开始/结束）。
- **[行 401]** `      prm.enter_subsection("Scenario");`：函数调用语句，触发对应计算或操作。
- **[行 402]** `      {`：作用域边界（代码块开始/结束）。
- **[行 403]** `        m_scenario = prm.get_integer("Scenario number");`：从参数文件读取并写入成员变量。
- **[行 404]** `        m_logfile_name = prm.get("Log file name");`：从参数文件读取并写入成员变量。
- **[行 405]** `        m_output_iteration_history = prm.get_bool("Output iteration history");`：从参数文件读取并写入成员变量。
- **[行 406]** `        m_plane_stress = prm.get_bool("Plane stress");`：从参数文件读取并写入成员变量。
- **[行 407]** `        m_type_nonlinear_solver = prm.get("Nonlinear solver type");`：从参数文件读取并写入成员变量。
- **[行 408]** `        m_type_line_search = prm.get("Line search type");`：从参数文件读取并写入成员变量。
- **[行 409]** `        m_type_linear_solver = prm.get("Linear solver type");`：从参数文件读取并写入成员变量。
- **[行 410]** `        m_type_preconditioner = prm.get("Preconditioner type for CG");`：从参数文件读取并写入成员变量。
- **[行 411]** `        m_CG_tolerace = prm.get_double("CG tolerance");`：从参数文件读取并写入成员变量。
- **[行 412]** `        m_refinement_strategy = prm.get("Mesh refinement strategy");`：从参数文件读取并写入成员变量。
- **[行 413]** `        m_LBFGS_m = prm.get_integer("LBFGS m");`：从参数文件读取并写入成员变量。
- **[行 414]** `        m_global_refine_times = prm.get_integer("Global refinement times");`：从参数文件读取并写入成员变量。
- **[行 415]** `        m_local_prerefine_times = prm.get_integer("Local prerefinement times");`：从参数文件读取并写入成员变量。
- **[行 416]** `        m_max_adaptive_refine_times = prm.get_integer("Max adaptive refinement times");`：从参数文件读取并写入成员变量。
- **[行 417]** `        m_max_allowed_refinement_level = prm.get_integer("Max allowed refinement level");`：从参数文件读取并写入成员变量。
- **[行 418]** `        m_phasefield_refine_threshold = prm.get_double("Phasefield refine threshold");`：从参数文件读取并写入成员变量。
- **[行 419]** `        m_allowed_max_h_l_ratio = prm.get_double("Allowed max hl ratio");`：从参数文件读取并写入成员变量。
- **[行 420]** `        m_total_material_regions = prm.get_integer("Material regions");`：从参数文件读取并写入成员变量。
- **[行 421]** `        m_material_file_name = prm.get("Material data file");`：从参数文件读取并写入成员变量。
- **[行 422]** `        m_reaction_force_face_id = prm.get_integer("Reaction force face ID");`：从参数文件读取并写入成员变量。
- **[行 423]** `      }`：作用域边界（代码块开始/结束）。
- **[行 424]** `      prm.leave_subsection();`：函数调用语句，触发对应计算或操作。
- **[行 425]** `    }`：作用域边界（代码块开始/结束）。
- **[行 426]** ``：空行，用于分隔逻辑块。
- **[行 427]** `    struct FESystem`：类型声明（类/结构体），封装数据与行为。
- **[行 428]** `    {`：作用域边界（代码块开始/结束）。
- **[行 429]** `      unsigned int m_poly_degree;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 430]** `      unsigned int m_quad_order;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 431]** ``：空行，用于分隔逻辑块。
- **[行 432]** `      static void declare_parameters(ParameterHandler &prm);`：函数调用语句，触发对应计算或操作。
- **[行 433]** ``：空行，用于分隔逻辑块。
- **[行 434]** `      void parse_parameters(ParameterHandler &prm);`：函数调用语句，触发对应计算或操作。
- **[行 435]** `    };`：作用域边界（代码块开始/结束）。
- **[行 436]** ``：空行，用于分隔逻辑块。
- **[行 437]** ``：空行，用于分隔逻辑块。
- **[行 438]** `    void FESystem::declare_parameters(ParameterHandler &prm)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 439]** `    {`：作用域边界（代码块开始/结束）。
- **[行 440]** `      prm.enter_subsection("Finite element system");`：函数调用语句，触发对应计算或操作。
- **[行 441]** `      {`：作用域边界（代码块开始/结束）。
- **[行 442]** `        prm.declare_entry("Polynomial degree",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 443]** `                          "1",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 444]** `                          Patterns::Integer(0),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 445]** `                          "Phase field polynomial order");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 446]** ``：空行，用于分隔逻辑块。
- **[行 447]** `        prm.declare_entry("Quadrature order",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 448]** `                          "2",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 449]** `                          Patterns::Integer(0),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 450]** `                          "Gauss quadrature order");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 451]** `      }`：作用域边界（代码块开始/结束）。
- **[行 452]** `      prm.leave_subsection();`：函数调用语句，触发对应计算或操作。
- **[行 453]** `    }`：作用域边界（代码块开始/结束）。
- **[行 454]** ``：空行，用于分隔逻辑块。
- **[行 455]** `    void FESystem::parse_parameters(ParameterHandler &prm)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 456]** `    {`：作用域边界（代码块开始/结束）。
- **[行 457]** `      prm.enter_subsection("Finite element system");`：函数调用语句，触发对应计算或操作。
- **[行 458]** `      {`：作用域边界（代码块开始/结束）。
- **[行 459]** `        m_poly_degree = prm.get_integer("Polynomial degree");`：从参数文件读取并写入成员变量。
- **[行 460]** `        m_quad_order  = prm.get_integer("Quadrature order");`：从参数文件读取并写入成员变量。
- **[行 461]** `      }`：作用域边界（代码块开始/结束）。
- **[行 462]** `      prm.leave_subsection();`：函数调用语句，触发对应计算或操作。
- **[行 463]** `    }`：作用域边界（代码块开始/结束）。
- **[行 464]** ``：空行，用于分隔逻辑块。
- **[行 465]** `    // body force (N/m^3)`：注释行，用于说明算法背景、假设或实现意图。
- **[行 466]** `    struct BodyForce`：类型声明（类/结构体），封装数据与行为。
- **[行 467]** `    {`：作用域边界（代码块开始/结束）。
- **[行 468]** `      double m_x_component;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 469]** `      double m_y_component;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 470]** `      double m_z_component;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 471]** ``：空行，用于分隔逻辑块。
- **[行 472]** `      static void declare_parameters(ParameterHandler &prm);`：函数调用语句，触发对应计算或操作。
- **[行 473]** ``：空行，用于分隔逻辑块。
- **[行 474]** `      void parse_parameters(ParameterHandler &prm);`：函数调用语句，触发对应计算或操作。
- **[行 475]** `    };`：作用域边界（代码块开始/结束）。
- **[行 476]** ``：空行，用于分隔逻辑块。
- **[行 477]** `    void BodyForce::declare_parameters(ParameterHandler &prm)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 478]** `    {`：作用域边界（代码块开始/结束）。
- **[行 479]** `      prm.enter_subsection("Body force");`：函数调用语句，触发对应计算或操作。
- **[行 480]** `      {`：作用域边界（代码块开始/结束）。
- **[行 481]** `        prm.declare_entry("Body force x component",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 482]** `			  "0.0",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 483]** `			  Patterns::Double(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 484]** `			  "Body force x-component (N/m^3)");`：函数调用语句，触发对应计算或操作。
- **[行 485]** ``：空行，用于分隔逻辑块。
- **[行 486]** `        prm.declare_entry("Body force y component",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 487]** `			  "0.0",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 488]** `			  Patterns::Double(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 489]** `			  "Body force y-component (N/m^3)");`：函数调用语句，触发对应计算或操作。
- **[行 490]** ``：空行，用于分隔逻辑块。
- **[行 491]** `        prm.declare_entry("Body force z component",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 492]** `			  "0.0",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 493]** `			  Patterns::Double(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 494]** `			  "Body force z-component (N/m^3)");`：函数调用语句，触发对应计算或操作。
- **[行 495]** `      }`：作用域边界（代码块开始/结束）。
- **[行 496]** `      prm.leave_subsection();`：函数调用语句，触发对应计算或操作。
- **[行 497]** `    }`：作用域边界（代码块开始/结束）。
- **[行 498]** ``：空行，用于分隔逻辑块。
- **[行 499]** `    void BodyForce::parse_parameters(ParameterHandler &prm)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 500]** `    {`：作用域边界（代码块开始/结束）。
- **[行 501]** `      prm.enter_subsection("Body force");`：函数调用语句，触发对应计算或操作。
- **[行 502]** `      {`：作用域边界（代码块开始/结束）。
- **[行 503]** `        m_x_component = prm.get_double("Body force x component");`：从参数文件读取并写入成员变量。
- **[行 504]** `        m_y_component = prm.get_double("Body force y component");`：从参数文件读取并写入成员变量。
- **[行 505]** `        m_z_component = prm.get_double("Body force z component");`：从参数文件读取并写入成员变量。
- **[行 506]** `      }`：作用域边界（代码块开始/结束）。
- **[行 507]** `      prm.leave_subsection();`：函数调用语句，触发对应计算或操作。
- **[行 508]** `    }`：作用域边界（代码块开始/结束）。
- **[行 509]** ``：空行，用于分隔逻辑块。
- **[行 510]** `    struct NonlinearSolver`：类型声明（类/结构体），封装数据与行为。
- **[行 511]** `    {`：作用域边界（代码块开始/结束）。
- **[行 512]** `      unsigned int m_max_iterations_BFGS;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 513]** `      bool m_relative_residual;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 514]** ``：空行，用于分隔逻辑块。
- **[行 515]** `      double       m_tol_u_residual;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 516]** `      double       m_tol_d_residual;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 517]** `      double       m_tol_u_incr;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 518]** `      double       m_tol_d_incr;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 519]** ``：空行，用于分隔逻辑块。
- **[行 520]** `      static void declare_parameters(ParameterHandler &prm);`：函数调用语句，触发对应计算或操作。
- **[行 521]** ``：空行，用于分隔逻辑块。
- **[行 522]** `      void parse_parameters(ParameterHandler &prm);`：函数调用语句，触发对应计算或操作。
- **[行 523]** `    };`：作用域边界（代码块开始/结束）。
- **[行 524]** ``：空行，用于分隔逻辑块。
- **[行 525]** `    void NonlinearSolver::declare_parameters(ParameterHandler &prm)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 526]** `    {`：作用域边界（代码块开始/结束）。
- **[行 527]** `      prm.enter_subsection("Nonlinear solver");`：函数调用语句，触发对应计算或操作。
- **[行 528]** `      {`：作用域边界（代码块开始/结束）。
- **[行 529]** `        prm.declare_entry("Max iterations BFGS",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 530]** `                          "20",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 531]** `                          Patterns::Integer(0),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 532]** `                          "Number of BFGS iterations allowed");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 533]** ``：空行，用于分隔逻辑块。
- **[行 534]** `        prm.declare_entry("Relative residual",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 535]** `			  "yes",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 536]** `                          Patterns::Selection("yes|no"),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 537]** `			  "Shall we use relative residual for convergence?");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 538]** ``：空行，用于分隔逻辑块。
- **[行 539]** `        prm.declare_entry("Tolerance displacement residual",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 540]** `                          "1.0e-9",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 541]** `                          Patterns::Double(0.0),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 542]** `                          "Displacement residual tolerance");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 543]** ``：空行，用于分隔逻辑块。
- **[行 544]** `        prm.declare_entry("Tolerance phasefield residual",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 545]** `                          "1.0e-9",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 546]** `                          Patterns::Double(0.0),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 547]** `                          "Phasefield residual tolerance");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 548]** ``：空行，用于分隔逻辑块。
- **[行 549]** `        prm.declare_entry("Tolerance displacement increment",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 550]** `                          "1.0e-9",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 551]** `                          Patterns::Double(0.0),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 552]** `                          "Displacement increment tolerance");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 553]** ``：空行，用于分隔逻辑块。
- **[行 554]** `        prm.declare_entry("Tolerance phasefield increment",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 555]** `                          "1.0e-9",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 556]** `                          Patterns::Double(0.0),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 557]** `                          "Phasefield increment tolerance");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 558]** `      }`：作用域边界（代码块开始/结束）。
- **[行 559]** `      prm.leave_subsection();`：函数调用语句，触发对应计算或操作。
- **[行 560]** `    }`：作用域边界（代码块开始/结束）。
- **[行 561]** ``：空行，用于分隔逻辑块。
- **[行 562]** `    void NonlinearSolver::parse_parameters(ParameterHandler &prm)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 563]** `    {`：作用域边界（代码块开始/结束）。
- **[行 564]** `      prm.enter_subsection("Nonlinear solver");`：函数调用语句，触发对应计算或操作。
- **[行 565]** `      {`：作用域边界（代码块开始/结束）。
- **[行 566]** `        m_max_iterations_BFGS = prm.get_integer("Max iterations BFGS");`：从参数文件读取并写入成员变量。
- **[行 567]** `        m_relative_residual = prm.get_bool("Relative residual");`：从参数文件读取并写入成员变量。
- **[行 568]** ``：空行，用于分隔逻辑块。
- **[行 569]** `        m_tol_u_residual           = prm.get_double("Tolerance displacement residual");`：从参数文件读取并写入成员变量。
- **[行 570]** `        m_tol_d_residual           = prm.get_double("Tolerance phasefield residual");`：从参数文件读取并写入成员变量。
- **[行 571]** `        m_tol_u_incr               = prm.get_double("Tolerance displacement increment");`：从参数文件读取并写入成员变量。
- **[行 572]** `        m_tol_d_incr               = prm.get_double("Tolerance phasefield increment");`：从参数文件读取并写入成员变量。
- **[行 573]** `      }`：作用域边界（代码块开始/结束）。
- **[行 574]** `      prm.leave_subsection();`：函数调用语句，触发对应计算或操作。
- **[行 575]** `    }`：作用域边界（代码块开始/结束）。
- **[行 576]** ``：空行，用于分隔逻辑块。
- **[行 577]** `    struct TimeInfo`：类型声明（类/结构体），封装数据与行为。
- **[行 578]** `    {`：作用域边界（代码块开始/结束）。
- **[行 579]** `      double m_end_time;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 580]** `      std::string m_time_file_name;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 581]** ``：空行，用于分隔逻辑块。
- **[行 582]** `      static void declare_parameters(ParameterHandler &prm);`：函数调用语句，触发对应计算或操作。
- **[行 583]** ``：空行，用于分隔逻辑块。
- **[行 584]** `      void parse_parameters(ParameterHandler &prm);`：函数调用语句，触发对应计算或操作。
- **[行 585]** `    };`：作用域边界（代码块开始/结束）。
- **[行 586]** ``：空行，用于分隔逻辑块。
- **[行 587]** `    void TimeInfo::declare_parameters(ParameterHandler &prm)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 588]** `    {`：作用域边界（代码块开始/结束）。
- **[行 589]** `      prm.enter_subsection("Time");`：函数调用语句，触发对应计算或操作。
- **[行 590]** `      {`：作用域边界（代码块开始/结束）。
- **[行 591]** `        prm.declare_entry("End time", "1", Patterns::Double(), "End time");`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 592]** ``：空行，用于分隔逻辑块。
- **[行 593]** `        prm.declare_entry("Time data file",`：注册输入参数条目（名称、默认值、校验规则、说明）。
- **[行 594]** `                          "1",`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 595]** `                          Patterns::FileName(Patterns::FileName::input),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 596]** `                          "Time data file");`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 597]** `      }`：作用域边界（代码块开始/结束）。
- **[行 598]** `      prm.leave_subsection();`：函数调用语句，触发对应计算或操作。
- **[行 599]** `    }`：作用域边界（代码块开始/结束）。
- **[行 600]** ``：空行，用于分隔逻辑块。
- **[行 601]** `    void TimeInfo::parse_parameters(ParameterHandler &prm)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 602]** `    {`：作用域边界（代码块开始/结束）。
- **[行 603]** `      prm.enter_subsection("Time");`：函数调用语句，触发对应计算或操作。
- **[行 604]** `      {`：作用域边界（代码块开始/结束）。
- **[行 605]** `        m_end_time = prm.get_double("End time");`：从参数文件读取并写入成员变量。
- **[行 606]** `        m_time_file_name = prm.get("Time data file");`：从参数文件读取并写入成员变量。
- **[行 607]** `      }`：作用域边界（代码块开始/结束）。
- **[行 608]** `      prm.leave_subsection();`：函数调用语句，触发对应计算或操作。
- **[行 609]** `    }`：作用域边界（代码块开始/结束）。
- **[行 610]** ``：空行，用于分隔逻辑块。
- **[行 611]** `    struct AllParameters : public Scenario,`：类型声明（类/结构体），封装数据与行为。
- **[行 612]** `	                   public FESystem,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 613]** `	                   public BodyForce,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 614]** `			   public NonlinearSolver,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 615]** `			   public TimeInfo`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 616]** `    {`：作用域边界（代码块开始/结束）。
- **[行 617]** `      AllParameters(const std::string &input_file);`：调用 C++ 标准库工具函数/容器接口。
- **[行 618]** ``：空行，用于分隔逻辑块。
- **[行 619]** `      static void declare_parameters(ParameterHandler &prm);`：函数调用语句，触发对应计算或操作。
- **[行 620]** ``：空行，用于分隔逻辑块。
- **[行 621]** `      void parse_parameters(ParameterHandler &prm);`：函数调用语句，触发对应计算或操作。
- **[行 622]** `    };`：作用域边界（代码块开始/结束）。
- **[行 623]** ``：空行，用于分隔逻辑块。
- **[行 624]** `    AllParameters::AllParameters(const std::string &input_file)`：调用 C++ 标准库工具函数/容器接口。
- **[行 625]** `    {`：作用域边界（代码块开始/结束）。
- **[行 626]** `      ParameterHandler prm;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 627]** `      declare_parameters(prm);`：函数调用语句，触发对应计算或操作。
- **[行 628]** `      prm.parse_input(input_file);`：函数调用语句，触发对应计算或操作。
- **[行 629]** `      parse_parameters(prm);`：函数调用语句，触发对应计算或操作。
- **[行 630]** `    }`：作用域边界（代码块开始/结束）。
- **[行 631]** ``：空行，用于分隔逻辑块。
- **[行 632]** `    void AllParameters::declare_parameters(ParameterHandler &prm)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 633]** `    {`：作用域边界（代码块开始/结束）。
- **[行 634]** `      Scenario::declare_parameters(prm);`：函数调用语句，触发对应计算或操作。
- **[行 635]** `      FESystem::declare_parameters(prm);`：函数调用语句，触发对应计算或操作。
- **[行 636]** `      BodyForce::declare_parameters(prm);`：函数调用语句，触发对应计算或操作。
- **[行 637]** `      NonlinearSolver::declare_parameters(prm);`：函数调用语句，触发对应计算或操作。
- **[行 638]** `      TimeInfo::declare_parameters(prm);`：函数调用语句，触发对应计算或操作。
- **[行 639]** `    }`：作用域边界（代码块开始/结束）。
- **[行 640]** ``：空行，用于分隔逻辑块。
- **[行 641]** `    void AllParameters::parse_parameters(ParameterHandler &prm)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 642]** `    {`：作用域边界（代码块开始/结束）。
- **[行 643]** `      Scenario::parse_parameters(prm);`：函数调用语句，触发对应计算或操作。
- **[行 644]** `      FESystem::parse_parameters(prm);`：函数调用语句，触发对应计算或操作。
- **[行 645]** `      BodyForce::parse_parameters(prm);`：函数调用语句，触发对应计算或操作。
- **[行 646]** `      NonlinearSolver::parse_parameters(prm);`：函数调用语句，触发对应计算或操作。
- **[行 647]** `      TimeInfo::parse_parameters(prm);`：函数调用语句，触发对应计算或操作。
- **[行 648]** `    }`：作用域边界（代码块开始/结束）。
- **[行 649]** `  } // namespace Parameters`：作用域边界（代码块开始/结束）。
- **[行 650]** ``：空行，用于分隔逻辑块。
- **[行 651]** `  class Time`：类型声明（类/结构体），封装数据与行为。
- **[行 652]** `  {`：作用域边界（代码块开始/结束）。
- **[行 653]** `  public:`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 654]** `    Time(const double time_end)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 655]** `      : m_timestep(0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 656]** `      , m_time_current(0.0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 657]** `      , m_time_end(time_end)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 658]** `      , m_delta_t(0.0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 659]** `      , m_magnitude(1.0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 660]** `    {}`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 661]** ``：空行，用于分隔逻辑块。
- **[行 662]** `    virtual ~Time() = default;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 663]** ``：空行，用于分隔逻辑块。
- **[行 664]** `    double current() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 665]** `    {`：作用域边界（代码块开始/结束）。
- **[行 666]** `      return m_time_current;`：返回当前函数结果。
- **[行 667]** `    }`：作用域边界（代码块开始/结束）。
- **[行 668]** `    double end() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 669]** `    {`：作用域边界（代码块开始/结束）。
- **[行 670]** `      return m_time_end;`：返回当前函数结果。
- **[行 671]** `    }`：作用域边界（代码块开始/结束）。
- **[行 672]** `    double get_delta_t() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 673]** `    {`：作用域边界（代码块开始/结束）。
- **[行 674]** `      return m_delta_t;`：返回当前函数结果。
- **[行 675]** `    }`：作用域边界（代码块开始/结束）。
- **[行 676]** `    double get_magnitude() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 677]** `    {`：作用域边界（代码块开始/结束）。
- **[行 678]** `      return m_magnitude;`：返回当前函数结果。
- **[行 679]** `    }`：作用域边界（代码块开始/结束）。
- **[行 680]** `    unsigned int get_timestep() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 681]** `    {`：作用域边界（代码块开始/结束）。
- **[行 682]** `      return m_timestep;`：返回当前函数结果。
- **[行 683]** `    }`：作用域边界（代码块开始/结束）。
- **[行 684]** `    void increment(std::vector<std::array<double, 4>> time_table)`：调用 C++ 标准库工具函数/容器接口。
- **[行 685]** `    {`：作用域边界（代码块开始/结束）。
- **[行 686]** `      double t_1, t_delta, t_magnitude;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 687]** `      for (auto & time_group : time_table)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 688]** `        {`：作用域边界（代码块开始/结束）。
- **[行 689]** `	  t_1 = time_group[1];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 690]** `	  t_delta = time_group[2];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 691]** `	  t_magnitude = time_group[3];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 692]** ``：空行，用于分隔逻辑块。
- **[行 693]** `	  if (m_time_current < t_1 - 1.0e-6*t_delta)`：条件分支：根据当前状态选择执行路径。
- **[行 694]** `	    {`：作用域边界（代码块开始/结束）。
- **[行 695]** `	      m_delta_t = t_delta;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 696]** `	      m_magnitude = t_magnitude;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 697]** `	      break;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 698]** `	    }`：作用域边界（代码块开始/结束）。
- **[行 699]** `        }`：作用域边界（代码块开始/结束）。
- **[行 700]** ``：空行，用于分隔逻辑块。
- **[行 701]** `      m_time_current += m_delta_t;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 702]** `      ++m_timestep;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 703]** `    }`：作用域边界（代码块开始/结束）。
- **[行 704]** ``：空行，用于分隔逻辑块。
- **[行 705]** `  private:`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 706]** `    unsigned int m_timestep;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 707]** `    double       m_time_current;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 708]** `    const double m_time_end;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 709]** `    double m_delta_t;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 710]** `    double m_magnitude;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 711]** `  };`：作用域边界（代码块开始/结束）。
- **[行 712]** ``：空行，用于分隔逻辑块。
- **[行 713]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 714]** `  class LinearIsotropicElasticityAdditiveSplit`：类型声明（类/结构体），封装数据与行为。
- **[行 715]** `  {`：作用域边界（代码块开始/结束）。
- **[行 716]** `  public:`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 717]** `    LinearIsotropicElasticityAdditiveSplit(const double lame_lambda,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 718]** `			                   const double lame_mu,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 719]** `				           const double residual_k,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 720]** `					   const double length_scale,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 721]** `					   const double viscosity,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 722]** `					   const double gc,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 723]** `					   const bool   plane_stress_flag)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 724]** `      : m_lame_lambda(lame_lambda)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 725]** `      , m_lame_mu(lame_mu)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 726]** `      , m_residual_k(residual_k)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 727]** `      , m_length_scale(length_scale)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 728]** `      , m_eta(viscosity)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 729]** `      , m_gc(gc)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 730]** `      , m_plane_stress(plane_stress_flag)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 731]** `      , m_phase_field_value(0.0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 732]** `      , m_grad_phasefield(Tensor<1, dim>())`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 733]** `      , m_strain(SymmetricTensor<2, dim>())`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 734]** `      , m_stress(SymmetricTensor<2, dim>())`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 735]** `      , m_stress_positive(SymmetricTensor<2, dim>())`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 736]** `      , m_mechanical_C(SymmetricTensor<4, dim>())`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 737]** `      , m_strain_energy_positive(0.0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 738]** `      , m_strain_energy_negative(0.0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 739]** `      , m_strain_energy_total(0.0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 740]** `      , m_crack_energy_dissipation(0.0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 741]** `    {`：作用域边界（代码块开始/结束）。
- **[行 742]** `      Assert(  ( lame_lambda / (2*(lame_lambda + lame_mu)) <= 0.5)`：运行期断言/检查，验证输入与状态合法性。
- **[行 743]** `	     & ( lame_lambda / (2*(lame_lambda + lame_mu)) >=-1.0),`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 744]** `	     ExcInternalError() );`：函数调用语句，触发对应计算或操作。
- **[行 745]** `    }`：作用域边界（代码块开始/结束）。
- **[行 746]** ``：空行，用于分隔逻辑块。
- **[行 747]** `    const SymmetricTensor<4, dim> & get_mechanical_C() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 748]** `    {`：作用域边界（代码块开始/结束）。
- **[行 749]** `      return m_mechanical_C;`：返回当前函数结果。
- **[行 750]** `    }`：作用域边界（代码块开始/结束）。
- **[行 751]** ``：空行，用于分隔逻辑块。
- **[行 752]** `    const SymmetricTensor<2, dim> & get_cauchy_stress() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 753]** `    {`：作用域边界（代码块开始/结束）。
- **[行 754]** `      return m_stress;`：返回当前函数结果。
- **[行 755]** `    }`：作用域边界（代码块开始/结束）。
- **[行 756]** ``：空行，用于分隔逻辑块。
- **[行 757]** `    const SymmetricTensor<2, dim> & get_cauchy_stress_positive() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 758]** `    {`：作用域边界（代码块开始/结束）。
- **[行 759]** `      return m_stress_positive;`：返回当前函数结果。
- **[行 760]** `    }`：作用域边界（代码块开始/结束）。
- **[行 761]** ``：空行，用于分隔逻辑块。
- **[行 762]** `    double get_positive_strain_energy() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 763]** `    {`：作用域边界（代码块开始/结束）。
- **[行 764]** `      return m_strain_energy_positive;`：返回当前函数结果。
- **[行 765]** `    }`：作用域边界（代码块开始/结束）。
- **[行 766]** ``：空行，用于分隔逻辑块。
- **[行 767]** `    double get_negative_strain_energy() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 768]** `    {`：作用域边界（代码块开始/结束）。
- **[行 769]** `      return m_strain_energy_negative;`：返回当前函数结果。
- **[行 770]** `    }`：作用域边界（代码块开始/结束）。
- **[行 771]** ``：空行，用于分隔逻辑块。
- **[行 772]** `    double get_total_strain_energy() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 773]** `    {`：作用域边界（代码块开始/结束）。
- **[行 774]** `      return m_strain_energy_total;`：返回当前函数结果。
- **[行 775]** `    }`：作用域边界（代码块开始/结束）。
- **[行 776]** ``：空行，用于分隔逻辑块。
- **[行 777]** `    double get_crack_energy_dissipation() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 778]** `    {`：作用域边界（代码块开始/结束）。
- **[行 779]** `      return m_crack_energy_dissipation;`：返回当前函数结果。
- **[行 780]** `    }`：作用域边界（代码块开始/结束）。
- **[行 781]** ``：空行，用于分隔逻辑块。
- **[行 782]** `    double get_phase_field_value() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 783]** `    {`：作用域边界（代码块开始/结束）。
- **[行 784]** `      return m_phase_field_value;`：返回当前函数结果。
- **[行 785]** `    }`：作用域边界（代码块开始/结束）。
- **[行 786]** ``：空行，用于分隔逻辑块。
- **[行 787]** `    const Tensor<1, dim> get_phase_field_gradient() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 788]** `    {`：作用域边界（代码块开始/结束）。
- **[行 789]** `      return m_grad_phasefield;`：返回当前函数结果。
- **[行 790]** `    }`：作用域边界（代码块开始/结束）。
- **[行 791]** ``：空行，用于分隔逻辑块。
- **[行 792]** `    void update_material_data(const SymmetricTensor<2, dim> & strain,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 793]** `			      const double phase_field_value,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 794]** `			      const Tensor<1, dim> & grad_phasefield,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 795]** `			      const double phase_field_value_previous_step,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 796]** `			      const double delta_time);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 797]** ``：空行，用于分隔逻辑块。
- **[行 798]** `  private:`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 799]** `    const double m_lame_lambda;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 800]** `    const double m_lame_mu;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 801]** `    const double m_residual_k;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 802]** `    const double m_length_scale;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 803]** `    const double m_eta;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 804]** `    const double m_gc;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 805]** `    const bool   m_plane_stress;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 806]** `    double m_phase_field_value;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 807]** `    Tensor<1, dim> m_grad_phasefield;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 808]** `    SymmetricTensor<2, dim> m_strain;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 809]** `    SymmetricTensor<2, dim> m_stress;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 810]** `    SymmetricTensor<2, dim> m_stress_positive;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 811]** `    SymmetricTensor<4, dim> m_mechanical_C;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 812]** `    double m_strain_energy_positive;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 813]** `    double m_strain_energy_negative;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 814]** `    double m_strain_energy_total;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 815]** `    double m_crack_energy_dissipation;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 816]** `  };`：作用域边界（代码块开始/结束）。
- **[行 817]** ``：空行，用于分隔逻辑块。
- **[行 818]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 819]** `  void LinearIsotropicElasticityAdditiveSplit<dim>::`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 820]** `   update_material_data(const SymmetricTensor<2, dim> & strain,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 821]** `			const double phase_field_value,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 822]** `			const Tensor<1, dim> & grad_phasefield,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 823]** `			const double phase_field_value_previous_step,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 824]** `			const double delta_time)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 825]** `  {`：作用域边界（代码块开始/结束）。
- **[行 826]** `    m_strain = strain;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 827]** `    m_phase_field_value = phase_field_value;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 828]** `    m_grad_phasefield = grad_phasefield;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 829]** `    Vector<double>              eigenvalues(dim);`：函数调用语句，触发对应计算或操作。
- **[行 830]** `    std::vector<Tensor<1, dim>> eigenvectors(dim);`：调用 C++ 标准库工具函数/容器接口。
- **[行 831]** `    usr_spectrum_decomposition::spectrum_decomposition<dim>(m_strain,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 832]** `  							      eigenvalues,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 833]** `  							      eigenvectors);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 834]** ``：空行，用于分隔逻辑块。
- **[行 835]** `    SymmetricTensor<2, dim> strain_positive, strain_negative;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 836]** `    strain_positive = usr_spectrum_decomposition::positive_tensor(eigenvalues, eigenvectors);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 837]** `    strain_negative = usr_spectrum_decomposition::negative_tensor(eigenvalues, eigenvectors);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 838]** ``：空行，用于分隔逻辑块。
- **[行 839]** `    SymmetricTensor<4, dim> projector_positive, projector_negative;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 840]** `    usr_spectrum_decomposition::positive_negative_projectors(eigenvalues,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 841]** `  							       eigenvectors,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 842]** `							       projector_positive,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 843]** `							       projector_negative);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 844]** ``：空行，用于分隔逻辑块。
- **[行 845]** `    SymmetricTensor<2, dim> stress_positive, stress_negative;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 846]** `    const double degradation = degradation_function(m_phase_field_value) + m_residual_k;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 847]** `    const double I_1 = trace(m_strain);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 848]** ``：空行，用于分隔逻辑块。
- **[行 849]** `    // 2D plane strain and 3D cases`：注释行，用于说明算法背景、假设或实现意图。
- **[行 850]** `    double my_lambda = m_lame_lambda;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 851]** ``：空行，用于分隔逻辑块。
- **[行 852]** `    // 2D plane stress case`：注释行，用于说明算法背景、假设或实现意图。
- **[行 853]** `    if (    dim == 2`：条件分支：根据当前状态选择执行路径。
- **[行 854]** `	   && m_plane_stress)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 855]** `      my_lambda = 2 * m_lame_mu * m_lame_lambda / (m_lame_lambda + 2 * m_lame_mu);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 856]** ``：空行，用于分隔逻辑块。
- **[行 857]** `    stress_positive = my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 858]** `                                    * Physics::Elasticity::StandardTensors<dim>::I`：注释行，用于说明算法背景、假设或实现意图。
- **[行 859]** `                    + 2 * m_lame_mu * strain_positive;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 860]** `    stress_negative = my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 861]** `                                    * Physics::Elasticity::StandardTensors<dim>::I`：注释行，用于说明算法背景、假设或实现意图。
- **[行 862]** `    		      + 2 * m_lame_mu * strain_negative;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 863]** ``：空行，用于分隔逻辑块。
- **[行 864]** `    m_stress = degradation * stress_positive + stress_negative;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 865]** `    m_stress_positive = stress_positive;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 866]** ``：空行，用于分隔逻辑块。
- **[行 867]** `    SymmetricTensor<4, dim> C_positive, C_negative;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 868]** `    C_positive = my_lambda * usr_spectrum_decomposition::heaviside_function(I_1)`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 869]** `                               * Physics::Elasticity::StandardTensors<dim>::IxI`：注释行，用于说明算法背景、假设或实现意图。
- **[行 870]** `		 + 2 * m_lame_mu * projector_positive;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 871]** `    C_negative = my_lambda * usr_spectrum_decomposition::heaviside_function(-I_1)`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 872]** `                               * Physics::Elasticity::StandardTensors<dim>::IxI`：注释行，用于说明算法背景、假设或实现意图。
- **[行 873]** `    		 + 2 * m_lame_mu * projector_negative;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 874]** `    m_mechanical_C = degradation * C_positive + C_negative;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 875]** ``：空行，用于分隔逻辑块。
- **[行 876]** `    m_strain_energy_positive = 0.5 * my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 877]** `                                                   * usr_spectrum_decomposition::positive_ramp_function(I_1)`：注释行，用于说明算法背景、假设或实现意图。
- **[行 878]** `                             + m_lame_mu * strain_positive * strain_positive;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 879]** ``：空行，用于分隔逻辑块。
- **[行 880]** `    m_strain_energy_negative = 0.5 * my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 881]** `                                                   * usr_spectrum_decomposition::negative_ramp_function(I_1)`：注释行，用于说明算法背景、假设或实现意图。
- **[行 882]** `                             + m_lame_mu * strain_negative * strain_negative;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 883]** ``：空行，用于分隔逻辑块。
- **[行 884]** `    m_strain_energy_total = degradation * m_strain_energy_positive + m_strain_energy_negative;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 885]** ``：空行，用于分隔逻辑块。
- **[行 886]** `    m_crack_energy_dissipation = m_gc * (  0.5 / m_length_scale * m_phase_field_value * m_phase_field_value`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 887]** `	                                   + 0.5 * m_length_scale * m_grad_phasefield * m_grad_phasefield)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 888]** `	                                   // the term due to viscosity regularization`：注释行，用于说明算法背景、假设或实现意图。
- **[行 889]** `	                                   + (m_phase_field_value - phase_field_value_previous_step)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 890]** `					   * (m_phase_field_value - phase_field_value_previous_step)`：注释行，用于说明算法背景、假设或实现意图。
- **[行 891]** `				           * 0.5 * m_eta / delta_time;`：注释行，用于说明算法背景、假设或实现意图。
- **[行 892]** `    //(void)delta_time;`：注释行，用于说明算法背景、假设或实现意图。
- **[行 893]** `    //(void)phase_field_value_previous_step;`：注释行，用于说明算法背景、假设或实现意图。
- **[行 894]** `  }`：作用域边界（代码块开始/结束）。
- **[行 895]** ``：空行，用于分隔逻辑块。
- **[行 896]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 897]** `  class PointHistory`：类型声明（类/结构体），封装数据与行为。
- **[行 898]** `  {`：作用域边界（代码块开始/结束）。
- **[行 899]** `  public:`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 900]** `    PointHistory()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 901]** `      : m_length_scale(0.0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 902]** `      , m_gc(0.0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 903]** `      , m_viscosity(0.0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 904]** `    {}`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 905]** ``：空行，用于分隔逻辑块。
- **[行 906]** `    virtual ~PointHistory() = default;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 907]** ``：空行，用于分隔逻辑块。
- **[行 908]** `    void setup_lqp(const double lame_lambda,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 909]** `		   const double lame_mu,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 910]** `		   const double length_scale,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 911]** `		   const double gc,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 912]** `		   const double viscosity,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 913]** `		   const double residual_k,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 914]** `		   const bool   plane_stress_flag)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 915]** `    {`：作用域边界（代码块开始/结束）。
- **[行 916]** `      m_material =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 917]** `              std::make_shared<LinearIsotropicElasticityAdditiveSplit<dim>>(lame_lambda,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 918]** `        	                                                            lame_mu,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 919]** `								            residual_k,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 920]** `									    length_scale,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 921]** `									    viscosity,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 922]** `									    gc,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 923]** `									    plane_stress_flag);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 924]** `      m_length_scale = length_scale;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 925]** `      m_gc = gc;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 926]** `      m_viscosity = viscosity;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 927]** ``：空行，用于分隔逻辑块。
- **[行 928]** `      update_field_values(SymmetricTensor<2, dim>(), 0.0, Tensor<1, dim>(), 0.0, 1.0);`：函数调用语句，触发对应计算或操作。
- **[行 929]** `    }`：作用域边界（代码块开始/结束）。
- **[行 930]** ``：空行，用于分隔逻辑块。
- **[行 931]** `    void update_field_values(const SymmetricTensor<2, dim> & strain,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 932]** `		             const double phase_field_value,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 933]** `			     const Tensor<1, dim> & grad_phasefield,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 934]** `			     const double phase_field_value_previous_step,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 935]** `			     const double delta_time)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 936]** `    {`：作用域边界（代码块开始/结束）。
- **[行 937]** `      m_material->update_material_data(strain, phase_field_value, grad_phasefield,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 938]** `				       phase_field_value_previous_step, delta_time);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 939]** `    }`：作用域边界（代码块开始/结束）。
- **[行 940]** ``：空行，用于分隔逻辑块。
- **[行 941]** `    double get_current_positive_strain_energy() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 942]** `    {`：作用域边界（代码块开始/结束）。
- **[行 943]** `      return m_material->get_positive_strain_energy();`：返回当前函数结果。
- **[行 944]** `    }`：作用域边界（代码块开始/结束）。
- **[行 945]** ``：空行，用于分隔逻辑块。
- **[行 946]** `    const SymmetricTensor<4, dim> & get_mechanical_C() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 947]** `    {`：作用域边界（代码块开始/结束）。
- **[行 948]** `      return m_material->get_mechanical_C();`：返回当前函数结果。
- **[行 949]** `    }`：作用域边界（代码块开始/结束）。
- **[行 950]** ``：空行，用于分隔逻辑块。
- **[行 951]** `    const SymmetricTensor<2, dim> & get_cauchy_stress() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 952]** `    {`：作用域边界（代码块开始/结束）。
- **[行 953]** `      return m_material->get_cauchy_stress();`：返回当前函数结果。
- **[行 954]** `    }`：作用域边界（代码块开始/结束）。
- **[行 955]** ``：空行，用于分隔逻辑块。
- **[行 956]** `    const SymmetricTensor<2, dim> & get_cauchy_stress_positive() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 957]** `    {`：作用域边界（代码块开始/结束）。
- **[行 958]** `      return m_material->get_cauchy_stress_positive();`：返回当前函数结果。
- **[行 959]** `    }`：作用域边界（代码块开始/结束）。
- **[行 960]** ``：空行，用于分隔逻辑块。
- **[行 961]** `    double get_total_strain_energy() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 962]** `    {`：作用域边界（代码块开始/结束）。
- **[行 963]** `      return m_material->get_total_strain_energy();`：返回当前函数结果。
- **[行 964]** `    }`：作用域边界（代码块开始/结束）。
- **[行 965]** ``：空行，用于分隔逻辑块。
- **[行 966]** `    double get_crack_energy_dissipation() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 967]** `    {`：作用域边界（代码块开始/结束）。
- **[行 968]** `      return m_material->get_crack_energy_dissipation();`：返回当前函数结果。
- **[行 969]** `    }`：作用域边界（代码块开始/结束）。
- **[行 970]** ``：空行，用于分隔逻辑块。
- **[行 971]** `    double get_phase_field_value() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 972]** `    {`：作用域边界（代码块开始/结束）。
- **[行 973]** `      return m_material->get_phase_field_value();`：返回当前函数结果。
- **[行 974]** `    }`：作用域边界（代码块开始/结束）。
- **[行 975]** ``：空行，用于分隔逻辑块。
- **[行 976]** `    const Tensor<1, dim> get_phase_field_gradient() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 977]** `    {`：作用域边界（代码块开始/结束）。
- **[行 978]** `      return m_material->get_phase_field_gradient();`：返回当前函数结果。
- **[行 979]** `    }`：作用域边界（代码块开始/结束）。
- **[行 980]** ``：空行，用于分隔逻辑块。
- **[行 981]** `    double get_length_scale() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 982]** `    {`：作用域边界（代码块开始/结束）。
- **[行 983]** `      return m_length_scale;`：返回当前函数结果。
- **[行 984]** `    }`：作用域边界（代码块开始/结束）。
- **[行 985]** ``：空行，用于分隔逻辑块。
- **[行 986]** `    double get_critical_energy_release_rate() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 987]** `    {`：作用域边界（代码块开始/结束）。
- **[行 988]** `      return m_gc;`：返回当前函数结果。
- **[行 989]** `    }`：作用域边界（代码块开始/结束）。
- **[行 990]** ``：空行，用于分隔逻辑块。
- **[行 991]** `    double get_viscosity() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 992]** `    {`：作用域边界（代码块开始/结束）。
- **[行 993]** `      return m_viscosity;`：返回当前函数结果。
- **[行 994]** `    }`：作用域边界（代码块开始/结束）。
- **[行 995]** `  private:`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 996]** `    std::shared_ptr<LinearIsotropicElasticityAdditiveSplit<dim>> m_material;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 997]** `    double m_length_scale;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 998]** `    double m_gc;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 999]** `    double m_viscosity;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1000]** `  };`：作用域边界（代码块开始/结束）。
- **[行 1001]** ``：空行，用于分隔逻辑块。
- **[行 1002]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 1003]** `  class PhaseFieldMonolithicSolve`：类型声明（类/结构体），封装数据与行为。
- **[行 1004]** `  {`：作用域边界（代码块开始/结束）。
- **[行 1005]** `  public:`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1006]** `    PhaseFieldMonolithicSolve(const std::string &input_file);`：调用 C++ 标准库工具函数/容器接口。
- **[行 1007]** ``：空行，用于分隔逻辑块。
- **[行 1008]** `    virtual ~PhaseFieldMonolithicSolve() = default;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1009]** `    void run();`：函数调用语句，触发对应计算或操作。
- **[行 1010]** ``：空行，用于分隔逻辑块。
- **[行 1011]** `  private:`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1012]** `    struct PerTaskData_ASM;`：类型声明（类/结构体），封装数据与行为。
- **[行 1013]** `    struct ScratchData_ASM;`：类型声明（类/结构体），封装数据与行为。
- **[行 1014]** ``：空行，用于分隔逻辑块。
- **[行 1015]** `    struct PerTaskData_ASM_RHS_BFGS;`：类型声明（类/结构体），封装数据与行为。
- **[行 1016]** `    struct ScratchData_ASM_RHS_BFGS;`：类型声明（类/结构体），封装数据与行为。
- **[行 1017]** ``：空行，用于分隔逻辑块。
- **[行 1018]** `    struct PerTaskData_UQPH;`：类型声明（类/结构体），封装数据与行为。
- **[行 1019]** `    struct ScratchData_UQPH;`：类型声明（类/结构体），封装数据与行为。
- **[行 1020]** ``：空行，用于分隔逻辑块。
- **[行 1021]** `    Parameters::AllParameters m_parameters;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1022]** `    Triangulation<dim> m_triangulation;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1023]** ``：空行，用于分隔逻辑块。
- **[行 1024]** `    CellDataStorage<typename Triangulation<dim>::cell_iterator,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1025]** `                    PointHistory<dim>>`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1026]** `      m_quadrature_point_history;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1027]** ``：空行，用于分隔逻辑块。
- **[行 1028]** `    Time                m_time;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1029]** `    std::ofstream m_logfile;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1030]** `    mutable TimerOutput m_timer;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1031]** ``：空行，用于分隔逻辑块。
- **[行 1032]** `    DoFHandler<dim>                  m_dof_handler;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1033]** `    FESystem<dim>                    m_fe;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1034]** `    const unsigned int               m_dofs_per_cell;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1035]** `    const FEValuesExtractors::Vector m_u_fe;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1036]** `    const FEValuesExtractors::Scalar m_d_fe;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1037]** ``：空行，用于分隔逻辑块。
- **[行 1038]** `    static const unsigned int m_n_blocks          = 2;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1039]** `    static const unsigned int m_n_components      = dim + 1;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1040]** `    static const unsigned int m_first_u_component = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1041]** `    static const unsigned int m_d_component       = dim;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1042]** ``：空行，用于分隔逻辑块。
- **[行 1043]** `    enum`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1044]** `    {`：作用域边界（代码块开始/结束）。
- **[行 1045]** `      m_u_dof = 0,`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1046]** `      m_d_dof = 1`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1047]** `    };`：作用域边界（代码块开始/结束）。
- **[行 1048]** ``：空行，用于分隔逻辑块。
- **[行 1049]** `    std::vector<types::global_dof_index> m_dofs_per_block;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1050]** ``：空行，用于分隔逻辑块。
- **[行 1051]** `    const QGauss<dim>     m_qf_cell;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1052]** `    const QGauss<dim - 1> m_qf_face;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1053]** `    const unsigned int    m_n_q_points;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1054]** ``：空行，用于分隔逻辑块。
- **[行 1055]** `    double m_vol_reference;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1056]** ``：空行，用于分隔逻辑块。
- **[行 1057]** `    AffineConstraints<double> m_constraints;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1058]** `    BlockSparsityPattern      m_sparsity_pattern;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1059]** `    BlockSparseMatrix<double> m_tangent_matrix;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1060]** `    BlockVector<double>       m_system_rhs;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1061]** `    BlockVector<double>       m_solution;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1062]** `    SparseDirectUMFPACK       m_A_direct;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1063]** `    // m_active_set_phasefield has 0 (inactive constraint)`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1064]** `    //                          or 1 (active constraint lower bound)`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1065]** `    //                          or 2 (active constraint upper bound)`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1066]** `    // In order to add active set into the VTK output, we have to declare`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1067]** `    // it as double, not int or unsigned int`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1068]** `    Vector<double> m_active_set_phasefield;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1069]** ``：空行，用于分隔逻辑块。
- **[行 1070]** `    std::map<unsigned int, std::vector<double>> m_material_data;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1071]** ``：空行，用于分隔逻辑块。
- **[行 1072]** `    std::vector<std::pair<double, std::vector<double>>> m_history_reaction_force;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1073]** `    std::vector<std::pair<double, std::array<double, 3>>> m_history_energy;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1074]** ``：空行，用于分隔逻辑块。
- **[行 1075]** ``：空行，用于分隔逻辑块。
- **[行 1076]** `    struct Errors`：类型声明（类/结构体），封装数据与行为。
- **[行 1077]** `    {`：作用域边界（代码块开始/结束）。
- **[行 1078]** `      Errors()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1079]** `        : m_norm(1.0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1080]** `        , m_u(1.0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1081]** `        , m_d(1.0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1082]** `      {}`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1083]** ``：空行，用于分隔逻辑块。
- **[行 1084]** `      void reset()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1085]** `      {`：作用域边界（代码块开始/结束）。
- **[行 1086]** `        m_norm = 1.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1087]** `        m_u    = 1.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1088]** `        m_d    = 1.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1089]** `      }`：作用域边界（代码块开始/结束）。
- **[行 1090]** ``：空行，用于分隔逻辑块。
- **[行 1091]** `      void normalize(const Errors &rhs)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1092]** `      {`：作用域边界（代码块开始/结束）。
- **[行 1093]** `        if (rhs.m_norm != 0.0)`：条件分支：根据当前状态选择执行路径。
- **[行 1094]** `          m_norm /= rhs.m_norm;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1095]** `        if (rhs.m_u != 0.0)`：条件分支：根据当前状态选择执行路径。
- **[行 1096]** `          m_u /= rhs.m_u;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1097]** `        if (rhs.m_d != 0.0)`：条件分支：根据当前状态选择执行路径。
- **[行 1098]** `          m_d /= rhs.m_d;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1099]** `      }`：作用域边界（代码块开始/结束）。
- **[行 1100]** ``：空行，用于分隔逻辑块。
- **[行 1101]** `      double m_norm, m_u, m_d;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1102]** `    };`：作用域边界（代码块开始/结束）。
- **[行 1103]** ``：空行，用于分隔逻辑块。
- **[行 1104]** `    Errors m_error_residual, m_error_residual_0, m_error_residual_norm, m_error_update,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1105]** `      m_error_update_0, m_error_update_norm;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1106]** ``：空行，用于分隔逻辑块。
- **[行 1107]** `    void get_error_residual(Errors &error_residual);`：函数调用语句，触发对应计算或操作。
- **[行 1108]** `    void get_error_residual_LBFGSB(Errors &error_residual,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1109]** `				   const BlockVector<double> & solution_delta);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1110]** ``：空行，用于分隔逻辑块。
- **[行 1111]** `    void get_error_update(const BlockVector<double> &newton_update,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1112]** `                          Errors & error_update);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1113]** ``：空行，用于分隔逻辑块。
- **[行 1114]** `    void make_grid();`：函数调用语句，触发对应计算或操作。
- **[行 1115]** `    void make_grid_case_1();`：函数调用语句，触发对应计算或操作。
- **[行 1116]** `    void make_grid_case_2();`：函数调用语句，触发对应计算或操作。
- **[行 1117]** `    void make_grid_case_3();`：函数调用语句，触发对应计算或操作。
- **[行 1118]** `    void make_grid_case_4();`：函数调用语句，触发对应计算或操作。
- **[行 1119]** `    void make_grid_case_5();`：函数调用语句，触发对应计算或操作。
- **[行 1120]** `    void make_grid_case_6();`：函数调用语句，触发对应计算或操作。
- **[行 1121]** `    void make_grid_case_7();`：函数调用语句，触发对应计算或操作。
- **[行 1122]** `    void make_grid_case_8();`：函数调用语句，触发对应计算或操作。
- **[行 1123]** `    void make_grid_case_9();`：函数调用语句，触发对应计算或操作。
- **[行 1124]** `    void make_grid_case_10();`：函数调用语句，触发对应计算或操作。
- **[行 1125]** `    void make_grid_case_11();`：函数调用语句，触发对应计算或操作。
- **[行 1126]** ``：空行，用于分隔逻辑块。
- **[行 1127]** `    void setup_system();`：函数调用语句，触发对应计算或操作。
- **[行 1128]** ``：空行，用于分隔逻辑块。
- **[行 1129]** `    void determine_component_extractors();`：函数调用语句，触发对应计算或操作。
- **[行 1130]** ``：空行，用于分隔逻辑块。
- **[行 1131]** `    void make_constraints(const unsigned int it_nr);`：函数调用语句，触发对应计算或操作。
- **[行 1132]** ``：空行，用于分隔逻辑块。
- **[行 1133]** `    void assemble_system_B0(const BlockVector<double> & solution_old);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 1134]** ``：空行，用于分隔逻辑块。
- **[行 1135]** `    void assemble_system_B0_one_cell(`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 1136]** `      const typename DoFHandler<dim>::active_cell_iterator &cell,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1137]** `      ScratchData_ASM &                                     scratch,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1138]** `      PerTaskData_ASM &                                     data) const;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1139]** ``：空行，用于分隔逻辑块。
- **[行 1140]** `    void assemble_system_rhs_BFGS_one_cell(`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 1141]** `      const typename DoFHandler<dim>::active_cell_iterator &cell,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1142]** `      ScratchData_ASM_RHS_BFGS &                           scratch,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1143]** `      PerTaskData_ASM_RHS_BFGS &                           data) const;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1144]** ``：空行，用于分隔逻辑块。
- **[行 1145]** `    void assemble_system_rhs_BFGS(const BlockVector<double> & solution_old,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 1146]** `				  BlockVector<double> & system_rhs);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1147]** ``：空行，用于分隔逻辑块。
- **[行 1148]** `    void assemble_system_rhs_BFGS_parallel(const BlockVector<double> & solution_old,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 1149]** `    				           BlockVector<double> & system_rhs);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1150]** ``：空行，用于分隔逻辑块。
- **[行 1151]** `    void solve_nonlinear_timestep_LBFGS(BlockVector<double> &solution_delta,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1152]** `				        BlockVector<double> & LBFGS_update_refine);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1153]** ``：空行，用于分隔逻辑块。
- **[行 1154]** `    void solve_nonlinear_timestep_LBFGS_B(BlockVector<double> &solution_delta,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1155]** `    				          BlockVector<double> & LBFGS_update_refine);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1156]** ``：空行，用于分隔逻辑块。
- **[行 1157]** `    void calculate_cauchy_point(const BlockSparseMatrix<double> & B0_matrix,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 1158]** `	                        const std::list<BlockVector<double>> & y_vector_list,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1159]** `				const std::list<BlockVector<double>> & b0xs_vector_list,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1160]** `				const FullMatrix<double> & M_matrix,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1161]** `				const BlockVector<double> & gradient_g,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1162]** `				const BlockVector<double> & solution_delta,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1163]** `				BlockVector<double> & solution_delta_cauchy_point);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1164]** ``：空行，用于分隔逻辑块。
- **[行 1165]** `    double line_search_stepsize_gradient_based(const BlockVector<double> & BFGS_p_vector,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 1166]** `					       const BlockVector<double> & solution_delta);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1167]** ``：空行，用于分隔逻辑块。
- **[行 1168]** `    double line_search_stepsize_strong_wolfe(const double phi_0,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 1169]** `				             const double phi_0_prime,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1170]** `				             const BlockVector<double> & BFGS_p_vector,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1171]** `				             const BlockVector<double> & solution_delta);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1172]** ``：空行，用于分隔逻辑块。
- **[行 1173]** `    double line_search_zoom_strong_wolfe(double phi_low, double phi_low_prime, double alpha_low,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 1174]** `					 double phi_high, double phi_high_prime, double alpha_high,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1175]** `					 double phi_0, double phi_0_prime, const BlockVector<double> & BFGS_p_vector,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1176]** `					 double c1, double c2, unsigned int max_iter,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1177]** `					 const BlockVector<double> & solution_delta);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1178]** ``：空行，用于分隔逻辑块。
- **[行 1179]** `    double line_search_interpolation_cubic(const double alpha_0, const double phi_0, const double phi_0_prime,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 1180]** `					   const double alpha_1, const double phi_1, const double phi_1_prime);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1181]** ``：空行，用于分隔逻辑块。
- **[行 1182]** `    std::pair<double, double> calculate_phi_and_phi_prime(const double alpha,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 1183]** `							  const BlockVector<double> & BFGS_p_vector,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1184]** `							  const BlockVector<double> & solution_delta);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1185]** ``：空行，用于分隔逻辑块。
- **[行 1186]** `    void LBFGS_B0(BlockVector<double> & LBFGS_r_vector,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1187]** `		  BlockVector<double> & LBFGS_q_vector);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1188]** ``：空行，用于分隔逻辑块。
- **[行 1189]** `    void output_results() const;`：函数调用语句，触发对应计算或操作。
- **[行 1190]** ``：空行，用于分隔逻辑块。
- **[行 1191]** `    void setup_qph();`：函数调用语句，触发对应计算或操作。
- **[行 1192]** ``：空行，用于分隔逻辑块。
- **[行 1193]** `    void update_qph_incremental(const BlockVector<double> &solution_delta,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 1194]** `				const BlockVector<double> &solution_old);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1195]** ``：空行，用于分隔逻辑块。
- **[行 1196]** `    void update_qph_incremental_one_cell(`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 1197]** `      const typename DoFHandler<dim>::active_cell_iterator &cell,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1198]** `      ScratchData_UQPH &                                    scratch,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1199]** `      PerTaskData_UQPH &                                    data);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1200]** ``：空行，用于分隔逻辑块。
- **[行 1201]** `    void copy_local_to_global_UQPH(const PerTaskData_UQPH & /*data*/)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1202]** `    {}`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1203]** ``：空行，用于分隔逻辑块。
- **[行 1204]** `    BlockVector<double>`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1205]** `    get_total_solution(const BlockVector<double> &solution_delta) const;`：函数调用语句，触发对应计算或操作。
- **[行 1206]** ``：空行，用于分隔逻辑块。
- **[行 1207]** `    // Should not make this function const`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1208]** `    void read_material_data(const std::string &data_file,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1209]** `			    const unsigned int total_material_regions);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1210]** ``：空行，用于分隔逻辑块。
- **[行 1211]** `    void read_time_data(const std::string &data_file,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1212]** `    		        std::vector<std::array<double, 4>> & time_table);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1213]** ``：空行，用于分隔逻辑块。
- **[行 1214]** `    void print_conv_header_LBFGS();`：函数调用语句，触发对应计算或操作。
- **[行 1215]** ``：空行，用于分隔逻辑块。
- **[行 1216]** `    void print_conv_header_LBFGSB();`：函数调用语句，触发对应计算或操作。
- **[行 1217]** ``：空行，用于分隔逻辑块。
- **[行 1218]** `    void print_parameter_information();`：函数调用语句，触发对应计算或操作。
- **[行 1219]** ``：空行，用于分隔逻辑块。
- **[行 1220]** `    void calculate_reaction_force(unsigned int face_ID);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 1221]** ``：空行，用于分隔逻辑块。
- **[行 1222]** `    void write_history_data();`：函数调用语句，触发对应计算或操作。
- **[行 1223]** ``：空行，用于分隔逻辑块。
- **[行 1224]** `    double calculate_energy_functional() const;`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 1225]** ``：空行，用于分隔逻辑块。
- **[行 1226]** `    std::pair<double, double> calculate_total_strain_energy_and_crack_energy_dissipation() const;`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 1227]** ``：空行，用于分隔逻辑块。
- **[行 1228]** `    bool local_refine_and_solution_transfer(BlockVector<double> & solution_delta,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1229]** `					    BlockVector<double> & LBFGS_update_refine);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1230]** ``：空行，用于分隔逻辑块。
- **[行 1231]** `    // L-BFGS-B subroutines`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1232]** `    void point_projection(BlockVector<double> & solution_delta);`：函数调用语句，触发对应计算或操作。
- **[行 1233]** ``：空行，用于分隔逻辑块。
- **[行 1234]** `    std::priority_queue< std::pair<double, unsigned int>,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1235]** `                         std::vector<std::pair<double, unsigned int>>,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1236]** `    		         std::greater<std::pair<double, unsigned int>> >`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1237]** `      calculate_break_points(const BlockVector<double> & solution_delta,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 1238]** `			     const BlockVector<double> & gradient_g,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1239]** `			     BlockVector<double> & gradient_d);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1240]** ``：空行，用于分隔逻辑块。
- **[行 1241]** `    double ebT_x_B0_x_v(const unsigned int b,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1242]** `			const BlockSparseMatrix<double> & B0_matrix,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1243]** `			const BlockVector<double> & v);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1244]** ``：空行，用于分隔逻辑块。
- **[行 1245]** `    void zT_x_vector(const BlockVector<double> & z,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1246]** `		     const BlockVector<double> & src_vector,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1247]** `		     BlockVector<double> & target_vector);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1248]** ``：空行，用于分隔逻辑块。
- **[行 1249]** `    void z_x_vector(const BlockVector<double> & z,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1250]** `		    const BlockVector<double> & src_vector,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1251]** `		    BlockVector<double> & target_vector);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1252]** ``：空行，用于分隔逻辑块。
- **[行 1253]** `    void zT_B0_z(const BlockVector<double> & z,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1254]** `		     BlockSparseMatrix<double> & B0_matrix);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1255]** ``：空行，用于分隔逻辑块。
- **[行 1256]** ``：空行，用于分隔逻辑块。
- **[行 1257]** `  }; // class PhaseFieldSplitSolve`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1258]** ``：空行，用于分隔逻辑块。
- **[行 1259]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 1260]** `  void PhaseFieldMonolithicSolve<dim>::zT_B0_z(const BlockVector<double> & z,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1261]** `					       BlockSparseMatrix<double> & B0_matrix)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1262]** `  {`：作用域边界（代码块开始/结束）。
- **[行 1263]** `    // block 1: displacement`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1264]** `    for (unsigned int i = 0; i < z.block(m_u_dof).size(); ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 1265]** `      {`：作用域边界（代码块开始/结束）。
- **[行 1266]** `	if (z.block(m_u_dof)[i] < 0)`：条件分支：根据当前状态选择执行路径。
- **[行 1267]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 1268]** `	    for (auto itr = B0_matrix.block(m_u_dof, m_u_dof).begin(i);`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 1269]** `		      itr != B0_matrix.block(m_u_dof, m_u_dof).end(i);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1270]** `		      ++itr)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1271]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 1272]** `		if (itr->column() != itr->row())`：条件分支：根据当前状态选择执行路径。
- **[行 1273]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 1274]** `		    B0_matrix.block(m_u_dof, m_u_dof).set(itr->row(),    itr->column(), 0.0);`：函数调用语句，触发对应计算或操作。
- **[行 1275]** `		    B0_matrix.block(m_u_dof, m_u_dof).set(itr->column(), itr->row(),    0.0);`：函数调用语句，触发对应计算或操作。
- **[行 1276]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 1277]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 1278]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 1279]** `      } // for (unsigned int i = 0; i < z.block(m_u_dof).size(); ++i)`：作用域边界（代码块开始/结束）。
- **[行 1280]** ``：空行，用于分隔逻辑块。
- **[行 1281]** `    // block 2: phasefield`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1282]** `    for (unsigned int i = 0; i < z.block(m_d_dof).size(); ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 1283]** `      {`：作用域边界（代码块开始/结束）。
- **[行 1284]** `	if (z.block(m_d_dof)[i] < 0)`：条件分支：根据当前状态选择执行路径。
- **[行 1285]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 1286]** `	    for (auto itr = B0_matrix.block(m_d_dof, m_d_dof).begin(i);`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 1287]** `		      itr != B0_matrix.block(m_d_dof, m_d_dof).end(i);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1288]** `		      ++itr)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1289]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 1290]** `		if (itr->column() != itr->row())`：条件分支：根据当前状态选择执行路径。
- **[行 1291]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 1292]** `		    B0_matrix.block(m_d_dof, m_d_dof).set(itr->row(),    itr->column(), 0.0);`：函数调用语句，触发对应计算或操作。
- **[行 1293]** `		    B0_matrix.block(m_d_dof, m_d_dof).set(itr->column(), itr->row(),    0.0);`：函数调用语句，触发对应计算或操作。
- **[行 1294]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 1295]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 1296]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 1297]** `      } // for (unsigned int i = 0; i < z.block(m_d_dof).size(); ++i)`：作用域边界（代码块开始/结束）。
- **[行 1298]** `  }`：作用域边界（代码块开始/结束）。
- **[行 1299]** ``：空行，用于分隔逻辑块。
- **[行 1300]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 1301]** `  void PhaseFieldMonolithicSolve<dim>::z_x_vector(const BlockVector<double> & z,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1302]** `						  const BlockVector<double> & src_vector,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1303]** `						  BlockVector<double> & target_vector)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1304]** `  {`：作用域边界（代码块开始/结束）。
- **[行 1305]** `    //We assume that the dimensions of all the block matrices are correct`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1306]** `    //block 1: displacement`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1307]** `    unsigned int target_vector_index = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1308]** `    for (unsigned int i = 0; i < z.block(m_u_dof).size(); ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 1309]** `      {`：作用域边界（代码块开始/结束）。
- **[行 1310]** `	if (z.block(m_u_dof)[i] > 0)`：条件分支：根据当前状态选择执行路径。
- **[行 1311]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 1312]** `	    target_vector.block(m_u_dof)[i] = src_vector.block(m_u_dof)[target_vector_index];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1313]** `	    ++target_vector_index;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1314]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 1315]** `	else`：条件分支的兜底路径。
- **[行 1316]** `	  target_vector.block(m_u_dof)[i] = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1317]** `      }`：作用域边界（代码块开始/结束）。
- **[行 1318]** ``：空行，用于分隔逻辑块。
- **[行 1319]** `    //block 2: phasefield`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1320]** `    target_vector_index = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1321]** `    for (unsigned int i = 0; i < z.block(m_d_dof).size(); ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 1322]** `      {`：作用域边界（代码块开始/结束）。
- **[行 1323]** `	if (z.block(m_d_dof)[i] > 0)`：条件分支：根据当前状态选择执行路径。
- **[行 1324]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 1325]** `	    target_vector.block(m_d_dof)[i] = src_vector.block(m_d_dof)[target_vector_index];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1326]** `	    ++target_vector_index;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1327]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 1328]** `	else`：条件分支的兜底路径。
- **[行 1329]** `	  target_vector.block(m_d_dof)[i] = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1330]** `      }`：作用域边界（代码块开始/结束）。
- **[行 1331]** `  }`：作用域边界（代码块开始/结束）。
- **[行 1332]** ``：空行，用于分隔逻辑块。
- **[行 1333]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 1334]** `  void PhaseFieldMonolithicSolve<dim>::zT_x_vector(const BlockVector<double> & z,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1335]** `						   const BlockVector<double> & src_vector,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1336]** `						   BlockVector<double> & target_vector)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1337]** `  {`：作用域边界（代码块开始/结束）。
- **[行 1338]** `    //We assume that the dimensions of all the block matrices are correct`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1339]** `    //block 1: displacement`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1340]** `    unsigned int target_vector_index = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1341]** `    for (unsigned int i = 0; i < z.block(m_u_dof).size(); ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 1342]** `      {`：作用域边界（代码块开始/结束）。
- **[行 1343]** `	if (z.block(m_u_dof)[i] > 0)`：条件分支：根据当前状态选择执行路径。
- **[行 1344]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 1345]** `	    target_vector.block(m_u_dof)[target_vector_index] = src_vector.block(m_u_dof)[i];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1346]** `	    ++target_vector_index;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1347]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 1348]** `      }`：作用域边界（代码块开始/结束）。
- **[行 1349]** ``：空行，用于分隔逻辑块。
- **[行 1350]** `    //block 2: phasefield`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1351]** `    target_vector_index = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1352]** `    for (unsigned int i = 0; i < z.block(m_d_dof).size(); ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 1353]** `      {`：作用域边界（代码块开始/结束）。
- **[行 1354]** `        if (z.block(m_d_dof)[i] > 0)`：条件分支：根据当前状态选择执行路径。
- **[行 1355]** `        {`：作用域边界（代码块开始/结束）。
- **[行 1356]** `	  target_vector.block(m_d_dof)[target_vector_index] = src_vector.block(m_d_dof)[i];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1357]** `	  ++target_vector_index;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1358]** `        }`：作用域边界（代码块开始/结束）。
- **[行 1359]** `      }`：作用域边界（代码块开始/结束）。
- **[行 1360]** `  }`：作用域边界（代码块开始/结束）。
- **[行 1361]** ``：空行，用于分隔逻辑块。
- **[行 1362]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 1363]** `  double PhaseFieldMonolithicSolve<dim>::ebT_x_B0_x_v(const unsigned int b,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1364]** `						      const BlockSparseMatrix<double> & B0_matrix,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1365]** `						      const BlockVector<double> & v)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1366]** `  {`：作用域边界（代码块开始/结束）。
- **[行 1367]** `    double row_sum = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1368]** `    for (auto itr = B0_matrix.block(m_d_dof, m_d_dof).begin(b);`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 1369]** `              itr != B0_matrix.block(m_d_dof, m_d_dof).end(b);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1370]** `              ++itr)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1371]** `      {`：作用域边界（代码块开始/结束）。
- **[行 1372]** `        row_sum += itr->value() * v.block(m_d_dof)(itr->column());`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1373]** `      }`：作用域边界（代码块开始/结束）。
- **[行 1374]** ``：空行，用于分隔逻辑块。
- **[行 1375]** `    return row_sum;`：返回当前函数结果。
- **[行 1376]** `  }`：作用域边界（代码块开始/结束）。
- **[行 1377]** ``：空行，用于分隔逻辑块。
- **[行 1378]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 1379]** `  void PhaseFieldMonolithicSolve<dim>::point_projection(BlockVector<double> & solution_delta)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1380]** `  {`：作用域边界（代码块开始/结束）。
- **[行 1381]** `    // Phase-field value cannot exceed 1.0`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1382]** `    const double upper_limit = 1.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1383]** ``：空行，用于分隔逻辑块。
- **[行 1384]** `    Vector<double> solution_phasefield_total(m_solution.block(m_d_dof));`：函数调用语句，触发对应计算或操作。
- **[行 1385]** `    solution_phasefield_total += solution_delta.block(m_d_dof);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1386]** ``：空行，用于分隔逻辑块。
- **[行 1387]** `    for (unsigned int i = 0; i < solution_phasefield_total.size(); ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 1388]** `      {`：作用域边界（代码块开始/结束）。
- **[行 1389]** `	if (solution_delta.block(m_d_dof)[i] < 0.0)`：条件分支：根据当前状态选择执行路径。
- **[行 1390]** `	  solution_delta.block(m_d_dof)[i] = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1391]** ``：空行，用于分隔逻辑块。
- **[行 1392]** `	if (solution_phasefield_total[i] > upper_limit)`：条件分支：根据当前状态选择执行路径。
- **[行 1393]** `	  solution_delta.block(m_d_dof)[i] = upper_limit - m_solution.block(m_d_dof)[i];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1394]** `      }`：作用域边界（代码块开始/结束）。
- **[行 1395]** `  }`：作用域边界（代码块开始/结束）。
- **[行 1396]** ``：空行，用于分隔逻辑块。
- **[行 1397]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 1398]** `  std::priority_queue< std::pair<double, unsigned int>,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1399]** `                       std::vector<std::pair<double, unsigned int>>,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1400]** `  		       std::greater<std::pair<double, unsigned int>> >`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1401]** `    PhaseFieldMonolithicSolve<dim>::calculate_break_points(const BlockVector<double> & solution_delta,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 1402]** `		       				           const BlockVector<double> & gradient_g,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1403]** `							   BlockVector<double> & gradient_d)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1404]** `  {`：作用域边界（代码块开始/结束）。
- **[行 1405]** `    // Creates a min heap of break points`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1406]** `    std::priority_queue< std::pair<double, unsigned int>,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1407]** `                         std::vector<std::pair<double, unsigned int>>,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1408]** `    		         std::greater<std::pair<double, unsigned int>> >`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1409]** `    break_points_sorted;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1410]** ``：空行，用于分隔逻辑块。
- **[行 1411]** `    double t = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1412]** ``：空行，用于分隔逻辑块。
- **[行 1413]** `    Vector<double> solution_phasefield_total(m_solution.block(m_d_dof));`：函数调用语句，触发对应计算或操作。
- **[行 1414]** `    solution_phasefield_total += solution_delta.block(m_d_dof);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1415]** ``：空行，用于分隔逻辑块。
- **[行 1416]** `    // upper bound is 1.0, lower bound is the solution at the previous step.`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1417]** `    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 1418]** `      {`：作用域边界（代码块开始/结束）。
- **[行 1419]** `	if (gradient_g.block(m_d_dof)[i] < 0)`：条件分支：根据当前状态选择执行路径。
- **[行 1420]** `	  t = (solution_phasefield_total[i] - 1.0 ) / gradient_g.block(m_d_dof)[i];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1421]** `	else if (gradient_g.block(m_d_dof)[i] > 0)`：多分支条件判断，处理备选情形。
- **[行 1422]** `	  t = solution_delta.block(m_d_dof)[i] / gradient_g.block(m_d_dof)[i];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1423]** `	else`：条件分支的兜底路径。
- **[行 1424]** `	  t = std::numeric_limits<double>::max();`：调用 C++ 标准库工具函数/容器接口。
- **[行 1425]** ``：空行，用于分隔逻辑块。
- **[行 1426]** `        //AssertThrow(t >= 0, ExcMessage("Break point has to be a non-negative t value"));`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1427]** ``：空行，用于分隔逻辑块。
- **[行 1428]** `        if (t > 0)`：条件分支：根据当前状态选择执行路径。
- **[行 1429]** `          {`：作用域边界（代码块开始/结束）。
- **[行 1430]** `	    break_points_sorted.push(std::make_pair(t, i));`：调用 C++ 标准库工具函数/容器接口。
- **[行 1431]** `          }`：作用域边界（代码块开始/结束）。
- **[行 1432]** `        else // if t == 0, i is in the active set`：条件分支的兜底路径。
- **[行 1433]** `          {`：作用域边界（代码块开始/结束）。
- **[行 1434]** `            gradient_d.block(m_d_dof)[i] = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1435]** `            if (gradient_g.block(m_d_dof)[i] > 0)`：条件分支：根据当前状态选择执行路径。
- **[行 1436]** `              m_active_set_phasefield(i) = 1; //lower bound`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1437]** `            else`：条件分支的兜底路径。
- **[行 1438]** `              m_active_set_phasefield(i) = 2; //upper bound`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1439]** `          }`：作用域边界（代码块开始/结束）。
- **[行 1440]** `      }`：作用域边界（代码块开始/结束）。
- **[行 1441]** ``：空行，用于分隔逻辑块。
- **[行 1442]** `    return break_points_sorted;`：返回当前函数结果。
- **[行 1443]** `  }`：作用域边界（代码块开始/结束）。
- **[行 1444]** ``：空行，用于分隔逻辑块。
- **[行 1445]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 1446]** `  void PhaseFieldMonolithicSolve<dim>::get_error_residual(Errors &error_residual)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1447]** `  {`：作用域边界（代码块开始/结束）。
- **[行 1448]** `    BlockVector<double> error_res(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 1449]** ``：空行，用于分隔逻辑块。
- **[行 1450]** `    for (unsigned int i = 0; i < m_dof_handler.n_dofs(); ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 1451]** `      if (!m_constraints.is_constrained(i))`：条件分支：根据当前状态选择执行路径。
- **[行 1452]** `        error_res(i) = m_system_rhs(i);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1453]** ``：空行，用于分隔逻辑块。
- **[行 1454]** `    error_residual.m_norm = error_res.l2_norm();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1455]** `    error_residual.m_u    = error_res.block(m_u_dof).l2_norm();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1456]** `    error_residual.m_d    = error_res.block(m_d_dof).l2_norm();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1457]** `  }`：作用域边界（代码块开始/结束）。
- **[行 1458]** ``：空行，用于分隔逻辑块。
- **[行 1459]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 1460]** `  void PhaseFieldMonolithicSolve<dim>::get_error_residual_LBFGSB(Errors &error_residual,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1461]** `								 const BlockVector<double> & solution_delta)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1462]** `  {`：作用域边界（代码块开始/结束）。
- **[行 1463]** `    // We use L_2 norm`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1464]** `    BlockVector<double> error_res(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 1465]** ``：空行，用于分隔逻辑块。
- **[行 1466]** `    // For displacement DOFs, except essential boundary conditions`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1467]** `    // and hanging-node constraints, there are no box constraints`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1468]** `    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 1469]** `      {`：作用域边界（代码块开始/结束）。
- **[行 1470]** `        if (!m_constraints.is_constrained(i))`：条件分支：根据当前状态选择执行路径。
- **[行 1471]** `	  error_res.block(m_u_dof)[i] = m_system_rhs.block(m_u_dof)[i];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1472]** `      }`：作用域边界（代码块开始/结束）。
- **[行 1473]** ``：空行，用于分隔逻辑块。
- **[行 1474]** `    // For phasefield DOFs, there are points with active box constraints`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1475]** `    // and points with inactive box constraints`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1476]** `    const double upper_limit = 1.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1477]** `    Vector<double> solution_phasefield_total(m_solution.block(m_d_dof));`：函数调用语句，触发对应计算或操作。
- **[行 1478]** `    solution_phasefield_total += solution_delta.block(m_d_dof);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1479]** ``：空行，用于分隔逻辑块。
- **[行 1480]** `    double trial_solution = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1481]** `    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 1482]** `      {`：作用域边界（代码块开始/结束）。
- **[行 1483]** `	// phasefield DOFs can still be constrained due to hanging-nodes`：注释行，用于说明算法背景、假设或实现意图。
- **[行 1484]** `        if (!m_constraints.is_constrained(i + m_dofs_per_block[m_u_dof]))`：条件分支：根据当前状态选择执行路径。
- **[行 1485]** `          {`：作用域边界（代码块开始/结束）。
- **[行 1486]** `            trial_solution = solution_phasefield_total(i) - m_system_rhs.block(m_d_dof)[i];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1487]** ``：空行，用于分隔逻辑块。
- **[行 1488]** `            if (trial_solution < m_solution.block(m_d_dof)[i])`：条件分支：根据当前状态选择执行路径。
- **[行 1489]** `              error_res.block(m_d_dof)[i] = m_solution.block(m_d_dof)[i] - solution_phasefield_total(i);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1490]** `            else if (trial_solution > upper_limit)`：多分支条件判断，处理备选情形。
- **[行 1491]** `              error_res.block(m_d_dof)[i] = upper_limit - solution_phasefield_total(i);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1492]** `            else`：条件分支的兜底路径。
- **[行 1493]** `              error_res.block(m_d_dof)[i] = (-m_system_rhs.block(m_d_dof)[i]);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1494]** `          }`：作用域边界（代码块开始/结束）。
- **[行 1495]** `      }`：作用域边界（代码块开始/结束）。
- **[行 1496]** ``：空行，用于分隔逻辑块。
- **[行 1497]** `    error_residual.m_norm = error_res.l2_norm();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1498]** `    error_residual.m_u    = error_res.block(m_u_dof).l2_norm();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1499]** `    error_residual.m_d    = error_res.block(m_d_dof).l2_norm();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1500]** `  }`：作用域边界（代码块开始/结束）。
- **[行 1501]** ``：空行，用于分隔逻辑块。
- **[行 1502]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 1503]** `  void PhaseFieldMonolithicSolve<dim>::get_error_update(const BlockVector<double> &newton_update,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1504]** `							Errors & error_update)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1505]** `  {`：作用域边界（代码块开始/结束）。
- **[行 1506]** `    BlockVector<double> error_ud(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 1507]** `    for (unsigned int i = 0; i < m_dof_handler.n_dofs(); ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 1508]** `      if (!m_constraints.is_constrained(i))`：条件分支：根据当前状态选择执行路径。
- **[行 1509]** `	error_ud(i) = newton_update(i);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1510]** ``：空行，用于分隔逻辑块。
- **[行 1511]** `    error_update.m_norm = error_ud.l2_norm();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1512]** `    error_update.m_u    = error_ud.block(m_u_dof).l2_norm();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1513]** `    error_update.m_d    = error_ud.block(m_d_dof).l2_norm();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1514]** `  }`：作用域边界（代码块开始/结束）。
- **[行 1515]** ``：空行，用于分隔逻辑块。
- **[行 1516]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 1517]** `  void PhaseFieldMonolithicSolve<dim>::read_material_data(const std::string &data_file,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1518]** `				                     const unsigned int total_material_regions)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1519]** `  {`：作用域边界（代码块开始/结束）。
- **[行 1520]** `    std::ifstream myfile (data_file);`：调用 C++ 标准库工具函数/容器接口。
- **[行 1521]** ``：空行，用于分隔逻辑块。
- **[行 1522]** `    double lame_lambda, lame_mu, length_scale, gc, viscosity, residual_k;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1523]** `    int material_region;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1524]** `    double poisson_ratio;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1525]** `    if (myfile.is_open())`：条件分支：根据当前状态选择执行路径。
- **[行 1526]** `      {`：作用域边界（代码块开始/结束）。
- **[行 1527]** `        m_logfile << "Reading material data file ..." << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 1528]** ``：空行，用于分隔逻辑块。
- **[行 1529]** `        while ( myfile >> material_region`：循环推进（通常用于时间步或迭代收敛过程）。
- **[行 1530]** `                       >> lame_lambda`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1531]** `		       >> lame_mu`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1532]** `		       >> length_scale`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1533]** `		       >> gc`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1534]** `		       >> viscosity`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1535]** `		       >> residual_k)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1536]** `          {`：作用域边界（代码块开始/结束）。
- **[行 1537]** `            m_material_data[material_region] = {lame_lambda,`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1538]** `        	                                lame_mu,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1539]** `						length_scale,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1540]** `						gc,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1541]** `						viscosity,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1542]** `                                                residual_k};`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1543]** `            poisson_ratio = lame_lambda / (2*(lame_lambda + lame_mu));`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1544]** `            Assert( (poisson_ratio <= 0.5)&(poisson_ratio >=-1.0) , ExcInternalError());`：运行期断言/检查，验证输入与状态合法性。
- **[行 1545]** ``：空行，用于分隔逻辑块。
- **[行 1546]** `            m_logfile << "\tRegion " << material_region << " : " << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 1547]** `            m_logfile << "\t\tLame lambda = " << lame_lambda << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 1548]** `            m_logfile << "\t\tLame mu = "  << lame_mu << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 1549]** `            m_logfile << "\t\tPoisson ratio = "  << poisson_ratio << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 1550]** `            m_logfile << "\t\tPhase field length scale (l) = " << length_scale << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 1551]** `            m_logfile << "\t\tCritical energy release rate (gc) = "  << gc << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 1552]** `            m_logfile << "\t\tViscosity for regularization (eta) = "  << viscosity << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 1553]** `            m_logfile << "\t\tResidual_k (k) = "  << residual_k << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 1554]** `          }`：作用域边界（代码块开始/结束）。
- **[行 1555]** ``：空行，用于分隔逻辑块。
- **[行 1556]** `        if (m_material_data.size() != total_material_regions)`：条件分支：根据当前状态选择执行路径。
- **[行 1557]** `          {`：作用域边界（代码块开始/结束）。
- **[行 1558]** `            m_logfile << "Material data file has " << m_material_data.size() << " rows. However, "`：写日志输出，记录当前计算状态与结果。
- **[行 1559]** `        	      << "the mesh has " << total_material_regions << " material regions."`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1560]** `		      << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1561]** `            Assert(m_material_data.size() == total_material_regions,`：运行期断言/检查，验证输入与状态合法性。
- **[行 1562]** `                       ExcDimensionMismatch(m_material_data.size(), total_material_regions));`：函数调用语句，触发对应计算或操作。
- **[行 1563]** `          }`：作用域边界（代码块开始/结束）。
- **[行 1564]** `        myfile.close();`：函数调用语句，触发对应计算或操作。
- **[行 1565]** `      }`：作用域边界（代码块开始/结束）。
- **[行 1566]** `    else`：条件分支的兜底路径。
- **[行 1567]** `      {`：作用域边界（代码块开始/结束）。
- **[行 1568]** `	m_logfile << "Material data file : " << data_file << " not exist!" << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 1569]** `	Assert(false, ExcMessage("Failed to read material data file"));`：运行期断言/检查，验证输入与状态合法性。
- **[行 1570]** `      }`：作用域边界（代码块开始/结束）。
- **[行 1571]** `  }`：作用域边界（代码块开始/结束）。
- **[行 1572]** ``：空行，用于分隔逻辑块。
- **[行 1573]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 1574]** `  void PhaseFieldMonolithicSolve<dim>::read_time_data(const std::string &data_file,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1575]** `				                 std::vector<std::array<double, 4>> & time_table)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1576]** `  {`：作用域边界（代码块开始/结束）。
- **[行 1577]** `    std::ifstream myfile (data_file);`：调用 C++ 标准库工具函数/容器接口。
- **[行 1578]** ``：空行，用于分隔逻辑块。
- **[行 1579]** `    double t_0, t_1, delta_t, t_magnitude;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1580]** ``：空行，用于分隔逻辑块。
- **[行 1581]** `    if (myfile.is_open())`：条件分支：根据当前状态选择执行路径。
- **[行 1582]** `      {`：作用域边界（代码块开始/结束）。
- **[行 1583]** `	m_logfile << "Reading time data file ..." << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 1584]** ``：空行，用于分隔逻辑块。
- **[行 1585]** `	while ( myfile >> t_0`：循环推进（通常用于时间步或迭代收敛过程）。
- **[行 1586]** `		       >> t_1`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1587]** `		       >> delta_t`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1588]** `		       >> t_magnitude)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1589]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 1590]** `	    Assert( t_0 < t_1,`：运行期断言/检查，验证输入与状态合法性。
- **[行 1591]** `		    ExcMessage("For each time pair, "`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1592]** `			       "the start time should be smaller than the end time"));`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1593]** `	    time_table.push_back({{t_0, t_1, delta_t, t_magnitude}});`：函数调用语句，触发对应计算或操作。
- **[行 1594]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 1595]** ``：空行，用于分隔逻辑块。
- **[行 1596]** `	Assert(std::fabs(t_1 - m_parameters.m_end_time) < 1.0e-9,`：运行期断言/检查，验证输入与状态合法性。
- **[行 1597]** `	       ExcMessage("End time in time table is inconsistent with input data in parameters.prm"));`：函数调用语句，触发对应计算或操作。
- **[行 1598]** ``：空行，用于分隔逻辑块。
- **[行 1599]** `	Assert(time_table.size() > 0,`：运行期断言/检查，验证输入与状态合法性。
- **[行 1600]** `	       ExcMessage("Time data file is empty."));`：函数调用语句，触发对应计算或操作。
- **[行 1601]** `	myfile.close();`：函数调用语句，触发对应计算或操作。
- **[行 1602]** `      }`：作用域边界（代码块开始/结束）。
- **[行 1603]** `    else`：条件分支的兜底路径。
- **[行 1604]** `      {`：作用域边界（代码块开始/结束）。
- **[行 1605]** `        m_logfile << "Time data file : " << data_file << " not exist!" << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 1606]** `        Assert(false, ExcMessage("Failed to read time data file"));`：运行期断言/检查，验证输入与状态合法性。
- **[行 1607]** `      }`：作用域边界（代码块开始/结束）。
- **[行 1608]** ``：空行，用于分隔逻辑块。
- **[行 1609]** `    for (auto & time_group : time_table)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 1610]** `      {`：作用域边界（代码块开始/结束）。
- **[行 1611]** `	m_logfile << "\t\t"`：写日志输出，记录当前计算状态与结果。
- **[行 1612]** `	          << time_group[0] << ",\t"`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1613]** `	          << time_group[1] << ",\t"`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1614]** `		  << time_group[2] << ",\t"`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1615]** `		  << time_group[3] << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1616]** `      }`：作用域边界（代码块开始/结束）。
- **[行 1617]** `  }`：作用域边界（代码块开始/结束）。
- **[行 1618]** ``：空行，用于分隔逻辑块。
- **[行 1619]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 1620]** `  void PhaseFieldMonolithicSolve<dim>::setup_qph()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1621]** `  {`：作用域边界（代码块开始/结束）。
- **[行 1622]** `    m_logfile << "\t\tSetting up quadrature point data ("`：写日志输出，记录当前计算状态与结果。
- **[行 1623]** `	      << m_n_q_points`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1624]** `	      << " points per cell)" << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1625]** ``：空行，用于分隔逻辑块。
- **[行 1626]** `    m_quadrature_point_history.clear();`：函数调用语句，触发对应计算或操作。
- **[行 1627]** `    for (auto const & cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 1628]** `      {`：作用域边界（代码块开始/结束）。
- **[行 1629]** `	m_quadrature_point_history.initialize(cell, m_n_q_points);`：函数调用语句，触发对应计算或操作。
- **[行 1630]** `      }`：作用域边界（代码块开始/结束）。
- **[行 1631]** ``：空行，用于分隔逻辑块。
- **[行 1632]** `    unsigned int material_id;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1633]** `    double lame_lambda = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1634]** `    double lame_mu = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1635]** `    double length_scale = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1636]** `    double gc = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1637]** `    double viscosity = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1638]** `    double residual_k = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1639]** ``：空行，用于分隔逻辑块。
- **[行 1640]** `    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 1641]** `      {`：作用域边界（代码块开始/结束）。
- **[行 1642]** `        material_id = cell->material_id();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1643]** `        if (m_material_data.find(material_id) != m_material_data.end())`：条件分支：根据当前状态选择执行路径。
- **[行 1644]** `          {`：作用域边界（代码块开始/结束）。
- **[行 1645]** `            lame_lambda                = m_material_data[material_id][0];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1646]** `            lame_mu                    = m_material_data[material_id][1];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1647]** `            length_scale               = m_material_data[material_id][2];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1648]** `            gc                         = m_material_data[material_id][3];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1649]** `            viscosity                  = m_material_data[material_id][4];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1650]** `            residual_k                 = m_material_data[material_id][5];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1651]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 1652]** `        else`：条件分支的兜底路径。
- **[行 1653]** `          {`：作用域边界（代码块开始/结束）。
- **[行 1654]** `            m_logfile << "Could not find material data for material id: " << material_id << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 1655]** `            AssertThrow(false, ExcMessage("Could not find material data for material id."));`：运行期断言/检查，验证输入与状态合法性。
- **[行 1656]** `          }`：作用域边界（代码块开始/结束）。
- **[行 1657]** ``：空行，用于分隔逻辑块。
- **[行 1658]** `        const std::vector<std::shared_ptr<PointHistory<dim>>> lqph =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1659]** `          m_quadrature_point_history.get_data(cell);`：函数调用语句，触发对应计算或操作。
- **[行 1660]** `        Assert(lqph.size() == m_n_q_points, ExcInternalError());`：运行期断言/检查，验证输入与状态合法性。
- **[行 1661]** ``：空行，用于分隔逻辑块。
- **[行 1662]** `        for (unsigned int q_point = 0; q_point < m_n_q_points; ++q_point)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 1663]** `          lqph[q_point]->setup_lqp(lame_lambda, lame_mu, length_scale,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1664]** `				   gc, viscosity, residual_k,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1665]** `				   m_parameters.m_plane_stress);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1666]** `      }`：作用域边界（代码块开始/结束）。
- **[行 1667]** `  }`：作用域边界（代码块开始/结束）。
- **[行 1668]** ``：空行，用于分隔逻辑块。
- **[行 1669]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 1670]** `  BlockVector<double> PhaseFieldMonolithicSolve<dim>::get_total_solution(`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1671]** `    const BlockVector<double> &solution_delta) const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1672]** `  {`：作用域边界（代码块开始/结束）。
- **[行 1673]** `    BlockVector<double> solution_total(m_solution);`：函数调用语句，触发对应计算或操作。
- **[行 1674]** `    solution_total += solution_delta;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1675]** `    return solution_total;`：返回当前函数结果。
- **[行 1676]** `  }`：作用域边界（代码块开始/结束）。
- **[行 1677]** ``：空行，用于分隔逻辑块。
- **[行 1678]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 1679]** `  void`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1680]** `  PhaseFieldMonolithicSolve<dim>::update_qph_incremental(const BlockVector<double> &solution_delta,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 1681]** `							 const BlockVector<double> &solution_old)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1682]** `  {`：作用域边界（代码块开始/结束）。
- **[行 1683]** `    m_timer.enter_subsection("Update QPH data");`：函数调用语句，触发对应计算或操作。
- **[行 1684]** ``：空行，用于分隔逻辑块。
- **[行 1685]** `    const BlockVector<double> solution_total(get_total_solution(solution_delta));`：函数调用语句，触发对应计算或操作。
- **[行 1686]** ``：空行，用于分隔逻辑块。
- **[行 1687]** `    const UpdateFlags uf_UQPH(update_values | update_gradients);`：函数调用语句，触发对应计算或操作。
- **[行 1688]** `    PerTaskData_UQPH  per_task_data_UQPH;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1689]** `    ScratchData_UQPH  scratch_data_UQPH(m_fe,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1690]** `					m_qf_cell,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1691]** `					uf_UQPH,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1692]** `					solution_total,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1693]** `					solution_old,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1694]** `					m_time.get_delta_t());`：函数调用语句，触发对应计算或操作。
- **[行 1695]** ``：空行，用于分隔逻辑块。
- **[行 1696]** `    auto worker = [this](const typename DoFHandler<dim>::active_cell_iterator &cell,`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1697]** `	                 ScratchData_UQPH & scratch,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1698]** `	                 PerTaskData_UQPH & data)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1699]** `      {`：作用域边界（代码块开始/结束）。
- **[行 1700]** `        this->update_qph_incremental_one_cell(cell, scratch, data);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 1701]** `      };`：作用域边界（代码块开始/结束）。
- **[行 1702]** ``：空行，用于分隔逻辑块。
- **[行 1703]** `    auto copier = [this](const PerTaskData_UQPH &data)`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1704]** `      {`：作用域边界（代码块开始/结束）。
- **[行 1705]** `        this->copy_local_to_global_UQPH(data);`：函数调用语句，触发对应计算或操作。
- **[行 1706]** `      };`：作用域边界（代码块开始/结束）。
- **[行 1707]** ``：空行，用于分隔逻辑块。
- **[行 1708]** `    WorkStream::run(`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1709]** `	m_dof_handler.begin_active(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1710]** `	m_dof_handler.end(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1711]** `	worker,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1712]** `	copier,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1713]** `	scratch_data_UQPH,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1714]** `	per_task_data_UQPH);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1715]** ``：空行，用于分隔逻辑块。
- **[行 1716]** `    m_timer.leave_subsection();`：函数调用语句，触发对应计算或操作。
- **[行 1717]** `  }`：作用域边界（代码块开始/结束）。
- **[行 1718]** ``：空行，用于分隔逻辑块。
- **[行 1719]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 1720]** `  struct PhaseFieldMonolithicSolve<dim>::PerTaskData_UQPH`：类型声明（类/结构体），封装数据与行为。
- **[行 1721]** `  {`：作用域边界（代码块开始/结束）。
- **[行 1722]** `    void reset()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1723]** `    {}`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1724]** `  };`：作用域边界（代码块开始/结束）。
- **[行 1725]** ``：空行，用于分隔逻辑块。
- **[行 1726]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 1727]** `  struct PhaseFieldMonolithicSolve<dim>::ScratchData_UQPH`：类型声明（类/结构体），封装数据与行为。
- **[行 1728]** `  {`：作用域边界（代码块开始/结束）。
- **[行 1729]** `    const BlockVector<double> & m_solution_UQPH;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1730]** ``：空行，用于分隔逻辑块。
- **[行 1731]** `    std::vector<SymmetricTensor<2, dim>> m_solution_symm_grads_u_cell;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1732]** `    std::vector<double>         m_solution_values_phasefield_cell;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1733]** `    std::vector<Tensor<1, dim>> m_solution_grad_phasefield_cell;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1734]** ``：空行，用于分隔逻辑块。
- **[行 1735]** `    FEValues<dim> m_fe_values;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1736]** ``：空行，用于分隔逻辑块。
- **[行 1737]** `    const BlockVector<double>&       m_solution_previous_step;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1738]** `    std::vector<double>              m_phasefield_previous_step_cell;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1739]** ``：空行，用于分隔逻辑块。
- **[行 1740]** `    const double                     m_delta_time;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1741]** ``：空行，用于分隔逻辑块。
- **[行 1742]** `    ScratchData_UQPH(const FiniteElement<dim> & fe_cell,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1743]** `                     const QGauss<dim> &        qf_cell,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1744]** `                     const UpdateFlags          uf_cell,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1745]** `                     const BlockVector<double> &solution_total,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1746]** `		     const BlockVector<double> &solution_old,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1747]** `		     const double delta_time)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1748]** `      : m_solution_UQPH(solution_total)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1749]** `      , m_solution_symm_grads_u_cell(qf_cell.size())`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1750]** `      , m_solution_values_phasefield_cell(qf_cell.size())`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1751]** `      , m_solution_grad_phasefield_cell(qf_cell.size())`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1752]** `      , m_fe_values(fe_cell, qf_cell, uf_cell)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1753]** `      , m_solution_previous_step(solution_old)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1754]** `      , m_phasefield_previous_step_cell(qf_cell.size())`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1755]** `      , m_delta_time(delta_time)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1756]** `    {}`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1757]** ``：空行，用于分隔逻辑块。
- **[行 1758]** `    ScratchData_UQPH(const ScratchData_UQPH &rhs)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1759]** `      : m_solution_UQPH(rhs.m_solution_UQPH)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1760]** `      , m_solution_symm_grads_u_cell(rhs.m_solution_symm_grads_u_cell)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1761]** `      , m_solution_values_phasefield_cell(rhs.m_solution_values_phasefield_cell)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1762]** `      , m_solution_grad_phasefield_cell(rhs.m_solution_grad_phasefield_cell)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1763]** `      , m_fe_values(rhs.m_fe_values.get_fe(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1764]** `                    rhs.m_fe_values.get_quadrature(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1765]** `                    rhs.m_fe_values.get_update_flags())`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1766]** `      , m_solution_previous_step(rhs.m_solution_previous_step)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1767]** `      , m_phasefield_previous_step_cell(rhs.m_phasefield_previous_step_cell)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1768]** `      , m_delta_time(rhs.m_delta_time)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1769]** `    {}`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1770]** ``：空行，用于分隔逻辑块。
- **[行 1771]** `    void reset()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1772]** `    {`：作用域边界（代码块开始/结束）。
- **[行 1773]** `      const unsigned int n_q_points = m_solution_symm_grads_u_cell.size();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1774]** `      for (unsigned int q = 0; q < n_q_points; ++q)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 1775]** `        {`：作用域边界（代码块开始/结束）。
- **[行 1776]** `          m_solution_symm_grads_u_cell[q]  = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1777]** `          m_solution_values_phasefield_cell[q] = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1778]** `          m_solution_grad_phasefield_cell[q] = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1779]** `          m_phasefield_previous_step_cell[q] = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1780]** `        }`：作用域边界（代码块开始/结束）。
- **[行 1781]** `    }`：作用域边界（代码块开始/结束）。
- **[行 1782]** `  };`：作用域边界（代码块开始/结束）。
- **[行 1783]** ``：空行，用于分隔逻辑块。
- **[行 1784]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 1785]** `  void PhaseFieldMonolithicSolve<dim>::update_qph_incremental_one_cell(`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 1786]** `    const typename DoFHandler<dim>::active_cell_iterator &cell,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1787]** `    ScratchData_UQPH & scratch,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1788]** `    PerTaskData_UQPH & /*data*/)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1789]** `  {`：作用域边界（代码块开始/结束）。
- **[行 1790]** `    scratch.reset();`：函数调用语句，触发对应计算或操作。
- **[行 1791]** ``：空行，用于分隔逻辑块。
- **[行 1792]** `    scratch.m_fe_values.reinit(cell);`：函数调用语句，触发对应计算或操作。
- **[行 1793]** ``：空行，用于分隔逻辑块。
- **[行 1794]** `    const std::vector<std::shared_ptr<PointHistory<dim>>> lqph =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1795]** `      m_quadrature_point_history.get_data(cell);`：函数调用语句，触发对应计算或操作。
- **[行 1796]** `    Assert(lqph.size() == m_n_q_points, ExcInternalError());`：运行期断言/检查，验证输入与状态合法性。
- **[行 1797]** ``：空行，用于分隔逻辑块。
- **[行 1798]** `    const FEValuesExtractors::Vector displacement(0);`：函数调用语句，触发对应计算或操作。
- **[行 1799]** ``：空行，用于分隔逻辑块。
- **[行 1800]** `    scratch.m_fe_values[m_u_fe].get_function_symmetric_gradients(`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1801]** `      scratch.m_solution_UQPH, scratch.m_solution_symm_grads_u_cell);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1802]** `    scratch.m_fe_values[m_d_fe].get_function_values(`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1803]** `      scratch.m_solution_UQPH, scratch.m_solution_values_phasefield_cell);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1804]** `    scratch.m_fe_values[m_d_fe].get_function_gradients(`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1805]** `      scratch.m_solution_UQPH, scratch.m_solution_grad_phasefield_cell);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1806]** ``：空行，用于分隔逻辑块。
- **[行 1807]** `    scratch.m_fe_values[m_d_fe].get_function_values(`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1808]** `      scratch.m_solution_previous_step, scratch.m_phasefield_previous_step_cell);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1809]** ``：空行，用于分隔逻辑块。
- **[行 1810]** `    for (const unsigned int q_point :`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 1811]** `         scratch.m_fe_values.quadrature_point_indices())`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1812]** `      lqph[q_point]->update_field_values(scratch.m_solution_symm_grads_u_cell[q_point],`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1813]** `                                         scratch.m_solution_values_phasefield_cell[q_point],`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1814]** `					 scratch.m_solution_grad_phasefield_cell[q_point],`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1815]** `					 scratch.m_phasefield_previous_step_cell[q_point],`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1816]** `					 scratch.m_delta_time);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1817]** `  }`：作用域边界（代码块开始/结束）。
- **[行 1818]** ``：空行，用于分隔逻辑块。
- **[行 1819]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 1820]** `  struct PhaseFieldMonolithicSolve<dim>::PerTaskData_ASM`：类型声明（类/结构体），封装数据与行为。
- **[行 1821]** `  {`：作用域边界（代码块开始/结束）。
- **[行 1822]** `    FullMatrix<double>                   m_cell_matrix;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1823]** `    Vector<double>                       m_cell_rhs;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1824]** `    std::vector<types::global_dof_index> m_local_dof_indices;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1825]** ``：空行，用于分隔逻辑块。
- **[行 1826]** `    PerTaskData_ASM(const unsigned int dofs_per_cell)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1827]** `      : m_cell_matrix(dofs_per_cell, dofs_per_cell)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1828]** `      , m_cell_rhs(dofs_per_cell)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1829]** `      , m_local_dof_indices(dofs_per_cell)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1830]** `    {}`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1831]** ``：空行，用于分隔逻辑块。
- **[行 1832]** `    void reset()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1833]** `    {`：作用域边界（代码块开始/结束）。
- **[行 1834]** `      m_cell_matrix = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1835]** `      m_cell_rhs    = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1836]** `    }`：作用域边界（代码块开始/结束）。
- **[行 1837]** `  };`：作用域边界（代码块开始/结束）。
- **[行 1838]** ``：空行，用于分隔逻辑块。
- **[行 1839]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 1840]** `  struct PhaseFieldMonolithicSolve<dim>::PerTaskData_ASM_RHS_BFGS`：类型声明（类/结构体），封装数据与行为。
- **[行 1841]** `  {`：作用域边界（代码块开始/结束）。
- **[行 1842]** `    Vector<double>                       m_cell_rhs;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1843]** `    std::vector<types::global_dof_index> m_local_dof_indices;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1844]** ``：空行，用于分隔逻辑块。
- **[行 1845]** `    PerTaskData_ASM_RHS_BFGS(const unsigned int dofs_per_cell)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1846]** `      : m_cell_rhs(dofs_per_cell)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1847]** `      , m_local_dof_indices(dofs_per_cell)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1848]** `    {}`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1849]** ``：空行，用于分隔逻辑块。
- **[行 1850]** `    void reset()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1851]** `    {`：作用域边界（代码块开始/结束）。
- **[行 1852]** `      m_cell_rhs    = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1853]** `    }`：作用域边界（代码块开始/结束）。
- **[行 1854]** `  };`：作用域边界（代码块开始/结束）。
- **[行 1855]** ``：空行，用于分隔逻辑块。
- **[行 1856]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 1857]** `  struct PhaseFieldMonolithicSolve<dim>::ScratchData_ASM`：类型声明（类/结构体），封装数据与行为。
- **[行 1858]** `  {`：作用域边界（代码块开始/结束）。
- **[行 1859]** `    FEValues<dim>     m_fe_values;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1860]** `    FEFaceValues<dim> m_fe_face_values;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1861]** ``：空行，用于分隔逻辑块。
- **[行 1862]** `    std::vector<std::vector<double>>                  m_Nx_phasefield;      // shape function values for phase-field`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1863]** `    std::vector<std::vector<Tensor<1, dim>>>          m_grad_Nx_phasefield; // gradient of shape function values for phase field`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1864]** ``：空行，用于分隔逻辑块。
- **[行 1865]** `    std::vector<std::vector<Tensor<1, dim>>>          m_Nx_disp;       // shape function values for displacement`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1866]** `    std::vector<std::vector<Tensor<2, dim>>>          m_grad_Nx_disp;  // gradient of shape function values for displacement`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1867]** `    std::vector<std::vector<SymmetricTensor<2, dim>>> m_symm_grad_Nx_disp;  // symmetric gradient of shape function values for displacement`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1868]** ``：空行，用于分隔逻辑块。
- **[行 1869]** `    const BlockVector<double>&       m_solution_previous_step;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1870]** `    std::vector<double>              m_phasefield_previous_step_cell;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1871]** ``：空行，用于分隔逻辑块。
- **[行 1872]** `    ScratchData_ASM(const FiniteElement<dim> & fe_cell,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1873]** `                    const QGauss<dim> &        qf_cell,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1874]** `                    const UpdateFlags          uf_cell,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1875]** `		    const QGauss<dim - 1> &    qf_face,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1876]** `		    const UpdateFlags          uf_face,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1877]** `		    const BlockVector<double>& solution_old)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1878]** `      : m_fe_values(fe_cell, qf_cell, uf_cell)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1879]** `      , m_fe_face_values(fe_cell, qf_face, uf_face)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1880]** `      , m_Nx_phasefield(qf_cell.size(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1881]** `	                std::vector<double>(fe_cell.n_dofs_per_cell()))`：调用 C++ 标准库工具函数/容器接口。
- **[行 1882]** `      , m_grad_Nx_phasefield(qf_cell.size(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1883]** `		             std::vector<Tensor<1, dim>>(fe_cell.n_dofs_per_cell()))`：调用 C++ 标准库工具函数/容器接口。
- **[行 1884]** `      , m_Nx_disp(qf_cell.size(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1885]** `		  std::vector<Tensor<1, dim>>(fe_cell.n_dofs_per_cell()))`：调用 C++ 标准库工具函数/容器接口。
- **[行 1886]** `      , m_grad_Nx_disp(qf_cell.size(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1887]** `                       std::vector<Tensor<2, dim>>(fe_cell.n_dofs_per_cell()))`：调用 C++ 标准库工具函数/容器接口。
- **[行 1888]** `      , m_symm_grad_Nx_disp(qf_cell.size(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1889]** `                            std::vector<SymmetricTensor<2, dim>>(fe_cell.n_dofs_per_cell()))`：调用 C++ 标准库工具函数/容器接口。
- **[行 1890]** `      , m_solution_previous_step(solution_old)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1891]** `      , m_phasefield_previous_step_cell(qf_cell.size())`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1892]** `    {}`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1893]** ``：空行，用于分隔逻辑块。
- **[行 1894]** `    ScratchData_ASM(const ScratchData_ASM &rhs)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1895]** `      : m_fe_values(rhs.m_fe_values.get_fe(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1896]** `                    rhs.m_fe_values.get_quadrature(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1897]** `                    rhs.m_fe_values.get_update_flags())`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1898]** `      , m_fe_face_values(rhs.m_fe_face_values.get_fe(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1899]** `	                 rhs.m_fe_face_values.get_quadrature(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1900]** `	                 rhs.m_fe_face_values.get_update_flags())`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1901]** `      , m_Nx_phasefield(rhs.m_Nx_phasefield)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1902]** `      , m_grad_Nx_phasefield(rhs.m_grad_Nx_phasefield)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1903]** `      , m_Nx_disp(rhs.m_Nx_disp)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1904]** `      , m_grad_Nx_disp(rhs.m_grad_Nx_disp)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1905]** `      , m_symm_grad_Nx_disp(rhs.m_symm_grad_Nx_disp)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1906]** `      , m_solution_previous_step(rhs.m_solution_previous_step)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1907]** `      , m_phasefield_previous_step_cell(rhs.m_phasefield_previous_step_cell)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1908]** `    {}`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1909]** ``：空行，用于分隔逻辑块。
- **[行 1910]** `    void reset()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1911]** `    {`：作用域边界（代码块开始/结束）。
- **[行 1912]** `      const unsigned int n_q_points      = m_Nx_phasefield.size();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1913]** `      const unsigned int n_dofs_per_cell = m_Nx_phasefield[0].size();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1914]** `      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 1915]** `        {`：作用域边界（代码块开始/结束）。
- **[行 1916]** `          Assert(m_Nx_phasefield[q_point].size() == n_dofs_per_cell,`：运行期断言/检查，验证输入与状态合法性。
- **[行 1917]** `		 ExcInternalError());`：函数调用语句，触发对应计算或操作。
- **[行 1918]** ``：空行，用于分隔逻辑块。
- **[行 1919]** `          Assert(m_grad_Nx_phasefield[q_point].size() == n_dofs_per_cell,`：运行期断言/检查，验证输入与状态合法性。
- **[行 1920]** `                 ExcInternalError());`：函数调用语句，触发对应计算或操作。
- **[行 1921]** ``：空行，用于分隔逻辑块。
- **[行 1922]** `          Assert(m_Nx_disp[q_point].size() == n_dofs_per_cell,`：运行期断言/检查，验证输入与状态合法性。
- **[行 1923]** `		 ExcInternalError());`：函数调用语句，触发对应计算或操作。
- **[行 1924]** ``：空行，用于分隔逻辑块。
- **[行 1925]** `          Assert(m_grad_Nx_disp[q_point].size() == n_dofs_per_cell,`：运行期断言/检查，验证输入与状态合法性。
- **[行 1926]** `                 ExcInternalError());`：函数调用语句，触发对应计算或操作。
- **[行 1927]** ``：空行，用于分隔逻辑块。
- **[行 1928]** `          Assert(m_symm_grad_Nx_disp[q_point].size() == n_dofs_per_cell,`：运行期断言/检查，验证输入与状态合法性。
- **[行 1929]** `                 ExcInternalError());`：函数调用语句，触发对应计算或操作。
- **[行 1930]** ``：空行，用于分隔逻辑块。
- **[行 1931]** `          m_phasefield_previous_step_cell[q_point] = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1932]** `          for (unsigned int k = 0; k < n_dofs_per_cell; ++k)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 1933]** `            {`：作用域边界（代码块开始/结束）。
- **[行 1934]** `              m_Nx_phasefield[q_point][k]           = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1935]** `              m_grad_Nx_phasefield[q_point][k]      = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1936]** `              m_Nx_disp[q_point][k]                 = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1937]** `              m_grad_Nx_disp[q_point][k]            = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1938]** `              m_symm_grad_Nx_disp[q_point][k]       = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 1939]** `            }`：作用域边界（代码块开始/结束）。
- **[行 1940]** `        }`：作用域边界（代码块开始/结束）。
- **[行 1941]** `    }`：作用域边界（代码块开始/结束）。
- **[行 1942]** `  };`：作用域边界（代码块开始/结束）。
- **[行 1943]** ``：空行，用于分隔逻辑块。
- **[行 1944]** ``：空行，用于分隔逻辑块。
- **[行 1945]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 1946]** `  struct PhaseFieldMonolithicSolve<dim>::ScratchData_ASM_RHS_BFGS`：类型声明（类/结构体），封装数据与行为。
- **[行 1947]** `  {`：作用域边界（代码块开始/结束）。
- **[行 1948]** `    FEValues<dim>     m_fe_values;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1949]** `    FEFaceValues<dim> m_fe_face_values;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1950]** ``：空行，用于分隔逻辑块。
- **[行 1951]** `    std::vector<std::vector<double>>                  m_Nx_phasefield;      // shape function values for phase-field`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1952]** `    std::vector<std::vector<Tensor<1, dim>>>          m_grad_Nx_phasefield; // gradient of shape function values for phase field`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1953]** ``：空行，用于分隔逻辑块。
- **[行 1954]** `    std::vector<std::vector<Tensor<1, dim>>>          m_Nx_disp;       // shape function values for displacement`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1955]** `    std::vector<std::vector<Tensor<2, dim>>>          m_grad_Nx_disp;  // gradient of shape function values for displacement`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1956]** `    std::vector<std::vector<SymmetricTensor<2, dim>>> m_symm_grad_Nx_disp;  // symmetric gradient of shape function values for displacement`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1957]** ``：空行，用于分隔逻辑块。
- **[行 1958]** `    const BlockVector<double>&       m_solution_previous_step;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1959]** `    std::vector<double>              m_phasefield_previous_step_cell;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1960]** ``：空行，用于分隔逻辑块。
- **[行 1961]** `    ScratchData_ASM_RHS_BFGS(const FiniteElement<dim> & fe_cell,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1962]** `                             const QGauss<dim> &        qf_cell,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1963]** `                             const UpdateFlags          uf_cell,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1964]** `		             const QGauss<dim - 1> &    qf_face,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1965]** `		             const UpdateFlags          uf_face,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1966]** `		             const BlockVector<double>& solution_old)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1967]** `      : m_fe_values(fe_cell, qf_cell, uf_cell)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1968]** `      , m_fe_face_values(fe_cell, qf_face, uf_face)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1969]** `      , m_Nx_phasefield(qf_cell.size(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1970]** `	                std::vector<double>(fe_cell.n_dofs_per_cell()))`：调用 C++ 标准库工具函数/容器接口。
- **[行 1971]** `      , m_grad_Nx_phasefield(qf_cell.size(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1972]** `		             std::vector<Tensor<1, dim>>(fe_cell.n_dofs_per_cell()))`：调用 C++ 标准库工具函数/容器接口。
- **[行 1973]** `      , m_Nx_disp(qf_cell.size(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1974]** `		  std::vector<Tensor<1, dim>>(fe_cell.n_dofs_per_cell()))`：调用 C++ 标准库工具函数/容器接口。
- **[行 1975]** `      , m_grad_Nx_disp(qf_cell.size(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1976]** `                       std::vector<Tensor<2, dim>>(fe_cell.n_dofs_per_cell()))`：调用 C++ 标准库工具函数/容器接口。
- **[行 1977]** `      , m_symm_grad_Nx_disp(qf_cell.size(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1978]** `                            std::vector<SymmetricTensor<2, dim>>(fe_cell.n_dofs_per_cell()))`：调用 C++ 标准库工具函数/容器接口。
- **[行 1979]** `      , m_solution_previous_step(solution_old)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1980]** `      , m_phasefield_previous_step_cell(qf_cell.size())`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1981]** `    {}`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1982]** ``：空行，用于分隔逻辑块。
- **[行 1983]** `    ScratchData_ASM_RHS_BFGS(const ScratchData_ASM_RHS_BFGS &rhs)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1984]** `      : m_fe_values(rhs.m_fe_values.get_fe(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1985]** `                    rhs.m_fe_values.get_quadrature(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1986]** `                    rhs.m_fe_values.get_update_flags())`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1987]** `      , m_fe_face_values(rhs.m_fe_face_values.get_fe(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1988]** `	                 rhs.m_fe_face_values.get_quadrature(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1989]** `	                 rhs.m_fe_face_values.get_update_flags())`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1990]** `      , m_Nx_phasefield(rhs.m_Nx_phasefield)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1991]** `      , m_grad_Nx_phasefield(rhs.m_grad_Nx_phasefield)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1992]** `      , m_Nx_disp(rhs.m_Nx_disp)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1993]** `      , m_grad_Nx_disp(rhs.m_grad_Nx_disp)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1994]** `      , m_symm_grad_Nx_disp(rhs.m_symm_grad_Nx_disp)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1995]** `      , m_solution_previous_step(rhs.m_solution_previous_step)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1996]** `      , m_phasefield_previous_step_cell(rhs.m_phasefield_previous_step_cell)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1997]** `    {}`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 1998]** ``：空行，用于分隔逻辑块。
- **[行 1999]** `    void reset()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2000]** `    {`：作用域边界（代码块开始/结束）。
- **[行 2001]** `      const unsigned int n_q_points      = m_Nx_phasefield.size();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2002]** `      const unsigned int n_dofs_per_cell = m_Nx_phasefield[0].size();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2003]** `      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2004]** `        {`：作用域边界（代码块开始/结束）。
- **[行 2005]** `          Assert(m_Nx_phasefield[q_point].size() == n_dofs_per_cell,`：运行期断言/检查，验证输入与状态合法性。
- **[行 2006]** `		 ExcInternalError());`：函数调用语句，触发对应计算或操作。
- **[行 2007]** ``：空行，用于分隔逻辑块。
- **[行 2008]** `          Assert(m_grad_Nx_phasefield[q_point].size() == n_dofs_per_cell,`：运行期断言/检查，验证输入与状态合法性。
- **[行 2009]** `                 ExcInternalError());`：函数调用语句，触发对应计算或操作。
- **[行 2010]** ``：空行，用于分隔逻辑块。
- **[行 2011]** `          Assert(m_Nx_disp[q_point].size() == n_dofs_per_cell,`：运行期断言/检查，验证输入与状态合法性。
- **[行 2012]** `		 ExcInternalError());`：函数调用语句，触发对应计算或操作。
- **[行 2013]** ``：空行，用于分隔逻辑块。
- **[行 2014]** `          Assert(m_grad_Nx_disp[q_point].size() == n_dofs_per_cell,`：运行期断言/检查，验证输入与状态合法性。
- **[行 2015]** `                 ExcInternalError());`：函数调用语句，触发对应计算或操作。
- **[行 2016]** ``：空行，用于分隔逻辑块。
- **[行 2017]** `          Assert(m_symm_grad_Nx_disp[q_point].size() == n_dofs_per_cell,`：运行期断言/检查，验证输入与状态合法性。
- **[行 2018]** `                 ExcInternalError());`：函数调用语句，触发对应计算或操作。
- **[行 2019]** ``：空行，用于分隔逻辑块。
- **[行 2020]** `          m_phasefield_previous_step_cell[q_point] = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2021]** `          for (unsigned int k = 0; k < n_dofs_per_cell; ++k)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2022]** `            {`：作用域边界（代码块开始/结束）。
- **[行 2023]** `              m_Nx_phasefield[q_point][k]           = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2024]** `              m_grad_Nx_phasefield[q_point][k]      = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2025]** `              m_Nx_disp[q_point][k]                 = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2026]** `              m_grad_Nx_disp[q_point][k]            = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2027]** `              m_symm_grad_Nx_disp[q_point][k]       = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2028]** `            }`：作用域边界（代码块开始/结束）。
- **[行 2029]** `        }`：作用域边界（代码块开始/结束）。
- **[行 2030]** `    }`：作用域边界（代码块开始/结束）。
- **[行 2031]** `  };`：作用域边界（代码块开始/结束）。
- **[行 2032]** ``：空行，用于分隔逻辑块。
- **[行 2033]** `  // constructor has no return type`：注释行，用于说明算法背景、假设或实现意图。
- **[行 2034]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 2035]** `  PhaseFieldMonolithicSolve<dim>::PhaseFieldMonolithicSolve(const std::string &input_file)`：调用 C++ 标准库工具函数/容器接口。
- **[行 2036]** `    : m_parameters(input_file)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2037]** `    , m_triangulation(Triangulation<dim>::maximum_smoothing)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2038]** `    , m_time(m_parameters.m_end_time)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2039]** `    , m_logfile(m_parameters.m_logfile_name)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2040]** `    , m_timer(m_logfile, TimerOutput::summary, TimerOutput::wall_times)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2041]** `    , m_dof_handler(m_triangulation)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2042]** `    , m_fe(FE_Q<dim>(m_parameters.m_poly_degree),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2043]** `	   dim, // displacement`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2044]** `	   FE_Q<dim>(m_parameters.m_poly_degree),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2045]** `	   1)   // phasefield`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2046]** `    , m_dofs_per_cell(m_fe.n_dofs_per_cell())`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2047]** `    , m_u_fe(m_first_u_component)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2048]** `    , m_d_fe(m_d_component)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2049]** `    , m_dofs_per_block(m_n_blocks)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2050]** `    , m_qf_cell(m_parameters.m_quad_order)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2051]** `    , m_qf_face(m_parameters.m_quad_order)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2052]** `    , m_n_q_points(m_qf_cell.size())`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2053]** `    , m_vol_reference(0.0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2054]** `  {}`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2055]** ``：空行，用于分隔逻辑块。
- **[行 2056]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 2057]** `  void PhaseFieldMonolithicSolve<dim>::make_grid()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2058]** `  {`：作用域边界（代码块开始/结束）。
- **[行 2059]** `    if (m_parameters.m_scenario == 1)`：条件分支：根据当前状态选择执行路径。
- **[行 2060]** `      make_grid_case_1();`：函数调用语句，触发对应计算或操作。
- **[行 2061]** `    else if (m_parameters.m_scenario == 2)`：多分支条件判断，处理备选情形。
- **[行 2062]** `      make_grid_case_2();`：函数调用语句，触发对应计算或操作。
- **[行 2063]** `    else if (m_parameters.m_scenario == 3)`：多分支条件判断，处理备选情形。
- **[行 2064]** `      make_grid_case_3();`：函数调用语句，触发对应计算或操作。
- **[行 2065]** `    else if (m_parameters.m_scenario == 4)`：多分支条件判断，处理备选情形。
- **[行 2066]** `      make_grid_case_4();`：函数调用语句，触发对应计算或操作。
- **[行 2067]** `    else if (m_parameters.m_scenario == 5)`：多分支条件判断，处理备选情形。
- **[行 2068]** `      make_grid_case_5();`：函数调用语句，触发对应计算或操作。
- **[行 2069]** `    else if (m_parameters.m_scenario == 6)`：多分支条件判断，处理备选情形。
- **[行 2070]** `      make_grid_case_6();`：函数调用语句，触发对应计算或操作。
- **[行 2071]** `    else if (m_parameters.m_scenario == 7)`：多分支条件判断，处理备选情形。
- **[行 2072]** `      make_grid_case_7();`：函数调用语句，触发对应计算或操作。
- **[行 2073]** `    else if (m_parameters.m_scenario == 8)`：多分支条件判断，处理备选情形。
- **[行 2074]** `      make_grid_case_8();`：函数调用语句，触发对应计算或操作。
- **[行 2075]** `    else if (m_parameters.m_scenario == 9)`：多分支条件判断，处理备选情形。
- **[行 2076]** `      make_grid_case_9();`：函数调用语句，触发对应计算或操作。
- **[行 2077]** `    else if (m_parameters.m_scenario == 10)`：多分支条件判断，处理备选情形。
- **[行 2078]** `      make_grid_case_10();`：函数调用语句，触发对应计算或操作。
- **[行 2079]** `    else if (m_parameters.m_scenario == 11)`：多分支条件判断，处理备选情形。
- **[行 2080]** `      make_grid_case_11();`：函数调用语句，触发对应计算或操作。
- **[行 2081]** `    else`：条件分支的兜底路径。
- **[行 2082]** `      Assert(false, ExcMessage("The scenario has not been implemented!"));`：运行期断言/检查，验证输入与状态合法性。
- **[行 2083]** ``：空行，用于分隔逻辑块。
- **[行 2084]** `    m_logfile << "\t\tTriangulation:"`：写日志输出，记录当前计算状态与结果。
- **[行 2085]** `              << "\n\t\t\tNumber of active cells: "`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2086]** `              << m_triangulation.n_active_cells()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2087]** `              << "\n\t\t\tNumber of used vertices: "`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2088]** `              << m_triangulation.n_used_vertices()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2089]** `	      << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2090]** ``：空行，用于分隔逻辑块。
- **[行 2091]** `    std::ofstream out("original_mesh.vtu");`：调用 C++ 标准库工具函数/容器接口。
- **[行 2092]** `    GridOut       grid_out;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2093]** `    grid_out.write_vtu(m_triangulation, out);`：函数调用语句，触发对应计算或操作。
- **[行 2094]** ``：空行，用于分隔逻辑块。
- **[行 2095]** `    m_vol_reference = GridTools::volume(m_triangulation);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2096]** `    m_logfile << "\t\tGrid:\n\t\t\tReference volume: " << m_vol_reference << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2097]** `  }`：作用域边界（代码块开始/结束）。
- **[行 2098]** ``：空行，用于分隔逻辑块。
- **[行 2099]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 2100]** `  void PhaseFieldMonolithicSolve<dim>::make_grid_case_1()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2101]** `  {`：作用域边界（代码块开始/结束）。
- **[行 2102]** `    for (unsigned int i = 0; i < 80; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2103]** `      m_logfile << "*";`：写日志输出，记录当前计算状态与结果。
- **[行 2104]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2105]** `    m_logfile << "\t\t\tSquare tension (unstructured)" << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2106]** `    for (unsigned int i = 0; i < 80; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2107]** `      m_logfile << "*";`：写日志输出，记录当前计算状态与结果。
- **[行 2108]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2109]** ``：空行，用于分隔逻辑块。
- **[行 2110]** `    AssertThrow(dim==2, ExcMessage("The dimension has to be 2D!"));`：运行期断言/检查，验证输入与状态合法性。
- **[行 2111]** ``：空行，用于分隔逻辑块。
- **[行 2112]** `    GridIn<dim> gridin;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2113]** `    gridin.attach_triangulation(m_triangulation);`：函数调用语句，触发对应计算或操作。
- **[行 2114]** `    std::ifstream f("square_tension_unstructured.msh");`：调用 C++ 标准库工具函数/容器接口。
- **[行 2115]** `    gridin.read_msh(f);`：函数调用语句，触发对应计算或操作。
- **[行 2116]** ``：空行，用于分隔逻辑块。
- **[行 2117]** `    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2118]** `      for (const auto &face : cell->face_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2119]** `	{`：作用域边界（代码块开始/结束）。
- **[行 2120]** `	  if (face->at_boundary() == true)`：条件分支：根据当前状态选择执行路径。
- **[行 2121]** `	    {`：作用域边界（代码块开始/结束）。
- **[行 2122]** `	      if (std::fabs(face->center()[1] + 0.5 ) < 1.0e-9 )`：条件分支：根据当前状态选择执行路径。
- **[行 2123]** `		face->set_boundary_id(0);`：函数调用语句，触发对应计算或操作。
- **[行 2124]** `	      else if (std::fabs(face->center()[1] - 0.5 ) < 1.0e-9)`：多分支条件判断，处理备选情形。
- **[行 2125]** `	        face->set_boundary_id(1);`：函数调用语句，触发对应计算或操作。
- **[行 2126]** `	      else`：条件分支的兜底路径。
- **[行 2127]** `	        face->set_boundary_id(2);`：函数调用语句，触发对应计算或操作。
- **[行 2128]** `	    }`：作用域边界（代码块开始/结束）。
- **[行 2129]** `	}`：作用域边界（代码块开始/结束）。
- **[行 2130]** ``：空行，用于分隔逻辑块。
- **[行 2131]** `    m_triangulation.refine_global(m_parameters.m_global_refine_times);`：函数调用语句，触发对应计算或操作。
- **[行 2132]** ``：空行，用于分隔逻辑块。
- **[行 2133]** `    if (m_parameters.m_refinement_strategy == "pre-refine")`：条件分支：根据当前状态选择执行路径。
- **[行 2134]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2135]** `	unsigned int material_id;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2136]** `	double length_scale;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2137]** `	for (unsigned int i = 0; i < m_parameters.m_local_prerefine_times; i++)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2138]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 2139]** `	    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2140]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 2141]** `		if (   std::fabs(cell->center()[1]) < 0.01`：条件分支：根据当前状态选择执行路径。
- **[行 2142]** `		    && cell->center()[0] > 0.495)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2143]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 2144]** `		    material_id = cell->material_id();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2145]** `		    length_scale = m_material_data[material_id][2];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2146]** `		    if (  std::sqrt(cell->measure())`：条件分支：根据当前状态选择执行路径。
- **[行 2147]** `			> length_scale * m_parameters.m_allowed_max_h_l_ratio )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2148]** `		      cell->set_refine_flag();`：函数调用语句，触发对应计算或操作。
- **[行 2149]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 2150]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 2151]** `	    m_triangulation.execute_coarsening_and_refinement();`：函数调用语句，触发对应计算或操作。
- **[行 2152]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 2153]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2154]** `    else if (m_parameters.m_refinement_strategy == "adaptive-refine")`：多分支条件判断，处理备选情形。
- **[行 2155]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2156]** `	unsigned int material_id;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2157]** `	double length_scale;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2158]** `	bool initiation_point_refine_unfinished = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2159]** `	while (initiation_point_refine_unfinished)`：循环推进（通常用于时间步或迭代收敛过程）。
- **[行 2160]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 2161]** `	    initiation_point_refine_unfinished = false;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2162]** `	    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2163]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 2164]** `		if (   std::fabs(cell->center()[1] - 0.0) < 0.05`：条件分支：根据当前状态选择执行路径。
- **[行 2165]** `		    && std::fabs(cell->center()[0] - 0.5) < 0.05)`：调用 C++ 标准库工具函数/容器接口。
- **[行 2166]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 2167]** `		    material_id = cell->material_id();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2168]** `		    length_scale = m_material_data[material_id][2];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2169]** `		    if (  std::sqrt(cell->measure())`：条件分支：根据当前状态选择执行路径。
- **[行 2170]** `			> length_scale * m_parameters.m_allowed_max_h_l_ratio )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2171]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 2172]** `		        cell->set_refine_flag();`：函数调用语句，触发对应计算或操作。
- **[行 2173]** `		        initiation_point_refine_unfinished = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2174]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 2175]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 2176]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 2177]** `	    m_triangulation.execute_coarsening_and_refinement();`：函数调用语句，触发对应计算或操作。
- **[行 2178]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 2179]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2180]** `    else`：条件分支的兜底路径。
- **[行 2181]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2182]** `	AssertThrow(false,`：运行期断言/检查，验证输入与状态合法性。
- **[行 2183]** `	            ExcMessage("Selected mesh refinement strategy not implemented!"));`：函数调用语句，触发对应计算或操作。
- **[行 2184]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2185]** `  }`：作用域边界（代码块开始/结束）。
- **[行 2186]** ``：空行，用于分隔逻辑块。
- **[行 2187]** ``：空行，用于分隔逻辑块。
- **[行 2188]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 2189]** `  void PhaseFieldMonolithicSolve<dim>::make_grid_case_2()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2190]** `  {`：作用域边界（代码块开始/结束）。
- **[行 2191]** `    for (unsigned int i = 0; i < 80; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2192]** `      m_logfile << "*";`：写日志输出，记录当前计算状态与结果。
- **[行 2193]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2194]** `    m_logfile << "\t\t\t\tSquare shear (unstructured)" << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2195]** `    for (unsigned int i = 0; i < 80; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2196]** `      m_logfile << "*";`：写日志输出，记录当前计算状态与结果。
- **[行 2197]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2198]** ``：空行，用于分隔逻辑块。
- **[行 2199]** `    AssertThrow(dim==2, ExcMessage("The dimension has to be 2D!"));`：运行期断言/检查，验证输入与状态合法性。
- **[行 2200]** ``：空行，用于分隔逻辑块。
- **[行 2201]** `    GridIn<dim> gridin;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2202]** `    gridin.attach_triangulation(m_triangulation);`：函数调用语句，触发对应计算或操作。
- **[行 2203]** `    std::ifstream f("square_shear_unstructured.msh");`：调用 C++ 标准库工具函数/容器接口。
- **[行 2204]** `    gridin.read_msh(f);`：函数调用语句，触发对应计算或操作。
- **[行 2205]** ``：空行，用于分隔逻辑块。
- **[行 2206]** `    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2207]** `      for (const auto &face : cell->face_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2208]** `	{`：作用域边界（代码块开始/结束）。
- **[行 2209]** `	  if (face->at_boundary() == true)`：条件分支：根据当前状态选择执行路径。
- **[行 2210]** `	    {`：作用域边界（代码块开始/结束）。
- **[行 2211]** `	      if (std::fabs(face->center()[1] + 0.5 ) < 1.0e-9 )`：条件分支：根据当前状态选择执行路径。
- **[行 2212]** `		face->set_boundary_id(0);`：函数调用语句，触发对应计算或操作。
- **[行 2213]** `	      else if (std::fabs(face->center()[1] - 0.5 ) < 1.0e-9)`：多分支条件判断，处理备选情形。
- **[行 2214]** `	        face->set_boundary_id(1);`：函数调用语句，触发对应计算或操作。
- **[行 2215]** `	      else if (   (std::fabs(face->center()[0] - 0.0 ) < 1.0e-9)`：多分支条件判断，处理备选情形。
- **[行 2216]** `		       || (std::fabs(face->center()[0] - 1.0 ) < 1.0e-9))`：调用 C++ 标准库工具函数/容器接口。
- **[行 2217]** `	        face->set_boundary_id(2);`：函数调用语句，触发对应计算或操作。
- **[行 2218]** `	      else`：条件分支的兜底路径。
- **[行 2219]** `	        face->set_boundary_id(3);`：函数调用语句，触发对应计算或操作。
- **[行 2220]** `	    }`：作用域边界（代码块开始/结束）。
- **[行 2221]** `	}`：作用域边界（代码块开始/结束）。
- **[行 2222]** ``：空行，用于分隔逻辑块。
- **[行 2223]** `    m_triangulation.refine_global(m_parameters.m_global_refine_times);`：函数调用语句，触发对应计算或操作。
- **[行 2224]** ``：空行，用于分隔逻辑块。
- **[行 2225]** `    if (m_parameters.m_refinement_strategy == "pre-refine")`：条件分支：根据当前状态选择执行路径。
- **[行 2226]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2227]** `	unsigned int material_id;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2228]** `	double length_scale;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2229]** `	for (unsigned int i = 0; i < m_parameters.m_local_prerefine_times; i++)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2230]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 2231]** `	    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2232]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 2233]** `		if (    (cell->center()[0] > 0.45)`：条件分支：根据当前状态选择执行路径。
- **[行 2234]** `		     && (cell->center()[1] < 0.05) )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2235]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 2236]** `		    material_id = cell->material_id();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2237]** `		    length_scale = m_material_data[material_id][2];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2238]** `		    if (  std::sqrt(cell->measure())`：条件分支：根据当前状态选择执行路径。
- **[行 2239]** `			> length_scale * m_parameters.m_allowed_max_h_l_ratio )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2240]** `		      cell->set_refine_flag();`：函数调用语句，触发对应计算或操作。
- **[行 2241]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 2242]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 2243]** `	    m_triangulation.execute_coarsening_and_refinement();`：函数调用语句，触发对应计算或操作。
- **[行 2244]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 2245]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2246]** `    else if (m_parameters.m_refinement_strategy == "adaptive-refine")`：多分支条件判断，处理备选情形。
- **[行 2247]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2248]** `	unsigned int material_id;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2249]** `	double length_scale;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2250]** `	bool initiation_point_refine_unfinished = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2251]** `	while (initiation_point_refine_unfinished)`：循环推进（通常用于时间步或迭代收敛过程）。
- **[行 2252]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 2253]** `	    initiation_point_refine_unfinished = false;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2254]** `	    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2255]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 2256]** `		if (    std::fabs(cell->center()[0] - 0.5) < 0.025`：条件分支：根据当前状态选择执行路径。
- **[行 2257]** `		     && cell->center()[1] < 0.0 && cell->center()[1] > -0.025)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2258]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 2259]** `		    material_id = cell->material_id();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2260]** `		    length_scale = m_material_data[material_id][2];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2261]** `		    if (  std::sqrt(cell->measure())`：条件分支：根据当前状态选择执行路径。
- **[行 2262]** `			> length_scale * m_parameters.m_allowed_max_h_l_ratio )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2263]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 2264]** `		        cell->set_refine_flag();`：函数调用语句，触发对应计算或操作。
- **[行 2265]** `		        initiation_point_refine_unfinished = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2266]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 2267]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 2268]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 2269]** `	    m_triangulation.execute_coarsening_and_refinement();`：函数调用语句，触发对应计算或操作。
- **[行 2270]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 2271]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2272]** `    else`：条件分支的兜底路径。
- **[行 2273]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2274]** `	AssertThrow(false,`：运行期断言/检查，验证输入与状态合法性。
- **[行 2275]** `	            ExcMessage("Selected mesh refinement strategy not implemented!"));`：函数调用语句，触发对应计算或操作。
- **[行 2276]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2277]** `  }`：作用域边界（代码块开始/结束）。
- **[行 2278]** ``：空行，用于分隔逻辑块。
- **[行 2279]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 2280]** `  void PhaseFieldMonolithicSolve<dim>::make_grid_case_3()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2281]** `  {`：作用域边界（代码块开始/结束）。
- **[行 2282]** `    for (unsigned int i = 0; i < 80; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2283]** `      m_logfile << "*";`：写日志输出，记录当前计算状态与结果。
- **[行 2284]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2285]** `    m_logfile << "\t\t\tSquare tension (structured)" << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2286]** `    for (unsigned int i = 0; i < 80; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2287]** `      m_logfile << "*";`：写日志输出，记录当前计算状态与结果。
- **[行 2288]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2289]** ``：空行，用于分隔逻辑块。
- **[行 2290]** `    AssertThrow(dim==2, ExcMessage("The dimension has to be 2D!"));`：运行期断言/检查，验证输入与状态合法性。
- **[行 2291]** ``：空行，用于分隔逻辑块。
- **[行 2292]** `    GridIn<dim> gridin;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2293]** `    gridin.attach_triangulation(m_triangulation);`：函数调用语句，触发对应计算或操作。
- **[行 2294]** `    std::ifstream f("square_tension_structured.msh");`：调用 C++ 标准库工具函数/容器接口。
- **[行 2295]** `    gridin.read_msh(f);`：函数调用语句，触发对应计算或操作。
- **[行 2296]** ``：空行，用于分隔逻辑块。
- **[行 2297]** `    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2298]** `      for (const auto &face : cell->face_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2299]** `	{`：作用域边界（代码块开始/结束）。
- **[行 2300]** `	  if (face->at_boundary() == true)`：条件分支：根据当前状态选择执行路径。
- **[行 2301]** `	    {`：作用域边界（代码块开始/结束）。
- **[行 2302]** `	      if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9 )`：条件分支：根据当前状态选择执行路径。
- **[行 2303]** `		face->set_boundary_id(0);`：函数调用语句，触发对应计算或操作。
- **[行 2304]** `	      else if (std::fabs(face->center()[1] - 1.0 ) < 1.0e-9)`：多分支条件判断，处理备选情形。
- **[行 2305]** `	        face->set_boundary_id(1);`：函数调用语句，触发对应计算或操作。
- **[行 2306]** `	      else`：条件分支的兜底路径。
- **[行 2307]** `	        face->set_boundary_id(2);`：函数调用语句，触发对应计算或操作。
- **[行 2308]** `	    }`：作用域边界（代码块开始/结束）。
- **[行 2309]** `	}`：作用域边界（代码块开始/结束）。
- **[行 2310]** ``：空行，用于分隔逻辑块。
- **[行 2311]** `    m_triangulation.refine_global(m_parameters.m_global_refine_times);`：函数调用语句，触发对应计算或操作。
- **[行 2312]** ``：空行，用于分隔逻辑块。
- **[行 2313]** `    if (m_parameters.m_refinement_strategy == "pre-refine")`：条件分支：根据当前状态选择执行路径。
- **[行 2314]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2315]** `	unsigned int material_id;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2316]** `	double length_scale;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2317]** `	for (unsigned int i = 0; i < m_parameters.m_local_prerefine_times; i++)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2318]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 2319]** `	    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2320]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 2321]** `		if (    (std::fabs(cell->center()[1] - 0.5) < 0.025)`：条件分支：根据当前状态选择执行路径。
- **[行 2322]** `		     && (cell->center()[0] > 0.475) )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2323]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 2324]** `		    material_id = cell->material_id();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2325]** `		    length_scale = m_material_data[material_id][2];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2326]** `		    if (  std::sqrt(cell->measure())`：条件分支：根据当前状态选择执行路径。
- **[行 2327]** `			> length_scale * m_parameters.m_allowed_max_h_l_ratio )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2328]** `		      cell->set_refine_flag();`：函数调用语句，触发对应计算或操作。
- **[行 2329]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 2330]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 2331]** `	    m_triangulation.execute_coarsening_and_refinement();`：函数调用语句，触发对应计算或操作。
- **[行 2332]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 2333]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2334]** `    else if (m_parameters.m_refinement_strategy == "adaptive-refine")`：多分支条件判断，处理备选情形。
- **[行 2335]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2336]** `	unsigned int material_id;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2337]** `	double length_scale;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2338]** `	bool initiation_point_refine_unfinished = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2339]** `	while (initiation_point_refine_unfinished)`：循环推进（通常用于时间步或迭代收敛过程）。
- **[行 2340]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 2341]** `	    initiation_point_refine_unfinished = false;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2342]** `	    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2343]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 2344]** `		if (    std::fabs(cell->center()[0] - 0.5) < 0.025`：条件分支：根据当前状态选择执行路径。
- **[行 2345]** `		     && std::fabs(cell->center()[1] - 0.5) < 0.025 )`：调用 C++ 标准库工具函数/容器接口。
- **[行 2346]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 2347]** `		    material_id = cell->material_id();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2348]** `		    length_scale = m_material_data[material_id][2];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2349]** `		    if (  std::sqrt(cell->measure())`：条件分支：根据当前状态选择执行路径。
- **[行 2350]** `			> length_scale * m_parameters.m_allowed_max_h_l_ratio )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2351]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 2352]** `		        cell->set_refine_flag();`：函数调用语句，触发对应计算或操作。
- **[行 2353]** `		        initiation_point_refine_unfinished = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2354]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 2355]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 2356]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 2357]** `	    m_triangulation.execute_coarsening_and_refinement();`：函数调用语句，触发对应计算或操作。
- **[行 2358]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 2359]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2360]** `    else`：条件分支的兜底路径。
- **[行 2361]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2362]** `	AssertThrow(false,`：运行期断言/检查，验证输入与状态合法性。
- **[行 2363]** `	            ExcMessage("Selected mesh refinement strategy not implemented!"));`：函数调用语句，触发对应计算或操作。
- **[行 2364]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2365]** `  }`：作用域边界（代码块开始/结束）。
- **[行 2366]** ``：空行，用于分隔逻辑块。
- **[行 2367]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 2368]** `  void PhaseFieldMonolithicSolve<dim>::make_grid_case_4()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2369]** `  {`：作用域边界（代码块开始/结束）。
- **[行 2370]** `    for (unsigned int i = 0; i < 80; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2371]** `      m_logfile << "*";`：写日志输出，记录当前计算状态与结果。
- **[行 2372]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2373]** `    m_logfile << "\t\t\t\tSquare shear (structured)" << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2374]** `    for (unsigned int i = 0; i < 80; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2375]** `      m_logfile << "*";`：写日志输出，记录当前计算状态与结果。
- **[行 2376]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2377]** ``：空行，用于分隔逻辑块。
- **[行 2378]** `    AssertThrow(dim==2, ExcMessage("The dimension has to be 2D!"));`：运行期断言/检查，验证输入与状态合法性。
- **[行 2379]** ``：空行，用于分隔逻辑块。
- **[行 2380]** `    GridIn<dim> gridin;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2381]** `    gridin.attach_triangulation(m_triangulation);`：函数调用语句，触发对应计算或操作。
- **[行 2382]** `    std::ifstream f("square_shear_structured.msh");`：调用 C++ 标准库工具函数/容器接口。
- **[行 2383]** `    gridin.read_msh(f);`：函数调用语句，触发对应计算或操作。
- **[行 2384]** ``：空行，用于分隔逻辑块。
- **[行 2385]** `    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2386]** `      for (const auto &face : cell->face_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2387]** `	{`：作用域边界（代码块开始/结束）。
- **[行 2388]** `	  if (face->at_boundary() == true)`：条件分支：根据当前状态选择执行路径。
- **[行 2389]** `	    {`：作用域边界（代码块开始/结束）。
- **[行 2390]** `	      if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9 )`：条件分支：根据当前状态选择执行路径。
- **[行 2391]** `		face->set_boundary_id(0);`：函数调用语句，触发对应计算或操作。
- **[行 2392]** `	      else if (std::fabs(face->center()[1] - 1.0 ) < 1.0e-9)`：多分支条件判断，处理备选情形。
- **[行 2393]** `	        face->set_boundary_id(1);`：函数调用语句，触发对应计算或操作。
- **[行 2394]** `	      else if (   (std::fabs(face->center()[0] - 0.0 ) < 1.0e-9)`：多分支条件判断，处理备选情形。
- **[行 2395]** `		       || (std::fabs(face->center()[0] - 1.0 ) < 1.0e-9))`：调用 C++ 标准库工具函数/容器接口。
- **[行 2396]** `	        face->set_boundary_id(2);`：函数调用语句，触发对应计算或操作。
- **[行 2397]** `	      else`：条件分支的兜底路径。
- **[行 2398]** `	        face->set_boundary_id(3);`：函数调用语句，触发对应计算或操作。
- **[行 2399]** `	    }`：作用域边界（代码块开始/结束）。
- **[行 2400]** `	}`：作用域边界（代码块开始/结束）。
- **[行 2401]** ``：空行，用于分隔逻辑块。
- **[行 2402]** `    m_triangulation.refine_global(m_parameters.m_global_refine_times);`：函数调用语句，触发对应计算或操作。
- **[行 2403]** ``：空行，用于分隔逻辑块。
- **[行 2404]** `    if (m_parameters.m_refinement_strategy == "pre-refine")`：条件分支：根据当前状态选择执行路径。
- **[行 2405]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2406]** `	unsigned int material_id;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2407]** `	double length_scale;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2408]** `	for (unsigned int i = 0; i < m_parameters.m_local_prerefine_times; i++)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2409]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 2410]** `	    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2411]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 2412]** `		if (   (   (cell->center()[0] > 0.475)`：条件分支：根据当前状态选择执行路径。
- **[行 2413]** `		        && (cell->center()[1] < 0.525) )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2414]** `		    //|| (    cell->center()[1] > 0.975)`：注释行，用于说明算法背景、假设或实现意图。
- **[行 2415]** `		   )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2416]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 2417]** `		    material_id = cell->material_id();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2418]** `		    length_scale = m_material_data[material_id][2];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2419]** `		    if (  std::sqrt(cell->measure())`：条件分支：根据当前状态选择执行路径。
- **[行 2420]** `			> length_scale * m_parameters.m_allowed_max_h_l_ratio )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2421]** `		      cell->set_refine_flag();`：函数调用语句，触发对应计算或操作。
- **[行 2422]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 2423]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 2424]** `	    m_triangulation.execute_coarsening_and_refinement();`：函数调用语句，触发对应计算或操作。
- **[行 2425]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 2426]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2427]** `    else if (m_parameters.m_refinement_strategy == "adaptive-refine")`：多分支条件判断，处理备选情形。
- **[行 2428]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2429]** `	unsigned int material_id;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2430]** `	double length_scale;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2431]** `	bool initiation_point_refine_unfinished = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2432]** `	while (initiation_point_refine_unfinished)`：循环推进（通常用于时间步或迭代收敛过程）。
- **[行 2433]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 2434]** `	    initiation_point_refine_unfinished = false;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2435]** `	    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2436]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 2437]** `		if (   (   std::fabs(cell->center()[0] - 0.5) < 0.025`：条件分支：根据当前状态选择执行路径。
- **[行 2438]** `		        && cell->center()[1] < 0.5 && cell->center()[1] > 0.475 )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2439]** `		    // we also need to refine the top edge, since there might be a conflict between`：注释行，用于说明算法背景、假设或实现意图。
- **[行 2440]** `		    // inhomogeneous boundary conditions and the hanging-node constraints at the`：注释行，用于说明算法背景、假设或实现意图。
- **[行 2441]** `		    // top edge`：注释行，用于说明算法背景、假设或实现意图。
- **[行 2442]** `		    || (   cell->center()[1] > 0.975 ) )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2443]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 2444]** `		    material_id = cell->material_id();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2445]** `		    length_scale = m_material_data[material_id][2];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2446]** `		    if (  std::sqrt(cell->measure())`：条件分支：根据当前状态选择执行路径。
- **[行 2447]** `			> length_scale * m_parameters.m_allowed_max_h_l_ratio )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2448]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 2449]** `		        cell->set_refine_flag();`：函数调用语句，触发对应计算或操作。
- **[行 2450]** `		        initiation_point_refine_unfinished = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2451]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 2452]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 2453]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 2454]** `	    m_triangulation.execute_coarsening_and_refinement();`：函数调用语句，触发对应计算或操作。
- **[行 2455]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 2456]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2457]** `    else`：条件分支的兜底路径。
- **[行 2458]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2459]** `	AssertThrow(false,`：运行期断言/检查，验证输入与状态合法性。
- **[行 2460]** `	            ExcMessage("Selected mesh refinement strategy not implemented!"));`：函数调用语句，触发对应计算或操作。
- **[行 2461]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2462]** `  }`：作用域边界（代码块开始/结束）。
- **[行 2463]** ``：空行，用于分隔逻辑块。
- **[行 2464]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 2465]** `  void PhaseFieldMonolithicSolve<dim>::make_grid_case_5()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2466]** `  {`：作用域边界（代码块开始/结束）。
- **[行 2467]** `    for (unsigned int i = 0; i < 80; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2468]** `      m_logfile << "*";`：写日志输出，记录当前计算状态与结果。
- **[行 2469]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2470]** `    m_logfile << "\t\t\t\tThree-point bending (structured)" << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2471]** `    for (unsigned int i = 0; i < 80; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2472]** `      m_logfile << "*";`：写日志输出，记录当前计算状态与结果。
- **[行 2473]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2474]** ``：空行，用于分隔逻辑块。
- **[行 2475]** `    AssertThrow(dim==2, ExcMessage("The dimension has to be 2D!"));`：运行期断言/检查，验证输入与状态合法性。
- **[行 2476]** ``：空行，用于分隔逻辑块。
- **[行 2477]** `    GridIn<dim> gridin;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2478]** `    gridin.attach_triangulation(m_triangulation);`：函数调用语句，触发对应计算或操作。
- **[行 2479]** `    std::ifstream f("three_point_bending_structured.msh");`：调用 C++ 标准库工具函数/容器接口。
- **[行 2480]** `    gridin.read_msh(f);`：函数调用语句，触发对应计算或操作。
- **[行 2481]** ``：空行，用于分隔逻辑块。
- **[行 2482]** `    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2483]** `      for (const auto &face : cell->face_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2484]** `	{`：作用域边界（代码块开始/结束）。
- **[行 2485]** `	  if (face->at_boundary() == true)`：条件分支：根据当前状态选择执行路径。
- **[行 2486]** `	    {`：作用域边界（代码块开始/结束）。
- **[行 2487]** `	      if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9 )`：条件分支：根据当前状态选择执行路径。
- **[行 2488]** `		face->set_boundary_id(0);`：函数调用语句，触发对应计算或操作。
- **[行 2489]** `	      else if (std::fabs(face->center()[1] - 2.0 ) < 1.0e-9)`：多分支条件判断，处理备选情形。
- **[行 2490]** `	        face->set_boundary_id(1);`：函数调用语句，触发对应计算或操作。
- **[行 2491]** `	      else`：条件分支的兜底路径。
- **[行 2492]** `	        face->set_boundary_id(2);`：函数调用语句，触发对应计算或操作。
- **[行 2493]** `	    }`：作用域边界（代码块开始/结束）。
- **[行 2494]** `	}`：作用域边界（代码块开始/结束）。
- **[行 2495]** ``：空行，用于分隔逻辑块。
- **[行 2496]** `    m_triangulation.refine_global(m_parameters.m_global_refine_times);`：函数调用语句，触发对应计算或操作。
- **[行 2497]** ``：空行，用于分隔逻辑块。
- **[行 2498]** `    if (m_parameters.m_refinement_strategy == "pre-refine")`：条件分支：根据当前状态选择执行路径。
- **[行 2499]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2500]** `	for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2501]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 2502]** `	    if (    std::fabs(cell->center()[0] - 4.0) < 0.075`：条件分支：根据当前状态选择执行路径。
- **[行 2503]** `		 && cell->center()[1] < 1.6)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2504]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 2505]** `		cell->set_refine_flag();`：函数调用语句，触发对应计算或操作。
- **[行 2506]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 2507]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 2508]** `	m_triangulation.execute_coarsening_and_refinement();`：函数调用语句，触发对应计算或操作。
- **[行 2509]** ``：空行，用于分隔逻辑块。
- **[行 2510]** `	unsigned int material_id;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2511]** `	double length_scale;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2512]** `	for (unsigned int i = 0; i < m_parameters.m_local_prerefine_times; i++)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2513]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 2514]** `	    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2515]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 2516]** `		if (    std::fabs(cell->center()[0] - 4.0) < 0.05`：条件分支：根据当前状态选择执行路径。
- **[行 2517]** `		     && cell->center()[1] < 1.6)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2518]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 2519]** `		    material_id = cell->material_id();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2520]** `		    length_scale = m_material_data[material_id][2];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2521]** `		    if (  std::sqrt(cell->measure())`：条件分支：根据当前状态选择执行路径。
- **[行 2522]** `			> length_scale * m_parameters.m_allowed_max_h_l_ratio )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2523]** `		      cell->set_refine_flag();`：函数调用语句，触发对应计算或操作。
- **[行 2524]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 2525]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 2526]** `	    m_triangulation.execute_coarsening_and_refinement();`：函数调用语句，触发对应计算或操作。
- **[行 2527]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 2528]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2529]** `    else if (m_parameters.m_refinement_strategy == "adaptive-refine")`：多分支条件判断，处理备选情形。
- **[行 2530]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2531]** `	unsigned int material_id;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2532]** `	double length_scale;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2533]** `	bool initiation_point_refine_unfinished = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2534]** `	while (initiation_point_refine_unfinished)`：循环推进（通常用于时间步或迭代收敛过程）。
- **[行 2535]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 2536]** `	    initiation_point_refine_unfinished = false;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2537]** `	    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2538]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 2539]** `		if (    std::fabs(cell->center()[0] - 4.0) < 0.075`：条件分支：根据当前状态选择执行路径。
- **[行 2540]** `		     && std::fabs(cell->center()[1] - 0.4) < 0.075 )`：调用 C++ 标准库工具函数/容器接口。
- **[行 2541]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 2542]** `		    material_id = cell->material_id();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2543]** `		    length_scale = m_material_data[material_id][2];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2544]** `		    if (  std::sqrt(cell->measure())`：条件分支：根据当前状态选择执行路径。
- **[行 2545]** `			> length_scale * m_parameters.m_allowed_max_h_l_ratio )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2546]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 2547]** `		        cell->set_refine_flag();`：函数调用语句，触发对应计算或操作。
- **[行 2548]** `		        initiation_point_refine_unfinished = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2549]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 2550]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 2551]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 2552]** `	    m_triangulation.execute_coarsening_and_refinement();`：函数调用语句，触发对应计算或操作。
- **[行 2553]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 2554]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2555]** `    else`：条件分支的兜底路径。
- **[行 2556]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2557]** `	AssertThrow(false,`：运行期断言/检查，验证输入与状态合法性。
- **[行 2558]** `	            ExcMessage("Selected mesh refinement strategy not implemented!"));`：函数调用语句，触发对应计算或操作。
- **[行 2559]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2560]** `  }`：作用域边界（代码块开始/结束）。
- **[行 2561]** ``：空行，用于分隔逻辑块。
- **[行 2562]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 2563]** `  void PhaseFieldMonolithicSolve<dim>::make_grid_case_6()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2564]** `  {`：作用域边界（代码块开始/结束）。
- **[行 2565]** `    AssertThrow(dim==3, ExcMessage("The dimension has to be 3D!"));`：运行期断言/检查，验证输入与状态合法性。
- **[行 2566]** ``：空行，用于分隔逻辑块。
- **[行 2567]** `    for (unsigned int i = 0; i < 80; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2568]** `      m_logfile << "*";`：写日志输出，记录当前计算状态与结果。
- **[行 2569]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2570]** `    m_logfile << "\t\t\t\tSphere inclusion (3D structured)" << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2571]** `    for (unsigned int i = 0; i < 80; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2572]** `      m_logfile << "*";`：写日志输出，记录当前计算状态与结果。
- **[行 2573]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2574]** ``：空行，用于分隔逻辑块。
- **[行 2575]** `    Triangulation<dim> tria_inner;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2576]** `    GridGenerator::hyper_ball(tria_inner, Point<dim>(), 0.5);`：函数调用语句，触发对应计算或操作。
- **[行 2577]** ``：空行，用于分隔逻辑块。
- **[行 2578]** `    Triangulation<dim> tria_outer;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2579]** `    GridGenerator::hyper_shell(`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2580]** `      tria_outer, Point<dim>(), 0.5, std::sqrt(dim), 2 * dim);`：调用 C++ 标准库工具函数/容器接口。
- **[行 2581]** ``：空行，用于分隔逻辑块。
- **[行 2582]** `    Triangulation<dim> tmp_triangulation;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2583]** ``：空行，用于分隔逻辑块。
- **[行 2584]** `    GridGenerator::merge_triangulations(tria_inner, tria_outer, tmp_triangulation);`：函数调用语句，触发对应计算或操作。
- **[行 2585]** ``：空行，用于分隔逻辑块。
- **[行 2586]** `    tmp_triangulation.reset_all_manifolds();`：函数调用语句，触发对应计算或操作。
- **[行 2587]** `    tmp_triangulation.set_all_manifold_ids(0);`：函数调用语句，触发对应计算或操作。
- **[行 2588]** ``：空行，用于分隔逻辑块。
- **[行 2589]** `    for (const auto &cell : tmp_triangulation.cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2590]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2591]** `        for (const auto &face : cell->face_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2592]** `          {`：作用域边界（代码块开始/结束）。
- **[行 2593]** `            bool face_at_sphere_boundary = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2594]** `            for (const auto v : face->vertex_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2595]** `              {`：作用域边界（代码块开始/结束）。
- **[行 2596]** `                if (std::abs(face->vertex(v).norm_square() - 0.25) > 1e-12)`：条件分支：根据当前状态选择执行路径。
- **[行 2597]** `                  face_at_sphere_boundary = false;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2598]** `              }`：作用域边界（代码块开始/结束）。
- **[行 2599]** `            if (face_at_sphere_boundary)`：条件分支：根据当前状态选择执行路径。
- **[行 2600]** `              face->set_all_manifold_ids(1);`：函数调用语句，触发对应计算或操作。
- **[行 2601]** `          }`：作用域边界（代码块开始/结束）。
- **[行 2602]** `        if (cell->center().norm_square() < 0.25)`：条件分支：根据当前状态选择执行路径。
- **[行 2603]** `          cell->set_material_id(1);`：函数调用语句，触发对应计算或操作。
- **[行 2604]** `        else`：条件分支的兜底路径。
- **[行 2605]** `          cell->set_material_id(0);`：函数调用语句，触发对应计算或操作。
- **[行 2606]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2607]** ``：空行，用于分隔逻辑块。
- **[行 2608]** `    tmp_triangulation.set_manifold(1, SphericalManifold<dim>());`：函数调用语句，触发对应计算或操作。
- **[行 2609]** ``：空行，用于分隔逻辑块。
- **[行 2610]** `    TransfiniteInterpolationManifold<dim> transfinite_manifold;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2611]** `    transfinite_manifold.initialize(tmp_triangulation);`：函数调用语句，触发对应计算或操作。
- **[行 2612]** `    tmp_triangulation.set_manifold(0, transfinite_manifold);`：函数调用语句，触发对应计算或操作。
- **[行 2613]** ``：空行，用于分隔逻辑块。
- **[行 2614]** `    tmp_triangulation.refine_global(m_parameters.m_global_refine_times);`：函数调用语句，触发对应计算或操作。
- **[行 2615]** ``：空行，用于分隔逻辑块。
- **[行 2616]** `    std::set<typename Triangulation< dim >::active_cell_iterator >`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2617]** `      cells_to_remove;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2618]** ``：空行，用于分隔逻辑块。
- **[行 2619]** `    for (const auto &cell : tmp_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2620]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2621]** `	if (   cell->center()[0] < 0.0`：条件分支：根据当前状态选择执行路径。
- **[行 2622]** `	    || cell->center()[1] < 0.0`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2623]** `	    || cell->center()[2] < 0.0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2624]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 2625]** `	    cells_to_remove.insert(cell);`：函数调用语句，触发对应计算或操作。
- **[行 2626]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 2627]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2628]** ``：空行，用于分隔逻辑块。
- **[行 2629]** `    GridGenerator::create_triangulation_with_removed_cells(tmp_triangulation,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2630]** `							   cells_to_remove,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2631]** `							   m_triangulation);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2632]** ``：空行，用于分隔逻辑块。
- **[行 2633]** `    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2634]** `      for (const auto &face : cell->face_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2635]** `	{`：作用域边界（代码块开始/结束）。
- **[行 2636]** `	  if (face->at_boundary() == true)`：条件分支：根据当前状态选择执行路径。
- **[行 2637]** `	    {`：作用域边界（代码块开始/结束）。
- **[行 2638]** `	      if (std::fabs(face->center()[0] - 0.0 ) < 1.0e-9 )`：条件分支：根据当前状态选择执行路径。
- **[行 2639]** `		face->set_boundary_id(0);`：函数调用语句，触发对应计算或操作。
- **[行 2640]** `	      else if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9)`：多分支条件判断，处理备选情形。
- **[行 2641]** `		face->set_boundary_id(1);`：函数调用语句，触发对应计算或操作。
- **[行 2642]** `	      else if (std::fabs(face->center()[2] - 0.0 ) < 1.0e-9)`：多分支条件判断，处理备选情形。
- **[行 2643]** `		face->set_boundary_id(2);`：函数调用语句，触发对应计算或操作。
- **[行 2644]** `	      else if (std::fabs(face->center()[2] - 1.0 ) < 1.0e-9)`：多分支条件判断，处理备选情形。
- **[行 2645]** `		face->set_boundary_id(3);`：函数调用语句，触发对应计算或操作。
- **[行 2646]** `	      else`：条件分支的兜底路径。
- **[行 2647]** `		face->set_boundary_id(4);`：函数调用语句，触发对应计算或操作。
- **[行 2648]** `	    }`：作用域边界（代码块开始/结束）。
- **[行 2649]** `	}`：作用域边界（代码块开始/结束）。
- **[行 2650]** ``：空行，用于分隔逻辑块。
- **[行 2651]** `    if (m_parameters.m_refinement_strategy == "adaptive-refine")`：条件分支：根据当前状态选择执行路径。
- **[行 2652]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2653]** `	unsigned int material_id;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2654]** `	double length_scale;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2655]** `	bool initiation_point_refine_unfinished = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2656]** `	while (initiation_point_refine_unfinished)`：循环推进（通常用于时间步或迭代收敛过程）。
- **[行 2657]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 2658]** `	    initiation_point_refine_unfinished = false;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2659]** `	    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2660]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 2661]** `		if (    cell->center()[2] > 0.525`：条件分支：根据当前状态选择执行路径。
- **[行 2662]** `		     && cell->center()[2] < 0.575`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2663]** `		     && cell->center()[0] < 0.05`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2664]** `		     && cell->center()[1] < 0.05 )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2665]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 2666]** `		    material_id = cell->material_id();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2667]** `		    length_scale = m_material_data[material_id][2];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2668]** `		    if (  std::cbrt(cell->measure())`：条件分支：根据当前状态选择执行路径。
- **[行 2669]** `			> length_scale * m_parameters.m_allowed_max_h_l_ratio )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2670]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 2671]** `			cell->set_refine_flag();`：函数调用语句，触发对应计算或操作。
- **[行 2672]** `			initiation_point_refine_unfinished = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2673]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 2674]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 2675]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 2676]** `	    m_triangulation.execute_coarsening_and_refinement();`：函数调用语句，触发对应计算或操作。
- **[行 2677]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 2678]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2679]** `    else`：条件分支的兜底路径。
- **[行 2680]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2681]** `	AssertThrow(false,`：运行期断言/检查，验证输入与状态合法性。
- **[行 2682]** `		    ExcMessage("Selected mesh refinement strategy not implemented!"));`：函数调用语句，触发对应计算或操作。
- **[行 2683]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2684]** `  }`：作用域边界（代码块开始/结束）。
- **[行 2685]** ``：空行，用于分隔逻辑块。
- **[行 2686]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 2687]** `  void PhaseFieldMonolithicSolve<dim>::make_grid_case_7()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2688]** `  {`：作用域边界（代码块开始/结束）。
- **[行 2689]** `    AssertThrow(dim==3, ExcMessage("The dimension has to be 3D!"));`：运行期断言/检查，验证输入与状态合法性。
- **[行 2690]** ``：空行，用于分隔逻辑块。
- **[行 2691]** `    for (unsigned int i = 0; i < 80; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2692]** `      m_logfile << "*";`：写日志输出，记录当前计算状态与结果。
- **[行 2693]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2694]** `    m_logfile << "\t\t\t\tSphere inclusion (3D structured version 2)" << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2695]** `    for (unsigned int i = 0; i < 80; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2696]** `      m_logfile << "*";`：写日志输出，记录当前计算状态与结果。
- **[行 2697]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2698]** ``：空行，用于分隔逻辑块。
- **[行 2699]** `    Triangulation<dim> tria_inner;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2700]** `    GridGenerator::hyper_ball(tria_inner, Point<dim>(), 0.49);`：函数调用语句，触发对应计算或操作。
- **[行 2701]** ``：空行，用于分隔逻辑块。
- **[行 2702]** `    Triangulation<dim> tria_outer;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2703]** `    GridGenerator::hyper_shell(`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2704]** `      tria_outer, Point<dim>(), 0.49, std::sqrt(dim)*0.5, 2 * dim);`：调用 C++ 标准库工具函数/容器接口。
- **[行 2705]** ``：空行，用于分隔逻辑块。
- **[行 2706]** `    Triangulation<dim> cube1;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2707]** `    GridGenerator::hyper_rectangle(cube1, Point<dim>(0, 0, 0.5), Point<dim>(1, 1, 1.5));`：函数调用语句，触发对应计算或操作。
- **[行 2708]** `    Triangulation<dim> cube2;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2709]** `    GridGenerator::hyper_rectangle(cube2, Point<dim>(0, 0.5, -0.5), Point<dim>(1, 1.5, 0.5));`：函数调用语句，触发对应计算或操作。
- **[行 2710]** `    Triangulation<dim> cube3;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2711]** `    GridGenerator::hyper_rectangle(cube3, Point<dim>(0.5, -0.5, -0.5), Point<dim>(1.5, 0.5, 0.5));`：函数调用语句，触发对应计算或操作。
- **[行 2712]** ``：空行，用于分隔逻辑块。
- **[行 2713]** `    Triangulation<dim> tmp_triangulation;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2714]** `    GridGenerator::merge_triangulations({&tria_inner, &tria_outer,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2715]** `                                         &cube1, &cube2, &cube3}, tmp_triangulation);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2716]** ``：空行，用于分隔逻辑块。
- **[行 2717]** `    tmp_triangulation.reset_all_manifolds();`：函数调用语句，触发对应计算或操作。
- **[行 2718]** `    tmp_triangulation.set_all_manifold_ids(0);`：函数调用语句，触发对应计算或操作。
- **[行 2719]** ``：空行，用于分隔逻辑块。
- **[行 2720]** `    for (const auto &cell : tmp_triangulation.cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2721]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2722]** `        for (const auto &face : cell->face_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2723]** `          {`：作用域边界（代码块开始/结束）。
- **[行 2724]** `            bool face_at_sphere_boundary = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2725]** `            for (const auto v : face->vertex_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2726]** `              {`：作用域边界（代码块开始/结束）。
- **[行 2727]** `                if (std::abs(face->vertex(v).norm_square() - 0.49 * 0.49) > 1e-12)`：条件分支：根据当前状态选择执行路径。
- **[行 2728]** `                  face_at_sphere_boundary = false;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2729]** `              }`：作用域边界（代码块开始/结束）。
- **[行 2730]** `            if (face_at_sphere_boundary)`：条件分支：根据当前状态选择执行路径。
- **[行 2731]** `              face->set_all_manifold_ids(1);`：函数调用语句，触发对应计算或操作。
- **[行 2732]** `          }`：作用域边界（代码块开始/结束）。
- **[行 2733]** `        if (cell->center().norm_square() < 0.1)`：条件分支：根据当前状态选择执行路径。
- **[行 2734]** `          cell->set_material_id(1);`：函数调用语句，触发对应计算或操作。
- **[行 2735]** `        else`：条件分支的兜底路径。
- **[行 2736]** `          cell->set_material_id(0);`：函数调用语句，触发对应计算或操作。
- **[行 2737]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2738]** ``：空行，用于分隔逻辑块。
- **[行 2739]** `    tmp_triangulation.set_manifold(1, SphericalManifold<dim>());`：函数调用语句，触发对应计算或操作。
- **[行 2740]** ``：空行，用于分隔逻辑块。
- **[行 2741]** `    TransfiniteInterpolationManifold<dim> transfinite_manifold;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2742]** `    transfinite_manifold.initialize(tmp_triangulation);`：函数调用语句，触发对应计算或操作。
- **[行 2743]** `    tmp_triangulation.set_manifold(0, transfinite_manifold);`：函数调用语句，触发对应计算或操作。
- **[行 2744]** ``：空行，用于分隔逻辑块。
- **[行 2745]** `    tmp_triangulation.refine_global(m_parameters.m_global_refine_times);`：函数调用语句，触发对应计算或操作。
- **[行 2746]** ``：空行，用于分隔逻辑块。
- **[行 2747]** `    std::set<typename Triangulation< dim >::active_cell_iterator >`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2748]** `      cells_to_remove;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2749]** ``：空行，用于分隔逻辑块。
- **[行 2750]** `    for (const auto &cell : tmp_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2751]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2752]** `	if (   cell->center()[0] < 0.0`：条件分支：根据当前状态选择执行路径。
- **[行 2753]** `	    || cell->center()[1] < 0.0`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2754]** `	    || cell->center()[2] < 0.0`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2755]** `	    || cell->center()[0] > 1.0`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2756]** `	    || cell->center()[1] > 1.0`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2757]** `	    || cell->center()[2] > 1.0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2758]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 2759]** `	    cells_to_remove.insert(cell);`：函数调用语句，触发对应计算或操作。
- **[行 2760]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 2761]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2762]** ``：空行，用于分隔逻辑块。
- **[行 2763]** `    GridGenerator::create_triangulation_with_removed_cells(tmp_triangulation,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2764]** `							   cells_to_remove,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2765]** `							   m_triangulation);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2766]** ``：空行，用于分隔逻辑块。
- **[行 2767]** `    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2768]** `      for (const auto &face : cell->face_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2769]** `	{`：作用域边界（代码块开始/结束）。
- **[行 2770]** `	  if (face->at_boundary() == true)`：条件分支：根据当前状态选择执行路径。
- **[行 2771]** `	    {`：作用域边界（代码块开始/结束）。
- **[行 2772]** `	      if (std::fabs(face->center()[0] - 0.0 ) < 1.0e-9 )`：条件分支：根据当前状态选择执行路径。
- **[行 2773]** `		face->set_boundary_id(0);`：函数调用语句，触发对应计算或操作。
- **[行 2774]** `	      else if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9)`：多分支条件判断，处理备选情形。
- **[行 2775]** `		face->set_boundary_id(1);`：函数调用语句，触发对应计算或操作。
- **[行 2776]** `	      else if (std::fabs(face->center()[2] - 0.0 ) < 1.0e-9)`：多分支条件判断，处理备选情形。
- **[行 2777]** `		face->set_boundary_id(2);`：函数调用语句，触发对应计算或操作。
- **[行 2778]** `	      else if (std::fabs(face->center()[2] - 1.0 ) < 1.0e-9)`：多分支条件判断，处理备选情形。
- **[行 2779]** `		face->set_boundary_id(3);`：函数调用语句，触发对应计算或操作。
- **[行 2780]** `	      else`：条件分支的兜底路径。
- **[行 2781]** `		face->set_boundary_id(4);`：函数调用语句，触发对应计算或操作。
- **[行 2782]** `	    }`：作用域边界（代码块开始/结束）。
- **[行 2783]** `	}`：作用域边界（代码块开始/结束）。
- **[行 2784]** ``：空行，用于分隔逻辑块。
- **[行 2785]** `    if (m_parameters.m_refinement_strategy == "adaptive-refine")`：条件分支：根据当前状态选择执行路径。
- **[行 2786]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2787]** `	unsigned int material_id;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2788]** `	double length_scale;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2789]** `	bool initiation_point_refine_unfinished = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2790]** `	while (initiation_point_refine_unfinished)`：循环推进（通常用于时间步或迭代收敛过程）。
- **[行 2791]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 2792]** `	    initiation_point_refine_unfinished = false;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2793]** `	    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2794]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 2795]** `		if (    cell->center()[2] > 0.505`：条件分支：根据当前状态选择执行路径。
- **[行 2796]** `		     && cell->center()[2] < 0.575`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2797]** `		     && cell->center()[0] < 0.05`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2798]** `		     && cell->center()[1] < 0.05 )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2799]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 2800]** `		    material_id = cell->material_id();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2801]** `		    length_scale = m_material_data[material_id][2];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2802]** `		    if (  std::cbrt(cell->measure())`：条件分支：根据当前状态选择执行路径。
- **[行 2803]** `			> length_scale * m_parameters.m_allowed_max_h_l_ratio )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2804]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 2805]** `			cell->set_refine_flag();`：函数调用语句，触发对应计算或操作。
- **[行 2806]** `			initiation_point_refine_unfinished = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2807]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 2808]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 2809]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 2810]** `	    m_triangulation.execute_coarsening_and_refinement();`：函数调用语句，触发对应计算或操作。
- **[行 2811]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 2812]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2813]** `    else`：条件分支的兜底路径。
- **[行 2814]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2815]** `	AssertThrow(false,`：运行期断言/检查，验证输入与状态合法性。
- **[行 2816]** `		    ExcMessage("Selected mesh refinement strategy not implemented!"));`：函数调用语句，触发对应计算或操作。
- **[行 2817]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2818]** `  }`：作用域边界（代码块开始/结束）。
- **[行 2819]** ``：空行，用于分隔逻辑块。
- **[行 2820]** ``：空行，用于分隔逻辑块。
- **[行 2821]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 2822]** `  void PhaseFieldMonolithicSolve<dim>::make_grid_case_8()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2823]** `  {`：作用域边界（代码块开始/结束）。
- **[行 2824]** `    AssertThrow(dim==3, ExcMessage("The dimension has to be 3D!"));`：运行期断言/检查，验证输入与状态合法性。
- **[行 2825]** ``：空行，用于分隔逻辑块。
- **[行 2826]** `    for (unsigned int i = 0; i < 80; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2827]** `      m_logfile << "*";`：写日志输出，记录当前计算状态与结果。
- **[行 2828]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2829]** `    m_logfile << "\t\t\t\tSphere inclusion (3D structured version 2 with barriers)" << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2830]** `    for (unsigned int i = 0; i < 80; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2831]** `      m_logfile << "*";`：写日志输出，记录当前计算状态与结果。
- **[行 2832]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2833]** ``：空行，用于分隔逻辑块。
- **[行 2834]** `    Triangulation<dim> tria_inner;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2835]** `    GridGenerator::hyper_ball(tria_inner, Point<dim>(), 0.49);`：函数调用语句，触发对应计算或操作。
- **[行 2836]** ``：空行，用于分隔逻辑块。
- **[行 2837]** `    Triangulation<dim> tria_outer;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2838]** `    GridGenerator::hyper_shell(`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2839]** `      tria_outer, Point<dim>(), 0.49, std::sqrt(dim)*0.5, 2 * dim);`：调用 C++ 标准库工具函数/容器接口。
- **[行 2840]** ``：空行，用于分隔逻辑块。
- **[行 2841]** `    Triangulation<dim> cube1;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2842]** `    GridGenerator::hyper_rectangle(cube1, Point<dim>(0, 0, 0.5), Point<dim>(1, 1, 1.5));`：函数调用语句，触发对应计算或操作。
- **[行 2843]** `    Triangulation<dim> cube2;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2844]** `    GridGenerator::hyper_rectangle(cube2, Point<dim>(0, 0.5, -0.5), Point<dim>(1, 1.5, 0.5));`：函数调用语句，触发对应计算或操作。
- **[行 2845]** `    Triangulation<dim> cube3;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2846]** `    GridGenerator::hyper_rectangle(cube3, Point<dim>(0.5, -0.5, -0.5), Point<dim>(1.5, 0.5, 0.5));`：函数调用语句，触发对应计算或操作。
- **[行 2847]** ``：空行，用于分隔逻辑块。
- **[行 2848]** `    Triangulation<dim> tmp_triangulation;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2849]** `    GridGenerator::merge_triangulations({&tria_inner, &tria_outer,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2850]** `                                         &cube1, &cube2, &cube3}, tmp_triangulation);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2851]** ``：空行，用于分隔逻辑块。
- **[行 2852]** `    tmp_triangulation.reset_all_manifolds();`：函数调用语句，触发对应计算或操作。
- **[行 2853]** `    tmp_triangulation.set_all_manifold_ids(0);`：函数调用语句，触发对应计算或操作。
- **[行 2854]** ``：空行，用于分隔逻辑块。
- **[行 2855]** `    for (const auto &cell : tmp_triangulation.cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2856]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2857]** `        for (const auto &face : cell->face_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2858]** `          {`：作用域边界（代码块开始/结束）。
- **[行 2859]** `            bool face_at_sphere_boundary = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2860]** `            for (const auto v : face->vertex_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2861]** `              {`：作用域边界（代码块开始/结束）。
- **[行 2862]** `                if (std::abs(face->vertex(v).norm_square() - 0.49 * 0.49) > 1e-12)`：条件分支：根据当前状态选择执行路径。
- **[行 2863]** `                  face_at_sphere_boundary = false;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2864]** `              }`：作用域边界（代码块开始/结束）。
- **[行 2865]** `            if (face_at_sphere_boundary)`：条件分支：根据当前状态选择执行路径。
- **[行 2866]** `              face->set_all_manifold_ids(1);`：函数调用语句，触发对应计算或操作。
- **[行 2867]** `          }`：作用域边界（代码块开始/结束）。
- **[行 2868]** `        if (cell->center().norm_square() < 0.1)`：条件分支：根据当前状态选择执行路径。
- **[行 2869]** `          cell->set_material_id(1);`：函数调用语句，触发对应计算或操作。
- **[行 2870]** `        else`：条件分支的兜底路径。
- **[行 2871]** `          cell->set_material_id(0);`：函数调用语句，触发对应计算或操作。
- **[行 2872]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2873]** ``：空行，用于分隔逻辑块。
- **[行 2874]** `    tmp_triangulation.set_manifold(1, SphericalManifold<dim>());`：函数调用语句，触发对应计算或操作。
- **[行 2875]** ``：空行，用于分隔逻辑块。
- **[行 2876]** `    TransfiniteInterpolationManifold<dim> transfinite_manifold;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2877]** `    transfinite_manifold.initialize(tmp_triangulation);`：函数调用语句，触发对应计算或操作。
- **[行 2878]** `    tmp_triangulation.set_manifold(0, transfinite_manifold);`：函数调用语句，触发对应计算或操作。
- **[行 2879]** ``：空行，用于分隔逻辑块。
- **[行 2880]** `    tmp_triangulation.refine_global(m_parameters.m_global_refine_times);`：函数调用语句，触发对应计算或操作。
- **[行 2881]** ``：空行，用于分隔逻辑块。
- **[行 2882]** `    // some extra barriers`：注释行，用于说明算法背景、假设或实现意图。
- **[行 2883]** `    for (const auto &cell : tmp_triangulation.cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2884]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2885]** `        if (    std::fabs(cell->center()[1] - 0.75) < 0.05`：条件分支：根据当前状态选择执行路径。
- **[行 2886]** `             && std::fabs(cell->center()[2] - 0.5625) < 0.05`：调用 C++ 标准库工具函数/容器接口。
- **[行 2887]** `             && std::fabs(cell->center()[0] - 0.0) < 0.2)`：调用 C++ 标准库工具函数/容器接口。
- **[行 2888]** `          cell->set_material_id(1);`：函数调用语句，触发对应计算或操作。
- **[行 2889]** ``：空行，用于分隔逻辑块。
- **[行 2890]** `        if (    std::fabs(cell->center()[1] - 0.0) < 0.2`：条件分支：根据当前状态选择执行路径。
- **[行 2891]** `             && std::fabs(cell->center()[2] - 0.5) < 0.1`：调用 C++ 标准库工具函数/容器接口。
- **[行 2892]** `             && std::fabs(cell->center()[0] - 0.75) < 0.05)`：调用 C++ 标准库工具函数/容器接口。
- **[行 2893]** `          cell->set_material_id(1);`：函数调用语句，触发对应计算或操作。
- **[行 2894]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2895]** ``：空行，用于分隔逻辑块。
- **[行 2896]** `    std::set<typename Triangulation< dim >::active_cell_iterator >`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2897]** `      cells_to_remove;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2898]** ``：空行，用于分隔逻辑块。
- **[行 2899]** `    for (const auto &cell : tmp_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2900]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2901]** `	if (   cell->center()[0] < 0.0`：条件分支：根据当前状态选择执行路径。
- **[行 2902]** `	    || cell->center()[1] < 0.0`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2903]** `	    || cell->center()[2] < 0.0`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2904]** `	    || cell->center()[0] > 1.0`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2905]** `	    || cell->center()[1] > 1.0`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2906]** `	    || cell->center()[2] > 1.0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2907]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 2908]** `	    cells_to_remove.insert(cell);`：函数调用语句，触发对应计算或操作。
- **[行 2909]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 2910]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2911]** ``：空行，用于分隔逻辑块。
- **[行 2912]** `    GridGenerator::create_triangulation_with_removed_cells(tmp_triangulation,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2913]** `							   cells_to_remove,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2914]** `							   m_triangulation);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2915]** ``：空行，用于分隔逻辑块。
- **[行 2916]** `    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2917]** `      for (const auto &face : cell->face_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2918]** `	{`：作用域边界（代码块开始/结束）。
- **[行 2919]** `	  if (face->at_boundary() == true)`：条件分支：根据当前状态选择执行路径。
- **[行 2920]** `	    {`：作用域边界（代码块开始/结束）。
- **[行 2921]** `	      if (std::fabs(face->center()[0] - 0.0 ) < 1.0e-9 )`：条件分支：根据当前状态选择执行路径。
- **[行 2922]** `		face->set_boundary_id(0);`：函数调用语句，触发对应计算或操作。
- **[行 2923]** `	      else if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9)`：多分支条件判断，处理备选情形。
- **[行 2924]** `		face->set_boundary_id(1);`：函数调用语句，触发对应计算或操作。
- **[行 2925]** `	      else if (std::fabs(face->center()[2] - 0.0 ) < 1.0e-9)`：多分支条件判断，处理备选情形。
- **[行 2926]** `		face->set_boundary_id(2);`：函数调用语句，触发对应计算或操作。
- **[行 2927]** `	      else if (std::fabs(face->center()[2] - 1.0 ) < 1.0e-9)`：多分支条件判断，处理备选情形。
- **[行 2928]** `		face->set_boundary_id(3);`：函数调用语句，触发对应计算或操作。
- **[行 2929]** `	      else`：条件分支的兜底路径。
- **[行 2930]** `		face->set_boundary_id(4);`：函数调用语句，触发对应计算或操作。
- **[行 2931]** `	    }`：作用域边界（代码块开始/结束）。
- **[行 2932]** `	}`：作用域边界（代码块开始/结束）。
- **[行 2933]** ``：空行，用于分隔逻辑块。
- **[行 2934]** `    if (m_parameters.m_refinement_strategy == "adaptive-refine")`：条件分支：根据当前状态选择执行路径。
- **[行 2935]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2936]** `	unsigned int material_id;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2937]** `	double length_scale;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2938]** `	bool initiation_point_refine_unfinished = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2939]** `	while (initiation_point_refine_unfinished)`：循环推进（通常用于时间步或迭代收敛过程）。
- **[行 2940]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 2941]** `	    initiation_point_refine_unfinished = false;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2942]** `	    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2943]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 2944]** `		if (    cell->center()[2] > 0.505`：条件分支：根据当前状态选择执行路径。
- **[行 2945]** `		     && cell->center()[2] < 0.575`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2946]** `		     && cell->center()[0] < 0.05`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2947]** `		     && cell->center()[1] < 0.05 )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2948]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 2949]** `		    material_id = cell->material_id();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2950]** `		    length_scale = m_material_data[material_id][2];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2951]** `		    if (  std::cbrt(cell->measure())`：条件分支：根据当前状态选择执行路径。
- **[行 2952]** `			> length_scale * m_parameters.m_allowed_max_h_l_ratio )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2953]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 2954]** `			cell->set_refine_flag();`：函数调用语句，触发对应计算或操作。
- **[行 2955]** `			initiation_point_refine_unfinished = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 2956]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 2957]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 2958]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 2959]** `	    m_triangulation.execute_coarsening_and_refinement();`：函数调用语句，触发对应计算或操作。
- **[行 2960]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 2961]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2962]** `    else`：条件分支的兜底路径。
- **[行 2963]** `      {`：作用域边界（代码块开始/结束）。
- **[行 2964]** `	AssertThrow(false,`：运行期断言/检查，验证输入与状态合法性。
- **[行 2965]** `		    ExcMessage("Selected mesh refinement strategy not implemented!"));`：函数调用语句，触发对应计算或操作。
- **[行 2966]** `      }`：作用域边界（代码块开始/结束）。
- **[行 2967]** `  }`：作用域边界（代码块开始/结束）。
- **[行 2968]** ``：空行，用于分隔逻辑块。
- **[行 2969]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 2970]** `  void PhaseFieldMonolithicSolve<dim>::make_grid_case_9()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2971]** `  {`：作用域边界（代码块开始/结束）。
- **[行 2972]** `    AssertThrow(dim==2, ExcMessage("The dimension has to be 2D!"));`：运行期断言/检查，验证输入与状态合法性。
- **[行 2973]** ``：空行，用于分隔逻辑块。
- **[行 2974]** `    for (unsigned int i = 0; i < 80; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2975]** `      m_logfile << "*";`：写日志输出，记录当前计算状态与结果。
- **[行 2976]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2977]** `    m_logfile << "\t\t\t\tL-shape bending (2D structured)" << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2978]** `    for (unsigned int i = 0; i < 80; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2979]** `      m_logfile << "*";`：写日志输出，记录当前计算状态与结果。
- **[行 2980]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 2981]** ``：空行，用于分隔逻辑块。
- **[行 2982]** `    GridIn<dim> gridin;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 2983]** `    gridin.attach_triangulation(m_triangulation);`：函数调用语句，触发对应计算或操作。
- **[行 2984]** `    std::ifstream f("L-Shape.msh");`：调用 C++ 标准库工具函数/容器接口。
- **[行 2985]** `    gridin.read_msh(f);`：函数调用语句，触发对应计算或操作。
- **[行 2986]** ``：空行，用于分隔逻辑块。
- **[行 2987]** `    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2988]** `      for (const auto &face : cell->face_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 2989]** `	{`：作用域边界（代码块开始/结束）。
- **[行 2990]** `	  if (face->at_boundary() == true)`：条件分支：根据当前状态选择执行路径。
- **[行 2991]** `	    {`：作用域边界（代码块开始/结束）。
- **[行 2992]** `	      if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9 )`：条件分支：根据当前状态选择执行路径。
- **[行 2993]** `		face->set_boundary_id(0);`：函数调用语句，触发对应计算或操作。
- **[行 2994]** `	      else`：条件分支的兜底路径。
- **[行 2995]** `	        face->set_boundary_id(1);`：函数调用语句，触发对应计算或操作。
- **[行 2996]** `	    }`：作用域边界（代码块开始/结束）。
- **[行 2997]** `	}`：作用域边界（代码块开始/结束）。
- **[行 2998]** ``：空行，用于分隔逻辑块。
- **[行 2999]** `    m_triangulation.refine_global(m_parameters.m_global_refine_times);`：函数调用语句，触发对应计算或操作。
- **[行 3000]** ``：空行，用于分隔逻辑块。
- **[行 3001]** `    if (m_parameters.m_refinement_strategy == "pre-refine")`：条件分支：根据当前状态选择执行路径。
- **[行 3002]** `      {`：作用域边界（代码块开始/结束）。
- **[行 3003]** `	unsigned int material_id;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3004]** `	double length_scale;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3005]** `	for (unsigned int i = 0; i < m_parameters.m_local_prerefine_times; i++)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3006]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 3007]** `	    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3008]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 3009]** `		if (    (cell->center()[1] > 242.0)`：条件分支：根据当前状态选择执行路径。
- **[行 3010]** `		     && (cell->center()[1] < 312.5)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3011]** `		     && (cell->center()[0] < 258.0) )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3012]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 3013]** `		    material_id = cell->material_id();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3014]** `		    length_scale = m_material_data[material_id][2];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3015]** `		    if (  std::sqrt(cell->measure())`：条件分支：根据当前状态选择执行路径。
- **[行 3016]** `			> length_scale * m_parameters.m_allowed_max_h_l_ratio )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3017]** `		      cell->set_refine_flag();`：函数调用语句，触发对应计算或操作。
- **[行 3018]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 3019]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 3020]** `	    m_triangulation.execute_coarsening_and_refinement();`：函数调用语句，触发对应计算或操作。
- **[行 3021]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 3022]** `      }`：作用域边界（代码块开始/结束）。
- **[行 3023]** `    else if (m_parameters.m_refinement_strategy == "adaptive-refine")`：多分支条件判断，处理备选情形。
- **[行 3024]** `      {`：作用域边界（代码块开始/结束）。
- **[行 3025]** `	unsigned int material_id;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3026]** `	double length_scale;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3027]** `	bool initiation_point_refine_unfinished = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3028]** `	while (initiation_point_refine_unfinished)`：循环推进（通常用于时间步或迭代收敛过程）。
- **[行 3029]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 3030]** `	    initiation_point_refine_unfinished = false;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3031]** `	    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3032]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 3033]** `		if (             (cell->center()[0] - 250) < 0.0`：条件分支：根据当前状态选择执行路径。
- **[行 3034]** `		     &&          (cell->center()[0] - 240) > 0.0`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3035]** `		     && std::fabs(cell->center()[1] - 250) < 10.0 )`：调用 C++ 标准库工具函数/容器接口。
- **[行 3036]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 3037]** `		    material_id = cell->material_id();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3038]** `		    length_scale = m_material_data[material_id][2];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3039]** `		    if (  std::sqrt(cell->measure())`：条件分支：根据当前状态选择执行路径。
- **[行 3040]** `			> length_scale * m_parameters.m_allowed_max_h_l_ratio )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3041]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 3042]** `		        cell->set_refine_flag();`：函数调用语句，触发对应计算或操作。
- **[行 3043]** `		        initiation_point_refine_unfinished = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3044]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 3045]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 3046]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 3047]** `	    m_triangulation.execute_coarsening_and_refinement();`：函数调用语句，触发对应计算或操作。
- **[行 3048]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 3049]** `      }`：作用域边界（代码块开始/结束）。
- **[行 3050]** `    else`：条件分支的兜底路径。
- **[行 3051]** `      {`：作用域边界（代码块开始/结束）。
- **[行 3052]** `	AssertThrow(false,`：运行期断言/检查，验证输入与状态合法性。
- **[行 3053]** `	            ExcMessage("Selected mesh refinement strategy not implemented!"));`：函数调用语句，触发对应计算或操作。
- **[行 3054]** `      }`：作用域边界（代码块开始/结束）。
- **[行 3055]** `  }`：作用域边界（代码块开始/结束）。
- **[行 3056]** ``：空行，用于分隔逻辑块。
- **[行 3057]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 3058]** `  void PhaseFieldMonolithicSolve<dim>::make_grid_case_10()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3059]** `  {`：作用域边界（代码块开始/结束）。
- **[行 3060]** `    AssertThrow(dim==3, ExcMessage("The dimension has to be 3D!"));`：运行期断言/检查，验证输入与状态合法性。
- **[行 3061]** ``：空行，用于分隔逻辑块。
- **[行 3062]** `    for (unsigned int i = 0; i < 80; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3063]** `      m_logfile << "*";`：写日志输出，记录当前计算状态与结果。
- **[行 3064]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 3065]** `    m_logfile << "\t\t\t\tL-shape bending (3D structured)" << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 3066]** `    for (unsigned int i = 0; i < 80; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3067]** `      m_logfile << "*";`：写日志输出，记录当前计算状态与结果。
- **[行 3068]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 3069]** ``：空行，用于分隔逻辑块。
- **[行 3070]** `    Triangulation<2> triangulation_2d;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3071]** ``：空行，用于分隔逻辑块。
- **[行 3072]** `    GridIn<2> gridin;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3073]** `    gridin.attach_triangulation(triangulation_2d);`：函数调用语句，触发对应计算或操作。
- **[行 3074]** `    std::ifstream f("L-Shape.msh");`：调用 C++ 标准库工具函数/容器接口。
- **[行 3075]** `    gridin.read_msh(f);`：函数调用语句，触发对应计算或操作。
- **[行 3076]** ``：空行，用于分隔逻辑块。
- **[行 3077]** `    const double thickness = 150.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3078]** `    const unsigned int n_layer = 11;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3079]** `    GridGenerator::extrude_triangulation(triangulation_2d, n_layer, thickness, m_triangulation);`：函数调用语句，触发对应计算或操作。
- **[行 3080]** ``：空行，用于分隔逻辑块。
- **[行 3081]** `    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3082]** `      for (const auto &face : cell->face_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3083]** `	{`：作用域边界（代码块开始/结束）。
- **[行 3084]** `	  if (face->at_boundary() == true)`：条件分支：根据当前状态选择执行路径。
- **[行 3085]** `	    {`：作用域边界（代码块开始/结束）。
- **[行 3086]** `	      if (std::fabs(face->center()[1] - 0.0 ) < 1.0e-9 )`：条件分支：根据当前状态选择执行路径。
- **[行 3087]** `		face->set_boundary_id(0);`：函数调用语句，触发对应计算或操作。
- **[行 3088]** `	      else`：条件分支的兜底路径。
- **[行 3089]** `	        face->set_boundary_id(1);`：函数调用语句，触发对应计算或操作。
- **[行 3090]** `	    }`：作用域边界（代码块开始/结束）。
- **[行 3091]** `	}`：作用域边界（代码块开始/结束）。
- **[行 3092]** ``：空行，用于分隔逻辑块。
- **[行 3093]** `    m_triangulation.refine_global(m_parameters.m_global_refine_times);`：函数调用语句，触发对应计算或操作。
- **[行 3094]** ``：空行，用于分隔逻辑块。
- **[行 3095]** `    if (m_parameters.m_refinement_strategy == "pre-refine")`：条件分支：根据当前状态选择执行路径。
- **[行 3096]** `      {`：作用域边界（代码块开始/结束）。
- **[行 3097]** `	unsigned int material_id;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3098]** `	double length_scale;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3099]** `	for (unsigned int i = 0; i < m_parameters.m_local_prerefine_times; i++)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3100]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 3101]** `	    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3102]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 3103]** `		if (    (std::fabs(cell->center()[1] - 250.0) < 10.0)`：条件分支：根据当前状态选择执行路径。
- **[行 3104]** `		     && (cell->center()[0] < 250.0) )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3105]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 3106]** `		    material_id = cell->material_id();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3107]** `		    length_scale = m_material_data[material_id][2];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3108]** `		    if (  std::cbrt(cell->measure())`：条件分支：根据当前状态选择执行路径。
- **[行 3109]** `			> length_scale * m_parameters.m_allowed_max_h_l_ratio )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3110]** `		      cell->set_refine_flag();`：函数调用语句，触发对应计算或操作。
- **[行 3111]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 3112]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 3113]** `	    m_triangulation.execute_coarsening_and_refinement();`：函数调用语句，触发对应计算或操作。
- **[行 3114]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 3115]** `      }`：作用域边界（代码块开始/结束）。
- **[行 3116]** `    else if (m_parameters.m_refinement_strategy == "adaptive-refine")`：多分支条件判断，处理备选情形。
- **[行 3117]** `      {`：作用域边界（代码块开始/结束）。
- **[行 3118]** `	unsigned int material_id;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3119]** `	double length_scale;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3120]** `	bool initiation_point_refine_unfinished = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3121]** `	while (initiation_point_refine_unfinished)`：循环推进（通常用于时间步或迭代收敛过程）。
- **[行 3122]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 3123]** `	    initiation_point_refine_unfinished = false;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3124]** `	    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3125]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 3126]** `		if (             (cell->center()[0] - 250) < 0.0`：条件分支：根据当前状态选择执行路径。
- **[行 3127]** `		     &&          (cell->center()[0] - 240) > 0.0`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3128]** `		     && std::fabs(cell->center()[1] - 250) < 10.0 )`：调用 C++ 标准库工具函数/容器接口。
- **[行 3129]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 3130]** `		    material_id = cell->material_id();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3131]** `		    length_scale = m_material_data[material_id][2];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3132]** `		    if (  std::cbrt(cell->measure())`：条件分支：根据当前状态选择执行路径。
- **[行 3133]** `			> length_scale * m_parameters.m_allowed_max_h_l_ratio )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3134]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 3135]** `		        cell->set_refine_flag();`：函数调用语句，触发对应计算或操作。
- **[行 3136]** `		        initiation_point_refine_unfinished = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3137]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 3138]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 3139]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 3140]** `	    m_triangulation.execute_coarsening_and_refinement();`：函数调用语句，触发对应计算或操作。
- **[行 3141]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 3142]** `      }`：作用域边界（代码块开始/结束）。
- **[行 3143]** `    else`：条件分支的兜底路径。
- **[行 3144]** `      {`：作用域边界（代码块开始/结束）。
- **[行 3145]** `	AssertThrow(false,`：运行期断言/检查，验证输入与状态合法性。
- **[行 3146]** `	            ExcMessage("Selected mesh refinement strategy not implemented!"));`：函数调用语句，触发对应计算或操作。
- **[行 3147]** `      }`：作用域边界（代码块开始/结束）。
- **[行 3148]** `  }`：作用域边界（代码块开始/结束）。
- **[行 3149]** ``：空行，用于分隔逻辑块。
- **[行 3150]** ``：空行，用于分隔逻辑块。
- **[行 3151]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 3152]** `  void PhaseFieldMonolithicSolve<dim>::make_grid_case_11()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3153]** `  {`：作用域边界（代码块开始/结束）。
- **[行 3154]** `    AssertThrow(dim==3, ExcMessage("The dimension has to be 3D!"));`：运行期断言/检查，验证输入与状态合法性。
- **[行 3155]** ``：空行，用于分隔逻辑块。
- **[行 3156]** `    for (unsigned int i = 0; i < 80; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3157]** `      m_logfile << "*";`：写日志输出，记录当前计算状态与结果。
- **[行 3158]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 3159]** `    m_logfile << "\t\t\t\tBrokenshire torsion (3D structured)" << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 3160]** `    for (unsigned int i = 0; i < 80; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3161]** `      m_logfile << "*";`：写日志输出，记录当前计算状态与结果。
- **[行 3162]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 3163]** ``：空行，用于分隔逻辑块。
- **[行 3164]** `    Triangulation<2> triangulation_2d;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3165]** ``：空行，用于分隔逻辑块。
- **[行 3166]** `    double const length = 200.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3167]** `    double const width = 50.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3168]** `    double const height = 50.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3169]** `    double const delta_L = 25.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3170]** `    double const tan_theta = delta_L / (0.5*width);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3171]** ``：空行，用于分隔逻辑块。
- **[行 3172]** `    std::vector<unsigned int> repetitions(2, 1);`：调用 C++ 标准库工具函数/容器接口。
- **[行 3173]** `    repetitions[0] = 20;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3174]** `    repetitions[1] = 5;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3175]** ``：空行，用于分隔逻辑块。
- **[行 3176]** `    Point<2> point1(0.0, 0.0);`：函数调用语句，触发对应计算或操作。
- **[行 3177]** `    Point<2> point2(length, width);`：函数调用语句，触发对应计算或操作。
- **[行 3178]** ``：空行，用于分隔逻辑块。
- **[行 3179]** `    GridGenerator::subdivided_hyper_rectangle(triangulation_2d,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3180]** `					      repetitions,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3181]** `					      point1,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3182]** `					      point2 );`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3183]** ``：空行，用于分隔逻辑块。
- **[行 3184]** `    typename Triangulation<2>::vertex_iterator vertex_ptr;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3185]** `    vertex_ptr = triangulation_2d.begin_active_vertex();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3186]** `    while (vertex_ptr != triangulation_2d.end_vertex())`：循环推进（通常用于时间步或迭代收敛过程）。
- **[行 3187]** `      {`：作用域边界（代码块开始/结束）。
- **[行 3188]** `	Point<2> & vertex_point = vertex_ptr->vertex();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3189]** ``：空行，用于分隔逻辑块。
- **[行 3190]** `	const double delta_x = (vertex_point(1) - 0.5*width) * tan_theta;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3191]** ``：空行，用于分隔逻辑块。
- **[行 3192]** `	if (std::fabs(vertex_point(0) - 0.5*length) < 1.0e-6)`：条件分支：根据当前状态选择执行路径。
- **[行 3193]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 3194]** `	    vertex_point(0) += delta_x;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3195]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 3196]** `	else if (std::fabs(vertex_point(0) + length/repetitions[0] - 0.5*length) < 1.0e-6)`：多分支条件判断，处理备选情形。
- **[行 3197]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 3198]** `	    vertex_point(0) += (delta_x + length/repetitions[0]*0.5);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3199]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 3200]** `	else if (std::fabs(vertex_point(0) - length/repetitions[0] - 0.5*length) < 1.0e-6)`：多分支条件判断，处理备选情形。
- **[行 3201]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 3202]** `	    vertex_point(0) += (delta_x - length/repetitions[0]*0.5);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3203]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 3204]** `	else if (vertex_point(0) < 0.5*length - length/repetitions[0] - 1.0e-6)`：多分支条件判断，处理备选情形。
- **[行 3205]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 3206]** `	    vertex_point(0) += (delta_x + length/repetitions[0]*0.5) * vertex_point(0)/(0.5*length - length/repetitions[0]);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3207]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 3208]** `	else if (vertex_point(0) > 0.5*length + length/repetitions[0] + 1.0e-6)`：多分支条件判断，处理备选情形。
- **[行 3209]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 3210]** `	    vertex_point(0) += (delta_x - length/repetitions[0]*0.5) * (length - vertex_point(0))/(0.5*length - length/repetitions[0]);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3211]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 3212]** ``：空行，用于分隔逻辑块。
- **[行 3213]** `	++vertex_ptr;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3214]** `      }`：作用域边界（代码块开始/结束）。
- **[行 3215]** ``：空行，用于分隔逻辑块。
- **[行 3216]** `    Triangulation<dim> tmp_triangulation;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3217]** `    const unsigned int n_layer = repetitions[1] + 1;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3218]** `    GridGenerator::extrude_triangulation(triangulation_2d, n_layer, height, tmp_triangulation);`：函数调用语句，触发对应计算或操作。
- **[行 3219]** ``：空行，用于分隔逻辑块。
- **[行 3220]** `    tmp_triangulation.refine_global(m_parameters.m_global_refine_times);`：函数调用语句，触发对应计算或操作。
- **[行 3221]** ``：空行，用于分隔逻辑块。
- **[行 3222]** `    std::set<typename Triangulation< dim >::active_cell_iterator >`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3223]** `      cells_to_remove;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3224]** ``：空行，用于分隔逻辑块。
- **[行 3225]** `    for (const auto &cell : tmp_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3226]** `      {`：作用域边界（代码块开始/结束）。
- **[行 3227]** `	if (    (std::fabs(cell->center()[0] - (cell->center()[1] - 0.5*width)*tan_theta - 0.5*length) < 2.5)`：条件分支：根据当前状态选择执行路径。
- **[行 3228]** `	     && cell->center()[2] > 0.5* height  )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3229]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 3230]** `	    cells_to_remove.insert(cell);`：函数调用语句，触发对应计算或操作。
- **[行 3231]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 3232]** `      }`：作用域边界（代码块开始/结束）。
- **[行 3233]** ``：空行，用于分隔逻辑块。
- **[行 3234]** `    GridGenerator::create_triangulation_with_removed_cells(tmp_triangulation,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3235]** `							   cells_to_remove,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3236]** `							   m_triangulation);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3237]** ``：空行，用于分隔逻辑块。
- **[行 3238]** `    if (m_parameters.m_refinement_strategy == "adaptive-refine")`：条件分支：根据当前状态选择执行路径。
- **[行 3239]** `      {`：作用域边界（代码块开始/结束）。
- **[行 3240]** `	unsigned int material_id;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3241]** `	double length_scale;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3242]** `	bool initiation_point_refine_unfinished = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3243]** `	while (initiation_point_refine_unfinished)`：循环推进（通常用于时间步或迭代收敛过程）。
- **[行 3244]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 3245]** `	    initiation_point_refine_unfinished = false;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3246]** `	    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3247]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 3248]** `		if (  (std::fabs(cell->center()[0] - (cell->center()[1] - 0.5*width)*tan_theta - 0.5*length) < 5.0)`：条件分支：根据当前状态选择执行路径。
- **[行 3249]** `		    && cell->center()[2] <= 0.5*height`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3250]** `		    && cell->center()[2] > 0.5*height - 5.0 )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3251]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 3252]** `		    material_id = cell->material_id();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3253]** `		    length_scale = m_material_data[material_id][2];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3254]** `		    if (  std::cbrt(cell->measure())`：条件分支：根据当前状态选择执行路径。
- **[行 3255]** `			> length_scale * m_parameters.m_allowed_max_h_l_ratio )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3256]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 3257]** `			cell->set_refine_flag();`：函数调用语句，触发对应计算或操作。
- **[行 3258]** `			initiation_point_refine_unfinished = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3259]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 3260]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 3261]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 3262]** `	    m_triangulation.execute_coarsening_and_refinement();`：函数调用语句，触发对应计算或操作。
- **[行 3263]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 3264]** `      }`：作用域边界（代码块开始/结束）。
- **[行 3265]** `    else`：条件分支的兜底路径。
- **[行 3266]** `      {`：作用域边界（代码块开始/结束）。
- **[行 3267]** `	AssertThrow(false,`：运行期断言/检查，验证输入与状态合法性。
- **[行 3268]** `		    ExcMessage("Selected mesh refinement strategy not implemented!"));`：函数调用语句，触发对应计算或操作。
- **[行 3269]** `      }`：作用域边界（代码块开始/结束）。
- **[行 3270]** ``：空行，用于分隔逻辑块。
- **[行 3271]** ``：空行，用于分隔逻辑块。
- **[行 3272]** `    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3273]** `      for (const auto &face : cell->face_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3274]** `	{`：作用域边界（代码块开始/结束）。
- **[行 3275]** `	  if (face->at_boundary() == true)`：条件分支：根据当前状态选择执行路径。
- **[行 3276]** `	    {`：作用域边界（代码块开始/结束）。
- **[行 3277]** `	      if (std::fabs(face->center()[0] - length) < 1.0e-6 )`：条件分支：根据当前状态选择执行路径。
- **[行 3278]** `		face->set_boundary_id(0);`：函数调用语句，触发对应计算或操作。
- **[行 3279]** `	      else if (std::fabs(face->center()[0] - 0.0) < 1.0e-6 )`：多分支条件判断，处理备选情形。
- **[行 3280]** `		face->set_boundary_id(1);`：函数调用语句，触发对应计算或操作。
- **[行 3281]** `	      else`：条件分支的兜底路径。
- **[行 3282]** `		face->set_boundary_id(2);`：函数调用语句，触发对应计算或操作。
- **[行 3283]** `	    }`：作用域边界（代码块开始/结束）。
- **[行 3284]** `	}`：作用域边界（代码块开始/结束）。
- **[行 3285]** `  }`：作用域边界（代码块开始/结束）。
- **[行 3286]** ``：空行，用于分隔逻辑块。
- **[行 3287]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 3288]** `  void PhaseFieldMonolithicSolve<dim>::setup_system()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3289]** `  {`：作用域边界（代码块开始/结束）。
- **[行 3290]** `    m_timer.enter_subsection("Setup system");`：函数调用语句，触发对应计算或操作。
- **[行 3291]** ``：空行，用于分隔逻辑块。
- **[行 3292]** `    std::vector<unsigned int> block_component(m_n_components,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3293]** `                                              m_u_dof); // displacement`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3294]** `    block_component[m_d_component] = m_d_dof;           // phasefield`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3295]** ``：空行，用于分隔逻辑块。
- **[行 3296]** `    m_dof_handler.distribute_dofs(m_fe);`：函数调用语句，触发对应计算或操作。
- **[行 3297]** `    DoFRenumbering::Cuthill_McKee(m_dof_handler);`：函数调用语句，触发对应计算或操作。
- **[行 3298]** `    DoFRenumbering::component_wise(m_dof_handler, block_component);`：函数调用语句，触发对应计算或操作。
- **[行 3299]** ``：空行，用于分隔逻辑块。
- **[行 3300]** `    m_constraints.clear();`：函数调用语句，触发对应计算或操作。
- **[行 3301]** `    DoFTools::make_hanging_node_constraints(m_dof_handler, m_constraints);`：函数调用语句，触发对应计算或操作。
- **[行 3302]** `    m_constraints.close();`：函数调用语句，触发对应计算或操作。
- **[行 3303]** ``：空行，用于分隔逻辑块。
- **[行 3304]** `    m_dofs_per_block =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3305]** `      DoFTools::count_dofs_per_fe_block(m_dof_handler, block_component);`：函数调用语句，触发对应计算或操作。
- **[行 3306]** ``：空行，用于分隔逻辑块。
- **[行 3307]** `    m_logfile << "\t\tTriangulation:"`：写日志输出，记录当前计算状态与结果。
- **[行 3308]** `              << "\n\t\t\t Number of active cells: "`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3309]** `              << m_triangulation.n_active_cells()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3310]** `              << "\n\t\t\t Number of used vertices: "`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3311]** `              << m_triangulation.n_used_vertices()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3312]** `              << "\n\t\t\t Number of active edges: "`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3313]** `              << m_triangulation.n_active_lines()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3314]** `              << "\n\t\t\t Number of active faces: "`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3315]** `              << m_triangulation.n_active_faces()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3316]** `              << "\n\t\t\t Number of degrees of freedom (total): "`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3317]** `	      << m_dof_handler.n_dofs()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3318]** `	      << "\n\t\t\t Number of degrees of freedom (disp): "`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3319]** `	      << m_dofs_per_block[m_u_dof]`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3320]** `	      << "\n\t\t\t Number of degrees of freedom (phasefield): "`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3321]** `	      << m_dofs_per_block[m_d_dof]`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3322]** `              << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3323]** ``：空行，用于分隔逻辑块。
- **[行 3324]** `    m_tangent_matrix.clear();`：函数调用语句，触发对应计算或操作。
- **[行 3325]** `    {`：作用域边界（代码块开始/结束）。
- **[行 3326]** `      BlockDynamicSparsityPattern dsp(m_dofs_per_block, m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 3327]** ``：空行，用于分隔逻辑块。
- **[行 3328]** `      Table<2, DoFTools::Coupling> coupling(m_n_components, m_n_components);`：函数调用语句，触发对应计算或操作。
- **[行 3329]** `      for (unsigned int ii = 0; ii < m_n_components; ++ii)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3330]** `        for (unsigned int jj = 0; jj < m_n_components; ++jj)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3331]** `          {`：作用域边界（代码块开始/结束）。
- **[行 3332]** `            if (   ((ii < m_d_component) && (jj == m_d_component))`：条件分支：根据当前状态选择执行路径。
- **[行 3333]** `                || ((ii == m_d_component) && (jj < m_d_component)) )`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3334]** `              coupling[ii][jj] = DoFTools::none;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3335]** `            else`：条件分支的兜底路径。
- **[行 3336]** `              coupling[ii][jj] = DoFTools::always;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3337]** `          }`：作用域边界（代码块开始/结束）。
- **[行 3338]** ``：空行，用于分隔逻辑块。
- **[行 3339]** `      DoFTools::make_sparsity_pattern(`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3340]** `        m_dof_handler, coupling, dsp, m_constraints, false);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3341]** `      m_sparsity_pattern.copy_from(dsp);`：函数调用语句，触发对应计算或操作。
- **[行 3342]** `    }`：作用域边界（代码块开始/结束）。
- **[行 3343]** ``：空行，用于分隔逻辑块。
- **[行 3344]** `    m_tangent_matrix.reinit(m_sparsity_pattern);`：函数调用语句，触发对应计算或操作。
- **[行 3345]** ``：空行，用于分隔逻辑块。
- **[行 3346]** `    m_system_rhs.reinit(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 3347]** `    m_solution.reinit(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 3348]** ``：空行，用于分隔逻辑块。
- **[行 3349]** `    m_active_set_phasefield.reinit(m_dofs_per_block[m_d_dof]);`：函数调用语句，触发对应计算或操作。
- **[行 3350]** ``：空行，用于分隔逻辑块。
- **[行 3351]** `    setup_qph();`：函数调用语句，触发对应计算或操作。
- **[行 3352]** ``：空行，用于分隔逻辑块。
- **[行 3353]** `    m_timer.leave_subsection();`：函数调用语句，触发对应计算或操作。
- **[行 3354]** `  }`：作用域边界（代码块开始/结束）。
- **[行 3355]** ``：空行，用于分隔逻辑块。
- **[行 3356]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 3357]** `  void PhaseFieldMonolithicSolve<dim>::make_constraints(const unsigned int it_nr)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3358]** `  {`：作用域边界（代码块开始/结束）。
- **[行 3359]** `    const bool apply_dirichlet_bc = (it_nr == 0);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3360]** ``：空行，用于分隔逻辑块。
- **[行 3361]** `    if (it_nr > 1)`：条件分支：根据当前状态选择执行路径。
- **[行 3362]** `      {`：作用域边界（代码块开始/结束）。
- **[行 3363]** `        return;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3364]** `      }`：作用域边界（代码块开始/结束）。
- **[行 3365]** ``：空行，用于分隔逻辑块。
- **[行 3366]** `    if (apply_dirichlet_bc)`：条件分支：根据当前状态选择执行路径。
- **[行 3367]** `      {`：作用域边界（代码块开始/结束）。
- **[行 3368]** `	m_constraints.clear();`：函数调用语句，触发对应计算或操作。
- **[行 3369]** `	DoFTools::make_hanging_node_constraints(m_dof_handler,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3370]** `						m_constraints);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3371]** ``：空行，用于分隔逻辑块。
- **[行 3372]** `	const FEValuesExtractors::Scalar x_displacement(0);`：函数调用语句，触发对应计算或操作。
- **[行 3373]** `	const FEValuesExtractors::Scalar y_displacement(1);`：函数调用语句，触发对应计算或操作。
- **[行 3374]** `	const FEValuesExtractors::Scalar z_displacement(2);`：函数调用语句，触发对应计算或操作。
- **[行 3375]** ``：空行，用于分隔逻辑块。
- **[行 3376]** `	const FEValuesExtractors::Vector displacements(0);`：函数调用语句，触发对应计算或操作。
- **[行 3377]** ``：空行，用于分隔逻辑块。
- **[行 3378]** `	if (   m_parameters.m_scenario == 1`：条件分支：根据当前状态选择执行路径。
- **[行 3379]** `	    || m_parameters.m_scenario == 3)`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3380]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 3381]** `	    // Dirichlet B,C. bottom surface`：注释行，用于说明算法背景、假设或实现意图。
- **[行 3382]** `	    const int boundary_id_bottom_surface = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3383]** `	    VectorTools::interpolate_boundary_values(m_dof_handler,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3384]** `						     boundary_id_bottom_surface,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3385]** `						     Functions::ZeroFunction<dim>(m_n_components),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3386]** `						     m_constraints,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3387]** `						     m_fe.component_mask(y_displacement));`：函数调用语句，触发对应计算或操作。
- **[行 3388]** ``：空行，用于分隔逻辑块。
- **[行 3389]** `	    typename Triangulation<dim>::active_vertex_iterator vertex_itr;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3390]** `	    vertex_itr = m_triangulation.begin_active_vertex();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3391]** `	    std::vector<types::global_dof_index> node_xy(m_fe.dofs_per_vertex);`：调用 C++ 标准库工具函数/容器接口。
- **[行 3392]** ``：空行，用于分隔逻辑块。
- **[行 3393]** `	    for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3394]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 3395]** `		if (   (std::fabs(vertex_itr->vertex()[0] - 0.0) < 1.0e-9)`：条件分支：根据当前状态选择执行路径。
- **[行 3396]** `		    && (std::fabs(vertex_itr->vertex()[1] - 0.0) < 1.0e-9) )`：调用 C++ 标准库工具函数/容器接口。
- **[行 3397]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 3398]** `		    node_xy = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3399]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 3400]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 3401]** `	    m_constraints.add_line(node_xy[0]);`：函数调用语句，触发对应计算或操作。
- **[行 3402]** `	    m_constraints.set_inhomogeneity(node_xy[0], 0.0);`：函数调用语句，触发对应计算或操作。
- **[行 3403]** ``：空行，用于分隔逻辑块。
- **[行 3404]** `	    m_constraints.add_line(node_xy[1]);`：函数调用语句，触发对应计算或操作。
- **[行 3405]** `	    m_constraints.set_inhomogeneity(node_xy[1], 0.0);`：函数调用语句，触发对应计算或操作。
- **[行 3406]** ``：空行，用于分隔逻辑块。
- **[行 3407]** `	    const int boundary_id_top_surface = 1;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3408]** `	    /*`：注释行，用于说明算法背景、假设或实现意图。
- **[行 3409]** `	    VectorTools::interpolate_boundary_values(m_dof_handler,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3410]** `						     boundary_id_top_surface,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3411]** `						     Functions::ZeroFunction<dim>(m_n_components),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3412]** `						     m_constraints,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3413]** `						     m_fe.component_mask(x_displacement));`：函数调用语句，触发对应计算或操作。
- **[行 3414]** `	    */`：注释行，用于说明算法背景、假设或实现意图。
- **[行 3415]** `            const double time_inc = m_time.get_delta_t();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3416]** `            double disp_magnitude = m_time.get_magnitude();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3417]** `	    VectorTools::interpolate_boundary_values(m_dof_handler,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3418]** `						     boundary_id_top_surface,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3419]** `						     Functions::ConstantFunction<dim>(`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3420]** `						       disp_magnitude*time_inc, m_n_components),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3421]** `						     m_constraints,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3422]** `						     m_fe.component_mask(y_displacement));`：函数调用语句，触发对应计算或操作。
- **[行 3423]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 3424]** `	else if (   m_parameters.m_scenario == 2`：多分支条件判断，处理备选情形。
- **[行 3425]** `	         || m_parameters.m_scenario == 4)`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3426]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 3427]** `	    // Dirichlet B,C. bottom surface`：注释行，用于说明算法背景、假设或实现意图。
- **[行 3428]** `	    const int boundary_id_bottom_surface = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3429]** `	    VectorTools::interpolate_boundary_values(m_dof_handler,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3430]** `						     boundary_id_bottom_surface,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3431]** `						     Functions::ZeroFunction<dim>(m_n_components),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3432]** `						     m_constraints,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3433]** `						     m_fe.component_mask(displacements));`：函数调用语句，触发对应计算或操作。
- **[行 3434]** ``：空行，用于分隔逻辑块。
- **[行 3435]** `	    const int boundary_id_top_surface = 1;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3436]** `	    VectorTools::interpolate_boundary_values(m_dof_handler,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3437]** `						     boundary_id_top_surface,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3438]** `						     Functions::ZeroFunction<dim>(m_n_components),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3439]** `						     m_constraints,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3440]** `						     m_fe.component_mask(y_displacement));`：函数调用语句，触发对应计算或操作。
- **[行 3441]** ``：空行，用于分隔逻辑块。
- **[行 3442]** `	    const double time_inc = m_time.get_delta_t();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3443]** `	    double disp_magnitude = m_time.get_magnitude();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3444]** `	    VectorTools::interpolate_boundary_values(m_dof_handler,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3445]** `						     boundary_id_top_surface,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3446]** `						     Functions::ConstantFunction<dim>(`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3447]** `						       disp_magnitude*time_inc, m_n_components),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3448]** `						     m_constraints,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3449]** `						     m_fe.component_mask(x_displacement));`：函数调用语句，触发对应计算或操作。
- **[行 3450]** ``：空行，用于分隔逻辑块。
- **[行 3451]** `	    const int boundary_id_side_surfaces = 2;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3452]** `	    VectorTools::interpolate_boundary_values(m_dof_handler,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3453]** `						     boundary_id_side_surfaces,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3454]** `						     Functions::ZeroFunction<dim>(m_n_components),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3455]** `						     m_constraints,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3456]** `						     m_fe.component_mask(y_displacement));`：函数调用语句，触发对应计算或操作。
- **[行 3457]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 3458]** `	else if (m_parameters.m_scenario == 5)`：多分支条件判断，处理备选情形。
- **[行 3459]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 3460]** `	    typename Triangulation<dim>::active_vertex_iterator vertex_itr;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3461]** `	    vertex_itr = m_triangulation.begin_active_vertex();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3462]** `	    std::vector<types::global_dof_index> node_bottomleft(m_fe.dofs_per_vertex);`：调用 C++ 标准库工具函数/容器接口。
- **[行 3463]** `	    std::vector<types::global_dof_index> node_bottomright(m_fe.dofs_per_vertex);`：调用 C++ 标准库工具函数/容器接口。
- **[行 3464]** `	    std::vector<types::global_dof_index> node_topcenter(m_fe.dofs_per_vertex);`：调用 C++ 标准库工具函数/容器接口。
- **[行 3465]** ``：空行，用于分隔逻辑块。
- **[行 3466]** `	    for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3467]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 3468]** `		if (   (std::fabs(vertex_itr->vertex()[0] - 0.0) < 1.0e-9)`：条件分支：根据当前状态选择执行路径。
- **[行 3469]** `		    && (std::fabs(vertex_itr->vertex()[1] - 0.0) < 1.0e-9) )`：调用 C++ 标准库工具函数/容器接口。
- **[行 3470]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 3471]** `		    node_bottomleft = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3472]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 3473]** `		if (   (std::fabs(vertex_itr->vertex()[0] - 8.0) < 1.0e-9)`：条件分支：根据当前状态选择执行路径。
- **[行 3474]** `		    && (std::fabs(vertex_itr->vertex()[1] - 0.0) < 1.0e-9) )`：调用 C++ 标准库工具函数/容器接口。
- **[行 3475]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 3476]** `		    node_bottomright = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3477]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 3478]** `		if (   (std::fabs(vertex_itr->vertex()[0] - 4.0) < 1.0e-9)`：条件分支：根据当前状态选择执行路径。
- **[行 3479]** `		    && (std::fabs(vertex_itr->vertex()[1] - 2.0) < 1.0e-9) )`：调用 C++ 标准库工具函数/容器接口。
- **[行 3480]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 3481]** `		    node_topcenter = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3482]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 3483]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 3484]** `	    // bottom-left node fixed in both x- and y-directions`：注释行，用于说明算法背景、假设或实现意图。
- **[行 3485]** `	    m_constraints.add_line(node_bottomleft[0]);`：函数调用语句，触发对应计算或操作。
- **[行 3486]** `	    m_constraints.set_inhomogeneity(node_bottomleft[0], 0.0);`：函数调用语句，触发对应计算或操作。
- **[行 3487]** ``：空行，用于分隔逻辑块。
- **[行 3488]** `	    m_constraints.add_line(node_bottomleft[1]);`：函数调用语句，触发对应计算或操作。
- **[行 3489]** `	    m_constraints.set_inhomogeneity(node_bottomleft[1], 0.0);`：函数调用语句，触发对应计算或操作。
- **[行 3490]** ``：空行，用于分隔逻辑块。
- **[行 3491]** `	    // bottom-right node only fixed in y-direction`：注释行，用于说明算法背景、假设或实现意图。
- **[行 3492]** `	    m_constraints.add_line(node_bottomright[1]);`：函数调用语句，触发对应计算或操作。
- **[行 3493]** `	    m_constraints.set_inhomogeneity(node_bottomright[1], 0.0);`：函数调用语句，触发对应计算或操作。
- **[行 3494]** ``：空行，用于分隔逻辑块。
- **[行 3495]** `	    // top-center node applied with y-displacement`：注释行，用于说明算法背景、假设或实现意图。
- **[行 3496]** `	    const double time_inc = m_time.get_delta_t();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3497]** `	    double disp_magnitude = m_time.get_magnitude();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3498]** ``：空行，用于分隔逻辑块。
- **[行 3499]** `	    m_constraints.add_line(node_topcenter[1]);`：函数调用语句，触发对应计算或操作。
- **[行 3500]** `	    m_constraints.set_inhomogeneity(node_topcenter[1], disp_magnitude*time_inc);`：函数调用语句，触发对应计算或操作。
- **[行 3501]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 3502]** `	else if (   m_parameters.m_scenario == 6`：多分支条件判断，处理备选情形。
- **[行 3503]** `	         || m_parameters.m_scenario == 7`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3504]** `		 || m_parameters.m_scenario == 8)`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3505]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 3506]** `	    const int x0_surface = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3507]** `	    VectorTools::interpolate_boundary_values(m_dof_handler,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3508]** `						     x0_surface,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3509]** `						     Functions::ZeroFunction<dim>(m_n_components),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3510]** `						     m_constraints,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3511]** `						     m_fe.component_mask(x_displacement));`：函数调用语句，触发对应计算或操作。
- **[行 3512]** `	    const int y0_surface = 1;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3513]** `	    VectorTools::interpolate_boundary_values(m_dof_handler,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3514]** `						     y0_surface,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3515]** `						     Functions::ZeroFunction<dim>(m_n_components),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3516]** `						     m_constraints,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3517]** `						     m_fe.component_mask(y_displacement));`：函数调用语句，触发对应计算或操作。
- **[行 3518]** `	    const int z0_surface = 2;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3519]** `	    VectorTools::interpolate_boundary_values(m_dof_handler,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3520]** `						     z0_surface,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3521]** `						     Functions::ZeroFunction<dim>(m_n_components),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3522]** `						     m_constraints,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3523]** `						     m_fe.component_mask(z_displacement));`：函数调用语句，触发对应计算或操作。
- **[行 3524]** ``：空行，用于分隔逻辑块。
- **[行 3525]** `	    const int z1_surface = 3;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3526]** `	    const double time_inc = m_time.get_delta_t();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3527]** `	    double disp_magnitude = 1.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3528]** `	    VectorTools::interpolate_boundary_values(m_dof_handler,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3529]** `						     z1_surface,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3530]** `						     Functions::ConstantFunction<dim>(`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3531]** `						       disp_magnitude*time_inc, m_n_components),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3532]** `						     m_constraints,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3533]** `						     m_fe.component_mask(z_displacement));`：函数调用语句，触发对应计算或操作。
- **[行 3534]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 3535]** `	else if (   m_parameters.m_scenario == 9`：多分支条件判断，处理备选情形。
- **[行 3536]** `	         || m_parameters.m_scenario == 10)`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3537]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 3538]** `	    // Dirichlet B,C. bottom surface`：注释行，用于说明算法背景、假设或实现意图。
- **[行 3539]** `	    const int boundary_id_bottom_surface = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3540]** `	    VectorTools::interpolate_boundary_values(m_dof_handler,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3541]** `						     boundary_id_bottom_surface,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3542]** `						     Functions::ZeroFunction<dim>(m_n_components),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3543]** `						     m_constraints,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3544]** `						     m_fe.component_mask(displacements));`：函数调用语句，触发对应计算或操作。
- **[行 3545]** ``：空行，用于分隔逻辑块。
- **[行 3546]** `	    typename Triangulation<dim>::active_vertex_iterator vertex_itr;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3547]** `	    vertex_itr = m_triangulation.begin_active_vertex();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3548]** `	    std::vector<types::global_dof_index> node_disp_control(m_fe.dofs_per_vertex);`：调用 C++ 标准库工具函数/容器接口。
- **[行 3549]** ``：空行，用于分隔逻辑块。
- **[行 3550]** `	    for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3551]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 3552]** `		if (   (std::fabs(vertex_itr->vertex()[0] - 470.0) < 1.0e-9)`：条件分支：根据当前状态选择执行路径。
- **[行 3553]** `		    && (std::fabs(vertex_itr->vertex()[1] - 250.0) < 1.0e-9) )`：调用 C++ 标准库工具函数/容器接口。
- **[行 3554]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 3555]** `		    node_disp_control = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3556]** `	            // node applied with y-displacement`：注释行，用于说明算法背景、假设或实现意图。
- **[行 3557]** `		    const double time_inc = m_time.get_delta_t();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3558]** `		    double disp_magnitude = m_time.get_magnitude();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3559]** ``：空行，用于分隔逻辑块。
- **[行 3560]** `		    m_constraints.add_line(node_disp_control[1]);`：函数调用语句，触发对应计算或操作。
- **[行 3561]** `		    m_constraints.set_inhomogeneity(node_disp_control[1], disp_magnitude*time_inc);`：函数调用语句，触发对应计算或操作。
- **[行 3562]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 3563]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 3564]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 3565]** `	else if (m_parameters.m_scenario == 11)`：多分支条件判断，处理备选情形。
- **[行 3566]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 3567]** `	    // Dirichlet B,C. right surface`：注释行，用于说明算法背景、假设或实现意图。
- **[行 3568]** `	    const int boundary_id_right_surface = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3569]** `	    VectorTools::interpolate_boundary_values(m_dof_handler,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3570]** `						     boundary_id_right_surface,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3571]** `						     Functions::ZeroFunction<dim>(m_n_components),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3572]** `						     m_constraints,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3573]** `						     m_fe.component_mask(displacements));`：函数调用语句，触发对应计算或操作。
- **[行 3574]** ``：空行，用于分隔逻辑块。
- **[行 3575]** `	    // Dirichlet B,C. left surface`：注释行，用于说明算法背景、假设或实现意图。
- **[行 3576]** `	    const int boundary_id_left_surface = 1;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3577]** `	    VectorTools::interpolate_boundary_values(m_dof_handler,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3578]** `						     boundary_id_left_surface,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3579]** `						     Functions::ZeroFunction<dim>(m_n_components),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3580]** `						     m_constraints,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3581]** `						     m_fe.component_mask(x_displacement));`：函数调用语句，触发对应计算或操作。
- **[行 3582]** ``：空行，用于分隔逻辑块。
- **[行 3583]** `	    typename Triangulation<dim>::active_vertex_iterator vertex_itr;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3584]** `	    vertex_itr = m_triangulation.begin_active_vertex();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3585]** `	    std::vector<types::global_dof_index> node_rotate(m_fe.dofs_per_vertex);`：调用 C++ 标准库工具函数/容器接口。
- **[行 3586]** `	    double node_dist = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3587]** `	    double disp_mag = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3588]** `	    double angle_theta = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3589]** `	    double disp_y = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3590]** `	    double disp_z = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3591]** ``：空行，用于分隔逻辑块。
- **[行 3592]** `	    for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3593]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 3594]** `		if (std::fabs(vertex_itr->vertex()[0] - 0.0) < 1.0e-9)`：条件分支：根据当前状态选择执行路径。
- **[行 3595]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 3596]** `		    node_rotate = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3597]** `		    node_dist = std::sqrt(  vertex_itr->vertex()[1] * vertex_itr->vertex()[1]`：调用 C++ 标准库工具函数/容器接口。
- **[行 3598]** `			                  + vertex_itr->vertex()[2] * vertex_itr->vertex()[2]);`：函数调用语句，触发对应计算或操作。
- **[行 3599]** ``：空行，用于分隔逻辑块。
- **[行 3600]** `		    angle_theta = m_time.get_delta_t() * m_time.get_magnitude();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3601]** `		    disp_mag = node_dist * std::tan(angle_theta);`：调用 C++ 标准库工具函数/容器接口。
- **[行 3602]** ``：空行，用于分隔逻辑块。
- **[行 3603]** `		    if (node_dist > 0)`：条件分支：根据当前状态选择执行路径。
- **[行 3604]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 3605]** `		        disp_y = vertex_itr->vertex()[2]/node_dist * disp_mag;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3606]** `		        disp_z = -vertex_itr->vertex()[1]/node_dist * disp_mag;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3607]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 3608]** `		    else`：条件分支的兜底路径。
- **[行 3609]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 3610]** `			disp_y = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3611]** `			disp_z = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3612]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 3613]** ``：空行，用于分隔逻辑块。
- **[行 3614]** `		    m_constraints.add_line(node_rotate[1]);`：函数调用语句，触发对应计算或操作。
- **[行 3615]** `		    m_constraints.set_inhomogeneity(node_rotate[1], disp_y);`：函数调用语句，触发对应计算或操作。
- **[行 3616]** ``：空行，用于分隔逻辑块。
- **[行 3617]** `		    m_constraints.add_line(node_rotate[2]);`：函数调用语句，触发对应计算或操作。
- **[行 3618]** `		    m_constraints.set_inhomogeneity(node_rotate[2], disp_z);`：函数调用语句，触发对应计算或操作。
- **[行 3619]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 3620]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 3621]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 3622]** `	else`：条件分支的兜底路径。
- **[行 3623]** `	  Assert(false, ExcMessage("The scenario has not been implemented!"));`：运行期断言/检查，验证输入与状态合法性。
- **[行 3624]** `      }`：作用域边界（代码块开始/结束）。
- **[行 3625]** `    else  // inhomogeneous constraints`：条件分支的兜底路径。
- **[行 3626]** `      {`：作用域边界（代码块开始/结束）。
- **[行 3627]** `        if (m_constraints.has_inhomogeneities())`：条件分支：根据当前状态选择执行路径。
- **[行 3628]** `          {`：作用域边界（代码块开始/结束）。
- **[行 3629]** `            AffineConstraints<double> homogeneous_constraints(m_constraints);`：函数调用语句，触发对应计算或操作。
- **[行 3630]** `            for (unsigned int dof = 0; dof != m_dof_handler.n_dofs(); ++dof)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3631]** `              if (homogeneous_constraints.is_inhomogeneously_constrained(dof))`：条件分支：根据当前状态选择执行路径。
- **[行 3632]** `                homogeneous_constraints.set_inhomogeneity(dof, 0.0);`：函数调用语句，触发对应计算或操作。
- **[行 3633]** `            m_constraints.clear();`：函数调用语句，触发对应计算或操作。
- **[行 3634]** `            m_constraints.copy_from(homogeneous_constraints);`：函数调用语句，触发对应计算或操作。
- **[行 3635]** `          }`：作用域边界（代码块开始/结束）。
- **[行 3636]** `      }`：作用域边界（代码块开始/结束）。
- **[行 3637]** `    m_constraints.close();`：函数调用语句，触发对应计算或操作。
- **[行 3638]** `  }`：作用域边界（代码块开始/结束）。
- **[行 3639]** ``：空行，用于分隔逻辑块。
- **[行 3640]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 3641]** `  void PhaseFieldMonolithicSolve<dim>::assemble_system_B0(const BlockVector<double> & solution_old)`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 3642]** `  {`：作用域边界（代码块开始/结束）。
- **[行 3643]** `    m_timer.enter_subsection("Assemble B0");`：函数调用语句，触发对应计算或操作。
- **[行 3644]** ``：空行，用于分隔逻辑块。
- **[行 3645]** `    m_tangent_matrix = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3646]** ``：空行，用于分隔逻辑块。
- **[行 3647]** `    const UpdateFlags uf_cell(update_values | update_gradients |`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3648]** `			      update_quadrature_points | update_JxW_values);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3649]** `    const UpdateFlags uf_face(update_values | update_normal_vectors |`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3650]** `                              update_JxW_values);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3651]** ``：空行，用于分隔逻辑块。
- **[行 3652]** `    PerTaskData_ASM per_task_data(m_fe.n_dofs_per_cell());`：函数调用语句，触发对应计算或操作。
- **[行 3653]** `    ScratchData_ASM scratch_data(m_fe, m_qf_cell, uf_cell, m_qf_face, uf_face, solution_old);`：函数调用语句，触发对应计算或操作。
- **[行 3654]** ``：空行，用于分隔逻辑块。
- **[行 3655]** `    auto worker =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3656]** `      [this](const typename DoFHandler<dim>::active_cell_iterator &cell,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3657]** `	     ScratchData_ASM & scratch,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3658]** `	     PerTaskData_ASM & data)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3659]** `      {`：作用域边界（代码块开始/结束）。
- **[行 3660]** `        this->assemble_system_B0_one_cell(cell, scratch, data);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 3661]** `      };`：作用域边界（代码块开始/结束）。
- **[行 3662]** ``：空行，用于分隔逻辑块。
- **[行 3663]** `    auto copier = [this](const PerTaskData_ASM &data)`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3664]** `      {`：作用域边界（代码块开始/结束）。
- **[行 3665]** `        this->m_constraints.distribute_local_to_global(data.m_cell_matrix,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3666]** `                                                       data.m_local_dof_indices,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3667]** `						       m_tangent_matrix);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3668]** `      };`：作用域边界（代码块开始/结束）。
- **[行 3669]** ``：空行，用于分隔逻辑块。
- **[行 3670]** `    WorkStream::run(`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3671]** `      m_dof_handler.active_cell_iterators(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3672]** `      worker,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3673]** `      copier,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3674]** `      scratch_data,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3675]** `      per_task_data);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3676]** ``：空行，用于分隔逻辑块。
- **[行 3677]** `    m_timer.leave_subsection();`：函数调用语句，触发对应计算或操作。
- **[行 3678]** `  }`：作用域边界（代码块开始/结束）。
- **[行 3679]** ``：空行，用于分隔逻辑块。
- **[行 3680]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 3681]** `  void PhaseFieldMonolithicSolve<dim>::assemble_system_rhs_BFGS_parallel(const BlockVector<double> & solution_old,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 3682]** `								         BlockVector<double> & system_rhs)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3683]** `  {`：作用域边界（代码块开始/结束）。
- **[行 3684]** `    m_timer.enter_subsection("Assemble RHS");`：函数调用语句，触发对应计算或操作。
- **[行 3685]** ``：空行，用于分隔逻辑块。
- **[行 3686]** `    system_rhs = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3687]** ``：空行，用于分隔逻辑块。
- **[行 3688]** `    const UpdateFlags uf_cell(update_values | update_gradients |`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3689]** `			      update_quadrature_points | update_JxW_values);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3690]** `    const UpdateFlags uf_face(update_values | update_normal_vectors |`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3691]** `			      update_JxW_values);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3692]** ``：空行，用于分隔逻辑块。
- **[行 3693]** `    PerTaskData_ASM_RHS_BFGS per_task_data(m_fe.n_dofs_per_cell());`：函数调用语句，触发对应计算或操作。
- **[行 3694]** `    ScratchData_ASM_RHS_BFGS scratch_data(m_fe, m_qf_cell, uf_cell, m_qf_face, uf_face, solution_old);`：函数调用语句，触发对应计算或操作。
- **[行 3695]** ``：空行，用于分隔逻辑块。
- **[行 3696]** `    auto worker =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3697]** `      [this](const typename DoFHandler<dim>::active_cell_iterator &cell,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3698]** `	     ScratchData_ASM_RHS_BFGS & scratch,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3699]** `	     PerTaskData_ASM_RHS_BFGS & data)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3700]** `      {`：作用域边界（代码块开始/结束）。
- **[行 3701]** `        this->assemble_system_rhs_BFGS_one_cell(cell, scratch, data);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 3702]** `      };`：作用域边界（代码块开始/结束）。
- **[行 3703]** ``：空行，用于分隔逻辑块。
- **[行 3704]** `    auto copier = [this, &system_rhs](const PerTaskData_ASM_RHS_BFGS &data)`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3705]** `      {`：作用域边界（代码块开始/结束）。
- **[行 3706]** `        this->m_constraints.distribute_local_to_global(data.m_cell_rhs,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3707]** `                                                       data.m_local_dof_indices,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3708]** `						       system_rhs);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3709]** `      };`：作用域边界（代码块开始/结束）。
- **[行 3710]** ``：空行，用于分隔逻辑块。
- **[行 3711]** `    WorkStream::run(`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3712]** `      m_dof_handler.active_cell_iterators(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3713]** `      worker,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3714]** `      copier,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3715]** `      scratch_data,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3716]** `      per_task_data);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3717]** ``：空行，用于分隔逻辑块。
- **[行 3718]** `    m_timer.leave_subsection();`：函数调用语句，触发对应计算或操作。
- **[行 3719]** `  }`：作用域边界（代码块开始/结束）。
- **[行 3720]** ``：空行，用于分隔逻辑块。
- **[行 3721]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 3722]** `  void PhaseFieldMonolithicSolve<dim>::assemble_system_rhs_BFGS_one_cell(`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 3723]** `      const typename DoFHandler<dim>::active_cell_iterator &cell,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3724]** `      ScratchData_ASM_RHS_BFGS & scratch,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3725]** `      PerTaskData_ASM_RHS_BFGS & data) const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3726]** `  {`：作用域边界（代码块开始/结束）。
- **[行 3727]** `    data.reset();`：函数调用语句，触发对应计算或操作。
- **[行 3728]** `    scratch.reset();`：函数调用语句，触发对应计算或操作。
- **[行 3729]** `    scratch.m_fe_values.reinit(cell);`：函数调用语句，触发对应计算或操作。
- **[行 3730]** `    cell->get_dof_indices(data.m_local_dof_indices);`：函数调用语句，触发对应计算或操作。
- **[行 3731]** ``：空行，用于分隔逻辑块。
- **[行 3732]** `    scratch.m_fe_values[m_d_fe].get_function_values(`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3733]** `      scratch.m_solution_previous_step, scratch.m_phasefield_previous_step_cell);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3734]** ``：空行，用于分隔逻辑块。
- **[行 3735]** `    const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3736]** `      m_quadrature_point_history.get_data(cell);`：函数调用语句，触发对应计算或操作。
- **[行 3737]** `    Assert(lqph.size() == m_n_q_points, ExcInternalError());`：运行期断言/检查，验证输入与状态合法性。
- **[行 3738]** ``：空行，用于分隔逻辑块。
- **[行 3739]** `    const double time_ramp = (m_time.current() / m_time.end());`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3740]** `    std::vector<Tensor<1, dim>> rhs_values(m_n_q_points);`：调用 C++ 标准库工具函数/容器接口。
- **[行 3741]** ``：空行，用于分隔逻辑块。
- **[行 3742]** `    right_hand_side(scratch.m_fe_values.get_quadrature_points(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3743]** `		    rhs_values,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3744]** `		    m_parameters.m_x_component*1.0,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3745]** `		    m_parameters.m_y_component*1.0,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3746]** `		    m_parameters.m_z_component*1.0);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3747]** ``：空行，用于分隔逻辑块。
- **[行 3748]** `    const double delta_time = m_time.get_delta_t();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3749]** ``：空行，用于分隔逻辑块。
- **[行 3750]** `    for (const unsigned int q_point : scratch.m_fe_values.quadrature_point_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3751]** `      {`：作用域边界（代码块开始/结束）。
- **[行 3752]** `        for (const unsigned int k : scratch.m_fe_values.dof_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3753]** `          {`：作用域边界（代码块开始/结束）。
- **[行 3754]** `            const unsigned int k_group = m_fe.system_to_base_index(k).first.first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3755]** ``：空行，用于分隔逻辑块。
- **[行 3756]** `            if (k_group == m_u_dof)`：条件分支：根据当前状态选择执行路径。
- **[行 3757]** `              {`：作用域边界（代码块开始/结束）。
- **[行 3758]** `                scratch.m_Nx_disp[q_point][k] =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3759]** `                  scratch.m_fe_values[m_u_fe].value(k, q_point);`：函数调用语句，触发对应计算或操作。
- **[行 3760]** `                scratch.m_grad_Nx_disp[q_point][k] =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3761]** `                  scratch.m_fe_values[m_u_fe].gradient(k, q_point);`：函数调用语句，触发对应计算或操作。
- **[行 3762]** `                scratch.m_symm_grad_Nx_disp[q_point][k] =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3763]** `                  symmetrize(scratch.m_grad_Nx_disp[q_point][k]);`：函数调用语句，触发对应计算或操作。
- **[行 3764]** `              }`：作用域边界（代码块开始/结束）。
- **[行 3765]** `            else if (k_group == m_d_dof)`：多分支条件判断，处理备选情形。
- **[行 3766]** `              {`：作用域边界（代码块开始/结束）。
- **[行 3767]** `		scratch.m_Nx_phasefield[q_point][k] =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3768]** `		  scratch.m_fe_values[m_d_fe].value(k, q_point);`：函数调用语句，触发对应计算或操作。
- **[行 3769]** `		scratch.m_grad_Nx_phasefield[q_point][k] =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3770]** `		  scratch.m_fe_values[m_d_fe].gradient(k, q_point);`：函数调用语句，触发对应计算或操作。
- **[行 3771]** `              }`：作用域边界（代码块开始/结束）。
- **[行 3772]** `            else`：条件分支的兜底路径。
- **[行 3773]** `              Assert(k_group <= m_d_dof, ExcInternalError());`：运行期断言/检查，验证输入与状态合法性。
- **[行 3774]** `          }`：作用域边界（代码块开始/结束）。
- **[行 3775]** `      }`：作用域边界（代码块开始/结束）。
- **[行 3776]** ``：空行，用于分隔逻辑块。
- **[行 3777]** `    for (const unsigned int q_point : scratch.m_fe_values.quadrature_point_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3778]** `      {`：作用域边界（代码块开始/结束）。
- **[行 3779]** `	const double length_scale            = lqph[q_point]->get_length_scale();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3780]** `	const double gc                      = lqph[q_point]->get_critical_energy_release_rate();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3781]** `	const double eta                     = lqph[q_point]->get_viscosity();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3782]** `	const double current_positive_strain_energy = lqph[q_point]->get_current_positive_strain_energy();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3783]** ``：空行，用于分隔逻辑块。
- **[行 3784]** `	const double phasefield_value        = lqph[q_point]->get_phase_field_value();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3785]** `	const Tensor<1, dim> phasefield_grad = lqph[q_point]->get_phase_field_gradient();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3786]** ``：空行，用于分隔逻辑块。
- **[行 3787]** `        const std::vector<double>         &      N_phasefield = scratch.m_Nx_phasefield[q_point];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3788]** `        const std::vector<Tensor<1, dim>> & grad_N_phasefield = scratch.m_grad_Nx_phasefield[q_point];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3789]** `        const double                old_phasefield = scratch.m_phasefield_previous_step_cell[q_point];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3790]** ``：空行，用于分隔逻辑块。
- **[行 3791]** `        const SymmetricTensor<2, dim> & cauchy_stress = lqph[q_point]->get_cauchy_stress();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3792]** ``：空行，用于分隔逻辑块。
- **[行 3793]** `        const std::vector<Tensor<1,dim>> & N_disp = scratch.m_Nx_disp[q_point];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3794]** `        const std::vector<SymmetricTensor<2, dim>> & symm_grad_N_disp =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3795]** `          scratch.m_symm_grad_Nx_disp[q_point];`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3796]** `        const double JxW = scratch.m_fe_values.JxW(q_point);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3797]** ``：空行，用于分隔逻辑块。
- **[行 3798]** `        SymmetricTensor<2, dim> symm_grad_Nx_i_x_C;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3799]** ``：空行，用于分隔逻辑块。
- **[行 3800]** `        for (const unsigned int i : scratch.m_fe_values.dof_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3801]** `          {`：作用域边界（代码块开始/结束）。
- **[行 3802]** `            const unsigned int i_group = m_fe.system_to_base_index(i).first.first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3803]** ``：空行，用于分隔逻辑块。
- **[行 3804]** `            if (i_group == m_u_dof)`：条件分支：根据当前状态选择执行路径。
- **[行 3805]** `              {`：作用域边界（代码块开始/结束）。
- **[行 3806]** `                data.m_cell_rhs(i) += (symm_grad_N_disp[i] * cauchy_stress) * JxW;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3807]** ``：空行，用于分隔逻辑块。
- **[行 3808]** `		// contributions from the body force to right-hand side`：注释行，用于说明算法背景、假设或实现意图。
- **[行 3809]** `		data.m_cell_rhs(i) -= N_disp[i] * rhs_values[q_point] * JxW;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3810]** `              }`：作用域边界（代码块开始/结束）。
- **[行 3811]** `            else if (i_group == m_d_dof)`：多分支条件判断，处理备选情形。
- **[行 3812]** `              {`：作用域边界（代码块开始/结束）。
- **[行 3813]** `    	        data.m_cell_rhs(i) += (    gc * length_scale * grad_N_phasefield[i] * phasefield_grad`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3814]** `    	                                +  (   gc / length_scale * phasefield_value`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3815]** `					     + eta / delta_time  * (phasefield_value - old_phasefield)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3816]** `					     + degradation_function_derivative(phasefield_value)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3817]** `					     * current_positive_strain_energy )`：注释行，用于说明算法背景、假设或实现意图。
- **[行 3818]** `					  * N_phasefield[i]`：注释行，用于说明算法背景、假设或实现意图。
- **[行 3819]** `				      ) * JxW;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3820]** `              }`：作用域边界（代码块开始/结束）。
- **[行 3821]** `            else`：条件分支的兜底路径。
- **[行 3822]** `              Assert(i_group <= m_d_dof, ExcInternalError());`：运行期断言/检查，验证输入与状态合法性。
- **[行 3823]** `          }  // i`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3824]** `      }  // q_point`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3825]** ``：空行，用于分隔逻辑块。
- **[行 3826]** `    // if there is surface pressure, this surface pressure always applied to the`：注释行，用于说明算法背景、假设或实现意图。
- **[行 3827]** `    // reference configuration`：注释行，用于说明算法背景、假设或实现意图。
- **[行 3828]** `    const unsigned int face_pressure_id = 100;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3829]** `    const double p0 = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3830]** ``：空行，用于分隔逻辑块。
- **[行 3831]** `    for (const auto &face : cell->face_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3832]** `      if (face->at_boundary() && face->boundary_id() == face_pressure_id)`：条件分支：根据当前状态选择执行路径。
- **[行 3833]** `        {`：作用域边界（代码块开始/结束）。
- **[行 3834]** `          scratch.m_fe_face_values.reinit(cell, face);`：函数调用语句，触发对应计算或操作。
- **[行 3835]** ``：空行，用于分隔逻辑块。
- **[行 3836]** `          for (const unsigned int f_q_point : scratch.m_fe_face_values.quadrature_point_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3837]** `            {`：作用域边界（代码块开始/结束）。
- **[行 3838]** `              const Tensor<1, dim> &N = scratch.m_fe_face_values.normal_vector(f_q_point);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3839]** ``：空行，用于分隔逻辑块。
- **[行 3840]** `              const double         pressure  = p0 * time_ramp;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3841]** `              const Tensor<1, dim> traction  = pressure * N;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3842]** ``：空行，用于分隔逻辑块。
- **[行 3843]** `              for (const unsigned int i : scratch.m_fe_values.dof_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3844]** `                {`：作用域边界（代码块开始/结束）。
- **[行 3845]** `                  const unsigned int i_group = m_fe.system_to_base_index(i).first.first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3846]** ``：空行，用于分隔逻辑块。
- **[行 3847]** `                  if (i_group == m_u_dof)`：条件分支：根据当前状态选择执行路径。
- **[行 3848]** `                    {`：作用域边界（代码块开始/结束）。
- **[行 3849]** `    		      const unsigned int component_i = m_fe.system_to_component_index(i).first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3850]** `    		      const double Ni = scratch.m_fe_face_values.shape_value(i, f_q_point);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3851]** `    		      const double JxW = scratch.m_fe_face_values.JxW(f_q_point);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3852]** `    		      data.m_cell_rhs(i) -= (Ni * traction[component_i]) * JxW;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3853]** `                    }`：作用域边界（代码块开始/结束）。
- **[行 3854]** `                }`：作用域边界（代码块开始/结束）。
- **[行 3855]** `            }`：作用域边界（代码块开始/结束）。
- **[行 3856]** `        }`：作用域边界（代码块开始/结束）。
- **[行 3857]** `  }`：作用域边界（代码块开始/结束）。
- **[行 3858]** ``：空行，用于分隔逻辑块。
- **[行 3859]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 3860]** `  void PhaseFieldMonolithicSolve<dim>::assemble_system_B0_one_cell(`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 3861]** `      const typename DoFHandler<dim>::active_cell_iterator &cell,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3862]** `      ScratchData_ASM & scratch,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3863]** `      PerTaskData_ASM & data) const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3864]** `  {`：作用域边界（代码块开始/结束）。
- **[行 3865]** `    data.reset();`：函数调用语句，触发对应计算或操作。
- **[行 3866]** `    scratch.reset();`：函数调用语句，触发对应计算或操作。
- **[行 3867]** `    scratch.m_fe_values.reinit(cell);`：函数调用语句，触发对应计算或操作。
- **[行 3868]** `    cell->get_dof_indices(data.m_local_dof_indices);`：函数调用语句，触发对应计算或操作。
- **[行 3869]** ``：空行，用于分隔逻辑块。
- **[行 3870]** `    scratch.m_fe_values[m_d_fe].get_function_values(`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3871]** `      scratch.m_solution_previous_step, scratch.m_phasefield_previous_step_cell);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3872]** ``：空行，用于分隔逻辑块。
- **[行 3873]** `    const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3874]** `      m_quadrature_point_history.get_data(cell);`：函数调用语句，触发对应计算或操作。
- **[行 3875]** `    Assert(lqph.size() == m_n_q_points, ExcInternalError());`：运行期断言/检查，验证输入与状态合法性。
- **[行 3876]** ``：空行，用于分隔逻辑块。
- **[行 3877]** `    const double delta_time = m_time.get_delta_t();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3878]** ``：空行，用于分隔逻辑块。
- **[行 3879]** `    for (const unsigned int q_point : scratch.m_fe_values.quadrature_point_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3880]** `      {`：作用域边界（代码块开始/结束）。
- **[行 3881]** `        for (const unsigned int k : scratch.m_fe_values.dof_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3882]** `          {`：作用域边界（代码块开始/结束）。
- **[行 3883]** `            const unsigned int k_group = m_fe.system_to_base_index(k).first.first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3884]** ``：空行，用于分隔逻辑块。
- **[行 3885]** `            if (k_group == m_u_dof)`：条件分支：根据当前状态选择执行路径。
- **[行 3886]** `              {`：作用域边界（代码块开始/结束）。
- **[行 3887]** `                scratch.m_Nx_disp[q_point][k] =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3888]** `                  scratch.m_fe_values[m_u_fe].value(k, q_point);`：函数调用语句，触发对应计算或操作。
- **[行 3889]** `                scratch.m_grad_Nx_disp[q_point][k] =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3890]** `                  scratch.m_fe_values[m_u_fe].gradient(k, q_point);`：函数调用语句，触发对应计算或操作。
- **[行 3891]** `                scratch.m_symm_grad_Nx_disp[q_point][k] =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3892]** `                  symmetrize(scratch.m_grad_Nx_disp[q_point][k]);`：函数调用语句，触发对应计算或操作。
- **[行 3893]** `              }`：作用域边界（代码块开始/结束）。
- **[行 3894]** `            else if (k_group == m_d_dof)`：多分支条件判断，处理备选情形。
- **[行 3895]** `              {`：作用域边界（代码块开始/结束）。
- **[行 3896]** `		scratch.m_Nx_phasefield[q_point][k] =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3897]** `		  scratch.m_fe_values[m_d_fe].value(k, q_point);`：函数调用语句，触发对应计算或操作。
- **[行 3898]** `		scratch.m_grad_Nx_phasefield[q_point][k] =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3899]** `		  scratch.m_fe_values[m_d_fe].gradient(k, q_point);`：函数调用语句，触发对应计算或操作。
- **[行 3900]** `              }`：作用域边界（代码块开始/结束）。
- **[行 3901]** `            else`：条件分支的兜底路径。
- **[行 3902]** `              Assert(k_group <= m_d_dof, ExcInternalError());`：运行期断言/检查，验证输入与状态合法性。
- **[行 3903]** `          }`：作用域边界（代码块开始/结束）。
- **[行 3904]** `      }`：作用域边界（代码块开始/结束）。
- **[行 3905]** ``：空行，用于分隔逻辑块。
- **[行 3906]** `    for (const unsigned int q_point : scratch.m_fe_values.quadrature_point_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3907]** `      {`：作用域边界（代码块开始/结束）。
- **[行 3908]** `	const double length_scale            = lqph[q_point]->get_length_scale();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3909]** `	const double gc                      = lqph[q_point]->get_critical_energy_release_rate();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3910]** `	const double eta                     = lqph[q_point]->get_viscosity();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3911]** `	const double current_positive_strain_energy = lqph[q_point]->get_current_positive_strain_energy();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3912]** ``：空行，用于分隔逻辑块。
- **[行 3913]** `	const double phasefield_value        = lqph[q_point]->get_phase_field_value();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3914]** ``：空行，用于分隔逻辑块。
- **[行 3915]** `        const std::vector<double>         &      N_phasefield = scratch.m_Nx_phasefield[q_point];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3916]** `        const std::vector<Tensor<1, dim>> & grad_N_phasefield = scratch.m_grad_Nx_phasefield[q_point];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3917]** ``：空行，用于分隔逻辑块。
- **[行 3918]** `        //const SymmetricTensor<2, dim> & cauchy_stress_positive = lqph[q_point]->get_cauchy_stress_positive();`：注释行，用于说明算法背景、假设或实现意图。
- **[行 3919]** `        const SymmetricTensor<4, dim> & mechanical_C  = lqph[q_point]->get_mechanical_C();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3920]** ``：空行，用于分隔逻辑块。
- **[行 3921]** `        const std::vector<SymmetricTensor<2, dim>> & symm_grad_N_disp =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3922]** `          scratch.m_symm_grad_Nx_disp[q_point];`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3923]** `        const double JxW = scratch.m_fe_values.JxW(q_point);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3924]** ``：空行，用于分隔逻辑块。
- **[行 3925]** `        SymmetricTensor<2, dim> symm_grad_Nx_i_x_C;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3926]** ``：空行，用于分隔逻辑块。
- **[行 3927]** `        for (const unsigned int i : scratch.m_fe_values.dof_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3928]** `          {`：作用域边界（代码块开始/结束）。
- **[行 3929]** `            const unsigned int i_group = m_fe.system_to_base_index(i).first.first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3930]** ``：空行，用于分隔逻辑块。
- **[行 3931]** `            if (i_group == m_u_dof)`：条件分支：根据当前状态选择执行路径。
- **[行 3932]** `              {`：作用域边界（代码块开始/结束）。
- **[行 3933]** `                symm_grad_Nx_i_x_C = symm_grad_N_disp[i] * mechanical_C;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3934]** `              }`：作用域边界（代码块开始/结束）。
- **[行 3935]** ``：空行，用于分隔逻辑块。
- **[行 3936]** `            for (const unsigned int j : scratch.m_fe_values.dof_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 3937]** `              {`：作用域边界（代码块开始/结束）。
- **[行 3938]** `                const unsigned int j_group = m_fe.system_to_base_index(j).first.first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3939]** ``：空行，用于分隔逻辑块。
- **[行 3940]** `                if ((i_group == j_group) && (i_group == m_u_dof))`：条件分支：根据当前状态选择执行路径。
- **[行 3941]** `                  {`：作用域边界（代码块开始/结束）。
- **[行 3942]** `                    data.m_cell_matrix(i, j) += symm_grad_Nx_i_x_C * symm_grad_N_disp[j] * JxW;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3943]** `                  }`：作用域边界（代码块开始/结束）。
- **[行 3944]** `                else if ((i_group == j_group) && (i_group == m_d_dof))`：多分支条件判断，处理备选情形。
- **[行 3945]** `                  {`：作用域边界（代码块开始/结束）。
- **[行 3946]** `                    data.m_cell_matrix(i, j) += (  (   gc/length_scale + eta/delta_time`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3947]** `                	                             + degradation_function_2nd_order_derivative(phasefield_value)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3948]** `						     * current_positive_strain_energy  )`：注释行，用于说明算法背景、假设或实现意图。
- **[行 3949]** `                	                          * N_phasefield[i] * N_phasefield[j]`：注释行，用于说明算法背景、假设或实现意图。
- **[行 3950]** `					          + gc * length_scale * grad_N_phasefield[i] * grad_N_phasefield[j]`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3951]** `					        ) * JxW;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3952]** `                  }`：作用域边界（代码块开始/结束）。
- **[行 3953]** `                else`：条件分支的兜底路径。
- **[行 3954]** `                  Assert((i_group <= m_d_dof) && (j_group <= m_d_dof),`：运行期断言/检查，验证输入与状态合法性。
- **[行 3955]** `                         ExcInternalError());`：函数调用语句，触发对应计算或操作。
- **[行 3956]** `              } // j`：作用域边界（代码块开始/结束）。
- **[行 3957]** `          }  // i`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3958]** `      }  // q_point`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3959]** `  }`：作用域边界（代码块开始/结束）。
- **[行 3960]** ``：空行，用于分隔逻辑块。
- **[行 3961]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 3962]** `  void PhaseFieldMonolithicSolve<dim>::assemble_system_rhs_BFGS(const BlockVector<double> & solution_old,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 3963]** `								BlockVector<double> & system_rhs)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3964]** `  {`：作用域边界（代码块开始/结束）。
- **[行 3965]** `    m_timer.enter_subsection("Assemble RHS");`：函数调用语句，触发对应计算或操作。
- **[行 3966]** ``：空行，用于分隔逻辑块。
- **[行 3967]** `    system_rhs = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3968]** ``：空行，用于分隔逻辑块。
- **[行 3969]** `    Vector<double> cell_rhs(m_dofs_per_cell);`：函数调用语句，触发对应计算或操作。
- **[行 3970]** `    std::vector<types::global_dof_index> local_dof_indices(m_dofs_per_cell);`：调用 C++ 标准库工具函数/容器接口。
- **[行 3971]** ``：空行，用于分隔逻辑块。
- **[行 3972]** `    const double time_ramp = (m_time.current() / m_time.end());`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3973]** `    const double delta_time = m_time.get_delta_t();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 3974]** ``：空行，用于分隔逻辑块。
- **[行 3975]** `    std::vector<Tensor<1, dim>> rhs_values(m_n_q_points);`：调用 C++ 标准库工具函数/容器接口。
- **[行 3976]** `    const UpdateFlags uf_cell(update_values | update_gradients |`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3977]** `			      update_quadrature_points | update_JxW_values);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3978]** `    const UpdateFlags uf_face(update_values | update_normal_vectors |`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3979]** `			      update_JxW_values);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3980]** ``：空行，用于分隔逻辑块。
- **[行 3981]** `    FEValues<dim> fe_values(m_fe, m_qf_cell, uf_cell);`：函数调用语句，触发对应计算或操作。
- **[行 3982]** `    FEFaceValues<dim> fe_face_values(m_fe, m_qf_face, uf_face);`：函数调用语句，触发对应计算或操作。
- **[行 3983]** ``：空行，用于分隔逻辑块。
- **[行 3984]** `    // shape function values for displacement field`：注释行，用于说明算法背景、假设或实现意图。
- **[行 3985]** `    std::vector<std::vector<Tensor<1, dim>>>`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3986]** `      Nx_disp(m_qf_cell.size(), std::vector<Tensor<1, dim>>(m_dofs_per_cell));`：调用 C++ 标准库工具函数/容器接口。
- **[行 3987]** `    std::vector<std::vector<Tensor<2, dim>>>`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3988]** `      grad_Nx_disp(m_qf_cell.size(), std::vector<Tensor<2, dim>>(m_dofs_per_cell));`：调用 C++ 标准库工具函数/容器接口。
- **[行 3989]** `    std::vector<std::vector<SymmetricTensor<2, dim>>>`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3990]** `      symm_grad_Nx_disp(m_qf_cell.size(), std::vector<SymmetricTensor<2, dim>>(m_dofs_per_cell));`：调用 C++ 标准库工具函数/容器接口。
- **[行 3991]** ``：空行，用于分隔逻辑块。
- **[行 3992]** `    // shape function values for phase field`：注释行，用于说明算法背景、假设或实现意图。
- **[行 3993]** `    std::vector<std::vector<double>>`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3994]** `      Nx_phasefield(m_qf_cell.size(), std::vector<double>(m_dofs_per_cell));`：调用 C++ 标准库工具函数/容器接口。
- **[行 3995]** `    std::vector<std::vector<Tensor<1, dim>>>`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 3996]** `      grad_Nx_phasefield(m_qf_cell.size(), std::vector<Tensor<1, dim>>(m_dofs_per_cell));`：调用 C++ 标准库工具函数/容器接口。
- **[行 3997]** ``：空行，用于分隔逻辑块。
- **[行 3998]** `    std::vector<double> phasefield_previous_step_cell(m_qf_cell.size());`：调用 C++ 标准库工具函数/容器接口。
- **[行 3999]** ``：空行，用于分隔逻辑块。
- **[行 4000]** `    for (const auto &cell : m_dof_handler.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4001]** `      {`：作用域边界（代码块开始/结束）。
- **[行 4002]** `	const std::vector<std::shared_ptr< PointHistory<dim>>> lqph =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4003]** `	  m_quadrature_point_history.get_data(cell);`：函数调用语句，触发对应计算或操作。
- **[行 4004]** `	Assert(lqph.size() == m_n_q_points, ExcInternalError());`：运行期断言/检查，验证输入与状态合法性。
- **[行 4005]** ``：空行，用于分隔逻辑块。
- **[行 4006]** `	cell_rhs = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4007]** `	fe_values.reinit(cell);`：函数调用语句，触发对应计算或操作。
- **[行 4008]** `	right_hand_side(fe_values.get_quadrature_points(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4009]** `			rhs_values,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4010]** `			m_parameters.m_x_component*time_ramp,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4011]** `			m_parameters.m_y_component*time_ramp,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4012]** `			m_parameters.m_z_component*time_ramp);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4013]** ``：空行，用于分隔逻辑块。
- **[行 4014]** `	fe_values[m_d_fe].get_function_values(`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4015]** `	    solution_old, phasefield_previous_step_cell);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4016]** ``：空行，用于分隔逻辑块。
- **[行 4017]** `	for (const unsigned int q_point : fe_values.quadrature_point_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4018]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 4019]** `	    for (const unsigned int k : fe_values.dof_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4020]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 4021]** `		const unsigned int k_group = m_fe.system_to_base_index(k).first.first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4022]** ``：空行，用于分隔逻辑块。
- **[行 4023]** `		if (k_group == m_u_dof)`：条件分支：根据当前状态选择执行路径。
- **[行 4024]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 4025]** `		    Nx_disp[q_point][k] = fe_values[m_u_fe].value(k, q_point);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4026]** `		    grad_Nx_disp[q_point][k] = fe_values[m_u_fe].gradient(k, q_point);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4027]** `		    symm_grad_Nx_disp[q_point][k] = symmetrize(grad_Nx_disp[q_point][k]);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4028]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 4029]** `		else if (k_group == m_d_dof)`：多分支条件判断，处理备选情形。
- **[行 4030]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 4031]** `		    Nx_phasefield[q_point][k] = fe_values[m_d_fe].value(k, q_point);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4032]** `		    grad_Nx_phasefield[q_point][k] = fe_values[m_d_fe].gradient(k, q_point);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4033]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 4034]** `		else`：条件分支的兜底路径。
- **[行 4035]** `		  Assert(k_group <= m_d_dof, ExcInternalError());`：运行期断言/检查，验证输入与状态合法性。
- **[行 4036]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 4037]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 4038]** ``：空行，用于分隔逻辑块。
- **[行 4039]** `	for (const unsigned int q_point : fe_values.quadrature_point_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4040]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 4041]** `	    const double length_scale            = lqph[q_point]->get_length_scale();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4042]** `	    const double gc                      = lqph[q_point]->get_critical_energy_release_rate();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4043]** `	    const double eta                     = lqph[q_point]->get_viscosity();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4044]** `	    const double current_positive_strain_energy = lqph[q_point]->get_current_positive_strain_energy();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4045]** ``：空行，用于分隔逻辑块。
- **[行 4046]** `	    const double phasefield_value        = lqph[q_point]->get_phase_field_value();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4047]** `	    const Tensor<1, dim> phasefield_grad = lqph[q_point]->get_phase_field_gradient();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4048]** ``：空行，用于分隔逻辑块。
- **[行 4049]** `	    const std::vector<double>         &      N_phasefield = Nx_phasefield[q_point];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4050]** `	    const std::vector<Tensor<1, dim>> & grad_N_phasefield = grad_Nx_phasefield[q_point];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4051]** `	    const double                old_phasefield = phasefield_previous_step_cell[q_point];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4052]** ``：空行，用于分隔逻辑块。
- **[行 4053]** `	    const SymmetricTensor<2, dim> & cauchy_stress = lqph[q_point]->get_cauchy_stress();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4054]** ``：空行，用于分隔逻辑块。
- **[行 4055]** `	    const std::vector<Tensor<1,dim>> & N = Nx_disp[q_point];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4056]** `	    const std::vector<SymmetricTensor<2, dim>> & symm_grad_N = symm_grad_Nx_disp[q_point];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4057]** `	    const double JxW = fe_values.JxW(q_point);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4058]** ``：空行，用于分隔逻辑块。
- **[行 4059]** `	    for (const unsigned int i : fe_values.dof_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4060]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 4061]** `		const unsigned int i_group = m_fe.system_to_base_index(i).first.first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4062]** ``：空行，用于分隔逻辑块。
- **[行 4063]** `		if (i_group == m_u_dof)`：条件分支：根据当前状态选择执行路径。
- **[行 4064]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 4065]** `		    cell_rhs(i) += (symm_grad_N[i] * cauchy_stress) * JxW;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4066]** `		    // contributions from the body force to right-hand side`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4067]** `		    cell_rhs(i) -= N[i] * rhs_values[q_point] * JxW;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4068]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 4069]** `		else if (i_group == m_d_dof)`：多分支条件判断，处理备选情形。
- **[行 4070]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 4071]** `		    cell_rhs(i) += (    gc * length_scale * grad_N_phasefield[i] * phasefield_grad`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4072]** `	    	                     +  (   gc / length_scale * phasefield_value`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4073]** `			                  + eta / delta_time  * (phasefield_value - old_phasefield)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4074]** `				          + degradation_function_derivative(phasefield_value)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4075]** `					  * current_positive_strain_energy )`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4076]** `				     * N_phasefield[i]`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4077]** `				   ) * JxW;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4078]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 4079]** `		else`：条件分支的兜底路径。
- **[行 4080]** `		  Assert(i_group <= m_d_dof, ExcInternalError());`：运行期断言/检查，验证输入与状态合法性。
- **[行 4081]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 4082]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 4083]** ``：空行，用于分隔逻辑块。
- **[行 4084]** `	// if there is surface pressure, this surface pressure always applied to the`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4085]** `	// reference configuration`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4086]** `	const unsigned int face_pressure_id = 100;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4087]** `	const double p0 = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4088]** ``：空行，用于分隔逻辑块。
- **[行 4089]** `	for (const auto &face : cell->face_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4090]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 4091]** `	    if (face->at_boundary() && face->boundary_id() == face_pressure_id)`：条件分支：根据当前状态选择执行路径。
- **[行 4092]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 4093]** `		fe_face_values.reinit(cell, face);`：函数调用语句，触发对应计算或操作。
- **[行 4094]** ``：空行，用于分隔逻辑块。
- **[行 4095]** `		for (const unsigned int f_q_point : fe_face_values.quadrature_point_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4096]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 4097]** `		    const Tensor<1, dim> &N = fe_face_values.normal_vector(f_q_point);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4098]** ``：空行，用于分隔逻辑块。
- **[行 4099]** `		    const double         pressure  = p0 * time_ramp;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4100]** `		    const Tensor<1, dim> traction  = pressure * N;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4101]** ``：空行，用于分隔逻辑块。
- **[行 4102]** `		    for (const unsigned int i : fe_values.dof_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4103]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 4104]** `			const unsigned int i_group = m_fe.system_to_base_index(i).first.first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4105]** ``：空行，用于分隔逻辑块。
- **[行 4106]** `			if (i_group == m_u_dof)`：条件分支：根据当前状态选择执行路径。
- **[行 4107]** `			  {`：作用域边界（代码块开始/结束）。
- **[行 4108]** `			    const unsigned int component_i = m_fe.system_to_component_index(i).first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4109]** `			    const double Ni = fe_face_values.shape_value(i, f_q_point);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4110]** `			    const double JxW = fe_face_values.JxW(f_q_point);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4111]** `			    cell_rhs(i) -= (Ni * traction[component_i]) * JxW;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4112]** `			  }`：作用域边界（代码块开始/结束）。
- **[行 4113]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 4114]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 4115]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 4116]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 4117]** ``：空行，用于分隔逻辑块。
- **[行 4118]** `	cell->get_dof_indices(local_dof_indices);`：函数调用语句，触发对应计算或操作。
- **[行 4119]** `	for (const unsigned int i : fe_values.dof_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4120]** `	  system_rhs(local_dof_indices[i]) += cell_rhs(i);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4121]** `      } // for (const auto &cell : m_dof_handler.active_cell_iterators())`：作用域边界（代码块开始/结束）。
- **[行 4122]** ``：空行，用于分隔逻辑块。
- **[行 4123]** `    m_timer.leave_subsection();`：函数调用语句，触发对应计算或操作。
- **[行 4124]** `  }`：作用域边界（代码块开始/结束）。
- **[行 4125]** ``：空行，用于分隔逻辑块。
- **[行 4126]** ``：空行，用于分隔逻辑块。
- **[行 4127]** ``：空行，用于分隔逻辑块。
- **[行 4128]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 4129]** `  double PhaseFieldMonolithicSolve<dim>::line_search_stepsize_gradient_based(const BlockVector<double> & BFGS_p_vector,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4130]** `				                                             const BlockVector<double> & solution_delta)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4131]** `  {`：作用域边界（代码块开始/结束）。
- **[行 4132]** `    BlockVector<double> g_old(m_system_rhs);`：函数调用语句，触发对应计算或操作。
- **[行 4133]** ``：空行，用于分隔逻辑块。
- **[行 4134]** `    // BFGS_p_vector is the search direction`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4135]** `    BlockVector<double> solution_delta_trial(solution_delta);`：函数调用语句，触发对应计算或操作。
- **[行 4136]** `    // take a full step size 1.0`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4137]** `    solution_delta_trial.add(1.0, BFGS_p_vector);`：函数调用语句，触发对应计算或操作。
- **[行 4138]** ``：空行，用于分隔逻辑块。
- **[行 4139]** `    update_qph_incremental(solution_delta_trial, m_solution);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4140]** ``：空行，用于分隔逻辑块。
- **[行 4141]** `    BlockVector<double> g_new(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 4142]** `    assemble_system_rhs_BFGS_parallel(m_solution, g_new);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4143]** ``：空行，用于分隔逻辑块。
- **[行 4144]** `    BlockVector<double> y_old(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 4145]** ``：空行，用于分隔逻辑块。
- **[行 4146]** `    y_old = g_new - g_old;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4147]** ``：空行，用于分隔逻辑块。
- **[行 4148]** `    double alpha = 1.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4149]** ``：空行，用于分隔逻辑块。
- **[行 4150]** `    double alpha_old = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4151]** ``：空行，用于分隔逻辑块。
- **[行 4152]** `    double delta_alpha_old = alpha - alpha_old;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4153]** ``：空行，用于分隔逻辑块。
- **[行 4154]** `    double delta_alpha_new;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4155]** ``：空行，用于分隔逻辑块。
- **[行 4156]** `    unsigned int ls_max = 10;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4157]** ``：空行，用于分隔逻辑块。
- **[行 4158]** `    unsigned int i = 1;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4159]** ``：空行，用于分隔逻辑块。
- **[行 4160]** `    for (; i <= ls_max; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4161]** `      {`：作用域边界（代码块开始/结束）。
- **[行 4162]** `	delta_alpha_new = -delta_alpha_old`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4163]** `	                * (g_new * BFGS_p_vector)/(y_old * BFGS_p_vector);`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4164]** `	alpha += delta_alpha_new;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4165]** ``：空行，用于分隔逻辑块。
- **[行 4166]** `	if (std::fabs(delta_alpha_new) < 1.0e-5)`：条件分支：根据当前状态选择执行路径。
- **[行 4167]** `	  break;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4168]** ``：空行，用于分隔逻辑块。
- **[行 4169]** `        if (i == ls_max)`：条件分支：根据当前状态选择执行路径。
- **[行 4170]** `          {`：作用域边界（代码块开始/结束）。
- **[行 4171]** `            alpha = 1.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4172]** `            break;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4173]** `          }`：作用域边界（代码块开始/结束）。
- **[行 4174]** ``：空行，用于分隔逻辑块。
- **[行 4175]** `        g_old = g_new;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4176]** ``：空行，用于分隔逻辑块。
- **[行 4177]** `        // BFGS_p_vector is the search direction`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4178]** `        solution_delta_trial = solution_delta;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4179]** `        solution_delta_trial.add(alpha, BFGS_p_vector);`：函数调用语句，触发对应计算或操作。
- **[行 4180]** `        update_qph_incremental(solution_delta_trial, m_solution);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4181]** `        assemble_system_rhs_BFGS_parallel(m_solution, g_new);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4182]** ``：空行，用于分隔逻辑块。
- **[行 4183]** `        y_old = g_new - g_old;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4184]** ``：空行，用于分隔逻辑块。
- **[行 4185]** `        delta_alpha_old = delta_alpha_new;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4186]** `      }`：作用域边界（代码块开始/结束）。
- **[行 4187]** ``：空行，用于分隔逻辑块。
- **[行 4188]** `    if (alpha < 1.0e-3)`：条件分支：根据当前状态选择执行路径。
- **[行 4189]** `      alpha = 1.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4190]** ``：空行，用于分隔逻辑块。
- **[行 4191]** `    //num_ls = i;`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4192]** `    return alpha;`：返回当前函数结果。
- **[行 4193]** `  }`：作用域边界（代码块开始/结束）。
- **[行 4194]** ``：空行，用于分隔逻辑块。
- **[行 4195]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 4196]** `  double PhaseFieldMonolithicSolve<dim>::line_search_stepsize_strong_wolfe(const double phi_0,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4197]** `				                                           const double phi_0_prime,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4198]** `				                                           const BlockVector<double> & BFGS_p_vector,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4199]** `				                                           const BlockVector<double> & solution_delta)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4200]** `  {`：作用域边界（代码块开始/结束）。
- **[行 4201]** `    //AssertThrow(phi_0_prime < 0,`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4202]** `    //            ExcMessage("The derivative of phi at alpha = 0 should be negative!"));`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4203]** ``：空行，用于分隔逻辑块。
- **[行 4204]** `    // Some line search parameters`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4205]** `    const double c1 = 0.0001;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4206]** `    const double c2 = 0.9;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4207]** `    const double alpha_max = 1.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4208]** `    const unsigned int max_iter = 20;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4209]** `    double alpha = 1.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4210]** ``：空行，用于分隔逻辑块。
- **[行 4211]** `    double phi_old = phi_0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4212]** `    double phi_prime_old = phi_0_prime;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4213]** `    double alpha_old = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4214]** ``：空行，用于分隔逻辑块。
- **[行 4215]** `    double phi, phi_prime;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4216]** ``：空行，用于分隔逻辑块。
- **[行 4217]** `    std::pair<double, double> current_phi_phi_prime;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4218]** ``：空行，用于分隔逻辑块。
- **[行 4219]** `    unsigned int i = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4220]** `    for (; i < max_iter; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4221]** `      {`：作用域边界（代码块开始/结束）。
- **[行 4222]** `	current_phi_phi_prime = calculate_phi_and_phi_prime(alpha, BFGS_p_vector, solution_delta);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4223]** `	phi = current_phi_phi_prime.first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4224]** `	phi_prime = current_phi_phi_prime.second;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4225]** ``：空行，用于分隔逻辑块。
- **[行 4226]** `	if (   ( phi > (phi_0 + c1 * alpha * phi_0_prime) )`：条件分支：根据当前状态选择执行路径。
- **[行 4227]** `	    || ( i > 0 && phi > phi_old ) )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4228]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 4229]** `	    return line_search_zoom_strong_wolfe(phi_old, phi_prime_old, alpha_old,`：返回当前函数结果。
- **[行 4230]** `						 phi,     phi_prime,     alpha,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4231]** `						 phi_0,   phi_0_prime,   BFGS_p_vector,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4232]** `						 c1,      c2,            max_iter, solution_delta);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4233]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 4234]** ``：空行，用于分隔逻辑块。
- **[行 4235]** `	if (std::fabs(phi_prime) <= c2 * std::fabs(phi_0_prime))`：条件分支：根据当前状态选择执行路径。
- **[行 4236]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 4237]** `	    return alpha;`：返回当前函数结果。
- **[行 4238]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 4239]** ``：空行，用于分隔逻辑块。
- **[行 4240]** `	if (phi_prime >= 0)`：条件分支：根据当前状态选择执行路径。
- **[行 4241]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 4242]** `	    return line_search_zoom_strong_wolfe(phi,     phi_prime,     alpha,`：返回当前函数结果。
- **[行 4243]** `						 phi_old, phi_prime_old, alpha_old,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4244]** `						 phi_0,   phi_0_prime,   BFGS_p_vector,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4245]** `						 c1,      c2,            max_iter, solution_delta);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4246]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 4247]** ``：空行，用于分隔逻辑块。
- **[行 4248]** `	phi_old = phi;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4249]** `	phi_prime_old = phi_prime;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4250]** `	alpha_old = alpha;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4251]** ``：空行，用于分隔逻辑块。
- **[行 4252]** `	alpha = std::min(0.6*alpha, alpha_max);`：调用 C++ 标准库工具函数/容器接口。
- **[行 4253]** `      }`：作用域边界（代码块开始/结束）。
- **[行 4254]** ``：空行，用于分隔逻辑块。
- **[行 4255]** `    return alpha;`：返回当前函数结果。
- **[行 4256]** `  }`：作用域边界（代码块开始/结束）。
- **[行 4257]** ``：空行，用于分隔逻辑块。
- **[行 4258]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 4259]** `  double PhaseFieldMonolithicSolve<dim>::`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4260]** `    line_search_zoom_strong_wolfe(double phi_low, double phi_low_prime, double alpha_low,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4261]** `				  double phi_high, double phi_high_prime, double alpha_high,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4262]** `				  double phi_0, double phi_0_prime, const BlockVector<double> & BFGS_p_vector,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4263]** `				  double c1, double c2, unsigned int max_iter, const BlockVector<double> & solution_delta)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4264]** `  {`：作用域边界（代码块开始/结束）。
- **[行 4265]** `    double alpha = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4266]** `    std::pair<double, double> current_phi_phi_prime;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4267]** `    double phi, phi_prime;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4268]** ``：空行，用于分隔逻辑块。
- **[行 4269]** `    unsigned int i = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4270]** `    for (; i < max_iter; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4271]** `      {`：作用域边界（代码块开始/结束）。
- **[行 4272]** `	// a simple bisection is faster than cubic interpolation`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4273]** `	alpha = 0.5 * (alpha_low + alpha_high);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4274]** `	//alpha = line_search_interpolation_cubic(alpha_low, phi_low, phi_low_prime,`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4275]** `	//					alpha_high, phi_high, phi_high_prime);`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4276]** `	current_phi_phi_prime = calculate_phi_and_phi_prime(alpha, BFGS_p_vector, solution_delta);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4277]** `	phi = current_phi_phi_prime.first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4278]** `	phi_prime = current_phi_phi_prime.second;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4279]** ``：空行，用于分隔逻辑块。
- **[行 4280]** `	if (   (phi > phi_0 + c1 * alpha * phi_0_prime)`：条件分支：根据当前状态选择执行路径。
- **[行 4281]** `	    || (phi > phi_low) )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4282]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 4283]** `	    alpha_high = alpha;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4284]** `	    phi_high = phi;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4285]** `	    phi_high_prime = phi_prime;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4286]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 4287]** `	else`：条件分支的兜底路径。
- **[行 4288]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 4289]** `	    if (std::fabs(phi_prime) <= c2 * std::fabs(phi_0_prime))`：条件分支：根据当前状态选择执行路径。
- **[行 4290]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 4291]** `		return alpha;`：返回当前函数结果。
- **[行 4292]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 4293]** ``：空行，用于分隔逻辑块。
- **[行 4294]** `	    if (phi_prime * (alpha_high - alpha_low) >= 0.0)`：条件分支：根据当前状态选择执行路径。
- **[行 4295]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 4296]** `		alpha_high = alpha_low;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4297]** `		phi_high_prime = phi_low_prime;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4298]** `		phi_high = phi_low;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4299]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 4300]** ``：空行，用于分隔逻辑块。
- **[行 4301]** `	    alpha_low = alpha;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4302]** `	    phi_low_prime = phi_prime;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4303]** `	    phi_low = phi;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4304]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 4305]** `      }`：作用域边界（代码块开始/结束）。
- **[行 4306]** ``：空行，用于分隔逻辑块。
- **[行 4307]** `    // avoid unused variable warnings from compiler`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4308]** `    (void)phi_high;`：函数调用语句，触发对应计算或操作。
- **[行 4309]** `    (void)phi_high_prime;`：函数调用语句，触发对应计算或操作。
- **[行 4310]** `    return alpha;`：返回当前函数结果。
- **[行 4311]** `  }`：作用域边界（代码块开始/结束）。
- **[行 4312]** ``：空行，用于分隔逻辑块。
- **[行 4313]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 4314]** `  double PhaseFieldMonolithicSolve<dim>::`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4315]** `    line_search_interpolation_cubic(const double alpha_0, const double phi_0, const double phi_0_prime,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4316]** `  			            const double alpha_1, const double phi_1, const double phi_1_prime)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4317]** `  {`：作用域边界（代码块开始/结束）。
- **[行 4318]** `    const double d1 = phi_0_prime + phi_1_prime - 3.0 * (phi_0 - phi_1) / (alpha_0 - alpha_1);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4319]** ``：空行，用于分隔逻辑块。
- **[行 4320]** `    const double temp = d1 * d1 - phi_0_prime * phi_1_prime;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4321]** ``：空行，用于分隔逻辑块。
- **[行 4322]** `    if (temp < 0.0)`：条件分支：根据当前状态选择执行路径。
- **[行 4323]** `      return 0.5 * (alpha_0 + alpha_1);`：返回当前函数结果。
- **[行 4324]** ``：空行，用于分隔逻辑块。
- **[行 4325]** `    int sign;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4326]** `    if (alpha_1 > alpha_0)`：条件分支：根据当前状态选择执行路径。
- **[行 4327]** `      sign = 1;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4328]** `    else`：条件分支的兜底路径。
- **[行 4329]** `      sign = -1;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4330]** ``：空行，用于分隔逻辑块。
- **[行 4331]** `    const double d2 = sign * std::sqrt(temp);`：调用 C++ 标准库工具函数/容器接口。
- **[行 4332]** ``：空行，用于分隔逻辑块。
- **[行 4333]** `    const double alpha = alpha_1 - (alpha_1 - alpha_0)`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4334]** `	               * (phi_1_prime + d2 - d1) / (phi_1_prime - phi_0_prime + 2*d2);`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4335]** ``：空行，用于分隔逻辑块。
- **[行 4336]** `    if (    (alpha_1 > alpha_0)`：条件分支：根据当前状态选择执行路径。
- **[行 4337]** `	 && (alpha > alpha_1 || alpha < alpha_0))`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4338]** `      return 0.5 * (alpha_0 + alpha_1);`：返回当前函数结果。
- **[行 4339]** ``：空行，用于分隔逻辑块。
- **[行 4340]** `    if (    (alpha_0 > alpha_1)`：条件分支：根据当前状态选择执行路径。
- **[行 4341]** `	 && (alpha > alpha_0 || alpha < alpha_1))`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4342]** `      return 0.5 * (alpha_0 + alpha_1);`：返回当前函数结果。
- **[行 4343]** ``：空行，用于分隔逻辑块。
- **[行 4344]** `    return alpha;`：返回当前函数结果。
- **[行 4345]** `  }`：作用域边界（代码块开始/结束）。
- **[行 4346]** ``：空行，用于分隔逻辑块。
- **[行 4347]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 4348]** `  std::pair<double, double> PhaseFieldMonolithicSolve<dim>::`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4349]** `    calculate_phi_and_phi_prime(const double alpha,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4350]** `				const BlockVector<double> & BFGS_p_vector,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4351]** `				const BlockVector<double> & solution_delta)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4352]** `  {`：作用域边界（代码块开始/结束）。
- **[行 4353]** `    // the first component is phi(alpha), the second component is phi_prime(alpha),`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4354]** `    std::pair<double, double> phi_values;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4355]** ``：空行，用于分隔逻辑块。
- **[行 4356]** `    BlockVector<double> solution_delta_trial(solution_delta);`：函数调用语句，触发对应计算或操作。
- **[行 4357]** `    solution_delta_trial.add(alpha, BFGS_p_vector);`：函数调用语句，触发对应计算或操作。
- **[行 4358]** ``：空行，用于分隔逻辑块。
- **[行 4359]** `    update_qph_incremental(solution_delta_trial, m_solution);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4360]** ``：空行，用于分隔逻辑块。
- **[行 4361]** `    BlockVector<double> system_rhs(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 4362]** `    assemble_system_rhs_BFGS_parallel(m_solution, system_rhs);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4363]** `    //m_constraints.condense(system_rhs);`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4364]** ``：空行，用于分隔逻辑块。
- **[行 4365]** `    phi_values.first = calculate_energy_functional();`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4366]** `    phi_values.second = system_rhs * BFGS_p_vector;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4367]** `    return phi_values;`：返回当前函数结果。
- **[行 4368]** `  }`：作用域边界（代码块开始/结束）。
- **[行 4369]** ``：空行，用于分隔逻辑块。
- **[行 4370]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 4371]** `  void PhaseFieldMonolithicSolve<dim>::LBFGS_B0(BlockVector<double> & LBFGS_r_vector,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4372]** `						BlockVector<double> & LBFGS_q_vector)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4373]** `  {`：作用域边界（代码块开始/结束）。
- **[行 4374]** `    m_timer.enter_subsection("Solve B0");`：函数调用语句，触发对应计算或操作。
- **[行 4375]** ``：空行，用于分隔逻辑块。
- **[行 4376]** `    assemble_system_B0(m_solution);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4377]** ``：空行，用于分隔逻辑块。
- **[行 4378]** `    if (m_parameters.m_type_linear_solver == "Direct")`：条件分支：根据当前状态选择执行路径。
- **[行 4379]** `      {`：作用域边界（代码块开始/结束）。
- **[行 4380]** `	SparseDirectUMFPACK A_direct;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4381]** `	A_direct.initialize(m_tangent_matrix);`：函数调用语句，触发对应计算或操作。
- **[行 4382]** `	A_direct.vmult(LBFGS_r_vector,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4383]** `		       LBFGS_q_vector);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4384]** `      }`：作用域边界（代码块开始/结束）。
- **[行 4385]** `    else if (m_parameters.m_type_linear_solver == "CG")`：多分支条件判断，处理备选情形。
- **[行 4386]** `      {`：作用域边界（代码块开始/结束）。
- **[行 4387]** `/*`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4388]** `	SolverControl            solver_control(1e6, 1e-9);`：函数调用语句，触发对应计算或操作。
- **[行 4389]** `	SolverCG<BlockVector<double>> cg(solver_control);`：函数调用语句，触发对应计算或操作。
- **[行 4390]** ``：空行，用于分隔逻辑块。
- **[行 4391]** `	PreconditionJacobi<BlockSparseMatrix<double>> preconditioner;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4392]** `	preconditioner.initialize(m_tangent_matrix, 1.0);`：函数调用语句，触发对应计算或操作。
- **[行 4393]** ``：空行，用于分隔逻辑块。
- **[行 4394]** `	cg.solve(m_tangent_matrix,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4395]** `		 LBFGS_r_vector,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4396]** `		 LBFGS_q_vector,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4397]** `		 preconditioner);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4398]** `*/`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4399]** `	SolverControl            solver_control_uu(1e6, 1e-9);`：函数调用语句，触发对应计算或操作。
- **[行 4400]** `	SolverCG<Vector<double>> cg_uu(solver_control_uu);`：函数调用语句，触发对应计算或操作。
- **[行 4401]** ``：空行，用于分隔逻辑块。
- **[行 4402]** `	PreconditionJacobi<SparseMatrix<double>> preconditioner_uu;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4403]** `	preconditioner_uu.initialize(m_tangent_matrix.block(m_u_dof, m_u_dof), 1.0);`：函数调用语句，触发对应计算或操作。
- **[行 4404]** `	cg_uu.solve(m_tangent_matrix.block(m_u_dof, m_u_dof),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4405]** `	            LBFGS_r_vector.block(m_u_dof),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4406]** `	            LBFGS_q_vector.block(m_u_dof),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4407]** `	            preconditioner_uu);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4408]** ``：空行，用于分隔逻辑块。
- **[行 4409]** `	SolverControl            solver_control_dd(1e6, 1e-15);`：函数调用语句，触发对应计算或操作。
- **[行 4410]** `	SolverCG<Vector<double>> cg_dd(solver_control_dd);`：函数调用语句，触发对应计算或操作。
- **[行 4411]** ``：空行，用于分隔逻辑块。
- **[行 4412]** `	PreconditionJacobi<SparseMatrix<double>> preconditioner_dd;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4413]** `	preconditioner_dd.initialize(m_tangent_matrix.block(m_d_dof, m_d_dof), 1.0);`：函数调用语句，触发对应计算或操作。
- **[行 4414]** `	cg_dd.solve(m_tangent_matrix.block(m_d_dof, m_d_dof),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4415]** `	            LBFGS_r_vector.block(m_d_dof),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4416]** `	            LBFGS_q_vector.block(m_d_dof),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4417]** `	            preconditioner_dd);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4418]** `      }`：作用域边界（代码块开始/结束）。
- **[行 4419]** `    else`：条件分支的兜底路径。
- **[行 4420]** `      {`：作用域边界（代码块开始/结束）。
- **[行 4421]** `	AssertThrow(false,`：运行期断言/检查，验证输入与状态合法性。
- **[行 4422]** `	            ExcMessage("Selected linear solver not implemented!"));`：函数调用语句，触发对应计算或操作。
- **[行 4423]** `      }`：作用域边界（代码块开始/结束）。
- **[行 4424]** ``：空行，用于分隔逻辑块。
- **[行 4425]** `    m_timer.leave_subsection();`：函数调用语句，触发对应计算或操作。
- **[行 4426]** `  }`：作用域边界（代码块开始/结束）。
- **[行 4427]** ``：空行，用于分隔逻辑块。
- **[行 4428]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 4429]** `  void PhaseFieldMonolithicSolve<dim>::print_conv_header_LBFGS()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4430]** `  {`：作用域边界（代码块开始/结束）。
- **[行 4431]** `    m_logfile << "\t\t" << "L-BFGS (warning: without phasefield irreversibility)" << std::endl;;`：写日志输出，记录当前计算状态与结果。
- **[行 4432]** `    static const unsigned int l_width = 100;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4433]** `    m_logfile << '\t' << '\t';`：写日志输出，记录当前计算状态与结果。
- **[行 4434]** `    for (unsigned int i = 0; i < l_width; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4435]** `      m_logfile << '_';`：写日志输出，记录当前计算状态与结果。
- **[行 4436]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 4437]** ``：空行，用于分隔逻辑块。
- **[行 4438]** `    m_logfile << "\t\t itr "`：写日志输出，记录当前计算状态与结果。
- **[行 4439]** `              << " |  LS-alpha     Energy      Res_Norm    "`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4440]** `              << " Res_u      Res_d    Inc_Norm   "`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4441]** `              << " Inc_u      Inc_d" << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4442]** ``：空行，用于分隔逻辑块。
- **[行 4443]** `    m_logfile << '\t' << '\t';`：写日志输出，记录当前计算状态与结果。
- **[行 4444]** `    for (unsigned int i = 0; i < l_width; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4445]** `      m_logfile << '_';`：写日志输出，记录当前计算状态与结果。
- **[行 4446]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 4447]** `  }`：作用域边界（代码块开始/结束）。
- **[行 4448]** ``：空行，用于分隔逻辑块。
- **[行 4449]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 4450]** `  void PhaseFieldMonolithicSolve<dim>::print_conv_header_LBFGSB()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4451]** `  {`：作用域边界（代码块开始/结束）。
- **[行 4452]** `    m_logfile << '\t' << "L-BFGS-B" << std::endl;;`：写日志输出，记录当前计算状态与结果。
- **[行 4453]** `    static const unsigned int l_width = 130;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4454]** `    m_logfile << '\t';`：写日志输出，记录当前计算状态与结果。
- **[行 4455]** `    for (unsigned int i = 0; i < l_width; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4456]** `      m_logfile << '_';`：写日志输出，记录当前计算状态与结果。
- **[行 4457]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 4458]** ``：空行，用于分隔逻辑块。
- **[行 4459]** `    m_logfile << "\t itr  "`：写日志输出，记录当前计算状态与结果。
- **[行 4460]** `              << " |    LB    UB    LUB    CG-itr    LS-alpha     Energy      Res_Norm    "`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4461]** `              << " Res_u      Res_d    Inc_Norm   "`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4462]** `              << " Inc_u      Inc_d" << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4463]** ``：空行，用于分隔逻辑块。
- **[行 4464]** `    m_logfile << '\t';`：写日志输出，记录当前计算状态与结果。
- **[行 4465]** `    for (unsigned int i = 0; i < l_width; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4466]** `      m_logfile << '_';`：写日志输出，记录当前计算状态与结果。
- **[行 4467]** `    m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 4468]** `  }`：作用域边界（代码块开始/结束）。
- **[行 4469]** ``：空行，用于分隔逻辑块。
- **[行 4470]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 4471]** `  void PhaseFieldMonolithicSolve<dim>::`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4472]** `  solve_nonlinear_timestep_LBFGS(BlockVector<double> & solution_delta,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4473]** `				 BlockVector<double> & LBFGS_update_refine)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4474]** `  {`：作用域边界（代码块开始/结束）。
- **[行 4475]** `    BlockVector<double> LBFGS_update(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 4476]** ``：空行，用于分隔逻辑块。
- **[行 4477]** `    LBFGS_update = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4478]** ``：空行，用于分隔逻辑块。
- **[行 4479]** `    m_error_residual.reset();`：函数调用语句，触发对应计算或操作。
- **[行 4480]** `    m_error_residual_0.reset();`：函数调用语句，触发对应计算或操作。
- **[行 4481]** `    m_error_residual_norm.reset();`：函数调用语句，触发对应计算或操作。
- **[行 4482]** `    m_error_update.reset();`：函数调用语句，触发对应计算或操作。
- **[行 4483]** `    m_error_update_0.reset();`：函数调用语句，触发对应计算或操作。
- **[行 4484]** `    m_error_update_norm.reset();`：函数调用语句，触发对应计算或操作。
- **[行 4485]** ``：空行，用于分隔逻辑块。
- **[行 4486]** `    if (m_parameters.m_output_iteration_history)`：条件分支：根据当前状态选择执行路径。
- **[行 4487]** `      print_conv_header_LBFGS();`：函数调用语句，触发对应计算或操作。
- **[行 4488]** ``：空行，用于分隔逻辑块。
- **[行 4489]** `    unsigned int LBFGS_iteration = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4490]** ``：空行，用于分隔逻辑块。
- **[行 4491]** `    BlockVector<double> LBFGS_r_vector(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 4492]** `    BlockVector<double> LBFGS_y_vector(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 4493]** `    BlockVector<double> LBFGS_q_vector(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 4494]** `    BlockVector<double> LBFGS_s_vector(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 4495]** `    std::list<std::pair< std::pair<BlockVector<double>,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4496]** `                                   BlockVector<double>>,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4497]** `                         double>> LBFGS_vector_list;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4498]** ``：空行，用于分隔逻辑块。
- **[行 4499]** `    const unsigned int LBFGS_m = m_parameters.m_LBFGS_m;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4500]** `    std::list<double> LBFGS_alpha_list;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4501]** ``：空行，用于分隔逻辑块。
- **[行 4502]** `    double line_search_parameter = 0.0;`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4503]** `    double LBFGS_beta = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4504]** `    double rho = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4505]** ``：空行，用于分隔逻辑块。
- **[行 4506]** `    for (; LBFGS_iteration < m_parameters.m_max_iterations_BFGS; ++LBFGS_iteration)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4507]** `      {`：作用域边界（代码块开始/结束）。
- **[行 4508]** `	if (m_parameters.m_output_iteration_history)`：条件分支：根据当前状态选择执行路径。
- **[行 4509]** `	  m_logfile << '\t' << '\t' << std::setw(4) << LBFGS_iteration << ' '`：写日志输出，记录当前计算状态与结果。
- **[行 4510]** `                    << std::flush;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4511]** ``：空行，用于分隔逻辑块。
- **[行 4512]** `        make_constraints(LBFGS_iteration);`：函数调用语句，触发对应计算或操作。
- **[行 4513]** ``：空行，用于分隔逻辑块。
- **[行 4514]** `        // At the first step, we simply distribute the inhomogeneous part of`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4515]** `        // the constraints`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4516]** `        if (LBFGS_iteration == 0)`：条件分支：根据当前状态选择执行路径。
- **[行 4517]** `          {`：作用域边界（代码块开始/结束）。
- **[行 4518]** `            // use the solution from the previous solve on the`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4519]** `            // refined mesh as initial guess`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4520]** `            LBFGS_update = LBFGS_update_refine;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4521]** ``：空行，用于分隔逻辑块。
- **[行 4522]** `            m_constraints.distribute(LBFGS_update);`：函数调用语句，触发对应计算或操作。
- **[行 4523]** `            solution_delta += LBFGS_update;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4524]** ``：空行，用于分隔逻辑块。
- **[行 4525]** `            update_qph_incremental(solution_delta, m_solution);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4526]** `            if (m_parameters.m_output_iteration_history)`：条件分支：根据当前状态选择执行路径。
- **[行 4527]** `              {`：作用域边界（代码块开始/结束）。
- **[行 4528]** `                m_logfile << " | " << std::flush;`：写日志输出，记录当前计算状态与结果。
- **[行 4529]** `                m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 4530]** `              }`：作用域边界（代码块开始/结束）。
- **[行 4531]** `            continue;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4532]** `          }`：作用域边界（代码块开始/结束）。
- **[行 4533]** `        else if (LBFGS_iteration == 1)`：多分支条件判断，处理备选情形。
- **[行 4534]** `          {`：作用域边界（代码块开始/结束）。
- **[行 4535]** `	    // Calculate the residual vector r. NOTICE that in the context of`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4536]** `	    // BFGS, this r is the gradient of the energy functional (objective function),`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4537]** `	    // NOT the negative gradient of the energy functional`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4538]** `	    assemble_system_rhs_BFGS_parallel(m_solution, m_system_rhs);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4539]** ``：空行，用于分隔逻辑块。
- **[行 4540]** `	    // We cannot simply zero out the dofs that are constrained, since we might`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4541]** `	    // have hanging node constraints. In this case, we need to modify the RHS`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4542]** `	    // as C^T * b, which C contains entries of 0.5 (x_3 = 0.5*x_1 + 0.5*x_2)`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4543]** `	    //for (unsigned int i = 0; i < m_dof_handler.n_dofs(); ++i)`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4544]** `	      //if (m_constraints.is_constrained(i))`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4545]** `		//m_system_rhs(i) = 0.0;`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4546]** ``：空行，用于分隔逻辑块。
- **[行 4547]** `	    // if m_constraints has inhomogeneity, we cannot call m_constraints.condense(m_system_rhs),`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4548]** `	    // since the m_system_matrix needs to be provided to modify the RHS properly. However, this`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4549]** `	    // error will not be detected in the release mode and only will be detected on the debug mode`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4550]** `	    // if we use assemble_system_rhs_BFGS_parallel, then condense() is not necessary`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4551]** `	    //m_constraints.condense(m_system_rhs);`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4552]** `          }`：作用域边界（代码块开始/结束）。
- **[行 4553]** ``：空行，用于分隔逻辑块。
- **[行 4554]** `        get_error_residual(m_error_residual);`：函数调用语句，触发对应计算或操作。
- **[行 4555]** `        if (LBFGS_iteration == 1)`：条件分支：根据当前状态选择执行路径。
- **[行 4556]** `          m_error_residual_0 = m_error_residual;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4557]** ``：空行，用于分隔逻辑块。
- **[行 4558]** `        m_error_residual_norm = m_error_residual;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4559]** `        // For three-point bending problem and 3D problem, we use absolute residual`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4560]** `        // for convergence test`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4561]** `        if (m_parameters.m_relative_residual)`：条件分支：根据当前状态选择执行路径。
- **[行 4562]** `          m_error_residual_norm.normalize(m_error_residual_0);`：函数调用语句，触发对应计算或操作。
- **[行 4563]** ``：空行，用于分隔逻辑块。
- **[行 4564]** `        if (LBFGS_iteration > 1 && m_error_update_norm.m_u <= m_parameters.m_tol_u_incr`：条件分支：根据当前状态选择执行路径。
- **[行 4565]** `                                && m_error_residual_norm.m_u <= m_parameters.m_tol_u_residual`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4566]** `			        && m_error_update_norm.m_d <= m_parameters.m_tol_d_incr`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4567]** `			        && m_error_residual_norm.m_d <= m_parameters.m_tol_d_residual`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4568]** `				)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4569]** `          {`：作用域边界（代码块开始/结束）。
- **[行 4570]** `            if (m_parameters.m_output_iteration_history)`：条件分支：根据当前状态选择执行路径。
- **[行 4571]** `              {`：作用域边界（代码块开始/结束）。
- **[行 4572]** `		m_logfile << " | ";`：写日志输出，记录当前计算状态与结果。
- **[行 4573]** `		m_logfile << " CONVERGED! " << std::fixed << std::setprecision(3) << std::setw(7)`：写日志输出，记录当前计算状态与结果。
- **[行 4574]** `			  << std::scientific`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4575]** `		      << "    ----    "`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4576]** `		      << "  " << m_error_residual_norm.m_norm`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4577]** `		      << "  " << m_error_residual_norm.m_u`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4578]** `		      << "  " << m_error_residual_norm.m_d`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4579]** `		      << "  " << m_error_update_norm.m_norm`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4580]** `		      << "  " << m_error_update_norm.m_u`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4581]** `		      << "  " << m_error_update_norm.m_d`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4582]** `		      << "  " << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4583]** ``：空行，用于分隔逻辑块。
- **[行 4584]** `		m_logfile << '\t' << '\t';`：写日志输出，记录当前计算状态与结果。
- **[行 4585]** `		for (unsigned int i = 0; i < 100; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4586]** `		  m_logfile << '_';`：写日志输出，记录当前计算状态与结果。
- **[行 4587]** `		m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 4588]** `              }`：作用域边界（代码块开始/结束）。
- **[行 4589]** ``：空行，用于分隔逻辑块。
- **[行 4590]** `            m_logfile << "\t\tConvergence is reached after "`：写日志输出，记录当前计算状态与结果。
- **[行 4591]** `        	      << LBFGS_iteration << " L-BFGS iterations."<< std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4592]** ``：空行，用于分隔逻辑块。
- **[行 4593]** `            m_logfile << "\t\tResidual information of convergence:" << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 4594]** ``：空行，用于分隔逻辑块。
- **[行 4595]** `            if (m_parameters.m_relative_residual)`：条件分支：根据当前状态选择执行路径。
- **[行 4596]** `              {`：作用域边界（代码块开始/结束）。
- **[行 4597]** `		m_logfile << "\t\t\tRelative residual of disp. equation: "`：写日志输出，记录当前计算状态与结果。
- **[行 4598]** `			  << m_error_residual_norm.m_u << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4599]** ``：空行，用于分隔逻辑块。
- **[行 4600]** `		m_logfile << "\t\t\tAbsolute residual of disp. equation: "`：写日志输出，记录当前计算状态与结果。
- **[行 4601]** `			  << m_error_residual_norm.m_u * m_error_residual_0.m_u << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4602]** ``：空行，用于分隔逻辑块。
- **[行 4603]** `		m_logfile << "\t\t\tRelative residual of phasefield equation: "`：写日志输出，记录当前计算状态与结果。
- **[行 4604]** `			  << m_error_residual_norm.m_d << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4605]** ``：空行，用于分隔逻辑块。
- **[行 4606]** `		m_logfile << "\t\t\tAbsolute residual of phasefield equation: "`：写日志输出，记录当前计算状态与结果。
- **[行 4607]** `			  << m_error_residual_norm.m_d * m_error_residual_0.m_d << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4608]** ``：空行，用于分隔逻辑块。
- **[行 4609]** `		m_logfile << "\t\t\tRelative increment of disp.: "`：写日志输出，记录当前计算状态与结果。
- **[行 4610]** `			  << m_error_update_norm.m_u << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4611]** ``：空行，用于分隔逻辑块。
- **[行 4612]** `		m_logfile << "\t\t\tAbsolute increment of disp.: "`：写日志输出，记录当前计算状态与结果。
- **[行 4613]** `			  << m_error_update_norm.m_u * m_error_update_0.m_u << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4614]** ``：空行，用于分隔逻辑块。
- **[行 4615]** `		m_logfile << "\t\t\tRelative increment of phasefield: "`：写日志输出，记录当前计算状态与结果。
- **[行 4616]** `			  << m_error_update_norm.m_d << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4617]** ``：空行，用于分隔逻辑块。
- **[行 4618]** `		m_logfile << "\t\t\tAbsolute increment of phasefield: "`：写日志输出，记录当前计算状态与结果。
- **[行 4619]** `			  << m_error_update_norm.m_d * m_error_update_0.m_d << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4620]** `              }`：作用域边界（代码块开始/结束）。
- **[行 4621]** `            else`：条件分支的兜底路径。
- **[行 4622]** `              {`：作用域边界（代码块开始/结束）。
- **[行 4623]** `		m_logfile << "\t\t\tAbsolute residual of disp. equation: "`：写日志输出，记录当前计算状态与结果。
- **[行 4624]** `			  << m_error_residual_norm.m_u << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4625]** ``：空行，用于分隔逻辑块。
- **[行 4626]** `		m_logfile << "\t\t\tAbsolute residual of phasefield equation: "`：写日志输出，记录当前计算状态与结果。
- **[行 4627]** `			  << m_error_residual_norm.m_d << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4628]** ``：空行，用于分隔逻辑块。
- **[行 4629]** `		m_logfile << "\t\t\tAbsolute increment of disp.: "`：写日志输出，记录当前计算状态与结果。
- **[行 4630]** `			  << m_error_update_norm.m_u << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4631]** ``：空行，用于分隔逻辑块。
- **[行 4632]** `		m_logfile << "\t\t\tAbsolute increment of phasefield: "`：写日志输出，记录当前计算状态与结果。
- **[行 4633]** `			  << m_error_update_norm.m_d << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4634]** `              }`：作用域边界（代码块开始/结束）。
- **[行 4635]** ``：空行，用于分隔逻辑块。
- **[行 4636]** `            break;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4637]** `          }`：作用域边界（代码块开始/结束）。
- **[行 4638]** ``：空行，用于分隔逻辑块。
- **[行 4639]** `        // LBFGS algorithm`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4640]** `        LBFGS_q_vector = m_system_rhs;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4641]** ``：空行，用于分隔逻辑块。
- **[行 4642]** `        LBFGS_alpha_list.clear();`：函数调用语句，触发对应计算或操作。
- **[行 4643]** `        for (auto itr = LBFGS_vector_list.begin(); itr != LBFGS_vector_list.end(); ++itr)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4644]** `          {`：作用域边界（代码块开始/结束）。
- **[行 4645]** `            LBFGS_s_vector = (itr->first).first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4646]** `            LBFGS_y_vector = (itr->first).second;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4647]** `            rho = itr->second;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4648]** ``：空行，用于分隔逻辑块。
- **[行 4649]** `            const double alpha = rho * (LBFGS_s_vector * LBFGS_q_vector);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4650]** `            LBFGS_alpha_list.push_back(alpha);`：函数调用语句，触发对应计算或操作。
- **[行 4651]** ``：空行，用于分隔逻辑块。
- **[行 4652]** `            LBFGS_q_vector.add(-alpha, LBFGS_y_vector);`：函数调用语句，触发对应计算或操作。
- **[行 4653]** `          }`：作用域边界（代码块开始/结束）。
- **[行 4654]** `/*`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4655]** `        double scale_gamma = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4656]** `        if (LBFGS_iteration == 1)`：条件分支：根据当前状态选择执行路径。
- **[行 4657]** `          {`：作用域边界（代码块开始/结束）。
- **[行 4658]** `            scale_gamma = 1.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4659]** `          }`：作用域边界（代码块开始/结束）。
- **[行 4660]** `        else`：条件分支的兜底路径。
- **[行 4661]** `          {`：作用域边界（代码块开始/结束）。
- **[行 4662]** `            LBFGS_s_vector = LBFGS_vector_list.front().first.first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4663]** `            LBFGS_y_vector = LBFGS_vector_list.front().first.second;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4664]** `            scale_gamma = (LBFGS_s_vector * LBFGS_y_vector)/(LBFGS_y_vector * LBFGS_y_vector);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4665]** `          }`：作用域边界（代码块开始/结束）。
- **[行 4666]** ``：空行，用于分隔逻辑块。
- **[行 4667]** `        LBFGS_q_vector *= scale_gamma;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4668]** `        LBFGS_r_vector = LBFGS_q_vector;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4669]** `*/`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4670]** `        LBFGS_B0(LBFGS_r_vector,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4671]** `		 LBFGS_q_vector);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4672]** ``：空行，用于分隔逻辑块。
- **[行 4673]** `        for (auto itr = LBFGS_vector_list.rbegin(); itr != LBFGS_vector_list.rend(); ++itr)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4674]** `          {`：作用域边界（代码块开始/结束）。
- **[行 4675]** `            LBFGS_s_vector = (itr->first).first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4676]** `            LBFGS_y_vector = (itr->first).second;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4677]** `            rho = itr->second;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4678]** ``：空行，用于分隔逻辑块。
- **[行 4679]** `            LBFGS_beta = rho * (LBFGS_y_vector * LBFGS_r_vector);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4680]** ``：空行，用于分隔逻辑块。
- **[行 4681]** `            const double alpha = LBFGS_alpha_list.back();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4682]** `            LBFGS_alpha_list.pop_back();`：函数调用语句，触发对应计算或操作。
- **[行 4683]** ``：空行，用于分隔逻辑块。
- **[行 4684]** `            LBFGS_r_vector.add(alpha - LBFGS_beta, LBFGS_s_vector);`：函数调用语句，触发对应计算或操作。
- **[行 4685]** `          }`：作用域边界（代码块开始/结束）。
- **[行 4686]** ``：空行，用于分隔逻辑块。
- **[行 4687]** `        LBFGS_r_vector *= -1.0; // this is the p_vector (search direction)`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4688]** ``：空行，用于分隔逻辑块。
- **[行 4689]** `        m_constraints.distribute(LBFGS_r_vector);`：函数调用语句，触发对应计算或操作。
- **[行 4690]** ``：空行，用于分隔逻辑块。
- **[行 4691]** `        // We need a line search algorithm to decide line_search_parameter`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4692]** ``：空行，用于分隔逻辑块。
- **[行 4693]** `        if(m_parameters.m_type_line_search == "StrongWolfe")`：条件分支：根据当前状态选择执行路径。
- **[行 4694]** `          {`：作用域边界（代码块开始/结束）。
- **[行 4695]** `            const double phi_0 = calculate_energy_functional();`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4696]** `            const double phi_0_prime = m_system_rhs * LBFGS_r_vector;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4697]** ``：空行，用于分隔逻辑块。
- **[行 4698]** `            line_search_parameter = line_search_stepsize_strong_wolfe(phi_0,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4699]** `		    				                      phi_0_prime,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4700]** `								      LBFGS_r_vector,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4701]** `						                      solution_delta);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4702]** `          }`：作用域边界（代码块开始/结束）。
- **[行 4703]** `        else if(m_parameters.m_type_line_search == "GradientBased")`：多分支条件判断，处理备选情形。
- **[行 4704]** `          {`：作用域边界（代码块开始/结束）。
- **[行 4705]** `	    // LBFGS_r_vector is the search direction`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4706]** `	    line_search_parameter = line_search_stepsize_gradient_based(LBFGS_r_vector,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4707]** `									solution_delta);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4708]** `          }`：作用域边界（代码块开始/结束）。
- **[行 4709]** `        else`：条件分支的兜底路径。
- **[行 4710]** `          {`：作用域边界（代码块开始/结束）。
- **[行 4711]** `            Assert(false, ExcMessage("An unknown line search method is called!"));`：运行期断言/检查，验证输入与状态合法性。
- **[行 4712]** `          }`：作用域边界（代码块开始/结束）。
- **[行 4713]** ``：空行，用于分隔逻辑块。
- **[行 4714]** `        LBFGS_r_vector *= line_search_parameter;`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4715]** `        LBFGS_update = LBFGS_r_vector;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4716]** ``：空行，用于分隔逻辑块。
- **[行 4717]** `        get_error_update(LBFGS_update, m_error_update);`：函数调用语句，触发对应计算或操作。
- **[行 4718]** `        if (LBFGS_iteration == 1)`：条件分支：根据当前状态选择执行路径。
- **[行 4719]** `          m_error_update_0 = m_error_update;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4720]** ``：空行，用于分隔逻辑块。
- **[行 4721]** `        m_error_update_norm = m_error_update;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4722]** `        // For three-point bending problem and the sphere inclusion problem,`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4723]** `        // we use absolute residual for convergence test`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4724]** `        if (m_parameters.m_relative_residual)`：条件分支：根据当前状态选择执行路径。
- **[行 4725]** `          m_error_update_norm.normalize(m_error_update_0);`：函数调用语句，触发对应计算或操作。
- **[行 4726]** ``：空行，用于分隔逻辑块。
- **[行 4727]** `        solution_delta += LBFGS_update;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4728]** `        update_qph_incremental(solution_delta, m_solution);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4729]** ``：空行，用于分隔逻辑块。
- **[行 4730]** `        LBFGS_y_vector = m_system_rhs;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4731]** `        LBFGS_y_vector *= -1.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4732]** `        assemble_system_rhs_BFGS_parallel(m_solution, m_system_rhs);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4733]** `        // if we use assemble_system_rhs_BFGS_parallel, then condense() is not necessary`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4734]** `        //m_constraints.condense(m_system_rhs);`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4735]** `        LBFGS_y_vector += m_system_rhs;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4736]** ``：空行，用于分隔逻辑块。
- **[行 4737]** `        LBFGS_s_vector = LBFGS_update;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4738]** ``：空行，用于分隔逻辑块。
- **[行 4739]** `        if (LBFGS_iteration > LBFGS_m)`：条件分支：根据当前状态选择执行路径。
- **[行 4740]** `          LBFGS_vector_list.pop_back();`：函数调用语句，触发对应计算或操作。
- **[行 4741]** ``：空行，用于分隔逻辑块。
- **[行 4742]** `        rho = 1.0 / (LBFGS_y_vector * LBFGS_s_vector);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4743]** ``：空行，用于分隔逻辑块。
- **[行 4744]** `        LBFGS_vector_list.push_front(std::make_pair(std::make_pair(LBFGS_s_vector,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4745]** `								   LBFGS_y_vector),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4746]** `						    rho));`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4747]** `        if (m_parameters.m_output_iteration_history)`：条件分支：根据当前状态选择执行路径。
- **[行 4748]** `          {`：作用域边界（代码块开始/结束）。
- **[行 4749]** `	    const double energy_functional = calculate_energy_functional();`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4750]** ``：空行，用于分隔逻辑块。
- **[行 4751]** `	    m_logfile << " | " << std::fixed << std::setprecision(3) << std::setw(1)`：写日志输出，记录当前计算状态与结果。
- **[行 4752]** `		      << std::scientific`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4753]** `		      << "" << line_search_parameter`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4754]** `		      << std::fixed << std::setprecision(6) << std::setw(1)`：调用 C++ 标准库工具函数/容器接口。
- **[行 4755]** `					<< std::scientific`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4756]** `		      << "  " << energy_functional`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4757]** `		      << std::fixed << std::setprecision(3) << std::setw(1)`：调用 C++ 标准库工具函数/容器接口。
- **[行 4758]** `					<< std::scientific`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4759]** `		      << "  " << m_error_residual_norm.m_norm`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4760]** `		      << "  " << m_error_residual_norm.m_u`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4761]** `		      << "  " << m_error_residual_norm.m_d`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4762]** `		      << "  " << m_error_update_norm.m_norm`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4763]** `		      << "  " << m_error_update_norm.m_u`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4764]** `		      << "  " << m_error_update_norm.m_d`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4765]** `		      << "  " << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4766]** `          }`：作用域边界（代码块开始/结束）。
- **[行 4767]** `      }`：作用域边界（代码块开始/结束）。
- **[行 4768]** ``：空行，用于分隔逻辑块。
- **[行 4769]** `    AssertThrow(LBFGS_iteration < m_parameters.m_max_iterations_BFGS,`：运行期断言/检查，验证输入与状态合法性。
- **[行 4770]** `                ExcMessage("No convergence in L-BFGS nonlinear solver!"));`：函数调用语句，触发对应计算或操作。
- **[行 4771]** `  }`：作用域边界（代码块开始/结束）。
- **[行 4772]** ``：空行，用于分隔逻辑块。
- **[行 4773]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 4774]** `  void PhaseFieldMonolithicSolve<dim>::`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4775]** `  calculate_cauchy_point(const BlockSparseMatrix<double> & B0_matrix,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4776]** `	                 const std::list<BlockVector<double>> & y_vector_list,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4777]** `		         const std::list<BlockVector<double>> & b0xs_vector_list,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4778]** `			 const FullMatrix<double> & M_matrix,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4779]** `			 const BlockVector<double> & gradient_g,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4780]** `			 const BlockVector<double> & solution_delta,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4781]** `			 BlockVector<double> & solution_delta_cauchy_point)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4782]** `  {`：作用域边界（代码块开始/结束）。
- **[行 4783]** `    m_timer.enter_subsection("Calculate Cauchy point");`：函数调用语句，触发对应计算或操作。
- **[行 4784]** ``：空行，用于分隔逻辑块。
- **[行 4785]** `    solution_delta_cauchy_point = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4786]** `    BlockVector<double> gradient_d(gradient_g);`：函数调用语句，触发对应计算或操作。
- **[行 4787]** `    gradient_d *= -1;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4788]** ``：空行，用于分隔逻辑块。
- **[行 4789]** `    const unsigned int list_size = y_vector_list.size();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4790]** `    const auto itr_y_begin    = y_vector_list.begin();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4791]** `    const auto itr_b0xs_begin = b0xs_vector_list.begin();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4792]** ``：空行，用于分隔逻辑块。
- **[行 4793]** `    // t_series only contains t > 0`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4794]** `    std::priority_queue< std::pair<double, unsigned int>,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4795]** `                         std::vector<std::pair<double, unsigned int>>,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4796]** `        		 std::greater<std::pair<double, unsigned int>> >`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4797]** `    t_series = calculate_break_points(solution_delta,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4798]** `    			              gradient_g,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4799]** `				      gradient_d);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4800]** ``：空行，用于分隔逻辑块。
- **[行 4801]** `    // m_active_set_phasefield contains 1 or 2 for active set and 0 for inactive set`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4802]** `    for (unsigned int i = 0; i < m_active_set_phasefield.size(); ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4803]** `      {`：作用域边界（代码块开始/结束）。
- **[行 4804]** `	if (m_active_set_phasefield(i) > 0.5)`：条件分支：根据当前状态选择执行路径。
- **[行 4805]** `	  solution_delta_cauchy_point.block(m_d_dof)[i] = solution_delta.block(m_d_dof)[i];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4806]** `      }`：作用域边界（代码块开始/结束）。
- **[行 4807]** ``：空行，用于分隔逻辑块。
- **[行 4808]** `    // p = W^T * d`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4809]** `    Vector<double> p(2 * list_size);`：函数调用语句，触发对应计算或操作。
- **[行 4810]** `    for (unsigned int i = 0; i < list_size; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4811]** `      {`：作用域边界（代码块开始/结束）。
- **[行 4812]** `        p(i)             = (*std::next(itr_y_begin,    i)) * gradient_d;`：调用 C++ 标准库工具函数/容器接口。
- **[行 4813]** `        p(i + list_size) = (*std::next(itr_b0xs_begin, i)) * gradient_d;`：调用 C++ 标准库工具函数/容器接口。
- **[行 4814]** `      }`：作用域边界（代码块开始/结束）。
- **[行 4815]** ``：空行，用于分隔逻辑块。
- **[行 4816]** `    Vector<double> c(2 * list_size);`：函数调用语句，触发对应计算或操作。
- **[行 4817]** `    c = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4818]** ``：空行，用于分隔逻辑块。
- **[行 4819]** `    double f_prime = -(gradient_d * gradient_d);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4820]** ``：空行，用于分隔逻辑块。
- **[行 4821]** `    // M * p`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4822]** `    Vector<double> Mp(2 * list_size);`：函数调用语句，触发对应计算或操作。
- **[行 4823]** `    if (list_size > 0)`：条件分支：根据当前状态选择执行路径。
- **[行 4824]** `      M_matrix.vmult(Mp, p);`：函数调用语句，触发对应计算或操作。
- **[行 4825]** ``：空行，用于分隔逻辑块。
- **[行 4826]** `    // B_0 * d`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4827]** `    BlockVector<double> B0_grandient_d(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 4828]** `    B0_matrix.vmult(B0_grandient_d, gradient_d);`：函数调用语句，触发对应计算或操作。
- **[行 4829]** ``：空行，用于分隔逻辑块。
- **[行 4830]** `    double f_prime_prime = gradient_d * B0_grandient_d;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4831]** `    if (list_size > 0)`：条件分支：根据当前状态选择执行路径。
- **[行 4832]** `      f_prime_prime -= (p * Mp);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4833]** ``：空行，用于分隔逻辑块。
- **[行 4834]** `    double delta_t_min = -f_prime / f_prime_prime;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4835]** ``：空行，用于分隔逻辑块。
- **[行 4836]** `    double t_old = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4837]** ``：空行，用于分隔逻辑块。
- **[行 4838]** `    std::pair<double, unsigned int> top_pair = t_series.top();`：调用 C++ 标准库工具函数/容器接口。
- **[行 4839]** `    double t = top_pair.first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4840]** `    unsigned int b = top_pair.second;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4841]** ``：空行，用于分隔逻辑块。
- **[行 4842]** `    double delta_t = t - t_old;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4843]** ``：空行，用于分隔逻辑块。
- **[行 4844]** `    BlockVector<double> z(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 4845]** `    z = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4846]** ``：空行，用于分隔逻辑块。
- **[行 4847]** `    // w_b = W^T * e_b`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4848]** `    Vector<double> w_b(2 * list_size);`：函数调用语句，触发对应计算或操作。
- **[行 4849]** ``：空行，用于分隔逻辑块。
- **[行 4850]** `    // w_b_T_x_M = w_b^T * M`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4851]** `    Vector<double> w_b_T_x_M(2 * list_size);`：函数调用语句，触发对应计算或操作。
- **[行 4852]** ``：空行，用于分隔逻辑块。
- **[行 4853]** `    Vector<double> temp_vector(m_dofs_per_block[m_d_dof]);`：函数调用语句，触发对应计算或操作。
- **[行 4854]** ``：空行，用于分隔逻辑块。
- **[行 4855]** `    while (delta_t_min >= delta_t)`：循环推进（通常用于时间步或迭代收敛过程）。
- **[行 4856]** `      {`：作用域边界（代码块开始/结束）。
- **[行 4857]** `	t_series.pop();`：函数调用语句，触发对应计算或操作。
- **[行 4858]** ``：空行，用于分隔逻辑块。
- **[行 4859]** `	if (gradient_d.block(m_d_dof)[b] > 0)`：条件分支：根据当前状态选择执行路径。
- **[行 4860]** `	  solution_delta_cauchy_point.block(m_d_dof)[b] = 1.0 - m_solution.block(m_d_dof)[b];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4861]** `	else if (gradient_d.block(m_d_dof)[b] < 0)`：多分支条件判断，处理备选情形。
- **[行 4862]** `	  solution_delta_cauchy_point.block(m_d_dof)[b] = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4863]** `	else`：条件分支的兜底路径。
- **[行 4864]** `	  AssertThrow(false,`：运行期断言/检查，验证输入与状态合法性。
- **[行 4865]** `	              ExcMessage("gradient_d(b) cannot be zero!"));`：函数调用语句，触发对应计算或操作。
- **[行 4866]** ``：空行，用于分隔逻辑块。
- **[行 4867]** `	if (gradient_d.block(m_d_dof)[b] < 0)`：条件分支：根据当前状态选择执行路径。
- **[行 4868]** `	  m_active_set_phasefield[b] = 1; //lower bound`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4869]** `	else`：条件分支的兜底路径。
- **[行 4870]** `	  m_active_set_phasefield[b] = 2; //upper bound`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4871]** ``：空行，用于分隔逻辑块。
- **[行 4872]** `        // z_{j} = z_{j-1} + delta_t_{j-1} * gradient_d_{j-1};`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4873]** `	z.sadd(1.0, delta_t, gradient_d);`：函数调用语句，触发对应计算或操作。
- **[行 4874]** ``：空行，用于分隔逻辑块。
- **[行 4875]** `	// c_{j} = c_{j-1} + delta_t_{j-1} * p_{j-1};`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4876]** `	if (list_size > 0)`：条件分支：根据当前状态选择执行路径。
- **[行 4877]** `	  c.sadd(1.0, delta_t, p);`：函数调用语句，触发对应计算或操作。
- **[行 4878]** ``：空行，用于分隔逻辑块。
- **[行 4879]** `        double temp_scalar = ebT_x_B0_x_v(b, B0_matrix, z);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4880]** ``：空行，用于分隔逻辑块。
- **[行 4881]** `        // w_b = W^T * e_b`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4882]** `        for (unsigned int i = 0; i < list_size; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4883]** `          {`：作用域边界（代码块开始/结束）。
- **[行 4884]** `            w_b(i)             = (*std::next(itr_y_begin,    i)).block(m_d_dof)[b];`：调用 C++ 标准库工具函数/容器接口。
- **[行 4885]** `            w_b(i + list_size) = (*std::next(itr_b0xs_begin, i)).block(m_d_dof)[b];`：调用 C++ 标准库工具函数/容器接口。
- **[行 4886]** `          }`：作用域边界（代码块开始/结束）。
- **[行 4887]** ``：空行，用于分隔逻辑块。
- **[行 4888]** `        if (list_size > 0)`：条件分支：根据当前状态选择执行路径。
- **[行 4889]** `          M_matrix.vmult(w_b_T_x_M, w_b);`：函数调用语句，触发对应计算或操作。
- **[行 4890]** ``：空行，用于分隔逻辑块。
- **[行 4891]** `	f_prime += delta_t * f_prime_prime`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4892]** `	         + gradient_g.block(m_d_dof)[b] * gradient_g.block(m_d_dof)[b]`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4893]** `	         + temp_scalar * gradient_g.block(m_d_dof)[b];`：函数调用语句，触发对应计算或操作。
- **[行 4894]** ``：空行，用于分隔逻辑块。
- **[行 4895]** `	if (list_size > 0)`：条件分支：根据当前状态选择执行路径。
- **[行 4896]** `	  f_prime -= (w_b_T_x_M * c) * gradient_g.block(m_d_dof)[b];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4897]** ``：空行，用于分隔逻辑块。
- **[行 4898]** `	temp_scalar = ebT_x_B0_x_v(b, B0_matrix, gradient_d);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4899]** ``：空行，用于分隔逻辑块。
- **[行 4900]** `	f_prime_prime += 2.0 * gradient_g.block(m_d_dof)[b] * temp_scalar`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4901]** `	               + gradient_g.block(m_d_dof)[b] * gradient_g.block(m_d_dof)[b]`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4902]** `		             * B0_matrix.block(m_d_dof, m_d_dof)(b, b);`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4903]** ``：空行，用于分隔逻辑块。
- **[行 4904]** `	if (list_size > 0)`：条件分支：根据当前状态选择执行路径。
- **[行 4905]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 4906]** `	    f_prime_prime -= 2.0 * gradient_g.block(m_d_dof)[b] * (w_b_T_x_M * p);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4907]** `	    f_prime_prime -= gradient_g.block(m_d_dof)[b] * gradient_g.block(m_d_dof)[b] * (w_b_T_x_M * w_b);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4908]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 4909]** ``：空行，用于分隔逻辑块。
- **[行 4910]** `	// p_{j} = p_{j-1} + g_b * w_b;`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4911]** `	if (list_size > 0)`：条件分支：根据当前状态选择执行路径。
- **[行 4912]** `	  p.sadd(1.0, gradient_g.block(m_d_dof)[b], w_b);`：函数调用语句，触发对应计算或操作。
- **[行 4913]** ``：空行，用于分隔逻辑块。
- **[行 4914]** `	gradient_d.block(m_d_dof)[b] = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4915]** ``：空行，用于分隔逻辑块。
- **[行 4916]** `	delta_t_min = -f_prime / f_prime_prime;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4917]** ``：空行，用于分隔逻辑块。
- **[行 4918]** `	t_old = t;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4919]** ``：空行，用于分隔逻辑块。
- **[行 4920]** `	top_pair = t_series.top();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4921]** `	t = top_pair.first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4922]** `	b = top_pair.second;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4923]** ``：空行，用于分隔逻辑块。
- **[行 4924]** `	delta_t = t - t_old;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4925]** `      }`：作用域边界（代码块开始/结束）。
- **[行 4926]** ``：空行，用于分隔逻辑块。
- **[行 4927]** `    if (delta_t_min < 0)`：条件分支：根据当前状态选择执行路径。
- **[行 4928]** `      delta_t_min = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4929]** ``：空行，用于分隔逻辑块。
- **[行 4930]** `    t_old += delta_t_min;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4931]** ``：空行，用于分隔逻辑块。
- **[行 4932]** `    for(unsigned int i = 0; i < m_active_set_phasefield.size(); ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4933]** `      {`：作用域边界（代码块开始/结束）。
- **[行 4934]** `	// inactive phasefield dof`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4935]** `	if (m_active_set_phasefield(i) < 0.5)`：条件分支：根据当前状态选择执行路径。
- **[行 4936]** `	  solution_delta_cauchy_point.block(m_d_dof)[i] = solution_delta.block(m_d_dof)[i]`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4937]** `						+ t_old * gradient_d.block(m_d_dof)[i];`：函数调用语句，触发对应计算或操作。
- **[行 4938]** `      }`：作用域边界（代码块开始/结束）。
- **[行 4939]** ``：空行，用于分隔逻辑块。
- **[行 4940]** `    // There are no active constraints in the displacement field`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4941]** `    solution_delta_cauchy_point.block(m_u_dof) = solution_delta.block(m_u_dof);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4942]** `    (solution_delta_cauchy_point.block(m_u_dof)).add(t_old, gradient_d.block(m_u_dof));`：函数调用语句，触发对应计算或操作。
- **[行 4943]** ``：空行，用于分隔逻辑块。
- **[行 4944]** `    // We need to make sure the solution_delta_cauchy_point satisfies the essential`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4945]** `    // boundary conditions and the hanging-node constraints`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4946]** `    m_constraints.distribute(solution_delta_cauchy_point);`：函数调用语句，触发对应计算或操作。
- **[行 4947]** ``：空行，用于分隔逻辑块。
- **[行 4948]** `    m_timer.leave_subsection();`：函数调用语句，触发对应计算或操作。
- **[行 4949]** `  }`：作用域边界（代码块开始/结束）。
- **[行 4950]** ``：空行，用于分隔逻辑块。
- **[行 4951]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 4952]** `  void PhaseFieldMonolithicSolve<dim>::`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4953]** `  solve_nonlinear_timestep_LBFGS_B(BlockVector<double> & solution_delta,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4954]** `				   BlockVector<double> & LBFGS_update_refine)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4955]** `  {`：作用域边界（代码块开始/结束）。
- **[行 4956]** `    BlockVector<double> LBFGS_update(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 4957]** `    BlockVector<double> solution_delta_cauchy_point(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 4958]** `    LBFGS_update = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4959]** ``：空行，用于分隔逻辑块。
- **[行 4960]** `    const unsigned int LBFGS_m = m_parameters.m_LBFGS_m;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4961]** ``：空行，用于分隔逻辑块。
- **[行 4962]** `    unsigned int LBFGS_iteration = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4963]** ``：空行，用于分隔逻辑块。
- **[行 4964]** `    m_error_residual.reset();`：函数调用语句，触发对应计算或操作。
- **[行 4965]** `    m_error_residual_0.reset();`：函数调用语句，触发对应计算或操作。
- **[行 4966]** `    m_error_residual_norm.reset();`：函数调用语句，触发对应计算或操作。
- **[行 4967]** `    m_error_update.reset();`：函数调用语句，触发对应计算或操作。
- **[行 4968]** `    m_error_update_0.reset();`：函数调用语句，触发对应计算或操作。
- **[行 4969]** `    m_error_update_norm.reset();`：函数调用语句，触发对应计算或操作。
- **[行 4970]** ``：空行，用于分隔逻辑块。
- **[行 4971]** `    if (m_parameters.m_output_iteration_history)`：条件分支：根据当前状态选择执行路径。
- **[行 4972]** `      print_conv_header_LBFGSB();`：函数调用语句，触发对应计算或操作。
- **[行 4973]** ``：空行，用于分隔逻辑块。
- **[行 4974]** `    BlockVector<double> LBFGS_s_vector(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 4975]** `    BlockVector<double> LBFGS_y_vector(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 4976]** `    BlockVector<double> free_dofs(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 4977]** `    BlockVector<double> b0xs_vector(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 4978]** ``：空行，用于分隔逻辑块。
- **[行 4979]** `    // all the list goes from k-m to k-1`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4980]** `    // the front is the oldest quantity,and the end is`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4981]** `    // newest quantity`：注释行，用于说明算法背景、假设或实现意图。
- **[行 4982]** `    std::list<BlockVector<double>> s_vector_list;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4983]** `    std::list<BlockVector<double>> y_vector_list;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4984]** `    std::list<double> s_dot_y_list;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4985]** `    std::list<BlockVector<double>> b0xs_vector_list;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 4986]** ``：空行，用于分隔逻辑块。
- **[行 4987]** `    double line_search_parameter = 0.0;`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 4988]** ``：空行，用于分隔逻辑块。
- **[行 4989]** `    unsigned int lower_bound_number_old = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4990]** `    unsigned int upper_bound_number_old = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4991]** `    unsigned int lowerupper_bound_number_old = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4992]** ``：空行，用于分隔逻辑块。
- **[行 4993]** `    unsigned int lower_bound_number_new = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4994]** `    unsigned int upper_bound_number_new = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4995]** `    unsigned int lowerupper_bound_number_new = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 4996]** ``：空行，用于分隔逻辑块。
- **[行 4997]** `    for (; LBFGS_iteration < m_parameters.m_max_iterations_BFGS; ++LBFGS_iteration)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 4998]** `      {`：作用域边界（代码块开始/结束）。
- **[行 4999]** `	if (m_parameters.m_output_iteration_history)`：条件分支：根据当前状态选择执行路径。
- **[行 5000]** `	  m_logfile << '\t' << std::setw(4) << LBFGS_iteration << ' '`：写日志输出，记录当前计算状态与结果。
- **[行 5001]** `                    << std::flush;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5002]** ``：空行，用于分隔逻辑块。
- **[行 5003]** `        make_constraints(LBFGS_iteration);`：函数调用语句，触发对应计算或操作。
- **[行 5004]** ``：空行，用于分隔逻辑块。
- **[行 5005]** `        // At the first step, we simply distribute the inhomogeneous part of`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5006]** `        // the constraints`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5007]** `        if (LBFGS_iteration == 0)`：条件分支：根据当前状态选择执行路径。
- **[行 5008]** `          {`：作用域边界（代码块开始/结束）。
- **[行 5009]** `            // use the solution from the previous solve on the`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5010]** `            // refined mesh as initial guess`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5011]** `            LBFGS_update = LBFGS_update_refine;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5012]** ``：空行，用于分隔逻辑块。
- **[行 5013]** `            m_constraints.distribute(LBFGS_update);`：函数调用语句，触发对应计算或操作。
- **[行 5014]** `            solution_delta += LBFGS_update;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5015]** `            update_qph_incremental(solution_delta, m_solution);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 5016]** `            assemble_system_rhs_BFGS_parallel(m_solution, m_system_rhs);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 5017]** `            m_logfile << "  | " << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 5018]** ``：空行，用于分隔逻辑块。
- **[行 5019]** `            continue;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5020]** `          }`：作用域边界（代码块开始/结束）。
- **[行 5021]** ``：空行，用于分隔逻辑块。
- **[行 5022]** `        get_error_residual_LBFGSB(m_error_residual,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5023]** `				  solution_delta);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5024]** ``：空行，用于分隔逻辑块。
- **[行 5025]** `        if (LBFGS_iteration == 1)`：条件分支：根据当前状态选择执行路径。
- **[行 5026]** `          m_error_residual_0 = m_error_residual;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5027]** ``：空行，用于分隔逻辑块。
- **[行 5028]** `        m_error_residual_norm = m_error_residual;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5029]** ``：空行，用于分隔逻辑块。
- **[行 5030]** `        if (m_parameters.m_relative_residual)`：条件分支：根据当前状态选择执行路径。
- **[行 5031]** `          m_error_residual_norm.normalize(m_error_residual_0);`：函数调用语句，触发对应计算或操作。
- **[行 5032]** ``：空行，用于分隔逻辑块。
- **[行 5033]** `        if (LBFGS_iteration > 1 && m_error_update_norm.m_u <= m_parameters.m_tol_u_incr`：条件分支：根据当前状态选择执行路径。
- **[行 5034]** `                                && m_error_residual_norm.m_u <= m_parameters.m_tol_u_residual`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5035]** `			        && m_error_update_norm.m_d <= m_parameters.m_tol_d_incr`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5036]** `			        && m_error_residual_norm.m_d <= m_parameters.m_tol_d_residual`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5037]** `				&& lower_bound_number_new == lower_bound_number_old`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5038]** `				&& upper_bound_number_new == upper_bound_number_old`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5039]** `				&& lowerupper_bound_number_new == lowerupper_bound_number_old)`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5040]** `          {`：作用域边界（代码块开始/结束）。
- **[行 5041]** `            if (m_parameters.m_output_iteration_history)`：条件分支：根据当前状态选择执行路径。
- **[行 5042]** `              {`：作用域边界（代码块开始/结束）。
- **[行 5043]** `		m_logfile << "  | ";`：写日志输出，记录当前计算状态与结果。
- **[行 5044]** `		m_logfile << " CONVERGED! " << std::fixed << std::setprecision(3) << std::setw(7)`：写日志输出，记录当前计算状态与结果。
- **[行 5045]** `			  << std::scientific`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5046]** `		      << "           ---      "`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5047]** `		      << "\t\t\t"`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5048]** `		      << "  " << m_error_residual_norm.m_norm`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5049]** `		      << "  " << m_error_residual_norm.m_u`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5050]** `		      << "  " << m_error_residual_norm.m_d`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5051]** `		      << "  " << m_error_update_norm.m_norm`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5052]** `		      << "  " << m_error_update_norm.m_u`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5053]** `		      << "  " << m_error_update_norm.m_d`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5054]** `		      << "  " << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5055]** ``：空行，用于分隔逻辑块。
- **[行 5056]** `		m_logfile << '\t';`：写日志输出，记录当前计算状态与结果。
- **[行 5057]** `		for (unsigned int i = 0; i < 130; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5058]** `		  m_logfile << '_';`：写日志输出，记录当前计算状态与结果。
- **[行 5059]** `		m_logfile << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 5060]** `              }`：作用域边界（代码块开始/结束）。
- **[行 5061]** ``：空行，用于分隔逻辑块。
- **[行 5062]** `            m_logfile << "\t\tThe current L-BFGS-B step converges in "`：写日志输出，记录当前计算状态与结果。
- **[行 5063]** `        	      << LBFGS_iteration`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5064]** `        	      << " iterations." << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5065]** `            m_logfile << "\t\tNumber of active lower bounds not changed." << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 5066]** `    	    m_logfile << "\t\tNumber of active upper bounds not changed." << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 5067]** `    	    m_logfile << "\t\tNumber of active lower-upper bounds not changed." << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 5068]** ``：空行，用于分隔逻辑块。
- **[行 5069]** `            if (m_parameters.m_relative_residual)`：条件分支：根据当前状态选择执行路径。
- **[行 5070]** `              {`：作用域边界（代码块开始/结束）。
- **[行 5071]** `		m_logfile << "\t\tProjected gradient of disp. (relative): "`：写日志输出，记录当前计算状态与结果。
- **[行 5072]** `			  << m_error_residual_norm.m_u << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5073]** ``：空行，用于分隔逻辑块。
- **[行 5074]** `		m_logfile << "\t\tProjected gradient of disp. (absolute): "`：写日志输出，记录当前计算状态与结果。
- **[行 5075]** `			  << m_error_residual_norm.m_u * m_error_residual_0.m_u << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5076]** ``：空行，用于分隔逻辑块。
- **[行 5077]** `		m_logfile << "\t\tProjected gradient of phasefield (relative): "`：写日志输出，记录当前计算状态与结果。
- **[行 5078]** `			  << m_error_residual_norm.m_d << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5079]** ``：空行，用于分隔逻辑块。
- **[行 5080]** `		m_logfile << "\t\tProjected gradient of phasefield (absolute): "`：写日志输出，记录当前计算状态与结果。
- **[行 5081]** `			  << m_error_residual_norm.m_d * m_error_residual_0.m_d << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5082]** ``：空行，用于分隔逻辑块。
- **[行 5083]** `		m_logfile << "\t\tRelative increment of disp.: "`：写日志输出，记录当前计算状态与结果。
- **[行 5084]** `			  << m_error_update_norm.m_u << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5085]** ``：空行，用于分隔逻辑块。
- **[行 5086]** `		m_logfile << "\t\tAbsolute increment of disp.: "`：写日志输出，记录当前计算状态与结果。
- **[行 5087]** `			  << m_error_update_norm.m_u * m_error_update_0.m_u << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5088]** ``：空行，用于分隔逻辑块。
- **[行 5089]** `		m_logfile << "\t\tRelative increment of phasefield: "`：写日志输出，记录当前计算状态与结果。
- **[行 5090]** `			  << m_error_update_norm.m_d << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5091]** ``：空行，用于分隔逻辑块。
- **[行 5092]** `		m_logfile << "\t\tAbsolute increment of phasefield: "`：写日志输出，记录当前计算状态与结果。
- **[行 5093]** `			  << m_error_update_norm.m_d * m_error_update_0.m_d << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5094]** `              }`：作用域边界（代码块开始/结束）。
- **[行 5095]** `            else`：条件分支的兜底路径。
- **[行 5096]** `              {`：作用域边界（代码块开始/结束）。
- **[行 5097]** `		m_logfile << "\t\tProjected gradient of disp. (absolute): "`：写日志输出，记录当前计算状态与结果。
- **[行 5098]** `			  << m_error_residual_norm.m_u << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5099]** ``：空行，用于分隔逻辑块。
- **[行 5100]** `		m_logfile << "\t\tProjected gradient of phasefield (absolute): "`：写日志输出，记录当前计算状态与结果。
- **[行 5101]** `			  << m_error_residual_norm.m_d << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5102]** ``：空行，用于分隔逻辑块。
- **[行 5103]** `		m_logfile << "\t\tAbsolute increment of disp.: "`：写日志输出，记录当前计算状态与结果。
- **[行 5104]** `			  << m_error_update_norm.m_u << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5105]** ``：空行，用于分隔逻辑块。
- **[行 5106]** `		m_logfile << "\t\tAbsolute increment of phasefield: "`：写日志输出，记录当前计算状态与结果。
- **[行 5107]** `			  << m_error_update_norm.m_d << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5108]** `              }`：作用域边界（代码块开始/结束）。
- **[行 5109]** ``：空行，用于分隔逻辑块。
- **[行 5110]** `            break;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5111]** `          }`：作用域边界（代码块开始/结束）。
- **[行 5112]** ``：空行，用于分隔逻辑块。
- **[行 5113]** `	// assemble the initial B_0 matrix at the k-th L-BFGS iteration`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5114]** `	// m_solution is the old solution from the previous converged step`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5115]** `	// it is needed only for the viscosity term`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5116]** `	// the output is m_tangent_matrix (B^0_k)`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5117]** `	assemble_system_B0(m_solution);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 5118]** ``：空行，用于分隔逻辑块。
- **[行 5119]** `	// B^0_k * s_vector has to be completely recalculated from scratch`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5120]** `	// at each L-BFGS iteration, since B^0_k is different`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5121]** `	b0xs_vector_list.clear();`：函数调用语句，触发对应计算或操作。
- **[行 5122]** `	for (auto itr = s_vector_list.begin(); itr != s_vector_list.end(); ++itr)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5123]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 5124]** `	    m_tangent_matrix.vmult(b0xs_vector, *itr);`：函数调用语句，触发对应计算或操作。
- **[行 5125]** `	    b0xs_vector_list.push_back(b0xs_vector);`：函数调用语句，触发对应计算或操作。
- **[行 5126]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 5127]** ``：空行，用于分隔逻辑块。
- **[行 5128]** `	// In the iteration LBFGS_iteration = 0, only the essential boundary conditions`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5129]** `	// are applied.`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5130]** `	// WHen LBFGS_iteration = 1, it is the first step of LBFGS update, and the`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5131]** `	// s_vector_list and y_vector_list are empty.`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5132]** `	// Since the pair of s and y will only be added to the list if s dot y > tol,`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5133]** `	// it is safer to decide the matrix dimension by the size of the list.`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5134]** `	const unsigned int list_size = s_vector_list.size();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5135]** `	const auto itr_s_begin    = s_vector_list.begin();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5136]** `	const auto itr_y_begin    = y_vector_list.begin();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5137]** `	const auto itr_b0xs_begin = b0xs_vector_list.begin();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5138]** `	const auto itr_s_dot_y_begin = s_dot_y_list.begin();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5139]** ``：空行，用于分隔逻辑块。
- **[行 5140]** `	FullMatrix<double> sTxBxs_matrix(list_size);`：函数调用语句，触发对应计算或操作。
- **[行 5141]** `	sTxBxs_matrix = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5142]** `	for (unsigned int i = 0; i < list_size; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5143]** `	  for (unsigned int j = 0; j <= i; ++j)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5144]** `	    {`：作用域边界（代码块开始/结束）。
- **[行 5145]** `	      sTxBxs_matrix(i, j) = (*std::next(itr_s_begin,    i))`：调用 C++ 标准库工具函数/容器接口。
- **[行 5146]** `		                  * (*std::next(itr_b0xs_begin, j));`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5147]** `	    }`：作用域边界（代码块开始/结束）。
- **[行 5148]** `	for (unsigned int i = 0; i < list_size; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5149]** `	  for (unsigned int j = i + 1; j < list_size; ++j)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5150]** `	    {`：作用域边界（代码块开始/结束）。
- **[行 5151]** `	      sTxBxs_matrix(i, j) = sTxBxs_matrix(j, i);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5152]** `	    }`：作用域边界（代码块开始/结束）。
- **[行 5153]** ``：空行，用于分隔逻辑块。
- **[行 5154]** `	FullMatrix<double> D_matrix(list_size);`：函数调用语句，触发对应计算或操作。
- **[行 5155]** `	D_matrix = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5156]** `	for (unsigned int i = 0; i < list_size; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5157]** `	  D_matrix(i, i) = (*std::next(itr_s_dot_y_begin, i));`：调用 C++ 标准库工具函数/容器接口。
- **[行 5158]** ``：空行，用于分隔逻辑块。
- **[行 5159]** `	FullMatrix<double> L_matrix(list_size);`：函数调用语句，触发对应计算或操作。
- **[行 5160]** `	L_matrix = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5161]** `	for (unsigned int i = 0; i < list_size; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5162]** `	  for (unsigned int j = 0; j < i; ++j)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5163]** `	    L_matrix(i, j) = (*std::next(itr_s_begin, i))`：调用 C++ 标准库工具函数/容器接口。
- **[行 5164]** `                           * (*std::next(itr_y_begin, j));`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5165]** ``：空行，用于分隔逻辑块。
- **[行 5166]** `	FullMatrix<double> M_matrix_inv(2 * list_size);`：函数调用语句，触发对应计算或操作。
- **[行 5167]** `	FullMatrix<double> M_matrix(2 * list_size);`：函数调用语句，触发对应计算或操作。
- **[行 5168]** ``：空行，用于分隔逻辑块。
- **[行 5169]** `	M_matrix_inv = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5170]** `	for (unsigned int i = 0; i < list_size; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5171]** `	  M_matrix_inv(i, i) = -D_matrix(i, i);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5172]** ``：空行，用于分隔逻辑块。
- **[行 5173]** `	for (unsigned int i = 0; i < list_size; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5174]** `          for (unsigned int j = 0; j < list_size; ++j)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5175]** `            {`：作用域边界（代码块开始/结束）。
- **[行 5176]** `              M_matrix_inv(i + list_size, j + list_size) = sTxBxs_matrix(i, j);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5177]** `              M_matrix_inv(i + list_size, j            ) = L_matrix(i, j);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5178]** `              M_matrix_inv(i            , j + list_size) = L_matrix(j, i);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5179]** `            }`：作用域边界（代码块开始/结束）。
- **[行 5180]** ``：空行，用于分隔逻辑块。
- **[行 5181]** `	if (!M_matrix_inv.empty())`：条件分支：根据当前状态选择执行路径。
- **[行 5182]** `	  M_matrix.invert(M_matrix_inv);`：函数调用语句，触发对应计算或操作。
- **[行 5183]** ``：空行，用于分隔逻辑块。
- **[行 5184]** `	m_active_set_phasefield = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5185]** `	calculate_cauchy_point(m_tangent_matrix,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 5186]** `			       y_vector_list,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5187]** `			       b0xs_vector_list,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5188]** `			       M_matrix,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5189]** `			       m_system_rhs,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5190]** `			       solution_delta,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5191]** `			       solution_delta_cauchy_point);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5192]** ``：空行，用于分隔逻辑块。
- **[行 5193]** `	// We need to find out which DOFs are free:`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5194]** `	// no essential boundary conditions, no hanging node constraints`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5195]** `	// no active box constraints`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5196]** `	unsigned int free_disp_number = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5197]** `	for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5198]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 5199]** `	    if (m_constraints.is_constrained(i))`：条件分支：根据当前状态选择执行路径。
- **[行 5200]** `	      free_dofs.block(m_u_dof)[i] = -1;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5201]** `	    else`：条件分支的兜底路径。
- **[行 5202]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 5203]** `	        free_dofs.block(m_u_dof)[i] = 1;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5204]** `	        ++free_disp_number;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5205]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 5206]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 5207]** ``：空行，用于分隔逻辑块。
- **[行 5208]** `	unsigned int free_phasefield_number = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5209]** `	for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5210]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 5211]** `	    if (   m_constraints.is_constrained(i + m_dofs_per_block[m_u_dof])`：条件分支：根据当前状态选择执行路径。
- **[行 5212]** `		|| m_active_set_phasefield(i) > 0.5)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5213]** `	      free_dofs.block(m_d_dof)[i] = -1;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5214]** `	    else`：条件分支的兜底路径。
- **[行 5215]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 5216]** `	        free_dofs.block(m_d_dof)[i] = 1;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5217]** `	        ++free_phasefield_number;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5218]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 5219]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 5220]** ``：空行，用于分隔逻辑块。
- **[行 5221]** `	// temp_vector_1 = x^c - x_k`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5222]** `	BlockVector<double> temp_vector_1(solution_delta_cauchy_point);`：函数调用语句，触发对应计算或操作。
- **[行 5223]** `	temp_vector_1 -= solution_delta;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5224]** ``：空行，用于分隔逻辑块。
- **[行 5225]** `	// temp_vector_2 = B_0 * (x^c - x_k)`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5226]** `	BlockVector<double> temp_vector_2(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 5227]** `	m_tangent_matrix.vmult(temp_vector_2, temp_vector_1);`：函数调用语句，触发对应计算或操作。
- **[行 5228]** ``：空行，用于分隔逻辑块。
- **[行 5229]** `	// temp_vector_3 = W^T * (x^c - x_k)`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5230]** `	Vector<double> temp_vector_3(2 * list_size);`：函数调用语句，触发对应计算或操作。
- **[行 5231]** `	for (unsigned int i = 0; i < list_size; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5232]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 5233]** `	    temp_vector_3(i)             = (*std::next(itr_y_begin,    i)) * temp_vector_1;`：调用 C++ 标准库工具函数/容器接口。
- **[行 5234]** `	    temp_vector_3(i + list_size) = (*std::next(itr_b0xs_begin, i)) * temp_vector_1;`：调用 C++ 标准库工具函数/容器接口。
- **[行 5235]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 5236]** ``：空行，用于分隔逻辑块。
- **[行 5237]** `	// temp_vector_4 = M * W^T * (x^c - x_k)`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5238]** `	Vector<double> temp_vector_4(2 * list_size);`：函数调用语句，触发对应计算或操作。
- **[行 5239]** `	if (list_size > 0)`：条件分支：根据当前状态选择执行路径。
- **[行 5240]** `	  M_matrix.vmult(temp_vector_4, temp_vector_3);`：函数调用语句，触发对应计算或操作。
- **[行 5241]** ``：空行，用于分隔逻辑块。
- **[行 5242]** `	// temp_vector_5 = W * M * W^T * (x^c - x_k)`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5243]** `	BlockVector<double> temp_vector_5(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 5244]** `	for (unsigned int i = 0; i < list_size; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5245]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 5246]** `	    temp_vector_5.add(temp_vector_4(i),             (*std::next(itr_y_begin,    i)));`：调用 C++ 标准库工具函数/容器接口。
- **[行 5247]** `	    temp_vector_5.add(temp_vector_4(i + list_size), (*std::next(itr_b0xs_begin, i)));`：调用 C++ 标准库工具函数/容器接口。
- **[行 5248]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 5249]** ``：空行，用于分隔逻辑块。
- **[行 5250]** `	// temp_vector_2 = B_0 * (x^c - x_k) - W * M * W^T * (x^c - x_k)`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5251]** `	if (list_size > 0)`：条件分支：根据当前状态选择执行路径。
- **[行 5252]** `	  temp_vector_2 -= temp_vector_5;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5253]** ``：空行，用于分隔逻辑块。
- **[行 5254]** `	// temp_vector_2 = g + B_0 * (x^c - x_k) - W * M * W^T * (x^c - x_k)`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5255]** `	temp_vector_2 += m_system_rhs;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5256]** ``：空行，用于分隔逻辑块。
- **[行 5257]** `	// temp_vector_2 = Z^T * [g + B_0 * (x^c - x_k) - W * M * W^T * (x^c - x_k)]`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5258]** `	for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5259]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 5260]** `	    if (free_dofs.block(m_u_dof)[i] < 0)`：条件分支：根据当前状态选择执行路径。
- **[行 5261]** `	      temp_vector_2.block(m_u_dof)[i] = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5262]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 5263]** ``：空行，用于分隔逻辑块。
- **[行 5264]** `	for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5265]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 5266]** `	    if (free_dofs.block(m_d_dof)[i] < 0)`：条件分支：根据当前状态选择执行路径。
- **[行 5267]** `	      temp_vector_2.block(m_d_dof)[i] = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5268]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 5269]** ``：空行，用于分隔逻辑块。
- **[行 5270]** `	BlockVector<double> rhs_vector(temp_vector_2);`：函数调用语句，触发对应计算或操作。
- **[行 5271]** `	rhs_vector *= -1;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5272]** ``：空行，用于分隔逻辑块。
- **[行 5273]** `	BlockVector<double> search_direction(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 5274]** `	search_direction = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5275]** ``：空行，用于分隔逻辑块。
- **[行 5276]** `	unsigned int cg_iterations = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5277]** ``：空行，用于分隔逻辑块。
- **[行 5278]** `	double alpha_backtrack = 1.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5279]** ``：空行，用于分隔逻辑块。
- **[行 5280]** `	if (m_parameters.m_type_linear_solver == "CG")`：条件分支：根据当前状态选择执行路径。
- **[行 5281]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 5282]** `	    m_timer.enter_subsection("Subspace CG solve (inverse operator)");`：函数调用语句，触发对应计算或操作。
- **[行 5283]** ``：空行，用于分隔逻辑块。
- **[行 5284]** `	    //const double rc_hat_norm = rhs_vector.l2_norm();`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5285]** `	    const double cg_tol = m_parameters.m_CG_tolerace; //std::min( 0.1, std::sqrt(rc_hat_norm) ) * rc_hat_norm;`：调用 C++ 标准库工具函数/容器接口。
- **[行 5286]** ``：空行，用于分隔逻辑块。
- **[行 5287]** `	    zT_B0_z(free_dofs, m_tangent_matrix);`：函数调用语句，触发对应计算或操作。
- **[行 5288]** ``：空行，用于分隔逻辑块。
- **[行 5289]** `	    const auto op_zT_B0_z = block_operator(m_tangent_matrix);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5290]** ``：空行，用于分隔逻辑块。
- **[行 5291]** `	    if (list_size > 0)`：条件分支：根据当前状态选择执行路径。
- **[行 5292]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 5293]** `		std::list<BlockVector<double>> zT_y_list;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5294]** `		BlockVector<double> zT_y_vector(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 5295]** `		for (auto itr = y_vector_list.begin(); itr != y_vector_list.end(); ++itr)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5296]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 5297]** `		    zT_y_vector = (*itr);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5298]** `		    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5299]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 5300]** `			if (free_dofs.block(m_u_dof)[i] < 0)`：条件分支：根据当前状态选择执行路径。
- **[行 5301]** `			  zT_y_vector.block(m_u_dof)[i] = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5302]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 5303]** ``：空行，用于分隔逻辑块。
- **[行 5304]** `		    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5305]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 5306]** `			if (free_dofs.block(m_d_dof)[i] < 0)`：条件分支：根据当前状态选择执行路径。
- **[行 5307]** `			  zT_y_vector.block(m_d_dof)[i] = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5308]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 5309]** ``：空行，用于分隔逻辑块。
- **[行 5310]** `		    zT_y_list.push_back(zT_y_vector);`：函数调用语句，触发对应计算或操作。
- **[行 5311]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 5312]** ``：空行，用于分隔逻辑块。
- **[行 5313]** `		std::list<BlockVector<double>> zT_b0xs_list;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5314]** `		BlockVector<double> zT_b0xs_vector(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 5315]** `		for (auto itr = b0xs_vector_list.begin(); itr != b0xs_vector_list.end(); ++itr)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5316]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 5317]** `		    zT_b0xs_vector = (*itr);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5318]** `		    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5319]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 5320]** `			if (free_dofs.block(m_u_dof)[i] < 0)`：条件分支：根据当前状态选择执行路径。
- **[行 5321]** `			  zT_b0xs_vector.block(m_u_dof)[i] = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5322]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 5323]** ``：空行，用于分隔逻辑块。
- **[行 5324]** `		    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5325]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 5326]** `			if (free_dofs.block(m_d_dof)[i] < 0)`：条件分支：根据当前状态选择执行路径。
- **[行 5327]** `			  zT_b0xs_vector.block(m_d_dof)[i] = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5328]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 5329]** ``：空行，用于分隔逻辑块。
- **[行 5330]** `		    zT_b0xs_list.push_back(zT_b0xs_vector);`：函数调用语句，触发对应计算或操作。
- **[行 5331]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 5332]** ``：空行，用于分隔逻辑块。
- **[行 5333]** `		const auto op_M_matrix = linear_operator(M_matrix);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5334]** ``：空行，用于分隔逻辑块。
- **[行 5335]** `		FullMatrix<double> zT_W_matrix_u(m_dofs_per_block[m_u_dof], 2*list_size);`：函数调用语句，触发对应计算或操作。
- **[行 5336]** `		unsigned int j = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5337]** `		for (auto itr = zT_y_list.begin(); itr != zT_y_list.end(); ++itr)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5338]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 5339]** `		    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5340]** `		      zT_W_matrix_u(i, j) = (*itr).block(m_u_dof)[i];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5341]** `		    ++j;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5342]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 5343]** `		j = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5344]** `		for (auto itr = zT_b0xs_list.begin(); itr != zT_b0xs_list.end(); ++itr)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5345]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 5346]** `		    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5347]** `		      zT_W_matrix_u(i, j + list_size) = (*itr).block(m_u_dof)[i];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5348]** `		    ++j;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5349]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 5350]** ``：空行，用于分隔逻辑块。
- **[行 5351]** `		FullMatrix<double> zT_W_matrix_d(m_dofs_per_block[m_d_dof], 2*list_size);`：函数调用语句，触发对应计算或操作。
- **[行 5352]** `		j = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5353]** `		for (auto itr = zT_y_list.begin(); itr != zT_y_list.end(); ++itr)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5354]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 5355]** `		    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5356]** `		      zT_W_matrix_d(i, j) = (*itr).block(m_d_dof)[i];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5357]** `		    ++j;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5358]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 5359]** `		j = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5360]** `		for (auto itr = zT_b0xs_list.begin(); itr != zT_b0xs_list.end(); ++itr)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5361]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 5362]** `		    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5363]** `		      zT_W_matrix_d(i, j + list_size) = (*itr).block(m_d_dof)[i];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5364]** `		    ++j;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5365]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 5366]** ``：空行，用于分隔逻辑块。
- **[行 5367]** `		const auto op_zT_W_matrix_u = linear_operator(zT_W_matrix_u);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5368]** `		const auto op_zT_W_matrix_d = linear_operator(zT_W_matrix_d);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5369]** ``：空行，用于分隔逻辑块。
- **[行 5370]** `		const auto op_uMuT = op_zT_W_matrix_u * op_M_matrix * transpose_operator(op_zT_W_matrix_u);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5371]** ``：空行，用于分隔逻辑块。
- **[行 5372]** `		const auto op_uMdT = op_zT_W_matrix_u * op_M_matrix * transpose_operator(op_zT_W_matrix_d);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5373]** ``：空行，用于分隔逻辑块。
- **[行 5374]** `		const auto op_dMuT = op_zT_W_matrix_d * op_M_matrix * transpose_operator(op_zT_W_matrix_u);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5375]** ``：空行，用于分隔逻辑块。
- **[行 5376]** `		const auto op_dMdT = op_zT_W_matrix_d * op_M_matrix * transpose_operator(op_zT_W_matrix_d);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5377]** ``：空行，用于分隔逻辑块。
- **[行 5378]** `		const auto op_zT_wMwT_z = block_operator<2, 2, BlockVector<double>>({op_uMuT, op_uMdT,`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5379]** `										     op_dMuT, op_dMdT});`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5380]** ``：空行，用于分隔逻辑块。
- **[行 5381]** `		const auto op_total = op_zT_B0_z - op_zT_wMwT_z;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5382]** ``：空行，用于分隔逻辑块。
- **[行 5383]** `		SolverControl            solver_control(1e5, cg_tol);`：函数调用语句，触发对应计算或操作。
- **[行 5384]** `		SolverCG<BlockVector<double>> cg(solver_control);`：函数调用语句，触发对应计算或操作。
- **[行 5385]** ``：空行，用于分隔逻辑块。
- **[行 5386]** `		if (m_parameters.m_type_preconditioner == "None")`：条件分支：根据当前状态选择执行路径。
- **[行 5387]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 5388]** `		    // somehow op_total_inv has to be made const, or the`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5389]** `		    // program will have compliation error`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5390]** `		    const auto op_total_inv = inverse_operator(op_total, cg);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5391]** `		    op_total_inv.vmult(search_direction, rhs_vector);`：函数调用语句，触发对应计算或操作。
- **[行 5392]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 5393]** `		else if (m_parameters.m_type_preconditioner == "Jacobi")`：多分支条件判断，处理备选情形。
- **[行 5394]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 5395]** `		    usr_Jacobi_preconditioner preconditioner(m_tangent_matrix);`：函数调用语句，触发对应计算或操作。
- **[行 5396]** ``：空行，用于分隔逻辑块。
- **[行 5397]** `		    const auto op_total_inv = inverse_operator(op_total, cg, preconditioner);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5398]** `		    op_total_inv.vmult(search_direction, rhs_vector);`：函数调用语句，触发对应计算或操作。
- **[行 5399]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 5400]** `		else if (m_parameters.m_type_preconditioner == "LU")`：多分支条件判断，处理备选情形。
- **[行 5401]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 5402]** `		    SparseDirectUMFPACK matrix_factorization;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5403]** `		    matrix_factorization.initialize(m_tangent_matrix);`：函数调用语句，触发对应计算或操作。
- **[行 5404]** ``：空行，用于分隔逻辑块。
- **[行 5405]** `		    usr_sparseLU_preconditioner preconditioner(matrix_factorization);`：函数调用语句，触发对应计算或操作。
- **[行 5406]** `		    const auto op_total_inv = inverse_operator(op_total, cg, preconditioner);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5407]** `		    op_total_inv.vmult(search_direction, rhs_vector);`：函数调用语句，触发对应计算或操作。
- **[行 5408]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 5409]** `		else if (m_parameters.m_type_preconditioner == "ILU")`：多分支条件判断，处理备选情形。
- **[行 5410]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 5411]** `		    SparseILU<double> SparseILU_disp;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5412]** `		    SparseILU_disp.initialize(m_tangent_matrix.block(0, 0));`：函数调用语句，触发对应计算或操作。
- **[行 5413]** `		    SparseILU<double> SparseILU_phasefield;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5414]** `		    SparseILU_phasefield.initialize(m_tangent_matrix.block(1, 1));`：函数调用语句，触发对应计算或操作。
- **[行 5415]** ``：空行，用于分隔逻辑块。
- **[行 5416]** `		    usr_sparseILU_preconditioner preconditioner(SparseILU_disp,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5417]** `								SparseILU_phasefield);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5418]** ``：空行，用于分隔逻辑块。
- **[行 5419]** `		    const auto op_total_inv = inverse_operator(op_total, cg, preconditioner);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5420]** `		    op_total_inv.vmult(search_direction, rhs_vector);`：函数调用语句，触发对应计算或操作。
- **[行 5421]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 5422]** `		else`：条件分支的兜底路径。
- **[行 5423]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 5424]** `		    AssertThrow(false, ExcMessage("Preconditioner type not implemented"));`：运行期断言/检查，验证输入与状态合法性。
- **[行 5425]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 5426]** ``：空行，用于分隔逻辑块。
- **[行 5427]** `		cg_iterations = solver_control.last_step();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5428]** `	      } // if (list_size > 0)`：作用域边界（代码块开始/结束）。
- **[行 5429]** `	    else`：条件分支的兜底路径。
- **[行 5430]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 5431]** `		const auto op_total = op_zT_B0_z;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5432]** `		SolverControl            solver_control(1e5, cg_tol);`：函数调用语句，触发对应计算或操作。
- **[行 5433]** `		SolverCG<BlockVector<double>> cg(solver_control);`：函数调用语句，触发对应计算或操作。
- **[行 5434]** ``：空行，用于分隔逻辑块。
- **[行 5435]** `		if (m_parameters.m_type_preconditioner == "None")`：条件分支：根据当前状态选择执行路径。
- **[行 5436]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 5437]** `		    // somehow op_total_inv has to be made const, or the`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5438]** `		    // program will have compliation error`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5439]** `		    const auto op_total_inv = inverse_operator(op_total, cg);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5440]** `		    op_total_inv.vmult(search_direction, rhs_vector);`：函数调用语句，触发对应计算或操作。
- **[行 5441]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 5442]** `		else if (m_parameters.m_type_preconditioner == "Jacobi")`：多分支条件判断，处理备选情形。
- **[行 5443]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 5444]** `		    usr_Jacobi_preconditioner preconditioner(m_tangent_matrix);`：函数调用语句，触发对应计算或操作。
- **[行 5445]** ``：空行，用于分隔逻辑块。
- **[行 5446]** `		    const auto op_total_inv = inverse_operator(op_total, cg, preconditioner);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5447]** `		    op_total_inv.vmult(search_direction, rhs_vector);`：函数调用语句，触发对应计算或操作。
- **[行 5448]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 5449]** `		else if (m_parameters.m_type_preconditioner == "LU")`：多分支条件判断，处理备选情形。
- **[行 5450]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 5451]** `		    SparseDirectUMFPACK matrix_factorization;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5452]** `		    matrix_factorization.initialize(m_tangent_matrix);`：函数调用语句，触发对应计算或操作。
- **[行 5453]** ``：空行，用于分隔逻辑块。
- **[行 5454]** `		    usr_sparseLU_preconditioner preconditioner(matrix_factorization);`：函数调用语句，触发对应计算或操作。
- **[行 5455]** `		    const auto op_total_inv = inverse_operator(op_total, cg, preconditioner);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5456]** `		    op_total_inv.vmult(search_direction, rhs_vector);`：函数调用语句，触发对应计算或操作。
- **[行 5457]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 5458]** `		else if (m_parameters.m_type_preconditioner == "ILU")`：多分支条件判断，处理备选情形。
- **[行 5459]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 5460]** `		    SparseILU<double> SparseILU_disp;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5461]** `		    SparseILU_disp.initialize(m_tangent_matrix.block(0, 0));`：函数调用语句，触发对应计算或操作。
- **[行 5462]** `		    SparseILU<double> SparseILU_phasefield;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5463]** `		    SparseILU_phasefield.initialize(m_tangent_matrix.block(1, 1));`：函数调用语句，触发对应计算或操作。
- **[行 5464]** ``：空行，用于分隔逻辑块。
- **[行 5465]** `		    usr_sparseILU_preconditioner preconditioner(SparseILU_disp,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5466]** `								SparseILU_phasefield);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5467]** ``：空行，用于分隔逻辑块。
- **[行 5468]** `		    const auto op_total_inv = inverse_operator(op_total, cg, preconditioner);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5469]** `		    op_total_inv.vmult(search_direction, rhs_vector);`：函数调用语句，触发对应计算或操作。
- **[行 5470]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 5471]** `		else`：条件分支的兜底路径。
- **[行 5472]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 5473]** `		    AssertThrow(false, ExcMessage("Preconditioner type not implemented"));`：运行期断言/检查，验证输入与状态合法性。
- **[行 5474]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 5475]** ``：空行，用于分隔逻辑块。
- **[行 5476]** `		cg_iterations = solver_control.last_step();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5477]** `	      } // // if (list_size == 0)`：作用域边界（代码块开始/结束）。
- **[行 5478]** ``：空行，用于分隔逻辑块。
- **[行 5479]** `            m_timer.leave_subsection();`：函数调用语句，触发对应计算或操作。
- **[行 5480]** `	  } // if (m_parameters.m_type_linear_solver == "CG")`：作用域边界（代码块开始/结束）。
- **[行 5481]** `	else if (m_parameters.m_type_linear_solver == "Direct")`：多分支条件判断，处理备选情形。
- **[行 5482]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 5483]** `	    m_timer.enter_subsection("Subspace direct solve (LU factorization)");`：函数调用语句，触发对应计算或操作。
- **[行 5484]** ``：空行，用于分隔逻辑块。
- **[行 5485]** `	    zT_B0_z(free_dofs, m_tangent_matrix);`：函数调用语句，触发对应计算或操作。
- **[行 5486]** ``：空行，用于分隔逻辑块。
- **[行 5487]** `	    SparseDirectUMFPACK zT_B0_z_inv;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5488]** `	    zT_B0_z_inv.initialize(m_tangent_matrix);`：函数调用语句，触发对应计算或操作。
- **[行 5489]** ``：空行，用于分隔逻辑块。
- **[行 5490]** `	    //SparseDirectUMFPACK zT_B0_z_inv_disp;`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5491]** `	    //zT_B0_z_inv_disp.initialize(m_tangent_matrix.block(m_u_dof, m_u_dof));`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5492]** ``：空行，用于分隔逻辑块。
- **[行 5493]** `	    //SparseDirectUMFPACK zT_B0_z_inv_phasefield;`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5494]** `	    //zT_B0_z_inv_phasefield.initialize(m_tangent_matrix.block(m_d_dof, m_d_dof));`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5495]** ``：空行，用于分隔逻辑块。
- **[行 5496]** `	    m_timer.leave_subsection();`：函数调用语句，触发对应计算或操作。
- **[行 5497]** ``：空行，用于分隔逻辑块。
- **[行 5498]** `	    m_timer.enter_subsection("Subspace direct solve (LU matrix-vector multiplication)");`：函数调用语句，触发对应计算或操作。
- **[行 5499]** ``：空行，用于分隔逻辑块。
- **[行 5500]** `	    zT_B0_z_inv.vmult(search_direction, rhs_vector);`：函数调用语句，触发对应计算或操作。
- **[行 5501]** `	    //zT_B0_z_inv_disp.vmult(search_direction.block(m_u_dof), rhs_vector.block(m_u_dof));`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5502]** `	    //zT_B0_z_inv_phasefield.vmult(search_direction.block(m_d_dof), rhs_vector.block(m_d_dof));`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5503]** ``：空行，用于分隔逻辑块。
- **[行 5504]** `	    BlockVector<double> update_vector(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 5505]** `	    update_vector = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5506]** `	    if (list_size > 0)`：条件分支：根据当前状态选择执行路径。
- **[行 5507]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 5508]** `		std::list<BlockVector<double>> zT_B0_z_inv_zT_y_list;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5509]** `		std::list<BlockVector<double>> zT_y_list;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5510]** `		BlockVector<double> zT_y_vector(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 5511]** `		BlockVector<double> zT_B0_z_inv_zT_y_vector(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 5512]** `		for (auto itr = y_vector_list.begin(); itr != y_vector_list.end(); ++itr)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5513]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 5514]** `		    zT_y_vector = (*itr);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5515]** `		    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5516]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 5517]** `			if (free_dofs.block(m_u_dof)[i] < 0)`：条件分支：根据当前状态选择执行路径。
- **[行 5518]** `			  zT_y_vector.block(m_u_dof)[i] = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5519]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 5520]** ``：空行，用于分隔逻辑块。
- **[行 5521]** `		    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5522]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 5523]** `			if (free_dofs.block(m_d_dof)[i] < 0)`：条件分支：根据当前状态选择执行路径。
- **[行 5524]** `			  zT_y_vector.block(m_d_dof)[i] = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5525]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 5526]** ``：空行，用于分隔逻辑块。
- **[行 5527]** `		    zT_y_list.push_back(zT_y_vector);`：函数调用语句，触发对应计算或操作。
- **[行 5528]** ``：空行，用于分隔逻辑块。
- **[行 5529]** `		    //zT_B0_z_inv_disp.vmult(zT_B0_z_inv_zT_y_vector.block(m_u_dof), zT_y_vector.block(m_u_dof));`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5530]** `		    //zT_B0_z_inv_phasefield.vmult(zT_B0_z_inv_zT_y_vector.block(m_d_dof), zT_y_vector.block(m_d_dof));`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5531]** `		    zT_B0_z_inv.vmult(zT_B0_z_inv_zT_y_vector, zT_y_vector);`：函数调用语句，触发对应计算或操作。
- **[行 5532]** ``：空行，用于分隔逻辑块。
- **[行 5533]** `		    zT_B0_z_inv_zT_y_list.push_back(zT_B0_z_inv_zT_y_vector);`：函数调用语句，触发对应计算或操作。
- **[行 5534]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 5535]** ``：空行，用于分隔逻辑块。
- **[行 5536]** `		std::list<BlockVector<double>> zT_B0_z_inv_zT_b0xs_list;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5537]** `		std::list<BlockVector<double>> zT_b0xs_list;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5538]** `		BlockVector<double> zT_b0xs_vector(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 5539]** `		BlockVector<double> zT_B0_z_inv_zT_b0xs_vector(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 5540]** `		for (auto itr = b0xs_vector_list.begin(); itr != b0xs_vector_list.end(); ++itr)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5541]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 5542]** `		    zT_b0xs_vector = (*itr);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5543]** `		    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5544]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 5545]** `			if (free_dofs.block(m_u_dof)[i] < 0)`：条件分支：根据当前状态选择执行路径。
- **[行 5546]** `			  zT_b0xs_vector.block(m_u_dof)[i] = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5547]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 5548]** ``：空行，用于分隔逻辑块。
- **[行 5549]** `		    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5550]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 5551]** `			if (free_dofs.block(m_d_dof)[i] < 0)`：条件分支：根据当前状态选择执行路径。
- **[行 5552]** `			  zT_b0xs_vector.block(m_d_dof)[i] = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5553]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 5554]** ``：空行，用于分隔逻辑块。
- **[行 5555]** `		    zT_b0xs_list.push_back(zT_b0xs_vector);`：函数调用语句，触发对应计算或操作。
- **[行 5556]** ``：空行，用于分隔逻辑块。
- **[行 5557]** `		    //zT_B0_z_inv_disp.vmult(zT_B0_z_inv_zT_b0xs_vector.block(m_u_dof), zT_b0xs_vector.block(m_u_dof));`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5558]** `		    //zT_B0_z_inv_phasefield.vmult(zT_B0_z_inv_zT_b0xs_vector.block(m_d_dof), zT_b0xs_vector.block(m_d_dof));`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5559]** `		    zT_B0_z_inv.vmult(zT_B0_z_inv_zT_b0xs_vector, zT_b0xs_vector);`：函数调用语句，触发对应计算或操作。
- **[行 5560]** ``：空行，用于分隔逻辑块。
- **[行 5561]** `		    zT_B0_z_inv_zT_b0xs_list.push_back(zT_B0_z_inv_zT_b0xs_vector);`：函数调用语句，触发对应计算或操作。
- **[行 5562]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 5563]** ``：空行，用于分隔逻辑块。
- **[行 5564]** `		FullMatrix<double> wT_z_zT_B0_z_inv_zT_w(2 * list_size);`：函数调用语句，触发对应计算或操作。
- **[行 5565]** `		const auto itr_zT_y_list_begin = zT_y_list.begin();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5566]** `		const auto itr_zT_b0xs_list_begin = zT_b0xs_list.begin();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5567]** `		const auto itr_zT_B0_z_inv_zT_y_list_begin = zT_B0_z_inv_zT_y_list.begin();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5568]** `		const auto itr_zT_B0_z_inv_zT_b0xs_list_begin = zT_B0_z_inv_zT_b0xs_list.begin();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5569]** `		for (unsigned int i = 0; i < list_size; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5570]** `		  for (unsigned int j = 0; j < list_size; ++j)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5571]** `		    {`：作用域边界（代码块开始/结束）。
- **[行 5572]** `		      wT_z_zT_B0_z_inv_zT_w(i          , j          ) = (*std::next(itr_zT_y_list_begin            , i))`：调用 C++ 标准库工具函数/容器接口。
- **[行 5573]** `								      * (*std::next(itr_zT_B0_z_inv_zT_y_list_begin, j));`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5574]** ``：空行，用于分隔逻辑块。
- **[行 5575]** `		      wT_z_zT_B0_z_inv_zT_w(i          , j+list_size) = (*std::next(itr_zT_y_list_begin               , i))`：调用 C++ 标准库工具函数/容器接口。
- **[行 5576]** `								      * (*std::next(itr_zT_B0_z_inv_zT_b0xs_list_begin, j));`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5577]** ``：空行，用于分隔逻辑块。
- **[行 5578]** `		      wT_z_zT_B0_z_inv_zT_w(i+list_size, j          ) = (*std::next(itr_zT_b0xs_list_begin         , i))`：调用 C++ 标准库工具函数/容器接口。
- **[行 5579]** `								      * (*std::next(itr_zT_B0_z_inv_zT_y_list_begin, j));`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5580]** ``：空行，用于分隔逻辑块。
- **[行 5581]** `		      wT_z_zT_B0_z_inv_zT_w(i+list_size, j+list_size) = (*std::next(itr_zT_b0xs_list_begin            , i))`：调用 C++ 标准库工具函数/容器接口。
- **[行 5582]** `								      * (*std::next(itr_zT_B0_z_inv_zT_b0xs_list_begin, j));`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5583]** `		    }`：作用域边界（代码块开始/结束）。
- **[行 5584]** ``：空行，用于分隔逻辑块。
- **[行 5585]** `		FullMatrix<double> temp_matrix(2 * list_size);`：函数调用语句，触发对应计算或操作。
- **[行 5586]** `		M_matrix.mmult(temp_matrix, wT_z_zT_B0_z_inv_zT_w);`：函数调用语句，触发对应计算或操作。
- **[行 5587]** ``：空行，用于分隔逻辑块。
- **[行 5588]** `		FullMatrix<double> middle_matrix(IdentityMatrix(2*list_size));`：函数调用语句，触发对应计算或操作。
- **[行 5589]** `		middle_matrix.add(-1.0, temp_matrix);`：函数调用语句，触发对应计算或操作。
- **[行 5590]** ``：空行，用于分隔逻辑块。
- **[行 5591]** `		FullMatrix<double> middle_matrix_inv(2 * list_size);`：函数调用语句，触发对应计算或操作。
- **[行 5592]** `		middle_matrix_inv.invert(middle_matrix);`：函数调用语句，触发对应计算或操作。
- **[行 5593]** ``：空行，用于分隔逻辑块。
- **[行 5594]** `		middle_matrix_inv.mmult(middle_matrix, M_matrix);`：函数调用语句，触发对应计算或操作。
- **[行 5595]** ``：空行，用于分隔逻辑块。
- **[行 5596]** `		Vector<double> wT_z_zT_B0_z_inv_rhs(2 * list_size);`：函数调用语句，触发对应计算或操作。
- **[行 5597]** `		for (unsigned int i = 0; i < list_size; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5598]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 5599]** `		    wT_z_zT_B0_z_inv_rhs(i            ) = (*std::next(itr_zT_B0_z_inv_zT_y_list_begin   , i)) * rhs_vector;`：调用 C++ 标准库工具函数/容器接口。
- **[行 5600]** `		    wT_z_zT_B0_z_inv_rhs(i + list_size) = (*std::next(itr_zT_B0_z_inv_zT_b0xs_list_begin, i)) * rhs_vector;`：调用 C++ 标准库工具函数/容器接口。
- **[行 5601]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 5602]** ``：空行，用于分隔逻辑块。
- **[行 5603]** `		Vector<double> middle_matrix_wT_z_zT_B0_z_inv_rhs(2 * list_size);`：函数调用语句，触发对应计算或操作。
- **[行 5604]** `		middle_matrix.vmult(middle_matrix_wT_z_zT_B0_z_inv_rhs,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5605]** `				    wT_z_zT_B0_z_inv_rhs);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5606]** ``：空行，用于分隔逻辑块。
- **[行 5607]** `		unsigned int index = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5608]** `		for (auto itr = zT_B0_z_inv_zT_y_list.begin(); itr != zT_B0_z_inv_zT_y_list.end(); ++itr)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5609]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 5610]** `		    update_vector.add(middle_matrix_wT_z_zT_B0_z_inv_rhs(index), *itr);`：函数调用语句，触发对应计算或操作。
- **[行 5611]** `		    ++index;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5612]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 5613]** `		for (auto itr = zT_B0_z_inv_zT_b0xs_list.begin(); itr != zT_B0_z_inv_zT_b0xs_list.end(); ++itr)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5614]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 5615]** `		    update_vector.add(middle_matrix_wT_z_zT_B0_z_inv_rhs(index), *itr);`：函数调用语句，触发对应计算或操作。
- **[行 5616]** `		    ++index;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5617]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 5618]** `	      } //	if (list_size > 0)`：作用域边界（代码块开始/结束）。
- **[行 5619]** ``：空行，用于分隔逻辑块。
- **[行 5620]** `	    search_direction += update_vector;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5621]** ``：空行，用于分隔逻辑块。
- **[行 5622]** `	    m_timer.leave_subsection();`：函数调用语句，触发对应计算或操作。
- **[行 5623]** `	  } // else if (m_parameters.m_type_linear_solver == "Direct")`：作用域边界（代码块开始/结束）。
- **[行 5624]** `	else`：条件分支的兜底路径。
- **[行 5625]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 5626]** `	    AssertThrow(false, ExcMessage("Linear solver type not implemented"));`：运行期断言/检查，验证输入与状态合法性。
- **[行 5627]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 5628]** ``：空行，用于分隔逻辑块。
- **[行 5629]** `	// We don't do backtrack yet. We will make sure phasefield`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5630]** `	// remains feasible later`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5631]** `	alpha_backtrack = 1.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5632]** `	search_direction *= alpha_backtrack;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5633]** ``：空行，用于分隔逻辑块。
- **[行 5634]** `	LBFGS_update.block(m_u_dof) = solution_delta_cauchy_point.block(m_u_dof);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5635]** `	LBFGS_update.block(m_u_dof) += search_direction.block(m_u_dof);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5636]** `	LBFGS_update.block(m_u_dof) -= solution_delta.block(m_u_dof);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5637]** ``：空行，用于分隔逻辑块。
- **[行 5638]** `	for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5639]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 5640]** `	    // phasefield active constraints`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5641]** `	    if (m_active_set_phasefield(i) > 0.5)`：条件分支：根据当前状态选择执行路径。
- **[行 5642]** `	      LBFGS_update.block(m_d_dof)[i] = solution_delta_cauchy_point.block(m_d_dof)[i]`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5643]** `					     - solution_delta.block(m_d_dof)[i];`：函数调用语句，触发对应计算或操作。
- **[行 5644]** `	    else`：条件分支的兜底路径。
- **[行 5645]** `	      LBFGS_update.block(m_d_dof)[i] = solution_delta_cauchy_point.block(m_d_dof)[i]`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5646]** `					     + search_direction.block(m_d_dof)[i]`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5647]** `					     - solution_delta.block(m_d_dof)[i];`：函数调用语句，触发对应计算或操作。
- **[行 5648]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 5649]** ``：空行，用于分隔逻辑块。
- **[行 5650]** `	// make sure the phasefield solutions are feasible`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5651]** `	for(unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5652]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 5653]** `	    if (solution_delta.block(m_d_dof)[i] + LBFGS_update.block(m_d_dof)[i] < 0.0)`：条件分支：根据当前状态选择执行路径。
- **[行 5654]** `	      LBFGS_update.block(m_d_dof)[i] = -solution_delta.block(m_d_dof)[i];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5655]** ``：空行，用于分隔逻辑块。
- **[行 5656]** `	    if (  solution_delta.block(m_d_dof)[i]`：条件分支：根据当前状态选择执行路径。
- **[行 5657]** `		+ m_solution.block(m_d_dof)[i]`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5658]** `		+ LBFGS_update.block(m_d_dof)[i] > 1.0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5659]** `	      LBFGS_update.block(m_d_dof)[i] = 1.0 - m_solution.block(m_d_dof)[i]`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5660]** `						   - solution_delta.block(m_d_dof)[i];`：函数调用语句，触发对应计算或操作。
- **[行 5661]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 5662]** ``：空行，用于分隔逻辑块。
- **[行 5663]** `	m_constraints.distribute(LBFGS_update);`：函数调用语句，触发对应计算或操作。
- **[行 5664]** ``：空行，用于分隔逻辑块。
- **[行 5665]** `	// We need a line search algorithm to decide line_search_parameter`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5666]** ``：空行，用于分隔逻辑块。
- **[行 5667]** `        if(m_parameters.m_type_line_search == "StrongWolfe")`：条件分支：根据当前状态选择执行路径。
- **[行 5668]** `          {`：作用域边界（代码块开始/结束）。
- **[行 5669]** `	    const double phi_0 = calculate_energy_functional();`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 5670]** `	    const double phi_0_prime = m_system_rhs * LBFGS_update;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5671]** ``：空行，用于分隔逻辑块。
- **[行 5672]** `	    line_search_parameter = line_search_stepsize_strong_wolfe(phi_0,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 5673]** `								      phi_0_prime,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5674]** `								      LBFGS_update,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5675]** `								      solution_delta);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5676]** `          }`：作用域边界（代码块开始/结束）。
- **[行 5677]** `        else if(m_parameters.m_type_line_search == "GradientBased")`：多分支条件判断，处理备选情形。
- **[行 5678]** `          {`：作用域边界（代码块开始/结束）。
- **[行 5679]** `	    // LBFGS_r_vector is the search direction`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5680]** `	    line_search_parameter = line_search_stepsize_gradient_based(LBFGS_update,`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 5681]** `									solution_delta);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5682]** `          }`：作用域边界（代码块开始/结束）。
- **[行 5683]** `        else`：条件分支的兜底路径。
- **[行 5684]** `          {`：作用域边界（代码块开始/结束）。
- **[行 5685]** `            Assert(false, ExcMessage("An unknown line search method is called!"));`：运行期断言/检查，验证输入与状态合法性。
- **[行 5686]** `          }`：作用域边界（代码块开始/结束）。
- **[行 5687]** ``：空行，用于分隔逻辑块。
- **[行 5688]** `	LBFGS_update *= line_search_parameter;`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 5689]** ``：空行，用于分隔逻辑块。
- **[行 5690]** `        get_error_update(LBFGS_update, m_error_update);`：函数调用语句，触发对应计算或操作。
- **[行 5691]** `        if (LBFGS_iteration == 1)`：条件分支：根据当前状态选择执行路径。
- **[行 5692]** `          m_error_update_0 = m_error_update;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5693]** ``：空行，用于分隔逻辑块。
- **[行 5694]** `        m_error_update_norm = m_error_update;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5695]** `        // For three-point bending problem and the sphere inclusion problem,`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5696]** `        // we use absolute residual for convergence test`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5697]** `        if (m_parameters.m_relative_residual)`：条件分支：根据当前状态选择执行路径。
- **[行 5698]** `          m_error_update_norm.normalize(m_error_update_0);`：函数调用语句，触发对应计算或操作。
- **[行 5699]** ``：空行，用于分隔逻辑块。
- **[行 5700]** `	solution_delta += LBFGS_update;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5701]** ``：空行，用于分隔逻辑块。
- **[行 5702]** `	update_qph_incremental(solution_delta, m_solution);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 5703]** ``：空行，用于分隔逻辑块。
- **[行 5704]** `        LBFGS_y_vector = m_system_rhs;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5705]** `        LBFGS_y_vector *= -1.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5706]** `        assemble_system_rhs_BFGS_parallel(m_solution, m_system_rhs);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 5707]** `        // if we use assemble_system_rhs_BFGS_parallel, then condense() is not necessary`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5708]** `        //m_constraints.condense(m_system_rhs);`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5709]** `        LBFGS_y_vector += m_system_rhs;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5710]** ``：空行，用于分隔逻辑块。
- **[行 5711]** `        LBFGS_s_vector = LBFGS_update;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5712]** ``：空行，用于分隔逻辑块。
- **[行 5713]** `	// s_vector_list, y_vector_list, s_dot_y_list only need to discard`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5714]** `	// the front (oldest) item and add the newest item to the end at`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5715]** `	// each L-BFGS iteration`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5716]** `	double s_dot_y = LBFGS_s_vector * LBFGS_y_vector;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5717]** `	if (s_dot_y > 1.0e-16 * LBFGS_y_vector.norm_sqr())`：条件分支：根据当前状态选择执行路径。
- **[行 5718]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 5719]** `	    if (list_size >= LBFGS_m)`：条件分支：根据当前状态选择执行路径。
- **[行 5720]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 5721]** `		s_vector_list.pop_front();`：函数调用语句，触发对应计算或操作。
- **[行 5722]** `		y_vector_list.pop_front();`：函数调用语句，触发对应计算或操作。
- **[行 5723]** `		s_dot_y_list.pop_front();`：函数调用语句，触发对应计算或操作。
- **[行 5724]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 5725]** ``：空行，用于分隔逻辑块。
- **[行 5726]** `	    s_vector_list.push_back(LBFGS_s_vector);`：函数调用语句，触发对应计算或操作。
- **[行 5727]** `	    y_vector_list.push_back(LBFGS_y_vector);`：函数调用语句，触发对应计算或操作。
- **[行 5728]** `	    s_dot_y_list.push_back(s_dot_y);`：函数调用语句，触发对应计算或操作。
- **[行 5729]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 5730]** ``：空行，用于分隔逻辑块。
- **[行 5731]** `	Vector<double> solution_phasefield_total(m_solution.block(m_d_dof));`：函数调用语句，触发对应计算或操作。
- **[行 5732]** `	solution_phasefield_total += solution_delta.block(m_d_dof);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5733]** ``：空行，用于分隔逻辑块。
- **[行 5734]** `	// Since line search parameter might be less than one, we need update`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5735]** `	// the phasefield active set status`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5736]** `	// upper bound is 1.0, lower bound is the solution at the previous step.`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5737]** `	unsigned int number_active_constraint_lower_bound = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5738]** `	unsigned int number_active_constraint_upper_bound = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5739]** `	unsigned int number_active_constraint_lowerupper_bound = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5740]** ``：空行，用于分隔逻辑块。
- **[行 5741]** `	for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5742]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 5743]** `	    if (   solution_delta.block(m_d_dof)[i] == 0.0`：条件分支：根据当前状态选择执行路径。
- **[行 5744]** `		&& solution_phasefield_total[i] == 1.0)`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5745]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 5746]** `		m_active_set_phasefield(i) = 3; //lower bound overlaps with upper bound`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5747]** `		++number_active_constraint_lowerupper_bound;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5748]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 5749]** `	    else if (   solution_delta.block(m_d_dof)[i] == 0.0`：多分支条件判断，处理备选情形。
- **[行 5750]** `		     && solution_phasefield_total[i] != 1.0)`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5751]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 5752]** `		m_active_set_phasefield(i) = 1; //lower bound`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5753]** `		++number_active_constraint_lower_bound;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5754]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 5755]** `	    else if (   solution_phasefield_total[i] == 1.0`：多分支条件判断，处理备选情形。
- **[行 5756]** `		     && solution_delta.block(m_d_dof)[i] != 0.0)`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5757]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 5758]** `	        m_active_set_phasefield(i) = 2; //upper bound`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5759]** `	        ++number_active_constraint_upper_bound;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5760]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 5761]** `	    else`：条件分支的兜底路径。
- **[行 5762]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 5763]** `	        m_active_set_phasefield(i) = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5764]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 5765]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 5766]** ``：空行，用于分隔逻辑块。
- **[行 5767]** `	lower_bound_number_old = lower_bound_number_new;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5768]** `	upper_bound_number_old = upper_bound_number_new;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5769]** `	lowerupper_bound_number_old = lowerupper_bound_number_new;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5770]** ``：空行，用于分隔逻辑块。
- **[行 5771]** `	lower_bound_number_new = number_active_constraint_lower_bound;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5772]** `	upper_bound_number_new = number_active_constraint_upper_bound;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5773]** `	lowerupper_bound_number_new = number_active_constraint_lowerupper_bound;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5774]** ``：空行，用于分隔逻辑块。
- **[行 5775]** `	if (m_parameters.m_output_iteration_history)`：条件分支：根据当前状态选择执行路径。
- **[行 5776]** `          {`：作用域边界（代码块开始/结束）。
- **[行 5777]** `	    const double energy_functional = calculate_energy_functional();`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 5778]** ``：空行，用于分隔逻辑块。
- **[行 5779]** `	    m_logfile << "  | "`：写日志输出，记录当前计算状态与结果。
- **[行 5780]** `		      << std::setw(6) << number_active_constraint_lower_bound`：调用 C++ 标准库工具函数/容器接口。
- **[行 5781]** `		      << std::setw(6) << number_active_constraint_upper_bound`：调用 C++ 标准库工具函数/容器接口。
- **[行 5782]** `		      << std::setw(6) << number_active_constraint_lowerupper_bound;`：调用 C++ 标准库工具函数/容器接口。
- **[行 5783]** `            if (m_parameters.m_type_linear_solver == "CG")`：条件分支：根据当前状态选择执行路径。
- **[行 5784]** `              m_logfile << std::setw(8) << cg_iterations;`：写日志输出，记录当前计算状态与结果。
- **[行 5785]** `            else`：条件分支的兜底路径。
- **[行 5786]** `              m_logfile << std::setw(8) << "---";`：写日志输出，记录当前计算状态与结果。
- **[行 5787]** `            m_logfile << "      "`：写日志输出，记录当前计算状态与结果。
- **[行 5788]** `        	      << std::fixed << std::setprecision(3) << std::scientific << line_search_parameter`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 5789]** `		      << std::fixed << std::setprecision(6) << std::setw(1) << std::scientific`：调用 C++ 标准库工具函数/容器接口。
- **[行 5790]** `		      << "  " << energy_functional`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5791]** `		      << std::fixed << std::setprecision(3) << std::setw(1)`：调用 C++ 标准库工具函数/容器接口。
- **[行 5792]** `					<< std::scientific`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5793]** `		      << "  " << m_error_residual_norm.m_norm`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5794]** `		      << "  " << m_error_residual_norm.m_u`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5795]** `		      << "  " << m_error_residual_norm.m_d`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5796]** `		      << "  " << m_error_update_norm.m_norm`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5797]** `		      << "  " << m_error_update_norm.m_u`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5798]** `		      << "  " << m_error_update_norm.m_d`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5799]** `		      << "  " << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5800]** `          }`：作用域边界（代码块开始/结束）。
- **[行 5801]** `      } // for (; LBFGS_iteration < m_parameters.m_max_iterations_BFGS; ++LBFGS_iteration)`：作用域边界（代码块开始/结束）。
- **[行 5802]** ``：空行，用于分隔逻辑块。
- **[行 5803]** `    AssertThrow(LBFGS_iteration < m_parameters.m_max_iterations_BFGS,`：运行期断言/检查，验证输入与状态合法性。
- **[行 5804]** `                ExcMessage("No convergence in L-BFGS-B nonlinear solver!"));`：函数调用语句，触发对应计算或操作。
- **[行 5805]** `  }`：作用域边界（代码块开始/结束）。
- **[行 5806]** ``：空行，用于分隔逻辑块。
- **[行 5807]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 5808]** `  void PhaseFieldMonolithicSolve<dim>::output_results() const`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5809]** `  {`：作用域边界（代码块开始/结束）。
- **[行 5810]** `    m_timer.enter_subsection("Output results");`：函数调用语句，触发对应计算或操作。
- **[行 5811]** ``：空行，用于分隔逻辑块。
- **[行 5812]** `    DataOut<dim> data_out;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5813]** ``：空行，用于分隔逻辑块。
- **[行 5814]** `    std::vector<DataComponentInterpretation::DataComponentInterpretation>`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5815]** `      data_component_interpretation(`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5816]** `        dim, DataComponentInterpretation::component_is_part_of_vector);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5817]** ``：空行，用于分隔逻辑块。
- **[行 5818]** `    data_component_interpretation.push_back(`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5819]** `      DataComponentInterpretation::component_is_scalar);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5820]** ``：空行，用于分隔逻辑块。
- **[行 5821]** `    std::vector<std::string> solution_name(dim, "displacement");`：调用 C++ 标准库工具函数/容器接口。
- **[行 5822]** `    solution_name.emplace_back("phasefield");`：函数调用语句，触发对应计算或操作。
- **[行 5823]** ``：空行，用于分隔逻辑块。
- **[行 5824]** `    data_out.attach_dof_handler(m_dof_handler);`：函数调用语句，触发对应计算或操作。
- **[行 5825]** `    data_out.add_data_vector(m_solution,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5826]** `                             solution_name,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5827]** `                             DataOut<dim>::type_dof_data,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5828]** `                             data_component_interpretation);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5829]** ``：空行，用于分隔逻辑块。
- **[行 5830]** `    // output phasefield active set status`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5831]** `    BlockVector<double> active_set_status(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 5832]** `    active_set_status.block(m_d_dof) = m_active_set_phasefield;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5833]** `    std::vector<DataComponentInterpretation::DataComponentInterpretation>`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5834]** `      data_component_interpretation_active_set(`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5835]** `        dim+1, DataComponentInterpretation::component_is_scalar);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5836]** `    std::vector<std::string> solution_name_active_set;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5837]** `    solution_name_active_set.emplace_back("disp_x_active_set");`：函数调用语句，触发对应计算或操作。
- **[行 5838]** `    solution_name_active_set.emplace_back("disp_y_active_set");`：函数调用语句，触发对应计算或操作。
- **[行 5839]** `    if (dim ==3)`：条件分支：根据当前状态选择执行路径。
- **[行 5840]** `      solution_name_active_set.emplace_back("disp_z_active_set");`：函数调用语句，触发对应计算或操作。
- **[行 5841]** `    solution_name_active_set.emplace_back("phasefield_active_set");`：函数调用语句，触发对应计算或操作。
- **[行 5842]** `    data_out.add_data_vector(active_set_status,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5843]** `			     solution_name_active_set,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5844]** `			     DataOut<dim>::type_dof_data,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5845]** `			     data_component_interpretation_active_set);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5846]** ``：空行，用于分隔逻辑块。
- **[行 5847]** `    Vector<double> cell_material_id(m_triangulation.n_active_cells());`：函数调用语句，触发对应计算或操作。
- **[行 5848]** `    // output material ID for each cell`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5849]** `    for (const auto &cell : m_triangulation.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5850]** `      {`：作用域边界（代码块开始/结束）。
- **[行 5851]** `	cell_material_id(cell->active_cell_index()) = cell->material_id();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5852]** `      }`：作用域边界（代码块开始/结束）。
- **[行 5853]** `    data_out.add_data_vector(cell_material_id, "materialID");`：函数调用语句，触发对应计算或操作。
- **[行 5854]** ``：空行，用于分隔逻辑块。
- **[行 5855]** `    // Stress L2 projection`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5856]** `    DoFHandler<dim> stresses_dof_handler_L2(m_triangulation);`：函数调用语句，触发对应计算或操作。
- **[行 5857]** `    FE_Q<dim>     stresses_fe_L2(m_parameters.m_poly_degree); //FE_Q element is continuous`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5858]** `    stresses_dof_handler_L2.distribute_dofs(stresses_fe_L2);`：函数调用语句，触发对应计算或操作。
- **[行 5859]** `    AffineConstraints<double> constraints;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5860]** `    constraints.clear();`：函数调用语句，触发对应计算或操作。
- **[行 5861]** `    DoFTools::make_hanging_node_constraints(stresses_dof_handler_L2, constraints);`：函数调用语句，触发对应计算或操作。
- **[行 5862]** `    constraints.close();`：函数调用语句，触发对应计算或操作。
- **[行 5863]** `    std::vector<DataComponentInterpretation::DataComponentInterpretation>`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5864]** `	  data_component_interpretation_stress(1,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5865]** `					       DataComponentInterpretation::component_is_scalar);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5866]** ``：空行，用于分隔逻辑块。
- **[行 5867]** `    for (unsigned int i = 0; i < dim; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5868]** `      for (unsigned int j = i; j < dim; ++j)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5869]** `	{`：作用域边界（代码块开始/结束）。
- **[行 5870]** `	  Vector<double> stress_field_L2;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5871]** `	  stress_field_L2.reinit(stresses_dof_handler_L2.n_dofs());`：函数调用语句，触发对应计算或操作。
- **[行 5872]** ``：空行，用于分隔逻辑块。
- **[行 5873]** `	  MappingQ<dim> mapping(m_parameters.m_poly_degree + 1);`：函数调用语句，触发对应计算或操作。
- **[行 5874]** `	  VectorTools::project(mapping,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5875]** `			       stresses_dof_handler_L2,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5876]** `			       constraints,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5877]** `			       m_qf_cell,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5878]** `			       [&] (const typename DoFHandler<dim>::active_cell_iterator & cell,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5879]** `				    const unsigned int q) -> double`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5880]** `			       {`：作用域边界（代码块开始/结束）。
- **[行 5881]** `				 return m_quadrature_point_history.get_data(cell)[q]->get_cauchy_stress()[i][j];`：返回当前函数结果。
- **[行 5882]** `			       },`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5883]** `			       stress_field_L2);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5884]** ``：空行，用于分隔逻辑块。
- **[行 5885]** `	  std::string stress_name = "Cauchy_stress_" + std::to_string(i+1) + std::to_string(j+1)`：调用 C++ 标准库工具函数/容器接口。
- **[行 5886]** `				  + "_L2";`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5887]** ``：空行，用于分隔逻辑块。
- **[行 5888]** `	  data_out.add_data_vector(stresses_dof_handler_L2,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5889]** `				   stress_field_L2,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5890]** `				   stress_name,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5891]** `				   data_component_interpretation_stress);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5892]** `	}`：作用域边界（代码块开始/结束）。
- **[行 5893]** ``：空行，用于分隔逻辑块。
- **[行 5894]** `    data_out.build_patches(m_parameters.m_poly_degree);`：函数调用语句，触发对应计算或操作。
- **[行 5895]** ``：空行，用于分隔逻辑块。
- **[行 5896]** `    std::ofstream output("Solution-" + std::to_string(dim) + "d-" +`：调用 C++ 标准库工具函数/容器接口。
- **[行 5897]** `			 Utilities::int_to_string(m_time.get_timestep(),4) + ".vtu");`：函数调用语句，触发对应计算或操作。
- **[行 5898]** ``：空行，用于分隔逻辑块。
- **[行 5899]** `    data_out.write_vtu(output);`：函数调用语句，触发对应计算或操作。
- **[行 5900]** `    m_timer.leave_subsection();`：函数调用语句，触发对应计算或操作。
- **[行 5901]** `  }`：作用域边界（代码块开始/结束）。
- **[行 5902]** ``：空行，用于分隔逻辑块。
- **[行 5903]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 5904]** `  void PhaseFieldMonolithicSolve<dim>::calculate_reaction_force(unsigned int face_ID)`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 5905]** `  {`：作用域边界（代码块开始/结束）。
- **[行 5906]** `    m_timer.enter_subsection("Calculate reaction force");`：函数调用语句，触发对应计算或操作。
- **[行 5907]** ``：空行，用于分隔逻辑块。
- **[行 5908]** `    BlockVector<double>       system_rhs;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5909]** `    system_rhs.reinit(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 5910]** ``：空行，用于分隔逻辑块。
- **[行 5911]** `    Vector<double> cell_rhs(m_dofs_per_cell);`：函数调用语句，触发对应计算或操作。
- **[行 5912]** `    std::vector<types::global_dof_index> local_dof_indices(m_dofs_per_cell);`：调用 C++ 标准库工具函数/容器接口。
- **[行 5913]** ``：空行，用于分隔逻辑块。
- **[行 5914]** `    const double time_ramp = (m_time.current() / m_time.end());`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5915]** `    std::vector<Tensor<1, dim>> rhs_values(m_n_q_points);`：调用 C++ 标准库工具函数/容器接口。
- **[行 5916]** `    const UpdateFlags uf_cell(update_values | update_gradients |`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5917]** `			      update_quadrature_points | update_JxW_values);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5918]** `    const UpdateFlags uf_face(update_values | update_normal_vectors |`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5919]** `                              update_JxW_values);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5920]** ``：空行，用于分隔逻辑块。
- **[行 5921]** `    FEValues<dim> fe_values(m_fe, m_qf_cell, uf_cell);`：函数调用语句，触发对应计算或操作。
- **[行 5922]** `    FEFaceValues<dim> fe_face_values(m_fe, m_qf_face, uf_face);`：函数调用语句，触发对应计算或操作。
- **[行 5923]** ``：空行，用于分隔逻辑块。
- **[行 5924]** `    // shape function values for displacement field`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5925]** `    std::vector<std::vector<Tensor<1, dim>>>`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5926]** `      Nx(m_qf_cell.size(), std::vector<Tensor<1, dim>>(m_dofs_per_cell));`：调用 C++ 标准库工具函数/容器接口。
- **[行 5927]** `    std::vector<std::vector<Tensor<2, dim>>>`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5928]** `      grad_Nx(m_qf_cell.size(), std::vector<Tensor<2, dim>>(m_dofs_per_cell));`：调用 C++ 标准库工具函数/容器接口。
- **[行 5929]** `    std::vector<std::vector<SymmetricTensor<2, dim>>>`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5930]** `      symm_grad_Nx(m_qf_cell.size(), std::vector<SymmetricTensor<2, dim>>(m_dofs_per_cell));`：调用 C++ 标准库工具函数/容器接口。
- **[行 5931]** ``：空行，用于分隔逻辑块。
- **[行 5932]** `    for (const auto &cell : m_dof_handler.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5933]** `      {`：作用域边界（代码块开始/结束）。
- **[行 5934]** `	// if calculate_reaction_force() is defined as const, then`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5935]** `	// we also need to put a const in std::shared_ptr,`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5936]** `	// that is, std::shared_ptr<const PointHistory<dim>>`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5937]** `	const std::vector<std::shared_ptr< PointHistory<dim>>> lqph =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5938]** `	  m_quadrature_point_history.get_data(cell);`：函数调用语句，触发对应计算或操作。
- **[行 5939]** `	Assert(lqph.size() == m_n_q_points, ExcInternalError());`：运行期断言/检查，验证输入与状态合法性。
- **[行 5940]** `        cell_rhs = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5941]** `        fe_values.reinit(cell);`：函数调用语句，触发对应计算或操作。
- **[行 5942]** `        right_hand_side(fe_values.get_quadrature_points(),`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5943]** `    		        rhs_values,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5944]** `    		        m_parameters.m_x_component*time_ramp,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5945]** `    		        m_parameters.m_y_component*time_ramp,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5946]** `    		        m_parameters.m_z_component*time_ramp);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 5947]** ``：空行，用于分隔逻辑块。
- **[行 5948]** `        for (const unsigned int q_point : fe_values.quadrature_point_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5949]** `          {`：作用域边界（代码块开始/结束）。
- **[行 5950]** `            for (const unsigned int k : fe_values.dof_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5951]** `              {`：作用域边界（代码块开始/结束）。
- **[行 5952]** `                const unsigned int k_group = m_fe.system_to_base_index(k).first.first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5953]** ``：空行，用于分隔逻辑块。
- **[行 5954]** `                if (k_group == m_u_dof)`：条件分支：根据当前状态选择执行路径。
- **[行 5955]** `                  {`：作用域边界（代码块开始/结束）。
- **[行 5956]** `    		    Nx[q_point][k] = fe_values[m_u_fe].value(k, q_point);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5957]** `    		    grad_Nx[q_point][k] = fe_values[m_u_fe].gradient(k, q_point);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5958]** `    		    symm_grad_Nx[q_point][k] = symmetrize(grad_Nx[q_point][k]);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5959]** `                  }`：作用域边界（代码块开始/结束）。
- **[行 5960]** `              }`：作用域边界（代码块开始/结束）。
- **[行 5961]** `          }`：作用域边界（代码块开始/结束）。
- **[行 5962]** ``：空行，用于分隔逻辑块。
- **[行 5963]** `        for (const unsigned int q_point : fe_values.quadrature_point_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5964]** `          {`：作用域边界（代码块开始/结束）。
- **[行 5965]** `            const SymmetricTensor<2, dim> & cauchy_stress = lqph[q_point]->get_cauchy_stress();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5966]** ``：空行，用于分隔逻辑块。
- **[行 5967]** `            const std::vector<Tensor<1,dim>> & N = Nx[q_point];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5968]** `            const std::vector<SymmetricTensor<2, dim>> & symm_grad_N = symm_grad_Nx[q_point];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5969]** `            const double JxW = fe_values.JxW(q_point);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5970]** ``：空行，用于分隔逻辑块。
- **[行 5971]** `            for (const unsigned int i : fe_values.dof_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5972]** `              {`：作用域边界（代码块开始/结束）。
- **[行 5973]** `                const unsigned int i_group = m_fe.system_to_base_index(i).first.first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5974]** ``：空行，用于分隔逻辑块。
- **[行 5975]** `                if (i_group == m_u_dof)`：条件分支：根据当前状态选择执行路径。
- **[行 5976]** `                  {`：作用域边界（代码块开始/结束）。
- **[行 5977]** `                    cell_rhs(i) -= (symm_grad_N[i] * cauchy_stress) * JxW;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5978]** `    		    // contributions from the body force to right-hand side`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5979]** `    		    cell_rhs(i) += N[i] * rhs_values[q_point] * JxW;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5980]** `                  }`：作用域边界（代码块开始/结束）。
- **[行 5981]** `              }`：作用域边界（代码块开始/结束）。
- **[行 5982]** `          }`：作用域边界（代码块开始/结束）。
- **[行 5983]** ``：空行，用于分隔逻辑块。
- **[行 5984]** `        // if there is surface pressure, this surface pressure always applied to the`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5985]** `        // reference configuration`：注释行，用于说明算法背景、假设或实现意图。
- **[行 5986]** `        const unsigned int face_pressure_id = 100;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5987]** `        const double p0 = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5988]** ``：空行，用于分隔逻辑块。
- **[行 5989]** `        for (const auto &face : cell->face_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5990]** `          {`：作用域边界（代码块开始/结束）。
- **[行 5991]** `	    if (face->at_boundary() && face->boundary_id() == face_pressure_id)`：条件分支：根据当前状态选择执行路径。
- **[行 5992]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 5993]** `		fe_face_values.reinit(cell, face);`：函数调用语句，触发对应计算或操作。
- **[行 5994]** ``：空行，用于分隔逻辑块。
- **[行 5995]** `		for (const unsigned int f_q_point : fe_face_values.quadrature_point_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 5996]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 5997]** `		    const Tensor<1, dim> &N = fe_face_values.normal_vector(f_q_point);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 5998]** ``：空行，用于分隔逻辑块。
- **[行 5999]** `		    const double         pressure  = p0 * time_ramp;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6000]** `		    const Tensor<1, dim> traction  = pressure * N;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6001]** ``：空行，用于分隔逻辑块。
- **[行 6002]** `		    for (const unsigned int i : fe_values.dof_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 6003]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 6004]** `			const unsigned int i_group = m_fe.system_to_base_index(i).first.first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6005]** ``：空行，用于分隔逻辑块。
- **[行 6006]** `			if (i_group == m_u_dof)`：条件分支：根据当前状态选择执行路径。
- **[行 6007]** `			  {`：作用域边界（代码块开始/结束）。
- **[行 6008]** `			    const unsigned int component_i = m_fe.system_to_component_index(i).first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6009]** `			    const double Ni = fe_face_values.shape_value(i, f_q_point);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6010]** `			    const double JxW = fe_face_values.JxW(f_q_point);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6011]** `			    cell_rhs(i) += (Ni * traction[component_i]) * JxW;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6012]** `			  }`：作用域边界（代码块开始/结束）。
- **[行 6013]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 6014]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 6015]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 6016]** `          }`：作用域边界（代码块开始/结束）。
- **[行 6017]** ``：空行，用于分隔逻辑块。
- **[行 6018]** `        cell->get_dof_indices(local_dof_indices);`：函数调用语句，触发对应计算或操作。
- **[行 6019]** `        for (const unsigned int i : fe_values.dof_indices())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 6020]** `          system_rhs(local_dof_indices[i]) += cell_rhs(i);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6021]** `      } // for (const auto &cell : m_dof_handler.active_cell_iterators())`：作用域边界（代码块开始/结束）。
- **[行 6022]** ``：空行，用于分隔逻辑块。
- **[行 6023]** `    // The difference between the above assembled system_rhs and m_system_rhs`：注释行，用于说明算法背景、假设或实现意图。
- **[行 6024]** `    // is that m_system_rhs is condensed by the m_constraints, which zero out`：注释行，用于说明算法背景、假设或实现意图。
- **[行 6025]** `    // the rhs values associated with the constrained DOFs and modify the rhs`：注释行，用于说明算法背景、假设或实现意图。
- **[行 6026]** `    // values associated with the unconstrained DOFs.`：注释行，用于说明算法背景、假设或实现意图。
- **[行 6027]** ``：空行，用于分隔逻辑块。
- **[行 6028]** `    std::vector< types::global_dof_index > mapping;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6029]** `    std::set<types::boundary_id> boundary_ids;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6030]** `    boundary_ids.insert(face_ID);`：函数调用语句，触发对应计算或操作。
- **[行 6031]** `    DoFTools::map_dof_to_boundary_indices(m_dof_handler,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6032]** `					  boundary_ids,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6033]** `					  mapping);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6034]** ``：空行，用于分隔逻辑块。
- **[行 6035]** `    std::vector<double> reaction_force(dim, 0.0);`：调用 C++ 标准库工具函数/容器接口。
- **[行 6036]** ``：空行，用于分隔逻辑块。
- **[行 6037]** `    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 6038]** `      {`：作用域边界（代码块开始/结束）。
- **[行 6039]** `	if (mapping[i] != numbers::invalid_dof_index)`：条件分支：根据当前状态选择执行路径。
- **[行 6040]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 6041]** `	    reaction_force[i % dim] += system_rhs.block(m_u_dof)(i);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6042]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 6043]** `      }`：作用域边界（代码块开始/结束）。
- **[行 6044]** ``：空行，用于分隔逻辑块。
- **[行 6045]** `    for (unsigned int i = 0; i < dim; i++)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 6046]** `      m_logfile << "\t\tReaction force in direction " << i << " on boundary ID " << face_ID`：写日志输出，记录当前计算状态与结果。
- **[行 6047]** `                << " = "`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6048]** `		<< std::fixed << std::setprecision(3) << std::setw(1)`：调用 C++ 标准库工具函数/容器接口。
- **[行 6049]** `                << std::scientific`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6050]** `		<< reaction_force[i] << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6051]** ``：空行，用于分隔逻辑块。
- **[行 6052]** `    std::pair<double, std::vector<double>> time_force;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6053]** `    time_force.first = m_time.current();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6054]** `    time_force.second = reaction_force;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6055]** `    m_history_reaction_force.push_back(time_force);`：函数调用语句，触发对应计算或操作。
- **[行 6056]** ``：空行，用于分隔逻辑块。
- **[行 6057]** `    m_timer.leave_subsection();`：函数调用语句，触发对应计算或操作。
- **[行 6058]** `  }`：作用域边界（代码块开始/结束）。
- **[行 6059]** ``：空行，用于分隔逻辑块。
- **[行 6060]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 6061]** `  void PhaseFieldMonolithicSolve<dim>::write_history_data()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6062]** `  {`：作用域边界（代码块开始/结束）。
- **[行 6063]** `    m_logfile << "\t\tWrite history data ... \n"<<std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6064]** ``：空行，用于分隔逻辑块。
- **[行 6065]** `    std::ofstream myfile_reaction_force ("Reaction_force.hist");`：调用 C++ 标准库工具函数/容器接口。
- **[行 6066]** `    if (myfile_reaction_force.is_open())`：条件分支：根据当前状态选择执行路径。
- **[行 6067]** `    {`：作用域边界（代码块开始/结束）。
- **[行 6068]** `      myfile_reaction_force << 0.0 << "\t";`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6069]** `      if (dim == 2)`：条件分支：根据当前状态选择执行路径。
- **[行 6070]** `	myfile_reaction_force << 0.0 << "\t"`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6071]** `	       << 0.0 << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6072]** `      if (dim == 3)`：条件分支：根据当前状态选择执行路径。
- **[行 6073]** `	myfile_reaction_force << 0.0 << "\t"`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6074]** `	       << 0.0 << "\t"`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6075]** `	       << 0.0 << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6076]** ``：空行，用于分隔逻辑块。
- **[行 6077]** `      for (auto const & time_force : m_history_reaction_force)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 6078]** `	{`：作用域边界（代码块开始/结束）。
- **[行 6079]** `	  myfile_reaction_force << time_force.first << "\t";`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6080]** `	  if (dim == 2)`：条件分支：根据当前状态选择执行路径。
- **[行 6081]** `	    myfile_reaction_force << time_force.second[0] << "\t"`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6082]** `	           << time_force.second[1] << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6083]** `	  if (dim == 3)`：条件分支：根据当前状态选择执行路径。
- **[行 6084]** `	    myfile_reaction_force << time_force.second[0] << "\t"`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6085]** `	           << time_force.second[1] << "\t"`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6086]** `		   << time_force.second[2] << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6087]** `	}`：作用域边界（代码块开始/结束）。
- **[行 6088]** `      myfile_reaction_force.close();`：函数调用语句，触发对应计算或操作。
- **[行 6089]** `    }`：作用域边界（代码块开始/结束）。
- **[行 6090]** `    else`：条件分支的兜底路径。
- **[行 6091]** `      m_logfile << "Unable to open file";`：写日志输出，记录当前计算状态与结果。
- **[行 6092]** ``：空行，用于分隔逻辑块。
- **[行 6093]** `    std::ofstream myfile_energy ("Energy.hist");`：调用 C++ 标准库工具函数/容器接口。
- **[行 6094]** `    if (myfile_energy.is_open())`：条件分支：根据当前状态选择执行路径。
- **[行 6095]** `    {`：作用域边界（代码块开始/结束）。
- **[行 6096]** `      myfile_energy << std::fixed << std::setprecision(10) << std::scientific`：调用 C++ 标准库工具函数/容器接口。
- **[行 6097]** `                    << 0.0 << "\t"`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6098]** `                    << 0.0 << "\t"`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6099]** `	            << 0.0 << "\t"`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6100]** `	            << 0.0 << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6101]** ``：空行，用于分隔逻辑块。
- **[行 6102]** `      for (auto const & time_energy : m_history_energy)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 6103]** `	{`：作用域边界（代码块开始/结束）。
- **[行 6104]** `	  myfile_energy << std::fixed << std::setprecision(10) << std::scientific`：调用 C++ 标准库工具函数/容器接口。
- **[行 6105]** `	                << time_energy.first     << "\t"`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6106]** `                        << time_energy.second[0] << "\t"`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6107]** `	                << time_energy.second[1] << "\t"`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6108]** `		        << time_energy.second[2] << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6109]** `	}`：作用域边界（代码块开始/结束）。
- **[行 6110]** `      myfile_energy.close();`：函数调用语句，触发对应计算或操作。
- **[行 6111]** `    }`：作用域边界（代码块开始/结束）。
- **[行 6112]** `    else`：条件分支的兜底路径。
- **[行 6113]** `      m_logfile << "Unable to open file";`：写日志输出，记录当前计算状态与结果。
- **[行 6114]** `  }`：作用域边界（代码块开始/结束）。
- **[行 6115]** ``：空行，用于分隔逻辑块。
- **[行 6116]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 6117]** `  double PhaseFieldMonolithicSolve<dim>::calculate_energy_functional() const`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 6118]** `  {`：作用域边界（代码块开始/结束）。
- **[行 6119]** `    double energy_functional = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6120]** ``：空行，用于分隔逻辑块。
- **[行 6121]** `    FEValues<dim> fe_values(m_fe, m_qf_cell, update_JxW_values);`：函数调用语句，触发对应计算或操作。
- **[行 6122]** ``：空行，用于分隔逻辑块。
- **[行 6123]** `    for (const auto &cell : m_dof_handler.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 6124]** `      {`：作用域边界（代码块开始/结束）。
- **[行 6125]** `        fe_values.reinit(cell);`：函数调用语句，触发对应计算或操作。
- **[行 6126]** ``：空行，用于分隔逻辑块。
- **[行 6127]** `        const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6128]** `          m_quadrature_point_history.get_data(cell);`：函数调用语句，触发对应计算或操作。
- **[行 6129]** `        Assert(lqph.size() == m_n_q_points, ExcInternalError());`：运行期断言/检查，验证输入与状态合法性。
- **[行 6130]** ``：空行，用于分隔逻辑块。
- **[行 6131]** `        for (unsigned int q_point = 0; q_point < m_n_q_points; ++q_point)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 6132]** `          {`：作用域边界（代码块开始/结束）。
- **[行 6133]** `            const double JxW = fe_values.JxW(q_point);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6134]** `            energy_functional += lqph[q_point]->get_total_strain_energy() * JxW;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6135]** `            energy_functional += lqph[q_point]->get_crack_energy_dissipation() * JxW;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6136]** `          }`：作用域边界（代码块开始/结束）。
- **[行 6137]** `      }`：作用域边界（代码块开始/结束）。
- **[行 6138]** ``：空行，用于分隔逻辑块。
- **[行 6139]** `    return energy_functional;`：返回当前函数结果。
- **[行 6140]** `  }`：作用域边界（代码块开始/结束）。
- **[行 6141]** ``：空行，用于分隔逻辑块。
- **[行 6142]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 6143]** `  std::pair<double, double>`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6144]** `    PhaseFieldMonolithicSolve<dim>::calculate_total_strain_energy_and_crack_energy_dissipation() const`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 6145]** `  {`：作用域边界（代码块开始/结束）。
- **[行 6146]** `    double total_strain_energy = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6147]** `    double crack_energy_dissipation = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6148]** ``：空行，用于分隔逻辑块。
- **[行 6149]** `    FEValues<dim> fe_values(m_fe, m_qf_cell, update_JxW_values);`：函数调用语句，触发对应计算或操作。
- **[行 6150]** ``：空行，用于分隔逻辑块。
- **[行 6151]** `    for (const auto &cell : m_dof_handler.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 6152]** `      {`：作用域边界（代码块开始/结束）。
- **[行 6153]** `        fe_values.reinit(cell);`：函数调用语句，触发对应计算或操作。
- **[行 6154]** ``：空行，用于分隔逻辑块。
- **[行 6155]** `        const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6156]** `          m_quadrature_point_history.get_data(cell);`：函数调用语句，触发对应计算或操作。
- **[行 6157]** `        Assert(lqph.size() == m_n_q_points, ExcInternalError());`：运行期断言/检查，验证输入与状态合法性。
- **[行 6158]** ``：空行，用于分隔逻辑块。
- **[行 6159]** `        for (unsigned int q_point = 0; q_point < m_n_q_points; ++q_point)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 6160]** `          {`：作用域边界（代码块开始/结束）。
- **[行 6161]** `            const double JxW = fe_values.JxW(q_point);`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6162]** `            total_strain_energy += lqph[q_point]->get_total_strain_energy() * JxW;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6163]** `            crack_energy_dissipation += lqph[q_point]->get_crack_energy_dissipation() * JxW;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6164]** `          }`：作用域边界（代码块开始/结束）。
- **[行 6165]** `      }`：作用域边界（代码块开始/结束）。
- **[行 6166]** ``：空行，用于分隔逻辑块。
- **[行 6167]** `    return std::make_pair(total_strain_energy, crack_energy_dissipation);`：返回当前函数结果。
- **[行 6168]** `  }`：作用域边界（代码块开始/结束）。
- **[行 6169]** ``：空行，用于分隔逻辑块。
- **[行 6170]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 6171]** `  bool PhaseFieldMonolithicSolve<dim>::local_refine_and_solution_transfer(BlockVector<double> & solution_delta,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6172]** `									  BlockVector<double> & LBFGS_update_refine)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6173]** `  {`：作用域边界（代码块开始/结束）。
- **[行 6174]** `    // This is the solution at (n+1) obtained from the old (coarse) mesh`：注释行，用于说明算法背景、假设或实现意图。
- **[行 6175]** `    BlockVector<double> solution_next_step(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 6176]** `    solution_next_step = m_solution + solution_delta;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6177]** `    bool mesh_is_same = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6178]** `    bool cell_refine_flag = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6179]** ``：空行，用于分隔逻辑块。
- **[行 6180]** `    unsigned int material_id;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6181]** `    double length_scale;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6182]** `    double cell_length;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6183]** `    while(cell_refine_flag)`：循环推进（通常用于时间步或迭代收敛过程）。
- **[行 6184]** `      {`：作用域边界（代码块开始/结束）。
- **[行 6185]** `	cell_refine_flag = false;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6186]** ``：空行，用于分隔逻辑块。
- **[行 6187]** `	std::vector<types::global_dof_index> local_dof_indices(m_fe.dofs_per_cell);`：调用 C++ 标准库工具函数/容器接口。
- **[行 6188]** `	for (const auto &cell : m_dof_handler.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 6189]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 6190]** `	    cell->get_dof_indices(local_dof_indices);`：函数调用语句，触发对应计算或操作。
- **[行 6191]** ``：空行，用于分隔逻辑块。
- **[行 6192]** `	    for (unsigned int i = 0; i< m_fe.dofs_per_cell; ++i)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 6193]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 6194]** `		const unsigned int comp_i = m_fe.system_to_component_index(i).first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6195]** `		if (comp_i == m_d_component) //phasefield component`：条件分支：根据当前状态选择执行路径。
- **[行 6196]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 6197]** `		    if (  solution_next_step(local_dof_indices[i])`：条件分支：根据当前状态选择执行路径。
- **[行 6198]** `			> m_parameters.m_phasefield_refine_threshold )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6199]** `		      {`：作用域边界（代码块开始/结束）。
- **[行 6200]** `			material_id = cell->material_id();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6201]** `	                length_scale = m_material_data[material_id][2];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6202]** `	                if (dim == 2)`：条件分支：根据当前状态选择执行路径。
- **[行 6203]** `	                  cell_length = std::sqrt(cell->measure());`：调用 C++ 标准库工具函数/容器接口。
- **[行 6204]** `	                else`：条件分支的兜底路径。
- **[行 6205]** `	                  cell_length = std::cbrt(cell->measure());`：调用 C++ 标准库工具函数/容器接口。
- **[行 6206]** `			if (  cell_length`：条件分支：根据当前状态选择执行路径。
- **[行 6207]** `			    > length_scale * m_parameters.m_allowed_max_h_l_ratio )`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6208]** `			  {`：作用域边界（代码块开始/结束）。
- **[行 6209]** `			    if (cell->level() < m_parameters.m_max_allowed_refinement_level)`：条件分支：根据当前状态选择执行路径。
- **[行 6210]** `			      {`：作用域边界（代码块开始/结束）。
- **[行 6211]** `			        cell->set_refine_flag();`：函数调用语句，触发对应计算或操作。
- **[行 6212]** `			        break;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6213]** `			      }`：作用域边界（代码块开始/结束）。
- **[行 6214]** `			  }`：作用域边界（代码块开始/结束）。
- **[行 6215]** `		      }`：作用域边界（代码块开始/结束）。
- **[行 6216]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 6217]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 6218]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 6219]** ``：空行，用于分隔逻辑块。
- **[行 6220]** `	for (const auto &cell : m_dof_handler.active_cell_iterators())`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 6221]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 6222]** `	    if (cell->refine_flag_set())`：条件分支：根据当前状态选择执行路径。
- **[行 6223]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 6224]** `		cell_refine_flag = true;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6225]** `		break;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6226]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 6227]** `	  }`：作用域边界（代码块开始/结束）。
- **[行 6228]** ``：空行，用于分隔逻辑块。
- **[行 6229]** `	// if any cell is refined, we need to project the solution`：注释行，用于说明算法背景、假设或实现意图。
- **[行 6230]** `	// to the newly refined mesh`：注释行，用于说明算法背景、假设或实现意图。
- **[行 6231]** `	if (cell_refine_flag)`：条件分支：根据当前状态选择执行路径。
- **[行 6232]** `	  {`：作用域边界（代码块开始/结束）。
- **[行 6233]** `	    mesh_is_same = false;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6234]** ``：空行，用于分隔逻辑块。
- **[行 6235]** `	    std::vector<BlockVector<double> > old_solutions(2);`：调用 C++ 标准库工具函数/容器接口。
- **[行 6236]** `	    old_solutions[0] = solution_next_step;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6237]** `	    old_solutions[1] = m_solution;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6238]** ``：空行，用于分隔逻辑块。
- **[行 6239]** `	    m_triangulation.prepare_coarsening_and_refinement();`：函数调用语句，触发对应计算或操作。
- **[行 6240]** `	    SolutionTransfer<dim, BlockVector<double>> solution_transfer(m_dof_handler);`：函数调用语句，触发对应计算或操作。
- **[行 6241]** `	    solution_transfer.prepare_for_coarsening_and_refinement(old_solutions);`：函数调用语句，触发对应计算或操作。
- **[行 6242]** `	    m_triangulation.execute_coarsening_and_refinement();`：函数调用语句，触发对应计算或操作。
- **[行 6243]** ``：空行，用于分隔逻辑块。
- **[行 6244]** `	    setup_system();`：函数调用语句，触发对应计算或操作。
- **[行 6245]** ``：空行，用于分隔逻辑块。
- **[行 6246]** `	    std::vector<BlockVector<double>> tmp_solutions(2);`：调用 C++ 标准库工具函数/容器接口。
- **[行 6247]** `	    tmp_solutions[0].reinit(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 6248]** `	    tmp_solutions[1].reinit(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 6249]** ``：空行，用于分隔逻辑块。
- **[行 6250]** `            #  if DEAL_II_VERSION_GTE(9, 7, 0)`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6251]** `	    solution_transfer.interpolate(tmp_solutions);`：函数调用语句，触发对应计算或操作。
- **[行 6252]** `	    #  else`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6253]** `	    // If an older version of dealII is used, for example, 9.4.0, interpolate()`：注释行，用于说明算法背景、假设或实现意图。
- **[行 6254]** `            // needs to use the following interface.`：注释行，用于说明算法背景、假设或实现意图。
- **[行 6255]** `            solution_transfer.interpolate(old_solutions, tmp_solutions);`：函数调用语句，触发对应计算或操作。
- **[行 6256]** `            #  endif`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6257]** `	    solution_next_step = tmp_solutions[0];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6258]** `	    m_solution = tmp_solutions[1];`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6259]** ``：空行，用于分隔逻辑块。
- **[行 6260]** `	    // make sure the projected solutions still satisfy`：注释行，用于说明算法背景、假设或实现意图。
- **[行 6261]** `	    // hanging node constraints`：注释行，用于说明算法背景、假设或实现意图。
- **[行 6262]** `	    m_constraints.distribute(solution_next_step);`：函数调用语句，触发对应计算或操作。
- **[行 6263]** `	    m_constraints.distribute(m_solution);`：函数调用语句，触发对应计算或操作。
- **[行 6264]** `	  } // if (cell_refine_flag)`：作用域边界（代码块开始/结束）。
- **[行 6265]** `      } // while(cell_refine_flag)`：作用域边界（代码块开始/结束）。
- **[行 6266]** ``：空行，用于分隔逻辑块。
- **[行 6267]** `    // calculate field variables for newly refined cells`：注释行，用于说明算法背景、假设或实现意图。
- **[行 6268]** `    if (!mesh_is_same)`：条件分支：根据当前状态选择执行路径。
- **[行 6269]** `      {`：作用域边界（代码块开始/结束）。
- **[行 6270]** `	BlockVector<double> temp_solution_delta(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 6271]** `	BlockVector<double> temp_previous_solution(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 6272]** `	temp_solution_delta = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6273]** `	temp_previous_solution = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6274]** `	update_qph_incremental(temp_solution_delta, temp_previous_solution);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 6275]** ``：空行，用于分隔逻辑块。
- **[行 6276]** `	// initial guess for the resolve on the refined mesh`：注释行，用于说明算法背景、假设或实现意图。
- **[行 6277]** `	LBFGS_update_refine = solution_next_step - m_solution;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6278]** `      }`：作用域边界（代码块开始/结束）。
- **[行 6279]** ``：空行，用于分隔逻辑块。
- **[行 6280]** `    return mesh_is_same;`：返回当前函数结果。
- **[行 6281]** `  }`：作用域边界（代码块开始/结束）。
- **[行 6282]** ``：空行，用于分隔逻辑块。
- **[行 6283]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 6284]** `  void PhaseFieldMonolithicSolve<dim>::print_parameter_information()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6285]** `  {`：作用域边界（代码块开始/结束）。
- **[行 6286]** `    if (m_parameters.m_type_nonlinear_solver == "LBFGS")`：条件分支：根据当前状态选择执行路径。
- **[行 6287]** `      {`：作用域边界（代码块开始/结束）。
- **[行 6288]** `	m_logfile << "WARNING: this version of LBFGS does not enforce"`：写日志输出，记录当前计算状态与结果。
- **[行 6289]** `	    " phase-field irreversibility." << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6290]** `	m_logfile << "It should only be used for demonstrating the importance of"`：写日志输出，记录当前计算状态与结果。
- **[行 6291]** `	    " inequality constraints." << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6292]** `	m_logfile << "The obtained result is meaningless!" << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6293]** `      }`：作用域边界（代码块开始/结束）。
- **[行 6294]** ``：空行，用于分隔逻辑块。
- **[行 6295]** `    m_logfile << "Scenario number = " << m_parameters.m_scenario << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6296]** `    m_logfile << "Log file = " << m_parameters.m_logfile_name << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6297]** `    m_logfile << "Write iteration history to log file? = " << std::boolalpha`：写日志输出，记录当前计算状态与结果。
- **[行 6298]** `	      << m_parameters.m_output_iteration_history << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6299]** ``：空行，用于分隔逻辑块。
- **[行 6300]** `    if (dim == 2)`：条件分支：根据当前状态选择执行路径。
- **[行 6301]** `      {`：作用域边界（代码块开始/结束）。
- **[行 6302]** `	if (m_parameters.m_plane_stress)`：条件分支：根据当前状态选择执行路径。
- **[行 6303]** `	  m_logfile << "2D plane-stress case" << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6304]** `	else`：条件分支的兜底路径。
- **[行 6305]** `	  m_logfile << "2D plane-strain case" << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6306]** `      }`：作用域边界（代码块开始/结束）。
- **[行 6307]** ``：空行，用于分隔逻辑块。
- **[行 6308]** `    m_logfile << "Nonlinear solver type = " << m_parameters.m_type_nonlinear_solver << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6309]** `    m_logfile << "Line search type = " << m_parameters.m_type_line_search << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6310]** `    m_logfile << "Linear solver type = " << m_parameters.m_type_linear_solver << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6311]** ``：空行，用于分隔逻辑块。
- **[行 6312]** `    if (m_parameters.m_type_linear_solver == "CG")`：条件分支：根据当前状态选择执行路径。
- **[行 6313]** `      {`：作用域边界（代码块开始/结束）。
- **[行 6314]** `        m_logfile << "Preconditioner type for CG = " << m_parameters.m_type_preconditioner << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6315]** `        m_logfile << "Convergence tolerance for CG iterations = " << m_parameters.m_CG_tolerace << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6316]** `      }`：作用域边界（代码块开始/结束）。
- **[行 6317]** ``：空行，用于分隔逻辑块。
- **[行 6318]** `    m_logfile << "Mesh refinement strategy = " << m_parameters.m_refinement_strategy << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6319]** ``：空行，用于分隔逻辑块。
- **[行 6320]** `    if (m_parameters.m_refinement_strategy == "adaptive-refine")`：条件分支：根据当前状态选择执行路径。
- **[行 6321]** `      {`：作用域边界（代码块开始/结束）。
- **[行 6322]** `	m_logfile << "\tMaximum adaptive refinement times allowed in each step = "`：写日志输出，记录当前计算状态与结果。
- **[行 6323]** `		  << m_parameters.m_max_adaptive_refine_times << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6324]** `	m_logfile << "\tMaximum allowed cell refinement level = "`：写日志输出，记录当前计算状态与结果。
- **[行 6325]** `		  << m_parameters.m_max_allowed_refinement_level << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6326]** `	m_logfile << "\tPhasefield-based refinement threshold value = "`：写日志输出，记录当前计算状态与结果。
- **[行 6327]** `		  << m_parameters.m_phasefield_refine_threshold << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6328]** `      }`：作用域边界（代码块开始/结束）。
- **[行 6329]** ``：空行，用于分隔逻辑块。
- **[行 6330]** `    m_logfile << "L-BFGS_m = " << m_parameters.m_LBFGS_m << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6331]** `    m_logfile << "Global refinement times = " << m_parameters.m_global_refine_times << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6332]** `    m_logfile << "Local prerefinement times = " <<m_parameters. m_local_prerefine_times << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6333]** `    m_logfile << "Allowed maximum h/l ratio = " << m_parameters.m_allowed_max_h_l_ratio << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6334]** `    m_logfile << "total number of material types = " << m_parameters.m_total_material_regions << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6335]** `    m_logfile << "material data file name = " << m_parameters.m_material_file_name << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6336]** `    if (m_parameters.m_reaction_force_face_id >= 0)`：条件分支：根据当前状态选择执行路径。
- **[行 6337]** `      m_logfile << "Calculate reaction forces on Face ID = " << m_parameters.m_reaction_force_face_id << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6338]** `    else`：条件分支的兜底路径。
- **[行 6339]** `      m_logfile << "No need to calculate reaction forces." << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6340]** ``：空行，用于分隔逻辑块。
- **[行 6341]** `    if (m_parameters.m_relative_residual)`：条件分支：根据当前状态选择执行路径。
- **[行 6342]** `      m_logfile << "Relative residual for convergence." << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6343]** `    else`：条件分支的兜底路径。
- **[行 6344]** `      m_logfile << "Absolute residual for convergence." << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6345]** ``：空行，用于分隔逻辑块。
- **[行 6346]** `    m_logfile << "Body force = (" << m_parameters.m_x_component << ", "`：写日志输出，记录当前计算状态与结果。
- **[行 6347]** `                                  << m_parameters.m_y_component << ", "`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6348]** `	                          << m_parameters.m_z_component << ") (N/m^3)"`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6349]** `				  << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6350]** ``：空行，用于分隔逻辑块。
- **[行 6351]** `    m_logfile << "End time = " << m_parameters.m_end_time << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6352]** `    m_logfile << "Time data file name = " << m_parameters.m_time_file_name << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6353]** `  }`：作用域边界（代码块开始/结束）。
- **[行 6354]** ``：空行，用于分隔逻辑块。
- **[行 6355]** `  template <int dim>`：模板声明，支持维度或类型的泛化实现。
- **[行 6356]** `  void PhaseFieldMonolithicSolve<dim>::run()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6357]** `  {`：作用域边界（代码块开始/结束）。
- **[行 6358]** `    print_parameter_information();`：函数调用语句，触发对应计算或操作。
- **[行 6359]** ``：空行，用于分隔逻辑块。
- **[行 6360]** `    read_material_data(m_parameters.m_material_file_name,`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6361]** `    		       m_parameters.m_total_material_regions);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6362]** ``：空行，用于分隔逻辑块。
- **[行 6363]** `    std::vector<std::array<double, 4>> time_table;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6364]** ``：空行，用于分隔逻辑块。
- **[行 6365]** `    read_time_data(m_parameters.m_time_file_name, time_table);`：函数调用语句，触发对应计算或操作。
- **[行 6366]** ``：空行，用于分隔逻辑块。
- **[行 6367]** `    make_grid();`：函数调用语句，触发对应计算或操作。
- **[行 6368]** `    setup_system();`：函数调用语句，触发对应计算或操作。
- **[行 6369]** `    output_results();`：函数调用语句，触发对应计算或操作。
- **[行 6370]** ``：空行，用于分隔逻辑块。
- **[行 6371]** `    m_time.increment(time_table);`：函数调用语句，触发对应计算或操作。
- **[行 6372]** ``：空行，用于分隔逻辑块。
- **[行 6373]** `    while(m_time.current() < m_time.end() + m_time.get_delta_t()*1.0e-6)`：循环推进（通常用于时间步或迭代收敛过程）。
- **[行 6374]** `      {`：作用域边界（代码块开始/结束）。
- **[行 6375]** `	m_logfile << std::endl`：写日志输出，记录当前计算状态与结果。
- **[行 6376]** `		  << "Timestep " << m_time.get_timestep() << " @ " << m_time.current()`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6377]** `		  << 's' << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6378]** ``：空行，用于分隔逻辑块。
- **[行 6379]** `        bool mesh_is_same = false;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6380]** ``：空行，用于分隔逻辑块。
- **[行 6381]** `        // initial guess for the resolve on the refined mesh`：注释行，用于说明算法背景、假设或实现意图。
- **[行 6382]** `	BlockVector<double> LBFGS_update_refine(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 6383]** `	LBFGS_update_refine = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6384]** ``：空行，用于分隔逻辑块。
- **[行 6385]** `        // local adaptive mesh refinement loop`：注释行，用于说明算法背景、假设或实现意图。
- **[行 6386]** `	unsigned int adp_refine_iteration = 0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6387]** `        for (; adp_refine_iteration < m_parameters.m_max_adaptive_refine_times + 1; ++adp_refine_iteration)`：循环遍历容器/自由度/单元，执行批量计算。
- **[行 6388]** `          {`：作用域边界（代码块开始/结束）。
- **[行 6389]** `	    if (m_parameters.m_refinement_strategy == "adaptive-refine")`：条件分支：根据当前状态选择执行路径。
- **[行 6390]** `	      m_logfile << "\tAdaptive refinement-"<< adp_refine_iteration << ": " << std::endl;`：写日志输出，记录当前计算状态与结果。
- **[行 6391]** ``：空行，用于分隔逻辑块。
- **[行 6392]** `	    BlockVector<double> solution_delta(m_dofs_per_block);`：函数调用语句，触发对应计算或操作。
- **[行 6393]** `	    solution_delta = 0.0;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6394]** ``：空行，用于分隔逻辑块。
- **[行 6395]** `            if (m_parameters.m_type_nonlinear_solver == "LBFGS")`：条件分支：根据当前状态选择执行路径。
- **[行 6396]** `	      solve_nonlinear_timestep_LBFGS(solution_delta, LBFGS_update_refine);`：函数调用语句，触发对应计算或操作。
- **[行 6397]** `            else if (m_parameters.m_type_nonlinear_solver == "LBFGSB")`：多分支条件判断，处理备选情形。
- **[行 6398]** `              solve_nonlinear_timestep_LBFGS_B(solution_delta, LBFGS_update_refine);`：函数调用语句，触发对应计算或操作。
- **[行 6399]** `	    else`：条件分支的兜底路径。
- **[行 6400]** `	      AssertThrow(false, ExcMessage("Nonlinear solver type not implemented"));`：运行期断言/检查，验证输入与状态合法性。
- **[行 6401]** ``：空行，用于分隔逻辑块。
- **[行 6402]** `	    if (m_parameters.m_refinement_strategy == "adaptive-refine")`：条件分支：根据当前状态选择执行路径。
- **[行 6403]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 6404]** ``：空行，用于分隔逻辑块。
- **[行 6405]** `		if (adp_refine_iteration == m_parameters.m_max_adaptive_refine_times)`：条件分支：根据当前状态选择执行路径。
- **[行 6406]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 6407]** `		    m_solution += solution_delta;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6408]** `		    break;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6409]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 6410]** ``：空行，用于分隔逻辑块。
- **[行 6411]** `		mesh_is_same = local_refine_and_solution_transfer(solution_delta,`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6412]** `								  LBFGS_update_refine);`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6413]** ``：空行，用于分隔逻辑块。
- **[行 6414]** `		if (mesh_is_same)`：条件分支：根据当前状态选择执行路径。
- **[行 6415]** `		  {`：作用域边界（代码块开始/结束）。
- **[行 6416]** `		    m_solution += solution_delta;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6417]** `		    break;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6418]** `		  }`：作用域边界（代码块开始/结束）。
- **[行 6419]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 6420]** `	    else if (m_parameters.m_refinement_strategy == "pre-refine")`：多分支条件判断，处理备选情形。
- **[行 6421]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 6422]** `		m_solution += solution_delta;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6423]** `	        break;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6424]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 6425]** `	    else`：条件分支的兜底路径。
- **[行 6426]** `	      {`：作用域边界（代码块开始/结束）。
- **[行 6427]** `		AssertThrow(false,`：运行期断言/检查，验证输入与状态合法性。
- **[行 6428]** `		            ExcMessage("Selected mesh refinement strategy not implemented!"));`：函数调用语句，触发对应计算或操作。
- **[行 6429]** `	      }`：作用域边界（代码块开始/结束）。
- **[行 6430]** `          } // for (; adp_refine_iteration < m_parameters.m_max_adaptive_refine_times; ++adp_refine_iteration)`：作用域边界（代码块开始/结束）。
- **[行 6431]** ``：空行，用于分隔逻辑块。
- **[行 6432]** `        //AssertThrow(adp_refine_iteration < m_parameters.m_max_adaptive_refine_times,`：注释行，用于说明算法背景、假设或实现意图。
- **[行 6433]** `        //            ExcMessage("Number of local adaptive mesh refinement exceeds allowed maximum times!"));`：注释行，用于说明算法背景、假设或实现意图。
- **[行 6434]** ``：空行，用于分隔逻辑块。
- **[行 6435]** `	// output vtk files every 10 steps if there are too`：注释行，用于说明算法背景、假设或实现意图。
- **[行 6436]** `	// many time steps`：注释行，用于说明算法背景、假设或实现意图。
- **[行 6437]** `	//if (m_time.get_timestep() % 10 == 0)`：注释行，用于说明算法背景、假设或实现意图。
- **[行 6438]** `        output_results();`：函数调用语句，触发对应计算或操作。
- **[行 6439]** ``：空行，用于分隔逻辑块。
- **[行 6440]** `	double energy_functional_current = calculate_energy_functional();`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 6441]** `	m_logfile << "\t\tEnergy functional (J) = " << std::fixed << std::setprecision(10) << std::scientific`：写日志输出，记录当前计算状态与结果。
- **[行 6442]** `	          << energy_functional_current << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6443]** ``：空行，用于分隔逻辑块。
- **[行 6444]** `	std::pair<double, double> energy_pair = calculate_total_strain_energy_and_crack_energy_dissipation();`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 6445]** `	m_logfile << "\t\tTotal strain energy (J) = " << std::fixed << std::setprecision(10) << std::scientific`：写日志输出，记录当前计算状态与结果。
- **[行 6446]** `		  << energy_pair.first << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6447]** `	m_logfile << "\t\tCrack energy dissipation (J) = " << std::fixed << std::setprecision(10) << std::scientific`：写日志输出，记录当前计算状态与结果。
- **[行 6448]** `		  << energy_pair.second << std::endl;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6449]** ``：空行，用于分隔逻辑块。
- **[行 6450]** `	std::pair<double, std::array<double, 3>> time_energy;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6451]** `	time_energy.first = m_time.current();`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6452]** `	time_energy.second[0] = energy_pair.first;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6453]** `	time_energy.second[1] = energy_pair.second;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6454]** `	time_energy.second[2] = energy_pair.first + energy_pair.second;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6455]** `	m_history_energy.push_back(time_energy);`：函数调用语句，触发对应计算或操作。
- **[行 6456]** ``：空行，用于分隔逻辑块。
- **[行 6457]** `	int face_ID = m_parameters.m_reaction_force_face_id;`：变量赋值/更新，准备后续计算或保存中间结果。
- **[行 6458]** `	if (face_ID >= 0)`：条件分支：根据当前状态选择执行路径。
- **[行 6459]** `	  calculate_reaction_force(face_ID);`：调用核心数值子程序（组装/更新/线搜索/后处理计算）。
- **[行 6460]** ``：空行，用于分隔逻辑块。
- **[行 6461]** `        write_history_data();`：函数调用语句，触发对应计算或操作。
- **[行 6462]** ``：空行，用于分隔逻辑块。
- **[行 6463]** `	m_time.increment(time_table);`：函数调用语句，触发对应计算或操作。
- **[行 6464]** `      } // while(m_time.current() < m_time.end() + m_time.get_delta_t()*1.0e-6)`：作用域边界（代码块开始/结束）。
- **[行 6465]** `  }`：作用域边界（代码块开始/结束）。
- **[行 6466]** `} // namespace PhaseField`：作用域边界（代码块开始/结束）。
- **[行 6467]** ``：空行，用于分隔逻辑块。
- **[行 6468]** `int main(int argc, char* argv[])`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6469]** `{`：作用域边界（代码块开始/结束）。
- **[行 6470]** `  using namespace dealii;`：实现语句：参与当前功能块的控制、数据或计算流程。
- **[行 6471]** ``：空行，用于分隔逻辑块。
- **[行 6472]** `  if (argc != 2)`：条件分支：根据当前状态选择执行路径。
- **[行 6473]** `    AssertThrow(false,`：运行期断言/检查，验证输入与状态合法性。
- **[行 6474]** `    		ExcMessage("The number of arguments provided to the program has to be 2!"));`：函数调用语句，触发对应计算或操作。
- **[行 6475]** ``：空行，用于分隔逻辑块。
- **[行 6476]** `  const unsigned int dim = std::stoi(argv[1]);`：调用 C++ 标准库工具函数/容器接口。
- **[行 6477]** `  if (dim == 2 )`：条件分支：根据当前状态选择执行路径。
- **[行 6478]** `    {`：作用域边界（代码块开始/结束）。
- **[行 6479]** `      PhaseField::PhaseFieldMonolithicSolve<2> problem_2D("parameters.prm");`：函数调用语句，触发对应计算或操作。
- **[行 6480]** `      problem_2D.run();`：函数调用语句，触发对应计算或操作。
- **[行 6481]** `    }`：作用域边界（代码块开始/结束）。
- **[行 6482]** `  else if (dim == 3)`：多分支条件判断，处理备选情形。
- **[行 6483]** `    {`：作用域边界（代码块开始/结束）。
- **[行 6484]** `      PhaseField::PhaseFieldMonolithicSolve<3> problem_3D("parameters.prm");`：函数调用语句，触发对应计算或操作。
- **[行 6485]** `      problem_3D.run();`：函数调用语句，触发对应计算或操作。
- **[行 6486]** `    }`：作用域边界（代码块开始/结束）。
- **[行 6487]** `  else`：条件分支的兜底路径。
- **[行 6488]** `    {`：作用域边界（代码块开始/结束）。
- **[行 6489]** `      AssertThrow(false,`：运行期断言/检查，验证输入与状态合法性。
- **[行 6490]** `                  ExcMessage("Dimension has to be either 2 or 3"));`：函数调用语句，触发对应计算或操作。
- **[行 6491]** `    }`：作用域边界（代码块开始/结束）。
- **[行 6492]** ``：空行，用于分隔逻辑块。
- **[行 6493]** `  return 0;`：返回当前函数结果。
- **[行 6494]** `}`：作用域边界（代码块开始/结束）。

## 5. 全部方法/函数索引（起止行号）

> 该索引覆盖 `main.cc` 中可独立定位的函数/方法定义（按起始行排序）。

| 起止行号 | 方法/函数签名 |
|---|---|
| 103~119 | `std::vector<types::global_dof_index> get_vertex_dofs( const typename Triangulation<dim>::active_vertex_iterator &vertex, const DoFHandler<dim> &dof_handler)` |
| 138~140 | `usr_Jacobi_preconditioner::usr_Jacobi_preconditioner(const BlockSparseMatrix<double> & S) : m_system_matrix(&S)` |
| 142~149 | `void usr_Jacobi_preconditioner::vmult(BlockVector<double> & dst, const BlockVector<double> & src) const` |
| 168~170 | `usr_sparseLU_preconditioner::usr_sparseLU_preconditioner(const SparseDirectUMFPACK & matrix_factorization) : m_matrix_LU(&matrix_factorization)` |
| 172~176 | `void usr_sparseLU_preconditioner::vmult(BlockVector<double> & dst, const BlockVector<double> & src) const` |
| 198~202 | `usr_sparseILU_preconditioner::usr_sparseILU_preconditioner(const SparseILU<double> & ILU_factorization_disp, const SparseILU<double> & ILU_factorization_phasefield) : m_ILU_factorization_disp(& ILU_factorization_disp) , m_ILU_factorization_phasefield(& ILU_factorization_phasefield)` |
| 204~214 | `void usr_sparseILU_preconditioner::vmult(BlockVector<double> & dst, const BlockVector<double> & src) const` |
| 218~242 | `void right_hand_side(const std::vector<Point<dim>> &points, std::vector<Tensor<1, dim>> &  values, const double fx, const double fy, const double fz)` |
| 244~247 | `double degradation_function(const double d)` |
| 249~252 | `double degradation_function_derivative(const double d)` |
| 254~258 | `double degradation_function_2nd_order_derivative(const double d)` |
| 289~397 | `void Scenario::declare_parameters(ParameterHandler &prm)` |
| 295~425 | `Patterns::Integer(0), "Geometry, loading and boundary conditions scenario");  prm.declare_entry("Log file name", "Output.log", Patterns::FileName(Patterns::FileName::input), "Name of the file for log");  prm.declare_entry("Output iteration history", "yes", Patterns::Selection("yes\|no"), "Shall we write iteration history to the log file?");  prm.declare_entry("Plane stress", "no", Patterns::Selection("yes\|no"), "If it is 2D, is it plane-stress?");  prm.declare_entry("Nonlinear solver type", "LBFGSB", Patterns::Selection("LBFGS\|LBFGSB"), "Type of solver used to solve the nonlinear system");  prm.declare_entry("Line search type", "GradientBased", Patterns::Selection("GradientBased\|StrongWolfe"), "Type of line search method, the gradient-based method " "should be preferred since it is generally faster");  prm.declare_entry("Linear solver type", "CG", Patterns::Selection("Direct\|CG"), "Type of solver used to solve the linear system");  prm.declare_entry("Preconditioner type for CG", "ILU", Patterns::Selection("None\|Jacobi\|LU\|ILU"), "Type of preconditioner used to solve the linear system");  prm.declare_entry("CG tolerance", "1.0e-6", Patterns::Double(0.0), "Convergence tolerance of CG iterations");  prm.declare_entry("Mesh refinement strategy", "adaptive-refine", Patterns::Selection("pre-refine\|adaptive-refine"), "Mesh refinement strategy: pre-refine or adaptive-refine");  prm.declare_entry("LBFGS m", "40", Patterns::Integer(0), "Number of vectors used for LBFGS");  prm.declare_entry("Global refinement times", "0", Patterns::Integer(0), "Global refinement times (across the entire domain)");  prm.declare_entry("Local prerefinement times", "0", Patterns::Integer(0), "Local pre-refinement times (assume crack path is known a priori), " "only refine along the crack path.");  prm.declare_entry("Max adaptive refinement times", "100", Patterns::Integer(0), "Maximum number of adaptive refinement times allowed in each step");  prm.declare_entry("Max allowed refinement level", "100", Patterns::Integer(0), "Maximum allowed cell refinement level");  prm.declare_entry("Phasefield refine threshold", "0.8", Patterns::Double(), "Phasefield-based refinement threshold value");  prm.declare_entry("Allowed max hl ratio", "0.25", Patterns::Double(), "Allowed maximum ratio between mesh size h and length scale l");  prm.declare_entry("Material regions", "1", Patterns::Integer(0), "Number of material regions");  prm.declare_entry("Material data file", "1", Patterns::FileName(Patterns::FileName::input), "Material data file");  prm.declare_entry("Reaction force face ID", "1", Patterns::Integer(), "Face id where reaction forces should be calculated " "(negative integer means not to calculate reaction force)"); } prm.leave_subsection(); }  void Scenario::parse_parameters(ParameterHandler &prm)` |
| 300~425 | `Patterns::FileName(Patterns::FileName::input), "Name of the file for log");  prm.declare_entry("Output iteration history", "yes", Patterns::Selection("yes\|no"), "Shall we write iteration history to the log file?");  prm.declare_entry("Plane stress", "no", Patterns::Selection("yes\|no"), "If it is 2D, is it plane-stress?");  prm.declare_entry("Nonlinear solver type", "LBFGSB", Patterns::Selection("LBFGS\|LBFGSB"), "Type of solver used to solve the nonlinear system");  prm.declare_entry("Line search type", "GradientBased", Patterns::Selection("GradientBased\|StrongWolfe"), "Type of line search method, the gradient-based method " "should be preferred since it is generally faster");  prm.declare_entry("Linear solver type", "CG", Patterns::Selection("Direct\|CG"), "Type of solver used to solve the linear system");  prm.declare_entry("Preconditioner type for CG", "ILU", Patterns::Selection("None\|Jacobi\|LU\|ILU"), "Type of preconditioner used to solve the linear system");  prm.declare_entry("CG tolerance", "1.0e-6", Patterns::Double(0.0), "Convergence tolerance of CG iterations");  prm.declare_entry("Mesh refinement strategy", "adaptive-refine", Patterns::Selection("pre-refine\|adaptive-refine"), "Mesh refinement strategy: pre-refine or adaptive-refine");  prm.declare_entry("LBFGS m", "40", Patterns::Integer(0), "Number of vectors used for LBFGS");  prm.declare_entry("Global refinement times", "0", Patterns::Integer(0), "Global refinement times (across the entire domain)");  prm.declare_entry("Local prerefinement times", "0", Patterns::Integer(0), "Local pre-refinement times (assume crack path is known a priori), " "only refine along the crack path.");  prm.declare_entry("Max adaptive refinement times", "100", Patterns::Integer(0), "Maximum number of adaptive refinement times allowed in each step");  prm.declare_entry("Max allowed refinement level", "100", Patterns::Integer(0), "Maximum allowed cell refinement level");  prm.declare_entry("Phasefield refine threshold", "0.8", Patterns::Double(), "Phasefield-based refinement threshold value");  prm.declare_entry("Allowed max hl ratio", "0.25", Patterns::Double(), "Allowed maximum ratio between mesh size h and length scale l");  prm.declare_entry("Material regions", "1", Patterns::Integer(0), "Number of material regions");  prm.declare_entry("Material data file", "1", Patterns::FileName(Patterns::FileName::input), "Material data file");  prm.declare_entry("Reaction force face ID", "1", Patterns::Integer(), "Face id where reaction forces should be calculated " "(negative integer means not to calculate reaction force)"); } prm.leave_subsection(); }  void Scenario::parse_parameters(ParameterHandler &prm)` |
| 305~425 | `Patterns::Selection("yes\|no"), "Shall we write iteration history to the log file?");  prm.declare_entry("Plane stress", "no", Patterns::Selection("yes\|no"), "If it is 2D, is it plane-stress?");  prm.declare_entry("Nonlinear solver type", "LBFGSB", Patterns::Selection("LBFGS\|LBFGSB"), "Type of solver used to solve the nonlinear system");  prm.declare_entry("Line search type", "GradientBased", Patterns::Selection("GradientBased\|StrongWolfe"), "Type of line search method, the gradient-based method " "should be preferred since it is generally faster");  prm.declare_entry("Linear solver type", "CG", Patterns::Selection("Direct\|CG"), "Type of solver used to solve the linear system");  prm.declare_entry("Preconditioner type for CG", "ILU", Patterns::Selection("None\|Jacobi\|LU\|ILU"), "Type of preconditioner used to solve the linear system");  prm.declare_entry("CG tolerance", "1.0e-6", Patterns::Double(0.0), "Convergence tolerance of CG iterations");  prm.declare_entry("Mesh refinement strategy", "adaptive-refine", Patterns::Selection("pre-refine\|adaptive-refine"), "Mesh refinement strategy: pre-refine or adaptive-refine");  prm.declare_entry("LBFGS m", "40", Patterns::Integer(0), "Number of vectors used for LBFGS");  prm.declare_entry("Global refinement times", "0", Patterns::Integer(0), "Global refinement times (across the entire domain)");  prm.declare_entry("Local prerefinement times", "0", Patterns::Integer(0), "Local pre-refinement times (assume crack path is known a priori), " "only refine along the crack path.");  prm.declare_entry("Max adaptive refinement times", "100", Patterns::Integer(0), "Maximum number of adaptive refinement times allowed in each step");  prm.declare_entry("Max allowed refinement level", "100", Patterns::Integer(0), "Maximum allowed cell refinement level");  prm.declare_entry("Phasefield refine threshold", "0.8", Patterns::Double(), "Phasefield-based refinement threshold value");  prm.declare_entry("Allowed max hl ratio", "0.25", Patterns::Double(), "Allowed maximum ratio between mesh size h and length scale l");  prm.declare_entry("Material regions", "1", Patterns::Integer(0), "Number of material regions");  prm.declare_entry("Material data file", "1", Patterns::FileName(Patterns::FileName::input), "Material data file");  prm.declare_entry("Reaction force face ID", "1", Patterns::Integer(), "Face id where reaction forces should be calculated " "(negative integer means not to calculate reaction force)"); } prm.leave_subsection(); }  void Scenario::parse_parameters(ParameterHandler &prm)` |
| 310~425 | `Patterns::Selection("yes\|no"), "If it is 2D, is it plane-stress?");  prm.declare_entry("Nonlinear solver type", "LBFGSB", Patterns::Selection("LBFGS\|LBFGSB"), "Type of solver used to solve the nonlinear system");  prm.declare_entry("Line search type", "GradientBased", Patterns::Selection("GradientBased\|StrongWolfe"), "Type of line search method, the gradient-based method " "should be preferred since it is generally faster");  prm.declare_entry("Linear solver type", "CG", Patterns::Selection("Direct\|CG"), "Type of solver used to solve the linear system");  prm.declare_entry("Preconditioner type for CG", "ILU", Patterns::Selection("None\|Jacobi\|LU\|ILU"), "Type of preconditioner used to solve the linear system");  prm.declare_entry("CG tolerance", "1.0e-6", Patterns::Double(0.0), "Convergence tolerance of CG iterations");  prm.declare_entry("Mesh refinement strategy", "adaptive-refine", Patterns::Selection("pre-refine\|adaptive-refine"), "Mesh refinement strategy: pre-refine or adaptive-refine");  prm.declare_entry("LBFGS m", "40", Patterns::Integer(0), "Number of vectors used for LBFGS");  prm.declare_entry("Global refinement times", "0", Patterns::Integer(0), "Global refinement times (across the entire domain)");  prm.declare_entry("Local prerefinement times", "0", Patterns::Integer(0), "Local pre-refinement times (assume crack path is known a priori), " "only refine along the crack path.");  prm.declare_entry("Max adaptive refinement times", "100", Patterns::Integer(0), "Maximum number of adaptive refinement times allowed in each step");  prm.declare_entry("Max allowed refinement level", "100", Patterns::Integer(0), "Maximum allowed cell refinement level");  prm.declare_entry("Phasefield refine threshold", "0.8", Patterns::Double(), "Phasefield-based refinement threshold value");  prm.declare_entry("Allowed max hl ratio", "0.25", Patterns::Double(), "Allowed maximum ratio between mesh size h and length scale l");  prm.declare_entry("Material regions", "1", Patterns::Integer(0), "Number of material regions");  prm.declare_entry("Material data file", "1", Patterns::FileName(Patterns::FileName::input), "Material data file");  prm.declare_entry("Reaction force face ID", "1", Patterns::Integer(), "Face id where reaction forces should be calculated " "(negative integer means not to calculate reaction force)"); } prm.leave_subsection(); }  void Scenario::parse_parameters(ParameterHandler &prm)` |
| 315~425 | `Patterns::Selection("LBFGS\|LBFGSB"), "Type of solver used to solve the nonlinear system");  prm.declare_entry("Line search type", "GradientBased", Patterns::Selection("GradientBased\|StrongWolfe"), "Type of line search method, the gradient-based method " "should be preferred since it is generally faster");  prm.declare_entry("Linear solver type", "CG", Patterns::Selection("Direct\|CG"), "Type of solver used to solve the linear system");  prm.declare_entry("Preconditioner type for CG", "ILU", Patterns::Selection("None\|Jacobi\|LU\|ILU"), "Type of preconditioner used to solve the linear system");  prm.declare_entry("CG tolerance", "1.0e-6", Patterns::Double(0.0), "Convergence tolerance of CG iterations");  prm.declare_entry("Mesh refinement strategy", "adaptive-refine", Patterns::Selection("pre-refine\|adaptive-refine"), "Mesh refinement strategy: pre-refine or adaptive-refine");  prm.declare_entry("LBFGS m", "40", Patterns::Integer(0), "Number of vectors used for LBFGS");  prm.declare_entry("Global refinement times", "0", Patterns::Integer(0), "Global refinement times (across the entire domain)");  prm.declare_entry("Local prerefinement times", "0", Patterns::Integer(0), "Local pre-refinement times (assume crack path is known a priori), " "only refine along the crack path.");  prm.declare_entry("Max adaptive refinement times", "100", Patterns::Integer(0), "Maximum number of adaptive refinement times allowed in each step");  prm.declare_entry("Max allowed refinement level", "100", Patterns::Integer(0), "Maximum allowed cell refinement level");  prm.declare_entry("Phasefield refine threshold", "0.8", Patterns::Double(), "Phasefield-based refinement threshold value");  prm.declare_entry("Allowed max hl ratio", "0.25", Patterns::Double(), "Allowed maximum ratio between mesh size h and length scale l");  prm.declare_entry("Material regions", "1", Patterns::Integer(0), "Number of material regions");  prm.declare_entry("Material data file", "1", Patterns::FileName(Patterns::FileName::input), "Material data file");  prm.declare_entry("Reaction force face ID", "1", Patterns::Integer(), "Face id where reaction forces should be calculated " "(negative integer means not to calculate reaction force)"); } prm.leave_subsection(); }  void Scenario::parse_parameters(ParameterHandler &prm)` |
| 320~425 | `Patterns::Selection("GradientBased\|StrongWolfe"), "Type of line search method, the gradient-based method " "should be preferred since it is generally faster");  prm.declare_entry("Linear solver type", "CG", Patterns::Selection("Direct\|CG"), "Type of solver used to solve the linear system");  prm.declare_entry("Preconditioner type for CG", "ILU", Patterns::Selection("None\|Jacobi\|LU\|ILU"), "Type of preconditioner used to solve the linear system");  prm.declare_entry("CG tolerance", "1.0e-6", Patterns::Double(0.0), "Convergence tolerance of CG iterations");  prm.declare_entry("Mesh refinement strategy", "adaptive-refine", Patterns::Selection("pre-refine\|adaptive-refine"), "Mesh refinement strategy: pre-refine or adaptive-refine");  prm.declare_entry("LBFGS m", "40", Patterns::Integer(0), "Number of vectors used for LBFGS");  prm.declare_entry("Global refinement times", "0", Patterns::Integer(0), "Global refinement times (across the entire domain)");  prm.declare_entry("Local prerefinement times", "0", Patterns::Integer(0), "Local pre-refinement times (assume crack path is known a priori), " "only refine along the crack path.");  prm.declare_entry("Max adaptive refinement times", "100", Patterns::Integer(0), "Maximum number of adaptive refinement times allowed in each step");  prm.declare_entry("Max allowed refinement level", "100", Patterns::Integer(0), "Maximum allowed cell refinement level");  prm.declare_entry("Phasefield refine threshold", "0.8", Patterns::Double(), "Phasefield-based refinement threshold value");  prm.declare_entry("Allowed max hl ratio", "0.25", Patterns::Double(), "Allowed maximum ratio between mesh size h and length scale l");  prm.declare_entry("Material regions", "1", Patterns::Integer(0), "Number of material regions");  prm.declare_entry("Material data file", "1", Patterns::FileName(Patterns::FileName::input), "Material data file");  prm.declare_entry("Reaction force face ID", "1", Patterns::Integer(), "Face id where reaction forces should be calculated " "(negative integer means not to calculate reaction force)"); } prm.leave_subsection(); }  void Scenario::parse_parameters(ParameterHandler &prm)` |
| 326~425 | `Patterns::Selection("Direct\|CG"), "Type of solver used to solve the linear system");  prm.declare_entry("Preconditioner type for CG", "ILU", Patterns::Selection("None\|Jacobi\|LU\|ILU"), "Type of preconditioner used to solve the linear system");  prm.declare_entry("CG tolerance", "1.0e-6", Patterns::Double(0.0), "Convergence tolerance of CG iterations");  prm.declare_entry("Mesh refinement strategy", "adaptive-refine", Patterns::Selection("pre-refine\|adaptive-refine"), "Mesh refinement strategy: pre-refine or adaptive-refine");  prm.declare_entry("LBFGS m", "40", Patterns::Integer(0), "Number of vectors used for LBFGS");  prm.declare_entry("Global refinement times", "0", Patterns::Integer(0), "Global refinement times (across the entire domain)");  prm.declare_entry("Local prerefinement times", "0", Patterns::Integer(0), "Local pre-refinement times (assume crack path is known a priori), " "only refine along the crack path.");  prm.declare_entry("Max adaptive refinement times", "100", Patterns::Integer(0), "Maximum number of adaptive refinement times allowed in each step");  prm.declare_entry("Max allowed refinement level", "100", Patterns::Integer(0), "Maximum allowed cell refinement level");  prm.declare_entry("Phasefield refine threshold", "0.8", Patterns::Double(), "Phasefield-based refinement threshold value");  prm.declare_entry("Allowed max hl ratio", "0.25", Patterns::Double(), "Allowed maximum ratio between mesh size h and length scale l");  prm.declare_entry("Material regions", "1", Patterns::Integer(0), "Number of material regions");  prm.declare_entry("Material data file", "1", Patterns::FileName(Patterns::FileName::input), "Material data file");  prm.declare_entry("Reaction force face ID", "1", Patterns::Integer(), "Face id where reaction forces should be calculated " "(negative integer means not to calculate reaction force)"); } prm.leave_subsection(); }  void Scenario::parse_parameters(ParameterHandler &prm)` |
| 331~425 | `Patterns::Selection("None\|Jacobi\|LU\|ILU"), "Type of preconditioner used to solve the linear system");  prm.declare_entry("CG tolerance", "1.0e-6", Patterns::Double(0.0), "Convergence tolerance of CG iterations");  prm.declare_entry("Mesh refinement strategy", "adaptive-refine", Patterns::Selection("pre-refine\|adaptive-refine"), "Mesh refinement strategy: pre-refine or adaptive-refine");  prm.declare_entry("LBFGS m", "40", Patterns::Integer(0), "Number of vectors used for LBFGS");  prm.declare_entry("Global refinement times", "0", Patterns::Integer(0), "Global refinement times (across the entire domain)");  prm.declare_entry("Local prerefinement times", "0", Patterns::Integer(0), "Local pre-refinement times (assume crack path is known a priori), " "only refine along the crack path.");  prm.declare_entry("Max adaptive refinement times", "100", Patterns::Integer(0), "Maximum number of adaptive refinement times allowed in each step");  prm.declare_entry("Max allowed refinement level", "100", Patterns::Integer(0), "Maximum allowed cell refinement level");  prm.declare_entry("Phasefield refine threshold", "0.8", Patterns::Double(), "Phasefield-based refinement threshold value");  prm.declare_entry("Allowed max hl ratio", "0.25", Patterns::Double(), "Allowed maximum ratio between mesh size h and length scale l");  prm.declare_entry("Material regions", "1", Patterns::Integer(0), "Number of material regions");  prm.declare_entry("Material data file", "1", Patterns::FileName(Patterns::FileName::input), "Material data file");  prm.declare_entry("Reaction force face ID", "1", Patterns::Integer(), "Face id where reaction forces should be calculated " "(negative integer means not to calculate reaction force)"); } prm.leave_subsection(); }  void Scenario::parse_parameters(ParameterHandler &prm)` |
| 336~425 | `Patterns::Double(0.0), "Convergence tolerance of CG iterations");  prm.declare_entry("Mesh refinement strategy", "adaptive-refine", Patterns::Selection("pre-refine\|adaptive-refine"), "Mesh refinement strategy: pre-refine or adaptive-refine");  prm.declare_entry("LBFGS m", "40", Patterns::Integer(0), "Number of vectors used for LBFGS");  prm.declare_entry("Global refinement times", "0", Patterns::Integer(0), "Global refinement times (across the entire domain)");  prm.declare_entry("Local prerefinement times", "0", Patterns::Integer(0), "Local pre-refinement times (assume crack path is known a priori), " "only refine along the crack path.");  prm.declare_entry("Max adaptive refinement times", "100", Patterns::Integer(0), "Maximum number of adaptive refinement times allowed in each step");  prm.declare_entry("Max allowed refinement level", "100", Patterns::Integer(0), "Maximum allowed cell refinement level");  prm.declare_entry("Phasefield refine threshold", "0.8", Patterns::Double(), "Phasefield-based refinement threshold value");  prm.declare_entry("Allowed max hl ratio", "0.25", Patterns::Double(), "Allowed maximum ratio between mesh size h and length scale l");  prm.declare_entry("Material regions", "1", Patterns::Integer(0), "Number of material regions");  prm.declare_entry("Material data file", "1", Patterns::FileName(Patterns::FileName::input), "Material data file");  prm.declare_entry("Reaction force face ID", "1", Patterns::Integer(), "Face id where reaction forces should be calculated " "(negative integer means not to calculate reaction force)"); } prm.leave_subsection(); }  void Scenario::parse_parameters(ParameterHandler &prm)` |
| 341~425 | `Patterns::Selection("pre-refine\|adaptive-refine"), "Mesh refinement strategy: pre-refine or adaptive-refine");  prm.declare_entry("LBFGS m", "40", Patterns::Integer(0), "Number of vectors used for LBFGS");  prm.declare_entry("Global refinement times", "0", Patterns::Integer(0), "Global refinement times (across the entire domain)");  prm.declare_entry("Local prerefinement times", "0", Patterns::Integer(0), "Local pre-refinement times (assume crack path is known a priori), " "only refine along the crack path.");  prm.declare_entry("Max adaptive refinement times", "100", Patterns::Integer(0), "Maximum number of adaptive refinement times allowed in each step");  prm.declare_entry("Max allowed refinement level", "100", Patterns::Integer(0), "Maximum allowed cell refinement level");  prm.declare_entry("Phasefield refine threshold", "0.8", Patterns::Double(), "Phasefield-based refinement threshold value");  prm.declare_entry("Allowed max hl ratio", "0.25", Patterns::Double(), "Allowed maximum ratio between mesh size h and length scale l");  prm.declare_entry("Material regions", "1", Patterns::Integer(0), "Number of material regions");  prm.declare_entry("Material data file", "1", Patterns::FileName(Patterns::FileName::input), "Material data file");  prm.declare_entry("Reaction force face ID", "1", Patterns::Integer(), "Face id where reaction forces should be calculated " "(negative integer means not to calculate reaction force)"); } prm.leave_subsection(); }  void Scenario::parse_parameters(ParameterHandler &prm)` |
| 346~425 | `Patterns::Integer(0), "Number of vectors used for LBFGS");  prm.declare_entry("Global refinement times", "0", Patterns::Integer(0), "Global refinement times (across the entire domain)");  prm.declare_entry("Local prerefinement times", "0", Patterns::Integer(0), "Local pre-refinement times (assume crack path is known a priori), " "only refine along the crack path.");  prm.declare_entry("Max adaptive refinement times", "100", Patterns::Integer(0), "Maximum number of adaptive refinement times allowed in each step");  prm.declare_entry("Max allowed refinement level", "100", Patterns::Integer(0), "Maximum allowed cell refinement level");  prm.declare_entry("Phasefield refine threshold", "0.8", Patterns::Double(), "Phasefield-based refinement threshold value");  prm.declare_entry("Allowed max hl ratio", "0.25", Patterns::Double(), "Allowed maximum ratio between mesh size h and length scale l");  prm.declare_entry("Material regions", "1", Patterns::Integer(0), "Number of material regions");  prm.declare_entry("Material data file", "1", Patterns::FileName(Patterns::FileName::input), "Material data file");  prm.declare_entry("Reaction force face ID", "1", Patterns::Integer(), "Face id where reaction forces should be calculated " "(negative integer means not to calculate reaction force)"); } prm.leave_subsection(); }  void Scenario::parse_parameters(ParameterHandler &prm)` |
| 351~425 | `Patterns::Integer(0), "Global refinement times (across the entire domain)");  prm.declare_entry("Local prerefinement times", "0", Patterns::Integer(0), "Local pre-refinement times (assume crack path is known a priori), " "only refine along the crack path.");  prm.declare_entry("Max adaptive refinement times", "100", Patterns::Integer(0), "Maximum number of adaptive refinement times allowed in each step");  prm.declare_entry("Max allowed refinement level", "100", Patterns::Integer(0), "Maximum allowed cell refinement level");  prm.declare_entry("Phasefield refine threshold", "0.8", Patterns::Double(), "Phasefield-based refinement threshold value");  prm.declare_entry("Allowed max hl ratio", "0.25", Patterns::Double(), "Allowed maximum ratio between mesh size h and length scale l");  prm.declare_entry("Material regions", "1", Patterns::Integer(0), "Number of material regions");  prm.declare_entry("Material data file", "1", Patterns::FileName(Patterns::FileName::input), "Material data file");  prm.declare_entry("Reaction force face ID", "1", Patterns::Integer(), "Face id where reaction forces should be calculated " "(negative integer means not to calculate reaction force)"); } prm.leave_subsection(); }  void Scenario::parse_parameters(ParameterHandler &prm)` |
| 356~425 | `Patterns::Integer(0), "Local pre-refinement times (assume crack path is known a priori), " "only refine along the crack path.");  prm.declare_entry("Max adaptive refinement times", "100", Patterns::Integer(0), "Maximum number of adaptive refinement times allowed in each step");  prm.declare_entry("Max allowed refinement level", "100", Patterns::Integer(0), "Maximum allowed cell refinement level");  prm.declare_entry("Phasefield refine threshold", "0.8", Patterns::Double(), "Phasefield-based refinement threshold value");  prm.declare_entry("Allowed max hl ratio", "0.25", Patterns::Double(), "Allowed maximum ratio between mesh size h and length scale l");  prm.declare_entry("Material regions", "1", Patterns::Integer(0), "Number of material regions");  prm.declare_entry("Material data file", "1", Patterns::FileName(Patterns::FileName::input), "Material data file");  prm.declare_entry("Reaction force face ID", "1", Patterns::Integer(), "Face id where reaction forces should be calculated " "(negative integer means not to calculate reaction force)"); } prm.leave_subsection(); }  void Scenario::parse_parameters(ParameterHandler &prm)` |
| 362~425 | `Patterns::Integer(0), "Maximum number of adaptive refinement times allowed in each step");  prm.declare_entry("Max allowed refinement level", "100", Patterns::Integer(0), "Maximum allowed cell refinement level");  prm.declare_entry("Phasefield refine threshold", "0.8", Patterns::Double(), "Phasefield-based refinement threshold value");  prm.declare_entry("Allowed max hl ratio", "0.25", Patterns::Double(), "Allowed maximum ratio between mesh size h and length scale l");  prm.declare_entry("Material regions", "1", Patterns::Integer(0), "Number of material regions");  prm.declare_entry("Material data file", "1", Patterns::FileName(Patterns::FileName::input), "Material data file");  prm.declare_entry("Reaction force face ID", "1", Patterns::Integer(), "Face id where reaction forces should be calculated " "(negative integer means not to calculate reaction force)"); } prm.leave_subsection(); }  void Scenario::parse_parameters(ParameterHandler &prm)` |
| 367~425 | `Patterns::Integer(0), "Maximum allowed cell refinement level");  prm.declare_entry("Phasefield refine threshold", "0.8", Patterns::Double(), "Phasefield-based refinement threshold value");  prm.declare_entry("Allowed max hl ratio", "0.25", Patterns::Double(), "Allowed maximum ratio between mesh size h and length scale l");  prm.declare_entry("Material regions", "1", Patterns::Integer(0), "Number of material regions");  prm.declare_entry("Material data file", "1", Patterns::FileName(Patterns::FileName::input), "Material data file");  prm.declare_entry("Reaction force face ID", "1", Patterns::Integer(), "Face id where reaction forces should be calculated " "(negative integer means not to calculate reaction force)"); } prm.leave_subsection(); }  void Scenario::parse_parameters(ParameterHandler &prm)` |
| 372~425 | `Patterns::Double(), "Phasefield-based refinement threshold value");  prm.declare_entry("Allowed max hl ratio", "0.25", Patterns::Double(), "Allowed maximum ratio between mesh size h and length scale l");  prm.declare_entry("Material regions", "1", Patterns::Integer(0), "Number of material regions");  prm.declare_entry("Material data file", "1", Patterns::FileName(Patterns::FileName::input), "Material data file");  prm.declare_entry("Reaction force face ID", "1", Patterns::Integer(), "Face id where reaction forces should be calculated " "(negative integer means not to calculate reaction force)"); } prm.leave_subsection(); }  void Scenario::parse_parameters(ParameterHandler &prm)` |
| 377~425 | `Patterns::Double(), "Allowed maximum ratio between mesh size h and length scale l");  prm.declare_entry("Material regions", "1", Patterns::Integer(0), "Number of material regions");  prm.declare_entry("Material data file", "1", Patterns::FileName(Patterns::FileName::input), "Material data file");  prm.declare_entry("Reaction force face ID", "1", Patterns::Integer(), "Face id where reaction forces should be calculated " "(negative integer means not to calculate reaction force)"); } prm.leave_subsection(); }  void Scenario::parse_parameters(ParameterHandler &prm)` |
| 382~425 | `Patterns::Integer(0), "Number of material regions");  prm.declare_entry("Material data file", "1", Patterns::FileName(Patterns::FileName::input), "Material data file");  prm.declare_entry("Reaction force face ID", "1", Patterns::Integer(), "Face id where reaction forces should be calculated " "(negative integer means not to calculate reaction force)"); } prm.leave_subsection(); }  void Scenario::parse_parameters(ParameterHandler &prm)` |
| 387~425 | `Patterns::FileName(Patterns::FileName::input), "Material data file");  prm.declare_entry("Reaction force face ID", "1", Patterns::Integer(), "Face id where reaction forces should be calculated " "(negative integer means not to calculate reaction force)"); } prm.leave_subsection(); }  void Scenario::parse_parameters(ParameterHandler &prm)` |
| 392~425 | `Patterns::Integer(), "Face id where reaction forces should be calculated " "(negative integer means not to calculate reaction force)"); } prm.leave_subsection(); }  void Scenario::parse_parameters(ParameterHandler &prm)` |
| 399~425 | `void Scenario::parse_parameters(ParameterHandler &prm)` |
| 438~453 | `void FESystem::declare_parameters(ParameterHandler &prm)` |
| 444~463 | `Patterns::Integer(0), "Phase field polynomial order");  prm.declare_entry("Quadrature order", "2", Patterns::Integer(0), "Gauss quadrature order"); } prm.leave_subsection(); }  void FESystem::parse_parameters(ParameterHandler &prm)` |
| 449~463 | `Patterns::Integer(0), "Gauss quadrature order"); } prm.leave_subsection(); }  void FESystem::parse_parameters(ParameterHandler &prm)` |
| 455~463 | `void FESystem::parse_parameters(ParameterHandler &prm)` |
| 477~497 | `void BodyForce::declare_parameters(ParameterHandler &prm)` |
| 483~508 | `Patterns::Double(), "Body force x-component (N/m^3)");  prm.declare_entry("Body force y component", "0.0", Patterns::Double(), "Body force y-component (N/m^3)");  prm.declare_entry("Body force z component", "0.0", Patterns::Double(), "Body force z-component (N/m^3)"); } prm.leave_subsection(); }  void BodyForce::parse_parameters(ParameterHandler &prm)` |
| 488~508 | `Patterns::Double(), "Body force y-component (N/m^3)");  prm.declare_entry("Body force z component", "0.0", Patterns::Double(), "Body force z-component (N/m^3)"); } prm.leave_subsection(); }  void BodyForce::parse_parameters(ParameterHandler &prm)` |
| 493~508 | `Patterns::Double(), "Body force z-component (N/m^3)"); } prm.leave_subsection(); }  void BodyForce::parse_parameters(ParameterHandler &prm)` |
| 499~508 | `void BodyForce::parse_parameters(ParameterHandler &prm)` |
| 525~560 | `void NonlinearSolver::declare_parameters(ParameterHandler &prm)` |
| 531~575 | `Patterns::Integer(0), "Number of BFGS iterations allowed");  prm.declare_entry("Relative residual", "yes", Patterns::Selection("yes\|no"), "Shall we use relative residual for convergence?");  prm.declare_entry("Tolerance displacement residual", "1.0e-9", Patterns::Double(0.0), "Displacement residual tolerance");  prm.declare_entry("Tolerance phasefield residual", "1.0e-9", Patterns::Double(0.0), "Phasefield residual tolerance");  prm.declare_entry("Tolerance displacement increment", "1.0e-9", Patterns::Double(0.0), "Displacement increment tolerance");  prm.declare_entry("Tolerance phasefield increment", "1.0e-9", Patterns::Double(0.0), "Phasefield increment tolerance"); } prm.leave_subsection(); }  void NonlinearSolver::parse_parameters(ParameterHandler &prm)` |
| 536~575 | `Patterns::Selection("yes\|no"), "Shall we use relative residual for convergence?");  prm.declare_entry("Tolerance displacement residual", "1.0e-9", Patterns::Double(0.0), "Displacement residual tolerance");  prm.declare_entry("Tolerance phasefield residual", "1.0e-9", Patterns::Double(0.0), "Phasefield residual tolerance");  prm.declare_entry("Tolerance displacement increment", "1.0e-9", Patterns::Double(0.0), "Displacement increment tolerance");  prm.declare_entry("Tolerance phasefield increment", "1.0e-9", Patterns::Double(0.0), "Phasefield increment tolerance"); } prm.leave_subsection(); }  void NonlinearSolver::parse_parameters(ParameterHandler &prm)` |
| 541~575 | `Patterns::Double(0.0), "Displacement residual tolerance");  prm.declare_entry("Tolerance phasefield residual", "1.0e-9", Patterns::Double(0.0), "Phasefield residual tolerance");  prm.declare_entry("Tolerance displacement increment", "1.0e-9", Patterns::Double(0.0), "Displacement increment tolerance");  prm.declare_entry("Tolerance phasefield increment", "1.0e-9", Patterns::Double(0.0), "Phasefield increment tolerance"); } prm.leave_subsection(); }  void NonlinearSolver::parse_parameters(ParameterHandler &prm)` |
| 546~575 | `Patterns::Double(0.0), "Phasefield residual tolerance");  prm.declare_entry("Tolerance displacement increment", "1.0e-9", Patterns::Double(0.0), "Displacement increment tolerance");  prm.declare_entry("Tolerance phasefield increment", "1.0e-9", Patterns::Double(0.0), "Phasefield increment tolerance"); } prm.leave_subsection(); }  void NonlinearSolver::parse_parameters(ParameterHandler &prm)` |
| 551~575 | `Patterns::Double(0.0), "Displacement increment tolerance");  prm.declare_entry("Tolerance phasefield increment", "1.0e-9", Patterns::Double(0.0), "Phasefield increment tolerance"); } prm.leave_subsection(); }  void NonlinearSolver::parse_parameters(ParameterHandler &prm)` |
| 556~575 | `Patterns::Double(0.0), "Phasefield increment tolerance"); } prm.leave_subsection(); }  void NonlinearSolver::parse_parameters(ParameterHandler &prm)` |
| 562~575 | `void NonlinearSolver::parse_parameters(ParameterHandler &prm)` |
| 587~599 | `void TimeInfo::declare_parameters(ParameterHandler &prm)` |
| 595~609 | `Patterns::FileName(Patterns::FileName::input), "Time data file"); } prm.leave_subsection(); }  void TimeInfo::parse_parameters(ParameterHandler &prm)` |
| 601~609 | `void TimeInfo::parse_parameters(ParameterHandler &prm)` |
| 624~630 | `AllParameters::AllParameters(const std::string &input_file)` |
| 632~639 | `void AllParameters::declare_parameters(ParameterHandler &prm)` |
| 641~648 | `void AllParameters::parse_parameters(ParameterHandler &prm)` |
| 840~1000 | `usr_spectrum_decomposition::positive_negative_projectors(eigenvalues, eigenvectors, projector_positive, projector_negative);  SymmetricTensor<2, dim> stress_positive, stress_negative; const double degradation = degradation_function(m_phase_field_value) + m_residual_k; const double I_1 = trace(m_strain);  // 2D plane strain and 3D cases double my_lambda = m_lame_lambda;  // 2D plane stress case if (    dim == 2 && m_plane_stress) my_lambda = 2 * m_lame_mu * m_lame_lambda / (m_lame_lambda + 2 * m_lame_mu);  stress_positive = my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1) * Physics::Elasticity::StandardTensors<dim>::I + 2 * m_lame_mu * strain_positive; stress_negative = my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1) * Physics::Elasticity::StandardTensors<dim>::I + 2 * m_lame_mu * strain_negative;  m_stress = degradation * stress_positive + stress_negative; m_stress_positive = stress_positive;  SymmetricTensor<4, dim> C_positive, C_negative; C_positive = my_lambda * usr_spectrum_decomposition::heaviside_function(I_1) * Physics::Elasticity::StandardTensors<dim>::IxI + 2 * m_lame_mu * projector_positive; C_negative = my_lambda * usr_spectrum_decomposition::heaviside_function(-I_1) * Physics::Elasticity::StandardTensors<dim>::IxI + 2 * m_lame_mu * projector_negative; m_mechanical_C = degradation * C_positive + C_negative;  m_strain_energy_positive = 0.5 * my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1) * usr_spectrum_decomposition::positive_ramp_function(I_1) + m_lame_mu * strain_positive * strain_positive;  m_strain_energy_negative = 0.5 * my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1) * usr_spectrum_decomposition::negative_ramp_function(I_1) + m_lame_mu * strain_negative * strain_negative;  m_strain_energy_total = degradation * m_strain_energy_positive + m_strain_energy_negative;  m_crack_energy_dissipation = m_gc * (  0.5 / m_length_scale * m_phase_field_value * m_phase_field_value + 0.5 * m_length_scale * m_grad_phasefield * m_grad_phasefield) // the term due to viscosity regularization + (m_phase_field_value - phase_field_value_previous_step) * (m_phase_field_value - phase_field_value_previous_step) * 0.5 * m_eta / delta_time; //(void)delta_time; //(void)phase_field_value_previous_step; }  template <int dim> class PointHistory` |
| 877~1000 | `* usr_spectrum_decomposition::positive_ramp_function(I_1) + m_lame_mu * strain_positive * strain_positive;  m_strain_energy_negative = 0.5 * my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1) * usr_spectrum_decomposition::negative_ramp_function(I_1) + m_lame_mu * strain_negative * strain_negative;  m_strain_energy_total = degradation * m_strain_energy_positive + m_strain_energy_negative;  m_crack_energy_dissipation = m_gc * (  0.5 / m_length_scale * m_phase_field_value * m_phase_field_value + 0.5 * m_length_scale * m_grad_phasefield * m_grad_phasefield) // the term due to viscosity regularization + (m_phase_field_value - phase_field_value_previous_step) * (m_phase_field_value - phase_field_value_previous_step) * 0.5 * m_eta / delta_time; //(void)delta_time; //(void)phase_field_value_previous_step; }  template <int dim> class PointHistory` |
| 881~1000 | `* usr_spectrum_decomposition::negative_ramp_function(I_1) + m_lame_mu * strain_negative * strain_negative;  m_strain_energy_total = degradation * m_strain_energy_positive + m_strain_energy_negative;  m_crack_energy_dissipation = m_gc * (  0.5 / m_length_scale * m_phase_field_value * m_phase_field_value + 0.5 * m_length_scale * m_grad_phasefield * m_grad_phasefield) // the term due to viscosity regularization + (m_phase_field_value - phase_field_value_previous_step) * (m_phase_field_value - phase_field_value_previous_step) * 0.5 * m_eta / delta_time; //(void)delta_time; //(void)phase_field_value_previous_step; }  template <int dim> class PointHistory` |
| 1260~1298 | `void PhaseFieldMonolithicSolve<dim>::zT_B0_z(const BlockVector<double> & z, BlockSparseMatrix<double> & B0_matrix)` |
| 1301~1331 | `void PhaseFieldMonolithicSolve<dim>::z_x_vector(const BlockVector<double> & z, const BlockVector<double> & src_vector, BlockVector<double> & target_vector)` |
| 1334~1360 | `void PhaseFieldMonolithicSolve<dim>::zT_x_vector(const BlockVector<double> & z, const BlockVector<double> & src_vector, BlockVector<double> & target_vector)` |
| 1363~1376 | `double PhaseFieldMonolithicSolve<dim>::ebT_x_B0_x_v(const unsigned int b, const BlockSparseMatrix<double> & B0_matrix, const BlockVector<double> & v)` |
| 1379~1395 | `void PhaseFieldMonolithicSolve<dim>::point_projection(BlockVector<double> & solution_delta)` |
| 1401~1443 | `PhaseFieldMonolithicSolve<dim>::calculate_break_points(const BlockVector<double> & solution_delta, const BlockVector<double> & gradient_g, BlockVector<double> & gradient_d)` |
| 1446~1457 | `void PhaseFieldMonolithicSolve<dim>::get_error_residual(Errors &error_residual)` |
| 1460~1500 | `void PhaseFieldMonolithicSolve<dim>::get_error_residual_LBFGSB(Errors &error_residual, const BlockVector<double> & solution_delta)` |
| 1503~1514 | `void PhaseFieldMonolithicSolve<dim>::get_error_update(const BlockVector<double> &newton_update, Errors & error_update)` |
| 1517~1571 | `void PhaseFieldMonolithicSolve<dim>::read_material_data(const std::string &data_file, const unsigned int total_material_regions)` |
| 1574~1617 | `void PhaseFieldMonolithicSolve<dim>::read_time_data(const std::string &data_file, std::vector<std::array<double, 4>> & time_table)` |
| 1620~1667 | `void PhaseFieldMonolithicSolve<dim>::setup_qph()` |
| 1670~1676 | `BlockVector<double> PhaseFieldMonolithicSolve<dim>::get_total_solution( const BlockVector<double> &solution_delta) const` |
| 1680~1717 | `PhaseFieldMonolithicSolve<dim>::update_qph_incremental(const BlockVector<double> &solution_delta, const BlockVector<double> &solution_old)` |
| 1708~1724 | `WorkStream::run( m_dof_handler.begin_active(), m_dof_handler.end(), worker, copier, scratch_data_UQPH, per_task_data_UQPH);  m_timer.leave_subsection(); }  template <int dim> struct PhaseFieldMonolithicSolve<dim>::PerTaskData_UQPH` |
| 1785~1817 | `void PhaseFieldMonolithicSolve<dim>::update_qph_incremental_one_cell( const typename DoFHandler<dim>::active_cell_iterator &cell, ScratchData_UQPH & scratch, PerTaskData_UQPH & /*data*/)` |
| 2035~2054 | `PhaseFieldMonolithicSolve<dim>::PhaseFieldMonolithicSolve(const std::string &input_file) : m_parameters(input_file) , m_triangulation(Triangulation<dim>::maximum_smoothing) , m_time(m_parameters.m_end_time) , m_logfile(m_parameters.m_logfile_name) , m_timer(m_logfile, TimerOutput::summary, TimerOutput::wall_times) , m_dof_handler(m_triangulation) , m_fe(FE_Q<dim>(m_parameters.m_poly_degree), dim, // displacement FE_Q<dim>(m_parameters.m_poly_degree), 1)   // phasefield , m_dofs_per_cell(m_fe.n_dofs_per_cell()) , m_u_fe(m_first_u_component) , m_d_fe(m_d_component) , m_dofs_per_block(m_n_blocks) , m_qf_cell(m_parameters.m_quad_order) , m_qf_face(m_parameters.m_quad_order) , m_n_q_points(m_qf_cell.size()) , m_vol_reference(0.0)` |
| 2057~2097 | `void PhaseFieldMonolithicSolve<dim>::make_grid()` |
| 2100~2185 | `void PhaseFieldMonolithicSolve<dim>::make_grid_case_1()` |
| 2165~2175 | `&& std::fabs(cell->center()[0] - 0.5) < 0.05)` |
| 2189~2277 | `void PhaseFieldMonolithicSolve<dim>::make_grid_case_2()` |
| 2280~2365 | `void PhaseFieldMonolithicSolve<dim>::make_grid_case_3()` |
| 2345~2355 | `&& std::fabs(cell->center()[1] - 0.5) < 0.025 )` |
| 2368~2462 | `void PhaseFieldMonolithicSolve<dim>::make_grid_case_4()` |
| 2465~2560 | `void PhaseFieldMonolithicSolve<dim>::make_grid_case_5()` |
| 2540~2550 | `&& std::fabs(cell->center()[1] - 0.4) < 0.075 )` |
| 2563~2684 | `void PhaseFieldMonolithicSolve<dim>::make_grid_case_6()` |
| 2579~2606 | `GridGenerator::hyper_shell( tria_outer, Point<dim>(), 0.5, std::sqrt(dim), 2 * dim);  Triangulation<dim> tmp_triangulation;  GridGenerator::merge_triangulations(tria_inner, tria_outer, tmp_triangulation);  tmp_triangulation.reset_all_manifolds(); tmp_triangulation.set_all_manifold_ids(0);  for (const auto &cell : tmp_triangulation.cell_iterators())` |
| 2629~2649 | `GridGenerator::create_triangulation_with_removed_cells(tmp_triangulation, cells_to_remove, m_triangulation);  for (const auto &cell : m_triangulation.active_cell_iterators()) for (const auto &face : cell->face_iterators())` |
| 2687~2818 | `void PhaseFieldMonolithicSolve<dim>::make_grid_case_7()` |
| 2703~2715 | `GridGenerator::hyper_shell( tria_outer, Point<dim>(), 0.49, std::sqrt(dim)*0.5, 2 * dim);  Triangulation<dim> cube1; GridGenerator::hyper_rectangle(cube1, Point<dim>(0, 0, 0.5), Point<dim>(1, 1, 1.5)); Triangulation<dim> cube2; GridGenerator::hyper_rectangle(cube2, Point<dim>(0, 0.5, -0.5), Point<dim>(1, 1.5, 0.5)); Triangulation<dim> cube3; GridGenerator::hyper_rectangle(cube3, Point<dim>(0.5, -0.5, -0.5), Point<dim>(1.5, 0.5, 0.5));  Triangulation<dim> tmp_triangulation; GridGenerator::merge_triangulations(` |
| 2714~2715 | `GridGenerator::merge_triangulations(` |
| 2763~2783 | `GridGenerator::create_triangulation_with_removed_cells(tmp_triangulation, cells_to_remove, m_triangulation);  for (const auto &cell : m_triangulation.active_cell_iterators()) for (const auto &face : cell->face_iterators())` |
| 2822~2967 | `void PhaseFieldMonolithicSolve<dim>::make_grid_case_8()` |
| 2838~2850 | `GridGenerator::hyper_shell( tria_outer, Point<dim>(), 0.49, std::sqrt(dim)*0.5, 2 * dim);  Triangulation<dim> cube1; GridGenerator::hyper_rectangle(cube1, Point<dim>(0, 0, 0.5), Point<dim>(1, 1, 1.5)); Triangulation<dim> cube2; GridGenerator::hyper_rectangle(cube2, Point<dim>(0, 0.5, -0.5), Point<dim>(1, 1.5, 0.5)); Triangulation<dim> cube3; GridGenerator::hyper_rectangle(cube3, Point<dim>(0.5, -0.5, -0.5), Point<dim>(1.5, 0.5, 0.5));  Triangulation<dim> tmp_triangulation; GridGenerator::merge_triangulations(` |
| 2849~2850 | `GridGenerator::merge_triangulations(` |
| 2886~2910 | `&& std::fabs(cell->center()[2] - 0.5625) < 0.05 && std::fabs(cell->center()[0] - 0.0) < 0.2) cell->set_material_id(1);  if (    std::fabs(cell->center()[1] - 0.0) < 0.2 && std::fabs(cell->center()[2] - 0.5) < 0.1 && std::fabs(cell->center()[0] - 0.75) < 0.05) cell->set_material_id(1); }  std::set<typename Triangulation< dim >::active_cell_iterator > cells_to_remove;  for (const auto &cell : tmp_triangulation.active_cell_iterators())` |
| 2887~2910 | `&& std::fabs(cell->center()[0] - 0.0) < 0.2) cell->set_material_id(1);  if (    std::fabs(cell->center()[1] - 0.0) < 0.2 && std::fabs(cell->center()[2] - 0.5) < 0.1 && std::fabs(cell->center()[0] - 0.75) < 0.05) cell->set_material_id(1); }  std::set<typename Triangulation< dim >::active_cell_iterator > cells_to_remove;  for (const auto &cell : tmp_triangulation.active_cell_iterators())` |
| 2891~2910 | `&& std::fabs(cell->center()[2] - 0.5) < 0.1 && std::fabs(cell->center()[0] - 0.75) < 0.05) cell->set_material_id(1); }  std::set<typename Triangulation< dim >::active_cell_iterator > cells_to_remove;  for (const auto &cell : tmp_triangulation.active_cell_iterators())` |
| 2892~2910 | `&& std::fabs(cell->center()[0] - 0.75) < 0.05) cell->set_material_id(1); }  std::set<typename Triangulation< dim >::active_cell_iterator > cells_to_remove;  for (const auto &cell : tmp_triangulation.active_cell_iterators())` |
| 2912~2932 | `GridGenerator::create_triangulation_with_removed_cells(tmp_triangulation, cells_to_remove, m_triangulation);  for (const auto &cell : m_triangulation.active_cell_iterators()) for (const auto &face : cell->face_iterators())` |
| 2970~3055 | `void PhaseFieldMonolithicSolve<dim>::make_grid_case_9()` |
| 3035~3045 | `&& std::fabs(cell->center()[1] - 250) < 10.0 )` |
| 3058~3148 | `void PhaseFieldMonolithicSolve<dim>::make_grid_case_10()` |
| 3128~3138 | `&& std::fabs(cell->center()[1] - 250) < 10.0 )` |
| 3152~3285 | `void PhaseFieldMonolithicSolve<dim>::make_grid_case_11()` |
| 3179~3214 | `GridGenerator::subdivided_hyper_rectangle(triangulation_2d, repetitions, point1, point2 );  typename Triangulation<2>::vertex_iterator vertex_ptr; vertex_ptr = triangulation_2d.begin_active_vertex(); while (vertex_ptr != triangulation_2d.end_vertex())` |
| 3234~3264 | `GridGenerator::create_triangulation_with_removed_cells(tmp_triangulation, cells_to_remove, m_triangulation);  if (m_parameters.m_refinement_strategy == "adaptive-refine")` |
| 3288~3354 | `void PhaseFieldMonolithicSolve<dim>::setup_system()` |
| 3339~3638 | `DoFTools::make_sparsity_pattern( m_dof_handler, coupling, dsp, m_constraints, false); m_sparsity_pattern.copy_from(dsp); }  m_tangent_matrix.reinit(m_sparsity_pattern);  m_system_rhs.reinit(m_dofs_per_block); m_solution.reinit(m_dofs_per_block);  m_active_set_phasefield.reinit(m_dofs_per_block[m_d_dof]);  setup_qph();  m_timer.leave_subsection(); }  template <int dim> void PhaseFieldMonolithicSolve<dim>::make_constraints(const unsigned int it_nr)` |
| 3357~3638 | `void PhaseFieldMonolithicSolve<dim>::make_constraints(const unsigned int it_nr)` |
| 3383~3400 | `VectorTools::interpolate_boundary_values(m_dof_handler, boundary_id_bottom_surface, Functions::ZeroFunction<dim>(m_n_components), m_constraints, m_fe.component_mask(y_displacement));  typename Triangulation<dim>::active_vertex_iterator vertex_itr; vertex_itr = m_triangulation.begin_active_vertex(); std::vector<types::global_dof_index> node_xy(m_fe.dofs_per_vertex);  for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)` |
| 3409~3457 | `VectorTools::interpolate_boundary_values(m_dof_handler, boundary_id_top_surface, Functions::ZeroFunction<dim>(m_n_components), m_constraints, m_fe.component_mask(x_displacement)); */ const double time_inc = m_time.get_delta_t(); double disp_magnitude = m_time.get_magnitude(); VectorTools::interpolate_boundary_values(m_dof_handler, boundary_id_top_surface, Functions::ConstantFunction<dim>( disp_magnitude*time_inc, m_n_components), m_constraints, m_fe.component_mask(y_displacement)); } else if (   m_parameters.m_scenario == 2 \|\| m_parameters.m_scenario == 4)` |
| 3417~3457 | `VectorTools::interpolate_boundary_values(m_dof_handler, boundary_id_top_surface, Functions::ConstantFunction<dim>( disp_magnitude*time_inc, m_n_components), m_constraints, m_fe.component_mask(y_displacement)); } else if (   m_parameters.m_scenario == 2 \|\| m_parameters.m_scenario == 4)` |
| 3429~3501 | `VectorTools::interpolate_boundary_values(m_dof_handler, boundary_id_bottom_surface, Functions::ZeroFunction<dim>(m_n_components), m_constraints, m_fe.component_mask(displacements));  const int boundary_id_top_surface = 1; VectorTools::interpolate_boundary_values(m_dof_handler, boundary_id_top_surface, Functions::ZeroFunction<dim>(m_n_components), m_constraints, m_fe.component_mask(y_displacement));  const double time_inc = m_time.get_delta_t(); double disp_magnitude = m_time.get_magnitude(); VectorTools::interpolate_boundary_values(m_dof_handler, boundary_id_top_surface, Functions::ConstantFunction<dim>( disp_magnitude*time_inc, m_n_components), m_constraints, m_fe.component_mask(x_displacement));  const int boundary_id_side_surfaces = 2; VectorTools::interpolate_boundary_values(m_dof_handler, boundary_id_side_surfaces, Functions::ZeroFunction<dim>(m_n_components), m_constraints, m_fe.component_mask(y_displacement)); } else if (m_parameters.m_scenario == 5)` |
| 3436~3501 | `VectorTools::interpolate_boundary_values(m_dof_handler, boundary_id_top_surface, Functions::ZeroFunction<dim>(m_n_components), m_constraints, m_fe.component_mask(y_displacement));  const double time_inc = m_time.get_delta_t(); double disp_magnitude = m_time.get_magnitude(); VectorTools::interpolate_boundary_values(m_dof_handler, boundary_id_top_surface, Functions::ConstantFunction<dim>( disp_magnitude*time_inc, m_n_components), m_constraints, m_fe.component_mask(x_displacement));  const int boundary_id_side_surfaces = 2; VectorTools::interpolate_boundary_values(m_dof_handler, boundary_id_side_surfaces, Functions::ZeroFunction<dim>(m_n_components), m_constraints, m_fe.component_mask(y_displacement)); } else if (m_parameters.m_scenario == 5)` |
| 3444~3501 | `VectorTools::interpolate_boundary_values(m_dof_handler, boundary_id_top_surface, Functions::ConstantFunction<dim>( disp_magnitude*time_inc, m_n_components), m_constraints, m_fe.component_mask(x_displacement));  const int boundary_id_side_surfaces = 2; VectorTools::interpolate_boundary_values(m_dof_handler, boundary_id_side_surfaces, Functions::ZeroFunction<dim>(m_n_components), m_constraints, m_fe.component_mask(y_displacement)); } else if (m_parameters.m_scenario == 5)` |
| 3452~3501 | `VectorTools::interpolate_boundary_values(m_dof_handler, boundary_id_side_surfaces, Functions::ZeroFunction<dim>(m_n_components), m_constraints, m_fe.component_mask(y_displacement)); } else if (m_parameters.m_scenario == 5)` |
| 3507~3564 | `VectorTools::interpolate_boundary_values(m_dof_handler, x0_surface, Functions::ZeroFunction<dim>(m_n_components), m_constraints, m_fe.component_mask(x_displacement)); const int y0_surface = 1; VectorTools::interpolate_boundary_values(m_dof_handler, y0_surface, Functions::ZeroFunction<dim>(m_n_components), m_constraints, m_fe.component_mask(y_displacement)); const int z0_surface = 2; VectorTools::interpolate_boundary_values(m_dof_handler, z0_surface, Functions::ZeroFunction<dim>(m_n_components), m_constraints, m_fe.component_mask(z_displacement));  const int z1_surface = 3; const double time_inc = m_time.get_delta_t(); double disp_magnitude = 1.0; VectorTools::interpolate_boundary_values(m_dof_handler, z1_surface, Functions::ConstantFunction<dim>( disp_magnitude*time_inc, m_n_components), m_constraints, m_fe.component_mask(z_displacement)); } else if (   m_parameters.m_scenario == 9 \|\| m_parameters.m_scenario == 10)` |
| 3513~3564 | `VectorTools::interpolate_boundary_values(m_dof_handler, y0_surface, Functions::ZeroFunction<dim>(m_n_components), m_constraints, m_fe.component_mask(y_displacement)); const int z0_surface = 2; VectorTools::interpolate_boundary_values(m_dof_handler, z0_surface, Functions::ZeroFunction<dim>(m_n_components), m_constraints, m_fe.component_mask(z_displacement));  const int z1_surface = 3; const double time_inc = m_time.get_delta_t(); double disp_magnitude = 1.0; VectorTools::interpolate_boundary_values(m_dof_handler, z1_surface, Functions::ConstantFunction<dim>( disp_magnitude*time_inc, m_n_components), m_constraints, m_fe.component_mask(z_displacement)); } else if (   m_parameters.m_scenario == 9 \|\| m_parameters.m_scenario == 10)` |
| 3519~3564 | `VectorTools::interpolate_boundary_values(m_dof_handler, z0_surface, Functions::ZeroFunction<dim>(m_n_components), m_constraints, m_fe.component_mask(z_displacement));  const int z1_surface = 3; const double time_inc = m_time.get_delta_t(); double disp_magnitude = 1.0; VectorTools::interpolate_boundary_values(m_dof_handler, z1_surface, Functions::ConstantFunction<dim>( disp_magnitude*time_inc, m_n_components), m_constraints, m_fe.component_mask(z_displacement)); } else if (   m_parameters.m_scenario == 9 \|\| m_parameters.m_scenario == 10)` |
| 3528~3564 | `VectorTools::interpolate_boundary_values(m_dof_handler, z1_surface, Functions::ConstantFunction<dim>( disp_magnitude*time_inc, m_n_components), m_constraints, m_fe.component_mask(z_displacement)); } else if (   m_parameters.m_scenario == 9 \|\| m_parameters.m_scenario == 10)` |
| 3540~3563 | `VectorTools::interpolate_boundary_values(m_dof_handler, boundary_id_bottom_surface, Functions::ZeroFunction<dim>(m_n_components), m_constraints, m_fe.component_mask(displacements));  typename Triangulation<dim>::active_vertex_iterator vertex_itr; vertex_itr = m_triangulation.begin_active_vertex(); std::vector<types::global_dof_index> node_disp_control(m_fe.dofs_per_vertex);  for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)` |
| 3569~3620 | `VectorTools::interpolate_boundary_values(m_dof_handler, boundary_id_right_surface, Functions::ZeroFunction<dim>(m_n_components), m_constraints, m_fe.component_mask(displacements));  // Dirichlet B,C. left surface const int boundary_id_left_surface = 1; VectorTools::interpolate_boundary_values(m_dof_handler, boundary_id_left_surface, Functions::ZeroFunction<dim>(m_n_components), m_constraints, m_fe.component_mask(x_displacement));  typename Triangulation<dim>::active_vertex_iterator vertex_itr; vertex_itr = m_triangulation.begin_active_vertex(); std::vector<types::global_dof_index> node_rotate(m_fe.dofs_per_vertex); double node_dist = 0.0; double disp_mag = 0.0; double angle_theta = 0.0; double disp_y = 0; double disp_z = 0;  for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)` |
| 3577~3620 | `VectorTools::interpolate_boundary_values(m_dof_handler, boundary_id_left_surface, Functions::ZeroFunction<dim>(m_n_components), m_constraints, m_fe.component_mask(x_displacement));  typename Triangulation<dim>::active_vertex_iterator vertex_itr; vertex_itr = m_triangulation.begin_active_vertex(); std::vector<types::global_dof_index> node_rotate(m_fe.dofs_per_vertex); double node_dist = 0.0; double disp_mag = 0.0; double angle_theta = 0.0; double disp_y = 0; double disp_z = 0;  for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)` |
| 3641~3678 | `void PhaseFieldMonolithicSolve<dim>::assemble_system_B0(const BlockVector<double> & solution_old)` |
| 3670~3719 | `WorkStream::run( m_dof_handler.active_cell_iterators(), worker, copier, scratch_data, per_task_data);  m_timer.leave_subsection(); }  template <int dim> void PhaseFieldMonolithicSolve<dim>::assemble_system_rhs_BFGS_parallel(const BlockVector<double> & solution_old, BlockVector<double> & system_rhs)` |
| 3681~3719 | `void PhaseFieldMonolithicSolve<dim>::assemble_system_rhs_BFGS_parallel(const BlockVector<double> & solution_old, BlockVector<double> & system_rhs)` |
| 3711~3857 | `WorkStream::run( m_dof_handler.active_cell_iterators(), worker, copier, scratch_data, per_task_data);  m_timer.leave_subsection(); }  template <int dim> void PhaseFieldMonolithicSolve<dim>::assemble_system_rhs_BFGS_one_cell( const typename DoFHandler<dim>::active_cell_iterator &cell, ScratchData_ASM_RHS_BFGS & scratch, PerTaskData_ASM_RHS_BFGS & data) const` |
| 3722~3857 | `void PhaseFieldMonolithicSolve<dim>::assemble_system_rhs_BFGS_one_cell( const typename DoFHandler<dim>::active_cell_iterator &cell, ScratchData_ASM_RHS_BFGS & scratch, PerTaskData_ASM_RHS_BFGS & data) const` |
| 3860~3959 | `void PhaseFieldMonolithicSolve<dim>::assemble_system_B0_one_cell( const typename DoFHandler<dim>::active_cell_iterator &cell, ScratchData_ASM & scratch, PerTaskData_ASM & data) const` |
| 3962~4124 | `void PhaseFieldMonolithicSolve<dim>::assemble_system_rhs_BFGS(const BlockVector<double> & solution_old, BlockVector<double> & system_rhs)` |
| 4129~4193 | `double PhaseFieldMonolithicSolve<dim>::line_search_stepsize_gradient_based(const BlockVector<double> & BFGS_p_vector, const BlockVector<double> & solution_delta)` |
| 4196~4256 | `double PhaseFieldMonolithicSolve<dim>::line_search_stepsize_strong_wolfe(const double phi_0, const double phi_0_prime, const BlockVector<double> & BFGS_p_vector, const BlockVector<double> & solution_delta)` |
| 4371~4426 | `void PhaseFieldMonolithicSolve<dim>::LBFGS_B0(BlockVector<double> & LBFGS_r_vector, BlockVector<double> & LBFGS_q_vector)` |
| 4429~4447 | `void PhaseFieldMonolithicSolve<dim>::print_conv_header_LBFGS()` |
| 4450~4468 | `void PhaseFieldMonolithicSolve<dim>::print_conv_header_LBFGSB()` |
| 4754~4949 | `<< std::fixed << std::setprecision(6) << std::setw(1) << std::scientific << "  " << energy_functional << std::fixed << std::setprecision(3) << std::setw(1) << std::scientific << "  " << m_error_residual_norm.m_norm << "  " << m_error_residual_norm.m_u << "  " << m_error_residual_norm.m_d << "  " << m_error_update_norm.m_norm << "  " << m_error_update_norm.m_u << "  " << m_error_update_norm.m_d << "  " << std::endl; } }  AssertThrow(LBFGS_iteration < m_parameters.m_max_iterations_BFGS, ExcMessage("No convergence in L-BFGS nonlinear solver!")); }  template <int dim> void PhaseFieldMonolithicSolve<dim>:: calculate_cauchy_point(const BlockSparseMatrix<double> & B0_matrix, const std::list<BlockVector<double>> & y_vector_list, const std::list<BlockVector<double>> & b0xs_vector_list, const FullMatrix<double> & M_matrix, const BlockVector<double> & gradient_g, const BlockVector<double> & solution_delta, BlockVector<double> & solution_delta_cauchy_point)` |
| 4757~4949 | `<< std::fixed << std::setprecision(3) << std::setw(1) << std::scientific << "  " << m_error_residual_norm.m_norm << "  " << m_error_residual_norm.m_u << "  " << m_error_residual_norm.m_d << "  " << m_error_update_norm.m_norm << "  " << m_error_update_norm.m_u << "  " << m_error_update_norm.m_d << "  " << std::endl; } }  AssertThrow(LBFGS_iteration < m_parameters.m_max_iterations_BFGS, ExcMessage("No convergence in L-BFGS nonlinear solver!")); }  template <int dim> void PhaseFieldMonolithicSolve<dim>:: calculate_cauchy_point(const BlockSparseMatrix<double> & B0_matrix, const std::list<BlockVector<double>> & y_vector_list, const std::list<BlockVector<double>> & b0xs_vector_list, const FullMatrix<double> & M_matrix, const BlockVector<double> & gradient_g, const BlockVector<double> & solution_delta, BlockVector<double> & solution_delta_cauchy_point)` |
| 5780~5901 | `<< std::setw(6) << number_active_constraint_lower_bound << std::setw(6) << number_active_constraint_upper_bound << std::setw(6) << number_active_constraint_lowerupper_bound; if (m_parameters.m_type_linear_solver == "CG") m_logfile << std::setw(8) << cg_iterations; else m_logfile << std::setw(8) << "---"; m_logfile << "      " << std::fixed << std::setprecision(3) << std::scientific << line_search_parameter << std::fixed << std::setprecision(6) << std::setw(1) << std::scientific << "  " << energy_functional << std::fixed << std::setprecision(3) << std::setw(1) << std::scientific << "  " << m_error_residual_norm.m_norm << "  " << m_error_residual_norm.m_u << "  " << m_error_residual_norm.m_d << "  " << m_error_update_norm.m_norm << "  " << m_error_update_norm.m_u << "  " << m_error_update_norm.m_d << "  " << std::endl; } } // for (; LBFGS_iteration < m_parameters.m_max_iterations_BFGS; ++LBFGS_iteration)  AssertThrow(LBFGS_iteration < m_parameters.m_max_iterations_BFGS, ExcMessage("No convergence in L-BFGS-B nonlinear solver!")); }  template <int dim> void PhaseFieldMonolithicSolve<dim>::output_results() const` |
| 5781~5901 | `<< std::setw(6) << number_active_constraint_upper_bound << std::setw(6) << number_active_constraint_lowerupper_bound; if (m_parameters.m_type_linear_solver == "CG") m_logfile << std::setw(8) << cg_iterations; else m_logfile << std::setw(8) << "---"; m_logfile << "      " << std::fixed << std::setprecision(3) << std::scientific << line_search_parameter << std::fixed << std::setprecision(6) << std::setw(1) << std::scientific << "  " << energy_functional << std::fixed << std::setprecision(3) << std::setw(1) << std::scientific << "  " << m_error_residual_norm.m_norm << "  " << m_error_residual_norm.m_u << "  " << m_error_residual_norm.m_d << "  " << m_error_update_norm.m_norm << "  " << m_error_update_norm.m_u << "  " << m_error_update_norm.m_d << "  " << std::endl; } } // for (; LBFGS_iteration < m_parameters.m_max_iterations_BFGS; ++LBFGS_iteration)  AssertThrow(LBFGS_iteration < m_parameters.m_max_iterations_BFGS, ExcMessage("No convergence in L-BFGS-B nonlinear solver!")); }  template <int dim> void PhaseFieldMonolithicSolve<dim>::output_results() const` |
| 5788~5901 | `<< std::fixed << std::setprecision(3) << std::scientific << line_search_parameter << std::fixed << std::setprecision(6) << std::setw(1) << std::scientific << "  " << energy_functional << std::fixed << std::setprecision(3) << std::setw(1) << std::scientific << "  " << m_error_residual_norm.m_norm << "  " << m_error_residual_norm.m_u << "  " << m_error_residual_norm.m_d << "  " << m_error_update_norm.m_norm << "  " << m_error_update_norm.m_u << "  " << m_error_update_norm.m_d << "  " << std::endl; } } // for (; LBFGS_iteration < m_parameters.m_max_iterations_BFGS; ++LBFGS_iteration)  AssertThrow(LBFGS_iteration < m_parameters.m_max_iterations_BFGS, ExcMessage("No convergence in L-BFGS-B nonlinear solver!")); }  template <int dim> void PhaseFieldMonolithicSolve<dim>::output_results() const` |
| 5789~5901 | `<< std::fixed << std::setprecision(6) << std::setw(1) << std::scientific << "  " << energy_functional << std::fixed << std::setprecision(3) << std::setw(1) << std::scientific << "  " << m_error_residual_norm.m_norm << "  " << m_error_residual_norm.m_u << "  " << m_error_residual_norm.m_d << "  " << m_error_update_norm.m_norm << "  " << m_error_update_norm.m_u << "  " << m_error_update_norm.m_d << "  " << std::endl; } } // for (; LBFGS_iteration < m_parameters.m_max_iterations_BFGS; ++LBFGS_iteration)  AssertThrow(LBFGS_iteration < m_parameters.m_max_iterations_BFGS, ExcMessage("No convergence in L-BFGS-B nonlinear solver!")); }  template <int dim> void PhaseFieldMonolithicSolve<dim>::output_results() const` |
| 5791~5901 | `<< std::fixed << std::setprecision(3) << std::setw(1) << std::scientific << "  " << m_error_residual_norm.m_norm << "  " << m_error_residual_norm.m_u << "  " << m_error_residual_norm.m_d << "  " << m_error_update_norm.m_norm << "  " << m_error_update_norm.m_u << "  " << m_error_update_norm.m_d << "  " << std::endl; } } // for (; LBFGS_iteration < m_parameters.m_max_iterations_BFGS; ++LBFGS_iteration)  AssertThrow(LBFGS_iteration < m_parameters.m_max_iterations_BFGS, ExcMessage("No convergence in L-BFGS-B nonlinear solver!")); }  template <int dim> void PhaseFieldMonolithicSolve<dim>::output_results() const` |
| 5808~5901 | `void PhaseFieldMonolithicSolve<dim>::output_results() const` |
| 5874~5882 | `VectorTools::project(mapping, stresses_dof_handler_L2, constraints, m_qf_cell, [&] (const typename DoFHandler<dim>::active_cell_iterator & cell, const unsigned int q) -> double` |
| 5904~6058 | `void PhaseFieldMonolithicSolve<dim>::calculate_reaction_force(unsigned int face_ID)` |
| 6031~6043 | `DoFTools::map_dof_to_boundary_indices(m_dof_handler, boundary_ids, mapping);  std::vector<double> reaction_force(dim, 0.0);  for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)` |
| 6048~6114 | `<< std::fixed << std::setprecision(3) << std::setw(1) << std::scientific << reaction_force[i] << std::endl;  std::pair<double, std::vector<double>> time_force; time_force.first = m_time.current(); time_force.second = reaction_force; m_history_reaction_force.push_back(time_force);  m_timer.leave_subsection(); }  template <int dim> void PhaseFieldMonolithicSolve<dim>::write_history_data()` |
| 6061~6114 | `void PhaseFieldMonolithicSolve<dim>::write_history_data()` |
| 6096~6109 | `myfile_energy << std::fixed << std::setprecision(10) << std::scientific << 0.0 << "\t" << 0.0 << "\t" << 0.0 << "\t" << 0.0 << std::endl;  for (auto const & time_energy : m_history_energy)` |
| 6104~6140 | `myfile_energy << std::fixed << std::setprecision(10) << std::scientific << time_energy.first     << "\t" << time_energy.second[0] << "\t" << time_energy.second[1] << "\t" << time_energy.second[2] << std::endl; } myfile_energy.close(); } else m_logfile << "Unable to open file"; }  template <int dim> double PhaseFieldMonolithicSolve<dim>::calculate_energy_functional() const` |
| 6117~6140 | `double PhaseFieldMonolithicSolve<dim>::calculate_energy_functional() const` |
| 6144~6168 | `PhaseFieldMonolithicSolve<dim>::calculate_total_strain_energy_and_crack_energy_dissipation() const` |
| 6171~6281 | `bool PhaseFieldMonolithicSolve<dim>::local_refine_and_solution_transfer(BlockVector<double> & solution_delta, BlockVector<double> & LBFGS_update_refine)` |
| 6284~6353 | `void PhaseFieldMonolithicSolve<dim>::print_parameter_information()` |
| 6356~6465 | `void PhaseFieldMonolithicSolve<dim>::run()` |
| 6468~6494 | `int main(int argc, char* argv[])` |