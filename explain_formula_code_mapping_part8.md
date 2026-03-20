## 031. ### 3.1. Algorithm overview（explain.md:L259-L261）

- 对应关系说明：本公式位于算法主循环层：投影、搜索方向、线搜索与收敛判据。
- 最底层代码链接：
  - [main.cc:L1379-L1442](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L1379-L1442)
  - [main.cc:L4953-L5060](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L4953-L5060)
  - [main.cc:L5658-L5691](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L5658-L5691)

### 代码片段 1（main.cc:L1379-L1442）
```cpp
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
```

### 代码片段 2（main.cc:L4953-L5060）
```cpp
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
```

### 代码片段 3（main.cc:L5658-L5691）
```cpp
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
```

### 公式（与 explain.md 一致）
\[
m_{k}(\pmb{x}) = \Pi_{k} + \pmb{r}_{k}^{\mathrm{T}}(\pmb{x} - \pmb{x}_{k}) + \frac{1}{2} (\pmb{x} - \pmb{x}_{k})^{\mathrm{T}}\mathbf{B}_{k}(\pmb{x} - \pmb{x}_{k}) \quad (15)
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 032. ### 3.1. Algorithm overview（explain.md:L267-L269）

- 对应关系说明：本公式位于算法主循环层：投影、搜索方向、线搜索与收敛判据。
- 最底层代码链接：
  - [main.cc:L1379-L1394](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L1379-L1394)
  - [main.cc:L1401-L1442](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L1401-L1442)
  - [main.cc:L4775-L4807](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L4775-L4807)

### 代码片段 1（main.cc:L1379-L1394）
```cpp
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
```

### 代码片段 2（main.cc:L1401-L1442）
```cpp
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
```

### 代码片段 3（main.cc:L4775-L4807）
```cpp
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
```

### 公式（与 explain.md 一致）
\[
\mathrm{Proj}_c(x_i,\mathrm{lb}_i,\mathrm{ub}_i) = \left\{ \begin{array}{ll}\mathrm{lb}_i & \mathrm{if}\; x_i < \mathrm{lb}_i,\\ x_i & \mathrm{if}\; \mathrm{lb}_i \leq x_i \leq \mathrm{ub}_i,\\ \mathrm{ub}_i & \mathrm{if}\; x_i > \mathrm{ub}_i, \end{array} \right. \quad (16)
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 033. ### 3.1. Algorithm overview（explain.md:L273-L275）

- 对应关系说明：本公式位于算法主循环层：投影、搜索方向、线搜索与收敛判据。
- 最底层代码链接：
  - [main.cc:L1379-L1394](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L1379-L1394)
  - [main.cc:L1401-L1442](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L1401-L1442)
  - [main.cc:L4775-L4807](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L4775-L4807)

### 代码片段 1（main.cc:L1379-L1394）
```cpp
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
```

### 代码片段 2（main.cc:L1401-L1442）
```cpp
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
```

### 代码片段 3（main.cc:L4775-L4807）
```cpp
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
```

### 公式（与 explain.md 一致）
\[
\pmb {x}(t) = \mathrm{Proj}_c(\pmb {x}_k - t\pmb {r}_k,\mathbf{lb},\mathbf{ub}),\quad t > 0, \quad (17)
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 034. ### 3.1. Algorithm overview（explain.md:L281-L283）

- 对应关系说明：本公式位于算法主循环层：投影、搜索方向、线搜索与收敛判据。
- 最底层代码链接：
  - [main.cc:L1379-L1442](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L1379-L1442)
  - [main.cc:L4953-L5060](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L4953-L5060)
  - [main.cc:L5658-L5691](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L5658-L5691)

### 代码片段 1（main.cc:L1379-L1442）
```cpp
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
```

### 代码片段 2（main.cc:L4953-L5060）
```cpp
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
```

### 代码片段 3（main.cc:L5658-L5691）
```cpp
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
```

### 公式（与 explain.md 一致）
\[
p_k(t) = m_k(\pmb{x}(t)) = \Pi_k + \pmb{r}_k^\mathrm{T}(\pmb{x}(t) - \pmb{x}_k) + \frac{1}{2} (\pmb{x}(t) - \pmb{x}_k)^\mathrm{T}\mathbf{B}_k(\pmb{x}(t) - \pmb{x}_k). \quad (18)
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 035. ### 3.1. Algorithm overview（explain.md:L289-L291）

- 对应关系说明：本公式位于算法主循环层：投影、搜索方向、线搜索与收敛判据。
- 最底层代码链接：
  - [main.cc:L1379-L1394](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L1379-L1394)
  - [main.cc:L1401-L1442](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L1401-L1442)
  - [main.cc:L4775-L4807](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L4775-L4807)

### 代码片段 1（main.cc:L1379-L1394）
```cpp
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
```

### 代码片段 2（main.cc:L1401-L1442）
```cpp
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
```

### 代码片段 3（main.cc:L4775-L4807）
```cpp
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
```

### 公式（与 explain.md 一致）
\[
\mathcal{A}(\pmb {x}) = \{i:x_i = \mathrm{lb}_i\} \cup \{i:x_i = \mathrm{ub}_i\}.
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。
