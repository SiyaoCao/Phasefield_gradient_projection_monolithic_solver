## 049. ### 3.1. Algorithm overview（explain.md:L375-L377）

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
(\mathbf{I}_{d_c})_{ii} = \begin{cases} 1 & \text{if DoF } i \text{ is constrained}, \\ 0 & \text{otherwise}, \end{cases}
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 050. ### 3.2. Compact representation of limited-memory BFGS matrix（explain.md:L384-L386）

- 对应关系说明：本公式位于L-BFGS紧凑存储层：s/y历史向量与M矩阵构造。
- 最底层代码链接：
  - [main.cc:L5119-L5182](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L5119-L5182)
  - [main.cc:L5713-L5728](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L5713-L5728)

### 代码片段 1（main.cc:L5119-L5182）
```cpp
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
```

### 代码片段 2（main.cc:L5713-L5728）
```cpp
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
```

### 公式（与 explain.md 一致）
\[
\mathbf{s}_k = \mathbf{x}_{k + 1} - \mathbf{x}_k,\quad \mathbf{y}_k = \mathbf{r}_{k + 1} - \mathbf{r}_k,
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 051. ### 3.2. Compact representation of limited-memory BFGS matrix（explain.md:L390-L392）

- 对应关系说明：本公式位于L-BFGS紧凑存储层：s/y历史向量与M矩阵构造。
- 最底层代码链接：
  - [main.cc:L5119-L5182](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L5119-L5182)
  - [main.cc:L5713-L5728](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L5713-L5728)

### 代码片段 1（main.cc:L5119-L5182）
```cpp
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
```

### 代码片段 2（main.cc:L5713-L5728）
```cpp
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
```

### 公式（与 explain.md 一致）
\[
\mathbf{B}_{k + 1} = \mathbf{B}_k - \frac{\mathbf{B}_k\mathbf{s}_k\mathbf{s}_k^{\mathrm{T}}\mathbf{B}_k}{\mathbf{s}_k^{\mathrm{T}}\mathbf{B}_k\mathbf{s}_k} + \frac{\mathbf{y}_k\mathbf{y}_k^{\mathrm{T}}}{\mathbf{y}_k^{\mathrm{T}}\mathbf{s}_k}.
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 052. ### 3.2. Compact representation of limited-memory BFGS matrix（explain.md:L396-L398）

- 对应关系说明：本公式位于L-BFGS紧凑存储层：s/y历史向量与M矩阵构造。
- 最底层代码链接：
  - [main.cc:L5119-L5182](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L5119-L5182)
  - [main.cc:L5713-L5728](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L5713-L5728)

### 代码片段 1（main.cc:L5119-L5182）
```cpp
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
```

### 代码片段 2（main.cc:L5713-L5728）
```cpp
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
```

### 公式（与 explain.md 一致）
\[
\mathbf{s}_k^{\mathrm{T}}\mathbf{y}_k > 0,
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 053. ### 3.2. Compact representation of limited-memory BFGS matrix（explain.md:L404-L406）

- 对应关系说明：本公式位于L-BFGS紧凑存储层：s/y历史向量与M矩阵构造。
- 最底层代码链接：
  - [main.cc:L5119-L5182](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L5119-L5182)

### 代码片段 1（main.cc:L5119-L5182）
```cpp
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
```

### 公式（与 explain.md 一致）
\[
\mathbf{S}_k = [\mathbf{s}_{k-m} \cdots \mathbf{s}_{k-1}]
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 054. ### 3.2. Compact representation of limited-memory BFGS matrix（explain.md:L408-L410）

- 对应关系说明：本公式位于L-BFGS紧凑存储层：s/y历史向量与M矩阵构造。
- 最底层代码链接：
  - [main.cc:L5119-L5182](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L5119-L5182)

### 代码片段 1（main.cc:L5119-L5182）
```cpp
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
```

### 公式（与 explain.md 一致）
\[
\mathbf{Y}_k = [\mathbf{y}_{k-m} \cdots \mathbf{y}_{k-1}].
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 055. ### 3.2. Compact representation of limited-memory BFGS matrix（explain.md:L414-L416）

- 对应关系说明：本公式位于L-BFGS紧凑存储层：s/y历史向量与M矩阵构造。
- 最底层代码链接：
  - [main.cc:L5119-L5182](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L5119-L5182)

### 代码片段 1（main.cc:L5119-L5182）
```cpp
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
```

### 公式（与 explain.md 一致）
\[
\mathbf{B}_k = \mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k, \quad (20)
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 056. ### 3.2. Compact representation of limited-memory BFGS matrix（explain.md:L420-L422）

- 对应关系说明：本公式位于L-BFGS紧凑存储层：s/y历史向量与M矩阵构造。
- 最底层代码链接：
  - [main.cc:L5119-L5182](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L5119-L5182)

### 代码片段 1（main.cc:L5119-L5182）
```cpp
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
```

### 公式（与 explain.md 一致）
\[
\mathbf{W}_k = [\mathbf{Y}_k \quad \mathbf{B}^0_k \mathbf{S}_k] \in \mathbb{R}^{n\times 2m}
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 057. ### 3.2. Compact representation of limited-memory BFGS matrix（explain.md:L425-L427）

- 对应关系说明：本公式位于L-BFGS紧凑存储层：s/y历史向量与M矩阵构造。
- 最底层代码链接：
  - [main.cc:L5119-L5182](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L5119-L5182)

### 代码片段 1（main.cc:L5119-L5182）
```cpp
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
```

### 公式（与 explain.md 一致）
\[
\mathbf{M}_k = \begin{bmatrix} -\mathbf{D}_k & \mathbf{L}^{\mathrm{T}}_k \\ \mathbf{L}_k & \mathbf{S}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{S}_k \end{bmatrix}^{-1} \in \mathbb{R}^{2m\times 2m}.
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。
