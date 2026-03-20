## 006. ## 1. Introduction（explain.md:L63-L65）

- 对应关系说明：该公式对应能量-残差-切线三层实现。
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
\mathrm{lb}_i\leq x_i\leq \mathrm{ub}_i,
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 007. ## 1. Introduction（explain.md:L69-L71）

- 对应关系说明：该公式对应能量-残差-切线三层实现。
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
\mathrm{lb}_i = d_i^{(n)}\leq d_i^{(n + 1)}\leq 1 = \mathrm{ub}_i,
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 008. ### 2.1. Phase-field formulation（explain.md:L87-L89）

- 对应关系说明：本公式位于本构层：应变谱分解 + 正负能量分裂 + 退化函数。
- 最底层代码链接：
  - [main.cc:L831-L874](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L831-L874)
  - [SpectrumDecomposition.h:L71-L143](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/SpectrumDecomposition.h#L71-L143)
  - [SpectrumDecomposition.cc:L19-L38](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/SpectrumDecomposition.cc#L19-L38)

### 代码片段 1（main.cc:L831-L874）
```cpp
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
```

### 代码片段 2（SpectrumDecomposition.h:L71-L143）
```cpp
  71:   template <int dim>
  72:   void positive_negative_projectors(Vector<double> const & eigenvalues,
  73:                                     std::vector<Tensor<1, dim>> const & eigenvectors,
  74: 			            SymmetricTensor<4, dim> & positive_projector,
  75: 				    SymmetricTensor<4, dim> & negative_projector)
  76:   {
  77:     Assert(dim <= 3,
  78: 	   ExcMessage("Project tensors only work for dim <= 3."));
  79: 
  80:     std::array<SymmetricTensor<2, dim>, dim> M;
  81:     for (int a = 0; a < dim; a++)
  82:       M[a] = symmetrize(outer_product(eigenvectors[a], eigenvectors[a]));
  83: 
  84:     std::array<SymmetricTensor<4, dim>, dim> Q;
  85:     for (int a = 0; a < dim; a++)
  86:       Q[a] = outer_product(M[a], M[a]);
  87: 
  88:     std::array<std::array<SymmetricTensor<4, dim>, dim>, dim> G;
  89:     for (int a = 0; a < dim; a++)
  90:       for (int b = 0; b < dim; b++)
  91: 	for (int i = 0; i < dim; i++)
  92: 	  for (int j = 0; j < dim; j++)
  93: 	    for (int k = 0; k < dim; k++)
  94:               for (int l = 0; l < dim; l++)
  95:         	G[a][b][i][j][k][l] = M[a][i][k] * M[b][j][l]
  96: 				    + M[a][i][l] * M[b][j][k];
  97: 
  98:     positive_projector = 0;
  99:     for (int a = 0; a < dim; a++)
 100:       {
 101: 	double lambda_a = eigenvalues[a];
 102: 	positive_projector += heaviside_function(lambda_a)
 103: 			    * Q[a];
 104: 	for (int b = 0; b < dim; b++)
 105: 	  {
 106: 	    if (b != a)
 107: 	      {
 108: 		double lambda_b = eigenvalues[b];
 109: 		double v_ab = 0.0;
 110: 		if (std::fabs(lambda_a - lambda_b) > 1.0e-12)
 111: 		  v_ab = (positive_ramp_function(lambda_a) - positive_ramp_function(lambda_b))
 112: 		       / (lambda_a - lambda_b);
 113: 		else
 114: 		  v_ab = 0.5 * (  heaviside_function(lambda_a)
 115: 		                + heaviside_function(lambda_b) );
 116: 		positive_projector += 0.5 * v_ab * 0.5 * (G[a][b] + G[b][a]);
 117: 	      }
 118: 	  }
 119:       }
 120: 
 121:     negative_projector = 0;
 122:     for (int a = 0; a < dim; a++)
 123:       {
 124: 	double lambda_a = eigenvalues[a];
 125: 	negative_projector += heaviside_function(-lambda_a)
 126: 			    * Q[a];
 127: 	for (int b = 0; b < dim; b++)
 128: 	  {
 129: 	    if (b != a)
 130: 	      {
 131: 		double lambda_b = eigenvalues[b];
 132: 		double v_ab = 0.0;
 133: 		if (std::fabs(lambda_a - lambda_b) > 1.0e-12)
 134: 		  v_ab = (negative_ramp_function(lambda_a) - negative_ramp_function(lambda_b))
 135: 		       / (lambda_a - lambda_b);
 136: 		else
 137: 		  v_ab = 0.5 * (  heaviside_function(-lambda_a)
 138: 		                + heaviside_function(-lambda_b) );
 139: 		negative_projector += 0.5 * v_ab * 0.5 * (G[a][b] + G[b][a]);
 140: 	      }
 141: 	  }
 142:       }
 143: 
```

### 代码片段 3（SpectrumDecomposition.cc:L19-L38）
```cpp
  19:   double positive_ramp_function(const double x)
  20:   {
  21:     return std::fmax(x, 0.0);
  22:   }
  23: 
  24:   double negative_ramp_function(const double x)
  25:   {
  26:     return std::fmin(x, 0.0);
  27:   }
  28: 
  29:   double heaviside_function(const double x)
  30:   {
  31:     if (std::fabs(x) < 1.0e-16)
  32:       return 0.5;
  33: 
  34:     if (x > 0)
  35:       return 1.0;
  36:     else
  37:       return 0.0;
  38:   }
```

### 公式（与 explain.md 一致）
\[
\psi(\pmb{\epsilon}, d) = [g(d) + k] \psi^{+}(\pmb{\epsilon}) + \psi^{-}(\pmb{\epsilon}), \quad (5)
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 009. ### 2.1. Phase-field formulation（explain.md:L93-L95）

- 对应关系说明：本公式位于本构层：应变谱分解 + 正负能量分裂 + 退化函数。
- 最底层代码链接：
  - [main.cc:L244-L258](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L244-L258)
  - [main.cc:L846-L847](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L846-L847)
  - [main.cc:L3813-L3818](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L3813-L3818)

### 代码片段 1（main.cc:L244-L258）
```cpp
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
```

### 代码片段 2（main.cc:L846-L847）
```cpp
 846:     const double degradation = degradation_function(m_phase_field_value) + m_residual_k;
 847:     const double I_1 = trace(m_strain);
```

### 代码片段 3（main.cc:L3813-L3818）
```cpp
3813:     	        data.m_cell_rhs(i) += (    gc * length_scale * grad_N_phasefield[i] * phasefield_grad
3814:     	                                +  (   gc / length_scale * phasefield_value
3815: 					     + eta / delta_time  * (phasefield_value - old_phasefield)
3816: 					     + degradation_function_derivative(phasefield_value)
3817: 					     * current_positive_strain_energy )
3818: 					  * N_phasefield[i]
```

### 公式（与 explain.md 一致）
\[
g(d) = (1 - d)^{2}. \quad (6)
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 010. ### 2.1. Phase-field formulation（explain.md:L101-L103）

- 对应关系说明：本公式位于本构层：应变谱分解 + 正负能量分裂 + 退化函数。
- 最底层代码链接：
  - [main.cc:L820-L891](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L820-L891)
  - [SpectrumDecomposition.h:L18-L143](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/SpectrumDecomposition.h#L18-L143)
  - [SpectrumDecomposition.cc:L19-L38](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/SpectrumDecomposition.cc#L19-L38)

### 代码片段 1（main.cc:L820-L891）
```cpp
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
```

### 代码片段 2（SpectrumDecomposition.h:L18-L143）
```cpp
  18:   double positive_ramp_function(const double x);
  19: 
  20:   double negative_ramp_function(const double x);
  21: 
  22:   double heaviside_function(const double x);
  23: 
  24:   // templated function has to be defined in the header file
  25:   // perform a spectrum decomposition of a symmetric tensor
  26:   // input: a symmetric tensor (SymmetricTensor<2, matrix_dimension>)
  27:   // output: eigenvalues  (Vector<double>)
  28:   //         eigenvectors (std::vector<Tensor<1, dim>>)
  29:   template <int dim>
  30:   void spectrum_decomposition(SymmetricTensor<2, dim> const & symmetric_tensor,
  31: 			      Vector<double> & myEigenvalues,
  32: 			      std::vector<Tensor<1, dim>> & myEigenvectors)
  33:   {
  34: 
  35:     const std::array< std::pair< double, Tensor< 1, dim > >, dim >
  36:       myEigenSystem = eigenvectors(symmetric_tensor);
  37: 
  38:     for (int i = 0; i < dim; i++)
  39:       {
  40:         myEigenvalues[i] = myEigenSystem[i].first;
  41:         myEigenvectors[i] = myEigenSystem[i].second;
  42:       }
  43:   }
  44: 
  45:   template <int dim>
  46:   SymmetricTensor<2, dim> positive_tensor(Vector<double> const & eigenvalues,
  47: 					  std::vector<Tensor<1, dim>> const & eigenvectors)
  48:   {
  49:     SymmetricTensor<2, dim> positive_part_tensor;
  50:     positive_part_tensor = 0;
  51:     for (int i = 0; i < dim; i++)
  52:       positive_part_tensor += positive_ramp_function(eigenvalues[i])
  53:                             * symmetrize(outer_product(eigenvectors[i],
  54:                                                        eigenvectors[i]));
  55:     return positive_part_tensor;
  56:   }
  57: 
  58:   template <int dim>
  59:   SymmetricTensor<2, dim> negative_tensor(Vector<double> const & eigenvalues,
  60: 					  std::vector<Tensor<1, dim>> const & eigenvectors)
  61:   {
  62:     SymmetricTensor<2, dim> negative_part_tensor;
  63:     negative_part_tensor = 0;
  64:     for (int i = 0; i < dim; i++)
  65:       negative_part_tensor += negative_ramp_function(eigenvalues[i])
  66:                             * symmetrize(outer_product(eigenvectors[i],
  67:                                                        eigenvectors[i]));
  68:     return negative_part_tensor;
  69:   }
  70: 
  71:   template <int dim>
  72:   void positive_negative_projectors(Vector<double> const & eigenvalues,
  73:                                     std::vector<Tensor<1, dim>> const & eigenvectors,
  74: 			            SymmetricTensor<4, dim> & positive_projector,
  75: 				    SymmetricTensor<4, dim> & negative_projector)
  76:   {
  77:     Assert(dim <= 3,
  78: 	   ExcMessage("Project tensors only work for dim <= 3."));
  79: 
  80:     std::array<SymmetricTensor<2, dim>, dim> M;
  81:     for (int a = 0; a < dim; a++)
  82:       M[a] = symmetrize(outer_product(eigenvectors[a], eigenvectors[a]));
  83: 
  84:     std::array<SymmetricTensor<4, dim>, dim> Q;
  85:     for (int a = 0; a < dim; a++)
  86:       Q[a] = outer_product(M[a], M[a]);
  87: 
  88:     std::array<std::array<SymmetricTensor<4, dim>, dim>, dim> G;
  89:     for (int a = 0; a < dim; a++)
  90:       for (int b = 0; b < dim; b++)
  91: 	for (int i = 0; i < dim; i++)
  92: 	  for (int j = 0; j < dim; j++)
  93: 	    for (int k = 0; k < dim; k++)
  94:               for (int l = 0; l < dim; l++)
  95:         	G[a][b][i][j][k][l] = M[a][i][k] * M[b][j][l]
  96: 				    + M[a][i][l] * M[b][j][k];
  97: 
  98:     positive_projector = 0;
  99:     for (int a = 0; a < dim; a++)
 100:       {
 101: 	double lambda_a = eigenvalues[a];
 102: 	positive_projector += heaviside_function(lambda_a)
 103: 			    * Q[a];
 104: 	for (int b = 0; b < dim; b++)
 105: 	  {
 106: 	    if (b != a)
 107: 	      {
 108: 		double lambda_b = eigenvalues[b];
 109: 		double v_ab = 0.0;
 110: 		if (std::fabs(lambda_a - lambda_b) > 1.0e-12)
 111: 		  v_ab = (positive_ramp_function(lambda_a) - positive_ramp_function(lambda_b))
 112: 		       / (lambda_a - lambda_b);
 113: 		else
 114: 		  v_ab = 0.5 * (  heaviside_function(lambda_a)
 115: 		                + heaviside_function(lambda_b) );
 116: 		positive_projector += 0.5 * v_ab * 0.5 * (G[a][b] + G[b][a]);
 117: 	      }
 118: 	  }
 119:       }
 120: 
 121:     negative_projector = 0;
 122:     for (int a = 0; a < dim; a++)
 123:       {
 124: 	double lambda_a = eigenvalues[a];
 125: 	negative_projector += heaviside_function(-lambda_a)
 126: 			    * Q[a];
 127: 	for (int b = 0; b < dim; b++)
 128: 	  {
 129: 	    if (b != a)
 130: 	      {
 131: 		double lambda_b = eigenvalues[b];
 132: 		double v_ab = 0.0;
 133: 		if (std::fabs(lambda_a - lambda_b) > 1.0e-12)
 134: 		  v_ab = (negative_ramp_function(lambda_a) - negative_ramp_function(lambda_b))
 135: 		       / (lambda_a - lambda_b);
 136: 		else
 137: 		  v_ab = 0.5 * (  heaviside_function(-lambda_a)
 138: 		                + heaviside_function(-lambda_b) );
 139: 		negative_projector += 0.5 * v_ab * 0.5 * (G[a][b] + G[b][a]);
 140: 	      }
 141: 	  }
 142:       }
 143: 
```

### 代码片段 3（SpectrumDecomposition.cc:L19-L38）
```cpp
  19:   double positive_ramp_function(const double x)
  20:   {
  21:     return std::fmax(x, 0.0);
  22:   }
  23: 
  24:   double negative_ramp_function(const double x)
  25:   {
  26:     return std::fmin(x, 0.0);
  27:   }
  28: 
  29:   double heaviside_function(const double x)
  30:   {
  31:     if (std::fabs(x) < 1.0e-16)
  32:       return 0.5;
  33: 
  34:     if (x > 0)
  35:       return 1.0;
  36:     else
  37:       return 0.0;
  38:   }
```

### 公式（与 explain.md 一致）
\[
\langle x\rangle_{+} = \frac{1}{2}(x + |x|), \quad \langle x\rangle_{-} = \frac{1}{2}(x - |x|), \quad H(x) = \begin{cases} 1 & x \geq 0, \\ 0 & x < 0. \end{cases}
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 011. ### 2.1. Phase-field formulation（explain.md:L107-L109）

- 对应关系说明：本公式位于本构层：应变谱分解 + 正负能量分裂 + 退化函数。
- 最底层代码链接：
  - [main.cc:L831-L874](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L831-L874)
  - [SpectrumDecomposition.h:L71-L143](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/SpectrumDecomposition.h#L71-L143)
  - [SpectrumDecomposition.cc:L19-L38](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/SpectrumDecomposition.cc#L19-L38)

### 代码片段 1（main.cc:L831-L874）
```cpp
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
```

### 代码片段 2（SpectrumDecomposition.h:L71-L143）
```cpp
  71:   template <int dim>
  72:   void positive_negative_projectors(Vector<double> const & eigenvalues,
  73:                                     std::vector<Tensor<1, dim>> const & eigenvectors,
  74: 			            SymmetricTensor<4, dim> & positive_projector,
  75: 				    SymmetricTensor<4, dim> & negative_projector)
  76:   {
  77:     Assert(dim <= 3,
  78: 	   ExcMessage("Project tensors only work for dim <= 3."));
  79: 
  80:     std::array<SymmetricTensor<2, dim>, dim> M;
  81:     for (int a = 0; a < dim; a++)
  82:       M[a] = symmetrize(outer_product(eigenvectors[a], eigenvectors[a]));
  83: 
  84:     std::array<SymmetricTensor<4, dim>, dim> Q;
  85:     for (int a = 0; a < dim; a++)
  86:       Q[a] = outer_product(M[a], M[a]);
  87: 
  88:     std::array<std::array<SymmetricTensor<4, dim>, dim>, dim> G;
  89:     for (int a = 0; a < dim; a++)
  90:       for (int b = 0; b < dim; b++)
  91: 	for (int i = 0; i < dim; i++)
  92: 	  for (int j = 0; j < dim; j++)
  93: 	    for (int k = 0; k < dim; k++)
  94:               for (int l = 0; l < dim; l++)
  95:         	G[a][b][i][j][k][l] = M[a][i][k] * M[b][j][l]
  96: 				    + M[a][i][l] * M[b][j][k];
  97: 
  98:     positive_projector = 0;
  99:     for (int a = 0; a < dim; a++)
 100:       {
 101: 	double lambda_a = eigenvalues[a];
 102: 	positive_projector += heaviside_function(lambda_a)
 103: 			    * Q[a];
 104: 	for (int b = 0; b < dim; b++)
 105: 	  {
 106: 	    if (b != a)
 107: 	      {
 108: 		double lambda_b = eigenvalues[b];
 109: 		double v_ab = 0.0;
 110: 		if (std::fabs(lambda_a - lambda_b) > 1.0e-12)
 111: 		  v_ab = (positive_ramp_function(lambda_a) - positive_ramp_function(lambda_b))
 112: 		       / (lambda_a - lambda_b);
 113: 		else
 114: 		  v_ab = 0.5 * (  heaviside_function(lambda_a)
 115: 		                + heaviside_function(lambda_b) );
 116: 		positive_projector += 0.5 * v_ab * 0.5 * (G[a][b] + G[b][a]);
 117: 	      }
 118: 	  }
 119:       }
 120: 
 121:     negative_projector = 0;
 122:     for (int a = 0; a < dim; a++)
 123:       {
 124: 	double lambda_a = eigenvalues[a];
 125: 	negative_projector += heaviside_function(-lambda_a)
 126: 			    * Q[a];
 127: 	for (int b = 0; b < dim; b++)
 128: 	  {
 129: 	    if (b != a)
 130: 	      {
 131: 		double lambda_b = eigenvalues[b];
 132: 		double v_ab = 0.0;
 133: 		if (std::fabs(lambda_a - lambda_b) > 1.0e-12)
 134: 		  v_ab = (negative_ramp_function(lambda_a) - negative_ramp_function(lambda_b))
 135: 		       / (lambda_a - lambda_b);
 136: 		else
 137: 		  v_ab = 0.5 * (  heaviside_function(-lambda_a)
 138: 		                + heaviside_function(-lambda_b) );
 139: 		negative_projector += 0.5 * v_ab * 0.5 * (G[a][b] + G[b][a]);
 140: 	      }
 141: 	  }
 142:       }
 143: 
```

### 代码片段 3（SpectrumDecomposition.cc:L19-L38）
```cpp
  19:   double positive_ramp_function(const double x)
  20:   {
  21:     return std::fmax(x, 0.0);
  22:   }
  23: 
  24:   double negative_ramp_function(const double x)
  25:   {
  26:     return std::fmin(x, 0.0);
  27:   }
  28: 
  29:   double heaviside_function(const double x)
  30:   {
  31:     if (std::fabs(x) < 1.0e-16)
  32:       return 0.5;
  33: 
  34:     if (x > 0)
  35:       return 1.0;
  36:     else
  37:       return 0.0;
  38:   }
```

### 公式（与 explain.md 一致）
\[
\pmb{\epsilon} = \sum_{\alpha} \epsilon_{\alpha} \mathbf{M}_{\alpha}, \quad \mathbf{M}_{\alpha} = \pmb{n}_{\alpha} \otimes \pmb{n}_{\alpha},
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。
