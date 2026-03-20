## 1. Introduction

### 未编号公式 #7 — 未编号公式：\mathrm{lb}_i = d_i^{(n)}\...
- 论文位置：[1. Introduction，`explain.md` 第 69 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L69)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\mathrm{lb}_i = d_i^{(n)}\leq d_i^{(n + 1)}\leq 1 = \mathrm{ub}_i,
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - 相场点投影（box 约束）：[main.cc:1379-1394](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L1379-L1394)
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
  - break points 计算与初始活跃集：[main.cc:1401-1443](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L1401-L1443)
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
1443:   }
```
  - 广义 Cauchy 点与活跃集更新：[main.cc:4774-4948](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4774-L4948)
```cpp
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
```

## 2. Phase-field formulation and finite element discretization > 2.1. Phase-field formulation

### Eq. (5) — 张拉-受压分裂能量密度
- 论文位置：[2. Phase-field formulation and finite element discretization > 2.1. Phase-field formulation，`explain.md` 第 87 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L87)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\psi(\pmb{\epsilon}, d) = [g(d) + k] \psi^{+}(\pmb{\epsilon}) + \psi^{-}(\pmb{\epsilon}), \quad (5)
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - 退化函数及其导数（最底层）：[main.cc:244-257](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L244-L257)
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
```
  - 相场残差中退化函数导数项：[main.cc:3813-3819](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3813-L3819)
```cpp
3813:     	        data.m_cell_rhs(i) += (    gc * length_scale * grad_N_phasefield[i] * phasefield_grad
3814:     	                                +  (   gc / length_scale * phasefield_value
3815: 					     + eta / delta_time  * (phasefield_value - old_phasefield)
3816: 					     + degradation_function_derivative(phasefield_value)
3817: 					     * current_positive_strain_energy )
3818: 					  * N_phasefield[i]
3819: 				      ) * JxW;
```
  - 相场切线中退化函数二阶导项：[main.cc:3946-3951](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3946-L3951)
```cpp
3946:                     data.m_cell_matrix(i, j) += (  (   gc/length_scale + eta/delta_time
3947:                 	                             + degradation_function_2nd_order_derivative(phasefield_value)
3948: 						     * current_positive_strain_energy  )
3949:                 	                          * N_phasefield[i] * N_phasefield[j]
3950: 					          + gc * length_scale * grad_N_phasefield[i] * grad_N_phasefield[j]
3951: 					        ) * JxW;
```
### Eq. (6) — 退化函数
- 论文位置：[2. Phase-field formulation and finite element discretization > 2.1. Phase-field formulation，`explain.md` 第 93 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L93)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
g(d) = (1 - d)^{2}. \quad (6)
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - 退化函数及其导数（最底层）：[main.cc:244-257](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L244-L257)
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
```
  - 相场残差中退化函数导数项：[main.cc:3813-3819](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3813-L3819)
```cpp
3813:     	        data.m_cell_rhs(i) += (    gc * length_scale * grad_N_phasefield[i] * phasefield_grad
3814:     	                                +  (   gc / length_scale * phasefield_value
3815: 					     + eta / delta_time  * (phasefield_value - old_phasefield)
3816: 					     + degradation_function_derivative(phasefield_value)
3817: 					     * current_positive_strain_energy )
3818: 					  * N_phasefield[i]
3819: 				      ) * JxW;
```
  - 相场切线中退化函数二阶导项：[main.cc:3946-3951](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3946-L3951)
```cpp
3946:                     data.m_cell_matrix(i, j) += (  (   gc/length_scale + eta/delta_time
3947:                 	                             + degradation_function_2nd_order_derivative(phasefield_value)
3948: 						     * current_positive_strain_energy  )
3949:                 	                          * N_phasefield[i] * N_phasefield[j]
3950: 					          + gc * length_scale * grad_N_phasefield[i] * grad_N_phasefield[j]
3951: 					        ) * JxW;
```
### 未编号公式 #10 — 未编号公式：\langle x\rangle_{+} = \fr...
- 论文位置：[2. Phase-field formulation and finite element discretization > 2.1. Phase-field formulation，`explain.md` 第 101 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L101)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\langle x\rangle_{+} = \frac{1}{2}(x + |x|), \quad \langle x\rangle_{-} = \frac{1}{2}(x - |x|), \quad H(x) = \begin{cases} 1 & x \geq 0, \\ 0 & x < 0. \end{cases}
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - Ramp/Heaviside 基础函数：[SpectrumDecomposition.cc:19-38](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.cc#L19-L38)
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
  - 应变正负分解：[SpectrumDecomposition.h:45-69](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L45-L69)
```cpp
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
```
  - 四阶投影张量 P+/P-：[SpectrumDecomposition.h:71-143](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L71-L143)
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
### 未编号公式 #11 — 未编号公式：\pmb{\epsilon} = \sum_{\al...
- 论文位置：[2. Phase-field formulation and finite element discretization > 2.1. Phase-field formulation，`explain.md` 第 107 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L107)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\pmb{\epsilon} = \sum_{\alpha} \epsilon_{\alpha} \mathbf{M}_{\alpha}, \quad \mathbf{M}_{\alpha} = \pmb{n}_{\alpha} \otimes \pmb{n}_{\alpha},
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - Ramp/Heaviside 基础函数：[SpectrumDecomposition.cc:19-38](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.cc#L19-L38)
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
  - 应变正负分解：[SpectrumDecomposition.h:45-69](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L45-L69)
```cpp
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
```
  - 四阶投影张量 P+/P-：[SpectrumDecomposition.h:71-143](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L71-L143)
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
### 未编号公式 #12 — 未编号公式：\pmb{\epsilon}^{+} = \sum_...
- 论文位置：[2. Phase-field formulation and finite element discretization > 2.1. Phase-field formulation，`explain.md` 第 113 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L113)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\pmb{\epsilon}^{+} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{+} \mathbf{M}_{\alpha}, \quad \pmb{\epsilon}^{-} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{-} \mathbf{M}_{\alpha}.
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - Ramp/Heaviside 基础函数：[SpectrumDecomposition.cc:19-38](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.cc#L19-L38)
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
  - 应变正负分解：[SpectrumDecomposition.h:45-69](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L45-L69)
```cpp
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
```
  - 四阶投影张量 P+/P-：[SpectrumDecomposition.h:71-143](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L71-L143)
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
### 未编号公式 #13 — 未编号公式：\psi^{+}(\pmb{\epsilon}) =...
- 论文位置：[2. Phase-field formulation and finite element discretization > 2.1. Phase-field formulation，`explain.md` 第 119 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L119)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\psi^{+}(\pmb{\epsilon}) = \frac{1}{2}\lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{+}^{2} + \mu \pmb{\epsilon}^{+} : \pmb{\epsilon}^{+}, \quad \psi^{-}(\pmb{\epsilon}) = \frac{1}{2} \lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{-}^{2} + \mu \pmb{\epsilon}^{-} : \pmb{\epsilon}^{-},
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - Ramp/Heaviside 基础函数：[SpectrumDecomposition.cc:19-38](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.cc#L19-L38)
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
  - 应变正负分解：[SpectrumDecomposition.h:45-69](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L45-L69)
```cpp
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
```
  - 四阶投影张量 P+/P-：[SpectrumDecomposition.h:71-143](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L71-L143)
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
### 未编号公式 #14 — 未编号公式：\pmb{\sigma} = \frac{\part...
- 论文位置：[2. Phase-field formulation and finite element discretization > 2.1. Phase-field formulation，`explain.md` 第 125 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L125)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\pmb{\sigma} = \frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial \pmb{\epsilon}} = [g(d) + k]\frac{\partial\psi^{+}(\pmb{\epsilon})}{\partial \pmb{\epsilon}} + \frac{\partial\psi^{-}(\pmb{\epsilon})}{\partial \pmb{\epsilon}} = [g(d) + k]\pmb{\sigma}^{+} + \pmb{\sigma}^{-},
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - 退化函数及其导数（最底层）：[main.cc:244-257](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L244-L257)
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
```
  - 相场残差中退化函数导数项：[main.cc:3813-3819](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3813-L3819)
```cpp
3813:     	        data.m_cell_rhs(i) += (    gc * length_scale * grad_N_phasefield[i] * phasefield_grad
3814:     	                                +  (   gc / length_scale * phasefield_value
3815: 					     + eta / delta_time  * (phasefield_value - old_phasefield)
3816: 					     + degradation_function_derivative(phasefield_value)
3817: 					     * current_positive_strain_energy )
3818: 					  * N_phasefield[i]
3819: 				      ) * JxW;
```
  - 相场切线中退化函数二阶导项：[main.cc:3946-3951](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3946-L3951)
```cpp
3946:                     data.m_cell_matrix(i, j) += (  (   gc/length_scale + eta/delta_time
3947:                 	                             + degradation_function_2nd_order_derivative(phasefield_value)
3948: 						     * current_positive_strain_energy  )
3949:                 	                          * N_phasefield[i] * N_phasefield[j]
3950: 					          + gc * length_scale * grad_N_phasefield[i] * grad_N_phasefield[j]
3951: 					        ) * JxW;
```
