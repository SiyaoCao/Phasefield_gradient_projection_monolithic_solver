## 3. Gradient projection based monolithic scheme > 3.4. Subspace minimization

### Eq. (26) — 子空间 box 约束
- 论文位置：[3. Gradient projection based monolithic scheme > 3.4. Subspace minimization，`explain.md` 第 579 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L579)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\mathrm{lb}_i \leq x_i \leq \mathrm{ub}_i, \quad \forall i \notin \mathcal{A}(\mathbf{x}^c). \quad (26)
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

## 3. Gradient projection based monolithic scheme > 3.4. Subspace minimization > 3.4.1. Direct matrix factorization for the primal approach

### 未编号公式 #82 — 未编号公式：\mathbf{x} = \mathbf{x}^c ...
- 论文位置：[3. Gradient projection based monolithic scheme > 3.4. Subspace minimization > 3.4.1. Direct matrix factorization for the primal approach，`explain.md` 第 589 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L589)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\mathbf{x} = \mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}},
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - 自由变量筛选与右端构造（Z^T[...]）：[main.cc:5193-5272](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5193-L5272)
```cpp
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
```
  - 子空间 CG 求解（原始法）：[main.cc:5280-5479](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5280-L5479)
```cpp
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
```
  - 子空间直接法 + SMW 等价实现：[main.cc:5481-5621](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5481-L5621)
```cpp
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
```
### 未编号公式 #83 — 未编号公式：\begin{array}{rl}  m_k(\ma...
- 论文位置：[3. Gradient projection based monolithic scheme > 3.4. Subspace minimization > 3.4.1. Direct matrix factorization for the primal approach，`explain.md` 第 595 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L595)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\begin{array}{rl} 
m_k(\mathbf{x}) &= \Pi_k + \mathbf{r}_k^{\mathrm{T}}(\mathbf{x} - \mathbf{x}_k) + \frac{1}{2}(\mathbf{x} - \mathbf{x}_k)^{\mathrm{T}}\mathbf{B}_k(\mathbf{x} - \mathbf{x}_k) \\
\implies \hat{m}_k(\hat{\mathbf{x}}) &= \Pi_k + \mathbf{r}_k^{\mathrm{T}}(\mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}} - \mathbf{x}_k) + \frac{1}{2}(\mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}} - \mathbf{x}_k)^{\mathrm{T}}\mathbf{B}_k(\mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}} - \mathbf{x}_k) \\
&= \Pi_k + \mathbf{r}_k^{\mathrm{T}}(\mathbf{x}^c - \mathbf{x}_k) + \frac{1}{2} (\mathbf{x}^c - \mathbf{x}_k)^{\mathrm{T}}\mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k) \\
&\quad + [\mathbf{r}_k + \mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k)]^{\mathrm{T}}\mathbf{Z}_k \hat{\mathbf{x}} + \frac{1}{2} \hat{\mathbf{x}}^{\mathrm{T}}\mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k \hat{\mathbf{x}} \\
&= \Pi_k + \mathbf{r}_k^{\mathrm{T}}(\mathbf{x}^c - \mathbf{x}_k) + \frac{1}{2} (\mathbf{x}^c - \mathbf{x}_k)^{\mathrm{T}}\mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k) \\
&\quad + [\mathbf{Z}^{\mathrm{T}}_k [\mathbf{r}_k + \mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k)]]^{\mathrm{T}} \hat{\mathbf{x}} + \frac{1}{2} \hat{\mathbf{x}}^{\mathrm{T}}\mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k \hat{\mathbf{x}}.
\end{array}
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - 裂纹耗散能密度（含梯度项）：[main.cc:886-891](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L886-L891)
```cpp
886:     m_crack_energy_dissipation = m_gc * (  0.5 / m_length_scale * m_phase_field_value * m_phase_field_value
887: 	                                   + 0.5 * m_length_scale * m_grad_phasefield * m_grad_phasefield)
888: 	                                   // the term due to viscosity regularization
889: 	                                   + (m_phase_field_value - phase_field_value_previous_step)
890: 					   * (m_phase_field_value - phase_field_value_previous_step)
891: 				           * 0.5 * m_eta / delta_time;
```
  - 总能量/裂纹能积分计算：[main.cc:6117-6167](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L6117-L6167)
```cpp
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
```
  - 紧凑 L-BFGS 矩阵构造（W/M）：[main.cc:5113-5183](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5113-L5183)
```cpp
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
```
