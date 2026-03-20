## 3. Gradient projection based monolithic scheme > 3.4. Subspace minimization > 3.4.3. Schur complement for the dual approach

### Eq. (31) — Schur 补方程
- 论文位置：[3. Gradient projection based monolithic scheme > 3.4. Subspace minimization > 3.4.3. Schur complement for the dual approach，`explain.md` 第 705 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L705)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
(\mathbf{Q}_k^{\mathrm{T}} \mathbf{B}_k^{-1} \mathbf{Q}_k) \lambda = -\mathbf{Q}_k^{\mathrm{T}} \mathbf{B}_k^{-1} \mathbf{r}_k - \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k). \quad (31)
\]
- 对应最底层计算代码（公式在前，代码在后）：
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
  - (s,y) 向量对更新与 limited-memory 维护：[main.cc:5713-5729](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5713-L5729)
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
5729: 	  }
```
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
### Eq. (32) — 从拉格朗日乘子回代
- 论文位置：[3. Gradient projection based monolithic scheme > 3.4. Subspace minimization > 3.4.3. Schur complement for the dual approach，`explain.md` 第 711 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L711)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\mathbf{B}_k \Delta \mathbf{x}_k = -(\mathbf{r}_k + \mathbf{Q}_k \lambda). \quad (32)
\]
- 对应最底层计算代码（公式在前，代码在后）：
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
  - (s,y) 向量对更新与 limited-memory 维护：[main.cc:5713-5729](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5713-L5729)
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
5729: 	  }
```
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
### 未编号公式 #103 — 未编号公式：\mathbf{p}_k = \mathbf{x}^...
- 论文位置：[3. Gradient projection based monolithic scheme > 3.4. Subspace minimization > 3.4.3. Schur complement for the dual approach，`explain.md` 第 719 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L719)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\mathbf{p}_k = \mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}} - \mathbf{x}_k.
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
### 未编号公式 #104 — 未编号公式：\mathbf{x}_{k + 1} = \math...
- 论文位置：[3. Gradient projection based monolithic scheme > 3.4. Subspace minimization > 3.4.3. Schur complement for the dual approach，`explain.md` 第 725 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L725)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\mathbf{x}_{k + 1} = \mathbf{x}_k + \alpha_k \mathbf{p}_k,
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - 强 Wolfe 线搜索主流程：[main.cc:4196-4256](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4196-L4256)
```cpp
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
```
  - zoom 子过程：[main.cc:4259-4311](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4259-L4311)
```cpp
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
```
  - phi/phi' 计算（能量+方向导）：[main.cc:4348-4367](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4348-L4367)
```cpp
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
```
