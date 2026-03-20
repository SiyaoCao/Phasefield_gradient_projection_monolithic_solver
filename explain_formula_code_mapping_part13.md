## 058. ### 3.2. Compact representation of limited-memory BFGS matrix（explain.md:L431-L433）

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
\mathbf{D}_k = \mathrm{diag}\{\mathbf{s}^{\mathrm{T}}_{k-m}\mathbf{y}_{k-m}, \ldots, \mathbf{s}^{\mathrm{T}}_{k-1}\mathbf{y}_{k-1}\} \in \mathbb{R}^{m\times m}
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 059. ### 3.2. Compact representation of limited-memory BFGS matrix（explain.md:L435-L437）

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
\mathbf{L}_k \in \mathbb{R}^{m\times m}, \quad (\mathbf{L}_k)_{i,j} = \begin{cases} \mathbf{s}^{\mathrm{T}}_{i+k-m-1}\mathbf{y}_{j+k-m-1} & \text{if } i > j, \\ 0 & \text{if } i \leq j. \end{cases}
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 060. ### 3.2. Compact representation of limited-memory BFGS matrix（explain.md:L441-L443）

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
\mathbf{B}^0_k = \hat{\mathbf{K}}^{(k)} = \begin{bmatrix} \mathbf{K}_{uu} & \mathbf{0} \\ \mathbf{0} & \mathbf{K}_{dd} \end{bmatrix}^{(k)}. \quad (21)
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 061. ### 3.3. Generalized Cauchy point（explain.md:L455-L457）

- 对应关系说明：本公式位于广义Cauchy点层：分段路径、断点和活动集更新。
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
t_i = \begin{cases} (x^0_i - \mathrm{ub}_i)/r_i & \text{if } r_i < 0, \\ (x^0_i - \mathrm{lb}_i)/r_i & \text{if } r_i > 0, \\ +\infty & \text{if } r_i = 0. \end{cases}
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 062. ### 3.3. Generalized Cauchy point（explain.md:L461-L463）

- 对应关系说明：本公式位于广义Cauchy点层：分段路径、断点和活动集更新。
- 最底层代码链接：
  - [main.cc:L1401-L1442](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L1401-L1442)
  - [main.cc:L4775-L4948](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L4775-L4948)

### 代码片段 1（main.cc:L1401-L1442）
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

### 代码片段 2（main.cc:L4775-L4948）
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

### 公式（与 explain.md 一致）
\[
x_i(t) = x_i^{0} - \min \{t,t_i\} r_i.
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 063. ### 3.3. Generalized Cauchy point（explain.md:L467-L469）

- 对应关系说明：本公式位于广义Cauchy点层：分段路径、断点和活动集更新。
- 最底层代码链接：
  - [main.cc:L1401-L1442](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L1401-L1442)
  - [main.cc:L4775-L4948](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L4775-L4948)

### 代码片段 1（main.cc:L1401-L1442）
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

### 代码片段 2（main.cc:L4775-L4948）
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

### 公式（与 explain.md 一致）
\[
\mathbf{x}(t) = \mathbf{x}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)},
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。
