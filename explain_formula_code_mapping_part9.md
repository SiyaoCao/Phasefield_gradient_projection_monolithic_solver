## 3. Gradient projection based monolithic scheme > 3.1. Algorithm overview

### 未编号公式 #47 — 未编号公式：\left(\mathbf{C}^{\mathrm{...
- 论文位置：[3. Gradient projection based monolithic scheme > 3.1. Algorithm overview，`explain.md` 第 364 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L364)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\left(\mathbf{C}^{\mathrm{T}}\mathbf{A}\mathbf{C} + \mathbf{I}_{d_{c}}\right)\hat{\mathbf{x}} = \mathbf{C}^{\mathrm{T}}(\mathbf{b} - \mathbf{A}\mathbf{k})
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
### 未编号公式 #48 — 未编号公式：\mathbf{x} = \mathbf{C}\ha...
- 论文位置：[3. Gradient projection based monolithic scheme > 3.1. Algorithm overview，`explain.md` 第 370 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L370)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\mathbf{x} = \mathbf{C}\hat{\mathbf{x}} + \mathbf{k}.
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
