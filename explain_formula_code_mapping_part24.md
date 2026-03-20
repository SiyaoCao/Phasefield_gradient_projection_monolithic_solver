## 088. #### 3.4.1. Direct matrix factorization for the primal approach（explain.md:L631-L633）

- 对应关系说明：本公式位于子空间最小化层：Z投影、受限自由变量求解与线性代数实现。
- 最底层代码链接：
  - [main.cc:L1259-L1363](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L1259-L1363)
  - [main.cc:L5221-L5635](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L5221-L5635)

### 代码片段 1（main.cc:L1259-L1363）
```cpp
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
```

### 代码片段 2（main.cc:L5221-L5635）
```cpp
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
```

### 公式（与 explain.md 一致）
\[
\hat{\mathbf{x}} = -\hat{\mathbf{B}}^{-1}_k \hat{\mathbf{r}}_k,
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 089. #### 3.4.1. Direct matrix factorization for the primal approach（explain.md:L634-L636）

- 对应关系说明：本公式位于子空间最小化层：Z投影、受限自由变量求解与线性代数实现。
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
\hat{\mathbf{B}}^{-1}_k = (\hat{\mathbf{B}}^0_k)^{-1} + (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \left[\mathbf{I} - \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \right]^{-1} \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1}, \quad (29)
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 090. #### 3.4.1. Direct matrix factorization for the primal approach（explain.md:L639-L641）

- 对应关系说明：本公式位于子空间最小化层：Z投影、受限自由变量求解与线性代数实现。
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
\mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k = \mathbf{Z}^{\mathrm{T}}_k [\mathbf{Y}_k \quad \mathbf{B}^0_k \mathbf{S}_k] = [\mathbf{Z}^{\mathrm{T}}_k \mathbf{y}_{k-m} \cdots \mathbf{Z}^{\mathrm{T}}_k \mathbf{y}_{k-1} \quad \mathbf{Z}^{\mathrm{T}}_k(\mathbf{B}^0_k \mathbf{s}_{k-m}) \cdots \mathbf{Z}^{\mathrm{T}}_k(\mathbf{B}^0_k \mathbf{s}_{k-1})].
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 091. #### 3.4.1. Direct matrix factorization for the primal approach（explain.md:L644-L646）

- 对应关系说明：本公式位于子空间最小化层：Z投影、受限自由变量求解与线性代数实现。
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
\mathbf{I} - \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \in \mathbb{R}^{2m\times 2m}
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。
