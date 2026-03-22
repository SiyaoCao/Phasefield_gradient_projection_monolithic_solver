# explain.md 与 main.cc 公式-代码对应（公式块 095）

- 所属章节：`3.4.3. Schur complement for the dual approach`
- explain.md 行号：`676`

论文公式：

\[
x_i = x^c_i, \quad \forall i \in \mathcal{A}_k(\mathbf{x}^c),
\]

对应 `main.cc` 代码：

```cpp
// main.cc:5193-5220
	// We need to find out which DOFs are free:
	// no essential boundary conditions, no hanging node constraints
	// no active box constraints
	unsigned int free_disp_number = 0;
	for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
	  {
	    if (m_constraints.is_constrained(i))
	      free_dofs.block(m_u_dof)[i] = -1;
	    else
	      {
	        free_dofs.block(m_u_dof)[i] = 1;
	        ++free_disp_number;
	      }
	  }

	unsigned int free_phasefield_number = 0;
	for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
	  {
	    if (   m_constraints.is_constrained(i + m_dofs_per_block[m_u_dof])
		|| m_active_set_phasefield(i) > 0.5)
	      free_dofs.block(m_d_dof)[i] = -1;
	    else
	      {
	        free_dofs.block(m_d_dof)[i] = 1;
	        ++free_phasefield_number;
	      }
	  }


// main.cc:5250-5272
	// temp_vector_2 = B_0 * (x^c - x_k) - W * M * W^T * (x^c - x_k)
	if (list_size > 0)
	  temp_vector_2 -= temp_vector_5;

	// temp_vector_2 = g + B_0 * (x^c - x_k) - W * M * W^T * (x^c - x_k)
	temp_vector_2 += m_system_rhs;

	// temp_vector_2 = Z^T * [g + B_0 * (x^c - x_k) - W * M * W^T * (x^c - x_k)]
	for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
	  {
	    if (free_dofs.block(m_u_dof)[i] < 0)
	      temp_vector_2.block(m_u_dof)[i] = 0;
	  }

	for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
	  {
	    if (free_dofs.block(m_d_dof)[i] < 0)
	      temp_vector_2.block(m_d_dof)[i] = 0;
	  }

	BlockVector<double> rhs_vector(temp_vector_2);
	rhs_vector *= -1;


// main.cc:5280-5382
	if (m_parameters.m_type_linear_solver == "CG")
	  {
	    m_timer.enter_subsection("Subspace CG solve (inverse operator)");

	    //const double rc_hat_norm = rhs_vector.l2_norm();
	    const double cg_tol = m_parameters.m_CG_tolerace; //std::min( 0.1, std::sqrt(rc_hat_norm) ) * rc_hat_norm;

	    zT_B0_z(free_dofs, m_tangent_matrix);

	    const auto op_zT_B0_z = block_operator(m_tangent_matrix);

	    if (list_size > 0)
	      {
		std::list<BlockVector<double>> zT_y_list;
		BlockVector<double> zT_y_vector(m_dofs_per_block);
		for (auto itr = y_vector_list.begin(); itr != y_vector_list.end(); ++itr)
		  {
		    zT_y_vector = (*itr);
		    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
		      {
			if (free_dofs.block(m_u_dof)[i] < 0)
			  zT_y_vector.block(m_u_dof)[i] = 0;
		      }

		    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
		      {
			if (free_dofs.block(m_d_dof)[i] < 0)
			  zT_y_vector.block(m_d_dof)[i] = 0;
		      }

		    zT_y_list.push_back(zT_y_vector);
		  }

		std::list<BlockVector<double>> zT_b0xs_list;
		BlockVector<double> zT_b0xs_vector(m_dofs_per_block);
		for (auto itr = b0xs_vector_list.begin(); itr != b0xs_vector_list.end(); ++itr)
		  {
		    zT_b0xs_vector = (*itr);
		    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
		      {
			if (free_dofs.block(m_u_dof)[i] < 0)
			  zT_b0xs_vector.block(m_u_dof)[i] = 0;
		      }

		    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
		      {
			if (free_dofs.block(m_d_dof)[i] < 0)
			  zT_b0xs_vector.block(m_d_dof)[i] = 0;
		      }

		    zT_b0xs_list.push_back(zT_b0xs_vector);
		  }

		const auto op_M_matrix = linear_operator(M_matrix);

		FullMatrix<double> zT_W_matrix_u(m_dofs_per_block[m_u_dof], 2*list_size);
		unsigned int j = 0;
		for (auto itr = zT_y_list.begin(); itr != zT_y_list.end(); ++itr)
		  {
		    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
		      zT_W_matrix_u(i, j) = (*itr).block(m_u_dof)[i];
		    ++j;
		  }
		j = 0;
		for (auto itr = zT_b0xs_list.begin(); itr != zT_b0xs_list.end(); ++itr)
		  {
		    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
		      zT_W_matrix_u(i, j + list_size) = (*itr).block(m_u_dof)[i];
		    ++j;
		  }

		FullMatrix<double> zT_W_matrix_d(m_dofs_per_block[m_d_dof], 2*list_size);
		j = 0;
		for (auto itr = zT_y_list.begin(); itr != zT_y_list.end(); ++itr)
		  {
		    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
		      zT_W_matrix_d(i, j) = (*itr).block(m_d_dof)[i];
		    ++j;
		  }
		j = 0;
		for (auto itr = zT_b0xs_list.begin(); itr != zT_b0xs_list.end(); ++itr)
		  {
		    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
		      zT_W_matrix_d(i, j + list_size) = (*itr).block(m_d_dof)[i];
		    ++j;
		  }

		const auto op_zT_W_matrix_u = linear_operator(zT_W_matrix_u);
		const auto op_zT_W_matrix_d = linear_operator(zT_W_matrix_d);

		const auto op_uMuT = op_zT_W_matrix_u * op_M_matrix * transpose_operator(op_zT_W_matrix_u);

		const auto op_uMdT = op_zT_W_matrix_u * op_M_matrix * transpose_operator(op_zT_W_matrix_d);

		const auto op_dMuT = op_zT_W_matrix_d * op_M_matrix * transpose_operator(op_zT_W_matrix_u);

		const auto op_dMdT = op_zT_W_matrix_d * op_M_matrix * transpose_operator(op_zT_W_matrix_d);

		const auto op_zT_wMwT_z = block_operator<2, 2, BlockVector<double>>({op_uMuT, op_uMdT,
										     op_dMuT, op_dMdT});

		const auto op_total = op_zT_B0_z - op_zT_wMwT_z;
```

简要说明：对应变量或算子在 main.cc 中有同名或等价离散实现。
