# explain.md 与 main.cc 公式-代码对应（公式块 083）

- 所属章节：`3.4.1. Direct matrix factorization for the primal approach`
- explain.md 行号：`595`

论文公式：

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

对应 `main.cc` 代码：

```cpp
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


// main.cc:5485-5501
	    zT_B0_z(free_dofs, m_tangent_matrix);

	    SparseDirectUMFPACK zT_B0_z_inv;
	    zT_B0_z_inv.initialize(m_tangent_matrix);

	    //SparseDirectUMFPACK zT_B0_z_inv_disp;
	    //zT_B0_z_inv_disp.initialize(m_tangent_matrix.block(m_u_dof, m_u_dof));

	    //SparseDirectUMFPACK zT_B0_z_inv_phasefield;
	    //zT_B0_z_inv_phasefield.initialize(m_tangent_matrix.block(m_d_dof, m_d_dof));

	    m_timer.leave_subsection();

	    m_timer.enter_subsection("Subspace direct solve (LU matrix-vector multiplication)");

	    zT_B0_z_inv.vmult(search_direction, rhs_vector);
	    //zT_B0_z_inv_disp.vmult(search_direction.block(m_u_dof), rhs_vector.block(m_u_dof));

// main.cc:5564-5621
		FullMatrix<double> wT_z_zT_B0_z_inv_zT_w(2 * list_size);
		const auto itr_zT_y_list_begin = zT_y_list.begin();
		const auto itr_zT_b0xs_list_begin = zT_b0xs_list.begin();
		const auto itr_zT_B0_z_inv_zT_y_list_begin = zT_B0_z_inv_zT_y_list.begin();
		const auto itr_zT_B0_z_inv_zT_b0xs_list_begin = zT_B0_z_inv_zT_b0xs_list.begin();
		for (unsigned int i = 0; i < list_size; ++i)
		  for (unsigned int j = 0; j < list_size; ++j)
		    {
		      wT_z_zT_B0_z_inv_zT_w(i          , j          ) = (*std::next(itr_zT_y_list_begin            , i))
								      * (*std::next(itr_zT_B0_z_inv_zT_y_list_begin, j));

		      wT_z_zT_B0_z_inv_zT_w(i          , j+list_size) = (*std::next(itr_zT_y_list_begin               , i))
								      * (*std::next(itr_zT_B0_z_inv_zT_b0xs_list_begin, j));

		      wT_z_zT_B0_z_inv_zT_w(i+list_size, j          ) = (*std::next(itr_zT_b0xs_list_begin         , i))
								      * (*std::next(itr_zT_B0_z_inv_zT_y_list_begin, j));

		      wT_z_zT_B0_z_inv_zT_w(i+list_size, j+list_size) = (*std::next(itr_zT_b0xs_list_begin            , i))
								      * (*std::next(itr_zT_B0_z_inv_zT_b0xs_list_begin, j));
		    }

		FullMatrix<double> temp_matrix(2 * list_size);
		M_matrix.mmult(temp_matrix, wT_z_zT_B0_z_inv_zT_w);

		FullMatrix<double> middle_matrix(IdentityMatrix(2*list_size));
		middle_matrix.add(-1.0, temp_matrix);

		FullMatrix<double> middle_matrix_inv(2 * list_size);
		middle_matrix_inv.invert(middle_matrix);

		middle_matrix_inv.mmult(middle_matrix, M_matrix);

		Vector<double> wT_z_zT_B0_z_inv_rhs(2 * list_size);
		for (unsigned int i = 0; i < list_size; ++i)
		  {
		    wT_z_zT_B0_z_inv_rhs(i            ) = (*std::next(itr_zT_B0_z_inv_zT_y_list_begin   , i)) * rhs_vector;
		    wT_z_zT_B0_z_inv_rhs(i + list_size) = (*std::next(itr_zT_B0_z_inv_zT_b0xs_list_begin, i)) * rhs_vector;
		  }

		Vector<double> middle_matrix_wT_z_zT_B0_z_inv_rhs(2 * list_size);
		middle_matrix.vmult(middle_matrix_wT_z_zT_B0_z_inv_rhs,
				    wT_z_zT_B0_z_inv_rhs);

		unsigned int index = 0;
		for (auto itr = zT_B0_z_inv_zT_y_list.begin(); itr != zT_B0_z_inv_zT_y_list.end(); ++itr)
		  {
		    update_vector.add(middle_matrix_wT_z_zT_B0_z_inv_rhs(index), *itr);
		    ++index;
		  }
		for (auto itr = zT_B0_z_inv_zT_b0xs_list.begin(); itr != zT_B0_z_inv_zT_b0xs_list.end(); ++itr)
		  {
		    update_vector.add(middle_matrix_wT_z_zT_B0_z_inv_rhs(index), *itr);
		    ++index;
		  }
	      } //	if (list_size > 0)

	    search_direction += update_vector;
```

简要说明：对应变量或算子在 main.cc 中有同名或等价离散实现。
