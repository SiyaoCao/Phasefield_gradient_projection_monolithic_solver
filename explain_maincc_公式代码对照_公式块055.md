# explain.md 与 main.cc 公式-代码对应（公式块 055）

- 所属章节：`3.2. Compact representation of limited-memory BFGS matrix`
- explain.md 行号：`414`

论文公式：

\[
\mathbf{B}_k = \mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k, \quad (20)
\]

对应 `main.cc` 代码：

```cpp
// main.cc:5113-5126
	// assemble the initial B_0 matrix at the k-th L-BFGS iteration
	// m_solution is the old solution from the previous converged step
	// it is needed only for the viscosity term
	// the output is m_tangent_matrix (B^0_k)
	assemble_system_B0(m_solution);

	// B^0_k * s_vector has to be completely recalculated from scratch
	// at each L-BFGS iteration, since B^0_k is different
	b0xs_vector_list.clear();
	for (auto itr = s_vector_list.begin(); itr != s_vector_list.end(); ++itr)
	  {
	    m_tangent_matrix.vmult(b0xs_vector, *itr);
	    b0xs_vector_list.push_back(b0xs_vector);
	  }

// main.cc:5166-5182
	FullMatrix<double> M_matrix_inv(2 * list_size);
	FullMatrix<double> M_matrix(2 * list_size);

	M_matrix_inv = 0;
	for (unsigned int i = 0; i < list_size; ++i)
	  M_matrix_inv(i, i) = -D_matrix(i, i);

	for (unsigned int i = 0; i < list_size; ++i)
          for (unsigned int j = 0; j < list_size; ++j)
            {
              M_matrix_inv(i + list_size, j + list_size) = sTxBxs_matrix(i, j);
              M_matrix_inv(i + list_size, j            ) = L_matrix(i, j);
              M_matrix_inv(i            , j + list_size) = L_matrix(j, i);
            }

	if (!M_matrix_inv.empty())
	  M_matrix.invert(M_matrix_inv);
```

简要说明：对应变量或算子在 main.cc 中有同名或等价离散实现。
