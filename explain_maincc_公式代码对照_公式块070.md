# explain.md 与 main.cc 公式-代码对应（公式块 070）

- 所属章节：`3.3. Generalized Cauchy point`
- explain.md 行号：`519`

论文公式：

\[
f'_j = \mathbf{r}^{\mathrm{T}}\mathbf{d}^{(j)} + \mathbf{d}^{(j)}{}^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j)} = f'_{j-1} + \Delta t^{(j-1)} f''_{j-1} + r^2_b + r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j)} \quad (23)
\]

对应 `main.cc` 代码：

```cpp
// main.cc:4819-4917
    double f_prime = -(gradient_d * gradient_d);

    // M * p
    Vector<double> Mp(2 * list_size);
    if (list_size > 0)
      M_matrix.vmult(Mp, p);

    // B_0 * d
    BlockVector<double> B0_grandient_d(m_dofs_per_block);
    B0_matrix.vmult(B0_grandient_d, gradient_d);

    double f_prime_prime = gradient_d * B0_grandient_d;
    if (list_size > 0)
      f_prime_prime -= (p * Mp);

    double delta_t_min = -f_prime / f_prime_prime;

    double t_old = 0.0;

    std::pair<double, unsigned int> top_pair = t_series.top();
    double t = top_pair.first;
    unsigned int b = top_pair.second;

    double delta_t = t - t_old;

    BlockVector<double> z(m_dofs_per_block);
    z = 0.0;

    // w_b = W^T * e_b
    Vector<double> w_b(2 * list_size);

    // w_b_T_x_M = w_b^T * M
    Vector<double> w_b_T_x_M(2 * list_size);

    Vector<double> temp_vector(m_dofs_per_block[m_d_dof]);

    while (delta_t_min >= delta_t)
      {
	t_series.pop();

	if (gradient_d.block(m_d_dof)[b] > 0)
	  solution_delta_cauchy_point.block(m_d_dof)[b] = 1.0 - m_solution.block(m_d_dof)[b];
	else if (gradient_d.block(m_d_dof)[b] < 0)
	  solution_delta_cauchy_point.block(m_d_dof)[b] = 0.0;
	else
	  AssertThrow(false,
	              ExcMessage("gradient_d(b) cannot be zero!"));

	if (gradient_d.block(m_d_dof)[b] < 0)
	  m_active_set_phasefield[b] = 1; //lower bound
	else
	  m_active_set_phasefield[b] = 2; //upper bound

        // z_{j} = z_{j-1} + delta_t_{j-1} * gradient_d_{j-1};
	z.sadd(1.0, delta_t, gradient_d);

	// c_{j} = c_{j-1} + delta_t_{j-1} * p_{j-1};
	if (list_size > 0)
	  c.sadd(1.0, delta_t, p);

        double temp_scalar = ebT_x_B0_x_v(b, B0_matrix, z);

        // w_b = W^T * e_b
        for (unsigned int i = 0; i < list_size; ++i)
          {
            w_b(i)             = (*std::next(itr_y_begin,    i)).block(m_d_dof)[b];
            w_b(i + list_size) = (*std::next(itr_b0xs_begin, i)).block(m_d_dof)[b];
          }

        if (list_size > 0)
          M_matrix.vmult(w_b_T_x_M, w_b);

	f_prime += delta_t * f_prime_prime
	         + gradient_g.block(m_d_dof)[b] * gradient_g.block(m_d_dof)[b]
	         + temp_scalar * gradient_g.block(m_d_dof)[b];

	if (list_size > 0)
	  f_prime -= (w_b_T_x_M * c) * gradient_g.block(m_d_dof)[b];

	temp_scalar = ebT_x_B0_x_v(b, B0_matrix, gradient_d);

	f_prime_prime += 2.0 * gradient_g.block(m_d_dof)[b] * temp_scalar
	               + gradient_g.block(m_d_dof)[b] * gradient_g.block(m_d_dof)[b]
		             * B0_matrix.block(m_d_dof, m_d_dof)(b, b);

	if (list_size > 0)
	  {
	    f_prime_prime -= 2.0 * gradient_g.block(m_d_dof)[b] * (w_b_T_x_M * p);
	    f_prime_prime -= gradient_g.block(m_d_dof)[b] * gradient_g.block(m_d_dof)[b] * (w_b_T_x_M * w_b);
	  }

	// p_{j} = p_{j-1} + g_b * w_b;
	if (list_size > 0)
	  p.sadd(1.0, gradient_g.block(m_d_dof)[b], w_b);

	gradient_d.block(m_d_dof)[b] = 0.0;

	delta_t_min = -f_prime / f_prime_prime;
```

简要说明：对应变量或算子在 main.cc 中有同名或等价离散实现。
