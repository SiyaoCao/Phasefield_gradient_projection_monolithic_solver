# explain.md 与 main.cc 公式-代码对应（公式块 062）

- 所属章节：`3.3. Generalized Cauchy point`
- explain.md 行号：`461`

论文公式：

\[
x_i(t) = x_i^{0} - \min \{t,t_i\} r_i.
\]

对应 `main.cc` 代码：

```cpp
// main.cc:1401-1442
    PhaseFieldMonolithicSolve<dim>::calculate_break_points(const BlockVector<double> & solution_delta,
		       				           const BlockVector<double> & gradient_g,
							   BlockVector<double> & gradient_d)
  {
    // Creates a min heap of break points
    std::priority_queue< std::pair<double, unsigned int>,
                         std::vector<std::pair<double, unsigned int>>,
    		         std::greater<std::pair<double, unsigned int>> >
    break_points_sorted;

    double t = 0.0;

    Vector<double> solution_phasefield_total(m_solution.block(m_d_dof));
    solution_phasefield_total += solution_delta.block(m_d_dof);

    // upper bound is 1.0, lower bound is the solution at the previous step.
    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
      {
	if (gradient_g.block(m_d_dof)[i] < 0)
	  t = (solution_phasefield_total[i] - 1.0 ) / gradient_g.block(m_d_dof)[i];
	else if (gradient_g.block(m_d_dof)[i] > 0)
	  t = solution_delta.block(m_d_dof)[i] / gradient_g.block(m_d_dof)[i];
	else
	  t = std::numeric_limits<double>::max();

        //AssertThrow(t >= 0, ExcMessage("Break point has to be a non-negative t value"));

        if (t > 0)
          {
	    break_points_sorted.push(std::make_pair(t, i));
          }
        else // if t == 0, i is in the active set
          {
            gradient_d.block(m_d_dof)[i] = 0;
            if (gradient_g.block(m_d_dof)[i] > 0)
              m_active_set_phasefield(i) = 1; //lower bound
            else
              m_active_set_phasefield(i) = 2; //upper bound
          }
      }

    return break_points_sorted;

// main.cc:4855-4943
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

	t_old = t;

	top_pair = t_series.top();
	t = top_pair.first;
	b = top_pair.second;

	delta_t = t - t_old;
      }

    if (delta_t_min < 0)
      delta_t_min = 0;

    t_old += delta_t_min;

    for(unsigned int i = 0; i < m_active_set_phasefield.size(); ++i)
      {
	// inactive phasefield dof
	if (m_active_set_phasefield(i) < 0.5)
	  solution_delta_cauchy_point.block(m_d_dof)[i] = solution_delta.block(m_d_dof)[i]
						+ t_old * gradient_d.block(m_d_dof)[i];
      }

    // There are no active constraints in the displacement field
    solution_delta_cauchy_point.block(m_u_dof) = solution_delta.block(m_u_dof);
    (solution_delta_cauchy_point.block(m_u_dof)).add(t_old, gradient_d.block(m_u_dof));
```

简要说明：对应变量或算子在 main.cc 中有同名或等价离散实现。
