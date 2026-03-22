# explain.md 与 main.cc 公式-代码对应（公式块 034）

- 所属章节：`3.1. Algorithm overview`
- explain.md 行号：`281`

论文公式：

\[
p_k(t) = m_k(\pmb{x}(t)) = \Pi_k + \pmb{r}_k^\mathrm{T}(\pmb{x}(t) - \pmb{x}_k) + \frac{1}{2} (\pmb{x}(t) - \pmb{x}_k)^\mathrm{T}\mathbf{B}_k(\pmb{x}(t) - \pmb{x}_k). \quad (18)
\]

对应 `main.cc` 代码：

```cpp
// main.cc:4775-4835
  calculate_cauchy_point(const BlockSparseMatrix<double> & B0_matrix,
	                 const std::list<BlockVector<double>> & y_vector_list,
		         const std::list<BlockVector<double>> & b0xs_vector_list,
			 const FullMatrix<double> & M_matrix,
			 const BlockVector<double> & gradient_g,
			 const BlockVector<double> & solution_delta,
			 BlockVector<double> & solution_delta_cauchy_point)
  {
    m_timer.enter_subsection("Calculate Cauchy point");

    solution_delta_cauchy_point = 0.0;
    BlockVector<double> gradient_d(gradient_g);
    gradient_d *= -1;

    const unsigned int list_size = y_vector_list.size();
    const auto itr_y_begin    = y_vector_list.begin();
    const auto itr_b0xs_begin = b0xs_vector_list.begin();

    // t_series only contains t > 0
    std::priority_queue< std::pair<double, unsigned int>,
                         std::vector<std::pair<double, unsigned int>>,
        		 std::greater<std::pair<double, unsigned int>> >
    t_series = calculate_break_points(solution_delta,
    			              gradient_g,
				      gradient_d);

    // m_active_set_phasefield contains 1 or 2 for active set and 0 for inactive set
    for (unsigned int i = 0; i < m_active_set_phasefield.size(); ++i)
      {
	if (m_active_set_phasefield(i) > 0.5)
	  solution_delta_cauchy_point.block(m_d_dof)[i] = solution_delta.block(m_d_dof)[i];
      }

    // p = W^T * d
    Vector<double> p(2 * list_size);
    for (unsigned int i = 0; i < list_size; ++i)
      {
        p(i)             = (*std::next(itr_y_begin,    i)) * gradient_d;
        p(i + list_size) = (*std::next(itr_b0xs_begin, i)) * gradient_d;
      }

    Vector<double> c(2 * list_size);
    c = 0.0;

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
```

简要说明：对应变量或算子在 main.cc 中有同名或等价离散实现。
