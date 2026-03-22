# explain.md 与 main.cc 公式-代码对应（公式块 031）

- 所属章节：`3.1. Algorithm overview`
- explain.md 行号：`259`

论文公式：

\[
m_{k}(\pmb{x}) = \Pi_{k} + \pmb{r}_{k}^{\mathrm{T}}(\pmb{x} - \pmb{x}_{k}) + \frac{1}{2} (\pmb{x} - \pmb{x}_{k})^{\mathrm{T}}\mathbf{B}_{k}(\pmb{x} - \pmb{x}_{k}) \quad (15)
\]

对应 `main.cc` 代码：

```cpp
// main.cc:4953-4977
  solve_nonlinear_timestep_LBFGS_B(BlockVector<double> & solution_delta,
				   BlockVector<double> & LBFGS_update_refine)
  {
    BlockVector<double> LBFGS_update(m_dofs_per_block);
    BlockVector<double> solution_delta_cauchy_point(m_dofs_per_block);
    LBFGS_update = 0.0;

    const unsigned int LBFGS_m = m_parameters.m_LBFGS_m;

    unsigned int LBFGS_iteration = 0;

    m_error_residual.reset();
    m_error_residual_0.reset();
    m_error_residual_norm.reset();
    m_error_update.reset();
    m_error_update_0.reset();
    m_error_update_norm.reset();

    if (m_parameters.m_output_iteration_history)
      print_conv_header_LBFGSB();

    BlockVector<double> LBFGS_s_vector(m_dofs_per_block);
    BlockVector<double> LBFGS_y_vector(m_dofs_per_block);
    BlockVector<double> free_dofs(m_dofs_per_block);
    BlockVector<double> b0xs_vector(m_dofs_per_block);

// main.cc:6117-6139
  double PhaseFieldMonolithicSolve<dim>::calculate_energy_functional() const
  {
    double energy_functional = 0.0;

    FEValues<dim> fe_values(m_fe, m_qf_cell, update_JxW_values);

    for (const auto &cell : m_dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);

        const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
          m_quadrature_point_history.get_data(cell);
        Assert(lqph.size() == m_n_q_points, ExcInternalError());

        for (unsigned int q_point = 0; q_point < m_n_q_points; ++q_point)
          {
            const double JxW = fe_values.JxW(q_point);
            energy_functional += lqph[q_point]->get_total_strain_energy() * JxW;
            energy_functional += lqph[q_point]->get_crack_energy_dissipation() * JxW;
          }
      }

    return energy_functional;
```

简要说明：对应变量或算子在 main.cc 中有同名或等价离散实现。
