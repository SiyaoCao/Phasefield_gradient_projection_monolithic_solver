# explain.md 与 main.cc 公式-代码对应（公式块 106）

- 所属章节：`4.1. Cyclic tension-compression test`
- explain.md 行号：`780`

论文公式：

\[
E_{\mathrm{crack}} = g_{\mathrm{c}}\Gamma_{l}(d) = g_{\mathrm{c}}\int_{\Omega}\gamma (d,\nabla d)\mathrm{d}\Omega = g_{\mathrm{c}}\int_{\Omega}\frac{1}{2l}\left(d^{2} + l^{2}\nabla d\cdot \nabla d\right)\mathrm{d}\Omega \quad (33)
\]

对应 `main.cc` 代码：

```cpp
// main.cc:886-887
    m_crack_energy_dissipation = m_gc * (  0.5 / m_length_scale * m_phase_field_value * m_phase_field_value
	                                   + 0.5 * m_length_scale * m_grad_phasefield * m_grad_phasefield)

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
