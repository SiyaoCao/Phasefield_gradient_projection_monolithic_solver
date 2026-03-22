# explain.md 与 main.cc 公式-代码对应（公式块 022）

- 所属章节：`2.2. Finite element discretization`
- explain.md 行号：`184`

论文公式：

\[
\begin{array}{rl} 
\Pi (\pmb {u}_{A},d_{A}) = &\int_{\Omega}\psi \left(\pmb{\epsilon}(\pmb {N}_{u_{A}}\pmb {u}_{A}),N_{d_{A}}d_{A}\right)\mathrm{d}\Omega \\ 
&+\int_{\Omega}\frac{g_{c}}{2l}\left((N_{d_{A}}d_{A})^{2} + l^{2}(\nabla N_{d_{A}}d_{A})\cdot (\nabla N_{d_{A}}d_{A})\right)\mathrm{d}\Omega \\ 
&-\int_{\Omega}\pmb {b}\cdot (\pmb {N}_{u_{A}}\pmb {u}_{A})\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot (\pmb {N}_{u_{A}}\pmb {u}_{A})\mathrm{d}\Gamma, 
\end{array} \quad (8)
\]

对应 `main.cc` 代码：

```cpp
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
