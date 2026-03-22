# explain.md 与 main.cc 公式-代码对应（公式块 025）

- 所属章节：`2.2. Finite element discretization`
- explain.md 行号：`210`

论文公式：

\[
\begin{array}{rl} 
\mathbf{K}_{u_{A}u_{B}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},\frac{\partial \pmb{\sigma}}{\partial \pmb{\epsilon}}:\nabla^{(s)}\pmb{N}_{u_{B}}\right),\qquad \mathbf{K}_{u_{A}d_{B}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},g^{\prime}(d)\pmb{\sigma}^{+}N_{d_{B}}\right),\\ 
\mathbf{K}_{d_{A}u_{B}} = \left(N_{d_{A}},g^{\prime}(d)\pmb{\sigma}^{+}:\nabla^{(s)}\pmb{N}_{u_{B}}\right),\qquad \mathbf{K}_{d_{A}d_{B}} = \left(N_{d_{A}},\left(\frac{g_{c}}{l} +g^{\prime \prime}(d)\psi^{+}\right)N_{d_{B}}\right) + \left(\nabla N_{d_{A}},g_{c}l\nabla N_{d_{B}}\right). 
\end{array} \quad (10)
\]

对应 `main.cc` 代码：

```cpp
// main.cc:3860-3952
  void PhaseFieldMonolithicSolve<dim>::assemble_system_B0_one_cell(
      const typename DoFHandler<dim>::active_cell_iterator &cell,
      ScratchData_ASM & scratch,
      PerTaskData_ASM & data) const
  {
    data.reset();
    scratch.reset();
    scratch.m_fe_values.reinit(cell);
    cell->get_dof_indices(data.m_local_dof_indices);

    scratch.m_fe_values[m_d_fe].get_function_values(
      scratch.m_solution_previous_step, scratch.m_phasefield_previous_step_cell);

    const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
      m_quadrature_point_history.get_data(cell);
    Assert(lqph.size() == m_n_q_points, ExcInternalError());

    const double delta_time = m_time.get_delta_t();

    for (const unsigned int q_point : scratch.m_fe_values.quadrature_point_indices())
      {
        for (const unsigned int k : scratch.m_fe_values.dof_indices())
          {
            const unsigned int k_group = m_fe.system_to_base_index(k).first.first;

            if (k_group == m_u_dof)
              {
                scratch.m_Nx_disp[q_point][k] =
                  scratch.m_fe_values[m_u_fe].value(k, q_point);
                scratch.m_grad_Nx_disp[q_point][k] =
                  scratch.m_fe_values[m_u_fe].gradient(k, q_point);
                scratch.m_symm_grad_Nx_disp[q_point][k] =
                  symmetrize(scratch.m_grad_Nx_disp[q_point][k]);
              }
            else if (k_group == m_d_dof)
              {
		scratch.m_Nx_phasefield[q_point][k] =
		  scratch.m_fe_values[m_d_fe].value(k, q_point);
		scratch.m_grad_Nx_phasefield[q_point][k] =
		  scratch.m_fe_values[m_d_fe].gradient(k, q_point);
              }
            else
              Assert(k_group <= m_d_dof, ExcInternalError());
          }
      }

    for (const unsigned int q_point : scratch.m_fe_values.quadrature_point_indices())
      {
	const double length_scale            = lqph[q_point]->get_length_scale();
	const double gc                      = lqph[q_point]->get_critical_energy_release_rate();
	const double eta                     = lqph[q_point]->get_viscosity();
	const double current_positive_strain_energy = lqph[q_point]->get_current_positive_strain_energy();

	const double phasefield_value        = lqph[q_point]->get_phase_field_value();

        const std::vector<double>         &      N_phasefield = scratch.m_Nx_phasefield[q_point];
        const std::vector<Tensor<1, dim>> & grad_N_phasefield = scratch.m_grad_Nx_phasefield[q_point];

        //const SymmetricTensor<2, dim> & cauchy_stress_positive = lqph[q_point]->get_cauchy_stress_positive();
        const SymmetricTensor<4, dim> & mechanical_C  = lqph[q_point]->get_mechanical_C();

        const std::vector<SymmetricTensor<2, dim>> & symm_grad_N_disp =
          scratch.m_symm_grad_Nx_disp[q_point];
        const double JxW = scratch.m_fe_values.JxW(q_point);

        SymmetricTensor<2, dim> symm_grad_Nx_i_x_C;

        for (const unsigned int i : scratch.m_fe_values.dof_indices())
          {
            const unsigned int i_group = m_fe.system_to_base_index(i).first.first;

            if (i_group == m_u_dof)
              {
                symm_grad_Nx_i_x_C = symm_grad_N_disp[i] * mechanical_C;
              }

            for (const unsigned int j : scratch.m_fe_values.dof_indices())
              {
                const unsigned int j_group = m_fe.system_to_base_index(j).first.first;

                if ((i_group == j_group) && (i_group == m_u_dof))
                  {
                    data.m_cell_matrix(i, j) += symm_grad_Nx_i_x_C * symm_grad_N_disp[j] * JxW;
                  }
                else if ((i_group == j_group) && (i_group == m_d_dof))
                  {
                    data.m_cell_matrix(i, j) += (  (   gc/length_scale + eta/delta_time
                	                             + degradation_function_2nd_order_derivative(phasefield_value)
						     * current_positive_strain_energy  )
                	                          * N_phasefield[i] * N_phasefield[j]
					          + gc * length_scale * grad_N_phasefield[i] * grad_N_phasefield[j]
					        ) * JxW;
                  }
```

简要说明：对应有限元离散中的块矩阵装配或约化矩阵表达。
