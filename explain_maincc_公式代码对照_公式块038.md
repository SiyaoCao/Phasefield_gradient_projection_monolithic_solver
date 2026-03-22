# explain.md 与 main.cc 公式-代码对应（公式块 038）

- 所属章节：`3.1. Algorithm overview`
- explain.md 行号：`309`

论文公式：

\[
\pmb{p}_k = \pmb{x}^* - \pmb{x}_k.
\]

对应 `main.cc` 代码：

```cpp
// main.cc:5629-5689
	// We don't do backtrack yet. We will make sure phasefield
	// remains feasible later
	alpha_backtrack = 1.0;
	search_direction *= alpha_backtrack;

	LBFGS_update.block(m_u_dof) = solution_delta_cauchy_point.block(m_u_dof);
	LBFGS_update.block(m_u_dof) += search_direction.block(m_u_dof);
	LBFGS_update.block(m_u_dof) -= solution_delta.block(m_u_dof);

	for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
	  {
	    // phasefield active constraints
	    if (m_active_set_phasefield(i) > 0.5)
	      LBFGS_update.block(m_d_dof)[i] = solution_delta_cauchy_point.block(m_d_dof)[i]
					     - solution_delta.block(m_d_dof)[i];
	    else
	      LBFGS_update.block(m_d_dof)[i] = solution_delta_cauchy_point.block(m_d_dof)[i]
					     + search_direction.block(m_d_dof)[i]
					     - solution_delta.block(m_d_dof)[i];
	  }

	// make sure the phasefield solutions are feasible
	for(unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
	  {
	    if (solution_delta.block(m_d_dof)[i] + LBFGS_update.block(m_d_dof)[i] < 0.0)
	      LBFGS_update.block(m_d_dof)[i] = -solution_delta.block(m_d_dof)[i];

	    if (  solution_delta.block(m_d_dof)[i]
		+ m_solution.block(m_d_dof)[i]
		+ LBFGS_update.block(m_d_dof)[i] > 1.0)
	      LBFGS_update.block(m_d_dof)[i] = 1.0 - m_solution.block(m_d_dof)[i]
						   - solution_delta.block(m_d_dof)[i];
	  }

	m_constraints.distribute(LBFGS_update);

	// We need a line search algorithm to decide line_search_parameter

        if(m_parameters.m_type_line_search == "StrongWolfe")
          {
	    const double phi_0 = calculate_energy_functional();
	    const double phi_0_prime = m_system_rhs * LBFGS_update;

	    line_search_parameter = line_search_stepsize_strong_wolfe(phi_0,
								      phi_0_prime,
								      LBFGS_update,
								      solution_delta);
          }
        else if(m_parameters.m_type_line_search == "GradientBased")
          {
	    // LBFGS_r_vector is the search direction
	    line_search_parameter = line_search_stepsize_gradient_based(LBFGS_update,
									solution_delta);
          }
        else
          {
            Assert(false, ExcMessage("An unknown line search method is called!"));
          }

	LBFGS_update *= line_search_parameter;
```

简要说明：对应变量或算子在 main.cc 中有同名或等价离散实现。
