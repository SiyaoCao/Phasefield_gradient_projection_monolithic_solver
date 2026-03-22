# explain.md 与 main.cc 公式-代码对应（公式块 044）

- 所属章节：`3.1. Algorithm overview`
- explain.md 行号：`344`

论文公式：

\[
\| \Delta \mathbf{x}_{k + 1}\|_{2} = \| \mathbf{x}_{k + 1} - \mathbf{x}_{k}\|_{2} < \mathrm{tol}.
\]

对应 `main.cc` 代码：

```cpp
// main.cc:5033-5039
        if (LBFGS_iteration > 1 && m_error_update_norm.m_u <= m_parameters.m_tol_u_incr
                                && m_error_residual_norm.m_u <= m_parameters.m_tol_u_residual
			        && m_error_update_norm.m_d <= m_parameters.m_tol_d_incr
			        && m_error_residual_norm.m_d <= m_parameters.m_tol_d_residual
				&& lower_bound_number_new == lower_bound_number_old
				&& upper_bound_number_new == upper_bound_number_old
				&& lowerupper_bound_number_new == lowerupper_bound_number_old)

// main.cc:5737-5773
	unsigned int number_active_constraint_lower_bound = 0;
	unsigned int number_active_constraint_upper_bound = 0;
	unsigned int number_active_constraint_lowerupper_bound = 0;

	for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
	  {
	    if (   solution_delta.block(m_d_dof)[i] == 0.0
		&& solution_phasefield_total[i] == 1.0)
	      {
		m_active_set_phasefield(i) = 3; //lower bound overlaps with upper bound
		++number_active_constraint_lowerupper_bound;
	      }
	    else if (   solution_delta.block(m_d_dof)[i] == 0.0
		     && solution_phasefield_total[i] != 1.0)
	      {
		m_active_set_phasefield(i) = 1; //lower bound
		++number_active_constraint_lower_bound;
	      }
	    else if (   solution_phasefield_total[i] == 1.0
		     && solution_delta.block(m_d_dof)[i] != 0.0)
	      {
	        m_active_set_phasefield(i) = 2; //upper bound
	        ++number_active_constraint_upper_bound;
	      }
	    else
	      {
	        m_active_set_phasefield(i) = 0;
	      }
	  }

	lower_bound_number_old = lower_bound_number_new;
	upper_bound_number_old = upper_bound_number_new;
	lowerupper_bound_number_old = lowerupper_bound_number_new;

	lower_bound_number_new = number_active_constraint_lower_bound;
	upper_bound_number_new = number_active_constraint_upper_bound;
	lowerupper_bound_number_new = number_active_constraint_lowerupper_bound;
```

简要说明：对应变量或算子在 main.cc 中有同名或等价离散实现。
