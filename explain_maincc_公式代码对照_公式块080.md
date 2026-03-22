# explain.md 与 main.cc 公式-代码对应（公式块 080）

- 所属章节：`3.4. Subspace minimization`
- explain.md 行号：`575`

论文公式：

\[
x_i = x^c_i, \quad \forall i \in \mathcal{A}(\mathbf{x}^c)
\]

对应 `main.cc` 代码：

```cpp
// main.cc:5193-5220
	// We need to find out which DOFs are free:
	// no essential boundary conditions, no hanging node constraints
	// no active box constraints
	unsigned int free_disp_number = 0;
	for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
	  {
	    if (m_constraints.is_constrained(i))
	      free_dofs.block(m_u_dof)[i] = -1;
	    else
	      {
	        free_dofs.block(m_u_dof)[i] = 1;
	        ++free_disp_number;
	      }
	  }

	unsigned int free_phasefield_number = 0;
	for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
	  {
	    if (   m_constraints.is_constrained(i + m_dofs_per_block[m_u_dof])
		|| m_active_set_phasefield(i) > 0.5)
	      free_dofs.block(m_d_dof)[i] = -1;
	    else
	      {
	        free_dofs.block(m_d_dof)[i] = 1;
	        ++free_phasefield_number;
	      }
	  }


// main.cc:5638-5647
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
```

简要说明：对应变量或算子在 main.cc 中有同名或等价离散实现。
