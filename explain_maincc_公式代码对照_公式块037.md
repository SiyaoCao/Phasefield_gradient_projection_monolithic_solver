# explain.md 与 main.cc 公式-代码对应（公式块 037）

- 所属章节：`3.1. Algorithm overview`
- explain.md 行号：`301`

论文公式：

\[
x_{i} = x_{i}^{c},\;\forall i\in \mathcal{A}_{k}(\pmb{x}^{c})\quad \text{and}\quad \mathrm{lb}_{i}\leq x_{i}\leq \mathrm{ub}_{i},\;\forall i\notin \mathcal{A}_{k}(\pmb{x}^{c}).
\]

对应 `main.cc` 代码：

```cpp
// main.cc:5731-5764
	Vector<double> solution_phasefield_total(m_solution.block(m_d_dof));
	solution_phasefield_total += solution_delta.block(m_d_dof);

	// Since line search parameter might be less than one, we need update
	// the phasefield active set status
	// upper bound is 1.0, lower bound is the solution at the previous step.
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
```

简要说明：对应变量或算子在 main.cc 中有同名或等价离散实现。
