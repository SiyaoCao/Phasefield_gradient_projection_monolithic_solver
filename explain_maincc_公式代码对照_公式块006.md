# explain.md 与 main.cc 公式-代码对应（公式块 006）

- 所属章节：`1. Introduction`
- explain.md 行号：`63`

论文公式：

\[
\mathrm{lb}_i\leq x_i\leq \mathrm{ub}_i,
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
```

简要说明：对应变量或算子在 main.cc 中有同名或等价离散实现。
