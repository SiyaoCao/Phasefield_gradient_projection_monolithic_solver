# explain.md 与 main.cc 公式-代码对应（公式块 045）

- 所属章节：`3.1. Algorithm overview`
- explain.md 行号：`350`

论文公式：

\[
\mathrm{Proj}_C(x_i - r_i,\mathrm{lb}_i,\mathrm{ub}_i) - x_i = \left\{ \begin{array}{ll}\mathrm{lb}_i - x_i & \mathrm{if}\; x_i - r_i < \mathrm{lb}_i,\\ - r_i & \mathrm{if}\; \mathrm{lb}_i \leq x_i - r_i \leq \mathrm{ub}_i,\\ \mathrm{ub}_i - x_i & \mathrm{if}\; x_i - r_i > \mathrm{ub}_i. \end{array} \right. \quad (19)
\]

对应 `main.cc` 代码：

```cpp
// main.cc:350-352
                          "0",
                          Patterns::Integer(0),
                          "Global refinement times (across the entire domain)");

// main.cc:1379-1394
  void PhaseFieldMonolithicSolve<dim>::point_projection(BlockVector<double> & solution_delta)
  {
    // Phase-field value cannot exceed 1.0
    const double upper_limit = 1.0;

    Vector<double> solution_phasefield_total(m_solution.block(m_d_dof));
    solution_phasefield_total += solution_delta.block(m_d_dof);

    for (unsigned int i = 0; i < solution_phasefield_total.size(); ++i)
      {
	if (solution_delta.block(m_d_dof)[i] < 0.0)
	  solution_delta.block(m_d_dof)[i] = 0.0;

	if (solution_phasefield_total[i] > upper_limit)
	  solution_delta.block(m_d_dof)[i] = upper_limit - m_solution.block(m_d_dof)[i];
      }
```

简要说明：对应箱约束投影、主动集识别与可行域修正。
