# explain.md 与 main.cc 公式-代码对应（公式块 026）

- 所属章节：`2.2. Finite element discretization`
- explain.md 行号：`219`

论文公式：

\[
d_A^{(n)}\leq d_A\leq 1,
\]

对应 `main.cc` 代码：

```cpp
// main.cc:1416-1423
    // upper bound is 1.0, lower bound is the solution at the previous step.
    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
      {
	if (gradient_g.block(m_d_dof)[i] < 0)
	  t = (solution_phasefield_total[i] - 1.0 ) / gradient_g.block(m_d_dof)[i];
	else if (gradient_g.block(m_d_dof)[i] > 0)
	  t = solution_delta.block(m_d_dof)[i] / gradient_g.block(m_d_dof)[i];
	else
```

简要说明：对应变量或算子在 main.cc 中有同名或等价离散实现。
