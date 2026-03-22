# explain.md 与 main.cc 公式-代码对应（公式块 048）

- 所属章节：`3.1. Algorithm overview`
- explain.md 行号：`370`

论文公式：

\[
\mathbf{x} = \mathbf{C}\hat{\mathbf{x}} + \mathbf{k}.
\]

对应 `main.cc` 代码：

```cpp
// main.cc:3625-3637
    else  // inhomogeneous constraints
      {
        if (m_constraints.has_inhomogeneities())
          {
            AffineConstraints<double> homogeneous_constraints(m_constraints);
            for (unsigned int dof = 0; dof != m_dof_handler.n_dofs(); ++dof)
              if (homogeneous_constraints.is_inhomogeneously_constrained(dof))
                homogeneous_constraints.set_inhomogeneity(dof, 0.0);
            m_constraints.clear();
            m_constraints.copy_from(homogeneous_constraints);
          }
      }
    m_constraints.close();

// main.cc:4512-4523
        make_constraints(LBFGS_iteration);

        // At the first step, we simply distribute the inhomogeneous part of
        // the constraints
        if (LBFGS_iteration == 0)
          {
            // use the solution from the previous solve on the
            // refined mesh as initial guess
            LBFGS_update = LBFGS_update_refine;

            m_constraints.distribute(LBFGS_update);
            solution_delta += LBFGS_update;
```

简要说明：对应变量或算子在 main.cc 中有同名或等价离散实现。
