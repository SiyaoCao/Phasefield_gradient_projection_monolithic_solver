# explain.md 与 main.cc 公式-代码对应（公式块 021）

- 所属章节：`2.2. Finite element discretization`
- explain.md 行号：`178`

论文公式：

\[
\delta \pmb {u} = \pmb{N}_{u_{A}}\delta \pmb{u}_{A}, \quad \text{and} \quad \delta d = N_{d_{A}}\delta d_{A},
\]

对应 `main.cc` 代码：

```cpp
// main.cc:3722-3731
  void PhaseFieldMonolithicSolve<dim>::assemble_system_rhs_BFGS_one_cell(
      const typename DoFHandler<dim>::active_cell_iterator &cell,
      ScratchData_ASM_RHS_BFGS & scratch,
      PerTaskData_ASM_RHS_BFGS & data) const
  {
    data.reset();
    scratch.reset();
    scratch.m_fe_values.reinit(cell);
    cell->get_dof_indices(data.m_local_dof_indices);
```

简要说明：对应变量或算子在 main.cc 中有同名或等价离散实现。
