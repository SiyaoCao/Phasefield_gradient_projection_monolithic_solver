# explain.md 与 main.cc 公式-代码对应（公式块 027）

- 所属章节：`2.2. Finite element discretization`
- explain.md 行号：`225`

论文公式：

\[
\hat{\mathbf{K}} = \left[ \begin{array}{cc}\mathbf{K}_{uu} & \mathbf{0}\\ \mathbf{0} & \mathbf{K}_{dd} \end{array} \right], \quad (11)
\]

对应 `main.cc` 代码：

```cpp
// main.cc:3641-3647
  void PhaseFieldMonolithicSolve<dim>::assemble_system_B0(const BlockVector<double> & solution_old)
  {
    m_timer.enter_subsection("Assemble B0");

    m_tangent_matrix = 0.0;

    const UpdateFlags uf_cell(update_values | update_gradients |

// main.cc:3940-3951
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
```

简要说明：对应有限元离散中的块矩阵装配或约化矩阵表达。
