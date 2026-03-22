# explain.md 与 main.cc 公式-代码对应（公式块 019）

- 所属章节：`2.1. Phase-field formulation`
- explain.md 行号：`162`

论文公式：

\[
\left\{ \begin{array}{ll} r_{\pmb{u}}(\pmb{u},d) = (\nabla^{(s)}\delta \pmb{u},\pmb{\sigma}) - (\delta \pmb{u},\pmb{b}) - (\delta \pmb{u},\pmb{t})_{\Gamma_{t}} = 0, \\ r_{d}(\pmb{u},d) = (\delta d,\frac{g_{c}}{l} d) + (\nabla \delta d,g_{c}l\nabla d) + (\delta d,g^{\prime}(d)\psi^{+}(\pmb{\epsilon})) = 0, \end{array} \right.
\]

对应 `main.cc` 代码：

```cpp
// main.cc:3800-3820
        for (const unsigned int i : scratch.m_fe_values.dof_indices())
          {
            const unsigned int i_group = m_fe.system_to_base_index(i).first.first;

            if (i_group == m_u_dof)
              {
                data.m_cell_rhs(i) += (symm_grad_N_disp[i] * cauchy_stress) * JxW;

		// contributions from the body force to right-hand side
		data.m_cell_rhs(i) -= N_disp[i] * rhs_values[q_point] * JxW;
              }
            else if (i_group == m_d_dof)
              {
    	        data.m_cell_rhs(i) += (    gc * length_scale * grad_N_phasefield[i] * phasefield_grad
    	                                +  (   gc / length_scale * phasefield_value
					     + eta / delta_time  * (phasefield_value - old_phasefield)
					     + degradation_function_derivative(phasefield_value)
					     * current_positive_strain_energy )
					  * N_phasefield[i]
				      ) * JxW;
              }
```

简要说明：对应变量或算子在 main.cc 中有同名或等价离散实现。
