# explain.md 与 main.cc 公式-代码对应（公式块 023）

- 所属章节：`2.2. Finite element discretization`
- explain.md 行号：`194`

论文公式：

\[
\begin{array}{rl} 
\pmb{r} = \nabla \Pi = (r_{\pmb{u}},r_{d})^{\mathrm{T}},\\ 
r_{\pmb{u}_{A}} = \frac{\partial\Pi}{\partial \pmb{u}_{A}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},\pmb{\sigma}\right) - \left(\pmb{N}_{u_{A}},\pmb{b}\right) - \left(\pmb{N}_{u_{A}},\pmb{t}\right)_{\Gamma_{t}},\\ 
r_{d_{A}} = \frac{\partial\Pi}{\partial d_{A}} = \left(N_{d_{A}},\frac{g_{c}}{l} d + g^{\prime}(d)\psi^{+}\right) + (\nabla N_{d_{A}},g_{c}l\nabla d). 
\end{array} \quad (9)
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
