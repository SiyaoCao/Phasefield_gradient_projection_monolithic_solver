# explain.md 与 main.cc 公式-代码对应（公式块 108）

- 所属章节：`5.2. Comparison of convergence behaviors`
- explain.md 行号：`937`

论文公式：

\[
\| \pmb {r}_u\| _2 < \mathrm{tol}, \quad \| \pmb {r}_d\| _2 < \mathrm{tol}, \quad \| \Delta \pmb {u}\| _2 < \mathrm{tol}, \quad \| \Delta d\| _2 < \mathrm{tol},
\]

对应 `main.cc` 代码：

```cpp
// main.cc:4564-4568
        if (LBFGS_iteration > 1 && m_error_update_norm.m_u <= m_parameters.m_tol_u_incr
                                && m_error_residual_norm.m_u <= m_parameters.m_tol_u_residual
			        && m_error_update_norm.m_d <= m_parameters.m_tol_d_incr
			        && m_error_residual_norm.m_d <= m_parameters.m_tol_d_residual
				)

// main.cc:5033-5039
        if (LBFGS_iteration > 1 && m_error_update_norm.m_u <= m_parameters.m_tol_u_incr
                                && m_error_residual_norm.m_u <= m_parameters.m_tol_u_residual
			        && m_error_update_norm.m_d <= m_parameters.m_tol_d_incr
			        && m_error_residual_norm.m_d <= m_parameters.m_tol_d_residual
				&& lower_bound_number_new == lower_bound_number_old
				&& upper_bound_number_new == upper_bound_number_old
				&& lowerupper_bound_number_new == lowerupper_bound_number_old)
```

简要说明：对应变量或算子在 main.cc 中有同名或等价离散实现。
