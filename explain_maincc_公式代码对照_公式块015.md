# explain.md 与 main.cc 公式-代码对应（公式块 015）

- 所属章节：`2.1. Phase-field formulation`
- explain.md 行号：`131`

论文公式：

\[
\pmb{\sigma}^{+} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{+}\mathbf{I} + 2\mu \pmb{\epsilon}^{+}, \quad \pmb{\sigma}^{-} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{-}\mathbf{I} + 2\mu \pmb{\epsilon}^{-}.
\]

对应 `main.cc` 代码：

```cpp
// main.cc:829-874
    Vector<double>              eigenvalues(dim);
    std::vector<Tensor<1, dim>> eigenvectors(dim);
    usr_spectrum_decomposition::spectrum_decomposition<dim>(m_strain,
  							      eigenvalues,
  							      eigenvectors);

    SymmetricTensor<2, dim> strain_positive, strain_negative;
    strain_positive = usr_spectrum_decomposition::positive_tensor(eigenvalues, eigenvectors);
    strain_negative = usr_spectrum_decomposition::negative_tensor(eigenvalues, eigenvectors);

    SymmetricTensor<4, dim> projector_positive, projector_negative;
    usr_spectrum_decomposition::positive_negative_projectors(eigenvalues,
  							       eigenvectors,
							       projector_positive,
							       projector_negative);

    SymmetricTensor<2, dim> stress_positive, stress_negative;
    const double degradation = degradation_function(m_phase_field_value) + m_residual_k;
    const double I_1 = trace(m_strain);

    // 2D plane strain and 3D cases
    double my_lambda = m_lame_lambda;

    // 2D plane stress case
    if (    dim == 2
	   && m_plane_stress)
      my_lambda = 2 * m_lame_mu * m_lame_lambda / (m_lame_lambda + 2 * m_lame_mu);

    stress_positive = my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)
                                    * Physics::Elasticity::StandardTensors<dim>::I
                    + 2 * m_lame_mu * strain_positive;
    stress_negative = my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)
                                    * Physics::Elasticity::StandardTensors<dim>::I
    		      + 2 * m_lame_mu * strain_negative;

    m_stress = degradation * stress_positive + stress_negative;
    m_stress_positive = stress_positive;

    SymmetricTensor<4, dim> C_positive, C_negative;
    C_positive = my_lambda * usr_spectrum_decomposition::heaviside_function(I_1)
                               * Physics::Elasticity::StandardTensors<dim>::IxI
		 + 2 * m_lame_mu * projector_positive;
    C_negative = my_lambda * usr_spectrum_decomposition::heaviside_function(-I_1)
                               * Physics::Elasticity::StandardTensors<dim>::IxI
    		 + 2 * m_lame_mu * projector_negative;
    m_mechanical_C = degradation * C_positive + C_negative;
```

简要说明：对应变量或算子在 main.cc 中有同名或等价离散实现。
