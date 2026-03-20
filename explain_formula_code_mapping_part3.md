## 012. ### 2.1. Phase-field formulation（explain.md:L113-L115）

- 对应关系说明：本公式位于本构层：应变谱分解 + 正负能量分裂 + 退化函数。
- 最底层代码链接：
  - [main.cc:L831-L874](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L831-L874)
  - [SpectrumDecomposition.h:L71-L143](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/SpectrumDecomposition.h#L71-L143)
  - [SpectrumDecomposition.cc:L19-L38](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/SpectrumDecomposition.cc#L19-L38)

### 代码片段 1（main.cc:L831-L874）
```cpp
 831:     usr_spectrum_decomposition::spectrum_decomposition<dim>(m_strain,
 832:   							      eigenvalues,
 833:   							      eigenvectors);
 834: 
 835:     SymmetricTensor<2, dim> strain_positive, strain_negative;
 836:     strain_positive = usr_spectrum_decomposition::positive_tensor(eigenvalues, eigenvectors);
 837:     strain_negative = usr_spectrum_decomposition::negative_tensor(eigenvalues, eigenvectors);
 838: 
 839:     SymmetricTensor<4, dim> projector_positive, projector_negative;
 840:     usr_spectrum_decomposition::positive_negative_projectors(eigenvalues,
 841:   							       eigenvectors,
 842: 							       projector_positive,
 843: 							       projector_negative);
 844: 
 845:     SymmetricTensor<2, dim> stress_positive, stress_negative;
 846:     const double degradation = degradation_function(m_phase_field_value) + m_residual_k;
 847:     const double I_1 = trace(m_strain);
 848: 
 849:     // 2D plane strain and 3D cases
 850:     double my_lambda = m_lame_lambda;
 851: 
 852:     // 2D plane stress case
 853:     if (    dim == 2
 854: 	   && m_plane_stress)
 855:       my_lambda = 2 * m_lame_mu * m_lame_lambda / (m_lame_lambda + 2 * m_lame_mu);
 856: 
 857:     stress_positive = my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)
 858:                                     * Physics::Elasticity::StandardTensors<dim>::I
 859:                     + 2 * m_lame_mu * strain_positive;
 860:     stress_negative = my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)
 861:                                     * Physics::Elasticity::StandardTensors<dim>::I
 862:     		      + 2 * m_lame_mu * strain_negative;
 863: 
 864:     m_stress = degradation * stress_positive + stress_negative;
 865:     m_stress_positive = stress_positive;
 866: 
 867:     SymmetricTensor<4, dim> C_positive, C_negative;
 868:     C_positive = my_lambda * usr_spectrum_decomposition::heaviside_function(I_1)
 869:                                * Physics::Elasticity::StandardTensors<dim>::IxI
 870: 		 + 2 * m_lame_mu * projector_positive;
 871:     C_negative = my_lambda * usr_spectrum_decomposition::heaviside_function(-I_1)
 872:                                * Physics::Elasticity::StandardTensors<dim>::IxI
 873:     		 + 2 * m_lame_mu * projector_negative;
 874:     m_mechanical_C = degradation * C_positive + C_negative;
```

### 代码片段 2（SpectrumDecomposition.h:L71-L143）
```cpp
  71:   template <int dim>
  72:   void positive_negative_projectors(Vector<double> const & eigenvalues,
  73:                                     std::vector<Tensor<1, dim>> const & eigenvectors,
  74: 			            SymmetricTensor<4, dim> & positive_projector,
  75: 				    SymmetricTensor<4, dim> & negative_projector)
  76:   {
  77:     Assert(dim <= 3,
  78: 	   ExcMessage("Project tensors only work for dim <= 3."));
  79: 
  80:     std::array<SymmetricTensor<2, dim>, dim> M;
  81:     for (int a = 0; a < dim; a++)
  82:       M[a] = symmetrize(outer_product(eigenvectors[a], eigenvectors[a]));
  83: 
  84:     std::array<SymmetricTensor<4, dim>, dim> Q;
  85:     for (int a = 0; a < dim; a++)
  86:       Q[a] = outer_product(M[a], M[a]);
  87: 
  88:     std::array<std::array<SymmetricTensor<4, dim>, dim>, dim> G;
  89:     for (int a = 0; a < dim; a++)
  90:       for (int b = 0; b < dim; b++)
  91: 	for (int i = 0; i < dim; i++)
  92: 	  for (int j = 0; j < dim; j++)
  93: 	    for (int k = 0; k < dim; k++)
  94:               for (int l = 0; l < dim; l++)
  95:         	G[a][b][i][j][k][l] = M[a][i][k] * M[b][j][l]
  96: 				    + M[a][i][l] * M[b][j][k];
  97: 
  98:     positive_projector = 0;
  99:     for (int a = 0; a < dim; a++)
 100:       {
 101: 	double lambda_a = eigenvalues[a];
 102: 	positive_projector += heaviside_function(lambda_a)
 103: 			    * Q[a];
 104: 	for (int b = 0; b < dim; b++)
 105: 	  {
 106: 	    if (b != a)
 107: 	      {
 108: 		double lambda_b = eigenvalues[b];
 109: 		double v_ab = 0.0;
 110: 		if (std::fabs(lambda_a - lambda_b) > 1.0e-12)
 111: 		  v_ab = (positive_ramp_function(lambda_a) - positive_ramp_function(lambda_b))
 112: 		       / (lambda_a - lambda_b);
 113: 		else
 114: 		  v_ab = 0.5 * (  heaviside_function(lambda_a)
 115: 		                + heaviside_function(lambda_b) );
 116: 		positive_projector += 0.5 * v_ab * 0.5 * (G[a][b] + G[b][a]);
 117: 	      }
 118: 	  }
 119:       }
 120: 
 121:     negative_projector = 0;
 122:     for (int a = 0; a < dim; a++)
 123:       {
 124: 	double lambda_a = eigenvalues[a];
 125: 	negative_projector += heaviside_function(-lambda_a)
 126: 			    * Q[a];
 127: 	for (int b = 0; b < dim; b++)
 128: 	  {
 129: 	    if (b != a)
 130: 	      {
 131: 		double lambda_b = eigenvalues[b];
 132: 		double v_ab = 0.0;
 133: 		if (std::fabs(lambda_a - lambda_b) > 1.0e-12)
 134: 		  v_ab = (negative_ramp_function(lambda_a) - negative_ramp_function(lambda_b))
 135: 		       / (lambda_a - lambda_b);
 136: 		else
 137: 		  v_ab = 0.5 * (  heaviside_function(-lambda_a)
 138: 		                + heaviside_function(-lambda_b) );
 139: 		negative_projector += 0.5 * v_ab * 0.5 * (G[a][b] + G[b][a]);
 140: 	      }
 141: 	  }
 142:       }
 143: 
```

### 代码片段 3（SpectrumDecomposition.cc:L19-L38）
```cpp
  19:   double positive_ramp_function(const double x)
  20:   {
  21:     return std::fmax(x, 0.0);
  22:   }
  23: 
  24:   double negative_ramp_function(const double x)
  25:   {
  26:     return std::fmin(x, 0.0);
  27:   }
  28: 
  29:   double heaviside_function(const double x)
  30:   {
  31:     if (std::fabs(x) < 1.0e-16)
  32:       return 0.5;
  33: 
  34:     if (x > 0)
  35:       return 1.0;
  36:     else
  37:       return 0.0;
  38:   }
```

### 公式（与 explain.md 一致）
\[
\pmb{\epsilon}^{+} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{+} \mathbf{M}_{\alpha}, \quad \pmb{\epsilon}^{-} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{-} \mathbf{M}_{\alpha}.
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 013. ### 2.1. Phase-field formulation（explain.md:L119-L121）

- 对应关系说明：本公式位于本构层：应变谱分解 + 正负能量分裂 + 退化函数。
- 最底层代码链接：
  - [main.cc:L831-L874](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L831-L874)
  - [SpectrumDecomposition.h:L71-L143](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/SpectrumDecomposition.h#L71-L143)
  - [SpectrumDecomposition.cc:L19-L38](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/SpectrumDecomposition.cc#L19-L38)

### 代码片段 1（main.cc:L831-L874）
```cpp
 831:     usr_spectrum_decomposition::spectrum_decomposition<dim>(m_strain,
 832:   							      eigenvalues,
 833:   							      eigenvectors);
 834: 
 835:     SymmetricTensor<2, dim> strain_positive, strain_negative;
 836:     strain_positive = usr_spectrum_decomposition::positive_tensor(eigenvalues, eigenvectors);
 837:     strain_negative = usr_spectrum_decomposition::negative_tensor(eigenvalues, eigenvectors);
 838: 
 839:     SymmetricTensor<4, dim> projector_positive, projector_negative;
 840:     usr_spectrum_decomposition::positive_negative_projectors(eigenvalues,
 841:   							       eigenvectors,
 842: 							       projector_positive,
 843: 							       projector_negative);
 844: 
 845:     SymmetricTensor<2, dim> stress_positive, stress_negative;
 846:     const double degradation = degradation_function(m_phase_field_value) + m_residual_k;
 847:     const double I_1 = trace(m_strain);
 848: 
 849:     // 2D plane strain and 3D cases
 850:     double my_lambda = m_lame_lambda;
 851: 
 852:     // 2D plane stress case
 853:     if (    dim == 2
 854: 	   && m_plane_stress)
 855:       my_lambda = 2 * m_lame_mu * m_lame_lambda / (m_lame_lambda + 2 * m_lame_mu);
 856: 
 857:     stress_positive = my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)
 858:                                     * Physics::Elasticity::StandardTensors<dim>::I
 859:                     + 2 * m_lame_mu * strain_positive;
 860:     stress_negative = my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)
 861:                                     * Physics::Elasticity::StandardTensors<dim>::I
 862:     		      + 2 * m_lame_mu * strain_negative;
 863: 
 864:     m_stress = degradation * stress_positive + stress_negative;
 865:     m_stress_positive = stress_positive;
 866: 
 867:     SymmetricTensor<4, dim> C_positive, C_negative;
 868:     C_positive = my_lambda * usr_spectrum_decomposition::heaviside_function(I_1)
 869:                                * Physics::Elasticity::StandardTensors<dim>::IxI
 870: 		 + 2 * m_lame_mu * projector_positive;
 871:     C_negative = my_lambda * usr_spectrum_decomposition::heaviside_function(-I_1)
 872:                                * Physics::Elasticity::StandardTensors<dim>::IxI
 873:     		 + 2 * m_lame_mu * projector_negative;
 874:     m_mechanical_C = degradation * C_positive + C_negative;
```

### 代码片段 2（SpectrumDecomposition.h:L71-L143）
```cpp
  71:   template <int dim>
  72:   void positive_negative_projectors(Vector<double> const & eigenvalues,
  73:                                     std::vector<Tensor<1, dim>> const & eigenvectors,
  74: 			            SymmetricTensor<4, dim> & positive_projector,
  75: 				    SymmetricTensor<4, dim> & negative_projector)
  76:   {
  77:     Assert(dim <= 3,
  78: 	   ExcMessage("Project tensors only work for dim <= 3."));
  79: 
  80:     std::array<SymmetricTensor<2, dim>, dim> M;
  81:     for (int a = 0; a < dim; a++)
  82:       M[a] = symmetrize(outer_product(eigenvectors[a], eigenvectors[a]));
  83: 
  84:     std::array<SymmetricTensor<4, dim>, dim> Q;
  85:     for (int a = 0; a < dim; a++)
  86:       Q[a] = outer_product(M[a], M[a]);
  87: 
  88:     std::array<std::array<SymmetricTensor<4, dim>, dim>, dim> G;
  89:     for (int a = 0; a < dim; a++)
  90:       for (int b = 0; b < dim; b++)
  91: 	for (int i = 0; i < dim; i++)
  92: 	  for (int j = 0; j < dim; j++)
  93: 	    for (int k = 0; k < dim; k++)
  94:               for (int l = 0; l < dim; l++)
  95:         	G[a][b][i][j][k][l] = M[a][i][k] * M[b][j][l]
  96: 				    + M[a][i][l] * M[b][j][k];
  97: 
  98:     positive_projector = 0;
  99:     for (int a = 0; a < dim; a++)
 100:       {
 101: 	double lambda_a = eigenvalues[a];
 102: 	positive_projector += heaviside_function(lambda_a)
 103: 			    * Q[a];
 104: 	for (int b = 0; b < dim; b++)
 105: 	  {
 106: 	    if (b != a)
 107: 	      {
 108: 		double lambda_b = eigenvalues[b];
 109: 		double v_ab = 0.0;
 110: 		if (std::fabs(lambda_a - lambda_b) > 1.0e-12)
 111: 		  v_ab = (positive_ramp_function(lambda_a) - positive_ramp_function(lambda_b))
 112: 		       / (lambda_a - lambda_b);
 113: 		else
 114: 		  v_ab = 0.5 * (  heaviside_function(lambda_a)
 115: 		                + heaviside_function(lambda_b) );
 116: 		positive_projector += 0.5 * v_ab * 0.5 * (G[a][b] + G[b][a]);
 117: 	      }
 118: 	  }
 119:       }
 120: 
 121:     negative_projector = 0;
 122:     for (int a = 0; a < dim; a++)
 123:       {
 124: 	double lambda_a = eigenvalues[a];
 125: 	negative_projector += heaviside_function(-lambda_a)
 126: 			    * Q[a];
 127: 	for (int b = 0; b < dim; b++)
 128: 	  {
 129: 	    if (b != a)
 130: 	      {
 131: 		double lambda_b = eigenvalues[b];
 132: 		double v_ab = 0.0;
 133: 		if (std::fabs(lambda_a - lambda_b) > 1.0e-12)
 134: 		  v_ab = (negative_ramp_function(lambda_a) - negative_ramp_function(lambda_b))
 135: 		       / (lambda_a - lambda_b);
 136: 		else
 137: 		  v_ab = 0.5 * (  heaviside_function(-lambda_a)
 138: 		                + heaviside_function(-lambda_b) );
 139: 		negative_projector += 0.5 * v_ab * 0.5 * (G[a][b] + G[b][a]);
 140: 	      }
 141: 	  }
 142:       }
 143: 
```

### 代码片段 3（SpectrumDecomposition.cc:L19-L38）
```cpp
  19:   double positive_ramp_function(const double x)
  20:   {
  21:     return std::fmax(x, 0.0);
  22:   }
  23: 
  24:   double negative_ramp_function(const double x)
  25:   {
  26:     return std::fmin(x, 0.0);
  27:   }
  28: 
  29:   double heaviside_function(const double x)
  30:   {
  31:     if (std::fabs(x) < 1.0e-16)
  32:       return 0.5;
  33: 
  34:     if (x > 0)
  35:       return 1.0;
  36:     else
  37:       return 0.0;
  38:   }
```

### 公式（与 explain.md 一致）
\[
\psi^{+}(\pmb{\epsilon}) = \frac{1}{2}\lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{+}^{2} + \mu \pmb{\epsilon}^{+} : \pmb{\epsilon}^{+}, \quad \psi^{-}(\pmb{\epsilon}) = \frac{1}{2} \lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{-}^{2} + \mu \pmb{\epsilon}^{-} : \pmb{\epsilon}^{-},
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 014. ### 2.1. Phase-field formulation（explain.md:L125-L127）

- 对应关系说明：本公式位于本构层：应变谱分解 + 正负能量分裂 + 退化函数。
- 最底层代码链接：
  - [main.cc:L831-L874](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L831-L874)
  - [SpectrumDecomposition.h:L71-L143](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/SpectrumDecomposition.h#L71-L143)
  - [SpectrumDecomposition.cc:L19-L38](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/SpectrumDecomposition.cc#L19-L38)

### 代码片段 1（main.cc:L831-L874）
```cpp
 831:     usr_spectrum_decomposition::spectrum_decomposition<dim>(m_strain,
 832:   							      eigenvalues,
 833:   							      eigenvectors);
 834: 
 835:     SymmetricTensor<2, dim> strain_positive, strain_negative;
 836:     strain_positive = usr_spectrum_decomposition::positive_tensor(eigenvalues, eigenvectors);
 837:     strain_negative = usr_spectrum_decomposition::negative_tensor(eigenvalues, eigenvectors);
 838: 
 839:     SymmetricTensor<4, dim> projector_positive, projector_negative;
 840:     usr_spectrum_decomposition::positive_negative_projectors(eigenvalues,
 841:   							       eigenvectors,
 842: 							       projector_positive,
 843: 							       projector_negative);
 844: 
 845:     SymmetricTensor<2, dim> stress_positive, stress_negative;
 846:     const double degradation = degradation_function(m_phase_field_value) + m_residual_k;
 847:     const double I_1 = trace(m_strain);
 848: 
 849:     // 2D plane strain and 3D cases
 850:     double my_lambda = m_lame_lambda;
 851: 
 852:     // 2D plane stress case
 853:     if (    dim == 2
 854: 	   && m_plane_stress)
 855:       my_lambda = 2 * m_lame_mu * m_lame_lambda / (m_lame_lambda + 2 * m_lame_mu);
 856: 
 857:     stress_positive = my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)
 858:                                     * Physics::Elasticity::StandardTensors<dim>::I
 859:                     + 2 * m_lame_mu * strain_positive;
 860:     stress_negative = my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)
 861:                                     * Physics::Elasticity::StandardTensors<dim>::I
 862:     		      + 2 * m_lame_mu * strain_negative;
 863: 
 864:     m_stress = degradation * stress_positive + stress_negative;
 865:     m_stress_positive = stress_positive;
 866: 
 867:     SymmetricTensor<4, dim> C_positive, C_negative;
 868:     C_positive = my_lambda * usr_spectrum_decomposition::heaviside_function(I_1)
 869:                                * Physics::Elasticity::StandardTensors<dim>::IxI
 870: 		 + 2 * m_lame_mu * projector_positive;
 871:     C_negative = my_lambda * usr_spectrum_decomposition::heaviside_function(-I_1)
 872:                                * Physics::Elasticity::StandardTensors<dim>::IxI
 873:     		 + 2 * m_lame_mu * projector_negative;
 874:     m_mechanical_C = degradation * C_positive + C_negative;
```

### 代码片段 2（SpectrumDecomposition.h:L71-L143）
```cpp
  71:   template <int dim>
  72:   void positive_negative_projectors(Vector<double> const & eigenvalues,
  73:                                     std::vector<Tensor<1, dim>> const & eigenvectors,
  74: 			            SymmetricTensor<4, dim> & positive_projector,
  75: 				    SymmetricTensor<4, dim> & negative_projector)
  76:   {
  77:     Assert(dim <= 3,
  78: 	   ExcMessage("Project tensors only work for dim <= 3."));
  79: 
  80:     std::array<SymmetricTensor<2, dim>, dim> M;
  81:     for (int a = 0; a < dim; a++)
  82:       M[a] = symmetrize(outer_product(eigenvectors[a], eigenvectors[a]));
  83: 
  84:     std::array<SymmetricTensor<4, dim>, dim> Q;
  85:     for (int a = 0; a < dim; a++)
  86:       Q[a] = outer_product(M[a], M[a]);
  87: 
  88:     std::array<std::array<SymmetricTensor<4, dim>, dim>, dim> G;
  89:     for (int a = 0; a < dim; a++)
  90:       for (int b = 0; b < dim; b++)
  91: 	for (int i = 0; i < dim; i++)
  92: 	  for (int j = 0; j < dim; j++)
  93: 	    for (int k = 0; k < dim; k++)
  94:               for (int l = 0; l < dim; l++)
  95:         	G[a][b][i][j][k][l] = M[a][i][k] * M[b][j][l]
  96: 				    + M[a][i][l] * M[b][j][k];
  97: 
  98:     positive_projector = 0;
  99:     for (int a = 0; a < dim; a++)
 100:       {
 101: 	double lambda_a = eigenvalues[a];
 102: 	positive_projector += heaviside_function(lambda_a)
 103: 			    * Q[a];
 104: 	for (int b = 0; b < dim; b++)
 105: 	  {
 106: 	    if (b != a)
 107: 	      {
 108: 		double lambda_b = eigenvalues[b];
 109: 		double v_ab = 0.0;
 110: 		if (std::fabs(lambda_a - lambda_b) > 1.0e-12)
 111: 		  v_ab = (positive_ramp_function(lambda_a) - positive_ramp_function(lambda_b))
 112: 		       / (lambda_a - lambda_b);
 113: 		else
 114: 		  v_ab = 0.5 * (  heaviside_function(lambda_a)
 115: 		                + heaviside_function(lambda_b) );
 116: 		positive_projector += 0.5 * v_ab * 0.5 * (G[a][b] + G[b][a]);
 117: 	      }
 118: 	  }
 119:       }
 120: 
 121:     negative_projector = 0;
 122:     for (int a = 0; a < dim; a++)
 123:       {
 124: 	double lambda_a = eigenvalues[a];
 125: 	negative_projector += heaviside_function(-lambda_a)
 126: 			    * Q[a];
 127: 	for (int b = 0; b < dim; b++)
 128: 	  {
 129: 	    if (b != a)
 130: 	      {
 131: 		double lambda_b = eigenvalues[b];
 132: 		double v_ab = 0.0;
 133: 		if (std::fabs(lambda_a - lambda_b) > 1.0e-12)
 134: 		  v_ab = (negative_ramp_function(lambda_a) - negative_ramp_function(lambda_b))
 135: 		       / (lambda_a - lambda_b);
 136: 		else
 137: 		  v_ab = 0.5 * (  heaviside_function(-lambda_a)
 138: 		                + heaviside_function(-lambda_b) );
 139: 		negative_projector += 0.5 * v_ab * 0.5 * (G[a][b] + G[b][a]);
 140: 	      }
 141: 	  }
 142:       }
 143: 
```

### 代码片段 3（SpectrumDecomposition.cc:L19-L38）
```cpp
  19:   double positive_ramp_function(const double x)
  20:   {
  21:     return std::fmax(x, 0.0);
  22:   }
  23: 
  24:   double negative_ramp_function(const double x)
  25:   {
  26:     return std::fmin(x, 0.0);
  27:   }
  28: 
  29:   double heaviside_function(const double x)
  30:   {
  31:     if (std::fabs(x) < 1.0e-16)
  32:       return 0.5;
  33: 
  34:     if (x > 0)
  35:       return 1.0;
  36:     else
  37:       return 0.0;
  38:   }
```

### 公式（与 explain.md 一致）
\[
\pmb{\sigma} = \frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial \pmb{\epsilon}} = [g(d) + k]\frac{\partial\psi^{+}(\pmb{\epsilon})}{\partial \pmb{\epsilon}} + \frac{\partial\psi^{-}(\pmb{\epsilon})}{\partial \pmb{\epsilon}} = [g(d) + k]\pmb{\sigma}^{+} + \pmb{\sigma}^{-},
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 015. ### 2.1. Phase-field formulation（explain.md:L131-L133）

- 对应关系说明：本公式位于本构层：应变谱分解 + 正负能量分裂 + 退化函数。
- 最底层代码链接：
  - [main.cc:L831-L874](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L831-L874)
  - [SpectrumDecomposition.h:L71-L143](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/SpectrumDecomposition.h#L71-L143)
  - [SpectrumDecomposition.cc:L19-L38](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/SpectrumDecomposition.cc#L19-L38)

### 代码片段 1（main.cc:L831-L874）
```cpp
 831:     usr_spectrum_decomposition::spectrum_decomposition<dim>(m_strain,
 832:   							      eigenvalues,
 833:   							      eigenvectors);
 834: 
 835:     SymmetricTensor<2, dim> strain_positive, strain_negative;
 836:     strain_positive = usr_spectrum_decomposition::positive_tensor(eigenvalues, eigenvectors);
 837:     strain_negative = usr_spectrum_decomposition::negative_tensor(eigenvalues, eigenvectors);
 838: 
 839:     SymmetricTensor<4, dim> projector_positive, projector_negative;
 840:     usr_spectrum_decomposition::positive_negative_projectors(eigenvalues,
 841:   							       eigenvectors,
 842: 							       projector_positive,
 843: 							       projector_negative);
 844: 
 845:     SymmetricTensor<2, dim> stress_positive, stress_negative;
 846:     const double degradation = degradation_function(m_phase_field_value) + m_residual_k;
 847:     const double I_1 = trace(m_strain);
 848: 
 849:     // 2D plane strain and 3D cases
 850:     double my_lambda = m_lame_lambda;
 851: 
 852:     // 2D plane stress case
 853:     if (    dim == 2
 854: 	   && m_plane_stress)
 855:       my_lambda = 2 * m_lame_mu * m_lame_lambda / (m_lame_lambda + 2 * m_lame_mu);
 856: 
 857:     stress_positive = my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)
 858:                                     * Physics::Elasticity::StandardTensors<dim>::I
 859:                     + 2 * m_lame_mu * strain_positive;
 860:     stress_negative = my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)
 861:                                     * Physics::Elasticity::StandardTensors<dim>::I
 862:     		      + 2 * m_lame_mu * strain_negative;
 863: 
 864:     m_stress = degradation * stress_positive + stress_negative;
 865:     m_stress_positive = stress_positive;
 866: 
 867:     SymmetricTensor<4, dim> C_positive, C_negative;
 868:     C_positive = my_lambda * usr_spectrum_decomposition::heaviside_function(I_1)
 869:                                * Physics::Elasticity::StandardTensors<dim>::IxI
 870: 		 + 2 * m_lame_mu * projector_positive;
 871:     C_negative = my_lambda * usr_spectrum_decomposition::heaviside_function(-I_1)
 872:                                * Physics::Elasticity::StandardTensors<dim>::IxI
 873:     		 + 2 * m_lame_mu * projector_negative;
 874:     m_mechanical_C = degradation * C_positive + C_negative;
```

### 代码片段 2（SpectrumDecomposition.h:L71-L143）
```cpp
  71:   template <int dim>
  72:   void positive_negative_projectors(Vector<double> const & eigenvalues,
  73:                                     std::vector<Tensor<1, dim>> const & eigenvectors,
  74: 			            SymmetricTensor<4, dim> & positive_projector,
  75: 				    SymmetricTensor<4, dim> & negative_projector)
  76:   {
  77:     Assert(dim <= 3,
  78: 	   ExcMessage("Project tensors only work for dim <= 3."));
  79: 
  80:     std::array<SymmetricTensor<2, dim>, dim> M;
  81:     for (int a = 0; a < dim; a++)
  82:       M[a] = symmetrize(outer_product(eigenvectors[a], eigenvectors[a]));
  83: 
  84:     std::array<SymmetricTensor<4, dim>, dim> Q;
  85:     for (int a = 0; a < dim; a++)
  86:       Q[a] = outer_product(M[a], M[a]);
  87: 
  88:     std::array<std::array<SymmetricTensor<4, dim>, dim>, dim> G;
  89:     for (int a = 0; a < dim; a++)
  90:       for (int b = 0; b < dim; b++)
  91: 	for (int i = 0; i < dim; i++)
  92: 	  for (int j = 0; j < dim; j++)
  93: 	    for (int k = 0; k < dim; k++)
  94:               for (int l = 0; l < dim; l++)
  95:         	G[a][b][i][j][k][l] = M[a][i][k] * M[b][j][l]
  96: 				    + M[a][i][l] * M[b][j][k];
  97: 
  98:     positive_projector = 0;
  99:     for (int a = 0; a < dim; a++)
 100:       {
 101: 	double lambda_a = eigenvalues[a];
 102: 	positive_projector += heaviside_function(lambda_a)
 103: 			    * Q[a];
 104: 	for (int b = 0; b < dim; b++)
 105: 	  {
 106: 	    if (b != a)
 107: 	      {
 108: 		double lambda_b = eigenvalues[b];
 109: 		double v_ab = 0.0;
 110: 		if (std::fabs(lambda_a - lambda_b) > 1.0e-12)
 111: 		  v_ab = (positive_ramp_function(lambda_a) - positive_ramp_function(lambda_b))
 112: 		       / (lambda_a - lambda_b);
 113: 		else
 114: 		  v_ab = 0.5 * (  heaviside_function(lambda_a)
 115: 		                + heaviside_function(lambda_b) );
 116: 		positive_projector += 0.5 * v_ab * 0.5 * (G[a][b] + G[b][a]);
 117: 	      }
 118: 	  }
 119:       }
 120: 
 121:     negative_projector = 0;
 122:     for (int a = 0; a < dim; a++)
 123:       {
 124: 	double lambda_a = eigenvalues[a];
 125: 	negative_projector += heaviside_function(-lambda_a)
 126: 			    * Q[a];
 127: 	for (int b = 0; b < dim; b++)
 128: 	  {
 129: 	    if (b != a)
 130: 	      {
 131: 		double lambda_b = eigenvalues[b];
 132: 		double v_ab = 0.0;
 133: 		if (std::fabs(lambda_a - lambda_b) > 1.0e-12)
 134: 		  v_ab = (negative_ramp_function(lambda_a) - negative_ramp_function(lambda_b))
 135: 		       / (lambda_a - lambda_b);
 136: 		else
 137: 		  v_ab = 0.5 * (  heaviside_function(-lambda_a)
 138: 		                + heaviside_function(-lambda_b) );
 139: 		negative_projector += 0.5 * v_ab * 0.5 * (G[a][b] + G[b][a]);
 140: 	      }
 141: 	  }
 142:       }
 143: 
```

### 代码片段 3（SpectrumDecomposition.cc:L19-L38）
```cpp
  19:   double positive_ramp_function(const double x)
  20:   {
  21:     return std::fmax(x, 0.0);
  22:   }
  23: 
  24:   double negative_ramp_function(const double x)
  25:   {
  26:     return std::fmin(x, 0.0);
  27:   }
  28: 
  29:   double heaviside_function(const double x)
  30:   {
  31:     if (std::fabs(x) < 1.0e-16)
  32:       return 0.5;
  33: 
  34:     if (x > 0)
  35:       return 1.0;
  36:     else
  37:       return 0.0;
  38:   }
```

### 公式（与 explain.md 一致）
\[
\pmb{\sigma}^{+} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{+}\mathbf{I} + 2\mu \pmb{\epsilon}^{+}, \quad \pmb{\sigma}^{-} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{-}\mathbf{I} + 2\mu \pmb{\epsilon}^{-}.
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 016. ### 2.1. Phase-field formulation（explain.md:L137-L139）

- 对应关系说明：本公式位于本构层：应变谱分解 + 正负能量分裂 + 退化函数。
- 最底层代码链接：
  - [main.cc:L831-L874](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L831-L874)
  - [SpectrumDecomposition.h:L71-L143](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/SpectrumDecomposition.h#L71-L143)
  - [SpectrumDecomposition.cc:L19-L38](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/SpectrumDecomposition.cc#L19-L38)

### 代码片段 1（main.cc:L831-L874）
```cpp
 831:     usr_spectrum_decomposition::spectrum_decomposition<dim>(m_strain,
 832:   							      eigenvalues,
 833:   							      eigenvectors);
 834: 
 835:     SymmetricTensor<2, dim> strain_positive, strain_negative;
 836:     strain_positive = usr_spectrum_decomposition::positive_tensor(eigenvalues, eigenvectors);
 837:     strain_negative = usr_spectrum_decomposition::negative_tensor(eigenvalues, eigenvectors);
 838: 
 839:     SymmetricTensor<4, dim> projector_positive, projector_negative;
 840:     usr_spectrum_decomposition::positive_negative_projectors(eigenvalues,
 841:   							       eigenvectors,
 842: 							       projector_positive,
 843: 							       projector_negative);
 844: 
 845:     SymmetricTensor<2, dim> stress_positive, stress_negative;
 846:     const double degradation = degradation_function(m_phase_field_value) + m_residual_k;
 847:     const double I_1 = trace(m_strain);
 848: 
 849:     // 2D plane strain and 3D cases
 850:     double my_lambda = m_lame_lambda;
 851: 
 852:     // 2D plane stress case
 853:     if (    dim == 2
 854: 	   && m_plane_stress)
 855:       my_lambda = 2 * m_lame_mu * m_lame_lambda / (m_lame_lambda + 2 * m_lame_mu);
 856: 
 857:     stress_positive = my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)
 858:                                     * Physics::Elasticity::StandardTensors<dim>::I
 859:                     + 2 * m_lame_mu * strain_positive;
 860:     stress_negative = my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)
 861:                                     * Physics::Elasticity::StandardTensors<dim>::I
 862:     		      + 2 * m_lame_mu * strain_negative;
 863: 
 864:     m_stress = degradation * stress_positive + stress_negative;
 865:     m_stress_positive = stress_positive;
 866: 
 867:     SymmetricTensor<4, dim> C_positive, C_negative;
 868:     C_positive = my_lambda * usr_spectrum_decomposition::heaviside_function(I_1)
 869:                                * Physics::Elasticity::StandardTensors<dim>::IxI
 870: 		 + 2 * m_lame_mu * projector_positive;
 871:     C_negative = my_lambda * usr_spectrum_decomposition::heaviside_function(-I_1)
 872:                                * Physics::Elasticity::StandardTensors<dim>::IxI
 873:     		 + 2 * m_lame_mu * projector_negative;
 874:     m_mechanical_C = degradation * C_positive + C_negative;
```

### 代码片段 2（SpectrumDecomposition.h:L71-L143）
```cpp
  71:   template <int dim>
  72:   void positive_negative_projectors(Vector<double> const & eigenvalues,
  73:                                     std::vector<Tensor<1, dim>> const & eigenvectors,
  74: 			            SymmetricTensor<4, dim> & positive_projector,
  75: 				    SymmetricTensor<4, dim> & negative_projector)
  76:   {
  77:     Assert(dim <= 3,
  78: 	   ExcMessage("Project tensors only work for dim <= 3."));
  79: 
  80:     std::array<SymmetricTensor<2, dim>, dim> M;
  81:     for (int a = 0; a < dim; a++)
  82:       M[a] = symmetrize(outer_product(eigenvectors[a], eigenvectors[a]));
  83: 
  84:     std::array<SymmetricTensor<4, dim>, dim> Q;
  85:     for (int a = 0; a < dim; a++)
  86:       Q[a] = outer_product(M[a], M[a]);
  87: 
  88:     std::array<std::array<SymmetricTensor<4, dim>, dim>, dim> G;
  89:     for (int a = 0; a < dim; a++)
  90:       for (int b = 0; b < dim; b++)
  91: 	for (int i = 0; i < dim; i++)
  92: 	  for (int j = 0; j < dim; j++)
  93: 	    for (int k = 0; k < dim; k++)
  94:               for (int l = 0; l < dim; l++)
  95:         	G[a][b][i][j][k][l] = M[a][i][k] * M[b][j][l]
  96: 				    + M[a][i][l] * M[b][j][k];
  97: 
  98:     positive_projector = 0;
  99:     for (int a = 0; a < dim; a++)
 100:       {
 101: 	double lambda_a = eigenvalues[a];
 102: 	positive_projector += heaviside_function(lambda_a)
 103: 			    * Q[a];
 104: 	for (int b = 0; b < dim; b++)
 105: 	  {
 106: 	    if (b != a)
 107: 	      {
 108: 		double lambda_b = eigenvalues[b];
 109: 		double v_ab = 0.0;
 110: 		if (std::fabs(lambda_a - lambda_b) > 1.0e-12)
 111: 		  v_ab = (positive_ramp_function(lambda_a) - positive_ramp_function(lambda_b))
 112: 		       / (lambda_a - lambda_b);
 113: 		else
 114: 		  v_ab = 0.5 * (  heaviside_function(lambda_a)
 115: 		                + heaviside_function(lambda_b) );
 116: 		positive_projector += 0.5 * v_ab * 0.5 * (G[a][b] + G[b][a]);
 117: 	      }
 118: 	  }
 119:       }
 120: 
 121:     negative_projector = 0;
 122:     for (int a = 0; a < dim; a++)
 123:       {
 124: 	double lambda_a = eigenvalues[a];
 125: 	negative_projector += heaviside_function(-lambda_a)
 126: 			    * Q[a];
 127: 	for (int b = 0; b < dim; b++)
 128: 	  {
 129: 	    if (b != a)
 130: 	      {
 131: 		double lambda_b = eigenvalues[b];
 132: 		double v_ab = 0.0;
 133: 		if (std::fabs(lambda_a - lambda_b) > 1.0e-12)
 134: 		  v_ab = (negative_ramp_function(lambda_a) - negative_ramp_function(lambda_b))
 135: 		       / (lambda_a - lambda_b);
 136: 		else
 137: 		  v_ab = 0.5 * (  heaviside_function(-lambda_a)
 138: 		                + heaviside_function(-lambda_b) );
 139: 		negative_projector += 0.5 * v_ab * 0.5 * (G[a][b] + G[b][a]);
 140: 	      }
 141: 	  }
 142:       }
 143: 
```

### 代码片段 3（SpectrumDecomposition.cc:L19-L38）
```cpp
  19:   double positive_ramp_function(const double x)
  20:   {
  21:     return std::fmax(x, 0.0);
  22:   }
  23: 
  24:   double negative_ramp_function(const double x)
  25:   {
  26:     return std::fmin(x, 0.0);
  27:   }
  28: 
  29:   double heaviside_function(const double x)
  30:   {
  31:     if (std::fabs(x) < 1.0e-16)
  32:       return 0.5;
  33: 
  34:     if (x > 0)
  35:       return 1.0;
  36:     else
  37:       return 0.0;
  38:   }
```

### 公式（与 explain.md 一致）
\[
\frac{\partial \pmb{\sigma}}{\partial \pmb{\epsilon}} = [g(d) + k]\frac{\partial \pmb{\sigma}^{+}}{\partial \pmb{\epsilon}} + \frac{\partial \pmb{\sigma}^{-}}{\partial \pmb{\epsilon}} = [g(d) + k]\left[\lambda H(\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes \mathbf{I} + 2\mu \mathbb{P}^{+}\right] + \left[\lambda H(-\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes \mathbf{I} + 2\mu \mathbb{P}^{-}\right].
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。
