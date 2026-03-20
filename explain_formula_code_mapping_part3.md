## 2. Phase-field formulation and finite element discretization > 2.1. Phase-field formulation

### 未编号公式 #15 — 未编号公式：\pmb{\sigma}^{+} = \lambda...
- 论文位置：[2. Phase-field formulation and finite element discretization > 2.1. Phase-field formulation，`explain.md` 第 131 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L131)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\pmb{\sigma}^{+} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{+}\mathbf{I} + 2\mu \pmb{\epsilon}^{+}, \quad \pmb{\sigma}^{-} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{-}\mathbf{I} + 2\mu \pmb{\epsilon}^{-}.
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - Ramp/Heaviside 基础函数：[SpectrumDecomposition.cc:19-38](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.cc#L19-L38)
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
  - 应变正负分解：[SpectrumDecomposition.h:45-69](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L45-L69)
```cpp
45:   template <int dim>
46:   SymmetricTensor<2, dim> positive_tensor(Vector<double> const & eigenvalues,
47: 					  std::vector<Tensor<1, dim>> const & eigenvectors)
48:   {
49:     SymmetricTensor<2, dim> positive_part_tensor;
50:     positive_part_tensor = 0;
51:     for (int i = 0; i < dim; i++)
52:       positive_part_tensor += positive_ramp_function(eigenvalues[i])
53:                             * symmetrize(outer_product(eigenvectors[i],
54:                                                        eigenvectors[i]));
55:     return positive_part_tensor;
56:   }
57: 
58:   template <int dim>
59:   SymmetricTensor<2, dim> negative_tensor(Vector<double> const & eigenvalues,
60: 					  std::vector<Tensor<1, dim>> const & eigenvectors)
61:   {
62:     SymmetricTensor<2, dim> negative_part_tensor;
63:     negative_part_tensor = 0;
64:     for (int i = 0; i < dim; i++)
65:       negative_part_tensor += negative_ramp_function(eigenvalues[i])
66:                             * symmetrize(outer_product(eigenvectors[i],
67:                                                        eigenvectors[i]));
68:     return negative_part_tensor;
69:   }
```
  - 四阶投影张量 P+/P-：[SpectrumDecomposition.h:71-143](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L71-L143)
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
### 未编号公式 #16 — 未编号公式：\frac{\partial \pmb{\sigma...
- 论文位置：[2. Phase-field formulation and finite element discretization > 2.1. Phase-field formulation，`explain.md` 第 137 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L137)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\frac{\partial \pmb{\sigma}}{\partial \pmb{\epsilon}} = [g(d) + k]\frac{\partial \pmb{\sigma}^{+}}{\partial \pmb{\epsilon}} + \frac{\partial \pmb{\sigma}^{-}}{\partial \pmb{\epsilon}} = [g(d) + k]\left[\lambda H(\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes \mathbf{I} + 2\mu \mathbb{P}^{+}\right] + \left[\lambda H(-\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes \mathbf{I} + 2\mu \mathbb{P}^{-}\right].
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - 退化函数及其导数（最底层）：[main.cc:244-257](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L244-L257)
```cpp
244:   double degradation_function(const double d)
245:   {
246:     return (1.0 - d) * (1.0 - d);
247:   }
248: 
249:   double degradation_function_derivative(const double d)
250:   {
251:     return 2.0 * (d - 1.0);
252:   }
253: 
254:   double degradation_function_2nd_order_derivative(const double d)
255:   {
256:     (void) d;
257:     return 2.0;
```
  - 相场残差中退化函数导数项：[main.cc:3813-3819](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3813-L3819)
```cpp
3813:     	        data.m_cell_rhs(i) += (    gc * length_scale * grad_N_phasefield[i] * phasefield_grad
3814:     	                                +  (   gc / length_scale * phasefield_value
3815: 					     + eta / delta_time  * (phasefield_value - old_phasefield)
3816: 					     + degradation_function_derivative(phasefield_value)
3817: 					     * current_positive_strain_energy )
3818: 					  * N_phasefield[i]
3819: 				      ) * JxW;
```
  - 相场切线中退化函数二阶导项：[main.cc:3946-3951](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3946-L3951)
```cpp
3946:                     data.m_cell_matrix(i, j) += (  (   gc/length_scale + eta/delta_time
3947:                 	                             + degradation_function_2nd_order_derivative(phasefield_value)
3948: 						     * current_positive_strain_energy  )
3949:                 	                          * N_phasefield[i] * N_phasefield[j]
3950: 					          + gc * length_scale * grad_N_phasefield[i] * grad_N_phasefield[j]
3951: 					        ) * JxW;
```
### 未编号公式 #17 — 未编号公式：\mathbb{P}^{+} = \frac{\pa...
- 论文位置：[2. Phase-field formulation and finite element discretization > 2.1. Phase-field formulation，`explain.md` 第 143 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L143)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\mathbb{P}^{+} = \frac{\partial \pmb{\epsilon}^{+}}{\partial \pmb{\epsilon}}, \quad \mathbb{P}^{-} = \frac{\partial \pmb{\epsilon}^{-}}{\partial \pmb{\epsilon}},
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - Ramp/Heaviside 基础函数：[SpectrumDecomposition.cc:19-38](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.cc#L19-L38)
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
  - 应变正负分解：[SpectrumDecomposition.h:45-69](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L45-L69)
```cpp
45:   template <int dim>
46:   SymmetricTensor<2, dim> positive_tensor(Vector<double> const & eigenvalues,
47: 					  std::vector<Tensor<1, dim>> const & eigenvectors)
48:   {
49:     SymmetricTensor<2, dim> positive_part_tensor;
50:     positive_part_tensor = 0;
51:     for (int i = 0; i < dim; i++)
52:       positive_part_tensor += positive_ramp_function(eigenvalues[i])
53:                             * symmetrize(outer_product(eigenvectors[i],
54:                                                        eigenvectors[i]));
55:     return positive_part_tensor;
56:   }
57: 
58:   template <int dim>
59:   SymmetricTensor<2, dim> negative_tensor(Vector<double> const & eigenvalues,
60: 					  std::vector<Tensor<1, dim>> const & eigenvectors)
61:   {
62:     SymmetricTensor<2, dim> negative_part_tensor;
63:     negative_part_tensor = 0;
64:     for (int i = 0; i < dim; i++)
65:       negative_part_tensor += negative_ramp_function(eigenvalues[i])
66:                             * symmetrize(outer_product(eigenvectors[i],
67:                                                        eigenvectors[i]));
68:     return negative_part_tensor;
69:   }
```
  - 四阶投影张量 P+/P-：[SpectrumDecomposition.h:71-143](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L71-L143)
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
### Eq. (7) — 总势能一阶变分/弱式来源
- 论文位置：[2. Phase-field formulation and finite element discretization > 2.1. Phase-field formulation，`explain.md` 第 151 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L151)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\begin{array}{rl} 
\delta \Pi (\pmb {u},d) &= D_{(\delta \pmb {u},\delta d)}\Pi (\pmb {u},d) = \left.\frac{\mathrm{d}}{\mathrm{d}\epsilon}\right|_{\epsilon = 0}\Pi (\pmb {u} + \epsilon \delta \pmb {u},d + \epsilon \delta d) \\ 
&= \int_{\Omega}\left(\frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial \pmb{\epsilon}}:\pmb{\epsilon}(\delta \pmb {u}) + \frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial d}\delta d\right)\mathrm{d}\Omega \\ 
&\quad +\int_{\Omega}\frac{g_{c}}{l}\left(d\delta d + l^{2}\nabla d\cdot \nabla \delta d\right)\mathrm{d}\Omega -\int_{\Omega}\pmb {b}\cdot \delta \pmb {u}\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot \delta \pmb {u}\mathrm{d}\Gamma \\ 
&= (\nabla^{(s)}\delta \pmb {u},\pmb {\sigma}) - (\delta \pmb {u},\pmb {b}) - (\delta \pmb {u},\pmb{t})_{\Gamma_{t}} + (\delta d,\frac{g_{c}}{l} d) + (\nabla \delta d,g_{c}l\nabla d) + (\delta d,g^{\prime}(d)\psi^{+}(\pmb{\epsilon})). 
\end{array} \quad (7)
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - 退化函数及其导数（最底层）：[main.cc:244-257](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L244-L257)
```cpp
244:   double degradation_function(const double d)
245:   {
246:     return (1.0 - d) * (1.0 - d);
247:   }
248: 
249:   double degradation_function_derivative(const double d)
250:   {
251:     return 2.0 * (d - 1.0);
252:   }
253: 
254:   double degradation_function_2nd_order_derivative(const double d)
255:   {
256:     (void) d;
257:     return 2.0;
```
  - 相场残差中退化函数导数项：[main.cc:3813-3819](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3813-L3819)
```cpp
3813:     	        data.m_cell_rhs(i) += (    gc * length_scale * grad_N_phasefield[i] * phasefield_grad
3814:     	                                +  (   gc / length_scale * phasefield_value
3815: 					     + eta / delta_time  * (phasefield_value - old_phasefield)
3816: 					     + degradation_function_derivative(phasefield_value)
3817: 					     * current_positive_strain_energy )
3818: 					  * N_phasefield[i]
3819: 				      ) * JxW;
```
  - 相场切线中退化函数二阶导项：[main.cc:3946-3951](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3946-L3951)
```cpp
3946:                     data.m_cell_matrix(i, j) += (  (   gc/length_scale + eta/delta_time
3947:                 	                             + degradation_function_2nd_order_derivative(phasefield_value)
3948: 						     * current_positive_strain_energy  )
3949:                 	                          * N_phasefield[i] * N_phasefield[j]
3950: 					          + gc * length_scale * grad_N_phasefield[i] * grad_N_phasefield[j]
3951: 					        ) * JxW;
```
### 未编号公式 #19 — 未编号公式：\left\{ \begin{array}{ll} ...
- 论文位置：[2. Phase-field formulation and finite element discretization > 2.1. Phase-field formulation，`explain.md` 第 162 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L162)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\left\{ \begin{array}{ll} r_{\pmb{u}}(\pmb{u},d) = (\nabla^{(s)}\delta \pmb{u},\pmb{\sigma}) - (\delta \pmb{u},\pmb{b}) - (\delta \pmb{u},\pmb{t})_{\Gamma_{t}} = 0, \\ r_{d}(\pmb{u},d) = (\delta d,\frac{g_{c}}{l} d) + (\nabla \delta d,g_{c}l\nabla d) + (\delta d,g^{\prime}(d)\psi^{+}(\pmb{\epsilon})) = 0, \end{array} \right.
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - 退化函数及其导数（最底层）：[main.cc:244-257](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L244-L257)
```cpp
244:   double degradation_function(const double d)
245:   {
246:     return (1.0 - d) * (1.0 - d);
247:   }
248: 
249:   double degradation_function_derivative(const double d)
250:   {
251:     return 2.0 * (d - 1.0);
252:   }
253: 
254:   double degradation_function_2nd_order_derivative(const double d)
255:   {
256:     (void) d;
257:     return 2.0;
```
  - 相场残差中退化函数导数项：[main.cc:3813-3819](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3813-L3819)
```cpp
3813:     	        data.m_cell_rhs(i) += (    gc * length_scale * grad_N_phasefield[i] * phasefield_grad
3814:     	                                +  (   gc / length_scale * phasefield_value
3815: 					     + eta / delta_time  * (phasefield_value - old_phasefield)
3816: 					     + degradation_function_derivative(phasefield_value)
3817: 					     * current_positive_strain_energy )
3818: 					  * N_phasefield[i]
3819: 				      ) * JxW;
```
  - 相场切线中退化函数二阶导项：[main.cc:3946-3951](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3946-L3951)
```cpp
3946:                     data.m_cell_matrix(i, j) += (  (   gc/length_scale + eta/delta_time
3947:                 	                             + degradation_function_2nd_order_derivative(phasefield_value)
3948: 						     * current_positive_strain_energy  )
3949:                 	                          * N_phasefield[i] * N_phasefield[j]
3950: 					          + gc * length_scale * grad_N_phasefield[i] * grad_N_phasefield[j]
3951: 					        ) * JxW;
```

## 2. Phase-field formulation and finite element discretization > 2.2. Finite element discretization

### 未编号公式 #20 — 未编号公式：\pmb {u} = \pmb{N}_{u_{A}}...
- 论文位置：[2. Phase-field formulation and finite element discretization > 2.2. Finite element discretization，`explain.md` 第 172 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L172)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\pmb {u} = \pmb{N}_{u_{A}}\pmb{u}_{A}, \quad \text{and} \quad d = N_{d_{A}}d_{A}.
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - 离散残差与能量相关核心装配：[main.cc:3681-3820](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3681-L3820)
```cpp
3681:   void PhaseFieldMonolithicSolve<dim>::assemble_system_rhs_BFGS_parallel(const BlockVector<double> & solution_old,
3682: 								         BlockVector<double> & system_rhs)
3683:   {
3684:     m_timer.enter_subsection("Assemble RHS");
3685: 
3686:     system_rhs = 0.0;
3687: 
3688:     const UpdateFlags uf_cell(update_values | update_gradients |
3689: 			      update_quadrature_points | update_JxW_values);
3690:     const UpdateFlags uf_face(update_values | update_normal_vectors |
3691: 			      update_JxW_values);
3692: 
3693:     PerTaskData_ASM_RHS_BFGS per_task_data(m_fe.n_dofs_per_cell());
3694:     ScratchData_ASM_RHS_BFGS scratch_data(m_fe, m_qf_cell, uf_cell, m_qf_face, uf_face, solution_old);
3695: 
3696:     auto worker =
3697:       [this](const typename DoFHandler<dim>::active_cell_iterator &cell,
3698: 	     ScratchData_ASM_RHS_BFGS & scratch,
3699: 	     PerTaskData_ASM_RHS_BFGS & data)
3700:       {
3701:         this->assemble_system_rhs_BFGS_one_cell(cell, scratch, data);
3702:       };
3703: 
3704:     auto copier = [this, &system_rhs](const PerTaskData_ASM_RHS_BFGS &data)
3705:       {
3706:         this->m_constraints.distribute_local_to_global(data.m_cell_rhs,
3707:                                                        data.m_local_dof_indices,
3708: 						       system_rhs);
3709:       };
3710: 
3711:     WorkStream::run(
3712:       m_dof_handler.active_cell_iterators(),
3713:       worker,
3714:       copier,
3715:       scratch_data,
3716:       per_task_data);
3717: 
3718:     m_timer.leave_subsection();
3719:   }
3720: 
3721:   template <int dim>
3722:   void PhaseFieldMonolithicSolve<dim>::assemble_system_rhs_BFGS_one_cell(
3723:       const typename DoFHandler<dim>::active_cell_iterator &cell,
3724:       ScratchData_ASM_RHS_BFGS & scratch,
3725:       PerTaskData_ASM_RHS_BFGS & data) const
3726:   {
3727:     data.reset();
3728:     scratch.reset();
3729:     scratch.m_fe_values.reinit(cell);
3730:     cell->get_dof_indices(data.m_local_dof_indices);
3731: 
3732:     scratch.m_fe_values[m_d_fe].get_function_values(
3733:       scratch.m_solution_previous_step, scratch.m_phasefield_previous_step_cell);
3734: 
3735:     const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
3736:       m_quadrature_point_history.get_data(cell);
3737:     Assert(lqph.size() == m_n_q_points, ExcInternalError());
3738: 
3739:     const double time_ramp = (m_time.current() / m_time.end());
3740:     std::vector<Tensor<1, dim>> rhs_values(m_n_q_points);
3741: 
3742:     right_hand_side(scratch.m_fe_values.get_quadrature_points(),
3743: 		    rhs_values,
3744: 		    m_parameters.m_x_component*1.0,
3745: 		    m_parameters.m_y_component*1.0,
3746: 		    m_parameters.m_z_component*1.0);
3747: 
3748:     const double delta_time = m_time.get_delta_t();
3749: 
3750:     for (const unsigned int q_point : scratch.m_fe_values.quadrature_point_indices())
3751:       {
3752:         for (const unsigned int k : scratch.m_fe_values.dof_indices())
3753:           {
3754:             const unsigned int k_group = m_fe.system_to_base_index(k).first.first;
3755: 
3756:             if (k_group == m_u_dof)
3757:               {
3758:                 scratch.m_Nx_disp[q_point][k] =
3759:                   scratch.m_fe_values[m_u_fe].value(k, q_point);
3760:                 scratch.m_grad_Nx_disp[q_point][k] =
3761:                   scratch.m_fe_values[m_u_fe].gradient(k, q_point);
3762:                 scratch.m_symm_grad_Nx_disp[q_point][k] =
3763:                   symmetrize(scratch.m_grad_Nx_disp[q_point][k]);
3764:               }
3765:             else if (k_group == m_d_dof)
3766:               {
3767: 		scratch.m_Nx_phasefield[q_point][k] =
3768: 		  scratch.m_fe_values[m_d_fe].value(k, q_point);
3769: 		scratch.m_grad_Nx_phasefield[q_point][k] =
3770: 		  scratch.m_fe_values[m_d_fe].gradient(k, q_point);
3771:               }
3772:             else
3773:               Assert(k_group <= m_d_dof, ExcInternalError());
3774:           }
3775:       }
3776: 
3777:     for (const unsigned int q_point : scratch.m_fe_values.quadrature_point_indices())
3778:       {
3779: 	const double length_scale            = lqph[q_point]->get_length_scale();
3780: 	const double gc                      = lqph[q_point]->get_critical_energy_release_rate();
3781: 	const double eta                     = lqph[q_point]->get_viscosity();
3782: 	const double current_positive_strain_energy = lqph[q_point]->get_current_positive_strain_energy();
3783: 
3784: 	const double phasefield_value        = lqph[q_point]->get_phase_field_value();
3785: 	const Tensor<1, dim> phasefield_grad = lqph[q_point]->get_phase_field_gradient();
3786: 
3787:         const std::vector<double>         &      N_phasefield = scratch.m_Nx_phasefield[q_point];
3788:         const std::vector<Tensor<1, dim>> & grad_N_phasefield = scratch.m_grad_Nx_phasefield[q_point];
3789:         const double                old_phasefield = scratch.m_phasefield_previous_step_cell[q_point];
3790: 
3791:         const SymmetricTensor<2, dim> & cauchy_stress = lqph[q_point]->get_cauchy_stress();
3792: 
3793:         const std::vector<Tensor<1,dim>> & N_disp = scratch.m_Nx_disp[q_point];
3794:         const std::vector<SymmetricTensor<2, dim>> & symm_grad_N_disp =
3795:           scratch.m_symm_grad_Nx_disp[q_point];
3796:         const double JxW = scratch.m_fe_values.JxW(q_point);
3797: 
3798:         SymmetricTensor<2, dim> symm_grad_Nx_i_x_C;
3799: 
3800:         for (const unsigned int i : scratch.m_fe_values.dof_indices())
3801:           {
3802:             const unsigned int i_group = m_fe.system_to_base_index(i).first.first;
3803: 
3804:             if (i_group == m_u_dof)
3805:               {
3806:                 data.m_cell_rhs(i) += (symm_grad_N_disp[i] * cauchy_stress) * JxW;
3807: 
3808: 		// contributions from the body force to right-hand side
3809: 		data.m_cell_rhs(i) -= N_disp[i] * rhs_values[q_point] * JxW;
3810:               }
3811:             else if (i_group == m_d_dof)
3812:               {
3813:     	        data.m_cell_rhs(i) += (    gc * length_scale * grad_N_phasefield[i] * phasefield_grad
3814:     	                                +  (   gc / length_scale * phasefield_value
3815: 					     + eta / delta_time  * (phasefield_value - old_phasefield)
3816: 					     + degradation_function_derivative(phasefield_value)
3817: 					     * current_positive_strain_energy )
3818: 					  * N_phasefield[i]
3819: 				      ) * JxW;
3820:               }
```
  - 目标能量积分计算：[main.cc:6117-6140](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L6117-L6140)
```cpp
6117:   double PhaseFieldMonolithicSolve<dim>::calculate_energy_functional() const
6118:   {
6119:     double energy_functional = 0.0;
6120: 
6121:     FEValues<dim> fe_values(m_fe, m_qf_cell, update_JxW_values);
6122: 
6123:     for (const auto &cell : m_dof_handler.active_cell_iterators())
6124:       {
6125:         fe_values.reinit(cell);
6126: 
6127:         const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
6128:           m_quadrature_point_history.get_data(cell);
6129:         Assert(lqph.size() == m_n_q_points, ExcInternalError());
6130: 
6131:         for (unsigned int q_point = 0; q_point < m_n_q_points; ++q_point)
6132:           {
6133:             const double JxW = fe_values.JxW(q_point);
6134:             energy_functional += lqph[q_point]->get_total_strain_energy() * JxW;
6135:             energy_functional += lqph[q_point]->get_crack_energy_dissipation() * JxW;
6136:           }
6137:       }
6138: 
6139:     return energy_functional;
6140:   }
```
### 未编号公式 #21 — 未编号公式：\delta \pmb {u} = \pmb{N}_...
- 论文位置：[2. Phase-field formulation and finite element discretization > 2.2. Finite element discretization，`explain.md` 第 178 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L178)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\delta \pmb {u} = \pmb{N}_{u_{A}}\delta \pmb{u}_{A}, \quad \text{and} \quad \delta d = N_{d_{A}}\delta d_{A},
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - 离散残差与能量相关核心装配：[main.cc:3681-3820](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3681-L3820)
```cpp
3681:   void PhaseFieldMonolithicSolve<dim>::assemble_system_rhs_BFGS_parallel(const BlockVector<double> & solution_old,
3682: 								         BlockVector<double> & system_rhs)
3683:   {
3684:     m_timer.enter_subsection("Assemble RHS");
3685: 
3686:     system_rhs = 0.0;
3687: 
3688:     const UpdateFlags uf_cell(update_values | update_gradients |
3689: 			      update_quadrature_points | update_JxW_values);
3690:     const UpdateFlags uf_face(update_values | update_normal_vectors |
3691: 			      update_JxW_values);
3692: 
3693:     PerTaskData_ASM_RHS_BFGS per_task_data(m_fe.n_dofs_per_cell());
3694:     ScratchData_ASM_RHS_BFGS scratch_data(m_fe, m_qf_cell, uf_cell, m_qf_face, uf_face, solution_old);
3695: 
3696:     auto worker =
3697:       [this](const typename DoFHandler<dim>::active_cell_iterator &cell,
3698: 	     ScratchData_ASM_RHS_BFGS & scratch,
3699: 	     PerTaskData_ASM_RHS_BFGS & data)
3700:       {
3701:         this->assemble_system_rhs_BFGS_one_cell(cell, scratch, data);
3702:       };
3703: 
3704:     auto copier = [this, &system_rhs](const PerTaskData_ASM_RHS_BFGS &data)
3705:       {
3706:         this->m_constraints.distribute_local_to_global(data.m_cell_rhs,
3707:                                                        data.m_local_dof_indices,
3708: 						       system_rhs);
3709:       };
3710: 
3711:     WorkStream::run(
3712:       m_dof_handler.active_cell_iterators(),
3713:       worker,
3714:       copier,
3715:       scratch_data,
3716:       per_task_data);
3717: 
3718:     m_timer.leave_subsection();
3719:   }
3720: 
3721:   template <int dim>
3722:   void PhaseFieldMonolithicSolve<dim>::assemble_system_rhs_BFGS_one_cell(
3723:       const typename DoFHandler<dim>::active_cell_iterator &cell,
3724:       ScratchData_ASM_RHS_BFGS & scratch,
3725:       PerTaskData_ASM_RHS_BFGS & data) const
3726:   {
3727:     data.reset();
3728:     scratch.reset();
3729:     scratch.m_fe_values.reinit(cell);
3730:     cell->get_dof_indices(data.m_local_dof_indices);
3731: 
3732:     scratch.m_fe_values[m_d_fe].get_function_values(
3733:       scratch.m_solution_previous_step, scratch.m_phasefield_previous_step_cell);
3734: 
3735:     const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
3736:       m_quadrature_point_history.get_data(cell);
3737:     Assert(lqph.size() == m_n_q_points, ExcInternalError());
3738: 
3739:     const double time_ramp = (m_time.current() / m_time.end());
3740:     std::vector<Tensor<1, dim>> rhs_values(m_n_q_points);
3741: 
3742:     right_hand_side(scratch.m_fe_values.get_quadrature_points(),
3743: 		    rhs_values,
3744: 		    m_parameters.m_x_component*1.0,
3745: 		    m_parameters.m_y_component*1.0,
3746: 		    m_parameters.m_z_component*1.0);
3747: 
3748:     const double delta_time = m_time.get_delta_t();
3749: 
3750:     for (const unsigned int q_point : scratch.m_fe_values.quadrature_point_indices())
3751:       {
3752:         for (const unsigned int k : scratch.m_fe_values.dof_indices())
3753:           {
3754:             const unsigned int k_group = m_fe.system_to_base_index(k).first.first;
3755: 
3756:             if (k_group == m_u_dof)
3757:               {
3758:                 scratch.m_Nx_disp[q_point][k] =
3759:                   scratch.m_fe_values[m_u_fe].value(k, q_point);
3760:                 scratch.m_grad_Nx_disp[q_point][k] =
3761:                   scratch.m_fe_values[m_u_fe].gradient(k, q_point);
3762:                 scratch.m_symm_grad_Nx_disp[q_point][k] =
3763:                   symmetrize(scratch.m_grad_Nx_disp[q_point][k]);
3764:               }
3765:             else if (k_group == m_d_dof)
3766:               {
3767: 		scratch.m_Nx_phasefield[q_point][k] =
3768: 		  scratch.m_fe_values[m_d_fe].value(k, q_point);
3769: 		scratch.m_grad_Nx_phasefield[q_point][k] =
3770: 		  scratch.m_fe_values[m_d_fe].gradient(k, q_point);
3771:               }
3772:             else
3773:               Assert(k_group <= m_d_dof, ExcInternalError());
3774:           }
3775:       }
3776: 
3777:     for (const unsigned int q_point : scratch.m_fe_values.quadrature_point_indices())
3778:       {
3779: 	const double length_scale            = lqph[q_point]->get_length_scale();
3780: 	const double gc                      = lqph[q_point]->get_critical_energy_release_rate();
3781: 	const double eta                     = lqph[q_point]->get_viscosity();
3782: 	const double current_positive_strain_energy = lqph[q_point]->get_current_positive_strain_energy();
3783: 
3784: 	const double phasefield_value        = lqph[q_point]->get_phase_field_value();
3785: 	const Tensor<1, dim> phasefield_grad = lqph[q_point]->get_phase_field_gradient();
3786: 
3787:         const std::vector<double>         &      N_phasefield = scratch.m_Nx_phasefield[q_point];
3788:         const std::vector<Tensor<1, dim>> & grad_N_phasefield = scratch.m_grad_Nx_phasefield[q_point];
3789:         const double                old_phasefield = scratch.m_phasefield_previous_step_cell[q_point];
3790: 
3791:         const SymmetricTensor<2, dim> & cauchy_stress = lqph[q_point]->get_cauchy_stress();
3792: 
3793:         const std::vector<Tensor<1,dim>> & N_disp = scratch.m_Nx_disp[q_point];
3794:         const std::vector<SymmetricTensor<2, dim>> & symm_grad_N_disp =
3795:           scratch.m_symm_grad_Nx_disp[q_point];
3796:         const double JxW = scratch.m_fe_values.JxW(q_point);
3797: 
3798:         SymmetricTensor<2, dim> symm_grad_Nx_i_x_C;
3799: 
3800:         for (const unsigned int i : scratch.m_fe_values.dof_indices())
3801:           {
3802:             const unsigned int i_group = m_fe.system_to_base_index(i).first.first;
3803: 
3804:             if (i_group == m_u_dof)
3805:               {
3806:                 data.m_cell_rhs(i) += (symm_grad_N_disp[i] * cauchy_stress) * JxW;
3807: 
3808: 		// contributions from the body force to right-hand side
3809: 		data.m_cell_rhs(i) -= N_disp[i] * rhs_values[q_point] * JxW;
3810:               }
3811:             else if (i_group == m_d_dof)
3812:               {
3813:     	        data.m_cell_rhs(i) += (    gc * length_scale * grad_N_phasefield[i] * phasefield_grad
3814:     	                                +  (   gc / length_scale * phasefield_value
3815: 					     + eta / delta_time  * (phasefield_value - old_phasefield)
3816: 					     + degradation_function_derivative(phasefield_value)
3817: 					     * current_positive_strain_energy )
3818: 					  * N_phasefield[i]
3819: 				      ) * JxW;
3820:               }
```
  - 目标能量积分计算：[main.cc:6117-6140](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L6117-L6140)
```cpp
6117:   double PhaseFieldMonolithicSolve<dim>::calculate_energy_functional() const
6118:   {
6119:     double energy_functional = 0.0;
6120: 
6121:     FEValues<dim> fe_values(m_fe, m_qf_cell, update_JxW_values);
6122: 
6123:     for (const auto &cell : m_dof_handler.active_cell_iterators())
6124:       {
6125:         fe_values.reinit(cell);
6126: 
6127:         const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
6128:           m_quadrature_point_history.get_data(cell);
6129:         Assert(lqph.size() == m_n_q_points, ExcInternalError());
6130: 
6131:         for (unsigned int q_point = 0; q_point < m_n_q_points; ++q_point)
6132:           {
6133:             const double JxW = fe_values.JxW(q_point);
6134:             energy_functional += lqph[q_point]->get_total_strain_energy() * JxW;
6135:             energy_functional += lqph[q_point]->get_crack_energy_dissipation() * JxW;
6136:           }
6137:       }
6138: 
6139:     return energy_functional;
6140:   }
```
### Eq. (8) — 有限元离散总能量
- 论文位置：[2. Phase-field formulation and finite element discretization > 2.2. Finite element discretization，`explain.md` 第 184 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L184)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\begin{array}{rl} 
\Pi (\pmb {u}_{A},d_{A}) = &\int_{\Omega}\psi \left(\pmb{\epsilon}(\pmb {N}_{u_{A}}\pmb {u}_{A}),N_{d_{A}}d_{A}\right)\mathrm{d}\Omega \\ 
&+\int_{\Omega}\frac{g_{c}}{2l}\left((N_{d_{A}}d_{A})^{2} + l^{2}(\nabla N_{d_{A}}d_{A})\cdot (\nabla N_{d_{A}}d_{A})\right)\mathrm{d}\Omega \\ 
&-\int_{\Omega}\pmb {b}\cdot (\pmb {N}_{u_{A}}\pmb {u}_{A})\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot (\pmb {N}_{u_{A}}\pmb {u}_{A})\mathrm{d}\Gamma, 
\end{array} \quad (8)
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - Ramp/Heaviside 基础函数：[SpectrumDecomposition.cc:19-38](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.cc#L19-L38)
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
  - 应变正负分解：[SpectrumDecomposition.h:45-69](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L45-L69)
```cpp
45:   template <int dim>
46:   SymmetricTensor<2, dim> positive_tensor(Vector<double> const & eigenvalues,
47: 					  std::vector<Tensor<1, dim>> const & eigenvectors)
48:   {
49:     SymmetricTensor<2, dim> positive_part_tensor;
50:     positive_part_tensor = 0;
51:     for (int i = 0; i < dim; i++)
52:       positive_part_tensor += positive_ramp_function(eigenvalues[i])
53:                             * symmetrize(outer_product(eigenvectors[i],
54:                                                        eigenvectors[i]));
55:     return positive_part_tensor;
56:   }
57: 
58:   template <int dim>
59:   SymmetricTensor<2, dim> negative_tensor(Vector<double> const & eigenvalues,
60: 					  std::vector<Tensor<1, dim>> const & eigenvectors)
61:   {
62:     SymmetricTensor<2, dim> negative_part_tensor;
63:     negative_part_tensor = 0;
64:     for (int i = 0; i < dim; i++)
65:       negative_part_tensor += negative_ramp_function(eigenvalues[i])
66:                             * symmetrize(outer_product(eigenvectors[i],
67:                                                        eigenvectors[i]));
68:     return negative_part_tensor;
69:   }
```
  - 四阶投影张量 P+/P-：[SpectrumDecomposition.h:71-143](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L71-L143)
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
### Eq. (9) — 离散梯度(残差)向量
- 论文位置：[2. Phase-field formulation and finite element discretization > 2.2. Finite element discretization，`explain.md` 第 194 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L194)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\begin{array}{rl} 
\pmb{r} = \nabla \Pi = (r_{\pmb{u}},r_{d})^{\mathrm{T}},\\ 
r_{\pmb{u}_{A}} = \frac{\partial\Pi}{\partial \pmb{u}_{A}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},\pmb{\sigma}\right) - \left(\pmb{N}_{u_{A}},\pmb{b}\right) - \left(\pmb{N}_{u_{A}},\pmb{t}\right)_{\Gamma_{t}},\\ 
r_{d_{A}} = \frac{\partial\Pi}{\partial d_{A}} = \left(N_{d_{A}},\frac{g_{c}}{l} d + g^{\prime}(d)\psi^{+}\right) + (\nabla N_{d_{A}},g_{c}l\nabla d). 
\end{array} \quad (9)
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - 退化函数及其导数（最底层）：[main.cc:244-257](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L244-L257)
```cpp
244:   double degradation_function(const double d)
245:   {
246:     return (1.0 - d) * (1.0 - d);
247:   }
248: 
249:   double degradation_function_derivative(const double d)
250:   {
251:     return 2.0 * (d - 1.0);
252:   }
253: 
254:   double degradation_function_2nd_order_derivative(const double d)
255:   {
256:     (void) d;
257:     return 2.0;
```
  - 相场残差中退化函数导数项：[main.cc:3813-3819](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3813-L3819)
```cpp
3813:     	        data.m_cell_rhs(i) += (    gc * length_scale * grad_N_phasefield[i] * phasefield_grad
3814:     	                                +  (   gc / length_scale * phasefield_value
3815: 					     + eta / delta_time  * (phasefield_value - old_phasefield)
3816: 					     + degradation_function_derivative(phasefield_value)
3817: 					     * current_positive_strain_energy )
3818: 					  * N_phasefield[i]
3819: 				      ) * JxW;
```
  - 相场切线中退化函数二阶导项：[main.cc:3946-3951](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3946-L3951)
```cpp
3946:                     data.m_cell_matrix(i, j) += (  (   gc/length_scale + eta/delta_time
3947:                 	                             + degradation_function_2nd_order_derivative(phasefield_value)
3948: 						     * current_positive_strain_energy  )
3949:                 	                          * N_phasefield[i] * N_phasefield[j]
3950: 					          + gc * length_scale * grad_N_phasefield[i] * grad_N_phasefield[j]
3951: 					        ) * JxW;
```
