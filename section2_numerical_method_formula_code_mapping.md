# explain.md 第 2 章（Numerical method）公式-代码逐段对应解读

> 说明：按你的要求，严格使用“**代码在前、公式在后**”格式；每个小节都给出**可访问链接**（GitHub 行号链接）+ **完整代码片段显示** + **公式渲染**。  
> 仓库链接基准分支使用稳定分支：`main`。

---

## 2.1 Phase-field formulation

### 2.1-1 公式（5）中的退化函数分裂项 \[ \psi(\boldsymbol\epsilon,d)=[g(d)+k]\psi^+(\boldsymbol\epsilon)+\psi^-(\boldsymbol\epsilon) \]

**对应最底层代码链接：**
- 应力/能量里真正应用 \[g(d)+k\] 的位置（核心执行点）  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L846-L885

**代码：**
```cpp
const double degradation = degradation_function(m_phase_field_value) + m_residual_k;
const double I_1 = trace(m_strain);

// ...

m_stress = degradation * stress_positive + stress_negative;
m_stress_positive = stress_positive;

// ...

m_mechanical_C = degradation * C_positive + C_negative;

m_strain_energy_total = degradation * m_strain_energy_positive + m_strain_energy_negative;
```

**公式（原文对应）：**
\[
\psi(\pmb{\epsilon}, d) = [g(d) + k] \psi^{+}(\pmb{\epsilon}) + \psi^{-}(\pmb{\epsilon}). \quad (5)
\]

---

### 2.1-2 公式（6） \[ g(d)=(1-d)^2 \]

**对应最底层代码链接：**
- 退化函数定义  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L244-L257

**代码：**
```cpp
double degradation_function(const double d)
{
  return (1.0 - d) * (1.0 - d);
}

double degradation_function_derivative(const double d)
{
  return 2.0 * (d - 1.0);
}

double degradation_function_2nd_order_derivative(const double d)
{
  (void) d;
  return 2.0;
}
```

**公式（原文对应）：**
\[
g(d) = (1 - d)^{2}. \quad (6)
\]

---

### 2.1-3 算子定义 \[\langle x\rangle_+,\langle x\rangle_-,H(x)\]

**对应最底层代码链接：**
- 正/负 ramp 与 Heaviside 的实现  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.cc#L19-L38

**代码：**
```cpp
double positive_ramp_function(const double x)
{
  return std::fmax(x, 0.0);
}

double negative_ramp_function(const double x)
{
  return std::fmin(x, 0.0);
}

double heaviside_function(const double x)
{
  if (std::fabs(x) < 1.0e-16)
    return 0.5;

  if (x > 0)
    return 1.0;
  else
    return 0.0;
}
```

**公式（原文对应）：**
\[
\langle x\rangle_{+} = \frac{1}{2}(x + |x|), \quad \langle x\rangle_{-} = \frac{1}{2}(x - |x|), \quad
H(x)=\begin{cases}1,&x\ge 0\\0,&x<0\end{cases}
\]

---

### 2.1-4 谱分解 \[\boldsymbol\epsilon=\sum_\alpha \epsilon_\alpha \mathbf M_\alpha\]

**对应最底层代码链接：**
- 应变谱分解入口  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L29-L43
- 在材料点调用谱分解  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L829-L833

**代码：**
```cpp
template <int dim>
void spectrum_decomposition(SymmetricTensor<2, dim> const & symmetric_tensor,
                            Vector<double> & myEigenvalues,
                            std::vector<Tensor<1, dim>> & myEigenvectors)
{
  const std::array< std::pair< double, Tensor< 1, dim > >, dim >
    myEigenSystem = eigenvectors(symmetric_tensor);

  for (int i = 0; i < dim; i++)
    {
      myEigenvalues[i] = myEigenSystem[i].first;
      myEigenvectors[i] = myEigenSystem[i].second;
    }
}

// in main.cc
Vector<double>              eigenvalues(dim);
std::vector<Tensor<1, dim>> eigenvectors(dim);
usr_spectrum_decomposition::spectrum_decomposition<dim>(m_strain,
                                                         eigenvalues,
                                                         eigenvectors);
```

**公式（原文对应）：**
\[
\pmb{\epsilon} = \sum_{\alpha} \epsilon_{\alpha} \mathbf{M}_{\alpha},\quad
\mathbf{M}_{\alpha}=\pmb n_{\alpha}\otimes\pmb n_{\alpha}.
\]

---

### 2.1-5 \[\boldsymbol\epsilon^+\] 与 \[\boldsymbol\epsilon^-\] 定义

**对应最底层代码链接：**
- 正/负应变张量构造  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L45-L69
- 材料点调用  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L835-L837

**代码：**
```cpp
template <int dim>
SymmetricTensor<2, dim> positive_tensor(Vector<double> const & eigenvalues,
                                        std::vector<Tensor<1, dim>> const & eigenvectors)
{
  SymmetricTensor<2, dim> positive_part_tensor;
  positive_part_tensor = 0;
  for (int i = 0; i < dim; i++)
    positive_part_tensor += positive_ramp_function(eigenvalues[i])
                          * symmetrize(outer_product(eigenvectors[i],
                                                     eigenvectors[i]));
  return positive_part_tensor;
}

template <int dim>
SymmetricTensor<2, dim> negative_tensor(Vector<double> const & eigenvalues,
                                        std::vector<Tensor<1, dim>> const & eigenvectors)
{
  SymmetricTensor<2, dim> negative_part_tensor;
  negative_part_tensor = 0;
  for (int i = 0; i < dim; i++)
    negative_part_tensor += negative_ramp_function(eigenvalues[i])
                          * symmetrize(outer_product(eigenvectors[i],
                                                     eigenvectors[i]));
  return negative_part_tensor;
}

// in main.cc
SymmetricTensor<2, dim> strain_positive, strain_negative;
strain_positive = usr_spectrum_decomposition::positive_tensor(eigenvalues, eigenvectors);
strain_negative = usr_spectrum_decomposition::negative_tensor(eigenvalues, eigenvectors);
```

**公式（原文对应）：**
\[
\pmb{\epsilon}^{+} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{+}\mathbf{M}_{\alpha},
\quad
\pmb{\epsilon}^{-} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{-}\mathbf{M}_{\alpha}.
\]

---

### 2.1-6 \[\psi^+(\boldsymbol\epsilon),\psi^-(\boldsymbol\epsilon)\] 的计算

**对应最底层代码链接：**
- 正负应变能密度（材料点）  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L876-L883

**代码：**
```cpp
m_strain_energy_positive = 0.5 * my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)
                                               * usr_spectrum_decomposition::positive_ramp_function(I_1)
                         + m_lame_mu * strain_positive * strain_positive;

m_strain_energy_negative = 0.5 * my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)
                                               * usr_spectrum_decomposition::negative_ramp_function(I_1)
                         + m_lame_mu * strain_negative * strain_negative;
```

**公式（原文对应）：**
\[
\psi^{+}(\pmb{\epsilon}) = \frac{1}{2}\lambda\langle\mathrm{tr}\pmb\epsilon\rangle_{+}^{2} + \mu\,\pmb\epsilon^{+}:\pmb\epsilon^{+},
\quad
\psi^{-}(\pmb{\epsilon}) = \frac{1}{2}\lambda\langle\mathrm{tr}\pmb\epsilon\rangle_{-}^{2} + \mu\,\pmb\epsilon^{-}:\pmb\epsilon^{-}.
\]

---

### 2.1-7 \[\boldsymbol\sigma=[g(d)+k]\boldsymbol\sigma^+ + \boldsymbol\sigma^-\]

**对应最底层代码链接：**
- 正负应力 + 合成总应力  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L857-L865

**代码：**
```cpp
stress_positive = my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)
                                * Physics::Elasticity::StandardTensors<dim>::I
                + 2 * m_lame_mu * strain_positive;
stress_negative = my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)
                                * Physics::Elasticity::StandardTensors<dim>::I
                + 2 * m_lame_mu * strain_negative;

m_stress = degradation * stress_positive + stress_negative;
m_stress_positive = stress_positive;
```

**公式（原文对应）：**
\[
\pmb\sigma
=\frac{\partial\psi(\pmb\epsilon(\pmb u),d)}{\partial\pmb\epsilon}
=[g(d)+k]\frac{\partial\psi^+}{\partial\pmb\epsilon}+\frac{\partial\psi^-}{\partial\pmb\epsilon}
=[g(d)+k]\pmb\sigma^+ + \pmb\sigma^-.
\]

---

### 2.1-8 \[\boldsymbol\sigma^+\] 与 \[\boldsymbol\sigma^-\] 显式表达

**对应最底层代码链接：**
- 与上式同一实现（同一最底层计算）  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L857-L862

**代码：**
```cpp
stress_positive = my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)
                                * Physics::Elasticity::StandardTensors<dim>::I
                + 2 * m_lame_mu * strain_positive;
stress_negative = my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)
                                * Physics::Elasticity::StandardTensors<dim>::I
                + 2 * m_lame_mu * strain_negative;
```

**公式（原文对应）：**
\[
\pmb\sigma^+ = \lambda\langle\mathrm{tr}\pmb\epsilon\rangle_+\mathbf I + 2\mu\pmb\epsilon^+,
\quad
\pmb\sigma^- = \lambda\langle\mathrm{tr}\pmb\epsilon\rangle_-\mathbf I + 2\mu\pmb\epsilon^-.
\]

---

### 2.1-9 切线模量 \[\partial\boldsymbol\sigma/\partial\boldsymbol\epsilon\]

**对应最底层代码链接：**
- \[\mathbf C^+\] / \[\mathbf C^-\] 与总切线 \[\mathbf C\]  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L867-L875
- 投影算子 \[\mathbb P^+,\mathbb P^-\]  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L72-L143

**代码：**
```cpp
// main.cc
SymmetricTensor<4, dim> C_positive, C_negative;
C_positive = my_lambda * usr_spectrum_decomposition::heaviside_function(I_1)
                           * Physics::Elasticity::StandardTensors<dim>::IxI
           + 2 * m_lame_mu * projector_positive;
C_negative = my_lambda * usr_spectrum_decomposition::heaviside_function(-I_1)
                           * Physics::Elasticity::StandardTensors<dim>::IxI
           + 2 * m_lame_mu * projector_negative;
m_mechanical_C = degradation * C_positive + C_negative;

// SpectrumDecomposition.h (projector 关键片段)
positive_projector += heaviside_function(lambda_a) * Q[a];
// ...
v_ab = (positive_ramp_function(lambda_a) - positive_ramp_function(lambda_b))
     / (lambda_a - lambda_b);
positive_projector += 0.5 * v_ab * 0.5 * (G[a][b] + G[b][a]);

negative_projector += heaviside_function(-lambda_a) * Q[a];
// ...
v_ab = (negative_ramp_function(lambda_a) - negative_ramp_function(lambda_b))
     / (lambda_a - lambda_b);
negative_projector += 0.5 * v_ab * 0.5 * (G[a][b] + G[b][a]);
```

**公式（原文对应）：**
\[
\frac{\partial\pmb\sigma}{\partial\pmb\epsilon}
=[g(d)+k]\left[\lambda H(\mathrm{tr}\pmb\epsilon)\mathbf I\otimes\mathbf I +2\mu\mathbb P^+\right]
+\left[\lambda H(-\mathrm{tr}\pmb\epsilon)\mathbf I\otimes\mathbf I +2\mu\mathbb P^-\right].
\]

---

### 2.1-10 投影算子定义 \[ \mathbb P^+=\partial\boldsymbol\epsilon^+/\partial\boldsymbol\epsilon,\;\mathbb P^-=\partial\boldsymbol\epsilon^-/\partial\boldsymbol\epsilon \]

**对应最底层代码链接：**
- 正负投影算子完整实现  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L72-L143

**代码：**
```cpp
template <int dim>
void positive_negative_projectors(Vector<double> const & eigenvalues,
                                  std::vector<Tensor<1, dim>> const & eigenvectors,
                                  SymmetricTensor<4, dim> & positive_projector,
                                  SymmetricTensor<4, dim> & negative_projector)
{
  // ... build M, Q, G ...

  positive_projector = 0;
  for (int a = 0; a < dim; a++)
    {
      double lambda_a = eigenvalues[a];
      positive_projector += heaviside_function(lambda_a) * Q[a];
      for (int b = 0; b < dim; b++)
        if (b != a)
          {
            double lambda_b = eigenvalues[b];
            double v_ab = 0.0;
            if (std::fabs(lambda_a - lambda_b) > 1.0e-12)
              v_ab = (positive_ramp_function(lambda_a) - positive_ramp_function(lambda_b))
                   / (lambda_a - lambda_b);
            else
              v_ab = 0.5 * (heaviside_function(lambda_a) + heaviside_function(lambda_b));
            positive_projector += 0.5 * v_ab * 0.5 * (G[a][b] + G[b][a]);
          }
    }

  negative_projector = 0;
  for (int a = 0; a < dim; a++)
    {
      double lambda_a = eigenvalues[a];
      negative_projector += heaviside_function(-lambda_a) * Q[a];
      for (int b = 0; b < dim; b++)
        if (b != a)
          {
            double lambda_b = eigenvalues[b];
            double v_ab = 0.0;
            if (std::fabs(lambda_a - lambda_b) > 1.0e-12)
              v_ab = (negative_ramp_function(lambda_a) - negative_ramp_function(lambda_b))
                   / (lambda_a - lambda_b);
            else
              v_ab = 0.5 * (heaviside_function(-lambda_a) + heaviside_function(-lambda_b));
            negative_projector += 0.5 * v_ab * 0.5 * (G[a][b] + G[b][a]);
          }
    }
}
```

**公式（原文对应）：**
\[
\mathbb P^{+}=\frac{\partial\pmb\epsilon^+}{\partial\pmb\epsilon},
\quad
\mathbb P^{-}=\frac{\partial\pmb\epsilon^-}{\partial\pmb\epsilon}.
\]

---

### 2.1-11 一阶变分公式（7）

**对应最底层代码链接：**
- 残量单元装配（\[\delta\Pi\] 对应离散残量）  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3804-L3820
- 四点场更新（\[\boldsymbol\epsilon(\mathbf u),d,\nabla d,d_n\]）  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L1800-L1816

**代码：**
```cpp
// update quadrature-point state from current FE solution
scratch.m_fe_values[m_u_fe].get_function_symmetric_gradients(
  scratch.m_solution_UQPH, scratch.m_solution_symm_grads_u_cell);
scratch.m_fe_values[m_d_fe].get_function_values(
  scratch.m_solution_UQPH, scratch.m_solution_values_phasefield_cell);
scratch.m_fe_values[m_d_fe].get_function_gradients(
  scratch.m_solution_UQPH, scratch.m_solution_grad_phasefield_cell);
scratch.m_fe_values[m_d_fe].get_function_values(
  scratch.m_solution_previous_step, scratch.m_phasefield_previous_step_cell);

lqph[q_point]->update_field_values(scratch.m_solution_symm_grads_u_cell[q_point],
                                   scratch.m_solution_values_phasefield_cell[q_point],
                                   scratch.m_solution_grad_phasefield_cell[q_point],
                                   scratch.m_phasefield_previous_step_cell[q_point],
                                   scratch.m_delta_time);

// residual assembly
if (i_group == m_u_dof)
  {
    data.m_cell_rhs(i) += (symm_grad_N_disp[i] * cauchy_stress) * JxW;
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

**公式（原文对应）：**
\[
\delta \Pi(\pmb u,d)=D_{(\delta\pmb u,\delta d)}\Pi(\pmb u,d)=\cdots \quad (7)
\]

---

### 2.1-12 弱式 \[r_{\mathbf u}=0,\;r_d=0\]

**对应最底层代码链接：**
- 位移残量项 + 体力项  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3804-L3810
- 相场残量项  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3811-L3820

**代码：**
```cpp
if (i_group == m_u_dof)
  {
    data.m_cell_rhs(i) += (symm_grad_N_disp[i] * cauchy_stress) * JxW;
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

**公式（原文对应）：**
\[
\begin{cases}
r_{\pmb u}(\pmb u,d)=(\nabla^{(s)}\delta\pmb u,\pmb\sigma)-(\delta\pmb u,\pmb b)-(\delta\pmb u,\pmb t)_{\Gamma_t}=0,\\
r_d(\pmb u,d)=(\delta d,\frac{g_c}{l}d)+(\nabla\delta d,g_cl\nabla d)+(\delta d,g'(d)\psi^+(\pmb\epsilon))=0.
\end{cases}
\]

---

## 2.2 Finite element discretization

### 2.2-1 插值公式 \[\mathbf u=\mathbf N_{u_A}\mathbf u_A,\; d=N_{d_A}d_A\]

**对应最底层代码链接：**
- 单元形函数值提取（离散场重建）  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3758-L3770

**代码：**
```cpp
if (k_group == m_u_dof)
  {
    scratch.m_Nx_disp[q_point][k] =
      scratch.m_fe_values[m_u_fe].value(k, q_point);
    scratch.m_grad_Nx_disp[q_point][k] =
      scratch.m_fe_values[m_u_fe].gradient(k, q_point);
    scratch.m_symm_grad_Nx_disp[q_point][k] =
      symmetrize(scratch.m_grad_Nx_disp[q_point][k]);
  }
else if (k_group == m_d_dof)
  {
    scratch.m_Nx_phasefield[q_point][k] =
      scratch.m_fe_values[m_d_fe].value(k, q_point);
    scratch.m_grad_Nx_phasefield[q_point][k] =
      scratch.m_fe_values[m_d_fe].gradient(k, q_point);
  }
```

**公式（原文对应）：**
\[
\pmb u=\pmb N_{u_A}\pmb u_A,\quad d=N_{d_A}d_A.
\]

---

### 2.2-2 变分插值 \[\delta\mathbf u=\mathbf N_{u_A}\delta\mathbf u_A,\;\delta d=N_{d_A}\delta d_A\]

**对应最底层代码链接：**
- 测试函数层面通过 \[i\] 号形函数参与残量装配  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3800-L3820

**代码：**
```cpp
for (const unsigned int i : scratch.m_fe_values.dof_indices())
  {
    const unsigned int i_group = m_fe.system_to_base_index(i).first.first;

    if (i_group == m_u_dof)
      {
        data.m_cell_rhs(i) += (symm_grad_N_disp[i] * cauchy_stress) * JxW;
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
  }
```

**公式（原文对应）：**
\[
\delta\pmb u=\pmb N_{u_A}\delta\pmb u_A,\quad \delta d=N_{d_A}\delta d_A.
\]

---

### 2.2-3 公式（8）离散总能量 \[\Pi(\mathbf u_A,d_A)\]

**对应最底层代码链接：**
- 全局能量积分（按单元/积分点累加）  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L6117-L6139
- 其中材料点能量定义（应变能 + 裂纹耗散）  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L884-L891

**代码：**
```cpp
// material-point energies
m_strain_energy_total = degradation * m_strain_energy_positive + m_strain_energy_negative;

m_crack_energy_dissipation = m_gc * (  0.5 / m_length_scale * m_phase_field_value * m_phase_field_value
                                   + 0.5 * m_length_scale * m_grad_phasefield * m_grad_phasefield)
                                   // the term due to viscosity regularization
                                   + (m_phase_field_value - phase_field_value_previous_step)
                                   * (m_phase_field_value - phase_field_value_previous_step)
                                   * 0.5 * m_eta / delta_time;

// global integration
double energy_functional = 0.0;
for (const auto &cell : m_dof_handler.active_cell_iterators())
  {
    fe_values.reinit(cell);
    const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
      m_quadrature_point_history.get_data(cell);

    for (unsigned int q_point = 0; q_point < m_n_q_points; ++q_point)
      {
        const double JxW = fe_values.JxW(q_point);
        energy_functional += lqph[q_point]->get_total_strain_energy() * JxW;
        energy_functional += lqph[q_point]->get_crack_energy_dissipation() * JxW;
      }
  }
```

**公式（原文对应）：**
\[
\Pi(\pmb u_A,d_A)=\int_\Omega\psi(\pmb\epsilon(\pmb N_{u_A}\pmb u_A),N_{d_A}d_A)\,\mathrm d\Omega
+\int_\Omega\frac{g_c}{2l}\Big((N_{d_A}d_A)^2+l^2(\nabla N_{d_A}d_A)\cdot(\nabla N_{d_A}d_A)\Big)\,\mathrm d\Omega
-\int_\Omega\pmb b\cdot(\pmb N_{u_A}\pmb u_A)\,\mathrm d\Omega
-\int_{\partial\Omega}\pmb t\cdot(\pmb N_{u_A}\pmb u_A)\,\mathrm d\Gamma.\quad (8)
\]

---

### 2.2-4 公式（9）梯度（残量） \[\mathbf r=\nabla\Pi=(r_{\mathbf u},r_d)^\mathrm T\]

**对应最底层代码链接：**
- 位移残量 \[r_{\mathbf u_A}\] 与相场残量 \[r_{d_A}\] 的离散装配  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3804-L3820

**代码：**
```cpp
if (i_group == m_u_dof)
  {
    data.m_cell_rhs(i) += (symm_grad_N_disp[i] * cauchy_stress) * JxW;
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

**公式（原文对应）：**
\[
\begin{aligned}
\pmb r&=\nabla\Pi=(r_{\pmb u},r_d)^\mathrm T,\\
r_{\pmb u_A}&=\frac{\partial\Pi}{\partial \pmb u_A}
=(\nabla^{(s)}\pmb N_{u_A},\pmb\sigma)-(\pmb N_{u_A},\pmb b)-(\pmb N_{u_A},\pmb t)_{\Gamma_t},\\
r_{d_A}&=\frac{\partial\Pi}{\partial d_A}
=\left(N_{d_A},\frac{g_c}{l}d+g'(d)\psi^+\right)+(\nabla N_{d_A},g_cl\nabla d).
\end{aligned}
\quad (9)
\]

---

### 2.2-5 公式（10）Hessian 分块 \[\mathbf K\]

**对应最底层代码链接：**
- 该仓库实际用于求解的是近似 \[\hat{\mathbf K}\]（见公式（11）），其单元级块装配在：  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3931-L3951

**代码：**
```cpp
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
  }
```

> 说明： 论文公式（10） 给出完整 \[\mathbf K_{uu},\mathbf K_{ud},\mathbf K_{du},\mathbf K_{dd}\]。 本仓库在求解流程中有意采用块对角近似（公式（11）），因此在最底层装配中只保留 \[\mathbf K_{uu}\] 与 \[\mathbf K_{dd}\]，并有意省略耦合块 \[\mathbf K_{ud},\mathbf K_{du}\]。

**公式（原文对应）：**
\[
\mathbf K=\nabla^2\Pi=
\begin{bmatrix}
\mathbf K_{uu} & \mathbf K_{ud}\\
\mathbf K_{du} & \mathbf K_{dd}
\end{bmatrix},\quad
\begin{aligned}
\mathbf K_{u_Au_B}&=(\nabla^{(s)}\pmb N_{u_A},\frac{\partial\pmb\sigma}{\partial\pmb\epsilon}:\nabla^{(s)}\pmb N_{u_B}),\\
\mathbf K_{u_Ad_B}&=(\nabla^{(s)}\pmb N_{u_A},g'(d)\pmb\sigma^+N_{d_B}),\\
\mathbf K_{d_Au_B}&=(N_{d_A},g'(d)\pmb\sigma^+:\nabla^{(s)}\pmb N_{u_B}),\\
\mathbf K_{d_Ad_B}&=\left(N_{d_A},\left(\frac{g_c}{l}+g''(d)\psi^+\right)N_{d_B}\right)+(\nabla N_{d_A},g_cl\nabla N_{d_B}).
\end{aligned}
\quad (10)
\]

---

### 2.2-6 不可逆盒约束 \[ d_A^{(n)}\le d_A\le 1 \]

**对应最底层代码链接：**
- 直接设置相场增量上下界（Cauchy 点阶段）  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4859-L4870
- 活跃集判定（下界/上界/重合）  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5731-L5759

**代码：**
```cpp
if (i_group == m_d_dof)
  {
    // lower bound is old solution => increment lower bound is 0
    if (solution_delta_cauchy_point.block(m_d_dof)[b] < 0.0)
      solution_delta_cauchy_point.block(m_d_dof)[b] = 0.0;

    // upper bound is 1.0
    if (solution_delta_cauchy_point.block(m_d_dof)[b] > 1.0 - m_solution.block(m_d_dof)[b])
      solution_delta_cauchy_point.block(m_d_dof)[b] = 1.0 - m_solution.block(m_d_dof)[b];

    if (std::fabs(solution_delta_cauchy_point.block(m_d_dof)[b]) < 1.0e-15)
      m_active_set_phasefield[b] = 1; // lower bound
    else if (std::fabs(solution_delta_cauchy_point.block(m_d_dof)[b]
                       - (1.0 - m_solution.block(m_d_dof)[b])) < 1.0e-15)
      m_active_set_phasefield[b] = 2; // upper bound
  }

// active-set check against total trial phasefield
Vector<double> solution_phasefield_total(m_solution.block(m_d_dof));
solution_phasefield_total += solution_delta.block(m_d_dof);

if (std::fabs(solution_phasefield_total(i) - m_solution.block(m_d_dof)[i]) < 1.0e-15
    && std::fabs(solution_phasefield_total(i) - 1.0) < 1.0e-15)
  {
    m_active_set_phasefield(i) = 3; // lower bound overlaps with upper bound
  }
else if (std::fabs(solution_phasefield_total(i) - m_solution.block(m_d_dof)[i]) < 1.0e-15)
  {
    m_active_set_phasefield(i) = 1; // lower bound
  }
else if (std::fabs(solution_phasefield_total(i) - 1.0) < 1.0e-15)
  {
    m_active_set_phasefield(i) = 2; // upper bound
  }
```

**公式（原文对应）：**
\[
d_A^{(n)} \le d_A \le 1.
\]

---

### 2.2-7 公式（11）块对角矩阵 \[\hat{\mathbf K}=\mathrm{diag}(\mathbf K_{uu},\mathbf K_{dd})\]

**对应最底层代码链接：**
- 只组装 \[uu\] 与 \[dd\]，不组装 \[ud\]/\[du\] 的单元矩阵逻辑（即块对角）  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3931-L3954

**代码：**
```cpp
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
  }
else
  Assert((i_group <= m_d_dof) && (j_group <= m_d_dof), ExcInternalError());
```

**公式（原文对应）：**
\[
\hat{\mathbf K}=\begin{bmatrix}
\mathbf K_{uu} & \mathbf 0\\
\mathbf 0 & \mathbf K_{dd}
\end{bmatrix}. \quad (11)
\]

---

## 补充：第 2 章公式到“最低层计算”调用链（便于审计）

**链接：**
- FE 解 \[\to\] 积分点场变量更新：  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L1798-L1816
- 积分点材料更新入口：  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L931-L939
- 材料本构/能量核心：  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L820-L891
- 单元残量装配：  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3722-L3857
- 单元块对角 Hessian 装配：  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3860-L3959
- 全局能量计算：  
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L6117-L6168

**代码（调用链关键片段）：**
```cpp
// FE solution -> q-point fields
scratch.m_fe_values[m_u_fe].get_function_symmetric_gradients(...);
scratch.m_fe_values[m_d_fe].get_function_values(...);
scratch.m_fe_values[m_d_fe].get_function_gradients(...);

// q-point material update
lqph[q_point]->update_field_values(...);

// constitutive update inside material
usr_spectrum_decomposition::spectrum_decomposition<dim>(m_strain, eigenvalues, eigenvectors);
strain_positive = usr_spectrum_decomposition::positive_tensor(eigenvalues, eigenvectors);
strain_negative = usr_spectrum_decomposition::negative_tensor(eigenvalues, eigenvectors);
usr_spectrum_decomposition::positive_negative_projectors(...);
m_stress = degradation * stress_positive + stress_negative;
m_mechanical_C = degradation * C_positive + C_negative;

// residual and matrix assembly
data.m_cell_rhs(i) += ...;
data.m_cell_matrix(i, j) += ...;

// global energy
energy_functional += lqph[q_point]->get_total_strain_energy() * JxW;
energy_functional += lqph[q_point]->get_crack_energy_dissipation() * JxW;
```

