# explain.md 第2章公式与最底层代码映射（代码在前、公式在后）

> 说明：
> 1. 本文档仅覆盖 `explain.md` 的第2章（`2.1` 与 `2.2`）中出现的全部公式与符号定义。
> 2. 每一小节严格按“代码在前、公式在后”的顺序。
> 3. 每个条目都给出可直接访问的 GitHub 链接（`blob/main` + 行号）。
> 4. 为保证可追溯性，代码片段直接展示，不省略关键计算行。

---

## 2.1 Phase-field formulation

### 2.1.1 公式(5)：应变能分裂与退化耦合

**代码（最底层计算）**

- 链接1（总应力与正负应变能耦合）：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L846-L885
- 链接2（退化函数定义）：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L244-L258

```cpp
// main.cc#L244-L258
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

```cpp
// main.cc#L846-L885
const double degradation = degradation_function(m_phase_field_value) + m_residual_k;
const double I_1 = trace(m_strain);

stress_positive = my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)
                                * Physics::Elasticity::StandardTensors<dim>::I
                + 2 * m_lame_mu * strain_positive;
stress_negative = my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)
                                * Physics::Elasticity::StandardTensors<dim>::I
                + 2 * m_lame_mu * strain_negative;

m_stress = degradation * stress_positive + stress_negative;

m_strain_energy_positive = 0.5 * my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)
                                               * usr_spectrum_decomposition::positive_ramp_function(I_1)
                         + m_lame_mu * strain_positive * strain_positive;

m_strain_energy_negative = 0.5 * my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)
                                               * usr_spectrum_decomposition::negative_ramp_function(I_1)
                         + m_lame_mu * strain_negative * strain_negative;

m_strain_energy_total = degradation * m_strain_energy_positive + m_strain_energy_negative;
```

**论文原公式（2.1, Eq.5）**

\[
\psi(\pmb{\epsilon}, d) = [g(d) + k] \, \psi^{+}(\pmb{\epsilon}) + \psi^{-}(\pmb{\epsilon}).
\]

**解读**

- 该式表达“仅拉伸能量退化、压缩能量不退化”的核心思想。
- 代码中 `degradation = degradation_function(d) + m_residual_k` 对应 `[g(d)+k]`，其中 `m_residual_k` 即文中小参数 `k`（防止全损伤时刚度退化为零导致病态）。
- `m_strain_energy_total = degradation * m_strain_energy_positive + m_strain_energy_negative` 是式(5)在积分点层面的直接离散实现。

---

### 2.1.2 公式(6)：退化函数

**代码（最底层计算）**

- 链接：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L244-L247

```cpp
// main.cc#L244-L247
double degradation_function(const double d)
{
  return (1.0 - d) * (1.0 - d);
}
```

**论文原公式（2.1, Eq.6）**

\[
g(d) = (1 - d)^2.
\]

**解读**

- `d=0` 时 `g(d)=1`，材料未退化；`d=1` 时 `g(d)=0`，拉伸承载能力降到最低。
- 代码逐字对应，不存在额外重标定。

---

### 2.1.3 符号定义：正负Ramp与Heaviside算子

**代码（最底层计算）**

- 链接：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.cc#L19-L38

```cpp
// SpectrumDecomposition.cc#L19-L38
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

**论文原公式（2.1，算子定义）**

\[
\langle x\rangle_{+} = \frac{1}{2}(x + |x|), \quad
\langle x\rangle_{-} = \frac{1}{2}(x - |x|), \quad
H(x) = \begin{cases}
1 & x \ge 0,\\
0 & x < 0.
\end{cases}
\]

**解读**

- `positive_ramp_function` 与 `negative_ramp_function` 等价于文中 \(\langle\cdot\rangle_+\)、\(\langle\cdot\rangle_-\) 的数值实现。
- `heaviside_function` 在 \(|x|<10^{-16}\) 时返回 `0.5`，这是常见的数值平滑处理，避免特征值重合或接近零时的不可导/不稳定。

---

### 2.1.4 谱分解公式：\(\epsilon = \sum_\alpha \epsilon_\alpha M_\alpha\)

**代码（最底层计算）**

- 链接1（谱分解函数）：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L29-L43
- 链接2（调用位置）：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L829-L833

```cpp
// SpectrumDecomposition.h#L29-L43
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
```

```cpp
// main.cc#L829-L833
Vector<double>              eigenvalues(dim);
std::vector<Tensor<1, dim>> eigenvectors(dim);
usr_spectrum_decomposition::spectrum_decomposition<dim>(m_strain,
                                                          eigenvalues,
                                                          eigenvectors);
```

**论文原公式（2.1，谱分解）**

\[
\pmb{\epsilon} = \sum_{\alpha} \epsilon_{\alpha} \mathbf{M}_{\alpha},
\quad
\mathbf{M}_{\alpha} = \pmb{n}_{\alpha} \otimes \pmb{n}_{\alpha}.
\]

**解读**

- `myEigenvalues[i]` 对应 \(\epsilon_\alpha\)，`myEigenvectors[i]` 对应 \(\pmb{n}_\alpha\)。
- 在后续 `positive_tensor/negative_tensor` 与 projector 计算里，\(\mathbf{M}_\alpha\) 由 `outer_product(eigenvectors[i], eigenvectors[i])` 构造。

---

### 2.1.5 正负应变张量：\(\epsilon^+\)、\(\epsilon^-\)

**代码（最底层计算）**

- 链接1（正负应变张量构造）：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L45-L69
- 链接2（调用位置）：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L835-L837

```cpp
// SpectrumDecomposition.h#L45-L69
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
```

```cpp
// main.cc#L835-L837
SymmetricTensor<2, dim> strain_positive, strain_negative;
strain_positive = usr_spectrum_decomposition::positive_tensor(eigenvalues, eigenvectors);
strain_negative = usr_spectrum_decomposition::negative_tensor(eigenvalues, eigenvectors);
```

**论文原公式（2.1）**

\[
\pmb{\epsilon}^{+} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{+} \mathbf{M}_{\alpha},
\quad
\pmb{\epsilon}^{-} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{-} \mathbf{M}_{\alpha}.
\]

**解读**

- 与论文完全同构：逐特征值施加正负Ramp，再以对应特征向量外积重组。
- 这是“拉压分裂”后续应力与切线构造的基础。

---

### 2.1.6 正负应变能：\(\psi^+\)、\(\psi^-\)

**代码（最底层计算）**

- 链接：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L876-L882

```cpp
// main.cc#L876-L882
m_strain_energy_positive = 0.5 * my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)
                                           * usr_spectrum_decomposition::positive_ramp_function(I_1)
                         + m_lame_mu * strain_positive * strain_positive;

m_strain_energy_negative = 0.5 * my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)
                                           * usr_spectrum_decomposition::negative_ramp_function(I_1)
                         + m_lame_mu * strain_negative * strain_negative;
```

**论文原公式（2.1）**

\[
\psi^{+}(\pmb{\epsilon}) = \frac{1}{2}\lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{+}^{2} + \mu \, \pmb{\epsilon}^{+} : \pmb{\epsilon}^{+},
\quad
\psi^{-}(\pmb{\epsilon}) = \frac{1}{2}\lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{-}^{2} + \mu \, \pmb{\epsilon}^{-} : \pmb{\epsilon}^{-}.
\]

**解读**

- `I_1 = trace(m_strain)` 对应 \(\mathrm{tr}\pmb{\epsilon}\)。
- `strain_positive * strain_positive` 对应双点积 \(\epsilon^+ : \epsilon^+\)（在 deal.II 的 `SymmetricTensor` 运算中即内积）。

---

### 2.1.7 应力分解：\(\sigma=[g(d)+k]\sigma^+ + \sigma^-\)

**代码（最底层计算）**

- 链接：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L857-L865

```cpp
// main.cc#L857-L865
stress_positive = my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)
                                * Physics::Elasticity::StandardTensors<dim>::I
                + 2 * m_lame_mu * strain_positive;
stress_negative = my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)
                                * Physics::Elasticity::StandardTensors<dim>::I
                + 2 * m_lame_mu * strain_negative;

m_stress = degradation * stress_positive + stress_negative;
m_stress_positive = stress_positive;
```

**论文原公式（2.1）**

\[
\pmb{\sigma} = \frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial \pmb{\epsilon}}
= [g(d)+k]\frac{\partial\psi^{+}(\pmb{\epsilon})}{\partial \pmb{\epsilon}} + \frac{\partial\psi^{-}(\pmb{\epsilon})}{\partial \pmb{\epsilon}}
= [g(d)+k]\pmb{\sigma}^{+} + \pmb{\sigma}^{-}.
\]

以及

\[
\pmb{\sigma}^{+} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{+}\mathbf{I} + 2\mu \pmb{\epsilon}^{+},
\quad
\pmb{\sigma}^{-} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{-}\mathbf{I} + 2\mu \pmb{\epsilon}^{-}.
\]

**解读**

- `stress_positive`、`stress_negative` 分别对应论文定义的 \(\sigma^+\)、\(\sigma^-\)。
- `m_stress` 是最终用于平衡方程残差组装的 Cauchy 应力。

---

### 2.1.8 切线模量：\(\partial\sigma/\partial\epsilon\)

**代码（最底层计算）**

- 链接1（projector构造）：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L71-L144
- 链接2（切线模量组合）：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L867-L874

```cpp
// SpectrumDecomposition.h#L71-L144（完整核心）
template <int dim>
void positive_negative_projectors(Vector<double> const & eigenvalues,
                                  std::vector<Tensor<1, dim>> const & eigenvectors,
                                  SymmetricTensor<4, dim> & positive_projector,
                                  SymmetricTensor<4, dim> & negative_projector)
{
  Assert(dim <= 3,
         ExcMessage("Project tensors only work for dim <= 3."));

  std::array<SymmetricTensor<2, dim>, dim> M;
  for (int a = 0; a < dim; a++)
    M[a] = symmetrize(outer_product(eigenvectors[a], eigenvectors[a]));

  std::array<SymmetricTensor<4, dim>, dim> Q;
  for (int a = 0; a < dim; a++)
    Q[a] = outer_product(M[a], M[a]);

  std::array<std::array<SymmetricTensor<4, dim>, dim>, dim> G;
  for (int a = 0; a < dim; a++)
    for (int b = 0; b < dim; b++)
      for (int i = 0; i < dim; i++)
        for (int j = 0; j < dim; j++)
          for (int k = 0; k < dim; k++)
            for (int l = 0; l < dim; l++)
              G[a][b][i][j][k][l] = M[a][i][k] * M[b][j][l]
                                   + M[a][i][l] * M[b][j][k];

  positive_projector = 0;
  for (int a = 0; a < dim; a++)
    {
      double lambda_a = eigenvalues[a];
      positive_projector += heaviside_function(lambda_a)
                          * Q[a];
      for (int b = 0; b < dim; b++)
        {
          if (b != a)
            {
              double lambda_b = eigenvalues[b];
              double v_ab = 0.0;
              if (std::fabs(lambda_a - lambda_b) > 1.0e-12)
                v_ab = (positive_ramp_function(lambda_a) - positive_ramp_function(lambda_b))
                     / (lambda_a - lambda_b);
              else
                v_ab = 0.5 * (  heaviside_function(lambda_a)
                              + heaviside_function(lambda_b) );
              positive_projector += 0.5 * v_ab * 0.5 * (G[a][b] + G[b][a]);
            }
        }
    }

  negative_projector = 0;
  for (int a = 0; a < dim; a++)
    {
      double lambda_a = eigenvalues[a];
      negative_projector += heaviside_function(-lambda_a)
                          * Q[a];
      for (int b = 0; b < dim; b++)
        {
          if (b != a)
            {
              double lambda_b = eigenvalues[b];
              double v_ab = 0.0;
              if (std::fabs(lambda_a - lambda_b) > 1.0e-12)
                v_ab = (negative_ramp_function(lambda_a) - negative_ramp_function(lambda_b))
                     / (lambda_a - lambda_b);
              else
                v_ab = 0.5 * (  heaviside_function(-lambda_a)
                              + heaviside_function(-lambda_b) );
              negative_projector += 0.5 * v_ab * 0.5 * (G[a][b] + G[b][a]);
            }
        }
    }
}
```

```cpp
// main.cc#L867-L874
SymmetricTensor<4, dim> C_positive, C_negative;
C_positive = my_lambda * usr_spectrum_decomposition::heaviside_function(I_1)
                           * Physics::Elasticity::StandardTensors<dim>::IxI
           + 2 * m_lame_mu * projector_positive;
C_negative = my_lambda * usr_spectrum_decomposition::heaviside_function(-I_1)
                           * Physics::Elasticity::StandardTensors<dim>::IxI
           + 2 * m_lame_mu * projector_negative;
m_mechanical_C = degradation * C_positive + C_negative;
```

**论文原公式（2.1）**

\[
\frac{\partial \pmb{\sigma}}{\partial \pmb{\epsilon}}
= [g(d)+k]\frac{\partial \pmb{\sigma}^{+}}{\partial \pmb{\epsilon}} + \frac{\partial \pmb{\sigma}^{-}}{\partial \pmb{\epsilon}}
= [g(d)+k]\left[\lambda H(\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes\mathbf{I} + 2\mu \mathbb{P}^{+}\right]
+ \left[\lambda H(-\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes\mathbf{I} + 2\mu \mathbb{P}^{-}\right].
\]

及

\[
\mathbb{P}^{+}=\frac{\partial \pmb{\epsilon}^{+}}{\partial\pmb{\epsilon}},
\quad
\mathbb{P}^{-}=\frac{\partial \pmb{\epsilon}^{-}}{\partial\pmb{\epsilon}}.
\]

**解读**

- `projector_positive/projector_negative` 是 \(\mathbb{P}^{+},\mathbb{P}^{-}\) 的直接数值实现。
- `m_mechanical_C` 是单元组装中 \(\mathbf{K}_{uu}\) 的核心材料切线，后续通过 `symm_grad_N * C * symm_grad_N` 进入刚度矩阵。

---

### 2.1.9 公式(7)：总能量一阶变分

**代码（最底层计算）**

- 链接（残差组装，直接对应 \(\delta\Pi\) 的两条方程）：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3804-L3820

```cpp
// main.cc#L3804-L3820
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

以及 Neumann 边界项：

- 链接：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3831-L3853

```cpp
// main.cc#L3831-L3853
for (const auto &face : cell->face_iterators())
  if (face->at_boundary() && face->boundary_id() == face_pressure_id)
    {
      scratch.m_fe_face_values.reinit(cell, face);

      for (const unsigned int f_q_point : scratch.m_fe_face_values.quadrature_point_indices())
        {
          const Tensor<1, dim> &N = scratch.m_fe_face_values.normal_vector(f_q_point);
          const double         pressure  = p0 * time_ramp;
          const Tensor<1, dim> traction  = pressure * N;

          for (const unsigned int i : scratch.m_fe_values.dof_indices())
            {
              const unsigned int i_group = m_fe.system_to_base_index(i).first.first;
              if (i_group == m_u_dof)
                {
                  const unsigned int component_i = m_fe.system_to_component_index(i).first;
                  const double Ni = scratch.m_fe_face_values.shape_value(i, f_q_point);
                  const double JxW = scratch.m_fe_face_values.JxW(f_q_point);
                  data.m_cell_rhs(i) -= (Ni * traction[component_i]) * JxW;
                }
            }
        }
    }
```

**论文原公式（2.1, Eq.7）**

\[
\begin{array}{rl}
\delta \Pi (\pmb {u},d)
&= D_{(\delta \pmb {u},\delta d)}\Pi (\pmb {u},d)
= \left.\frac{\mathrm{d}}{\mathrm{d}\epsilon}\right|_{\epsilon = 0}\Pi (\pmb {u} + \epsilon \delta \pmb {u},d + \epsilon \delta d) \\
&= \int_{\Omega}\left(\frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial \pmb{\epsilon}}:\pmb{\epsilon}(\delta \pmb {u}) + \frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial d}\delta d\right)\mathrm{d}\Omega \\
&\quad +\int_{\Omega}\frac{g_{c}}{l}\left(d\delta d + l^{2}\nabla d\cdot \nabla \delta d\right)\mathrm{d}\Omega
-\int_{\Omega}\pmb {b}\cdot \delta \pmb {u}\mathrm{d}\Omega
-\int_{\partial \Omega}\pmb {t}\cdot \delta \pmb {u}\mathrm{d}\Gamma \\
&= (\nabla^{(s)}\delta \pmb {u},\pmb {\sigma}) - (\delta \pmb {u},\pmb {b}) - (\delta \pmb {u},\pmb{t})_{\Gamma_{t}}
+ (\delta d,\frac{g_{c}}{l} d)
+ (\nabla \delta d,g_{c}l\nabla d)
+ (\delta d,g^{\prime}(d)\psi^{+}(\pmb{\epsilon})).
\end{array}
\]

**解读**

- 代码中位移方程、相场方程的单元右端累加项，正是 \(\delta\Pi\) 对应的离散后残差。
- 代码还加入了粘性正则项 `eta/delta_time*(d-d_old)`，这是对文中基本模型的扩展（有助于数值稳定）。

---

### 2.1.10 弱式残差：\(r_u=0\), \(r_d=0\)

**代码（最底层计算）**

- 链接：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3804-L3820

```cpp
// main.cc#L3804-L3820（同上）
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

**论文原公式（2.1，弱式）**

\[
\left\{
\begin{array}{ll}
r_{\pmb{u}}(\pmb{u},d)
= (\nabla^{(s)}\delta \pmb{u},\pmb{\sigma}) - (\delta \pmb{u},\pmb{b}) - (\delta \pmb{u},\pmb{t})_{\Gamma_{t}} = 0, \\
r_{d}(\pmb{u},d)
= (\delta d,\frac{g_{c}}{l} d) + (\nabla \delta d,g_{c}l\nabla d) + (\delta d,g^{\prime}(d)\psi^{+}(\pmb{\epsilon})) = 0.
\end{array}
\right.
\]

**解读**

- `i_group == m_u_dof` 的分支就是 \(r_u\)；`i_group == m_d_dof` 分支就是 \(r_d\)。
- 单元贡献通过 `JxW` 积分权重累计到全局右端，构成离散非线性方程组残差。

---

## 2.2 Finite element discretization

### 2.2.1 插值表达：\(\mathbf{u}=\mathbf{N}_{u_A}\mathbf{u}_A\), \(d=N_{d_A}d_A\)

**代码（最底层计算）**

- 链接：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3756-L3770

```cpp
// main.cc#L3756-L3770
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

**论文原公式（2.2）**

\[
\pmb {u} = \pmb{N}_{u_{A}}\pmb{u}_{A},
\quad
\text{and}
\quad
d = N_{d_{A}}d_{A}.
\]

**解读**

- 代码从 `FEValues::value/gradient` 获取形函数值与梯度，等价实现离散插值。
- 在单元层，`scratch.m_Nx_disp`、`scratch.m_Nx_phasefield` 就是上式中的 \(\pmb N_{u_A}\) 与 \(N_{d_A}\)。

---

### 2.2.2 变分插值：\(\delta\mathbf{u}\), \(\delta d\)

**代码（最底层计算）**

- 链接（测试函数在装配中的使用）：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3804-L3819

```cpp
// main.cc#L3804-L3819
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

**论文原公式（2.2）**

\[
\delta \pmb {u} = \pmb{N}_{u_{A}}\delta \pmb{u}_{A},
\quad
\text{and}
\quad
\delta d = N_{d_{A}}\delta d_{A}.
\]

**解读**

- 弱式离散中，`N_disp[i]`、`symm_grad_N_disp[i]`、`N_phasefield[i]`、`grad_N_phasefield[i]` 即试函数/检验函数离散基。
- 残差每一项都由这些形函数乘物理量积分得到。

---

### 2.2.3 公式(8)：离散总能量泛函

**代码（最底层计算）**

- 链接1（积分点总能量定义）：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L884-L891
- 链接2（域积分累加为全局能量）：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L6117-L6139

```cpp
// main.cc#L884-L891
m_strain_energy_total = degradation * m_strain_energy_positive + m_strain_energy_negative;

m_crack_energy_dissipation = m_gc * (  0.5 / m_length_scale * m_phase_field_value * m_phase_field_value
                                   + 0.5 * m_length_scale * m_grad_phasefield * m_grad_phasefield)
                                   // the term due to viscosity regularization
                                   + (m_phase_field_value - phase_field_value_previous_step)
                                   * (m_phase_field_value - phase_field_value_previous_step)
                                   * 0.5 * m_eta / delta_time;
```

```cpp
// main.cc#L6117-L6139
double PhaseFieldMonolithicSolve<dim>::calculate_energy_functional() const
{
  double energy_functional = 0.0;

  FEValues<dim> fe_values(m_fe, m_qf_cell, update_JxW_values);

  for (const auto &cell : m_dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);

      const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
        m_quadrature_point_history.get_data(cell);
      Assert(lqph.size() == m_n_q_points, ExcInternalError());

      for (unsigned int q_point = 0; q_point < m_n_q_points; ++q_point)
        {
          const double JxW = fe_values.JxW(q_point);
          energy_functional += lqph[q_point]->get_total_strain_energy() * JxW;
          energy_functional += lqph[q_point]->get_crack_energy_dissipation() * JxW;
        }
    }

  return energy_functional;
}
```

**论文原公式（2.2, Eq.8）**

\[
\begin{array}{rl}
\Pi (\pmb {u}_{A},d_{A}) = &\int_{\Omega}\psi \left(\pmb{\epsilon}(\pmb {N}_{u_{A}}\pmb {u}_{A}),N_{d_{A}}d_{A}\right)\mathrm{d}\Omega \\
&+\int_{\Omega}\frac{g_{c}}{2l}\left((N_{d_{A}}d_{A})^{2} + l^{2}(\nabla N_{d_{A}}d_{A})\cdot (\nabla N_{d_{A}}d_{A})\right)\mathrm{d}\Omega \\
&-\int_{\Omega}\pmb {b}\cdot (\pmb {N}_{u_{A}}\pmb {u}_{A})\mathrm{d}\Omega
-\int_{\partial \Omega}\pmb {t}\cdot (\pmb {N}_{u_{A}}\pmb {u}_{A})\mathrm{d}\Gamma.
\end{array}
\]

**解读**

- `get_total_strain_energy()` 对应 \(\psi\) 项（已包含拉压分裂与退化）。
- `get_crack_energy_dissipation()` 对应 \(g_c/(2l)(d^2+l^2|\nabla d|^2)\) 裂纹表面近似能；代码另加粘性项。
- 体力与表面力在残差里显式出现（见 2.1.9），与能量泛函外载项一致。

---

### 2.2.4 公式(9)：离散梯度/残差

**代码（最底层计算）**

- 链接：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3804-L3820

```cpp
// main.cc#L3804-L3820
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

**论文原公式（2.2, Eq.9）**

\[
\begin{array}{rl}
\pmb{r} = \nabla \Pi = (r_{\pmb{u}},r_{d})^{\mathrm{T}},\\
r_{\pmb{u}_{A}} = \frac{\partial\Pi}{\partial \pmb{u}_{A}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},\pmb{\sigma}\right) - \left(\pmb{N}_{u_{A}},\pmb{b}\right) - \left(\pmb{N}_{u_{A}},\pmb{t}\right)_{\Gamma_{t}},\\
r_{d_{A}} = \frac{\partial\Pi}{\partial d_{A}} = \left(N_{d_{A}},\frac{g_{c}}{l} d + g^{\prime}(d)\psi^{+}\right) + (\nabla N_{d_{A}},g_{c}l\nabla d).
\end{array}
\]

**解读**

- 上述代码分量级实现了 \(r_{u_A}\) 与 \(r_{d_A}\) 的每一积分项。
- `degradation_function_derivative(phasefield_value) * current_positive_strain_energy` 即 \(g'(d)\psi^+\) 的离散值。

---

### 2.2.5 Hessian分块结构：\(\mathbf{K}=\nabla^2\Pi\)

**代码（最底层计算）**

- 链接（单元矩阵装配主循环）：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3927-L3955

```cpp
// main.cc#L3927-L3955
for (const unsigned int i : scratch.m_fe_values.dof_indices())
  {
    const unsigned int i_group = m_fe.system_to_base_index(i).first.first;

    if (i_group == m_u_dof)
      {
        symm_grad_Nx_i_x_C = symm_grad_N_disp[i] * mechanical_C;
      }

    for (const unsigned int j : scratch.m_fe_values.dof_indices())
      {
        const unsigned int j_group = m_fe.system_to_base_index(j).first.first;

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
          Assert((i_group <= m_d_dof) && (j_group <= m_d_dof),
                 ExcInternalError());
      } // j
  }  // i
```

**论文原公式（2.2）**

\[
\mathbf{K} = \nabla^{2}\Pi =
\left[
\begin{array}{ll}
\mathbf{K}_{uu} & \mathbf{K}_{ud}\\
\mathbf{K}_{du} & \mathbf{K}_{dd}
\end{array}
\right].
\]

**解读**

- 该段代码明确按 `u` 与 `d` 自由度分组组装刚度，体现论文的分块结构。
- 当前实现主要显式构造了对角块贡献（`uu` 与 `dd`），并通过其它流程与近似策略服务于BFGS/L-BFGS框架。

---

### 2.2.6 公式(10)：各分块显式表达

**代码（最底层计算）**

- 链接1（\(K_{uu}\), \(K_{dd}\)）：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3940-L3951
- 链接2（机械切线 \(\partial\sigma/\partial\epsilon\) 输入）：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3919-L3934

```cpp
// main.cc#L3919-L3934
const SymmetricTensor<4, dim> & mechanical_C  = lqph[q_point]->get_mechanical_C();
const std::vector<SymmetricTensor<2, dim>> & symm_grad_N_disp =
  scratch.m_symm_grad_Nx_disp[q_point];

SymmetricTensor<2, dim> symm_grad_Nx_i_x_C;

for (const unsigned int i : scratch.m_fe_values.dof_indices())
  {
    const unsigned int i_group = m_fe.system_to_base_index(i).first.first;

    if (i_group == m_u_dof)
      {
        symm_grad_Nx_i_x_C = symm_grad_N_disp[i] * mechanical_C;
      }
```

```cpp
// main.cc#L3940-L3951
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

**论文原公式（2.2, Eq.10）**

\[
\begin{array}{rl}
\mathbf{K}_{u_{A}u_{B}} &= \left(\nabla^{(s)}\pmb{N}_{u_{A}},\frac{\partial \pmb{\sigma}}{\partial \pmb{\epsilon}}:\nabla^{(s)}\pmb{N}_{u_{B}}\right), \qquad
\mathbf{K}_{u_{A}d_{B}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},g^{\prime}(d)\pmb{\sigma}^{+}N_{d_{B}}\right),\\
\mathbf{K}_{d_{A}u_{B}} &= \left(N_{d_{A}},g^{\prime}(d)\pmb{\sigma}^{+}:\nabla^{(s)}\pmb{N}_{u_{B}}\right), \qquad
\mathbf{K}_{d_{A}d_{B}} = \left(N_{d_{A}},\left(\frac{g_{c}}{l} +g^{\prime \prime}(d)\psi^{+}\right)N_{d_{B}}\right)
+ \left(\nabla N_{d_{A}},g_{c}l\nabla N_{d_{B}}\right).
\end{array}
\]

**解读**

- `uu` 块：`symm_grad_N * mechanical_C * symm_grad_N` 直接对应 \(\mathbf{K}_{u_Au_B}\)。
- `dd` 块：`gc/l + g''(d)psi+` 与梯度项 `gc*l*gradN_i*gradN_j` 一一对应 \(\mathbf{K}_{d_Ad_B}\)。
- 论文中的 `ud/du` 混合块在该实现中并未在此装配函数显式形成完整稀疏块，而是通过求解策略与近似 Hessian 结构处理（与后续章节的 BFGS/\(\hat{K}\) 思路一致）。

---

### 2.2.7 节点盒约束：\(d_A^{(n)} \le d_A \le 1\)

**代码（最底层计算）**

- 链接（下界来自上一步、上界为1.0的投影实现）：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L1379-L1394

```cpp
// main.cc#L1379-L1394
void PhaseFieldMonolithicSolve<dim>::point_projection(BlockVector<double> & solution_delta)
{
  // Phase-field value cannot exceed 1.0
  const double upper_limit = 1.0;

  Vector<double> solution_phasefield_total(m_solution.block(m_d_dof));
  solution_phasefield_total += solution_delta.block(m_d_dof);

  for (unsigned int i = 0; i < solution_phasefield_total.size(); ++i)
    {
      if (solution_delta.block(m_d_dof)[i] < 0.0)
        solution_delta.block(m_d_dof)[i] = 0.0;

      if (solution_phasefield_total[i] > upper_limit)
        solution_delta.block(m_d_dof)[i] = upper_limit - m_solution.block(m_d_dof)[i];
    }
}
```

- 链接（迭代后依据边界更新活动集状态）：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5731-L5764

```cpp
// main.cc#L5731-L5764（活动集判定：下界=上一时刻值， 上界=1）
Vector<double> solution_phasefield_total(m_solution.block(m_d_dof));
solution_phasefield_total += solution_delta.block(m_d_dof);

for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
  {
    if (   solution_delta.block(m_d_dof)[i] == 0.0
        && solution_phasefield_total[i] == 1.0)
      {
        m_active_set_phasefield(i) = 3; //lower bound overlaps with upper bound
      }
    else if (   solution_delta.block(m_d_dof)[i] == 0.0
             && solution_phasefield_total[i] != 1.0)
      {
        m_active_set_phasefield(i) = 1; //lower bound
      }
    else if (   solution_phasefield_total[i] == 1.0
             && solution_delta.block(m_d_dof)[i] != 0.0)
      {
        m_active_set_phasefield(i) = 2; //upper bound
      }
    else
      {
        m_active_set_phasefield(i) = 0;
      }
  }
```

**论文原公式（2.2）**

\[
d_A^{(n)}\leq d_A\leq 1.
\]

**解读**

- 下界取上一步解 `m_solution_previous_step`，体现不可逆条件；上界固定为 `1.0`。
- 该约束随后在梯度投影与线搜索路径中持续生效。

---

### 2.2.8 公式(11)：对角块矩阵 \(\hat{\mathbf K}\)

**代码（最底层计算）**

- 链接（只保留 `uu` 与 `dd` 组装逻辑）：
