  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3940-L3951

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

**论文原公式（2.2, Eq.11）**

\[
\hat{\mathbf{K}} =
\left[
\begin{array}{cc}
\mathbf{K}_{uu} & \mathbf{0}\\
\mathbf{0} & \mathbf{K}_{dd}
\end{array}
\right].
\]

**解读**
- 代码装配阶段显式区分 `u-u` 与 `d-d` 两个对角块，和 Eq.(11) 的结构一致。
- 这类块对角结构在文中用于后续单调迭代、子空间最小化与预条件设计。

---

## 第2章符号级补充映射（避免遗漏）

> 下面补齐第2章里出现但容易被忽略的符号/量，仍保持“代码在前、公式在后”。

### A. 一阶不变量 \(I_1 = \mathrm{tr}(\epsilon)\)

**代码**

- 链接：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L847

```cpp
const double I_1 = trace(m_strain);
```

**公式**

\[
I_1 = \mathrm{tr}(\pmb\epsilon).
\]

**解读**

- 该标量是拉压判别（Ramp/Heaviside）和体积应力贡献的核心输入。

---

### B. \(g'(d)\) 与 \(g''(d)\)

**代码**

- 链接：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L249-L258

```cpp
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

**公式**

\[
g(d) = (1-d)^2,
\quad
g'(d) = 2(d-1),
\quad
g''(d)=2.
\]

**解读**

- 与相场残差中的 \(g'(d)\psi^+\) 项、刚度中的 \(g''(d)\psi^+\) 项直接对应。

---

### C. 裂纹表面密度项 \(\gamma(d,\nabla d)\) 的离散核

**代码**

- 链接：
  https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L886-L887

```cpp
m_crack_energy_dissipation = m_gc * (  0.5 / m_length_scale * m_phase_field_value * m_phase_field_value
                                   + 0.5 * m_length_scale * m_grad_phasefield * m_grad_phasefield)
```

**公式**

\[
\gamma(d,\nabla d) = \frac{1}{2l}\left(d^2 + l^2\nabla d\cdot\nabla d\right).
\]

**解读**

- 代码中 `0.5/l*d^2 + 0.5*l*|grad d|^2` 与上式逐项一致，再乘 `g_c` 形成裂纹耗散能密度。

---

## 结论：第2章“公式→代码”最短路径总结

1. **材料本构与拉压分裂**：`main.cc::LinearIsotropicElasticityAdditiveSplit::update_material_data` + `SpectrumDecomposition.*`。
2. **退化函数与导数**：`main.cc` 顶层 `degradation_function*` 三个函数。
3. **弱式残差/梯度**：`assemble_system_rhs_BFGS_one_cell`。
4. **Hessian与块结构**：`assemble_system_B0_one_cell`。
5. **能量泛函积分**：`calculate_energy_functional`。
6. **不可逆盒约束**：相场上下界初始化 + 迭代投影截断。

以上覆盖 `explain.md` 第2章全部出现的公式与关键符号，并给出了最底层可定位的实现代码与链接。
