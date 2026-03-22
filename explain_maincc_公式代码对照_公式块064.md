# explain.md 与 main.cc 公式-代码对应（公式块 064）

- 所属章节：`3.3. Generalized Cauchy point`
- explain.md 行号：`473`

论文公式：

\[
\begin{array}{rl} 
p(t) &= m(\mathbf{x}(t)) = \Pi + \mathbf{r}^{\mathrm{T}}(\mathbf{x}(t) - \mathbf{x}^0) + \frac{1}{2} (\mathbf{x}(t) - \mathbf{x}^0)^{\mathrm{T}}\mathbf{B}(\mathbf{x}(t) - \mathbf{x}^0) \\
&= \Pi + \mathbf{r}^{\mathrm{T}}(\mathbf{z}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)}) + \frac{1}{2} (\mathbf{z}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)})^{\mathrm{T}}\mathbf{B}(\mathbf{z}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)}) \\
&= \left(\Pi + \mathbf{r}^{\mathrm{T}}\mathbf{z}^{(j-1)} + \frac{1}{2} \mathbf{z}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j-1)}\right) \\
&\quad + \left(\mathbf{r}^{\mathrm{T}}\mathbf{d}^{(j-1)} + \mathbf{d}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j-1)}\right)\Delta t + \frac{1}{2}\left(\mathbf{d}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j-1)}\right)\Delta t^{2} \\
&= f_{j-1} + f_{j-1}^{\prime}\Delta t + \frac{1}{2} f_{j-1}^{\prime \prime}\Delta t^{2} = \hat{p} (\Delta t), 
\end{array} \quad (22)
\]

对应 `main.cc` 代码：

```cpp
// main.cc:4819-4835
    double f_prime = -(gradient_d * gradient_d);

    // M * p
    Vector<double> Mp(2 * list_size);
    if (list_size > 0)
      M_matrix.vmult(Mp, p);

    // B_0 * d
    BlockVector<double> B0_grandient_d(m_dofs_per_block);
    B0_matrix.vmult(B0_grandient_d, gradient_d);

    double f_prime_prime = gradient_d * B0_grandient_d;
    if (list_size > 0)
      f_prime_prime -= (p * Mp);

    double delta_t_min = -f_prime / f_prime_prime;


// main.cc:4891-4917
	f_prime += delta_t * f_prime_prime
	         + gradient_g.block(m_d_dof)[b] * gradient_g.block(m_d_dof)[b]
	         + temp_scalar * gradient_g.block(m_d_dof)[b];

	if (list_size > 0)
	  f_prime -= (w_b_T_x_M * c) * gradient_g.block(m_d_dof)[b];

	temp_scalar = ebT_x_B0_x_v(b, B0_matrix, gradient_d);

	f_prime_prime += 2.0 * gradient_g.block(m_d_dof)[b] * temp_scalar
	               + gradient_g.block(m_d_dof)[b] * gradient_g.block(m_d_dof)[b]
		             * B0_matrix.block(m_d_dof, m_d_dof)(b, b);

	if (list_size > 0)
	  {
	    f_prime_prime -= 2.0 * gradient_g.block(m_d_dof)[b] * (w_b_T_x_M * p);
	    f_prime_prime -= gradient_g.block(m_d_dof)[b] * gradient_g.block(m_d_dof)[b] * (w_b_T_x_M * w_b);
	  }

	// p_{j} = p_{j-1} + g_b * w_b;
	if (list_size > 0)
	  p.sadd(1.0, gradient_g.block(m_d_dof)[b], w_b);

	gradient_d.block(m_d_dof)[b] = 0.0;

	delta_t_min = -f_prime / f_prime_prime;
```

简要说明：对应变量或算子在 main.cc 中有同名或等价离散实现。
