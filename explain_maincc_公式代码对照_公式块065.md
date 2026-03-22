# explain.md 与 main.cc 公式-代码对应（公式块 065）

- 所属章节：`3.3. Generalized Cauchy point`
- explain.md 行号：`487`

论文公式：

\[
f_{j-1}^{\prime \prime} = \mathbf{d}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j-1)} > 0.
\]

对应 `main.cc` 代码：

```cpp
// main.cc:4830-4835
    double f_prime_prime = gradient_d * B0_grandient_d;
    if (list_size > 0)
      f_prime_prime -= (p * Mp);

    double delta_t_min = -f_prime / f_prime_prime;


// main.cc:4916-4924
	delta_t_min = -f_prime / f_prime_prime;

	t_old = t;

	top_pair = t_series.top();
	t = top_pair.first;
	b = top_pair.second;

	delta_t = t - t_old;
```

简要说明：对应变量或算子在 main.cc 中有同名或等价离散实现。
