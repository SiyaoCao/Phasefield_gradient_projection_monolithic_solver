# explain.md 与 main.cc 公式-代码对应（公式块 041）

- 所属章节：`3.1. Algorithm overview`
- explain.md 行号：`327`

论文公式：

\[
\left|\pmb{r}_{k + 1}^\mathrm{T}\pmb {p}_k\right| = \left|\pmb {r}(\pmb{x}_k + \alpha_k\pmb {p}_k)^\mathrm{T}\pmb {p}_k\right| \leq c_2\left|\pmb{r}_k^\mathrm{T}\pmb {p}_k\right|.
\]

对应 `main.cc` 代码：

```cpp
// main.cc:4196-4253
  double PhaseFieldMonolithicSolve<dim>::line_search_stepsize_strong_wolfe(const double phi_0,
				                                           const double phi_0_prime,
				                                           const BlockVector<double> & BFGS_p_vector,
				                                           const BlockVector<double> & solution_delta)
  {
    //AssertThrow(phi_0_prime < 0,
    //            ExcMessage("The derivative of phi at alpha = 0 should be negative!"));

    // Some line search parameters
    const double c1 = 0.0001;
    const double c2 = 0.9;
    const double alpha_max = 1.0;
    const unsigned int max_iter = 20;
    double alpha = 1.0;

    double phi_old = phi_0;
    double phi_prime_old = phi_0_prime;
    double alpha_old = 0.0;

    double phi, phi_prime;

    std::pair<double, double> current_phi_phi_prime;

    unsigned int i = 0;
    for (; i < max_iter; ++i)
      {
	current_phi_phi_prime = calculate_phi_and_phi_prime(alpha, BFGS_p_vector, solution_delta);
	phi = current_phi_phi_prime.first;
	phi_prime = current_phi_phi_prime.second;

	if (   ( phi > (phi_0 + c1 * alpha * phi_0_prime) )
	    || ( i > 0 && phi > phi_old ) )
	  {
	    return line_search_zoom_strong_wolfe(phi_old, phi_prime_old, alpha_old,
						 phi,     phi_prime,     alpha,
						 phi_0,   phi_0_prime,   BFGS_p_vector,
						 c1,      c2,            max_iter, solution_delta);
	  }

	if (std::fabs(phi_prime) <= c2 * std::fabs(phi_0_prime))
	  {
	    return alpha;
	  }

	if (phi_prime >= 0)
	  {
	    return line_search_zoom_strong_wolfe(phi,     phi_prime,     alpha,
						 phi_old, phi_prime_old, alpha_old,
						 phi_0,   phi_0_prime,   BFGS_p_vector,
						 c1,      c2,            max_iter, solution_delta);
	  }

	phi_old = phi;
	phi_prime_old = phi_prime;
	alpha_old = alpha;

	alpha = std::min(0.6*alpha, alpha_max);
      }

// main.cc:4260-4311
    line_search_zoom_strong_wolfe(double phi_low, double phi_low_prime, double alpha_low,
				  double phi_high, double phi_high_prime, double alpha_high,
				  double phi_0, double phi_0_prime, const BlockVector<double> & BFGS_p_vector,
				  double c1, double c2, unsigned int max_iter, const BlockVector<double> & solution_delta)
  {
    double alpha = 0;
    std::pair<double, double> current_phi_phi_prime;
    double phi, phi_prime;

    unsigned int i = 0;
    for (; i < max_iter; ++i)
      {
	// a simple bisection is faster than cubic interpolation
	alpha = 0.5 * (alpha_low + alpha_high);
	//alpha = line_search_interpolation_cubic(alpha_low, phi_low, phi_low_prime,
	//					alpha_high, phi_high, phi_high_prime);
	current_phi_phi_prime = calculate_phi_and_phi_prime(alpha, BFGS_p_vector, solution_delta);
	phi = current_phi_phi_prime.first;
	phi_prime = current_phi_phi_prime.second;

	if (   (phi > phi_0 + c1 * alpha * phi_0_prime)
	    || (phi > phi_low) )
	  {
	    alpha_high = alpha;
	    phi_high = phi;
	    phi_high_prime = phi_prime;
	  }
	else
	  {
	    if (std::fabs(phi_prime) <= c2 * std::fabs(phi_0_prime))
	      {
		return alpha;
	      }

	    if (phi_prime * (alpha_high - alpha_low) >= 0.0)
	      {
		alpha_high = alpha_low;
		phi_high_prime = phi_low_prime;
		phi_high = phi_low;
	      }

	    alpha_low = alpha;
	    phi_low_prime = phi_prime;
	    phi_low = phi;
	  }
      }

    // avoid unused variable warnings from compiler
    (void)phi_high;
    (void)phi_high_prime;
    return alpha;
  }
```

简要说明：对应线搜索（Armijo + 曲率条件）的数值实现。
