# Formula-to-Code Mapping (Part 2): Phase-Field Fracture Monolithic Solver

Continuation of `formula_code_mapping_1.md`.
This file covers Sections 3.3–3.4 (Cauchy point and subspace minimization), Section 3.5
(line-search Wolfe conditions, L-BFGS two-loop recursion), Eq. (33) (crack dissipation energy),
and the torsion boundary conditions for Section 4 numerical example.

---

## 3.3. Computation of the Cauchy Point

### Overall Function Signature

[`main.cc#L4773-L4782`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4773)

```cpp
  template <int dim>
  void PhaseFieldMonolithicSolve<dim>::
  calculate_cauchy_point(const BlockSparseMatrix<double> & B0_matrix,
	                 const std::list<BlockVector<double>> & y_vector_list,
		         const std::list<BlockVector<double>> & b0xs_vector_list,
			 const FullMatrix<double> & M_matrix,
			 const BlockVector<double> & gradient_g,
			 const BlockVector<double> & solution_delta,
			 BlockVector<double> & solution_delta_cauchy_point)
```

The Cauchy point is the minimizer of the quadratic model
[ m_k(\mathbf{x}) = \Pi_k + \mathbf{r}_k^T(\mathbf{x}-\mathbf{x}_k) + \tfrac{1}{2}(\mathbf{x}-\mathbf{x}_k)^T\mathbf{B}_k(\mathbf{x}-\mathbf{x}_k) ]
along the piecewise-linear path [ \mathbf{x}(t) = \mathrm{Proj}_C(\mathbf{x}_k - t\mathbf{r}_k) ].

---

### Break-Point Calculation (Eq. 23 related)

For each component, the break point [ t_i ] is the smallest positive [ t ] at which the projected
path hits a bound.  For phase-field DOFs ( [ \mathrm{lb}_i = d_i^{(n)} ], [ \mathrm{ub}_i = 1 ] ):

[ t_i = \begin{cases}(x_i^0 - \mathrm{ub}_i)/r_i & r_i < 0,\\ (x_i^0 - \mathrm{lb}_i)/r_i & r_i > 0,\\ +\infty & r_i = 0.\end{cases} ]

For displacement DOFs there are no box constraints, so [ t_i = +\infty ] for all [ u ] DOFs.

**Code — break-point computation:**

[`main.cc#L1386-L1450`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L1386)

```cpp
  template <int dim>
  std::priority_queue< std::pair<double, unsigned int>,
                       std::vector<std::pair<double, unsigned int>>,
                       std::greater<std::pair<double, unsigned int>> >
  PhaseFieldMonolithicSolve<dim>::calculate_break_points(const BlockVector<double> & solution_delta,
                                                         const BlockVector<double> & gradient_g,
                                                         BlockVector<double> & gradient_d)
  {
    std::priority_queue< std::pair<double, unsigned int>,
                         std::vector<std::pair<double, unsigned int>>,
                         std::greater<std::pair<double, unsigned int>> >
    break_points_sorted;

    // sentinel value used so that while() loop always has a valid top()
    break_points_sorted.push(std::make_pair(std::numeric_limits<double>::max(),
                                            std::numeric_limits<unsigned int>::max()));

    // displacement DOFs: no box constraints, they never hit a boundary in finite time
    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
      gradient_d.block(m_u_dof)[i] = -gradient_g.block(m_u_dof)[i];

    double t = 0.0;
    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
      {
	if (gradient_g.block(m_d_dof)[i] < 0)
	  t = (solution_phasefield_total[i] - 1.0 ) / gradient_g.block(m_d_dof)[i];
	else if (gradient_g.block(m_d_dof)[i] > 0)
	  t = solution_delta.block(m_d_dof)[i] / gradient_g.block(m_d_dof)[i];
	else
	  t = std::numeric_limits<double>::max();

        if (t > 0)
          {
	    break_points_sorted.push(std::make_pair(t, i));
          }
        else
          {
            gradient_d.block(m_d_dof)[i] = 0;
            if (gradient_g.block(m_d_dof)[i] > 0)
              m_active_set_phasefield(i) = 1; //lower bound
            else
              m_active_set_phasefield(i) = 2; //upper bound
          }
      }
    return break_points_sorted;
  }
```

- Phase-field `gradient_g` = [ \mathbf{r} ] (gradient of [ \Pi ]).
- `gradient_d = -gradient_g` = [ -\mathbf{r} ] (steepest-descent direction [ \mathbf{d} ]).
- If [ t \le 0 ]: the DOF is already on the boundary; it enters the active set immediately.
- Break points are sorted smallest-first using `std::priority_queue` with `std::greater`.

---

### Eq. (23): Cauchy Point Formulas — Initialization

Before the first break point, the model derivative and second derivative are:

[ f'(0) = -\mathbf{d}^T\mathbf{r} = -\mathbf{d}\cdot\mathbf{r} ]
[ f''(0) = \mathbf{d}^T\mathbf{B}_k\mathbf{d} = \mathbf{d}^T\mathbf{B}^0_k\mathbf{d} - (\mathbf{p}^T\mathbf{M}_k)\mathbf{p},\quad \mathbf{p} = \mathbf{W}_k^T\mathbf{d} ]

**Code:**

[`main.cc#L4808-L4834`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4808)

```cpp
    // p = W^T * d
    Vector<double> p(2 * list_size);
    for (unsigned int i = 0; i < list_size; ++i)
      {
        p(i)             = (*std::next(itr_y_begin,    i)) * gradient_d;
        p(i + list_size) = (*std::next(itr_b0xs_begin, i)) * gradient_d;
      }

    Vector<double> c(2 * list_size);
    c = 0.0;

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
```

- `p(i) = y_i * d` and `p(i+m) = b0xs_i * d` → [ \mathbf{p} = \mathbf{W}_k^T\mathbf{d} ]
- `f_prime = -(d * d)` → [ f'(t_0^-) = -\|\mathbf{d}\|^2 = -\mathbf{d}^T\mathbf{r} ]
  (because [ \mathbf{d} = -\mathbf{r} ], so [ -\mathbf{d}\cdot\mathbf{r} = \|\mathbf{d}\|^2 ];
  then sign convention: `f_prime` is stored as [ -\mathbf{d}^T\mathbf{r} ])
- `f_prime_prime = d * B0*d - p * Mp` → [ f''(0) = \mathbf{d}^T\mathbf{B}^0_k\mathbf{d} - \mathbf{p}^T\mathbf{M}_k\mathbf{p} ]
- `delta_t_min = -f_prime / f_prime_prime` → [ \hat{t} = -f'(0)/f''(0) ]

If [ \hat{t} \ge t_{j+1} - t_j ] (first segment), step through break points; else stop.

---

### Eq. (23): Update After Each Break Point (Segment [ j \to j+1 ])

When the minimizer lies beyond break point [ t_j ] (component [ b ] hits its bound), update:

[ f'_j = f'_{j-1} + \Delta t_{j-1} f''_{j-1} + g_b^2 + g_b\,\mathbf{e}_b^T\mathbf{B}^0_k\mathbf{z}_{j-1} - g_b(\mathbf{w}_b^T\mathbf{M}_k\mathbf{c}_{j-1}) ]
[ f''_j = f''_{j-1} + 2g_b\,\mathbf{e}_b^T\mathbf{B}^0_k\mathbf{d}_{j-1} + g_b^2(\mathbf{B}^0_k)_{bb} - 2g_b(\mathbf{w}_b^T\mathbf{M}_k\mathbf{p}_{j-1}) - g_b^2(\mathbf{w}_b^T\mathbf{M}_k\mathbf{w}_b) ]

where [ \mathbf{z}_j = \mathbf{z}_{j-1} + \Delta t_{j-1}\mathbf{d}_{j-1} ] and [ \mathbf{c}_j = \mathbf{c}_{j-1} + \Delta t_{j-1}\mathbf{p}_{j-1} ].

**Code:**

[`main.cc#L4872-L4916`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4872)

```cpp
        // z_{j} = z_{j-1} + delta_t_{j-1} * gradient_d_{j-1};
	z.sadd(1.0, delta_t, gradient_d);

	// c_{j} = c_{j-1} + delta_t_{j-1} * p_{j-1};
	if (list_size > 0)
	  c.sadd(1.0, delta_t, p);

        double temp_scalar = ebT_x_B0_x_v(b, B0_matrix, z);

        // w_b = W^T * e_b
        for (unsigned int i = 0; i < list_size; ++i)
          {
            w_b(i)             = (*std::next(itr_y_begin,    i)).block(m_d_dof)[b];
            w_b(i + list_size) = (*std::next(itr_b0xs_begin, i)).block(m_d_dof)[b];
          }

        if (list_size > 0)
          M_matrix.vmult(w_b_T_x_M, w_b);

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

- `z.sadd(1, delta_t, gradient_d)` → [ \mathbf{z}_j = \mathbf{z}_{j-1} + \Delta t_{j-1}\mathbf{d}_{j-1} ]
- `c.sadd(1, delta_t, p)` → [ \mathbf{c}_j = \mathbf{c}_{j-1} + \Delta t_{j-1}\mathbf{p}_{j-1} ]
- `w_b` → [ \mathbf{w}_b = \mathbf{W}_k^T\mathbf{e}_b ]
- `w_b_T_x_M = M_matrix * w_b` → [ \mathbf{M}_k\mathbf{w}_b ]
- `f_prime += ...` → update of [ f'_j ]
- `f_prime_prime += ...` → update of [ f''_j ]
- `p.sadd(1, g_b, w_b)` → [ \mathbf{p}_j = \mathbf{p}_{j-1} + g_b\mathbf{w}_b ]
- `gradient_d[b] = 0` → zeroes out component [ b ] (it is now fixed at its bound)
- `delta_t_min = -f_prime / f_prime_prime` → [ \hat{t}_j = t_{j-1} - f'_{j}/f''_{j} ]

---

### Cauchy Point — Final Position

After the while-loop finds the segment where [ \hat{t} ] falls, the Cauchy point is set:

[`main.cc#L4927-L4942`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4927)

```cpp
    if (delta_t_min < 0)
      delta_t_min = 0;

    t_old += delta_t_min;

    for(unsigned int i = 0; i < m_active_set_phasefield.size(); ++i)
      {
	// inactive phasefield dof
	if (m_active_set_phasefield(i) < 0.5)
	  solution_delta_cauchy_point.block(m_d_dof)[i] = solution_delta.block(m_d_dof)[i]
						+ t_old * gradient_d.block(m_d_dof)[i];
      }

    // There are no active constraints in the displacement field
    solution_delta_cauchy_point.block(m_u_dof) = solution_delta.block(m_u_dof);
    (solution_delta_cauchy_point.block(m_u_dof)).add(t_old, gradient_d.block(m_u_dof));
```

For inactive DOFs: [ x_i^c = x_i^0 + t^*\, d_i^j ] (the final search distance is `t_old`).
For active DOFs set at their bound during the loop:
[ x_b^c = \mathrm{ub}_b ] or [ x_b^c = \mathrm{lb}_b ] (already set earlier in the loop).

---

## 3.4. Subspace Minimization

### Eq. (24): Subspace Minimization Problem

[ \min_{\hat{\mathbf{x}}} \frac{1}{2}\hat{\mathbf{x}}^T(\mathbf{Z}^T\mathbf{B}_k\mathbf{Z})\hat{\mathbf{x}} - (\mathbf{Z}^T\mathbf{r}^c)^T\hat{\mathbf{x}} ]

where [ \mathbf{Z} ] is the restriction to free (inactive) DOFs,
[ \mathbf{r}^c = \mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k) + \mathbf{r}_k ] is the reduced gradient at the Cauchy point,
and [ \hat{\mathbf{x}} = \mathbf{Z}^T(\mathbf{x} - \mathbf{x}^c) ] is the free-variable update.

The reduced operator [ \mathbf{Z}^T\mathbf{B}_k\mathbf{Z} = \mathbf{Z}^T\mathbf{B}^0_k\mathbf{Z} - \mathbf{Z}^T\mathbf{W}_k\mathbf{M}_k\mathbf{W}_k^T\mathbf{Z} ].

**Code — compute r^c (RHS of subspace problem):**

[`main.cc#L5221-L5272`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5221)

```cpp
	// temp_vector_1 = x^c - x_k
	BlockVector<double> temp_vector_1(solution_delta_cauchy_point);
	temp_vector_1 -= solution_delta;

	// temp_vector_2 = B_0 * (x^c - x_k)
	BlockVector<double> temp_vector_2(m_dofs_per_block);
	m_tangent_matrix.vmult(temp_vector_2, temp_vector_1);

	// temp_vector_3 = W^T * (x^c - x_k)
	Vector<double> temp_vector_3(2 * list_size);
	for (unsigned int i = 0; i < list_size; ++i)
	  {
	    temp_vector_3(i)             = (*std::next(itr_y_begin,    i)) * temp_vector_1;
	    temp_vector_3(i + list_size) = (*std::next(itr_b0xs_begin, i)) * temp_vector_1;
	  }

	// temp_vector_4 = M * W^T * (x^c - x_k)
	Vector<double> temp_vector_4(2 * list_size);
	if (list_size > 0)
	  M_matrix.vmult(temp_vector_4, temp_vector_3);

	// temp_vector_5 = W * M * W^T * (x^c - x_k)
	BlockVector<double> temp_vector_5(m_dofs_per_block);
	for (unsigned int i = 0; i < list_size; ++i)
	  {
	    temp_vector_5.add(temp_vector_4(i),             (*std::next(itr_y_begin,    i)));
	    temp_vector_5.add(temp_vector_4(i + list_size), (*std::next(itr_b0xs_begin, i)));
	  }

	// temp_vector_2 = B_0 * (x^c - x_k) - W * M * W^T * (x^c - x_k) = B_k * (x^c - x_k)
	if (list_size > 0)
	  temp_vector_2 -= temp_vector_5;

	// rhs_vector = -(B_k * (x^c - x_k) + r_k) for inactive DOFs only
	BlockVector<double> rhs_vector(temp_vector_2);
	rhs_vector *= -1;
```

After zeroing out constrained/active DOFs:
- `temp_vector_2 = B_k * (x^c - x_k)` → [ \mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k) ]
- Adding [ \mathbf{r}_k ] (from `m_system_rhs`) gives [ \mathbf{r}^c ]
- `rhs_vector = -r^c` (negative RHS for the linear system)

---

### Primal CG Approach (Section 3.4)

The subspace problem is solved approximately by CG in the free-variable subspace.
The operator applied is [ \mathbf{Z}^T\mathbf{B}_k\mathbf{Z} ] using the compact L-BFGS formula.

**Code — CG path:**

[`main.cc#L5280-L5428`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5280)

```cpp
	if (m_parameters.m_type_linear_solver == "CG")
	  {
	    m_timer.enter_subsection("Subspace CG solve (inverse operator)");

	    const double cg_tol = m_parameters.m_CG_tolerace;

	    zT_B0_z(free_dofs, m_tangent_matrix);

	    const auto op_zT_B0_z = block_operator(m_tangent_matrix);

	    if (list_size > 0)
	      {
		// Build zT_W_matrix (Z^T * W projected to free DOFs)
		// ...
		const auto op_zT_wMwT_z = block_operator<2, 2, BlockVector<double>>({op_uMuT, op_uMdT,
										     op_dMuT, op_dMdT});

		const auto op_total = op_zT_B0_z - op_zT_wMwT_z;

		SolverControl            solver_control(1e5, cg_tol);
		SolverCG<BlockVector<double>> cg(solver_control);

		const auto op_total_inv = inverse_operator(op_total, cg);
		op_total_inv.vmult(search_direction, rhs_vector);
	      }
	    else
	      {
		const auto op_total = op_zT_B0_z;
		// ... CG solve with B0 alone
	      }
	  }
```

- `zT_B0_z(free_dofs, m_tangent_matrix)` → zeroes rows/cols of constrained DOFs in [ \mathbf{B}^0_k ], giving [ \mathbf{Z}^T\mathbf{B}^0_k\mathbf{Z} ]
- `op_zT_wMwT_z` → [ \mathbf{Z}^T\mathbf{W}_k\mathbf{M}_k\mathbf{W}_k^T\mathbf{Z} ] (built from `zT_y_list`, `zT_b0xs_list`)
- `op_total = op_zT_B0_z - op_zT_wMwT_z` → [ \mathbf{Z}^T\mathbf{B}_k\mathbf{Z} ]
- CG solves: [ (\mathbf{Z}^T\mathbf{B}_k\mathbf{Z})\hat{\mathbf{x}} = -\mathbf{Z}^T\mathbf{r}^c ]

---

### Direct Subspace Solve (Section 3.4 — alternative)

For the Direct solver option, the Woodbury identity is applied explicitly:

[ (\mathbf{Z}^T\mathbf{B}_k\mathbf{Z})^{-1} = (\mathbf{Z}^T\mathbf{B}^0_k\mathbf{Z})^{-1} + (\mathbf{Z}^T\mathbf{B}^0_k\mathbf{Z})^{-1}\mathbf{Z}^T\mathbf{W}_k (\mathbf{I} - \mathbf{M}_k\mathbf{W}_k^T\mathbf{Z}(\mathbf{Z}^T\mathbf{B}^0_k\mathbf{Z})^{-1}\mathbf{Z}^T\mathbf{W}_k)^{-1}\mathbf{M}_k\mathbf{W}_k^T\mathbf{Z}(\mathbf{Z}^T\mathbf{B}^0_k\mathbf{Z})^{-1} ]

**Code:**

[`main.cc#L5481-L5620`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5481)

```cpp
	else if (m_parameters.m_type_linear_solver == "Direct")
	  {
	    m_timer.enter_subsection("Subspace direct solve (LU factorization)");
	    zT_B0_z(free_dofs, m_tangent_matrix);
	    SparseDirectUMFPACK zT_B0_z_inv;
	    zT_B0_z_inv.initialize(m_tangent_matrix);

	    m_timer.leave_subsection();
	    m_timer.enter_subsection("Subspace direct solve (LU matrix-vector multiplication)");

	    zT_B0_z_inv.vmult(search_direction, rhs_vector);

	    BlockVector<double> update_vector(m_dofs_per_block);
	    update_vector = 0;
	    if (list_size > 0)
	      {
		// Compute (Z^T B0 Z)^{-1} Z^T W columns (y and b0xs)
		// Build W^T Z (Z^T B0 Z)^{-1} Z^T W matrix (2m x 2m)
		FullMatrix<double> wT_z_zT_B0_z_inv_zT_w(2 * list_size);
		// ...

		// middle_matrix = (I - M * W^T Z (Z^T B0 Z)^{-1} Z^T W)^{-1} * M
		FullMatrix<double> temp_matrix(2 * list_size);
		M_matrix.mmult(temp_matrix, wT_z_zT_B0_z_inv_zT_w);
		FullMatrix<double> middle_matrix(IdentityMatrix(2*list_size));
		middle_matrix.add(-1.0, temp_matrix);
		FullMatrix<double> middle_matrix_inv(2 * list_size);
		middle_matrix_inv.invert(middle_matrix);
		middle_matrix_inv.mmult(middle_matrix, M_matrix);

		// update_vector = Z B0^{-1} Z^T W * middle * W^T Z B0^{-1} Z^T * rhs
		// ...
	      }

	    search_direction += update_vector;
	  }
```

This implements the Woodbury identity in three stages:
1. [ (\mathbf{Z}^T\mathbf{B}^0_k\mathbf{Z})^{-1}\mathbf{r}^c ] via `vmult`.
2. [ (\mathbf{I} - \mathbf{M}_k\mathbf{W}_k^T\mathbf{Z}(\mathbf{Z}^T\mathbf{B}^0_k\mathbf{Z})^{-1}\mathbf{Z}^T\mathbf{W}_k)^{-1}\mathbf{M}_k ] built as `middle_matrix`.
3. `search_direction += update_vector` adds the Woodbury correction.

---

### Combining Cauchy and Subspace Solution

After the subspace solve, the total update is:

[ \mathbf{x}_{k+1} = \mathbf{x}^c + \mathbf{Z}\hat{\mathbf{x}} ]

**Code:**

[`main.cc#L5634-L5648`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5634)

```cpp
	LBFGS_update.block(m_u_dof) = solution_delta_cauchy_point.block(m_u_dof);
	LBFGS_update.block(m_u_dof) += search_direction.block(m_u_dof);
	LBFGS_update.block(m_u_dof) -= solution_delta.block(m_u_dof);

	for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
	  {
	    // phasefield active constraints
	    if (m_active_set_phasefield(i) > 0.5)
	      LBFGS_update.block(m_d_dof)[i] = solution_delta_cauchy_point.block(m_d_dof)[i]
					     - solution_delta.block(m_d_dof)[i];
	    else
	      LBFGS_update.block(m_d_dof)[i] = solution_delta_cauchy_point.block(m_d_dof)[i]
					     + search_direction.block(m_d_dof)[i]
					     - solution_delta.block(m_d_dof)[i];
	  }
```

For active DOFs: `LBFGS_update[i] = x^c[i] - x_k[i]` (move to bound).
For inactive DOFs: `LBFGS_update[i] = x^c[i] + search_direction[i] - x_k[i]` (add subspace correction).

---

## 3.5. L-BFGS Two-Loop Recursion (unconstrained version, `solve_nonlinear_timestep_LBFGS`)

### L-BFGS Search Direction (Two-Loop)

[ \mathbf{p}_k = -\mathbf{H}_k\mathbf{r}_k ]

using the Nocedal–Wright two-loop recursion where [ \mathbf{H}_k = \mathbf{B}_k^{-1} ].

**Code — loop 1 (backward):**

[`main.cc#L4639-L4671`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4639)

```cpp
        // LBFGS algorithm
        LBFGS_q_vector = m_system_rhs;

        LBFGS_alpha_list.clear();
        for (auto itr = LBFGS_vector_list.begin(); itr != LBFGS_vector_list.end(); ++itr)
          {
            LBFGS_s_vector = (itr->first).first;
            LBFGS_y_vector = (itr->first).second;
            rho = itr->second;

            const double alpha = rho * (LBFGS_s_vector * LBFGS_q_vector);
            LBFGS_alpha_list.push_back(alpha);

            LBFGS_q_vector.add(-alpha, LBFGS_y_vector);
          }

        LBFGS_B0(LBFGS_r_vector,
		 LBFGS_q_vector);

        for (auto itr = LBFGS_vector_list.rbegin(); itr != LBFGS_vector_list.rend(); ++itr)
          {
            LBFGS_s_vector = (itr->first).first;
            LBFGS_y_vector = (itr->first).second;
            rho = itr->second;

            LBFGS_beta = rho * (LBFGS_y_vector * LBFGS_r_vector);

            const double alpha = LBFGS_alpha_list.back();
            LBFGS_alpha_list.pop_back();

            LBFGS_r_vector.add(alpha - LBFGS_beta, LBFGS_s_vector);
          }

        LBFGS_r_vector *= -1.0; // this is the p_vector (search direction)
```

Standard two-loop recursion with [ \rho_i = 1/(\mathbf{y}_i^T\mathbf{s}_i) ],
[ \alpha_i = \rho_i\,\mathbf{s}_i^T\mathbf{q} ], [ \beta_i = \rho_i\,\mathbf{y}_i^T\mathbf{r} ].
Initial Hessian [ \mathbf{H}^0_k ] is applied via `LBFGS_B0(r, q)` which computes
[ \mathbf{r} = (\mathbf{B}^0_k)^{-1}\mathbf{q} ] (block-diagonal linear solve).
Final `*= -1.0` gives the search direction [ \mathbf{p}_k = -\mathbf{H}_k\mathbf{r}_k ].

---

## 3.5. Line Search with Wolfe Conditions

### Eq. (25): Strong Wolfe Conditions

[ \phi(\alpha) \leq \phi(0) + c_1\alpha\phi'(0) \quad \text{(sufficient decrease, Armijo)} ]
[ |\phi'(\alpha)| \leq c_2\,|\phi'(0)| \quad \text{(curvature condition)} ]

with [ 0 < c_1 < c_2 < 1 ].

**Code:**

[`main.cc#L4196-L4256`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4196)

```cpp
  double PhaseFieldMonolithicSolve<dim>::line_search_stepsize_strong_wolfe(const double phi_0,
				                                           const double phi_0_prime,
				                                           const BlockVector<double> & BFGS_p_vector,
				                                           const BlockVector<double> & solution_delta)
  {
    const double c1 = 0.0001;
    const double c2 = 0.9;
    const double alpha_max = 1.0;
    const unsigned int max_iter = 20;
    double alpha = 1.0;

    double phi_old = phi_0;
    double phi_prime_old = phi_0_prime;
    double alpha_old = 0.0;

    double phi, phi_prime;

    unsigned int i = 0;
    for (; i < max_iter; ++i)
      {
	current_phi_phi_prime = calculate_phi_and_phi_prime(alpha, BFGS_p_vector, solution_delta);
	phi = current_phi_phi_prime.first;
	phi_prime = current_phi_phi_prime.second;

	if (   ( phi > (phi_0 + c1 * alpha * phi_0_prime) )
	    || ( i > 0 && phi > phi_old ) )
	  {
	    return line_search_zoom_strong_wolfe(...);
	  }

	if (std::fabs(phi_prime) <= c2 * std::fabs(phi_0_prime))
	  {
	    return alpha;
	  }

	if (phi_prime >= 0)
	  {
	    return line_search_zoom_strong_wolfe(...);
	  }

	phi_old = phi;
	phi_prime_old = phi_prime;
	alpha_old = alpha;

	alpha = std::min(0.6*alpha, alpha_max);
      }

    return alpha;
  }
```

- `c1 = 0.0001` → [ c_1 ] (Armijo constant)
- `c2 = 0.9` → [ c_2 ] (curvature constant)
- `phi > phi_0 + c1 * alpha * phi_0_prime` → Armijo condition violated → zoom
- `fabs(phi_prime) <= c2 * fabs(phi_0_prime)` → strong curvature condition → accept
- `phi_prime >= 0` → overshot → zoom

---

### [ \phi(\alpha) ] and [ \phi'(\alpha) ] Computation

[ \phi(\alpha) = \Pi(\mathbf{x}_k + \alpha\mathbf{p}_k), \quad \phi'(\alpha) = \nabla\Pi(\mathbf{x}_k + \alpha\mathbf{p}_k)\cdot\mathbf{p}_k ]

**Code:**

[`main.cc#L4347-L4368`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4347)

```cpp
  std::pair<double, double> PhaseFieldMonolithicSolve<dim>::
    calculate_phi_and_phi_prime(const double alpha,
				const BlockVector<double> & BFGS_p_vector,
				const BlockVector<double> & solution_delta)
  {
    std::pair<double, double> phi_values;

    BlockVector<double> solution_delta_trial(solution_delta);
    solution_delta_trial.add(alpha, BFGS_p_vector);

    update_qph_incremental(solution_delta_trial, m_solution);

    BlockVector<double> system_rhs(m_dofs_per_block);
    assemble_system_rhs_BFGS_parallel(m_solution, system_rhs);

    phi_values.first = calculate_energy_functional();
    phi_values.second = system_rhs * BFGS_p_vector;
    return phi_values;
  }
```

- `solution_delta_trial = x_k + alpha * p` → [ \mathbf{x}_k + \alpha\mathbf{p}_k ]
- `calculate_energy_functional()` → [ \phi(\alpha) = \Pi(\mathbf{x}_k + \alpha\mathbf{p}_k) ]
- `system_rhs * BFGS_p_vector` → [ \phi'(\alpha) = \mathbf{r}^T\mathbf{p}_k ] (note: RHS = gradient)

---

### Zoom Algorithm (Bisection)

[`main.cc#L4258-L4311`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4258)

```cpp
  double PhaseFieldMonolithicSolve<dim>::
    line_search_zoom_strong_wolfe(double phi_low, double phi_low_prime, double alpha_low,
				  double phi_high, double phi_high_prime, double alpha_high,
				  double phi_0, double phi_0_prime, const BlockVector<double> & BFGS_p_vector,
				  double c1, double c2, unsigned int max_iter, const BlockVector<double> & solution_delta)
  {
    double alpha = 0;
    unsigned int i = 0;
    for (; i < max_iter; ++i)
      {
	// a simple bisection is faster than cubic interpolation
	alpha = 0.5 * (alpha_low + alpha_high);
	current_phi_phi_prime = calculate_phi_and_phi_prime(alpha, BFGS_p_vector, solution_delta);
	phi = current_phi_phi_prime.first;
	phi_prime = current_phi_phi_prime.second;

	if (   (phi > phi_0 + c1 * alpha * phi_0_prime)
	    || (phi > phi_low) )
	  {
	    alpha_high = alpha;
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
	      }

	    alpha_low = alpha;
	  }
      }
    return alpha;
  }
```

Bisection maintains an interval [ [\alpha_\mathrm{low}, \alpha_\mathrm{high}] ] containing a
Wolfe-condition satisfying point, bisecting at each step.
A commented-out cubic interpolation variant (`line_search_interpolation_cubic`) is also available.

---

### Cubic Interpolation (Optional — commented out)

[ \alpha^* = \alpha_1 - (\alpha_1 - \alpha_0)\frac{\phi'_1 + d_2 - d_1}{\phi'_1 - \phi'_0 + 2d_2}, \quad d_1 = \phi'_0 + \phi'_1 - 3\frac{\phi_0 - \phi_1}{\alpha_0 - \alpha_1},\quad d_2 = \sqrt{d_1^2 - \phi'_0\phi'_1} ]

**Code:**

[`main.cc#L4313-L4345`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4313)

```cpp
  double PhaseFieldMonolithicSolve<dim>::
    line_search_interpolation_cubic(const double alpha_0, const double phi_0, const double phi_0_prime,
  			            const double alpha_1, const double phi_1, const double phi_1_prime)
  {
    const double d1 = phi_0_prime + phi_1_prime - 3.0 * (phi_0 - phi_1) / (alpha_0 - alpha_1);

    const double temp = d1 * d1 - phi_0_prime * phi_1_prime;

    if (temp < 0.0)
      return 0.5 * (alpha_0 + alpha_1);

    int sign;
    if (alpha_1 > alpha_0)
      sign = 1;
    else
      sign = -1;

    const double d2 = sign * std::sqrt(temp);

    const double alpha = alpha_1 - (alpha_1 - alpha_0)
	               * (phi_1_prime + d2 - d1) / (phi_1_prime - phi_0_prime + 2*d2);

    if (    (alpha_1 > alpha_0)
	 && (alpha > alpha_1 || alpha < alpha_0))
      return 0.5 * (alpha_0 + alpha_1);

    if (    (alpha_0 > alpha_1)
	 && (alpha > alpha_0 || alpha < alpha_1))
      return 0.5 * (alpha_0 + alpha_1);

    return alpha;
  }
```

`d1` → [ d_1 ], `d2` → [ d_2 ], `alpha` → [ \alpha^* ].
Falls back to bisection if the interpolant is outside [ [\alpha_0, \alpha_1] ].

---

## 4. Section 4: Numerical Examples

### Eq. (33): Crack Dissipation Energy

[ \Gamma = g_c\int_\Omega\gamma(d,\nabla d)\,\mathrm{d}\Omega = g_c\int_\Omega\left(\frac{d^2}{2l} + \frac{l}{2}|\nabla d|^2\right)\mathrm{d}\Omega ]

**Code — per-quadrature-point contribution (stored in m_crack_energy_dissipation, already scaled by g_c):**

[`main.cc#L886-L891`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L886)

```cpp
    m_crack_energy_dissipation = m_gc * (  0.5 / m_length_scale * m_phase_field_value * m_phase_field_value
	                                   + 0.5 * m_length_scale * m_grad_phasefield * m_grad_phasefield)
	                                   // the term due to viscosity regularization
	                                   + (m_phase_field_value - phase_field_value_previous_step)
				           * (m_phase_field_value - phase_field_value_previous_step)
				           * 0.5 * m_eta / delta_time;
```

- [ g_c\frac{d^2}{2l} ] → `m_gc * 0.5 / m_length_scale * d * d`
- [ g_c\frac{l}{2}|\nabla d|^2 ] → `m_gc * 0.5 * m_length_scale * grad_d * grad_d`
- Viscosity term [ \frac{\eta}{2\Delta t}(d - d^{(n)})^2 ] is also included here (from Eq. 13).

**Code — global integral over [ \Omega ]:**

[`main.cc#L3064-L3087`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3064)

```cpp
    double energy = 0.0;
    for (auto & cell : m_dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          // ...
          for (unsigned int q_point = 0; q_point < m_n_q_points; ++q_point)
            {
              const double JxW = lqph[q_point]->get_JxW();
              cell_energy += (  lqph[q_point]->get_strain_energy_total()
                              + lqph[q_point]->get_crack_energy_dissipation() ) * JxW;
            }
          energy += cell_energy;
        }
    energy = Utilities::MPI::sum(energy, mpi_communicator);
    return energy;
```

`get_crack_energy_dissipation() * JxW` integrated over all cells gives
[ g_c\int_\Omega\gamma\,\mathrm{d}\Omega = g_c\Gamma_l(d) ] (the second term of Eq. 1).

---

### Torsion Boundary Conditions (3D Torsion Scenario)

For the 3D torsion test (scenario 11), nodes on the left face at [ x = 0 ] are given rotational
displacements about the [ x ]-axis:

[ u_y = z\tan\theta, \quad u_z = -y\tan\theta, \quad \tan\theta = \Omega\,\Delta t ]

where [ \Omega ] is the prescribed rotation rate and [ \Delta t ] is the time increment.

**Code:**

[`main.cc#L3565-L3620`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3565)

```cpp
	else if (m_parameters.m_scenario == 11)
	  {
	    // Dirichlet B.C. right surface: all displacements = 0
	    const int boundary_id_right_surface = 0;
	    VectorTools::interpolate_boundary_values(m_dof_handler,
						     boundary_id_right_surface,
						     Functions::ZeroFunction<dim>(m_n_components),
						     m_constraints,
						     m_fe.component_mask(displacements));

	    // Dirichlet B.C. left surface: x-displacement = 0
	    const int boundary_id_left_surface = 1;
	    VectorTools::interpolate_boundary_values(m_dof_handler,
						     boundary_id_left_surface,
						     Functions::ZeroFunction<dim>(m_n_components),
						     m_constraints,
						     m_fe.component_mask(x_displacement));

	    typename Triangulation<dim>::active_vertex_iterator vertex_itr;
	    vertex_itr = m_triangulation.begin_active_vertex();
	    std::vector<types::global_dof_index> node_rotate(m_fe.dofs_per_vertex);
	    double node_dist = 0.0;
	    double disp_mag = 0.0;
	    double angle_theta = 0.0;
	    double disp_y = 0;
	    double disp_z = 0;

	    for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)
	      {
		if (std::fabs(vertex_itr->vertex()[0] - 0.0) < 1.0e-9)
		  {
		    node_rotate = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);
		    node_dist = std::sqrt(  vertex_itr->vertex()[1] * vertex_itr->vertex()[1]
			                  + vertex_itr->vertex()[2] * vertex_itr->vertex()[2]);

		    angle_theta = m_time.get_delta_t() * m_time.get_magnitude();
		    disp_mag = node_dist * std::tan(angle_theta);

		    if (node_dist > 0)
		      {
		        disp_y = vertex_itr->vertex()[2]/node_dist * disp_mag;
		        disp_z = -vertex_itr->vertex()[1]/node_dist * disp_mag;
		      }
		    else
		      {
			disp_y = 0.0;
			disp_z = 0.0;
		      }

		    m_constraints.add_line(node_rotate[1]);
		    m_constraints.set_inhomogeneity(node_rotate[1], disp_y);

		    m_constraints.add_line(node_rotate[2]);
		    m_constraints.set_inhomogeneity(node_rotate[2], disp_z);
		  }
	      }
	  }
```

Mapping to formulas:
- `angle_theta = delta_t * magnitude` → [ \theta = \Omega\,\Delta t ]
- `disp_mag = node_dist * tan(angle_theta)` → [ |\mathbf{u}| = r\tan\theta ]
- `disp_y = z / node_dist * disp_mag` → [ u_y = \frac{z}{r}\cdot r\tan\theta = z\tan\theta ]
- `disp_z = -y / node_dist * disp_mag` → [ u_z = -\frac{y}{r}\cdot r\tan\theta = -y\tan\theta ]

Geometric interpretation: the displacement vector is tangent to the circle of radius [ r = \sqrt{y^2+z^2} ],
perpendicular to the radial direction, producing pure torsion about the [ x ]-axis.

---

## Summary Table

| Formula | Description | Code Location |
|---------|-------------|---------------|
| Eq. (1) | Total energy [ \Pi ] | [`main.cc#L3064`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3064) |
| Eq. (2) | Crack surface density [ \Gamma_l ] | [`main.cc#L886`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L886) |
| Eq. (3) | Minimization problem | [`main.cc#L4997`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4997) |
| Eq. (4) | Irreversibility constraints | [`main.cc#L5730`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5730) |
| Box constraints | [ \mathrm{lb} \le x \le \mathrm{ub} ] | [`main.cc#L1398`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L1398) |
| Eq. (5) | [ \psi = [g(d)+k]\psi^+ + \psi^- ] | [`main.cc#L884`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L884) |
| Eq. (6) | [ g(d) = (1-d)^2 ] | [`main.cc#L244`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L244) |
| [ \langle\cdot\rangle_\pm ], [ H ] | Ramp / Heaviside functions | [`SpectrumDecomposition.cc#L19`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.cc#L19) |
| Spectral decomp | [ \boldsymbol{\epsilon} = \sum_\alpha\epsilon_\alpha\mathbf{M}_\alpha ] | [`SpectrumDecomposition.h#L29`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L29) |
| [ \boldsymbol{\epsilon}^\pm ] | Positive/negative strain | [`SpectrumDecomposition.h#L45`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L45) |
| Eq. (7) | [ \psi^+ ] | [`main.cc#L876`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L876) |
| Eq. (8) | [ \psi^- ] | [`main.cc#L880`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L880) |
| Eq. (9) | [ \boldsymbol{\sigma}^\pm ] | [`main.cc#L857`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L857) |
| Eq. (10) | [ \boldsymbol{\sigma} = [g(d)+k]\boldsymbol{\sigma}^+ + \boldsymbol{\sigma}^- ] | [`main.cc#L864`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L864) |
| Eq. (11) | Elasticity tensor C | [`main.cc#L867`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L867) + [`SpectrumDecomposition.h#L71`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L71) |
| Eq. (12) | Viscosity term [ \eta\dot{d} ] | [`main.cc#L3813`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3813) |
| Eq. (13) | Viscosity energy | [`main.cc#L886`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L886) |
| Eq. (14) | FE discretization of u | [`main.cc#L3793`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3793) |
| Eq. (15) | FE discretization of d | [`main.cc#L3787`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3787) |
| Eq. (16) | Residual u | [`main.cc#L3804`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3804) |
| Eq. (17) | Residual d | [`main.cc#L3811`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3811) |
| Eq. (18) | Hessian K_uu, K_dd | [`main.cc#L3940`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3940) |
| Eq. (19) | Projected gradient criterion | [`main.cc#L1460`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L1460) |
| Eq. (20) | Compact L-BFGS [ \mathbf{B}_k = \mathbf{B}^0_k - \mathbf{WMW}^T ] | [`main.cc#L5221`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5221) |
| Eq. (21) | W matrix [ [\mathbf{Y}\;\mathbf{B}^0\mathbf{S}] ] | [`main.cc#L5119`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5119) |
| Eq. (22) | M^{-1} | [`main.cc#L5140`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5140) |
| Eq. (23) | Cauchy point formulas | [`main.cc#L4808`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4808) + [`main.cc#L4872`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4872) |
| Break points [ t_j^i ] | Bound intersection times | [`main.cc#L1386`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L1386) |
| Eq. (24) | Subspace minimization | [`main.cc#L5221`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5221) + [`main.cc#L5280`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5280) |
| Primal CG | CG in free subspace | [`main.cc#L5280`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5280) |
| Eq. (25) | Wolfe conditions | [`main.cc#L4196`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4196) |
| Eq. (33) | Crack dissipation energy | [`main.cc#L886`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L886) |
| Torsion BCs | [ u_y = z\tan\theta,\; u_z=-y\tan\theta ] | [`main.cc#L3565`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3565) |
