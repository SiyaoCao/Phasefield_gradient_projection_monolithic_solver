# Formula-to-Code Mapping (Part 1): Phase-Field Fracture Monolithic Solver

**Repository:** https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver

This document maps every formula in `explain.md` to the corresponding implementation in
`main.cc`, `SpectrumDecomposition.h`, `SpectrumDecomposition.cc`, and `Utilities.h`.
Continued in `formula_code_mapping_2.md` for Sections 3.1 (continued), 3.3, 3.4, and 4.

---

## 1. Introduction

### Eq. (1): Total Energy Functional

[ \Pi(\mathbf{u}, d) = \int_\Omega \psi(\boldsymbol{\epsilon}(\mathbf{u}), d)\,\mathrm{d}\Omega + g_c\Gamma_l(d) - \int_\Omega \mathbf{b}\cdot\mathbf{u}\,\mathrm{d}\Omega - \int_{\partial\Omega} \mathbf{t}\cdot\mathbf{u}\,\mathrm{d}\Gamma ]

The total energy is assembled and evaluated in `calculate_energy_functional()`.

[`main.cc#L3064-L3087`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3064)

```cpp
  template <int dim>
  double PhaseFieldMonolithicSolve<dim>::calculate_energy_functional() const
  {
    double energy = 0.0;
    for (auto & cell : m_dof_handler.active_cell_iterators())
      if (cell->is_locally_owned())
        {
          const std::vector<std::shared_ptr<PointHistory<dim>>> lqph =
            m_quadrature_point_history.get_data(cell);
          Assert(lqph.size() == m_n_q_points, ExcInternalError());

          double cell_energy = 0.0;
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
  }
```

- `get_strain_energy_total()` returns [ [g(d)+k]\psi^+(\boldsymbol{\epsilon}) + \psi^-(\boldsymbol{\epsilon}) ] (the elastic portion).
- `get_crack_energy_dissipation()` returns [ g_c\gamma(d,\nabla d) ] at each quadrature point (fracture energy density).
- The body force term [ -\int_\Omega\mathbf{b}\cdot\mathbf{u}\,\mathrm{d}\Omega ] and traction term are assembled into the RHS during gradient computation (see Eq. (9) mapping below).

---

### Eq. (2): Crack Surface Density

[ \Gamma_l(d) = \int_\Omega \gamma(d,\nabla d)\,\mathrm{d}\Omega = \int_\Omega \frac{1}{2l}\left(d^2 + l^2\nabla d\cdot\nabla d\right)\mathrm{d}\Omega ]

**Code — crack energy dissipation stored per quadrature point:**

[`main.cc#L886-L891`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L886)

```cpp
    m_crack_energy_dissipation = m_gc * (  0.5 / m_length_scale * m_phase_field_value * m_phase_field_value
	                                   + 0.5 * m_length_scale * m_grad_phasefield * m_grad_phasefield)
	                                   // the term due to viscosity regularization
	                                   + (m_phase_field_value - phase_field_value_previous_step)
				           * (m_phase_field_value - phase_field_value_previous_step)
				           * 0.5 * m_eta / delta_time;
```

The first two terms implement [ g_c\left(\frac{d^2}{2l} + \frac{l}{2}|\nabla d|^2\right) ], which
is [ g_c\gamma(d,\nabla d) ]. The last two lines add the optional viscosity regularization (see Eqs. 12–13).

---

### Eq. (3): Constrained Minimization Problem

[ (\mathbf{u}_{n+1}, d_{n+1}) = \arg\min\;\Pi(\mathbf{u}, d) ]

subject to the inequality constraints in Eq. (4).

**Code — outer L-BFGS-B solver loop:**

[`main.cc#L4997-L5000`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4997)

```cpp
    for (; LBFGS_iteration < m_parameters.m_max_iterations_BFGS; ++LBFGS_iteration)
      {
	if (m_parameters.m_output_iteration_history)
	  m_logfile << '\t' << std::setw(4) << LBFGS_iteration << ' '
                    << std::flush;
```

This loop in `solve_nonlinear_timestep_LBFGS_B` repeatedly minimizes the quadratic model
of [ \Pi ] subject to box constraints until convergence.

---

### Eq. (4): Inequality (Irreversibility) Constraints

[ 0 \leq d_n \leq d_{n+1} \leq 1 ]

**Code — active-set detection per node:**

[`main.cc#L5730-L5760`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5730)

```cpp
	for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
	  {
	    if (   solution_delta.block(m_d_dof)[i] == 0.0
		&& solution_phasefield_total[i] == 1.0)
	      {
		m_active_set_phasefield(i) = 3; //lower bound overlaps with upper bound
		++number_active_constraint_lowerupper_bound;
	      }
	    else if (   solution_delta.block(m_d_dof)[i] == 0.0
		     && solution_phasefield_total[i] != 1.0)
	      {
		m_active_set_phasefield(i) = 1; //lower bound
		++number_active_constraint_lower_bound;
	      }
	    else if (   solution_phasefield_total[i] == 1.0
		     && solution_delta.block(m_d_dof)[i] != 0.0)
	      {
	        m_active_set_phasefield(i) = 2; //upper bound
	        ++number_active_constraint_upper_bound;
	      }
```

At the lower bound: [ d_i = d_i^{(n)} ] (status `1`; encoded as `solution_delta[i] == 0`).
At the upper bound: [ d_i = 1 ] (status `2`).
When both bounds coincide: [ d_i^{(n)} = d_i = 1 ] (status `3`).

---

### Box Constraints (General Form)

[ \mathrm{lb}_i \leq x_i \leq \mathrm{ub}_i ]

**Code — break-point computation uses ub and lb directly:**

[`main.cc#L1398-L1443`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L1398)

```cpp
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
        else // if t == 0, i is in the active set
          {
            gradient_d.block(m_d_dof)[i] = 0;
            if (gradient_g.block(m_d_dof)[i] > 0)
              m_active_set_phasefield(i) = 1; //lower bound
            else
              m_active_set_phasefield(i) = 2; //upper bound
          }
      }
```

- `ub_i = 1.0` (hard-coded upper bound for phase-field).
- `lb_i = m_solution.block(m_d_dof)[i]` (phase-field from previous time step).
- When [ r_i < 0 ]: [ t_i = (x_i^0 - 1.0)/r_i ] (distance to upper bound).
- When [ r_i > 0 ]: [ t_i = \text{solution\_delta}[i]/r_i ] (distance to lower bound in delta-space).

---

### Phase-Field Box Constraints

[ \mathrm{lb}_i = d_i^{(n)} \leq d_i^{(n+1)} \leq 1 = \mathrm{ub}_i ]

**Code — feasibility clamp after subspace minimization:**

[`main.cc#L5650-L5661`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5650)

```cpp
	// make sure the phasefield solutions are feasible
	for(unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
	  {
	    if (solution_delta.block(m_d_dof)[i] + LBFGS_update.block(m_d_dof)[i] < 0.0)
	      LBFGS_update.block(m_d_dof)[i] = -solution_delta.block(m_d_dof)[i];

	    if (  solution_delta.block(m_d_dof)[i]
		+ m_solution.block(m_d_dof)[i]
		+ LBFGS_update.block(m_d_dof)[i] > 1.0)
	      LBFGS_update.block(m_d_dof)[i] = 1.0 - m_solution.block(m_d_dof)[i]
					   - solution_delta.block(m_d_dof)[i];
	  }
```

After the CG/direct subspace solve, the result is clamped so that
[ \mathrm{lb}_i \leq x_i^{\text{new}} \leq \mathrm{ub}_i ] holds for every phase-field DOF.

---

## 2. Phase-Field Formulation and Finite Element Discretization

### 2.1. Phase-Field Formulation

#### Eq. (5): Strain Energy Density with Degradation

[ \psi(\boldsymbol{\epsilon}, d) = [g(d) + k]\psi^+(\boldsymbol{\epsilon}) + \psi^-(\boldsymbol{\epsilon}) ]

**Code:**

[`main.cc#L846-L884`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L846)

```cpp
    const double degradation = degradation_function(m_phase_field_value) + m_residual_k;
    // ...
    m_strain_energy_total = degradation * m_strain_energy_positive + m_strain_energy_negative;
```

`degradation = g(d) + k` where `k = m_residual_k` is a small numerical constant.
`m_strain_energy_positive` = [ \psi^+ ], `m_strain_energy_negative` = [ \psi^- ].

---

#### Eq. (6): Degradation Function

[ g(d) = (1-d)^2 ]

**Code:**

[`main.cc#L244-L258`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L244)

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

- `degradation_function(d)` → [ g(d) = (1-d)^2 ]
- `degradation_function_derivative(d)` → [ g'(d) = 2(d-1) = -2(1-d) ]
- `degradation_function_2nd_order_derivative(d)` → [ g''(d) = 2 ]

---

#### Operators [ \langle x\rangle_+ ], [ \langle x\rangle_- ], and [ H(x) ]

[ \langle x\rangle_+ = \max(x,\,0), \quad \langle x\rangle_- = \min(x,\,0), \quad H(x) = \begin{cases}1 & x > 0,\\ 1/2 & x = 0,\\ 0 & x < 0.\end{cases} ]

**Code:**

[`SpectrumDecomposition.cc#L19-L38`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.cc#L19)

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

- `positive_ramp_function(x)` = `std::fmax(x, 0)` → [ \langle x\rangle_+ ]
- `negative_ramp_function(x)` = `std::fmin(x, 0)` → [ \langle x\rangle_- ]
- `heaviside_function(x)` → [ H(x) ] with [ H(0) = 0.5 ] for symmetry at repeated eigenvalues

---

#### Spectral Decomposition [ \boldsymbol{\epsilon} = \sum_\alpha \epsilon_\alpha \mathbf{M}_\alpha ]

[ \boldsymbol{\epsilon} = \sum_{\alpha=1}^{\dim} \epsilon_\alpha \,\mathbf{M}_\alpha, \quad \mathbf{M}_\alpha = \mathbf{n}_\alpha \otimes \mathbf{n}_\alpha ]

**Code:**

[`SpectrumDecomposition.h#L29-L43`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L29)

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
```

`myEigenvalues[alpha]` = [ \epsilon_\alpha ] (eigenvalue), `myEigenvectors[alpha]` = [ \mathbf{n}_\alpha ] (eigenvector). Their outer product [ \mathbf{M}_\alpha = \mathbf{n}_\alpha\otimes\mathbf{n}_\alpha ] is used in `positive_tensor` and `negative_tensor`.

**Call site in material update:**

[`main.cc#L831-L837`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L831)

```cpp
    usr_spectrum_decomposition::spectrum_decomposition<dim>(m_strain,
  							      eigenvalues,
  							      eigenvectors);

    SymmetricTensor<2, dim> strain_positive, strain_negative;
    strain_positive = usr_spectrum_decomposition::positive_tensor(eigenvalues, eigenvectors);
    strain_negative = usr_spectrum_decomposition::negative_tensor(eigenvalues, eigenvectors);
```

---

#### Positive and Negative Strain Tensors

[ \boldsymbol{\epsilon}^+ = \sum_\alpha \langle\epsilon_\alpha\rangle_+ \mathbf{M}_\alpha, \quad \boldsymbol{\epsilon}^- = \sum_\alpha \langle\epsilon_\alpha\rangle_- \mathbf{M}_\alpha ]

**Code:**

[`SpectrumDecomposition.h#L45-L69`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L45)

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
```

- `positive_ramp_function(eigenvalues[i])` gives [ \langle\epsilon_\alpha\rangle_+ ]
- `outer_product(eigenvectors[i], eigenvectors[i])` gives [ \mathbf{M}_\alpha ]
- The loop sums to [ \boldsymbol{\epsilon}^+ ] (or [ \boldsymbol{\epsilon}^- ]) exactly as in the formula.

---

#### Eq. (7): Positive Strain Energy Density

[ \psi^+(\boldsymbol{\epsilon}) = \tfrac{1}{2}\lambda\langle\mathrm{tr}\,\boldsymbol{\epsilon}\rangle_+^2 + \mu\,\boldsymbol{\epsilon}^+:\boldsymbol{\epsilon}^+ ]

**Code:**

[`main.cc#L876-L878`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L876)

```cpp
    m_strain_energy_positive = 0.5 * my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)
                                                   * usr_spectrum_decomposition::positive_ramp_function(I_1)
                             + m_lame_mu * strain_positive * strain_positive;
```

- `0.5 * lambda * pos_ramp(I_1)^2` → [ \tfrac{1}{2}\lambda\langle\mathrm{tr}\,\boldsymbol{\epsilon}\rangle_+^2 ]
- `mu * strain_positive * strain_positive` → [ \mu\,\boldsymbol{\epsilon}^+:\boldsymbol{\epsilon}^+ ]

where `I_1 = trace(m_strain)` = [ \mathrm{tr}\,\boldsymbol{\epsilon} ].

---

#### Eq. (8): Negative Strain Energy Density

[ \psi^-(\boldsymbol{\epsilon}) = \tfrac{1}{2}\lambda\langle\mathrm{tr}\,\boldsymbol{\epsilon}\rangle_-^2 + \mu\,\boldsymbol{\epsilon}^-:\boldsymbol{\epsilon}^- ]

**Code:**

[`main.cc#L880-L882`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L880)

```cpp
    m_strain_energy_negative = 0.5 * my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)
                                                   * usr_spectrum_decomposition::negative_ramp_function(I_1)
                             + m_lame_mu * strain_negative * strain_negative;
```

- `0.5 * lambda * neg_ramp(I_1)^2` → [ \tfrac{1}{2}\lambda\langle\mathrm{tr}\,\boldsymbol{\epsilon}\rangle_-^2 ]
- `mu * strain_negative * strain_negative` → [ \mu\,\boldsymbol{\epsilon}^-:\boldsymbol{\epsilon}^- ]

---

#### Eq. (9): Cauchy Stress Components [ \boldsymbol{\sigma}^+ ] and [ \boldsymbol{\sigma}^- ]

[ \boldsymbol{\sigma}^+ = \lambda\langle\mathrm{tr}\,\boldsymbol{\epsilon}\rangle_+\mathbf{I} + 2\mu\,\boldsymbol{\epsilon}^+, \quad \boldsymbol{\sigma}^- = \lambda\langle\mathrm{tr}\,\boldsymbol{\epsilon}\rangle_-\mathbf{I} + 2\mu\,\boldsymbol{\epsilon}^- ]

**Code:**

[`main.cc#L857-L862`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L857)

```cpp
    stress_positive = my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)
                                    * Physics::Elasticity::StandardTensors<dim>::I
                    + 2 * m_lame_mu * strain_positive;
    stress_negative = my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)
                                    * Physics::Elasticity::StandardTensors<dim>::I
    		      + 2 * m_lame_mu * strain_negative;
```

- `lambda * pos_ramp(I_1) * I` → [ \lambda\langle\mathrm{tr}\,\boldsymbol{\epsilon}\rangle_+\mathbf{I} ]
- `2 * mu * strain_positive` → [ 2\mu\,\boldsymbol{\epsilon}^+ ]

---

#### Eq. (10): Cauchy Stress with Degradation

[ \boldsymbol{\sigma} = \frac{\partial\psi}{\partial\boldsymbol{\epsilon}} = [g(d)+k]\boldsymbol{\sigma}^+ + \boldsymbol{\sigma}^- ]

**Code:**

[`main.cc#L864`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L864)

```cpp
    m_stress = degradation * stress_positive + stress_negative;
```

`degradation = g(d) + k` multiplies the positive stress part only, implementing the
tension–compression asymmetry exactly as in Eq. (10).

---

#### Eq. (11): Material Tangent Modulus (Elasticity Tensor C)

[ \mathbf{C} = \frac{\partial\boldsymbol{\sigma}}{\partial\boldsymbol{\epsilon}} = [g(d)+k]\left[\lambda H(\mathrm{tr}\,\boldsymbol{\epsilon})\,\mathbf{I}\otimes\mathbf{I} + 2\mu\,\mathbb{P}^+\right] + \left[\lambda H(-\mathrm{tr}\,\boldsymbol{\epsilon})\,\mathbf{I}\otimes\mathbf{I} + 2\mu\,\mathbb{P}^-\right] ]

with [ \mathbb{P}^+ = \frac{\partial\boldsymbol{\epsilon}^+}{\partial\boldsymbol{\epsilon}} ] and [ \mathbb{P}^- = \frac{\partial\boldsymbol{\epsilon}^-}{\partial\boldsymbol{\epsilon}} ].

**Code — C assembly:**

[`main.cc#L867-L874`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L867)

```cpp
    SymmetricTensor<4, dim> C_positive, C_negative;
    C_positive = my_lambda * usr_spectrum_decomposition::heaviside_function(I_1)
                               * Physics::Elasticity::StandardTensors<dim>::IxI
		 + 2 * m_lame_mu * projector_positive;
    C_negative = my_lambda * usr_spectrum_decomposition::heaviside_function(-I_1)
                               * Physics::Elasticity::StandardTensors<dim>::IxI
    		 + 2 * m_lame_mu * projector_negative;
    m_mechanical_C = degradation * C_positive + C_negative;
```

- `heaviside_function(I_1)` → [ H(\mathrm{tr}\,\boldsymbol{\epsilon}) ]
- `heaviside_function(-I_1)` → [ H(-\mathrm{tr}\,\boldsymbol{\epsilon}) ]
- `IxI` → [ \mathbf{I}\otimes\mathbf{I} ]
- `projector_positive` → [ \mathbb{P}^+ ] (from `positive_negative_projectors`)
- `projector_negative` → [ \mathbb{P}^- ]

**Code — fourth-order projection tensors** [ \mathbb{P}^+ ] **and** [ \mathbb{P}^- ]:

[`SpectrumDecomposition.h#L71-L155`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/SpectrumDecomposition.h#L71)

```cpp
  template <int dim>
  void positive_negative_projectors(Vector<double> const & eigenvalues,
                                    std::vector<Tensor<1, dim>> const & eigenvectors,
			            SymmetricTensor<4, dim> & positive_projector,
				    SymmetricTensor<4, dim> & negative_projector)
  {
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
```

[ \mathbb{P}^+ ] is assembled using the spectral formula:
[ (\mathbb{P}^+)_{ijkl} = \sum_a H(\epsilon_a)(M_a)_{ij}(M_a)_{kl} + \sum_{a\neq b} \frac{\langle\epsilon_a\rangle_+ - \langle\epsilon_b\rangle_+}{\epsilon_a - \epsilon_b} \cdot \tfrac{1}{2}(G_{ab}+G_{ba})_{ijkl} ]

The limit [ L'H\hat{o}pital ] (when [ |\epsilon_a - \epsilon_b| < 10^{-12} ]):
[ v_{ab} \to \tfrac{1}{2}(H(\epsilon_a) + H(\epsilon_b)) ]

---

### Phase-Field Viscosity Regularization (Eqs. 12–13)

**Eq. (12) — viscosity term added to phase-field balance:**

[ \eta\dot{d} \approx \frac{\eta}{\Delta t}(d^{(n+1)} - d^{(n)}) ]

**Eq. (13) — total energy with viscosity:**

[ \Pi_\eta(\mathbf{u}, d) = \Pi(\mathbf{u}, d) + \frac{\eta}{2\Delta t}\int_\Omega(d - d^{(n)})^2\,\mathrm{d}\Omega ]

**Code — viscosity in gradient (RHS):**

[`main.cc#L3813-L3819`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3813)

```cpp
    	        data.m_cell_rhs(i) += (    gc * length_scale * grad_N_phasefield[i] * phasefield_grad
    	                                +  (   gc / length_scale * phasefield_value
					     + eta / delta_time  * (phasefield_value - old_phasefield)
					     + degradation_function_derivative(phasefield_value)
					     * current_positive_strain_energy )
					  * N_phasefield[i]
				      ) * JxW;
```

`eta / delta_time * (phasefield_value - old_phasefield)` is [ \frac{\eta}{\Delta t}(d - d^{(n)}) ].

**Code — viscosity in Hessian (K_dd):**

[`main.cc#L3946-L3958`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3946)

```cpp
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

`eta/delta_time` contributes [ \frac{\eta}{\Delta t} ] to the second derivative of the viscosity
term with respect to [ d ].

---

### 2.2. Finite Element Discretization

#### Eq. (14): FE Discretization of Displacement

[ \mathbf{u}(\mathbf{x}) = \mathbf{N}_{u_A}(\mathbf{x})\,\mathbf{u}_A ]

**Code:**

[`main.cc#L3793-L3795`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3793)

```cpp
        const std::vector<Tensor<1,dim>> & N_disp = scratch.m_Nx_disp[q_point];
        const std::vector<SymmetricTensor<2, dim>> & symm_grad_N_disp =
          scratch.m_symm_grad_Nx_disp[q_point];
```

`N_disp[i]` is [ \mathbf{N}_{u_A} ] (vector-valued shape function for node/DOF [ A ]).
`symm_grad_N_disp[i]` is [ \nabla^{(s)}\mathbf{N}_{u_A} ] (symmetrized gradient).

---

#### Eq. (15): FE Discretization of Phase-Field

[ d(\mathbf{x}) = N_{d_A}(\mathbf{x})\,d_A ]

**Code:**

[`main.cc#L3787-L3789`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3787)

```cpp
        const std::vector<double>         &      N_phasefield = scratch.m_Nx_phasefield[q_point];
        const std::vector<Tensor<1, dim>> & grad_N_phasefield = scratch.m_grad_Nx_phasefield[q_point];
        const double                old_phasefield = scratch.m_phasefield_previous_step_cell[q_point];
```

`N_phasefield[i]` is [ N_{d_A} ] (scalar shape function), `grad_N_phasefield[i]` is [ \nabla N_{d_A} ].

---

#### Eq. (16): Discrete Residual/Gradient — Displacement Part

[ r_{\mathbf{u}_A} = \frac{\partial\Pi}{\partial\mathbf{u}_A} = (\nabla^{(s)}\mathbf{N}_{u_A},\,\boldsymbol{\sigma}) - (\mathbf{N}_{u_A},\,\mathbf{b}) - (\mathbf{N}_{u_A},\,\mathbf{t})_{\Gamma_t} ]

**Code:**

[`main.cc#L3804-L3810`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3804)

```cpp
            if (i_group == m_u_dof)
              {
                data.m_cell_rhs(i) += (symm_grad_N_disp[i] * cauchy_stress) * JxW;

		// contributions from the body force to right-hand side
		data.m_cell_rhs(i) -= N_disp[i] * rhs_values[q_point] * JxW;
              }
```

- `symm_grad_N_disp[i] * cauchy_stress` → [ \nabla^{(s)}\mathbf{N}_{u_A}:\boldsymbol{\sigma} ]
- `N_disp[i] * rhs_values[q_point]` → [ \mathbf{N}_{u_A}\cdot\mathbf{b} ]
- Neumann (traction) term is applied via constraints/boundary integration.

---

#### Eq. (17): Discrete Residual/Gradient — Phase-Field Part

[ r_{d_A} = \frac{\partial\Pi}{\partial d_A} = \left(N_{d_A},\,\frac{g_c}{l}d + g'(d)\psi^+\right) + (\nabla N_{d_A},\,g_c l\,\nabla d) ]

**Code:**

[`main.cc#L3811-L3820`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3811)

```cpp
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

- `gc * length_scale * grad_N * phasefield_grad` → [ g_c l\,\nabla N_{d_A}\cdot\nabla d ]
- `gc / length_scale * phasefield_value * N` → [ \frac{g_c}{l}d\,N_{d_A} ]
- `degradation_function_derivative(d) * psi_plus * N` → [ g'(d)\psi^+\,N_{d_A} ]
- `eta / delta_time * (d - d_old) * N` → viscosity term [ \frac{\eta}{\Delta t}(d-d^{(n)}) N_{d_A} ]

---

#### Eq. (18): Discrete Hessian (Stiffness) Matrix

**K_uu block:**

[ K_{u_A u_B} = \left(\nabla^{(s)}\mathbf{N}_{u_A},\,\mathbf{C}:\nabla^{(s)}\mathbf{N}_{u_B}\right) ]

**K_dd block:**

[ K_{d_A d_B} = \left(N_{d_A},\,\left(\frac{g_c}{l} + g''(d)\psi^+\right)N_{d_B}\right) + (\nabla N_{d_A},\,g_c l\,\nabla N_{d_B}) ]

**Code — B0 matrix (only diagonal blocks assembled):**

[`main.cc#L3940-L3964`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3940)

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
                  Assert((i_group <= m_d_dof) && (j_group <= m_d_dof),
                         ExcInternalError());
```

For `K_uu`: `symm_grad_Nx_i_x_C = C : symm_grad_N[i]`, so
`symm_grad_Nx_i_x_C * symm_grad_N_disp[j]` = [ \nabla^{(s)}\mathbf{N}_{u_A}:\mathbf{C}:\nabla^{(s)}\mathbf{N}_{u_B} ].
For `K_dd`: each term maps directly to the formula terms.
Note the off-diagonal blocks [ K_{ud} ] and [ K_{du} ] are **not** assembled in B0 (they are set to zero in the block structure).

---

## 3. Gradient Projection Based Monolithic Scheme

### 3.1. Algorithm Overview

#### Eq. (19): Projected Gradient (Convergence Criterion)

[ P_i = \mathrm{Proj}_{[\mathrm{lb}_i,\,\mathrm{ub}_i]}(x_i - r_i) - x_i = \begin{cases}\mathrm{lb}_i - x_i & x_i - r_i < \mathrm{lb}_i,\\ -r_i & \mathrm{lb}_i \leq x_i - r_i \leq \mathrm{ub}_i,\\ \mathrm{ub}_i - x_i & x_i - r_i > \mathrm{ub}_i.\end{cases} ]

**Code:**

[`main.cc#L1460-L1500`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L1460)

```cpp
  template <int dim>
  void PhaseFieldMonolithicSolve<dim>::get_error_residual_LBFGSB(Errors &error_residual,
								 const BlockVector<double> & solution_delta)
  {
    BlockVector<double> error_res(m_dofs_per_block);

    for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
      {
        if (!m_constraints.is_constrained(i))
	  error_res.block(m_u_dof)[i] = m_system_rhs.block(m_u_dof)[i];
      }

    const double upper_limit = 1.0;
    Vector<double> solution_phasefield_total(m_solution.block(m_d_dof));
    solution_phasefield_total += solution_delta.block(m_d_dof);

    double trial_solution = 0.0;
    for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
      {
        if (!m_constraints.is_constrained(i + m_dofs_per_block[m_u_dof]))
          {
            trial_solution = solution_phasefield_total(i) - m_system_rhs.block(m_d_dof)[i];

            if (trial_solution < m_solution.block(m_d_dof)[i])
              error_res.block(m_d_dof)[i] = m_solution.block(m_d_dof)[i] - solution_phasefield_total(i);
            else if (trial_solution > upper_limit)
              error_res.block(m_d_dof)[i] = upper_limit - solution_phasefield_total(i);
            else
              error_res.block(m_d_dof)[i] = (-m_system_rhs.block(m_d_dof)[i]);
          }
      }

    error_residual.m_norm = error_res.l2_norm();
    error_residual.m_u    = error_res.block(m_u_dof).l2_norm();
    error_residual.m_d    = error_res.block(m_d_dof).l2_norm();
  }
```

`trial_solution = x_i - r_i`. Three cases:
- `< lb_i`: `error_res[i] = lb_i - x_i` → [ \mathrm{lb}_i - x_i ]
- `> ub_i`: `error_res[i] = ub_i - x_i` → [ \mathrm{ub}_i - x_i ]
- else: `error_res[i] = -r_i` → [ -r_i ]

This matches Eq. (19) exactly. The [ l_2 ]-norm [ \|P\|_2 ] is used as the convergence test.

---

#### Convergence Criteria (Algorithm, Section 3.1)

1. Active set unchanged: [ \mathcal{A}_{k+1} = \mathcal{A}_k ]
2. Projected gradient: [ \|\mathbf{P}(\mathbf{x}_{k+1})\|_2 < \mathrm{tol\_residual} ]
3. Solution increment: [ \|\mathbf{x}_{k+1} - \mathbf{x}_k\|_2 < \mathrm{tol\_incr} ]

**Code:**

[`main.cc#L5033-L5040`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5033)

```cpp
        if (LBFGS_iteration > 1 && m_error_update_norm.m_u <= m_parameters.m_tol_u_incr
                                && m_error_residual_norm.m_u <= m_parameters.m_tol_u_residual
			        && m_error_update_norm.m_d <= m_parameters.m_tol_d_incr
			        && m_error_residual_norm.m_d <= m_parameters.m_tol_d_residual
				&& lower_bound_number_new == lower_bound_number_old
				&& upper_bound_number_new == upper_bound_number_old
				&& lowerupper_bound_number_new == lowerupper_bound_number_old)
```

All three criteria must hold simultaneously before declaring convergence.

---

### 3.2. Compact Representation of Limited-Memory BFGS Matrix

#### BFGS Update Vectors

[ \mathbf{s}_k = \mathbf{x}_{k+1} - \mathbf{x}_k, \quad \mathbf{y}_k = \mathbf{r}_{k+1} - \mathbf{r}_k ]

**Code — s and y computed and stored:**

[`main.cc#L5704-L5728`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5704)

```cpp
        LBFGS_y_vector = m_system_rhs;
        LBFGS_y_vector *= -1.0;
        assemble_system_rhs_BFGS_parallel(m_solution, m_system_rhs);
        LBFGS_y_vector += m_system_rhs;

        LBFGS_s_vector = LBFGS_update;

	double s_dot_y = LBFGS_s_vector * LBFGS_y_vector;
	if (s_dot_y > 1.0e-16 * LBFGS_y_vector.norm_sqr())
	  {
	    if (list_size >= LBFGS_m)
	      {
		s_vector_list.pop_front();
		y_vector_list.pop_front();
		s_dot_y_list.pop_front();
	      }

	    s_vector_list.push_back(LBFGS_s_vector);
	    y_vector_list.push_back(LBFGS_y_vector);
	    s_dot_y_list.push_back(s_dot_y);
	  }
```

- `LBFGS_y_vector = r_{k+1} - r_k` = [ \mathbf{y}_k ]
- `LBFGS_s_vector = LBFGS_update` (step) = [ \mathbf{s}_k ]
- Curvature condition [ \mathbf{s}_k^T\mathbf{y}_k > 0 ] is enforced before adding to list.
- `std::list` is a bounded FIFO: oldest entry popped when `list_size >= m`.

---

#### BFGS Matrix Update

[ \mathbf{B}_{k+1} = \mathbf{B}_k - \frac{\mathbf{B}_k\mathbf{s}_k\mathbf{s}_k^T\mathbf{B}_k}{\mathbf{s}_k^T\mathbf{B}_k\mathbf{s}_k} + \frac{\mathbf{y}_k\mathbf{y}_k^T}{\mathbf{y}_k^T\mathbf{s}_k} ]

This update is **never formed explicitly**; only the compact representation (Eq. 20) is used.

---

#### Curvature Condition

[ \mathbf{s}_k^T\mathbf{y}_k > 0 ]

**Code:** `if (s_dot_y > 1.0e-16 * LBFGS_y_vector.norm_sqr())` at
[`main.cc#L5712`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5712).

---

#### Correction Matrices S and Y

[ \mathbf{S}_k = [\mathbf{s}_{k-m}\;\cdots\;\mathbf{s}_{k-1}], \quad \mathbf{Y}_k = [\mathbf{y}_{k-m}\;\cdots\;\mathbf{y}_{k-1}] ]

**Code:**

[`main.cc#L4982-L4985`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4982)

```cpp
    std::list<BlockVector<double>> s_vector_list;
    std::list<BlockVector<double>> y_vector_list;
    std::list<double> s_dot_y_list;
    std::list<BlockVector<double>> b0xs_vector_list;
```

`s_vector_list` and `y_vector_list` store the columns of [ \mathbf{S}_k ] and [ \mathbf{Y}_k ]
as linked lists with at most `m` entries.

---

#### Eq. (20): Compact L-BFGS Representation

[ \mathbf{B}_k = \mathbf{B}_k^0 - \mathbf{W}_k\mathbf{M}_k\mathbf{W}_k^T ]

The matrix-vector product [ \mathbf{B}_k\mathbf{v} ] is computed implicitly
(see full code block quoted in Part 2 at [`main.cc#L5221-L5258`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5221)).

**Code — key multiply steps:**

[`main.cc#L5221-L5258`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5221)

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

	// temp_vector_2 = B_0 * (x^c - x_k) - W * M * W^T * (x^c - x_k)
	if (list_size > 0)
	  temp_vector_2 -= temp_vector_5;
```

Step by step: (1) [ \mathbf{B}^0_k\mathbf{v} ], (2) [ \mathbf{W}_k^T\mathbf{v} ],
(3) [ \mathbf{M}_k\mathbf{W}_k^T\mathbf{v} ], (4) [ \mathbf{W}_k\mathbf{M}_k\mathbf{W}_k^T\mathbf{v} ],
(5) subtract → [ \mathbf{B}_k\mathbf{v} ].

---

#### Eq. (21): W Matrix Definition

[ \mathbf{W}_k = [\mathbf{Y}_k \quad \mathbf{B}_k^0\mathbf{S}_k] \in \mathbb{R}^{n\times 2m} ]

**Code — W stored as two separate lists (never as one dense matrix):**

[`main.cc#L5119-L5126`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5119)

```cpp
	// B^0_k * s_vector has to be completely recalculated from scratch
	// at each L-BFGS iteration, since B^0_k is different
	b0xs_vector_list.clear();
	for (auto itr = s_vector_list.begin(); itr != s_vector_list.end(); ++itr)
	  {
	    m_tangent_matrix.vmult(b0xs_vector, *itr);
	    b0xs_vector_list.push_back(b0xs_vector);
	  }
```

- `y_vector_list` → columns of [ \mathbf{Y}_k ] (first [ m ] columns of [ \mathbf{W}_k ])
- `b0xs_vector_list` → [ \mathbf{B}^0_k\mathbf{s}_i ] (last [ m ] columns of [ \mathbf{W}_k ])
- Because [ \mathbf{B}^0_k ] is reassembled every iteration, [ \mathbf{B}^0_k\mathbf{S}_k ] is recomputed each time.

---

#### Eq. (22): M^{-1} Matrix Definition

[ \mathbf{M}_k^{-1} = \begin{bmatrix} -\mathbf{D}_k & \mathbf{L}_k^T \\ \mathbf{L}_k & \mathbf{S}_k^T\mathbf{B}_k^0\mathbf{S}_k \end{bmatrix} ]

where [ (\mathbf{D}_k)_{ii} = \mathbf{s}_i^T\mathbf{y}_i ] and
[ (\mathbf{L}_k)_{ij} = \mathbf{s}_i^T\mathbf{y}_j ] for [ i > j ],  0 otherwise.

**Code:**

[`main.cc#L5140-L5182`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5140)

```cpp
	FullMatrix<double> sTxBxs_matrix(list_size);
	sTxBxs_matrix = 0;
	for (unsigned int i = 0; i < list_size; ++i)
	  for (unsigned int j = 0; j <= i; ++j)
	    {
	      sTxBxs_matrix(i, j) = (*std::next(itr_s_begin,    i))
		                  * (*std::next(itr_b0xs_begin, j));
	    }
	for (unsigned int i = 0; i < list_size; ++i)
	  for (unsigned int j = i + 1; j < list_size; ++j)
	    {
	      sTxBxs_matrix(i, j) = sTxBxs_matrix(j, i);
	    }

	FullMatrix<double> D_matrix(list_size);
	D_matrix = 0;
	for (unsigned int i = 0; i < list_size; ++i)
	  D_matrix(i, i) = (*std::next(itr_s_dot_y_begin, i));

	FullMatrix<double> L_matrix(list_size);
	L_matrix = 0;
	for (unsigned int i = 0; i < list_size; ++i)
	  for (unsigned int j = 0; j < i; ++j)
	    L_matrix(i, j) = (*std::next(itr_s_begin, i))
                           * (*std::next(itr_y_begin, j));

	FullMatrix<double> M_matrix_inv(2 * list_size);
	FullMatrix<double> M_matrix(2 * list_size);

	M_matrix_inv = 0;
	for (unsigned int i = 0; i < list_size; ++i)
	  M_matrix_inv(i, i) = -D_matrix(i, i);

	for (unsigned int i = 0; i < list_size; ++i)
          for (unsigned int j = 0; j < list_size; ++j)
            {
              M_matrix_inv(i + list_size, j + list_size) = sTxBxs_matrix(i, j);
              M_matrix_inv(i + list_size, j            ) = L_matrix(i, j);
              M_matrix_inv(i            , j + list_size) = L_matrix(j, i);
            }

	if (!M_matrix_inv.empty())
	  M_matrix.invert(M_matrix_inv);
```

Block-by-block construction:
- `D_matrix(i,i) = s_dot_y_list[i]` → [ \mathbf{D}_k = \mathrm{diag}\{\mathbf{s}_i^T\mathbf{y}_i\} ]
- `L_matrix(i,j) = s_i * y_j` ([ i > j ]) → [ (\mathbf{L}_k)_{ij} ]
- `sTxBxs_matrix(i,j)` → [ \mathbf{S}_k^T\mathbf{B}_k^0\mathbf{S}_k ]
- `M_matrix_inv` assembles the [ 2m\times 2m ] block matrix [ \mathbf{M}_k^{-1} ]
- `M_matrix.invert(M_matrix_inv)` computes [ \mathbf{M}_k = (\mathbf{M}_k^{-1})^{-1} ]

---

#### Initial B0 Matrix (Eq. 11 / Eq. 21 in paper)

[ \mathbf{B}_k^0 = \hat{\mathbf{K}}^{(k)} = \begin{bmatrix}\mathbf{K}_{uu} & \mathbf{0} \\ \mathbf{0} & \mathbf{K}_{dd}\end{bmatrix}^{(k)} ]

**Code — assembly call:**

[`main.cc#L5113-L5117`](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5113)

```cpp
	// assemble the initial B_0 matrix at the k-th L-BFGS iteration
	// m_solution is the old solution from the previous converged step
	// it is needed only for the viscosity term
	// the output is m_tangent_matrix (B^0_k)
	assemble_system_B0(m_solution);
```

`assemble_system_B0` builds only the diagonal blocks [ \mathbf{K}_{uu} ] and [ \mathbf{K}_{dd} ];
the off-diagonal blocks are identically zero (not assembled), giving the exact block-diagonal
structure of [ \hat{\mathbf{K}} ].

---

*Continued in `formula_code_mapping_2.md` — covering Sections 3.3 Cauchy point,*
*3.4 subspace minimization, line-search Wolfe conditions (Eq. 25), Eq. (33) crack*
*dissipation energy, and the torsion boundary conditions.*
