## 3. Gradient projection based monolithic scheme > 3.1. Algorithm overview

### 未编号公式 #38 — 未编号公式：\pmb{p}_k = \pmb{x}^* - \p...
- 论文位置：[3. Gradient projection based monolithic scheme > 3.1. Algorithm overview，`explain.md` 第 309 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L309)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\pmb{p}_k = \pmb{x}^* - \pmb{x}_k.
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
### 未编号公式 #39 — 未编号公式：\pmb{x}_{k + 1} = \pmb{x}_...
- 论文位置：[3. Gradient projection based monolithic scheme > 3.1. Algorithm overview，`explain.md` 第 315 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L315)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\pmb{x}_{k + 1} = \pmb{x}_k + \alpha_k\pmb{p}_k.
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - 强 Wolfe 线搜索主流程：[main.cc:4196-4256](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4196-L4256)
```cpp
4196:   double PhaseFieldMonolithicSolve<dim>::line_search_stepsize_strong_wolfe(const double phi_0,
4197: 				                                           const double phi_0_prime,
4198: 				                                           const BlockVector<double> & BFGS_p_vector,
4199: 				                                           const BlockVector<double> & solution_delta)
4200:   {
4201:     //AssertThrow(phi_0_prime < 0,
4202:     //            ExcMessage("The derivative of phi at alpha = 0 should be negative!"));
4203: 
4204:     // Some line search parameters
4205:     const double c1 = 0.0001;
4206:     const double c2 = 0.9;
4207:     const double alpha_max = 1.0;
4208:     const unsigned int max_iter = 20;
4209:     double alpha = 1.0;
4210: 
4211:     double phi_old = phi_0;
4212:     double phi_prime_old = phi_0_prime;
4213:     double alpha_old = 0.0;
4214: 
4215:     double phi, phi_prime;
4216: 
4217:     std::pair<double, double> current_phi_phi_prime;
4218: 
4219:     unsigned int i = 0;
4220:     for (; i < max_iter; ++i)
4221:       {
4222: 	current_phi_phi_prime = calculate_phi_and_phi_prime(alpha, BFGS_p_vector, solution_delta);
4223: 	phi = current_phi_phi_prime.first;
4224: 	phi_prime = current_phi_phi_prime.second;
4225: 
4226: 	if (   ( phi > (phi_0 + c1 * alpha * phi_0_prime) )
4227: 	    || ( i > 0 && phi > phi_old ) )
4228: 	  {
4229: 	    return line_search_zoom_strong_wolfe(phi_old, phi_prime_old, alpha_old,
4230: 						 phi,     phi_prime,     alpha,
4231: 						 phi_0,   phi_0_prime,   BFGS_p_vector,
4232: 						 c1,      c2,            max_iter, solution_delta);
4233: 	  }
4234: 
4235: 	if (std::fabs(phi_prime) <= c2 * std::fabs(phi_0_prime))
4236: 	  {
4237: 	    return alpha;
4238: 	  }
4239: 
4240: 	if (phi_prime >= 0)
4241: 	  {
4242: 	    return line_search_zoom_strong_wolfe(phi,     phi_prime,     alpha,
4243: 						 phi_old, phi_prime_old, alpha_old,
4244: 						 phi_0,   phi_0_prime,   BFGS_p_vector,
4245: 						 c1,      c2,            max_iter, solution_delta);
4246: 	  }
4247: 
4248: 	phi_old = phi;
4249: 	phi_prime_old = phi_prime;
4250: 	alpha_old = alpha;
4251: 
4252: 	alpha = std::min(0.6*alpha, alpha_max);
4253:       }
4254: 
4255:     return alpha;
4256:   }
```
  - zoom 子过程：[main.cc:4259-4311](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4259-L4311)
```cpp
4259:   double PhaseFieldMonolithicSolve<dim>::
4260:     line_search_zoom_strong_wolfe(double phi_low, double phi_low_prime, double alpha_low,
4261: 				  double phi_high, double phi_high_prime, double alpha_high,
4262: 				  double phi_0, double phi_0_prime, const BlockVector<double> & BFGS_p_vector,
4263: 				  double c1, double c2, unsigned int max_iter, const BlockVector<double> & solution_delta)
4264:   {
4265:     double alpha = 0;
4266:     std::pair<double, double> current_phi_phi_prime;
4267:     double phi, phi_prime;
4268: 
4269:     unsigned int i = 0;
4270:     for (; i < max_iter; ++i)
4271:       {
4272: 	// a simple bisection is faster than cubic interpolation
4273: 	alpha = 0.5 * (alpha_low + alpha_high);
4274: 	//alpha = line_search_interpolation_cubic(alpha_low, phi_low, phi_low_prime,
4275: 	//					alpha_high, phi_high, phi_high_prime);
4276: 	current_phi_phi_prime = calculate_phi_and_phi_prime(alpha, BFGS_p_vector, solution_delta);
4277: 	phi = current_phi_phi_prime.first;
4278: 	phi_prime = current_phi_phi_prime.second;
4279: 
4280: 	if (   (phi > phi_0 + c1 * alpha * phi_0_prime)
4281: 	    || (phi > phi_low) )
4282: 	  {
4283: 	    alpha_high = alpha;
4284: 	    phi_high = phi;
4285: 	    phi_high_prime = phi_prime;
4286: 	  }
4287: 	else
4288: 	  {
4289: 	    if (std::fabs(phi_prime) <= c2 * std::fabs(phi_0_prime))
4290: 	      {
4291: 		return alpha;
4292: 	      }
4293: 
4294: 	    if (phi_prime * (alpha_high - alpha_low) >= 0.0)
4295: 	      {
4296: 		alpha_high = alpha_low;
4297: 		phi_high_prime = phi_low_prime;
4298: 		phi_high = phi_low;
4299: 	      }
4300: 
4301: 	    alpha_low = alpha;
4302: 	    phi_low_prime = phi_prime;
4303: 	    phi_low = phi;
4304: 	  }
4305:       }
4306: 
4307:     // avoid unused variable warnings from compiler
4308:     (void)phi_high;
4309:     (void)phi_high_prime;
4310:     return alpha;
4311:   }
```
  - phi/phi' 计算（能量+方向导）：[main.cc:4348-4367](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4348-L4367)
```cpp
4348:   std::pair<double, double> PhaseFieldMonolithicSolve<dim>::
4349:     calculate_phi_and_phi_prime(const double alpha,
4350: 				const BlockVector<double> & BFGS_p_vector,
4351: 				const BlockVector<double> & solution_delta)
4352:   {
4353:     // the first component is phi(alpha), the second component is phi_prime(alpha),
4354:     std::pair<double, double> phi_values;
4355: 
4356:     BlockVector<double> solution_delta_trial(solution_delta);
4357:     solution_delta_trial.add(alpha, BFGS_p_vector);
4358: 
4359:     update_qph_incremental(solution_delta_trial, m_solution);
4360: 
4361:     BlockVector<double> system_rhs(m_dofs_per_block);
4362:     assemble_system_rhs_BFGS_parallel(m_solution, system_rhs);
4363:     //m_constraints.condense(system_rhs);
4364: 
4365:     phi_values.first = calculate_energy_functional();
4366:     phi_values.second = system_rhs * BFGS_p_vector;
4367:     return phi_values;
```
### 未编号公式 #40 — 未编号公式：\Pi_{k + 1} = \Pi (\pmb{x}...
- 论文位置：[3. Gradient projection based monolithic scheme > 3.1. Algorithm overview，`explain.md` 第 321 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L321)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\Pi_{k + 1} = \Pi (\pmb{x}_k + \alpha_k\pmb {p}_k) \leq \Pi_k + c_1\alpha_k\pmb {p}_k^\mathrm{T}\pmb{r}_k
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - 裂纹耗散能密度（含梯度项）：[main.cc:886-891](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L886-L891)
```cpp
886:     m_crack_energy_dissipation = m_gc * (  0.5 / m_length_scale * m_phase_field_value * m_phase_field_value
887: 	                                   + 0.5 * m_length_scale * m_grad_phasefield * m_grad_phasefield)
888: 	                                   // the term due to viscosity regularization
889: 	                                   + (m_phase_field_value - phase_field_value_previous_step)
890: 					   * (m_phase_field_value - phase_field_value_previous_step)
891: 				           * 0.5 * m_eta / delta_time;
```
  - 总能量/裂纹能积分计算：[main.cc:6117-6167](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L6117-L6167)
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
6141: 
6142:   template <int dim>
6143:   std::pair<double, double>
6144:     PhaseFieldMonolithicSolve<dim>::calculate_total_strain_energy_and_crack_energy_dissipation() const
6145:   {
6146:     double total_strain_energy = 0.0;
6147:     double crack_energy_dissipation = 0.0;
6148: 
6149:     FEValues<dim> fe_values(m_fe, m_qf_cell, update_JxW_values);
6150: 
6151:     for (const auto &cell : m_dof_handler.active_cell_iterators())
6152:       {
6153:         fe_values.reinit(cell);
6154: 
6155:         const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
6156:           m_quadrature_point_history.get_data(cell);
6157:         Assert(lqph.size() == m_n_q_points, ExcInternalError());
6158: 
6159:         for (unsigned int q_point = 0; q_point < m_n_q_points; ++q_point)
6160:           {
6161:             const double JxW = fe_values.JxW(q_point);
6162:             total_strain_energy += lqph[q_point]->get_total_strain_energy() * JxW;
6163:             crack_energy_dissipation += lqph[q_point]->get_crack_energy_dissipation() * JxW;
6164:           }
6165:       }
6166: 
6167:     return std::make_pair(total_strain_energy, crack_energy_dissipation);
```
  - 强 Wolfe 线搜索主流程：[main.cc:4196-4256](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4196-L4256)
```cpp
4196:   double PhaseFieldMonolithicSolve<dim>::line_search_stepsize_strong_wolfe(const double phi_0,
4197: 				                                           const double phi_0_prime,
4198: 				                                           const BlockVector<double> & BFGS_p_vector,
4199: 				                                           const BlockVector<double> & solution_delta)
4200:   {
4201:     //AssertThrow(phi_0_prime < 0,
4202:     //            ExcMessage("The derivative of phi at alpha = 0 should be negative!"));
4203: 
4204:     // Some line search parameters
4205:     const double c1 = 0.0001;
4206:     const double c2 = 0.9;
4207:     const double alpha_max = 1.0;
4208:     const unsigned int max_iter = 20;
4209:     double alpha = 1.0;
4210: 
4211:     double phi_old = phi_0;
4212:     double phi_prime_old = phi_0_prime;
4213:     double alpha_old = 0.0;
4214: 
4215:     double phi, phi_prime;
4216: 
4217:     std::pair<double, double> current_phi_phi_prime;
4218: 
4219:     unsigned int i = 0;
4220:     for (; i < max_iter; ++i)
4221:       {
4222: 	current_phi_phi_prime = calculate_phi_and_phi_prime(alpha, BFGS_p_vector, solution_delta);
4223: 	phi = current_phi_phi_prime.first;
4224: 	phi_prime = current_phi_phi_prime.second;
4225: 
4226: 	if (   ( phi > (phi_0 + c1 * alpha * phi_0_prime) )
4227: 	    || ( i > 0 && phi > phi_old ) )
4228: 	  {
4229: 	    return line_search_zoom_strong_wolfe(phi_old, phi_prime_old, alpha_old,
4230: 						 phi,     phi_prime,     alpha,
4231: 						 phi_0,   phi_0_prime,   BFGS_p_vector,
4232: 						 c1,      c2,            max_iter, solution_delta);
4233: 	  }
4234: 
4235: 	if (std::fabs(phi_prime) <= c2 * std::fabs(phi_0_prime))
4236: 	  {
4237: 	    return alpha;
4238: 	  }
4239: 
4240: 	if (phi_prime >= 0)
4241: 	  {
4242: 	    return line_search_zoom_strong_wolfe(phi,     phi_prime,     alpha,
4243: 						 phi_old, phi_prime_old, alpha_old,
4244: 						 phi_0,   phi_0_prime,   BFGS_p_vector,
4245: 						 c1,      c2,            max_iter, solution_delta);
4246: 	  }
4247: 
4248: 	phi_old = phi;
4249: 	phi_prime_old = phi_prime;
4250: 	alpha_old = alpha;
4251: 
4252: 	alpha = std::min(0.6*alpha, alpha_max);
4253:       }
4254: 
4255:     return alpha;
4256:   }
```
### 未编号公式 #41 — 未编号公式：\left|\pmb{r}_{k + 1}^\mat...
- 论文位置：[3. Gradient projection based monolithic scheme > 3.1. Algorithm overview，`explain.md` 第 327 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L327)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\left|\pmb{r}_{k + 1}^\mathrm{T}\pmb {p}_k\right| = \left|\pmb {r}(\pmb{x}_k + \alpha_k\pmb {p}_k)^\mathrm{T}\pmb {p}_k\right| \leq c_2\left|\pmb{r}_k^\mathrm{T}\pmb {p}_k\right|.
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - 强 Wolfe 线搜索主流程：[main.cc:4196-4256](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4196-L4256)
```cpp
4196:   double PhaseFieldMonolithicSolve<dim>::line_search_stepsize_strong_wolfe(const double phi_0,
4197: 				                                           const double phi_0_prime,
4198: 				                                           const BlockVector<double> & BFGS_p_vector,
4199: 				                                           const BlockVector<double> & solution_delta)
4200:   {
4201:     //AssertThrow(phi_0_prime < 0,
4202:     //            ExcMessage("The derivative of phi at alpha = 0 should be negative!"));
4203: 
4204:     // Some line search parameters
4205:     const double c1 = 0.0001;
4206:     const double c2 = 0.9;
4207:     const double alpha_max = 1.0;
4208:     const unsigned int max_iter = 20;
4209:     double alpha = 1.0;
4210: 
4211:     double phi_old = phi_0;
4212:     double phi_prime_old = phi_0_prime;
4213:     double alpha_old = 0.0;
4214: 
4215:     double phi, phi_prime;
4216: 
4217:     std::pair<double, double> current_phi_phi_prime;
4218: 
4219:     unsigned int i = 0;
4220:     for (; i < max_iter; ++i)
4221:       {
4222: 	current_phi_phi_prime = calculate_phi_and_phi_prime(alpha, BFGS_p_vector, solution_delta);
4223: 	phi = current_phi_phi_prime.first;
4224: 	phi_prime = current_phi_phi_prime.second;
4225: 
4226: 	if (   ( phi > (phi_0 + c1 * alpha * phi_0_prime) )
4227: 	    || ( i > 0 && phi > phi_old ) )
4228: 	  {
4229: 	    return line_search_zoom_strong_wolfe(phi_old, phi_prime_old, alpha_old,
4230: 						 phi,     phi_prime,     alpha,
4231: 						 phi_0,   phi_0_prime,   BFGS_p_vector,
4232: 						 c1,      c2,            max_iter, solution_delta);
4233: 	  }
4234: 
4235: 	if (std::fabs(phi_prime) <= c2 * std::fabs(phi_0_prime))
4236: 	  {
4237: 	    return alpha;
4238: 	  }
4239: 
4240: 	if (phi_prime >= 0)
4241: 	  {
4242: 	    return line_search_zoom_strong_wolfe(phi,     phi_prime,     alpha,
4243: 						 phi_old, phi_prime_old, alpha_old,
4244: 						 phi_0,   phi_0_prime,   BFGS_p_vector,
4245: 						 c1,      c2,            max_iter, solution_delta);
4246: 	  }
4247: 
4248: 	phi_old = phi;
4249: 	phi_prime_old = phi_prime;
4250: 	alpha_old = alpha;
4251: 
4252: 	alpha = std::min(0.6*alpha, alpha_max);
4253:       }
4254: 
4255:     return alpha;
4256:   }
```
  - zoom 子过程：[main.cc:4259-4311](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4259-L4311)
```cpp
4259:   double PhaseFieldMonolithicSolve<dim>::
4260:     line_search_zoom_strong_wolfe(double phi_low, double phi_low_prime, double alpha_low,
4261: 				  double phi_high, double phi_high_prime, double alpha_high,
4262: 				  double phi_0, double phi_0_prime, const BlockVector<double> & BFGS_p_vector,
4263: 				  double c1, double c2, unsigned int max_iter, const BlockVector<double> & solution_delta)
4264:   {
4265:     double alpha = 0;
4266:     std::pair<double, double> current_phi_phi_prime;
4267:     double phi, phi_prime;
4268: 
4269:     unsigned int i = 0;
4270:     for (; i < max_iter; ++i)
4271:       {
4272: 	// a simple bisection is faster than cubic interpolation
4273: 	alpha = 0.5 * (alpha_low + alpha_high);
4274: 	//alpha = line_search_interpolation_cubic(alpha_low, phi_low, phi_low_prime,
4275: 	//					alpha_high, phi_high, phi_high_prime);
4276: 	current_phi_phi_prime = calculate_phi_and_phi_prime(alpha, BFGS_p_vector, solution_delta);
4277: 	phi = current_phi_phi_prime.first;
4278: 	phi_prime = current_phi_phi_prime.second;
4279: 
4280: 	if (   (phi > phi_0 + c1 * alpha * phi_0_prime)
4281: 	    || (phi > phi_low) )
4282: 	  {
4283: 	    alpha_high = alpha;
4284: 	    phi_high = phi;
4285: 	    phi_high_prime = phi_prime;
4286: 	  }
4287: 	else
4288: 	  {
4289: 	    if (std::fabs(phi_prime) <= c2 * std::fabs(phi_0_prime))
4290: 	      {
4291: 		return alpha;
4292: 	      }
4293: 
4294: 	    if (phi_prime * (alpha_high - alpha_low) >= 0.0)
4295: 	      {
4296: 		alpha_high = alpha_low;
4297: 		phi_high_prime = phi_low_prime;
4298: 		phi_high = phi_low;
4299: 	      }
4300: 
4301: 	    alpha_low = alpha;
4302: 	    phi_low_prime = phi_prime;
4303: 	    phi_low = phi;
4304: 	  }
4305:       }
4306: 
4307:     // avoid unused variable warnings from compiler
4308:     (void)phi_high;
4309:     (void)phi_high_prime;
4310:     return alpha;
4311:   }
```
  - phi/phi' 计算（能量+方向导）：[main.cc:4348-4367](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4348-L4367)
```cpp
4348:   std::pair<double, double> PhaseFieldMonolithicSolve<dim>::
4349:     calculate_phi_and_phi_prime(const double alpha,
4350: 				const BlockVector<double> & BFGS_p_vector,
4351: 				const BlockVector<double> & solution_delta)
4352:   {
4353:     // the first component is phi(alpha), the second component is phi_prime(alpha),
4354:     std::pair<double, double> phi_values;
4355: 
4356:     BlockVector<double> solution_delta_trial(solution_delta);
4357:     solution_delta_trial.add(alpha, BFGS_p_vector);
4358: 
4359:     update_qph_incremental(solution_delta_trial, m_solution);
4360: 
4361:     BlockVector<double> system_rhs(m_dofs_per_block);
4362:     assemble_system_rhs_BFGS_parallel(m_solution, system_rhs);
4363:     //m_constraints.condense(system_rhs);
4364: 
4365:     phi_values.first = calculate_energy_functional();
4366:     phi_values.second = system_rhs * BFGS_p_vector;
4367:     return phi_values;
```
### 未编号公式 #42 — 未编号公式：\mathcal{A}_{k + 1}(\pmb{x...
- 论文位置：[3. Gradient projection based monolithic scheme > 3.1. Algorithm overview，`explain.md` 第 336 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L336)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\mathcal{A}_{k + 1}(\pmb{x}^c) = \mathcal{A}_k(\pmb{x}^c).
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - 相场点投影（box 约束）：[main.cc:1379-1394](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L1379-L1394)
```cpp
1379:   void PhaseFieldMonolithicSolve<dim>::point_projection(BlockVector<double> & solution_delta)
1380:   {
1381:     // Phase-field value cannot exceed 1.0
1382:     const double upper_limit = 1.0;
1383: 
1384:     Vector<double> solution_phasefield_total(m_solution.block(m_d_dof));
1385:     solution_phasefield_total += solution_delta.block(m_d_dof);
1386: 
1387:     for (unsigned int i = 0; i < solution_phasefield_total.size(); ++i)
1388:       {
1389: 	if (solution_delta.block(m_d_dof)[i] < 0.0)
1390: 	  solution_delta.block(m_d_dof)[i] = 0.0;
1391: 
1392: 	if (solution_phasefield_total[i] > upper_limit)
1393: 	  solution_delta.block(m_d_dof)[i] = upper_limit - m_solution.block(m_d_dof)[i];
1394:       }
```
  - break points 计算与初始活跃集：[main.cc:1401-1443](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L1401-L1443)
```cpp
1401:     PhaseFieldMonolithicSolve<dim>::calculate_break_points(const BlockVector<double> & solution_delta,
1402: 		       				           const BlockVector<double> & gradient_g,
1403: 							   BlockVector<double> & gradient_d)
1404:   {
1405:     // Creates a min heap of break points
1406:     std::priority_queue< std::pair<double, unsigned int>,
1407:                          std::vector<std::pair<double, unsigned int>>,
1408:     		         std::greater<std::pair<double, unsigned int>> >
1409:     break_points_sorted;
1410: 
1411:     double t = 0.0;
1412: 
1413:     Vector<double> solution_phasefield_total(m_solution.block(m_d_dof));
1414:     solution_phasefield_total += solution_delta.block(m_d_dof);
1415: 
1416:     // upper bound is 1.0, lower bound is the solution at the previous step.
1417:     for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
1418:       {
1419: 	if (gradient_g.block(m_d_dof)[i] < 0)
1420: 	  t = (solution_phasefield_total[i] - 1.0 ) / gradient_g.block(m_d_dof)[i];
1421: 	else if (gradient_g.block(m_d_dof)[i] > 0)
1422: 	  t = solution_delta.block(m_d_dof)[i] / gradient_g.block(m_d_dof)[i];
1423: 	else
1424: 	  t = std::numeric_limits<double>::max();
1425: 
1426:         //AssertThrow(t >= 0, ExcMessage("Break point has to be a non-negative t value"));
1427: 
1428:         if (t > 0)
1429:           {
1430: 	    break_points_sorted.push(std::make_pair(t, i));
1431:           }
1432:         else // if t == 0, i is in the active set
1433:           {
1434:             gradient_d.block(m_d_dof)[i] = 0;
1435:             if (gradient_g.block(m_d_dof)[i] > 0)
1436:               m_active_set_phasefield(i) = 1; //lower bound
1437:             else
1438:               m_active_set_phasefield(i) = 2; //upper bound
1439:           }
1440:       }
1441: 
1442:     return break_points_sorted;
1443:   }
```
  - 广义 Cauchy 点与活跃集更新：[main.cc:4774-4948](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L4774-L4948)
```cpp
4774:   void PhaseFieldMonolithicSolve<dim>::
4775:   calculate_cauchy_point(const BlockSparseMatrix<double> & B0_matrix,
4776: 	                 const std::list<BlockVector<double>> & y_vector_list,
4777: 		         const std::list<BlockVector<double>> & b0xs_vector_list,
4778: 			 const FullMatrix<double> & M_matrix,
4779: 			 const BlockVector<double> & gradient_g,
4780: 			 const BlockVector<double> & solution_delta,
4781: 			 BlockVector<double> & solution_delta_cauchy_point)
4782:   {
4783:     m_timer.enter_subsection("Calculate Cauchy point");
4784: 
4785:     solution_delta_cauchy_point = 0.0;
4786:     BlockVector<double> gradient_d(gradient_g);
4787:     gradient_d *= -1;
4788: 
4789:     const unsigned int list_size = y_vector_list.size();
4790:     const auto itr_y_begin    = y_vector_list.begin();
4791:     const auto itr_b0xs_begin = b0xs_vector_list.begin();
4792: 
4793:     // t_series only contains t > 0
4794:     std::priority_queue< std::pair<double, unsigned int>,
4795:                          std::vector<std::pair<double, unsigned int>>,
4796:         		 std::greater<std::pair<double, unsigned int>> >
4797:     t_series = calculate_break_points(solution_delta,
4798:     			              gradient_g,
4799: 				      gradient_d);
4800: 
4801:     // m_active_set_phasefield contains 1 or 2 for active set and 0 for inactive set
4802:     for (unsigned int i = 0; i < m_active_set_phasefield.size(); ++i)
4803:       {
4804: 	if (m_active_set_phasefield(i) > 0.5)
4805: 	  solution_delta_cauchy_point.block(m_d_dof)[i] = solution_delta.block(m_d_dof)[i];
4806:       }
4807: 
4808:     // p = W^T * d
4809:     Vector<double> p(2 * list_size);
4810:     for (unsigned int i = 0; i < list_size; ++i)
4811:       {
4812:         p(i)             = (*std::next(itr_y_begin,    i)) * gradient_d;
4813:         p(i + list_size) = (*std::next(itr_b0xs_begin, i)) * gradient_d;
4814:       }
4815: 
4816:     Vector<double> c(2 * list_size);
4817:     c = 0.0;
4818: 
4819:     double f_prime = -(gradient_d * gradient_d);
4820: 
4821:     // M * p
4822:     Vector<double> Mp(2 * list_size);
4823:     if (list_size > 0)
4824:       M_matrix.vmult(Mp, p);
4825: 
4826:     // B_0 * d
4827:     BlockVector<double> B0_grandient_d(m_dofs_per_block);
4828:     B0_matrix.vmult(B0_grandient_d, gradient_d);
4829: 
4830:     double f_prime_prime = gradient_d * B0_grandient_d;
4831:     if (list_size > 0)
4832:       f_prime_prime -= (p * Mp);
4833: 
4834:     double delta_t_min = -f_prime / f_prime_prime;
4835: 
4836:     double t_old = 0.0;
4837: 
4838:     std::pair<double, unsigned int> top_pair = t_series.top();
4839:     double t = top_pair.first;
4840:     unsigned int b = top_pair.second;
4841: 
4842:     double delta_t = t - t_old;
4843: 
4844:     BlockVector<double> z(m_dofs_per_block);
4845:     z = 0.0;
4846: 
4847:     // w_b = W^T * e_b
4848:     Vector<double> w_b(2 * list_size);
4849: 
4850:     // w_b_T_x_M = w_b^T * M
4851:     Vector<double> w_b_T_x_M(2 * list_size);
4852: 
4853:     Vector<double> temp_vector(m_dofs_per_block[m_d_dof]);
4854: 
4855:     while (delta_t_min >= delta_t)
4856:       {
4857: 	t_series.pop();
4858: 
4859: 	if (gradient_d.block(m_d_dof)[b] > 0)
4860: 	  solution_delta_cauchy_point.block(m_d_dof)[b] = 1.0 - m_solution.block(m_d_dof)[b];
4861: 	else if (gradient_d.block(m_d_dof)[b] < 0)
4862: 	  solution_delta_cauchy_point.block(m_d_dof)[b] = 0.0;
4863: 	else
4864: 	  AssertThrow(false,
4865: 	              ExcMessage("gradient_d(b) cannot be zero!"));
4866: 
4867: 	if (gradient_d.block(m_d_dof)[b] < 0)
4868: 	  m_active_set_phasefield[b] = 1; //lower bound
4869: 	else
4870: 	  m_active_set_phasefield[b] = 2; //upper bound
4871: 
4872:         // z_{j} = z_{j-1} + delta_t_{j-1} * gradient_d_{j-1};
4873: 	z.sadd(1.0, delta_t, gradient_d);
4874: 
4875: 	// c_{j} = c_{j-1} + delta_t_{j-1} * p_{j-1};
4876: 	if (list_size > 0)
4877: 	  c.sadd(1.0, delta_t, p);
4878: 
4879:         double temp_scalar = ebT_x_B0_x_v(b, B0_matrix, z);
4880: 
4881:         // w_b = W^T * e_b
4882:         for (unsigned int i = 0; i < list_size; ++i)
4883:           {
4884:             w_b(i)             = (*std::next(itr_y_begin,    i)).block(m_d_dof)[b];
4885:             w_b(i + list_size) = (*std::next(itr_b0xs_begin, i)).block(m_d_dof)[b];
4886:           }
4887: 
4888:         if (list_size > 0)
4889:           M_matrix.vmult(w_b_T_x_M, w_b);
4890: 
4891: 	f_prime += delta_t * f_prime_prime
4892: 	         + gradient_g.block(m_d_dof)[b] * gradient_g.block(m_d_dof)[b]
4893: 	         + temp_scalar * gradient_g.block(m_d_dof)[b];
4894: 
4895: 	if (list_size > 0)
4896: 	  f_prime -= (w_b_T_x_M * c) * gradient_g.block(m_d_dof)[b];
4897: 
4898: 	temp_scalar = ebT_x_B0_x_v(b, B0_matrix, gradient_d);
4899: 
4900: 	f_prime_prime += 2.0 * gradient_g.block(m_d_dof)[b] * temp_scalar
4901: 	               + gradient_g.block(m_d_dof)[b] * gradient_g.block(m_d_dof)[b]
4902: 		             * B0_matrix.block(m_d_dof, m_d_dof)(b, b);
4903: 
4904: 	if (list_size > 0)
4905: 	  {
4906: 	    f_prime_prime -= 2.0 * gradient_g.block(m_d_dof)[b] * (w_b_T_x_M * p);
4907: 	    f_prime_prime -= gradient_g.block(m_d_dof)[b] * gradient_g.block(m_d_dof)[b] * (w_b_T_x_M * w_b);
4908: 	  }
4909: 
4910: 	// p_{j} = p_{j-1} + g_b * w_b;
4911: 	if (list_size > 0)
4912: 	  p.sadd(1.0, gradient_g.block(m_d_dof)[b], w_b);
4913: 
4914: 	gradient_d.block(m_d_dof)[b] = 0.0;
4915: 
4916: 	delta_t_min = -f_prime / f_prime_prime;
4917: 
4918: 	t_old = t;
4919: 
4920: 	top_pair = t_series.top();
4921: 	t = top_pair.first;
4922: 	b = top_pair.second;
4923: 
4924: 	delta_t = t - t_old;
4925:       }
4926: 
4927:     if (delta_t_min < 0)
4928:       delta_t_min = 0;
4929: 
4930:     t_old += delta_t_min;
4931: 
4932:     for(unsigned int i = 0; i < m_active_set_phasefield.size(); ++i)
4933:       {
4934: 	// inactive phasefield dof
4935: 	if (m_active_set_phasefield(i) < 0.5)
4936: 	  solution_delta_cauchy_point.block(m_d_dof)[i] = solution_delta.block(m_d_dof)[i]
4937: 						+ t_old * gradient_d.block(m_d_dof)[i];
4938:       }
4939: 
4940:     // There are no active constraints in the displacement field
4941:     solution_delta_cauchy_point.block(m_u_dof) = solution_delta.block(m_u_dof);
4942:     (solution_delta_cauchy_point.block(m_u_dof)).add(t_old, gradient_d.block(m_u_dof));
4943: 
4944:     // We need to make sure the solution_delta_cauchy_point satisfies the essential
4945:     // boundary conditions and the hanging-node constraints
4946:     m_constraints.distribute(solution_delta_cauchy_point);
4947: 
4948:     m_timer.leave_subsection();
```
