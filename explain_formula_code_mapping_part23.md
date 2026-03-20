## 3. Gradient projection based monolithic scheme > 3.4. Subspace minimization > 3.4.3. Schur complement for the dual approach

### 未编号公式 #105 — 未编号公式：\mathbf{lb} \leq \mathbf{x...
- 论文位置：[3. Gradient projection based monolithic scheme > 3.4. Subspace minimization > 3.4.3. Schur complement for the dual approach，`explain.md` 第 732 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L732)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\mathbf{lb} \leq \mathbf{x}_k + \alpha_k \mathbf{p}_k \leq \mathbf{ub}.
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

## 4. Numerical examples > 4.1. Cyclic tension-compression test

### Eq. (33) — 裂纹耗散能
- 论文位置：[4. Numerical examples > 4.1. Cyclic tension-compression test，`explain.md` 第 780 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L780)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
E_{\mathrm{crack}} = g_{\mathrm{c}}\Gamma_{l}(d) = g_{\mathrm{c}}\int_{\Omega}\gamma (d,\nabla d)\mathrm{d}\Omega = g_{\mathrm{c}}\int_{\Omega}\frac{1}{2l}\left(d^{2} + l^{2}\nabla d\cdot \nabla d\right)\mathrm{d}\Omega \quad (33)
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

## 4. Numerical examples > 4.4. Three-dimensional torsion test

### 未编号公式 #107 — 未编号公式：u_{y} = z \tan t, \quad u_...
- 论文位置：[4. Numerical examples > 4.4. Three-dimensional torsion test，`explain.md` 第 878 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L878)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
u_{y} = z \tan t, \quad u_{z} = -y \tan t,
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - 3D 扭转边界条件离散施加：[main.cc:3592-3618](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3592-L3618)
```cpp
3592: 	    for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)
3593: 	      {
3594: 		if (std::fabs(vertex_itr->vertex()[0] - 0.0) < 1.0e-9)
3595: 		  {
3596: 		    node_rotate = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);
3597: 		    node_dist = std::sqrt(  vertex_itr->vertex()[1] * vertex_itr->vertex()[1]
3598: 			                  + vertex_itr->vertex()[2] * vertex_itr->vertex()[2]);
3599: 
3600: 		    angle_theta = m_time.get_delta_t() * m_time.get_magnitude();
3601: 		    disp_mag = node_dist * std::tan(angle_theta);
3602: 
3603: 		    if (node_dist > 0)
3604: 		      {
3605: 		        disp_y = vertex_itr->vertex()[2]/node_dist * disp_mag;
3606: 		        disp_z = -vertex_itr->vertex()[1]/node_dist * disp_mag;
3607: 		      }
3608: 		    else
3609: 		      {
3610: 			disp_y = 0.0;
3611: 			disp_z = 0.0;
3612: 		      }
3613: 
3614: 		    m_constraints.add_line(node_rotate[1]);
3615: 		    m_constraints.set_inhomogeneity(node_rotate[1], disp_y);
3616: 
3617: 		    m_constraints.add_line(node_rotate[2]);
3618: 		    m_constraints.set_inhomogeneity(node_rotate[2], disp_z);
```

## 5. Wall-clock time and convergence behavior > 5.2. Comparison of convergence behaviors

### 未编号公式 #108 — 未编号公式：\| \pmb {r}_u\| _2 < \math...
- 论文位置：[5. Wall-clock time and convergence behavior > 5.2. Comparison of convergence behaviors，`explain.md` 第 937 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L937)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\| \pmb {r}_u\| _2 < \mathrm{tol}, \quad \| \pmb {r}_d\| _2 < \mathrm{tol}, \quad \| \Delta \pmb {u}\| _2 < \mathrm{tol}, \quad \| \Delta d\| _2 < \mathrm{tol},
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - LBFGS-B 投影残差范数计算：[main.cc:1460-1500](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L1460-L1500)
```cpp
1460:   void PhaseFieldMonolithicSolve<dim>::get_error_residual_LBFGSB(Errors &error_residual,
1461: 								 const BlockVector<double> & solution_delta)
1462:   {
1463:     // We use L_2 norm
1464:     BlockVector<double> error_res(m_dofs_per_block);
1465: 
1466:     // For displacement DOFs, except essential boundary conditions
1467:     // and hanging-node constraints, there are no box constraints
1468:     for (unsigned int i = 0; i < m_dofs_per_block[m_u_dof]; ++i)
1469:       {
1470:         if (!m_constraints.is_constrained(i))
1471: 	  error_res.block(m_u_dof)[i] = m_system_rhs.block(m_u_dof)[i];
1472:       }
1473: 
1474:     // For phasefield DOFs, there are points with active box constraints
1475:     // and points with inactive box constraints
1476:     const double upper_limit = 1.0;
1477:     Vector<double> solution_phasefield_total(m_solution.block(m_d_dof));
1478:     solution_phasefield_total += solution_delta.block(m_d_dof);
1479: 
1480:     double trial_solution = 0.0;
1481:     for (unsigned int i = 0; i < m_dofs_per_block[m_d_dof]; ++i)
1482:       {
1483: 	// phasefield DOFs can still be constrained due to hanging-nodes
1484:         if (!m_constraints.is_constrained(i + m_dofs_per_block[m_u_dof]))
1485:           {
1486:             trial_solution = solution_phasefield_total(i) - m_system_rhs.block(m_d_dof)[i];
1487: 
1488:             if (trial_solution < m_solution.block(m_d_dof)[i])
1489:               error_res.block(m_d_dof)[i] = m_solution.block(m_d_dof)[i] - solution_phasefield_total(i);
1490:             else if (trial_solution > upper_limit)
1491:               error_res.block(m_d_dof)[i] = upper_limit - solution_phasefield_total(i);
1492:             else
1493:               error_res.block(m_d_dof)[i] = (-m_system_rhs.block(m_d_dof)[i]);
1494:           }
1495:       }
1496: 
1497:     error_residual.m_norm = error_res.l2_norm();
1498:     error_residual.m_u    = error_res.block(m_u_dof).l2_norm();
1499:     error_residual.m_d    = error_res.block(m_d_dof).l2_norm();
1500:   }
```
  - LBFGS-B 收敛判据（含活跃集稳定）：[main.cc:5033-5040](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L5033-L5040)
```cpp
5033:         if (LBFGS_iteration > 1 && m_error_update_norm.m_u <= m_parameters.m_tol_u_incr
5034:                                 && m_error_residual_norm.m_u <= m_parameters.m_tol_u_residual
5035: 			        && m_error_update_norm.m_d <= m_parameters.m_tol_d_incr
5036: 			        && m_error_residual_norm.m_d <= m_parameters.m_tol_d_residual
5037: 				&& lower_bound_number_new == lower_bound_number_old
5038: 				&& upper_bound_number_new == upper_bound_number_old
5039: 				&& lowerupper_bound_number_new == lowerupper_bound_number_old)
5040:           {
```

## Appendix

### 未编号公式 #109 — 未编号公式：\hat{\mathbf{A}} = \mathbf...
- 论文位置：[Appendix，`explain.md` 第 973 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L973)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\hat{\mathbf{A}} = \mathbf{A} + \mathbf{U}\mathbf{V}^{\mathrm{T}},
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
### 未编号公式 #110 — 未编号公式：\hat{\mathbf{A}}^{-1} = \m...
- 论文位置：[Appendix，`explain.md` 第 979 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L979)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\hat{\mathbf{A}}^{-1} = \mathbf{A}^{-1} - \mathbf{A}^{-1}\mathbf{U}\left(\mathbf{I} + \mathbf{V}^{\mathrm{T}}\mathbf{A}^{-1}\mathbf{U}\right)^{-1}\mathbf{V}^{\mathrm{T}}\mathbf{A}^{-1}.
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
