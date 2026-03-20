## 027. ### 2.2. Finite element discretization（explain.md:L225-L227）

- 对应关系说明：本公式位于离散层：单元残差与切线子块的组装。
- 最底层代码链接：
  - [main.cc:L3722-L3820](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L3722-L3820)
  - [main.cc:L3860-L3951](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L3860-L3951)
  - [main.cc:L3641-L3676](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L3641-L3676)

### 代码片段 1（main.cc:L3722-L3820）
```cpp
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

### 代码片段 2（main.cc:L3860-L3951）
```cpp
3860:   void PhaseFieldMonolithicSolve<dim>::assemble_system_B0_one_cell(
3861:       const typename DoFHandler<dim>::active_cell_iterator &cell,
3862:       ScratchData_ASM & scratch,
3863:       PerTaskData_ASM & data) const
3864:   {
3865:     data.reset();
3866:     scratch.reset();
3867:     scratch.m_fe_values.reinit(cell);
3868:     cell->get_dof_indices(data.m_local_dof_indices);
3869: 
3870:     scratch.m_fe_values[m_d_fe].get_function_values(
3871:       scratch.m_solution_previous_step, scratch.m_phasefield_previous_step_cell);
3872: 
3873:     const std::vector<std::shared_ptr<const PointHistory<dim>>> lqph =
3874:       m_quadrature_point_history.get_data(cell);
3875:     Assert(lqph.size() == m_n_q_points, ExcInternalError());
3876: 
3877:     const double delta_time = m_time.get_delta_t();
3878: 
3879:     for (const unsigned int q_point : scratch.m_fe_values.quadrature_point_indices())
3880:       {
3881:         for (const unsigned int k : scratch.m_fe_values.dof_indices())
3882:           {
3883:             const unsigned int k_group = m_fe.system_to_base_index(k).first.first;
3884: 
3885:             if (k_group == m_u_dof)
3886:               {
3887:                 scratch.m_Nx_disp[q_point][k] =
3888:                   scratch.m_fe_values[m_u_fe].value(k, q_point);
3889:                 scratch.m_grad_Nx_disp[q_point][k] =
3890:                   scratch.m_fe_values[m_u_fe].gradient(k, q_point);
3891:                 scratch.m_symm_grad_Nx_disp[q_point][k] =
3892:                   symmetrize(scratch.m_grad_Nx_disp[q_point][k]);
3893:               }
3894:             else if (k_group == m_d_dof)
3895:               {
3896: 		scratch.m_Nx_phasefield[q_point][k] =
3897: 		  scratch.m_fe_values[m_d_fe].value(k, q_point);
3898: 		scratch.m_grad_Nx_phasefield[q_point][k] =
3899: 		  scratch.m_fe_values[m_d_fe].gradient(k, q_point);
3900:               }
3901:             else
3902:               Assert(k_group <= m_d_dof, ExcInternalError());
3903:           }
3904:       }
3905: 
3906:     for (const unsigned int q_point : scratch.m_fe_values.quadrature_point_indices())
3907:       {
3908: 	const double length_scale            = lqph[q_point]->get_length_scale();
3909: 	const double gc                      = lqph[q_point]->get_critical_energy_release_rate();
3910: 	const double eta                     = lqph[q_point]->get_viscosity();
3911: 	const double current_positive_strain_energy = lqph[q_point]->get_current_positive_strain_energy();
3912: 
3913: 	const double phasefield_value        = lqph[q_point]->get_phase_field_value();
3914: 
3915:         const std::vector<double>         &      N_phasefield = scratch.m_Nx_phasefield[q_point];
3916:         const std::vector<Tensor<1, dim>> & grad_N_phasefield = scratch.m_grad_Nx_phasefield[q_point];
3917: 
3918:         //const SymmetricTensor<2, dim> & cauchy_stress_positive = lqph[q_point]->get_cauchy_stress_positive();
3919:         const SymmetricTensor<4, dim> & mechanical_C  = lqph[q_point]->get_mechanical_C();
3920: 
3921:         const std::vector<SymmetricTensor<2, dim>> & symm_grad_N_disp =
3922:           scratch.m_symm_grad_Nx_disp[q_point];
3923:         const double JxW = scratch.m_fe_values.JxW(q_point);
3924: 
3925:         SymmetricTensor<2, dim> symm_grad_Nx_i_x_C;
3926: 
3927:         for (const unsigned int i : scratch.m_fe_values.dof_indices())
3928:           {
3929:             const unsigned int i_group = m_fe.system_to_base_index(i).first.first;
3930: 
3931:             if (i_group == m_u_dof)
3932:               {
3933:                 symm_grad_Nx_i_x_C = symm_grad_N_disp[i] * mechanical_C;
3934:               }
3935: 
3936:             for (const unsigned int j : scratch.m_fe_values.dof_indices())
3937:               {
3938:                 const unsigned int j_group = m_fe.system_to_base_index(j).first.first;
3939: 
3940:                 if ((i_group == j_group) && (i_group == m_u_dof))
3941:                   {
3942:                     data.m_cell_matrix(i, j) += symm_grad_Nx_i_x_C * symm_grad_N_disp[j] * JxW;
3943:                   }
3944:                 else if ((i_group == j_group) && (i_group == m_d_dof))
3945:                   {
3946:                     data.m_cell_matrix(i, j) += (  (   gc/length_scale + eta/delta_time
3947:                 	                             + degradation_function_2nd_order_derivative(phasefield_value)
3948: 						     * current_positive_strain_energy  )
3949:                 	                          * N_phasefield[i] * N_phasefield[j]
3950: 					          + gc * length_scale * grad_N_phasefield[i] * grad_N_phasefield[j]
3951: 					        ) * JxW;
```

### 代码片段 3（main.cc:L3641-L3676）
```cpp
3641:   void PhaseFieldMonolithicSolve<dim>::assemble_system_B0(const BlockVector<double> & solution_old)
3642:   {
3643:     m_timer.enter_subsection("Assemble B0");
3644: 
3645:     m_tangent_matrix = 0.0;
3646: 
3647:     const UpdateFlags uf_cell(update_values | update_gradients |
3648: 			      update_quadrature_points | update_JxW_values);
3649:     const UpdateFlags uf_face(update_values | update_normal_vectors |
3650:                               update_JxW_values);
3651: 
3652:     PerTaskData_ASM per_task_data(m_fe.n_dofs_per_cell());
3653:     ScratchData_ASM scratch_data(m_fe, m_qf_cell, uf_cell, m_qf_face, uf_face, solution_old);
3654: 
3655:     auto worker =
3656:       [this](const typename DoFHandler<dim>::active_cell_iterator &cell,
3657: 	     ScratchData_ASM & scratch,
3658: 	     PerTaskData_ASM & data)
3659:       {
3660:         this->assemble_system_B0_one_cell(cell, scratch, data);
3661:       };
3662: 
3663:     auto copier = [this](const PerTaskData_ASM &data)
3664:       {
3665:         this->m_constraints.distribute_local_to_global(data.m_cell_matrix,
3666:                                                        data.m_local_dof_indices,
3667: 						       m_tangent_matrix);
3668:       };
3669: 
3670:     WorkStream::run(
3671:       m_dof_handler.active_cell_iterators(),
3672:       worker,
3673:       copier,
3674:       scratch_data,
3675:       per_task_data);
3676: 
```

### 公式（与 explain.md 一致）
\[
\hat{\mathbf{K}} = \left[ \begin{array}{cc}\mathbf{K}_{uu} & \mathbf{0}\\ \mathbf{0} & \mathbf{K}_{dd} \end{array} \right], \quad (11)
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 028. ### 3.1. Algorithm overview（explain.md:L239-L241）

- 对应关系说明：本公式位于算法主循环层：投影、搜索方向、线搜索与收敛判据。
- 最底层代码链接：
  - [main.cc:L1379-L1442](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L1379-L1442)
  - [main.cc:L4953-L5060](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L4953-L5060)
  - [main.cc:L5658-L5691](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L5658-L5691)

### 代码片段 1（main.cc:L1379-L1442）
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
1395:   }
1396: 
1397:   template <int dim>
1398:   std::priority_queue< std::pair<double, unsigned int>,
1399:                        std::vector<std::pair<double, unsigned int>>,
1400:   		       std::greater<std::pair<double, unsigned int>> >
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
```

### 代码片段 2（main.cc:L4953-L5060）
```cpp
4953:   solve_nonlinear_timestep_LBFGS_B(BlockVector<double> & solution_delta,
4954: 				   BlockVector<double> & LBFGS_update_refine)
4955:   {
4956:     BlockVector<double> LBFGS_update(m_dofs_per_block);
4957:     BlockVector<double> solution_delta_cauchy_point(m_dofs_per_block);
4958:     LBFGS_update = 0.0;
4959: 
4960:     const unsigned int LBFGS_m = m_parameters.m_LBFGS_m;
4961: 
4962:     unsigned int LBFGS_iteration = 0;
4963: 
4964:     m_error_residual.reset();
4965:     m_error_residual_0.reset();
4966:     m_error_residual_norm.reset();
4967:     m_error_update.reset();
4968:     m_error_update_0.reset();
4969:     m_error_update_norm.reset();
4970: 
4971:     if (m_parameters.m_output_iteration_history)
4972:       print_conv_header_LBFGSB();
4973: 
4974:     BlockVector<double> LBFGS_s_vector(m_dofs_per_block);
4975:     BlockVector<double> LBFGS_y_vector(m_dofs_per_block);
4976:     BlockVector<double> free_dofs(m_dofs_per_block);
4977:     BlockVector<double> b0xs_vector(m_dofs_per_block);
4978: 
4979:     // all the list goes from k-m to k-1
4980:     // the front is the oldest quantity,and the end is
4981:     // newest quantity
4982:     std::list<BlockVector<double>> s_vector_list;
4983:     std::list<BlockVector<double>> y_vector_list;
4984:     std::list<double> s_dot_y_list;
4985:     std::list<BlockVector<double>> b0xs_vector_list;
4986: 
4987:     double line_search_parameter = 0.0;
4988: 
4989:     unsigned int lower_bound_number_old = 0;
4990:     unsigned int upper_bound_number_old = 0;
4991:     unsigned int lowerupper_bound_number_old = 0;
4992: 
4993:     unsigned int lower_bound_number_new = 0;
4994:     unsigned int upper_bound_number_new = 0;
4995:     unsigned int lowerupper_bound_number_new = 0;
4996: 
4997:     for (; LBFGS_iteration < m_parameters.m_max_iterations_BFGS; ++LBFGS_iteration)
4998:       {
4999: 	if (m_parameters.m_output_iteration_history)
5000: 	  m_logfile << '\t' << std::setw(4) << LBFGS_iteration << ' '
5001:                     << std::flush;
5002: 
5003:         make_constraints(LBFGS_iteration);
5004: 
5005:         // At the first step, we simply distribute the inhomogeneous part of
5006:         // the constraints
5007:         if (LBFGS_iteration == 0)
5008:           {
5009:             // use the solution from the previous solve on the
5010:             // refined mesh as initial guess
5011:             LBFGS_update = LBFGS_update_refine;
5012: 
5013:             m_constraints.distribute(LBFGS_update);
5014:             solution_delta += LBFGS_update;
5015:             update_qph_incremental(solution_delta, m_solution);
5016:             assemble_system_rhs_BFGS_parallel(m_solution, m_system_rhs);
5017:             m_logfile << "  | " << std::endl;
5018: 
5019:             continue;
5020:           }
5021: 
5022:         get_error_residual_LBFGSB(m_error_residual,
5023: 				  solution_delta);
5024: 
5025:         if (LBFGS_iteration == 1)
5026:           m_error_residual_0 = m_error_residual;
5027: 
5028:         m_error_residual_norm = m_error_residual;
5029: 
5030:         if (m_parameters.m_relative_residual)
5031:           m_error_residual_norm.normalize(m_error_residual_0);
5032: 
5033:         if (LBFGS_iteration > 1 && m_error_update_norm.m_u <= m_parameters.m_tol_u_incr
5034:                                 && m_error_residual_norm.m_u <= m_parameters.m_tol_u_residual
5035: 			        && m_error_update_norm.m_d <= m_parameters.m_tol_d_incr
5036: 			        && m_error_residual_norm.m_d <= m_parameters.m_tol_d_residual
5037: 				&& lower_bound_number_new == lower_bound_number_old
5038: 				&& upper_bound_number_new == upper_bound_number_old
5039: 				&& lowerupper_bound_number_new == lowerupper_bound_number_old)
5040:           {
5041:             if (m_parameters.m_output_iteration_history)
5042:               {
5043: 		m_logfile << "  | ";
5044: 		m_logfile << " CONVERGED! " << std::fixed << std::setprecision(3) << std::setw(7)
5045: 			  << std::scientific
5046: 		      << "           ---      "
5047: 		      << "\t\t\t"
5048: 		      << "  " << m_error_residual_norm.m_norm
5049: 		      << "  " << m_error_residual_norm.m_u
5050: 		      << "  " << m_error_residual_norm.m_d
5051: 		      << "  " << m_error_update_norm.m_norm
5052: 		      << "  " << m_error_update_norm.m_u
5053: 		      << "  " << m_error_update_norm.m_d
5054: 		      << "  " << std::endl;
5055: 
5056: 		m_logfile << '\t';
5057: 		for (unsigned int i = 0; i < 130; ++i)
5058: 		  m_logfile << '_';
5059: 		m_logfile << std::endl;
5060:               }
```

### 代码片段 3（main.cc:L5658-L5691）
```cpp
5658: 		+ LBFGS_update.block(m_d_dof)[i] > 1.0)
5659: 	      LBFGS_update.block(m_d_dof)[i] = 1.0 - m_solution.block(m_d_dof)[i]
5660: 						   - solution_delta.block(m_d_dof)[i];
5661: 	  }
5662: 
5663: 	m_constraints.distribute(LBFGS_update);
5664: 
5665: 	// We need a line search algorithm to decide line_search_parameter
5666: 
5667:         if(m_parameters.m_type_line_search == "StrongWolfe")
5668:           {
5669: 	    const double phi_0 = calculate_energy_functional();
5670: 	    const double phi_0_prime = m_system_rhs * LBFGS_update;
5671: 
5672: 	    line_search_parameter = line_search_stepsize_strong_wolfe(phi_0,
5673: 								      phi_0_prime,
5674: 								      LBFGS_update,
5675: 								      solution_delta);
5676:           }
5677:         else if(m_parameters.m_type_line_search == "GradientBased")
5678:           {
5679: 	    // LBFGS_r_vector is the search direction
5680: 	    line_search_parameter = line_search_stepsize_gradient_based(LBFGS_update,
5681: 									solution_delta);
5682:           }
5683:         else
5684:           {
5685:             Assert(false, ExcMessage("An unknown line search method is called!"));
5686:           }
5687: 
5688: 	LBFGS_update *= line_search_parameter;
5689: 
5690:         get_error_update(LBFGS_update, m_error_update);
5691:         if (LBFGS_iteration == 1)
```

### 公式（与 explain.md 一致）
\[
\left(\pmb{u}_A,d_A\right) = \arg \min \Pi (\pmb {u}_A,d_A) \quad (12)
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 029. ### 3.1. Algorithm overview（explain.md:L245-L247）

- 对应关系说明：本公式位于算法主循环层：投影、搜索方向、线搜索与收敛判据。
- 最底层代码链接：
  - [main.cc:L1379-L1442](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L1379-L1442)
  - [main.cc:L4953-L5060](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L4953-L5060)
  - [main.cc:L5658-L5691](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L5658-L5691)

### 代码片段 1（main.cc:L1379-L1442）
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
1395:   }
1396: 
1397:   template <int dim>
1398:   std::priority_queue< std::pair<double, unsigned int>,
1399:                        std::vector<std::pair<double, unsigned int>>,
1400:   		       std::greater<std::pair<double, unsigned int>> >
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
```

### 代码片段 2（main.cc:L4953-L5060）
```cpp
4953:   solve_nonlinear_timestep_LBFGS_B(BlockVector<double> & solution_delta,
4954: 				   BlockVector<double> & LBFGS_update_refine)
4955:   {
4956:     BlockVector<double> LBFGS_update(m_dofs_per_block);
4957:     BlockVector<double> solution_delta_cauchy_point(m_dofs_per_block);
4958:     LBFGS_update = 0.0;
4959: 
4960:     const unsigned int LBFGS_m = m_parameters.m_LBFGS_m;
4961: 
4962:     unsigned int LBFGS_iteration = 0;
4963: 
4964:     m_error_residual.reset();
4965:     m_error_residual_0.reset();
4966:     m_error_residual_norm.reset();
4967:     m_error_update.reset();
4968:     m_error_update_0.reset();
4969:     m_error_update_norm.reset();
4970: 
4971:     if (m_parameters.m_output_iteration_history)
4972:       print_conv_header_LBFGSB();
4973: 
4974:     BlockVector<double> LBFGS_s_vector(m_dofs_per_block);
4975:     BlockVector<double> LBFGS_y_vector(m_dofs_per_block);
4976:     BlockVector<double> free_dofs(m_dofs_per_block);
4977:     BlockVector<double> b0xs_vector(m_dofs_per_block);
4978: 
4979:     // all the list goes from k-m to k-1
4980:     // the front is the oldest quantity,and the end is
4981:     // newest quantity
4982:     std::list<BlockVector<double>> s_vector_list;
4983:     std::list<BlockVector<double>> y_vector_list;
4984:     std::list<double> s_dot_y_list;
4985:     std::list<BlockVector<double>> b0xs_vector_list;
4986: 
4987:     double line_search_parameter = 0.0;
4988: 
4989:     unsigned int lower_bound_number_old = 0;
4990:     unsigned int upper_bound_number_old = 0;
4991:     unsigned int lowerupper_bound_number_old = 0;
4992: 
4993:     unsigned int lower_bound_number_new = 0;
4994:     unsigned int upper_bound_number_new = 0;
4995:     unsigned int lowerupper_bound_number_new = 0;
4996: 
4997:     for (; LBFGS_iteration < m_parameters.m_max_iterations_BFGS; ++LBFGS_iteration)
4998:       {
4999: 	if (m_parameters.m_output_iteration_history)
5000: 	  m_logfile << '\t' << std::setw(4) << LBFGS_iteration << ' '
5001:                     << std::flush;
5002: 
5003:         make_constraints(LBFGS_iteration);
5004: 
5005:         // At the first step, we simply distribute the inhomogeneous part of
5006:         // the constraints
5007:         if (LBFGS_iteration == 0)
5008:           {
5009:             // use the solution from the previous solve on the
5010:             // refined mesh as initial guess
5011:             LBFGS_update = LBFGS_update_refine;
5012: 
5013:             m_constraints.distribute(LBFGS_update);
5014:             solution_delta += LBFGS_update;
5015:             update_qph_incremental(solution_delta, m_solution);
5016:             assemble_system_rhs_BFGS_parallel(m_solution, m_system_rhs);
5017:             m_logfile << "  | " << std::endl;
5018: 
5019:             continue;
5020:           }
5021: 
5022:         get_error_residual_LBFGSB(m_error_residual,
5023: 				  solution_delta);
5024: 
5025:         if (LBFGS_iteration == 1)
5026:           m_error_residual_0 = m_error_residual;
5027: 
5028:         m_error_residual_norm = m_error_residual;
5029: 
5030:         if (m_parameters.m_relative_residual)
5031:           m_error_residual_norm.normalize(m_error_residual_0);
5032: 
5033:         if (LBFGS_iteration > 1 && m_error_update_norm.m_u <= m_parameters.m_tol_u_incr
5034:                                 && m_error_residual_norm.m_u <= m_parameters.m_tol_u_residual
5035: 			        && m_error_update_norm.m_d <= m_parameters.m_tol_d_incr
5036: 			        && m_error_residual_norm.m_d <= m_parameters.m_tol_d_residual
5037: 				&& lower_bound_number_new == lower_bound_number_old
5038: 				&& upper_bound_number_new == upper_bound_number_old
5039: 				&& lowerupper_bound_number_new == lowerupper_bound_number_old)
5040:           {
5041:             if (m_parameters.m_output_iteration_history)
5042:               {
5043: 		m_logfile << "  | ";
5044: 		m_logfile << " CONVERGED! " << std::fixed << std::setprecision(3) << std::setw(7)
5045: 			  << std::scientific
5046: 		      << "           ---      "
5047: 		      << "\t\t\t"
5048: 		      << "  " << m_error_residual_norm.m_norm
5049: 		      << "  " << m_error_residual_norm.m_u
5050: 		      << "  " << m_error_residual_norm.m_d
5051: 		      << "  " << m_error_update_norm.m_norm
5052: 		      << "  " << m_error_update_norm.m_u
5053: 		      << "  " << m_error_update_norm.m_d
5054: 		      << "  " << std::endl;
5055: 
5056: 		m_logfile << '\t';
5057: 		for (unsigned int i = 0; i < 130; ++i)
5058: 		  m_logfile << '_';
5059: 		m_logfile << std::endl;
5060:               }
```

### 代码片段 3（main.cc:L5658-L5691）
```cpp
5658: 		+ LBFGS_update.block(m_d_dof)[i] > 1.0)
5659: 	      LBFGS_update.block(m_d_dof)[i] = 1.0 - m_solution.block(m_d_dof)[i]
5660: 						   - solution_delta.block(m_d_dof)[i];
5661: 	  }
5662: 
5663: 	m_constraints.distribute(LBFGS_update);
5664: 
5665: 	// We need a line search algorithm to decide line_search_parameter
5666: 
5667:         if(m_parameters.m_type_line_search == "StrongWolfe")
5668:           {
5669: 	    const double phi_0 = calculate_energy_functional();
5670: 	    const double phi_0_prime = m_system_rhs * LBFGS_update;
5671: 
5672: 	    line_search_parameter = line_search_stepsize_strong_wolfe(phi_0,
5673: 								      phi_0_prime,
5674: 								      LBFGS_update,
5675: 								      solution_delta);
5676:           }
5677:         else if(m_parameters.m_type_line_search == "GradientBased")
5678:           {
5679: 	    // LBFGS_r_vector is the search direction
5680: 	    line_search_parameter = line_search_stepsize_gradient_based(LBFGS_update,
5681: 									solution_delta);
5682:           }
5683:         else
5684:           {
5685:             Assert(false, ExcMessage("An unknown line search method is called!"));
5686:           }
5687: 
5688: 	LBFGS_update *= line_search_parameter;
5689: 
5690:         get_error_update(LBFGS_update, m_error_update);
5691:         if (LBFGS_iteration == 1)
```

### 公式（与 explain.md 一致）
\[
d_A^{(n)}\leq d_A\leq 1 \quad (13)
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 030. ### 3.1. Algorithm overview（explain.md:L251-L253）

- 对应关系说明：本公式位于算法主循环层：投影、搜索方向、线搜索与收敛判据。
- 最底层代码链接：
  - [main.cc:L1379-L1442](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L1379-L1442)
  - [main.cc:L4953-L5060](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L4953-L5060)
  - [main.cc:L5658-L5691](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L5658-L5691)

### 代码片段 1（main.cc:L1379-L1442）
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
1395:   }
1396: 
1397:   template <int dim>
1398:   std::priority_queue< std::pair<double, unsigned int>,
1399:                        std::vector<std::pair<double, unsigned int>>,
1400:   		       std::greater<std::pair<double, unsigned int>> >
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
```

### 代码片段 2（main.cc:L4953-L5060）
```cpp
4953:   solve_nonlinear_timestep_LBFGS_B(BlockVector<double> & solution_delta,
4954: 				   BlockVector<double> & LBFGS_update_refine)
4955:   {
4956:     BlockVector<double> LBFGS_update(m_dofs_per_block);
4957:     BlockVector<double> solution_delta_cauchy_point(m_dofs_per_block);
4958:     LBFGS_update = 0.0;
4959: 
4960:     const unsigned int LBFGS_m = m_parameters.m_LBFGS_m;
4961: 
4962:     unsigned int LBFGS_iteration = 0;
4963: 
4964:     m_error_residual.reset();
4965:     m_error_residual_0.reset();
4966:     m_error_residual_norm.reset();
4967:     m_error_update.reset();
4968:     m_error_update_0.reset();
4969:     m_error_update_norm.reset();
4970: 
4971:     if (m_parameters.m_output_iteration_history)
4972:       print_conv_header_LBFGSB();
4973: 
4974:     BlockVector<double> LBFGS_s_vector(m_dofs_per_block);
4975:     BlockVector<double> LBFGS_y_vector(m_dofs_per_block);
4976:     BlockVector<double> free_dofs(m_dofs_per_block);
4977:     BlockVector<double> b0xs_vector(m_dofs_per_block);
4978: 
4979:     // all the list goes from k-m to k-1
4980:     // the front is the oldest quantity,and the end is
4981:     // newest quantity
4982:     std::list<BlockVector<double>> s_vector_list;
4983:     std::list<BlockVector<double>> y_vector_list;
4984:     std::list<double> s_dot_y_list;
4985:     std::list<BlockVector<double>> b0xs_vector_list;
4986: 
4987:     double line_search_parameter = 0.0;
4988: 
4989:     unsigned int lower_bound_number_old = 0;
4990:     unsigned int upper_bound_number_old = 0;
4991:     unsigned int lowerupper_bound_number_old = 0;
4992: 
4993:     unsigned int lower_bound_number_new = 0;
4994:     unsigned int upper_bound_number_new = 0;
4995:     unsigned int lowerupper_bound_number_new = 0;
4996: 
4997:     for (; LBFGS_iteration < m_parameters.m_max_iterations_BFGS; ++LBFGS_iteration)
4998:       {
4999: 	if (m_parameters.m_output_iteration_history)
5000: 	  m_logfile << '\t' << std::setw(4) << LBFGS_iteration << ' '
5001:                     << std::flush;
5002: 
5003:         make_constraints(LBFGS_iteration);
5004: 
5005:         // At the first step, we simply distribute the inhomogeneous part of
5006:         // the constraints
5007:         if (LBFGS_iteration == 0)
5008:           {
5009:             // use the solution from the previous solve on the
5010:             // refined mesh as initial guess
5011:             LBFGS_update = LBFGS_update_refine;
5012: 
5013:             m_constraints.distribute(LBFGS_update);
5014:             solution_delta += LBFGS_update;
5015:             update_qph_incremental(solution_delta, m_solution);
5016:             assemble_system_rhs_BFGS_parallel(m_solution, m_system_rhs);
5017:             m_logfile << "  | " << std::endl;
5018: 
5019:             continue;
5020:           }
5021: 
5022:         get_error_residual_LBFGSB(m_error_residual,
5023: 				  solution_delta);
5024: 
5025:         if (LBFGS_iteration == 1)
5026:           m_error_residual_0 = m_error_residual;
5027: 
5028:         m_error_residual_norm = m_error_residual;
5029: 
5030:         if (m_parameters.m_relative_residual)
5031:           m_error_residual_norm.normalize(m_error_residual_0);
5032: 
5033:         if (LBFGS_iteration > 1 && m_error_update_norm.m_u <= m_parameters.m_tol_u_incr
5034:                                 && m_error_residual_norm.m_u <= m_parameters.m_tol_u_residual
5035: 			        && m_error_update_norm.m_d <= m_parameters.m_tol_d_incr
5036: 			        && m_error_residual_norm.m_d <= m_parameters.m_tol_d_residual
5037: 				&& lower_bound_number_new == lower_bound_number_old
5038: 				&& upper_bound_number_new == upper_bound_number_old
5039: 				&& lowerupper_bound_number_new == lowerupper_bound_number_old)
5040:           {
5041:             if (m_parameters.m_output_iteration_history)
5042:               {
5043: 		m_logfile << "  | ";
5044: 		m_logfile << " CONVERGED! " << std::fixed << std::setprecision(3) << std::setw(7)
5045: 			  << std::scientific
5046: 		      << "           ---      "
5047: 		      << "\t\t\t"
5048: 		      << "  " << m_error_residual_norm.m_norm
5049: 		      << "  " << m_error_residual_norm.m_u
5050: 		      << "  " << m_error_residual_norm.m_d
5051: 		      << "  " << m_error_update_norm.m_norm
5052: 		      << "  " << m_error_update_norm.m_u
5053: 		      << "  " << m_error_update_norm.m_d
5054: 		      << "  " << std::endl;
5055: 
5056: 		m_logfile << '\t';
5057: 		for (unsigned int i = 0; i < 130; ++i)
5058: 		  m_logfile << '_';
5059: 		m_logfile << std::endl;
5060:               }
```

### 代码片段 3（main.cc:L5658-L5691）
```cpp
5658: 		+ LBFGS_update.block(m_d_dof)[i] > 1.0)
5659: 	      LBFGS_update.block(m_d_dof)[i] = 1.0 - m_solution.block(m_d_dof)[i]
5660: 						   - solution_delta.block(m_d_dof)[i];
5661: 	  }
5662: 
5663: 	m_constraints.distribute(LBFGS_update);
5664: 
5665: 	// We need a line search algorithm to decide line_search_parameter
5666: 
5667:         if(m_parameters.m_type_line_search == "StrongWolfe")
5668:           {
5669: 	    const double phi_0 = calculate_energy_functional();
5670: 	    const double phi_0_prime = m_system_rhs * LBFGS_update;
5671: 
5672: 	    line_search_parameter = line_search_stepsize_strong_wolfe(phi_0,
5673: 								      phi_0_prime,
5674: 								      LBFGS_update,
5675: 								      solution_delta);
5676:           }
5677:         else if(m_parameters.m_type_line_search == "GradientBased")
5678:           {
5679: 	    // LBFGS_r_vector is the search direction
5680: 	    line_search_parameter = line_search_stepsize_gradient_based(LBFGS_update,
5681: 									solution_delta);
5682:           }
5683:         else
5684:           {
5685:             Assert(false, ExcMessage("An unknown line search method is called!"));
5686:           }
5687: 
5688: 	LBFGS_update *= line_search_parameter;
5689: 
5690:         get_error_update(LBFGS_update, m_error_update);
5691:         if (LBFGS_iteration == 1)
```

### 公式（与 explain.md 一致）
\[
\left[\pmb {u}_A,d_A\right]^{\mathrm{T}} = \mathbf{C}\left[\pmb {u}_A,d_A\right]^{\mathrm{T}} + \pmb {k}. \quad (14)
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。
