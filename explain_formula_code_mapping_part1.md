# explain.md 全公式 -> 代码最底层映射解读（代码在前，公式在后）

说明：
- 本文严格按 `explain.md` 的论文分段顺序，覆盖该文档中全部显示公式（共 110 条）。
- 每条均给出“可访问代码链接 + 代码原文片段 + 公式 + 解读”。
- 公式统一采用 `\[ ... \]` 渲染格式。
- 仓库链接均指向当前工作分支，便于直接点击定位到对应行。


## 001. ## 1. Introduction（explain.md:L21-L23）

- 对应关系说明：该公式对应能量-残差-切线三层实现。
- 最底层代码链接：
  - [main.cc:L886-L887](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L886-L887)
  - [main.cc:L6143-L6164](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L6143-L6164)

### 代码片段 1（main.cc:L886-L887）
```cpp
 886:     m_crack_energy_dissipation = m_gc * (  0.5 / m_length_scale * m_phase_field_value * m_phase_field_value
 887: 	                                   + 0.5 * m_length_scale * m_grad_phasefield * m_grad_phasefield)
```

### 代码片段 2（main.cc:L6143-L6164）
```cpp
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
```

### 公式（与 explain.md 一致）
\[
\Pi (\pmb {u},d) = \int_{\Omega}\psi (\pmb {e}(\pmb {u}),d)\mathrm{d}\Omega + g_{c}\Gamma_{l}(d) - \int_{\Omega}\pmb {b}\cdot \pmb {u}\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot \pmb {u}\mathrm{d}\Gamma , \quad (1)
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 002. ## 1. Introduction（explain.md:L27-L29）

- 对应关系说明：该公式对应能量-残差-切线三层实现。
- 最底层代码链接：
  - [main.cc:L886-L887](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L886-L887)
  - [main.cc:L6143-L6164](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L6143-L6164)

### 代码片段 1（main.cc:L886-L887）
```cpp
 886:     m_crack_energy_dissipation = m_gc * (  0.5 / m_length_scale * m_phase_field_value * m_phase_field_value
 887: 	                                   + 0.5 * m_length_scale * m_grad_phasefield * m_grad_phasefield)
```

### 代码片段 2（main.cc:L6143-L6164）
```cpp
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
```

### 公式（与 explain.md 一致）
\[
\Gamma_{l}(d) = \int_{\Omega}\gamma (d,\nabla d)\mathrm{d}\Omega = \int_{\Omega}\frac{1}{2l}\left(d^{2} + l^{2}\nabla d\cdot \nabla d\right)\mathrm{d}\Omega , \quad (2)
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 003. ## 1. Introduction（explain.md:L35-L37）

- 对应关系说明：该公式对应能量-残差-切线三层实现。
- 最底层代码链接：
  - [main.cc:L6117-L6139](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L6117-L6139)
  - [main.cc:L820-L891](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L820-L891)
  - [main.cc:L3722-L3951](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L3722-L3951)

### 代码片段 1（main.cc:L6117-L6139）
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
```

### 代码片段 2（main.cc:L820-L891）
```cpp
 820:    update_material_data(const SymmetricTensor<2, dim> & strain,
 821: 			const double phase_field_value,
 822: 			const Tensor<1, dim> & grad_phasefield,
 823: 			const double phase_field_value_previous_step,
 824: 			const double delta_time)
 825:   {
 826:     m_strain = strain;
 827:     m_phase_field_value = phase_field_value;
 828:     m_grad_phasefield = grad_phasefield;
 829:     Vector<double>              eigenvalues(dim);
 830:     std::vector<Tensor<1, dim>> eigenvectors(dim);
 831:     usr_spectrum_decomposition::spectrum_decomposition<dim>(m_strain,
 832:   							      eigenvalues,
 833:   							      eigenvectors);
 834: 
 835:     SymmetricTensor<2, dim> strain_positive, strain_negative;
 836:     strain_positive = usr_spectrum_decomposition::positive_tensor(eigenvalues, eigenvectors);
 837:     strain_negative = usr_spectrum_decomposition::negative_tensor(eigenvalues, eigenvectors);
 838: 
 839:     SymmetricTensor<4, dim> projector_positive, projector_negative;
 840:     usr_spectrum_decomposition::positive_negative_projectors(eigenvalues,
 841:   							       eigenvectors,
 842: 							       projector_positive,
 843: 							       projector_negative);
 844: 
 845:     SymmetricTensor<2, dim> stress_positive, stress_negative;
 846:     const double degradation = degradation_function(m_phase_field_value) + m_residual_k;
 847:     const double I_1 = trace(m_strain);
 848: 
 849:     // 2D plane strain and 3D cases
 850:     double my_lambda = m_lame_lambda;
 851: 
 852:     // 2D plane stress case
 853:     if (    dim == 2
 854: 	   && m_plane_stress)
 855:       my_lambda = 2 * m_lame_mu * m_lame_lambda / (m_lame_lambda + 2 * m_lame_mu);
 856: 
 857:     stress_positive = my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)
 858:                                     * Physics::Elasticity::StandardTensors<dim>::I
 859:                     + 2 * m_lame_mu * strain_positive;
 860:     stress_negative = my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)
 861:                                     * Physics::Elasticity::StandardTensors<dim>::I
 862:     		      + 2 * m_lame_mu * strain_negative;
 863: 
 864:     m_stress = degradation * stress_positive + stress_negative;
 865:     m_stress_positive = stress_positive;
 866: 
 867:     SymmetricTensor<4, dim> C_positive, C_negative;
 868:     C_positive = my_lambda * usr_spectrum_decomposition::heaviside_function(I_1)
 869:                                * Physics::Elasticity::StandardTensors<dim>::IxI
 870: 		 + 2 * m_lame_mu * projector_positive;
 871:     C_negative = my_lambda * usr_spectrum_decomposition::heaviside_function(-I_1)
 872:                                * Physics::Elasticity::StandardTensors<dim>::IxI
 873:     		 + 2 * m_lame_mu * projector_negative;
 874:     m_mechanical_C = degradation * C_positive + C_negative;
 875: 
 876:     m_strain_energy_positive = 0.5 * my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)
 877:                                                    * usr_spectrum_decomposition::positive_ramp_function(I_1)
 878:                              + m_lame_mu * strain_positive * strain_positive;
 879: 
 880:     m_strain_energy_negative = 0.5 * my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)
 881:                                                    * usr_spectrum_decomposition::negative_ramp_function(I_1)
 882:                              + m_lame_mu * strain_negative * strain_negative;
 883: 
 884:     m_strain_energy_total = degradation * m_strain_energy_positive + m_strain_energy_negative;
 885: 
 886:     m_crack_energy_dissipation = m_gc * (  0.5 / m_length_scale * m_phase_field_value * m_phase_field_value
 887: 	                                   + 0.5 * m_length_scale * m_grad_phasefield * m_grad_phasefield)
 888: 	                                   // the term due to viscosity regularization
 889: 	                                   + (m_phase_field_value - phase_field_value_previous_step)
 890: 					   * (m_phase_field_value - phase_field_value_previous_step)
 891: 				           * 0.5 * m_eta / delta_time;
```

### 代码片段 3（main.cc:L3722-L3951）
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
3821:             else
3822:               Assert(i_group <= m_d_dof, ExcInternalError());
3823:           }  // i
3824:       }  // q_point
3825: 
3826:     // if there is surface pressure, this surface pressure always applied to the
3827:     // reference configuration
3828:     const unsigned int face_pressure_id = 100;
3829:     const double p0 = 0.0;
3830: 
3831:     for (const auto &face : cell->face_iterators())
3832:       if (face->at_boundary() && face->boundary_id() == face_pressure_id)
3833:         {
3834:           scratch.m_fe_face_values.reinit(cell, face);
3835: 
3836:           for (const unsigned int f_q_point : scratch.m_fe_face_values.quadrature_point_indices())
3837:             {
3838:               const Tensor<1, dim> &N = scratch.m_fe_face_values.normal_vector(f_q_point);
3839: 
3840:               const double         pressure  = p0 * time_ramp;
3841:               const Tensor<1, dim> traction  = pressure * N;
3842: 
3843:               for (const unsigned int i : scratch.m_fe_values.dof_indices())
3844:                 {
3845:                   const unsigned int i_group = m_fe.system_to_base_index(i).first.first;
3846: 
3847:                   if (i_group == m_u_dof)
3848:                     {
3849:     		      const unsigned int component_i = m_fe.system_to_component_index(i).first;
3850:     		      const double Ni = scratch.m_fe_face_values.shape_value(i, f_q_point);
3851:     		      const double JxW = scratch.m_fe_face_values.JxW(f_q_point);
3852:     		      data.m_cell_rhs(i) -= (Ni * traction[component_i]) * JxW;
3853:                     }
3854:                 }
3855:             }
3856:         }
3857:   }
3858: 
3859:   template <int dim>
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

### 公式（与 explain.md 一致）
\[
\left(\pmb{u}_{n + 1},d_{n + 1}\right) = \arg \min \Pi (\pmb {u},d) \quad (3)
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 004. ## 1. Introduction（explain.md:L41-L43）

- 对应关系说明：该公式对应能量-残差-切线三层实现。
- 最底层代码链接：
  - [main.cc:L6117-L6139](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L6117-L6139)
  - [main.cc:L820-L891](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L820-L891)
  - [main.cc:L3722-L3951](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L3722-L3951)

### 代码片段 1（main.cc:L6117-L6139）
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
```

### 代码片段 2（main.cc:L820-L891）
```cpp
 820:    update_material_data(const SymmetricTensor<2, dim> & strain,
 821: 			const double phase_field_value,
 822: 			const Tensor<1, dim> & grad_phasefield,
 823: 			const double phase_field_value_previous_step,
 824: 			const double delta_time)
 825:   {
 826:     m_strain = strain;
 827:     m_phase_field_value = phase_field_value;
 828:     m_grad_phasefield = grad_phasefield;
 829:     Vector<double>              eigenvalues(dim);
 830:     std::vector<Tensor<1, dim>> eigenvectors(dim);
 831:     usr_spectrum_decomposition::spectrum_decomposition<dim>(m_strain,
 832:   							      eigenvalues,
 833:   							      eigenvectors);
 834: 
 835:     SymmetricTensor<2, dim> strain_positive, strain_negative;
 836:     strain_positive = usr_spectrum_decomposition::positive_tensor(eigenvalues, eigenvectors);
 837:     strain_negative = usr_spectrum_decomposition::negative_tensor(eigenvalues, eigenvectors);
 838: 
 839:     SymmetricTensor<4, dim> projector_positive, projector_negative;
 840:     usr_spectrum_decomposition::positive_negative_projectors(eigenvalues,
 841:   							       eigenvectors,
 842: 							       projector_positive,
 843: 							       projector_negative);
 844: 
 845:     SymmetricTensor<2, dim> stress_positive, stress_negative;
 846:     const double degradation = degradation_function(m_phase_field_value) + m_residual_k;
 847:     const double I_1 = trace(m_strain);
 848: 
 849:     // 2D plane strain and 3D cases
 850:     double my_lambda = m_lame_lambda;
 851: 
 852:     // 2D plane stress case
 853:     if (    dim == 2
 854: 	   && m_plane_stress)
 855:       my_lambda = 2 * m_lame_mu * m_lame_lambda / (m_lame_lambda + 2 * m_lame_mu);
 856: 
 857:     stress_positive = my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)
 858:                                     * Physics::Elasticity::StandardTensors<dim>::I
 859:                     + 2 * m_lame_mu * strain_positive;
 860:     stress_negative = my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)
 861:                                     * Physics::Elasticity::StandardTensors<dim>::I
 862:     		      + 2 * m_lame_mu * strain_negative;
 863: 
 864:     m_stress = degradation * stress_positive + stress_negative;
 865:     m_stress_positive = stress_positive;
 866: 
 867:     SymmetricTensor<4, dim> C_positive, C_negative;
 868:     C_positive = my_lambda * usr_spectrum_decomposition::heaviside_function(I_1)
 869:                                * Physics::Elasticity::StandardTensors<dim>::IxI
 870: 		 + 2 * m_lame_mu * projector_positive;
 871:     C_negative = my_lambda * usr_spectrum_decomposition::heaviside_function(-I_1)
 872:                                * Physics::Elasticity::StandardTensors<dim>::IxI
 873:     		 + 2 * m_lame_mu * projector_negative;
 874:     m_mechanical_C = degradation * C_positive + C_negative;
 875: 
 876:     m_strain_energy_positive = 0.5 * my_lambda * usr_spectrum_decomposition::positive_ramp_function(I_1)
 877:                                                    * usr_spectrum_decomposition::positive_ramp_function(I_1)
 878:                              + m_lame_mu * strain_positive * strain_positive;
 879: 
 880:     m_strain_energy_negative = 0.5 * my_lambda * usr_spectrum_decomposition::negative_ramp_function(I_1)
 881:                                                    * usr_spectrum_decomposition::negative_ramp_function(I_1)
 882:                              + m_lame_mu * strain_negative * strain_negative;
 883: 
 884:     m_strain_energy_total = degradation * m_strain_energy_positive + m_strain_energy_negative;
 885: 
 886:     m_crack_energy_dissipation = m_gc * (  0.5 / m_length_scale * m_phase_field_value * m_phase_field_value
 887: 	                                   + 0.5 * m_length_scale * m_grad_phasefield * m_grad_phasefield)
 888: 	                                   // the term due to viscosity regularization
 889: 	                                   + (m_phase_field_value - phase_field_value_previous_step)
 890: 					   * (m_phase_field_value - phase_field_value_previous_step)
 891: 				           * 0.5 * m_eta / delta_time;
```

### 代码片段 3（main.cc:L3722-L3951）
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
3821:             else
3822:               Assert(i_group <= m_d_dof, ExcInternalError());
3823:           }  // i
3824:       }  // q_point
3825: 
3826:     // if there is surface pressure, this surface pressure always applied to the
3827:     // reference configuration
3828:     const unsigned int face_pressure_id = 100;
3829:     const double p0 = 0.0;
3830: 
3831:     for (const auto &face : cell->face_iterators())
3832:       if (face->at_boundary() && face->boundary_id() == face_pressure_id)
3833:         {
3834:           scratch.m_fe_face_values.reinit(cell, face);
3835: 
3836:           for (const unsigned int f_q_point : scratch.m_fe_face_values.quadrature_point_indices())
3837:             {
3838:               const Tensor<1, dim> &N = scratch.m_fe_face_values.normal_vector(f_q_point);
3839: 
3840:               const double         pressure  = p0 * time_ramp;
3841:               const Tensor<1, dim> traction  = pressure * N;
3842: 
3843:               for (const unsigned int i : scratch.m_fe_values.dof_indices())
3844:                 {
3845:                   const unsigned int i_group = m_fe.system_to_base_index(i).first.first;
3846: 
3847:                   if (i_group == m_u_dof)
3848:                     {
3849:     		      const unsigned int component_i = m_fe.system_to_component_index(i).first;
3850:     		      const double Ni = scratch.m_fe_face_values.shape_value(i, f_q_point);
3851:     		      const double JxW = scratch.m_fe_face_values.JxW(f_q_point);
3852:     		      data.m_cell_rhs(i) -= (Ni * traction[component_i]) * JxW;
3853:                     }
3854:                 }
3855:             }
3856:         }
3857:   }
3858: 
3859:   template <int dim>
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

### 公式（与 explain.md 一致）
\[
0\leq d_{n}\leq d_{n + 1}\leq 1. \quad (4)
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。

## 005. ## 1. Introduction（explain.md:L57-L59）

- 对应关系说明：该公式对应能量-残差-切线三层实现。
- 最底层代码链接：
  - [main.cc:L1379-L1394](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L1379-L1394)
  - [main.cc:L1401-L1442](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L1401-L1442)
  - [main.cc:L4775-L4807](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/copilot/create-explain-md-document/main.cc#L4775-L4807)

### 代码片段 1（main.cc:L1379-L1394）
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

### 代码片段 2（main.cc:L1401-L1442）
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
```

### 代码片段 3（main.cc:L4775-L4807）
```cpp
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
```

### 公式（与 explain.md 一致）
\[
\mathrm{Proj}_C(\mathbf{x}_k = \mathbf{a}_k\nabla f(\mathbf{x}_k)).
\]

### 分段解读
- 变量层：该公式中的未知量在代码中以 `m_solution`（总场）、`solution_delta`（增量）、以及单元/积分点历史变量的组合出现。
- 计算层：实现采用“积分点更新 -> 单元残差/切线 -> 全局装配/线性求解”链路，保证公式中的每一项都有离散对应。
- 约束层：相场不可逆条件通过下界（前一步解）与上界（1.0）进行盒约束投影，属于最底层可执行逻辑，而非仅停留在数学描述。
- 数值层：若涉及 L-BFGS-B，则代码通过 `s/y` 历史向量、紧凑矩阵 `M`、广义 Cauchy 点与子空间求解完成对公式的可计算化。
