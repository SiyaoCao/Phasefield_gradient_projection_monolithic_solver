## 3. Gradient projection based monolithic scheme > 3.1. Algorithm overview

### Eq. (14) — 线性约束(Dirichlet/悬挂节点)
- 论文位置：[3. Gradient projection based monolithic scheme > 3.1. Algorithm overview，`explain.md` 第 251 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L251)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\left[\pmb {u}_A,d_A\right]^{\mathrm{T}} = \mathbf{C}\left[\pmb {u}_A,d_A\right]^{\mathrm{T}} + \pmb {k}. \quad (14)
\]
- 对应最底层计算代码（公式在前，代码在后）：
  - Dirichlet/悬挂节点线性约束构造：[main.cc:3357-3638](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3357-L3638)
```cpp
3357:   void PhaseFieldMonolithicSolve<dim>::make_constraints(const unsigned int it_nr)
3358:   {
3359:     const bool apply_dirichlet_bc = (it_nr == 0);
3360: 
3361:     if (it_nr > 1)
3362:       {
3363:         return;
3364:       }
3365: 
3366:     if (apply_dirichlet_bc)
3367:       {
3368: 	m_constraints.clear();
3369: 	DoFTools::make_hanging_node_constraints(m_dof_handler,
3370: 						m_constraints);
3371: 
3372: 	const FEValuesExtractors::Scalar x_displacement(0);
3373: 	const FEValuesExtractors::Scalar y_displacement(1);
3374: 	const FEValuesExtractors::Scalar z_displacement(2);
3375: 
3376: 	const FEValuesExtractors::Vector displacements(0);
3377: 
3378: 	if (   m_parameters.m_scenario == 1
3379: 	    || m_parameters.m_scenario == 3)
3380: 	  {
3381: 	    // Dirichlet B,C. bottom surface
3382: 	    const int boundary_id_bottom_surface = 0;
3383: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3384: 						     boundary_id_bottom_surface,
3385: 						     Functions::ZeroFunction<dim>(m_n_components),
3386: 						     m_constraints,
3387: 						     m_fe.component_mask(y_displacement));
3388: 
3389: 	    typename Triangulation<dim>::active_vertex_iterator vertex_itr;
3390: 	    vertex_itr = m_triangulation.begin_active_vertex();
3391: 	    std::vector<types::global_dof_index> node_xy(m_fe.dofs_per_vertex);
3392: 
3393: 	    for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)
3394: 	      {
3395: 		if (   (std::fabs(vertex_itr->vertex()[0] - 0.0) < 1.0e-9)
3396: 		    && (std::fabs(vertex_itr->vertex()[1] - 0.0) < 1.0e-9) )
3397: 		  {
3398: 		    node_xy = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);
3399: 		  }
3400: 	      }
3401: 	    m_constraints.add_line(node_xy[0]);
3402: 	    m_constraints.set_inhomogeneity(node_xy[0], 0.0);
3403: 
3404: 	    m_constraints.add_line(node_xy[1]);
3405: 	    m_constraints.set_inhomogeneity(node_xy[1], 0.0);
3406: 
3407: 	    const int boundary_id_top_surface = 1;
3408: 	    /*
3409: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3410: 						     boundary_id_top_surface,
3411: 						     Functions::ZeroFunction<dim>(m_n_components),
3412: 						     m_constraints,
3413: 						     m_fe.component_mask(x_displacement));
3414: 	    */
3415:             const double time_inc = m_time.get_delta_t();
3416:             double disp_magnitude = m_time.get_magnitude();
3417: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3418: 						     boundary_id_top_surface,
3419: 						     Functions::ConstantFunction<dim>(
3420: 						       disp_magnitude*time_inc, m_n_components),
3421: 						     m_constraints,
3422: 						     m_fe.component_mask(y_displacement));
3423: 	  }
3424: 	else if (   m_parameters.m_scenario == 2
3425: 	         || m_parameters.m_scenario == 4)
3426: 	  {
3427: 	    // Dirichlet B,C. bottom surface
3428: 	    const int boundary_id_bottom_surface = 0;
3429: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3430: 						     boundary_id_bottom_surface,
3431: 						     Functions::ZeroFunction<dim>(m_n_components),
3432: 						     m_constraints,
3433: 						     m_fe.component_mask(displacements));
3434: 
3435: 	    const int boundary_id_top_surface = 1;
3436: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3437: 						     boundary_id_top_surface,
3438: 						     Functions::ZeroFunction<dim>(m_n_components),
3439: 						     m_constraints,
3440: 						     m_fe.component_mask(y_displacement));
3441: 
3442: 	    const double time_inc = m_time.get_delta_t();
3443: 	    double disp_magnitude = m_time.get_magnitude();
3444: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3445: 						     boundary_id_top_surface,
3446: 						     Functions::ConstantFunction<dim>(
3447: 						       disp_magnitude*time_inc, m_n_components),
3448: 						     m_constraints,
3449: 						     m_fe.component_mask(x_displacement));
3450: 
3451: 	    const int boundary_id_side_surfaces = 2;
3452: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3453: 						     boundary_id_side_surfaces,
3454: 						     Functions::ZeroFunction<dim>(m_n_components),
3455: 						     m_constraints,
3456: 						     m_fe.component_mask(y_displacement));
3457: 	  }
3458: 	else if (m_parameters.m_scenario == 5)
3459: 	  {
3460: 	    typename Triangulation<dim>::active_vertex_iterator vertex_itr;
3461: 	    vertex_itr = m_triangulation.begin_active_vertex();
3462: 	    std::vector<types::global_dof_index> node_bottomleft(m_fe.dofs_per_vertex);
3463: 	    std::vector<types::global_dof_index> node_bottomright(m_fe.dofs_per_vertex);
3464: 	    std::vector<types::global_dof_index> node_topcenter(m_fe.dofs_per_vertex);
3465: 
3466: 	    for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)
3467: 	      {
3468: 		if (   (std::fabs(vertex_itr->vertex()[0] - 0.0) < 1.0e-9)
3469: 		    && (std::fabs(vertex_itr->vertex()[1] - 0.0) < 1.0e-9) )
3470: 		  {
3471: 		    node_bottomleft = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);
3472: 		  }
3473: 		if (   (std::fabs(vertex_itr->vertex()[0] - 8.0) < 1.0e-9)
3474: 		    && (std::fabs(vertex_itr->vertex()[1] - 0.0) < 1.0e-9) )
3475: 		  {
3476: 		    node_bottomright = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);
3477: 		  }
3478: 		if (   (std::fabs(vertex_itr->vertex()[0] - 4.0) < 1.0e-9)
3479: 		    && (std::fabs(vertex_itr->vertex()[1] - 2.0) < 1.0e-9) )
3480: 		  {
3481: 		    node_topcenter = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);
3482: 		  }
3483: 	      }
3484: 	    // bottom-left node fixed in both x- and y-directions
3485: 	    m_constraints.add_line(node_bottomleft[0]);
3486: 	    m_constraints.set_inhomogeneity(node_bottomleft[0], 0.0);
3487: 
3488: 	    m_constraints.add_line(node_bottomleft[1]);
3489: 	    m_constraints.set_inhomogeneity(node_bottomleft[1], 0.0);
3490: 
3491: 	    // bottom-right node only fixed in y-direction
3492: 	    m_constraints.add_line(node_bottomright[1]);
3493: 	    m_constraints.set_inhomogeneity(node_bottomright[1], 0.0);
3494: 
3495: 	    // top-center node applied with y-displacement
3496: 	    const double time_inc = m_time.get_delta_t();
3497: 	    double disp_magnitude = m_time.get_magnitude();
3498: 
3499: 	    m_constraints.add_line(node_topcenter[1]);
3500: 	    m_constraints.set_inhomogeneity(node_topcenter[1], disp_magnitude*time_inc);
3501: 	  }
3502: 	else if (   m_parameters.m_scenario == 6
3503: 	         || m_parameters.m_scenario == 7
3504: 		 || m_parameters.m_scenario == 8)
3505: 	  {
3506: 	    const int x0_surface = 0;
3507: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3508: 						     x0_surface,
3509: 						     Functions::ZeroFunction<dim>(m_n_components),
3510: 						     m_constraints,
3511: 						     m_fe.component_mask(x_displacement));
3512: 	    const int y0_surface = 1;
3513: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3514: 						     y0_surface,
3515: 						     Functions::ZeroFunction<dim>(m_n_components),
3516: 						     m_constraints,
3517: 						     m_fe.component_mask(y_displacement));
3518: 	    const int z0_surface = 2;
3519: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3520: 						     z0_surface,
3521: 						     Functions::ZeroFunction<dim>(m_n_components),
3522: 						     m_constraints,
3523: 						     m_fe.component_mask(z_displacement));
3524: 
3525: 	    const int z1_surface = 3;
3526: 	    const double time_inc = m_time.get_delta_t();
3527: 	    double disp_magnitude = 1.0;
3528: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3529: 						     z1_surface,
3530: 						     Functions::ConstantFunction<dim>(
3531: 						       disp_magnitude*time_inc, m_n_components),
3532: 						     m_constraints,
3533: 						     m_fe.component_mask(z_displacement));
3534: 	  }
3535: 	else if (   m_parameters.m_scenario == 9
3536: 	         || m_parameters.m_scenario == 10)
3537: 	  {
3538: 	    // Dirichlet B,C. bottom surface
3539: 	    const int boundary_id_bottom_surface = 0;
3540: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3541: 						     boundary_id_bottom_surface,
3542: 						     Functions::ZeroFunction<dim>(m_n_components),
3543: 						     m_constraints,
3544: 						     m_fe.component_mask(displacements));
3545: 
3546: 	    typename Triangulation<dim>::active_vertex_iterator vertex_itr;
3547: 	    vertex_itr = m_triangulation.begin_active_vertex();
3548: 	    std::vector<types::global_dof_index> node_disp_control(m_fe.dofs_per_vertex);
3549: 
3550: 	    for (; vertex_itr != m_triangulation.end_vertex(); ++vertex_itr)
3551: 	      {
3552: 		if (   (std::fabs(vertex_itr->vertex()[0] - 470.0) < 1.0e-9)
3553: 		    && (std::fabs(vertex_itr->vertex()[1] - 250.0) < 1.0e-9) )
3554: 		  {
3555: 		    node_disp_control = usr_utilities::get_vertex_dofs(vertex_itr, m_dof_handler);
3556: 	            // node applied with y-displacement
3557: 		    const double time_inc = m_time.get_delta_t();
3558: 		    double disp_magnitude = m_time.get_magnitude();
3559: 
3560: 		    m_constraints.add_line(node_disp_control[1]);
3561: 		    m_constraints.set_inhomogeneity(node_disp_control[1], disp_magnitude*time_inc);
3562: 		  }
3563: 	      }
3564: 	  }
3565: 	else if (m_parameters.m_scenario == 11)
3566: 	  {
3567: 	    // Dirichlet B,C. right surface
3568: 	    const int boundary_id_right_surface = 0;
3569: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3570: 						     boundary_id_right_surface,
3571: 						     Functions::ZeroFunction<dim>(m_n_components),
3572: 						     m_constraints,
3573: 						     m_fe.component_mask(displacements));
3574: 
3575: 	    // Dirichlet B,C. left surface
3576: 	    const int boundary_id_left_surface = 1;
3577: 	    VectorTools::interpolate_boundary_values(m_dof_handler,
3578: 						     boundary_id_left_surface,
3579: 						     Functions::ZeroFunction<dim>(m_n_components),
3580: 						     m_constraints,
3581: 						     m_fe.component_mask(x_displacement));
3582: 
3583: 	    typename Triangulation<dim>::active_vertex_iterator vertex_itr;
3584: 	    vertex_itr = m_triangulation.begin_active_vertex();
3585: 	    std::vector<types::global_dof_index> node_rotate(m_fe.dofs_per_vertex);
3586: 	    double node_dist = 0.0;
3587: 	    double disp_mag = 0.0;
3588: 	    double angle_theta = 0.0;
3589: 	    double disp_y = 0;
3590: 	    double disp_z = 0;
3591: 
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
3619: 		  }
3620: 	      }
3621: 	  }
3622: 	else
3623: 	  Assert(false, ExcMessage("The scenario has not been implemented!"));
3624:       }
3625:     else  // inhomogeneous constraints
3626:       {
3627:         if (m_constraints.has_inhomogeneities())
3628:           {
3629:             AffineConstraints<double> homogeneous_constraints(m_constraints);
3630:             for (unsigned int dof = 0; dof != m_dof_handler.n_dofs(); ++dof)
3631:               if (homogeneous_constraints.is_inhomogeneously_constrained(dof))
3632:                 homogeneous_constraints.set_inhomogeneity(dof, 0.0);
3633:             m_constraints.clear();
3634:             m_constraints.copy_from(homogeneous_constraints);
3635:           }
3636:       }
3637:     m_constraints.close();
3638:   }
```
  - C^T(·) 的全局分配实现：[main.cc:3663-3668](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/main.cc#L3663-L3668)
```cpp
3663:     auto copier = [this](const PerTaskData_ASM &data)
3664:       {
3665:         this->m_constraints.distribute_local_to_global(data.m_cell_matrix,
3666:                                                        data.m_local_dof_indices,
3667: 						       m_tangent_matrix);
3668:       };
```
### Eq. (15) — BFGS 二次模型
- 论文位置：[3. Gradient projection based monolithic scheme > 3.1. Algorithm overview，`explain.md` 第 259 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L259)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
m_{k}(\pmb{x}) = \Pi_{k} + \pmb{r}_{k}^{\mathrm{T}}(\pmb{x} - \pmb{x}_{k}) + \frac{1}{2} (\pmb{x} - \pmb{x}_{k})^{\mathrm{T}}\mathbf{B}_{k}(\pmb{x} - \pmb{x}_{k}) \quad (15)
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
### Eq. (16) — 分量投影算子
- 论文位置：[3. Gradient projection based monolithic scheme > 3.1. Algorithm overview，`explain.md` 第 267 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L267)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\mathrm{Proj}_c(x_i,\mathrm{lb}_i,\mathrm{ub}_i) = \left\{ \begin{array}{ll}\mathrm{lb}_i & \mathrm{if}\; x_i < \mathrm{lb}_i,\\ x_i & \mathrm{if}\; \mathrm{lb}_i \leq x_i \leq \mathrm{ub}_i,\\ \mathrm{ub}_i & \mathrm{if}\; x_i > \mathrm{ub}_i, \end{array} \right. \quad (16)
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
### Eq. (17) — 投影梯度路径
- 论文位置：[3. Gradient projection based monolithic scheme > 3.1. Algorithm overview，`explain.md` 第 273 行附近](https://github.com/SiyaoCao/Phasefield_gradient_projection_monolithic_solver/blob/main/explain.md#L273)
- 公式解读：该式用于描述本章节对应的变分量、离散量或约束/算法步骤，代码中通过下述函数在单元积分、向量更新、约束投影或线搜索阶段执行。
- 公式（按要求渲染）：
\[
\pmb {x}(t) = \mathrm{Proj}_c(\pmb {x}_k - t\pmb {r}_k,\mathbf{lb},\mathbf{ub}),\quad t > 0, \quad (17)
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
