# explain.md 与 main.cc 公式-代码对应（公式块 107）

- 所属章节：`4.4. Three-dimensional torsion test`
- explain.md 行号：`878`

论文公式：

\[
u_{y} = z \tan t, \quad u_{z} = -y \tan t,
\]

对应 `main.cc` 代码：

```cpp
// main.cc:3599-3618

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
```

简要说明：对应三维扭转边界位移函数。
