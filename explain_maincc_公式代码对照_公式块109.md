# explain.md 与 main.cc 公式-代码对应（公式块 109）

- 所属章节：`Appendix`
- explain.md 行号：`973`

论文公式：

\[
\hat{\mathbf{A}} = \mathbf{A} + \mathbf{U}\mathbf{V}^{\mathrm{T}},
\]

对应 `main.cc` 代码：

```cpp
// main.cc:5564-5621
		FullMatrix<double> wT_z_zT_B0_z_inv_zT_w(2 * list_size);
		const auto itr_zT_y_list_begin = zT_y_list.begin();
		const auto itr_zT_b0xs_list_begin = zT_b0xs_list.begin();
		const auto itr_zT_B0_z_inv_zT_y_list_begin = zT_B0_z_inv_zT_y_list.begin();
		const auto itr_zT_B0_z_inv_zT_b0xs_list_begin = zT_B0_z_inv_zT_b0xs_list.begin();
		for (unsigned int i = 0; i < list_size; ++i)
		  for (unsigned int j = 0; j < list_size; ++j)
		    {
		      wT_z_zT_B0_z_inv_zT_w(i          , j          ) = (*std::next(itr_zT_y_list_begin            , i))
								      * (*std::next(itr_zT_B0_z_inv_zT_y_list_begin, j));

		      wT_z_zT_B0_z_inv_zT_w(i          , j+list_size) = (*std::next(itr_zT_y_list_begin               , i))
								      * (*std::next(itr_zT_B0_z_inv_zT_b0xs_list_begin, j));

		      wT_z_zT_B0_z_inv_zT_w(i+list_size, j          ) = (*std::next(itr_zT_b0xs_list_begin         , i))
								      * (*std::next(itr_zT_B0_z_inv_zT_y_list_begin, j));

		      wT_z_zT_B0_z_inv_zT_w(i+list_size, j+list_size) = (*std::next(itr_zT_b0xs_list_begin            , i))
								      * (*std::next(itr_zT_B0_z_inv_zT_b0xs_list_begin, j));
		    }

		FullMatrix<double> temp_matrix(2 * list_size);
		M_matrix.mmult(temp_matrix, wT_z_zT_B0_z_inv_zT_w);

		FullMatrix<double> middle_matrix(IdentityMatrix(2*list_size));
		middle_matrix.add(-1.0, temp_matrix);

		FullMatrix<double> middle_matrix_inv(2 * list_size);
		middle_matrix_inv.invert(middle_matrix);

		middle_matrix_inv.mmult(middle_matrix, M_matrix);

		Vector<double> wT_z_zT_B0_z_inv_rhs(2 * list_size);
		for (unsigned int i = 0; i < list_size; ++i)
		  {
		    wT_z_zT_B0_z_inv_rhs(i            ) = (*std::next(itr_zT_B0_z_inv_zT_y_list_begin   , i)) * rhs_vector;
		    wT_z_zT_B0_z_inv_rhs(i + list_size) = (*std::next(itr_zT_B0_z_inv_zT_b0xs_list_begin, i)) * rhs_vector;
		  }

		Vector<double> middle_matrix_wT_z_zT_B0_z_inv_rhs(2 * list_size);
		middle_matrix.vmult(middle_matrix_wT_z_zT_B0_z_inv_rhs,
				    wT_z_zT_B0_z_inv_rhs);

		unsigned int index = 0;
		for (auto itr = zT_B0_z_inv_zT_y_list.begin(); itr != zT_B0_z_inv_zT_y_list.end(); ++itr)
		  {
		    update_vector.add(middle_matrix_wT_z_zT_B0_z_inv_rhs(index), *itr);
		    ++index;
		  }
		for (auto itr = zT_B0_z_inv_zT_b0xs_list.begin(); itr != zT_B0_z_inv_zT_b0xs_list.end(); ++itr)
		  {
		    update_vector.add(middle_matrix_wT_z_zT_B0_z_inv_rhs(index), *itr);
		    ++index;
		  }
	      } //	if (list_size > 0)

	    search_direction += update_vector;
```

简要说明：对应附录中的低秩修正逆关系（代码中以矩阵恒等式实现）。
