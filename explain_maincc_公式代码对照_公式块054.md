# explain.md 与 main.cc 公式-代码对应（公式块 054）

- 所属章节：`3.2. Compact representation of limited-memory BFGS matrix`
- explain.md 行号：`408`

论文公式：

\[
\mathbf{Y}_k = [\mathbf{y}_{k-m} \cdots \mathbf{y}_{k-1}].
\]

对应 `main.cc` 代码：

```cpp
// main.cc:5713-5729
	// s_vector_list, y_vector_list, s_dot_y_list only need to discard
	// the front (oldest) item and add the newest item to the end at
	// each L-BFGS iteration
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

简要说明：对应变量或算子在 main.cc 中有同名或等价离散实现。
