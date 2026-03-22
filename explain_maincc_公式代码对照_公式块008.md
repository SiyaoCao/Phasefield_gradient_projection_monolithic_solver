# explain.md 与 main.cc 公式-代码对应（公式块 008）

- 所属章节：`2.1. Phase-field formulation`
- explain.md 行号：`87`

论文公式：

\[
\psi(\pmb{\epsilon}, d) = [g(d) + k] \psi^{+}(\pmb{\epsilon}) + \psi^{-}(\pmb{\epsilon}), \quad (5)
\]

对应 `main.cc` 代码：

```cpp
// main.cc:846-847
    const double degradation = degradation_function(m_phase_field_value) + m_residual_k;
    const double I_1 = trace(m_strain);

// main.cc:864-864
    m_stress = degradation * stress_positive + stress_negative;
```

简要说明：对应变量或算子在 main.cc 中有同名或等价离散实现。
