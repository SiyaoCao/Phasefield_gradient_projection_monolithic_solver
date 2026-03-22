# explain.md 与 main.cc 公式-代码对应（公式块 009）

- 所属章节：`2.1. Phase-field formulation`
- explain.md 行号：`93`

论文公式：

\[
g(d) = (1 - d)^{2}. \quad (6)
\]

对应 `main.cc` 代码：

```cpp
// main.cc:244-247
  double degradation_function(const double d)
  {
    return (1.0 - d) * (1.0 - d);
  }
```

简要说明：对应变量或算子在 main.cc 中有同名或等价离散实现。
