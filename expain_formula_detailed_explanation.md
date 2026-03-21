# explain.md 论文公式逐一详解（含推导与内在联系）

> 说明：本文档对 `/home/runner/work/Phasefield_gradient_projection_monolithic_solver/Phasefield_gradient_projection_monolithic_solver/explain.md` 中检测到的全部公式块进行逐一解释。每条都包含：出现位置、完整表达、公式含义、必要推导过程、与其他公式的内在联系。
>
> 统一记号约定（均按你要求用 `[ ... ]` 方式渲染）：
> - 位移场：[ \mathbf{u} ]，离散自由度向量：[ \mathbf{u}_A ]。
> - 相场/损伤变量：[ d ]，历史下界：[ d^{(n)} ]。
> - 总势能：[ \Pi ]，残量：[ \mathbf{r} ]，切线/近似 Hessian：[ \mathbf{K} ] 或 [ \mathbf{B}_k ]。
> - 盒约束上下界：[ \mathbf{lb}, \mathbf{ub} ]。

## 全局逻辑主线（先给“公式网络图”）

1. **物理建模层**：由总势能泛函 [ \Pi(\mathbf{u},d) ]（公式(1)）与相场裂纹表面密度 [ \Gamma_l(d) ]（公式(2)）构成最优化目标。
2. **不可逆约束层**：通过 [ 0\le d_n\le d_{n+1}\le 1 ]（公式(4)）转化为离散盒约束 [ d_A^{(n)}\le d_A\le 1 ]，形成有约束最优化问题。
3. **离散与线性化层**：由变分残量 [ r_u,r_d ]（公式(7)）到有限元离散残量/切线（公式(9)(10)），再到块对角初值矩阵（公式(11)(21)）。
4. **算法层**：构造二次模型 [ m_k ]（公式(15)），先求广义 Cauchy 点（公式(17)(18)(22)–(24)），再做子空间最小化（公式(25)–(32)）。
5. **工程量与附录层**：示例中的能量与边界位移关系（公式(33)等）用于结果解释；附录 Sherman–Morrison–Woodbury 公式支撑快速逆更新。


## 逐条公式解释


### 公式 001

- **出现位置**：`explain.md` 第 21-23 行；章节：[## 1. Introduction]；编号标注：(1)。
- **公式完整表达**：
[
\Pi (\pmb {u},d) = \int_{\Omega}\psi (\pmb {e}(\pmb {u}),d)\mathrm{d}\Omega + g_{c}\Gamma_{l}(d) - \int_{\Omega}\pmb {b}\cdot \pmb {u}\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot \pmb {u}\mathrm{d}\Gamma , \quad (1)
]
- **公式含义**：该式处于问题定义层：用于给出相场断裂的目标泛函、约束或其优化表述。核心是把力学平衡与断裂演化统一为对 [ \Pi ] 的最小化。
- **必要推导过程**：
- 由 Griffith 断裂能思想：总势能 = 体能 + 裂纹表面能 - 外力势。
- 用相场正则化把离散裂纹面积替换为体积分 [ \Gamma_l(d) ]，得到可微近似。
- 施加不可逆条件后形成带盒约束的时序最小化问题。
- **与其他公式的内在联系**：
- 与公式(1)–(4)共同定义“带不可逆约束的能量最小化”总问题。
- 该层结果直接进入第 2 节本构与第 3 节数值优化。

### 公式 002

- **出现位置**：`explain.md` 第 27-29 行；章节：[## 1. Introduction]；编号标注：(2)。
- **公式完整表达**：
[
\Gamma_{l}(d) = \int_{\Omega}\gamma (d,\nabla d)\mathrm{d}\Omega = \int_{\Omega}\frac{1}{2l}\left(d^{2} + l^{2}\nabla d\cdot \nabla d\right)\mathrm{d}\Omega , \quad (2)
]
- **公式含义**：该式处于问题定义层：用于给出相场断裂的目标泛函、约束或其优化表述。核心是把力学平衡与断裂演化统一为对 [ \Pi ] 的最小化。
- **必要推导过程**：
- 由 Griffith 断裂能思想：总势能 = 体能 + 裂纹表面能 - 外力势。
- 用相场正则化把离散裂纹面积替换为体积分 [ \Gamma_l(d) ]，得到可微近似。
- 施加不可逆条件后形成带盒约束的时序最小化问题。
- **与其他公式的内在联系**：
- 与公式(1)–(4)共同定义“带不可逆约束的能量最小化”总问题。
- 该层结果直接进入第 2 节本构与第 3 节数值优化。

### 公式 003

- **出现位置**：`explain.md` 第 35-37 行；章节：[## 1. Introduction]；编号标注：(3)。
- **公式完整表达**：
[
\left(\pmb{u}_{n + 1},d_{n + 1}\right) = \arg \min \Pi (\pmb {u},d) \quad (3)
]
- **公式含义**：该式处于问题定义层：用于给出相场断裂的目标泛函、约束或其优化表述。核心是把力学平衡与断裂演化统一为对 [ \Pi ] 的最小化。
- **必要推导过程**：
- 由 Griffith 断裂能思想：总势能 = 体能 + 裂纹表面能 - 外力势。
- 用相场正则化把离散裂纹面积替换为体积分 [ \Gamma_l(d) ]，得到可微近似。
- 施加不可逆条件后形成带盒约束的时序最小化问题。
- **与其他公式的内在联系**：
- 与公式(1)–(4)共同定义“带不可逆约束的能量最小化”总问题。
- 该层结果直接进入第 2 节本构与第 3 节数值优化。

### 公式 004

- **出现位置**：`explain.md` 第 41-43 行；章节：[## 1. Introduction]；编号标注：(4)。
- **公式完整表达**：
[
0\leq d_{n}\leq d_{n + 1}\leq 1. \quad (4)
]
- **公式含义**：该式处于问题定义层：用于给出相场断裂的目标泛函、约束或其优化表述。核心是把力学平衡与断裂演化统一为对 [ \Pi ] 的最小化。
- **必要推导过程**：
- 由 Griffith 断裂能思想：总势能 = 体能 + 裂纹表面能 - 外力势。
- 用相场正则化把离散裂纹面积替换为体积分 [ \Gamma_l(d) ]，得到可微近似。
- 施加不可逆条件后形成带盒约束的时序最小化问题。
- **与其他公式的内在联系**：
- 与公式(1)–(4)共同定义“带不可逆约束的能量最小化”总问题。
- 该层结果直接进入第 2 节本构与第 3 节数值优化。

### 公式 005

- **出现位置**：`explain.md` 第 57-59 行；章节：[## 1. Introduction]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\mathrm{Proj}_C(\mathbf{x}_k = \mathbf{a}_k\nabla f(\mathbf{x}_k)).
]
- **公式含义**：该式处于问题定义层：用于给出相场断裂的目标泛函、约束或其优化表述。核心是把力学平衡与断裂演化统一为对 [ \Pi ] 的最小化。
- **必要推导过程**：
- 由 Griffith 断裂能思想：总势能 = 体能 + 裂纹表面能 - 外力势。
- 用相场正则化把离散裂纹面积替换为体积分 [ \Gamma_l(d) ]，得到可微近似。
- 施加不可逆条件后形成带盒约束的时序最小化问题。
- **与其他公式的内在联系**：
- 与公式(1)–(4)共同定义“带不可逆约束的能量最小化”总问题。
- 该层结果直接进入第 2 节本构与第 3 节数值优化。

### 公式 006

- **出现位置**：`explain.md` 第 63-65 行；章节：[## 1. Introduction]；编号标注：(4)。
- **公式完整表达**：
[
\mathrm{lb}_i\leq x_i\leq \mathrm{ub}_i,
]
- **公式含义**：该式处于问题定义层：用于给出相场断裂的目标泛函、约束或其优化表述。核心是把力学平衡与断裂演化统一为对 [ \Pi ] 的最小化。
- **必要推导过程**：
- 由 Griffith 断裂能思想：总势能 = 体能 + 裂纹表面能 - 外力势。
- 用相场正则化把离散裂纹面积替换为体积分 [ \Gamma_l(d) ]，得到可微近似。
- 施加不可逆条件后形成带盒约束的时序最小化问题。
- **与其他公式的内在联系**：
- 与公式(1)–(4)共同定义“带不可逆约束的能量最小化”总问题。
- 该层结果直接进入第 2 节本构与第 3 节数值优化。

### 公式 007

- **出现位置**：`explain.md` 第 69-71 行；章节：[## 1. Introduction]；编号标注：(4)。
- **公式完整表达**：
[
\mathrm{lb}_i = d_i^{(n)}\leq d_i^{(n + 1)}\leq 1 = \mathrm{ub}_i,
]
- **公式含义**：该式处于问题定义层：用于给出相场断裂的目标泛函、约束或其优化表述。核心是把力学平衡与断裂演化统一为对 [ \Pi ] 的最小化。
- **必要推导过程**：
- 由 Griffith 断裂能思想：总势能 = 体能 + 裂纹表面能 - 外力势。
- 用相场正则化把离散裂纹面积替换为体积分 [ \Gamma_l(d) ]，得到可微近似。
- 施加不可逆条件后形成带盒约束的时序最小化问题。
- **与其他公式的内在联系**：
- 与公式(1)–(4)共同定义“带不可逆约束的能量最小化”总问题。
- 该层结果直接进入第 2 节本构与第 3 节数值优化。

### 公式 008

- **出现位置**：`explain.md` 第 87-89 行；章节：[### 2.1. Phase-field formulation]；编号标注：(1)、(5)。
- **公式完整表达**：
[
\psi(\pmb{\epsilon}, d) = [g(d) + k] \psi^{+}(\pmb{\epsilon}) + \psi^{-}(\pmb{\epsilon}), \quad (5)
]
- **公式含义**：该式处于本构层：定义应变能张拉/压缩分裂、退化函数、应力与切线算子，是“材料响应如何进入能量与残量”的关键桥梁。
- **必要推导过程**：
- 从小应变假设 [ \boldsymbol{\epsilon}=\nabla^{(s)}\mathbf{u} ] 出发，对主应变做谱分解。
- 通过 Macaulay 括号把张拉/压缩能分离，令退化函数只作用于张拉部分，避免压缩下虚假损伤。
- 对能量对 [ \boldsymbol{\epsilon} ] 求导得应力，再求导得一致切线。
- **与其他公式的内在联系**：
- 向前接公式(1)(2)：把能量密度 [ \psi ] 具体化。
- 向后接公式(7)(9)(10)：为变分残量与切线矩阵提供闭式表达。

### 公式 009

- **出现位置**：`explain.md` 第 93-95 行；章节：[### 2.1. Phase-field formulation]；编号标注：(6)。
- **公式完整表达**：
[
g(d) = (1 - d)^{2}. \quad (6)
]
- **公式含义**：该式处于本构层：定义应变能张拉/压缩分裂、退化函数、应力与切线算子，是“材料响应如何进入能量与残量”的关键桥梁。
- **必要推导过程**：
- 从小应变假设 [ \boldsymbol{\epsilon}=\nabla^{(s)}\mathbf{u} ] 出发，对主应变做谱分解。
- 通过 Macaulay 括号把张拉/压缩能分离，令退化函数只作用于张拉部分，避免压缩下虚假损伤。
- 对能量对 [ \boldsymbol{\epsilon} ] 求导得应力，再求导得一致切线。
- **与其他公式的内在联系**：
- 向前接公式(1)(2)：把能量密度 [ \psi ] 具体化。
- 向后接公式(7)(9)(10)：为变分残量与切线矩阵提供闭式表达。

### 公式 010

- **出现位置**：`explain.md` 第 101-103 行；章节：[### 2.1. Phase-field formulation]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\langle x\rangle_{+} = \frac{1}{2}(x + |x|), \quad \langle x\rangle_{-} = \frac{1}{2}(x - |x|), \quad H(x) = \begin{cases} 1 & x \geq 0, \\ 0 & x < 0. \end{cases}
]
- **公式含义**：该式处于本构层：定义应变能张拉/压缩分裂、退化函数、应力与切线算子，是“材料响应如何进入能量与残量”的关键桥梁。
- **必要推导过程**：
- 从小应变假设 [ \boldsymbol{\epsilon}=\nabla^{(s)}\mathbf{u} ] 出发，对主应变做谱分解。
- 通过 Macaulay 括号把张拉/压缩能分离，令退化函数只作用于张拉部分，避免压缩下虚假损伤。
- 对能量对 [ \boldsymbol{\epsilon} ] 求导得应力，再求导得一致切线。
- **与其他公式的内在联系**：
- 向前接公式(1)(2)：把能量密度 [ \psi ] 具体化。
- 向后接公式(7)(9)(10)：为变分残量与切线矩阵提供闭式表达。

### 公式 011

- **出现位置**：`explain.md` 第 107-109 行；章节：[### 2.1. Phase-field formulation]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\pmb{\epsilon} = \sum_{\alpha} \epsilon_{\alpha} \mathbf{M}_{\alpha}, \quad \mathbf{M}_{\alpha} = \pmb{n}_{\alpha} \otimes \pmb{n}_{\alpha},
]
- **公式含义**：该式处于本构层：定义应变能张拉/压缩分裂、退化函数、应力与切线算子，是“材料响应如何进入能量与残量”的关键桥梁。
- **必要推导过程**：
- 从小应变假设 [ \boldsymbol{\epsilon}=\nabla^{(s)}\mathbf{u} ] 出发，对主应变做谱分解。
- 通过 Macaulay 括号把张拉/压缩能分离，令退化函数只作用于张拉部分，避免压缩下虚假损伤。
- 对能量对 [ \boldsymbol{\epsilon} ] 求导得应力，再求导得一致切线。
- **与其他公式的内在联系**：
- 向前接公式(1)(2)：把能量密度 [ \psi ] 具体化。
- 向后接公式(7)(9)(10)：为变分残量与切线矩阵提供闭式表达。

### 公式 012

- **出现位置**：`explain.md` 第 113-115 行；章节：[### 2.1. Phase-field formulation]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\pmb{\epsilon}^{+} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{+} \mathbf{M}_{\alpha}, \quad \pmb{\epsilon}^{-} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{-} \mathbf{M}_{\alpha}.
]
- **公式含义**：该式处于本构层：定义应变能张拉/压缩分裂、退化函数、应力与切线算子，是“材料响应如何进入能量与残量”的关键桥梁。
- **必要推导过程**：
- 从小应变假设 [ \boldsymbol{\epsilon}=\nabla^{(s)}\mathbf{u} ] 出发，对主应变做谱分解。
- 通过 Macaulay 括号把张拉/压缩能分离，令退化函数只作用于张拉部分，避免压缩下虚假损伤。
- 对能量对 [ \boldsymbol{\epsilon} ] 求导得应力，再求导得一致切线。
- **与其他公式的内在联系**：
- 向前接公式(1)(2)：把能量密度 [ \psi ] 具体化。
- 向后接公式(7)(9)(10)：为变分残量与切线矩阵提供闭式表达。

### 公式 013

- **出现位置**：`explain.md` 第 119-121 行；章节：[### 2.1. Phase-field formulation]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\psi^{+}(\pmb{\epsilon}) = \frac{1}{2}\lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{+}^{2} + \mu \pmb{\epsilon}^{+} : \pmb{\epsilon}^{+}, \quad \psi^{-}(\pmb{\epsilon}) = \frac{1}{2} \lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{-}^{2} + \mu \pmb{\epsilon}^{-} : \pmb{\epsilon}^{-},
]
- **公式含义**：该式处于本构层：定义应变能张拉/压缩分裂、退化函数、应力与切线算子，是“材料响应如何进入能量与残量”的关键桥梁。
- **必要推导过程**：
- 从小应变假设 [ \boldsymbol{\epsilon}=\nabla^{(s)}\mathbf{u} ] 出发，对主应变做谱分解。
- 通过 Macaulay 括号把张拉/压缩能分离，令退化函数只作用于张拉部分，避免压缩下虚假损伤。
- 对能量对 [ \boldsymbol{\epsilon} ] 求导得应力，再求导得一致切线。
- **与其他公式的内在联系**：
- 向前接公式(1)(2)：把能量密度 [ \psi ] 具体化。
- 向后接公式(7)(9)(10)：为变分残量与切线矩阵提供闭式表达。

### 公式 014

- **出现位置**：`explain.md` 第 125-127 行；章节：[### 2.1. Phase-field formulation]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\pmb{\sigma} = \frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial \pmb{\epsilon}} = [g(d) + k]\frac{\partial\psi^{+}(\pmb{\epsilon})}{\partial \pmb{\epsilon}} + \frac{\partial\psi^{-}(\pmb{\epsilon})}{\partial \pmb{\epsilon}} = [g(d) + k]\pmb{\sigma}^{+} + \pmb{\sigma}^{-},
]
- **公式含义**：该式处于本构层：定义应变能张拉/压缩分裂、退化函数、应力与切线算子，是“材料响应如何进入能量与残量”的关键桥梁。
- **必要推导过程**：
- 从小应变假设 [ \boldsymbol{\epsilon}=\nabla^{(s)}\mathbf{u} ] 出发，对主应变做谱分解。
- 通过 Macaulay 括号把张拉/压缩能分离，令退化函数只作用于张拉部分，避免压缩下虚假损伤。
- 对能量对 [ \boldsymbol{\epsilon} ] 求导得应力，再求导得一致切线。
- **与其他公式的内在联系**：
- 向前接公式(1)(2)：把能量密度 [ \psi ] 具体化。
- 向后接公式(7)(9)(10)：为变分残量与切线矩阵提供闭式表达。

### 公式 015

- **出现位置**：`explain.md` 第 131-133 行；章节：[### 2.1. Phase-field formulation]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\pmb{\sigma}^{+} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{+}\mathbf{I} + 2\mu \pmb{\epsilon}^{+}, \quad \pmb{\sigma}^{-} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{-}\mathbf{I} + 2\mu \pmb{\epsilon}^{-}.
]
- **公式含义**：该式处于本构层：定义应变能张拉/压缩分裂、退化函数、应力与切线算子，是“材料响应如何进入能量与残量”的关键桥梁。
- **必要推导过程**：
- 从小应变假设 [ \boldsymbol{\epsilon}=\nabla^{(s)}\mathbf{u} ] 出发，对主应变做谱分解。
- 通过 Macaulay 括号把张拉/压缩能分离，令退化函数只作用于张拉部分，避免压缩下虚假损伤。
- 对能量对 [ \boldsymbol{\epsilon} ] 求导得应力，再求导得一致切线。
- **与其他公式的内在联系**：
- 向前接公式(1)(2)：把能量密度 [ \psi ] 具体化。
- 向后接公式(7)(9)(10)：为变分残量与切线矩阵提供闭式表达。

### 公式 016

- **出现位置**：`explain.md` 第 137-139 行；章节：[### 2.1. Phase-field formulation]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\frac{\partial \pmb{\sigma}}{\partial \pmb{\epsilon}} = [g(d) + k]\frac{\partial \pmb{\sigma}^{+}}{\partial \pmb{\epsilon}} + \frac{\partial \pmb{\sigma}^{-}}{\partial \pmb{\epsilon}} = [g(d) + k]\left[\lambda H(\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes \mathbf{I} + 2\mu \mathbb{P}^{+}\right] + \left[\lambda H(-\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes \mathbf{I} + 2\mu \mathbb{P}^{-}\right].
]
- **公式含义**：该式处于本构层：定义应变能张拉/压缩分裂、退化函数、应力与切线算子，是“材料响应如何进入能量与残量”的关键桥梁。
- **必要推导过程**：
- 从小应变假设 [ \boldsymbol{\epsilon}=\nabla^{(s)}\mathbf{u} ] 出发，对主应变做谱分解。
- 通过 Macaulay 括号把张拉/压缩能分离，令退化函数只作用于张拉部分，避免压缩下虚假损伤。
- 对能量对 [ \boldsymbol{\epsilon} ] 求导得应力，再求导得一致切线。
- **与其他公式的内在联系**：
- 向前接公式(1)(2)：把能量密度 [ \psi ] 具体化。
- 向后接公式(7)(9)(10)：为变分残量与切线矩阵提供闭式表达。

### 公式 017

- **出现位置**：`explain.md` 第 143-145 行；章节：[### 2.1. Phase-field formulation]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\mathbb{P}^{+} = \frac{\partial \pmb{\epsilon}^{+}}{\partial \pmb{\epsilon}}, \quad \mathbb{P}^{-} = \frac{\partial \pmb{\epsilon}^{-}}{\partial \pmb{\epsilon}},
]
- **公式含义**：该式处于本构层：定义应变能张拉/压缩分裂、退化函数、应力与切线算子，是“材料响应如何进入能量与残量”的关键桥梁。
- **必要推导过程**：
- 从小应变假设 [ \boldsymbol{\epsilon}=\nabla^{(s)}\mathbf{u} ] 出发，对主应变做谱分解。
- 通过 Macaulay 括号把张拉/压缩能分离，令退化函数只作用于张拉部分，避免压缩下虚假损伤。
- 对能量对 [ \boldsymbol{\epsilon} ] 求导得应力，再求导得一致切线。
- **与其他公式的内在联系**：
- 向前接公式(1)(2)：把能量密度 [ \psi ] 具体化。
- 向后接公式(7)(9)(10)：为变分残量与切线矩阵提供闭式表达。

### 公式 018

- **出现位置**：`explain.md` 第 151-158 行；章节：[### 2.1. Phase-field formulation]；编号标注：(7)。
- **公式完整表达**：
[
\begin{array}{rl} 
\delta \Pi (\pmb {u},d) &= D_{(\delta \pmb {u},\delta d)}\Pi (\pmb {u},d) = \left.\frac{\mathrm{d}}{\mathrm{d}\epsilon}\right|_{\epsilon = 0}\Pi (\pmb {u} + \epsilon \delta \pmb {u},d + \epsilon \delta d) \\ 
&= \int_{\Omega}\left(\frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial \pmb{\epsilon}}:\pmb{\epsilon}(\delta \pmb {u}) + \frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial d}\delta d\right)\mathrm{d}\Omega \\ 
&\quad +\int_{\Omega}\frac{g_{c}}{l}\left(d\delta d + l^{2}\nabla d\cdot \nabla \delta d\right)\mathrm{d}\Omega -\int_{\Omega}\pmb {b}\cdot \delta \pmb {u}\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot \delta \pmb {u}\mathrm{d}\Gamma \\ 
&= (\nabla^{(s)}\delta \pmb {u},\pmb {\sigma}) - (\delta \pmb {u},\pmb {b}) - (\delta \pmb {u},\pmb{t})_{\Gamma_{t}} + (\delta d,\frac{g_{c}}{l} d) + (\nabla \delta d,g_{c}l\nabla d) + (\delta d,g^{\prime}(d)\psi^{+}(\pmb{\epsilon})). 
\end{array} \quad (7)
]
- **公式含义**：该式处于本构层：定义应变能张拉/压缩分裂、退化函数、应力与切线算子，是“材料响应如何进入能量与残量”的关键桥梁。
- **必要推导过程**：
- 从小应变假设 [ \boldsymbol{\epsilon}=\nabla^{(s)}\mathbf{u} ] 出发，对主应变做谱分解。
- 通过 Macaulay 括号把张拉/压缩能分离，令退化函数只作用于张拉部分，避免压缩下虚假损伤。
- 对能量对 [ \boldsymbol{\epsilon} ] 求导得应力，再求导得一致切线。
- **与其他公式的内在联系**：
- 向前接公式(1)(2)：把能量密度 [ \psi ] 具体化。
- 向后接公式(7)(9)(10)：为变分残量与切线矩阵提供闭式表达。

### 公式 019

- **出现位置**：`explain.md` 第 162-164 行；章节：[### 2.1. Phase-field formulation]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\left\{ \begin{array}{ll} r_{\pmb{u}}(\pmb{u},d) = (\nabla^{(s)}\delta \pmb{u},\pmb{\sigma}) - (\delta \pmb{u},\pmb{b}) - (\delta \pmb{u},\pmb{t})_{\Gamma_{t}} = 0, \\ r_{d}(\pmb{u},d) = (\delta d,\frac{g_{c}}{l} d) + (\nabla \delta d,g_{c}l\nabla d) + (\delta d,g^{\prime}(d)\psi^{+}(\pmb{\epsilon})) = 0, \end{array} \right.
]
- **公式含义**：该式处于本构层：定义应变能张拉/压缩分裂、退化函数、应力与切线算子，是“材料响应如何进入能量与残量”的关键桥梁。
- **必要推导过程**：
- 从小应变假设 [ \boldsymbol{\epsilon}=\nabla^{(s)}\mathbf{u} ] 出发，对主应变做谱分解。
- 通过 Macaulay 括号把张拉/压缩能分离，令退化函数只作用于张拉部分，避免压缩下虚假损伤。
- 对能量对 [ \boldsymbol{\epsilon} ] 求导得应力，再求导得一致切线。
- **与其他公式的内在联系**：
- 向前接公式(1)(2)：把能量密度 [ \psi ] 具体化。
- 向后接公式(7)(9)(10)：为变分残量与切线矩阵提供闭式表达。

### 公式 020

- **出现位置**：`explain.md` 第 172-174 行；章节：[### 2.2. Finite element discretization]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\pmb {u} = \pmb{N}_{u_{A}}\pmb{u}_{A}, \quad \text{and} \quad d = N_{d_{A}}d_{A}.
]
- **公式含义**：该式处于有限元离散层：把连续变量 [ \mathbf{u},d ] 投影到离散自由度 [ \mathbf{u}_A,d_A ]，并由此得到离散残量与雅可比/切线矩阵。
- **必要推导过程**：
- 采用形函数插值 [ \mathbf{u}=\mathbf{N}_{u_A}\mathbf{u}_A,\ d=N_{d_A}d_A ]。
- 把插值代入连续泛函后，对节点自由度求偏导得到离散残量。
- 再对残量求导得到切线子块 [ \mathbf{K}_{uu},\mathbf{K}_{ud},\mathbf{K}_{du},\mathbf{K}_{dd} ]。
- **与其他公式的内在联系**：
- 向前接第 2.1 节本构关系。
- 向后接第 3 节优化算法：离散残量 [ \mathbf{r} ] 与切线 [ \mathbf{K} ] 是每次迭代输入。

### 公式 021

- **出现位置**：`explain.md` 第 178-180 行；章节：[### 2.2. Finite element discretization]；编号标注：(1)。
- **公式完整表达**：
[
\delta \pmb {u} = \pmb{N}_{u_{A}}\delta \pmb{u}_{A}, \quad \text{and} \quad \delta d = N_{d_{A}}\delta d_{A},
]
- **公式含义**：该式处于有限元离散层：把连续变量 [ \mathbf{u},d ] 投影到离散自由度 [ \mathbf{u}_A,d_A ]，并由此得到离散残量与雅可比/切线矩阵。
- **必要推导过程**：
- 采用形函数插值 [ \mathbf{u}=\mathbf{N}_{u_A}\mathbf{u}_A,\ d=N_{d_A}d_A ]。
- 把插值代入连续泛函后，对节点自由度求偏导得到离散残量。
- 再对残量求导得到切线子块 [ \mathbf{K}_{uu},\mathbf{K}_{ud},\mathbf{K}_{du},\mathbf{K}_{dd} ]。
- **与其他公式的内在联系**：
- 向前接第 2.1 节本构关系。
- 向后接第 3 节优化算法：离散残量 [ \mathbf{r} ] 与切线 [ \mathbf{K} ] 是每次迭代输入。

### 公式 022

- **出现位置**：`explain.md` 第 184-190 行；章节：[### 2.2. Finite element discretization]；编号标注：(1)、(8)。
- **公式完整表达**：
[
\begin{array}{rl} 
\Pi (\pmb {u}_{A},d_{A}) = &\int_{\Omega}\psi \left(\pmb{\epsilon}(\pmb {N}_{u_{A}}\pmb {u}_{A}),N_{d_{A}}d_{A}\right)\mathrm{d}\Omega \\ 
&+\int_{\Omega}\frac{g_{c}}{2l}\left((N_{d_{A}}d_{A})^{2} + l^{2}(\nabla N_{d_{A}}d_{A})\cdot (\nabla N_{d_{A}}d_{A})\right)\mathrm{d}\Omega \\ 
&-\int_{\Omega}\pmb {b}\cdot (\pmb {N}_{u_{A}}\pmb {u}_{A})\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot (\pmb {N}_{u_{A}}\pmb {u}_{A})\mathrm{d}\Gamma, 
\end{array} \quad (8)
]
- **公式含义**：该式处于有限元离散层：把连续变量 [ \mathbf{u},d ] 投影到离散自由度 [ \mathbf{u}_A,d_A ]，并由此得到离散残量与雅可比/切线矩阵。
- **必要推导过程**：
- 采用形函数插值 [ \mathbf{u}=\mathbf{N}_{u_A}\mathbf{u}_A,\ d=N_{d_A}d_A ]。
- 把插值代入连续泛函后，对节点自由度求偏导得到离散残量。
- 再对残量求导得到切线子块 [ \mathbf{K}_{uu},\mathbf{K}_{ud},\mathbf{K}_{du},\mathbf{K}_{dd} ]。
- **与其他公式的内在联系**：
- 向前接第 2.1 节本构关系。
- 向后接第 3 节优化算法：离散残量 [ \mathbf{r} ] 与切线 [ \mathbf{K} ] 是每次迭代输入。

### 公式 023

- **出现位置**：`explain.md` 第 194-200 行；章节：[### 2.2. Finite element discretization]；编号标注：(9)。
- **公式完整表达**：
[
\begin{array}{rl} 
\pmb{r} = \nabla \Pi = (r_{\pmb{u}},r_{d})^{\mathrm{T}},\\ 
r_{\pmb{u}_{A}} = \frac{\partial\Pi}{\partial \pmb{u}_{A}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},\pmb{\sigma}\right) - \left(\pmb{N}_{u_{A}},\pmb{b}\right) - \left(\pmb{N}_{u_{A}},\pmb{t}\right)_{\Gamma_{t}},\\ 
r_{d_{A}} = \frac{\partial\Pi}{\partial d_{A}} = \left(N_{d_{A}},\frac{g_{c}}{l} d + g^{\prime}(d)\psi^{+}\right) + (\nabla N_{d_{A}},g_{c}l\nabla d). 
\end{array} \quad (9)
]
- **公式含义**：该式处于有限元离散层：把连续变量 [ \mathbf{u},d ] 投影到离散自由度 [ \mathbf{u}_A,d_A ]，并由此得到离散残量与雅可比/切线矩阵。
- **必要推导过程**：
- 采用形函数插值 [ \mathbf{u}=\mathbf{N}_{u_A}\mathbf{u}_A,\ d=N_{d_A}d_A ]。
- 把插值代入连续泛函后，对节点自由度求偏导得到离散残量。
- 再对残量求导得到切线子块 [ \mathbf{K}_{uu},\mathbf{K}_{ud},\mathbf{K}_{du},\mathbf{K}_{dd} ]。
- **与其他公式的内在联系**：
- 向前接第 2.1 节本构关系。
- 向后接第 3 节优化算法：离散残量 [ \mathbf{r} ] 与切线 [ \mathbf{K} ] 是每次迭代输入。

### 公式 024

- **出现位置**：`explain.md` 第 204-206 行；章节：[### 2.2. Finite element discretization]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\mathbf{K} = \nabla^{2}\Pi = \left[ \begin{array}{ll}\mathbf{K}_{uu} & \mathbf{K}_{ud}\\ \mathbf{K}_{du} & \mathbf{K}_{dd} \end{array} \right],
]
- **公式含义**：该式处于有限元离散层：把连续变量 [ \mathbf{u},d ] 投影到离散自由度 [ \mathbf{u}_A,d_A ]，并由此得到离散残量与雅可比/切线矩阵。
- **必要推导过程**：
- 采用形函数插值 [ \mathbf{u}=\mathbf{N}_{u_A}\mathbf{u}_A,\ d=N_{d_A}d_A ]。
- 把插值代入连续泛函后，对节点自由度求偏导得到离散残量。
- 再对残量求导得到切线子块 [ \mathbf{K}_{uu},\mathbf{K}_{ud},\mathbf{K}_{du},\mathbf{K}_{dd} ]。
- **与其他公式的内在联系**：
- 向前接第 2.1 节本构关系。
- 向后接第 3 节优化算法：离散残量 [ \mathbf{r} ] 与切线 [ \mathbf{K} ] 是每次迭代输入。

### 公式 025

- **出现位置**：`explain.md` 第 210-215 行；章节：[### 2.2. Finite element discretization]；编号标注：(10)。
- **公式完整表达**：
[
\begin{array}{rl} 
\mathbf{K}_{u_{A}u_{B}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},\frac{\partial \pmb{\sigma}}{\partial \pmb{\epsilon}}:\nabla^{(s)}\pmb{N}_{u_{B}}\right),\qquad \mathbf{K}_{u_{A}d_{B}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},g^{\prime}(d)\pmb{\sigma}^{+}N_{d_{B}}\right),\\ 
\mathbf{K}_{d_{A}u_{B}} = \left(N_{d_{A}},g^{\prime}(d)\pmb{\sigma}^{+}:\nabla^{(s)}\pmb{N}_{u_{B}}\right),\qquad \mathbf{K}_{d_{A}d_{B}} = \left(N_{d_{A}},\left(\frac{g_{c}}{l} +g^{\prime \prime}(d)\psi^{+}\right)N_{d_{B}}\right) + \left(\nabla N_{d_{A}},g_{c}l\nabla N_{d_{B}}\right). 
\end{array} \quad (10)
]
- **公式含义**：该式处于有限元离散层：把连续变量 [ \mathbf{u},d ] 投影到离散自由度 [ \mathbf{u}_A,d_A ]，并由此得到离散残量与雅可比/切线矩阵。
- **必要推导过程**：
- 采用形函数插值 [ \mathbf{u}=\mathbf{N}_{u_A}\mathbf{u}_A,\ d=N_{d_A}d_A ]。
- 把插值代入连续泛函后，对节点自由度求偏导得到离散残量。
- 再对残量求导得到切线子块 [ \mathbf{K}_{uu},\mathbf{K}_{ud},\mathbf{K}_{du},\mathbf{K}_{dd} ]。
- **与其他公式的内在联系**：
- 向前接第 2.1 节本构关系。
- 向后接第 3 节优化算法：离散残量 [ \mathbf{r} ] 与切线 [ \mathbf{K} ] 是每次迭代输入。

### 公式 026

- **出现位置**：`explain.md` 第 219-221 行；章节：[### 2.2. Finite element discretization]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
d_A^{(n)}\leq d_A\leq 1,
]
- **公式含义**：该式处于有限元离散层：把连续变量 [ \mathbf{u},d ] 投影到离散自由度 [ \mathbf{u}_A,d_A ]，并由此得到离散残量与雅可比/切线矩阵。
- **必要推导过程**：
- 采用形函数插值 [ \mathbf{u}=\mathbf{N}_{u_A}\mathbf{u}_A,\ d=N_{d_A}d_A ]。
- 把插值代入连续泛函后，对节点自由度求偏导得到离散残量。
- 再对残量求导得到切线子块 [ \mathbf{K}_{uu},\mathbf{K}_{ud},\mathbf{K}_{du},\mathbf{K}_{dd} ]。
- **与其他公式的内在联系**：
- 向前接第 2.1 节本构关系。
- 向后接第 3 节优化算法：离散残量 [ \mathbf{r} ] 与切线 [ \mathbf{K} ] 是每次迭代输入。

### 公式 027

- **出现位置**：`explain.md` 第 225-227 行；章节：[### 2.2. Finite element discretization]；编号标注：(11)。
- **公式完整表达**：
[
\hat{\mathbf{K}} = \left[ \begin{array}{cc}\mathbf{K}_{uu} & \mathbf{0}\\ \mathbf{0} & \mathbf{K}_{dd} \end{array} \right], \quad (11)
]
- **公式含义**：该式处于有限元离散层：把连续变量 [ \mathbf{u},d ] 投影到离散自由度 [ \mathbf{u}_A,d_A ]，并由此得到离散残量与雅可比/切线矩阵。
- **必要推导过程**：
- 采用形函数插值 [ \mathbf{u}=\mathbf{N}_{u_A}\mathbf{u}_A,\ d=N_{d_A}d_A ]。
- 把插值代入连续泛函后，对节点自由度求偏导得到离散残量。
- 再对残量求导得到切线子块 [ \mathbf{K}_{uu},\mathbf{K}_{ud},\mathbf{K}_{du},\mathbf{K}_{dd} ]。
- **与其他公式的内在联系**：
- 向前接第 2.1 节本构关系。
- 向后接第 3 节优化算法：离散残量 [ \mathbf{r} ] 与切线 [ \mathbf{K} ] 是每次迭代输入。

### 公式 028

- **出现位置**：`explain.md` 第 239-241 行；章节：[### 3.1. Algorithm overview]；编号标注：(12)。
- **公式完整表达**：
[
\left(\pmb{u}_A,d_A\right) = \arg \min \Pi (\pmb {u}_A,d_A) \quad (12)
]
- **公式含义**：该式处于求解器总览层：描述约束优化模型、投影算子、线搜索、收敛判据与受限自由度处理，是单调收敛与可行性保持的核心。
- **必要推导过程**：
- 在第 [ k ] 次迭代，用二次模型 [ m_k ] 近似原函数 [ \Pi ]。
- 先做投影得到可行下降方向（广义 Cauchy 点），再在自由子空间内二次优化。
- 用 Wolfe/Armijo 类线搜索确定 [ \alpha_k ]，并以投影残量与步长范数作为停止准则。
- **与其他公式的内在联系**：
- 与公式(15)–(19)形成“二次模型 + 投影 + 活动集 + 线搜索 + 停止判据”闭环。
- 为第 3.2–3.4 节的 L-BFGS、GCP 与子空间法提供框架。

### 公式 029

- **出现位置**：`explain.md` 第 245-247 行；章节：[### 3.1. Algorithm overview]；编号标注：(13)。
- **公式完整表达**：
[
d_A^{(n)}\leq d_A\leq 1 \quad (13)
]
- **公式含义**：该式处于求解器总览层：描述约束优化模型、投影算子、线搜索、收敛判据与受限自由度处理，是单调收敛与可行性保持的核心。
- **必要推导过程**：
- 在第 [ k ] 次迭代，用二次模型 [ m_k ] 近似原函数 [ \Pi ]。
- 先做投影得到可行下降方向（广义 Cauchy 点），再在自由子空间内二次优化。
- 用 Wolfe/Armijo 类线搜索确定 [ \alpha_k ]，并以投影残量与步长范数作为停止准则。
- **与其他公式的内在联系**：
- 与公式(15)–(19)形成“二次模型 + 投影 + 活动集 + 线搜索 + 停止判据”闭环。
- 为第 3.2–3.4 节的 L-BFGS、GCP 与子空间法提供框架。

### 公式 030

- **出现位置**：`explain.md` 第 251-253 行；章节：[### 3.1. Algorithm overview]；编号标注：(14)、(8)、(9)。
- **公式完整表达**：
[
\left[\pmb {u}_A,d_A\right]^{\mathrm{T}} = \mathbf{C}\left[\pmb {u}_A,d_A\right]^{\mathrm{T}} + \pmb {k}. \quad (14)
]
- **公式含义**：该式处于求解器总览层：描述约束优化模型、投影算子、线搜索、收敛判据与受限自由度处理，是单调收敛与可行性保持的核心。
- **必要推导过程**：
- 在第 [ k ] 次迭代，用二次模型 [ m_k ] 近似原函数 [ \Pi ]。
- 先做投影得到可行下降方向（广义 Cauchy 点），再在自由子空间内二次优化。
- 用 Wolfe/Armijo 类线搜索确定 [ \alpha_k ]，并以投影残量与步长范数作为停止准则。
- **与其他公式的内在联系**：
- 与公式(15)–(19)形成“二次模型 + 投影 + 活动集 + 线搜索 + 停止判据”闭环。
- 为第 3.2–3.4 节的 L-BFGS、GCP 与子空间法提供框架。

### 公式 031

- **出现位置**：`explain.md` 第 259-261 行；章节：[### 3.1. Algorithm overview]；编号标注：(12)、(15)。
- **公式完整表达**：
[
m_{k}(\pmb{x}) = \Pi_{k} + \pmb{r}_{k}^{\mathrm{T}}(\pmb{x} - \pmb{x}_{k}) + \frac{1}{2} (\pmb{x} - \pmb{x}_{k})^{\mathrm{T}}\mathbf{B}_{k}(\pmb{x} - \pmb{x}_{k}) \quad (15)
]
- **公式含义**：该式处于求解器总览层：描述约束优化模型、投影算子、线搜索、收敛判据与受限自由度处理，是单调收敛与可行性保持的核心。
- **必要推导过程**：
- 在第 [ k ] 次迭代，用二次模型 [ m_k ] 近似原函数 [ \Pi ]。
- 先做投影得到可行下降方向（广义 Cauchy 点），再在自由子空间内二次优化。
- 用 Wolfe/Armijo 类线搜索确定 [ \alpha_k ]，并以投影残量与步长范数作为停止准则。
- **与其他公式的内在联系**：
- 与公式(15)–(19)形成“二次模型 + 投影 + 活动集 + 线搜索 + 停止判据”闭环。
- 为第 3.2–3.4 节的 L-BFGS、GCP 与子空间法提供框架。

### 公式 032

- **出现位置**：`explain.md` 第 267-269 行；章节：[### 3.1. Algorithm overview]；编号标注：(15)、(13)、(16)。
- **公式完整表达**：
[
\mathrm{Proj}_c(x_i,\mathrm{lb}_i,\mathrm{ub}_i) = \left\{ \begin{array}{ll}\mathrm{lb}_i & \mathrm{if}\; x_i < \mathrm{lb}_i,\\ x_i & \mathrm{if}\; \mathrm{lb}_i \leq x_i \leq \mathrm{ub}_i,\\ \mathrm{ub}_i & \mathrm{if}\; x_i > \mathrm{ub}_i, \end{array} \right. \quad (16)
]
- **公式含义**：该式处于求解器总览层：描述约束优化模型、投影算子、线搜索、收敛判据与受限自由度处理，是单调收敛与可行性保持的核心。
- **必要推导过程**：
- 在第 [ k ] 次迭代，用二次模型 [ m_k ] 近似原函数 [ \Pi ]。
- 先做投影得到可行下降方向（广义 Cauchy 点），再在自由子空间内二次优化。
- 用 Wolfe/Armijo 类线搜索确定 [ \alpha_k ]，并以投影残量与步长范数作为停止准则。
- **与其他公式的内在联系**：
- 与公式(15)–(19)形成“二次模型 + 投影 + 活动集 + 线搜索 + 停止判据”闭环。
- 为第 3.2–3.4 节的 L-BFGS、GCP 与子空间法提供框架。

### 公式 033

- **出现位置**：`explain.md` 第 273-275 行；章节：[### 3.1. Algorithm overview]；编号标注：(17)。
- **公式完整表达**：
[
\pmb {x}(t) = \mathrm{Proj}_c(\pmb {x}_k - t\pmb {r}_k,\mathbf{lb},\mathbf{ub}),\quad t > 0, \quad (17)
]
- **公式含义**：该式处于求解器总览层：描述约束优化模型、投影算子、线搜索、收敛判据与受限自由度处理，是单调收敛与可行性保持的核心。
- **必要推导过程**：
- 在第 [ k ] 次迭代，用二次模型 [ m_k ] 近似原函数 [ \Pi ]。
- 先做投影得到可行下降方向（广义 Cauchy 点），再在自由子空间内二次优化。
- 用 Wolfe/Armijo 类线搜索确定 [ \alpha_k ]，并以投影残量与步长范数作为停止准则。
- **与其他公式的内在联系**：
- 与公式(15)–(19)形成“二次模型 + 投影 + 活动集 + 线搜索 + 停止判据”闭环。
- 为第 3.2–3.4 节的 L-BFGS、GCP 与子空间法提供框架。

### 公式 034

- **出现位置**：`explain.md` 第 281-283 行；章节：[### 3.1. Algorithm overview]；编号标注：(15)、(18)。
- **公式完整表达**：
[
p_k(t) = m_k(\pmb{x}(t)) = \Pi_k + \pmb{r}_k^\mathrm{T}(\pmb{x}(t) - \pmb{x}_k) + \frac{1}{2} (\pmb{x}(t) - \pmb{x}_k)^\mathrm{T}\mathbf{B}_k(\pmb{x}(t) - \pmb{x}_k). \quad (18)
]
- **公式含义**：该式处于求解器总览层：描述约束优化模型、投影算子、线搜索、收敛判据与受限自由度处理，是单调收敛与可行性保持的核心。
- **必要推导过程**：
- 在第 [ k ] 次迭代，用二次模型 [ m_k ] 近似原函数 [ \Pi ]。
- 先做投影得到可行下降方向（广义 Cauchy 点），再在自由子空间内二次优化。
- 用 Wolfe/Armijo 类线搜索确定 [ \alpha_k ]，并以投影残量与步长范数作为停止准则。
- **与其他公式的内在联系**：
- 与公式(15)–(19)形成“二次模型 + 投影 + 活动集 + 线搜索 + 停止判据”闭环。
- 为第 3.2–3.4 节的 L-BFGS、GCP 与子空间法提供框架。

### 公式 035

- **出现位置**：`explain.md` 第 289-291 行；章节：[### 3.1. Algorithm overview]；编号标注：(13)、(15)。
- **公式完整表达**：
[
\mathcal{A}(\pmb {x}) = \{i:x_i = \mathrm{lb}_i\} \cup \{i:x_i = \mathrm{ub}_i\}.
]
- **公式含义**：该式处于求解器总览层：描述约束优化模型、投影算子、线搜索、收敛判据与受限自由度处理，是单调收敛与可行性保持的核心。
- **必要推导过程**：
- 在第 [ k ] 次迭代，用二次模型 [ m_k ] 近似原函数 [ \Pi ]。
- 先做投影得到可行下降方向（广义 Cauchy 点），再在自由子空间内二次优化。
- 用 Wolfe/Armijo 类线搜索确定 [ \alpha_k ]，并以投影残量与步长范数作为停止准则。
- **与其他公式的内在联系**：
- 与公式(15)–(19)形成“二次模型 + 投影 + 活动集 + 线搜索 + 停止判据”闭环。
- 为第 3.2–3.4 节的 L-BFGS、GCP 与子空间法提供框架。

### 公式 036

- **出现位置**：`explain.md` 第 295-297 行；章节：[### 3.1. Algorithm overview]；编号标注：(15)。
- **公式完整表达**：
[
\pmb{x}^* = \arg \min m_k(\pmb {x})
]
- **公式含义**：该式处于求解器总览层：描述约束优化模型、投影算子、线搜索、收敛判据与受限自由度处理，是单调收敛与可行性保持的核心。
- **必要推导过程**：
- 在第 [ k ] 次迭代，用二次模型 [ m_k ] 近似原函数 [ \Pi ]。
- 先做投影得到可行下降方向（广义 Cauchy 点），再在自由子空间内二次优化。
- 用 Wolfe/Armijo 类线搜索确定 [ \alpha_k ]，并以投影残量与步长范数作为停止准则。
- **与其他公式的内在联系**：
- 与公式(15)–(19)形成“二次模型 + 投影 + 活动集 + 线搜索 + 停止判据”闭环。
- 为第 3.2–3.4 节的 L-BFGS、GCP 与子空间法提供框架。

### 公式 037

- **出现位置**：`explain.md` 第 301-303 行；章节：[### 3.1. Algorithm overview]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
x_{i} = x_{i}^{c},\;\forall i\in \mathcal{A}_{k}(\pmb{x}^{c})\quad \text{and}\quad \mathrm{lb}_{i}\leq x_{i}\leq \mathrm{ub}_{i},\;\forall i\notin \mathcal{A}_{k}(\pmb{x}^{c}).
]
- **公式含义**：该式处于求解器总览层：描述约束优化模型、投影算子、线搜索、收敛判据与受限自由度处理，是单调收敛与可行性保持的核心。
- **必要推导过程**：
- 在第 [ k ] 次迭代，用二次模型 [ m_k ] 近似原函数 [ \Pi ]。
- 先做投影得到可行下降方向（广义 Cauchy 点），再在自由子空间内二次优化。
- 用 Wolfe/Armijo 类线搜索确定 [ \alpha_k ]，并以投影残量与步长范数作为停止准则。
- **与其他公式的内在联系**：
- 与公式(15)–(19)形成“二次模型 + 投影 + 活动集 + 线搜索 + 停止判据”闭环。
- 为第 3.2–3.4 节的 L-BFGS、GCP 与子空间法提供框架。

### 公式 038

- **出现位置**：`explain.md` 第 309-311 行；章节：[### 3.1. Algorithm overview]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\pmb{p}_k = \pmb{x}^* - \pmb{x}_k.
]
- **公式含义**：该式处于求解器总览层：描述约束优化模型、投影算子、线搜索、收敛判据与受限自由度处理，是单调收敛与可行性保持的核心。
- **必要推导过程**：
- 在第 [ k ] 次迭代，用二次模型 [ m_k ] 近似原函数 [ \Pi ]。
- 先做投影得到可行下降方向（广义 Cauchy 点），再在自由子空间内二次优化。
- 用 Wolfe/Armijo 类线搜索确定 [ \alpha_k ]，并以投影残量与步长范数作为停止准则。
- **与其他公式的内在联系**：
- 与公式(15)–(19)形成“二次模型 + 投影 + 活动集 + 线搜索 + 停止判据”闭环。
- 为第 3.2–3.4 节的 L-BFGS、GCP 与子空间法提供框架。

### 公式 039

- **出现位置**：`explain.md` 第 315-317 行；章节：[### 3.1. Algorithm overview]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\pmb{x}_{k + 1} = \pmb{x}_k + \alpha_k\pmb{p}_k.
]
- **公式含义**：该式处于求解器总览层：描述约束优化模型、投影算子、线搜索、收敛判据与受限自由度处理，是单调收敛与可行性保持的核心。
- **必要推导过程**：
- 在第 [ k ] 次迭代，用二次模型 [ m_k ] 近似原函数 [ \Pi ]。
- 先做投影得到可行下降方向（广义 Cauchy 点），再在自由子空间内二次优化。
- 用 Wolfe/Armijo 类线搜索确定 [ \alpha_k ]，并以投影残量与步长范数作为停止准则。
- **与其他公式的内在联系**：
- 与公式(15)–(19)形成“二次模型 + 投影 + 活动集 + 线搜索 + 停止判据”闭环。
- 为第 3.2–3.4 节的 L-BFGS、GCP 与子空间法提供框架。

### 公式 040

- **出现位置**：`explain.md` 第 321-323 行；章节：[### 3.1. Algorithm overview]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\Pi_{k + 1} = \Pi (\pmb{x}_k + \alpha_k\pmb {p}_k) \leq \Pi_k + c_1\alpha_k\pmb {p}_k^\mathrm{T}\pmb{r}_k
]
- **公式含义**：该式处于求解器总览层：描述约束优化模型、投影算子、线搜索、收敛判据与受限自由度处理，是单调收敛与可行性保持的核心。
- **必要推导过程**：
- 在第 [ k ] 次迭代，用二次模型 [ m_k ] 近似原函数 [ \Pi ]。
- 先做投影得到可行下降方向（广义 Cauchy 点），再在自由子空间内二次优化。
- 用 Wolfe/Armijo 类线搜索确定 [ \alpha_k ]，并以投影残量与步长范数作为停止准则。
- **与其他公式的内在联系**：
- 与公式(15)–(19)形成“二次模型 + 投影 + 活动集 + 线搜索 + 停止判据”闭环。
- 为第 3.2–3.4 节的 L-BFGS、GCP 与子空间法提供框架。

### 公式 041

- **出现位置**：`explain.md` 第 327-329 行；章节：[### 3.1. Algorithm overview]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\left|\pmb{r}_{k + 1}^\mathrm{T}\pmb {p}_k\right| = \left|\pmb {r}(\pmb{x}_k + \alpha_k\pmb {p}_k)^\mathrm{T}\pmb {p}_k\right| \leq c_2\left|\pmb{r}_k^\mathrm{T}\pmb {p}_k\right|.
]
- **公式含义**：该式处于求解器总览层：描述约束优化模型、投影算子、线搜索、收敛判据与受限自由度处理，是单调收敛与可行性保持的核心。
- **必要推导过程**：
- 在第 [ k ] 次迭代，用二次模型 [ m_k ] 近似原函数 [ \Pi ]。
- 先做投影得到可行下降方向（广义 Cauchy 点），再在自由子空间内二次优化。
- 用 Wolfe/Armijo 类线搜索确定 [ \alpha_k ]，并以投影残量与步长范数作为停止准则。
- **与其他公式的内在联系**：
- 与公式(15)–(19)形成“二次模型 + 投影 + 活动集 + 线搜索 + 停止判据”闭环。
- 为第 3.2–3.4 节的 L-BFGS、GCP 与子空间法提供框架。

### 公式 042

- **出现位置**：`explain.md` 第 336-338 行；章节：[### 3.1. Algorithm overview]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\mathcal{A}_{k + 1}(\pmb{x}^c) = \mathcal{A}_k(\pmb{x}^c).
]
- **公式含义**：该式处于求解器总览层：描述约束优化模型、投影算子、线搜索、收敛判据与受限自由度处理，是单调收敛与可行性保持的核心。
- **必要推导过程**：
- 在第 [ k ] 次迭代，用二次模型 [ m_k ] 近似原函数 [ \Pi ]。
- 先做投影得到可行下降方向（广义 Cauchy 点），再在自由子空间内二次优化。
- 用 Wolfe/Armijo 类线搜索确定 [ \alpha_k ]，并以投影残量与步长范数作为停止准则。
- **与其他公式的内在联系**：
- 与公式(15)–(19)形成“二次模型 + 投影 + 活动集 + 线搜索 + 停止判据”闭环。
- 为第 3.2–3.4 节的 L-BFGS、GCP 与子空间法提供框架。

### 公式 043

- **出现位置**：`explain.md` 第 340-342 行；章节：[### 3.1. Algorithm overview]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\| \mathrm{Proj}_C(\mathbf{x}_{k + 1} - \mathbf{r}_{k + 1},\mathbf{lb},\mathbf{ub}) - \mathbf{x}_{k + 1}\|_{2} < \mathrm{tol}.
]
- **公式含义**：该式处于求解器总览层：描述约束优化模型、投影算子、线搜索、收敛判据与受限自由度处理，是单调收敛与可行性保持的核心。
- **必要推导过程**：
- 在第 [ k ] 次迭代，用二次模型 [ m_k ] 近似原函数 [ \Pi ]。
- 先做投影得到可行下降方向（广义 Cauchy 点），再在自由子空间内二次优化。
- 用 Wolfe/Armijo 类线搜索确定 [ \alpha_k ]，并以投影残量与步长范数作为停止准则。
- **与其他公式的内在联系**：
- 与公式(15)–(19)形成“二次模型 + 投影 + 活动集 + 线搜索 + 停止判据”闭环。
- 为第 3.2–3.4 节的 L-BFGS、GCP 与子空间法提供框架。

### 公式 044

- **出现位置**：`explain.md` 第 344-346 行；章节：[### 3.1. Algorithm overview]；编号标注：(16)。
- **公式完整表达**：
[
\| \Delta \mathbf{x}_{k + 1}\|_{2} = \| \mathbf{x}_{k + 1} - \mathbf{x}_{k}\|_{2} < \mathrm{tol}.
]
- **公式含义**：该式处于求解器总览层：描述约束优化模型、投影算子、线搜索、收敛判据与受限自由度处理，是单调收敛与可行性保持的核心。
- **必要推导过程**：
- 在第 [ k ] 次迭代，用二次模型 [ m_k ] 近似原函数 [ \Pi ]。
- 先做投影得到可行下降方向（广义 Cauchy 点），再在自由子空间内二次优化。
- 用 Wolfe/Armijo 类线搜索确定 [ \alpha_k ]，并以投影残量与步长范数作为停止准则。
- **与其他公式的内在联系**：
- 与公式(15)–(19)形成“二次模型 + 投影 + 活动集 + 线搜索 + 停止判据”闭环。
- 为第 3.2–3.4 节的 L-BFGS、GCP 与子空间法提供框架。

### 公式 045

- **出现位置**：`explain.md` 第 350-352 行；章节：[### 3.1. Algorithm overview]；编号标注：(16)、(19)。
- **公式完整表达**：
[
\mathrm{Proj}_C(x_i - r_i,\mathrm{lb}_i,\mathrm{ub}_i) - x_i = \left\{ \begin{array}{ll}\mathrm{lb}_i - x_i & \mathrm{if}\; x_i - r_i < \mathrm{lb}_i,\\ - r_i & \mathrm{if}\; \mathrm{lb}_i \leq x_i - r_i \leq \mathrm{ub}_i,\\ \mathrm{ub}_i - x_i & \mathrm{if}\; x_i - r_i > \mathrm{ub}_i. \end{array} \right. \quad (19)
]
- **公式含义**：该式处于求解器总览层：描述约束优化模型、投影算子、线搜索、收敛判据与受限自由度处理，是单调收敛与可行性保持的核心。
- **必要推导过程**：
- 在第 [ k ] 次迭代，用二次模型 [ m_k ] 近似原函数 [ \Pi ]。
- 先做投影得到可行下降方向（广义 Cauchy 点），再在自由子空间内二次优化。
- 用 Wolfe/Armijo 类线搜索确定 [ \alpha_k ]，并以投影残量与步长范数作为停止准则。
- **与其他公式的内在联系**：
- 与公式(15)–(19)形成“二次模型 + 投影 + 活动集 + 线搜索 + 停止判据”闭环。
- 为第 3.2–3.4 节的 L-BFGS、GCP 与子空间法提供框架。

### 公式 046

- **出现位置**：`explain.md` 第 358-360 行；章节：[### 3.1. Algorithm overview]；编号标注：(14)。
- **公式完整表达**：
[
\mathbf{x} = \mathbf{C}\mathbf{x} + \mathbf{k}.
]
- **公式含义**：该式处于求解器总览层：描述约束优化模型、投影算子、线搜索、收敛判据与受限自由度处理，是单调收敛与可行性保持的核心。
- **必要推导过程**：
- 在第 [ k ] 次迭代，用二次模型 [ m_k ] 近似原函数 [ \Pi ]。
- 先做投影得到可行下降方向（广义 Cauchy 点），再在自由子空间内二次优化。
- 用 Wolfe/Armijo 类线搜索确定 [ \alpha_k ]，并以投影残量与步长范数作为停止准则。
- **与其他公式的内在联系**：
- 与公式(15)–(19)形成“二次模型 + 投影 + 活动集 + 线搜索 + 停止判据”闭环。
- 为第 3.2–3.4 节的 L-BFGS、GCP 与子空间法提供框架。

### 公式 047

- **出现位置**：`explain.md` 第 364-366 行；章节：[### 3.1. Algorithm overview]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\left(\mathbf{C}^{\mathrm{T}}\mathbf{A}\mathbf{C} + \mathbf{I}_{d_{c}}\right)\hat{\mathbf{x}} = \mathbf{C}^{\mathrm{T}}(\mathbf{b} - \mathbf{A}\mathbf{k})
]
- **公式含义**：该式处于求解器总览层：描述约束优化模型、投影算子、线搜索、收敛判据与受限自由度处理，是单调收敛与可行性保持的核心。
- **必要推导过程**：
- 在第 [ k ] 次迭代，用二次模型 [ m_k ] 近似原函数 [ \Pi ]。
- 先做投影得到可行下降方向（广义 Cauchy 点），再在自由子空间内二次优化。
- 用 Wolfe/Armijo 类线搜索确定 [ \alpha_k ]，并以投影残量与步长范数作为停止准则。
- **与其他公式的内在联系**：
- 与公式(15)–(19)形成“二次模型 + 投影 + 活动集 + 线搜索 + 停止判据”闭环。
- 为第 3.2–3.4 节的 L-BFGS、GCP 与子空间法提供框架。

### 公式 048

- **出现位置**：`explain.md` 第 370-372 行；章节：[### 3.1. Algorithm overview]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\mathbf{x} = \mathbf{C}\hat{\mathbf{x}} + \mathbf{k}.
]
- **公式含义**：该式处于求解器总览层：描述约束优化模型、投影算子、线搜索、收敛判据与受限自由度处理，是单调收敛与可行性保持的核心。
- **必要推导过程**：
- 在第 [ k ] 次迭代，用二次模型 [ m_k ] 近似原函数 [ \Pi ]。
- 先做投影得到可行下降方向（广义 Cauchy 点），再在自由子空间内二次优化。
- 用 Wolfe/Armijo 类线搜索确定 [ \alpha_k ]，并以投影残量与步长范数作为停止准则。
- **与其他公式的内在联系**：
- 与公式(15)–(19)形成“二次模型 + 投影 + 活动集 + 线搜索 + 停止判据”闭环。
- 为第 3.2–3.4 节的 L-BFGS、GCP 与子空间法提供框架。

### 公式 049

- **出现位置**：`explain.md` 第 375-377 行；章节：[### 3.1. Algorithm overview]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
(\mathbf{I}_{d_c})_{ii} = \begin{cases} 1 & \text{if DoF } i \text{ is constrained}, \\ 0 & \text{otherwise}, \end{cases}
]
- **公式含义**：该式处于求解器总览层：描述约束优化模型、投影算子、线搜索、收敛判据与受限自由度处理，是单调收敛与可行性保持的核心。
- **必要推导过程**：
- 在第 [ k ] 次迭代，用二次模型 [ m_k ] 近似原函数 [ \Pi ]。
- 先做投影得到可行下降方向（广义 Cauchy 点），再在自由子空间内二次优化。
- 用 Wolfe/Armijo 类线搜索确定 [ \alpha_k ]，并以投影残量与步长范数作为停止准则。
- **与其他公式的内在联系**：
- 与公式(15)–(19)形成“二次模型 + 投影 + 活动集 + 线搜索 + 停止判据”闭环。
- 为第 3.2–3.4 节的 L-BFGS、GCP 与子空间法提供框架。

### 公式 050

- **出现位置**：`explain.md` 第 384-386 行；章节：[### 3.2. Compact representation of limited-memory BFGS matrix]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\mathbf{s}_k = \mathbf{x}_{k + 1} - \mathbf{x}_k,\quad \mathbf{y}_k = \mathbf{r}_{k + 1} - \mathbf{r}_k,
]
- **公式含义**：该式处于 L-BFGS 紧凑表示层：定义 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对、矩阵 [ \mathbf{W}_k,\mathbf{M}_k ]，用于低秩更新近似 Hessian。
- **必要推导过程**：
- 由 BFGS 割线条件 [ \mathbf{B}_{k+1}\mathbf{s}_k=\mathbf{y}_k ] 与最小改变量原则得到标准更新式。
- 限内存版本仅保留最近 [ m ] 组 [ (\mathbf{s},\mathbf{y}) ]，写成紧凑低秩形式。
- 紧凑式便于与附录中的低秩逆公式联用，降低每步复杂度。
- **与其他公式的内在联系**：
- 与公式(20)(21)一起把 [ \mathbf{B}_k ] 表示为“块对角初值 + 低秩修正”。
- 直接服务于第 3.3/3.4 节中 [ f'_j,f''_j ]、子空间线性系统与 Schur 补。

### 公式 051

- **出现位置**：`explain.md` 第 390-392 行；章节：[### 3.2. Compact representation of limited-memory BFGS matrix]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\mathbf{B}_{k + 1} = \mathbf{B}_k - \frac{\mathbf{B}_k\mathbf{s}_k\mathbf{s}_k^{\mathrm{T}}\mathbf{B}_k}{\mathbf{s}_k^{\mathrm{T}}\mathbf{B}_k\mathbf{s}_k} + \frac{\mathbf{y}_k\mathbf{y}_k^{\mathrm{T}}}{\mathbf{y}_k^{\mathrm{T}}\mathbf{s}_k}.
]
- **公式含义**：该式处于 L-BFGS 紧凑表示层：定义 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对、矩阵 [ \mathbf{W}_k,\mathbf{M}_k ]，用于低秩更新近似 Hessian。
- **必要推导过程**：
- 由 BFGS 割线条件 [ \mathbf{B}_{k+1}\mathbf{s}_k=\mathbf{y}_k ] 与最小改变量原则得到标准更新式。
- 限内存版本仅保留最近 [ m ] 组 [ (\mathbf{s},\mathbf{y}) ]，写成紧凑低秩形式。
- 紧凑式便于与附录中的低秩逆公式联用，降低每步复杂度。
- **与其他公式的内在联系**：
- 与公式(20)(21)一起把 [ \mathbf{B}_k ] 表示为“块对角初值 + 低秩修正”。
- 直接服务于第 3.3/3.4 节中 [ f'_j,f''_j ]、子空间线性系统与 Schur 补。

### 公式 052

- **出现位置**：`explain.md` 第 396-398 行；章节：[### 3.2. Compact representation of limited-memory BFGS matrix]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\mathbf{s}_k^{\mathrm{T}}\mathbf{y}_k > 0,
]
- **公式含义**：该式处于 L-BFGS 紧凑表示层：定义 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对、矩阵 [ \mathbf{W}_k,\mathbf{M}_k ]，用于低秩更新近似 Hessian。
- **必要推导过程**：
- 由 BFGS 割线条件 [ \mathbf{B}_{k+1}\mathbf{s}_k=\mathbf{y}_k ] 与最小改变量原则得到标准更新式。
- 限内存版本仅保留最近 [ m ] 组 [ (\mathbf{s},\mathbf{y}) ]，写成紧凑低秩形式。
- 紧凑式便于与附录中的低秩逆公式联用，降低每步复杂度。
- **与其他公式的内在联系**：
- 与公式(20)(21)一起把 [ \mathbf{B}_k ] 表示为“块对角初值 + 低秩修正”。
- 直接服务于第 3.3/3.4 节中 [ f'_j,f''_j ]、子空间线性系统与 Schur 补。

### 公式 053

- **出现位置**：`explain.md` 第 404-406 行；章节：[### 3.2. Compact representation of limited-memory BFGS matrix]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\mathbf{S}_k = [\mathbf{s}_{k-m} \cdots \mathbf{s}_{k-1}]
]
- **公式含义**：该式处于 L-BFGS 紧凑表示层：定义 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对、矩阵 [ \mathbf{W}_k,\mathbf{M}_k ]，用于低秩更新近似 Hessian。
- **必要推导过程**：
- 由 BFGS 割线条件 [ \mathbf{B}_{k+1}\mathbf{s}_k=\mathbf{y}_k ] 与最小改变量原则得到标准更新式。
- 限内存版本仅保留最近 [ m ] 组 [ (\mathbf{s},\mathbf{y}) ]，写成紧凑低秩形式。
- 紧凑式便于与附录中的低秩逆公式联用，降低每步复杂度。
- **与其他公式的内在联系**：
- 与公式(20)(21)一起把 [ \mathbf{B}_k ] 表示为“块对角初值 + 低秩修正”。
- 直接服务于第 3.3/3.4 节中 [ f'_j,f''_j ]、子空间线性系统与 Schur 补。

### 公式 054

- **出现位置**：`explain.md` 第 408-410 行；章节：[### 3.2. Compact representation of limited-memory BFGS matrix]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\mathbf{Y}_k = [\mathbf{y}_{k-m} \cdots \mathbf{y}_{k-1}].
]
- **公式含义**：该式处于 L-BFGS 紧凑表示层：定义 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对、矩阵 [ \mathbf{W}_k,\mathbf{M}_k ]，用于低秩更新近似 Hessian。
- **必要推导过程**：
- 由 BFGS 割线条件 [ \mathbf{B}_{k+1}\mathbf{s}_k=\mathbf{y}_k ] 与最小改变量原则得到标准更新式。
- 限内存版本仅保留最近 [ m ] 组 [ (\mathbf{s},\mathbf{y}) ]，写成紧凑低秩形式。
- 紧凑式便于与附录中的低秩逆公式联用，降低每步复杂度。
- **与其他公式的内在联系**：
- 与公式(20)(21)一起把 [ \mathbf{B}_k ] 表示为“块对角初值 + 低秩修正”。
- 直接服务于第 3.3/3.4 节中 [ f'_j,f''_j ]、子空间线性系统与 Schur 补。

### 公式 055

- **出现位置**：`explain.md` 第 414-416 行；章节：[### 3.2. Compact representation of limited-memory BFGS matrix]；编号标注：(20)。
- **公式完整表达**：
[
\mathbf{B}_k = \mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k, \quad (20)
]
- **公式含义**：该式处于 L-BFGS 紧凑表示层：定义 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对、矩阵 [ \mathbf{W}_k,\mathbf{M}_k ]，用于低秩更新近似 Hessian。
- **必要推导过程**：
- 由 BFGS 割线条件 [ \mathbf{B}_{k+1}\mathbf{s}_k=\mathbf{y}_k ] 与最小改变量原则得到标准更新式。
- 限内存版本仅保留最近 [ m ] 组 [ (\mathbf{s},\mathbf{y}) ]，写成紧凑低秩形式。
- 紧凑式便于与附录中的低秩逆公式联用，降低每步复杂度。
- **与其他公式的内在联系**：
- 与公式(20)(21)一起把 [ \mathbf{B}_k ] 表示为“块对角初值 + 低秩修正”。
- 直接服务于第 3.3/3.4 节中 [ f'_j,f''_j ]、子空间线性系统与 Schur 补。

### 公式 056

- **出现位置**：`explain.md` 第 420-422 行；章节：[### 3.2. Compact representation of limited-memory BFGS matrix]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\mathbf{W}_k = [\mathbf{Y}_k \quad \mathbf{B}^0_k \mathbf{S}_k] \in \mathbb{R}^{n\times 2m}
]
- **公式含义**：该式处于 L-BFGS 紧凑表示层：定义 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对、矩阵 [ \mathbf{W}_k,\mathbf{M}_k ]，用于低秩更新近似 Hessian。
- **必要推导过程**：
- 由 BFGS 割线条件 [ \mathbf{B}_{k+1}\mathbf{s}_k=\mathbf{y}_k ] 与最小改变量原则得到标准更新式。
- 限内存版本仅保留最近 [ m ] 组 [ (\mathbf{s},\mathbf{y}) ]，写成紧凑低秩形式。
- 紧凑式便于与附录中的低秩逆公式联用，降低每步复杂度。
- **与其他公式的内在联系**：
- 与公式(20)(21)一起把 [ \mathbf{B}_k ] 表示为“块对角初值 + 低秩修正”。
- 直接服务于第 3.3/3.4 节中 [ f'_j,f''_j ]、子空间线性系统与 Schur 补。

### 公式 057

- **出现位置**：`explain.md` 第 425-427 行；章节：[### 3.2. Compact representation of limited-memory BFGS matrix]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\mathbf{M}_k = \begin{bmatrix} -\mathbf{D}_k & \mathbf{L}^{\mathrm{T}}_k \\ \mathbf{L}_k & \mathbf{S}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{S}_k \end{bmatrix}^{-1} \in \mathbb{R}^{2m\times 2m}.
]
- **公式含义**：该式处于 L-BFGS 紧凑表示层：定义 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对、矩阵 [ \mathbf{W}_k,\mathbf{M}_k ]，用于低秩更新近似 Hessian。
- **必要推导过程**：
- 由 BFGS 割线条件 [ \mathbf{B}_{k+1}\mathbf{s}_k=\mathbf{y}_k ] 与最小改变量原则得到标准更新式。
- 限内存版本仅保留最近 [ m ] 组 [ (\mathbf{s},\mathbf{y}) ]，写成紧凑低秩形式。
- 紧凑式便于与附录中的低秩逆公式联用，降低每步复杂度。
- **与其他公式的内在联系**：
- 与公式(20)(21)一起把 [ \mathbf{B}_k ] 表示为“块对角初值 + 低秩修正”。
- 直接服务于第 3.3/3.4 节中 [ f'_j,f''_j ]、子空间线性系统与 Schur 补。

### 公式 058

- **出现位置**：`explain.md` 第 431-433 行；章节：[### 3.2. Compact representation of limited-memory BFGS matrix]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\mathbf{D}_k = \mathrm{diag}\{\mathbf{s}^{\mathrm{T}}_{k-m}\mathbf{y}_{k-m}, \ldots, \mathbf{s}^{\mathrm{T}}_{k-1}\mathbf{y}_{k-1}\} \in \mathbb{R}^{m\times m}
]
- **公式含义**：该式处于 L-BFGS 紧凑表示层：定义 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对、矩阵 [ \mathbf{W}_k,\mathbf{M}_k ]，用于低秩更新近似 Hessian。
- **必要推导过程**：
- 由 BFGS 割线条件 [ \mathbf{B}_{k+1}\mathbf{s}_k=\mathbf{y}_k ] 与最小改变量原则得到标准更新式。
- 限内存版本仅保留最近 [ m ] 组 [ (\mathbf{s},\mathbf{y}) ]，写成紧凑低秩形式。
- 紧凑式便于与附录中的低秩逆公式联用，降低每步复杂度。
- **与其他公式的内在联系**：
- 与公式(20)(21)一起把 [ \mathbf{B}_k ] 表示为“块对角初值 + 低秩修正”。
- 直接服务于第 3.3/3.4 节中 [ f'_j,f''_j ]、子空间线性系统与 Schur 补。

### 公式 059

- **出现位置**：`explain.md` 第 435-437 行；章节：[### 3.2. Compact representation of limited-memory BFGS matrix]；编号标注：(11)。
- **公式完整表达**：
[
\mathbf{L}_k \in \mathbb{R}^{m\times m}, \quad (\mathbf{L}_k)_{i,j} = \begin{cases} \mathbf{s}^{\mathrm{T}}_{i+k-m-1}\mathbf{y}_{j+k-m-1} & \text{if } i > j, \\ 0 & \text{if } i \leq j. \end{cases}
]
- **公式含义**：该式处于 L-BFGS 紧凑表示层：定义 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对、矩阵 [ \mathbf{W}_k,\mathbf{M}_k ]，用于低秩更新近似 Hessian。
- **必要推导过程**：
- 由 BFGS 割线条件 [ \mathbf{B}_{k+1}\mathbf{s}_k=\mathbf{y}_k ] 与最小改变量原则得到标准更新式。
- 限内存版本仅保留最近 [ m ] 组 [ (\mathbf{s},\mathbf{y}) ]，写成紧凑低秩形式。
- 紧凑式便于与附录中的低秩逆公式联用，降低每步复杂度。
- **与其他公式的内在联系**：
- 与公式(20)(21)一起把 [ \mathbf{B}_k ] 表示为“块对角初值 + 低秩修正”。
- 直接服务于第 3.3/3.4 节中 [ f'_j,f''_j ]、子空间线性系统与 Schur 补。

### 公式 060

- **出现位置**：`explain.md` 第 441-443 行；章节：[### 3.2. Compact representation of limited-memory BFGS matrix]；编号标注：(11)、(21)、(20)。
- **公式完整表达**：
[
\mathbf{B}^0_k = \hat{\mathbf{K}}^{(k)} = \begin{bmatrix} \mathbf{K}_{uu} & \mathbf{0} \\ \mathbf{0} & \mathbf{K}_{dd} \end{bmatrix}^{(k)}. \quad (21)
]
- **公式含义**：该式处于 L-BFGS 紧凑表示层：定义 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对、矩阵 [ \mathbf{W}_k,\mathbf{M}_k ]，用于低秩更新近似 Hessian。
- **必要推导过程**：
- 由 BFGS 割线条件 [ \mathbf{B}_{k+1}\mathbf{s}_k=\mathbf{y}_k ] 与最小改变量原则得到标准更新式。
- 限内存版本仅保留最近 [ m ] 组 [ (\mathbf{s},\mathbf{y}) ]，写成紧凑低秩形式。
- 紧凑式便于与附录中的低秩逆公式联用，降低每步复杂度。
- **与其他公式的内在联系**：
- 与公式(20)(21)一起把 [ \mathbf{B}_k ] 表示为“块对角初值 + 低秩修正”。
- 直接服务于第 3.3/3.4 节中 [ f'_j,f''_j ]、子空间线性系统与 Schur 补。

### 公式 061

- **出现位置**：`explain.md` 第 455-457 行；章节：[### 3.3. Generalized Cauchy point]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
t_i = \begin{cases} (x^0_i - \mathrm{ub}_i)/r_i & \text{if } r_i < 0, \\ (x^0_i - \mathrm{lb}_i)/r_i & \text{if } r_i > 0, \\ +\infty & \text{if } r_i = 0. \end{cases}
]
- **公式含义**：该式处于广义 Cauchy 点（GCP）层：沿投影梯度折线搜索，动态识别活动集并更新一维二次模型导数，保证约束可行下降。
- **必要推导过程**：
- 沿投影梯度轨迹分段前进：每遇到一个变量触边界，就更新活动集。
- 在每段上模型值是关于段参数 [ \Delta t ] 的一元二次函数，可解析求极小点。
- 若段内极小点不可行，则走到下一个断点继续，直到得到广义 Cauchy 点。
- **与其他公式的内在联系**：
- 以第 3.2 节近似 Hessian 为曲率信息，构造可行下降的起点 [ \mathbf{x}^c ]。
- 其输出活动集 [ \mathcal{A}(\mathbf{x}^c) ] 直接传入第 3.4 节子空间最小化。

### 公式 062

- **出现位置**：`explain.md` 第 461-463 行；章节：[### 3.3. Generalized Cauchy point]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
x_i(t) = x_i^{0} - \min \{t,t_i\} r_i.
]
- **公式含义**：该式处于广义 Cauchy 点（GCP）层：沿投影梯度折线搜索，动态识别活动集并更新一维二次模型导数，保证约束可行下降。
- **必要推导过程**：
- 沿投影梯度轨迹分段前进：每遇到一个变量触边界，就更新活动集。
- 在每段上模型值是关于段参数 [ \Delta t ] 的一元二次函数，可解析求极小点。
- 若段内极小点不可行，则走到下一个断点继续，直到得到广义 Cauchy 点。
- **与其他公式的内在联系**：
- 以第 3.2 节近似 Hessian 为曲率信息，构造可行下降的起点 [ \mathbf{x}^c ]。
- 其输出活动集 [ \mathcal{A}(\mathbf{x}^c) ] 直接传入第 3.4 节子空间最小化。

### 公式 063

- **出现位置**：`explain.md` 第 467-469 行；章节：[### 3.3. Generalized Cauchy point]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\mathbf{x}(t) = \mathbf{x}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)},
]
- **公式含义**：该式处于广义 Cauchy 点（GCP）层：沿投影梯度折线搜索，动态识别活动集并更新一维二次模型导数，保证约束可行下降。
- **必要推导过程**：
- 沿投影梯度轨迹分段前进：每遇到一个变量触边界，就更新活动集。
- 在每段上模型值是关于段参数 [ \Delta t ] 的一元二次函数，可解析求极小点。
- 若段内极小点不可行，则走到下一个断点继续，直到得到广义 Cauchy 点。
- **与其他公式的内在联系**：
- 以第 3.2 节近似 Hessian 为曲率信息，构造可行下降的起点 [ \mathbf{x}^c ]。
- 其输出活动集 [ \mathcal{A}(\mathbf{x}^c) ] 直接传入第 3.4 节子空间最小化。

### 公式 064

- **出现位置**：`explain.md` 第 473-481 行；章节：[### 3.3. Generalized Cauchy point]；编号标注：(22)。
- **公式完整表达**：
[
\begin{array}{rl} 
p(t) &= m(\mathbf{x}(t)) = \Pi + \mathbf{r}^{\mathrm{T}}(\mathbf{x}(t) - \mathbf{x}^0) + \frac{1}{2} (\mathbf{x}(t) - \mathbf{x}^0)^{\mathrm{T}}\mathbf{B}(\mathbf{x}(t) - \mathbf{x}^0) \\
&= \Pi + \mathbf{r}^{\mathrm{T}}(\mathbf{z}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)}) + \frac{1}{2} (\mathbf{z}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)})^{\mathrm{T}}\mathbf{B}(\mathbf{z}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)}) \\
&= \left(\Pi + \mathbf{r}^{\mathrm{T}}\mathbf{z}^{(j-1)} + \frac{1}{2} \mathbf{z}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j-1)}\right) \\
&\quad + \left(\mathbf{r}^{\mathrm{T}}\mathbf{d}^{(j-1)} + \mathbf{d}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j-1)}\right)\Delta t + \frac{1}{2}\left(\mathbf{d}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j-1)}\right)\Delta t^{2} \\
&= f_{j-1} + f_{j-1}^{\prime}\Delta t + \frac{1}{2} f_{j-1}^{\prime \prime}\Delta t^{2} = \hat{p} (\Delta t), 
\end{array} \quad (22)
]
- **公式含义**：该式处于广义 Cauchy 点（GCP）层：沿投影梯度折线搜索，动态识别活动集并更新一维二次模型导数，保证约束可行下降。
- **必要推导过程**：
- 沿投影梯度轨迹分段前进：每遇到一个变量触边界，就更新活动集。
- 在每段上模型值是关于段参数 [ \Delta t ] 的一元二次函数，可解析求极小点。
- 若段内极小点不可行，则走到下一个断点继续，直到得到广义 Cauchy 点。
- **与其他公式的内在联系**：
- 以第 3.2 节近似 Hessian 为曲率信息，构造可行下降的起点 [ \mathbf{x}^c ]。
- 其输出活动集 [ \mathcal{A}(\mathbf{x}^c) ] 直接传入第 3.4 节子空间最小化。

### 公式 065

- **出现位置**：`explain.md` 第 487-489 行；章节：[### 3.3. Generalized Cauchy point]；编号标注：(18)。
- **公式完整表达**：
[
f_{j-1}^{\prime \prime} = \mathbf{d}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j-1)} > 0.
]
- **公式含义**：该式处于广义 Cauchy 点（GCP）层：沿投影梯度折线搜索，动态识别活动集并更新一维二次模型导数，保证约束可行下降。
- **必要推导过程**：
- 沿投影梯度轨迹分段前进：每遇到一个变量触边界，就更新活动集。
- 在每段上模型值是关于段参数 [ \Delta t ] 的一元二次函数，可解析求极小点。
- 若段内极小点不可行，则走到下一个断点继续，直到得到广义 Cauchy 点。
- **与其他公式的内在联系**：
- 以第 3.2 节近似 Hessian 为曲率信息，构造可行下降的起点 [ \mathbf{x}^c ]。
- 其输出活动集 [ \mathcal{A}(\mathbf{x}^c) ] 直接传入第 3.4 节子空间最小化。

### 公式 066

- **出现位置**：`explain.md` 第 493-495 行；章节：[### 3.3. Generalized Cauchy point]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\hat{p}'(\Delta t) = 0 \Rightarrow \Delta t^* = -f_{j-1}' / f_{j-1}^{\prime \prime}.
]
- **公式含义**：该式处于广义 Cauchy 点（GCP）层：沿投影梯度折线搜索，动态识别活动集并更新一维二次模型导数，保证约束可行下降。
- **必要推导过程**：
- 沿投影梯度轨迹分段前进：每遇到一个变量触边界，就更新活动集。
- 在每段上模型值是关于段参数 [ \Delta t ] 的一元二次函数，可解析求极小点。
- 若段内极小点不可行，则走到下一个断点继续，直到得到广义 Cauchy 点。
- **与其他公式的内在联系**：
- 以第 3.2 节近似 Hessian 为曲率信息，构造可行下降的起点 [ \mathbf{x}^c ]。
- 其输出活动集 [ \mathcal{A}(\mathbf{x}^c) ] 直接传入第 3.4 节子空间最小化。

### 公式 067

- **出现位置**：`explain.md` 第 504-506 行；章节：[### 3.3. Generalized Cauchy point]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\mathcal{A}(\mathbf{x}^c) = \{i : t_i \leq t^{(j-1)}, i = 1, \ldots, n\}.
]
- **公式含义**：该式处于广义 Cauchy 点（GCP）层：沿投影梯度折线搜索，动态识别活动集并更新一维二次模型导数，保证约束可行下降。
- **必要推导过程**：
- 沿投影梯度轨迹分段前进：每遇到一个变量触边界，就更新活动集。
- 在每段上模型值是关于段参数 [ \Delta t ] 的一元二次函数，可解析求极小点。
- 若段内极小点不可行，则走到下一个断点继续，直到得到广义 Cauchy 点。
- **与其他公式的内在联系**：
- 以第 3.2 节近似 Hessian 为曲率信息，构造可行下降的起点 [ \mathbf{x}^c ]。
- 其输出活动集 [ \mathcal{A}(\mathbf{x}^c) ] 直接传入第 3.4 节子空间最小化。

### 公式 068

- **出现位置**：`explain.md` 第 509-511 行；章节：[### 3.3. Generalized Cauchy point]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\Delta t^{(j-1)} = t^{(j)} - t^{(j-1)}, \quad \mathbf{x}^{(j)} = \mathbf{x}^{(j-1)} + \Delta t^{(j-1)}\mathbf{d}^{(j-1)}.
]
- **公式含义**：该式处于广义 Cauchy 点（GCP）层：沿投影梯度折线搜索，动态识别活动集并更新一维二次模型导数，保证约束可行下降。
- **必要推导过程**：
- 沿投影梯度轨迹分段前进：每遇到一个变量触边界，就更新活动集。
- 在每段上模型值是关于段参数 [ \Delta t ] 的一元二次函数，可解析求极小点。
- 若段内极小点不可行，则走到下一个断点继续，直到得到广义 Cauchy 点。
- **与其他公式的内在联系**：
- 以第 3.2 节近似 Hessian 为曲率信息，构造可行下降的起点 [ \mathbf{x}^c ]。
- 其输出活动集 [ \mathcal{A}(\mathbf{x}^c) ] 直接传入第 3.4 节子空间最小化。

### 公式 069

- **出现位置**：`explain.md` 第 514-516 行；章节：[### 3.3. Generalized Cauchy point]；编号标注：(22)。
- **公式完整表达**：
[
\mathbf{d}^{(j)} = \mathbf{d}^{(j-1)} + r_b \mathbf{e}_b,
]
- **公式含义**：该式处于广义 Cauchy 点（GCP）层：沿投影梯度折线搜索，动态识别活动集并更新一维二次模型导数，保证约束可行下降。
- **必要推导过程**：
- 沿投影梯度轨迹分段前进：每遇到一个变量触边界，就更新活动集。
- 在每段上模型值是关于段参数 [ \Delta t ] 的一元二次函数，可解析求极小点。
- 若段内极小点不可行，则走到下一个断点继续，直到得到广义 Cauchy 点。
- **与其他公式的内在联系**：
- 以第 3.2 节近似 Hessian 为曲率信息，构造可行下降的起点 [ \mathbf{x}^c ]。
- 其输出活动集 [ \mathcal{A}(\mathbf{x}^c) ] 直接传入第 3.4 节子空间最小化。

### 公式 070

- **出现位置**：`explain.md` 第 519-521 行；章节：[### 3.3. Generalized Cauchy point]；编号标注：(22)、(23)。
- **公式完整表达**：
[
f'_j = \mathbf{r}^{\mathrm{T}}\mathbf{d}^{(j)} + \mathbf{d}^{(j)}{}^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j)} = f'_{j-1} + \Delta t^{(j-1)} f''_{j-1} + r^2_b + r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j)} \quad (23)
]
- **公式含义**：该式处于广义 Cauchy 点（GCP）层：沿投影梯度折线搜索，动态识别活动集并更新一维二次模型导数，保证约束可行下降。
- **必要推导过程**：
- 沿投影梯度轨迹分段前进：每遇到一个变量触边界，就更新活动集。
- 在每段上模型值是关于段参数 [ \Delta t ] 的一元二次函数，可解析求极小点。
- 若段内极小点不可行，则走到下一个断点继续，直到得到广义 Cauchy 点。
- **与其他公式的内在联系**：
- 以第 3.2 节近似 Hessian 为曲率信息，构造可行下降的起点 [ \mathbf{x}^c ]。
- 其输出活动集 [ \mathcal{A}(\mathbf{x}^c) ] 直接传入第 3.4 节子空间最小化。

### 公式 071

- **出现位置**：`explain.md` 第 523-525 行；章节：[### 3.3. Generalized Cauchy point]；编号标注：(24)、(20)。
- **公式完整表达**：
[
f''_j = \mathbf{d}^{(j)}{}^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j)} = f''_{j-1} + 2r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j-1)} + r^2_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{e}_b. \quad (24)
]
- **公式含义**：该式处于广义 Cauchy 点（GCP）层：沿投影梯度折线搜索，动态识别活动集并更新一维二次模型导数，保证约束可行下降。
- **必要推导过程**：
- 沿投影梯度轨迹分段前进：每遇到一个变量触边界，就更新活动集。
- 在每段上模型值是关于段参数 [ \Delta t ] 的一元二次函数，可解析求极小点。
- 若段内极小点不可行，则走到下一个断点继续，直到得到广义 Cauchy 点。
- **与其他公式的内在联系**：
- 以第 3.2 节近似 Hessian 为曲率信息，构造可行下降的起点 [ \mathbf{x}^c ]。
- 其输出活动集 [ \mathcal{A}(\mathbf{x}^c) ] 直接传入第 3.4 节子空间最小化。

### 公式 072

- **出现位置**：`explain.md` 第 528-530 行；章节：[### 3.3. Generalized Cauchy point]；编号标注：(20)。
- **公式完整表达**：
[
\mathbf{B} = \mathbf{B}^0 - \mathbf{W}\mathbf{M}\mathbf{W}^{\mathrm{T}}.
]
- **公式含义**：该式处于广义 Cauchy 点（GCP）层：沿投影梯度折线搜索，动态识别活动集并更新一维二次模型导数，保证约束可行下降。
- **必要推导过程**：
- 沿投影梯度轨迹分段前进：每遇到一个变量触边界，就更新活动集。
- 在每段上模型值是关于段参数 [ \Delta t ] 的一元二次函数，可解析求极小点。
- 若段内极小点不可行，则走到下一个断点继续，直到得到广义 Cauchy 点。
- **与其他公式的内在联系**：
- 以第 3.2 节近似 Hessian 为曲率信息，构造可行下降的起点 [ \mathbf{x}^c ]。
- 其输出活动集 [ \mathcal{A}(\mathbf{x}^c) ] 直接传入第 3.4 节子空间最小化。

### 公式 073

- **出现位置**：`explain.md` 第 532-534 行；章节：[### 3.3. Generalized Cauchy point]；编号标注：(23)、(24)。
- **公式完整表达**：
[
\mathbf{w}_b = \mathbf{W}^{\mathrm{T}}\mathbf{e}_b, \quad \mathbf{p}^{(j)} = \mathbf{W}^{\mathrm{T}}\mathbf{d}^{(j)} = \mathbf{W}^{\mathrm{T}}(\mathbf{d}^{(j-1)} + r_b \mathbf{e}_b) = \mathbf{p}^{(j-1)} + r_b \mathbf{w}_b, \quad \mathbf{c}^{(j)} = \mathbf{W}^{\mathrm{T}}\mathbf{z}^{(j)} = \mathbf{W}^{\mathrm{T}}(\mathbf{z}^{(j-1)} + \Delta t^{(j-1)}\mathbf{d}^{(j-1)}) = \mathbf{c}^{(j-1)} + \Delta t^{(j-1)}\mathbf{p}^{(j-1)}.
]
- **公式含义**：该式处于广义 Cauchy 点（GCP）层：沿投影梯度折线搜索，动态识别活动集并更新一维二次模型导数，保证约束可行下降。
- **必要推导过程**：
- 沿投影梯度轨迹分段前进：每遇到一个变量触边界，就更新活动集。
- 在每段上模型值是关于段参数 [ \Delta t ] 的一元二次函数，可解析求极小点。
- 若段内极小点不可行，则走到下一个断点继续，直到得到广义 Cauchy 点。
- **与其他公式的内在联系**：
- 以第 3.2 节近似 Hessian 为曲率信息，构造可行下降的起点 [ \mathbf{x}^c ]。
- 其输出活动集 [ \mathcal{A}(\mathbf{x}^c) ] 直接传入第 3.4 节子空间最小化。

### 公式 074

- **出现位置**：`explain.md` 第 537-539 行；章节：[### 3.3. Generalized Cauchy point]；编号标注：(23)、(24)。
- **公式完整表达**：
[
f'_j = f'_{j-1} + \Delta t^{(j-1)} f''_{j-1} + r^2_b + r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{z}^{(j)} - r_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{c}^{(j)}
]
- **公式含义**：该式处于广义 Cauchy 点（GCP）层：沿投影梯度折线搜索，动态识别活动集并更新一维二次模型导数，保证约束可行下降。
- **必要推导过程**：
- 沿投影梯度轨迹分段前进：每遇到一个变量触边界，就更新活动集。
- 在每段上模型值是关于段参数 [ \Delta t ] 的一元二次函数，可解析求极小点。
- 若段内极小点不可行，则走到下一个断点继续，直到得到广义 Cauchy 点。
- **与其他公式的内在联系**：
- 以第 3.2 节近似 Hessian 为曲率信息，构造可行下降的起点 [ \mathbf{x}^c ]。
- 其输出活动集 [ \mathcal{A}(\mathbf{x}^c) ] 直接传入第 3.4 节子空间最小化。

### 公式 075

- **出现位置**：`explain.md` 第 541-543 行；章节：[### 3.3. Generalized Cauchy point]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
f''_j = f''_{j-1} + 2r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{d}^{(j-1)} - 2r_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{p}^{(j-1)} + r^2_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{e}_b - r^2_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{w}_b.
]
- **公式含义**：该式处于广义 Cauchy 点（GCP）层：沿投影梯度折线搜索，动态识别活动集并更新一维二次模型导数，保证约束可行下降。
- **必要推导过程**：
- 沿投影梯度轨迹分段前进：每遇到一个变量触边界，就更新活动集。
- 在每段上模型值是关于段参数 [ \Delta t ] 的一元二次函数，可解析求极小点。
- 若段内极小点不可行，则走到下一个断点继续，直到得到广义 Cauchy 点。
- **与其他公式的内在联系**：
- 以第 3.2 节近似 Hessian 为曲率信息，构造可行下降的起点 [ \mathbf{x}^c ]。
- 其输出活动集 [ \mathcal{A}(\mathbf{x}^c) ] 直接传入第 3.4 节子空间最小化。

### 公式 076

- **出现位置**：`explain.md` 第 546-548 行；章节：[### 3.3. Generalized Cauchy point]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\mathcal{A}(\mathbf{x}^c) = \{i : t_i \leq t^{(j)}, i = 1, \ldots, n\}.
]
- **公式含义**：该式处于广义 Cauchy 点（GCP）层：沿投影梯度折线搜索，动态识别活动集并更新一维二次模型导数，保证约束可行下降。
- **必要推导过程**：
- 沿投影梯度轨迹分段前进：每遇到一个变量触边界，就更新活动集。
- 在每段上模型值是关于段参数 [ \Delta t ] 的一元二次函数，可解析求极小点。
- 若段内极小点不可行，则走到下一个断点继续，直到得到广义 Cauchy 点。
- **与其他公式的内在联系**：
- 以第 3.2 节近似 Hessian 为曲率信息，构造可行下降的起点 [ \mathbf{x}^c ]。
- 其输出活动集 [ \mathcal{A}(\mathbf{x}^c) ] 直接传入第 3.4 节子空间最小化。

### 公式 077

- **出现位置**：`explain.md` 第 556-558 行；章节：[### 3.4. Subspace minimization]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\mathcal{A}(\mathbf{x}^c) = \{i : x^c_i = \mathrm{lb}_i\} \cup \{i : x^c_i = \mathrm{ub}_i\} = \{i : t_i \leq t^{(j)}, i = 1, \ldots, n\}.
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 078

- **出现位置**：`explain.md` 第 561-563 行；章节：[### 3.4. Subspace minimization]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
x_i = x^c_i = \begin{cases} \mathrm{ub}_i & \text{if } r_i < 0 \\ \mathrm{lb}_i & \text{if } r_i > 0 \end{cases}, \quad \forall i \in \mathcal{A}(\mathbf{x}^c).
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 079

- **出现位置**：`explain.md` 第 569-571 行；章节：[### 3.4. Subspace minimization]；编号标注：(25)。
- **公式完整表达**：
[
\mathbf{x}^* = \arg \min m_k(\mathbf{x}) \quad (25)
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 080

- **出现位置**：`explain.md` 第 575-577 行；章节：[### 3.4. Subspace minimization]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
x_i = x^c_i, \quad \forall i \in \mathcal{A}(\mathbf{x}^c)
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 081

- **出现位置**：`explain.md` 第 579-581 行；章节：[### 3.4. Subspace minimization]；编号标注：(26)、(15)。
- **公式完整表达**：
[
\mathrm{lb}_i \leq x_i \leq \mathrm{ub}_i, \quad \forall i \notin \mathcal{A}(\mathbf{x}^c). \quad (26)
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 082

- **出现位置**：`explain.md` 第 589-591 行；章节：[#### 3.4.1. Direct matrix factorization for the primal approach]；编号标注：(15)。
- **公式完整表达**：
[
\mathbf{x} = \mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}},
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 083

- **出现位置**：`explain.md` 第 595-604 行；章节：[#### 3.4.1. Direct matrix factorization for the primal approach]；编号标注：(15)。
- **公式完整表达**：
[
\begin{array}{rl} 
m_k(\mathbf{x}) &= \Pi_k + \mathbf{r}_k^{\mathrm{T}}(\mathbf{x} - \mathbf{x}_k) + \frac{1}{2}(\mathbf{x} - \mathbf{x}_k)^{\mathrm{T}}\mathbf{B}_k(\mathbf{x} - \mathbf{x}_k) \\
\implies \hat{m}_k(\hat{\mathbf{x}}) &= \Pi_k + \mathbf{r}_k^{\mathrm{T}}(\mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}} - \mathbf{x}_k) + \frac{1}{2}(\mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}} - \mathbf{x}_k)^{\mathrm{T}}\mathbf{B}_k(\mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}} - \mathbf{x}_k) \\
&= \Pi_k + \mathbf{r}_k^{\mathrm{T}}(\mathbf{x}^c - \mathbf{x}_k) + \frac{1}{2} (\mathbf{x}^c - \mathbf{x}_k)^{\mathrm{T}}\mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k) \\
&\quad + [\mathbf{r}_k + \mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k)]^{\mathrm{T}}\mathbf{Z}_k \hat{\mathbf{x}} + \frac{1}{2} \hat{\mathbf{x}}^{\mathrm{T}}\mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k \hat{\mathbf{x}} \\
&= \Pi_k + \mathbf{r}_k^{\mathrm{T}}(\mathbf{x}^c - \mathbf{x}_k) + \frac{1}{2} (\mathbf{x}^c - \mathbf{x}_k)^{\mathrm{T}}\mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k) \\
&\quad + [\mathbf{Z}^{\mathrm{T}}_k [\mathbf{r}_k + \mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k)]]^{\mathrm{T}} \hat{\mathbf{x}} + \frac{1}{2} \hat{\mathbf{x}}^{\mathrm{T}}\mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k \hat{\mathbf{x}}.
\end{array}
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 084

- **出现位置**：`explain.md` 第 608-610 行；章节：[#### 3.4.1. Direct matrix factorization for the primal approach]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\hat{m}'_k(\hat{\mathbf{x}}) = 0 \implies \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k \hat{\mathbf{x}} = -\mathbf{Z}^{\mathrm{T}}_k [\mathbf{r}_k + \mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k)].
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 085

- **出现位置**：`explain.md` 第 614-616 行；章节：[#### 3.4.1. Direct matrix factorization for the primal approach]；编号标注：(27)。
- **公式完整表达**：
[
\hat{\mathbf{B}}_k \hat{\mathbf{x}} = -\hat{\mathbf{r}}_k. \quad (27)
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 086

- **出现位置**：`explain.md` 第 619-621 行；章节：[#### 3.4.1. Direct matrix factorization for the primal approach]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\mathbf{B}_k = \mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k.
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 087

- **出现位置**：`explain.md` 第 623-625 行；章节：[#### 3.4.1. Direct matrix factorization for the primal approach]；编号标注：(28)、(27)。
- **公式完整表达**：
[
\hat{\mathbf{B}}_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k = \mathbf{Z}^{\mathrm{T}}_k (\mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k) \mathbf{Z}_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{Z}_k - \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k. \quad (28)
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 088

- **出现位置**：`explain.md` 第 631-633 行；章节：[#### 3.4.1. Direct matrix factorization for the primal approach]；编号标注：(29)。
- **公式完整表达**：
[
\hat{\mathbf{x}} = -\hat{\mathbf{B}}^{-1}_k \hat{\mathbf{r}}_k,
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 089

- **出现位置**：`explain.md` 第 634-636 行；章节：[#### 3.4.1. Direct matrix factorization for the primal approach]；编号标注：(29)。
- **公式完整表达**：
[
\hat{\mathbf{B}}^{-1}_k = (\hat{\mathbf{B}}^0_k)^{-1} + (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \left[\mathbf{I} - \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \right]^{-1} \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1}, \quad (29)
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 090

- **出现位置**：`explain.md` 第 639-641 行；章节：[#### 3.4.1. Direct matrix factorization for the primal approach]；编号标注：(29)。
- **公式完整表达**：
[
\mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k = \mathbf{Z}^{\mathrm{T}}_k [\mathbf{Y}_k \quad \mathbf{B}^0_k \mathbf{S}_k] = [\mathbf{Z}^{\mathrm{T}}_k \mathbf{y}_{k-m} \cdots \mathbf{Z}^{\mathrm{T}}_k \mathbf{y}_{k-1} \quad \mathbf{Z}^{\mathrm{T}}_k(\mathbf{B}^0_k \mathbf{s}_{k-m}) \cdots \mathbf{Z}^{\mathrm{T}}_k(\mathbf{B}^0_k \mathbf{s}_{k-1})].
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 091

- **出现位置**：`explain.md` 第 644-646 行；章节：[#### 3.4.1. Direct matrix factorization for the primal approach]；编号标注：(29)。
- **公式完整表达**：
[
\mathbf{I} - \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \in \mathbb{R}^{2m\times 2m}
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 092

- **出现位置**：`explain.md` 第 655-657 行；章节：[#### 3.4.2. Conjugate gradient method for the primal approach]；编号标注：(18)。
- **公式完整表达**：
[
\mathbf{P}^{-1} = \mathrm{ILU}(\hat{\mathbf{B}}^0_k) = \mathrm{ILU}(\mathbf{Z}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{Z}_k).
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 093

- **出现位置**：`explain.md` 第 665-667 行；章节：[#### 3.4.3. Schur complement for the dual approach]；编号标注：(25)、(26)、(15)。
- **公式完整表达**：
[
\mathbf{x} = \mathbf{x}_k + \Delta \mathbf{x}_k.
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 094

- **出现位置**：`explain.md` 第 671-673 行；章节：[#### 3.4.3. Schur complement for the dual approach]；编号标注：(15)。
- **公式完整表达**：
[
m_k(\mathbf{x}) = \Pi_k + \mathbf{r}_k^{\mathrm{T}} \Delta \mathbf{x}_k + \frac{1}{2} \Delta \mathbf{x}_k^{\mathrm{T}} \mathbf{B}_k \Delta \mathbf{x}_k.
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 095

- **出现位置**：`explain.md` 第 676-678 行；章节：[#### 3.4.3. Schur complement for the dual approach]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
x_i = x^c_i, \quad \forall i \in \mathcal{A}_k(\mathbf{x}^c),
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 096

- **出现位置**：`explain.md` 第 680-682 行；章节：[#### 3.4.3. Schur complement for the dual approach]；编号标注：(25)、(26)。
- **公式完整表达**：
[
\mathbf{Q}^{\mathrm{T}}_k \Delta \mathbf{x}_k = \mathbf{Q}^{\mathrm{T}}_k (\mathbf{x} - \mathbf{x}_k) = \mathbf{Q}^{\mathrm{T}}_k (\mathbf{x}^c - \mathbf{x}_k).
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 097

- **出现位置**：`explain.md` 第 685-687 行；章节：[#### 3.4.3. Schur complement for the dual approach]；编号标注：(25)、(26)。
- **公式完整表达**：
[
\min \left\{ \Pi_k + \mathbf{r}_k^{\mathrm{T}} \Delta \mathbf{x}_k + \frac{1}{2} \Delta \mathbf{x}_k^{\mathrm{T}} \mathbf{B}_k \Delta \mathbf{x}_k \right\}
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 098

- **出现位置**：`explain.md` 第 689-691 行；章节：[#### 3.4.3. Schur complement for the dual approach]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\mathbf{Q}_k^{\mathrm{T}} \Delta \mathbf{x}_k = \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k)
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 099

- **出现位置**：`explain.md` 第 693-695 行；章节：[#### 3.4.3. Schur complement for the dual approach]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\mathbf{lb} - \mathbf{x}_k \leq \Delta \mathbf{x}_k \leq \mathbf{ub} - \mathbf{x}_k.
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 100

- **出现位置**：`explain.md` 第 699-701 行；章节：[#### 3.4.3. Schur complement for the dual approach]；编号标注：(30)。
- **公式完整表达**：
[
\begin{bmatrix} \mathbf{B}_k & \mathbf{Q}_k \\ \mathbf{Q}_k^{\mathrm{T}} & \mathbf{0} \end{bmatrix} \begin{bmatrix} \Delta \mathbf{x}_k \\ \lambda \end{bmatrix} = \begin{bmatrix} -\mathbf{r}_k \\ \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k) \end{bmatrix}. \quad (30)
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 101

- **出现位置**：`explain.md` 第 705-707 行；章节：[#### 3.4.3. Schur complement for the dual approach]；编号标注：(31)。
- **公式完整表达**：
[
(\mathbf{Q}_k^{\mathrm{T}} \mathbf{B}_k^{-1} \mathbf{Q}_k) \lambda = -\mathbf{Q}_k^{\mathrm{T}} \mathbf{B}_k^{-1} \mathbf{r}_k - \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k). \quad (31)
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 102

- **出现位置**：`explain.md` 第 711-713 行；章节：[#### 3.4.3. Schur complement for the dual approach]；编号标注：(32)、(31)。
- **公式完整表达**：
[
\mathbf{B}_k \Delta \mathbf{x}_k = -(\mathbf{r}_k + \mathbf{Q}_k \lambda). \quad (32)
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 103

- **出现位置**：`explain.md` 第 719-721 行；章节：[#### 3.4.3. Schur complement for the dual approach]；编号标注：(27)。
- **公式完整表达**：
[
\mathbf{p}_k = \mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}} - \mathbf{x}_k.
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 104

- **出现位置**：`explain.md` 第 725-727 行；章节：[#### 3.4.3. Schur complement for the dual approach]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\mathbf{x}_{k + 1} = \mathbf{x}_k + \alpha_k \mathbf{p}_k,
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 105

- **出现位置**：`explain.md` 第 732-734 行；章节：[#### 3.4.3. Schur complement for the dual approach]；编号标注：(18)。
- **公式完整表达**：
[
\mathbf{lb} \leq \mathbf{x}_k + \alpha_k \mathbf{p}_k \leq \mathbf{ub}.
]
- **公式含义**：该式处于子空间最小化层：固定活动集上的变量，在自由子空间内解降维二次子问题；可采用原始法或对偶 Schur 补法。
- **必要推导过程**：
- 固定活动集变量后，引入基矩阵 [ \mathbf{Z}_k ]（自由子空间）或约束矩阵 [ \mathbf{Q}_k ]（对偶）。
- 原始法：解降维线性系统 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ]。
- 对偶法：由 KKT 方程消元得到 Schur 补，再回代求步长。
- **与其他公式的内在联系**：
- 向前承接 GCP 给出的活动集与可行点。
- 向后输出搜索方向 [ \mathbf{p}_k ] 与更新 [ \mathbf{x}_{k+1} ]，回到第 3.1 节线搜索与收敛判据。

### 公式 106

- **出现位置**：`explain.md` 第 780-782 行；章节：[### 4.1. Cyclic tension-compression test]；编号标注：(33)。
- **公式完整表达**：
[
E_{\mathrm{crack}} = g_{\mathrm{c}}\Gamma_{l}(d) = g_{\mathrm{c}}\int_{\Omega}\gamma (d,\nabla d)\mathrm{d}\Omega = g_{\mathrm{c}}\int_{\Omega}\frac{1}{2l}\left(d^{2} + l^{2}\nabla d\cdot \nabla d\right)\mathrm{d}\Omega \quad (33)
]
- **公式含义**：该式处于数值算例层：定义载荷路径、裂纹能或边界位移函数，用于把算法输出与物理实验量关联。
- **必要推导过程**：
- 将通用模型中的能量/位移量映射为具体算例控制量（循环载荷、扭转角等）。
- 这些式子通常是边界条件或后处理量定义，不改变求解框架但决定结果解读维度。
- **与其他公式的内在联系**：
- 向前承接前述统一模型与算法。
- 向后把解转化为可比较的工程指标（裂纹能、荷载-位移、扭转响应）。

### 公式 107

- **出现位置**：`explain.md` 第 878-880 行；章节：[### 4.4. Three-dimensional torsion test]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
u_{y} = z \tan t, \quad u_{z} = -y \tan t,
]
- **公式含义**：该式处于数值算例层：定义载荷路径、裂纹能或边界位移函数，用于把算法输出与物理实验量关联。
- **必要推导过程**：
- 将通用模型中的能量/位移量映射为具体算例控制量（循环载荷、扭转角等）。
- 这些式子通常是边界条件或后处理量定义，不改变求解框架但决定结果解读维度。
- **与其他公式的内在联系**：
- 向前承接前述统一模型与算法。
- 向后把解转化为可比较的工程指标（裂纹能、荷载-位移、扭转响应）。

### 公式 108

- **出现位置**：`explain.md` 第 937-939 行；章节：[### 5.2. Comparison of convergence behaviors]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\| \pmb {r}_u\| _2 < \mathrm{tol}, \quad \| \pmb {r}_d\| _2 < \mathrm{tol}, \quad \| \Delta \pmb {u}\| _2 < \mathrm{tol}, \quad \| \Delta d\| _2 < \mathrm{tol},
]
- **公式含义**：该式处于收敛评估层：用残量范数与增量范数判定非线性迭代终止，体现算法稳定性与计算成本。
- **必要推导过程**：
- 由非线性方程求解的必要条件 [ \mathbf{r}\to 0 ] 与增量收敛 [ \Delta\mathbf{x}\to 0 ] 共同构造停止准则。
- 双重准则可避免“残量小但变量仍震荡”或“步长小但未平衡”的误判。
- **与其他公式的内在联系**：
- 与第 3.1 节停止准则呼应，用于比较不同线性求解器/预条件器的鲁棒性。

### 公式 109

- **出现位置**：`explain.md` 第 973-975 行；章节：[## Appendix]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\hat{\mathbf{A}} = \mathbf{A} + \mathbf{U}\mathbf{V}^{\mathrm{T}},
]
- **公式含义**：该式处于线性代数附录层：利用低秩修正逆矩阵公式，把大系统求逆转化为小维系统求逆，提高计算效率。
- **必要推导过程**：
- 把矩阵写成“基矩阵 + 低秩修正”形式。
- 使用 Sherman–Morrison–Woodbury 恒等式，将高维逆运算改写为低维矩阵 [ \mathbf{I}+\mathbf{V}^T\mathbf{A}^{-1}\mathbf{U} ] 的求逆。
- **与其他公式的内在联系**：
- 为第 3 节中反复出现的“低秩修正矩阵逆”提供理论工具。
- 与 L-BFGS 紧凑表示天然耦合，解释为何算法可高效实现。

### 公式 110

- **出现位置**：`explain.md` 第 979-981 行；章节：[## Appendix]；编号标注：未显式编号（上下文公式）。
- **公式完整表达**：
[
\hat{\mathbf{A}}^{-1} = \mathbf{A}^{-1} - \mathbf{A}^{-1}\mathbf{U}\left(\mathbf{I} + \mathbf{V}^{\mathrm{T}}\mathbf{A}^{-1}\mathbf{U}\right)^{-1}\mathbf{V}^{\mathrm{T}}\mathbf{A}^{-1}.
]
- **公式含义**：该式处于线性代数附录层：利用低秩修正逆矩阵公式，把大系统求逆转化为小维系统求逆，提高计算效率。
- **必要推导过程**：
- 把矩阵写成“基矩阵 + 低秩修正”形式。
- 使用 Sherman–Morrison–Woodbury 恒等式，将高维逆运算改写为低维矩阵 [ \mathbf{I}+\mathbf{V}^T\mathbf{A}^{-1}\mathbf{U} ] 的求逆。
- **与其他公式的内在联系**：
- 为第 3 节中反复出现的“低秩修正矩阵逆”提供理论工具。
- 与 L-BFGS 紧凑表示天然耦合，解释为何算法可高效实现。

## 关键编号公式之间的主干联系（跨章节汇总）

1. **建模主干**：
   - [ (1) ] 定义总势能；[ (2) ] 给出裂纹相场正则项；[ (3)(4) ] 给出增量最小化与不可逆约束。  
   - 这四式共同定义“物理正确 + 数值可解”的起点。

2. **离散主干**：
   - [ (7) ] 一阶变分得到残量；[ (8)(9)(10) ] 离散化并得到切线矩阵；[ (11) ] 给出块对角近似。  
   - 该链条把连续偏微分问题转为代数优化问题。

3. **优化主干**：
   - [ (12)(13) ] 盒约束最小化；[ (15) ] 局部二次模型；[ (16)(17)(18)(19) ] 投影、Cauchy 点与停止判据。  
   - 对应“可行性 + 下降性 + 收敛性”的核心机制。

4. **L-BFGS 与高效线代主干**：
   - [ (20)(21) ] 紧凑 L-BFGS 表示与初值矩阵；[ (22)(23)(24) ] GCP 段内导数更新；[ (27)(28)(29) ] 子空间求解与逆更新。  
   - 附录 SMW 公式解释了上述逆运算为何可在低维空间高效实现。

5. **KKT/对偶主干**：
   - [ (30)(31)(32) ] 从带约束二次子问题到 KKT 与 Schur 补系统，再回代得到步长。  
   - 与原始子空间法互补，提高了大规模问题的求解灵活性。

6. **算例与验证主干**：
   - [ (33) ] 等式把裂纹表面密度积分转为可观测裂纹能。  
   - 与载荷位移边界公式共同构成“算法结果 → 物理指标”的解释闭环。

## 备注

- 本文档严格覆盖 `explain.md` 中提取到的全部 110 个公式块。
- 若你希望，我可以在下一轮把“每条公式的推导”进一步扩写为**逐步代数展开版**（例如把 [ (7) ]、[ (10) ]、[ (20) ]、[ (30) ] 分解成逐行推导）。
