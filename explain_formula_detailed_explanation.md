# explain.md 论文全部公式逐一解释（位置 + 完整表达 + 推导 + 联系）

> 对象文件：`explain.md`

> 说明：以下按公式出现顺序逐一解释，每条包含：出现位置、完整表达、含义、必要推导、内在联系。

> 渲染说明：按需求统一采用 `[ ... ]` 形式呈现行内符号与公式（这是一种文档风格约定）。

## 逐条公式解释


### 公式 001

- **出现位置**：`explain.md` 第 21-23 行；章节：[ ## 1. Introduction ]；文中编号：[ (1) ]。
- **公式完整表达**：
[
\Pi (\pmb {u},d) = \int_{\Omega}\psi (\pmb {e}(\pmb {u}),d)\mathrm{d}\Omega + g_{c}\Gamma_{l}(d) - \int_{\Omega}\pmb {b}\cdot \pmb {u}\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot \pmb {u}\mathrm{d}\Gamma , \quad (1)
]
- **公式含义**：
  - 定义总势能目标 [ \Pi ]，把内能、断裂能与外载统一到同一优化框架。
  - 给出裂纹正则项 [ \Gamma_l ] 或裂纹密度 [ \gamma ]，控制裂纹带宽与能量耗散。
  - 邻近语境（简要）：前句“During the past decade, the phase-field method becomes a popular technique for modeling fracture propagation, particularly in brittle materials, due to its cap…”；后句“where \(\pmb {u}(\pmb {x})\) represents the unknown vector displacement field, \(d(\pmb {x})\) represents the unknown scalar phase-field, \(\pmb{b}\) is the bo…”。
- **必要推导过程**：
  1. 从最小势能原理写出内能减外力势，并加入断裂表面能项。
  2. 采用相场正则化 [ \gamma=\frac{1}{2l}(d^2+l^2|\nabla d|^2) ]，对域积分得到 [ \Gamma_l ]。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。
  - 与文中编号公式互引：[ (1) ]。

### 公式 002

- **出现位置**：`explain.md` 第 27-29 行；章节：[ ## 1. Introduction ]；文中编号：[ (2) ]。
- **公式完整表达**：
[
\Gamma_{l}(d) = \int_{\Omega}\gamma (d,\nabla d)\mathrm{d}\Omega = \int_{\Omega}\frac{1}{2l}\left(d^{2} + l^{2}\nabla d\cdot \nabla d\right)\mathrm{d}\Omega , \quad (2)
]
- **公式含义**：
  - 给出裂纹正则项 [ \Gamma_l ] 或裂纹密度 [ \gamma ]，控制裂纹带宽与能量耗散。
  - 邻近语境（简要）：前句“where \(\pmb {u}(\pmb {x})\) represents the unknown vector displacement field, \(d(\pmb {x})\) represents the unknown scalar phase-field, \(\pmb{b}\) is the bo…”；后句“where \(\gamma (d,\nabla d)\) is considered as the crack surface density function, and \(l\) is the phase-field length-scale parameter.”。
- **必要推导过程**：
  1. 采用相场正则化 [ \gamma=\frac{1}{2l}(d^2+l^2|\nabla d|^2) ]，对域积分得到 [ \Gamma_l ]。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。
  - 与文中编号公式互引：[ (2) ]。

### 公式 003

- **出现位置**：`explain.md` 第 35-37 行；章节：[ ## 1. Introduction ]；文中编号：[ (3) ]。
- **公式完整表达**：
[
\left(\pmb{u}_{n + 1},d_{n + 1}\right) = \arg \min \Pi (\pmb {u},d) \quad (3)
]
- **公式含义**：
  - 定义总势能目标 [ \Pi ]，把内能、断裂能与外载统一到同一优化框架。
  - 将求解过程明确为带约束极小化问题。
  - 邻近语境（简要）：前句“For a quasi-static problem, let \(t\) represent the pseudo time instead of the real time. The pseudo time \(t\) enters the problem as the load step through tim…”；后句“subject to the inequality constraints”。
- **必要推导过程**：
  1. 从最小势能原理写出内能减外力势，并加入断裂表面能项。
  2. 把增量步未知量作为优化变量，配合约束写成极小化问题。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。
  - 与文中编号公式互引：[ (3) ]。

### 公式 004

- **出现位置**：`explain.md` 第 41-43 行；章节：[ ## 1. Introduction ]；文中编号：[ (4) ]。
- **公式完整表达**：
[
0\leq d_{n}\leq d_{n + 1}\leq 1. \quad (4)
]
- **公式含义**：
  - 用于盒约束可行性维护，防止相场变量违反不可逆与上界约束。
  - 邻近语境（简要）：前句“subject to the inequality constraints”；后句“The above inequality constraints represent the following considerations. First, the phase-field cannot decrease in order to ensure the thermodynamic consistenc…”。
- **必要推导过程**：
  1. 按分量截断规则实施投影：低于下界取下界，高于上界取上界。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与文中编号公式互引：[ (4) ]。

### 公式 005

- **出现位置**：`explain.md` 第 57-59 行；章节：[ ## 1. Introduction ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathrm{Proj}_C(\mathbf{x}_k = \mathbf{a}_k\nabla f(\mathbf{x}_k)).
]
- **公式含义**：
  - 用于盒约束可行性维护，防止相场变量违反不可逆与上界约束。
  - 邻近语境（简要）：前句“The gradient projection method [39-42] is a special type of active set methods that can be used to solve the inequality-constrained optimization problem. This …”；后句“For a general convex feasible region, for example, as shown in Fig. 1(a), the projection operation \(\mathrm{Proj}_C(\cdot)\) incurs significant computational …”。
- **必要推导过程**：
  1. 按分量截断规则实施投影：低于下界取下界，高于上界取上界。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 006

- **出现位置**：`explain.md` 第 63-65 行；章节：[ ## 1. Introduction ]；文中编号：[ (4) ]。
- **公式完整表达**：
[
\mathrm{lb}_i\leq x_i\leq \mathrm{ub}_i,
]
- **公式含义**：
  - 用于盒约束可行性维护，防止相场变量违反不可逆与上界约束。
  - 邻近语境（简要）：前句“For a general convex feasible region, for example, as shown in Fig. 1(a), the projection operation \(\mathrm{Proj}_C(\cdot)\) incurs significant computational …”；后句“where \(x_{i}\) represents the \(i\) th component of the unknown vector \(\mathbf{x}\) \(\mathrm{lb}_i\) and \(\mathrm{ub}_i\) represent its corresponding lowe…”。
- **必要推导过程**：
  1. 按分量截断规则实施投影：低于下界取下界，高于上界取上界。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与文中编号公式互引：[ (4) ]。

### 公式 007

- **出现位置**：`explain.md` 第 69-71 行；章节：[ ## 1. Introduction ]；文中编号：[ (4) ]。
- **公式完整表达**：
[
\mathrm{lb}_i = d_i^{(n)}\leq d_i^{(n + 1)}\leq 1 = \mathrm{ub}_i,
]
- **公式含义**：
  - 用于盒约束可行性维护，防止相场变量违反不可逆与上界约束。
  - 邻近语境（简要）：前句“where \(x_{i}\) represents the \(i\) th component of the unknown vector \(\mathbf{x}\) \(\mathrm{lb}_i\) and \(\mathrm{ub}_i\) represent its corresponding lowe…”；后句“where \(d_i^{(n)}\) is known at the beginning of the current time (load) step, and \(d_i^{(n + 1)}\) is the unknown phase-field value that needs to be solved.”。
- **必要推导过程**：
  1. 按分量截断规则实施投影：低于下界取下界，高于上界取上界。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与文中编号公式互引：[ (4) ]。

### 公式 008

- **出现位置**：`explain.md` 第 87-89 行；章节：[ ### 2.1. Phase-field formulation ]；文中编号：[ (1) ]、[ (5) ]。
- **公式完整表达**：
[
\psi(\pmb{\epsilon}, d) = [g(d) + k] \psi^{+}(\pmb{\epsilon}) + \psi^{-}(\pmb{\epsilon}), \quad (5)
]
- **公式含义**：
  - 给出退化函数，决定损伤对张拉能的衰减规律。
  - 邻近语境（简要）：前句“In order to consider the fracture tension–compression asymmetry, the strain energy density function used in the total energy functional in Eq. (1) is additivel…”；后句“where \(\psi^{+}\) is the positive strain energy, \(\psi^{-}\) is the negative strain energy, \(k\) is a small non-negative number, and \(g(d)\) represents the…”。
- **必要推导过程**：
  1. 将退化函数乘于张拉能并保留微小残余刚度项，防止病态奇异。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。
  - 与文中编号公式互引：[ (1) ]、[ (5) ]。

### 公式 009

- **出现位置**：`explain.md` 第 93-95 行；章节：[ ### 2.1. Phase-field formulation ]；文中编号：[ (6) ]。
- **公式完整表达**：
[
g(d) = (1 - d)^{2}. \quad (6)
]
- **公式含义**：
  - 给出退化函数，决定损伤对张拉能的衰减规律。
  - 邻近语境（简要）：前句“where \(\psi^{+}\) is the positive strain energy, \(\psi^{-}\) is the negative strain energy, \(k\) is a small non-negative number, and \(g(d)\) represents the…”；后句“In this work, the small non-negative number \(k\) is set as zero (\(k = 0\)) in all the numerical examples.”。
- **必要推导过程**：
  1. 将退化函数乘于张拉能并保留微小残余刚度项，防止病态奇异。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。
  - 与文中编号公式互引：[ (6) ]。

### 公式 010

- **出现位置**：`explain.md` 第 101-103 行；章节：[ ### 2.1. Phase-field formulation ]；文中编号：未显式编号。
- **公式完整表达**：
[
\langle x\rangle_{+} = \frac{1}{2}(x + |x|), \quad \langle x\rangle_{-} = \frac{1}{2}(x - |x|), \quad H(x) = \begin{cases} 1 & x \geq 0, \\ 0 & x < 0. \end{cases}
]
- **公式含义**：
  - 实现张压分裂，避免压缩态非物理裂纹扩展。
  - 邻近语境（简要）：前句“The following operators are introduced to describe the constitutive relationship based on the additive decomposition of the strain energy density function,”；后句“The spectrum decomposition of the strain tensor \(\pmb{\epsilon}\) is expressed as”。
- **必要推导过程**：
  1. 先做谱分解，再用 Macaulay 括号分离正负主应变贡献。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。

### 公式 011

- **出现位置**：`explain.md` 第 107-109 行；章节：[ ### 2.1. Phase-field formulation ]；文中编号：未显式编号。
- **公式完整表达**：
[
\pmb{\epsilon} = \sum_{\alpha} \epsilon_{\alpha} \mathbf{M}_{\alpha}, \quad \mathbf{M}_{\alpha} = \pmb{n}_{\alpha} \otimes \pmb{n}_{\alpha},
]
- **公式含义**：
  - 实现张压分裂，避免压缩态非物理裂纹扩展。
  - 邻近语境（简要）：前句“The spectrum decomposition of the strain tensor \(\pmb{\epsilon}\) is expressed as”；后句“where \(\epsilon_{\alpha}\) and \(\pmb{n}_{\alpha}\) represent a pair of eigenvalue and eigenvector. The positive and negative parts of the strain tensor are d…”。
- **必要推导过程**：
  1. 先做谱分解，再用 Macaulay 括号分离正负主应变贡献。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。

### 公式 012

- **出现位置**：`explain.md` 第 113-115 行；章节：[ ### 2.1. Phase-field formulation ]；文中编号：未显式编号。
- **公式完整表达**：
[
\pmb{\epsilon}^{+} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{+} \mathbf{M}_{\alpha}, \quad \pmb{\epsilon}^{-} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{-} \mathbf{M}_{\alpha}.
]
- **公式含义**：
  - 实现张压分裂，避免压缩态非物理裂纹扩展。
  - 邻近语境（简要）：前句“where \(\epsilon_{\alpha}\) and \(\pmb{n}_{\alpha}\) represent a pair of eigenvalue and eigenvector. The positive and negative parts of the strain tensor are d…”；后句“Using the above definitions, the positive and negative parts of the strain energy are expressed as”。
- **必要推导过程**：
  1. 先做谱分解，再用 Macaulay 括号分离正负主应变贡献。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。

### 公式 013

- **出现位置**：`explain.md` 第 119-121 行；章节：[ ### 2.1. Phase-field formulation ]；文中编号：未显式编号。
- **公式完整表达**：
[
\psi^{+}(\pmb{\epsilon}) = \frac{1}{2}\lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{+}^{2} + \mu \pmb{\epsilon}^{+} : \pmb{\epsilon}^{+}, \quad \psi^{-}(\pmb{\epsilon}) = \frac{1}{2} \lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{-}^{2} + \mu \pmb{\epsilon}^{-} : \pmb{\epsilon}^{-},
]
- **公式含义**：
  - 实现张压分裂，避免压缩态非物理裂纹扩展。
  - 属于子空间最小化/KKT-Scher 补链，生成最终搜索方向。
  - 邻近语境（简要）：前句“Using the above definitions, the positive and negative parts of the strain energy are expressed as”；后句“where \(\lambda\) and \(\mu\) are the Lame parameters, and \(\mathrm{tr}\pmb{\epsilon}\) is the trace of the strain tensor. The stress tensor \(\pmb{\sigma}\) …”。
- **必要推导过程**：
  1. 先做谱分解，再用 Macaulay 括号分离正负主应变贡献。
  2. 固定活动变量后，解原始降维系统或对偶 Schur 补系统并回代。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 014

- **出现位置**：`explain.md` 第 125-127 行；章节：[ ### 2.1. Phase-field formulation ]；文中编号：未显式编号。
- **公式完整表达**：
[
\pmb{\sigma} = \frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial \pmb{\epsilon}} = [g(d) + k]\frac{\partial\psi^{+}(\pmb{\epsilon})}{\partial \pmb{\epsilon}} + \frac{\partial\psi^{-}(\pmb{\epsilon})}{\partial \pmb{\epsilon}} = [g(d) + k]\pmb{\sigma}^{+} + \pmb{\sigma}^{-},
]
- **公式含义**：
  - 给出退化函数，决定损伤对张拉能的衰减规律。
  - 给出一致切线/雅可比信息，影响牛顿或准牛顿迭代收敛效率。
  - 邻近语境（简要）：前句“where \(\lambda\) and \(\mu\) are the Lame parameters, and \(\mathrm{tr}\pmb{\epsilon}\) is the trace of the strain tensor. The stress tensor \(\pmb{\sigma}\) …”；后句“where”。
- **必要推导过程**：
  1. 将退化函数乘于张拉能并保留微小残余刚度项，防止病态奇异。
  2. 对能量对应变求一阶导得到应力。
  3. 再对应力/残量求导得到一致切线矩阵或投影算子导数。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。

### 公式 015

- **出现位置**：`explain.md` 第 131-133 行；章节：[ ### 2.1. Phase-field formulation ]；文中编号：未显式编号。
- **公式完整表达**：
[
\pmb{\sigma}^{+} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{+}\mathbf{I} + 2\mu \pmb{\epsilon}^{+}, \quad \pmb{\sigma}^{-} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{-}\mathbf{I} + 2\mu \pmb{\epsilon}^{-}.
]
- **公式含义**：
  - 实现张压分裂，避免压缩态非物理裂纹扩展。
  - 由能量导出应力，是平衡方程残量的核心材料项。
  - 属于子空间最小化/KKT-Scher 补链，生成最终搜索方向。
  - 邻近语境（简要）：前句“where”；后句“Due to the crack tension-compression asymmetry, the phase-field degradation function is exclusively applied to the positive part of the strain energy. Therefor…”。
- **必要推导过程**：
  1. 先做谱分解，再用 Macaulay 括号分离正负主应变贡献。
  2. 对能量对应变求一阶导得到应力。
  3. 固定活动变量后，解原始降维系统或对偶 Schur 补系统并回代。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 016

- **出现位置**：`explain.md` 第 137-139 行；章节：[ ### 2.1. Phase-field formulation ]；文中编号：未显式编号。
- **公式完整表达**：
[
\frac{\partial \pmb{\sigma}}{\partial \pmb{\epsilon}} = [g(d) + k]\frac{\partial \pmb{\sigma}^{+}}{\partial \pmb{\epsilon}} + \frac{\partial \pmb{\sigma}^{-}}{\partial \pmb{\epsilon}} = [g(d) + k]\left[\lambda H(\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes \mathbf{I} + 2\mu \mathbb{P}^{+}\right] + \left[\lambda H(-\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes \mathbf{I} + 2\mu \mathbb{P}^{-}\right].
]
- **公式含义**：
  - 给出退化函数，决定损伤对张拉能的衰减规律。
  - 实现张压分裂，避免压缩态非物理裂纹扩展。
  - 给出一致切线/雅可比信息，影响牛顿或准牛顿迭代收敛效率。
  - 属于子空间最小化/KKT-Scher 补链，生成最终搜索方向。
  - 邻近语境（简要）：前句“Due to the crack tension-compression asymmetry, the phase-field degradation function is exclusively applied to the positive part of the strain energy. Therefor…”；后句“Inside the above tangent modulus, the two fourth-order projection tensors \(\mathbb{P}^{+}\) and \(\mathbb{P}^{-}\) are defined as”。
- **必要推导过程**：
  1. 将退化函数乘于张拉能并保留微小残余刚度项，防止病态奇异。
  2. 先做谱分解，再用 Macaulay 括号分离正负主应变贡献。
  3. 对能量对应变求一阶导得到应力。
  4. 再对应力/残量求导得到一致切线矩阵或投影算子导数。
  5. 固定活动变量后，解原始降维系统或对偶 Schur 补系统并回代。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 017

- **出现位置**：`explain.md` 第 143-145 行；章节：[ ### 2.1. Phase-field formulation ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbb{P}^{+} = \frac{\partial \pmb{\epsilon}^{+}}{\partial \pmb{\epsilon}}, \quad \mathbb{P}^{-} = \frac{\partial \pmb{\epsilon}^{-}}{\partial \pmb{\epsilon}},
]
- **公式含义**：
  - 给出一致切线/雅可比信息，影响牛顿或准牛顿迭代收敛效率。
  - 邻近语境（简要）：前句“Inside the above tangent modulus, the two fourth-order projection tensors \(\mathbb{P}^{+}\) and \(\mathbb{P}^{-}\) are defined as”；后句“the specific expressions of which can be found in [26].”。
- **必要推导过程**：
  1. 再对应力/残量求导得到一致切线矩阵或投影算子导数。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。

### 公式 018

- **出现位置**：`explain.md` 第 151-158 行；章节：[ ### 2.1. Phase-field formulation ]；文中编号：[ (7) ]。
- **公式完整表达**：
[
\begin{array}{rl} 
\delta \Pi (\pmb {u},d) &= D_{(\delta \pmb {u},\delta d)}\Pi (\pmb {u},d) = \left.\frac{\mathrm{d}}{\mathrm{d}\epsilon}\right|_{\epsilon = 0}\Pi (\pmb {u} + \epsilon \delta \pmb {u},d + \epsilon \delta d) \\ 
&= \int_{\Omega}\left(\frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial \pmb{\epsilon}}:\pmb{\epsilon}(\delta \pmb {u}) + \frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial d}\delta d\right)\mathrm{d}\Omega \\ 
&\quad +\int_{\Omega}\frac{g_{c}}{l}\left(d\delta d + l^{2}\nabla d\cdot \nabla \delta d\right)\mathrm{d}\Omega -\int_{\Omega}\pmb {b}\cdot \delta \pmb {u}\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot \delta \pmb {u}\mathrm{d}\Gamma \\ 
&= (\nabla^{(s)}\delta \pmb {u},\pmb {\sigma}) - (\delta \pmb {u},\pmb {b}) - (\delta \pmb {u},\pmb{t})_{\Gamma_{t}} + (\delta d,\frac{g_{c}}{l} d) + (\nabla \delta d,g_{c}l\nabla d) + (\delta d,g^{\prime}(d)\psi^{+}(\pmb{\epsilon})). 
\end{array} \quad (7)
]
- **公式含义**：
  - 定义总势能目标 [ \Pi ]，把内能、断裂能与外载统一到同一优化框架。
  - 由能量导出应力，是平衡方程残量的核心材料项。
  - 邻近语境（简要）：前句“Using the directional derivative, the first variation of the energy functional is written as”；后句“Therefore, the weak form of the phase-field formulation is expressed as”。
- **必要推导过程**：
  1. 从最小势能原理写出内能减外力势，并加入断裂表面能项。
  2. 对能量对应变求一阶导得到应力。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。
  - 与文中编号公式互引：[ (7) ]。

### 公式 019

- **出现位置**：`explain.md` 第 162-164 行；章节：[ ### 2.1. Phase-field formulation ]；文中编号：未显式编号。
- **公式完整表达**：
[
\left\{ \begin{array}{ll} r_{\pmb{u}}(\pmb{u},d) = (\nabla^{(s)}\delta \pmb{u},\pmb{\sigma}) - (\delta \pmb{u},\pmb{b}) - (\delta \pmb{u},\pmb{t})_{\Gamma_{t}} = 0, \\ r_{d}(\pmb{u},d) = (\delta d,\frac{g_{c}}{l} d) + (\nabla \delta d,g_{c}l\nabla d) + (\delta d,g^{\prime}(d)\psi^{+}(\pmb{\epsilon})) = 0, \end{array} \right.
]
- **公式含义**：
  - 由能量导出应力，是平衡方程残量的核心材料项。
  - 定义残量或梯度，残量趋零对应离散平衡。
  - 邻近语境（简要）：前句“Therefore, the weak form of the phase-field formulation is expressed as”；后句“subject to the phase-field irreversibility condition.”。
- **必要推导过程**：
  1. 对能量对应变求一阶导得到应力。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。
  - 与离散求解链相连：输出 [ \mathbf{r} ] 与 [ \mathbf{K}/\mathbf{B}_k ] 供迭代器使用。

### 公式 020

- **出现位置**：`explain.md` 第 172-174 行；章节：[ ### 2.2. Finite element discretization ]；文中编号：未显式编号。
- **公式完整表达**：
[
\pmb {u} = \pmb{N}_{u_{A}}\pmb{u}_{A}, \quad \text{and} \quad d = N_{d_{A}}d_{A}.
]
- **公式含义**：
  - 把连续场映射到离散自由度，建立有限元代数系统。
  - 邻近语境（简要）：前句“Let \(\mathbf{N}_{u_{A}}\) represent the vector-valued shape function for the displacement field \(\pmb{u}\) at node \(A\) and \(N_{d_{A}}\) represent the scal…”；后句“Correspondingly, the displacement variation and the phase-field variation can be expressed as”。
- **必要推导过程**：
  1. 以形函数插值代入弱式并组装，得到离散残量与矩阵块。
- **与其他公式的内在联系**：
  - 与离散求解链相连：输出 [ \mathbf{r} ] 与 [ \mathbf{K}/\mathbf{B}_k ] 供迭代器使用。

### 公式 021

- **出现位置**：`explain.md` 第 178-180 行；章节：[ ### 2.2. Finite element discretization ]；文中编号：[ (1) ]。
- **公式完整表达**：
[
\delta \pmb {u} = \pmb{N}_{u_{A}}\delta \pmb{u}_{A}, \quad \text{and} \quad \delta d = N_{d_{A}}\delta d_{A},
]
- **公式含义**：
  - 把连续场映射到离散自由度，建立有限元代数系统。
  - 邻近语境（简要）：前句“Correspondingly, the displacement variation and the phase-field variation can be expressed as”；后句“where the Einstein summation is used. Plug the above expressions into the total energy defined in Eq. (1), the objective function for the inequality-constraine…”。
- **必要推导过程**：
  1. 以形函数插值代入弱式并组装，得到离散残量与矩阵块。
- **与其他公式的内在联系**：
  - 与离散求解链相连：输出 [ \mathbf{r} ] 与 [ \mathbf{K}/\mathbf{B}_k ] 供迭代器使用。
  - 与文中编号公式互引：[ (1) ]。

### 公式 022

- **出现位置**：`explain.md` 第 184-190 行；章节：[ ### 2.2. Finite element discretization ]；文中编号：[ (1) ]、[ (8) ]。
- **公式完整表达**：
[
\begin{array}{rl} 
\Pi (\pmb {u}_{A},d_{A}) = &\int_{\Omega}\psi \left(\pmb{\epsilon}(\pmb {N}_{u_{A}}\pmb {u}_{A}),N_{d_{A}}d_{A}\right)\mathrm{d}\Omega \\ 
&+\int_{\Omega}\frac{g_{c}}{2l}\left((N_{d_{A}}d_{A})^{2} + l^{2}(\nabla N_{d_{A}}d_{A})\cdot (\nabla N_{d_{A}}d_{A})\right)\mathrm{d}\Omega \\ 
&-\int_{\Omega}\pmb {b}\cdot (\pmb {N}_{u_{A}}\pmb {u}_{A})\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot (\pmb {N}_{u_{A}}\pmb {u}_{A})\mathrm{d}\Gamma, 
\end{array} \quad (8)
]
- **公式含义**：
  - 定义总势能目标 [ \Pi ]，把内能、断裂能与外载统一到同一优化框架。
  - 把连续场映射到离散自由度，建立有限元代数系统。
  - 邻近语境（简要）：前句“where the Einstein summation is used. Plug the above expressions into the total energy defined in Eq. (1), the objective function for the inequality-constraine…”；后句“The gradient of the discretized energy functional is derived as”。
- **必要推导过程**：
  1. 从最小势能原理写出内能减外力势，并加入断裂表面能项。
  2. 以形函数插值代入弱式并组装，得到离散残量与矩阵块。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。
  - 与离散求解链相连：输出 [ \mathbf{r} ] 与 [ \mathbf{K}/\mathbf{B}_k ] 供迭代器使用。
  - 与文中编号公式互引：[ (1) ]、[ (8) ]。

### 公式 023

- **出现位置**：`explain.md` 第 194-200 行；章节：[ ### 2.2. Finite element discretization ]；文中编号：[ (9) ]。
- **公式完整表达**：
[
\begin{array}{rl} 
\pmb{r} = \nabla \Pi = (r_{\pmb{u}},r_{d})^{\mathrm{T}},\\ 
r_{\pmb{u}_{A}} = \frac{\partial\Pi}{\partial \pmb{u}_{A}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},\pmb{\sigma}\right) - \left(\pmb{N}_{u_{A}},\pmb{b}\right) - \left(\pmb{N}_{u_{A}},\pmb{t}\right)_{\Gamma_{t}},\\ 
r_{d_{A}} = \frac{\partial\Pi}{\partial d_{A}} = \left(N_{d_{A}},\frac{g_{c}}{l} d + g^{\prime}(d)\psi^{+}\right) + (\nabla N_{d_{A}},g_{c}l\nabla d). 
\end{array} \quad (9)
]
- **公式含义**：
  - 由能量导出应力，是平衡方程残量的核心材料项。
  - 把连续场映射到离散自由度，建立有限元代数系统。
  - 定义残量或梯度，残量趋零对应离散平衡。
  - 邻近语境（简要）：前句“The gradient of the discretized energy functional is derived as”；后句“The Hessian matrix of the total energy after the finite element discretization is”。
- **必要推导过程**：
  1. 对能量对应变求一阶导得到应力。
  2. 以形函数插值代入弱式并组装，得到离散残量与矩阵块。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。
  - 与离散求解链相连：输出 [ \mathbf{r} ] 与 [ \mathbf{K}/\mathbf{B}_k ] 供迭代器使用。
  - 与文中编号公式互引：[ (9) ]。

### 公式 024

- **出现位置**：`explain.md` 第 204-206 行；章节：[ ### 2.2. Finite element discretization ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{K} = \nabla^{2}\Pi = \left[ \begin{array}{ll}\mathbf{K}_{uu} & \mathbf{K}_{ud}\\ \mathbf{K}_{du} & \mathbf{K}_{dd} \end{array} \right],
]
- **公式含义**：
  - 给出一致切线/雅可比信息，影响牛顿或准牛顿迭代收敛效率。
  - 邻近语境（简要）：前句“The Hessian matrix of the total energy after the finite element discretization is”；后句“where”。
- **必要推导过程**：
  1. 再对应力/残量求导得到一致切线矩阵或投影算子导数。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。

### 公式 025

- **出现位置**：`explain.md` 第 210-215 行；章节：[ ### 2.2. Finite element discretization ]；文中编号：[ (10) ]。
- **公式完整表达**：
[
\begin{array}{rl} 
\mathbf{K}_{u_{A}u_{B}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},\frac{\partial \pmb{\sigma}}{\partial \pmb{\epsilon}}:\nabla^{(s)}\pmb{N}_{u_{B}}\right),\qquad \mathbf{K}_{u_{A}d_{B}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},g^{\prime}(d)\pmb{\sigma}^{+}N_{d_{B}}\right),\\ 
\mathbf{K}_{d_{A}u_{B}} = \left(N_{d_{A}},g^{\prime}(d)\pmb{\sigma}^{+}:\nabla^{(s)}\pmb{N}_{u_{B}}\right),\qquad \mathbf{K}_{d_{A}d_{B}} = \left(N_{d_{A}},\left(\frac{g_{c}}{l} +g^{\prime \prime}(d)\psi^{+}\right)N_{d_{B}}\right) + \left(\nabla N_{d_{A}},g_{c}l\nabla N_{d_{B}}\right). 
\end{array} \quad (10)
]
- **公式含义**：
  - 给出一致切线/雅可比信息，影响牛顿或准牛顿迭代收敛效率。
  - 把连续场映射到离散自由度，建立有限元代数系统。
  - 邻近语境（简要）：前句“where”；后句“In the above block matrices, the subscripts \(A\) and \(B\) represent the finite element nodal number, respectively. During a typical time (load) step \([t_n,t…”。
- **必要推导过程**：
  1. 对能量对应变求一阶导得到应力。
  2. 再对应力/残量求导得到一致切线矩阵或投影算子导数。
  3. 以形函数插值代入弱式并组装，得到离散残量与矩阵块。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。
  - 与离散求解链相连：输出 [ \mathbf{r} ] 与 [ \mathbf{K}/\mathbf{B}_k ] 供迭代器使用。
  - 与文中编号公式互引：[ (10) ]。

### 公式 026

- **出现位置**：`explain.md` 第 219-221 行；章节：[ ### 2.2. Finite element discretization ]；文中编号：未显式编号。
- **公式完整表达**：
[
d_A^{(n)}\leq d_A\leq 1,
]
- **公式含义**：
  - 用于盒约束可行性维护，防止相场变量违反不可逆与上界约束。
  - 邻近语境（简要）：前句“In the above block matrices, the subscripts \(A\) and \(B\) represent the finite element nodal number, respectively. During a typical time (load) step \([t_n,t…”；后句“where the subscript \(A\) represents the finite element nodal number, and \(d_A^{(n)}\) is known at the beginning of the current time step. Lastly, a diagonal …”。
- **必要推导过程**：
  1. 按分量截断规则实施投影：低于下界取下界，高于上界取上界。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 027

- **出现位置**：`explain.md` 第 225-227 行；章节：[ ### 2.2. Finite element discretization ]；文中编号：[ (11) ]。
- **公式完整表达**：
[
\hat{\mathbf{K}} = \left[ \begin{array}{cc}\mathbf{K}_{uu} & \mathbf{0}\\ \mathbf{0} & \mathbf{K}_{dd} \end{array} \right], \quad (11)
]
- **公式含义**：
  - 给出一致切线/雅可比信息，影响牛顿或准牛顿迭代收敛效率。
  - 邻近语境（简要）：前句“where the subscript \(A\) represents the finite element nodal number, and \(d_A^{(n)}\) is known at the beginning of the current time step. Lastly, a diagonal …”；后句“which is needed in the monolithic scheme presented in Section 3.”。
- **必要推导过程**：
  1. 再对应力/残量求导得到一致切线矩阵或投影算子导数。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。
  - 与文中编号公式互引：[ (11) ]。

### 公式 028

- **出现位置**：`explain.md` 第 239-241 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (12) ]。
- **公式完整表达**：
[
\left(\pmb{u}_A,d_A\right) = \arg \min \Pi (\pmb {u}_A,d_A) \quad (12)
]
- **公式含义**：
  - 定义总势能目标 [ \Pi ]，把内能、断裂能与外载统一到同一优化框架。
  - 将求解过程明确为带约束极小化问题。
  - 邻近语境（简要）：前句“For the quasi-static phase-field fracture propagation problem, the pseudo time \(t\) is used to represent the actual load step. Inside the time (load) step \([…”；后句“subject to the inequality constraints”。
- **必要推导过程**：
  1. 从最小势能原理写出内能减外力势，并加入断裂表面能项。
  2. 把增量步未知量作为优化变量，配合约束写成极小化问题。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。
  - 与文中编号公式互引：[ (12) ]。

### 公式 029

- **出现位置**：`explain.md` 第 245-247 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (13) ]。
- **公式完整表达**：
[
d_A^{(n)}\leq d_A\leq 1 \quad (13)
]
- **公式含义**：
  - 用于盒约束可行性维护，防止相场变量违反不可逆与上界约束。
  - 邻近语境（简要）：前句“subject to the inequality constraints”；后句“and the linear constraints in the following form”。
- **必要推导过程**：
  1. 按分量截断规则实施投影：低于下界取下界，高于上界取上界。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与文中编号公式互引：[ (13) ]。

### 公式 030

- **出现位置**：`explain.md` 第 251-253 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (14) ]、[ (8) ]、[ (9) ]。
- **公式完整表达**：
[
\left[\pmb {u}_A,d_A\right]^{\mathrm{T}} = \mathbf{C}\left[\pmb {u}_A,d_A\right]^{\mathrm{T}} + \pmb {k}. \quad (14)
]
- **公式含义**：
  - 该式承担局部定义或过渡作用，连接前后推导步骤。
  - 邻近语境（简要）：前句“and the linear constraints in the following form”；后句“In the above constrained minimization problem, \(\{\pmb {u}_A,d_A\}\) represents the finite element nodal solution of the displacement field and the phase-fiel…”。
- **必要推导过程**：
  1. 由上一式变量定义代入并做代数整理可得本式。
- **与其他公式的内在联系**：
  - 与文中编号公式互引：[ (14) ]、[ (8) ]、[ (9) ]。

### 公式 031

- **出现位置**：`explain.md` 第 259-261 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (12) ]、[ (15) ]。
- **公式完整表达**：
[
m_{k}(\pmb{x}) = \Pi_{k} + \pmb{r}_{k}^{\mathrm{T}}(\pmb{x} - \pmb{x}_{k}) + \frac{1}{2} (\pmb{x} - \pmb{x}_{k})^{\mathrm{T}}\mathbf{B}_{k}(\pmb{x} - \pmb{x}_{k}) \quad (15)
]
- **公式含义**：
  - 该式承担局部定义或过渡作用，连接前后推导步骤。
  - 邻近语境（简要）：前句“To overcome the convergence difficulties arising from the non-convexity of the objective function in Eq. (12), the BFGS method is adopted. As a type of quasi-N…”；后句“where the subscript \(k\) represents the \(k\) th iteration, \(\mathbf{x}\) represents the unknowns \(\{\pmb{u},d\}\), and \(\pmb{r}_k = \nabla \Pi_k\) represe…”。
- **必要推导过程**：
  1. 由上一式变量定义代入并做代数整理可得本式。
- **与其他公式的内在联系**：
  - 与文中编号公式互引：[ (12) ]、[ (15) ]。

### 公式 032

- **出现位置**：`explain.md` 第 267-269 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (15) ]、[ (13) ]、[ (16) ]。
- **公式完整表达**：
[
\mathrm{Proj}_c(x_i,\mathrm{lb}_i,\mathrm{ub}_i) = \left\{ \begin{array}{ll}\mathrm{lb}_i & \mathrm{if}\; x_i < \mathrm{lb}_i,\\ x_i & \mathrm{if}\; \mathrm{lb}_i \leq x_i \leq \mathrm{ub}_i,\\ \mathrm{ub}_i & \mathrm{if}\; x_i > \mathrm{ub}_i, \end{array} \right. \quad (16)
]
- **公式含义**：
  - 用于盒约束可行性维护，防止相场变量违反不可逆与上界约束。
  - 邻近语境（简要）：前句“In order to consider the phase-field irreversibility condition, the quadratic model in Eq. (15) needs to satisfy the box constraints shown in Eq. (13). Let \(\…”；后句“where the feasible region \(C\) formed by the phase-field inequality constraints is obviously convex. Based on the above projection operator, at the \(k\) th i…”。
- **必要推导过程**：
  1. 按分量截断规则实施投影：低于下界取下界，高于上界取上界。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与文中编号公式互引：[ (15) ]、[ (13) ]、[ (16) ]。

### 公式 033

- **出现位置**：`explain.md` 第 273-275 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (17) ]。
- **公式完整表达**：
[
\pmb {x}(t) = \mathrm{Proj}_c(\pmb {x}_k - t\pmb {r}_k,\mathbf{lb},\mathbf{ub}),\quad t > 0, \quad (17)
]
- **公式含义**：
  - 用于盒约束可行性维护，防止相场变量违反不可逆与上界约束。
  - 邻近语境（简要）：前句“where the feasible region \(C\) formed by the phase-field inequality constraints is obviously convex. Based on the above projection operator, at the \(k\) th i…”；后句“in which the projection operator is applied to the vector \(\pmb{x}_k - t\pmb{r}_k\) component wise.”。
- **必要推导过程**：
  1. 按分量截断规则实施投影：低于下界取下界，高于上界取上界。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与文中编号公式互引：[ (17) ]。

### 公式 034

- **出现位置**：`explain.md` 第 281-283 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (15) ]、[ (18) ]。
- **公式完整表达**：
[
p_k(t) = m_k(\pmb{x}(t)) = \Pi_k + \pmb{r}_k^\mathrm{T}(\pmb{x}(t) - \pmb{x}_k) + \frac{1}{2} (\pmb{x}(t) - \pmb{x}_k)^\mathrm{T}\mathbf{B}_k(\pmb{x}(t) - \pmb{x}_k). \quad (18)
]
- **公式含义**：
  - 该式承担局部定义或过渡作用，连接前后推导步骤。
  - 邻近语境（简要）：前句“Following the piecewise linear path formed by the projected gradient, the quadratic model in the BFGS method, as shown in Eq. (15), is transformed into a univa…”；后句“The first local minimizer \(\pmb{x}^c\) of the above model is called the generalized Cauchy point, the calculation of which is presented in Section 3.3.”。
- **必要推导过程**：
  1. 由上一式变量定义代入并做代数整理可得本式。
- **与其他公式的内在联系**：
  - 与文中编号公式互引：[ (15) ]、[ (18) ]。

### 公式 035

- **出现位置**：`explain.md` 第 289-291 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (13) ]、[ (15) ]。
- **公式完整表达**：
[
\mathcal{A}(\pmb {x}) = \{i:x_i = \mathrm{lb}_i\} \cup \{i:x_i = \mathrm{ub}_i\}.
]
- **公式含义**：
  - 用于盒约束可行性维护，防止相场变量违反不可逆与上界约束。
  - 定义活动集，用于后续自由子空间构造。
  - 邻近语境（简要）：前句“For an arbitrary solution state \(\pmb{x}\) the active set \(\mathcal{A}(\pmb {x})\) of the inequality constraints shown in Eq. (13) include those components w…”；后句“Based on the generalized Cauchy point \(\pmb{x}^c\), the variables in the active set \(\mathcal{A}_k(\pmb{x}^c)\) (located at the boundary of the box constrain…”。
- **必要推导过程**：
  1. 按分量截断规则实施投影：低于下界取下界，高于上界取上界。
  2. 根据触边条件更新活动集，并在其补集上继续优化。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与文中编号公式互引：[ (13) ]、[ (15) ]。

### 公式 036

- **出现位置**：`explain.md` 第 295-297 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (15) ]。
- **公式完整表达**：
[
\pmb{x}^* = \arg \min m_k(\pmb {x})
]
- **公式含义**：
  - 将求解过程明确为带约束极小化问题。
  - 邻近语境（简要）：前句“Based on the generalized Cauchy point \(\pmb{x}^c\), the variables in the active set \(\mathcal{A}_k(\pmb{x}^c)\) (located at the boundary of the box constrain…”；后句“subject to”。
- **必要推导过程**：
  1. 把增量步未知量作为优化变量，配合约束写成极小化问题。
- **与其他公式的内在联系**：
  - 与文中编号公式互引：[ (15) ]。

### 公式 037

- **出现位置**：`explain.md` 第 301-303 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：未显式编号。
- **公式完整表达**：
[
x_{i} = x_{i}^{c},\;\forall i\in \mathcal{A}_{k}(\pmb{x}^{c})\quad \text{and}\quad \mathrm{lb}_{i}\leq x_{i}\leq \mathrm{ub}_{i},\;\forall i\notin \mathcal{A}_{k}(\pmb{x}^{c}).
]
- **公式含义**：
  - 用于盒约束可行性维护，防止相场变量违反不可逆与上界约束。
  - 定义活动集，用于后续自由子空间构造。
  - 邻近语境（简要）：前句“subject to”；后句“Firstly ignoring the inequality constraints on the free variables, the above minimization can be solved using either the primal approach on the subspace of the…”。
- **必要推导过程**：
  1. 按分量截断规则实施投影：低于下界取下界，高于上界取上界。
  2. 根据触边条件更新活动集，并在其补集上继续优化。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 038

- **出现位置**：`explain.md` 第 309-311 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：未显式编号。
- **公式完整表达**：
[
\pmb{p}_k = \pmb{x}^* - \pmb{x}_k.
]
- **公式含义**：
  - 该式承担局部定义或过渡作用，连接前后推导步骤。
  - 邻近语境（简要）：前句“Based on the solution \(\pmb{x}^*\) of the above minimization, a search direction \(\pmb{p}_k\) is defined as”；后句“The search direction \(\pmb{p}_k\) is a descent direction of the objective function \(\Pi (\pmb{x})\), see [41] for the proof. The new solution \(\pmb{x}_{k + …”。
- **必要推导过程**：
  1. 由上一式变量定义代入并做代数整理可得本式。
- **与其他公式的内在联系**：
  - 在局部推导中承前启后，为下一式提供变量或算子。

### 公式 039

- **出现位置**：`explain.md` 第 315-317 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：未显式编号。
- **公式完整表达**：
[
\pmb{x}_{k + 1} = \pmb{x}_k + \alpha_k\pmb{p}_k.
]
- **公式含义**：
  - 定义线搜索更新与步长准则。
  - 邻近语境（简要）：前句“The search direction \(\pmb{p}_k\) is a descent direction of the objective function \(\Pi (\pmb{x})\), see [41] for the proof. The new solution \(\pmb{x}_{k + …”；后句“The positive step length \(\alpha_k > 0\) is determined using the line search based on the strong Wolfe conditions, which include the sufficient decrease condi…”。
- **必要推导过程**：
  1. 给定方向后使用 Armijo/Wolfe 条件搜索 [ \alpha_k ]。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 040

- **出现位置**：`explain.md` 第 321-323 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：未显式编号。
- **公式完整表达**：
[
\Pi_{k + 1} = \Pi (\pmb{x}_k + \alpha_k\pmb {p}_k) \leq \Pi_k + c_1\alpha_k\pmb {p}_k^\mathrm{T}\pmb{r}_k
]
- **公式含义**：
  - 定义总势能目标 [ \Pi ]，把内能、断裂能与外载统一到同一优化框架。
  - 用于盒约束可行性维护，防止相场变量违反不可逆与上界约束。
  - 定义线搜索更新与步长准则。
  - 邻近语境（简要）：前句“The positive step length \(\alpha_k > 0\) is determined using the line search based on the strong Wolfe conditions, which include the sufficient decrease condi…”；后句“and the curvature condition”。
- **必要推导过程**：
  1. 从最小势能原理写出内能减外力势，并加入断裂表面能项。
  2. 按分量截断规则实施投影：低于下界取下界，高于上界取上界。
  3. 给定方向后使用 Armijo/Wolfe 条件搜索 [ \alpha_k ]。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 041

- **出现位置**：`explain.md` 第 327-329 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：未显式编号。
- **公式完整表达**：
[
\left|\pmb{r}_{k + 1}^\mathrm{T}\pmb {p}_k\right| = \left|\pmb {r}(\pmb{x}_k + \alpha_k\pmb {p}_k)^\mathrm{T}\pmb {p}_k\right| \leq c_2\left|\pmb{r}_k^\mathrm{T}\pmb {p}_k\right|.
]
- **公式含义**：
  - 用于盒约束可行性维护，防止相场变量违反不可逆与上界约束。
  - 定义线搜索更新与步长准则。
  - 邻近语境（简要）：前句“and the curvature condition”；后句“To find a step length \(\alpha_k\) satisfying the strong Wolfe conditions, this work uses the implementation proposed by More and Thuente [52] with the paramet…”。
- **必要推导过程**：
  1. 按分量截断规则实施投影：低于下界取下界，高于上界取上界。
  2. 给定方向后使用 Armijo/Wolfe 条件搜索 [ \alpha_k ]。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 042

- **出现位置**：`explain.md` 第 336-338 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathcal{A}_{k + 1}(\pmb{x}^c) = \mathcal{A}_k(\pmb{x}^c).
]
- **公式含义**：
  - 定义活动集，用于后续自由子空间构造。
  - 邻近语境（简要）：前句“1. The active set remains unchanged between two consecutive iterations:”；后句“2. The \(l_{2}\)-norm of the projected gradient is smaller than the prescribed tolerance:”。
- **必要推导过程**：
  1. 根据触边条件更新活动集，并在其补集上继续优化。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 043

- **出现位置**：`explain.md` 第 340-342 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：未显式编号。
- **公式完整表达**：
[
\| \mathrm{Proj}_C(\mathbf{x}_{k + 1} - \mathbf{r}_{k + 1},\mathbf{lb},\mathbf{ub}) - \mathbf{x}_{k + 1}\|_{2} < \mathrm{tol}.
]
- **公式含义**：
  - 用于盒约束可行性维护，防止相场变量违反不可逆与上界约束。
  - 定义残量或梯度，残量趋零对应离散平衡。
  - 邻近语境（简要）：前句“2. The \(l_{2}\)-norm of the projected gradient is smaller than the prescribed tolerance:”；后句“3. The \(l_{2}\)-norm of the solution increment is smaller than the prescribed tolerance:”。
- **必要推导过程**：
  1. 按分量截断规则实施投影：低于下界取下界，高于上界取上界。
- **与其他公式的内在联系**：
  - 与离散求解链相连：输出 [ \mathbf{r} ] 与 [ \mathbf{K}/\mathbf{B}_k ] 供迭代器使用。
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 044

- **出现位置**：`explain.md` 第 344-346 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (16) ]。
- **公式完整表达**：
[
\| \Delta \mathbf{x}_{k + 1}\|_{2} = \| \mathbf{x}_{k + 1} - \mathbf{x}_{k}\|_{2} < \mathrm{tol}.
]
- **公式含义**：
  - 该式承担局部定义或过渡作用，连接前后推导步骤。
  - 邻近语境（简要）：前句“3. The \(l_{2}\)-norm of the solution increment is smaller than the prescribed tolerance:”；后句“Recall that the projection operator is defined in Eq. (16). For the \(i\) th component of the projected gradient,”。
- **必要推导过程**：
  1. 由上一式变量定义代入并做代数整理可得本式。
- **与其他公式的内在联系**：
  - 与文中编号公式互引：[ (16) ]。

### 公式 045

- **出现位置**：`explain.md` 第 350-352 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (16) ]、[ (19) ]。
- **公式完整表达**：
[
\mathrm{Proj}_C(x_i - r_i,\mathrm{lb}_i,\mathrm{ub}_i) - x_i = \left\{ \begin{array}{ll}\mathrm{lb}_i - x_i & \mathrm{if}\; x_i - r_i < \mathrm{lb}_i,\\ - r_i & \mathrm{if}\; \mathrm{lb}_i \leq x_i - r_i \leq \mathrm{ub}_i,\\ \mathrm{ub}_i - x_i & \mathrm{if}\; x_i - r_i > \mathrm{ub}_i. \end{array} \right. \quad (19)
]
- **公式含义**：
  - 用于盒约束可行性维护，防止相场变量违反不可逆与上界约束。
  - 邻近语境（简要）：前句“Recall that the projection operator is defined in Eq. (16). For the \(i\) th component of the projected gradient,”；后句“The second convergence criterion based on the projected gradient means that if the \(i\) th component of the projected gradient is within the feasible region (…”。
- **必要推导过程**：
  1. 按分量截断规则实施投影：低于下界取下界，高于上界取上界。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与文中编号公式互引：[ (16) ]、[ (19) ]。

### 公式 046

- **出现位置**：`explain.md` 第 358-360 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (14) ]。
- **公式完整表达**：
[
\mathbf{x} = \mathbf{C}\mathbf{x} + \mathbf{k}.
]
- **公式含义**：
  - 该式承担局部定义或过渡作用，连接前后推导步骤。
  - 邻近语境（简要）：前句“**Comment 1.** The linear constraints shown in Eq. (14) include the essential (Dirichlet) boundary conditions imposed on the domain boundary and the hanging-no…”；后句“For a linear system \(\mathbf{A}\mathbf{x} = \mathbf{b}\) with the above set of constraints, we can instead solve the following modified linear system [53]”。
- **必要推导过程**：
  1. 由上一式变量定义代入并做代数整理可得本式。
- **与其他公式的内在联系**：
  - 与文中编号公式互引：[ (14) ]。

### 公式 047

- **出现位置**：`explain.md` 第 364-366 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：未显式编号。
- **公式完整表达**：
[
\left(\mathbf{C}^{\mathrm{T}}\mathbf{A}\mathbf{C} + \mathbf{I}_{d_{c}}\right)\hat{\mathbf{x}} = \mathbf{C}^{\mathrm{T}}(\mathbf{b} - \mathbf{A}\mathbf{k})
]
- **公式含义**：
  - 该式承担局部定义或过渡作用，连接前后推导步骤。
  - 邻近语境（简要）：前句“For a linear system \(\mathbf{A}\mathbf{x} = \mathbf{b}\) with the above set of constraints, we can instead solve the following modified linear system [53]”；后句“and then recover the true solution \(\mathbf{x}\) as”。
- **必要推导过程**：
  1. 由上一式变量定义代入并做代数整理可得本式。
- **与其他公式的内在联系**：
  - 在局部推导中承前启后，为下一式提供变量或算子。

### 公式 048

- **出现位置**：`explain.md` 第 370-372 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{x} = \mathbf{C}\hat{\mathbf{x}} + \mathbf{k}.
]
- **公式含义**：
  - 该式承担局部定义或过渡作用，连接前后推导步骤。
  - 邻近语境（简要）：前句“and then recover the true solution \(\mathbf{x}\) as”；后句“In the modified linear system, the matrix \(\mathbf{I}_{d_c}\) is defined as”。
- **必要推导过程**：
  1. 由上一式变量定义代入并做代数整理可得本式。
- **与其他公式的内在联系**：
  - 在局部推导中承前启后，为下一式提供变量或算子。

### 公式 049

- **出现位置**：`explain.md` 第 375-377 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：未显式编号。
- **公式完整表达**：
[
(\mathbf{I}_{d_c})_{ii} = \begin{cases} 1 & \text{if DoF } i \text{ is constrained}, \\ 0 & \text{otherwise}, \end{cases}
]
- **公式含义**：
  - 该式承担局部定义或过渡作用，连接前后推导步骤。
  - 邻近语境（简要）：前句“In the modified linear system, the matrix \(\mathbf{I}_{d_c}\) is defined as”；后句“where \(\mathcal{T}\) represents the set of the degrees of freedom at the constrained nodes, including the hanging nodes and the nodes prescribed with essentia…”。
- **必要推导过程**：
  1. 由上一式变量定义代入并做代数整理可得本式。
- **与其他公式的内在联系**：
  - 在局部推导中承前启后，为下一式提供变量或算子。

### 公式 050

- **出现位置**：`explain.md` 第 384-386 行；章节：[ ### 3.2. Compact representation of limited-memory BFGS matrix ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{s}_k = \mathbf{x}_{k + 1} - \mathbf{x}_k,\quad \mathbf{y}_k = \mathbf{r}_{k + 1} - \mathbf{r}_k,
]
- **公式含义**：
  - 定义残量或梯度，残量趋零对应离散平衡。
  - 属于 L-BFGS 紧凑更新链，利用历史向量近似曲率。
  - 邻近语境（简要）：前句“The detailed approach of the standard BFGS method can be found in many textbooks in the field of numerical optimization, for instance, see the classical textbo…”；后句“and the BFGS matrix is updated as”。
- **必要推导过程**：
  1. 由 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对构造紧凑低秩更新。
- **与其他公式的内在联系**：
  - 与离散求解链相连：输出 [ \mathbf{r} ] 与 [ \mathbf{K}/\mathbf{B}_k ] 供迭代器使用。
  - 与低秩线代链相连：L-BFGS 紧凑更新与 SMW 逆更新互补。

### 公式 051

- **出现位置**：`explain.md` 第 390-392 行；章节：[ ### 3.2. Compact representation of limited-memory BFGS matrix ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{B}_{k + 1} = \mathbf{B}_k - \frac{\mathbf{B}_k\mathbf{s}_k\mathbf{s}_k^{\mathrm{T}}\mathbf{B}_k}{\mathbf{s}_k^{\mathrm{T}}\mathbf{B}_k\mathbf{s}_k} + \frac{\mathbf{y}_k\mathbf{y}_k^{\mathrm{T}}}{\mathbf{y}_k^{\mathrm{T}}\mathbf{s}_k}.
]
- **公式含义**：
  - 属于 L-BFGS 紧凑更新链，利用历史向量近似曲率。
  - 邻近语境（简要）：前句“and the BFGS matrix is updated as”；后句“If the old BFGS matrix \(\mathbf{B}_k\) is positive definite and the vector pairs \(\mathbf{s}_k,\mathbf{y}_k\) satisfy the curvature condition,”。
- **必要推导过程**：
  1. 由 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对构造紧凑低秩更新。
- **与其他公式的内在联系**：
  - 与低秩线代链相连：L-BFGS 紧凑更新与 SMW 逆更新互补。

### 公式 052

- **出现位置**：`explain.md` 第 396-398 行；章节：[ ### 3.2. Compact representation of limited-memory BFGS matrix ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{s}_k^{\mathrm{T}}\mathbf{y}_k > 0,
]
- **公式含义**：
  - 属于 L-BFGS 紧凑更新链，利用历史向量近似曲率。
  - 邻近语境（简要）：前句“If the old BFGS matrix \(\mathbf{B}_k\) is positive definite and the vector pairs \(\mathbf{s}_k,\mathbf{y}_k\) satisfy the curvature condition,”；后句“then the updated BFGS matrix \(\mathbf{B}_{k + 1}\) is also positive definite [45]. Since \(\mathbf{s}_k\) and \(\mathbf{y}_k\) are two vectors with \(n\) comp…”。
- **必要推导过程**：
  1. 由 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对构造紧凑低秩更新。
- **与其他公式的内在联系**：
  - 与低秩线代链相连：L-BFGS 紧凑更新与 SMW 逆更新互补。

### 公式 053

- **出现位置**：`explain.md` 第 404-406 行；章节：[ ### 3.2. Compact representation of limited-memory BFGS matrix ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{S}_k = [\mathbf{s}_{k-m} \cdots \mathbf{s}_{k-1}]
]
- **公式含义**：
  - 该式承担局部定义或过渡作用，连接前后推导步骤。
  - 邻近语境（简要）：前句“In this work, the compact representation of the limited-memory BFGS (L-BFGS) matrix, originally proposed by Byrd et al. [29], is adopted to avoid the storage o…”；后句“and”。
- **必要推导过程**：
  1. 由上一式变量定义代入并做代数整理可得本式。
- **与其他公式的内在联系**：
  - 在局部推导中承前启后，为下一式提供变量或算子。

### 公式 054

- **出现位置**：`explain.md` 第 408-410 行；章节：[ ### 3.2. Compact representation of limited-memory BFGS matrix ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{Y}_k = [\mathbf{y}_{k-m} \cdots \mathbf{y}_{k-1}].
]
- **公式含义**：
  - 该式承担局部定义或过渡作用，连接前后推导步骤。
  - 邻近语境（简要）：前句“and”；后句“Then, the L-BFGS matrix can be represented in the following compact form,”。
- **必要推导过程**：
  1. 由上一式变量定义代入并做代数整理可得本式。
- **与其他公式的内在联系**：
  - 在局部推导中承前启后，为下一式提供变量或算子。

### 公式 055

- **出现位置**：`explain.md` 第 414-416 行；章节：[ ### 3.2. Compact representation of limited-memory BFGS matrix ]；文中编号：[ (20) ]。
- **公式完整表达**：
[
\mathbf{B}_k = \mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k, \quad (20)
]
- **公式含义**：
  - 属于 L-BFGS 紧凑更新链，利用历史向量近似曲率。
  - 邻近语境（简要）：前句“Then, the L-BFGS matrix can be represented in the following compact form,”；后句“where”。
- **必要推导过程**：
  1. 由 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对构造紧凑低秩更新。
- **与其他公式的内在联系**：
  - 与低秩线代链相连：L-BFGS 紧凑更新与 SMW 逆更新互补。
  - 与文中编号公式互引：[ (20) ]。

### 公式 056

- **出现位置**：`explain.md` 第 420-422 行；章节：[ ### 3.2. Compact representation of limited-memory BFGS matrix ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{W}_k = [\mathbf{Y}_k \quad \mathbf{B}^0_k \mathbf{S}_k] \in \mathbb{R}^{n\times 2m}
]
- **公式含义**：
  - 属于 L-BFGS 紧凑更新链，利用历史向量近似曲率。
  - 邻近语境（简要）：前句“where”；后句“and”。
- **必要推导过程**：
  1. 由 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对构造紧凑低秩更新。
- **与其他公式的内在联系**：
  - 与低秩线代链相连：L-BFGS 紧凑更新与 SMW 逆更新互补。

### 公式 057

- **出现位置**：`explain.md` 第 425-427 行；章节：[ ### 3.2. Compact representation of limited-memory BFGS matrix ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{M}_k = \begin{bmatrix} -\mathbf{D}_k & \mathbf{L}^{\mathrm{T}}_k \\ \mathbf{L}_k & \mathbf{S}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{S}_k \end{bmatrix}^{-1} \in \mathbb{R}^{2m\times 2m}.
]
- **公式含义**：
  - 属于 L-BFGS 紧凑更新链，利用历史向量近似曲率。
  - 邻近语境（简要）：前句“and”；后句“The matrix \(\mathbf{M}_k\) is a \(2m\times 2m\) matrix with the block matrices defined as”。
- **必要推导过程**：
  1. 由 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对构造紧凑低秩更新。
- **与其他公式的内在联系**：
  - 与低秩线代链相连：L-BFGS 紧凑更新与 SMW 逆更新互补。

### 公式 058

- **出现位置**：`explain.md` 第 431-433 行；章节：[ ### 3.2. Compact representation of limited-memory BFGS matrix ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{D}_k = \mathrm{diag}\{\mathbf{s}^{\mathrm{T}}_{k-m}\mathbf{y}_{k-m}, \ldots, \mathbf{s}^{\mathrm{T}}_{k-1}\mathbf{y}_{k-1}\} \in \mathbb{R}^{m\times m}
]
- **公式含义**：
  - 属于 L-BFGS 紧凑更新链，利用历史向量近似曲率。
  - 邻近语境（简要）：前句“The matrix \(\mathbf{M}_k\) is a \(2m\times 2m\) matrix with the block matrices defined as”；后句“and”。
- **必要推导过程**：
  1. 由 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对构造紧凑低秩更新。
- **与其他公式的内在联系**：
  - 与低秩线代链相连：L-BFGS 紧凑更新与 SMW 逆更新互补。

### 公式 059

- **出现位置**：`explain.md` 第 435-437 行；章节：[ ### 3.2. Compact representation of limited-memory BFGS matrix ]；文中编号：[ (11) ]。
- **公式完整表达**：
[
\mathbf{L}_k \in \mathbb{R}^{m\times m}, \quad (\mathbf{L}_k)_{i,j} = \begin{cases} \mathbf{s}^{\mathrm{T}}_{i+k-m-1}\mathbf{y}_{j+k-m-1} & \text{if } i > j, \\ 0 & \text{if } i \leq j. \end{cases}
]
- **公式含义**：
  - 用于盒约束可行性维护，防止相场变量违反不可逆与上界约束。
  - 属于 L-BFGS 紧凑更新链，利用历史向量近似曲率。
  - 邻近语境（简要）：前句“and”；后句“At the \(k\) th iteration of the phase-field monolithic scheme, the initial BFGS matrix \(\mathbf{B}^0_k\) takes the diagonal block matrix defined in Eq. (11),…”。
- **必要推导过程**：
  1. 按分量截断规则实施投影：低于下界取下界，高于上界取上界。
  2. 由 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对构造紧凑低秩更新。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与低秩线代链相连：L-BFGS 紧凑更新与 SMW 逆更新互补。
  - 与文中编号公式互引：[ (11) ]。

### 公式 060

- **出现位置**：`explain.md` 第 441-443 行；章节：[ ### 3.2. Compact representation of limited-memory BFGS matrix ]；文中编号：[ (11) ]、[ (21) ]、[ (20) ]。
- **公式完整表达**：
[
\mathbf{B}^0_k = \hat{\mathbf{K}}^{(k)} = \begin{bmatrix} \mathbf{K}_{uu} & \mathbf{0} \\ \mathbf{0} & \mathbf{K}_{dd} \end{bmatrix}^{(k)}. \quad (21)
]
- **公式含义**：
  - 给出一致切线/雅可比信息，影响牛顿或准牛顿迭代收敛效率。
  - 邻近语境（简要）：前句“At the \(k\) th iteration of the phase-field monolithic scheme, the initial BFGS matrix \(\mathbf{B}^0_k\) takes the diagonal block matrix defined in Eq. (11),…”；后句“The matrices \(\mathbf{K}_{uu}\) and \(\mathbf{K}_{dd}\) are both assembled from the finite element procedure, and therefore, are sparse. Also, both \(\mathbf{…”。
- **必要推导过程**：
  1. 再对应力/残量求导得到一致切线矩阵或投影算子导数。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。
  - 与文中编号公式互引：[ (11) ]、[ (21) ]、[ (20) ]。

### 公式 061

- **出现位置**：`explain.md` 第 455-457 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：未显式编号。
- **公式完整表达**：
[
t_i = \begin{cases} (x^0_i - \mathrm{ub}_i)/r_i & \text{if } r_i < 0, \\ (x^0_i - \mathrm{lb}_i)/r_i & \text{if } r_i > 0, \\ +\infty & \text{if } r_i = 0. \end{cases}
]
- **公式含义**：
  - 用于盒约束可行性维护，防止相场变量违反不可逆与上界约束。
  - 属于广义 Cauchy 点计算链，分段推进并更新一维导数。
  - 邻近语境（简要）：前句“Once the L-BFGS matrix is formed using the compact representation, the next step is to locate the generalized Cauchy point, which is defined as the first local…”；后句“Note that if there is no upper or lower bound on the \(i\) th direction, then set \(\mathrm{ub}_i = +\infty\) or \(\mathrm{lb}_i = -\infty\). After all the bre…”。
- **必要推导过程**：
  1. 按分量截断规则实施投影：低于下界取下界，高于上界取上界。
  2. 计算断点 [ t_i ]，在分段二次模型上求段内极小或推进到下一断点。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 062

- **出现位置**：`explain.md` 第 461-463 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：未显式编号。
- **公式完整表达**：
[
x_i(t) = x_i^{0} - \min \{t,t_i\} r_i.
]
- **公式含义**：
  - 属于广义 Cauchy 点计算链，分段推进并更新一维导数。
  - 邻近语境（简要）：前句“Note that if there is no upper or lower bound on the \(i\) th direction, then set \(\mathrm{ub}_i = +\infty\) or \(\mathrm{lb}_i = -\infty\). After all the bre…”；后句“Inside an interval \([t^{(j-1)}, t^{(j)}]\) along the piecewise linear path formed by the projected gradient, the solution vector can be expressed as a functio…”。
- **必要推导过程**：
  1. 计算断点 [ t_i ]，在分段二次模型上求段内极小或推进到下一断点。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 063

- **出现位置**：`explain.md` 第 467-469 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{x}(t) = \mathbf{x}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)},
]
- **公式含义**：
  - 属于广义 Cauchy 点计算链，分段推进并更新一维导数。
  - 邻近语境（简要）：前句“Inside an interval \([t^{(j-1)}, t^{(j)}]\) along the piecewise linear path formed by the projected gradient, the solution vector can be expressed as a functio…”；后句“where \(\mathbf{x}^{(j-1)} = \mathbf{x}(t^{(j-1)})\), \(\Delta t = t - t^{(j-1)}\), and \(\mathbf{d}^{(j-1)}\) is the direction of the projected gradient in th…”。
- **必要推导过程**：
  1. 计算断点 [ t_i ]，在分段二次模型上求段内极小或推进到下一断点。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 064

- **出现位置**：`explain.md` 第 473-481 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：[ (22) ]。
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
- **公式含义**：
  - 定义残量或梯度，残量趋零对应离散平衡。
  - 属于广义 Cauchy 点计算链，分段推进并更新一维导数。
  - 邻近语境（简要）：前句“where \(\mathbf{x}^{(j-1)} = \mathbf{x}(t^{(j-1)})\), \(\Delta t = t - t^{(j-1)}\), and \(\mathbf{d}^{(j-1)}\) is the direction of the projected gradient in th…”；后句“where \(f_{j-1}, f_{j-1}^{\prime}\) and \(f_{j-1}^{\prime \prime}\) represent the coefficients of the quadratic polynomial \(\hat{p}(\Delta t)\).”。
- **必要推导过程**：
  1. 计算断点 [ t_i ]，在分段二次模型上求段内极小或推进到下一断点。
- **与其他公式的内在联系**：
  - 与离散求解链相连：输出 [ \mathbf{r} ] 与 [ \mathbf{K}/\mathbf{B}_k ] 供迭代器使用。
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与文中编号公式互引：[ (22) ]。

### 公式 065

- **出现位置**：`explain.md` 第 487-489 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：[ (18) ]。
- **公式完整表达**：
[
f_{j-1}^{\prime \prime} = \mathbf{d}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j-1)} > 0.
]
- **公式含义**：
  - 该式承担局部定义或过渡作用，连接前后推导步骤。
  - 邻近语境（简要）：前句“As discussed in Section 3.2, the updated L-BFGS matrix \(\mathbf{B}\) is positive definite if the vector-pairs \(\{\mathbf{s}_i,\mathbf{y}_i\}_{i=k-m}^{k-1}\) …”；后句“The quadratic polynomial \(\hat{p} (\Delta t)\) reaches to its minimum value at the critical point \(\Delta t^*\)”。
- **必要推导过程**：
  1. 由上一式变量定义代入并做代数整理可得本式。
- **与其他公式的内在联系**：
  - 与文中编号公式互引：[ (18) ]。

### 公式 066

- **出现位置**：`explain.md` 第 493-495 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：未显式编号。
- **公式完整表达**：
[
\hat{p}'(\Delta t) = 0 \Rightarrow \Delta t^* = -f_{j-1}' / f_{j-1}^{\prime \prime}.
]
- **公式含义**：
  - 属于广义 Cauchy 点计算链，分段推进并更新一维导数。
  - 邻近语境（简要）：前句“The quadratic polynomial \(\hat{p} (\Delta t)\) reaches to its minimum value at the critical point \(\Delta t^*\)”；后句“Depending on the value of \(\Delta t^*\), there are the following three scenarios:”。
- **必要推导过程**：
  1. 计算断点 [ t_i ]，在分段二次模型上求段内极小或推进到下一断点。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 067

- **出现位置**：`explain.md` 第 504-506 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathcal{A}(\mathbf{x}^c) = \{i : t_i \leq t^{(j-1)}, i = 1, \ldots, n\}.
]
- **公式含义**：
  - 用于盒约束可行性维护，防止相场变量违反不可逆与上界约束。
  - 定义活动集，用于后续自由子空间构造。
  - 属于广义 Cauchy 点计算链，分段推进并更新一维导数。
  - 邻近语境（简要）：前句“For the first two scenarios, the set of active constraints include those components whose break point \(t\)-value is no larger than \(t^{(j-1)}\), that is,”；后句“During the next interval \([t^{(j)}, t^{(j+1)}]\), it is not necessary to calculate the values of \(f'_j\) and \(f''_j\) from scratch. Rather, they can be obta…”。
- **必要推导过程**：
  1. 按分量截断规则实施投影：低于下界取下界，高于上界取上界。
  2. 根据触边条件更新活动集，并在其补集上继续优化。
  3. 计算断点 [ t_i ]，在分段二次模型上求段内极小或推进到下一断点。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 068

- **出现位置**：`explain.md` 第 509-511 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：未显式编号。
- **公式完整表达**：
[
\Delta t^{(j-1)} = t^{(j)} - t^{(j-1)}, \quad \mathbf{x}^{(j)} = \mathbf{x}^{(j-1)} + \Delta t^{(j-1)}\mathbf{d}^{(j-1)}.
]
- **公式含义**：
  - 属于广义 Cauchy 点计算链，分段推进并更新一维导数。
  - 邻近语境（简要）：前句“During the next interval \([t^{(j)}, t^{(j+1)}]\), it is not necessary to calculate the values of \(f'_j\) and \(f''_j\) from scratch. Rather, they can be obta…”；后句“Moving from the previous interval \([t^{(j-1)}, t^{(j)}]\) to the new interval \([t^{(j)}, t^{(j+1)}]\), assume that the constraint associated with the \(b\) t…”。
- **必要推导过程**：
  1. 计算断点 [ t_i ]，在分段二次模型上求段内极小或推进到下一断点。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 069

- **出现位置**：`explain.md` 第 514-516 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：[ (22) ]。
- **公式完整表达**：
[
\mathbf{d}^{(j)} = \mathbf{d}^{(j-1)} + r_b \mathbf{e}_b,
]
- **公式含义**：
  - 该式承担局部定义或过渡作用，连接前后推导步骤。
  - 邻近语境（简要）：前句“Moving from the previous interval \([t^{(j-1)}, t^{(j)}]\) to the new interval \([t^{(j)}, t^{(j+1)}]\), assume that the constraint associated with the \(b\) t…”；后句“where \(\mathbf{e}_b\) is the \(b\) th unit vector, and \(r_b\) is the \(b\) th component of the gradient vector \(\mathbf{r}\). Similar to Eq. (22), the new v…”。
- **必要推导过程**：
  1. 由上一式变量定义代入并做代数整理可得本式。
- **与其他公式的内在联系**：
  - 与文中编号公式互引：[ (22) ]。

### 公式 070

- **出现位置**：`explain.md` 第 519-521 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：[ (22) ]、[ (23) ]。
- **公式完整表达**：
[
f'_j = \mathbf{r}^{\mathrm{T}}\mathbf{d}^{(j)} + \mathbf{d}^{(j)}{}^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j)} = f'_{j-1} + \Delta t^{(j-1)} f''_{j-1} + r^2_b + r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j)} \quad (23)
]
- **公式含义**：
  - 定义残量或梯度，残量趋零对应离散平衡。
  - 属于广义 Cauchy 点计算链，分段推进并更新一维导数。
  - 邻近语境（简要）：前句“where \(\mathbf{e}_b\) is the \(b\) th unit vector, and \(r_b\) is the \(b\) th component of the gradient vector \(\mathbf{r}\). Similar to Eq. (22), the new v…”；后句“and”。
- **必要推导过程**：
  1. 计算断点 [ t_i ]，在分段二次模型上求段内极小或推进到下一断点。
- **与其他公式的内在联系**：
  - 与离散求解链相连：输出 [ \mathbf{r} ] 与 [ \mathbf{K}/\mathbf{B}_k ] 供迭代器使用。
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与文中编号公式互引：[ (22) ]、[ (23) ]。

### 公式 071

- **出现位置**：`explain.md` 第 523-525 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：[ (24) ]、[ (20) ]。
- **公式完整表达**：
[
f''_j = \mathbf{d}^{(j)}{}^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j)} = f''_{j-1} + 2r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j-1)} + r^2_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{e}_b. \quad (24)
]
- **公式含义**：
  - 属于广义 Cauchy 点计算链，分段推进并更新一维导数。
  - 邻近语境（简要）：前句“and”；后句“Based on the compact representation of the L-BFGS matrix, as shown in Eq. (20),”。
- **必要推导过程**：
  1. 计算断点 [ t_i ]，在分段二次模型上求段内极小或推进到下一断点。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与文中编号公式互引：[ (24) ]、[ (20) ]。

### 公式 072

- **出现位置**：`explain.md` 第 528-530 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：[ (20) ]。
- **公式完整表达**：
[
\mathbf{B} = \mathbf{B}^0 - \mathbf{W}\mathbf{M}\mathbf{W}^{\mathrm{T}}.
]
- **公式含义**：
  - 该式承担局部定义或过渡作用，连接前后推导步骤。
  - 邻近语境（简要）：前句“Based on the compact representation of the L-BFGS matrix, as shown in Eq. (20),”；后句“Let”。
- **必要推导过程**：
  1. 由上一式变量定义代入并做代数整理可得本式。
- **与其他公式的内在联系**：
  - 与文中编号公式互引：[ (20) ]。

### 公式 073

- **出现位置**：`explain.md` 第 532-534 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：[ (23) ]、[ (24) ]。
- **公式完整表达**：
[
\mathbf{w}_b = \mathbf{W}^{\mathrm{T}}\mathbf{e}_b, \quad \mathbf{p}^{(j)} = \mathbf{W}^{\mathrm{T}}\mathbf{d}^{(j)} = \mathbf{W}^{\mathrm{T}}(\mathbf{d}^{(j-1)} + r_b \mathbf{e}_b) = \mathbf{p}^{(j-1)} + r_b \mathbf{w}_b, \quad \mathbf{c}^{(j)} = \mathbf{W}^{\mathrm{T}}\mathbf{z}^{(j)} = \mathbf{W}^{\mathrm{T}}(\mathbf{z}^{(j-1)} + \Delta t^{(j-1)}\mathbf{d}^{(j-1)}) = \mathbf{c}^{(j-1)} + \Delta t^{(j-1)}\mathbf{p}^{(j-1)}.
]
- **公式含义**：
  - 属于广义 Cauchy 点计算链，分段推进并更新一维导数。
  - 邻近语境（简要）：前句“Let”；后句“Plug all the above relationships into Eqs. (23) and (24), the coefficients of the quadratic polynomial \(\hat{p}(\Delta t)\) for the interval \([t^{(j)}, t^{(j…”。
- **必要推导过程**：
  1. 计算断点 [ t_i ]，在分段二次模型上求段内极小或推进到下一断点。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与文中编号公式互引：[ (23) ]、[ (24) ]。

### 公式 074

- **出现位置**：`explain.md` 第 537-539 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：[ (23) ]、[ (24) ]。
- **公式完整表达**：
[
f'_j = f'_{j-1} + \Delta t^{(j-1)} f''_{j-1} + r^2_b + r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{z}^{(j)} - r_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{c}^{(j)}
]
- **公式含义**：
  - 属于广义 Cauchy 点计算链，分段推进并更新一维导数。
  - 邻近语境（简要）：前句“Plug all the above relationships into Eqs. (23) and (24), the coefficients of the quadratic polynomial \(\hat{p}(\Delta t)\) for the interval \([t^{(j)}, t^{(j…”；后句“and”。
- **必要推导过程**：
  1. 计算断点 [ t_i ]，在分段二次模型上求段内极小或推进到下一断点。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与文中编号公式互引：[ (23) ]、[ (24) ]。

### 公式 075

- **出现位置**：`explain.md` 第 541-543 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：未显式编号。
- **公式完整表达**：
[
f''_j = f''_{j-1} + 2r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{d}^{(j-1)} - 2r_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{p}^{(j-1)} + r^2_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{e}_b - r^2_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{w}_b.
]
- **公式含义**：
  - 属于广义 Cauchy 点计算链，分段推进并更新一维导数。
  - 邻近语境（简要）：前句“and”；后句“Similar to the steps in the previous interval, the critical point of the quadratic polynomial \(\hat{p}(\Delta t)\) is calculated as \(\Delta t^* = -f'_j / f''…”。
- **必要推导过程**：
  1. 计算断点 [ t_i ]，在分段二次模型上求段内极小或推进到下一断点。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 076

- **出现位置**：`explain.md` 第 546-548 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathcal{A}(\mathbf{x}^c) = \{i : t_i \leq t^{(j)}, i = 1, \ldots, n\}.
]
- **公式含义**：
  - 用于盒约束可行性维护，防止相场变量违反不可逆与上界约束。
  - 定义活动集，用于后续自由子空间构造。
  - 属于广义 Cauchy 点计算链，分段推进并更新一维导数。
  - 邻近语境（简要）：前句“Similar to the steps in the previous interval, the critical point of the quadratic polynomial \(\hat{p}(\Delta t)\) is calculated as \(\Delta t^* = -f'_j / f''…”；后句“If \(\Delta t^* \geq t^{(j+1)} - t^{(j)}\), the search of the generalized Cauchy point moves onto the next interval (line segment) of the projected gradient. T…”。
- **必要推导过程**：
  1. 按分量截断规则实施投影：低于下界取下界，高于上界取上界。
  2. 根据触边条件更新活动集，并在其补集上继续优化。
  3. 计算断点 [ t_i ]，在分段二次模型上求段内极小或推进到下一断点。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 077

- **出现位置**：`explain.md` 第 556-558 行；章节：[ ### 3.4. Subspace minimization ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathcal{A}(\mathbf{x}^c) = \{i : x^c_i = \mathrm{lb}_i\} \cup \{i : x^c_i = \mathrm{ub}_i\} = \{i : t_i \leq t^{(j)}, i = 1, \ldots, n\}.
]
- **公式含义**：
  - 用于盒约束可行性维护，防止相场变量违反不可逆与上界约束。
  - 定义活动集，用于后续自由子空间构造。
  - 属于广义 Cauchy 点计算链，分段推进并更新一维导数。
  - 邻近语境（简要）：前句“Once the generalized Cauchy point \(\mathbf{x}^c\) is located inside the interval \([\mathbf{x}(t^{(j)}), \mathbf{x}(t^{(j+1)})]\) on the piecewise linear path…”；后句“The components of the solution vector \(\mathbf{x}\) inside the active set are also known, each of which takes either the lower bound or the upper bound value,”。
- **必要推导过程**：
  1. 按分量截断规则实施投影：低于下界取下界，高于上界取上界。
  2. 根据触边条件更新活动集，并在其补集上继续优化。
  3. 计算断点 [ t_i ]，在分段二次模型上求段内极小或推进到下一断点。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 078

- **出现位置**：`explain.md` 第 561-563 行；章节：[ ### 3.4. Subspace minimization ]；文中编号：未显式编号。
- **公式完整表达**：
[
x_i = x^c_i = \begin{cases} \mathrm{ub}_i & \text{if } r_i < 0 \\ \mathrm{lb}_i & \text{if } r_i > 0 \end{cases}, \quad \forall i \in \mathcal{A}(\mathbf{x}^c).
]
- **公式含义**：
  - 用于盒约束可行性维护，防止相场变量违反不可逆与上界约束。
  - 定义活动集，用于后续自由子空间构造。
  - 邻近语境（简要）：前句“The components of the solution vector \(\mathbf{x}\) inside the active set are also known, each of which takes either the lower bound or the upper bound value,”；后句“Assume that the size of the active set \(\mathcal{A}(\mathbf{x}^c)\) is \(q\), that is, there are totally \(q\) components of the generalized Cauchy point \(\m…”。
- **必要推导过程**：
  1. 按分量截断规则实施投影：低于下界取下界，高于上界取上界。
  2. 根据触边条件更新活动集，并在其补集上继续优化。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 079

- **出现位置**：`explain.md` 第 569-571 行；章节：[ ### 3.4. Subspace minimization ]；文中编号：[ (25) ]。
- **公式完整表达**：
[
\mathbf{x}^* = \arg \min m_k(\mathbf{x}) \quad (25)
]
- **公式含义**：
  - 将求解过程明确为带约束极小化问题。
  - 邻近语境（简要）：前句“During the \(k\) th L-BFGS iteration, since at the generalized Cauchy point \(\mathbf{x}^c\), the variables contained in the active set of the box constraints …”；后句“subject to”。
- **必要推导过程**：
  1. 把增量步未知量作为优化变量，配合约束写成极小化问题。
- **与其他公式的内在联系**：
  - 与文中编号公式互引：[ (25) ]。

### 公式 080

- **出现位置**：`explain.md` 第 575-577 行；章节：[ ### 3.4. Subspace minimization ]；文中编号：未显式编号。
- **公式完整表达**：
[
x_i = x^c_i, \quad \forall i \in \mathcal{A}(\mathbf{x}^c)
]
- **公式含义**：
  - 定义活动集，用于后续自由子空间构造。
  - 邻近语境（简要）：前句“subject to”；后句“and”。
- **必要推导过程**：
  1. 根据触边条件更新活动集，并在其补集上继续优化。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 081

- **出现位置**：`explain.md` 第 579-581 行；章节：[ ### 3.4. Subspace minimization ]；文中编号：[ (26) ]、[ (15) ]。
- **公式完整表达**：
[
\mathrm{lb}_i \leq x_i \leq \mathrm{ub}_i, \quad \forall i \notin \mathcal{A}(\mathbf{x}^c). \quad (26)
]
- **公式含义**：
  - 用于盒约束可行性维护，防止相场变量违反不可逆与上界约束。
  - 定义活动集，用于后续自由子空间构造。
  - 邻近语境（简要）：前句“and”；后句“Recall that the objective function \(m_k(\mathbf{x})\) is the quadratic function defined in Eq. (15). There are generally two approaches to solve the above min…”。
- **必要推导过程**：
  1. 按分量截断规则实施投影：低于下界取下界，高于上界取上界。
  2. 根据触边条件更新活动集，并在其补集上继续优化。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与文中编号公式互引：[ (26) ]、[ (15) ]。

### 公式 082

- **出现位置**：`explain.md` 第 589-591 行；章节：[ #### 3.4.1. Direct matrix factorization for the primal approach ]；文中编号：[ (15) ]。
- **公式完整表达**：
[
\mathbf{x} = \mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}},
]
- **公式含义**：
  - 属于子空间最小化/KKT-Scher 补链，生成最终搜索方向。
  - 邻近语境（简要）：前句“At the \(k\) th L-BFGS iteration, let \(\mathbf{x}^c \in \mathbb{R}^n\) represent the generalized Cauchy point, and \(\mathbf{Z}_k \in \mathbb{R}^{n\times (n-q…”；后句“where \(\hat{\mathbf{x}} \in \mathbb{R}^{n-q}\) is the vector of the free variables. The quadratic model in Eq. (15) can be transformed as”。
- **必要推导过程**：
  1. 固定活动变量后，解原始降维系统或对偶 Schur 补系统并回代。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与文中编号公式互引：[ (15) ]。

### 公式 083

- **出现位置**：`explain.md` 第 595-604 行；章节：[ #### 3.4.1. Direct matrix factorization for the primal approach ]；文中编号：[ (15) ]。
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
- **公式含义**：
  - 定义残量或梯度，残量趋零对应离散平衡。
  - 属于子空间最小化/KKT-Scher 补链，生成最终搜索方向。
  - 邻近语境（简要）：前句“where \(\hat{\mathbf{x}} \in \mathbb{R}^{n-q}\) is the vector of the free variables. The quadratic model in Eq. (15) can be transformed as”；后句“Since the L-BFGS matrix \(\mathbf{B}_k\) is positive definite, the above quadratic function \(\hat{m}_k(\hat{\mathbf{x}})\) of the free variables \(\hat{\mathb…”。
- **必要推导过程**：
  1. 固定活动变量后，解原始降维系统或对偶 Schur 补系统并回代。
- **与其他公式的内在联系**：
  - 与离散求解链相连：输出 [ \mathbf{r} ] 与 [ \mathbf{K}/\mathbf{B}_k ] 供迭代器使用。
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与文中编号公式互引：[ (15) ]。

### 公式 084

- **出现位置**：`explain.md` 第 608-610 行；章节：[ #### 3.4.1. Direct matrix factorization for the primal approach ]；文中编号：未显式编号。
- **公式完整表达**：
[
\hat{m}'_k(\hat{\mathbf{x}}) = 0 \implies \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k \hat{\mathbf{x}} = -\mathbf{Z}^{\mathrm{T}}_k [\mathbf{r}_k + \mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k)].
]
- **公式含义**：
  - 定义残量或梯度，残量趋零对应离散平衡。
  - 属于子空间最小化/KKT-Scher 补链，生成最终搜索方向。
  - 邻近语境（简要）：前句“Since the L-BFGS matrix \(\mathbf{B}_k\) is positive definite, the above quadratic function \(\hat{m}_k(\hat{\mathbf{x}})\) of the free variables \(\hat{\mathb…”；后句“Let \(\hat{\mathbf{B}}_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k\) and \(\hat{\mathbf{r}}_k = \mathbf{Z}^{\mathrm{T}}_k [\mathbf{r}_k + \mathbf{B…”。
- **必要推导过程**：
  1. 固定活动变量后，解原始降维系统或对偶 Schur 补系统并回代。
- **与其他公式的内在联系**：
  - 与离散求解链相连：输出 [ \mathbf{r} ] 与 [ \mathbf{K}/\mathbf{B}_k ] 供迭代器使用。
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 085

- **出现位置**：`explain.md` 第 614-616 行；章节：[ #### 3.4.1. Direct matrix factorization for the primal approach ]；文中编号：[ (27) ]。
- **公式完整表达**：
[
\hat{\mathbf{B}}_k \hat{\mathbf{x}} = -\hat{\mathbf{r}}_k. \quad (27)
]
- **公式含义**：
  - 定义残量或梯度，残量趋零对应离散平衡。
  - 属于子空间最小化/KKT-Scher 补链，生成最终搜索方向。
  - 邻近语境（简要）：前句“Let \(\hat{\mathbf{B}}_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k\) and \(\hat{\mathbf{r}}_k = \mathbf{Z}^{\mathrm{T}}_k [\mathbf{r}_k + \mathbf{B…”；后句“Recall that \(\mathbf{B}_k\) is not directly available, since it is not stored in the memory component wise. Rather, \(\mathbf{B}_k\) is only available in the …”。
- **必要推导过程**：
  1. 固定活动变量后，解原始降维系统或对偶 Schur 补系统并回代。
- **与其他公式的内在联系**：
  - 与离散求解链相连：输出 [ \mathbf{r} ] 与 [ \mathbf{K}/\mathbf{B}_k ] 供迭代器使用。
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与文中编号公式互引：[ (27) ]。

### 公式 086

- **出现位置**：`explain.md` 第 619-621 行；章节：[ #### 3.4.1. Direct matrix factorization for the primal approach ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{B}_k = \mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k.
]
- **公式含义**：
  - 属于 L-BFGS 紧凑更新链，利用历史向量近似曲率。
  - 邻近语境（简要）：前句“Recall that \(\mathbf{B}_k\) is not directly available, since it is not stored in the memory component wise. Rather, \(\mathbf{B}_k\) is only available in the …”；后句“Therefore,”。
- **必要推导过程**：
  1. 由 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对构造紧凑低秩更新。
- **与其他公式的内在联系**：
  - 与低秩线代链相连：L-BFGS 紧凑更新与 SMW 逆更新互补。

### 公式 087

- **出现位置**：`explain.md` 第 623-625 行；章节：[ #### 3.4.1. Direct matrix factorization for the primal approach ]；文中编号：[ (28) ]、[ (27) ]。
- **公式完整表达**：
[
\hat{\mathbf{B}}_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k = \mathbf{Z}^{\mathrm{T}}_k (\mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k) \mathbf{Z}_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{Z}_k - \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k. \quad (28)
]
- **公式含义**：
  - 属于 L-BFGS 紧凑更新链，利用历史向量近似曲率。
  - 属于子空间最小化/KKT-Scher 补链，生成最终搜索方向。
  - 邻近语境（简要）：前句“Therefore,”；后句“Similar to the L-BFGS matrix \(\mathbf{B}_k\), the reduced matrix \(\hat{\mathbf{B}}_k\) is not stored in the memory component wise. As a result, the direct so…”。
- **必要推导过程**：
  1. 由 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对构造紧凑低秩更新。
  2. 固定活动变量后，解原始降维系统或对偶 Schur 补系统并回代。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与低秩线代链相连：L-BFGS 紧凑更新与 SMW 逆更新互补。
  - 与文中编号公式互引：[ (28) ]、[ (27) ]。

### 公式 088

- **出现位置**：`explain.md` 第 631-633 行；章节：[ #### 3.4.1. Direct matrix factorization for the primal approach ]；文中编号：[ (29) ]。
- **公式完整表达**：
[
\hat{\mathbf{x}} = -\hat{\mathbf{B}}^{-1}_k \hat{\mathbf{r}}_k,
]
- **公式含义**：
  - 定义残量或梯度，残量趋零对应离散平衡。
  - 邻近语境（简要）：前句“Using the Sherman–Morrison–Woodbury formula [45], as stated in the Appendix, the vector of the free variables can be solved as”；后句“\hat{\mathbf{B}}^{-1}_k = (\hat{\mathbf{B}}^0_k)^{-1} + (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \left[\mathbf{I} - \mathbf{M}_k \mat…”。
- **必要推导过程**：
  1. 由上一式变量定义代入并做代数整理可得本式。
- **与其他公式的内在联系**：
  - 与离散求解链相连：输出 [ \mathbf{r} ] 与 [ \mathbf{K}/\mathbf{B}_k ] 供迭代器使用。
  - 与文中编号公式互引：[ (29) ]。

### 公式 089

- **出现位置**：`explain.md` 第 634-636 行；章节：[ #### 3.4.1. Direct matrix factorization for the primal approach ]；文中编号：[ (29) ]。
- **公式完整表达**：
[
\hat{\mathbf{B}}^{-1}_k = (\hat{\mathbf{B}}^0_k)^{-1} + (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \left[\mathbf{I} - \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \right]^{-1} \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1}, \quad (29)
]
- **公式含义**：
  - 属于 L-BFGS 紧凑更新链，利用历史向量近似曲率。
  - 属于子空间最小化/KKT-Scher 补链，生成最终搜索方向。
  - 邻近语境（简要）：前句“\hat{\mathbf{x}} = -\hat{\mathbf{B}}^{-1}_k \hat{\mathbf{r}}_k,”；后句“where”。
- **必要推导过程**：
  1. 由 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对构造紧凑低秩更新。
  2. 固定活动变量后，解原始降维系统或对偶 Schur 补系统并回代。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与低秩线代链相连：L-BFGS 紧凑更新与 SMW 逆更新互补。
  - 与文中编号公式互引：[ (29) ]。

### 公式 090

- **出现位置**：`explain.md` 第 639-641 行；章节：[ #### 3.4.1. Direct matrix factorization for the primal approach ]；文中编号：[ (29) ]。
- **公式完整表达**：
[
\mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k = \mathbf{Z}^{\mathrm{T}}_k [\mathbf{Y}_k \quad \mathbf{B}^0_k \mathbf{S}_k] = [\mathbf{Z}^{\mathrm{T}}_k \mathbf{y}_{k-m} \cdots \mathbf{Z}^{\mathrm{T}}_k \mathbf{y}_{k-1} \quad \mathbf{Z}^{\mathrm{T}}_k(\mathbf{B}^0_k \mathbf{s}_{k-m}) \cdots \mathbf{Z}^{\mathrm{T}}_k(\mathbf{B}^0_k \mathbf{s}_{k-1})].
]
- **公式含义**：
  - 属于 L-BFGS 紧凑更新链，利用历史向量近似曲率。
  - 邻近语境（简要）：前句“where”；后句“In order to calculate the inverse of \(\hat{\mathbf{B}}_k\), the inverse of \(\hat{\mathbf{B}}^0_k\) needs to be repeatedly applied to the column vectors in \(…”。
- **必要推导过程**：
  1. 由 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对构造紧凑低秩更新。
- **与其他公式的内在联系**：
  - 与低秩线代链相连：L-BFGS 紧凑更新与 SMW 逆更新互补。
  - 与文中编号公式互引：[ (29) ]。

### 公式 091

- **出现位置**：`explain.md` 第 644-646 行；章节：[ #### 3.4.1. Direct matrix factorization for the primal approach ]；文中编号：[ (29) ]。
- **公式完整表达**：
[
\mathbf{I} - \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \in \mathbb{R}^{2m\times 2m}
]
- **公式含义**：
  - 属于 L-BFGS 紧凑更新链，利用历史向量近似曲率。
  - 属于子空间最小化/KKT-Scher 补链，生成最终搜索方向。
  - 邻近语境（简要）：前句“In order to calculate the inverse of \(\hat{\mathbf{B}}_k\), the inverse of \(\hat{\mathbf{B}}^0_k\) needs to be repeatedly applied to the column vectors in \(…”；后句“is a \(2m\times 2m\) matrix, therefore, its inverse can be directly calculated with a negligible computational cost.”。
- **必要推导过程**：
  1. 由 [ \mathbf{s}_k,\mathbf{y}_k ] 历史对构造紧凑低秩更新。
  2. 固定活动变量后，解原始降维系统或对偶 Schur 补系统并回代。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与低秩线代链相连：L-BFGS 紧凑更新与 SMW 逆更新互补。
  - 与文中编号公式互引：[ (29) ]。

### 公式 092

- **出现位置**：`explain.md` 第 655-657 行；章节：[ #### 3.4.2. Conjugate gradient method for the primal approach ]；文中编号：[ (18) ]。
- **公式完整表达**：
[
\mathbf{P}^{-1} = \mathrm{ILU}(\hat{\mathbf{B}}^0_k) = \mathrm{ILU}(\mathbf{Z}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{Z}_k).
]
- **公式含义**：
  - 属于子空间最小化/KKT-Scher 补链，生成最终搜索方向。
  - 邻近语境（简要）：前句“It is well known that an effective preconditioner could significantly improve the performance of an iterative linear solver such as the conjugate gradient meth…”；后句“**Comment 4.** The reduced matrix \(\hat{\mathbf{B}}^0_k\) is positive definite since the curvature condition in Eq. (18) is satisfied due to the line search b…”。
- **必要推导过程**：
  1. 固定活动变量后，解原始降维系统或对偶 Schur 补系统并回代。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与文中编号公式互引：[ (18) ]。

### 公式 093

- **出现位置**：`explain.md` 第 665-667 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：[ (25) ]、[ (26) ]、[ (15) ]。
- **公式完整表达**：
[
\mathbf{x} = \mathbf{x}_k + \Delta \mathbf{x}_k.
]
- **公式含义**：
  - 该式承担局部定义或过渡作用，连接前后推导步骤。
  - 邻近语境（简要）：前句“Instead of formulating the minimization problem in Eq. (25) subject to the constraints in Eq. (26) as a subspace minimization, an alternative approach is to fo…”；后句“Plug the above equation into the quadratic model shown in Eq. (15),”。
- **必要推导过程**：
  1. 由上一式变量定义代入并做代数整理可得本式。
- **与其他公式的内在联系**：
  - 与文中编号公式互引：[ (25) ]、[ (26) ]、[ (15) ]。

### 公式 094

- **出现位置**：`explain.md` 第 671-673 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：[ (15) ]。
- **公式完整表达**：
[
m_k(\mathbf{x}) = \Pi_k + \mathbf{r}_k^{\mathrm{T}} \Delta \mathbf{x}_k + \frac{1}{2} \Delta \mathbf{x}_k^{\mathrm{T}} \mathbf{B}_k \Delta \mathbf{x}_k.
]
- **公式含义**：
  - 定义残量或梯度，残量趋零对应离散平衡。
  - 邻近语境（简要）：前句“Plug the above equation into the quadratic model shown in Eq. (15),”；后句“Also, at the generalized Cauchy point \(\mathbf{x}^c\), for the variables in the active set of the box constraints,”。
- **必要推导过程**：
  1. 由上一式变量定义代入并做代数整理可得本式。
- **与其他公式的内在联系**：
  - 与离散求解链相连：输出 [ \mathbf{r} ] 与 [ \mathbf{K}/\mathbf{B}_k ] 供迭代器使用。
  - 与文中编号公式互引：[ (15) ]。

### 公式 095

- **出现位置**：`explain.md` 第 676-678 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：未显式编号。
- **公式完整表达**：
[
x_i = x^c_i, \quad \forall i \in \mathcal{A}_k(\mathbf{x}^c),
]
- **公式含义**：
  - 定义活动集，用于后续自由子空间构造。
  - 邻近语境（简要）：前句“Also, at the generalized Cauchy point \(\mathbf{x}^c\), for the variables in the active set of the box constraints,”；后句“where \(x^c_i\) takes the upper bound value \(\mathrm{ub}_i\) or the lower bound value \(\mathrm{lb}_i\) of the box constraints. Recall that the matrix \(\math…”。
- **必要推导过程**：
  1. 根据触边条件更新活动集，并在其补集上继续优化。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 096

- **出现位置**：`explain.md` 第 680-682 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：[ (25) ]、[ (26) ]。
- **公式完整表达**：
[
\mathbf{Q}^{\mathrm{T}}_k \Delta \mathbf{x}_k = \mathbf{Q}^{\mathrm{T}}_k (\mathbf{x} - \mathbf{x}_k) = \mathbf{Q}^{\mathrm{T}}_k (\mathbf{x}^c - \mathbf{x}_k).
]
- **公式含义**：
  - 该式承担局部定义或过渡作用，连接前后推导步骤。
  - 邻近语境（简要）：前句“where \(x^c_i\) takes the upper bound value \(\mathrm{ub}_i\) or the lower bound value \(\mathrm{lb}_i\) of the box constraints. Recall that the matrix \(\math…”；后句“The minimization problem in Eq. (25) and the corresponding constraints in Eq. (26) can be rewritten as”。
- **必要推导过程**：
  1. 由上一式变量定义代入并做代数整理可得本式。
- **与其他公式的内在联系**：
  - 与文中编号公式互引：[ (25) ]、[ (26) ]。

### 公式 097

- **出现位置**：`explain.md` 第 685-687 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：[ (25) ]、[ (26) ]。
- **公式完整表达**：
[
\min \left\{ \Pi_k + \mathbf{r}_k^{\mathrm{T}} \Delta \mathbf{x}_k + \frac{1}{2} \Delta \mathbf{x}_k^{\mathrm{T}} \mathbf{B}_k \Delta \mathbf{x}_k \right\}
]
- **公式含义**：
  - 定义残量或梯度，残量趋零对应离散平衡。
  - 邻近语境（简要）：前句“The minimization problem in Eq. (25) and the corresponding constraints in Eq. (26) can be rewritten as”；后句“subject to”。
- **必要推导过程**：
  1. 由上一式变量定义代入并做代数整理可得本式。
- **与其他公式的内在联系**：
  - 与离散求解链相连：输出 [ \mathbf{r} ] 与 [ \mathbf{K}/\mathbf{B}_k ] 供迭代器使用。
  - 与文中编号公式互引：[ (25) ]、[ (26) ]。

### 公式 098

- **出现位置**：`explain.md` 第 689-691 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{Q}_k^{\mathrm{T}} \Delta \mathbf{x}_k = \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k)
]
- **公式含义**：
  - 属于子空间最小化/KKT-Scher 补链，生成最终搜索方向。
  - 邻近语境（简要）：前句“subject to”；后句“and the box constraints”。
- **必要推导过程**：
  1. 固定活动变量后，解原始降维系统或对偶 Schur 补系统并回代。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 099

- **出现位置**：`explain.md` 第 693-695 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{lb} - \mathbf{x}_k \leq \Delta \mathbf{x}_k \leq \mathbf{ub} - \mathbf{x}_k.
]
- **公式含义**：
  - 用于盒约束可行性维护，防止相场变量违反不可逆与上界约束。
  - 邻近语境（简要）：前句“and the box constraints”；后句“For now, ignore the box constraints and let \(\lambda\) represent the Lagrange multiplier for the equality constraints. The optimality condition of the above c…”。
- **必要推导过程**：
  1. 按分量截断规则实施投影：低于下界取下界，高于上界取上界。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 100

- **出现位置**：`explain.md` 第 699-701 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：[ (30) ]。
- **公式完整表达**：
[
\begin{bmatrix} \mathbf{B}_k & \mathbf{Q}_k \\ \mathbf{Q}_k^{\mathrm{T}} & \mathbf{0} \end{bmatrix} \begin{bmatrix} \Delta \mathbf{x}_k \\ \lambda \end{bmatrix} = \begin{bmatrix} -\mathbf{r}_k \\ \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k) \end{bmatrix}. \quad (30)
]
- **公式含义**：
  - 定义残量或梯度，残量趋零对应离散平衡。
  - 属于子空间最小化/KKT-Scher 补链，生成最终搜索方向。
  - 邻近语境（简要）：前句“For now, ignore the box constraints and let \(\lambda\) represent the Lagrange multiplier for the equality constraints. The optimality condition of the above c…”；后句“A popular method to solve the above linear system is to use the Schur complement approach by firstly solving \(\lambda\) from”。
- **必要推导过程**：
  1. 固定活动变量后，解原始降维系统或对偶 Schur 补系统并回代。
- **与其他公式的内在联系**：
  - 与离散求解链相连：输出 [ \mathbf{r} ] 与 [ \mathbf{K}/\mathbf{B}_k ] 供迭代器使用。
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与文中编号公式互引：[ (30) ]。

### 公式 101

- **出现位置**：`explain.md` 第 705-707 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：[ (31) ]。
- **公式完整表达**：
[
(\mathbf{Q}_k^{\mathrm{T}} \mathbf{B}_k^{-1} \mathbf{Q}_k) \lambda = -\mathbf{Q}_k^{\mathrm{T}} \mathbf{B}_k^{-1} \mathbf{r}_k - \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k). \quad (31)
]
- **公式含义**：
  - 定义残量或梯度，残量趋零对应离散平衡。
  - 属于子空间最小化/KKT-Scher 补链，生成最终搜索方向。
  - 邻近语境（简要）：前句“A popular method to solve the above linear system is to use the Schur complement approach by firstly solving \(\lambda\) from”；后句“Then, solve \(\Delta \mathbf{x}_k\) from”。
- **必要推导过程**：
  1. 固定活动变量后，解原始降维系统或对偶 Schur 补系统并回代。
- **与其他公式的内在联系**：
  - 与离散求解链相连：输出 [ \mathbf{r} ] 与 [ \mathbf{K}/\mathbf{B}_k ] 供迭代器使用。
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与文中编号公式互引：[ (31) ]。

### 公式 102

- **出现位置**：`explain.md` 第 711-713 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：[ (32) ]、[ (31) ]。
- **公式完整表达**：
[
\mathbf{B}_k \Delta \mathbf{x}_k = -(\mathbf{r}_k + \mathbf{Q}_k \lambda). \quad (32)
]
- **公式含义**：
  - 定义残量或梯度，残量趋零对应离散平衡。
  - 属于子空间最小化/KKT-Scher 补链，生成最终搜索方向。
  - 邻近语境（简要）：前句“Then, solve \(\Delta \mathbf{x}_k\) from”；后句“Notice that \(\mathbf{B}_k^{-1}\) is contained in Eq. (31). Recall that \(\mathbf{B}_k\) only has the compact representation form. Therefore, similar to the pr…”。
- **必要推导过程**：
  1. 固定活动变量后，解原始降维系统或对偶 Schur 补系统并回代。
- **与其他公式的内在联系**：
  - 与离散求解链相连：输出 [ \mathbf{r} ] 与 [ \mathbf{K}/\mathbf{B}_k ] 供迭代器使用。
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与文中编号公式互引：[ (32) ]、[ (31) ]。

### 公式 103

- **出现位置**：`explain.md` 第 719-721 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：[ (27) ]。
- **公式完整表达**：
[
\mathbf{p}_k = \mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}} - \mathbf{x}_k.
]
- **公式含义**：
  - 定义线搜索更新与步长准则。
  - 属于子空间最小化/KKT-Scher 补链，生成最终搜索方向。
  - 邻近语境（简要）：前句“Once the free variables \(\hat{\mathbf{x}}\) at the generalized Cauchy point \(\mathbf{x}^c\) are solved from the reduced linear system shown in Eq. (27), the …”；后句“The updated solution is obtained as”。
- **必要推导过程**：
  1. 给定方向后使用 Armijo/Wolfe 条件搜索 [ \alpha_k ]。
  2. 固定活动变量后，解原始降维系统或对偶 Schur 补系统并回代。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与文中编号公式互引：[ (27) ]。

### 公式 104

- **出现位置**：`explain.md` 第 725-727 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{x}_{k + 1} = \mathbf{x}_k + \alpha_k \mathbf{p}_k,
]
- **公式含义**：
  - 定义线搜索更新与步长准则。
  - 邻近语境（简要）：前句“The updated solution is obtained as”；后句“where the step length parameter \(\alpha_k\) is determined by simultaneously satisfying the following two conditions:”。
- **必要推导过程**：
  1. 给定方向后使用 Armijo/Wolfe 条件搜索 [ \alpha_k ]。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。

### 公式 105

- **出现位置**：`explain.md` 第 732-734 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：[ (18) ]。
- **公式完整表达**：
[
\mathbf{lb} \leq \mathbf{x}_k + \alpha_k \mathbf{p}_k \leq \mathbf{ub}.
]
- **公式含义**：
  - 用于盒约束可行性维护，防止相场变量违反不可逆与上界约束。
  - 定义线搜索更新与步长准则。
  - 邻近语境（简要）：前句“1. The updated solution \(\mathbf{x}_{k + 1}\) should stay feasible, that is,”；后句“2. The positive step length \(\alpha_k\) should satisfy the strong Wolfe conditions such that the curvature condition in Eq. (18) is satisfied to ensure the po…”。
- **必要推导过程**：
  1. 按分量截断规则实施投影：低于下界取下界，高于上界取上界。
  2. 给定方向后使用 Armijo/Wolfe 条件搜索 [ \alpha_k ]。
- **与其他公式的内在联系**：
  - 与约束优化链相连：可行投影 → 活动集识别 → 子空间/线搜索更新。
  - 与文中编号公式互引：[ (18) ]。

### 公式 106

- **出现位置**：`explain.md` 第 780-782 行；章节：[ ### 4.1. Cyclic tension-compression test ]；文中编号：[ (33) ]。
- **公式完整表达**：
[
E_{\mathrm{crack}} = g_{\mathrm{c}}\Gamma_{l}(d) = g_{\mathrm{c}}\int_{\Omega}\gamma (d,\nabla d)\mathrm{d}\Omega = g_{\mathrm{c}}\int_{\Omega}\frac{1}{2l}\left(d^{2} + l^{2}\nabla d\cdot \nabla d\right)\mathrm{d}\Omega \quad (33)
]
- **公式含义**：
  - 给出裂纹正则项 [ \Gamma_l ] 或裂纹密度 [ \gamma ]，控制裂纹带宽与能量耗散。
  - 定义算例中的加载或能量观测量，连接数值结果与物理解释。
  - 邻近语境（简要）：前句“Fig. 4 shows the phase-field distribution in the unit square under the vertical cyclic load at various pseudo time steps. When the vertical displacement increa…”；后句“with and without enforcing the phase-field irreversibility condition. When the irreversibility condition is properly enforced by the proposed L-BFGS-B scheme, …”。
- **必要推导过程**：
  1. 采用相场正则化 [ \gamma=\frac{1}{2l}(d^2+l^2|\nabla d|^2) ]，对域积分得到 [ \Gamma_l ]。
  2. 将控制参数代入边界位移/能量积分表达，形成可比工程指标。
- **与其他公式的内在联系**：
  - 与本构—能量链相连：上接总势能定义，下接残量与切线离散。
  - 与算例验证链相连：将算法输出映射为可观测物理量。
  - 与文中编号公式互引：[ (33) ]。

### 公式 107

- **出现位置**：`explain.md` 第 878-880 行；章节：[ ### 4.4. Three-dimensional torsion test ]；文中编号：未显式编号。
- **公式完整表达**：
[
u_{y} = z \tan t, \quad u_{z} = -y \tan t,
]
- **公式含义**：
  - 定义算例中的加载或能量观测量，连接数值结果与物理解释。
  - 邻近语境（简要）：前句“In the last example, a 3D beam shown in Fig. 17 undergoes a torsional load. The beam has a square cross section of \(50 \times 50 \mathrm{mm}\) and a length of…”；后句“where \(t\) is the pseudo time representing the load step. The initial mesh is pre-refined near the bottom region of the notch where the crack is expected to i…”。
- **必要推导过程**：
  1. 将控制参数代入边界位移/能量积分表达，形成可比工程指标。
- **与其他公式的内在联系**：
  - 与算例验证链相连：将算法输出映射为可观测物理量。

### 公式 108

- **出现位置**：`explain.md` 第 937-939 行；章节：[ ### 5.2. Comparison of convergence behaviors ]；文中编号：未显式编号。
- **公式完整表达**：
[
\| \pmb {r}_u\| _2 < \mathrm{tol}, \quad \| \pmb {r}_d\| _2 < \mathrm{tol}, \quad \| \Delta \pmb {u}\| _2 < \mathrm{tol}, \quad \| \Delta d\| _2 < \mathrm{tol},
]
- **公式含义**：
  - 该式承担局部定义或过渡作用，连接前后推导步骤。
  - 邻近语境（简要）：前句“For the staggered scheme using the alternate minimization, inside one iteration, the nonlinear displacement sub-problem is firstly solved via the Newton-Raphso…”；后句“where \(\pmb{r}_{u}\) and \(\pmb{r}_{d}\) represent the residuals of the displacement and the phase-field sub-problems, and \(\Delta \pmb {u}\) and \(\Delta d\…”。
- **必要推导过程**：
  1. 由上一式变量定义代入并做代数整理可得本式。
- **与其他公式的内在联系**：
  - 在局部推导中承前启后，为下一式提供变量或算子。

### 公式 109

- **出现位置**：`explain.md` 第 973-975 行；章节：[ ## Appendix ]；文中编号：未显式编号。
- **公式完整表达**：
[
\hat{\mathbf{A}} = \mathbf{A} + \mathbf{U}\mathbf{V}^{\mathrm{T}},
]
- **公式含义**：
  - 给出低秩逆更新公式，是高效线性代数实现的理论基础。
  - 邻近语境（简要）：前句“Let \(\mathbf{U}\in \mathbb{R}^{n\times p}\) and \(\mathbf{V}\in \mathbb{R}^{n\times p}\) represent two matrices for some \(p\) between 1 and \(n\). If”；后句“according to the Sherman–Morrison–Woodbury formula [45], \(\hat{\mathbf{A}}\) is non-singular if and only if \((\mathbf{I} + \mathbf{V}^{\mathrm{T}}\mathbf{A}^…”。
- **必要推导过程**：
  1. 将矩阵写为“基矩阵+低秩修正”，应用 SMW 公式化简求逆。
- **与其他公式的内在联系**：
  - 与低秩线代链相连：L-BFGS 紧凑更新与 SMW 逆更新互补。

### 公式 110

- **出现位置**：`explain.md` 第 979-981 行；章节：[ ## Appendix ]；文中编号：未显式编号。
- **公式完整表达**：
[
\hat{\mathbf{A}}^{-1} = \mathbf{A}^{-1} - \mathbf{A}^{-1}\mathbf{U}\left(\mathbf{I} + \mathbf{V}^{\mathrm{T}}\mathbf{A}^{-1}\mathbf{U}\right)^{-1}\mathbf{V}^{\mathrm{T}}\mathbf{A}^{-1}.
]
- **公式含义**：
  - 给出低秩逆更新公式，是高效线性代数实现的理论基础。
  - 邻近语境（简要）：前句“according to the Sherman–Morrison–Woodbury formula [45], \(\hat{\mathbf{A}}\) is non-singular if and only if \((\mathbf{I} + \mathbf{V}^{\mathrm{T}}\mathbf{A}^…”；后句“## Data availability”。
- **必要推导过程**：
  1. 将矩阵写为“基矩阵+低秩修正”，应用 SMW 公式化简求逆。
- **与其他公式的内在联系**：
  - 与低秩线代链相连：L-BFGS 紧凑更新与 SMW 逆更新互补。