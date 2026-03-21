# explain.md 第2章公式逐一详解（完整推导与物理含义）

> 说明：本文档严格对应 `explain.md` 的 **2. Phase-field formulation and finite element discretization**（即 2.1 与 2.2）中出现的全部公式。  
> 组织方式遵循论文分段；每条公式都先完整展示，再给出逐步解释、必要推导、符号含义与物理解释。  
> 为满足格式要求，公式统一采用如下渲染形式：
>
> \[
> \text{公式内容}
> \]

---

## 2.1 Phase-field formulation

### 公式(5)：张拉-受压不对称下的能量分解

\[
\psi(\pmb{\epsilon}, d) = [g(d) + k] \, \psi^{+}(\pmb{\epsilon}) + \psi^{-}(\pmb{\epsilon}), \quad (5)
\]

**逐项含义**

- \[ \psi(\pmb{\epsilon}, d) \]：退化后材料在当前应变 \[\pmb{\epsilon}\] 与相场 \[d\] 下的总应变能密度。
- \[ \psi^{+}(\pmb{\epsilon}) \]：正应变能（主导裂纹扩展的张拉部分）。
- \[ \psi^{-}(\pmb{\epsilon}) \]：负应变能（受压部分，通常不直接驱动开裂）。
- \[ g(d) \]：相场退化函数。
- \[ k \]：小的非负残余刚度参数（避免完全退化导致数值奇异；文中算例取 0）。

**为何这样分解**

1. 脆性断裂中，裂纹在张拉状态下扩展明显，而纯压缩状态不应产生同等损伤驱动力。
2. 因此只对 \[\psi^{+}\] 乘退化系数 \([g(d)+k]\)，\[\psi^{-}\] 保持不退化或弱退化，从而体现张压不对称。

**物理解读**

- 当 \[d\to 1\]（完全破坏）时，张拉承载能力趋于显著降低；
- 受压部分 \[\psi^{-}\] 仍保留，使裂纹面在压缩闭合时具有合理响应。

---

### 公式(6)：退化函数定义

\[
g(d) = (1 - d)^{2}. \quad (6)
\]

**逐项含义**

- \[d\in[0,1]\]：相场变量，\[d=0\] 未损伤，\[d=1\] 完全损伤。
- \[g(0)=1\]：完好材料不退化。
- \[g(1)=0\]：完全损伤后张拉刚度退化到零（若 \[k=0\]）。

**必要推导（导数）**

\[
g'(d) = -2(1-d), \qquad g''(d)=2.
\]

这些导数在后续残量与 Hessian（公式(9)、(10)）中直接出现。

**物理解读**

- 二次函数在 \[d\] 方向平滑、单调递减，兼顾数值稳定与物理合理性。

---

### 公式(2.1-a)：Macaulay 括号与 Heaviside 函数

\[
\langle x\rangle_{+} = \frac{1}{2}(x + |x|), \quad \langle x\rangle_{-} = \frac{1}{2}(x - |x|), \quad H(x) = \begin{cases} 1 & x \geq 0, \\ 0 & x < 0. \end{cases}
\]

**逐项含义**

- \[\langle x\rangle_{+}=\max(x,0)\]，\[\langle x\rangle_{-}=\min(x,0)\]。
- \[H(x)\] 为阶跃函数，后续用于切线模量中区分拉压区间。

**必要推导**

由绝对值分段定义可得：

- 若 \[x\ge 0\]，\[|x|=x\]，所以 \[\langle x\rangle_+=x,\;\langle x\rangle_-=0\]；
- 若 \[x<0\]，\[|x|=-x\]，所以 \[\langle x\rangle_+=0,\;\langle x\rangle_-=x\]。

并且（在可导区间）

\[
\frac{\mathrm d}{\mathrm dx}\langle x\rangle_{+}=H(x), \qquad \frac{\mathrm d}{\mathrm dx}\langle x\rangle_{-}=H(-x).
\]

这正是后续出现 \[H(\mathrm{tr}\,\pmb\epsilon)\] 与 \[H(-\mathrm{tr}\,\pmb\epsilon)\] 的来源。

---

### 公式(2.1-b)：应变谱分解

\[
\pmb{\epsilon} = \sum_{\alpha} \epsilon_{\alpha} \mathbf{M}_{\alpha}, \quad \mathbf{M}_{\alpha} = \pmb{n}_{\alpha} \otimes \pmb{n}_{\alpha},
\]

**逐项含义**

- \[\epsilon_\alpha\]：主应变（特征值）。
- \[\pmb n_\alpha\]：对应特征向量（主方向）。
- \[\mathbf M_\alpha\]：谱投影张量，将某方向分量投影出来。

**必要推导思路**

对对称小应变张量 \[\pmb\epsilon\]，可做正交特征分解：

\[
\pmb\epsilon = \mathbf Q\,\mathrm{diag}(\epsilon_1,\epsilon_2,\epsilon_3)\,\mathbf Q^\mathrm T
\]

等价写成谱和式即上式。该表示使“仅退化正主应变方向”可被直接表达。

---

### 公式(2.1-c)：正负应变张量

\[
\pmb{\epsilon}^{+} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{+} \mathbf{M}_{\alpha}, \quad \pmb{\epsilon}^{-} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{-} \mathbf{M}_{\alpha}.
\]

**逐项含义**

- \[\pmb\epsilon^+\] 聚合全部“拉伸主应变”贡献；
- \[\pmb\epsilon^-\] 聚合全部“压缩主应变”贡献。

**直接性质**

\[
\pmb\epsilon = \pmb\epsilon^+ + \pmb\epsilon^-.
\]

这是后续能量加性分裂（张拉/压缩）在张量层面的基础。

---

### 公式(2.1-d)：正负应变能密度

\[
\psi^{+}(\pmb{\epsilon}) = \frac{1}{2}\lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{+}^{2} + \mu \, \pmb{\epsilon}^{+} : \pmb{\epsilon}^{+}, \quad \psi^{-}(\pmb{\epsilon}) = \frac{1}{2} \lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{-}^{2} + \mu \, \pmb{\epsilon}^{-} : \pmb{\epsilon}^{-},
\]

**逐项含义**

- \[\lambda,\mu\]：Lamé 常数。
- \[\mathrm{tr}\pmb\epsilon\]：体应变（体积变化）指标。
- \[\pmb\epsilon^{\pm}:\pmb\epsilon^{\pm}\]：对应偏/总应变的二次型能量项。

**必要推导（与线弹性能量关系）**

经典线弹性能量

\[
\psi_0(\pmb\epsilon)=\frac{\lambda}{2}(\mathrm{tr}\pmb\epsilon)^2+\mu\,\pmb\epsilon:\pmb\epsilon
\]

将其按张拉/压缩分量替换为

\[
\mathrm{tr}\pmb\epsilon \rightarrow \langle\mathrm{tr}\pmb\epsilon\rangle_+\;\text{与}\;\langle\mathrm{tr}\pmb\epsilon\rangle_-
\]

及

\[
\pmb\epsilon \rightarrow \pmb\epsilon^+\;\text{与}\;\pmb\epsilon^-
\]

即可得到上式，保证 \[\psi^+\ge 0\] 且主要代表开裂驱动力。

---

### 公式(2.1-e)：总应力表达

\[
\pmb{\sigma} = \frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial \pmb{\epsilon}} = [g(d) + k]\frac{\partial\psi^{+}(\pmb{\epsilon})}{\partial \pmb{\epsilon}} + \frac{\partial\psi^{-}(\pmb{\epsilon})}{\partial \pmb{\epsilon}} = [g(d) + k]\pmb{\sigma}^{+} + \pmb{\sigma}^{-},
\]

**必要推导**

由公式(5)

\[
\psi = [g(d)+k]\psi^+(\pmb\epsilon)+\psi^-(\pmb\epsilon)
\]

对 \[\pmb\epsilon\] 求偏导，且 \[d\] 固定：

\[
\frac{\partial\psi}{\partial\pmb\epsilon} = [g(d)+k]\frac{\partial\psi^+}{\partial\pmb\epsilon} + \frac{\partial\psi^-}{\partial\pmb\epsilon}.
\]

定义

\[
\pmb\sigma^+ := \frac{\partial\psi^+}{\partial\pmb\epsilon},\qquad
\pmb\sigma^- := \frac{\partial\psi^-}{\partial\pmb\epsilon},
\]

即得上式。

**物理解读**

- 张拉应力部分被 \[g(d)+k\] 缩放；
- 压缩应力部分不缩放（或弱缩放），体现裂纹闭合效应。

---

### 公式(2.1-f)：正负应力显式表达

\[
\pmb{\sigma}^{+} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{+}\mathbf{I} + 2\mu \pmb{\epsilon}^{+}, \quad \pmb{\sigma}^{-} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{-}\mathbf{I} + 2\mu \pmb{\epsilon}^{-}.
\]

**必要推导**

由公式(2.1-d)：

\[
\psi^+ = \frac{\lambda}{2}\langle\mathrm{tr}\pmb\epsilon\rangle_+^2 + \mu\,\pmb\epsilon^+:\pmb\epsilon^+.
\]

对 \[\pmb\epsilon\] 求导：

1. 第一项导数
\[
\frac{\partial}{\partial\pmb\epsilon}\Big(\frac{\lambda}{2}\langle\mathrm{tr}\pmb\epsilon\rangle_+^2\Big)
=\lambda\langle\mathrm{tr}\pmb\epsilon\rangle_+\,\mathbf I.
\]
2. 第二项导数给出
\[
\frac{\partial}{\partial\pmb\epsilon}(\mu\,\pmb\epsilon^+:\pmb\epsilon^+)=2\mu\,\pmb\epsilon^+.
\]

组合得 \[\pmb\sigma^+\]；\[\pmb\sigma^-\] 同理。

---

### 公式(2.1-g)：材料切线模量

\[
\frac{\partial \pmb{\sigma}}{\partial \pmb{\epsilon}} = [g(d) + k]\frac{\partial \pmb{\sigma}^{+}}{\partial \pmb{\epsilon}} + \frac{\partial \pmb{\sigma}^{-}}{\partial \pmb{\epsilon}} = [g(d) + k]\left[\lambda H(\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes \mathbf{I} + 2\mu \mathbb{P}^{+}\right] + \left[\lambda H(-\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes \mathbf{I} + 2\mu \mathbb{P}^{-}\right].
\]

**逐步推导框架**

1. 从上一式
\[
\pmb\sigma=[g(d)+k]\pmb\sigma^+ + \pmb\sigma^-.
\]
对 \[\pmb\epsilon\] 求导（\[d\] 固定）：
\[
\frac{\partial\pmb\sigma}{\partial\pmb\epsilon}=[g(d)+k]\frac{\partial\pmb\sigma^+}{\partial\pmb\epsilon}+\frac{\partial\pmb\sigma^-}{\partial\pmb\epsilon}.
\]

2. 对 \[\pmb\sigma^+=\lambda\langle\mathrm{tr}\pmb\epsilon\rangle_+\mathbf I+2\mu\pmb\epsilon^+\] 求导：
\[
\frac{\partial\pmb\sigma^+}{\partial\pmb\epsilon}=
\lambda H(\mathrm{tr}\pmb\epsilon)\,\mathbf I\otimes\mathbf I +2\mu\,\frac{\partial\pmb\epsilon^+}{\partial\pmb\epsilon}.
\]

3. 定义
\[
\mathbb P^+ := \frac{\partial\pmb\epsilon^+}{\partial\pmb\epsilon},\qquad
\mathbb P^- := \frac{\partial\pmb\epsilon^-}{\partial\pmb\epsilon},
\]
并代入即可。

**物理意义**

- 切线模量是牛顿/拟牛顿迭代的核心“刚度”；
- 其中 \[H(\cdot)\] 与 \[\mathbb P^{\pm}\] 体现了拉压切换与主方向投影，导致模型非线性。

---

### 公式(2.1-h)：四阶投影算子定义

\[
\mathbb{P}^{+} = \frac{\partial \pmb{\epsilon}^{+}}{\partial \pmb{\epsilon}}, \quad \mathbb{P}^{-} = \frac{\partial \pmb{\epsilon}^{-}}{\partial \pmb{\epsilon}}.
\]

**含义**

- \[\mathbb P^+\]：把应变增量映射到“正应变子空间”的线性算子；
- \[\mathbb P^-\]：对应负应变子空间。

**解读**

在主值重根附近其解析式较复杂（文中引用文献给出）。数值实现通常用谱分解安全处理，保证切线一致性。

---

### 公式(7)：总势能一阶变分（方向导数形式）

\[
\begin{array}{rl}
\delta \Pi (\pmb {u},d) &= D_{(\delta \pmb {u},\delta d)}\Pi (\pmb {u},d) = \left.\frac{\mathrm{d}}{\mathrm{d}\epsilon}\right|_{\epsilon = 0}\Pi (\pmb {u} + \epsilon \delta \pmb {u},d + \epsilon \delta d) \\
&= \int_{\Omega}\left(\frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial \pmb{\epsilon}}:\pmb{\epsilon}(\delta \pmb {u}) + \frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial d}\delta d\right)\mathrm{d}\Omega \\
&\quad +\int_{\Omega}\frac{g_{c}}{l}\left(d\delta d + l^{2}\nabla d\cdot \nabla \delta d\right)\mathrm{d}\Omega -\int_{\Omega}\pmb {b}\cdot \delta \pmb {u}\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot \delta \pmb {u}\mathrm{d}\Gamma \\
&= (\nabla^{(s)}\delta \pmb {u},\pmb {\sigma}) - (\delta \pmb {u},\pmb {b}) - (\delta \pmb {u},\pmb{t})_{\Gamma_{t}} + (\delta d,\frac{g_{c}}{l} d) + (\nabla \delta d,g_{c}l\nabla d) + (\delta d,g^{\prime}(d)\psi^{+}(\pmb{\epsilon})).
\end{array} \quad (7)
\]

**详细推导步骤**

1. 定义方向导数（Gateaux 导数）：
\[
\delta\Pi = \left.\frac{\mathrm d}{\mathrm d\epsilon}\right|_{0}\Pi(\pmb u+\epsilon\delta\pmb u,d+\epsilon\delta d).
\]

2. 对内部能项 \[\int_\Omega \psi(\pmb\epsilon(\pmb u),d)\,\mathrm d\Omega\] 用链式法则：
\[
\delta\psi = \frac{\partial\psi}{\partial\pmb\epsilon}:\delta\pmb\epsilon + \frac{\partial\psi}{\partial d}\,\delta d,
\quad \delta\pmb\epsilon=\nabla^{(s)}\delta\pmb u.
\]

3. 对裂纹正则项
\[
\int_\Omega \frac{g_c}{2l}\left(d^2+l^2|\nabla d|^2\right)\mathrm d\Omega
\]
变分得
\[
\int_\Omega \frac{g_c}{l}\left(d\,\delta d+l^2\nabla d\cdot\nabla\delta d\right)\mathrm d\Omega.
\]

4. 外力势能项线性，直接得到
\[
-\int_\Omega\pmb b\cdot\delta\pmb u\,\mathrm d\Omega
-\int_{\partial\Omega}\pmb t\cdot\delta\pmb u\,\mathrm d\Gamma.
\]

5. 使用内积记号
\[
(a,b):=\int_\Omega a\,b\,\mathrm d\Omega,
\quad (a,b)_{\Gamma_t}:=\int_{\Gamma_t}a\,b\,\mathrm d\Gamma
\]
可写成最后一行紧凑形式。

**物理意义**

- \[\delta\Pi=0\] 是平衡条件（驻值条件）；
- 位移测试方向给力学平衡；相场测试方向给损伤演化平衡。

---

### 公式(2.1-i)：弱式残量方程组

\[
\left\{ \begin{array}{ll} r_{\pmb{u}}(\pmb{u},d) = (\nabla^{(s)}\delta \pmb{u},\pmb{\sigma}) - (\delta \pmb{u},\pmb{b}) - (\delta \pmb{u},\pmb{t})_{\Gamma_{t}} = 0, \\ r_{d}(\pmb{u},d) = (\delta d,\frac{g_{c}}{l} d) + (\nabla \delta d,g_{c}l\nabla d) + (\delta d,g^{\prime}(d)\psi^{+}(\pmb{\epsilon})) = 0, \end{array} \right.
\]

**来源**

把公式(7)中与 \[\delta\pmb u\]、\[\delta d\] 分别对应的项分组，即得两组弱式残量。

**含义**

- \[r_{\pmb u}=0\]：动量平衡的弱形式；
- \[r_d=0\]：相场控制方程的弱形式；
- 还需附加不可逆约束 \[d^{(n)}\le d\le 1\]。

---

## 2.2 Finite element discretization

### 公式(2.2-a)：位移与相场的有限元插值

\[
\pmb {u} = \pmb{N}_{u_{A}}\pmb{u}_{A}, \quad \text{and} \quad d = N_{d_{A}}d_{A}.
\]

**逐项含义**

- \[\pmb N_{u_A}\]：结点 \[A\] 处位移向量形函数。
- \[N_{d_A}\]：结点 \[A\] 处相场标量形函数。
- \[\pmb u_A,d_A\]：对应结点自由度。
- 采用 Einstein 求和约定，对重复下标 \[A\] 自动求和。

**意义**

将无限维场变量近似到有限维结点自由度空间，是离散化第一步。

---

### 公式(2.2-b)：测试函数（变分）插值

\[
\delta \pmb {u} = \pmb{N}_{u_{A}}\delta \pmb{u}_{A}, \quad \text{and} \quad \delta d = N_{d_{A}}\delta d_{A},
\]

**含义**

- 与主变量同一离散空间（Galerkin 思想）；
- 保证离散变分方程与离散能量导数一致。

---

### 公式(8)：离散后的目标函数

\[
\begin{array}{rl}
\Pi (\pmb {u}_{A},d_{A}) = &\int_{\Omega}\psi \left(\pmb{\epsilon}(\pmb {N}_{u_{A}}\pmb {u}_{A}),N_{d_{A}}d_{A}\right)\mathrm{d}\Omega \\
&+\int_{\Omega}\frac{g_{c}}{2l}\left((N_{d_{A}}d_{A})^{2} + l^{2}(\nabla N_{d_{A}}d_{A})\cdot (\nabla N_{d_{A}}d_{A})\right)\mathrm{d}\Omega \\
&-\int_{\Omega}\pmb {b}\cdot (\pmb {N}_{u_{A}}\pmb {u}_{A})\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot (\pmb {N}_{u_{A}}\pmb {u}_{A})\mathrm{d}\Gamma,
\end{array} \quad (8)
\]

**推导来源**

把公式(2.2-a) 的插值关系代入连续总势能（论文 Eq.(1)）即可得到该离散目标函数。

**含义**

- 第一项：离散弹性能；
- 第二项：相场裂纹表面正则化能；
- 后两项：体力与边界力外力势。

该式是后续优化问题 \[\min\Pi\] 的核心。

---

### 公式(9)：离散能量梯度（残量向量）

\[
\begin{array}{rl}
\pmb{r} = \nabla \Pi = (r_{\pmb{u}},r_{d})^{\mathrm{T}},\\
r_{\pmb{u}_{A}} = \frac{\partial\Pi}{\partial \pmb{u}_{A}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},\pmb{\sigma}\right) - \left(\pmb{N}_{u_{A}},\pmb{b}\right) - \left(\pmb{N}_{u_{A}},\pmb{t}\right)_{\Gamma_{t}},\\
r_{d_{A}} = \frac{\partial\Pi}{\partial d_{A}} = \left(N_{d_{A}},\frac{g_{c}}{l} d + g^{\prime}(d)\psi^{+}\right) + (\nabla N_{d_{A}},g_{c}l\nabla d).
\end{array} \quad (9)
\]

**必要推导**

1. 按自由度求偏导：
\[
r_{\pmb u_A}=\frac{\partial\Pi}{\partial\pmb u_A},\qquad r_{d_A}=\frac{\partial\Pi}{\partial d_A}.
\]

2. 位移部分：
\[
\frac{\partial\pmb\epsilon(\pmb N_{u_B}\pmb u_B)}{\partial\pmb u_A}=\nabla^{(s)}\pmb N_{u_A}
\]
代入能量一阶导并减去外力项，得到 \[r_{\pmb u_A}\]。

3. 相场部分：

- 来自 \[\psi\] 的相场导数项
\[
\frac{\partial\psi}{\partial d}=g'(d)\psi^+,
\]
- 来自 \[\frac{g_c}{2l}d^2\] 的导数给 \[\frac{g_c}{l}d\]；
- 来自梯度项给 \[(\nabla N_{d_A}, g_c l\nabla d)\]。

组合即得 \[r_{d_A}\]。

**数值意义**

- \[\pmb r=\mathbf 0\] 即离散平衡方程；
- 优化/牛顿迭代中，\[\pmb r\] 是一阶最优性残差。

---

### 公式(2.2-c)：Hessian 分块结构

\[
\mathbf{K} = \nabla^{2}\Pi = \left[ \begin{array}{ll}\mathbf{K}_{uu} & \mathbf{K}_{ud}\\ \mathbf{K}_{du} & \mathbf{K}_{dd} \end{array} \right],
\]

**含义**

- \[\mathbf K\] 是离散能量二阶导（切线刚度/雅可比）。
- 四个分块分别对应
  - 位移-位移耦合 \[uu\]；
  - 位移-相场耦合 \[ud\]；
  - 相场-位移耦合 \[du\]；
  - 相场-相场 \[dd\]。

**物理解释**

耦合块 \[\mathbf K_{ud},\mathbf K_{du}\] 表示“损伤影响力学刚度”与“力学状态反馈损伤驱动力”的双向耦合。

---

### 公式(10)：Hessian 各分块显式表达

\[
\begin{array}{rl}
\mathbf{K}_{u_{A}u_{B}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},\frac{\partial \pmb{\sigma}}{\partial \pmb{\epsilon}}:\nabla^{(s)}\pmb{N}_{u_{B}}\right),\qquad \mathbf{K}_{u_{A}d_{B}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},g^{\prime}(d)\pmb{\sigma}^{+}N_{d_{B}}\right),\\
\mathbf{K}_{d_{A}u_{B}} = \left(N_{d_{A}},g^{\prime}(d)\pmb{\sigma}^{+}:\nabla^{(s)}\pmb{N}_{u_{B}}\right),\qquad \mathbf{K}_{d_{A}d_{B}} = \left(N_{d_{A}},\left(\frac{g_{c}}{l} +g^{\prime \prime}(d)\psi^{+}\right)N_{d_{B}}\right) + \left(\nabla N_{d_{A}},g_{c}l\nabla N_{d_{B}}\right).
\end{array} \quad (10)
\]

**逐块推导要点**

1. \[\mathbf K_{u_Au_B}=\partial r_{u_A}/\partial u_B\]

由公式(9)中 \[r_{u_A}\] 对 \[u_B\] 求导，核心是
\[
\frac{\partial\pmb\sigma}{\partial\pmb\epsilon}:\nabla^{(s)}\pmb N_{u_B}.
\]

2. \[\mathbf K_{u_Ad_B}=\partial r_{u_A}/\partial d_B\]

由于 \[\pmb\sigma=[g(d)+k]\pmb\sigma^+ + \pmb\sigma^-\]，对 \[d\] 求导仅第一项有贡献：
\[
\frac{\partial\pmb\sigma}{\partial d}=g'(d)\pmb\sigma^+.
\]
再乘 \[N_{d_B}\] 得上式。

3. \[\mathbf K_{d_Au_B}=\partial r_{d_A}/\partial u_B\]

\[r_{d_A}\] 中含 \[\psi^+(\pmb\epsilon)\]，对 \[u_B\] 求导：
\[
\frac{\partial\psi^+}{\partial\pmb\epsilon}:\nabla^{(s)}\pmb N_{u_B}=\pmb\sigma^+:\nabla^{(s)}\pmb N_{u_B}.
\]
再乘前因子 \[g'(d)N_{d_A}\]。

4. \[\mathbf K_{d_Ad_B}=\partial r_{d_A}/\partial d_B\]

- \[\frac{g_c}{l}d\] 对 \[d_B\] 导数给 \[\frac{g_c}{l}N_{d_B}\]；
- \[g'(d)\psi^+\] 对 \[d_B\] 导数给 \[g''(d)\psi^+N_{d_B}\]（此处 \[\psi^+\] 对 \[d\] 不显含）；
- 梯度项导数给 \[(\nabla N_{d_A},g_c l\nabla N_{d_B})\]。

**数值意义**

该矩阵是牛顿/拟牛顿与预条件器构造的重要依据，也是单元刚度装配的直接模板。

---

### 公式(2.2-d)：相场不可逆箱约束（结点级）

\[
d_A^{(n)}\leq d_A\leq 1,
\]

**含义**

- 下界 \[d_A^{(n)}\]：当前时间步起点的已收敛损伤，不允许“愈合”；
- 上界 \[1\]：相场变量的物理上限。

**物理解释**

裂纹不可逆：一旦损伤形成，不应在后续加载/卸载中自发恢复到更小损伤。

---

### 公式(11)：仅保留对角耦合块的块对角矩阵

\[
\hat{\mathbf{K}} = \left[ \begin{array}{cc}\mathbf{K}_{uu} & \mathbf{0}\\ \mathbf{0} & \mathbf{K}_{dd} \end{array} \right], \quad (11)
\]

**含义与用途**

- 将完整 Hessian 的耦合块 \[\mathbf K_{ud},\mathbf K_{du}\] 暂时忽略，保留两个主对角块；
- 在后续单体化（monolithic）算法中常用于构造简化预条件、初始近似或子问题矩阵。

**数值直觉**

- 优点：更易求解、存储更友好；
- 代价：忽略耦合可能降低单步逼近精度，但整体可通过迭代修正。

---

## 2章公式之间的逻辑闭环（总览）

为避免“只看单式不看联系”，这里把 2 章全部公式按求解流程串起来：

1. **本构与退化定义**：公式(5)、(6)、(2.1-a)~(2.1-d)给出能量分解与拉压区分。  
2. **由能量到应力/切线**：公式(2.1-e)~(2.1-h)给出 \[\pmb\sigma\] 与 \[\partial\pmb\sigma/\partial\pmb\epsilon\]。  
3. **由总势能到弱式**：公式(7)、(2.1-i)给出连续层面的平衡方程。  
4. **有限元离散**：公式(2.2-a)、(2.2-b)、(8)把问题变成有限维优化。  
5. **一阶与二阶离散导数**：公式(9)、(2.2-c)、(10)形成残量与切线矩阵。  
6. **不可逆约束与算法矩阵**：公式(2.2-d)、(11)进入后续梯度投影单体化求解框架。

---

## 符号表（按第2章出现顺序）

- \[\Omega\]：计算域；\[\partial\Omega\]：边界；\[\Gamma_t\]：受力边界。  
- \[\pmb u\]：位移场；\[d\]：相场（损伤）变量。  
- \[\delta\pmb u,\delta d\]：对应测试函数/变分方向。  
- \[\pmb\epsilon(\pmb u)=\nabla^{(s)}\pmb u\]：小应变张量。  
- \[\psi,\psi^+,\psi^-\]：总/正/负应变能密度。  
- \[g(d)\]：退化函数；\[g'(d),g''(d)\]：一、二阶导。  
- \[\pmb\sigma,\pmb\sigma^+,\pmb\sigma^-\]：总/正/负应力。  
- \[\lambda,\mu\]：Lamé 参数。  
- \[\mathbb P^+,\mathbb P^-\]：正负应变投影四阶张量。  
- \[g_c\]：断裂韧度；\[l\]：相场长度尺度。  
- \[\pmb b\]：体力；\[\pmb t\]：边界力。  
- \[\pmb r\]：离散残量；\[\mathbf K\]：离散 Hessian（切线刚度）。  
- \[A,B\]：有限元结点编号。  
- \[d_A^{(n)}\]：时间步 \[n\] 已收敛相场值（下界）。

---

## 完整性核对（第2章无遗漏）

本文件已覆盖 `explain.md` 第2章中全部公式类型：

- 编号公式：\[(5)\]、\[(6)\]、\[(7)\]、\[(8)\]、\[(9)\]、\[(10)\]、\[(11)\]；
- 非编号但关键定义式：算子定义、谱分解、\[\pmb\epsilon^{\pm}\]、\[\psi^{\pm}\]、\[\pmb\sigma^{\pm}\]、切线模量式、\[\mathbb P^{\pm}\] 定义、弱式系统、插值与变分插值、Hessian 分块、不可逆约束。

