# explain.md 第二章节全部公式逐一解释（完整版）

## 范围说明

- 本文档只解读 `explain.md` 中 **第 2 章节**（`## 2. Phase-field formulation and finite element discretization`）出现的全部公式与符号。
- 每个条目均按以下顺序组织：
  1. 先完整显示公式（可渲染）；
  2. 解释各部分符号含义；
  3. 给出公式整体含义；
  4. 给出必要推导过程。
- 公式渲染统一采用 `\[ ... \]`（块公式）与 `\( ... \)`（行内公式）。

---

## 2.1 Phase-field formulation

### 公式 2.1-1（原文式 (5)）

\[
\psi(\pmb{\epsilon}, d) = [g(d) + k] \, \psi^{+}(\pmb{\epsilon}) + \psi^{-}(\pmb{\epsilon}), \quad (5)
\]

**各部分含义：**

- \(\psi(\pmb{\epsilon}, d)\)：总应变能密度（受损伤变量 \(d\) 影响）。
- \(\pmb{\epsilon}\)：小应变张量，通常 \(\pmb{\epsilon}=\nabla^{(s)}\pmb{u}\)。
- \(d\)：相场/损伤变量，\(d=0\) 未损伤，\(d=1\) 完全损伤。
- \(\psi^+(\pmb{\epsilon})\)：拉伸（正）部分应变能密度。
- \(\psi^-(\pmb{\epsilon})\)：压缩（负）部分应变能密度。
- \(g(d)\)：退化函数（损伤对刚度或能量的削弱）。
- \(k\)：小的非负参数（避免数值奇异，本文示例取 \(k=0\)）。

**公式含义：**

该式表达“拉压不对称”断裂建模思想：只让拉伸能量被损伤退化，压缩能量保持不退化，从而避免裂缝在压缩下不合理扩展。

**必要推导过程：**

1. 从完整弹性能 \(\psi_0(\pmb\epsilon)\) 出发，先做谱分解并拆成正负两部分：
   \(\psi_0=\psi^+ + \psi^-\)。
2. 引入损伤变量 \(d\) 后，仅对 \(\psi^+\) 乘以退化因子 \(g(d)\)（再加 \(k\) 稳定项）：
   \(\psi = [g(d)+k]\psi^+ + \psi^-\)。
3. 这一步是建模假设（constitutive assumption），其目的不是纯数学推导唯一得到，而是满足断裂的物理合理性（张拉开裂、压缩闭合）。

---

### 公式 2.1-2（原文式 (6)）

\[
g(d) = (1 - d)^{2}. \quad (6)
\]

**各部分含义：**

- \(g(d)\)：退化函数。
- \(d\in[0,1]\)：损伤变量。

**公式含义：**

二次型退化函数满足典型要求：
- \(g(0)=1\)：无损时不退化；
- \(g(1)=0\)：全损时拉伸承载能力消失；
- 连续可导，便于牛顿/拟牛顿线性化。

**必要推导过程：**

这属于“函数选型”而非由前式唯一推导。常见约束：
\[
g(0)=1,\quad g(1)=0,\quad g'(d)\le 0.
\]
最简单且光滑的多项式选型即 \((1-d)^2\)。

---

### 公式 2.1-3（正负括号与 Heaviside 算子）

\[
\langle x\rangle_{+} = \frac{1}{2}(x + |x|), \quad
\langle x\rangle_{-} = \frac{1}{2}(x - |x|), \quad
H(x) = \begin{cases} 1 & x \geq 0, \\ 0 & x < 0. \end{cases}
\]

**各部分含义：**

- \(\langle x\rangle_+\)：取 \(x\) 的正部分。
- \(\langle x\rangle_-\)：取 \(x\) 的负部分。
- \(H(x)\)：Heaviside 阶跃函数，在线性化中判断“正/负激活区间”。

**公式含义：**

这是拉压分解的基础算子：将标量拆成正负两部分，便于作用到应变主值及迹上。

**必要推导过程：**

利用绝对值分段定义：
- 若 \(x\ge0\)，\(|x|=x\)，则 \(\langle x\rangle_+=x,\langle x\rangle_-=0\)；
- 若 \(x<0\)，\(|x|=-x\)，则 \(\langle x\rangle_+=0,\langle x\rangle_-=x\)。
故可统一写成上述代数形式。

---

### 公式 2.1-4（应变谱分解）

\[
\pmb{\epsilon} = \sum_{\alpha} \epsilon_{\alpha} \mathbf{M}_{\alpha}, \quad
\mathbf{M}_{\alpha} = \pmb{n}_{\alpha} \otimes \pmb{n}_{\alpha},
\]

**各部分含义：**

- \(\epsilon_\alpha\)：应变张量的第 \(\alpha\) 个主应变（特征值）。
- \(\pmb n_\alpha\)：对应主方向特征向量。
- \(\mathbf M_\alpha\)：谱投影张量（二阶）。

**公式含义：**

将对称应变张量写成主值-主方向基的和，为后续“只退化拉伸主应变”提供可操作表达。

**必要推导过程：**

对称张量 \(\pmb\epsilon\) 可正交对角化：
\[
\pmb\epsilon = \mathbf Q\,\mathrm{diag}(\epsilon_1,\epsilon_2,\epsilon_3)\,\mathbf Q^T.
\]
令 \(\pmb n_\alpha\) 为 \(\mathbf Q\) 列向量，即得谱分解形式。

---

### 公式 2.1-5（正负应变张量）

\[
\pmb{\epsilon}^{+} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{+} \mathbf{M}_{\alpha}, \quad
\pmb{\epsilon}^{-} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{-} \mathbf{M}_{\alpha}.
\]

**各部分含义：**

- \(\pmb\epsilon^+\)：由正主应变组装的“拉伸应变部分”。
- \(\pmb\epsilon^-\)：由负主应变组装的“压缩应变部分”。

**公式含义：**

在主方向上逐项截断（positive/negative part），把应变张量拆成两部分，供能量与应力的拉压分裂构造使用。

**必要推导过程：**

由上一个谱分解：
\[
\pmb\epsilon = \sum_\alpha \epsilon_\alpha \mathbf M_\alpha
= \sum_\alpha \left(\langle\epsilon_\alpha\rangle_+ + \langle\epsilon_\alpha\rangle_-\right)\mathbf M_\alpha
= \pmb\epsilon^+ + \pmb\epsilon^-.
\]

---

### 公式 2.1-6（正负应变能密度）

\[
\psi^{+}(\pmb{\epsilon}) = \frac{1}{2}\lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{+}^{2} + \mu \, \pmb{\epsilon}^{+} : \pmb{\epsilon}^{+}, \quad
\psi^{-}(\pmb{\epsilon}) = \frac{1}{2} \lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{-}^{2} + \mu \, \pmb{\epsilon}^{-} : \pmb{\epsilon}^{-},
\]

**各部分含义：**

- \(\lambda,\mu\)：Lamé 常数。
- \(\mathrm{tr}\pmb\epsilon\)：体积应变（迹）。
- \(\mathbf A:\mathbf B\)：双点积。

**公式含义：**

把线弹性能 \(\frac{\lambda}{2}(\mathrm{tr}\epsilon)^2+\mu\epsilon:\epsilon\) 拆为拉伸部分与压缩部分，使得损伤只作用于拉伸能。

**必要推导过程：**

1. 经典各向同性线弹性能：
\[
\psi_0(\epsilon)=\frac{\lambda}{2}(\mathrm{tr}\epsilon)^2+\mu\,\epsilon:\epsilon.
\]
2. 将体积项和偏量项分别做正负拆分：
\((\mathrm{tr}\epsilon)^2\to \langle\mathrm{tr}\epsilon\rangle_+^2 + \langle\mathrm{tr}\epsilon\rangle_-^2\)，
\(\epsilon\to \epsilon^+,\epsilon^-\)。
3. 定义上式中的 \(\psi^+,\psi^-\)，并满足 \(\psi_0=\psi^+ + \psi^-\)（在该分裂定义下）。

---

### 公式 2.1-7（总应力分裂表达）

\[
\pmb{\sigma} = \frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial \pmb{\epsilon}}
= [g(d) + k]\frac{\partial\psi^{+}(\pmb{\epsilon})}{\partial \pmb{\epsilon}} + \frac{\partial\psi^{-}(\pmb{\epsilon})}{\partial \pmb{\epsilon}}
= [g(d) + k]\pmb{\sigma}^{+} + \pmb{\sigma}^{-},
\]

**各部分含义：**

- \(\sigma\)：Cauchy 应力（小变形下可直接由能量对应变导出）。
- \(\sigma^+=\partial\psi^+/\partial\epsilon\)，\(\sigma^- = \partial\psi^-/\partial\epsilon\)。

**公式含义：**

应力同样按拉压两部分分裂；拉伸应力乘以退化因子，压缩应力不退化。

**必要推导过程：**

把公式 (5) 对 \(\epsilon\) 求偏导：
\[
\frac{\partial\psi}{\partial\epsilon}
=\frac{\partial\left(([g(d)+k]\psi^+)+\psi^-\right)}{\partial\epsilon}
=[g(d)+k]\frac{\partial\psi^+}{\partial\epsilon}+\frac{\partial\psi^-}{\partial\epsilon},
\]
因为 \(g(d)+k\) 对 \(\epsilon\) 视作常数（\(d\) 在此偏导中固定）。

---

### 公式 2.1-8（正负应力显式形式）

\[
\pmb{\sigma}^{+} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{+}\mathbf{I} + 2\mu \, \pmb{\epsilon}^{+}, \quad
\pmb{\sigma}^{-} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{-}\mathbf{I} + 2\mu \, \pmb{\epsilon}^{-}.
\]

**各部分含义：**

- \(\mathbf I\)：二阶单位张量。
- 第一项：体积响应；第二项：偏量/剪切响应。

**公式含义：**

正负应力分别由正负体积应变与正负应变张量构成，是线弹性本构在拉压分裂框架下的对应形式。

**必要推导过程：**

对上一节的 \(\psi^+,\psi^-\) 分别求导：
\[
\frac{\partial}{\partial\epsilon}\left(\frac{\lambda}{2}\langle\mathrm{tr}\epsilon\rangle_\pm^2\right)
= \lambda\langle\mathrm{tr}\epsilon\rangle_\pm\,\mathbf I,
\]
\[
\frac{\partial}{\partial\epsilon}(\mu\,\epsilon^\pm:\epsilon^\pm)
=2\mu\,\epsilon^\pm,
\]
组合即可得到 \(\sigma^\pm\)。

---

### 公式 2.1-9（材料切线模量）

\[
\frac{\partial \pmb{\sigma}}{\partial \pmb{\epsilon}}
= [g(d) + k]\frac{\partial \pmb{\sigma}^{+}}{\partial \pmb{\epsilon}} + \frac{\partial \pmb{\sigma}^{-}}{\partial \pmb{\epsilon}}
= [g(d) + k]\left[\lambda H(\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes \mathbf{I} + 2\mu \mathbb{P}^{+}\right]
+ \left[\lambda H(-\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes \mathbf{I} + 2\mu \mathbb{P}^{-}\right].
\]

**各部分含义：**

- \(\partial\sigma/\partial\epsilon\)：四阶切线刚度（consistent tangent）。
- \(\mathbf I\otimes\mathbf I\)：体积投影相关四阶算子。
- \(\mathbb P^+,\mathbb P^-\)：正/负应变张量对应变的导数（四阶投影张量）。
- \(H(\cdot)\)：决定当前点处于拉伸还是压缩主导区间。

**公式含义：**

该式是牛顿或拟牛顿迭代组装刚度矩阵时的核心本构线性化；由于拉压分裂与退化耦合，本构是分段非线性的。

**必要推导过程：**

1. 从 \(\sigma=[g(d)+k]\sigma^+ + \sigma^-\) 出发，对 \(\epsilon\) 求导。
2. \([g(d)+k]\) 与 \(\epsilon\) 无关（该偏导下），可直接提出。
3. 对 \(\sigma^\pm\) 求导时：
   - \(\partial\langle\mathrm{tr}\epsilon\rangle_+/\partial\epsilon = H(\mathrm{tr}\epsilon)\,\mathbf I\)；
   - \(\partial\langle\mathrm{tr}\epsilon\rangle_-/\partial\epsilon = H(-\mathrm{tr}\epsilon)\,\mathbf I\)；
   - \(\partial\epsilon^\pm/\partial\epsilon = \mathbb P^\pm\)。
4. 代回并整理得到最终表达。

---

### 公式 2.1-10（投影张量定义）

\[
\mathbb{P}^{+} = \frac{\partial \pmb{\epsilon}^{+}}{\partial \pmb{\epsilon}}, \quad
\mathbb{P}^{-} = \frac{\partial \pmb{\epsilon}^{-}}{\partial \pmb{\epsilon}},
\]

**各部分含义：**

- \(\mathbb P^\pm\)：四阶张量，描述“总应变微扰”如何投影到“正/负应变部分微扰”。

**公式含义：**

它们是分裂本构里最关键的线性化对象，决定切线刚度在谱分裂下的精确形式。

**必要推导过程：**

由定义直接得出：
\(\epsilon^\pm=\epsilon^\pm(\epsilon)\) 是 \(\epsilon\) 的函数，线性化时
\(\mathrm d\epsilon^\pm = \mathbb P^\pm : \mathrm d\epsilon\)，故有上式定义。

---

### 公式 2.1-11（原文式 (7)：能量一阶变分）

\[
\begin{array}{rl}
\delta \Pi (\pmb {u},d) &= D_{(\delta \pmb {u},\delta d)}\Pi (\pmb {u},d) = \left.\frac{\mathrm{d}}{\mathrm{d}\epsilon}\right|_{\epsilon = 0}\Pi (\pmb {u} + \epsilon \delta \pmb {u},d + \epsilon \delta d) \\
&= \int_{\Omega}\left(\frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial \pmb{\epsilon}}:\pmb{\epsilon}(\delta \pmb {u}) + \frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial d}\delta d\right)\mathrm{d}\Omega \\
&\quad +\int_{\Omega}\frac{g_{c}}{l}\left(d\delta d + l^{2}\nabla d\cdot \nabla \delta d\right)\mathrm{d}\Omega -\int_{\Omega}\pmb {b}\cdot \delta \pmb {u}\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot \delta \pmb {u}\mathrm{d}\Gamma \\
&= (\nabla^{(s)}\delta \pmb {u},\pmb {\sigma}) - (\delta \pmb {u},\pmb {b}) - (\delta \pmb {u},\pmb{t})_{\Gamma_{t}} + (\delta d,\frac{g_{c}}{l} d) + (\nabla \delta d,g_{c}l\nabla d) + (\delta d,g^{\prime}(d)\psi^{+}(\pmb{\epsilon})).
\end{array} \quad (7)
\]

**各部分含义：**

- \(\delta\Pi\)：总势能在方向 \((\delta\mathbf u,\delta d)\) 上的一阶方向导数。
- \((\cdot,\cdot)\)：\(L^2\) 内积，\((\cdot,\cdot)_{\Gamma_t}\) 为边界内积。
- \(g_c\)：断裂韧度（临界能量释放率）。
- \(l\)：相场长度尺度。

**公式含义：**

该式给出弱式残量来源：令 \(\delta\Pi=0\) 就得到位移平衡方程残量和相场演化方程残量。

**必要推导过程：**

1. 从总能量函数 \(\Pi(u,d)\)（文中式 (1)+(2)）出发，对扰动 \(u+\epsilon\delta u,\ d+\epsilon\delta d\) 求导。
2. 使用链式法则：
   - \(\mathrm d\psi = (\partial\psi/\partial\epsilon):\epsilon(\delta u) + (\partial\psi/\partial d)\delta d\)；
   - 其中 \(\partial\psi/\partial\epsilon=\sigma\)。
3. 裂纹表面密度项 \(\frac{g_c}{2l}(d^2+l^2|\nabla d|^2)\) 的变分：
\[
\delta\left(\frac{g_c}{2l}d^2\right)=\frac{g_c}{l}d\,\delta d,
\quad
\delta\left(\frac{g_c l}{2}|\nabla d|^2\right)=g_c l\,\nabla d\cdot\nabla\delta d.
\]
4. 外力势项对 \(d\) 无贡献，仅对 \(u\) 贡献负号项。
5. 把积分表达改写成内积记号，得到式 (7) 的最后一行。

---

### 公式 2.1-12（弱式方程组）

\[
\left\{ \begin{array}{ll}
 r_{\pmb{u}}(\pmb{u},d) = (\nabla^{(s)}\delta \pmb{u},\pmb{\sigma}) - (\delta \pmb{u},\pmb{b}) - (\delta \pmb{u},\pmb{t})_{\Gamma_{t}} = 0, \\
 r_{d}(\pmb{u},d) = (\delta d,\frac{g_{c}}{l} d) + (\nabla \delta d,g_{c}l\nabla d) + (\delta d,g^{\prime}(d)\psi^{+}(\pmb{\epsilon})) = 0,
\end{array} \right.
\]

**各部分含义：**

- \(r_u\)：位移方程残量。
- \(r_d\)：相场方程残量。
- 测试函数分别是 \(\delta\mathbf u\) 与 \(\delta d\)。

**公式含义：**

这是变分原理 \(\delta\Pi=0\) 对任意测试函数成立的等价弱式系统。

**必要推导过程：**

由上一个式 (7) 的项按 \(\delta\mathbf u\) 与 \(\delta d\) 分组：
\[
\delta\Pi = r_u(\delta\mathbf u) + r_d(\delta d).
\]
因为测试函数任意，分别令两组系数泛函为零，即得两条弱式方程。

---

## 2.2 Finite element discretization

### 公式 2.2-1（位移与相场插值）

\[
\pmb {u} = \pmb{N}_{u_{A}}\pmb{u}_{A}, \quad \text{and} \quad d = N_{d_{A}}d_{A}.
\]

**各部分含义：**

- \(\pmb N_{u_A}\)：位移场向量形函数矩阵。
- \(N_{d_A}\)：相场标量形函数。
- \(\pmb u_A, d_A\)：节点自由度。
- 默认使用 Einstein 求和（重复下标求和）。

**公式含义：**

把连续场投影到有限元离散空间，得到以节点自由度为未知量的近似表示。

**必要推导过程：**

有限元近似基本假设：在单元上用形函数线性组合近似未知场。向量场与标量场形式分别如上。

---

### 公式 2.2-2（变分场插值）

\[
\delta \pmb {u} = \pmb{N}_{u_{A}}\delta \pmb{u}_{A}, \quad \text{and} \quad \delta d = N_{d_{A}}\delta d_{A},
\]

**各部分含义：**

- \(\delta \pmb u_A,\delta d_A\)：节点层面的虚位移/虚相场增量。

**公式含义：**

测试函数与试函数取同一离散空间（Galerkin 框架）。

**必要推导过程：**

由 Galerkin 法，变分空间与近似空间一致，故对变分量使用同形函数展开。

---

### 公式 2.2-3（原文式 (8)：离散后总势能）

\[
\begin{array}{rl}
\Pi (\pmb {u}_{A},d_{A}) = &\int_{\Omega}\psi \left(\pmb{\epsilon}(\pmb {N}_{u_{A}}\pmb {u}_{A}),N_{d_{A}}d_{A}\right)\mathrm{d}\Omega \\
&+\int_{\Omega}\frac{g_{c}}{2l}\left((N_{d_{A}}d_{A})^{2} + l^{2}(\nabla N_{d_{A}}d_{A})\cdot (\nabla N_{d_{A}}d_{A})\right)\mathrm{d}\Omega \\
&-\int_{\Omega}\pmb {b}\cdot (\pmb {N}_{u_{A}}\pmb {u}_{A})\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot (\pmb {N}_{u_{A}}\pmb {u}_{A})\mathrm{d}\Gamma,
\end{array} \quad (8)
\]

**各部分含义：**

- 该式已完全转为节点未知量 \(\{\pmb u_A,d_A\}\) 的函数。
- 四项依次为：弹性能、裂纹表面正则化能、体力势、表面力势。

**公式含义：**

把连续泛函离散成有限维优化目标，为后续求梯度（残量）和 Hessian（切线刚度）奠定基础。

**必要推导过程：**

1. 将 \(u,d\) 与 \(\delta u,\delta d\) 的形函数展开代入连续能量式。
2. 用 \(\epsilon(u)=\nabla^{(s)}u\) 把应变写成节点自由度的线性映射。
3. 保持积分形式（实际实现中再用数值积分求值）。

---

### 公式 2.2-4（原文式 (9)：离散梯度/残量）

\[
\begin{array}{rl}
\pmb{r} = \nabla \Pi = (r_{\pmb{u}},r_{d})^{\mathrm{T}},\\
r_{\pmb{u}_{A}} = \frac{\partial\Pi}{\partial \pmb{u}_{A}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},\pmb{\sigma}\right) - \left(\pmb{N}_{u_{A}},\pmb{b}\right) - \left(\pmb{N}_{u_{A}},\pmb{t}\right)_{\Gamma_{t}},\\
r_{d_{A}} = \frac{\partial\Pi}{\partial d_{A}} = \left(N_{d_{A}},\frac{g_{c}}{l} d + g^{\prime}(d)\psi^{+}\right) + (\nabla N_{d_{A}},g_{c}l\nabla d).
\end{array} \quad (9)
\]

**各部分含义：**

- \(\mathbf r\)：离散残量向量（优化中即目标函数梯度）。
- \(r_{u_A}\)：第 \(A\) 个位移自由度对应分量。
- \(r_{d_A}\)：第 \(A\) 个相场自由度对应分量。

**公式含义：**

离散平衡条件就是 \(\mathbf r=\mathbf 0\)。这也是非线性迭代的主方程。

**必要推导过程：**

1. 对离散总能量 \(\Pi(\mathbf u_A,d_A)\) 分别对 \(\mathbf u_A\) 与 \(d_A\) 求偏导。
2. 运用链式法则和内积导数规则：
   - \(\partial\epsilon/\partial\mathbf u_A = \nabla^{(s)}\mathbf N_{u_A}\)；
   - \(\partial d/\partial d_A = N_{d_A}\)。
3. 将连续弱式中测试函数替换成对应形函数，即得到离散残量表达。

---

### 公式 2.2-5（Hessian 分块结构）

\[
\mathbf{K} = \nabla^{2}\Pi = \left[ \begin{array}{ll}\mathbf{K}_{uu} & \mathbf{K}_{ud}\\ \mathbf{K}_{du} & \mathbf{K}_{dd} \end{array} \right],
\]

**各部分含义：**

- \(\mathbf K\)：离散 Hessian/切线刚度矩阵。
- \(\mathbf K_{uu}\)：位移-位移块。
- \(\mathbf K_{dd}\)：相场-相场块。
- \(\mathbf K_{ud},\mathbf K_{du}\)：机-损耦合块。

**公式含义：**

反映单元组装后的耦合结构：位移与相场是强耦合非线性系统。

**必要推导过程：**

由定义 \(\mathbf K=\partial\mathbf r/\partial\mathbf x\)（\(\mathbf x=[\mathbf u;\mathbf d]\)），按变量块划分即可得到该分块矩阵。

---

### 公式 2.2-6（原文式 (10)：各块矩阵显式表达）

\[
\begin{array}{rl}
\mathbf{K}_{u_{A}u_{B}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},\frac{\partial \pmb{\sigma}}{\partial \pmb{\epsilon}}:\nabla^{(s)}\pmb{N}_{u_{B}}\right),\qquad
\mathbf{K}_{u_{A}d_{B}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},g^{\prime}(d)\pmb{\sigma}^{+}N_{d_{B}}\right),\\
\mathbf{K}_{d_{A}u_{B}} = \left(N_{d_{A}},g^{\prime}(d)\pmb{\sigma}^{+}:\nabla^{(s)}\pmb{N}_{u_{B}}\right),\qquad
\mathbf{K}_{d_{A}d_{B}} = \left(N_{d_{A}},\left(\frac{g_{c}}{l} +g^{\prime \prime}(d)\psi^{+}\right)N_{d_{B}}\right) + \left(\nabla N_{d_{A}},g_{c}l\nabla N_{d_{B}}\right).
\end{array} \quad (10)
\]

**各部分含义：**

- 下标 \(A,B\)：节点编号。
- \(g'(d),g''(d)\)：退化函数的一阶、二阶导数。
- \(\partial\sigma/\partial\epsilon\)：一致切线本构。

**公式含义：**

该式给出牛顿线性化后单元刚度组装所需的全部块项，是算法实现最直接的离散表达。

**必要推导过程：**

1. 对式 (9) 的两个残量分别再求导：
   - \(K_{uu}=\partial r_u/\partial u\)，
   - \(K_{ud}=\partial r_u/\partial d\)，
   - \(K_{du}=\partial r_d/\partial u\)，
   - \(K_{dd}=\partial r_d/\partial d\)。
2. 例如 \(K_{ud}\)：
\[
\frac{\partial}{\partial d_B}(\nabla^sN_{u_A},\sigma)
= (\nabla^sN_{u_A},\frac{\partial\sigma}{\partial d}N_{d_B})
= (\nabla^sN_{u_A},g'(d)\sigma^+N_{d_B}),
\]
因为 \(\sigma=[g(d)+k]\sigma^+ + \sigma^-\Rightarrow \partial\sigma/\partial d = g'(d)\sigma^+\)。
3. \(K_{dd}\) 来自对
\((N_{d_A},\frac{g_c}{l}d + g'(d)\psi^+) + (\nabla N_{d_A},g_cl\nabla d)\)
再对 \(d_B\) 求导，分别得到质量样项、反应样项与扩散样项。

---

### 公式 2.2-7（相场不可逆盒约束：节点形式）

\[
d_A^{(n)}\leq d_A\leq 1,
\]

**各部分含义：**

- \(d_A^{(n)}\)：上一步（已收敛）损伤值，当前步已知。
- \(d_A\)：当前步待求损伤值。
- 上界 \(1\)：完全损伤极限。

**公式含义：**

把连续不可逆条件 \(d_n\le d_{n+1}\le1\) 离散到每个节点，形成典型 box constraints，正是后续梯度投影/L-BFGS-B 的约束形式。

**必要推导过程：**

从连续约束在节点插值点（或节点自由度）直接离散：
\[
0\le d_n(\mathbf x)\le d(\mathbf x)\le1
\quad\Longrightarrow\quad
0\le d_A^{(n)}\le d_A\le1.
\]
在本文实现中下界主要取 \(d_A^{(n)}\)。

---

### 公式 2.2-8（原文式 (11)：对角块矩阵）

\[
\hat{\mathbf{K}} = \left[ \begin{array}{cc}\mathbf{K}_{uu} & \mathbf{0}\\ \mathbf{0} & \mathbf{K}_{dd} \end{array} \right], \quad (11)
\]

**各部分含义：**

- \(\hat{\mathbf K}\)：去除耦合块后的块对角近似矩阵。
- 与完整 \(\mathbf K\) 对比，省略 \(K_{ud},K_{du}\)。

**公式含义：**

这是为第 3 节算法（尤其子空间/预条件相关步骤）准备的结构化近似矩阵，保留主块、降低处理成本。

**必要推导过程：**

不是从代数恒等式推得，而是算法上的“定义与选取”：
从完整 Hessian
\(\mathbf K=\begin{bmatrix}K_{uu}&K_{ud}\\K_{du}&K_{dd}\end{bmatrix}\)
提取主对角块形成
\(\hat K\)，用于后续数值流程。

---

## 补充：第 2 章节中出现的关键行内符号/关系（逐条解释）

> 下列内容虽多为行内表达，但同样属于第 2 章公式信息的一部分，故一并完整列出。

### 行内式 A：\(\mathbf{V} = H_0^1(\Omega)\), \(\mathbf{W}=H^1(\Omega)\)

\[
\mathbf{V} = H_0^1(\Omega), \qquad \mathbf{W}=H^1(\Omega).
\]

- 含义： 位移场函数空间含齐次 Dirichlet 条件；相场函数空间为一般 \(H^1\) 空间。
- 推导背景：由弱式对梯度可积性的要求确定（至少一阶弱导数平方可积）。

### 行内式 B：\(\pmb{\epsilon}=\nabla^{(s)}\pmb{u}\)

\[
\pmb{\epsilon}=\nabla^{(s)}\pmb{u}=\frac{1}{2}(\nabla\pmb u+\nabla\pmb u^T).
\]

- 含义：小变形应变定义。
- 推导背景：小应变假设下 Green-Lagrange 应变线性化得到。

### 行内式 C：\(\partial\psi/\partial d = g'(d)\psi^+\)

\[
\frac{\partial\psi}{\partial d} = g'(d)\psi^+(\epsilon), \qquad (k\text{ 对 }d\text{ 不依赖}).
\]

- 含义：相场方程中的“反应项”驱动力来自拉伸能。
- 推导： 由 \(\psi=([g(d)+k]\psi^+) + \psi^-\) 对 \(d\) 求导，\(\psi^\pm\) 在此不显含 \(d\)。

### 行内式 D：当 \(g(d)=(1-d)^2\) 时

\[
g'(d)=-2(1-d), \qquad g''(d)=2.
\]

- 含义：在式 (9)(10) 中实际可代入的显式导数。
- 推导： 直接对二次函数求导。

### 行内式 E：Einstein 求和约定

\[
\pmb u = \pmb N_{u_A}\pmb u_A \equiv \sum_A \pmb N_{u_A}\pmb u_A,
\qquad
d = N_{d_A}d_A \equiv \sum_A N_{d_A}d_A.
\]

- 含义：重复下标代表对节点索引求和，简化记号。
- 推导背景： 张量/有限元文献的标准缩并记法。

---

## 与第 2 章结构一致的总结

1. **2.1 小结**：先给出拉压分裂本构（式 (5) 到式 (7) 及相关定义），得到两条弱式残量。
2. **2.2 小结**：将弱式投影到有限元空间，得到离散目标函数（式 (8)）、梯度（式 (9)）和 Hessian（式 (10)），并给出节点盒约束与块对角矩阵（式 (11) 与节点约束式）。
3. **算法意义**：第 2 章给出的所有量正是第 3 章梯度投影 + L-BFGS-B 单调迭代所需的“目标-梯度-（近似）Hessian-约束”四件套。

---

## 完整性核对清单（第 2 章公式）

- [x] 式 (5)
- [x] 式 (6)
- [x] 正负括号与 Heaviside 定义
- [x] 应变谱分解与投影张量 \(\mathbf M_\alpha\)
- [x] \(\epsilon^+,\epsilon^-\) 定义
- [x] \(\psi^+,\psi^-\) 定义
- [x] \(\sigma\)、\(\sigma^+\)、\(\sigma^-\)
- [x] 切线模量 \(\partial\sigma/\partial\epsilon\)
- [x] \(\mathbb P^+,\mathbb P^-\) 定义
- [x] 式 (7)
- [x] 弱式残量方程组
- [x] 有限元插值 \(u,d\)
- [x] 变分插值 \(\delta u,\delta d\)
- [x] 式 (8)
- [x] 式 (9)
- [x] Hessian 分块结构
- [x] 式 (10)
- [x] 节点盒约束 \(d_A^{(n)}\le d_A\le1\)
- [x] 式 (11)
- [x] 第 2 章关键行内关系补充
