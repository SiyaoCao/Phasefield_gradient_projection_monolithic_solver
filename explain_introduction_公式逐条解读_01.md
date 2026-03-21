# explain.md 第1章（Introduction）公式逐一解读（含必要推导）

> 说明：本文件仅解读 `/home/runner/work/Phasefield_gradient_projection_monolithic_solver/Phasefield_gradient_projection_monolithic_solver/explain.md` 中 `## 1. Introduction` 章节出现的全部公式与符号表达。  
> 组织方式：严格按章节叙述顺序展开；每条先给出完整公式，再给出含义、推导和解读。  
> 公式渲染：统一使用 `[ ... ]`（即 `\[ ... \]`）形式；行内符号也单独列出并解释。

---

## 1) 总势能泛函（式(1)）

\[
\Pi (\pmb {u},d) = \int_{\Omega}\psi (\pmb {e}(\pmb {u}),d)\mathrm{d}\Omega + g_{c}\Gamma_{l}(d) - \int_{\Omega}\pmb {b}\cdot \pmb {u}\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot \pmb {u}\mathrm{d}\Gamma , \quad (1)
\]

### 含义
- 这是脆性断裂相场模型在准静态条件下的**总势能**（待最小化目标函数）。
- 前两项是系统“储能/耗能”贡献：
  - 弹性应变能（受相场损伤变量影响）；
  - 裂纹表面能（通过相场正则化近似）。
- 后两项是外力势能（体力、边界力）做功项，带负号。

### 逐项解释
1. \[ \int_{\Omega}\psi (\pmb {e}(\pmb {u}),d)\,\mathrm{d}\Omega \]
   - 体积分形式的应变能。
   - \[\psi\] 依赖应变 \[\pmb e\] 与相场 \[d\]，表示损伤导致材料刚度退化。

2. \[ g_c\,\Gamma_l(d) \]
   - 裂纹扩展耗能。
   - \[g_c\] 是临界能量释放率（单位面积断裂能）。
   - \[\Gamma_l(d)\] 是裂纹表面积（在相场框架中的弥散近似）。

3. \[ -\int_{\Omega}\pmb b\cdot\pmb u\,\mathrm d\Omega \]
   - 体力做功对应的势能项。

4. \[ -\int_{\partial\Omega}\pmb t\cdot\pmb u\,\mathrm d\Gamma \]
   - 边界牵引力做功对应的势能项。

### 必要推导（从“内能-外功”结构得到式(1)）
- 准静态体系总势能一般写为
  \[
  \Pi = \mathcal U - \mathcal W_{\text{ext}}
  \]
- 其中
  \[
  \mathcal U = \int_\Omega \psi\,\mathrm d\Omega + g_c\Gamma_l(d),
  \quad
  \mathcal W_{\text{ext}} = \int_\Omega \pmb b\cdot\pmb u\,\mathrm d\Omega + \int_{\partial\Omega}\pmb t\cdot\pmb u\,\mathrm d\Gamma
  \]
- 代入可得
  \[
  \Pi = \int_\Omega \psi\,\mathrm d\Omega + g_c\Gamma_l(d)-\int_\Omega \pmb b\cdot\pmb u\,\mathrm d\Omega-\int_{\partial\Omega}\pmb t\cdot\pmb u\,\mathrm d\Gamma
  \]
- 即式(1)。

### 与求解器的关系
- 数值算法本质上是在每个载荷步求解
  \[
  \min_{\pmb u,d}\Pi(\pmb u,d)
  \]
- 非凸性、约束可行性（不可逆）和计算规模，构成后续算法设计三大核心难点。

---

## 2) 引入的未知场与符号（式(1)后紧随定义）

### 2.1 位移场
\[
\pmb u(\pmb x)
\]
- 含义：位置 \[\pmb x\] 处的向量位移未知量。

### 2.2 相场（损伤）变量
\[
d(\pmb x)
\]
- 含义：标量损伤变量。
- 常见物理解释：
  - \[d=0\]：完好；
  - \[d=1\]：完全断裂；
  - \[0<d<1\]：弥散裂纹带。

### 2.3 体力、牵引力
\[
\pmb b,\quad \pmb t
\]
- \[\pmb b\]：体力密度；\[\pmb t\]：边界牵引力。

### 2.4 小变形应变张量
\[
\pmb e=\nabla^{(s)}\pmb u
\]
- 其中对称梯度定义为
  \[
  \nabla^{(s)}\pmb u=\frac{1}{2}\left(\nabla\pmb u+(\nabla\pmb u)^\mathsf T\right)
  \]
- 含义：在线弹性小变形框架下，应变由位移对称梯度给出。

### 2.5 其余关键参数
\[
\psi,\quad g_c,\quad \Gamma_l
\]
- \[\psi\]：应变能密度；\[g_c\]：断裂韧度参数；\[\Gamma_l\]：弥散裂纹面积泛函。

---

## 3) 裂纹表面积正则化近似（式(2)）

\[
\Gamma_{l}(d) = \int_{\Omega}\gamma (d,\nabla d)\mathrm{d}\Omega = \int_{\Omega}\frac{1}{2l}\left(d^{2} + l^{2}\nabla d\cdot \nabla d\right)\mathrm{d}\Omega , \quad (2)
\]

### 含义
- 该式是 Ambrosio–Tortorelli 类型近似：用平滑相场 \[d\] 近似尖锐裂纹面。
- 被积函数
  \[
  \gamma(d,\nabla d)=\frac{1}{2l}(d^2+l^2|\nabla d|^2)
  \]
  称为裂纹表面密度函数。

### 参数 \[l\] 的作用
\[
l
\]
- \[l\] 是长度尺度参数，控制裂纹弥散带宽度。
- \[l\] 小：裂纹更尖锐、网格要求更高；\[l\] 大：裂纹更平滑、计算更稳但几何精度下降。

### 必要推导（1D 截面理解）
考虑一维理想化情形（裂纹法向坐标 \[x\]）：
\[
\Gamma_l \sim \int \frac{1}{2l}\left(d^2+l^2(d')^2\right)\,\mathrm dx
\]
对被积函数使用配方法可见其最小化倾向于平衡两部分：
\[
\frac{1}{2l}d^2 \quad\text{与}\quad \frac{l}{2}(d')^2
\]
- 第一项惩罚“到处损伤”；
- 第二项惩罚“损伤梯度过陡”。

因此最优 \[d\] 会在有限宽度内平滑过渡，形成“弥散裂纹带”，其宽度量级由 \[l\] 控制。此即式(2)的核心物理-数值意义。

---

## 4) 准静态伪时间与边界载荷表达（Introduction 中行内公式）

### 4.1 伪时间
\[
t
\]
- 含义：不是动力学真实时间，而是载荷步参数。

### 4.2 位移型边界条件
\[
\pmb u=\hat{\pmb u}(t)
\]
- 含义：随载荷步变化的 Dirichlet 条件。

### 4.3 压力型边界条件
\[
p=\hat p(t)
\]
- 含义：随载荷步变化的 Neumann 条件。

### 4.4 当前时间步区间
\[
[t_n,t_{n+1}]
\]
- 含义：从第 \[n\] 步到 \[n+1\] 步的离散载荷增量区间。

---

## 5) 函数空间与离散步未知量（Introduction 中行内公式）

### 5.1 上一步已知解
\[
(\pmb u_n,d_n)\in \mathbf V\times\mathbf W
\]

### 5.2 函数空间定义
\[
\mathbf V=\mathbf H_0^1(\Omega),\quad \mathbf W=\mathbf H^1(\Omega)
\]

### 5.3 当前步待求解
\[
(\pmb u_{n+1},d_{n+1})\in \mathbf V\times\mathbf W
\]

### 含义解读
- \[\mathbf H_0^1(\Omega)\]：位移空间（满足齐次本质边界条件的 Sobolev 空间）。
- \[\mathbf H^1(\Omega)\]：相场空间（函数及一阶导可平方可积）。
- 这是典型的弱形式有限元设定，保证变分问题与离散求解兼容。

---

## 6) 当前载荷步的最小化问题（式(3)）

\[
(\pmb u_{n+1},d_{n+1})=\arg\min\Pi(\pmb u,d)\quad (3)
\]

### 含义
- 在当前步内，求一组 \[(\pmb u,d)\] 使总势能最小。
- 该式将断裂演化转化为“增量变分最小化”。

### 必要推导（从变分到离散优化）
1. 连续层面：
   \[
   \delta\Pi(\pmb u,d;\delta\pmb u,\delta d)=0
   \]
   给出 Euler–Lagrange 方程（平衡方程 + 相场方程）。
2. 离散层面（有限元）：
   \[
   \Pi(\mathbf z),\quad \mathbf z=[\mathbf u;\mathbf d]
   \]
   转化为有限维优化：
   \[
   \min_{\mathbf z}\Pi(\mathbf z)
   \]
3. 加入不可逆约束后，就是式(3)+(4)的约束优化问题。

### 实际数值意义
- 由于 \[\Pi\] 非凸，直接 Newton 容易失稳；这也是文中选择 BFGS / L-BFGS-B 路线的重要动机。

---

## 7) 相场不可逆与盒约束（式(4)）

\[
0\le d_n\le d_{n+1}\le 1.\quad (4)
\]

### 含义
- 约束 1：\[d_{n+1}\ge d_n\]，损伤不可恢复（裂纹不自愈）。
- 约束 2：\[0\le d_{n+1}\le 1\]，相场物理有界。

### 必要推导（由物理约束到数学不等式）
1. 物理上定义 \[d\in[0,1]\]；
2. 不可逆性要求 \[\dot d\ge 0\]；
3. 载荷步离散（后向差分观点）可写成
   \[
   \frac{d_{n+1}-d_n}{\Delta t}\ge 0
   \Rightarrow d_{n+1}\ge d_n
   \]
4. 合并上下界得到
   \[
   0\le d_n\le d_{n+1}\le 1
   \]
即式(4)。

### 额外相关初值表达（Introduction 行内）
\[
t=0
\]
\[
d_0(\pmb x)=0
\]
- 含义：初始时刻默认全域无损伤。

---

## 8) 图1对应的梯度投影核心表达（Introduction 行内公式）

文中先引入目标函数与可行域：
\[
f(x)
\]
\[
C
\]
\[
\nabla f(x)
\]

并给出投影操作写法（下式为论文原文记法，属于示意性写法）：
\[
\mathrm{Proj}_C(\mathbf{x}_k = \mathbf{a}_k\nabla f(\mathbf{x}_k)).
\]

为避免歧义，也可将其理解为“对某个与梯度方向相关的中间量 \\(\\mathbf y_k\\) 做投影”，即
\[
\mathbf x_{k+1}=\mathrm{Proj}_C(\mathbf y_k),\quad \mathbf y_k=\mathbf x_k-\alpha_k\nabla f(\mathbf x_k).
\]

### 含义与规范写法说明
- 该行表达的是“把梯度相关下降步投影到可行域 \[C\]”。
- 优化文献中更常见写法是
  \[
  \mathbf x_{k+1}=\mathrm{Proj}_C\left(\mathbf x_k-\alpha_k\nabla f(\mathbf x_k)\right)
  \]
  其中 \[\alpha_k>0\] 为步长。
- 本文 Introduction 的写法可理解为“投影算子作用于梯度步相关表达式”的示意写法。

### 投影的数学定义
\[
\mathrm{Proj}_C(\mathbf y)=\arg\min_{\mathbf z\in C}\|\mathbf z-\mathbf y\|_2
\]
- 当 \[C\] 为一般凸集，投影可能较贵；
- 当 \[C\] 为盒约束，投影可分量截断，代价极低。

---

## 9) 盒约束分量形式（Introduction 独立公式）

\[
\mathrm{lb}_i\le x_i\le \mathrm{ub}_i,
\]

### 含义
- 对未知向量第 \[i\] 个分量 \[x_i\] 给上下界。
- \[\mathrm{lb}_i\]、\[\mathrm{ub}_i\] 分别是下界与上界。

### 对应投影闭式（必要推导）
对任意分量 \[y_i\]，投影到区间 \[[\mathrm{lb}_i,\mathrm{ub}_i]\] 为
\[
(\mathrm{Proj}_C(\mathbf y))_i=\min\left(\max(y_i,\mathrm{lb}_i),\mathrm{ub}_i\right)
\]
推导要点：一维区间上的欧氏投影就是“区间截断（clipping）”。

这就是 L-BFGS-B 高效处理 box constraints 的根本原因之一。

---

## 10) 相场离散后盒约束实例化（Introduction 独立公式）

\[
\mathrm{lb}_i=d_i^{(n)}\le d_i^{(n+1)}\le 1=\mathrm{ub}_i,
\]

### 含义
- 在当前载荷步，上一时刻节点相场 \[d_i^{(n)}\] 已知，因此它自然成为当前步下界；
- 上界取 \[1\]，符合相场定义域。

### 与式(4)的离散对应关系
- 将式(4)对每个自由度 \[i\] 写成分量形式即得上式。
- 这一步把“连续约束”落到“代数未知向量约束”，从而可直接用 L-BFGS-B 类算法处理。

---

## 11) Introduction 中出现的关键符号/不等式补充逐条解释

以下虽多为行内表达，但均在该章承担明确数学含义：

### 11.1 未知与已知变量分层
\[
(\pmb u_n,d_n)\ \text{已知},\qquad (\pmb u_{n+1},d_{n+1})\ \text{未知}
\]
- 意义：增量求解框架，历史步结果作为当前步约束和初值信息来源。

### 11.2 物理边界值
\[
d=0\quad (\text{undamaged}),\qquad d=1\quad (\text{fully damaged})
\]
- 意义：相场的两端物理态定义，决定了约束上、下界解释。

### 11.3 不可逆核心不等式
\[
d_{n+1}\ge d_n
\]
- 意义：裂纹演化单调，不允许“愈合回退”。

### 11.4 空间域与边界域积分符号
\[
\int_\Omega(\cdot)\,\mathrm d\Omega,\qquad \int_{\partial\Omega}(\cdot)\,\mathrm d\Gamma
\]
- 体积分与边界积分分别对应体力功与边界功，亦对应有限元中单元积分与边界面积分。

---

## 12) 将 Introduction 全部公式串联成统一优化问题（总结性重写）

根据第1章中所有公式，可将每一载荷步的数学问题写为：

\[
\min_{\pmb u,d}\ \Pi(\pmb u,d)
\]

其中

\[
\Pi (\pmb {u},d) = \int_{\Omega}\psi (\nabla^{(s)}\pmb {u},d)\,\mathrm{d}\Omega
+ g_{c}\int_{\Omega}\frac{1}{2l}\left(d^{2} + l^{2}\nabla d\cdot \nabla d\right)\,\mathrm{d}\Omega
- \int_{\Omega}\pmb {b}\cdot \pmb {u}\,\mathrm{d}\Omega
- \int_{\partial \Omega}\pmb {t}\cdot \pmb {u}\,\mathrm{d}\Gamma
\]

并满足

\[
\pmb u\in \mathbf H_0^1(\Omega),\qquad d\in \mathbf H^1(\Omega)
\]
\[
0\le d_n\le d_{n+1}\le 1
\]

其离散到节点分量后等价为盒约束：

\[
\mathrm{lb}_i=d_i^{(n)}\le d_i^{(n+1)}\le 1=\mathrm{ub}_i
\]

因此可采用“投影 + 准牛顿”范式：

\[
\mathbf x_{k+1}=\mathrm{Proj}_C\left(\mathbf x_k-\alpha_k\nabla f(\mathbf x_k)\right)
\]

这正是 Introduction 里提出后续 L-BFGS-B 思路的数学落脚点。

---

## 13) 对第1章公式的工程实现含义（面向本仓库求解器）

1. \[\Pi\] 决定了残差与切线（或近似 Hessian）构造来源。  
2. \[\Gamma_l\] 的 \[|\nabla d|^2\] 项决定相场方程具有扩散型特征。  
3. \[0\le d_n\le d_{n+1}\le 1\] 是实现中必须逐自由度严格维护的硬约束。  
4. 盒约束形式 \[\mathrm{lb}_i\le x_i\le\mathrm{ub}_i\] 使投影步骤可做廉价截断。  
5. 这也是本文将“不可逆约束处理”与“非凸稳健迭代”统一到 L-BFGS-B 框架的关键逻辑。

---

## 14) 逐条核对清单（确保第1章公式无遗漏）

已覆盖公式/符号如下：

1. \[\Pi(\pmb u,d)\]（式(1)）  
2. \[\Gamma_l(d)=\int\gamma=\int\frac{1}{2l}(d^2+l^2\nabla d\cdot\nabla d)\]（式(2)）  
3. \[\pmb u=\hat{\pmb u}(t)\]  
4. \[p=\hat p(t)\]  
5. \[[t_n,t_{n+1}]\]  
6. \[(\pmb u_n,d_n)\in\mathbf V\times\mathbf W\]  
7. \[\mathbf V=\mathbf H_0^1(\Omega),\ \mathbf W=\mathbf H^1(\Omega)\]  
8. \[(\pmb u_{n+1},d_{n+1})\in\mathbf V\times\mathbf W\]  
9. \[(\pmb u_{n+1},d_{n+1})=\arg\min\Pi(\pmb u,d)\]（式(3)）  
10. \[0\le d_n\le d_{n+1}\le 1\]（式(4)）  
11. \[t=0\]  
12. \[d_0(\pmb x)=0\]  
13. \[f(x),\ C,\ \nabla f(x)\]  
14. \[\mathrm{Proj}_C(\mathbf{x}_k = \mathbf{a}_k\nabla f(\mathbf{x}_k))\]（论文原文中的非标准示意写法；规范实现见第8节给出的标准投影迭代表达）  
15. \[\mathrm{lb}_i\le x_i\le\mathrm{ub}_i\]  
16. \[\mathrm{lb}_i=d_i^{(n)}\le d_i^{(n+1)}\le 1=\mathrm{ub}_i\]

核对结论：`1. Introduction` 中出现的公式与符号表达均已逐条解释并给出必要推导。
