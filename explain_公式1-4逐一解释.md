# explain.md 公式(1)–(4)逐一解释（含必要推导）

> 说明：本文档严格结合 `explain.md` 中第 1 节与第 2 节开头的原始叙述，对公式 (1)–(4) 逐条解释。每一条均先完整展示公式，再给出符号含义、物理含义与必要推导步骤。为保证可渲染性，行间公式统一使用 `\[ ... \]`，行内符号统一使用 `\( ... \)`。

---

## 1. 公式 (1)：断裂体系总势能泛函

\[
\Pi (\mathbf{u},d) = \int_{\Omega}\psi \big(\mathbf{e}(\mathbf{u}),d\big)\,\mathrm{d}\Omega + g_{c}\,\Gamma_{l}(d) - \int_{\Omega}\mathbf{b}\cdot \mathbf{u}\,\mathrm{d}\Omega -\int_{\partial \Omega}\mathbf{t}\cdot \mathbf{u}\,\mathrm{d}\Gamma, \quad (1)
\]

### 1.1 各部分与符号含义

- \(\Pi(\mathbf{u},d)\)：系统总势能（待最小化的目标泛函）。
- \(\mathbf{u}(\mathbf{x})\)：位移场（向量未知量）。
- \(d(\mathbf{x})\)：相场/损伤场（标量未知量），通常 \(d=0\) 表示完好，\(d=1\) 表示完全损伤。
- \(\Omega\)：计算域（体域）。
- \(\partial\Omega\)：边界。
- \(\psi\big(\mathbf{e}(\mathbf{u}),d\big)\)：受损材料应变能密度。
- \(\mathbf{e}(\mathbf{u})=\nabla^{(s)}\mathbf{u}\)：小变形应变张量（对称梯度）。
- \(g_c\)：临界能量释放率（断裂韧度，单位通常为 \(\mathrm{N/m}\)）。
- \(\Gamma_l(d)\)：裂纹表面积（在相场框架下的正则化近似，见公式 (2)）。
- \(\mathbf{b}\)：体力密度（如重力项）。
- \(\mathbf{t}\)：Neumann 边界上的外力（面力/牵引力）。
- \(-\int_{\Omega}\mathbf{b}\cdot\mathbf{u}\,\mathrm{d}\Omega\)：体力外载势能（负号体现“外力作功降低系统势能”）。
- \(-\int_{\partial\Omega}\mathbf{t}\cdot\mathbf{u}\,\mathrm{d}\Gamma\)：边界外载势能。

### 1.2 公式整体物理含义

公式 (1) 表达的是 Griffith 型断裂思想在相场方法中的能量形式：

1. 第一项 \(\int_{\Omega}\psi\,\mathrm{d}\Omega\) 是材料体内的弹性能（受损后会退化）；
2. 第二项 \(g_c\Gamma_l(d)\) 是“生成新裂纹面”所需的断裂能；
3. 后两项是外载做功（体力与边界力）。

在准静态问题中，系统在每个载荷步趋向于选择使 \(\Pi\) 最小的 \((\mathbf{u},d)\)。因此，裂纹扩展不是“手动指定路径”，而是由能量竞争自动决定。

### 1.3 必要推导（从经典 Griffith 到相场势能）

#### 第一步：经典脆性断裂能量最小化思想

经典 Griffith 思想可写为

\[
\Pi_{\text{Griffith}} = \int_{\Omega\setminus\Gamma} \psi_0\big(\mathbf{e}(\mathbf{u})\big)\,\mathrm{d}\Omega + g_c\,|\Gamma| - W_{\text{ext}},
\]

其中 \(|\Gamma|\) 为真实裂纹面测度，\(W_{\text{ext}}\) 为外载做功。

#### 第二步：用相场变量 \(d\) 正则化裂纹几何

将“尖锐裂纹面” \(\Gamma\) 用扩散带表示，引入裂纹表面近似 \(\Gamma_l(d)\)（公式 (2)）；同时将完好材料弹性能 \(\psi_0\) 改写为受损形式 \(\psi(\mathbf{e}(\mathbf{u}),d)\)（常见做法是乘以退化函数 \(g(d)\)）。

于是得到

\[
\Pi(\mathbf{u},d)=\int_{\Omega}\psi\big(\mathbf{e}(\mathbf{u}),d\big)\,\mathrm{d}\Omega + g_c\Gamma_l(d) - W_{\text{ext}}.
\]

#### 第三步：展开外载做功

\[
W_{\text{ext}}=\int_{\Omega}\mathbf{b}\cdot\mathbf{u}\,\mathrm{d}\Omega + \int_{\partial\Omega}\mathbf{t}\cdot\mathbf{u}\,\mathrm{d}\Gamma.
\]

代回后即得到 explain.md 中公式 (1) 的完整形式。

---

## 2. 公式 (2)：裂纹表面积正则化（AT2 型）

\[
\Gamma_{l}(d) = \int_{\Omega}\gamma (d,\nabla d)\,\mathrm{d}\Omega = \int_{\Omega}\frac{1}{2l}\left(d^{2} + l^{2}\nabla d\cdot \nabla d\right)\,\mathrm{d}\Omega, \quad (2)
\]

### 2.1 各部分与符号含义

- \(\Gamma_l(d)\)：真实裂纹面积的相场近似。
- \(\gamma(d,\nabla d)\)：裂纹表面密度函数（单位体积上的“裂纹面贡献”）。
- \(l\)：相场长度尺度参数，控制损伤带宽度。
- \(\dfrac{1}{2l}d^2\)：局部损伤“势”项，鼓励 \(d\) 在无裂处接近 0。
- \(\dfrac{l}{2}|\nabla d|^2\)（等价于上式第二部分）：梯度正则项，抑制 \(d\) 的非物理振荡，控制裂纹扩散宽度。

### 2.2 公式整体物理含义

公式 (2) 的作用是把“几何上的裂纹面测度”转化为“定义在整个域上的体积分”。

- 当 \(l\) 越小，损伤带越窄，\(\Gamma_l(d)\) 越接近尖锐裂纹面积；
- 当 \(l\) 较大，裂纹被更宽地“抹平”，数值更稳定但几何更扩散。

这使得裂纹扩展问题可以在标准有限元网格上求解，而不必显式跟踪裂纹面。

### 2.3 必要推导（为何是该二次型）

#### 第一步：将表面能转写为体积分近似

相场理论中希望有

\[
\lim_{l\to 0}\Gamma_l(d)=|\Gamma|,
\]

即 \(\Gamma_l\) 在 \(l\to 0\) 时逼近真实裂纹面测度。常用构造是

\[
\Gamma_l(d)=\int_\Omega \left(\frac{w(d)}{l}+l\,|\nabla d|^2\right)\,\mathrm{d}\Omega,
\]

其中 \(w(d)\) 是损伤耗散势。

#### 第二步：AT2 选择

AT2 模型取

\[
w(d)=\frac{1}{2}d^2,
\]

于是

\[
\Gamma_l(d)=\int_\Omega\left(\frac{1}{2l}d^2+\frac{l}{2}|\nabla d|^2\right)\,\mathrm{d}\Omega
=\int_{\Omega}\frac{1}{2l}\left(d^2+l^2\nabla d\cdot\nabla d\right)\,\mathrm{d}\Omega,
\]

即 explain.md 的公式 (2)。

#### 第三步：欧拉方程层面的直观解释

对 \(d\) 做变分会出现

\[
\frac{1}{l}d-l\,\Delta d,
\]

其中 \(-l\Delta d\) 使损伤空间分布平滑，\(\frac{1}{l}d\) 则将损伤限制在局部，从而形成具有有限宽度 \(\mathcal{O}(l)\) 的裂纹过渡带。

---

## 3. 公式 (3)：当前载荷步的约束最小化问题

\[
\left(\mathbf{u}_{n + 1},d_{n + 1}\right) = \arg \min \Pi (\mathbf{u},d), \quad (3)
\]

### 3.1 各部分与符号含义

- \(n\)：当前时间（或伪时间/载荷步）编号。
- \((\mathbf{u}_n,d_n)\)：上一步已收敛解。
- \((\mathbf{u}_{n+1},d_{n+1})\)：本步要求解的未知场。
- \(\arg\min\)：给出使目标函数最小的自变量。
- 目标函数即公式 (1) 的 \(\Pi(\mathbf{u},d)\)。

### 3.2 公式整体物理含义

公式 (3) 表示：在每个载荷增量内，结构响应通过“最小势能原理”确定。由于 \(d\) 被纳入未知量，这里是位移场与损伤场的单体化（monolithic）联合求解。

### 3.3 必要推导（从变分原理到离散优化）

#### 第一步：连续问题的稳定平衡条件

准静态、等温条件下，稳定解满足在可行集合内的势能最小：

\[
(\mathbf{u},d)\in \mathcal{K},\quad \Pi(\mathbf{u},d)\le \Pi(\tilde{\mathbf{u}},\tilde d),\ \forall (\tilde{\mathbf{u}},\tilde d)\in\mathcal{K}.
\]

其中 \(\mathcal{K}\) 包含位移边界条件和损伤不可逆约束。

#### 第二步：时间离散到第 \(n+1\) 步

把已知历史 \((\mathbf{u}_n,d_n)\) 作为本步约束的下界来源，可得本步优化问题

\[
(\mathbf{u}_{n+1},d_{n+1})=\underset{(\mathbf{u},d)\in\mathcal{K}_{n+1}}{\operatorname{argmin}}\ \Pi(\mathbf{u},d),
\]

这与公式 (3) 完全一致（约束细节由公式 (4) 给出）。

#### 第三步：离散有限元后的向量形式

离散后可写成

\[
\mathbf{x}_{n+1}=\arg\min\ \Pi_h(\mathbf{x}),\quad \mathbf{x}:=\begin{bmatrix}\mathbf{u}_h\\ d_h\end{bmatrix},
\]

并配合分量级 box 约束（由公式 (4) 离散得到）。这正是后续 explain.md 第 3 节中梯度投影/活动集方法处理的对象。

---

## 4. 公式 (4)：相场不可逆与上下界约束

\[
0\leq d_{n}\leq d_{n + 1}\leq 1. \quad (4)
\]

### 4.1 各部分与符号含义

- \(d_n\)：上一步损伤值（已知）。
- \(d_{n+1}\)：当前步损伤值（未知）。
- 下界 \(0\)：无损伤状态。
- 上界 \(1\)：完全损伤状态。
- 中间不等式 \(d_n\le d_{n+1}\)：损伤不可逆（裂纹不愈合）。

### 4.2 公式整体物理含义

该式同时编码三层物理约束：

1. **状态可行性**：\(d\in[0,1]\)。
2. **热力学一致性**：损伤只能累积，不能回退。
3. **数值算法友好性**：在离散后天然形成分量独立的 box 约束，适配梯度投影/L-BFGS-B。

### 4.3 必要推导（从耗散不等式到 box 约束）

#### 第一步：不可逆条件

脆性断裂中常用内变量不可逆性写作

\[
\dot d\ge 0.
\]

#### 第二步：时间离散

在一步增量 \([t_n,t_{n+1}]\) 内做后向离散，可得

\[
\frac{d_{n+1}-d_n}{\Delta t}\ge 0 \ \Longrightarrow\ d_{n+1}\ge d_n.
\]

#### 第三步：物理边界

定义上要求

\[
0\le d\le 1.
\]

对第 \(n\) 与 \(n+1\) 两步写出即

\[
0\le d_n\le d_{n+1}\le 1,
\]

得到公式 (4)。

#### 第四步：有限元节点上的分量形式（算法直接使用）

对每个节点自由度 \(i\) 有

\[
\mathrm{lb}_i = d_i^{(n)},\qquad d_i^{(n+1)}\in[\mathrm{lb}_i,\mathrm{ub}_i],\qquad \mathrm{ub}_i=1,
\]

即

\[
 d_i^{(n)} \le d_i^{(n+1)} \le 1,
\]

这正是 explain.md 后续将约束写成 box constraints 的来源。

---

## 5. 四个公式之间的逻辑关系（按论文叙述顺序）

1. 公式 (1) 给出总势能目标泛函 \(\Pi\)。
2. 公式 (2) 具体定义了公式 (1) 中断裂能项 \(\Gamma_l(d)\) 的正则化表达。
3. 公式 (3) 规定在每个载荷步通过最小化 \(\Pi\) 来求 \((\mathbf{u}_{n+1},d_{n+1})\)。
4. 公式 (4) 给出该最小化问题必须满足的相场不等式约束。

因此，(1)+(2) 决定“优化目标”，(4) 决定“可行域”，(3) 给出“求解任务”的数学形式。

---

## 6. 便于对照 explain.md 的紧凑版符号表

- \(\mathbf{u}\)：位移向量场。
- \(d\)：相场/损伤标量场。
- \(\mathbf{e}(\mathbf{u})=\nabla^{(s)}\mathbf{u}\)：小应变。
- \(\psi\)：应变能密度（含损伤退化）。
- \(g_c\)：断裂韧度（单位裂纹面能量）。
- \(\Gamma_l(d)\)：裂纹面近似。
- \(\gamma(d,\nabla d)\)：裂纹密度函数。
- \(l\)：相场长度尺度。
- \(\mathbf{b},\mathbf{t}\)：体力与边界牵引。
- \(\Omega,\partial\Omega\)：域与边界。
- \(n,n+1\)：前后两个载荷步编号。

---

## 7. 与用户格式要求的一致性说明

本文档中的所有行间公式均采用如下可渲染形式：

\[
\text{示例：}\quad \mathbf{u}_{b} \approx \alpha(\chi)\,\big(\mathbf{u}_{\mathrm{pres}}-\mathbf{u}_{f}\big)+\mathbf{u}_{\mathrm{pres}}.
\]

即统一使用 `\[ ... \]` 形式展示公式；行内符号统一使用 `\( ... \)` 形式，以保证在 Markdown 数学渲染器中的一致显示效果。
