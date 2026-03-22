# explain.md 第三章节全部公式逐一解释（完整版）

> 范围：仅对应 `explain.md` 中 `## 3. Gradient projection based monolithic scheme`（含 3.1、3.2、3.3、3.4 与 Comment 1-4）里出现的全部显示公式。  
> 说明：每个条目都按“公式完整显示 → 推导/来源 → 各部分含义 → 公式整体含义”组织；公式统一使用 `\[ ... \]` 形式渲染。

---

## 3.1 Algorithm overview

### 公式 3.1-1（Eq.12）
\[
\left(\pmb{u}_A,d_A\right) = \arg \min \Pi (\pmb {u}_A,d_A)
\]

- 推导/来源：这是离散后相场断裂问题的基本优化模型，直接由“总势能最小原理”写成离散节点变量上的最小化问题。  
- 各部分含义：\(\pmb{u}_A\) 是节点位移自由度，\(d_A\) 是节点相场自由度，\(\Pi\) 是离散总势能泛函，\(\arg\min\) 表示取使能量最小的变量。  
- 公式整体含义：在每个载荷步内，通过联合求解位移和损伤相场，使系统达到当前约束下的稳定（或局部稳定）能量态。

### 公式 3.1-2（Eq.13）
\[
d_A^{(n)}\leq d_A\leq 1
\]

- 推导/来源：来自相场不可逆条件与相场取值区间定义。不可逆要求当前步损伤不小于前一步，物理上相场上界为 1。  
- 各部分含义：\(d_A^{(n)}\) 为上一步节点相场，\(d_A\) 为当前步相场，1 表示完全损伤上限。  
- 公式整体含义：该盒约束同时编码“不可自愈”和“损伤变量有界”。

### 公式 3.1-3（Eq.14）
\[
\left[\pmb {u}_A,d_A\right]^{\mathrm{T}} = \mathbf{C}\left[\pmb {u}_A,d_A\right]^{\mathrm{T}} + \pmb {k}.
\]

- 推导/来源：将 Dirichlet 边界与悬挂节点约束统一写成线性仿射约束。  
- 各部分含义：\(\mathbf{C}\) 为约束传递矩阵，\(\pmb{k}\) 为非齐次项。  
- 公式整体含义：除盒约束外，还必须满足有限元线性约束（边界条件、网格约束）。

### 公式 3.1-4（Eq.15）
\[
m_{k}(\pmb{x}) = \Pi_{k} + \pmb{r}_{k}^{\mathrm{T}}(\pmb{x} - \pmb{x}_{k}) + \frac{1}{2} (\pmb{x} - \pmb{x}_{k})^{\mathrm{T}}\mathbf{B}_{k}(\pmb{x} - \pmb{x}_{k})
\]

- 推导/来源：对目标函数在迭代点 \(\pmb{x}_k\) 做二次模型近似（零阶+一阶+二阶近似）。  
- 各部分含义：\(\Pi_k\) 为当前能量值，\(\pmb{r}_k=\nabla\Pi_k\) 为梯度，\(\mathbf{B}_k\) 为 BFGS 近似 Hessian。  
- 公式整体含义：把难解的非线性能量最小化转成可控的“局部二次子问题”。

### 公式 3.1-5（Eq.16）
\[
\mathrm{Proj}_c(x_i,\mathrm{lb}_i,\mathrm{ub}_i) = \left\{ \begin{array}{ll}\mathrm{lb}_i & \mathrm{if}\; x_i < \mathrm{lb}_i,\\ x_i & \mathrm{if}\; \mathrm{lb}_i \leq x_i \leq \mathrm{ub}_i,\\ \mathrm{ub}_i & \mathrm{if}\; x_i > \mathrm{ub}_i, \end{array} \right.
\]

- 推导/来源：凸盒约束集合上的欧几里得投影算子，按分量独立截断得到。  
- 各部分含义：\(\mathrm{lb}_i,\mathrm{ub}_i\) 分别是第 \(i\) 个分量下上界。  
- 公式整体含义：任何越界分量都被强制拉回可行区间。

### 公式 3.1-6（Eq.17）
\[
\pmb {x}(t) = \mathrm{Proj}_c(\pmb {x}_k - t\pmb {r}_k,\mathbf{lb},\mathbf{ub}),\quad t > 0
\]

- 推导/来源：先沿负梯度走 \(\pmb{x}_k-t\pmb r_k\)，再投影到可行域。  
- 各部分含义：\(t\) 是路径参数，\(\mathbf{lb},\mathbf{ub}\) 为盒约束向量。  
- 公式整体含义：给出“投影梯度路径”，用于后续 Cauchy 点搜索。

### 公式 3.1-7（Eq.18）
\[
p_k(t) = m_k(\pmb{x}(t)) = \Pi_k + \pmb{r}_k^\mathrm{T}(\pmb{x}(t) - \pmb{x}_k) + \frac{1}{2} (\pmb{x}(t) - \pmb{x}_k)^\mathrm{T}\mathbf{B}_k(\pmb{x}(t) - \pmb{x}_k)
\]

- 推导/来源：把 Eq.15 的二次模型沿 Eq.17 的路径进行复合。  
- 各部分含义：\(p_k(t)\) 是单变量分段二次函数。  
- 公式整体含义：将多维约束优化沿路径降维为一维搜索。

### 公式 3.1-8（活动集定义）
\[
\mathcal{A}(\pmb {x}) = \{i:x_i = \mathrm{lb}_i\} \cup \{i:x_i = \mathrm{ub}_i\}
\]

- 推导/来源：活动约束定义。  
- 各部分含义：集合中索引是“当前正好贴边界”的分量。  
- 公式整体含义：活动集决定哪些变量在子空间最小化中固定。

### 公式 3.1-9（子空间问题目标）
\[
\pmb{x}^* = \arg \min m_k(\pmb {x})
\]

- 推导/来源：在 Cauchy 点确定活动集后，继续最小化二次模型。  
- 各部分含义：\(\pmb{x}^*\) 为子空间最小化解。  
- 公式整体含义：从“首次可行下降点”继续走向更优局部解。

### 公式 3.1-10（子空间约束）
\[
x_{i} = x_{i}^{c},\;\forall i\in \mathcal{A}_{k}(\pmb{x}^{c})\quad \text{and}\quad \mathrm{lb}_{i}\leq x_{i}\leq \mathrm{ub}_{i},\;\forall i\notin \mathcal{A}_{k}(\pmb{x}^{c})
\]

- 推导/来源：活动变量固定，自由变量保留盒约束。  
- 各部分含义：\(x_i^c\) 为广义 Cauchy 点上的分量值。  
- 公式整体含义：把原约束分解成“固定变量 + 自由变量有界”。

### 公式 3.1-11（搜索方向）
\[
\pmb{p}_k = \pmb{x}^* - \pmb{x}_k
\]

- 推导/来源：由当前点到子空间解的位移定义方向。  
- 公式整体含义：用于线搜索的候选下降方向。

### 公式 3.1-12（迭代更新）
\[
\pmb{x}_{k + 1} = \pmb{x}_k + \alpha_k\pmb{p}_k
\]

- 推导/来源：标准线搜索更新。  
- 各部分含义：\(\alpha_k\) 为步长。  
- 公式整体含义：在可行下降方向上推进到下一迭代点。

### 公式 3.1-13（强 Wolfe：充分下降）
\[
\Pi_{k + 1} = \Pi (\pmb{x}_k + \alpha_k\pmb {p}_k) \leq \Pi_k + c_1\alpha_k\pmb {p}_k^\mathrm{T}\pmb{r}_k
\]

- 推导/来源：线搜索充分下降条件。  
- 各部分含义：\(c_1\in(0,1)\) 为小常数。  
- 公式整体含义：新点能量必须有足够下降。

### 公式 3.1-14（强 Wolfe：曲率）
\[
\left|\pmb{r}_{k + 1}^\mathrm{T}\pmb {p}_k\right| = \left|\pmb {r}(\pmb{x}_k + \alpha_k\pmb {p}_k)^\mathrm{T}\pmb {p}_k\right| \leq c_2\left|\pmb{r}_k^\mathrm{T}\pmb {p}_k\right|
\]

- 推导/来源：线搜索曲率条件。  
- 各部分含义：\(c_2\in(c_1,1)\)。  
- 公式整体含义：避免步长过小并保障拟牛顿更新曲率信息质量。

### 公式 3.1-15（收敛准则 1）
\[
\mathcal{A}_{k + 1}(\pmb{x}^c) = \mathcal{A}_k(\pmb{x}^c)
\]

- 推导/来源：活动集稳定性判据。  
- 公式整体含义：边界状态不再变化，表明约束结构趋于稳定。

### 公式 3.1-16（收敛准则 2）
\[
\| \mathrm{Proj}_C(\mathbf{x}_{k + 1} - \mathbf{r}_{k + 1},\mathbf{lb},\mathbf{ub}) - \mathbf{x}_{k + 1}\|_{2} < \mathrm{tol}
\]

- 推导/来源：投影梯度范数判据。  
- 公式整体含义：在约束意义下的一阶最优性残差足够小。

### 公式 3.1-17（收敛准则 3）
\[
\| \Delta \mathbf{x}_{k + 1}\|_{2} = \| \mathbf{x}_{k + 1} - \mathbf{x}_{k}\|_{2} < \mathrm{tol}
\]

- 推导/来源：增量判据。  
- 公式整体含义：解更新幅度足够小，迭代趋于停滞/收敛。

### 公式 3.1-18（Eq.19）
\[
\mathrm{Proj}_C(x_i - r_i,\mathrm{lb}_i,\mathrm{ub}_i) - x_i = \left\{ \begin{array}{ll}\mathrm{lb}_i - x_i & \mathrm{if}\; x_i - r_i < \mathrm{lb}_i,\\ - r_i & \mathrm{if}\; \mathrm{lb}_i \leq x_i - r_i \leq \mathrm{ub}_i,\\ \mathrm{ub}_i - x_i & \mathrm{if}\; x_i - r_i > \mathrm{ub}_i. \end{array} \right.
\]

- 推导/来源：把 Eq.16 代入 \(x_i-r_i\) 并减 \(x_i\) 逐段展开。  
- 各部分含义：中间段等于负梯度，越界段等于到边界距离。  
- 公式整体含义：统一刻画“最优性误差 + 约束违背距离”。

### 公式 3.1-19（Comment 1）
\[
\mathbf{x} = \mathbf{C}\mathbf{x} + \mathbf{k}
\]

- 推导/来源：Eq.14 的紧凑写法。  
- 公式整体含义：总自由度满足线性仿射约束。

### 公式 3.1-20（Comment 1 约束改写系统）
\[
\left(\mathbf{C}^{\mathrm{T}}\mathbf{A}\mathbf{C} + \mathbf{I}_{d_{c}}\right)\hat{\mathbf{x}} = \mathbf{C}^{\mathrm{T}}(\mathbf{b} - \mathbf{A}\mathbf{k})
\]

- 推导/来源：将 \(\mathbf{x}=\mathbf{C}\hat{\mathbf{x}}+\mathbf{k}\) 代入线性系统 \(\mathbf{A}\mathbf{x}=\mathbf{b}\)，并对受约束自由度加单位对角稳定项。  
- 公式整体含义：把原受约束系统转成可直接求解的改写系统。

### 公式 3.1-21（Comment 1 解恢复）
\[
\mathbf{x} = \mathbf{C}\hat{\mathbf{x}} + \mathbf{k}
\]

- 推导/来源：由变量替换反代得到。  
- 公式整体含义：先求辅助变量，再恢复真实自由度解。

### 公式 3.1-22（Comment 1 指示矩阵）
\[
(\mathbf{I}_{d_c})_{ii} = \begin{cases} 1 & \text{if DoF } i \text{ is constrained}, \\ 0 & \text{otherwise}, \end{cases}
\]

- 推导/来源：按受约束 DoF 构造对角指示矩阵。  
- 公式整体含义：只在受约束分量上施加单位稳定项。

---

## 3.2 Compact representation of limited-memory BFGS matrix

### 公式 3.2-1（位移量与梯度差）
\[
\mathbf{s}_k = \mathbf{x}_{k + 1} - \mathbf{x}_k,\quad \mathbf{y}_k = \mathbf{r}_{k + 1} - \mathbf{r}_k
\]

- 推导/来源：拟牛顿更新的标准差分向量定义。  
- 公式整体含义：\(\mathbf{s}_k,\mathbf{y}_k\) 携带局部曲率信息。

### 公式 3.2-2（BFGS 更新）
\[
\mathbf{B}_{k + 1} = \mathbf{B}_k - \frac{\mathbf{B}_k\mathbf{s}_k\mathbf{s}_k^{\mathrm{T}}\mathbf{B}_k}{\mathbf{s}_k^{\mathrm{T}}\mathbf{B}_k\mathbf{s}_k} + \frac{\mathbf{y}_k\mathbf{y}_k^{\mathrm{T}}}{\mathbf{y}_k^{\mathrm{T}}\mathbf{s}_k}
\]

- 推导/来源：满足割线条件与对称正定保持的经典 BFGS 公式。  
- 公式整体含义：用两次秩一修正更新 Hessian 近似。

### 公式 3.2-3（曲率条件）
\[
\mathbf{s}_k^{\mathrm{T}}\mathbf{y}_k > 0
\]

- 推导/来源：来自强 Wolfe 线搜索保证。  
- 公式整体含义：确保更新后 \(\mathbf{B}_{k+1}\) 保持正定。

### 公式 3.2-4（S 校正矩阵）
\[
\mathbf{S}_k = [\mathbf{s}_{k-m} \cdots \mathbf{s}_{k-1}]
\]

- 推导/来源：把最近 \(m\) 个 \(\mathbf{s}\) 向量拼成列块矩阵。  
- 公式整体含义：有限记忆窗口内位移历史。

### 公式 3.2-5（Y 校正矩阵）
\[
\mathbf{Y}_k = [\mathbf{y}_{k-m} \cdots \mathbf{y}_{k-1}]
\]

- 推导/来源：把最近 \(m\) 个 \(\mathbf{y}\) 向量拼成列块矩阵。  
- 公式整体含义：有限记忆窗口内梯度变化历史。

### 公式 3.2-6（Eq.20，紧凑 L-BFGS）
\[
\mathbf{B}_k = \mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k
\]

- 推导/来源：Byrd 等给出的紧凑表示，把“稠密更新”写成低秩修正。  
- 公式整体含义：不显式存稠密 \(\mathbf{B}_k\)，只存小矩阵和窄矩阵。

### 公式 3.2-7（W 定义）
\[
\mathbf{W}_k = [\mathbf{Y}_k \quad \mathbf{B}^0_k \mathbf{S}_k] \in \mathbb{R}^{n\times 2m}
\]

- 推导/来源：紧凑表示中的基矩阵定义。  
- 公式整体含义：把历史梯度差和预乘块统一到一个列空间。

### 公式 3.2-8（M 定义）
\[
\mathbf{M}_k = \begin{bmatrix} -\mathbf{D}_k & \mathbf{L}^{\mathrm{T}}_k \\ \mathbf{L}_k & \mathbf{S}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{S}_k \end{bmatrix}^{-1} \in \mathbb{R}^{2m\times 2m}
\]

- 推导/来源：由 BFGS 多步更新等价变换得到的小规模核矩阵。  
- 公式整体含义：所有历史曲率耦合信息浓缩在 \(2m\times2m\) 系统里。

### 公式 3.2-9（D 定义）
\[
\mathbf{D}_k = \mathrm{diag}\{\mathbf{s}^{\mathrm{T}}_{k-m}\mathbf{y}_{k-m}, \ldots, \mathbf{s}^{\mathrm{T}}_{k-1}\mathbf{y}_{k-1}\} \in \mathbb{R}^{m\times m}
\]

- 推导/来源：取各历史对的曲率内积做对角阵。  
- 公式整体含义：记录每对向量的局部曲率强度。

### 公式 3.2-10（L 定义）
\[
\mathbf{L}_k \in \mathbb{R}^{m\times m}, \quad (\mathbf{L}_k)_{i,j} = \begin{cases} \mathbf{s}^{\mathrm{T}}_{i+k-m-1}\mathbf{y}_{j+k-m-1} & \text{if } i > j, \\ 0 & \text{if } i \leq j. \end{cases}
\]

- 推导/来源：按索引构造严格下三角历史耦合矩阵。  
- 公式整体含义：保存不同历史步之间的非对角曲率耦合。

### 公式 3.2-11（Eq.21，初始块对角矩阵）
\[
\mathbf{B}^0_k = \hat{\mathbf{K}}^{(k)} = \begin{bmatrix} \mathbf{K}_{uu} & \mathbf{0} \\ \mathbf{0} & \mathbf{K}_{dd} \end{bmatrix}^{(k)}
\]

- 推导/来源：用单场切线块对角矩阵作为初始 Hessian 近似。  
- 各部分含义：\(\mathbf{K}_{uu}\) 位移子块，\(\mathbf{K}_{dd}\) 相场子块。  
- 公式整体含义：把有限元可装配、稀疏且正定的结构注入 L-BFGS 初值。

---

## 3.3 Generalized Cauchy point

### 公式 3.3-1（断点定义）
\[
t_i = \begin{cases} (x^0_i - \mathrm{ub}_i)/r_i & \text{if } r_i < 0, \\ (x^0_i - \mathrm{lb}_i)/r_i & \text{if } r_i > 0, \\ +\infty & \text{if } r_i = 0. \end{cases}
\]

- 推导/来源：沿 \(x_i^0-tr_i\) 到达边界时解 \(t\) 得到。  
- 公式整体含义：\(t_i\) 表示第 \(i\) 分量首次触及约束的“时间”。

### 公式 3.3-2（分量路径）
\[
x_i(t) = x_i^{0} - \min \{t,t_i\} r_i
\]

- 推导/来源：未触界时按梯度走，触界后冻结，故有 \(\min\)。  
- 公式整体含义：投影路径按分量的统一表达。

### 公式 3.3-3（区间线性段）
\[
\mathbf{x}(t) = \mathbf{x}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)}
\]

- 推导/来源：在相邻断点之间活动集不变，路径是线段。  
- 公式整体含义：把分段投影路径局部线性化。

### 公式 3.3-4（Eq.22，区间二次多项式展开）
\[
\begin{array}{rl}
p(t) &= m(\mathbf{x}(t)) = \Pi + \mathbf{r}^{\mathrm{T}}(\mathbf{x}(t) - \mathbf{x}^0) + \frac{1}{2} (\mathbf{x}(t) - \mathbf{x}^0)^{\mathrm{T}}\mathbf{B}(\mathbf{x}(t) - \mathbf{x}^0) \\
&= \Pi + \mathbf{r}^{\mathrm{T}}(\mathbf{z}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)}) + \frac{1}{2} (\mathbf{z}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)})^{\mathrm{T}}\mathbf{B}(\mathbf{z}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)}) \\
&= \left(\Pi + \mathbf{r}^{\mathrm{T}}\mathbf{z}^{(j-1)} + \frac{1}{2} \mathbf{z}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j-1)}\right) \\
&\quad + \left(\mathbf{r}^{\mathrm{T}}\mathbf{d}^{(j-1)} + \mathbf{d}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j-1)}\right)\Delta t + \frac{1}{2}\left(\mathbf{d}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j-1)}\right)\Delta t^{2} \\
&= f_{j-1} + f_{j-1}^{\prime}\Delta t + \frac{1}{2} f_{j-1}^{\prime \prime}\Delta t^{2} = \hat{p} (\Delta t)
\end{array}
\]

- 推导/来源：将 \(\mathbf{x}(t)=\mathbf{x}^{(j-1)}+\Delta t\mathbf{d}^{(j-1)}\) 代入二次模型并按 \(\Delta t\) 收集常数、一次、二次项。  
- 各部分含义：\(f_{j-1},f'_{j-1},f''_{j-1}\) 是该区间一维二次函数系数。  
- 公式整体含义：Cauchy 点搜索被转化为“每段上最小化一元二次函数”。

### 公式 3.3-5（二次项正性）
\[
f_{j-1}^{\prime \prime} = \mathbf{d}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j-1)} > 0
\]

- 推导/来源：\(\mathbf{B}\) 正定，任意非零方向二次型为正。  
- 公式整体含义：每个区间上 \(\hat p\) 是开口向上的凸抛物线。

### 公式 3.3-6（临界点）
\[
\hat{p}'(\Delta t) = 0 \Rightarrow \Delta t^* = -f_{j-1}' / f_{j-1}^{\prime \prime}
\]

- 推导/来源：一元二次函数导数置零。  
- 公式整体含义：给出当前线段上的候选最小点。

### 公式 3.3-7（活动集判定）
\[
\mathcal{A}(\mathbf{x}^c) = \{i : t_i \leq t^{(j-1)}, i = 1, \ldots, n\}
\]

- 推导/来源：断点不大于当前段起点的分量已触边界。  
- 公式整体含义：快速由断点索引恢复活动集。

### 公式 3.3-8（区间推进）
\[
\Delta t^{(j-1)} = t^{(j)} - t^{(j-1)}, \quad \mathbf{x}^{(j)} = \mathbf{x}^{(j-1)} + \Delta t^{(j-1)}\mathbf{d}^{(j-1)}
\]

- 推导/来源：线段长度与端点更新定义。  
- 公式整体含义：从上一段推进到下一段的状态传递。

### 公式 3.3-9（方向更新）
\[
\mathbf{d}^{(j)} = \mathbf{d}^{(j-1)} + r_b \mathbf{e}_b
\]

- 推导/来源：第 \(b\) 分量在新断点被固定后，投影方向分量变化。  
- 公式整体含义：活动集变化只引起单分量方向修正。

### 公式 3.3-10（Eq.23）
\[
f'_j = \mathbf{r}^{\mathrm{T}}\mathbf{d}^{(j)} + \mathbf{d}^{(j)}{}^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j)} = f'_{j-1} + \Delta t^{(j-1)} f''_{j-1} + r^2_b + r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j)}
\]

- 推导/来源：将 \(\mathbf{d}^{(j)}\)、\(\mathbf{z}^{(j)}\) 的递推关系代入 \(f'_j\) 定义并展开。  
- 公式整体含义：不用从零计算即可增量更新一次项系数。

### 公式 3.3-11（Eq.24）
\[
f''_j = \mathbf{d}^{(j)}{}^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j)} = f''_{j-1} + 2r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j-1)} + r^2_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{e}_b
\]

- 推导/来源：把 \(\mathbf{d}^{(j)}=\mathbf{d}^{(j-1)}+r_b\mathbf e_b\) 代入二次型并展开。  
- 公式整体含义：二次项系数也可低成本递推更新。

### 公式 3.3-12（紧凑 B 重申）
\[
\mathbf{B} = \mathbf{B}^0 - \mathbf{W}\mathbf{M}\mathbf{W}^{\mathrm{T}}
\]

- 推导/来源：Eq.20 在本节记号下重写。  
- 公式整体含义：后续把 \(\mathbf{e}_b^T\mathbf{B}(\cdot)\) 运算拆成稀疏部分与小矩阵部分。

### 公式 3.3-13（辅助向量递推）
\[
\mathbf{w}_b = \mathbf{W}^{\mathrm{T}}\mathbf{e}_b, \quad \mathbf{p}^{(j)} = \mathbf{W}^{\mathrm{T}}\mathbf{d}^{(j)} = \mathbf{p}^{(j-1)} + r_b \mathbf{w}_b, \quad \mathbf{c}^{(j)} = \mathbf{W}^{\mathrm{T}}\mathbf{z}^{(j)} = \mathbf{c}^{(j-1)} + \Delta t^{(j-1)}\mathbf{p}^{(j-1)}
\]

- 推导/来源：分别用 \(\mathbf d\)、\(\mathbf z\) 的递推关系左乘 \(\mathbf W^T\)。  
- 公式整体含义：把大维度向量映射到 \(2m\) 维辅助空间递推。

### 公式 3.3-14（改写后的 \(f'_j\)）
\[
f'_j = f'_{j-1} + \Delta t^{(j-1)} f''_{j-1} + r^2_b + r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{z}^{(j)} - r_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{c}^{(j)}
\]

- 推导/来源：把 Eq.23 中 \(\mathbf B\) 用紧凑式替换并用 3.3-13 的变量重写。  
- 公式整体含义：显式避免构造稠密 \(\mathbf B\) 仍可更新系数。

### 公式 3.3-15（改写后的 \(f''_j\)）
\[
f''_j = f''_{j-1} + 2r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{d}^{(j-1)} - 2r_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{p}^{(j-1)} + r^2_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{e}_b - r^2_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{w}_b
\]

- 推导/来源：同上，对 Eq.24 进行紧凑分解代入。  
- 公式整体含义：二次项更新也只需稀疏行操作与 \(2m\) 维小矩阵运算。

### 公式 3.3-16（下一段活动集）
\[
\mathcal{A}(\mathbf{x}^c) = \{i : t_i \leq t^{(j)}, i = 1, \ldots, n\}
\]

- 推导/来源：当搜索推进到第 \(j\) 段后，活动条件同步前移。  
- 公式整体含义：Cauchy 点落入当前段时的活动集判定规则。

---

## 3.4 Subspace minimization

### 公式 3.4-1（活动集等价表示）
\[
\mathcal{A}(\mathbf{x}^c) = \{i : x^c_i = \mathrm{lb}_i\} \cup \{i : x^c_i = \mathrm{ub}_i\} = \{i : t_i \leq t^{(j)}, i = 1, \ldots, n\}
\]

- 推导/来源：由“贴边界”等价于“断点已发生”得到。  
- 公式整体含义：给出几何定义与断点定义之间的一致性。

### 公式 3.4-2（活动变量取值）
\[
x_i = x^c_i = \begin{cases} \mathrm{ub}_i & \text{if } r_i < 0 \\ \mathrm{lb}_i & \text{if } r_i > 0 \end{cases}, \quad \forall i \in \mathcal{A}(\mathbf{x}^c)
\]

- 推导/来源：沿负梯度移动时，\(r_i<0\) 会触上界，\(r_i>0\) 会触下界。  
- 公式整体含义：活动变量已知且固定，不再参与自由求解。

### 公式 3.4-3（Eq.25）
\[
\mathbf{x}^* = \arg \min m_k(\mathbf{x})
\]

- 推导/来源：在活动集已知后继续最小化二次模型。  
- 公式整体含义：定义子空间最优化目标解。

### 公式 3.4-4（活动变量固定约束）
\[
x_i = x^c_i, \quad \forall i \in \mathcal{A}(\mathbf{x}^c)
\]

- 推导/来源：活动集方法核心约束。  
- 公式整体含义：约束降维的第一步。

### 公式 3.4-5（Eq.26，自由变量盒约束）
\[
\mathrm{lb}_i \leq x_i \leq \mathrm{ub}_i, \quad \forall i \notin \mathcal{A}(\mathbf{x}^c)
\]

- 推导/来源：未激活变量仍需满足原盒约束。  
- 公式整体含义：约束降维的第二步。

### 公式 3.4-6（自由变量参数化）
\[
\mathbf{x} = \mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}}
\]

- 推导/来源：\(\mathbf Z_k\) 张成自由子空间，用 \(\hat{\mathbf x}\) 表示自由坐标。  
- 公式整体含义：把原变量问题转成低维自由变量问题。

### 公式 3.4-7（子空间二次模型展开）
\[
\begin{array}{rl}
m_k(\mathbf{x}) &= \Pi_k + \mathbf{r}_k^{\mathrm{T}}(\mathbf{x} - \mathbf{x}_k) + \frac{1}{2}(\mathbf{x} - \mathbf{x}_k)^{\mathrm{T}}\mathbf{B}_k(\mathbf{x} - \mathbf{x}_k) \\
\implies \hat{m}_k(\hat{\mathbf{x}}) &= \Pi_k + \mathbf{r}_k^{\mathrm{T}}(\mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}} - \mathbf{x}_k) + \frac{1}{2}(\mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}} - \mathbf{x}_k)^{\mathrm{T}}\mathbf{B}_k(\mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}} - \mathbf{x}_k) \\
&= \Pi_k + \mathbf{r}_k^{\mathrm{T}}(\mathbf{x}^c - \mathbf{x}_k) + \frac{1}{2} (\mathbf{x}^c - \mathbf{x}_k)^{\mathrm{T}}\mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k) \\
&\quad + [\mathbf{r}_k + \mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k)]^{\mathrm{T}}\mathbf{Z}_k \hat{\mathbf{x}} + \frac{1}{2} \hat{\mathbf{x}}^{\mathrm{T}}\mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k \hat{\mathbf{x}} \\
&= \Pi_k + \mathbf{r}_k^{\mathrm{T}}(\mathbf{x}^c - \mathbf{x}_k) + \frac{1}{2} (\mathbf{x}^c - \mathbf{x}_k)^{\mathrm{T}}\mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k) \\
&\quad + [\mathbf{Z}^{\mathrm{T}}_k [\mathbf{r}_k + \mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k)]]^{\mathrm{T}} \hat{\mathbf{x}} + \frac{1}{2} \hat{\mathbf{x}}^{\mathrm{T}}\mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k \hat{\mathbf{x}}
\end{array}
\]

- 推导/来源：将 3.4-6 代入 Eq.15，逐项展开并整理成对 \(\hat{\mathbf x}\) 的二次型。  
- 公式整体含义：明确了子空间问题的常数项、线性项、二次项。

### 公式 3.4-8（一阶最优条件）
\[
\hat{m}'_k(\hat{\mathbf{x}}) = 0 \implies \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k \hat{\mathbf{x}} = -\mathbf{Z}^{\mathrm{T}}_k [\mathbf{r}_k + \mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k)]
\]

- 推导/来源：对 3.4-7 对 \(\hat{\mathbf x}\) 求导并置零。  
- 公式整体含义：得到自由变量满足的线性方程。

### 公式 3.4-9（Eq.27）
\[
\hat{\mathbf{B}}_k \hat{\mathbf{x}} = -\hat{\mathbf{r}}_k
\]

- 推导/来源：将 3.4-8 左右两项分别记为 \(\hat{\mathbf B}_k\)、\(\hat{\mathbf r}_k\)。  
- 公式整体含义：子空间最小化核心代数系统。

### 公式 3.4-10（紧凑 B 重申）
\[
\mathbf{B}_k = \mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k
\]

- 推导/来源：Eq.20 再次用于 3.4.1。  
- 公式整体含义：说明 \(\hat{\mathbf B}_k\) 也应按紧凑结构计算。

### 公式 3.4-11（Eq.28）
\[
\hat{\mathbf{B}}_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k = \mathbf{Z}^{\mathrm{T}}_k (\mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k) \mathbf{Z}_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{Z}_k - \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k
\]

- 推导/来源：直接把 3.4-10 代入 \(\hat{\mathbf B}_k\) 定义。  
- 公式整体含义：子空间矩阵由“稀疏主项 + 低秩修正”组成。

### 公式 3.4-12（自由变量解）
\[
\hat{\mathbf{x}} = -\hat{\mathbf{B}}^{-1}_k \hat{\mathbf{r}}_k
\]

- 推导/来源：由 Eq.27 直接求逆形式。  
- 公式整体含义：显式表达自由变量解依赖于 \(\hat{\mathbf B}_k^{-1}\)。

### 公式 3.4-13（Eq.29，SMW 逆）
\[
\hat{\mathbf{B}}^{-1}_k = (\hat{\mathbf{B}}^0_k)^{-1} + (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \left[\mathbf{I} - \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \right]^{-1} \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1}
\]

- 推导/来源：把 3.4-11 识别为“可逆矩阵 + 低秩修正”，应用 Sherman–Morrison–Woodbury 公式。  
- 公式整体含义：将大规模逆运算分解成若干次稀疏求解和一个小矩阵求逆。

### 公式 3.4-14（\(\mathbf Z^T\mathbf W\) 展开）
\[
\mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k = \mathbf{Z}^{\mathrm{T}}_k [\mathbf{Y}_k \quad \mathbf{B}^0_k \mathbf{S}_k] = [\mathbf{Z}^{\mathrm{T}}_k \mathbf{y}_{k-m} \cdots \mathbf{Z}^{\mathrm{T}}_k \mathbf{y}_{k-1} \quad \mathbf{Z}^{\mathrm{T}}_k(\mathbf{B}^0_k \mathbf{s}_{k-m}) \cdots \mathbf{Z}^{\mathrm{T}}_k(\mathbf{B}^0_k \mathbf{s}_{k-1})]
\]

- 推导/来源：按列块乘法直接展开。  
- 公式整体含义：说明需要反复稀疏解算的向量列具体来自哪里。

### 公式 3.4-15（小矩阵维度）
\[
\mathbf{I} - \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \in \mathbb{R}^{2m\times 2m}
\]

- 推导/来源：由各乘积维度链直接得到。  
- 公式整体含义：最终昂贵逆只发生在小尺度 \(2m\times2m\) 上。

### 公式 3.4-16（预条件子）
\[
\mathbf{P}^{-1} = \mathrm{ILU}(\hat{\mathbf{B}}^0_k) = \mathrm{ILU}(\mathbf{Z}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{Z}_k)
\]

- 推导/来源：3.4.2 中以 \(\hat{\mathbf B}^0_k\) 的 ILU 作为 \(\hat{\mathbf B}_k\) 的近似逆。  
- 公式整体含义：提升 CG 收敛效率且构造代价可控。

### 公式 3.4-17（增量表示）
\[
\mathbf{x} = \mathbf{x}_k + \Delta \mathbf{x}_k
\]

- 推导/来源：3.4.3 对偶/KKT 推导的变量替换。  
- 公式整体含义：把未知量从绝对坐标转为迭代增量。

### 公式 3.4-18（二次模型增量形式）
\[
m_k(\mathbf{x}) = \Pi_k + \mathbf{r}_k^{\mathrm{T}} \Delta \mathbf{x}_k + \frac{1}{2} \Delta \mathbf{x}_k^{\mathrm{T}} \mathbf{B}_k \Delta \mathbf{x}_k
\]

- 推导/来源：将 3.4-17 代入 Eq.15。  
- 公式整体含义：KKT 形式下的标准二次目标函数。

### 公式 3.4-19（活动约束重述）
\[
x_i = x^c_i, \quad \forall i \in \mathcal{A}_k(\mathbf{x}^c)
\]

- 推导/来源：对偶推导中再次显式写出固定变量条件。  
- 公式整体含义：等式约束来源。

### 公式 3.4-20（Q 投影约束）
\[
\mathbf{Q}^{\mathrm{T}}_k \Delta \mathbf{x}_k = \mathbf{Q}^{\mathrm{T}}_k (\mathbf{x} - \mathbf{x}_k) = \mathbf{Q}^{\mathrm{T}}_k (\mathbf{x}^c - \mathbf{x}_k)
\]

- 推导/来源：将活动分量等式约束写成矩阵形式。  
- 公式整体含义：仅在固定变量子空间上施加约束方程。

### 公式 3.4-21（对偶目标）
\[
\min \left\{ \Pi_k + \mathbf{r}_k^{\mathrm{T}} \Delta \mathbf{x}_k + \frac{1}{2} \Delta \mathbf{x}_k^{\mathrm{T}} \mathbf{B}_k \Delta \mathbf{x}_k \right\}
\]

- 推导/来源：3.4-18 直接作为约束优化目标。  
- 公式整体含义：形成带等式/不等式约束的二次规划。

### 公式 3.4-22（对偶等式约束）
\[
\mathbf{Q}_k^{\mathrm{T}} \Delta \mathbf{x}_k = \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k)
\]

- 推导/来源：同 3.4-20。  
- 公式整体含义：固定变量条件在增量坐标中的线性等式。

### 公式 3.4-23（增量盒约束）
\[
\mathbf{lb} - \mathbf{x}_k \leq \Delta \mathbf{x}_k \leq \mathbf{ub} - \mathbf{x}_k
\]

- 推导/来源：由原变量盒约束减去 \(\mathbf x_k\) 得到。  
- 公式整体含义：增量解也必须保持可行。

### 公式 3.4-24（Eq.30，KKT 系统）
\[
\begin{bmatrix} \mathbf{B}_k & \mathbf{Q}_k \\ \mathbf{Q}_k^{\mathrm{T}} & \mathbf{0} \end{bmatrix} \begin{bmatrix} \Delta \mathbf{x}_k \\ \lambda \end{bmatrix} = \begin{bmatrix} -\mathbf{r}_k \\ \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k) \end{bmatrix}
\]

- 推导/来源：构造拉格朗日函数并写一阶最优条件。  
- 各部分含义：\(\lambda\) 为拉格朗日乘子。  
- 公式整体含义：对偶法下的鞍点线性系统。

### 公式 3.4-25（Eq.31，Schur 补）
\[
(\mathbf{Q}_k^{\mathrm{T}} \mathbf{B}_k^{-1} \mathbf{Q}_k) \lambda = -\mathbf{Q}_k^{\mathrm{T}} \mathbf{B}_k^{-1} \mathbf{r}_k - \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k)
\]

- 推导/来源：由 Eq.30 消去 \(\Delta\mathbf x_k\) 得到 Schur 补方程。  
- 公式整体含义：先求乘子再回代求增量。

### 公式 3.4-26（Eq.32，回代）
\[
\mathbf{B}_k \Delta \mathbf{x}_k = -(\mathbf{r}_k + \mathbf{Q}_k \lambda)
\]

- 推导/来源：Eq.30 第一行回代。  
- 公式整体含义：乘子已知后恢复主变量增量。

### 公式 3.4-27（搜索方向重述）
\[
\mathbf{p}_k = \mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}} - \mathbf{x}_k
\]

- 推导/来源：由子空间解回到全空间后与当前点作差。  
- 公式整体含义：综合“活动变量固定 + 自由变量修正”得到方向。

### 公式 3.4-28（状态更新重述）
\[
\mathbf{x}_{k + 1} = \mathbf{x}_k + \alpha_k \mathbf{p}_k
\]

- 推导/来源：标准线搜索步。  
- 公式整体含义：生成下一次 L-BFGS-B 迭代状态。

### 公式 3.4-29（可行性步长限制）
\[
\mathbf{lb} \leq \mathbf{x}_k + \alpha_k \mathbf{p}_k \leq \mathbf{ub}
\]

- 推导/来源：把更新式代入盒约束。  
- 公式整体含义：线搜索步长必须保证新解不越界。

---

## 公式间主线总结（第三章）

1. 以 Eq.12–Eq.14 建立“能量最小化 + 盒约束 + 线性约束”的完整离散问题。  
2. 以 Eq.15 的二次模型替代原问题局部行为，并通过 Eq.16–Eq.19 的投影与最优性刻画处理不可逆约束。  
3. 以 Eq.20–Eq.21（及其组成块）实现 L-BFGS 紧凑存储，避免稠密 Hessian 内存灾难。  
4. 以 Eq.22–Eq.24 及其紧凑改写，在分段路径上高效定位广义 Cauchy 点和活动集。  
5. 以 Eq.25–Eq.29 在自由子空间求解修正量，并可用 CG+ILU（3.4-16）或 SMW（3.4-13）实现。  
6. 以 Eq.30–Eq.32 说明对偶 Schur 补可行但实现成本高，从而支撑本文采用原始子空间路线。

