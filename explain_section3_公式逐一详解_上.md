# explain.md 第3章公式逐一详解（上）

> 范围：对应 `explain.md` 中 `## 3. Gradient projection based monolithic scheme` 的前半部分（3.1 与 3.2）。  
> 说明：每个公式均先完整显示，再给出“含义 + 必要推导 + 详细解读”。文中符号尽量统一用 `\[ \cdot \]` 的数学渲染格式表示。

---

## 3.1 Algorithm overview（公式逐一解读）

### 公式 3.1（原式(12)）
\[
\left(\pmb{u}_A,d_A\right) = \arg \min \Pi (\pmb {u}_A,d_A) \quad (12)
\]

- 含义：在当前载荷步内，把离散后的位移自由度 \[ \pmb{u}_A \] 与相场自由度 \[ d_A \] 作为联合未知量，求使总势能 \[ \Pi \] 最小的解。
- 必要推导：相场断裂在准静态条件下满足能量最小原理。连续问题离散后，泛函 \[ \Pi(\pmb{u},d) \] 变成有限维函数 \[ \Pi(\pmb{u}_A,d_A) \]，于是得到标准优化问题 \[ \arg\min \].
- 详细解读：这一步把“求解偏微分方程”转写为“受约束最优化”。后续 L-BFGS-B、活动集、线搜索都围绕该目标函数的下降展开。

### 公式 3.2（原式(13)）
\[
d_A^{(n)}\leq d_A\leq 1 \quad (13)
\]

- 含义：相场不可逆约束与物理上界。下界 \[ d_A^{(n)} \] 是上一步历史值，上界为 \[1\]（完全损伤）。
- 必要推导：不可逆性要求 \[ \dot d \ge 0 \Rightarrow d^{(n+1)} \ge d^{(n)} \]；相场定义域是 \[ [0,1] \]，故当前步综合为箱约束。
- 详细解读：这是第3章核心约束。算法上用投影算子把每次迭代都截断到可行域，避免“裂纹愈合”数值误差。

### 公式 3.3（原式(14)）
\[
\left[\pmb {u}_A,d_A\right]^{\mathrm{T}} = \mathbf{C}\left[\pmb {u}_A,d_A\right]^{\mathrm{T}} + \pmb {k}. \quad (14)
\]

- 含义：线性约束统一表达（Dirichlet 边界 + hanging node 约束）。
- 必要推导：有限元约束消元可写成仿射关系 \[ \mathbf{x}=\mathbf{C}\mathbf{x}+\mathbf{k} \]，其中 \[ \mathbf{x}=[\pmb{u}_A,d_A]^T \].
- 详细解读：\(\mathbf{C}\) 负责“从独立自由度映射到全自由度”，\(\mathbf{k}\) 负责非齐次边值贡献。该结构允许在迭代中稳定施加线性约束而不破坏优化框架。

### 公式 3.4（原式(15)）
\[
m_{k}(\pmb{x}) = \Pi_{k} + \pmb{r}_{k}^{\mathrm{T}}(\pmb{x} - \pmb{x}_{k}) + \frac{1}{2} (\pmb{x} - \pmb{x}_{k})^{\mathrm{T}}\mathbf{B}_{k}(\pmb{x} - \pmb{x}_{k}) \quad (15)
\]

- 含义：第 \[k\] 次迭代的局部二次模型。
- 必要推导：对 \[ \Pi(\mathbf{x}) \] 在 \[ \mathbf{x}_k \] 处作二阶近似：常数项 + 一阶项 + 二次项。牛顿法用真 Hessian，这里用拟牛顿矩阵 \[ \mathbf{B}_k \] 近似。
- 详细解读：该模型是 L-BFGS-B 的“工作模型”。后续广义 Cauchy 点与子空间最小化都是在这个模型上求解，兼顾效率和收敛稳健性。

### 公式 3.5（原式(16)）
\[
\mathrm{Proj}_c(x_i,\mathrm{lb}_i,\mathrm{ub}_i) = \left\{ \begin{array}{ll}\mathrm{lb}_i & \mathrm{if}\; x_i < \mathrm{lb}_i,\\ x_i & \mathrm{if}\; \mathrm{lb}_i \leq x_i \leq \mathrm{ub}_i,\\ \mathrm{ub}_i & \mathrm{if}\; x_i > \mathrm{ub}_i, \end{array} \right. \quad (16)
\]

- 含义：分量级箱约束投影。
- 必要推导：欧氏投影 \[ \Pi_{[\mathrm{lb}_i,\mathrm{ub}_i]}(x_i) \] 的闭式形式正是“低于下界取下界，高于上界取上界，中间不变”。
- 详细解读：计算代价 \[ \mathcal{O}(1) \] / 分量，极适合大规模问题。它保证每次更新后变量不离开可行域。

### 公式 3.6（原式(17)）
\[
\pmb {x}(t) = \mathrm{Proj}_c(\pmb {x}_k - t\pmb {r}_k,\mathbf{lb},\mathbf{ub}),\quad t > 0, \quad (17)
\]

- 含义：沿负梯度方向前进并实时投影得到的“分段线性轨迹”。
- 必要推导：无约束最速下降轨迹是 \[ \mathbf{x}_k-t\mathbf{r}_k \]；有箱约束时，需对该点逐分量投影。
- 详细解读：随着 \[ t \] 增大，不同分量会在不同“断点”撞到边界并被冻结，轨迹因此转折，构成后续广义 Cauchy 点搜索路径。

### 公式 3.7（原式(18)）
\[
p_k(t) = m_k(\pmb{x}(t)) = \Pi_k + \pmb{r}_k^\mathrm{T}(\pmb{x}(t) - \pmb{x}_k) + \frac{1}{2} (\pmb{x}(t) - \pmb{x}_k)^\mathrm{T}\mathbf{B}_k(\pmb{x}(t) - \pmb{x}_k). \quad (18)
\]

- 含义：把多变量二次模型限制到路径 \[ \mathbf{x}(t) \] 上，得到一维分段二次函数。
- 必要推导：复合函数 \[ p_k(t)=m_k(\mathbf{x}(t)) \]；由于 \[ \mathbf{x}(t) \] 分段线性，故 \[ p_k(t) \] 分段二次。
- 详细解读：一维化后可高效定位“第一个局部极小点”（广义 Cauchy 点），避免直接在全空间处理活动集切换。

### 公式 3.8（活动集定义）
\[
\mathcal{A}(\pmb {x}) = \{i:x_i = \mathrm{lb}_i\} \cup \{i:x_i = \mathrm{ub}_i\}.
\]

- 含义：处于箱约束边界的分量索引集合。
- 必要推导：约束优化中，等式成立的约束称活动约束。箱约束恰有上下两类边界，因此取并集。
- 详细解读：活动集把变量分成“固定变量（active）”和“自由变量（free）”，这是子空间最小化的核心分解。

### 公式 3.9（子问题目标）
\[
\pmb{x}^* = \arg \min m_k(\pmb {x})
\]

- 含义：在当前二次模型上求近似最优点 \[ \mathbf{x}^* \]。
- 必要推导：L-BFGS-B 每次迭代并不直接最小化原函数，而是先解局部二次近似子问题。
- 详细解读：\(\mathbf{x}^*\) 与广义 Cauchy 点共同定义搜索方向，影响线搜索效率和全局收敛表现。

### 公式 3.10（子问题约束）
\[
x_{i} = x_{i}^{c},\;\forall i\in \mathcal{A}_{k}(\pmb{x}^{c})\quad \text{and}\quad \mathrm{lb}_{i}\leq x_{i}\leq \mathrm{ub}_{i},\;\forall i\notin \mathcal{A}_{k}(\pmb{x}^{c}).
\]

- 含义：活动变量固定在 \[ \mathbf{x}^c \]；自由变量仍受箱约束。
- 必要推导：活动集法经典策略：先固定活动边界，再在自由子空间最小化。
- 详细解读：这样可避免越界并降低问题维度，数值上等价于在可行面上做二次规划。

### 公式 3.11（搜索方向）
\[
\pmb{p}_k = \pmb{x}^* - \pmb{x}_k.
\]

- 含义：由当前点到子问题解的位移向量。
- 必要推导：优化迭代统一写成 \[ \mathbf{x}_{k+1}=\mathbf{x}_k+\alpha_k\mathbf{p}_k \]，因此先定义方向。
- 详细解读：理论上该方向是下降方向（结合模型正定性与活动集处理），使线搜索能找到有效步长。

### 公式 3.12（迭代更新）
\[
\pmb{x}_{k + 1} = \pmb{x}_k + \alpha_k\pmb{p}_k.
\]

- 含义：线搜索更新式。
- 必要推导：一维参数化 \[ \phi(\alpha)=\Pi(\mathbf{x}_k+\alpha\mathbf{p}_k) \] 后，选择满足 Wolfe 条件的 \[ \alpha_k \].
- 详细解读：\(\alpha_k\) 平衡“下降幅度”和“曲率信息可靠性”，过大可能越界或震荡，过小会导致慢收敛。

### 公式 3.13（强 Wolfe：充分下降）
\[
\Pi_{k + 1} = \Pi (\pmb{x}_k + \alpha_k\pmb {p}_k) \leq \Pi_k + c_1\alpha_k\pmb {p}_k^\mathrm{T}\pmb{r}_k
\]

- 含义：Armijo 条件，保证函数值确实下降。
- 必要推导：由一阶近似 \[ \Pi(\mathbf{x}_k+\alpha\mathbf{p}_k)\approx\Pi_k+\alpha\nabla\Pi_k^T\mathbf{p}_k \] 得到可接受下降阈值。
- 详细解读：\(c_1\in(0,1)\) 很小（文中 \(10^{-4}\)），代表“至少达到一小部分线性预测下降”。

### 公式 3.14（强 Wolfe：曲率条件）
\[
\left|\pmb{r}_{k + 1}^\mathrm{T}\pmb {p}_k\right| = \left|\pmb {r}(\pmb{x}_k + \alpha_k\pmb {p}_k)^\mathrm{T}\pmb {p}_k\right| \leq c_2\left|\pmb{r}_k^\mathrm{T}\pmb {p}_k\right|.
\]

- 含义：要求新点沿方向导数绝对值明显减小。
- 必要推导：若导数衰减不足，说明尚未靠近该方向上的“合适”驻点；强 Wolfe 用绝对值增强稳健性。
- 详细解读：该条件和 BFGS 更新正定性密切相关，可避免步长过短导致曲率对 \[ \mathbf{y}_k \] 估计失真。

### 公式 3.15（收敛判据1）
\[
\mathcal{A}_{k + 1}(\pmb{x}^c) = \mathcal{A}_k(\pmb{x}^c).
\]

- 含义：活动集在相邻迭代间不再变化。
- 必要推导：约束识别稳定意味着“哪几个分量贴边界”已确定。
- 详细解读：该判据是受约束问题特有停止指标，和梯度范数联合使用更可靠。

### 公式 3.16（收敛判据2）
\[
\| \mathrm{Proj}_C(\mathbf{x}_{k + 1} - \mathbf{r}_{k + 1},\mathbf{lb},\mathbf{ub}) - \mathbf{x}_{k + 1}\|_{2} < \mathrm{tol}.
\]

- 含义：投影梯度范数足够小。
- 必要推导：无约束时驻点条件是 \[ \nabla\Pi=0 \]；有箱约束时等价推广为“投影后位移接近 0”。
- 详细解读：它同时检测站位性与可行性，是 L-BFGS-B 常用一阶最优性度量。

### 公式 3.17（收敛判据3）
\[
\| \Delta \mathbf{x}_{k + 1}\|_{2} = \| \mathbf{x}_{k + 1} - \mathbf{x}_{k}\|_{2} < \mathrm{tol}.
\]

- 含义：迭代增量很小。
- 必要推导：若解更新趋近于零，说明算法已到稳定邻域。
- 详细解读：与判据2配合可避免“梯度小但步长异常”或“步长小但未最优”的误判。

### 公式 3.18（原式(19)）
\[
\mathrm{Proj}_C(x_i - r_i,\mathrm{lb}_i,\mathrm{ub}_i) - x_i = \left\{ \begin{array}{ll}\mathrm{lb}_i - x_i & \mathrm{if}\; x_i - r_i < \mathrm{lb}_i,\\ - r_i & \mathrm{if}\; \mathrm{lb}_i \leq x_i - r_i \leq \mathrm{ub}_i,\\ \mathrm{ub}_i - x_i & \mathrm{if}\; x_i - r_i > \mathrm{ub}_i. \end{array} \right. \quad (19)
\]

- 含义：投影梯度在单分量上的显式表达。
- 必要推导：将式(16)代入 \[ \mathrm{Proj}_C(x_i-r_i)-x_i \]，按三段区间逐段展开即可。
- 详细解读：若试探点仍在区间内，则该量就是 \[ -r_i \]；若越界，则变成到边界的距离，精确刻画约束违反方向。

### 公式 3.19（Comment 1 约束形式）
\[
\mathbf{x} = \mathbf{C}\mathbf{x} + \mathbf{k}.
\]

- 含义：把所有线性约束写成统一仿射方程。
- 必要推导：由有限元约束矩阵装配可得该表达；它和式(14)同构，只是记号更紧凑。
- 详细解读：后续线性系统改写、解回代都依赖这个结构。

### 公式 3.20（约束后线性系统）
\[
\left(\mathbf{C}^{\mathrm{T}}\mathbf{A}\mathbf{C} + \mathbf{I}_{d_{c}}\right)\hat{\mathbf{x}} = \mathbf{C}^{\mathrm{T}}(\mathbf{b} - \mathbf{A}\mathbf{k})
\]

- 含义：在受约束空间中求解等价系统。
- 必要推导：代入 \[ \mathbf{x}=\mathbf{C}\hat{\mathbf{x}}+\mathbf{k} \] 到 \[ \mathbf{A}\mathbf{x}=\mathbf{b} \]，左乘 \[ \mathbf{C}^T \] 并添加 \[ \mathbf{I}_{d_c} \] 以处理受限自由度。
- 详细解读：该形式在工程实现中常用于把“约束处理”与“线性求解器”解耦。

### 公式 3.21（回代恢复）
\[
\mathbf{x} = \mathbf{C}\hat{\mathbf{x}} + \mathbf{k}.
\]

- 含义：由辅助变量 \[ \hat{\mathbf{x}} \] 还原真实自由度解。
- 必要推导：直接来自仿射约束定义。
- 详细解读：这一步保证最终解自动满足边界条件与 hanging-node 关系。

### 公式 3.22（约束指示对角阵）
\[
(\mathbf{I}_{d_c})_{ii} = \begin{cases} 1 & \text{if DoF } i \text{ is constrained}, \\ 0 & \text{otherwise}, \end{cases}
\]

- 含义：对受约束自由度打“1”的对角掩码矩阵。
- 必要推导：以指示函数构造稀疏对角矩阵。
- 详细解读：它在实现中用于稳定数值系统并清晰标记约束 DoF。

---

## 3.2 Compact representation of limited-memory BFGS matrix（公式逐一解读）

### 公式 3.23（位移与梯度差分）
\[
\mathbf{s}_k = \mathbf{x}_{k + 1} - \mathbf{x}_k,\quad \mathbf{y}_k = \mathbf{r}_{k + 1} - \mathbf{r}_k,
\]

- 含义：拟牛顿更新的两组核心向量对。
- 必要推导：\(\mathbf{s}_k\) 表示步长位移，\(\mathbf{y}_k\) 表示梯度变化，二者共同编码局部曲率。
- 详细解读：L-BFGS 只保存最近若干 \[ (\mathbf{s}_k,\mathbf{y}_k) \]，以低存储逼近 Hessian 信息。

### 公式 3.24（标准 BFGS 更新）
\[
\mathbf{B}_{k + 1} = \mathbf{B}_k - \frac{\mathbf{B}_k\mathbf{s}_k\mathbf{s}_k^{\mathrm{T}}\mathbf{B}_k}{\mathbf{s}_k^{\mathrm{T}}\mathbf{B}_k\mathbf{s}_k} + \frac{\mathbf{y}_k\mathbf{y}_k^{\mathrm{T}}}{\mathbf{y}_k^{\mathrm{T}}\mathbf{s}_k}.
\]

- 含义：用两次秩一修正更新 Hessian 近似。
- 必要推导：满足拟牛顿割线条件 \[ \mathbf{B}_{k+1}\mathbf{s}_k=\mathbf{y}_k \] 且尽量接近 \[ \mathbf{B}_k \]（在特定矩阵度量下）可导出该式。
- 详细解读：第一项沿 \[ \mathbf{s}_k \] 方向移除旧曲率，第二项注入新曲率。

### 公式 3.25（曲率条件）
\[
\mathbf{s}_k^{\mathrm{T}}\mathbf{y}_k > 0,
\]

- 含义：保证 BFGS 更新后正定性。
- 必要推导：BFGS 正定性理论要求分母 \[ \mathbf{y}_k^T\mathbf{s}_k \] 为正。
- 详细解读：与线搜索（强 Wolfe）联动，通常可自然满足该条件。

### 公式 3.26（S 矩阵）
\[
\mathbf{S}_k = [\mathbf{s}_{k-m} \cdots \mathbf{s}_{k-1}]
\]

- 含义：保存最近 \[m\] 个位移向量列块。
- 必要推导：L-BFGS 截断历史，只取近邻曲率信息。
- 详细解读：维度 \[ n\times m \]，其中 \[ n \] 是总自由度，\(m\ll n\)。

### 公式 3.27（Y 矩阵）
\[
\mathbf{Y}_k = [\mathbf{y}_{k-m} \cdots \mathbf{y}_{k-1}].
\]

- 含义：保存最近 \[m\] 个梯度增量列块。
- 必要推导：与 \[ \mathbf{S}_k \] 成对使用构建紧凑表示。
- 详细解读：\((\mathbf{S}_k,\mathbf{Y}_k)\) 是 L-BFGS“记忆体”。

### 公式 3.28（原式(20) 紧凑表达）
\[
\mathbf{B}_k = \mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k, \quad (20)
\]

- 含义：不显式组装稠密 \[ \mathbf{B}_k \] 的等价表达。
- 必要推导：Byrd 等人将多次 BFGS 更新代数重排，得到“初始矩阵 + 低秩修正”形式。
- 详细解读：关键收益是内存从 \[ \mathcal{O}(n^2) \] 降到 \[ \mathcal{O}(nm+m^2) \]。

### 公式 3.29（W 矩阵）
\[
\mathbf{W}_k = [\mathbf{Y}_k \quad \mathbf{B}^0_k \mathbf{S}_k] \in \mathbb{R}^{n\times 2m}
\]

- 含义：把历史梯度增量与初始矩阵作用后的位移增量并排拼接。
- 必要推导：紧凑表达中低秩项由这两个子空间共同张成。
- 详细解读：后续矩阵-向量乘只需 \[ \mathbf{W}_k \] 与 \[ \mathbf{M}_k \]，无需构造 \[ \mathbf{B}_k \]。

### 公式 3.30（M 矩阵）
\[
\mathbf{M}_k = \begin{bmatrix} -\mathbf{D}_k & \mathbf{L}^{\mathrm{T}}_k \\ \mathbf{L}_k & \mathbf{S}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{S}_k \end{bmatrix}^{-1} \in \mathbb{R}^{2m\times 2m}.
\]

- 含义：低维（\(2m\times2m\)）核心矩阵。
- 必要推导：由 BFGS 多步更新的块代数关系整理并取逆得到。
- 详细解读：尽管 \[ \mathbf{B}_k \] 维度高达 \[ n\times n \]，关键难点被压缩到小矩阵运算。

### 公式 3.31（D 块）
\[
\mathbf{D}_k = \mathrm{diag}\{\mathbf{s}^{\mathrm{T}}_{k-m}\mathbf{y}_{k-m}, \ldots, \mathbf{s}^{\mathrm{T}}_{k-1}\mathbf{y}_{k-1}\} \in \mathbb{R}^{m\times m}
\]

- 含义：每个历史对的曲率标量构成的对角阵。
- 必要推导：取各 \[ \mathbf{s}_i^T\mathbf{y}_i \] 放置在对角线上。
- 详细解读：该块直接反映每次迭代沿步长方向采样到的局部二阶信息。

### 公式 3.32（L 块）
\[
\mathbf{L}_k \in \mathbb{R}^{m\times m}, \quad (\mathbf{L}_k)_{i,j} = \begin{cases} \mathbf{s}^{\mathrm{T}}_{i+k-m-1}\mathbf{y}_{j+k-m-1} & \text{if } i > j, \\ 0 & \text{if } i \leq j. \end{cases}
\]

- 含义：严格下三角块，记录“跨历史对”的互相关曲率。
- 必要推导：按索引顺序把 \[ \mathbf{s}_i^T\mathbf{y}_j \] 的下三角部分提取出来。
- 详细解读：与 \[ \mathbf{D}_k \] 一起构成紧凑表达所需全部二阶统计量。

### 公式 3.33（原式(21) 初始矩阵）
\[
\mathbf{B}^0_k = \hat{\mathbf{K}}^{(k)} = \begin{bmatrix} \mathbf{K}_{uu} & \mathbf{0} \\ \mathbf{0} & \mathbf{K}_{dd} \end{bmatrix}^{(k)}. \quad (21)
\]

- 含义：L-BFGS 的初始近似采用分块对角刚度矩阵。
- 必要推导：在本问题中位移与相场可构造各自正定切线刚度块，取块对角作为稳定初值。
- 详细解读：这一步非常关键：\(\mathbf{B}^0_k\) 既保留有限元稀疏结构，又提供正确尺度，使低秩校正在物理上更可信。

---

## 本文件小结

- 本文件覆盖了第3章中 3.1 与 3.2 的全部展示公式（含 Comment 1 的实现公式）。
- 下一文件将继续 3.3 与 3.4 的全部公式，并保持同样“先公式、后推导与解读”的结构。
