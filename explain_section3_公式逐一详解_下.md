# explain.md 第3章公式逐一详解（下）

> 范围：对应 `explain.md` 中第3章后半部分（3.3 与 3.4）全部展示公式。  
> 说明：每个公式均完整显示在前，并给出含义、必要推导与详细解读。

---

## 3.3 Generalized Cauchy point（公式逐一解读）

### 公式 3.34（断点定义）
\[
t_i = \begin{cases} (x^0_i - \mathrm{ub}_i)/r_i & \text{if } r_i < 0, \\ (x^0_i - \mathrm{lb}_i)/r_i & \text{if } r_i > 0, \\ +\infty & \text{if } r_i = 0. \end{cases}
\]

- 含义：第 \[i\] 个分量沿投影前下降方向首次触碰边界所需“时间”。
- 推导：无约束轨迹是 \[ x_i(t)=x_i^0-tr_i \]。令其等于上/下界并解 \[ t \]，按 \[ r_i \] 符号分段得到该式；\(r_i=0\) 永不触边，故 \[ +\infty \]。
- 解读：所有 \[ t_i \] 排序后决定活动集切换顺序。

### 公式 3.35（单分量投影轨迹）
\[
x_i(t) = x_i^{0} - \min \{t,t_i\} r_i.
\]

- 含义：在断点前线性移动，断点后冻结在边界。
- 推导：\(t<t_i\) 时是原下降；\(t\ge t_i\) 时增量固定为 \[ t_i r_i \]。
- 解读：这一“截断线性”结构使整体轨迹分段线性。

### 公式 3.36（区间内向量轨迹）
\[
\mathbf{x}(t) = \mathbf{x}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)},
\]

- 含义：在第 \[j\] 个区间内，轨迹是仿射直线。
- 推导：当活动集不变时，方向常量 \[ \mathbf{d}^{(j-1)} \] 不变，故线性参数化成立。
- 解读：把多次折线段分解后，每段可做标准一元二次最小化。

### 公式 3.37（区间内一元二次模型）
\[
\begin{array}{rl} 
p(t) &= m(\mathbf{x}(t)) = \Pi + \mathbf{r}^{\mathrm{T}}(\mathbf{x}(t) - \mathbf{x}^0) + \frac{1}{2} (\mathbf{x}(t) - \mathbf{x}^0)^{\mathrm{T}}\mathbf{B}(\mathbf{x}(t) - \mathbf{x}^0)\\
&= f_{j-1} + f'_{j-1}\Delta t + \frac{1}{2} f''_{j-1}\Delta t^2,
\end{array}
\]

- 含义：在固定区间内，模型可写为标准二次多项式。
- 推导：把式(3.36)代入二次型 \[ m(\mathbf{x}) \]，按 \[ \Delta t \] 收集常数/一次/二次项。
- 解读：只需计算 \[ f'_{j-1},f''_{j-1} \] 即可判断该段极小点位置。

### 公式 3.38（二阶导正性）
\[
f_{j-1}^{\prime \prime} = \mathbf{d}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j-1)} > 0.
\]

- 含义：区间模型在该方向上是凸的。
- 推导：\(\mathbf{B}\) 正定（或在相关子空间正定）时，任意非零方向二次型为正。
- 解读：确保每段最多一个局部极小点。

### 公式 3.39（区间极小点）
\[
\hat{p}'(\Delta t) = 0 \Rightarrow \Delta t^* = -f_{j-1}' / f_{j-1}^{\prime \prime}.
\]

- 含义：二次函数驻点闭式解。
- 推导：\(\hat p'(\Delta t)=f'_{j-1}+f''_{j-1}\Delta t\)，令其为零得到。
- 解读：若 \[ \Delta t^* \] 落在当前区间，则得到广义 Cauchy 点；否则继续下一区间。

### 公式 3.40（活动集阶段判定）
\[
\mathcal{A}(\mathbf{x}^c) = \{i : t_i \leq t^{(j-1)}, i = 1, \ldots, n\}.
\]

- 含义：在时刻 \[ t^{(j-1)} \] 前已触边分量构成活动集。
- 推导：由断点定义直接得到。
- 解读：该集合随 \[ t \] 单调增加。

### 公式 3.41（跨断点更新）
\[
\Delta t^{(j-1)} = t^{(j)} - t^{(j-1)}, \quad \mathbf{x}^{(j)} = \mathbf{x}^{(j-1)} + \Delta t^{(j-1)}\mathbf{d}^{(j-1)}.
\]

- 含义：从一个断点推进到下一个断点。
- 推导：区间线性关系直接积分得到。
- 解读：这是 piecewise 搜索的主循环更新式。

### 公式 3.42（方向修正）
\[
\mathbf{d}^{(j)} = \mathbf{d}^{(j-1)} + r_b \mathbf{e}_b,
\]

- 含义：当第 \[b\] 个分量变为活动后，将其方向分量置零的等价写法。
- 推导：原方向该分量为 \[ -r_b \]，加 \[ r_b \] 后得到0。
- 解读：\(\mathbf{e}_b\) 是标准基向量，用于低成本增量更新。

### 公式 3.43（一次导递推）
\[
f'_j = \mathbf{r}^{\mathrm{T}}\mathbf{d}^{(j)} + \mathbf{d}^{(j)}{}^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j)} = f'_{j-1} + \Delta t^{(j-1)} f''_{j-1} + r^2_b + r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j)}.
\]

- 含义：跨断点后新区间一阶系数的递推更新。
- 推导：将 \[ \mathbf{d}^{(j)} \] 与 \[ \mathbf{z}^{(j)} \] 的增量表达代入并整理。
- 解读：避免每段从头计算，显著降低复杂度。

### 公式 3.44（二次导递推）
\[
f''_j = \mathbf{d}^{(j)}{}^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j)} = f''_{j-1} + 2r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j-1)} + r^2_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{e}_b.
\]

- 含义：跨断点后曲率项的递推公式。
- 推导：把 \[ \mathbf{d}^{(j)}=\mathbf{d}^{(j-1)}+r_b\mathbf{e}_b \] 代入二次型并展开。
- 解读：对应“一个分量被冻结”带来的曲率重分配。

### 公式 3.45（紧凑B再写）
\[
\mathbf{B} = \mathbf{B}^0 - \mathbf{W}\mathbf{M}\mathbf{W}^{\mathrm{T}}.
\]

- 含义：将3.2结论用于3.3中快速计算 \[ \mathbf{B} \] 作用。
- 推导：直接来自紧凑表示式(20)。
- 解读：后续所有 \[ \mathbf{e}_b^T\mathbf{B}\cdot \] 均可改写为稀疏+低秩计算。

### 公式 3.46（辅助向量）
\[
\mathbf{w}_b = \mathbf{W}^{\mathrm{T}}\mathbf{e}_b, \quad \mathbf{p}^{(j)} = \mathbf{W}^{\mathrm{T}}\mathbf{d}^{(j)} = \mathbf{W}^{\mathrm{T}}(\mathbf{d}^{(j-1)} + r_b \mathbf{e}_b).
\]

- 含义：把标准基与方向向量投影到低维 \[2m\] 子空间。
- 推导：由线性映射 \[ \mathbf{W}^T \] 直接定义。
- 解读：后续与 \[ \mathbf{M} \] 配合，实现低维快速更新。

### 公式 3.47（紧凑形式的一阶导递推）
\[
f'_j = f'_{j-1} + \Delta t^{(j-1)} f''_{j-1} + r^2_b + r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{z}^{(j)} - r_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{c}^{(j)}
\]

- 含义：用 \[ \mathbf{B}^0,\mathbf{W},\mathbf{M} \] 替代显式 \[ \mathbf{B} \] 的递推版。
- 推导：把式3.45代入式3.43中 \[ \mathbf{B} \] 相关项并重排。
- 解读：核心好处是避免任何 \[ n\times n \] 稠密矩阵操作。

### 公式 3.48（紧凑形式的二阶导递推）
\[
f''_j = f''_{j-1} + 2r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{d}^{(j-1)} - 2r_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{p}^{(j-1)} + r^2_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{e}_b - r^2_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M}\mathbf{w}_b.
\]

- 含义：曲率递推的紧凑低秩实现。
- 推导：同理，将式3.45代入式3.44并将项按 \[ \mathbf{B}^0 \] 与低秩修正拆分。
- 解读：对大型有限元系统，这是性能关键式。

### 公式 3.49（活动集更新）
\[
\mathcal{A}(\mathbf{x}^c) = \{i : t_i \leq t^{(j)}, i = 1, \ldots, n\}.
\]

- 含义：推进到新断点后更新活动集。
- 推导：阈值时间从 \[ t^{(j-1)} \] 更新为 \[ t^{(j)} \]。
- 解读：活动集随区间推进逐步扩展直到找到 \[ \mathbf{x}^c \]。

### 公式 3.50（活动集等价定义）
\[
\mathcal{A}(\mathbf{x}^c) = \{i : x^c_i = \mathrm{lb}_i\} \cup \{i : x^c_i = \mathrm{ub}_i\} = \{i : t_i \leq t^{(j)}, i = 1, \ldots, n\}.
\]

- 含义：几何定义与断点定义一致。
- 推导：由 \[ x_i(t) \] 与 \[ t_i \] 的构造可得一一对应。
- 解读：实现上既可按值判定，也可按断点判定。

### 公式 3.51（活动变量取值）
\[
x_i = x^c_i = \begin{cases} \mathrm{ub}_i & \text{if } r_i < 0 \\ \mathrm{lb}_i & \text{if } r_i > 0 \end{cases}, \quad \forall i \in \mathcal{A}(\mathbf{x}^c).
\]

- 含义：活动变量最终固定在哪个边界由梯度符号决定。
- 推导：若 \[ r_i<0 \]，无约束方向使 \[ x_i \] 增大，先碰上界；反之先碰下界。
- 解读：这是“广义 Cauchy 点 = 首个局部极小 + 边界锁定”的具体落点规则。

---

## 3.4 Subspace minimization（公式逐一解读）

### 公式 3.52（原式(25)）
\[
\mathbf{x}^* = \arg \min m_k(\mathbf{x}) \quad (25)
\]

- 含义：在 Cauchy 点识别后，进一步做子空间二次最小化。
- 推导：L-BFGS-B 标准两阶段：先 Cauchy 点，再子空间修正。
- 解读：该步用于提升仅沿投影梯度方向搜索的精度。

### 公式 3.53（活动变量等式约束）
\[
x_i = x^c_i, \quad \forall i \in \mathcal{A}(\mathbf{x}^c)
\]

- 含义：活动变量固定。
- 推导：活动集法的可行面约束。
- 解读：把问题维度降到自由变量子空间。

### 公式 3.54（原式(26) 自由变量箱约束）
\[
\mathrm{lb}_i \leq x_i \leq \mathrm{ub}_i, \quad \forall i \notin \mathcal{A}(\mathbf{x}^c). \quad (26)
\]

- 含义：自由变量仍须在箱内。
- 推导：原始全局约束在子空间中的继承。
- 解读：实践中常先解无界子问题再截断回箱内。

### 公式 3.55（自由变量参数化）
\[
\mathbf{x} = \mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}},
\]

- 含义：用基矩阵 \[ \mathbf{Z}_k \] 表示自由子空间坐标。
- 推导：\(\mathbf{Z}_k\) 的列对应自由变量单位基；活动变量对应分量为0。
- 解读：把原问题转为低维变量 \[ \hat{\mathbf{x}} \] 问题。

### 公式 3.56（代入后的子空间模型）
\[
\begin{array}{rl}
m_k(\mathbf{x}) &= \Pi_k + \mathbf{r}_k^{\mathrm{T}}(\mathbf{x} - \mathbf{x}_k) + \frac{1}{2}(\mathbf{x} - \mathbf{x}_k)^{\mathrm{T}}\mathbf{B}_k(\mathbf{x} - \mathbf{x}_k) \\
&= \Pi_k + \mathbf{r}_k^{\mathrm{T}}(\mathbf{x}^c - \mathbf{x}_k + \mathbf{Z}_k\hat{\mathbf{x}}) + \frac{1}{2}(\mathbf{x}^c - \mathbf{x}_k + \mathbf{Z}_k\hat{\mathbf{x}})^{\mathrm{T}}\mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k + \mathbf{Z}_k\hat{\mathbf{x}}),
\end{array}
\]

- 含义：将目标显式改写为 \[ \hat{\mathbf{x}} \] 的二次函数。
- 推导：直接代入式3.55并展开。
- 解读：后续一阶最优条件可得到线性方程组。

### 公式 3.57（子空间一阶条件）
\[
\hat{m}'_k(\hat{\mathbf{x}}) = 0 \implies \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k \hat{\mathbf{x}} = -\mathbf{Z}^{\mathrm{T}}_k [\mathbf{r}_k + \mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k)].
\]

- 含义：自由子空间的正规方程。
- 推导：对式3.56对 \[ \hat{\mathbf{x}} \] 求导并令零。
- 解读：这是 primal 方案的核心线性系统。

### 公式 3.58（原式(27)）
\[
\hat{\mathbf{B}}_k \hat{\mathbf{x}} = -\hat{\mathbf{r}}_k. \quad (27)
\]

- 含义：对式3.57做记号压缩。
- 推导：定义 \[ \hat{\mathbf{B}}_k=\mathbf{Z}_k^T\mathbf{B}_k\mathbf{Z}_k,\ \hat{\mathbf{r}}_k=\mathbf{Z}_k^T[\cdots] \] 即得。
- 解读：便于调用 PCG 等迭代解法。

### 公式 3.59（B_k 再次紧凑表达）
\[
\mathbf{B}_k = \mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k.
\]

- 含义：把3.2紧凑形式移入子空间系统。
- 推导：直接代入定义。
- 解读：是快速算子实现基础。

### 公式 3.60（子空间矩阵展开）
\[
\hat{\mathbf{B}}_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k = \mathbf{Z}^{\mathrm{T}}_k (\mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k) \mathbf{Z}_k.
\]

- 含义：\(\hat{\mathbf{B}}_k\) 由稀疏主项与低秩修正组成。
- 推导：式3.59左右乘 \[ \mathbf{Z}_k^T,\mathbf{Z}_k \]。
- 解读：这使子空间乘法仍可分解为“稀疏乘 + 小矩阵乘”。

### 公式 3.61（显式解）
\[
\hat{\mathbf{x}} = -\hat{\mathbf{B}}^{-1}_k \hat{\mathbf{r}}_k,
\]

- 含义：若可显式求逆时的闭式。
- 推导：线性系统直接解。
- 解读：实际大规模中通常不显式求逆，而用迭代法。

### 公式 3.62（逆矩阵表达）
\[
\hat{\mathbf{B}}^{-1}_k = (\hat{\mathbf{B}}^0_k)^{-1} + (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \left[\mathbf{I} - \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k\right]^{-1} \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1},
\]

- 含义：子空间逆的低秩修正形式（Woodbury 思想）。
- 推导：对 \[ \hat{\mathbf{B}}_k=\hat{\mathbf{B}}_k^0-\mathbf{U}\mathbf{M}\mathbf{U}^T \] 应用矩阵求逆引理。
- 解读：把大维求逆问题转为 \[2m\] 维小系统求逆。

### 公式 3.63（Z^T W 展开）
\[
\mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k = \mathbf{Z}^{\mathrm{T}}_k [\mathbf{Y}_k \quad \mathbf{B}^0_k \mathbf{S}_k] = [\mathbf{Z}^{\mathrm{T}}_k \mathbf{y}_{k-m} \cdots \mathbf{Z}^{\mathrm{T}}_k \mathbf{y}_{k-1} \quad \mathbf{Z}^{\mathrm{T}}_k(\mathbf{B}^0_k \mathbf{s}_{k-m}) \cdots \mathbf{Z}^{\mathrm{T}}_k(\mathbf{B}^0_k \mathbf{s}_{k-1})].
\]

- 含义：子空间与历史曲率信息的耦合矩阵。
- 推导：按列块线性分配展开。
- 解读：便于逐列增量更新，不必整体重组。

### 公式 3.64（小矩阵维度）
\[
\mathbf{I} - \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \in \mathbb{R}^{2m\times 2m}
\]

- 含义：关键求逆对象仅是 \[ 2m \times 2m \]。
- 推导：根据各矩阵维度相乘可得。
- 解读：说明算法规模主要受 \[ m \] 影响而非 \[ n \]。

### 公式 3.65（预条件子）
\[
\mathbf{P}^{-1} = \mathrm{ILU}(\hat{\mathbf{B}}^0_k) = \mathrm{ILU}(\mathbf{Z}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{Z}_k).
\]

- 含义：以初始子空间矩阵的不完全LU作为预条件。
- 推导：\(\hat{\mathbf{B}}^0_k\) 稀疏且正定性好，适合 ILU 构造。
- 解读：显著改善 PCG 收敛速度。

### 公式 3.66（增量变量）
\[
\mathbf{x} = \mathbf{x}_k + \Delta \mathbf{x}_k.
\]

- 含义：dual 推导中改用增量未知量。
- 推导：标准变量替换便于写 KKT 系统。
- 解读：使线性项和约束写法更对称。

### 公式 3.67（增量形式目标）
\[
m_k(\mathbf{x}) = \Pi_k + \mathbf{r}_k^{\mathrm{T}} \Delta \mathbf{x}_k + \frac{1}{2} \Delta \mathbf{x}_k^{\mathrm{T}} \mathbf{B}_k \Delta \mathbf{x}_k.
\]

- 含义：以 \[ \Delta\mathbf{x}_k \] 为变量的二次模型。
- 推导：将式3.66代入原模型并抵消常量项。
- 解读：更适合构造拉格朗日函数。

### 公式 3.68（活动等式约束）
\[
x_i = x^c_i, \quad \forall i \in \mathcal{A}_k(\mathbf{x}^c),
\]

- 含义：dual 方案下的等式约束。
- 推导：与 primal 一致，只是求解路径不同。
- 解读：通过乘子显式施加。

### 公式 3.69（Q 矩阵约束写法）
\[
\mathbf{Q}^{\mathrm{T}}_k \Delta \mathbf{x}_k = \mathbf{Q}^{\mathrm{T}}_k (\mathbf{x} - \mathbf{x}_k) = \mathbf{Q}^{\mathrm{T}}_k (\mathbf{x}^c - \mathbf{x}_k).
\]

- 含义：把活动约束压缩成矩阵等式。
- 推导：\(\mathbf{Q}_k\) 的列挑选活动分量单位基，故左侧提取对应分量增量。
- 解读：便于统一写成 KKT 块系统。

### 公式 3.70（dual 子问题）
\[
\min \left\{ \Pi_k + \mathbf{r}_k^{\mathrm{T}} \Delta \mathbf{x}_k + \frac{1}{2} \Delta \mathbf{x}_k^{\mathrm{T}} \mathbf{B}_k \Delta \mathbf{x}_k \right\}
\]

- 含义：dual 视角下同一二次目标。
- 推导：目标不变，仅约束处理方式改变。
- 解读：后续用拉格朗日乘子构造联立方程。

### 公式 3.71（dual 等式约束）
\[
\mathbf{Q}_k^{\mathrm{T}} \Delta \mathbf{x}_k = \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k)
\]

- 含义：式3.69的简写。
- 推导：直接沿用。
- 解读：活动分量增量被固定为已知值。

### 公式 3.72（增量箱约束）
\[
\mathbf{lb} - \mathbf{x}_k \leq \Delta \mathbf{x}_k \leq \mathbf{ub} - \mathbf{x}_k.
\]

- 含义：原箱约束转成对增量的上下界。
- 推导：由 \[ \mathbf{x}=\mathbf{x}_k+\Delta\mathbf{x}_k \] 直接平移。
- 解读：实际实现中常在求得候选步后再做截断投影。

### 公式 3.73（KKT块系统）
\[
\begin{bmatrix} \mathbf{B}_k & \mathbf{Q}_k \\ \mathbf{Q}_k^{\mathrm{T}} & \mathbf{0} \end{bmatrix} \begin{bmatrix} \Delta \mathbf{x}_k \\ \lambda \end{bmatrix} = \begin{bmatrix} -\mathbf{r}_k \\ \mathbf{Q}^{\mathrm{T}}_k (\mathbf{x}^c - \mathbf{x}_k) \end{bmatrix}.
\]

- 含义：一阶最优条件与等式约束组成的 saddle-point 系统。
- 推导：构造拉格朗日函数 \[ \mathcal{L}=m_k+\lambda^T(\mathbf{Q}^T\Delta\mathbf{x}-\mathbf{q}) \]，对 \[ \Delta\mathbf{x},\lambda \] 求导为零。
- 解读：dual 方法数值稳定性取决于该系统求解策略。

### 公式 3.74（Schur 补方程）
\[
(\mathbf{Q}_k^{\mathrm{T}} \mathbf{B}_k^{-1} \mathbf{Q}_k) \lambda = -\mathbf{Q}_k^{\mathrm{T}} \mathbf{B}_k^{-1} \mathbf{r}_k - \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k).
\]

- 含义：消去 \[ \Delta\mathbf{x}_k \] 后得到乘子方程。
- 推导：由式3.73第一行得 \[ \Delta\mathbf{x}_k=-\mathbf{B}^{-1}(\mathbf{r}+\mathbf{Q}\lambda) \]，代入第二行。
- 解读：维度只与活动约束数相关，适于约束数较少情形。

### 公式 3.75（原式(32)）
\[
\mathbf{B}_k \Delta \mathbf{x}_k = -(\mathbf{r}_k + \mathbf{Q}_k \lambda). \quad (32)
\]

- 含义：已知乘子后回代求增量。
- 推导：即 KKT 第一行重写。
- 解读：可继续利用 L-BFGS 紧凑算子求解。

### 公式 3.76（最终方向）
\[
\mathbf{p}_k = \mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}} - \mathbf{x}_k.
\]

- 含义：由 Cauchy 点 + 子空间修正构成最终搜索方向。
- 推导：\(\mathbf{x}^*=\mathbf{x}^c+\mathbf{Z}_k\hat{\mathbf{x}}\) 后减去当前点。
- 解读：这是 L-BFGS-B 与纯投影梯度法的关键区别（多了子空间二次修正）。

### 公式 3.77（线搜索更新）
\[
\mathbf{x}_{k + 1} = \mathbf{x}_k + \alpha_k \mathbf{p}_k,
\]

- 含义：按步长沿最终方向更新。
- 推导：标准线搜索框架。
- 解读：\(\alpha_k\) 仍由强 Wolfe 条件确定。

### 公式 3.78（步长可行性）
\[
\mathbf{lb} \leq \mathbf{x}_k + \alpha_k \mathbf{p}_k \leq \mathbf{ub}.
\]

- 含义：线搜索时必须保持箱约束可行。
- 推导：原始约束对更新式直接施加。
- 解读：实践中可通过可行步长上界裁剪 + Wolfe 搜索联合实现。

---

## 本文件小结

- 本文件覆盖了第3章中 3.3 与 3.4 的全部展示公式，并保持逐式解释。
- 与上篇合并后，已完成对 `explain.md` 第3章所有可见公式的完整解读与必要推导说明。
