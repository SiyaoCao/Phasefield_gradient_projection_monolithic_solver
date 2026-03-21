# explain.md 公式逐一解释、必要推导与内在联系

> 说明：本文档按 `explain.md` 中公式出现顺序编号为 F1–F110（对应我对原文全部显示公式的顺序提取）。
> 
> 每条都包含四部分：**公式完整显示**、**含义**、**必要推导/变形**、**与其他公式联系（被引用公式可在对应 Fx 小节看到完整表达）**。

---

## 0. 总体逻辑主线（先看这一节）

整篇论文公式可归纳为 6 条连续链路：

1. **物理建模链**：F1–F4，定义相场断裂总势能与不可逆约束。
2. **本构与变分链**：F8–F19，从能量分裂到应力、切线模量，再到弱式残量。
3. **有限元离散链**：F20–F27，把连续问题离散成残量-切线系统。
4. **约束优化链**：F28–F49，把不可逆性写成盒约束并通过投影/活动集处理。
5. **L-BFGS-B 数值链**：F50–F105，构造近似 Hessian、广义 Cauchy 点、子空间最小化、线搜索更新。
6. **算例与附录链**：F106–F110，裂纹能监测、扭转边界、收敛准则、SMW 公式。

因此，内在关系是：

\[
\text{能量泛函 }(F1,F2) \Rightarrow \text{约束最小化 }(F3,F4) \Rightarrow \text{变分残量/切线 }(F18\text{--}F25)
\Rightarrow \text{离散有约束优化 }(F28\text{--}F49)
\Rightarrow \text{L-BFGS-B 高效求解 }(F50\text{--}F105).
\]

---

## 1. F1–F7：问题定义与约束结构

### F1（Eq.1）总势能泛函
- 公式（完整显示）：

\[
\Pi (\pmb {u},d) = \int_{\Omega}\psi (\pmb {e}(\pmb {u}),d)\mathrm{d}\Omega + g_{c}\Gamma_{l}(d) - \int_{\Omega}\pmb {b}\cdot \pmb {u}\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot \pmb {u}\mathrm{d}\Gamma , \quad (1)
\]

- 含义：总势能 = 弹性能 + 断裂表面能 - 外力做功。
- 推导：由 Griffith 断裂能思想和小应变弹性势能叠加得到；相场变量 \(d\) 进入 \(\psi(\epsilon,d)\) 与裂纹项。
- 联系：F2 给出其中 \(\Gamma_l\) 的具体正则化形式；F3/F4 以 F1 为目标函数建立约束最小化。

### F2（Eq.2）相场裂纹表面近似
- 公式（完整显示）：

\[
\Gamma_{l}(d) = \int_{\Omega}\gamma (d,\nabla d)\mathrm{d}\Omega = \int_{\Omega}\frac{1}{2l}\left(d^{2} + l^{2}\nabla d\cdot \nabla d\right)\mathrm{d}\Omega , \quad (2)
\]

- 含义：用体积分近似真实裂纹面积。
- 推导：Ambrosio–Tortorelli 正则化：\(\gamma=\frac1{2l}(d^2+l^2|\nabla d|^2)\)。
- 联系：代回 F1 得完整目标函数；后续 F106 直接复用该表达监测裂纹能。

### F3（Eq.3）时步最小化问题
- 公式（完整显示）：

\[
\left(\pmb{u}_{n + 1},d_{n + 1}\right) = \arg \min \Pi (\pmb {u},d) \quad (3)
\]

- 含义：每一加载步求 \((u_{n+1},d_{n+1})\) 使总势能最小。
- 推导：准静态假设下惯性忽略，增量平衡等价于势能最小。
- 联系：约束条件由 F4 给出；离散后得到 F28/F29。

### F4（Eq.4）不可逆与有界性
- 公式（完整显示）：

\[
0\leq d_{n}\leq d_{n + 1}\leq 1. \quad (4)
\]

- 含义：\(d\) 不能下降，且始终在 [0,1]。
- 推导：热力学不可逆 + 相场变量物理定义（0 完好，1 完全损伤）。
- 联系：离散后成为盒约束 F6/F7/F26/F29，是 L-BFGS-B 的核心前提。

### F5 投影算子原型
- 公式（完整显示）：

\[
\mathrm{Proj}_C(\mathbf{x}_k = \mathbf{a}_k\nabla f(\mathbf{x}_k)).
\]

- 含义：沿下降方向后投影回可行域 \(C\)。
- 推导：有约束一阶法通式 \(x^+=\mathrm{Proj}_C(x-a\nabla f)\)（原文式子写法略简写）。
- 联系：当 \(C\) 是盒约束时，投影退化为分量截断（F32）。

### F6 盒约束分量形式
- 公式（完整显示）：

\[
\mathrm{lb}_i\leq x_i\leq \mathrm{ub}_i,
\]

- 含义：每个自由度独立上下界。
- 推导：由约束集合 \(C=\prod_i[lb_i,ub_i]\) 得分量不等式。
- 联系：F32 的投影是对 F6 逐分量应用。

### F7 相场离散自由度的具体上下界
- 公式（完整显示）：

\[
\mathrm{lb}_i = d_i^{(n)}\leq d_i^{(n + 1)}\leq 1 = \mathrm{ub}_i,
\]

- 含义：\(lb_i=d_i^{(n)},\ ub_i=1\)。
- 推导：由 F4 在离散节点 i 上直接得到。
- 联系：把物理不可逆约束转成标准优化盒约束（连接物理与算法的关键桥梁）。

---

## 2. F8–F19：本构分裂、应力与弱式残量

### F8（Eq.5）能量正负分裂退化
- 公式（完整显示）：

\[
\psi(\pmb{\epsilon}, d) = [g(d) + k] \psi^{+}(\pmb{\epsilon}) + \psi^{-}(\pmb{\epsilon}), \quad (5)
\]

- 含义：只退化拉伸正能 \(\psi^+\)，压缩 \(\psi^-\) 不退化。
- 推导：脆性断裂主要受拉开裂，采用 \([g(d)+k]\psi^+ + \psi^-\)。
- 联系：F9 给出退化函数；F14/F25 中出现 \(g'(d),g''(d)\)。

### F9（Eq.6）退化函数
- 公式（完整显示）：

\[
g(d) = (1 - d)^{2}. \quad (6)
\]

- 含义：\(g(d)=(1-d)^2\) 单调下降且 \(g(0)=1,g(1)=0\)。
- 推导：常用二次退化，保证连续可导。
- 联系：导数 \(g'=-2(1-d), g''=2\) 在 F18/F25 进入残量与切线。

### F10 正负括号与 Heaviside
- 公式（完整显示）：

\[
\langle x\rangle_{+} = \frac{1}{2}(x + |x|), \quad \langle x\rangle_{-} = \frac{1}{2}(x - |x|), \quad H(x) = \begin{cases} 1 & x \geq 0, \\ 0 & x < 0. \end{cases}
\]

- 含义：定义张压分解与分段导数符号函数。
- 推导：\(x=\langle x\rangle_+ + \langle x\rangle_-\)，\(H\) 用于不可导点的分段表达。
- 联系：F13–F16 都依赖这些算子。

### F11 应变谱分解
- 公式（完整显示）：

\[
\pmb{\epsilon} = \sum_{\alpha} \epsilon_{\alpha} \mathbf{M}_{\alpha}, \quad \mathbf{M}_{\alpha} = \pmb{n}_{\alpha} \otimes \pmb{n}_{\alpha},
\]

- 含义：按主应变基分解 \(\epsilon\)。
- 推导：对称二阶张量可正交特征分解。
- 联系：F12/F13 通过特征值正负部分构造能量分裂。

### F12 正负应变张量
- 公式（完整显示）：

\[
\pmb{\epsilon}^{+} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{+} \mathbf{M}_{\alpha}, \quad \pmb{\epsilon}^{-} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{-} \mathbf{M}_{\alpha}.
\]

- 含义：仅对特征值做 \(\langle\cdot\rangle_\pm\) 后重组。
- 推导：由 F11 分量替换得到。
- 联系：直接代入 F13（正负应变能）与 F15（正负应力）。

### F13 正负应变能密度
- 公式（完整显示）：

\[
\psi^{+}(\pmb{\epsilon}) = \frac{1}{2}\lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{+}^{2} + \mu \pmb{\epsilon}^{+} : \pmb{\epsilon}^{+}, \quad \psi^{-}(\pmb{\epsilon}) = \frac{1}{2} \lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{-}^{2} + \mu \pmb{\epsilon}^{-} : \pmb{\epsilon}^{-},
\]

- 含义：拉伸/压缩能分别按 Lamé 常数构造。
- 推导：把 \(\epsilon^\pm\)、\(\langle\mathrm{tr}\epsilon\rangle_\pm\) 代入线弹性二次型。
- 联系：F14 对 F13 对 \(\epsilon\) 求导得应力。

### F14 总应力表达
- 公式（完整显示）：

\[
\pmb{\sigma} = \frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial \pmb{\epsilon}} = [g(d) + k]\frac{\partial\psi^{+}(\pmb{\epsilon})}{\partial \pmb{\epsilon}} + \frac{\partial\psi^{-}(\pmb{\epsilon})}{\partial \pmb{\epsilon}} = [g(d) + k]\pmb{\sigma}^{+} + \pmb{\sigma}^{-},
\]

- 含义：\(\sigma=\partial\psi/\partial\epsilon=[g(d)+k]\sigma^+ + \sigma^-\)。
- 推导：链式法则 + F8 的线性组合结构。
- 联系：F15 是 \(\sigma^\pm\) 显式式；F18/F23/F25 用于残量和切线。

### F15 正负应力
- 公式（完整显示）：

\[
\pmb{\sigma}^{+} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{+}\mathbf{I} + 2\mu \pmb{\epsilon}^{+}, \quad \pmb{\sigma}^{-} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{-}\mathbf{I} + 2\mu \pmb{\epsilon}^{-}.
\]

- 含义：\(\sigma^\pm = \lambda\langle\mathrm{tr}\epsilon\rangle_\pm I +2\mu\epsilon^\pm\)。
- 推导：对 F13 分别对 \(\epsilon\) 求导。
- 联系：代回 F14 得总应力；其导数进入 F16。

### F16 切线刚度（四阶张量）
- 公式（完整显示）：

\[
\frac{\partial \pmb{\sigma}}{\partial \pmb{\epsilon}} = [g(d) + k]\frac{\partial \pmb{\sigma}^{+}}{\partial \pmb{\epsilon}} + \frac{\partial \pmb{\sigma}^{-}}{\partial \pmb{\epsilon}} = [g(d) + k]\left[\lambda H(\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes \mathbf{I} + 2\mu \mathbb{P}^{+}\right] + \left[\lambda H(-\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes \mathbf{I} + 2\mu \mathbb{P}^{-}\right].
\]

- 含义：\(\partial\sigma/\partial\epsilon\) 的正负分裂表达。
- 推导：对 F14 再求导，\(\partial\langle x\rangle_+/\partial x=H(x)\)，得到 \(I\otimes I\) 与投影算子项。
- 联系：F25 的 \(K_{uu}\) 直接使用该切线。

### F17 投影算子定义
- 公式（完整显示）：

\[
\mathbb{P}^{+} = \frac{\partial \pmb{\epsilon}^{+}}{\partial \pmb{\epsilon}}, \quad \mathbb{P}^{-} = \frac{\partial \pmb{\epsilon}^{-}}{\partial \pmb{\epsilon}},
\]

- 含义：\(\mathbb P^\pm=\partial\epsilon^\pm/\partial\epsilon\)。
- 推导：由谱分解的 Fréchet 导数定义。
- 联系：为 F16 提供紧凑写法。

### F18（Eq.7）总势能一阶变分
- 公式（完整显示）：

\[
\begin{array}{rl} 
\delta \Pi (\pmb {u},d) &= D_{(\delta \pmb {u},\delta d)}\Pi (\pmb {u},d) = \left.\frac{\mathrm{d}}{\mathrm{d}\epsilon}\right|_{\epsilon = 0}\Pi (\pmb {u} + \epsilon \delta \pmb {u},d + \epsilon \delta d) \\ 
&= \int_{\Omega}\left(\frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial \pmb{\epsilon}}:\pmb{\epsilon}(\delta \pmb {u}) + \frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial d}\delta d\right)\mathrm{d}\Omega \\ 
&\quad +\int_{\Omega}\frac{g_{c}}{l}\left(d\delta d + l^{2}\nabla d\cdot \nabla \delta d\right)\mathrm{d}\Omega -\int_{\Omega}\pmb {b}\cdot \delta \pmb {u}\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot \delta \pmb {u}\mathrm{d}\Gamma \\ 
&= (\nabla^{(s)}\delta \pmb {u},\pmb {\sigma}) - (\delta \pmb {u},\pmb {b}) - (\delta \pmb {u},\pmb{t})_{\Gamma_{t}} + (\delta d,\frac{g_{c}}{l} d) + (\nabla \delta d,g_{c}l\nabla d) + (\delta d,g^{\prime}(d)\psi^{+}(\pmb{\epsilon})). 
\end{array} \quad (7)
\]

- 含义：得到弱式平衡方程的统一表达。
- 推导：对 F1 中各项分别 Gateaux 导数：
  - 弹性能项给 \(\sigma:\epsilon(\delta u)\) 与 \(\partial\psi/\partial d\,\delta d\)
  - 裂纹正则项给 \((g_c/l)d\delta d + g_cl\nabla d\cdot\nabla\delta d\)
  - 外载项给负号做功项。
- 联系：将 \(\delta u\)、\(\delta d\) 系数分别置零得到 F19。

### F19 位移与相场残量方程
- 公式（完整显示）：

\[
\left\{ \begin{array}{ll} r_{\pmb{u}}(\pmb{u},d) = (\nabla^{(s)}\delta \pmb{u},\pmb{\sigma}) - (\delta \pmb{u},\pmb{b}) - (\delta \pmb{u},\pmb{t})_{\Gamma_{t}} = 0, \\ r_{d}(\pmb{u},d) = (\delta d,\frac{g_{c}}{l} d) + (\nabla \delta d,g_{c}l\nabla d) + (\delta d,g^{\prime}(d)\psi^{+}(\pmb{\epsilon})) = 0, \end{array} \right.
\]

- 含义：弱式双场耦合方程组 \(r_u=0,r_d=0\)。
- 推导：由 F18 对任意试函数成立的基本引理得到。
- 联系：F23 是其有限元离散版本；F108 用于收敛判据。

---

## 3. F20–F30：有限元离散与约束化

### F20 离散插值
- 公式（完整显示）：

\[
\pmb {u} = \pmb{N}_{u_{A}}\pmb{u}_{A}, \quad \text{and} \quad d = N_{d_{A}}d_{A}.
\]

- 含义：\(u,d\) 用形函数展开。
- 推导：标准有限元近似。
- 联系：代入 F1 得 F22。

### F21 变分插值
- 公式（完整显示）：

\[
\delta \pmb {u} = \pmb{N}_{u_{A}}\delta \pmb{u}_{A}, \quad \text{and} \quad \delta d = N_{d_{A}}\delta d_{A},
\]

- 含义：试函数也在同一近似空间展开。
- 推导：Galerkin 方法。
- 联系：用于从 F18 得离散残量 F23。

### F22（Eq.8）离散总势能
- 公式（完整显示）：

\[
\begin{array}{rl} 
\Pi (\pmb {u}_{A},d_{A}) = &\int_{\Omega}\psi \left(\pmb{\epsilon}(\pmb {N}_{u_{A}}\pmb {u}_{A}),N_{d_{A}}d_{A}\right)\mathrm{d}\Omega \\ 
&+\int_{\Omega}\frac{g_{c}}{2l}\left((N_{d_{A}}d_{A})^{2} + l^{2}(\nabla N_{d_{A}}d_{A})\cdot (\nabla N_{d_{A}}d_{A})\right)\mathrm{d}\Omega \\ 
&-\int_{\Omega}\pmb {b}\cdot (\pmb {N}_{u_{A}}\pmb {u}_{A})\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot (\pmb {N}_{u_{A}}\pmb {u}_{A})\mathrm{d}\Gamma, 
\end{array} \quad (8)
\]

- 含义：把连续泛函写成节点未知量函数。
- 推导：F20 代入 F1/F2，积分域不变，未知转为 \(u_A,d_A\)。
- 联系：对其求梯度与 Hessian 得 F23/F24。

### F23（Eq.9）离散梯度/残量
- 公式（完整显示）：

\[
\begin{array}{rl} 
\pmb{r} = \nabla \Pi = (r_{\pmb{u}},r_{d})^{\mathrm{T}},\\ 
r_{\pmb{u}_{A}} = \frac{\partial\Pi}{\partial \pmb{u}_{A}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},\pmb{\sigma}\right) - \left(\pmb{N}_{u_{A}},\pmb{b}\right) - \left(\pmb{N}_{u_{A}},\pmb{t}\right)_{\Gamma_{t}},\\ 
r_{d_{A}} = \frac{\partial\Pi}{\partial d_{A}} = \left(N_{d_{A}},\frac{g_{c}}{l} d + g^{\prime}(d)\psi^{+}\right) + (\nabla N_{d_{A}},g_{c}l\nabla d). 
\end{array} \quad (9)
\]

- 含义：\(r=\nabla\Pi=(r_u,r_d)^T\)。
- 推导：对 F22 分别对 \(u_A,d_A\) 偏导。
- 联系：F24/F25 是再求导；F31 中 \(r_k\) 用于二次模型。

### F24 离散 Hessian 分块
- 公式（完整显示）：

\[
\mathbf{K} = \nabla^{2}\Pi = \left[ \begin{array}{ll}\mathbf{K}_{uu} & \mathbf{K}_{ud}\\ \mathbf{K}_{du} & \mathbf{K}_{dd} \end{array} \right],
\]

- 含义：\(K=\nabla^2\Pi\) 的 2×2 块结构。
- 推导：对 F23 再偏导。
- 联系：F25 给各块表达；F27/F60 取其块对角近似。

### F25（Eq.10）各块切线刚度
- 公式（完整显示）：

\[
\begin{array}{rl} 
\mathbf{K}_{u_{A}u_{B}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},\frac{\partial \pmb{\sigma}}{\partial \pmb{\epsilon}}:\nabla^{(s)}\pmb{N}_{u_{B}}\right),\qquad \mathbf{K}_{u_{A}d_{B}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},g^{\prime}(d)\pmb{\sigma}^{+}N_{d_{B}}\right),\\ 
\mathbf{K}_{d_{A}u_{B}} = \left(N_{d_{A}},g^{\prime}(d)\pmb{\sigma}^{+}:\nabla^{(s)}\pmb{N}_{u_{B}}\right),\qquad \mathbf{K}_{d_{A}d_{B}} = \left(N_{d_{A}},\left(\frac{g_{c}}{l} +g^{\prime \prime}(d)\psi^{+}\right)N_{d_{B}}\right) + \left(\nabla N_{d_{A}},g_{c}l\nabla N_{d_{B}}\right). 
\end{array} \quad (10)
\]

- 含义：给出 \(K_{uu},K_{ud},K_{du},K_{dd}\) 显式积分形式。
- 推导：分别对残量项线性化，利用 F14/F16 与 \(g',g''\)。
- 联系：数值实现中的牛顿/拟牛顿基础；F27 仅保留块对角。

### F26 离散相场不可逆约束
- 公式（完整显示）：

\[
d_A^{(n)}\leq d_A\leq 1,
\]

- 含义：节点层面的 \(d_A^{(n)}\le d_A\le1\)。
- 推导：F4 离散化。
- 联系：和 F29 等价；用于定义 lb/ub。

### F27（Eq.11）块对角近似 Hessian
- 公式（完整显示）：

\[
\hat{\mathbf{K}} = \left[ \begin{array}{cc}\mathbf{K}_{uu} & \mathbf{0}\\ \mathbf{0} & \mathbf{K}_{dd} \end{array} \right], \quad (11)
\]

- 含义：忽略 \(K_{ud},K_{du}\) 的耦合项。
- 推导：为提高稳健性与效率采用近似 \(\hat K\)。
- 联系：作为 L-BFGS 初始矩阵 F60/F55 的 \(B_k^0\)。

### F28（Eq.12）离散约束最小化
- 公式（完整显示）：

\[
\left(\pmb{u}_A,d_A\right) = \arg \min \Pi (\pmb {u}_A,d_A) \quad (12)
\]

- 含义：离散未知量层面的目标函数最小化。
- 推导：F3 的离散版本。
- 联系：约束由 F29 给出，数值求解由第 3 节公式链完成。

### F29（Eq.13）离散盒约束
- 公式（完整显示）：

\[
d_A^{(n)}\leq d_A\leq 1 \quad (13)
\]

- 含义：和 F26 相同，强调优化角度。
- 推导：直接继承不可逆物理约束。
- 联系：构成 L-BFGS-B 的 B（bound）。

### F30（Eq.14）约束消元形式
- 公式（完整显示）：

\[
\left[\pmb {u}_A,d_A\right]^{\mathrm{T}} = \mathbf{C}\left[\pmb {u}_A,d_A\right]^{\mathrm{T}} + \pmb {k}. \quad (14)
\]

- 含义：\(x=Cx+k\) 统一处理 Dirichlet/固定分量。
- 推导：约束 DoF 的值写入 \(k\)，自由 DoF 由投影矩阵 \(C\) 保留。
- 联系：F46–F49 给线性系统层面的同类处理。

---

## 4. F31–F49：投影路径、活动集与约束线性代数

### F31（Eq.15）局部二次模型
- 公式（完整显示）：

\[
m_{k}(\pmb{x}) = \Pi_{k} + \pmb{r}_{k}^{\mathrm{T}}(\pmb{x} - \pmb{x}_{k}) + \frac{1}{2} (\pmb{x} - \pmb{x}_{k})^{\mathrm{T}}\mathbf{B}_{k}(\pmb{x} - \pmb{x}_{k}) \quad (15)
\]

- 含义：在 \(x_k\) 处二阶近似 \(\Pi\)。
- 推导：泰勒展开到二阶：常数 + 一阶梯度 + 二次项。
- 联系：F34/F83/F94 都是该模型在不同子空间的重写。

### F32（Eq.16）盒投影分量公式
- 公式（完整显示）：

\[
\mathrm{Proj}_c(x_i,\mathrm{lb}_i,\mathrm{ub}_i) = \left\{ \begin{array}{ll}\mathrm{lb}_i & \mathrm{if}\; x_i < \mathrm{lb}_i,\\ x_i & \mathrm{if}\; \mathrm{lb}_i \leq x_i \leq \mathrm{ub}_i,\\ \mathrm{ub}_i & \mathrm{if}\; x_i > \mathrm{ub}_i, \end{array} \right. \quad (16)
\]

- 含义：超下界取下界，超上界取上界，否则不变。
- 推导：\(\min_{y\in[lb_i,ub_i]}|y-x_i|\) 的闭式解。
- 联系：F33 逐向量应用此投影。

### F33（Eq.17）投影路径
- 公式（完整显示）：

\[
\pmb {x}(t) = \mathrm{Proj}_c(\pmb {x}_k - t\pmb {r}_k,\mathbf{lb},\mathbf{ub}),\quad t > 0, \quad (17)
\]

- 含义：沿负梯度 \(-r_k\) 走步长 \(t\) 后投影。
- 推导：约束最速下降路径定义。
- 联系：F61/F62 给出分量“撞界时间”与分段线性轨迹。

### F34（Eq.18）路径上一维模型
- 公式（完整显示）：

\[
p_k(t) = m_k(\pmb{x}(t)) = \Pi_k + \pmb{r}_k^\mathrm{T}(\pmb{x}(t) - \pmb{x}_k) + \frac{1}{2} (\pmb{x}(t) - \pmb{x}_k)^\mathrm{T}\mathbf{B}_k(\pmb{x}(t) - \pmb{x}_k). \quad (18)
\]

- 含义：把多维二次模型限制到 \(x(t)\) 上。
- 推导：把 F33 代入 F31。
- 联系：广义 Cauchy 点通过最小化该一维分段二次函数得到（见 F64–F66）。

### F35 活动集定义
- 公式（完整显示）：

\[
\mathcal{A}(\pmb {x}) = \{i:x_i = \mathrm{lb}_i\} \cup \{i:x_i = \mathrm{ub}_i\}.
\]

- 含义：在边界上的分量集合。
- 推导：KKT 中活跃不等式约束概念。
- 联系：F37/F42/F67/F76/F77/F80/F95 都围绕活动集更新。

### F36 子问题最优点定义
- 公式（完整显示）：

\[
\pmb{x}^* = \arg \min m_k(\pmb {x})
\]

- 含义：\(x^*=\arg\min m_k(x)\)。
- 推导：在当前迭代二次近似上求最优。
- 联系：实际用“两步法”：先 F33 得 \(x^c\)，再 F79–F85 子空间修正。

### F37 活动集固定后的可行域
- 公式（完整显示）：

\[
x_{i} = x_{i}^{c},\;\forall i\in \mathcal{A}_{k}(\pmb{x}^{c})\quad \text{and}\quad \mathrm{lb}_{i}\leq x_{i}\leq \mathrm{ub}_{i},\;\forall i\notin \mathcal{A}_{k}(\pmb{x}^{c}).
\]

- 含义：活动变量固定为 \(x^c\)，其余保持盒约束。
- 推导：active-set 思想：边界变量冻结，内点变量再优化。
- 联系：F79–F82 是该约束下的数学化。

### F38 搜索方向
- 公式（完整显示）：

\[
\pmb{p}_k = \pmb{x}^* - \pmb{x}_k.
\]

- 含义：\(p_k=x^*-x_k\)。
- 推导：从当前点指向子问题近似最优点。
- 联系：F39/F40/F41 线搜索沿该方向进行。

### F39 迭代更新
- 公式（完整显示）：

\[
\pmb{x}_{k + 1} = \pmb{x}_k + \alpha_k\pmb{p}_k.
\]

- 含义：\(x_{k+1}=x_k+\alpha_k p_k\)。
- 推导：一维步长策略。
- 联系：\(\alpha_k\) 由 Wolfe 条件 F40/F41 选取。

### F40 Armijo 条件
- 公式（完整显示）：

\[
\Pi_{k + 1} = \Pi (\pmb{x}_k + \alpha_k\pmb {p}_k) \leq \Pi_k + c_1\alpha_k\pmb {p}_k^\mathrm{T}\pmb{r}_k
\]

- 含义：保证足够下降。
- 推导：一阶下降上界。
- 联系：与 F41 共同形成强 Wolfe 条件。

### F41 曲率条件
- 公式（完整显示）：

\[
\left|\pmb{r}_{k + 1}^\mathrm{T}\pmb {p}_k\right| = \left|\pmb {r}(\pmb{x}_k + \alpha_k\pmb {p}_k)^\mathrm{T}\pmb {p}_k\right| \leq c_2\left|\pmb{r}_k^\mathrm{T}\pmb {p}_k\right|.
\]

- 含义：避免步长过小，保证曲率信息有效。
- 推导：限制新梯度在方向 \(p_k\) 上的投影。
- 联系：满足该条件有助于后续 BFGS 正定更新（关联 F52）。

### F42 活动集稳定判据
- 公式（完整显示）：

\[
\mathcal{A}_{k + 1}(\pmb{x}^c) = \mathcal{A}_k(\pmb{x}^c).
\]

- 含义：若活动集不再变化可进入下一阶段/终止判定。
- 推导：active-set 算法常见收敛信号。
- 联系：与 F43/F44 组成终止标准。

### F43 投影梯度范数判据
- 公式（完整显示）：

\[
\| \mathrm{Proj}_C(\mathbf{x}_{k + 1} - \mathbf{r}_{k + 1},\mathbf{lb},\mathbf{ub}) - \mathbf{x}_{k + 1}\|_{2} < \mathrm{tol}.
\]

- 含义：\(\|\mathrm{Proj}(x-r)-x\|\) 小表示近似满足约束一阶最优条件。
- 推导：盒约束下投影梯度是 KKT 残量。
- 联系：F45 给其分量解释。

### F44 变量增量判据
- 公式（完整显示）：

\[
\| \Delta \mathbf{x}_{k + 1}\|_{2} = \| \mathbf{x}_{k + 1} - \mathbf{x}_{k}\|_{2} < \mathrm{tol}.
\]

- 含义：步长足够小则可停止。
- 推导：迭代稳定性判据。
- 联系：与 F43 互补：一个看最优性，一个看变化量。

### F45（Eq.19）投影梯度分量展开
- 公式（完整显示）：

\[
\mathrm{Proj}_C(x_i - r_i,\mathrm{lb}_i,\mathrm{ub}_i) - x_i = \left\{ \begin{array}{ll}\mathrm{lb}_i - x_i & \mathrm{if}\; x_i - r_i < \mathrm{lb}_i,\\ - r_i & \mathrm{if}\; \mathrm{lb}_i \leq x_i - r_i \leq \mathrm{ub}_i,\\ \mathrm{ub}_i - x_i & \mathrm{if}\; x_i - r_i > \mathrm{ub}_i. \end{array} \right. \quad (19)
\]

- 含义：三种状态（撞下界/自由/撞上界）对应不同残量。
- 推导：把 F32 代入 F43 的分量表达。
- 联系：清晰揭示活动集判定机制（自由变量时残量就是 \(-r_i\)）。

### F46 线性约束编码
- 公式（完整显示）：

\[
\mathbf{x} = \mathbf{C}\mathbf{x} + \mathbf{k}.
\]

- 含义：\(x=Cx+k\) 的线性系统级重述。
- 推导：固定自由度重排写成投影矩阵形式。
- 联系：F47/F48 给求解形式。

### F47 约束后线性系统
- 公式（完整显示）：

\[
\left(\mathbf{C}^{\mathrm{T}}\mathbf{A}\mathbf{C} + \mathbf{I}_{d_{c}}\right)\hat{\mathbf{x}} = \mathbf{C}^{\mathrm{T}}(\mathbf{b} - \mathbf{A}\mathbf{k})
\]

- 含义：在受约束空间求 \(\hat x\) 的方程。
- 推导：将 F46 代入 \(Ax=b\)，左乘 \(C^T\) 并加入约束单位对角稳定项。
- 联系：F49 定义该单位对角矩阵含义。

### F48 回代恢复原变量
- 公式（完整显示）：

\[
\mathbf{x} = \mathbf{C}\hat{\mathbf{x}} + \mathbf{k}.
\]

- 含义：先求 \(\hat x\) 再映射回 \(x\)。
- 推导：直接由 F46 改写。
- 联系：与 F47 成对出现。

### F49 约束指示矩阵
- 公式（完整显示）：

\[
(\mathbf{I}_{d_c})_{ii} = \begin{cases} 1 & \text{if DoF } i \text{ is constrained}, \\ 0 & \text{otherwise}, \end{cases}
\]

- 含义：\(I_{d_c}\) 在受约束自由度对角为 1。
- 推导：用于在约束行上强制可逆与稳定。
- 联系：F47 中不可缺少。

---

## 5. F50–F60：BFGS 与紧凑 L-BFGS 表示

### F50 位移对与梯度差
- 公式（完整显示）：

\[
\mathbf{s}_k = \mathbf{x}_{k + 1} - \mathbf{x}_k,\quad \mathbf{y}_k = \mathbf{r}_{k + 1} - \mathbf{r}_k,
\]

- 含义：\(s_k=x_{k+1}-x_k,\ y_k=r_{k+1}-r_k\)。
- 推导：拟牛顿标准历史对。
- 联系：F51 更新公式直接使用。

### F51 BFGS 更新
- 公式（完整显示）：

\[
\mathbf{B}_{k + 1} = \mathbf{B}_k - \frac{\mathbf{B}_k\mathbf{s}_k\mathbf{s}_k^{\mathrm{T}}\mathbf{B}_k}{\mathbf{s}_k^{\mathrm{T}}\mathbf{B}_k\mathbf{s}_k} + \frac{\mathbf{y}_k\mathbf{y}_k^{\mathrm{T}}}{\mathbf{y}_k^{\mathrm{T}}\mathbf{s}_k}.
\]

- 含义：用秩二更新改进 \(B_k\)。
- 推导：满足割线条件 \(B_{k+1}s_k=y_k\) 且最小改变量的经典结果。
- 联系：F52 是其正定性条件。

### F52 曲率条件
- 公式（完整显示）：

\[
\mathbf{s}_k^{\mathrm{T}}\mathbf{y}_k > 0,
\]

- 含义：\(s_k^Ty_k>0\) 保证更新后矩阵保持正定。
- 推导：BFGS 理论结论。
- 联系：线搜索满足 Wolfe（F40/F41）通常可促成该条件。

### F53 历史步矩阵 \(S_k\)
- 公式（完整显示）：

\[
\mathbf{S}_k = [\mathbf{s}_{k-m} \cdots \mathbf{s}_{k-1}]
\]

- 含义：存储最近 m 个 \(s\) 向量。
- 推导：limited-memory 思想。
- 联系：与 F54 一起构成紧凑表示输入。

### F54 历史梯度差矩阵 \(Y_k\)
- 公式（完整显示）：

\[
\mathbf{Y}_k = [\mathbf{y}_{k-m} \cdots \mathbf{y}_{k-1}].
\]

- 含义：存储最近 m 个 \(y\) 向量。
- 推导：同上。
- 联系：F56/F57 使用。

### F55（Eq.20）紧凑表示
- 公式（完整显示）：

\[
\mathbf{B}_k = \mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k, \quad (20)
\]

- 含义：\(B_k=B_k^0-W_kM_kW_k^T\)。
- 推导：将多次 BFGS 秩二更新合并成低秩修正形式。
- 联系：F72/F86 在后续子空间最小化再次使用该结构。

### F56 矩阵 \(W_k\)
- 公式（完整显示）：

\[
\mathbf{W}_k = [\mathbf{Y}_k \quad \mathbf{B}^0_k \mathbf{S}_k] \in \mathbb{R}^{n\times 2m}
\]

- 含义：由 \(Y_k\) 和 \(B_k^0S_k\) 拼接。
- 推导：紧凑表示的基矩阵定义。
- 联系：F55/F57/F90。

### F57 矩阵 \(M_k\)
- 公式（完整显示）：

\[
\mathbf{M}_k = \begin{bmatrix} -\mathbf{D}_k & \mathbf{L}^{\mathrm{T}}_k \\ \mathbf{L}_k & \mathbf{S}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{S}_k \end{bmatrix}^{-1} \in \mathbb{R}^{2m\times 2m}.
\]

- 含义：由 \(D_k,L_k,S_k^TB_k^0S_k\) 组成的小规模矩阵逆。
- 推导：由 BFGS 累积更新代数整理得到。
- 联系：维度仅 \(2m\times2m\)，是高效关键。

### F58 对角矩阵 \(D_k\)
- 公式（完整显示）：

\[
\mathbf{D}_k = \mathrm{diag}\{\mathbf{s}^{\mathrm{T}}_{k-m}\mathbf{y}_{k-m}, \ldots, \mathbf{s}^{\mathrm{T}}_{k-1}\mathbf{y}_{k-1}\} \in \mathbb{R}^{m\times m}
\]

- 含义：对角元为 \(s_i^Ty_i\)。
- 推导：来自 BFGS 更新中的曲率标量。
- 联系：进入 F57。

### F59 下三角矩阵 \(L_k\)
- 公式（完整显示）：

\[
\mathbf{L}_k \in \mathbb{R}^{m\times m}, \quad (\mathbf{L}_k)_{i,j} = \begin{cases} \mathbf{s}^{\mathrm{T}}_{i+k-m-1}\mathbf{y}_{j+k-m-1} & \text{if } i > j, \\ 0 & \text{if } i \leq j. \end{cases}
\]

- 含义：存储历史向量的交叉内积（严格下三角）。
- 推导：由时间顺序配对的 \(s^Ty\) 填充。
- 联系：进入 F57。

### F60（Eq.21）初始近似矩阵
- 公式（完整显示）：

\[
\mathbf{B}^0_k = \hat{\mathbf{K}}^{(k)} = \begin{bmatrix} \mathbf{K}_{uu} & \mathbf{0} \\ \mathbf{0} & \mathbf{K}_{dd} \end{bmatrix}^{(k)}. \quad (21)
\]

- 含义：取块对角 \(\hat K\) 作为 \(B_k^0\)。
- 推导：由 F27 在当前迭代 \(k\) 取值。
- 联系：连接有限元切线与 L-BFGS 初始度量。

---

## 6. F61–F78：广义 Cauchy 点与活动集更新推导

### F61 撞界时间
- 公式（完整显示）：

\[
t_i = \begin{cases} (x^0_i - \mathrm{ub}_i)/r_i & \text{if } r_i < 0, \\ (x^0_i - \mathrm{lb}_i)/r_i & \text{if } r_i > 0, \\ +\infty & \text{if } r_i = 0. \end{cases}
\]

- 含义：第 i 分量沿 \(-r\) 前进何时触碰边界。
- 推导：解 \(x_i^0-tr_i=ub_i\) 或 \(lb_i\) 得分段公式。
- 联系：F62/F67/F76/F77 基于该时间构造活动集。

### F62 分量轨迹
- 公式（完整显示）：

\[
x_i(t) = x_i^{0} - \min \{t,t_i\} r_i.
\]

- 含义：\(x_i(t)=x_i^0-\min(t,t_i)r_i\)，触界后冻结。
- 推导：把 F61 引入投影路径可得。
- 联系：组成整体路径 F33 的分量表达。

### F63 分段线性段表达
- 公式（完整显示）：

\[
\mathbf{x}(t) = \mathbf{x}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)},
\]

- 含义：在第 \(j-1\) 段内 \(x(t)=x^{(j-1)}+\Delta t d^{(j-1)}\)。
- 推导：活动集固定时方向不变，故线性。
- 联系：代入 F31 得段内二次函数 F64。

### F64（Eq.22）段内一维二次模型
- 公式（完整显示）：

\[
\begin{array}{rl} 
p(t) &= m(\mathbf{x}(t)) = \Pi + \mathbf{r}^{\mathrm{T}}(\mathbf{x}(t) - \mathbf{x}^0) + \frac{1}{2} (\mathbf{x}(t) - \mathbf{x}^0)^{\mathrm{T}}\mathbf{B}(\mathbf{x}(t) - \mathbf{x}^0) \\
&= \Pi + \mathbf{r}^{\mathrm{T}}(\mathbf{z}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)}) + \frac{1}{2} (\mathbf{z}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)})^{\mathrm{T}}\mathbf{B}(\mathbf{z}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)}) \\
&= \left(\Pi + \mathbf{r}^{\mathrm{T}}\mathbf{z}^{(j-1)} + \frac{1}{2} \mathbf{z}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j-1)}\right) \\
&\quad + \left(\mathbf{r}^{\mathrm{T}}\mathbf{d}^{(j-1)} + \mathbf{d}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j-1)}\right)\Delta t + \frac{1}{2}\left(\mathbf{d}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j-1)}\right)\Delta t^{2} \\
&= f_{j-1} + f_{j-1}^{\prime}\Delta t + \frac{1}{2} f_{j-1}^{\prime \prime}\Delta t^{2} = \hat{p} (\Delta t), 
\end{array} \quad (22)
\]

- 含义：\(p(t)\) 在每段可写为 \(f+f'\Delta t+\frac12f''\Delta t^2\)。
- 推导：将 F63 代入二次模型并按 \(\Delta t\) 收集同类项。
- 联系：F65/F66 由该二次式直接得到。

### F65 二阶系数正性
- 公式（完整显示）：

\[
f_{j-1}^{\prime \prime} = \mathbf{d}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j-1)} > 0.
\]

- 含义：\(f''=d^TBd>0\) 表示该段凸。
- 推导：若 \(B\) 在子空间正定则成立。
- 联系：保证 F66 的极小点唯一。

### F66 段内极小步长
- 公式（完整显示）：

\[
\hat{p}'(\Delta t) = 0 \Rightarrow \Delta t^* = -f_{j-1}' / f_{j-1}^{\prime \prime}.
\]

- 含义：\(\Delta t^*=-f'/f''\)。
- 推导：令 F64 导数为零。
- 联系：与区间端点比较决定是否在本段取 Cauchy 点。

### F67 当前段活动集
- 公式（完整显示）：

\[
\mathcal{A}(\mathbf{x}^c) = \{i : t_i \leq t^{(j-1)}, i = 1, \ldots, n\}.
\]

- 含义：已满足 \(t_i\le t^{(j-1)}\) 的分量属于活动集。
- 推导：时间阈值定义。
- 联系：下一段更新见 F68/F76。

### F68 段推进
- 公式（完整显示）：

\[
\Delta t^{(j-1)} = t^{(j)} - t^{(j-1)}, \quad \mathbf{x}^{(j)} = \mathbf{x}^{(j-1)} + \Delta t^{(j-1)}\mathbf{d}^{(j-1)}.
\]

- 含义：更新时间与位置。
- 推导：离散步进公式。
- 联系：供 F70/F71 的递推量更新。

### F69 方向更新
- 公式（完整显示）：

\[
\mathbf{d}^{(j)} = \mathbf{d}^{(j-1)} + r_b \mathbf{e}_b,
\]

- 含义：当新分量 b 撞界时，从方向中去除该分量（通过 \(+r_b e_b\)）。
- 推导：使边界分量后续变化为 0。
- 联系：带来导数递推 F70/F71。

### F70（Eq.23）一阶导递推
- 公式（完整显示）：

\[
f'_j = \mathbf{r}^{\mathrm{T}}\mathbf{d}^{(j)} + \mathbf{d}^{(j)}{}^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j)} = f'_{j-1} + \Delta t^{(j-1)} f''_{j-1} + r^2_b + r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j)} \quad (23)
\]

- 含义：跨段更新 \(f'_j\) 的显式递推式。
- 推导：将 F68/F69 代入 \(f'=r^Td+d^TBz\) 展开整理。
- 联系：用于快速扫描分段最小点。

### F71（Eq.24）二阶导递推
- 公式（完整显示）：

\[
f''_j = \mathbf{d}^{(j)}{}^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j)} = f''_{j-1} + 2r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j-1)} + r^2_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{e}_b. \quad (24)
\]

- 含义：跨段更新 \(f''_j\)。
- 推导：对 \(d^{(j)}=d^{(j-1)}+r_be_b\) 代入 \(d^TBd\) 展开。
- 联系：与 F70 成对，避免每段全量重算。

### F72 紧凑表示重申
- 公式（完整显示）：

\[
\mathbf{B} = \mathbf{B}^0 - \mathbf{W}\mathbf{M}\mathbf{W}^{\mathrm{T}}.
\]

- 含义：\(B=B^0-WMW^T\) 用于把 F70/F71 的矩阵乘降维。
- 推导：直接用 F55 的当前迭代形式。
- 联系：引出 F73–F75 的低秩加速公式。

### F73 辅助向量递推
- 公式（完整显示）：

\[
\mathbf{w}_b = \mathbf{W}^{\mathrm{T}}\mathbf{e}_b, \quad \mathbf{p}^{(j)} = \mathbf{W}^{\mathrm{T}}\mathbf{d}^{(j)} = \mathbf{W}^{\mathrm{T}}(\mathbf{d}^{(j-1)} + r_b \mathbf{e}_b) = \mathbf{p}^{(j-1)} + r_b \mathbf{w}_b, \quad \mathbf{c}^{(j)} = \mathbf{W}^{\mathrm{T}}\mathbf{z}^{(j)} = \mathbf{W}^{\mathrm{T}}(\mathbf{z}^{(j-1)} + \Delta t^{(j-1)}\mathbf{d}^{(j-1)}) = \mathbf{c}^{(j-1)} + \Delta t^{(j-1)}\mathbf{p}^{(j-1)}.
\]

- 含义：定义 \(w_b,p^{(j)},c^{(j)}\) 便于低秩项复用。
- 推导：按 F69/F68 在线性映射 \(W^T\) 下的递推。
- 联系：代入 F74/F75。

### F74 低秩化的一阶导更新
- 公式（完整显示）：

\[
f'_j = f'_{j-1} + \Delta t^{(j-1)} f''_{j-1} + r^2_b + r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{z}^{(j)} - r_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{c}^{(j)}
\]

- 含义：把 F70 中 \(B\) 项改写为 \(B^0\) 与小矩阵 \(M\) 项。
- 推导：用 F72 展开 \(e_b^TBz\)。
- 联系：降低计算复杂度。

### F75 低秩化的二阶导更新
- 公式（完整显示）：

\[
f''_j = f''_{j-1} + 2r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{d}^{(j-1)} - 2r_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{p}^{(j-1)} + r^2_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{e}_b - r^2_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{w}_b.
\]

- 含义：把 F71 同样改写为 \(B^0\) + 低秩修正。
- 推导：用 F72 展开 \(e_b^TBd, e_b^TBe_b\)。
- 联系：与 F74 配合完成高效 Cauchy 点计算。

### F76 新时间下活动集
- 公式（完整显示）：

\[
\mathcal{A}(\mathbf{x}^c) = \{i : t_i \leq t^{(j)}, i = 1, \ldots, n\}.
\]

- 含义：推进到 \(t^{(j)}\) 后的活动集。
- 推导：F67 的时间更新版。
- 联系：与 F77 等价表述。

### F77 活动集等价定义
- 公式（完整显示）：

\[
\mathcal{A}(\mathbf{x}^c) = \{i : x^c_i = \mathrm{lb}_i\} \cup \{i : x^c_i = \mathrm{ub}_i\} = \{i : t_i \leq t^{(j)}, i = 1, \ldots, n\}.
\]

- 含义：边界集合与时间阈值集合等价。
- 推导：由 F62 可知 \(t_i\le t^{(j)}\iff x_i^c\) 在边界。
- 联系：连接几何定义（边界）与路径定义（撞界时间）。

### F78 活动分量边界取值
- 公式（完整显示）：

\[
x_i = x^c_i = \begin{cases} \mathrm{ub}_i & \text{if } r_i < 0 \\ \mathrm{lb}_i & \text{if } r_i > 0 \end{cases}, \quad \forall i \in \mathcal{A}(\mathbf{x}^c).
\]

- 含义：按梯度符号决定落在上界或下界。
- 推导：若 \(r_i<0\)，\(-r_i\) 方向增大 \(x_i\) 会先碰上界；反之碰下界。
- 联系：给出广义 Cauchy 点中活动变量的具体值。

---

## 7. F79–F105：子空间最小化（原始/对偶）与线搜索

### F79（Eq.25）子空间目标
- 公式（完整显示）：

\[
\mathbf{x}^* = \arg \min m_k(\mathbf{x}) \quad (25)
\]

- 含义：在活动集固定后继续最小化二次模型。
- 推导：承接 F37。
- 联系：约束细化见 F80/F81。

### F80 活动变量固定
- 公式（完整显示）：

\[
x_i = x^c_i, \quad \forall i \in \mathcal{A}(\mathbf{x}^c)
\]

- 含义：活动变量不再优化。
- 推导：active-set 原理。
- 联系：配合 F81 定义自由变量域。

### F81（Eq.26）自由变量盒约束
- 公式（完整显示）：

\[
\mathrm{lb}_i \leq x_i \leq \mathrm{ub}_i, \quad \forall i \notin \mathcal{A}(\mathbf{x}^c). \quad (26)
\]

- 含义：非活动变量仍需满足上下界。
- 推导：由原盒约束继承。
- 联系：F82 用零空间参数化进一步消元。

### F82 零空间参数化
- 公式（完整显示）：

\[
\mathbf{x} = \mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}},
\]

- 含义：\(x=x^c+Z_k\hat x\) 只在自由子空间移动。
- 推导：\(Z_k\) 选取自由变量基，活动变量对应行为 0。
- 联系：代入 F31 得 F83。

### F83 子空间模型展开
- 公式（完整显示）：

\[
\begin{array}{rl} 
m_k(\mathbf{x}) &= \Pi_k + \mathbf{r}_k^{\mathrm{T}}(\mathbf{x} - \mathbf{x}_k) + \frac{1}{2}(\mathbf{x} - \mathbf{x}_k)^{\mathrm{T}}\mathbf{B}_k(\mathbf{x} - \mathbf{x}_k) \\
\implies \hat{m}_k(\hat{\mathbf{x}}) &= \Pi_k + \mathbf{r}_k^{\mathrm{T}}(\mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}} - \mathbf{x}_k) + \frac{1}{2}(\mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}} - \mathbf{x}_k)^{\mathrm{T}}\mathbf{B}_k(\mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}} - \mathbf{x}_k) \\
&= \Pi_k + \mathbf{r}_k^{\mathrm{T}}(\mathbf{x}^c - \mathbf{x}_k) + \frac{1}{2} (\mathbf{x}^c - \mathbf{x}_k)^{\mathrm{T}}\mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k) \\
&\quad + [\mathbf{r}_k + \mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k)]^{\mathrm{T}}\mathbf{Z}_k \hat{\mathbf{x}} + \frac{1}{2} \hat{\mathbf{x}}^{\mathrm{T}}\mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k \hat{\mathbf{x}} \\
&= \Pi_k + \mathbf{r}_k^{\mathrm{T}}(\mathbf{x}^c - \mathbf{x}_k) + \frac{1}{2} (\mathbf{x}^c - \mathbf{x}_k)^{\mathrm{T}}\mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k) \\
&\quad + [\mathbf{Z}^{\mathrm{T}}_k [\mathbf{r}_k + \mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k)]]^{\mathrm{T}} \hat{\mathbf{x}} + \frac{1}{2} \hat{\mathbf{x}}^{\mathrm{T}}\mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k \hat{\mathbf{x}}.
\end{array}
\]

- 含义：把原模型重写为 \(\hat x\) 的二次函数。
- 推导：将 F82 代入 F31 并展开整理。
- 联系：其驻点方程是 F84。

### F84 子空间一阶最优条件
- 公式（完整显示）：

\[
\hat{m}'_k(\hat{\mathbf{x}}) = 0 \implies \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k \hat{\mathbf{x}} = -\mathbf{Z}^{\mathrm{T}}_k [\mathbf{r}_k + \mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k)].
\]

- 含义：\(Z^TBZ\hat x=-Z^T[r+B(x^c-x_k)]\)。
- 推导：对 F83 对 \(\hat x\) 求导为零。
- 联系：简写即 F85。

### F85（Eq.27）约化线性方程
- 公式（完整显示）：

\[
\hat{\mathbf{B}}_k \hat{\mathbf{x}} = -\hat{\mathbf{r}}_k. \quad (27)
\]

- 含义：\(\hat B_k\hat x=-\hat r_k\)。
- 推导：定义 \(\hat B_k,\hat r_k\) 后的紧凑形式。
- 联系：求解效率由 F87–F92 提升。

### F86 紧凑 Hessian 再代入
- 公式（完整显示）：

\[
\mathbf{B}_k = \mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k.
\]

- 含义：将 \(B_k=B_k^0-W_kM_kW_k^T\) 代入约化系统准备使用 SMW。
- 推导：直接代入。
- 联系：F87/F89。

### F87（Eq.28）约化矩阵分解
- 公式（完整显示）：

\[
\hat{\mathbf{B}}_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k = \mathbf{Z}^{\mathrm{T}}_k (\mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k) \mathbf{Z}_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{Z}_k - \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k. \quad (28)
\]

- 含义：\(\hat B_k=Z^TB_k^0Z - Z^TW_kM_kW_k^TZ\)。
- 推导：由 F86 左右乘 \(Z^T, Z\)。
- 联系：是低秩修正结构，可用 F89 求逆。

### F88 约化解形式
- 公式（完整显示）：

\[
\hat{\mathbf{x}} = -\hat{\mathbf{B}}^{-1}_k \hat{\mathbf{r}}_k,
\]

- 含义：\(\hat x=-\hat B_k^{-1}\hat r_k\)。
- 推导：线性方程显式解。
- 联系：关键在高效算 \(\hat B_k^{-1}\)（F89）。

### F89（Eq.29）SMW 型逆公式
- 公式（完整显示）：

\[
\hat{\mathbf{B}}^{-1}_k = (\hat{\mathbf{B}}^0_k)^{-1} + (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \left[\mathbf{I} - \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \right]^{-1} \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1}, \quad (29)
\]

- 含义：把大矩阵逆转化为 \(2m\) 维小矩阵逆。
- 推导：对 F87 的“基矩阵 + 低秩修正”应用 Sherman-Morrison-Woodbury。
- 联系：附录 F109/F110 给出理论依据。

### F90 \(Z^TW\) 的具体拼接
- 公式（完整显示）：

\[
\mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k = \mathbf{Z}^{\mathrm{T}}_k [\mathbf{Y}_k \quad \mathbf{B}^0_k \mathbf{S}_k] = [\mathbf{Z}^{\mathrm{T}}_k \mathbf{y}_{k-m} \cdots \mathbf{Z}^{\mathrm{T}}_k \mathbf{y}_{k-1} \quad \mathbf{Z}^{\mathrm{T}}_k(\mathbf{B}^0_k \mathbf{s}_{k-m}) \cdots \mathbf{Z}^{\mathrm{T}}_k(\mathbf{B}^0_k \mathbf{s}_{k-1})].
\]

- 含义：明确需预计算的块向量。
- 推导：由 F56 的列拼接直接左乘 \(Z^T\)。
- 联系：用于 F89 中小矩阵构造。

### F91 小矩阵维度说明
- 公式（完整显示）：

\[
\mathbf{I} - \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \in \mathbb{R}^{2m\times 2m}
\]

- 含义：需要求逆的仅是 \(2m\times2m\) 矩阵。
- 推导：由 \(W\) 列数为 \(2m\) 得出。
- 联系：说明 L-BFGS-B 的内存/计算优势。

### F92 预条件子
- 公式（完整显示）：

\[
\mathbf{P}^{-1} = \mathrm{ILU}(\hat{\mathbf{B}}^0_k) = \mathrm{ILU}(\mathbf{Z}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{Z}_k).
\]

- 含义：用 \(\hat B_k^0\) 的 ILU 作为迭代解算器预条件。
- 推导：近似主导刚度块，易构造且有效。
- 联系：用于求解 F85（原始子空间法）。

### F93 增量变量替换
- 公式（完整显示）：

\[
\mathbf{x} = \mathbf{x}_k + \Delta \mathbf{x}_k.
\]

- 含义：\(x=x_k+\Delta x_k\)。
- 推导：标准增量化。
- 联系：导出对偶子问题 F94–F102。

### F94 增量下二次模型
- 公式（完整显示）：

\[
m_k(\mathbf{x}) = \Pi_k + \mathbf{r}_k^{\mathrm{T}} \Delta \mathbf{x}_k + \frac{1}{2} \Delta \mathbf{x}_k^{\mathrm{T}} \mathbf{B}_k \Delta \mathbf{x}_k.
\]

- 含义：把 F31 改写成 \(\Delta x\) 变量。
- 推导：代入 F93。
- 联系：F97 作为其优化目标。

### F95 活动约束（增量形式）
- 公式（完整显示）：

\[
x_i = x^c_i, \quad \forall i \in \mathcal{A}_k(\mathbf{x}^c),
\]

- 含义：仍固定活动变量到 Cauchy 点。
- 推导：F80 的等价写法。
- 联系：用线性约束矩阵 \(Q_k\) 写成 F96/F98。

### F96 线性等式约束
- 公式（完整显示）：

\[
\mathbf{Q}^{\mathrm{T}}_k \Delta \mathbf{x}_k = \mathbf{Q}^{\mathrm{T}}_k (\mathbf{x} - \mathbf{x}_k) = \mathbf{Q}^{\mathrm{T}}_k (\mathbf{x}^c - \mathbf{x}_k).
\]

- 含义：\(Q_k^T\Delta x_k=Q_k^T(x^c-x_k)\)。
- 推导：把 F95 写成矩阵形式。
- 联系：与 F97 组成等式约束二次规划。

### F97 对偶法目标函数
- 公式（完整显示）：

\[
\min \left\{ \Pi_k + \mathbf{r}_k^{\mathrm{T}} \Delta \mathbf{x}_k + \frac{1}{2} \Delta \mathbf{x}_k^{\mathrm{T}} \mathbf{B}_k \Delta \mathbf{x}_k \right\}
\]

- 含义：最小化增量二次模型。
- 推导：即 F94 的目标。
- 联系：约束由 F98/F99 给出。

### F98 对偶法等式约束重述
- 公式（完整显示）：

\[
\mathbf{Q}_k^{\mathrm{T}} \Delta \mathbf{x}_k = \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k)
\]

- 含义：与 F96 同义（符号一致化）。
- 推导：同上。
- 联系：用于构造 KKT 系统 F100。

### F99 对偶法盒约束
- 公式（完整显示）：

\[
\mathbf{lb} - \mathbf{x}_k \leq \Delta \mathbf{x}_k \leq \mathbf{ub} - \mathbf{x}_k.
\]

- 含义：增量变量的上下界。
- 推导：由 \(lb\le x_k+\Delta x\le ub\) 变形。
- 联系：与活动集策略共同保证可行性。

### F100（Eq.30）KKT 线性系统
- 公式（完整显示）：

\[
\begin{bmatrix} \mathbf{B}_k & \mathbf{Q}_k \\ \mathbf{Q}_k^{\mathrm{T}} & \mathbf{0} \end{bmatrix} \begin{bmatrix} \Delta \mathbf{x}_k \\ \lambda \end{bmatrix} = \begin{bmatrix} -\mathbf{r}_k \\ \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k) \end{bmatrix}. \quad (30)
\]

- 含义：增量与拉格朗日乘子的联立方程。
- 推导：对 F97+F98 构造拉格朗日函数并求驻点。
- 联系：消元得到 F101/F102。

### F101（Eq.31）乘子方程
- 公式（完整显示）：

\[
(\mathbf{Q}_k^{\mathrm{T}} \mathbf{B}_k^{-1} \mathbf{Q}_k) \lambda = -\mathbf{Q}_k^{\mathrm{T}} \mathbf{B}_k^{-1} \mathbf{r}_k - \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k). \quad (31)
\]

- 含义：先解 \(\lambda\) 的 Schur 补系统。
- 推导：从 F100 消去 \(\Delta x\) 得到。
- 联系：再由 F102 回算 \(\Delta x\)。

### F102（Eq.32）增量回代方程
- 公式（完整显示）：

\[
\mathbf{B}_k \Delta \mathbf{x}_k = -(\mathbf{r}_k + \mathbf{Q}_k \lambda). \quad (32)
\]

- 含义：\(B_k\Delta x_k=-(r_k+Q_k\lambda)\)。
- 推导：F100 第一行或 F101 回代。
- 联系：完成对偶子空间最小化。

### F103 最终搜索方向
- 公式（完整显示）：

\[
\mathbf{p}_k = \mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}} - \mathbf{x}_k.
\]

- 含义：\(p_k=x^c+Z_k\hat x-x_k\)。
- 推导：Cauchy 点 + 子空间修正量。
- 联系：用于全局线搜索更新 F104。

### F104 主迭代更新
- 公式（完整显示）：

\[
\mathbf{x}_{k + 1} = \mathbf{x}_k + \alpha_k \mathbf{p}_k,
\]

- 含义：\(x_{k+1}=x_k+\alpha_k p_k\)。
- 推导：同 F39，但此处是完整 L-BFGS-B 方向。
- 联系：步长受可行性 F105 与 Wolfe 条件 F40/F41 约束。

### F105 步长可行域条件
- 公式（完整显示）：

\[
\mathbf{lb} \leq \mathbf{x}_k + \alpha_k \mathbf{p}_k \leq \mathbf{ub}.
\]

- 含义：更新后仍在盒约束内。
- 推导：直接将 F104 代入盒约束。
- 联系：与投影策略共同保证不可逆约束绝不被破坏。

---

## 8. F106–F110：算例监测、收敛判据与附录理论

### F106（Eq.33）裂纹能
- 公式（完整显示）：

\[
E_{\mathrm{crack}} = g_{\mathrm{c}}\Gamma_{l}(d) = g_{\mathrm{c}}\int_{\Omega}\gamma (d,\nabla d)\mathrm{d}\Omega = g_{\mathrm{c}}\int_{\Omega}\frac{1}{2l}\left(d^{2} + l^{2}\nabla d\cdot \nabla d\right)\mathrm{d}\Omega \quad (33)
\]

- 含义：\(E_{crack}=g_c\Gamma_l\)，用于后处理与能量演化分析。
- 推导：由 F1 中断裂能项直接抽取，并代入 F2。
- 联系：用于数值算例比较裂纹发展阶段。

### F107 扭转边界位移
- 公式（完整显示）：

\[
u_{y} = z \tan t, \quad u_{z} = -y \tan t,
\]

- 含义：给 3D 扭转例中的几何旋转位移边界。
- 推导：截面绕 x 轴小角度旋转的位移关系。
- 联系：为 F3 的每一载荷步提供外载/位移输入。

### F108 非线性收敛判据
- 公式（完整显示）：

\[
\| \pmb {r}_u\| _2 < \mathrm{tol}, \quad \| \pmb {r}_d\| _2 < \mathrm{tol}, \quad \| \Delta \pmb {u}\| _2 < \mathrm{tol}, \quad \| \Delta d\| _2 < \mathrm{tol},
\]

- 含义：位移残量、相场残量及两者增量均小于阈值。
- 推导：双场耦合问题常用“四条件”停止标准。
- 联系：和 F43/F44（优化视角）互相呼应。

### F109 附录低秩更新结构
- 公式（完整显示）：

\[
\hat{\mathbf{A}} = \mathbf{A} + \mathbf{U}\mathbf{V}^{\mathrm{T}},
\]

- 含义：\(\hat A=A+UV^T\) 的基本形式。
- 推导：为应用 SMW 做准备。
- 联系：对应正文 F89 的矩阵结构。

### F110 Sherman-Morrison-Woodbury 公式
- 公式（完整显示）：

\[
\hat{\mathbf{A}}^{-1} = \mathbf{A}^{-1} - \mathbf{A}^{-1}\mathbf{U}\left(\mathbf{I} + \mathbf{V}^{\mathrm{T}}\mathbf{A}^{-1}\mathbf{U}\right)^{-1}\mathbf{V}^{\mathrm{T}}\mathbf{A}^{-1}.
\]

- 含义：给出 \((A+UV^T)^{-1}\) 的闭式表达。
- 推导：由分块矩阵求逆或 Woodbury 恒等式可证。
- 联系：是 F89 高效求逆的理论基础。

---

## 9. 关键推导串联（把“逐式”连成“整体”）

### 串联 A：物理模型到离散优化
1. F1+F2 定义能量；
2. F3+F4 给约束最小化；
3. F18/F19 得弱式；
4. F20–F25 离散成 \(r,K\)；
5. F28/F29 转为离散盒约束优化。

### 串联 B：不可逆约束到投影算法
1. F4 \(\Rightarrow\) F7/F29（节点上下界）；
2. F32 定义投影，F33 给投影路径；
3. F61–F78 找广义 Cauchy 点与活动集；
4. F79–F105 在自由子空间继续最小化并更新。

### 串联 C：大规模效率来源
1. F51 需要 Hessian 近似更新；
2. F55–F57 把全量矩阵写成“块对角 + 低秩修正”；
3. F89+F110 把大系统逆转化为小系统逆；
4. F91 指出小系统仅 \(2m\times2m\)，因此适合大规模 FE。

---

## 10. 对“公式之间内在联系”的一句话总结

`explain.md` 的所有公式本质上是在完成同一件事：

> **把“相场断裂的不可逆物理约束（F4）”严格、稳定、可扩展地嵌入到“有限元离散后的非凸能量最小化（F22/F28）”中，并通过 L-BFGS-B 的投影 + 活动集 + 低秩更新机制（F31–F105）实现高效求解。**


---

## 11. 内在联系（按“完整公式”直接对照）

> 本节把“哪几个公式互相关联”用**完整公式并排**写出，便于不跳转地直接对照。

### 11.1 总势能 → 变分弱式（F1, F2, F18, F19）

- 起点（总势能）：

\[
\Pi (\pmb {u},d) = \int_{\Omega}\psi (\pmb {e}(\pmb {u}),d)\mathrm{d}\Omega + g_{c}\Gamma_{l}(d) - \int_{\Omega}\pmb {b}\cdot \pmb {u}\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot \pmb {u}\mathrm{d}\Gamma , \quad (1)
\]

- 裂纹项定义：

\[
\Gamma_{l}(d) = \int_{\Omega}\gamma (d,\nabla d)\mathrm{d}\Omega = \int_{\Omega}\frac{1}{2l}\left(d^{2} + l^{2}\nabla d\cdot \nabla d\right)\mathrm{d}\Omega , \quad (2)
\]

- 对 F1 作一阶变分并代入 F2，得到：

\[
\begin{array}{rl} 
\delta \Pi (\pmb {u},d) &= D_{(\delta \pmb {u},\delta d)}\Pi (\pmb {u},d) = \left.\frac{\mathrm{d}}{\mathrm{d}\epsilon}\right|_{\epsilon = 0}\Pi (\pmb {u} + \epsilon \delta \pmb {u},d + \epsilon \delta d) \\
&= \int_{\Omega}\left(\frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial \pmb{\epsilon}}:\pmb{\epsilon}(\delta \pmb {u}) + \frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial d}\delta d\right)\mathrm{d}\Omega \\
&\quad +\int_{\Omega}\frac{g_{c}}{l}\left(d\delta d + l^{2}\nabla d\cdot \nabla \delta d\right)\mathrm{d}\Omega -\int_{\Omega}\pmb {b}\cdot \delta \pmb {u}\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot \delta \pmb {u}\mathrm{d}\Gamma \\
&= (\nabla^{(s)}\delta \pmb {u},\pmb {\sigma}) - (\delta \pmb {u},\pmb {b}) - (\delta \pmb {u},\pmb{t})_{\Gamma_{t}} + (\delta d,\frac{g_{c}}{l} d) + (\nabla \delta d,g_{c}l\nabla d) + (\delta d,g^{\prime}(d)\psi^{+}(\pmb{\epsilon})). 
\end{array} \quad (7)
\]

- 对任意试函数令变分系数为零，得到残量方程：

\[
\left\{ \begin{array}{ll} r_{\pmb{u}}(\pmb{u},d) = (\nabla^{(s)}\delta \pmb{u},\pmb{\sigma}) - (\delta \pmb{u},\pmb{b}) - (\delta \pmb{u},\pmb{t})_{\Gamma_{t}} = 0, \\ r_{d}(\pmb{u},d) = (\delta d,\frac{g_{c}}{l} d) + (\nabla \delta d,g_{c}l\nabla d) + (\delta d,g^{\prime}(d)\psi^{+}(\pmb{\epsilon})) = 0, \end{array} \right.
\]

**内在联系**：F1/F2 是能量层定义，F18 是其微分，F19 是平衡条件（即变分驻值条件）的方程化。

### 11.2 不可逆约束 → 盒约束投影（F4, F7, F32, F33, F45）

- 不可逆物理约束：

\[
0\leq d_{n}\leq d_{n + 1}\leq 1. \quad (4)
\]

- 离散后每个自由度边界：

\[
\mathrm{lb}_i = d_i^{(n)}\leq d_i^{(n + 1)}\leq 1 = \mathrm{ub}_i,
\]

- 盒投影算子：

\[
\mathrm{Proj}_c(x_i,\mathrm{lb}_i,\mathrm{ub}_i) = \left\{ \begin{array}{ll}\mathrm{lb}_i & \mathrm{if}\; x_i < \mathrm{lb}_i,\\ x_i & \mathrm{if}\; \mathrm{lb}_i \leq x_i \leq \mathrm{ub}_i,\\ \mathrm{ub}_i & \mathrm{if}\; x_i > \mathrm{ub}_i, \end{array} \right. \quad (16)
\]

- 投影路径：

\[
\pmb {x}(t) = \mathrm{Proj}_c(\pmb {x}_k - t\pmb {r}_k,\mathbf{lb},\mathbf{ub}),\quad t > 0, \quad (17)
\]

- 投影梯度分量表达：

\[
\mathrm{Proj}_C(x_i - r_i,\mathrm{lb}_i,\mathrm{ub}_i) - x_i = \left\{ \begin{array}{ll}\mathrm{lb}_i - x_i & \mathrm{if}\; x_i - r_i < \mathrm{lb}_i,\\ - r_i & \mathrm{if}\; \mathrm{lb}_i \leq x_i - r_i \leq \mathrm{ub}_i,\\ \mathrm{ub}_i - x_i & \mathrm{if}\; x_i - r_i > \mathrm{ub}_i. \end{array} \right. \quad (19)
\]

**内在联系**：F4 给出物理不可逆，F7 把它转成节点盒约束，F32/F33/F45 则给出严格保持可行性的数值实现。

### 11.3 L-BFGS 紧凑表示 → SMW 快速求逆（F55, F87, F89, F110）

- L-BFGS 紧凑表示：

\[
\mathbf{B}_k = \mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k, \quad (20)
\]

- 子空间矩阵分解：

\[
\hat{\mathbf{B}}_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{Z}_k - \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k. \quad (28)
\]

- 由 SMW 得到逆：

\[
\hat{\mathbf{B}}^{-1}_k = (\hat{\mathbf{B}}^0_k)^{-1} + (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \left[\mathbf{I} - \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \right]^{-1} \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1}, \quad (29)
\]

- 附录理论依据：

\[
\hat{\mathbf{A}}^{-1} = \mathbf{A}^{-1} - \mathbf{A}^{-1}\mathbf{U}\left(\mathbf{I} + \mathbf{V}^{\mathrm{T}}\mathbf{A}^{-1}\mathbf{U}\right)^{-1}\mathbf{V}^{\mathrm{T}}\mathbf{A}^{-1}.
\]

**内在联系**：F55/F87 把问题写成“基矩阵 + 低秩修正”，F89 直接套用 F110，从而把大规模逆运算降到小规模矩阵逆。
