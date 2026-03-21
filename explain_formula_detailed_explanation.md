# explain.md 论文全部公式逐一解释（位置 + 完整表达 + 推导 + 联系）

> 对象文件：`/home/runner/work/Phasefield_gradient_projection_monolithic_solver/Phasefield_gradient_projection_monolithic_solver/explain.md`
>
> 说明：以下按公式出现顺序逐一解释，每条包含：
> 1) 公式出现位置；2) 公式完整表达；3) 公式含义；4) 必要推导过程；5) 与其他公式的内在联系。
>
> 统一渲染约定：行内符号与公式均使用 `[ ... ]` 形式，例如：[ \mathbf{u} ]、[ \Pi ]、[ \mathbf{u}_b \approx \alpha(\chi)(\mathbf{u}_{pres}-\mathbf{u}_f)+\mathbf{u}_{pres} ]。

## 全文主线（总览）

- **第1层（物理）**：由 [ \Pi ]、[ \Gamma_l ] 与不可逆约束给出“要最小化什么”。
- **第2层（离散）**：由插值、残量、切线给出“离散后怎么写”。
- **第3层（优化）**：由投影、L-BFGS、GCP、子空间/KKT 给出“数值上怎么解”。
- **第4层（算例）**：由能量与边界位移关系给出“结果怎么解释”。
- **第5层（附录）**：由 SMW 公式给出“为什么线代实现高效”。


## 逐条公式解释


### 公式 001

- **出现位置**：`explain.md` 第 21-23 行；章节：[ ## 1. Introduction ]；文中编号：[ (1) ]。
- **公式完整表达**：
[
\Pi (\pmb {u},d) = \int_{\Omega}\psi (\pmb {e}(\pmb {u}),d)\mathrm{d}\Omega + g_{c}\Gamma_{l}(d) - \int_{\Omega}\pmb {b}\cdot \pmb {u}\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot \pmb {u}\mathrm{d}\Gamma , \quad (1)
]
- **公式含义**：
  - 该式把系统目标定义为总势能 [ \Pi ]：内能、裂纹能与外力势共同决定平衡/演化。
  - 该式给出相场裂纹表面密度正则项，利用长度尺度 [ l ] 把尖锐裂纹转化为可微弥散带。
  - 本式位于章节 [## 1. Introduction]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：## 1. Introduction；During the past decade, the phase-field method becomes a popular technique for modeling fracture propagation, particularly in brittle materials, due to its capability of naturally handling complex crack geometry and propagation paths such as branching and merging [1-10]. Comparing with the traditional level-set based tracking strategy [11-13] that attempts to represent the discontinuity using an auxiliary scalar field, the phase-field method relies on a variational approach [14-16] and regularizes the sharp crack geometry in a diffusive manner. According to Francfort and Marigo [14], the energy functional of a fractured quasi-static elastic solid system is expressed as；where \(\pmb {u}(\pmb {x})\) represents the unknown vector displacement field, \(d(\pmb {x})\) represents the unknown scalar phase-field, \(\pmb{b}\) is the body force, \(\pmb{t}\) is the traction load, \(\pmb {e} = \nabla^{(s)}\pmb {u}\) is the small deformation linear strain tensor, \(\psi\) represents the strain energy density function, \(g_{c}\) is the critical energy release rate, and \(\Gamma_{l}\) is an approximation of the crack surface area. Particularly, the approximated crack surface \(\Gamma_{l}\) is defined as [1]；\Gamma_{l}(d) = \int_{\Omega}\gamma (d,\nabla d)\mathrm{d}\Omega = \int_{\Omega}\frac{1}{2l}\left(d^{2} + l^{2}\nabla d\cdot \nabla d\right)\mathrm{d}\Omega , \quad (2)
- **必要推导过程**：
  1. 从广义最小势能原理出发：内能项与外力势项相减得到平衡判据。
  2. 把裂纹面能以相场形式加入，得到可变分、可离散的目标泛函。
  3. 由相场正则化设定裂纹密度 [ \gamma(d,\nabla d)=\frac{1}{2l}(d^2+l^2|\nabla d|^2) ]。
  4. 对全域积分后得到近似裂纹表面积函数 [ \Gamma_l(d) ]。
- **与其他公式的内在联系**：
  - 向后连接到变分残量 [ r_u,r_d ] 与离散切线矩阵，是全篇物理起点。
  - 与编号公式 [ (1) ] 在文中直接互引。

### 公式 002

- **出现位置**：`explain.md` 第 27-29 行；章节：[ ## 1. Introduction ]；文中编号：[ (2) ]。
- **公式完整表达**：
[
\Gamma_{l}(d) = \int_{\Omega}\gamma (d,\nabla d)\mathrm{d}\Omega = \int_{\Omega}\frac{1}{2l}\left(d^{2} + l^{2}\nabla d\cdot \nabla d\right)\mathrm{d}\Omega , \quad (2)
]
- **公式含义**：
  - 该式给出相场裂纹表面密度正则项，利用长度尺度 [ l ] 把尖锐裂纹转化为可微弥散带。
  - 本式位于章节 [## 1. Introduction]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\Pi (\pmb {u},d) = \int_{\Omega}\psi (\pmb {e}(\pmb {u}),d)\mathrm{d}\Omega + g_{c}\Gamma_{l}(d) - \int_{\Omega}\pmb {b}\cdot \pmb {u}\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot \pmb {u}\mathrm{d}\Gamma , \quad (1)；where \(\pmb {u}(\pmb {x})\) represents the unknown vector displacement field, \(d(\pmb {x})\) represents the unknown scalar phase-field, \(\pmb{b}\) is the body force, \(\pmb{t}\) is the traction load, \(\pmb {e} = \nabla^{(s)}\pmb {u}\) is the small deformation linear strain tensor, \(\psi\) represents the strain energy density function, \(g_{c}\) is the critical energy release rate, and \(\Gamma_{l}\) is an approximation of the crack surface area. Particularly, the approximated crack surface \(\Gamma_{l}\) is defined as [1]；where \(\gamma (d,\nabla d)\) is considered as the crack surface density function, and \(l\) is the phase-field length-scale parameter.；For a quasi-static problem, let \(t\) represent the pseudo time instead of the real time. The pseudo time \(t\) enters the problem as the load step through time-dependent boundary conditions, for instance, the displacement \(\pmb {u} = \hat{\pmb{u}} (t)\) prescribed as the Dirichlet boundary condition or the pressure \(p = \hat{p} (t)\) prescribed as the Neumann boundary condition. During a typical time (load) step \(\left[t_{n},t_{n + 1}\right]\), let \(\left(\pmb{u}_{n},d_{n}\right)\in \mathbf{V}\times \mathbf{W}\) represent the solution from the previous time step, where \(\mathbf{V} = \mathbf{H}_{0}^{1}(\Omega)\) and \(\mathbf{W} = \mathbf{H}^{1}(\Omega)\). The primary unknown fields \(\left(\pmb{u}_{n + 1},d_{n + 1}\right)\in \mathbf{V}\times \mathbf{W}\) can be obtained from the following minimization [4,17]:
- **必要推导过程**：
  1. 由相场正则化设定裂纹密度 [ \gamma(d,\nabla d)=\frac{1}{2l}(d^2+l^2|\nabla d|^2) ]。
  2. 对全域积分后得到近似裂纹表面积函数 [ \Gamma_l(d) ]。
  3. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  4. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 向后连接到变分残量 [ r_u,r_d ] 与离散切线矩阵，是全篇物理起点。
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 与编号公式 [ (2) ] 在文中直接互引。

### 公式 003

- **出现位置**：`explain.md` 第 35-37 行；章节：[ ## 1. Introduction ]；文中编号：[ (3) ]。
- **公式完整表达**：
[
\left(\pmb{u}_{n + 1},d_{n + 1}\right) = \arg \min \Pi (\pmb {u},d) \quad (3)
]
- **公式含义**：
  - 该式把系统目标定义为总势能 [ \Pi ]：内能、裂纹能与外力势共同决定平衡/演化。
  - 该式把当前步未知量写成优化问题的极小点，明确“求解 = 约束最小化”。
  - 本式位于章节 [## 1. Introduction]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：where \(\gamma (d,\nabla d)\) is considered as the crack surface density function, and \(l\) is the phase-field length-scale parameter.；For a quasi-static problem, let \(t\) represent the pseudo time instead of the real time. The pseudo time \(t\) enters the problem as the load step through time-dependent boundary conditions, for instance, the displacement \(\pmb {u} = \hat{\pmb{u}} (t)\) prescribed as the Dirichlet boundary condition or the pressure \(p = \hat{p} (t)\) prescribed as the Neumann boundary condition. During a typical time (load) step \(\left[t_{n},t_{n + 1}\right]\), let \(\left(\pmb{u}_{n},d_{n}\right)\in \mathbf{V}\times \mathbf{W}\) represent the solution from the previous time step, where \(\mathbf{V} = \mathbf{H}_{0}^{1}(\Omega)\) and \(\mathbf{W} = \mathbf{H}^{1}(\Omega)\). The primary unknown fields \(\left(\pmb{u}_{n + 1},d_{n + 1}\right)\in \mathbf{V}\times \mathbf{W}\) can be obtained from the following minimization [4,17]:；subject to the inequality constraints；0\leq d_{n}\leq d_{n + 1}\leq 1. \quad (4)
- **必要推导过程**：
  1. 从广义最小势能原理出发：内能项与外力势项相减得到平衡判据。
  2. 把裂纹面能以相场形式加入，得到可变分、可离散的目标泛函。
  3. 把增量步未知量组合为优化变量向量。
  4. 施加不可逆约束后，写成带盒约束的极小化问题。
  5. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  6. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 向后连接到变分残量 [ r_u,r_d ] 与离散切线矩阵，是全篇物理起点。
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 与编号公式 [ (3) ] 在文中直接互引。

### 公式 004

- **出现位置**：`explain.md` 第 41-43 行；章节：[ ## 1. Introduction ]；文中编号：[ (4) ]。
- **公式完整表达**：
[
0\leq d_{n}\leq d_{n + 1}\leq 1. \quad (4)
]
- **公式含义**：
  - 该式在当前章节承担定义/变换/约束中的一个局部环节，用于连接前后公式。
  - 本式位于章节 [## 1. Introduction]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\left(\pmb{u}_{n + 1},d_{n + 1}\right) = \arg \min \Pi (\pmb {u},d) \quad (3)；subject to the inequality constraints；The above inequality constraints represent the following considerations. First, the phase-field cannot decrease in order to ensure the thermodynamic consistency. Second, the phase-field should be between \(d = 0\) (the undamaged state) and \(d = 1\) (the fully damaged state). The initial damage state \(d_{0}\) of the computational domain \(\Omega\) is known. For convenience, initially \((t = 0)\) the computational domain \(\Omega\) is assumed to be in the undamaged state, that is, \(d_{0}(\pmb{x}) = 0\).；Applying the phase-field approach to model fracture propagation encounters the following three challenges. First, it is well known that the energy functional shown in Eq. (1) is non-convex [15,16,18]. As a result, the Newton-based monolithic scheme usually faces convergence difficulties. Therefore, carefully designed numerical techniques are required to solve the nonlinear system. Second, cracks do not self-heal after formation. The phase-field (damage) irreversibility needs to be enforced for the thermodynamic consistency, which renders the phase-field formulation into a constrained minimization with inequality constraints, as shown in Eqs. (3) and (4). Third, the phase-field approach to represent fracture is intrinsically expensive, since the phase-field length-scale needs to be resolved around the crack region. Discretizing the entire computational domain with the same level of high mesh resolution is impractical, particularly for 3D problems. The adoption of adaptive mesh refinement technique is necessary. In the literature, a large body of work has been devoted to overcome the aforementioned three challenges.
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 与编号公式 [ (4) ] 在文中直接互引。

### 公式 005

- **出现位置**：`explain.md` 第 57-59 行；章节：[ ## 1. Introduction ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathrm{Proj}_C(\mathbf{x}_k = \mathbf{a}_k\nabla f(\mathbf{x}_k)).
]
- **公式含义**：
  - 该式引入投影算子，保证每次迭代都留在盒约束可行域内。
  - 本式位于章节 [## 1. Introduction]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：The phase-field crack formulation is computationally expensive since highly refined meshes around the crack path are needed to represent the sharp crack in a diffusive way. Globally refining the entire domain to achieve the required mesh resolution around the crack region is expensive for 2D problems and simply impractical for 3D problems. Several research efforts are devoted to develop the adaptive mesh refinement technique [4,36-38]. Particularly, Heister et al. [4] proposed a predictor-corrector local mesh adaptivity scheme. This scheme does not rely on the knowledge of the crack path a priori. Rather, it solely uses the phase-field solution as the local refinement criterion and performs the adaptive mesh refinement through multiple prediction-correction cycles. This scheme is appealing due to its computational efficiency, since the phase-field is part of the primary unknowns solved from the nonlinear system and no extra quantities are required.；The gradient projection method [39-42] is a special type of active set methods that can be used to solve the inequality-constrained optimization problem. This method has the advantage of allowing the set of active constraints to change rapidly during each iteration. Fig. 1 illustrates the basic idea of the gradient projection method. Let \(f(x)\) represent the objective function to be minimized. Let \(C\) represent the feasible region that is convex and formed by a series of inequality constraints. The gradient projection method typically involves an operation that projects the descent direction \(\nabla f(x)\) onto the convex feasible region \(C\), that is,；For a general convex feasible region, for example, as shown in Fig. 1(a), the projection operation \(\mathrm{Proj}_C(\cdot)\) incurs significant computational cost. Therefore, this method is not competitive comparing with other inequality-constrained minimization techniques such as the interior-point method. However, for inequality constraints that are imposed as bounds on the unknown vector \(\mathbf{x}\) as shown in Fig. 1(b), the gradient projection operation \(\mathrm{Proj}_C(\cdot)\) becomes trivial, making the gradient projection method extremely appealing [41,43,44]. The bound constraints, also known as the box constraints, can be expressed in the following component form,；\mathrm{lb}_i\leq x_i\leq \mathrm{ub}_i,
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。

### 公式 006

- **出现位置**：`explain.md` 第 63-65 行；章节：[ ## 1. Introduction ]；文中编号：[ (4) ]。
- **公式完整表达**：
[
\mathrm{lb}_i\leq x_i\leq \mathrm{ub}_i,
]
- **公式含义**：
  - 该式在当前章节承担定义/变换/约束中的一个局部环节，用于连接前后公式。
  - 本式位于章节 [## 1. Introduction]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathrm{Proj}_C(\mathbf{x}_k = \mathbf{a}_k\nabla f(\mathbf{x}_k)).；For a general convex feasible region, for example, as shown in Fig. 1(a), the projection operation \(\mathrm{Proj}_C(\cdot)\) incurs significant computational cost. Therefore, this method is not competitive comparing with other inequality-constrained minimization techniques such as the interior-point method. However, for inequality constraints that are imposed as bounds on the unknown vector \(\mathbf{x}\) as shown in Fig. 1(b), the gradient projection operation \(\mathrm{Proj}_C(\cdot)\) becomes trivial, making the gradient projection method extremely appealing [41,43,44]. The bound constraints, also known as the box constraints, can be expressed in the following component form,；where \(x_{i}\) represents the \(i\) th component of the unknown vector \(\mathbf{x}\) \(\mathrm{lb}_i\) and \(\mathrm{ub}_i\) represent its corresponding lower and upper bounds, respectively. Comparing the above box constraints with the inequality constraints shown in Eq. (4) imposed on the minimization of the phase-field energy functional, we can easily observe that after the spatial and temporal discretizations, the lower and upper bounds imposed on the discretized phase-field (nodal value) become；\mathrm{lb}_i = d_i^{(n)}\leq d_i^{(n + 1)}\leq 1 = \mathrm{ub}_i,
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 与编号公式 [ (4) ] 在文中直接互引。

### 公式 007

- **出现位置**：`explain.md` 第 69-71 行；章节：[ ## 1. Introduction ]；文中编号：[ (4) ]。
- **公式完整表达**：
[
\mathrm{lb}_i = d_i^{(n)}\leq d_i^{(n + 1)}\leq 1 = \mathrm{ub}_i,
]
- **公式含义**：
  - 该式在当前章节承担定义/变换/约束中的一个局部环节，用于连接前后公式。
  - 本式位于章节 [## 1. Introduction]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathrm{lb}_i\leq x_i\leq \mathrm{ub}_i,；where \(x_{i}\) represents the \(i\) th component of the unknown vector \(\mathbf{x}\) \(\mathrm{lb}_i\) and \(\mathrm{ub}_i\) represent its corresponding lower and upper bounds, respectively. Comparing the above box constraints with the inequality constraints shown in Eq. (4) imposed on the minimization of the phase-field energy functional, we can easily observe that after the spatial and temporal discretizations, the lower and upper bounds imposed on the discretized phase-field (nodal value) become；where \(d_i^{(n)}\) is known at the beginning of the current time (load) step, and \(d_i^{(n + 1)}\) is the unknown phase-field value that needs to be solved.；In this paper, a phase-field monolithic scheme based on the gradient projection method is proposed to enforce the phase-field irreversibility condition as box constraints. Comparing with the history variable approach, the gradient projection method preserves the variational structure of the energy functional and ensures that the obtained phase-field value is between 0 and 1. Unlike the penalty method that requires the proper adjustment of the penalty parameter, the proposed method does not involve extra algorithmic parameters. Furthermore, this method is combined with the limited-memory BFGS (L-BFGS) method to overcome the convergence difficulties arising from the non-convex energy functional. In the field of numerical optimization, the combination of the gradient projection method and the limited-memory BFGS method forms the so-called L-BFGS-B method [41,45], in which “L” stands for limited-memory and “B” stands for box constraints. Comparing with the existing phase-field monolithic schemes in the literature, such as the primal–dual active set method [4], the augmented Lagrangian method [6], and the interior-point method [19], the proposed monolithic scheme based on the L-BFGS-B method can rigorously and efficiently enforce the phase-field irreversibility. Moreover, the BFGS quasi-Newton approach can naturally handle non-convex optimizations without the necessity of introducing any heuristic parameters. The proposed monolithic scheme is further combined with the predictor–corrector adaptive mesh refinement technique [4]. Therefore, this paper presents an integrated monolithic scheme that is able to overcome the convergence difficulties associated with the non-convex energy functional, rigorously and efficiently enforce the phase-field irreversibility, and alleviate the computational cost through the adaptive mesh refinement.
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 与编号公式 [ (4) ] 在文中直接互引。

### 公式 008

- **出现位置**：`explain.md` 第 87-89 行；章节：[ ### 2.1. Phase-field formulation ]；文中编号：[ (1) ]、[ (5) ]。
- **公式完整表达**：
[
\psi(\pmb{\epsilon}, d) = [g(d) + k] \psi^{+}(\pmb{\epsilon}) + \psi^{-}(\pmb{\epsilon}), \quad (5)
]
- **公式含义**：
  - 该式指定退化函数 [ g(d) ]，控制损伤对张拉能的削弱强度。
  - 本式位于章节 [### 2.1. Phase-field formulation]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：### 2.1. Phase-field formulation；In order to consider the fracture tension–compression asymmetry, the strain energy density function used in the total energy functional in Eq. (1) is additively decomposed into two parts [1],；where \(\psi^{+}\) is the positive strain energy, \(\psi^{-}\) is the negative strain energy, \(k\) is a small non-negative number, and \(g(d)\) represents the phase-field degradation function and adopts the following form,；g(d) = (1 - d)^{2}. \quad (6)
- **必要推导过程**：
  1. 对完整应变能做张拉/压缩分裂。
  2. 令 [ g(d) ] 只衰减张拉部分，压缩部分保持不退化。
- **与其他公式的内在联系**：
  - 向前承接总势能密度定义，向后为离散残量/切线提供材料响应项。
  - 与编号公式 [ (1) ]、[ (5) ] 在文中直接互引。

### 公式 009

- **出现位置**：`explain.md` 第 93-95 行；章节：[ ### 2.1. Phase-field formulation ]；文中编号：[ (6) ]。
- **公式完整表达**：
[
g(d) = (1 - d)^{2}. \quad (6)
]
- **公式含义**：
  - 该式指定退化函数 [ g(d) ]，控制损伤对张拉能的削弱强度。
  - 本式位于章节 [### 2.1. Phase-field formulation]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\psi(\pmb{\epsilon}, d) = [g(d) + k] \psi^{+}(\pmb{\epsilon}) + \psi^{-}(\pmb{\epsilon}), \quad (5)；where \(\psi^{+}\) is the positive strain energy, \(\psi^{-}\) is the negative strain energy, \(k\) is a small non-negative number, and \(g(d)\) represents the phase-field degradation function and adopts the following form,；In this work, the small non-negative number \(k\) is set as zero (\(k = 0\)) in all the numerical examples.；The following operators are introduced to describe the constitutive relationship based on the additive decomposition of the strain energy density function,
- **必要推导过程**：
  1. 对完整应变能做张拉/压缩分裂。
  2. 令 [ g(d) ] 只衰减张拉部分，压缩部分保持不退化。
- **与其他公式的内在联系**：
  - 向前承接总势能密度定义，向后为离散残量/切线提供材料响应项。
  - 与编号公式 [ (6) ] 在文中直接互引。

### 公式 010

- **出现位置**：`explain.md` 第 101-103 行；章节：[ ### 2.1. Phase-field formulation ]；文中编号：未显式编号。
- **公式完整表达**：
[
\langle x\rangle_{+} = \frac{1}{2}(x + |x|), \quad \langle x\rangle_{-} = \frac{1}{2}(x - |x|), \quad H(x) = \begin{cases} 1 & x \geq 0, \\ 0 & x < 0. \end{cases}
]
- **公式含义**：
  - 该式通过正负括号函数分离张拉/压缩贡献，避免压缩驱动裂纹。
  - 本式位于章节 [### 2.1. Phase-field formulation]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：In this work, the small non-negative number \(k\) is set as zero (\(k = 0\)) in all the numerical examples.；The following operators are introduced to describe the constitutive relationship based on the additive decomposition of the strain energy density function,；The spectrum decomposition of the strain tensor \(\pmb{\epsilon}\) is expressed as；\pmb{\epsilon} = \sum_{\alpha} \epsilon_{\alpha} \mathbf{M}_{\alpha}, \quad \mathbf{M}_{\alpha} = \pmb{n}_{\alpha} \otimes \pmb{n}_{\alpha},
- **必要推导过程**：
  1. 由 [ \langle x \rangle_\pm=\frac12(x\pm|x|) ] 定义正负分解。
  2. 结合 Heaviside 函数可得到分段可导表达。
- **与其他公式的内在联系**：
  - 向前承接总势能密度定义，向后为离散残量/切线提供材料响应项。

### 公式 011

- **出现位置**：`explain.md` 第 107-109 行；章节：[ ### 2.1. Phase-field formulation ]；文中编号：未显式编号。
- **公式完整表达**：
[
\pmb{\epsilon} = \sum_{\alpha} \epsilon_{\alpha} \mathbf{M}_{\alpha}, \quad \mathbf{M}_{\alpha} = \pmb{n}_{\alpha} \otimes \pmb{n}_{\alpha},
]
- **公式含义**：
  - 该式使用主值谱分解，把应变张量投影到主方向后再做能量分裂。
  - 本式位于章节 [### 2.1. Phase-field formulation]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\langle x\rangle_{+} = \frac{1}{2}(x + |x|), \quad \langle x\rangle_{-} = \frac{1}{2}(x - |x|), \quad H(x) = \begin{cases} 1 & x \geq 0, \\ 0 & x < 0. \end{cases}；The spectrum decomposition of the strain tensor \(\pmb{\epsilon}\) is expressed as；where \(\epsilon_{\alpha}\) and \(\pmb{n}_{\alpha}\) represent a pair of eigenvalue and eigenvector. The positive and negative parts of the strain tensor are defined as,；\pmb{\epsilon}^{+} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{+} \mathbf{M}_{\alpha}, \quad \pmb{\epsilon}^{-} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{-} \mathbf{M}_{\alpha}.
- **必要推导过程**：
  1. 根据本式中的等号关系执行代入、求导或约束投影，可从上一式直接得到当前式。
- **与其他公式的内在联系**：
  - 向前承接总势能密度定义，向后为离散残量/切线提供材料响应项。

### 公式 012

- **出现位置**：`explain.md` 第 113-115 行；章节：[ ### 2.1. Phase-field formulation ]；文中编号：未显式编号。
- **公式完整表达**：
[
\pmb{\epsilon}^{+} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{+} \mathbf{M}_{\alpha}, \quad \pmb{\epsilon}^{-} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{-} \mathbf{M}_{\alpha}.
]
- **公式含义**：
  - 该式通过正负括号函数分离张拉/压缩贡献，避免压缩驱动裂纹。
  - 该式使用主值谱分解，把应变张量投影到主方向后再做能量分裂。
  - 本式位于章节 [### 2.1. Phase-field formulation]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\pmb{\epsilon} = \sum_{\alpha} \epsilon_{\alpha} \mathbf{M}_{\alpha}, \quad \mathbf{M}_{\alpha} = \pmb{n}_{\alpha} \otimes \pmb{n}_{\alpha},；where \(\epsilon_{\alpha}\) and \(\pmb{n}_{\alpha}\) represent a pair of eigenvalue and eigenvector. The positive and negative parts of the strain tensor are defined as,；Using the above definitions, the positive and negative parts of the strain energy are expressed as；\psi^{+}(\pmb{\epsilon}) = \frac{1}{2}\lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{+}^{2} + \mu \pmb{\epsilon}^{+} : \pmb{\epsilon}^{+}, \quad \psi^{-}(\pmb{\epsilon}) = \frac{1}{2} \lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{-}^{2} + \mu \pmb{\epsilon}^{-} : \pmb{\epsilon}^{-},
- **必要推导过程**：
  1. 由 [ \langle x \rangle_\pm=\frac12(x\pm|x|) ] 定义正负分解。
  2. 结合 Heaviside 函数可得到分段可导表达。
- **与其他公式的内在联系**：
  - 向前承接总势能密度定义，向后为离散残量/切线提供材料响应项。

### 公式 013

- **出现位置**：`explain.md` 第 119-121 行；章节：[ ### 2.1. Phase-field formulation ]；文中编号：未显式编号。
- **公式完整表达**：
[
\psi^{+}(\pmb{\epsilon}) = \frac{1}{2}\lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{+}^{2} + \mu \pmb{\epsilon}^{+} : \pmb{\epsilon}^{+}, \quad \psi^{-}(\pmb{\epsilon}) = \frac{1}{2} \lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{-}^{2} + \mu \pmb{\epsilon}^{-} : \pmb{\epsilon}^{-},
]
- **公式含义**：
  - 该式通过正负括号函数分离张拉/压缩贡献，避免压缩驱动裂纹。
  - 该式属于对偶/KKT 路线：通过 Schur 补先求乘子再回代主变量。
  - 本式位于章节 [### 2.1. Phase-field formulation]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\pmb{\epsilon}^{+} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{+} \mathbf{M}_{\alpha}, \quad \pmb{\epsilon}^{-} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{-} \mathbf{M}_{\alpha}.；Using the above definitions, the positive and negative parts of the strain energy are expressed as；where \(\lambda\) and \(\mu\) are the Lame parameters, and \(\mathrm{tr}\pmb{\epsilon}\) is the trace of the strain tensor. The stress tensor \(\pmb{\sigma}\) is derived as；\pmb{\sigma} = \frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial \pmb{\epsilon}} = [g(d) + k]\frac{\partial\psi^{+}(\pmb{\epsilon})}{\partial \pmb{\epsilon}} + \frac{\partial\psi^{-}(\pmb{\epsilon})}{\partial \pmb{\epsilon}} = [g(d) + k]\pmb{\sigma}^{+} + \pmb{\sigma}^{-},
- **必要推导过程**：
  1. 由 [ \langle x \rangle_\pm=\frac12(x\pm|x|) ] 定义正负分解。
  2. 结合 Heaviside 函数可得到分段可导表达。
  3. 写出 KKT 增广系统。
  4. 对主变量消元得到 Schur 补方程求乘子 [ \lambda ]，再回代 [ \Delta\mathbf{x}_k ]。
- **与其他公式的内在联系**：
  - 向前承接总势能密度定义，向后为离散残量/切线提供材料响应项。
  - 产生最终步长方向 [ \mathbf{p}_k ]，再进入线搜索更新 [ \mathbf{x}_{k+1} ]。

### 公式 014

- **出现位置**：`explain.md` 第 125-127 行；章节：[ ### 2.1. Phase-field formulation ]；文中编号：未显式编号。
- **公式完整表达**：
[
\pmb{\sigma} = \frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial \pmb{\epsilon}} = [g(d) + k]\frac{\partial\psi^{+}(\pmb{\epsilon})}{\partial \pmb{\epsilon}} + \frac{\partial\psi^{-}(\pmb{\epsilon})}{\partial \pmb{\epsilon}} = [g(d) + k]\pmb{\sigma}^{+} + \pmb{\sigma}^{-},
]
- **公式含义**：
  - 该式指定退化函数 [ g(d) ]，控制损伤对张拉能的削弱强度。
  - 该式给出一致切线（雅可比/Hessian 子块），直接决定牛顿或准牛顿收敛速度。
  - 本式位于章节 [### 2.1. Phase-field formulation]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\psi^{+}(\pmb{\epsilon}) = \frac{1}{2}\lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{+}^{2} + \mu \pmb{\epsilon}^{+} : \pmb{\epsilon}^{+}, \quad \psi^{-}(\pmb{\epsilon}) = \frac{1}{2} \lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{-}^{2} + \mu \pmb{\epsilon}^{-} : \pmb{\epsilon}^{-},；where \(\lambda\) and \(\mu\) are the Lame parameters, and \(\mathrm{tr}\pmb{\epsilon}\) is the trace of the strain tensor. The stress tensor \(\pmb{\sigma}\) is derived as；where；\pmb{\sigma}^{+} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{+}\mathbf{I} + 2\mu \pmb{\epsilon}^{+}, \quad \pmb{\sigma}^{-} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{-}\mathbf{I} + 2\mu \pmb{\epsilon}^{-}.
- **必要推导过程**：
  1. 对完整应变能做张拉/压缩分裂。
  2. 令 [ g(d) ] 只衰减张拉部分，压缩部分保持不退化。
  3. 对能量密度对 [ \boldsymbol{\epsilon} ] 求一阶导，得到应力。
  4. 若含退化函数，则链式法则带来 [ g(d)+k ] 缩放因子。
  5. 继续对应力对 [ \boldsymbol{\epsilon} ] 或对 [ d ] 求导。
  6. 得到一致切线，供牛顿/准牛顿迭代构造雅可比。
- **与其他公式的内在联系**：
  - 向前承接总势能密度定义，向后为离散残量/切线提供材料响应项。

### 公式 015

- **出现位置**：`explain.md` 第 131-133 行；章节：[ ### 2.1. Phase-field formulation ]；文中编号：未显式编号。
- **公式完整表达**：
[
\pmb{\sigma}^{+} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{+}\mathbf{I} + 2\mu \pmb{\epsilon}^{+}, \quad \pmb{\sigma}^{-} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{-}\mathbf{I} + 2\mu \pmb{\epsilon}^{-}.
]
- **公式含义**：
  - 该式通过正负括号函数分离张拉/压缩贡献，避免压缩驱动裂纹。
  - 该式把能量密度对应变求导，得到一致应力表达。
  - 该式属于对偶/KKT 路线：通过 Schur 补先求乘子再回代主变量。
  - 本式位于章节 [### 2.1. Phase-field formulation]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\pmb{\sigma} = \frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial \pmb{\epsilon}} = [g(d) + k]\frac{\partial\psi^{+}(\pmb{\epsilon})}{\partial \pmb{\epsilon}} + \frac{\partial\psi^{-}(\pmb{\epsilon})}{\partial \pmb{\epsilon}} = [g(d) + k]\pmb{\sigma}^{+} + \pmb{\sigma}^{-},；where；Due to the crack tension-compression asymmetry, the phase-field degradation function is exclusively applied to the positive part of the strain energy. Therefore, the stress-strain relationship is nonlinear, and the material tangent modulus is written as；\frac{\partial \pmb{\sigma}}{\partial \pmb{\epsilon}} = [g(d) + k]\frac{\partial \pmb{\sigma}^{+}}{\partial \pmb{\epsilon}} + \frac{\partial \pmb{\sigma}^{-}}{\partial \pmb{\epsilon}} = [g(d) + k]\left[\lambda H(\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes \mathbf{I} + 2\mu \mathbb{P}^{+}\right] + \left[\lambda H(-\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes \mathbf{I} + 2\mu \mathbb{P}^{-}\right].
- **必要推导过程**：
  1. 由 [ \langle x \rangle_\pm=\frac12(x\pm|x|) ] 定义正负分解。
  2. 结合 Heaviside 函数可得到分段可导表达。
  3. 对能量密度对 [ \boldsymbol{\epsilon} ] 求一阶导，得到应力。
  4. 若含退化函数，则链式法则带来 [ g(d)+k ] 缩放因子。
  5. 写出 KKT 增广系统。
  6. 对主变量消元得到 Schur 补方程求乘子 [ \lambda ]，再回代 [ \Delta\mathbf{x}_k ]。
- **与其他公式的内在联系**：
  - 向前承接总势能密度定义，向后为离散残量/切线提供材料响应项。
  - 产生最终步长方向 [ \mathbf{p}_k ]，再进入线搜索更新 [ \mathbf{x}_{k+1} ]。

### 公式 016

- **出现位置**：`explain.md` 第 137-139 行；章节：[ ### 2.1. Phase-field formulation ]；文中编号：未显式编号。
- **公式完整表达**：
[
\frac{\partial \pmb{\sigma}}{\partial \pmb{\epsilon}} = [g(d) + k]\frac{\partial \pmb{\sigma}^{+}}{\partial \pmb{\epsilon}} + \frac{\partial \pmb{\sigma}^{-}}{\partial \pmb{\epsilon}} = [g(d) + k]\left[\lambda H(\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes \mathbf{I} + 2\mu \mathbb{P}^{+}\right] + \left[\lambda H(-\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes \mathbf{I} + 2\mu \mathbb{P}^{-}\right].
]
- **公式含义**：
  - 该式指定退化函数 [ g(d) ]，控制损伤对张拉能的削弱强度。
  - 该式使用主值谱分解，把应变张量投影到主方向后再做能量分裂。
  - 该式给出一致切线（雅可比/Hessian 子块），直接决定牛顿或准牛顿收敛速度。
  - 该式属于对偶/KKT 路线：通过 Schur 补先求乘子再回代主变量。
  - 本式位于章节 [### 2.1. Phase-field formulation]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\pmb{\sigma}^{+} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{+}\mathbf{I} + 2\mu \pmb{\epsilon}^{+}, \quad \pmb{\sigma}^{-} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{-}\mathbf{I} + 2\mu \pmb{\epsilon}^{-}.；Due to the crack tension-compression asymmetry, the phase-field degradation function is exclusively applied to the positive part of the strain energy. Therefore, the stress-strain relationship is nonlinear, and the material tangent modulus is written as；Inside the above tangent modulus, the two fourth-order projection tensors \(\mathbb{P}^{+}\) and \(\mathbb{P}^{-}\) are defined as；\mathbb{P}^{+} = \frac{\partial \pmb{\epsilon}^{+}}{\partial \pmb{\epsilon}}, \quad \mathbb{P}^{-} = \frac{\partial \pmb{\epsilon}^{-}}{\partial \pmb{\epsilon}},
- **必要推导过程**：
  1. 对完整应变能做张拉/压缩分裂。
  2. 令 [ g(d) ] 只衰减张拉部分，压缩部分保持不退化。
  3. 对能量密度对 [ \boldsymbol{\epsilon} ] 求一阶导，得到应力。
  4. 若含退化函数，则链式法则带来 [ g(d)+k ] 缩放因子。
  5. 继续对应力对 [ \boldsymbol{\epsilon} ] 或对 [ d ] 求导。
  6. 得到一致切线，供牛顿/准牛顿迭代构造雅可比。
  7. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  8. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
  9. 写出 KKT 增广系统。
  10. 对主变量消元得到 Schur 补方程求乘子 [ \lambda ]，再回代 [ \Delta\mathbf{x}_k ]。
- **与其他公式的内在联系**：
  - 向前承接总势能密度定义，向后为离散残量/切线提供材料响应项。
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 产生最终步长方向 [ \mathbf{p}_k ]，再进入线搜索更新 [ \mathbf{x}_{k+1} ]。

### 公式 017

- **出现位置**：`explain.md` 第 143-145 行；章节：[ ### 2.1. Phase-field formulation ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbb{P}^{+} = \frac{\partial \pmb{\epsilon}^{+}}{\partial \pmb{\epsilon}}, \quad \mathbb{P}^{-} = \frac{\partial \pmb{\epsilon}^{-}}{\partial \pmb{\epsilon}},
]
- **公式含义**：
  - 该式给出一致切线（雅可比/Hessian 子块），直接决定牛顿或准牛顿收敛速度。
  - 本式位于章节 [### 2.1. Phase-field formulation]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\frac{\partial \pmb{\sigma}}{\partial \pmb{\epsilon}} = [g(d) + k]\frac{\partial \pmb{\sigma}^{+}}{\partial \pmb{\epsilon}} + \frac{\partial \pmb{\sigma}^{-}}{\partial \pmb{\epsilon}} = [g(d) + k]\left[\lambda H(\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes \mathbf{I} + 2\mu \mathbb{P}^{+}\right] + \left[\lambda H(-\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes \mathbf{I} + 2\mu \mathbb{P}^{-}\right].；Inside the above tangent modulus, the two fourth-order projection tensors \(\mathbb{P}^{+}\) and \(\mathbb{P}^{-}\) are defined as；the specific expressions of which can be found in [26].；Using the directional derivative, the first variation of the energy functional is written as
- **必要推导过程**：
  1. 继续对应力对 [ \boldsymbol{\epsilon} ] 或对 [ d ] 求导。
  2. 得到一致切线，供牛顿/准牛顿迭代构造雅可比。
- **与其他公式的内在联系**：
  - 该式在局部推导链中承接前式并提供下一式所需变量定义。

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
  - 该式把系统目标定义为总势能 [ \Pi ]：内能、裂纹能与外力势共同决定平衡/演化。
  - 该式把能量密度对应变求导，得到一致应力表达。
  - 本式位于章节 [### 2.1. Phase-field formulation]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：the specific expressions of which can be found in [26].；Using the directional derivative, the first variation of the energy functional is written as；Therefore, the weak form of the phase-field formulation is expressed as；\left\{ \begin{array}{ll} r_{\pmb{u}}(\pmb{u},d) = (\nabla^{(s)}\delta \pmb{u},\pmb{\sigma}) - (\delta \pmb{u},\pmb{b}) - (\delta \pmb{u},\pmb{t})_{\Gamma_{t}} = 0, \\ r_{d}(\pmb{u},d) = (\delta d,\frac{g_{c}}{l} d) + (\nabla \delta d,g_{c}l\nabla d) + (\delta d,g^{\prime}(d)\psi^{+}(\pmb{\epsilon})) = 0, \end{array} \right.
- **必要推导过程**：
  1. 从广义最小势能原理出发：内能项与外力势项相减得到平衡判据。
  2. 把裂纹面能以相场形式加入，得到可变分、可离散的目标泛函。
  3. 对能量密度对 [ \boldsymbol{\epsilon} ] 求一阶导，得到应力。
  4. 若含退化函数，则链式法则带来 [ g(d)+k ] 缩放因子。
  5. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  6. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 向后连接到变分残量 [ r_u,r_d ] 与离散切线矩阵，是全篇物理起点。
  - 向前承接总势能密度定义，向后为离散残量/切线提供材料响应项。
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 与编号公式 [ (7) ] 在文中直接互引。

### 公式 019

- **出现位置**：`explain.md` 第 162-164 行；章节：[ ### 2.1. Phase-field formulation ]；文中编号：未显式编号。
- **公式完整表达**：
[
\left\{ \begin{array}{ll} r_{\pmb{u}}(\pmb{u},d) = (\nabla^{(s)}\delta \pmb{u},\pmb{\sigma}) - (\delta \pmb{u},\pmb{b}) - (\delta \pmb{u},\pmb{t})_{\Gamma_{t}} = 0, \\ r_{d}(\pmb{u},d) = (\delta d,\frac{g_{c}}{l} d) + (\nabla \delta d,g_{c}l\nabla d) + (\delta d,g^{\prime}(d)\psi^{+}(\pmb{\epsilon})) = 0, \end{array} \right.
]
- **公式含义**：
  - 该式把能量密度对应变求导，得到一致应力表达。
  - 该式定义残量向量，残量为零对应离散平衡/最优性条件成立。
  - 本式位于章节 [### 2.1. Phase-field formulation]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\end{array} \quad (7)；Therefore, the weak form of the phase-field formulation is expressed as；subject to the phase-field irreversibility condition.；### 2.2. Finite element discretization
- **必要推导过程**：
  1. 对能量密度对 [ \boldsymbol{\epsilon} ] 求一阶导，得到应力。
  2. 若含退化函数，则链式法则带来 [ g(d)+k ] 缩放因子。
  3. 对离散势能对节点自由度求偏导得到残量。
  4. 残量组装后形成非线性代数方程组。
  5. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  6. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 向前承接总势能密度定义，向后为离散残量/切线提供材料响应项。
  - 把连续模型桥接到离散代数系统，直接作为第 3 节优化器输入。
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。

### 公式 020

- **出现位置**：`explain.md` 第 172-174 行；章节：[ ### 2.2. Finite element discretization ]；文中编号：未显式编号。
- **公式完整表达**：
[
\pmb {u} = \pmb{N}_{u_{A}}\pmb{u}_{A}, \quad \text{and} \quad d = N_{d_{A}}d_{A}.
]
- **公式含义**：
  - 该式是有限元插值关系：把场变量映射到节点自由度上。
  - 本式位于章节 [### 2.2. Finite element discretization]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：### 2.2. Finite element discretization；Let \(\mathbf{N}_{u_{A}}\) represent the vector-valued shape function for the displacement field \(\pmb{u}\) at node \(A\) and \(N_{d_{A}}\) represent the scalar-valued shape function for the phase-field \(d\) at node \(A\). Based on the shape functions \(\{N_{u_{A}},N_{d_{A}}\}\) and the nodal values \(\{\pmb{u}_{A},d_{A}\}\), the displacement field and the phase-field can be expressed as；Correspondingly, the displacement variation and the phase-field variation can be expressed as；\delta \pmb {u} = \pmb{N}_{u_{A}}\delta \pmb{u}_{A}, \quad \text{and} \quad \delta d = N_{d_{A}}\delta d_{A},
- **必要推导过程**：
  1. 以形函数展开场变量与其变分。
  2. 代入弱式后可把连续积分转换为离散矩阵-向量形式。
- **与其他公式的内在联系**：
  - 把连续模型桥接到离散代数系统，直接作为第 3 节优化器输入。

### 公式 021

- **出现位置**：`explain.md` 第 178-180 行；章节：[ ### 2.2. Finite element discretization ]；文中编号：[ (1) ]。
- **公式完整表达**：
[
\delta \pmb {u} = \pmb{N}_{u_{A}}\delta \pmb{u}_{A}, \quad \text{and} \quad \delta d = N_{d_{A}}\delta d_{A},
]
- **公式含义**：
  - 该式是有限元插值关系：把场变量映射到节点自由度上。
  - 本式位于章节 [### 2.2. Finite element discretization]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\pmb {u} = \pmb{N}_{u_{A}}\pmb{u}_{A}, \quad \text{and} \quad d = N_{d_{A}}d_{A}.；Correspondingly, the displacement variation and the phase-field variation can be expressed as；where the Einstein summation is used. Plug the above expressions into the total energy defined in Eq. (1), the objective function for the inequality-constrained minimization after the finite element discretization is expressed as；\begin{array}{rl}
- **必要推导过程**：
  1. 以形函数展开场变量与其变分。
  2. 代入弱式后可把连续积分转换为离散矩阵-向量形式。
- **与其他公式的内在联系**：
  - 把连续模型桥接到离散代数系统，直接作为第 3 节优化器输入。
  - 与编号公式 [ (1) ] 在文中直接互引。

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
  - 该式把系统目标定义为总势能 [ \Pi ]：内能、裂纹能与外力势共同决定平衡/演化。
  - 该式是有限元插值关系：把场变量映射到节点自由度上。
  - 本式位于章节 [### 2.2. Finite element discretization]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\delta \pmb {u} = \pmb{N}_{u_{A}}\delta \pmb{u}_{A}, \quad \text{and} \quad \delta d = N_{d_{A}}\delta d_{A},；where the Einstein summation is used. Plug the above expressions into the total energy defined in Eq. (1), the objective function for the inequality-constrained minimization after the finite element discretization is expressed as；The gradient of the discretized energy functional is derived as；\begin{array}{rl}
- **必要推导过程**：
  1. 从广义最小势能原理出发：内能项与外力势项相减得到平衡判据。
  2. 把裂纹面能以相场形式加入，得到可变分、可离散的目标泛函。
  3. 以形函数展开场变量与其变分。
  4. 代入弱式后可把连续积分转换为离散矩阵-向量形式。
  5. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  6. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 向后连接到变分残量 [ r_u,r_d ] 与离散切线矩阵，是全篇物理起点。
  - 把连续模型桥接到离散代数系统，直接作为第 3 节优化器输入。
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 与编号公式 [ (1) ]、[ (8) ] 在文中直接互引。

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
  - 该式把能量密度对应变求导，得到一致应力表达。
  - 该式是有限元插值关系：把场变量映射到节点自由度上。
  - 该式定义残量向量，残量为零对应离散平衡/最优性条件成立。
  - 本式位于章节 [### 2.2. Finite element discretization]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\end{array} \quad (8)；The gradient of the discretized energy functional is derived as；The Hessian matrix of the total energy after the finite element discretization is；\mathbf{K} = \nabla^{2}\Pi = \left[ \begin{array}{ll}\mathbf{K}_{uu} & \mathbf{K}_{ud}\\ \mathbf{K}_{du} & \mathbf{K}_{dd} \end{array} \right],
- **必要推导过程**：
  1. 对能量密度对 [ \boldsymbol{\epsilon} ] 求一阶导，得到应力。
  2. 若含退化函数，则链式法则带来 [ g(d)+k ] 缩放因子。
  3. 以形函数展开场变量与其变分。
  4. 代入弱式后可把连续积分转换为离散矩阵-向量形式。
  5. 对离散势能对节点自由度求偏导得到残量。
  6. 残量组装后形成非线性代数方程组。
  7. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  8. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 向前承接总势能密度定义，向后为离散残量/切线提供材料响应项。
  - 把连续模型桥接到离散代数系统，直接作为第 3 节优化器输入。
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 与编号公式 [ (9) ] 在文中直接互引。

### 公式 024

- **出现位置**：`explain.md` 第 204-206 行；章节：[ ### 2.2. Finite element discretization ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{K} = \nabla^{2}\Pi = \left[ \begin{array}{ll}\mathbf{K}_{uu} & \mathbf{K}_{ud}\\ \mathbf{K}_{du} & \mathbf{K}_{dd} \end{array} \right],
]
- **公式含义**：
  - 该式给出一致切线（雅可比/Hessian 子块），直接决定牛顿或准牛顿收敛速度。
  - 该式给出二阶导矩阵的块结构，揭示位移-相场耦合路径。
  - 本式位于章节 [### 2.2. Finite element discretization]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\end{array} \quad (9)；The Hessian matrix of the total energy after the finite element discretization is；where；\begin{array}{rl}
- **必要推导过程**：
  1. 继续对应力对 [ \boldsymbol{\epsilon} ] 或对 [ d ] 求导。
  2. 得到一致切线，供牛顿/准牛顿迭代构造雅可比。
  3. 对残量再求导，得到切线矩阵块。
  4. 按 [ uu,ud,du,dd ] 分块体现耦合与对称结构。
  5. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  6. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 把连续模型桥接到离散代数系统，直接作为第 3 节优化器输入。
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。

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
  - 该式给出一致切线（雅可比/Hessian 子块），直接决定牛顿或准牛顿收敛速度。
  - 该式是有限元插值关系：把场变量映射到节点自由度上。
  - 本式位于章节 [### 2.2. Finite element discretization]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{K} = \nabla^{2}\Pi = \left[ \begin{array}{ll}\mathbf{K}_{uu} & \mathbf{K}_{ud}\\ \mathbf{K}_{du} & \mathbf{K}_{dd} \end{array} \right],；where；In the above block matrices, the subscripts \(A\) and \(B\) represent the finite element nodal number, respectively. During a typical time (load) step \([t_n,t_{n + 1}]\), the phase-field irreversibility condition can be expressed as the following box constraints applied at individual finite element node, that is,；d_A^{(n)}\leq d_A\leq 1,
- **必要推导过程**：
  1. 对能量密度对 [ \boldsymbol{\epsilon} ] 求一阶导，得到应力。
  2. 若含退化函数，则链式法则带来 [ g(d)+k ] 缩放因子。
  3. 继续对应力对 [ \boldsymbol{\epsilon} ] 或对 [ d ] 求导。
  4. 得到一致切线，供牛顿/准牛顿迭代构造雅可比。
  5. 以形函数展开场变量与其变分。
  6. 代入弱式后可把连续积分转换为离散矩阵-向量形式。
  7. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  8. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 向前承接总势能密度定义，向后为离散残量/切线提供材料响应项。
  - 把连续模型桥接到离散代数系统，直接作为第 3 节优化器输入。
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 与编号公式 [ (10) ] 在文中直接互引。

### 公式 026

- **出现位置**：`explain.md` 第 219-221 行；章节：[ ### 2.2. Finite element discretization ]；文中编号：未显式编号。
- **公式完整表达**：
[
d_A^{(n)}\leq d_A\leq 1,
]
- **公式含义**：
  - 该式在当前章节承担定义/变换/约束中的一个局部环节，用于连接前后公式。
  - 本式位于章节 [### 2.2. Finite element discretization]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\end{array} \quad (10)；In the above block matrices, the subscripts \(A\) and \(B\) represent the finite element nodal number, respectively. During a typical time (load) step \([t_n,t_{n + 1}]\), the phase-field irreversibility condition can be expressed as the following box constraints applied at individual finite element node, that is,；where the subscript \(A\) represents the finite element nodal number, and \(d_A^{(n)}\) is known at the beginning of the current time step. Lastly, a diagonal block matrix is defined as,；\hat{\mathbf{K}} = \left[ \begin{array}{cc}\mathbf{K}_{uu} & \mathbf{0}\\ \mathbf{0} & \mathbf{K}_{dd} \end{array} \right], \quad (11)
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。

### 公式 027

- **出现位置**：`explain.md` 第 225-227 行；章节：[ ### 2.2. Finite element discretization ]；文中编号：[ (11) ]。
- **公式完整表达**：
[
\hat{\mathbf{K}} = \left[ \begin{array}{cc}\mathbf{K}_{uu} & \mathbf{0}\\ \mathbf{0} & \mathbf{K}_{dd} \end{array} \right], \quad (11)
]
- **公式含义**：
  - 该式给出一致切线（雅可比/Hessian 子块），直接决定牛顿或准牛顿收敛速度。
  - 该式给出二阶导矩阵的块结构，揭示位移-相场耦合路径。
  - 本式位于章节 [### 2.2. Finite element discretization]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：d_A^{(n)}\leq d_A\leq 1,；where the subscript \(A\) represents the finite element nodal number, and \(d_A^{(n)}\) is known at the beginning of the current time step. Lastly, a diagonal block matrix is defined as,；which is needed in the monolithic scheme presented in Section 3.；## 3. Gradient projection based monolithic scheme
- **必要推导过程**：
  1. 继续对应力对 [ \boldsymbol{\epsilon} ] 或对 [ d ] 求导。
  2. 得到一致切线，供牛顿/准牛顿迭代构造雅可比。
  3. 对残量再求导，得到切线矩阵块。
  4. 按 [ uu,ud,du,dd ] 分块体现耦合与对称结构。
  5. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  6. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 把连续模型桥接到离散代数系统，直接作为第 3 节优化器输入。
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 与编号公式 [ (11) ] 在文中直接互引。

### 公式 028

- **出现位置**：`explain.md` 第 239-241 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (12) ]。
- **公式完整表达**：
[
\left(\pmb{u}_A,d_A\right) = \arg \min \Pi (\pmb {u}_A,d_A) \quad (12)
]
- **公式含义**：
  - 该式把系统目标定义为总势能 [ \Pi ]：内能、裂纹能与外力势共同决定平衡/演化。
  - 该式把当前步未知量写成优化问题的极小点，明确“求解 = 约束最小化”。
  - 本式位于章节 [### 3.1. Algorithm overview]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：### 3.1. Algorithm overview；For the quasi-static phase-field fracture propagation problem, the pseudo time \(t\) is used to represent the actual load step. Inside the time (load) step \([t_n,t_{n + 1}]\), the solution \((\cdot)^{(n)}\) at the beginning is known. The goal is to determine the solution \((\cdot)^{(n + 1)}\) at the end of the current time step by solving the following constrained minimization after the finite element discretization,；subject to the inequality constraints；d_A^{(n)}\leq d_A\leq 1 \quad (13)
- **必要推导过程**：
  1. 从广义最小势能原理出发：内能项与外力势项相减得到平衡判据。
  2. 把裂纹面能以相场形式加入，得到可变分、可离散的目标泛函。
  3. 把增量步未知量组合为优化变量向量。
  4. 施加不可逆约束后，写成带盒约束的极小化问题。
  5. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  6. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 向后连接到变分残量 [ r_u,r_d ] 与离散切线矩阵，是全篇物理起点。
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 与编号公式 [ (12) ] 在文中直接互引。

### 公式 029

- **出现位置**：`explain.md` 第 245-247 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (13) ]。
- **公式完整表达**：
[
d_A^{(n)}\leq d_A\leq 1 \quad (13)
]
- **公式含义**：
  - 该式在当前章节承担定义/变换/约束中的一个局部环节，用于连接前后公式。
  - 本式位于章节 [### 3.1. Algorithm overview]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\left(\pmb{u}_A,d_A\right) = \arg \min \Pi (\pmb {u}_A,d_A) \quad (12)；subject to the inequality constraints；and the linear constraints in the following form；\left[\pmb {u}_A,d_A\right]^{\mathrm{T}} = \mathbf{C}\left[\pmb {u}_A,d_A\right]^{\mathrm{T}} + \pmb {k}. \quad (14)
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 与编号公式 [ (13) ] 在文中直接互引。

### 公式 030

- **出现位置**：`explain.md` 第 251-253 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (14) ]、[ (8) ]、[ (9) ]。
- **公式完整表达**：
[
\left[\pmb {u}_A,d_A\right]^{\mathrm{T}} = \mathbf{C}\left[\pmb {u}_A,d_A\right]^{\mathrm{T}} + \pmb {k}. \quad (14)
]
- **公式含义**：
  - 该式在当前章节承担定义/变换/约束中的一个局部环节，用于连接前后公式。
  - 本式位于章节 [### 3.1. Algorithm overview]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：d_A^{(n)}\leq d_A\leq 1 \quad (13)；and the linear constraints in the following form；In the above constrained minimization problem, \(\{\pmb {u}_A,d_A\}\) represents the finite element nodal solution of the displacement field and the phase-field. The objective function is shown in Eq. (8), and the gradient of the objective function is shown in Eq. (9). The linear constraints shown in Eq. (14) include the essential (Dirichlet) boundary conditions imposed on the domain boundary and the hanging-node constraints if the adaptive mesh refinement is adopted. The constraint coefficient matrix \(\mathbf{C}\) and the inhomogeneous constraints \(\mathbf{k}\) can be obtained from the essential (Dirichlet) boundary conditions and the adaptive mesh refinement process, for instance, see [26].；To overcome the convergence difficulties arising from the non-convexity of the objective function in Eq. (12), the BFGS method is adopted. As a type of quasi-Newton method, the BFGS method is originally based on the idea from Davidon [46,47] and further developed by Broyden, Fletcher, Goldfarb, and Shanno [48-51]. The basic idea is to construct a series of quadratic models to approximate the original objective function,
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 与编号公式 [ (14) ]、[ (8) ]、[ (9) ] 在文中直接互引。

### 公式 031

- **出现位置**：`explain.md` 第 259-261 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (12) ]、[ (15) ]。
- **公式完整表达**：
[
m_{k}(\pmb{x}) = \Pi_{k} + \pmb{r}_{k}^{\mathrm{T}}(\pmb{x} - \pmb{x}_{k}) + \frac{1}{2} (\pmb{x} - \pmb{x}_{k})^{\mathrm{T}}\mathbf{B}_{k}(\pmb{x} - \pmb{x}_{k}) \quad (15)
]
- **公式含义**：
  - 该式在当前章节承担定义/变换/约束中的一个局部环节，用于连接前后公式。
  - 本式位于章节 [### 3.1. Algorithm overview]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：In the above constrained minimization problem, \(\{\pmb {u}_A,d_A\}\) represents the finite element nodal solution of the displacement field and the phase-field. The objective function is shown in Eq. (8), and the gradient of the objective function is shown in Eq. (9). The linear constraints shown in Eq. (14) include the essential (Dirichlet) boundary conditions imposed on the domain boundary and the hanging-node constraints if the adaptive mesh refinement is adopted. The constraint coefficient matrix \(\mathbf{C}\) and the inhomogeneous constraints \(\mathbf{k}\) can be obtained from the essential (Dirichlet) boundary conditions and the adaptive mesh refinement process, for instance, see [26].；To overcome the convergence difficulties arising from the non-convexity of the objective function in Eq. (12), the BFGS method is adopted. As a type of quasi-Newton method, the BFGS method is originally based on the idea from Davidon [46,47] and further developed by Broyden, Fletcher, Goldfarb, and Shanno [48-51]. The basic idea is to construct a series of quadratic models to approximate the original objective function,；where the subscript \(k\) represents the \(k\) th iteration, \(\mathbf{x}\) represents the unknowns \(\{\pmb{u},d\}\), and \(\pmb{r}_k = \nabla \Pi_k\) represents the gradient of the objective function. The BFGS matrix \(\mathbf{B}_k\) is constructed from two rank-one updates on top of the previous matrix \(\mathbf{B}_{k - 1}\) and is positive definite [45], regardless whether the original objective function is convex or not. However, in the context of the finite element method, the BFGS matrix is fully dense and loses the underlying sparsity pattern associated with the finite element spatial discretization. Therefore, the conventional BFGS method cannot be directly applied to large-scale finite element simulations due to the memory required to store the fully dense matrix. This limitation motivated the development of the “limited-memory” BFGS (L-BFGS) method [27-29], which is further discussed in Section 3.2.；In order to consider the phase-field irreversibility condition, the quadratic model in Eq. (15) needs to satisfy the box constraints shown in Eq. (13). Let \(\mathbf{lb}_A = d_A^{(n)}\) and \(\mathbf{ub}_A = 1.0\) represent the lower and upper bounds of the phase-field at node \(A\) during time step \([t_n,t_{n + 1}]\). The gradient projection method [39-42], which is a type of active set method, can be used to effectively and efficiently enforce the box constraints. Define the projection operator \(\mathrm{Proj}_c(\cdot)\) as
- **必要推导过程**：
  1. 根据本式中的等号关系执行代入、求导或约束投影，可从上一式直接得到当前式。
- **与其他公式的内在联系**：
  - 与编号公式 [ (12) ]、[ (15) ] 在文中直接互引。

### 公式 032

- **出现位置**：`explain.md` 第 267-269 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (15) ]、[ (13) ]、[ (16) ]。
- **公式完整表达**：
[
\mathrm{Proj}_c(x_i,\mathrm{lb}_i,\mathrm{ub}_i) = \left\{ \begin{array}{ll}\mathrm{lb}_i & \mathrm{if}\; x_i < \mathrm{lb}_i,\\ x_i & \mathrm{if}\; \mathrm{lb}_i \leq x_i \leq \mathrm{ub}_i,\\ \mathrm{ub}_i & \mathrm{if}\; x_i > \mathrm{ub}_i, \end{array} \right. \quad (16)
]
- **公式含义**：
  - 该式引入投影算子，保证每次迭代都留在盒约束可行域内。
  - 本式位于章节 [### 3.1. Algorithm overview]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：where the subscript \(k\) represents the \(k\) th iteration, \(\mathbf{x}\) represents the unknowns \(\{\pmb{u},d\}\), and \(\pmb{r}_k = \nabla \Pi_k\) represents the gradient of the objective function. The BFGS matrix \(\mathbf{B}_k\) is constructed from two rank-one updates on top of the previous matrix \(\mathbf{B}_{k - 1}\) and is positive definite [45], regardless whether the original objective function is convex or not. However, in the context of the finite element method, the BFGS matrix is fully dense and loses the underlying sparsity pattern associated with the finite element spatial discretization. Therefore, the conventional BFGS method cannot be directly applied to large-scale finite element simulations due to the memory required to store the fully dense matrix. This limitation motivated the development of the “limited-memory” BFGS (L-BFGS) method [27-29], which is further discussed in Section 3.2.；In order to consider the phase-field irreversibility condition, the quadratic model in Eq. (15) needs to satisfy the box constraints shown in Eq. (13). Let \(\mathbf{lb}_A = d_A^{(n)}\) and \(\mathbf{ub}_A = 1.0\) represent the lower and upper bounds of the phase-field at node \(A\) during time step \([t_n,t_{n + 1}]\). The gradient projection method [39-42], which is a type of active set method, can be used to effectively and efficiently enforce the box constraints. Define the projection operator \(\mathrm{Proj}_c(\cdot)\) as；where the feasible region \(C\) formed by the phase-field inequality constraints is obviously convex. Based on the above projection operator, at the \(k\) th iteration, the steepest descent direction \(\pmb {r}_k = \nabla \Pi_k\) can be projected onto the feasible region with a negligible computational cost. The projected gradient forms a piecewise linear path, for instance, see Fig. 1(b), and can be written as；\pmb {x}(t) = \mathrm{Proj}_c(\pmb {x}_k - t\pmb {r}_k,\mathbf{lb},\mathbf{ub}),\quad t > 0, \quad (17)
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 与编号公式 [ (15) ]、[ (13) ]、[ (16) ] 在文中直接互引。

### 公式 033

- **出现位置**：`explain.md` 第 273-275 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (17) ]。
- **公式完整表达**：
[
\pmb {x}(t) = \mathrm{Proj}_c(\pmb {x}_k - t\pmb {r}_k,\mathbf{lb},\mathbf{ub}),\quad t > 0, \quad (17)
]
- **公式含义**：
  - 该式引入投影算子，保证每次迭代都留在盒约束可行域内。
  - 本式位于章节 [### 3.1. Algorithm overview]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathrm{Proj}_c(x_i,\mathrm{lb}_i,\mathrm{ub}_i) = \left\{ \begin{array}{ll}\mathrm{lb}_i & \mathrm{if}\; x_i < \mathrm{lb}_i,\\ x_i & \mathrm{if}\; \mathrm{lb}_i \leq x_i \leq \mathrm{ub}_i,\\ \mathrm{ub}_i & \mathrm{if}\; x_i > \mathrm{ub}_i, \end{array} \right. \quad (16)；where the feasible region \(C\) formed by the phase-field inequality constraints is obviously convex. Based on the above projection operator, at the \(k\) th iteration, the steepest descent direction \(\pmb {r}_k = \nabla \Pi_k\) can be projected onto the feasible region with a negligible computational cost. The projected gradient forms a piecewise linear path, for instance, see Fig. 1(b), and can be written as；in which the projection operator is applied to the vector \(\pmb{x}_k - t\pmb{r}_k\) component wise.；Following the piecewise linear path formed by the projected gradient, the quadratic model in the BFGS method, as shown in Eq. (15), is transformed into a univariate piecewise quadratic model,
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 与编号公式 [ (17) ] 在文中直接互引。

### 公式 034

- **出现位置**：`explain.md` 第 281-283 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (15) ]、[ (18) ]。
- **公式完整表达**：
[
p_k(t) = m_k(\pmb{x}(t)) = \Pi_k + \pmb{r}_k^\mathrm{T}(\pmb{x}(t) - \pmb{x}_k) + \frac{1}{2} (\pmb{x}(t) - \pmb{x}_k)^\mathrm{T}\mathbf{B}_k(\pmb{x}(t) - \pmb{x}_k). \quad (18)
]
- **公式含义**：
  - 该式在当前章节承担定义/变换/约束中的一个局部环节，用于连接前后公式。
  - 本式位于章节 [### 3.1. Algorithm overview]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：in which the projection operator is applied to the vector \(\pmb{x}_k - t\pmb{r}_k\) component wise.；Following the piecewise linear path formed by the projected gradient, the quadratic model in the BFGS method, as shown in Eq. (15), is transformed into a univariate piecewise quadratic model,；The first local minimizer \(\pmb{x}^c\) of the above model is called the generalized Cauchy point, the calculation of which is presented in Section 3.3.；For an arbitrary solution state \(\pmb{x}\) the active set \(\mathcal{A}(\pmb {x})\) of the inequality constraints shown in Eq. (13) include those components whose values are at either the lower bound or the upper bound. Note that the active set \(\mathcal{A}(\pmb {x})\) actually only contains component index, that is,
- **必要推导过程**：
  1. 根据本式中的等号关系执行代入、求导或约束投影，可从上一式直接得到当前式。
- **与其他公式的内在联系**：
  - 与编号公式 [ (15) ]、[ (18) ] 在文中直接互引。

### 公式 035

- **出现位置**：`explain.md` 第 289-291 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (13) ]、[ (15) ]。
- **公式完整表达**：
[
\mathcal{A}(\pmb {x}) = \{i:x_i = \mathrm{lb}_i\} \cup \{i:x_i = \mathrm{ub}_i\}.
]
- **公式含义**：
  - 该式刻画活动集（被约束卡住的自由度），是后续降维求解的关键。
  - 本式位于章节 [### 3.1. Algorithm overview]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：The first local minimizer \(\pmb{x}^c\) of the above model is called the generalized Cauchy point, the calculation of which is presented in Section 3.3.；For an arbitrary solution state \(\pmb{x}\) the active set \(\mathcal{A}(\pmb {x})\) of the inequality constraints shown in Eq. (13) include those components whose values are at either the lower bound or the upper bound. Note that the active set \(\mathcal{A}(\pmb {x})\) actually only contains component index, that is,；Based on the generalized Cauchy point \(\pmb{x}^c\), the variables in the active set \(\mathcal{A}_k(\pmb{x}^c)\) (located at the boundary of the box constraints) are held fixed. The solution \(\pmb{x}^*\) is obtained by minimizing the quadratic model shown in Eq. (15) in the subspace of the free variables (located inside the box constraints) subject to the box constraints, that is,；\pmb{x}^* = \arg \min m_k(\pmb {x})
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
  3. 根据变量是否触碰上下界识别活动集。
  4. 固定活动变量后，仅在自由子空间继续优化。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 与编号公式 [ (13) ]、[ (15) ] 在文中直接互引。

### 公式 036

- **出现位置**：`explain.md` 第 295-297 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (15) ]。
- **公式完整表达**：
[
\pmb{x}^* = \arg \min m_k(\pmb {x})
]
- **公式含义**：
  - 该式把当前步未知量写成优化问题的极小点，明确“求解 = 约束最小化”。
  - 本式位于章节 [### 3.1. Algorithm overview]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathcal{A}(\pmb {x}) = \{i:x_i = \mathrm{lb}_i\} \cup \{i:x_i = \mathrm{ub}_i\}.；Based on the generalized Cauchy point \(\pmb{x}^c\), the variables in the active set \(\mathcal{A}_k(\pmb{x}^c)\) (located at the boundary of the box constraints) are held fixed. The solution \(\pmb{x}^*\) is obtained by minimizing the quadratic model shown in Eq. (15) in the subspace of the free variables (located inside the box constraints) subject to the box constraints, that is,；subject to；x_{i} = x_{i}^{c},\;\forall i\in \mathcal{A}_{k}(\pmb{x}^{c})\quad \text{and}\quad \mathrm{lb}_{i}\leq x_{i}\leq \mathrm{ub}_{i},\;\forall i\notin \mathcal{A}_{k}(\pmb{x}^{c}).
- **必要推导过程**：
  1. 把增量步未知量组合为优化变量向量。
  2. 施加不可逆约束后，写成带盒约束的极小化问题。
- **与其他公式的内在联系**：
  - 与编号公式 [ (15) ] 在文中直接互引。

### 公式 037

- **出现位置**：`explain.md` 第 301-303 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：未显式编号。
- **公式完整表达**：
[
x_{i} = x_{i}^{c},\;\forall i\in \mathcal{A}_{k}(\pmb{x}^{c})\quad \text{and}\quad \mathrm{lb}_{i}\leq x_{i}\leq \mathrm{ub}_{i},\;\forall i\notin \mathcal{A}_{k}(\pmb{x}^{c}).
]
- **公式含义**：
  - 该式刻画活动集（被约束卡住的自由度），是后续降维求解的关键。
  - 本式位于章节 [### 3.1. Algorithm overview]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\pmb{x}^* = \arg \min m_k(\pmb {x})；subject to；Firstly ignoring the inequality constraints on the free variables, the above minimization can be solved using either the primal approach on the subspace of the free variables or the dual approach enforcing the active constraints via the Lagrange multipliers. Then, \(\pmb{x}^*\) is obtained by truncating the above solution to satisfy the box constraints. Section 3.4 discusses in detail the numerical techniques to solve the subspace minimization problem.；Based on the solution \(\pmb{x}^*\) of the above minimization, a search direction \(\pmb{p}_k\) is defined as
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
  3. 根据变量是否触碰上下界识别活动集。
  4. 固定活动变量后，仅在自由子空间继续优化。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。

### 公式 038

- **出现位置**：`explain.md` 第 309-311 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：未显式编号。
- **公式完整表达**：
[
\pmb{p}_k = \pmb{x}^* - \pmb{x}_k.
]
- **公式含义**：
  - 该式在当前章节承担定义/变换/约束中的一个局部环节，用于连接前后公式。
  - 本式位于章节 [### 3.1. Algorithm overview]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：Firstly ignoring the inequality constraints on the free variables, the above minimization can be solved using either the primal approach on the subspace of the free variables or the dual approach enforcing the active constraints via the Lagrange multipliers. Then, \(\pmb{x}^*\) is obtained by truncating the above solution to satisfy the box constraints. Section 3.4 discusses in detail the numerical techniques to solve the subspace minimization problem.；Based on the solution \(\pmb{x}^*\) of the above minimization, a search direction \(\pmb{p}_k\) is defined as；The search direction \(\pmb{p}_k\) is a descent direction of the objective function \(\Pi (\pmb{x})\), see [41] for the proof. The new solution \(\pmb{x}_{k + 1}\) is obtained as；\pmb{x}_{k + 1} = \pmb{x}_k + \alpha_k\pmb{p}_k.
- **必要推导过程**：
  1. 根据本式中的等号关系执行代入、求导或约束投影，可从上一式直接得到当前式。
- **与其他公式的内在联系**：
  - 该式在局部推导链中承接前式并提供下一式所需变量定义。

### 公式 039

- **出现位置**：`explain.md` 第 315-317 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：未显式编号。
- **公式完整表达**：
[
\pmb{x}_{k + 1} = \pmb{x}_k + \alpha_k\pmb{p}_k.
]
- **公式含义**：
  - 该式控制步长 [ \alpha_k ] 的选取，在下降性与曲率条件之间平衡。
  - 本式位于章节 [### 3.1. Algorithm overview]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\pmb{p}_k = \pmb{x}^* - \pmb{x}_k.；The search direction \(\pmb{p}_k\) is a descent direction of the objective function \(\Pi (\pmb{x})\), see [41] for the proof. The new solution \(\pmb{x}_{k + 1}\) is obtained as；The positive step length \(\alpha_k > 0\) is determined using the line search based on the strong Wolfe conditions, which include the sufficient decrease condition；\Pi_{k + 1} = \Pi (\pmb{x}_k + \alpha_k\pmb {p}_k) \leq \Pi_k + c_1\alpha_k\pmb {p}_k^\mathrm{T}\pmb{r}_k
- **必要推导过程**：
  1. 先给定可行方向 [ \mathbf{p}_k ]。
  2. 再用 Armijo/Wolfe 条件选 [ \alpha_k ]，兼顾下降与曲率。
- **与其他公式的内在联系**：
  - 与搜索方向方程联立，保证每步更新既可行又下降。

### 公式 040

- **出现位置**：`explain.md` 第 321-323 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：未显式编号。
- **公式完整表达**：
[
\Pi_{k + 1} = \Pi (\pmb{x}_k + \alpha_k\pmb {p}_k) \leq \Pi_k + c_1\alpha_k\pmb {p}_k^\mathrm{T}\pmb{r}_k
]
- **公式含义**：
  - 该式把系统目标定义为总势能 [ \Pi ]：内能、裂纹能与外力势共同决定平衡/演化。
  - 该式控制步长 [ \alpha_k ] 的选取，在下降性与曲率条件之间平衡。
  - 本式位于章节 [### 3.1. Algorithm overview]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\pmb{x}_{k + 1} = \pmb{x}_k + \alpha_k\pmb{p}_k.；The positive step length \(\alpha_k > 0\) is determined using the line search based on the strong Wolfe conditions, which include the sufficient decrease condition；and the curvature condition；\left|\pmb{r}_{k + 1}^\mathrm{T}\pmb {p}_k\right| = \left|\pmb {r}(\pmb{x}_k + \alpha_k\pmb {p}_k)^\mathrm{T}\pmb {p}_k\right| \leq c_2\left|\pmb{r}_k^\mathrm{T}\pmb {p}_k\right|.
- **必要推导过程**：
  1. 从广义最小势能原理出发：内能项与外力势项相减得到平衡判据。
  2. 把裂纹面能以相场形式加入，得到可变分、可离散的目标泛函。
  3. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  4. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
  5. 先给定可行方向 [ \mathbf{p}_k ]。
  6. 再用 Armijo/Wolfe 条件选 [ \alpha_k ]，兼顾下降与曲率。
- **与其他公式的内在联系**：
  - 向后连接到变分残量 [ r_u,r_d ] 与离散切线矩阵，是全篇物理起点。
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 与搜索方向方程联立，保证每步更新既可行又下降。

### 公式 041

- **出现位置**：`explain.md` 第 327-329 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：未显式编号。
- **公式完整表达**：
[
\left|\pmb{r}_{k + 1}^\mathrm{T}\pmb {p}_k\right| = \left|\pmb {r}(\pmb{x}_k + \alpha_k\pmb {p}_k)^\mathrm{T}\pmb {p}_k\right| \leq c_2\left|\pmb{r}_k^\mathrm{T}\pmb {p}_k\right|.
]
- **公式含义**：
  - 该式控制步长 [ \alpha_k ] 的选取，在下降性与曲率条件之间平衡。
  - 本式位于章节 [### 3.1. Algorithm overview]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\Pi_{k + 1} = \Pi (\pmb{x}_k + \alpha_k\pmb {p}_k) \leq \Pi_k + c_1\alpha_k\pmb {p}_k^\mathrm{T}\pmb{r}_k；and the curvature condition；To find a step length \(\alpha_k\) satisfying the strong Wolfe conditions, this work uses the implementation proposed by More and Thuente [52] with the parameters \(c_1 = 10^{-4}\) and \(c_2 = 0.9\).；Once the solution \(\pmb{x}_{k + 1}\) is obtained, a new iteration starts by updating the L-BFGS matrix \(\mathbf{B}_{k + 1}\), calculating the generalized Cauchy point \(\pmb{x}^c\) from the univariate piecewise quadratic model \(p_{k + 1}(t)\), updating the active set \(\mathcal{A}_{k + 1}(\pmb{x}^c)\), solving the subspace minimization to obtain the search direction \(\pmb{p}_{k + 1}\), and performing the line search based on the strong Wolfe conditions. The nonlinear iterations stop when the following convergence criteria are met simultaneously:
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
  3. 先给定可行方向 [ \mathbf{p}_k ]。
  4. 再用 Armijo/Wolfe 条件选 [ \alpha_k ]，兼顾下降与曲率。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 与搜索方向方程联立，保证每步更新既可行又下降。

### 公式 042

- **出现位置**：`explain.md` 第 336-338 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathcal{A}_{k + 1}(\pmb{x}^c) = \mathcal{A}_k(\pmb{x}^c).
]
- **公式含义**：
  - 该式刻画活动集（被约束卡住的自由度），是后续降维求解的关键。
  - 本式位于章节 [### 3.1. Algorithm overview]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：Once the solution \(\pmb{x}_{k + 1}\) is obtained, a new iteration starts by updating the L-BFGS matrix \(\mathbf{B}_{k + 1}\), calculating the generalized Cauchy point \(\pmb{x}^c\) from the univariate piecewise quadratic model \(p_{k + 1}(t)\), updating the active set \(\mathcal{A}_{k + 1}(\pmb{x}^c)\), solving the subspace minimization to obtain the search direction \(\pmb{p}_{k + 1}\), and performing the line search based on the strong Wolfe conditions. The nonlinear iterations stop when the following convergence criteria are met simultaneously:；1. The active set remains unchanged between two consecutive iterations:；2. The \(l_{2}\)-norm of the projected gradient is smaller than the prescribed tolerance:；\| \mathrm{Proj}_C(\mathbf{x}_{k + 1} - \mathbf{r}_{k + 1},\mathbf{lb},\mathbf{ub}) - \mathbf{x}_{k + 1}\|_{2} < \mathrm{tol}.
- **必要推导过程**：
  1. 根据变量是否触碰上下界识别活动集。
  2. 固定活动变量后，仅在自由子空间继续优化。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。

### 公式 043

- **出现位置**：`explain.md` 第 340-342 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：未显式编号。
- **公式完整表达**：
[
\| \mathrm{Proj}_C(\mathbf{x}_{k + 1} - \mathbf{r}_{k + 1},\mathbf{lb},\mathbf{ub}) - \mathbf{x}_{k + 1}\|_{2} < \mathrm{tol}.
]
- **公式含义**：
  - 该式引入投影算子，保证每次迭代都留在盒约束可行域内。
  - 该式定义残量向量，残量为零对应离散平衡/最优性条件成立。
  - 本式位于章节 [### 3.1. Algorithm overview]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathcal{A}_{k + 1}(\pmb{x}^c) = \mathcal{A}_k(\pmb{x}^c).；2. The \(l_{2}\)-norm of the projected gradient is smaller than the prescribed tolerance:；3. The \(l_{2}\)-norm of the solution increment is smaller than the prescribed tolerance:；\| \Delta \mathbf{x}_{k + 1}\|_{2} = \| \mathbf{x}_{k + 1} - \mathbf{x}_{k}\|_{2} < \mathrm{tol}.
- **必要推导过程**：
  1. 对离散势能对节点自由度求偏导得到残量。
  2. 残量组装后形成非线性代数方程组。
  3. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  4. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 把连续模型桥接到离散代数系统，直接作为第 3 节优化器输入。
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。

### 公式 044

- **出现位置**：`explain.md` 第 344-346 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (16) ]。
- **公式完整表达**：
[
\| \Delta \mathbf{x}_{k + 1}\|_{2} = \| \mathbf{x}_{k + 1} - \mathbf{x}_{k}\|_{2} < \mathrm{tol}.
]
- **公式含义**：
  - 该式在当前章节承担定义/变换/约束中的一个局部环节，用于连接前后公式。
  - 本式位于章节 [### 3.1. Algorithm overview]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\| \mathrm{Proj}_C(\mathbf{x}_{k + 1} - \mathbf{r}_{k + 1},\mathbf{lb},\mathbf{ub}) - \mathbf{x}_{k + 1}\|_{2} < \mathrm{tol}.；3. The \(l_{2}\)-norm of the solution increment is smaller than the prescribed tolerance:；Recall that the projection operator is defined in Eq. (16). For the \(i\) th component of the projected gradient,；\mathrm{Proj}_C(x_i - r_i,\mathrm{lb}_i,\mathrm{ub}_i) - x_i = \left\{ \begin{array}{ll}\mathrm{lb}_i - x_i & \mathrm{if}\; x_i - r_i < \mathrm{lb}_i,\\ - r_i & \mathrm{if}\; \mathrm{lb}_i \leq x_i - r_i \leq \mathrm{ub}_i,\\ \mathrm{ub}_i - x_i & \mathrm{if}\; x_i - r_i > \mathrm{ub}_i. \end{array} \right. \quad (19)
- **必要推导过程**：
  1. 根据本式中的等号关系执行代入、求导或约束投影，可从上一式直接得到当前式。
- **与其他公式的内在联系**：
  - 与编号公式 [ (16) ] 在文中直接互引。

### 公式 045

- **出现位置**：`explain.md` 第 350-352 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (16) ]、[ (19) ]。
- **公式完整表达**：
[
\mathrm{Proj}_C(x_i - r_i,\mathrm{lb}_i,\mathrm{ub}_i) - x_i = \left\{ \begin{array}{ll}\mathrm{lb}_i - x_i & \mathrm{if}\; x_i - r_i < \mathrm{lb}_i,\\ - r_i & \mathrm{if}\; \mathrm{lb}_i \leq x_i - r_i \leq \mathrm{ub}_i,\\ \mathrm{ub}_i - x_i & \mathrm{if}\; x_i - r_i > \mathrm{ub}_i. \end{array} \right. \quad (19)
]
- **公式含义**：
  - 该式引入投影算子，保证每次迭代都留在盒约束可行域内。
  - 本式位于章节 [### 3.1. Algorithm overview]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\| \Delta \mathbf{x}_{k + 1}\|_{2} = \| \mathbf{x}_{k + 1} - \mathbf{x}_{k}\|_{2} < \mathrm{tol}.；Recall that the projection operator is defined in Eq. (16). For the \(i\) th component of the projected gradient,；The second convergence criterion based on the projected gradient means that if the \(i\) th component of the projected gradient is within the feasible region (the box constraints), the \(l_{2}\)-norm uses the gradient value \(r_i\) itself. If the \(i\) th component of the projected gradient is outside the feasible region, the \(l_{2}\)-norm uses the distance between the solution component \(x_i\) and the boundary of the box constraints \(\mathrm{lb}_i\) or \(\mathrm{ub}_i\). This distance measures how well the inequality constraint is enforced.；**Comment 1.** The linear constraints shown in Eq. (14) include the essential (Dirichlet) boundary conditions imposed on the domain boundary and the hanging-node constraints if the adaptive mesh refinement is applied. Let \(\mathbf{x} = \{\mathbf{u},d\}\) represent the unknowns. Then, the linear constraints can be written in the following form
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 与编号公式 [ (16) ]、[ (19) ] 在文中直接互引。

### 公式 046

- **出现位置**：`explain.md` 第 358-360 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：[ (14) ]。
- **公式完整表达**：
[
\mathbf{x} = \mathbf{C}\mathbf{x} + \mathbf{k}.
]
- **公式含义**：
  - 该式在当前章节承担定义/变换/约束中的一个局部环节，用于连接前后公式。
  - 本式位于章节 [### 3.1. Algorithm overview]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：The second convergence criterion based on the projected gradient means that if the \(i\) th component of the projected gradient is within the feasible region (the box constraints), the \(l_{2}\)-norm uses the gradient value \(r_i\) itself. If the \(i\) th component of the projected gradient is outside the feasible region, the \(l_{2}\)-norm uses the distance between the solution component \(x_i\) and the boundary of the box constraints \(\mathrm{lb}_i\) or \(\mathrm{ub}_i\). This distance measures how well the inequality constraint is enforced.；**Comment 1.** The linear constraints shown in Eq. (14) include the essential (Dirichlet) boundary conditions imposed on the domain boundary and the hanging-node constraints if the adaptive mesh refinement is applied. Let \(\mathbf{x} = \{\mathbf{u},d\}\) represent the unknowns. Then, the linear constraints can be written in the following form；For a linear system \(\mathbf{A}\mathbf{x} = \mathbf{b}\) with the above set of constraints, we can instead solve the following modified linear system [53]；\left(\mathbf{C}^{\mathrm{T}}\mathbf{A}\mathbf{C} + \mathbf{I}_{d_{c}}\right)\hat{\mathbf{x}} = \mathbf{C}^{\mathrm{T}}(\mathbf{b} - \mathbf{A}\mathbf{k})
- **必要推导过程**：
  1. 根据本式中的等号关系执行代入、求导或约束投影，可从上一式直接得到当前式。
- **与其他公式的内在联系**：
  - 与编号公式 [ (14) ] 在文中直接互引。

### 公式 047

- **出现位置**：`explain.md` 第 364-366 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：未显式编号。
- **公式完整表达**：
[
\left(\mathbf{C}^{\mathrm{T}}\mathbf{A}\mathbf{C} + \mathbf{I}_{d_{c}}\right)\hat{\mathbf{x}} = \mathbf{C}^{\mathrm{T}}(\mathbf{b} - \mathbf{A}\mathbf{k})
]
- **公式含义**：
  - 该式在当前章节承担定义/变换/约束中的一个局部环节，用于连接前后公式。
  - 本式位于章节 [### 3.1. Algorithm overview]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{x} = \mathbf{C}\mathbf{x} + \mathbf{k}.；For a linear system \(\mathbf{A}\mathbf{x} = \mathbf{b}\) with the above set of constraints, we can instead solve the following modified linear system [53]；and then recover the true solution \(\mathbf{x}\) as；\mathbf{x} = \mathbf{C}\hat{\mathbf{x}} + \mathbf{k}.
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。

### 公式 048

- **出现位置**：`explain.md` 第 370-372 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{x} = \mathbf{C}\hat{\mathbf{x}} + \mathbf{k}.
]
- **公式含义**：
  - 该式在当前章节承担定义/变换/约束中的一个局部环节，用于连接前后公式。
  - 本式位于章节 [### 3.1. Algorithm overview]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\left(\mathbf{C}^{\mathrm{T}}\mathbf{A}\mathbf{C} + \mathbf{I}_{d_{c}}\right)\hat{\mathbf{x}} = \mathbf{C}^{\mathrm{T}}(\mathbf{b} - \mathbf{A}\mathbf{k})；and then recover the true solution \(\mathbf{x}\) as；In the modified linear system, the matrix \(\mathbf{I}_{d_c}\) is defined as；(\mathbf{I}_{d_c})_{ii} = \begin{cases} 1 & \text{if DoF } i \text{ is constrained}, \\ 0 & \text{otherwise}, \end{cases}
- **必要推导过程**：
  1. 根据本式中的等号关系执行代入、求导或约束投影，可从上一式直接得到当前式。
- **与其他公式的内在联系**：
  - 该式在局部推导链中承接前式并提供下一式所需变量定义。

### 公式 049

- **出现位置**：`explain.md` 第 375-377 行；章节：[ ### 3.1. Algorithm overview ]；文中编号：未显式编号。
- **公式完整表达**：
[
(\mathbf{I}_{d_c})_{ii} = \begin{cases} 1 & \text{if DoF } i \text{ is constrained}, \\ 0 & \text{otherwise}, \end{cases}
]
- **公式含义**：
  - 该式在当前章节承担定义/变换/约束中的一个局部环节，用于连接前后公式。
  - 本式位于章节 [### 3.1. Algorithm overview]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{x} = \mathbf{C}\hat{\mathbf{x}} + \mathbf{k}.；In the modified linear system, the matrix \(\mathbf{I}_{d_c}\) is defined as；where \(\mathcal{T}\) represents the set of the degrees of freedom at the constrained nodes, including the hanging nodes and the nodes prescribed with essential boundary conditions. For a detailed discussion about how to apply the above linear constraints, see [26].；### 3.2. Compact representation of limited-memory BFGS matrix
- **必要推导过程**：
  1. 根据本式中的等号关系执行代入、求导或约束投影，可从上一式直接得到当前式。
- **与其他公式的内在联系**：
  - 该式在局部推导链中承接前式并提供下一式所需变量定义。

### 公式 050

- **出现位置**：`explain.md` 第 384-386 行；章节：[ ### 3.2. Compact representation of limited-memory BFGS matrix ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{s}_k = \mathbf{x}_{k + 1} - \mathbf{x}_k,\quad \mathbf{y}_k = \mathbf{r}_{k + 1} - \mathbf{r}_k,
]
- **公式含义**：
  - 该式定义残量向量，残量为零对应离散平衡/最优性条件成立。
  - 该式属于 L-BFGS 曲率信息更新链，用少量历史向量近似全 Hessian。
  - 本式位于章节 [### 3.2. Compact representation of limited-memory BFGS matrix]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：### 3.2. Compact representation of limited-memory BFGS matrix；The detailed approach of the standard BFGS method can be found in many textbooks in the field of numerical optimization, for instance, see the classical textbook by Nocedal and Wright [45]. Let \(\mathbf{x} = \{\mathbf{u}_A,d_A\}\) represent the solutions of the displacement field and the phase-field after the finite element discretization. Recall that \(\Pi (\mathbf{x})\) represents the objective function (in this case the discretized total energy) to be minimized, and \(\mathbf{r} = \nabla \Pi (\mathbf{x})\) represents its gradient. At the \(k\) th BFGS iteration, the following two vectors are defined:；and the BFGS matrix is updated as；\mathbf{B}_{k + 1} = \mathbf{B}_k - \frac{\mathbf{B}_k\mathbf{s}_k\mathbf{s}_k^{\mathrm{T}}\mathbf{B}_k}{\mathbf{s}_k^{\mathrm{T}}\mathbf{B}_k\mathbf{s}_k} + \frac{\mathbf{y}_k\mathbf{y}_k^{\mathrm{T}}}{\mathbf{y}_k^{\mathrm{T}}\mathbf{s}_k}.
- **必要推导过程**：
  1. 对离散势能对节点自由度求偏导得到残量。
  2. 残量组装后形成非线性代数方程组。
  3. 由迭代差分定义 [ \mathbf{s}_k=\mathbf{x}_{k+1}-\mathbf{x}_k ] 与 [ \mathbf{y}_k=\mathbf{r}_{k+1}-\mathbf{r}_k ]。
  4. 以割线条件保证近似曲率与真实曲率一致。
- **与其他公式的内在联系**：
  - 把连续模型桥接到离散代数系统，直接作为第 3 节优化器输入。
  - 为 GCP 段内导数与子空间系统提供近似曲率信息。

### 公式 051

- **出现位置**：`explain.md` 第 390-392 行；章节：[ ### 3.2. Compact representation of limited-memory BFGS matrix ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{B}_{k + 1} = \mathbf{B}_k - \frac{\mathbf{B}_k\mathbf{s}_k\mathbf{s}_k^{\mathrm{T}}\mathbf{B}_k}{\mathbf{s}_k^{\mathrm{T}}\mathbf{B}_k\mathbf{s}_k} + \frac{\mathbf{y}_k\mathbf{y}_k^{\mathrm{T}}}{\mathbf{y}_k^{\mathrm{T}}\mathbf{s}_k}.
]
- **公式含义**：
  - 该式属于 L-BFGS 曲率信息更新链，用少量历史向量近似全 Hessian。
  - 本式位于章节 [### 3.2. Compact representation of limited-memory BFGS matrix]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{s}_k = \mathbf{x}_{k + 1} - \mathbf{x}_k,\quad \mathbf{y}_k = \mathbf{r}_{k + 1} - \mathbf{r}_k,；and the BFGS matrix is updated as；If the old BFGS matrix \(\mathbf{B}_k\) is positive definite and the vector pairs \(\mathbf{s}_k,\mathbf{y}_k\) satisfy the curvature condition,；\mathbf{s}_k^{\mathrm{T}}\mathbf{y}_k > 0,
- **必要推导过程**：
  1. 由迭代差分定义 [ \mathbf{s}_k=\mathbf{x}_{k+1}-\mathbf{x}_k ] 与 [ \mathbf{y}_k=\mathbf{r}_{k+1}-\mathbf{r}_k ]。
  2. 以割线条件保证近似曲率与真实曲率一致。
- **与其他公式的内在联系**：
  - 为 GCP 段内导数与子空间系统提供近似曲率信息。

### 公式 052

- **出现位置**：`explain.md` 第 396-398 行；章节：[ ### 3.2. Compact representation of limited-memory BFGS matrix ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{s}_k^{\mathrm{T}}\mathbf{y}_k > 0,
]
- **公式含义**：
  - 该式属于 L-BFGS 曲率信息更新链，用少量历史向量近似全 Hessian。
  - 本式位于章节 [### 3.2. Compact representation of limited-memory BFGS matrix]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{B}_{k + 1} = \mathbf{B}_k - \frac{\mathbf{B}_k\mathbf{s}_k\mathbf{s}_k^{\mathrm{T}}\mathbf{B}_k}{\mathbf{s}_k^{\mathrm{T}}\mathbf{B}_k\mathbf{s}_k} + \frac{\mathbf{y}_k\mathbf{y}_k^{\mathrm{T}}}{\mathbf{y}_k^{\mathrm{T}}\mathbf{s}_k}.；If the old BFGS matrix \(\mathbf{B}_k\) is positive definite and the vector pairs \(\mathbf{s}_k,\mathbf{y}_k\) satisfy the curvature condition,；then the updated BFGS matrix \(\mathbf{B}_{k + 1}\) is also positive definite [45]. Since \(\mathbf{s}_k\) and \(\mathbf{y}_k\) are two vectors with \(n\) components, where \(n\) represents the total number of degrees of freedom of the discretized system, each of the two terms \(\mathbf{B}_k\mathbf{s}_k\mathbf{s}_k^{\mathrm{T}}\mathbf{B}_k\) and \(\mathbf{y}_k\mathbf{y}_k^{\mathrm{T}}\) is a rank-one fully dense \(n\times n\) matrix, respectively. Therefore, in the BFGS matrix update, regardless whether \(\mathbf{B}_k\) is a sparse matrix or not, the updated matrix \(\mathbf{B}_{k + 1}\) is a fully dense matrix. Storing fully dense matrices is impractical due to the required memory. As a result, the conventional BFGS method cannot be directly used in large-scale finite element simulations. This limitation is the motivation for the development of the “limited-memory” feature for the BFGS method in the literature of numerical optimization, see [27-29].；In this work, the compact representation of the limited-memory BFGS (L-BFGS) matrix, originally proposed by Byrd et al. [29], is adopted to avoid the storage of a fully dense matrix. At the \(k\) th iteration, the BFGS matrix is essentially composed of an initial matrix \(\mathbf{B}^0_k\) and a series of vector-pairs \(\{\mathbf{s}_i, \mathbf{y}_i\}_{i=1}^{k-1}\). These vector-pairs contain the curvature information of the objective function. According to Nocedal and Wright [45], “the main idea of this method (L-BFGS) is to use curvature information from only the most recent iterations to construct the Hessian approximation. Curvature information from earlier iterations, which is less likely to be relevant to the actual behavior of the Hessian at the current iteration, is discarded in the interest of saving storage”. Instead of storing the entire history of the vector pairs, only the most recent \(m\) vector-pairs \(\{\mathbf{s}_i, \mathbf{y}_i\}_{i=k-m}^{k-1}\) need to be stored, where \(m\) is a small number. Obviously, if \(k < m\), all the vector-pairs are stored. According to the compact representation proposed by Byrd et al. [29], the following two \(n\times m\) correction matrices are defined as
- **必要推导过程**：
  1. 由迭代差分定义 [ \mathbf{s}_k=\mathbf{x}_{k+1}-\mathbf{x}_k ] 与 [ \mathbf{y}_k=\mathbf{r}_{k+1}-\mathbf{r}_k ]。
  2. 以割线条件保证近似曲率与真实曲率一致。
- **与其他公式的内在联系**：
  - 为 GCP 段内导数与子空间系统提供近似曲率信息。

### 公式 053

- **出现位置**：`explain.md` 第 404-406 行；章节：[ ### 3.2. Compact representation of limited-memory BFGS matrix ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{S}_k = [\mathbf{s}_{k-m} \cdots \mathbf{s}_{k-1}]
]
- **公式含义**：
  - 该式在当前章节承担定义/变换/约束中的一个局部环节，用于连接前后公式。
  - 本式位于章节 [### 3.2. Compact representation of limited-memory BFGS matrix]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：then the updated BFGS matrix \(\mathbf{B}_{k + 1}\) is also positive definite [45]. Since \(\mathbf{s}_k\) and \(\mathbf{y}_k\) are two vectors with \(n\) components, where \(n\) represents the total number of degrees of freedom of the discretized system, each of the two terms \(\mathbf{B}_k\mathbf{s}_k\mathbf{s}_k^{\mathrm{T}}\mathbf{B}_k\) and \(\mathbf{y}_k\mathbf{y}_k^{\mathrm{T}}\) is a rank-one fully dense \(n\times n\) matrix, respectively. Therefore, in the BFGS matrix update, regardless whether \(\mathbf{B}_k\) is a sparse matrix or not, the updated matrix \(\mathbf{B}_{k + 1}\) is a fully dense matrix. Storing fully dense matrices is impractical due to the required memory. As a result, the conventional BFGS method cannot be directly used in large-scale finite element simulations. This limitation is the motivation for the development of the “limited-memory” feature for the BFGS method in the literature of numerical optimization, see [27-29].；In this work, the compact representation of the limited-memory BFGS (L-BFGS) matrix, originally proposed by Byrd et al. [29], is adopted to avoid the storage of a fully dense matrix. At the \(k\) th iteration, the BFGS matrix is essentially composed of an initial matrix \(\mathbf{B}^0_k\) and a series of vector-pairs \(\{\mathbf{s}_i, \mathbf{y}_i\}_{i=1}^{k-1}\). These vector-pairs contain the curvature information of the objective function. According to Nocedal and Wright [45], “the main idea of this method (L-BFGS) is to use curvature information from only the most recent iterations to construct the Hessian approximation. Curvature information from earlier iterations, which is less likely to be relevant to the actual behavior of the Hessian at the current iteration, is discarded in the interest of saving storage”. Instead of storing the entire history of the vector pairs, only the most recent \(m\) vector-pairs \(\{\mathbf{s}_i, \mathbf{y}_i\}_{i=k-m}^{k-1}\) need to be stored, where \(m\) is a small number. Obviously, if \(k < m\), all the vector-pairs are stored. According to the compact representation proposed by Byrd et al. [29], the following two \(n\times m\) correction matrices are defined as；and；\mathbf{Y}_k = [\mathbf{y}_{k-m} \cdots \mathbf{y}_{k-1}].
- **必要推导过程**：
  1. 根据本式中的等号关系执行代入、求导或约束投影，可从上一式直接得到当前式。
- **与其他公式的内在联系**：
  - 该式在局部推导链中承接前式并提供下一式所需变量定义。

### 公式 054

- **出现位置**：`explain.md` 第 408-410 行；章节：[ ### 3.2. Compact representation of limited-memory BFGS matrix ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{Y}_k = [\mathbf{y}_{k-m} \cdots \mathbf{y}_{k-1}].
]
- **公式含义**：
  - 该式在当前章节承担定义/变换/约束中的一个局部环节，用于连接前后公式。
  - 本式位于章节 [### 3.2. Compact representation of limited-memory BFGS matrix]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{S}_k = [\mathbf{s}_{k-m} \cdots \mathbf{s}_{k-1}]；and；Then, the L-BFGS matrix can be represented in the following compact form,；\mathbf{B}_k = \mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k, \quad (20)
- **必要推导过程**：
  1. 根据本式中的等号关系执行代入、求导或约束投影，可从上一式直接得到当前式。
- **与其他公式的内在联系**：
  - 该式在局部推导链中承接前式并提供下一式所需变量定义。

### 公式 055

- **出现位置**：`explain.md` 第 414-416 行；章节：[ ### 3.2. Compact representation of limited-memory BFGS matrix ]；文中编号：[ (20) ]。
- **公式完整表达**：
[
\mathbf{B}_k = \mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k, \quad (20)
]
- **公式含义**：
  - 该式属于 L-BFGS 曲率信息更新链，用少量历史向量近似全 Hessian。
  - 本式位于章节 [### 3.2. Compact representation of limited-memory BFGS matrix]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{Y}_k = [\mathbf{y}_{k-m} \cdots \mathbf{y}_{k-1}].；Then, the L-BFGS matrix can be represented in the following compact form,；where；\mathbf{W}_k = [\mathbf{Y}_k \quad \mathbf{B}^0_k \mathbf{S}_k] \in \mathbb{R}^{n\times 2m}
- **必要推导过程**：
  1. 把多次 BFGS 更新改写为紧凑低秩形式 [ \mathbf{B}_k=\mathbf{B}^0_k-\mathbf{W}_k\mathbf{M}_k\mathbf{W}_k^T ]。
  2. 只需维护小矩阵 [ \mathbf{M}_k ] 与历史向量列块，降低存储和计算。
- **与其他公式的内在联系**：
  - 为 GCP 段内导数与子空间系统提供近似曲率信息。
  - 与编号公式 [ (20) ] 在文中直接互引。

### 公式 056

- **出现位置**：`explain.md` 第 420-422 行；章节：[ ### 3.2. Compact representation of limited-memory BFGS matrix ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{W}_k = [\mathbf{Y}_k \quad \mathbf{B}^0_k \mathbf{S}_k] \in \mathbb{R}^{n\times 2m}
]
- **公式含义**：
  - 该式属于 L-BFGS 曲率信息更新链，用少量历史向量近似全 Hessian。
  - 本式位于章节 [### 3.2. Compact representation of limited-memory BFGS matrix]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{B}_k = \mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k, \quad (20)；where；and；\mathbf{M}_k = \begin{bmatrix} -\mathbf{D}_k & \mathbf{L}^{\mathrm{T}}_k \\ \mathbf{L}_k & \mathbf{S}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{S}_k \end{bmatrix}^{-1} \in \mathbb{R}^{2m\times 2m}.
- **必要推导过程**：
  1. 把多次 BFGS 更新改写为紧凑低秩形式 [ \mathbf{B}_k=\mathbf{B}^0_k-\mathbf{W}_k\mathbf{M}_k\mathbf{W}_k^T ]。
  2. 只需维护小矩阵 [ \mathbf{M}_k ] 与历史向量列块，降低存储和计算。
- **与其他公式的内在联系**：
  - 为 GCP 段内导数与子空间系统提供近似曲率信息。

### 公式 057

- **出现位置**：`explain.md` 第 425-427 行；章节：[ ### 3.2. Compact representation of limited-memory BFGS matrix ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{M}_k = \begin{bmatrix} -\mathbf{D}_k & \mathbf{L}^{\mathrm{T}}_k \\ \mathbf{L}_k & \mathbf{S}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{S}_k \end{bmatrix}^{-1} \in \mathbb{R}^{2m\times 2m}.
]
- **公式含义**：
  - 该式属于 L-BFGS 曲率信息更新链，用少量历史向量近似全 Hessian。
  - 本式位于章节 [### 3.2. Compact representation of limited-memory BFGS matrix]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{W}_k = [\mathbf{Y}_k \quad \mathbf{B}^0_k \mathbf{S}_k] \in \mathbb{R}^{n\times 2m}；and；The matrix \(\mathbf{M}_k\) is a \(2m\times 2m\) matrix with the block matrices defined as；\mathbf{D}_k = \mathrm{diag}\{\mathbf{s}^{\mathrm{T}}_{k-m}\mathbf{y}_{k-m}, \ldots, \mathbf{s}^{\mathrm{T}}_{k-1}\mathbf{y}_{k-1}\} \in \mathbb{R}^{m\times m}
- **必要推导过程**：
  1. 把多次 BFGS 更新改写为紧凑低秩形式 [ \mathbf{B}_k=\mathbf{B}^0_k-\mathbf{W}_k\mathbf{M}_k\mathbf{W}_k^T ]。
  2. 只需维护小矩阵 [ \mathbf{M}_k ] 与历史向量列块，降低存储和计算。
- **与其他公式的内在联系**：
  - 为 GCP 段内导数与子空间系统提供近似曲率信息。

### 公式 058

- **出现位置**：`explain.md` 第 431-433 行；章节：[ ### 3.2. Compact representation of limited-memory BFGS matrix ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{D}_k = \mathrm{diag}\{\mathbf{s}^{\mathrm{T}}_{k-m}\mathbf{y}_{k-m}, \ldots, \mathbf{s}^{\mathrm{T}}_{k-1}\mathbf{y}_{k-1}\} \in \mathbb{R}^{m\times m}
]
- **公式含义**：
  - 该式属于 L-BFGS 曲率信息更新链，用少量历史向量近似全 Hessian。
  - 本式位于章节 [### 3.2. Compact representation of limited-memory BFGS matrix]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{M}_k = \begin{bmatrix} -\mathbf{D}_k & \mathbf{L}^{\mathrm{T}}_k \\ \mathbf{L}_k & \mathbf{S}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{S}_k \end{bmatrix}^{-1} \in \mathbb{R}^{2m\times 2m}.；The matrix \(\mathbf{M}_k\) is a \(2m\times 2m\) matrix with the block matrices defined as；and；\mathbf{L}_k \in \mathbb{R}^{m\times m}, \quad (\mathbf{L}_k)_{i,j} = \begin{cases} \mathbf{s}^{\mathrm{T}}_{i+k-m-1}\mathbf{y}_{j+k-m-1} & \text{if } i > j, \\ 0 & \text{if } i \leq j. \end{cases}
- **必要推导过程**：
  1. 把多次 BFGS 更新改写为紧凑低秩形式 [ \mathbf{B}_k=\mathbf{B}^0_k-\mathbf{W}_k\mathbf{M}_k\mathbf{W}_k^T ]。
  2. 只需维护小矩阵 [ \mathbf{M}_k ] 与历史向量列块，降低存储和计算。
- **与其他公式的内在联系**：
  - 为 GCP 段内导数与子空间系统提供近似曲率信息。

### 公式 059

- **出现位置**：`explain.md` 第 435-437 行；章节：[ ### 3.2. Compact representation of limited-memory BFGS matrix ]；文中编号：[ (11) ]。
- **公式完整表达**：
[
\mathbf{L}_k \in \mathbb{R}^{m\times m}, \quad (\mathbf{L}_k)_{i,j} = \begin{cases} \mathbf{s}^{\mathrm{T}}_{i+k-m-1}\mathbf{y}_{j+k-m-1} & \text{if } i > j, \\ 0 & \text{if } i \leq j. \end{cases}
]
- **公式含义**：
  - 该式属于 L-BFGS 曲率信息更新链，用少量历史向量近似全 Hessian。
  - 本式位于章节 [### 3.2. Compact representation of limited-memory BFGS matrix]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{D}_k = \mathrm{diag}\{\mathbf{s}^{\mathrm{T}}_{k-m}\mathbf{y}_{k-m}, \ldots, \mathbf{s}^{\mathrm{T}}_{k-1}\mathbf{y}_{k-1}\} \in \mathbb{R}^{m\times m}；and；At the \(k\) th iteration of the phase-field monolithic scheme, the initial BFGS matrix \(\mathbf{B}^0_k\) takes the diagonal block matrix defined in Eq. (11), that is,；\mathbf{B}^0_k = \hat{\mathbf{K}}^{(k)} = \begin{bmatrix} \mathbf{K}_{uu} & \mathbf{0} \\ \mathbf{0} & \mathbf{K}_{dd} \end{bmatrix}^{(k)}. \quad (21)
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
  3. 把多次 BFGS 更新改写为紧凑低秩形式 [ \mathbf{B}_k=\mathbf{B}^0_k-\mathbf{W}_k\mathbf{M}_k\mathbf{W}_k^T ]。
  4. 只需维护小矩阵 [ \mathbf{M}_k ] 与历史向量列块，降低存储和计算。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 为 GCP 段内导数与子空间系统提供近似曲率信息。
  - 与编号公式 [ (11) ] 在文中直接互引。

### 公式 060

- **出现位置**：`explain.md` 第 441-443 行；章节：[ ### 3.2. Compact representation of limited-memory BFGS matrix ]；文中编号：[ (11) ]、[ (21) ]、[ (20) ]。
- **公式完整表达**：
[
\mathbf{B}^0_k = \hat{\mathbf{K}}^{(k)} = \begin{bmatrix} \mathbf{K}_{uu} & \mathbf{0} \\ \mathbf{0} & \mathbf{K}_{dd} \end{bmatrix}^{(k)}. \quad (21)
]
- **公式含义**：
  - 该式给出一致切线（雅可比/Hessian 子块），直接决定牛顿或准牛顿收敛速度。
  - 该式给出二阶导矩阵的块结构，揭示位移-相场耦合路径。
  - 该式属于 L-BFGS 曲率信息更新链，用少量历史向量近似全 Hessian。
  - 本式位于章节 [### 3.2. Compact representation of limited-memory BFGS matrix]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{L}_k \in \mathbb{R}^{m\times m}, \quad (\mathbf{L}_k)_{i,j} = \begin{cases} \mathbf{s}^{\mathrm{T}}_{i+k-m-1}\mathbf{y}_{j+k-m-1} & \text{if } i > j, \\ 0 & \text{if } i \leq j. \end{cases}；At the \(k\) th iteration of the phase-field monolithic scheme, the initial BFGS matrix \(\mathbf{B}^0_k\) takes the diagonal block matrix defined in Eq. (11), that is,；The matrices \(\mathbf{K}_{uu}\) and \(\mathbf{K}_{dd}\) are both assembled from the finite element procedure, and therefore, are sparse. Also, both \(\mathbf{K}_{uu}\) and \(\mathbf{K}_{dd}\) are positive definite. Therefore, \(\mathbf{B}^0_k\) is sparse and positive definite, and can be stored in the memory without any issue. On the other hand, the second term \(\mathbf{W}_k\mathbf{M}_k\mathbf{W}^{\mathrm{T}}_k\) in Eq. (20) results a fully dense \(n\times n\) matrix and is never calculated or stored as a whole. Rather, only the matrices \(\mathbf{W}_k \in \mathbb{R}^{n\times 2m}\) and \(\mathbf{M}_k \in \mathbb{R}^{2m\times 2m}\) are stored individually. Recall that \(m \ll n\) is only a small number and typically ranges from 10 to 40 [45].；In summary, according to the compact representation of the L-BFGS matrix, at the \(k\) th iteration, only the following matrices need to be stored in the memory, including the sparse and positive definite matrix \(\mathbf{B}^0_k = \hat{\mathbf{K}}\), the \(2m\times 2m\) matrix \(\mathbf{M}_k\), and the \(n\times 2m\) matrix \(\mathbf{W}_k\). At the end of the \(k\) th iteration, the oldest vector-pair \(\{\mathbf{s}_{k-m}, \mathbf{y}_{k-m}\}\) is discarded, and the newest pair \(\{\mathbf{s}_k, \mathbf{y}_k\}\) is calculated and added into the list.
- **必要推导过程**：
  1. 继续对应力对 [ \boldsymbol{\epsilon} ] 或对 [ d ] 求导。
  2. 得到一致切线，供牛顿/准牛顿迭代构造雅可比。
  3. 对残量再求导，得到切线矩阵块。
  4. 按 [ uu,ud,du,dd ] 分块体现耦合与对称结构。
  5. 把多次 BFGS 更新改写为紧凑低秩形式 [ \mathbf{B}_k=\mathbf{B}^0_k-\mathbf{W}_k\mathbf{M}_k\mathbf{W}_k^T ]。
  6. 只需维护小矩阵 [ \mathbf{M}_k ] 与历史向量列块，降低存储和计算。
- **与其他公式的内在联系**：
  - 把连续模型桥接到离散代数系统，直接作为第 3 节优化器输入。
  - 为 GCP 段内导数与子空间系统提供近似曲率信息。
  - 与编号公式 [ (11) ]、[ (21) ]、[ (20) ] 在文中直接互引。

### 公式 061

- **出现位置**：`explain.md` 第 455-457 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：未显式编号。
- **公式完整表达**：
[
t_i = \begin{cases} (x^0_i - \mathrm{ub}_i)/r_i & \text{if } r_i < 0, \\ (x^0_i - \mathrm{lb}_i)/r_i & \text{if } r_i > 0, \\ +\infty & \text{if } r_i = 0. \end{cases}
]
- **公式含义**：
  - 该式属于广义 Cauchy 点构造过程，用分段一维二次模型追踪可行下降轨迹。
  - 本式位于章节 [### 3.3. Generalized Cauchy point]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：### 3.3. Generalized Cauchy point；Once the L-BFGS matrix is formed using the compact representation, the next step is to locate the generalized Cauchy point, which is defined as the first local minimizer along the piecewise linear path formed by the projected gradient. The calculation of the generalized Cauchy point mainly follows the steps detailed in [41]. However, several places are modified to consider the impact of the initial BFGS matrix \(\mathbf{B}^0 = \hat{\mathbf{K}}\) on the calculation. At the \(k\) th nonlinear iteration, we use \(\mathbf{x}^0\) to represent the solution at the beginning of the current iteration. For notation convenience, we drop the iteration index \(k\) in all the vectors and matrices such that \(\mathbf{r}= \mathbf{r}_k\) and \(\mathbf{B}= \mathbf{B}_k\). Instead, the subscript index is used to represent the component of a vector. For instance, \(r_i\) represents the \(i\) th component of the gradient vector \(\mathbf{r}\). Also, the concept of “break points” is introduced to represent the intersections between the projected gradient and the box constraints. The schematic in Fig. 2 illustrates the path of the projected gradient and the break points. For the \(i\) th coordinate direction, the break point along the projected gradient can be represented by a positive scalar as；Note that if there is no upper or lower bound on the \(i\) th direction, then set \(\mathrm{ub}_i = +\infty\) or \(\mathrm{lb}_i = -\infty\). After all the break points are computed for each direction \(\{t_i\}_{i = 1}^n\), they need to be reordered to form an increasing order set \(\{t^{(j)}: t^{(j)} \leq t^{(j+1)}\}_{j=1}^n\). Note that the subscript \(i\) represents the \(i\) th vector component, while the superscript \(j\) represents the \(j\) th element in the ordered set. Since the break points need to be sorted multiple times during each iteration, it is important to efficiently perform the sort operation using the appropriate data structure and algorithm, for instance, the heap sort algorithm with a complexity of \(\mathcal{O}(n\log n)\). Using the break points, the piecewise linear path formed by the projected gradient \(\mathrm{Proj}_C(\mathbf{x}^0 - t\mathbf{r},\mathbf{lb},\mathbf{ub})\) can be expressed in the following component form；x_i(t) = x_i^{0} - \min \{t,t_i\} r_i.
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
  3. 对每个自由度计算触边断点 [ t_i ]。
  4. 按断点顺序推进分段路径 [ \mathbf{x}(t) ]，逐步扩大活动集。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 输出广义 Cauchy 点与活动集，作为子空间最小化的初值与约束结构。

### 公式 062

- **出现位置**：`explain.md` 第 461-463 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：未显式编号。
- **公式完整表达**：
[
x_i(t) = x_i^{0} - \min \{t,t_i\} r_i.
]
- **公式含义**：
  - 该式属于广义 Cauchy 点构造过程，用分段一维二次模型追踪可行下降轨迹。
  - 本式位于章节 [### 3.3. Generalized Cauchy point]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：t_i = \begin{cases} (x^0_i - \mathrm{ub}_i)/r_i & \text{if } r_i < 0, \\ (x^0_i - \mathrm{lb}_i)/r_i & \text{if } r_i > 0, \\ +\infty & \text{if } r_i = 0. \end{cases}；Note that if there is no upper or lower bound on the \(i\) th direction, then set \(\mathrm{ub}_i = +\infty\) or \(\mathrm{lb}_i = -\infty\). After all the break points are computed for each direction \(\{t_i\}_{i = 1}^n\), they need to be reordered to form an increasing order set \(\{t^{(j)}: t^{(j)} \leq t^{(j+1)}\}_{j=1}^n\). Note that the subscript \(i\) represents the \(i\) th vector component, while the superscript \(j\) represents the \(j\) th element in the ordered set. Since the break points need to be sorted multiple times during each iteration, it is important to efficiently perform the sort operation using the appropriate data structure and algorithm, for instance, the heap sort algorithm with a complexity of \(\mathcal{O}(n\log n)\). Using the break points, the piecewise linear path formed by the projected gradient \(\mathrm{Proj}_C(\mathbf{x}^0 - t\mathbf{r},\mathbf{lb},\mathbf{ub})\) can be expressed in the following component form；Inside an interval \([t^{(j-1)}, t^{(j)}]\) along the piecewise linear path formed by the projected gradient, the solution vector can be expressed as a function of the scalar \(t\)；\mathbf{x}(t) = \mathbf{x}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)},
- **必要推导过程**：
  1. 对每个自由度计算触边断点 [ t_i ]。
  2. 按断点顺序推进分段路径 [ \mathbf{x}(t) ]，逐步扩大活动集。
- **与其他公式的内在联系**：
  - 输出广义 Cauchy 点与活动集，作为子空间最小化的初值与约束结构。

### 公式 063

- **出现位置**：`explain.md` 第 467-469 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{x}(t) = \mathbf{x}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)},
]
- **公式含义**：
  - 该式属于广义 Cauchy 点构造过程，用分段一维二次模型追踪可行下降轨迹。
  - 本式位于章节 [### 3.3. Generalized Cauchy point]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：x_i(t) = x_i^{0} - \min \{t,t_i\} r_i.；Inside an interval \([t^{(j-1)}, t^{(j)}]\) along the piecewise linear path formed by the projected gradient, the solution vector can be expressed as a function of the scalar \(t\)；where \(\mathbf{x}^{(j-1)} = \mathbf{x}(t^{(j-1)})\), \(\Delta t = t - t^{(j-1)}\), and \(\mathbf{d}^{(j-1)}\) is the direction of the projected gradient in that interval. Let \(\mathbf{z}^{(j-1)} = \mathbf{x}^{(j-1)} - \mathbf{x}^0\), and the univariate quadratic model on this interval can be expressed as；\begin{array}{rl}
- **必要推导过程**：
  1. 对每个自由度计算触边断点 [ t_i ]。
  2. 按断点顺序推进分段路径 [ \mathbf{x}(t) ]，逐步扩大活动集。
- **与其他公式的内在联系**：
  - 输出广义 Cauchy 点与活动集，作为子空间最小化的初值与约束结构。

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
  - 该式定义残量向量，残量为零对应离散平衡/最优性条件成立。
  - 该式属于广义 Cauchy 点构造过程，用分段一维二次模型追踪可行下降轨迹。
  - 本式位于章节 [### 3.3. Generalized Cauchy point]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{x}(t) = \mathbf{x}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)},；where \(\mathbf{x}^{(j-1)} = \mathbf{x}(t^{(j-1)})\), \(\Delta t = t - t^{(j-1)}\), and \(\mathbf{d}^{(j-1)}\) is the direction of the projected gradient in that interval. Let \(\mathbf{z}^{(j-1)} = \mathbf{x}^{(j-1)} - \mathbf{x}^0\), and the univariate quadratic model on this interval can be expressed as；where \(f_{j-1}, f_{j-1}^{\prime}\) and \(f_{j-1}^{\prime \prime}\) represent the coefficients of the quadratic polynomial \(\hat{p}(\Delta t)\).；As discussed in Section 3.2, the updated L-BFGS matrix \(\mathbf{B}\) is positive definite if the vector-pairs \(\{\mathbf{s}_i,\mathbf{y}_i\}_{i=k-m}^{k-1}\) satisfy the curvature condition shown in Eq. (18). Therefore, the coefficient \(f_{j-1}^{\prime \prime}\) is positive, that is,
- **必要推导过程**：
  1. 对离散势能对节点自由度求偏导得到残量。
  2. 残量组装后形成非线性代数方程组。
  3. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  4. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
  5. 对每个自由度计算触边断点 [ t_i ]。
  6. 按断点顺序推进分段路径 [ \mathbf{x}(t) ]，逐步扩大活动集。
- **与其他公式的内在联系**：
  - 把连续模型桥接到离散代数系统，直接作为第 3 节优化器输入。
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 输出广义 Cauchy 点与活动集，作为子空间最小化的初值与约束结构。
  - 与编号公式 [ (22) ] 在文中直接互引。

### 公式 065

- **出现位置**：`explain.md` 第 487-489 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：[ (18) ]。
- **公式完整表达**：
[
f_{j-1}^{\prime \prime} = \mathbf{d}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j-1)} > 0.
]
- **公式含义**：
  - 该式在当前章节承担定义/变换/约束中的一个局部环节，用于连接前后公式。
  - 本式位于章节 [### 3.3. Generalized Cauchy point]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：where \(f_{j-1}, f_{j-1}^{\prime}\) and \(f_{j-1}^{\prime \prime}\) represent the coefficients of the quadratic polynomial \(\hat{p}(\Delta t)\).；As discussed in Section 3.2, the updated L-BFGS matrix \(\mathbf{B}\) is positive definite if the vector-pairs \(\{\mathbf{s}_i,\mathbf{y}_i\}_{i=k-m}^{k-1}\) satisfy the curvature condition shown in Eq. (18). Therefore, the coefficient \(f_{j-1}^{\prime \prime}\) is positive, that is,；The quadratic polynomial \(\hat{p} (\Delta t)\) reaches to its minimum value at the critical point \(\Delta t^*\)；\hat{p}'(\Delta t) = 0 \Rightarrow \Delta t^* = -f_{j-1}' / f_{j-1}^{\prime \prime}.
- **必要推导过程**：
  1. 根据本式中的等号关系执行代入、求导或约束投影，可从上一式直接得到当前式。
- **与其他公式的内在联系**：
  - 与编号公式 [ (18) ] 在文中直接互引。

### 公式 066

- **出现位置**：`explain.md` 第 493-495 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：未显式编号。
- **公式完整表达**：
[
\hat{p}'(\Delta t) = 0 \Rightarrow \Delta t^* = -f_{j-1}' / f_{j-1}^{\prime \prime}.
]
- **公式含义**：
  - 该式属于广义 Cauchy 点构造过程，用分段一维二次模型追踪可行下降轨迹。
  - 本式位于章节 [### 3.3. Generalized Cauchy point]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：f_{j-1}^{\prime \prime} = \mathbf{d}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j-1)} > 0.；The quadratic polynomial \(\hat{p} (\Delta t)\) reaches to its minimum value at the critical point \(\Delta t^*\)；Depending on the value of \(\Delta t^*\), there are the following three scenarios:；1. If \(\Delta t^* \in [0, t^{(j)} - t^{(j-1)})\), the generalized Cauchy point \(\mathbf{x}^c = \mathbf{x}(t^{(j-1)} + \Delta t^*)\).
- **必要推导过程**：
  1. 对每个自由度计算触边断点 [ t_i ]。
  2. 按断点顺序推进分段路径 [ \mathbf{x}(t) ]，逐步扩大活动集。
  3. 在当前分段上写出一元二次模型及一、二阶导。
  4. 若解析极小点位于段内则终止，否则进入下一段。
- **与其他公式的内在联系**：
  - 输出广义 Cauchy 点与活动集，作为子空间最小化的初值与约束结构。

### 公式 067

- **出现位置**：`explain.md` 第 504-506 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathcal{A}(\mathbf{x}^c) = \{i : t_i \leq t^{(j-1)}, i = 1, \ldots, n\}.
]
- **公式含义**：
  - 该式刻画活动集（被约束卡住的自由度），是后续降维求解的关键。
  - 该式属于广义 Cauchy 点构造过程，用分段一维二次模型追踪可行下降轨迹。
  - 本式位于章节 [### 3.3. Generalized Cauchy point]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：3. If \(\Delta t^* \geq t^{(j)} - t^{(j-1)}\), the search of the generalized Cauchy point moves onto the next interval (line segment) of the projected gradient \([t^{(j)}, t^{(j+1)}]\).；For the first two scenarios, the set of active constraints include those components whose break point \(t\)-value is no larger than \(t^{(j-1)}\), that is,；During the next interval \([t^{(j)}, t^{(j+1)}]\), it is not necessary to calculate the values of \(f'_j\) and \(f''_j\) from scratch. Rather, they can be obtained from updating the old values of \(f'_{j-1}\) and \(f''_{j-1}\), which could significantly lower the computational cost. Let；\Delta t^{(j-1)} = t^{(j)} - t^{(j-1)}, \quad \mathbf{x}^{(j)} = \mathbf{x}^{(j-1)} + \Delta t^{(j-1)}\mathbf{d}^{(j-1)}.
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
  3. 根据变量是否触碰上下界识别活动集。
  4. 固定活动变量后，仅在自由子空间继续优化。
  5. 对每个自由度计算触边断点 [ t_i ]。
  6. 按断点顺序推进分段路径 [ \mathbf{x}(t) ]，逐步扩大活动集。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 输出广义 Cauchy 点与活动集，作为子空间最小化的初值与约束结构。

### 公式 068

- **出现位置**：`explain.md` 第 509-511 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：未显式编号。
- **公式完整表达**：
[
\Delta t^{(j-1)} = t^{(j)} - t^{(j-1)}, \quad \mathbf{x}^{(j)} = \mathbf{x}^{(j-1)} + \Delta t^{(j-1)}\mathbf{d}^{(j-1)}.
]
- **公式含义**：
  - 该式属于广义 Cauchy 点构造过程，用分段一维二次模型追踪可行下降轨迹。
  - 本式位于章节 [### 3.3. Generalized Cauchy point]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathcal{A}(\mathbf{x}^c) = \{i : t_i \leq t^{(j-1)}, i = 1, \ldots, n\}.；During the next interval \([t^{(j)}, t^{(j+1)}]\), it is not necessary to calculate the values of \(f'_j\) and \(f''_j\) from scratch. Rather, they can be obtained from updating the old values of \(f'_{j-1}\) and \(f''_{j-1}\), which could significantly lower the computational cost. Let；Moving from the previous interval \([t^{(j-1)}, t^{(j)}]\) to the new interval \([t^{(j)}, t^{(j+1)}]\), assume that the constraint associated with the \(b\) th vector component becomes active, that is, \(t_b = t^{(j)}\). The direction of the projected gradient in the current interval becomes；\mathbf{d}^{(j)} = \mathbf{d}^{(j-1)} + r_b \mathbf{e}_b,
- **必要推导过程**：
  1. 对每个自由度计算触边断点 [ t_i ]。
  2. 按断点顺序推进分段路径 [ \mathbf{x}(t) ]，逐步扩大活动集。
- **与其他公式的内在联系**：
  - 输出广义 Cauchy 点与活动集，作为子空间最小化的初值与约束结构。

### 公式 069

- **出现位置**：`explain.md` 第 514-516 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：[ (22) ]。
- **公式完整表达**：
[
\mathbf{d}^{(j)} = \mathbf{d}^{(j-1)} + r_b \mathbf{e}_b,
]
- **公式含义**：
  - 该式在当前章节承担定义/变换/约束中的一个局部环节，用于连接前后公式。
  - 本式位于章节 [### 3.3. Generalized Cauchy point]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\Delta t^{(j-1)} = t^{(j)} - t^{(j-1)}, \quad \mathbf{x}^{(j)} = \mathbf{x}^{(j-1)} + \Delta t^{(j-1)}\mathbf{d}^{(j-1)}.；Moving from the previous interval \([t^{(j-1)}, t^{(j)}]\) to the new interval \([t^{(j)}, t^{(j+1)}]\), assume that the constraint associated with the \(b\) th vector component becomes active, that is, \(t_b = t^{(j)}\). The direction of the projected gradient in the current interval becomes；where \(\mathbf{e}_b\) is the \(b\) th unit vector, and \(r_b\) is the \(b\) th component of the gradient vector \(\mathbf{r}\). Similar to Eq. (22), the new values \(f'_j\) and \(f''_j\) are written as；f'_j = \mathbf{r}^{\mathrm{T}}\mathbf{d}^{(j)} + \mathbf{d}^{(j)}{}^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j)} = f'_{j-1} + \Delta t^{(j-1)} f''_{j-1} + r^2_b + r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j)} \quad (23)
- **必要推导过程**：
  1. 根据本式中的等号关系执行代入、求导或约束投影，可从上一式直接得到当前式。
- **与其他公式的内在联系**：
  - 与编号公式 [ (22) ] 在文中直接互引。

### 公式 070

- **出现位置**：`explain.md` 第 519-521 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：[ (22) ]、[ (23) ]。
- **公式完整表达**：
[
f'_j = \mathbf{r}^{\mathrm{T}}\mathbf{d}^{(j)} + \mathbf{d}^{(j)}{}^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j)} = f'_{j-1} + \Delta t^{(j-1)} f''_{j-1} + r^2_b + r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j)} \quad (23)
]
- **公式含义**：
  - 该式定义残量向量，残量为零对应离散平衡/最优性条件成立。
  - 该式属于广义 Cauchy 点构造过程，用分段一维二次模型追踪可行下降轨迹。
  - 本式位于章节 [### 3.3. Generalized Cauchy point]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{d}^{(j)} = \mathbf{d}^{(j-1)} + r_b \mathbf{e}_b,；where \(\mathbf{e}_b\) is the \(b\) th unit vector, and \(r_b\) is the \(b\) th component of the gradient vector \(\mathbf{r}\). Similar to Eq. (22), the new values \(f'_j\) and \(f''_j\) are written as；and；f''_j = \mathbf{d}^{(j)}{}^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j)} = f''_{j-1} + 2r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j-1)} + r^2_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{e}_b. \quad (24)
- **必要推导过程**：
  1. 对离散势能对节点自由度求偏导得到残量。
  2. 残量组装后形成非线性代数方程组。
  3. 对每个自由度计算触边断点 [ t_i ]。
  4. 按断点顺序推进分段路径 [ \mathbf{x}(t) ]，逐步扩大活动集。
  5. 在当前分段上写出一元二次模型及一、二阶导。
  6. 若解析极小点位于段内则终止，否则进入下一段。
- **与其他公式的内在联系**：
  - 把连续模型桥接到离散代数系统，直接作为第 3 节优化器输入。
  - 输出广义 Cauchy 点与活动集，作为子空间最小化的初值与约束结构。
  - 与编号公式 [ (22) ]、[ (23) ] 在文中直接互引。

### 公式 071

- **出现位置**：`explain.md` 第 523-525 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：[ (24) ]、[ (20) ]。
- **公式完整表达**：
[
f''_j = \mathbf{d}^{(j)}{}^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j)} = f''_{j-1} + 2r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j-1)} + r^2_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{e}_b. \quad (24)
]
- **公式含义**：
  - 该式属于广义 Cauchy 点构造过程，用分段一维二次模型追踪可行下降轨迹。
  - 本式位于章节 [### 3.3. Generalized Cauchy point]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：f'_j = \mathbf{r}^{\mathrm{T}}\mathbf{d}^{(j)} + \mathbf{d}^{(j)}{}^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j)} = f'_{j-1} + \Delta t^{(j-1)} f''_{j-1} + r^2_b + r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j)} \quad (23)；and；Based on the compact representation of the L-BFGS matrix, as shown in Eq. (20),；\mathbf{B} = \mathbf{B}^0 - \mathbf{W}\mathbf{M}\mathbf{W}^{\mathrm{T}}.
- **必要推导过程**：
  1. 在当前分段上写出一元二次模型及一、二阶导。
  2. 若解析极小点位于段内则终止，否则进入下一段。
- **与其他公式的内在联系**：
  - 输出广义 Cauchy 点与活动集，作为子空间最小化的初值与约束结构。
  - 与编号公式 [ (24) ]、[ (20) ] 在文中直接互引。

### 公式 072

- **出现位置**：`explain.md` 第 528-530 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：[ (20) ]。
- **公式完整表达**：
[
\mathbf{B} = \mathbf{B}^0 - \mathbf{W}\mathbf{M}\mathbf{W}^{\mathrm{T}}.
]
- **公式含义**：
  - 该式在当前章节承担定义/变换/约束中的一个局部环节，用于连接前后公式。
  - 本式位于章节 [### 3.3. Generalized Cauchy point]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：f''_j = \mathbf{d}^{(j)}{}^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j)} = f''_{j-1} + 2r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j-1)} + r^2_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{e}_b. \quad (24)；Based on the compact representation of the L-BFGS matrix, as shown in Eq. (20),；Let；\mathbf{w}_b = \mathbf{W}^{\mathrm{T}}\mathbf{e}_b, \quad \mathbf{p}^{(j)} = \mathbf{W}^{\mathrm{T}}\mathbf{d}^{(j)} = \mathbf{W}^{\mathrm{T}}(\mathbf{d}^{(j-1)} + r_b \mathbf{e}_b) = \mathbf{p}^{(j-1)} + r_b \mathbf{w}_b, \quad \mathbf{c}^{(j)} = \mathbf{W}^{\mathrm{T}}\mathbf{z}^{(j)} = \mathbf{W}^{\mathrm{T}}(\mathbf{z}^{(j-1)} + \Delta t^{(j-1)}\mathbf{d}^{(j-1)}) = \mathbf{c}^{(j-1)} + \Delta t^{(j-1)}\mathbf{p}^{(j-1)}.
- **必要推导过程**：
  1. 根据本式中的等号关系执行代入、求导或约束投影，可从上一式直接得到当前式。
- **与其他公式的内在联系**：
  - 与编号公式 [ (20) ] 在文中直接互引。

### 公式 073

- **出现位置**：`explain.md` 第 532-534 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：[ (23) ]、[ (24) ]。
- **公式完整表达**：
[
\mathbf{w}_b = \mathbf{W}^{\mathrm{T}}\mathbf{e}_b, \quad \mathbf{p}^{(j)} = \mathbf{W}^{\mathrm{T}}\mathbf{d}^{(j)} = \mathbf{W}^{\mathrm{T}}(\mathbf{d}^{(j-1)} + r_b \mathbf{e}_b) = \mathbf{p}^{(j-1)} + r_b \mathbf{w}_b, \quad \mathbf{c}^{(j)} = \mathbf{W}^{\mathrm{T}}\mathbf{z}^{(j)} = \mathbf{W}^{\mathrm{T}}(\mathbf{z}^{(j-1)} + \Delta t^{(j-1)}\mathbf{d}^{(j-1)}) = \mathbf{c}^{(j-1)} + \Delta t^{(j-1)}\mathbf{p}^{(j-1)}.
]
- **公式含义**：
  - 该式属于广义 Cauchy 点构造过程，用分段一维二次模型追踪可行下降轨迹。
  - 本式位于章节 [### 3.3. Generalized Cauchy point]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{B} = \mathbf{B}^0 - \mathbf{W}\mathbf{M}\mathbf{W}^{\mathrm{T}}.；Let；Plug all the above relationships into Eqs. (23) and (24), the coefficients of the quadratic polynomial \(\hat{p}(\Delta t)\) for the interval \([t^{(j)}, t^{(j+1)}]\) are written as；f'_j = f'_{j-1} + \Delta t^{(j-1)} f''_{j-1} + r^2_b + r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{z}^{(j)} - r_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{c}^{(j)}
- **必要推导过程**：
  1. 对每个自由度计算触边断点 [ t_i ]。
  2. 按断点顺序推进分段路径 [ \mathbf{x}(t) ]，逐步扩大活动集。
- **与其他公式的内在联系**：
  - 输出广义 Cauchy 点与活动集，作为子空间最小化的初值与约束结构。
  - 与编号公式 [ (23) ]、[ (24) ] 在文中直接互引。

### 公式 074

- **出现位置**：`explain.md` 第 537-539 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：[ (23) ]、[ (24) ]。
- **公式完整表达**：
[
f'_j = f'_{j-1} + \Delta t^{(j-1)} f''_{j-1} + r^2_b + r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{z}^{(j)} - r_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{c}^{(j)}
]
- **公式含义**：
  - 该式属于广义 Cauchy 点构造过程，用分段一维二次模型追踪可行下降轨迹。
  - 本式位于章节 [### 3.3. Generalized Cauchy point]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{w}_b = \mathbf{W}^{\mathrm{T}}\mathbf{e}_b, \quad \mathbf{p}^{(j)} = \mathbf{W}^{\mathrm{T}}\mathbf{d}^{(j)} = \mathbf{W}^{\mathrm{T}}(\mathbf{d}^{(j-1)} + r_b \mathbf{e}_b) = \mathbf{p}^{(j-1)} + r_b \mathbf{w}_b, \quad \mathbf{c}^{(j)} = \mathbf{W}^{\mathrm{T}}\mathbf{z}^{(j)} = \mathbf{W}^{\mathrm{T}}(\mathbf{z}^{(j-1)} + \Delta t^{(j-1)}\mathbf{d}^{(j-1)}) = \mathbf{c}^{(j-1)} + \Delta t^{(j-1)}\mathbf{p}^{(j-1)}.；Plug all the above relationships into Eqs. (23) and (24), the coefficients of the quadratic polynomial \(\hat{p}(\Delta t)\) for the interval \([t^{(j)}, t^{(j+1)}]\) are written as；and；f''_j = f''_{j-1} + 2r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{d}^{(j-1)} - 2r_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{p}^{(j-1)} + r^2_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{e}_b - r^2_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{w}_b.
- **必要推导过程**：
  1. 对每个自由度计算触边断点 [ t_i ]。
  2. 按断点顺序推进分段路径 [ \mathbf{x}(t) ]，逐步扩大活动集。
  3. 在当前分段上写出一元二次模型及一、二阶导。
  4. 若解析极小点位于段内则终止，否则进入下一段。
- **与其他公式的内在联系**：
  - 输出广义 Cauchy 点与活动集，作为子空间最小化的初值与约束结构。
  - 与编号公式 [ (23) ]、[ (24) ] 在文中直接互引。

### 公式 075

- **出现位置**：`explain.md` 第 541-543 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：未显式编号。
- **公式完整表达**：
[
f''_j = f''_{j-1} + 2r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{d}^{(j-1)} - 2r_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{p}^{(j-1)} + r^2_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{e}_b - r^2_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{w}_b.
]
- **公式含义**：
  - 该式属于广义 Cauchy 点构造过程，用分段一维二次模型追踪可行下降轨迹。
  - 本式位于章节 [### 3.3. Generalized Cauchy point]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：f'_j = f'_{j-1} + \Delta t^{(j-1)} f''_{j-1} + r^2_b + r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{z}^{(j)} - r_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{c}^{(j)}；and；Similar to the steps in the previous interval, the critical point of the quadratic polynomial \(\hat{p}(\Delta t)\) is calculated as \(\Delta t^* = -f'_j / f''_j\). If \(\Delta t^* \in [0, t^{(j+1)} - t^{(j)})\), the generalized Cauchy point \(\mathbf{x}^c = \mathbf{x}(t^{(j)} + \Delta t^*)\). If \(\Delta t^* < 0\), the generalized Cauchy point \(\mathbf{x}^c = \mathbf{x}(t^{(j)})\). In the above two scenarios, the set of active constraints include those components whose break point \(t\)-value is no larger than \(t^{(j)}\), that is,；\mathcal{A}(\mathbf{x}^c) = \{i : t_i \leq t^{(j)}, i = 1, \ldots, n\}.
- **必要推导过程**：
  1. 在当前分段上写出一元二次模型及一、二阶导。
  2. 若解析极小点位于段内则终止，否则进入下一段。
- **与其他公式的内在联系**：
  - 输出广义 Cauchy 点与活动集，作为子空间最小化的初值与约束结构。

### 公式 076

- **出现位置**：`explain.md` 第 546-548 行；章节：[ ### 3.3. Generalized Cauchy point ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathcal{A}(\mathbf{x}^c) = \{i : t_i \leq t^{(j)}, i = 1, \ldots, n\}.
]
- **公式含义**：
  - 该式刻画活动集（被约束卡住的自由度），是后续降维求解的关键。
  - 该式属于广义 Cauchy 点构造过程，用分段一维二次模型追踪可行下降轨迹。
  - 本式位于章节 [### 3.3. Generalized Cauchy point]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：f''_j = f''_{j-1} + 2r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{d}^{(j-1)} - 2r_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{p}^{(j-1)} + r^2_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{e}_b - r^2_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{w}_b.；Similar to the steps in the previous interval, the critical point of the quadratic polynomial \(\hat{p}(\Delta t)\) is calculated as \(\Delta t^* = -f'_j / f''_j\). If \(\Delta t^* \in [0, t^{(j+1)} - t^{(j)})\), the generalized Cauchy point \(\mathbf{x}^c = \mathbf{x}(t^{(j)} + \Delta t^*)\). If \(\Delta t^* < 0\), the generalized Cauchy point \(\mathbf{x}^c = \mathbf{x}(t^{(j)})\). In the above two scenarios, the set of active constraints include those components whose break point \(t\)-value is no larger than \(t^{(j)}\), that is,；If \(\Delta t^* \geq t^{(j+1)} - t^{(j)}\), the search of the generalized Cauchy point moves onto the next interval (line segment) of the projected gradient. The above processes stop until the critical point \(\Delta t^*\) is smaller than the interval length formed by two neighboring break points.；**Comment 3.** The calculation of the generalized Cauchy point is extremely efficient. Recall that \(\mathbf{B}^0 = \hat{\mathbf{K}}\) is a sparse matrix assembled from the standard finite element procedure, and \(\mathbf{e}_b\) is the \(b\) th unit vector. Assume the \(b\) th row of the sparse matrix \(\mathbf{B}^0\) has \(n_b\) non-zero entries, where \(n_b\) is way smaller than the degrees of freedom \(n\) of the discretized system. When only the multiplication is counted, the operation counts of \(\mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{d}^{(j-1)}\) and \(\mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{z}^{(j)}\) are both \(\mathcal{O}(n_b)\), and the operation count of \(\mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{e}_b\) is \(\mathcal{O}(1)\). Also, since \(\mathbf{M} \in \mathbb{R}^{2m\times 2m}\), and \(\mathbf{w}_b, \mathbf{c}^{(j)}, \mathbf{p}^{(j-1)} \in \mathbb{R}^{2m}\), the operation counts of \(\mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{c}^{(j)}\), \(\mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{p}^{(j-1)}\), and \(\mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{w}_b\) are all \(\mathcal{O}(4m^2)\), where \(m\) is the small number of vector-pairs stored for the compact representation of the L-BFGS matrix.
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
  3. 根据变量是否触碰上下界识别活动集。
  4. 固定活动变量后，仅在自由子空间继续优化。
  5. 对每个自由度计算触边断点 [ t_i ]。
  6. 按断点顺序推进分段路径 [ \mathbf{x}(t) ]，逐步扩大活动集。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 输出广义 Cauchy 点与活动集，作为子空间最小化的初值与约束结构。

### 公式 077

- **出现位置**：`explain.md` 第 556-558 行；章节：[ ### 3.4. Subspace minimization ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathcal{A}(\mathbf{x}^c) = \{i : x^c_i = \mathrm{lb}_i\} \cup \{i : x^c_i = \mathrm{ub}_i\} = \{i : t_i \leq t^{(j)}, i = 1, \ldots, n\}.
]
- **公式含义**：
  - 该式刻画活动集（被约束卡住的自由度），是后续降维求解的关键。
  - 该式属于广义 Cauchy 点构造过程，用分段一维二次模型追踪可行下降轨迹。
  - 本式位于章节 [### 3.4. Subspace minimization]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：### 3.4. Subspace minimization；Once the generalized Cauchy point \(\mathbf{x}^c\) is located inside the interval \([\mathbf{x}(t^{(j)}), \mathbf{x}(t^{(j+1)})]\) on the piecewise linear path formed by the projected gradient, the active set of the box constraints is known, that is,；The components of the solution vector \(\mathbf{x}\) inside the active set are also known, each of which takes either the lower bound or the upper bound value,；x_i = x^c_i = \begin{cases} \mathrm{ub}_i & \text{if } r_i < 0 \\ \mathrm{lb}_i & \text{if } r_i > 0 \end{cases}, \quad \forall i \in \mathcal{A}(\mathbf{x}^c).
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
  3. 根据变量是否触碰上下界识别活动集。
  4. 固定活动变量后，仅在自由子空间继续优化。
  5. 对每个自由度计算触边断点 [ t_i ]。
  6. 按断点顺序推进分段路径 [ \mathbf{x}(t) ]，逐步扩大活动集。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 输出广义 Cauchy 点与活动集，作为子空间最小化的初值与约束结构。

### 公式 078

- **出现位置**：`explain.md` 第 561-563 行；章节：[ ### 3.4. Subspace minimization ]；文中编号：未显式编号。
- **公式完整表达**：
[
x_i = x^c_i = \begin{cases} \mathrm{ub}_i & \text{if } r_i < 0 \\ \mathrm{lb}_i & \text{if } r_i > 0 \end{cases}, \quad \forall i \in \mathcal{A}(\mathbf{x}^c).
]
- **公式含义**：
  - 该式刻画活动集（被约束卡住的自由度），是后续降维求解的关键。
  - 本式位于章节 [### 3.4. Subspace minimization]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathcal{A}(\mathbf{x}^c) = \{i : x^c_i = \mathrm{lb}_i\} \cup \{i : x^c_i = \mathrm{ub}_i\} = \{i : t_i \leq t^{(j)}, i = 1, \ldots, n\}.；The components of the solution vector \(\mathbf{x}\) inside the active set are also known, each of which takes either the lower bound or the upper bound value,；Assume that the size of the active set \(\mathcal{A}(\mathbf{x}^c)\) is \(q\), that is, there are totally \(q\) components of the generalized Cauchy point \(\mathbf{x}^c\) located at the boundary of the box constraints. These \(q\) components are fixed in the subsequent process. Assume that the \(i\) th component of the vector \(\mathbf{x}^c\) is fixed (at the boundary of the box constraints), let \(\mathbf{e}_i\) represent the corresponding unit vector. Then, a matrix \(\mathbf{Q} = \{\mathbf{e}_i\} \in \mathbb{R}^{n\times q}\) can be defined such that its columns are these unit vectors spanning the subspace of the fixed variables at \(\mathbf{x}^c\). Similarly, assume that the \(j\) th component of the vector \(\mathbf{x}^c\) is located inside the box constraints (the corresponding constraint is inactive), let \(\mathbf{e}_j\) represent the corresponding unit vector. Then, a matrix \(\mathbf{Z} \in \mathbb{R}^{n\times (n-q)}\) can be defined such that its columns are these unit vectors spanning the subspace of the free variables at \(\mathbf{x}^c\).；During the \(k\) th L-BFGS iteration, since at the generalized Cauchy point \(\mathbf{x}^c\), the variables contained in the active set of the box constraints \(\mathcal{A}(\mathbf{x}^c)\) are known, the goal is to solve the free variables through the minimization while keeping the components inside the active set of the box constraints fixed, that is,
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
  3. 根据变量是否触碰上下界识别活动集。
  4. 固定活动变量后，仅在自由子空间继续优化。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。

### 公式 079

- **出现位置**：`explain.md` 第 569-571 行；章节：[ ### 3.4. Subspace minimization ]；文中编号：[ (25) ]。
- **公式完整表达**：
[
\mathbf{x}^* = \arg \min m_k(\mathbf{x}) \quad (25)
]
- **公式含义**：
  - 该式把当前步未知量写成优化问题的极小点，明确“求解 = 约束最小化”。
  - 本式位于章节 [### 3.4. Subspace minimization]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：Assume that the size of the active set \(\mathcal{A}(\mathbf{x}^c)\) is \(q\), that is, there are totally \(q\) components of the generalized Cauchy point \(\mathbf{x}^c\) located at the boundary of the box constraints. These \(q\) components are fixed in the subsequent process. Assume that the \(i\) th component of the vector \(\mathbf{x}^c\) is fixed (at the boundary of the box constraints), let \(\mathbf{e}_i\) represent the corresponding unit vector. Then, a matrix \(\mathbf{Q} = \{\mathbf{e}_i\} \in \mathbb{R}^{n\times q}\) can be defined such that its columns are these unit vectors spanning the subspace of the fixed variables at \(\mathbf{x}^c\). Similarly, assume that the \(j\) th component of the vector \(\mathbf{x}^c\) is located inside the box constraints (the corresponding constraint is inactive), let \(\mathbf{e}_j\) represent the corresponding unit vector. Then, a matrix \(\mathbf{Z} \in \mathbb{R}^{n\times (n-q)}\) can be defined such that its columns are these unit vectors spanning the subspace of the free variables at \(\mathbf{x}^c\).；During the \(k\) th L-BFGS iteration, since at the generalized Cauchy point \(\mathbf{x}^c\), the variables contained in the active set of the box constraints \(\mathcal{A}(\mathbf{x}^c)\) are known, the goal is to solve the free variables through the minimization while keeping the components inside the active set of the box constraints fixed, that is,；subject to；x_i = x^c_i, \quad \forall i \in \mathcal{A}(\mathbf{x}^c)
- **必要推导过程**：
  1. 把增量步未知量组合为优化变量向量。
  2. 施加不可逆约束后，写成带盒约束的极小化问题。
- **与其他公式的内在联系**：
  - 与编号公式 [ (25) ] 在文中直接互引。

### 公式 080

- **出现位置**：`explain.md` 第 575-577 行；章节：[ ### 3.4. Subspace minimization ]；文中编号：未显式编号。
- **公式完整表达**：
[
x_i = x^c_i, \quad \forall i \in \mathcal{A}(\mathbf{x}^c)
]
- **公式含义**：
  - 该式刻画活动集（被约束卡住的自由度），是后续降维求解的关键。
  - 本式位于章节 [### 3.4. Subspace minimization]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{x}^* = \arg \min m_k(\mathbf{x}) \quad (25)；subject to；and；\mathrm{lb}_i \leq x_i \leq \mathrm{ub}_i, \quad \forall i \notin \mathcal{A}(\mathbf{x}^c). \quad (26)
- **必要推导过程**：
  1. 根据变量是否触碰上下界识别活动集。
  2. 固定活动变量后，仅在自由子空间继续优化。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。

### 公式 081

- **出现位置**：`explain.md` 第 579-581 行；章节：[ ### 3.4. Subspace minimization ]；文中编号：[ (26) ]、[ (15) ]。
- **公式完整表达**：
[
\mathrm{lb}_i \leq x_i \leq \mathrm{ub}_i, \quad \forall i \notin \mathcal{A}(\mathbf{x}^c). \quad (26)
]
- **公式含义**：
  - 该式刻画活动集（被约束卡住的自由度），是后续降维求解的关键。
  - 本式位于章节 [### 3.4. Subspace minimization]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：x_i = x^c_i, \quad \forall i \in \mathcal{A}(\mathbf{x}^c)；and；Recall that the objective function \(m_k(\mathbf{x})\) is the quadratic function defined in Eq. (15). There are generally two approaches to solve the above minimization, the primal approach that directly solves the free variables through the subspace minimization, or the dual approach that enforces the fixed variables via the Lagrange multipliers. Sections 3.4.1 and 3.4.2 present the primal approach via a direct matrix factorization and a conjugate gradient method, respectively. Section 3.4.3 lays out the Schur complement method for the dual approach and discusses the difficulties of using this method to solve the current problem.；#### 3.4.1. Direct matrix factorization for the primal approach
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
  3. 根据变量是否触碰上下界识别活动集。
  4. 固定活动变量后，仅在自由子空间继续优化。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 与编号公式 [ (26) ]、[ (15) ] 在文中直接互引。

### 公式 082

- **出现位置**：`explain.md` 第 589-591 行；章节：[ #### 3.4.1. Direct matrix factorization for the primal approach ]；文中编号：[ (15) ]。
- **公式完整表达**：
[
\mathbf{x} = \mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}},
]
- **公式含义**：
  - 该式属于原始子空间法：固定活动变量，仅在自由变量子空间解线性系统。
  - 本式位于章节 [#### 3.4.1. Direct matrix factorization for the primal approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：#### 3.4.1. Direct matrix factorization for the primal approach；At the \(k\) th L-BFGS iteration, let \(\mathbf{x}^c \in \mathbb{R}^n\) represent the generalized Cauchy point, and \(\mathbf{Z}_k \in \mathbb{R}^{n\times (n-q)}\) represent the matrix whose columns span the subspace of the free variables at \(\mathbf{x}^c\). Assume that \(q\) components of the generalized Cauchy point \(\mathbf{x}^c\) are at the boundary of the box constraints. Then, the unknown vector \(\mathbf{x}\) can be expressed in the following form,；where \(\hat{\mathbf{x}} \in \mathbb{R}^{n-q}\) is the vector of the free variables. The quadratic model in Eq. (15) can be transformed as；\begin{array}{rl}
- **必要推导过程**：
  1. 构造自由子空间基 [ \mathbf{Z}_k ] 并降维。
  2. 解 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ] 后回代到全空间。
- **与其他公式的内在联系**：
  - 产生最终步长方向 [ \mathbf{p}_k ]，再进入线搜索更新 [ \mathbf{x}_{k+1} ]。
  - 与编号公式 [ (15) ] 在文中直接互引。

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
  - 该式定义残量向量，残量为零对应离散平衡/最优性条件成立。
  - 该式属于原始子空间法：固定活动变量，仅在自由变量子空间解线性系统。
  - 本式位于章节 [#### 3.4.1. Direct matrix factorization for the primal approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{x} = \mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}},；where \(\hat{\mathbf{x}} \in \mathbb{R}^{n-q}\) is the vector of the free variables. The quadratic model in Eq. (15) can be transformed as；Since the L-BFGS matrix \(\mathbf{B}_k\) is positive definite, the above quadratic function \(\hat{m}_k(\hat{\mathbf{x}})\) of the free variables \(\hat{\mathbf{x}}\) reaches to its minimum at the critical point,；\hat{m}'_k(\hat{\mathbf{x}}) = 0 \implies \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k \hat{\mathbf{x}} = -\mathbf{Z}^{\mathrm{T}}_k [\mathbf{r}_k + \mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k)].
- **必要推导过程**：
  1. 对离散势能对节点自由度求偏导得到残量。
  2. 残量组装后形成非线性代数方程组。
  3. 构造自由子空间基 [ \mathbf{Z}_k ] 并降维。
  4. 解 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ] 后回代到全空间。
- **与其他公式的内在联系**：
  - 把连续模型桥接到离散代数系统，直接作为第 3 节优化器输入。
  - 产生最终步长方向 [ \mathbf{p}_k ]，再进入线搜索更新 [ \mathbf{x}_{k+1} ]。
  - 与编号公式 [ (15) ] 在文中直接互引。

### 公式 084

- **出现位置**：`explain.md` 第 608-610 行；章节：[ #### 3.4.1. Direct matrix factorization for the primal approach ]；文中编号：未显式编号。
- **公式完整表达**：
[
\hat{m}'_k(\hat{\mathbf{x}}) = 0 \implies \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k \hat{\mathbf{x}} = -\mathbf{Z}^{\mathrm{T}}_k [\mathbf{r}_k + \mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k)].
]
- **公式含义**：
  - 该式定义残量向量，残量为零对应离散平衡/最优性条件成立。
  - 该式属于原始子空间法：固定活动变量，仅在自由变量子空间解线性系统。
  - 本式位于章节 [#### 3.4.1. Direct matrix factorization for the primal approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\end{array}；Since the L-BFGS matrix \(\mathbf{B}_k\) is positive definite, the above quadratic function \(\hat{m}_k(\hat{\mathbf{x}})\) of the free variables \(\hat{\mathbf{x}}\) reaches to its minimum at the critical point,；Let \(\hat{\mathbf{B}}_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k\) and \(\hat{\mathbf{r}}_k = \mathbf{Z}^{\mathrm{T}}_k [\mathbf{r}_k + \mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k)]\), the reduced linear system of the subspace minimization becomes；\hat{\mathbf{B}}_k \hat{\mathbf{x}} = -\hat{\mathbf{r}}_k. \quad (27)
- **必要推导过程**：
  1. 对离散势能对节点自由度求偏导得到残量。
  2. 残量组装后形成非线性代数方程组。
  3. 构造自由子空间基 [ \mathbf{Z}_k ] 并降维。
  4. 解 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ] 后回代到全空间。
- **与其他公式的内在联系**：
  - 把连续模型桥接到离散代数系统，直接作为第 3 节优化器输入。
  - 产生最终步长方向 [ \mathbf{p}_k ]，再进入线搜索更新 [ \mathbf{x}_{k+1} ]。

### 公式 085

- **出现位置**：`explain.md` 第 614-616 行；章节：[ #### 3.4.1. Direct matrix factorization for the primal approach ]；文中编号：[ (27) ]。
- **公式完整表达**：
[
\hat{\mathbf{B}}_k \hat{\mathbf{x}} = -\hat{\mathbf{r}}_k. \quad (27)
]
- **公式含义**：
  - 该式定义残量向量，残量为零对应离散平衡/最优性条件成立。
  - 该式属于原始子空间法：固定活动变量，仅在自由变量子空间解线性系统。
  - 本式位于章节 [#### 3.4.1. Direct matrix factorization for the primal approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\hat{m}'_k(\hat{\mathbf{x}}) = 0 \implies \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k \hat{\mathbf{x}} = -\mathbf{Z}^{\mathrm{T}}_k [\mathbf{r}_k + \mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k)].；Let \(\hat{\mathbf{B}}_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k\) and \(\hat{\mathbf{r}}_k = \mathbf{Z}^{\mathrm{T}}_k [\mathbf{r}_k + \mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k)]\), the reduced linear system of the subspace minimization becomes；Recall that \(\mathbf{B}_k\) is not directly available, since it is not stored in the memory component wise. Rather, \(\mathbf{B}_k\) is only available in the compact representation form as discussed in Section 3.2, that is,；\mathbf{B}_k = \mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k.
- **必要推导过程**：
  1. 对离散势能对节点自由度求偏导得到残量。
  2. 残量组装后形成非线性代数方程组。
  3. 构造自由子空间基 [ \mathbf{Z}_k ] 并降维。
  4. 解 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ] 后回代到全空间。
- **与其他公式的内在联系**：
  - 把连续模型桥接到离散代数系统，直接作为第 3 节优化器输入。
  - 产生最终步长方向 [ \mathbf{p}_k ]，再进入线搜索更新 [ \mathbf{x}_{k+1} ]。
  - 与编号公式 [ (27) ] 在文中直接互引。

### 公式 086

- **出现位置**：`explain.md` 第 619-621 行；章节：[ #### 3.4.1. Direct matrix factorization for the primal approach ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{B}_k = \mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k.
]
- **公式含义**：
  - 该式属于 L-BFGS 曲率信息更新链，用少量历史向量近似全 Hessian。
  - 本式位于章节 [#### 3.4.1. Direct matrix factorization for the primal approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\hat{\mathbf{B}}_k \hat{\mathbf{x}} = -\hat{\mathbf{r}}_k. \quad (27)；Recall that \(\mathbf{B}_k\) is not directly available, since it is not stored in the memory component wise. Rather, \(\mathbf{B}_k\) is only available in the compact representation form as discussed in Section 3.2, that is,；Therefore,；\hat{\mathbf{B}}_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k = \mathbf{Z}^{\mathrm{T}}_k (\mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k) \mathbf{Z}_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{Z}_k - \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k. \quad (28)
- **必要推导过程**：
  1. 把多次 BFGS 更新改写为紧凑低秩形式 [ \mathbf{B}_k=\mathbf{B}^0_k-\mathbf{W}_k\mathbf{M}_k\mathbf{W}_k^T ]。
  2. 只需维护小矩阵 [ \mathbf{M}_k ] 与历史向量列块，降低存储和计算。
- **与其他公式的内在联系**：
  - 为 GCP 段内导数与子空间系统提供近似曲率信息。

### 公式 087

- **出现位置**：`explain.md` 第 623-625 行；章节：[ #### 3.4.1. Direct matrix factorization for the primal approach ]；文中编号：[ (28) ]、[ (27) ]。
- **公式完整表达**：
[
\hat{\mathbf{B}}_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k = \mathbf{Z}^{\mathrm{T}}_k (\mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k) \mathbf{Z}_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{Z}_k - \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k. \quad (28)
]
- **公式含义**：
  - 该式属于 L-BFGS 曲率信息更新链，用少量历史向量近似全 Hessian。
  - 该式属于原始子空间法：固定活动变量，仅在自由变量子空间解线性系统。
  - 本式位于章节 [#### 3.4.1. Direct matrix factorization for the primal approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{B}_k = \mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k.；Therefore,；Similar to the L-BFGS matrix \(\mathbf{B}_k\), the reduced matrix \(\hat{\mathbf{B}}_k\) is not stored in the memory component wise. As a result, the direct solvers based on the matrix factorization, such as the LU decomposition and the Cholesky decomposition, cannot be directly used to solve the reduced system in Eq. (27). However, \(\mathbf{B}^0_k\) is sparse (recall that it is assembled from the finite element procedure). Therefore, the reduced matrix \(\hat{\mathbf{B}}^0_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{Z}_k\) is also sparse and can be factorized using these techniques.；Using the Sherman–Morrison–Woodbury formula [45], as stated in the Appendix, the vector of the free variables can be solved as
- **必要推导过程**：
  1. 把多次 BFGS 更新改写为紧凑低秩形式 [ \mathbf{B}_k=\mathbf{B}^0_k-\mathbf{W}_k\mathbf{M}_k\mathbf{W}_k^T ]。
  2. 只需维护小矩阵 [ \mathbf{M}_k ] 与历史向量列块，降低存储和计算。
  3. 构造自由子空间基 [ \mathbf{Z}_k ] 并降维。
  4. 解 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ] 后回代到全空间。
- **与其他公式的内在联系**：
  - 为 GCP 段内导数与子空间系统提供近似曲率信息。
  - 产生最终步长方向 [ \mathbf{p}_k ]，再进入线搜索更新 [ \mathbf{x}_{k+1} ]。
  - 与编号公式 [ (28) ]、[ (27) ] 在文中直接互引。

### 公式 088

- **出现位置**：`explain.md` 第 631-633 行；章节：[ #### 3.4.1. Direct matrix factorization for the primal approach ]；文中编号：[ (29) ]。
- **公式完整表达**：
[
\hat{\mathbf{x}} = -\hat{\mathbf{B}}^{-1}_k \hat{\mathbf{r}}_k,
]
- **公式含义**：
  - 该式定义残量向量，残量为零对应离散平衡/最优性条件成立。
  - 该式属于原始子空间法：固定活动变量，仅在自由变量子空间解线性系统。
  - 本式位于章节 [#### 3.4.1. Direct matrix factorization for the primal approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：Similar to the L-BFGS matrix \(\mathbf{B}_k\), the reduced matrix \(\hat{\mathbf{B}}_k\) is not stored in the memory component wise. As a result, the direct solvers based on the matrix factorization, such as the LU decomposition and the Cholesky decomposition, cannot be directly used to solve the reduced system in Eq. (27). However, \(\mathbf{B}^0_k\) is sparse (recall that it is assembled from the finite element procedure). Therefore, the reduced matrix \(\hat{\mathbf{B}}^0_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{Z}_k\) is also sparse and can be factorized using these techniques.；Using the Sherman–Morrison–Woodbury formula [45], as stated in the Appendix, the vector of the free variables can be solved as；\hat{\mathbf{B}}^{-1}_k = (\hat{\mathbf{B}}^0_k)^{-1} + (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \left[\mathbf{I} - \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \right]^{-1} \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1}, \quad (29)；where
- **必要推导过程**：
  1. 对离散势能对节点自由度求偏导得到残量。
  2. 残量组装后形成非线性代数方程组。
  3. 构造自由子空间基 [ \mathbf{Z}_k ] 并降维。
  4. 解 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ] 后回代到全空间。
- **与其他公式的内在联系**：
  - 把连续模型桥接到离散代数系统，直接作为第 3 节优化器输入。
  - 产生最终步长方向 [ \mathbf{p}_k ]，再进入线搜索更新 [ \mathbf{x}_{k+1} ]。
  - 与编号公式 [ (29) ] 在文中直接互引。

### 公式 089

- **出现位置**：`explain.md` 第 634-636 行；章节：[ #### 3.4.1. Direct matrix factorization for the primal approach ]；文中编号：[ (29) ]。
- **公式完整表达**：
[
\hat{\mathbf{B}}^{-1}_k = (\hat{\mathbf{B}}^0_k)^{-1} + (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \left[\mathbf{I} - \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \right]^{-1} \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1}, \quad (29)
]
- **公式含义**：
  - 该式属于 L-BFGS 曲率信息更新链，用少量历史向量近似全 Hessian。
  - 该式属于原始子空间法：固定活动变量，仅在自由变量子空间解线性系统。
  - 本式位于章节 [#### 3.4.1. Direct matrix factorization for the primal approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：Using the Sherman–Morrison–Woodbury formula [45], as stated in the Appendix, the vector of the free variables can be solved as；\hat{\mathbf{x}} = -\hat{\mathbf{B}}^{-1}_k \hat{\mathbf{r}}_k,；where；\mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k = \mathbf{Z}^{\mathrm{T}}_k [\mathbf{Y}_k \quad \mathbf{B}^0_k \mathbf{S}_k] = [\mathbf{Z}^{\mathrm{T}}_k \mathbf{y}_{k-m} \cdots \mathbf{Z}^{\mathrm{T}}_k \mathbf{y}_{k-1} \quad \mathbf{Z}^{\mathrm{T}}_k(\mathbf{B}^0_k \mathbf{s}_{k-m}) \cdots \mathbf{Z}^{\mathrm{T}}_k(\mathbf{B}^0_k \mathbf{s}_{k-1})].
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
  3. 把多次 BFGS 更新改写为紧凑低秩形式 [ \mathbf{B}_k=\mathbf{B}^0_k-\mathbf{W}_k\mathbf{M}_k\mathbf{W}_k^T ]。
  4. 只需维护小矩阵 [ \mathbf{M}_k ] 与历史向量列块，降低存储和计算。
  5. 构造自由子空间基 [ \mathbf{Z}_k ] 并降维。
  6. 解 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ] 后回代到全空间。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 为 GCP 段内导数与子空间系统提供近似曲率信息。
  - 产生最终步长方向 [ \mathbf{p}_k ]，再进入线搜索更新 [ \mathbf{x}_{k+1} ]。
  - 与编号公式 [ (29) ] 在文中直接互引。

### 公式 090

- **出现位置**：`explain.md` 第 639-641 行；章节：[ #### 3.4.1. Direct matrix factorization for the primal approach ]；文中编号：[ (29) ]。
- **公式完整表达**：
[
\mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k = \mathbf{Z}^{\mathrm{T}}_k [\mathbf{Y}_k \quad \mathbf{B}^0_k \mathbf{S}_k] = [\mathbf{Z}^{\mathrm{T}}_k \mathbf{y}_{k-m} \cdots \mathbf{Z}^{\mathrm{T}}_k \mathbf{y}_{k-1} \quad \mathbf{Z}^{\mathrm{T}}_k(\mathbf{B}^0_k \mathbf{s}_{k-m}) \cdots \mathbf{Z}^{\mathrm{T}}_k(\mathbf{B}^0_k \mathbf{s}_{k-1})].
]
- **公式含义**：
  - 该式属于 L-BFGS 曲率信息更新链，用少量历史向量近似全 Hessian。
  - 本式位于章节 [#### 3.4.1. Direct matrix factorization for the primal approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\hat{\mathbf{B}}^{-1}_k = (\hat{\mathbf{B}}^0_k)^{-1} + (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \left[\mathbf{I} - \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \right]^{-1} \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1}, \quad (29)；where；In order to calculate the inverse of \(\hat{\mathbf{B}}_k\), the inverse of \(\hat{\mathbf{B}}^0_k\) needs to be repeatedly applied to the column vectors in \(\mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k\). This can be achieved through the sparse LU decomposition of \(\hat{\mathbf{B}}^0_k\). Also, notice that in Eq. (29), the following term；\mathbf{I} - \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \in \mathbb{R}^{2m\times 2m}
- **必要推导过程**：
  1. 把多次 BFGS 更新改写为紧凑低秩形式 [ \mathbf{B}_k=\mathbf{B}^0_k-\mathbf{W}_k\mathbf{M}_k\mathbf{W}_k^T ]。
  2. 只需维护小矩阵 [ \mathbf{M}_k ] 与历史向量列块，降低存储和计算。
- **与其他公式的内在联系**：
  - 为 GCP 段内导数与子空间系统提供近似曲率信息。
  - 与编号公式 [ (29) ] 在文中直接互引。

### 公式 091

- **出现位置**：`explain.md` 第 644-646 行；章节：[ #### 3.4.1. Direct matrix factorization for the primal approach ]；文中编号：[ (29) ]。
- **公式完整表达**：
[
\mathbf{I} - \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \in \mathbb{R}^{2m\times 2m}
]
- **公式含义**：
  - 该式属于 L-BFGS 曲率信息更新链，用少量历史向量近似全 Hessian。
  - 该式属于原始子空间法：固定活动变量，仅在自由变量子空间解线性系统。
  - 本式位于章节 [#### 3.4.1. Direct matrix factorization for the primal approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k = \mathbf{Z}^{\mathrm{T}}_k [\mathbf{Y}_k \quad \mathbf{B}^0_k \mathbf{S}_k] = [\mathbf{Z}^{\mathrm{T}}_k \mathbf{y}_{k-m} \cdots \mathbf{Z}^{\mathrm{T}}_k \mathbf{y}_{k-1} \quad \mathbf{Z}^{\mathrm{T}}_k(\mathbf{B}^0_k \mathbf{s}_{k-m}) \cdots \mathbf{Z}^{\mathrm{T}}_k(\mathbf{B}^0_k \mathbf{s}_{k-1})].；In order to calculate the inverse of \(\hat{\mathbf{B}}_k\), the inverse of \(\hat{\mathbf{B}}^0_k\) needs to be repeatedly applied to the column vectors in \(\mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k\). This can be achieved through the sparse LU decomposition of \(\hat{\mathbf{B}}^0_k\). Also, notice that in Eq. (29), the following term；is a \(2m\times 2m\) matrix, therefore, its inverse can be directly calculated with a negligible computational cost.；#### 3.4.2. Conjugate gradient method for the primal approach
- **必要推导过程**：
  1. 把多次 BFGS 更新改写为紧凑低秩形式 [ \mathbf{B}_k=\mathbf{B}^0_k-\mathbf{W}_k\mathbf{M}_k\mathbf{W}_k^T ]。
  2. 只需维护小矩阵 [ \mathbf{M}_k ] 与历史向量列块，降低存储和计算。
  3. 构造自由子空间基 [ \mathbf{Z}_k ] 并降维。
  4. 解 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ] 后回代到全空间。
- **与其他公式的内在联系**：
  - 为 GCP 段内导数与子空间系统提供近似曲率信息。
  - 产生最终步长方向 [ \mathbf{p}_k ]，再进入线搜索更新 [ \mathbf{x}_{k+1} ]。
  - 与编号公式 [ (29) ] 在文中直接互引。

### 公式 092

- **出现位置**：`explain.md` 第 655-657 行；章节：[ #### 3.4.2. Conjugate gradient method for the primal approach ]；文中编号：[ (18) ]。
- **公式完整表达**：
[
\mathbf{P}^{-1} = \mathrm{ILU}(\hat{\mathbf{B}}^0_k) = \mathrm{ILU}(\mathbf{Z}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{Z}_k).
]
- **公式含义**：
  - 该式属于 L-BFGS 曲率信息更新链，用少量历史向量近似全 Hessian。
  - 该式属于原始子空间法：固定活动变量，仅在自由变量子空间解线性系统。
  - 本式位于章节 [#### 3.4.2. Conjugate gradient method for the primal approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：As mentioned before, \(\mathbf{B}_k\) is not stored in the memory component wise and is only available in the compact representation form shown in Eq. (20). Therefore, the reduced matrix \(\hat{\mathbf{B}}_k\) is not directly available and only has the compact form shown in Eq. (28). Recall that \(\mathbf{B}_k\) is positive definite as long as the vector-pairs \(\{\mathbf{s}_i, \mathbf{y}_i\}_{i=k-m}^{k-1}\) satisfy the curvature condition Eq. (18). Therefore, \(\hat{\mathbf{B}}_k\) is also positive definite, making the conjugate gradient method an ideal candidate to solve the reduced linear system Eq. (27). Moreover, as a type of iterative solver, the conjugate gradient method does not require the component form of the reduced matrix \(\hat{\mathbf{B}}_k\). Rather, it only needs to know how to perform the matrix–vector multiplication between the matrix \(\hat{\mathbf{B}}_k\) and an arbitrary vector. Indeed, this matrix–vector multiplication is well defined according to Eq. (28), in which the matrix form of \(\hat{\mathbf{B}}^0_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{Z}_k\), \(\mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k\), and \(\mathbf{M}_k\) are all available. Therefore, besides the direct matrix factorization approach provided in Section 3.4.1, the conjugate gradient method provides another option to solve the reduced linear system in order to obtain the free variables \(\hat{\mathbf{x}}\) at the generalized Cauchy point \(\mathbf{x}^c\).；It is well known that an effective preconditioner could significantly improve the performance of an iterative linear solver such as the conjugate gradient method. A good preconditioner should balance between the computational cost incurred for its own construction and the performance boost provided to the iterative solver. Let \(\mathbf{P}^{-1}\) represent the preconditioner. Considering the similarity between \(\hat{\mathbf{B}}^0_k\) and \(\hat{\mathbf{B}}_k\), in this work, the incomplete LU decomposition of \(\hat{\mathbf{B}}^0_k\) is chosen as the preconditioner for the conjugate gradient method to solve the reduced linear system, that is,；**Comment 4.** The reduced matrix \(\hat{\mathbf{B}}^0_k\) is positive definite since the curvature condition in Eq. (18) is satisfied due to the line search based on the strong Wolfe conditions. Therefore, the Cholesky decomposition is more efficient than the LU decomposition in the direct matrix factorization approach. Similarly, in the conjugate gradient method, the preconditioner based on the incomplete Cholesky decomposition is less expensive to construct compared with the counterpart based on the incomplete LU decomposition. However, since the UMFPACK used in the deal.II library is built on the direct sparse LU factorization, and the incomplete LU preconditioner is also readily available, these techniques are chosen here mainly for implementation convenience.；#### 3.4.3. Schur complement for the dual approach
- **必要推导过程**：
  1. 把多次 BFGS 更新改写为紧凑低秩形式 [ \mathbf{B}_k=\mathbf{B}^0_k-\mathbf{W}_k\mathbf{M}_k\mathbf{W}_k^T ]。
  2. 只需维护小矩阵 [ \mathbf{M}_k ] 与历史向量列块，降低存储和计算。
  3. 构造自由子空间基 [ \mathbf{Z}_k ] 并降维。
  4. 解 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ] 后回代到全空间。
- **与其他公式的内在联系**：
  - 为 GCP 段内导数与子空间系统提供近似曲率信息。
  - 产生最终步长方向 [ \mathbf{p}_k ]，再进入线搜索更新 [ \mathbf{x}_{k+1} ]。
  - 与编号公式 [ (18) ] 在文中直接互引。

### 公式 093

- **出现位置**：`explain.md` 第 665-667 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：[ (25) ]、[ (26) ]、[ (15) ]。
- **公式完整表达**：
[
\mathbf{x} = \mathbf{x}_k + \Delta \mathbf{x}_k.
]
- **公式含义**：
  - 该式在当前章节承担定义/变换/约束中的一个局部环节，用于连接前后公式。
  - 本式位于章节 [#### 3.4.3. Schur complement for the dual approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：#### 3.4.3. Schur complement for the dual approach；Instead of formulating the minimization problem in Eq. (25) subject to the constraints in Eq. (26) as a subspace minimization, an alternative approach is to formulate the above problem via a dual approach using the Lagrange multipliers. Recall that \(\mathbf{x}_k\) is known at the beginning of the current iteration, and the generalized Cauchy point \(\mathbf{x}^c\) is already located on the piecewise linear path formed by the projected gradient. The unknown vector \(\mathbf{x}\) is expressed as；Plug the above equation into the quadratic model shown in Eq. (15),；m_k(\mathbf{x}) = \Pi_k + \mathbf{r}_k^{\mathrm{T}} \Delta \mathbf{x}_k + \frac{1}{2} \Delta \mathbf{x}_k^{\mathrm{T}} \mathbf{B}_k \Delta \mathbf{x}_k.
- **必要推导过程**：
  1. 根据本式中的等号关系执行代入、求导或约束投影，可从上一式直接得到当前式。
- **与其他公式的内在联系**：
  - 与编号公式 [ (25) ]、[ (26) ]、[ (15) ] 在文中直接互引。

### 公式 094

- **出现位置**：`explain.md` 第 671-673 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：[ (15) ]。
- **公式完整表达**：
[
m_k(\mathbf{x}) = \Pi_k + \mathbf{r}_k^{\mathrm{T}} \Delta \mathbf{x}_k + \frac{1}{2} \Delta \mathbf{x}_k^{\mathrm{T}} \mathbf{B}_k \Delta \mathbf{x}_k.
]
- **公式含义**：
  - 该式定义残量向量，残量为零对应离散平衡/最优性条件成立。
  - 本式位于章节 [#### 3.4.3. Schur complement for the dual approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{x} = \mathbf{x}_k + \Delta \mathbf{x}_k.；Plug the above equation into the quadratic model shown in Eq. (15),；Also, at the generalized Cauchy point \(\mathbf{x}^c\), for the variables in the active set of the box constraints,；x_i = x^c_i, \quad \forall i \in \mathcal{A}_k(\mathbf{x}^c),
- **必要推导过程**：
  1. 对离散势能对节点自由度求偏导得到残量。
  2. 残量组装后形成非线性代数方程组。
- **与其他公式的内在联系**：
  - 把连续模型桥接到离散代数系统，直接作为第 3 节优化器输入。
  - 与编号公式 [ (15) ] 在文中直接互引。

### 公式 095

- **出现位置**：`explain.md` 第 676-678 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：未显式编号。
- **公式完整表达**：
[
x_i = x^c_i, \quad \forall i \in \mathcal{A}_k(\mathbf{x}^c),
]
- **公式含义**：
  - 该式刻画活动集（被约束卡住的自由度），是后续降维求解的关键。
  - 本式位于章节 [#### 3.4.3. Schur complement for the dual approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：m_k(\mathbf{x}) = \Pi_k + \mathbf{r}_k^{\mathrm{T}} \Delta \mathbf{x}_k + \frac{1}{2} \Delta \mathbf{x}_k^{\mathrm{T}} \mathbf{B}_k \Delta \mathbf{x}_k.；Also, at the generalized Cauchy point \(\mathbf{x}^c\), for the variables in the active set of the box constraints,；where \(x^c_i\) takes the upper bound value \(\mathrm{ub}_i\) or the lower bound value \(\mathrm{lb}_i\) of the box constraints. Recall that the matrix \(\mathbf{Q}_k = \{\mathbf{e}_i\} \in \mathbb{R}^{n\times q}\) is defined such that its columns are the unit vectors spanning the subspace of the fixed variables at \(\mathbf{x}^c\). Therefore,；\mathbf{Q}^{\mathrm{T}}_k \Delta \mathbf{x}_k = \mathbf{Q}^{\mathrm{T}}_k (\mathbf{x} - \mathbf{x}_k) = \mathbf{Q}^{\mathrm{T}}_k (\mathbf{x}^c - \mathbf{x}_k).
- **必要推导过程**：
  1. 根据变量是否触碰上下界识别活动集。
  2. 固定活动变量后，仅在自由子空间继续优化。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。

### 公式 096

- **出现位置**：`explain.md` 第 680-682 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：[ (25) ]、[ (26) ]。
- **公式完整表达**：
[
\mathbf{Q}^{\mathrm{T}}_k \Delta \mathbf{x}_k = \mathbf{Q}^{\mathrm{T}}_k (\mathbf{x} - \mathbf{x}_k) = \mathbf{Q}^{\mathrm{T}}_k (\mathbf{x}^c - \mathbf{x}_k).
]
- **公式含义**：
  - 该式在当前章节承担定义/变换/约束中的一个局部环节，用于连接前后公式。
  - 本式位于章节 [#### 3.4.3. Schur complement for the dual approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：x_i = x^c_i, \quad \forall i \in \mathcal{A}_k(\mathbf{x}^c),；where \(x^c_i\) takes the upper bound value \(\mathrm{ub}_i\) or the lower bound value \(\mathrm{lb}_i\) of the box constraints. Recall that the matrix \(\mathbf{Q}_k = \{\mathbf{e}_i\} \in \mathbb{R}^{n\times q}\) is defined such that its columns are the unit vectors spanning the subspace of the fixed variables at \(\mathbf{x}^c\). Therefore,；The minimization problem in Eq. (25) and the corresponding constraints in Eq. (26) can be rewritten as；\min \left\{ \Pi_k + \mathbf{r}_k^{\mathrm{T}} \Delta \mathbf{x}_k + \frac{1}{2} \Delta \mathbf{x}_k^{\mathrm{T}} \mathbf{B}_k \Delta \mathbf{x}_k \right\}
- **必要推导过程**：
  1. 根据本式中的等号关系执行代入、求导或约束投影，可从上一式直接得到当前式。
- **与其他公式的内在联系**：
  - 与编号公式 [ (25) ]、[ (26) ] 在文中直接互引。

### 公式 097

- **出现位置**：`explain.md` 第 685-687 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：[ (25) ]、[ (26) ]。
- **公式完整表达**：
[
\min \left\{ \Pi_k + \mathbf{r}_k^{\mathrm{T}} \Delta \mathbf{x}_k + \frac{1}{2} \Delta \mathbf{x}_k^{\mathrm{T}} \mathbf{B}_k \Delta \mathbf{x}_k \right\}
]
- **公式含义**：
  - 该式定义残量向量，残量为零对应离散平衡/最优性条件成立。
  - 本式位于章节 [#### 3.4.3. Schur complement for the dual approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{Q}^{\mathrm{T}}_k \Delta \mathbf{x}_k = \mathbf{Q}^{\mathrm{T}}_k (\mathbf{x} - \mathbf{x}_k) = \mathbf{Q}^{\mathrm{T}}_k (\mathbf{x}^c - \mathbf{x}_k).；The minimization problem in Eq. (25) and the corresponding constraints in Eq. (26) can be rewritten as；subject to；\mathbf{Q}_k^{\mathrm{T}} \Delta \mathbf{x}_k = \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k)
- **必要推导过程**：
  1. 对离散势能对节点自由度求偏导得到残量。
  2. 残量组装后形成非线性代数方程组。
  3. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  4. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 把连续模型桥接到离散代数系统，直接作为第 3 节优化器输入。
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 与编号公式 [ (25) ]、[ (26) ] 在文中直接互引。

### 公式 098

- **出现位置**：`explain.md` 第 689-691 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{Q}_k^{\mathrm{T}} \Delta \mathbf{x}_k = \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k)
]
- **公式含义**：
  - 该式属于对偶/KKT 路线：通过 Schur 补先求乘子再回代主变量。
  - 本式位于章节 [#### 3.4.3. Schur complement for the dual approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\min \left\{ \Pi_k + \mathbf{r}_k^{\mathrm{T}} \Delta \mathbf{x}_k + \frac{1}{2} \Delta \mathbf{x}_k^{\mathrm{T}} \mathbf{B}_k \Delta \mathbf{x}_k \right\}；subject to；and the box constraints；\mathbf{lb} - \mathbf{x}_k \leq \Delta \mathbf{x}_k \leq \mathbf{ub} - \mathbf{x}_k.
- **必要推导过程**：
  1. 写出 KKT 增广系统。
  2. 对主变量消元得到 Schur 补方程求乘子 [ \lambda ]，再回代 [ \Delta\mathbf{x}_k ]。
- **与其他公式的内在联系**：
  - 产生最终步长方向 [ \mathbf{p}_k ]，再进入线搜索更新 [ \mathbf{x}_{k+1} ]。

### 公式 099

- **出现位置**：`explain.md` 第 693-695 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{lb} - \mathbf{x}_k \leq \Delta \mathbf{x}_k \leq \mathbf{ub} - \mathbf{x}_k.
]
- **公式含义**：
  - 该式在当前章节承担定义/变换/约束中的一个局部环节，用于连接前后公式。
  - 本式位于章节 [#### 3.4.3. Schur complement for the dual approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{Q}_k^{\mathrm{T}} \Delta \mathbf{x}_k = \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k)；and the box constraints；For now, ignore the box constraints and let \(\lambda\) represent the Lagrange multiplier for the equality constraints. The optimality condition of the above constrained minimization is；\begin{bmatrix} \mathbf{B}_k & \mathbf{Q}_k \\ \mathbf{Q}_k^{\mathrm{T}} & \mathbf{0} \end{bmatrix} \begin{bmatrix} \Delta \mathbf{x}_k \\ \lambda \end{bmatrix} = \begin{bmatrix} -\mathbf{r}_k \\ \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k) \end{bmatrix}. \quad (30)
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。

### 公式 100

- **出现位置**：`explain.md` 第 699-701 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：[ (30) ]。
- **公式完整表达**：
[
\begin{bmatrix} \mathbf{B}_k & \mathbf{Q}_k \\ \mathbf{Q}_k^{\mathrm{T}} & \mathbf{0} \end{bmatrix} \begin{bmatrix} \Delta \mathbf{x}_k \\ \lambda \end{bmatrix} = \begin{bmatrix} -\mathbf{r}_k \\ \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k) \end{bmatrix}. \quad (30)
]
- **公式含义**：
  - 该式定义残量向量，残量为零对应离散平衡/最优性条件成立。
  - 该式属于对偶/KKT 路线：通过 Schur 补先求乘子再回代主变量。
  - 本式位于章节 [#### 3.4.3. Schur complement for the dual approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{lb} - \mathbf{x}_k \leq \Delta \mathbf{x}_k \leq \mathbf{ub} - \mathbf{x}_k.；For now, ignore the box constraints and let \(\lambda\) represent the Lagrange multiplier for the equality constraints. The optimality condition of the above constrained minimization is；A popular method to solve the above linear system is to use the Schur complement approach by firstly solving \(\lambda\) from；(\mathbf{Q}_k^{\mathrm{T}} \mathbf{B}_k^{-1} \mathbf{Q}_k) \lambda = -\mathbf{Q}_k^{\mathrm{T}} \mathbf{B}_k^{-1} \mathbf{r}_k - \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k). \quad (31)
- **必要推导过程**：
  1. 对离散势能对节点自由度求偏导得到残量。
  2. 残量组装后形成非线性代数方程组。
  3. 写出 KKT 增广系统。
  4. 对主变量消元得到 Schur 补方程求乘子 [ \lambda ]，再回代 [ \Delta\mathbf{x}_k ]。
- **与其他公式的内在联系**：
  - 把连续模型桥接到离散代数系统，直接作为第 3 节优化器输入。
  - 产生最终步长方向 [ \mathbf{p}_k ]，再进入线搜索更新 [ \mathbf{x}_{k+1} ]。
  - 与编号公式 [ (30) ] 在文中直接互引。

### 公式 101

- **出现位置**：`explain.md` 第 705-707 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：[ (31) ]。
- **公式完整表达**：
[
(\mathbf{Q}_k^{\mathrm{T}} \mathbf{B}_k^{-1} \mathbf{Q}_k) \lambda = -\mathbf{Q}_k^{\mathrm{T}} \mathbf{B}_k^{-1} \mathbf{r}_k - \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k). \quad (31)
]
- **公式含义**：
  - 该式定义残量向量，残量为零对应离散平衡/最优性条件成立。
  - 该式属于对偶/KKT 路线：通过 Schur 补先求乘子再回代主变量。
  - 本式位于章节 [#### 3.4.3. Schur complement for the dual approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\begin{bmatrix} \mathbf{B}_k & \mathbf{Q}_k \\ \mathbf{Q}_k^{\mathrm{T}} & \mathbf{0} \end{bmatrix} \begin{bmatrix} \Delta \mathbf{x}_k \\ \lambda \end{bmatrix} = \begin{bmatrix} -\mathbf{r}_k \\ \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k) \end{bmatrix}. \quad (30)；A popular method to solve the above linear system is to use the Schur complement approach by firstly solving \(\lambda\) from；Then, solve \(\Delta \mathbf{x}_k\) from；\mathbf{B}_k \Delta \mathbf{x}_k = -(\mathbf{r}_k + \mathbf{Q}_k \lambda). \quad (32)
- **必要推导过程**：
  1. 对离散势能对节点自由度求偏导得到残量。
  2. 残量组装后形成非线性代数方程组。
  3. 写出 KKT 增广系统。
  4. 对主变量消元得到 Schur 补方程求乘子 [ \lambda ]，再回代 [ \Delta\mathbf{x}_k ]。
- **与其他公式的内在联系**：
  - 把连续模型桥接到离散代数系统，直接作为第 3 节优化器输入。
  - 产生最终步长方向 [ \mathbf{p}_k ]，再进入线搜索更新 [ \mathbf{x}_{k+1} ]。
  - 与编号公式 [ (31) ] 在文中直接互引。

### 公式 102

- **出现位置**：`explain.md` 第 711-713 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：[ (32) ]、[ (31) ]。
- **公式完整表达**：
[
\mathbf{B}_k \Delta \mathbf{x}_k = -(\mathbf{r}_k + \mathbf{Q}_k \lambda). \quad (32)
]
- **公式含义**：
  - 该式定义残量向量，残量为零对应离散平衡/最优性条件成立。
  - 该式属于对偶/KKT 路线：通过 Schur 补先求乘子再回代主变量。
  - 本式位于章节 [#### 3.4.3. Schur complement for the dual approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：(\mathbf{Q}_k^{\mathrm{T}} \mathbf{B}_k^{-1} \mathbf{Q}_k) \lambda = -\mathbf{Q}_k^{\mathrm{T}} \mathbf{B}_k^{-1} \mathbf{r}_k - \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k). \quad (31)；Then, solve \(\Delta \mathbf{x}_k\) from；Notice that \(\mathbf{B}_k^{-1}\) is contained in Eq. (31). Recall that \(\mathbf{B}_k\) only has the compact representation form. Therefore, similar to the primal approach, either the direct matrix factorization based on the Sherman–Morrison–Woodbury formula or the conjugate gradient method needs to be applied to obtain \(\mathbf{B}_k^{-1}\). To make the situation more complex, on top of the matrix factorization or the conjugate gradient method required for \(\mathbf{B}_k^{-1}\), \((\mathbf{Q}_k^{\mathrm{T}} \mathbf{B}_k^{-1} \mathbf{Q}_k)^{-1}\) has to be formed either via another direct linear solver or iterative linear solver to solve for the Lagrange multiplier \(\lambda\). Therefore, the Schur complement approach, which generally works well to solve a Karush–Kuhn–Tucker (KKT) system, is not practical for the current problem due to the lack of the component form of the L-BFGS matrix \(\mathbf{B}_k\).；Once the free variables \(\hat{\mathbf{x}}\) at the generalized Cauchy point \(\mathbf{x}^c\) are solved from the reduced linear system shown in Eq. (27), the search direction at the current \(k\) th L-BFGS iteration can be defined as
- **必要推导过程**：
  1. 对离散势能对节点自由度求偏导得到残量。
  2. 残量组装后形成非线性代数方程组。
  3. 写出 KKT 增广系统。
  4. 对主变量消元得到 Schur 补方程求乘子 [ \lambda ]，再回代 [ \Delta\mathbf{x}_k ]。
- **与其他公式的内在联系**：
  - 把连续模型桥接到离散代数系统，直接作为第 3 节优化器输入。
  - 产生最终步长方向 [ \mathbf{p}_k ]，再进入线搜索更新 [ \mathbf{x}_{k+1} ]。
  - 与编号公式 [ (32) ]、[ (31) ] 在文中直接互引。

### 公式 103

- **出现位置**：`explain.md` 第 719-721 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：[ (27) ]。
- **公式完整表达**：
[
\mathbf{p}_k = \mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}} - \mathbf{x}_k.
]
- **公式含义**：
  - 该式控制步长 [ \alpha_k ] 的选取，在下降性与曲率条件之间平衡。
  - 该式属于原始子空间法：固定活动变量，仅在自由变量子空间解线性系统。
  - 本式位于章节 [#### 3.4.3. Schur complement for the dual approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：Notice that \(\mathbf{B}_k^{-1}\) is contained in Eq. (31). Recall that \(\mathbf{B}_k\) only has the compact representation form. Therefore, similar to the primal approach, either the direct matrix factorization based on the Sherman–Morrison–Woodbury formula or the conjugate gradient method needs to be applied to obtain \(\mathbf{B}_k^{-1}\). To make the situation more complex, on top of the matrix factorization or the conjugate gradient method required for \(\mathbf{B}_k^{-1}\), \((\mathbf{Q}_k^{\mathrm{T}} \mathbf{B}_k^{-1} \mathbf{Q}_k)^{-1}\) has to be formed either via another direct linear solver or iterative linear solver to solve for the Lagrange multiplier \(\lambda\). Therefore, the Schur complement approach, which generally works well to solve a Karush–Kuhn–Tucker (KKT) system, is not practical for the current problem due to the lack of the component form of the L-BFGS matrix \(\mathbf{B}_k\).；Once the free variables \(\hat{\mathbf{x}}\) at the generalized Cauchy point \(\mathbf{x}^c\) are solved from the reduced linear system shown in Eq. (27), the search direction at the current \(k\) th L-BFGS iteration can be defined as；The updated solution is obtained as；\mathbf{x}_{k + 1} = \mathbf{x}_k + \alpha_k \mathbf{p}_k,
- **必要推导过程**：
  1. 先给定可行方向 [ \mathbf{p}_k ]。
  2. 再用 Armijo/Wolfe 条件选 [ \alpha_k ]，兼顾下降与曲率。
  3. 构造自由子空间基 [ \mathbf{Z}_k ] 并降维。
  4. 解 [ \hat{\mathbf{B}}_k\hat{\mathbf{x}}=-\hat{\mathbf{r}}_k ] 后回代到全空间。
- **与其他公式的内在联系**：
  - 与搜索方向方程联立，保证每步更新既可行又下降。
  - 产生最终步长方向 [ \mathbf{p}_k ]，再进入线搜索更新 [ \mathbf{x}_{k+1} ]。
  - 与编号公式 [ (27) ] 在文中直接互引。

### 公式 104

- **出现位置**：`explain.md` 第 725-727 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：未显式编号。
- **公式完整表达**：
[
\mathbf{x}_{k + 1} = \mathbf{x}_k + \alpha_k \mathbf{p}_k,
]
- **公式含义**：
  - 该式控制步长 [ \alpha_k ] 的选取，在下降性与曲率条件之间平衡。
  - 本式位于章节 [#### 3.4.3. Schur complement for the dual approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\mathbf{p}_k = \mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}} - \mathbf{x}_k.；The updated solution is obtained as；where the step length parameter \(\alpha_k\) is determined by simultaneously satisfying the following two conditions:；1. The updated solution \(\mathbf{x}_{k + 1}\) should stay feasible, that is,
- **必要推导过程**：
  1. 先给定可行方向 [ \mathbf{p}_k ]。
  2. 再用 Armijo/Wolfe 条件选 [ \alpha_k ]，兼顾下降与曲率。
- **与其他公式的内在联系**：
  - 与搜索方向方程联立，保证每步更新既可行又下降。

### 公式 105

- **出现位置**：`explain.md` 第 732-734 行；章节：[ #### 3.4.3. Schur complement for the dual approach ]；文中编号：[ (18) ]。
- **公式完整表达**：
[
\mathbf{lb} \leq \mathbf{x}_k + \alpha_k \mathbf{p}_k \leq \mathbf{ub}.
]
- **公式含义**：
  - 该式控制步长 [ \alpha_k ] 的选取，在下降性与曲率条件之间平衡。
  - 本式位于章节 [#### 3.4.3. Schur complement for the dual approach]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：where the step length parameter \(\alpha_k\) is determined by simultaneously satisfying the following two conditions:；1. The updated solution \(\mathbf{x}_{k + 1}\) should stay feasible, that is,；2. The positive step length \(\alpha_k\) should satisfy the strong Wolfe conditions such that the curvature condition in Eq. (18) is satisfied to ensure the positive-definiteness of the updated L-BFGS matrix \(\mathbf{B}_{k + 1}\).；Then, a new vector-pair \(\{\mathbf{s}_k, \mathbf{y}_k\}\) is computed. The oldest vector-pair \(\{\mathbf{s}_{k-m}, \mathbf{y}_{k-m}\}\) is removed from the vector-pair list \(\{\mathbf{s}_i, \mathbf{y}_i\}_{i=k-m+1}^{k-1}\) and the newly computed vector-pair \(\{\mathbf{s}_k, \mathbf{y}_k\}\) is inserted to form the new vector-pair list \(\{\mathbf{s}_i, \mathbf{y}_i\}_{i=k-m+1}^{k}\) for the next L-BFGS iteration. Algorithm 1 summaries the major steps of the gradient projection L-BFGS-B scheme.
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
  3. 先给定可行方向 [ \mathbf{p}_k ]。
  4. 再用 Armijo/Wolfe 条件选 [ \alpha_k ]，兼顾下降与曲率。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 与搜索方向方程联立，保证每步更新既可行又下降。
  - 与编号公式 [ (18) ] 在文中直接互引。

### 公式 106

- **出现位置**：`explain.md` 第 780-782 行；章节：[ ### 4.1. Cyclic tension-compression test ]；文中编号：[ (33) ]。
- **公式完整表达**：
[
E_{\mathrm{crack}} = g_{\mathrm{c}}\Gamma_{l}(d) = g_{\mathrm{c}}\int_{\Omega}\gamma (d,\nabla d)\mathrm{d}\Omega = g_{\mathrm{c}}\int_{\Omega}\frac{1}{2l}\left(d^{2} + l^{2}\nabla d\cdot \nabla d\right)\mathrm{d}\Omega \quad (33)
]
- **公式含义**：
  - 该式给出相场裂纹表面密度正则项，利用长度尺度 [ l ] 把尖锐裂纹转化为可微弥散带。
  - 该式定义算例中的边界控制量或后处理能量，是数值结果物理解释的直接接口。
  - 本式位于章节 [### 4.1. Cyclic tension-compression test]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：The 2D example shown in Fig. 3(a) is widely adopted in the literature of the phase-field fracture simulation to test the robustness of the numerical method. Under the displacement-controlled load shown in Fig. 3(b) that is vertically applied at the top edge of the unit square, the preexisting crack is expected to propagate from the domain center to the right edge in a unstable fashion. Since the crack propagation is known a priori, a preferred mesh shown in Fig. 3(a) is sufficient and no adaptive mesh refinement is necessary. The material parameters include the Lamé parameters \(\lambda = 121.15 \mathrm{kN/mm^2}\) and \(\mu = 80.77 \mathrm{kN/mm^2}\), the critical energy release rate \(g_c = 2.7 \times 10^{-3} \mathrm{kN/mm}\), and the small positive parameter \(k = 0.0\) in the degradation function \(g(d)\). The phase-field length-scale parameter is chosen as \(l = 0.0075 \mathrm{mm}\), and the mesh size in the preferred region is \(h = 0.00265 \mathrm{mm}\) such that the ratio is \(h/l \approx 1/3\).；Fig. 4 shows the phase-field distribution in the unit square under the vertical cyclic load at various pseudo time steps. When the vertical displacement increases from \(u_{y} = 5.9 \times 10^{-3} \mathrm{mm}\) to \(u_{y} = 6.0 \times 10^{-3} \mathrm{mm}\), the crack abruptly propagates from the preexisting crack tip at the domain center to the right edge of the domain. Moreover, after the crack is fully developed, the phase-field distribution remains unchanged in the subsequent load step, and the phase-field nodal value never decreases during the unloading and compression stages due to the enforced irreversibility condition. Also, the phase-field value is strictly in the range of \([0, 1]\). Fig. 5(a) shows the relationship between the vertical reaction force at the bottom edge and the pseudo time. Fig. 5(b) shows the relationship between the vertical reaction force at the bottom edge and the vertical displacement at the top edge. The sample loses the tensile load bearing capacity after the crack is fully developed at \(t = 6.0 \times 10^{-3}\). Subsequently, the crack never self-heals during the unloading stage due to the enforced phase-field (damage) irreversibility. The sample starts to be compressed after \(t = 14.0 \times 10^{-3}\). Due to the tension-compression asymmetry, the sample still retains the fully elastic behavior under compression. Fig. 5(c) compares the crack dissipation energy calculated as；with and without enforcing the phase-field irreversibility condition. When the irreversibility condition is properly enforced by the proposed L-BFGS-B scheme, the crack dissipation energy remains constant after the crack is fully developed. When the irreversibility condition is not enforced, the crack dissipation energy decreases to zero during the unloading stage, which obviously violates the thermodynamic consistency.；Fig. 6 shows the total energy of the system, the \(l_{2}\)-norm of the projected gradient of the displacement field \(\| \mathrm{Proj}_{c}^{\pmb{u}}\|_{2}\) and the phase-field \(\| \mathrm{Proj}_{c}^{d}\|_{2}\) according to Eq. (19), and the number of the active constraints during the L-BFGS-B iterations at the critical load step when the crack unstably propagates. Due to the line search method based on the strong Wolfe conditions, the total energy of the system consistently decreases during the iterations and eventually reaches to a plateau. The convergence criteria are satisfied when the \(l_{2}\)-norms of the projected gradient and the solution increments are reduced to the prescribed tolerance \((10^{-6})\) and the active constraints remain unchanged during two consecutive iterations.
- **必要推导过程**：
  1. 由相场正则化设定裂纹密度 [ \gamma(d,\nabla d)=\frac{1}{2l}(d^2+l^2|\nabla d|^2) ]。
  2. 对全域积分后得到近似裂纹表面积函数 [ \Gamma_l(d) ]。
  3. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  4. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
  5. 将控制参数（如扭转角、位移幅值）代入边界位移表达。
  6. 在后处理中代入能量积分定义得到可比物理量。
- **与其他公式的内在联系**：
  - 向后连接到变分残量 [ r_u,r_d ] 与离散切线矩阵，是全篇物理起点。
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 把算法解映射到可观测物理指标，用于与实验/文献曲线比较。
  - 与编号公式 [ (33) ] 在文中直接互引。

### 公式 107

- **出现位置**：`explain.md` 第 878-880 行；章节：[ ### 4.4. Three-dimensional torsion test ]；文中编号：未显式编号。
- **公式完整表达**：
[
u_{y} = z \tan t, \quad u_{z} = -y \tan t,
]
- **公式含义**：
  - 该式定义算例中的边界控制量或后处理能量，是数值结果物理解释的直接接口。
  - 本式位于章节 [### 4.4. Three-dimensional torsion test]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：### 4.4. Three-dimensional torsion test；In the last example, a 3D beam shown in Fig. 17 undergoes a torsional load. The beam has a square cross section of \(50 \times 50 \mathrm{mm}\) and a length of \(200 \mathrm{mm}\). In the middle of the beam, there is a preexisting \(45^{\circ}\) notch. The width of the notch opening is \(5 \mathrm{mm}\), and the depth of the notch is \(25 \mathrm{mm}\). The right end of the beam is fixed in all three directions. A displacement controlled rotation is applied at the left end of the beam so that it undergoes torsion. At an arbitrary node \((0,y,z)\) located at the left surface of the beam, the Dirichlet boundary conditions are applied as；where \(t\) is the pseudo time representing the load step. The initial mesh is pre-refined near the bottom region of the notch where the crack is expected to initialize. The material parameters include the Lamé parameters \(\lambda = 9.72 \mathrm{kN/mm^2}\) and \(\mu = 14.58 \mathrm{kN/mm^2}\), the critical energy release rate \(g_{c} = 1.1 \times 10^{-4} \mathrm{kN/mm}\), and the small positive parameter \(k = 0.0\) in the degradation function \(g(d)\). The phase-field length-scale parameter is chosen as \(l = 1.0 \mathrm{mm}\).；Fig. 18 shows different views of the crack surface (the phase-field value \(d \geq 0.8\)) and the corresponding adaptively refined mesh when the rotation angle (pseudo time) \(t = 8.0 \times 10^{-3}\). Obviously, the crack surface is non-planar under the torsional load. With the adaptive mesh refinement technique, the total number of DoFs is still more than one million. From this example, we can see that adaptive mesh refinement is absolutely necessary for any practical 3D phase-field fracture simulations to control the computational cost. Furthermore, Fig. 19 shows the constraint status of the phase-field degree of freedom (DoF) at each finite element node. The lower bounds of the box constraints prevent the phase-field from decreasing, while the upper bounds of the box constraints prevent the phase-field from going over 1.0.
- **必要推导过程**：
  1. 将控制参数（如扭转角、位移幅值）代入边界位移表达。
  2. 在后处理中代入能量积分定义得到可比物理量。
- **与其他公式的内在联系**：
  - 把算法解映射到可观测物理指标，用于与实验/文献曲线比较。

### 公式 108

- **出现位置**：`explain.md` 第 937-939 行；章节：[ ### 5.2. Comparison of convergence behaviors ]；文中编号：未显式编号。
- **公式完整表达**：
[
\| \pmb {r}_u\| _2 < \mathrm{tol}, \quad \| \pmb {r}_d\| _2 < \mathrm{tol}, \quad \| \Delta \pmb {u}\| _2 < \mathrm{tol}, \quad \| \Delta d\| _2 < \mathrm{tol},
]
- **公式含义**：
  - 该式在当前章节承担定义/变换/约束中的一个局部环节，用于连接前后公式。
  - 本式位于章节 [### 5.2. Comparison of convergence behaviors]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：All the three numerical tests go through a complete loading-unloading cycle, and the required number of iterations is reported for each load step. For the tensile test, during the loading stage, the vertical displacement \(u_{y}\) first increases to \(5.0 \times 10^{-3} \mathrm{mm}\) with the step size \(1.0 \times 10^{-3} \mathrm{mm}\). Then, \(u_{y}\) further increases to \(7.0 \times 10^{-3} \mathrm{mm}\) with the step size \(1.0 \times 10^{-4} \mathrm{mm}\). Afterwards, \(u_{y}\) decreases to zero with the step size \(-1.0 \times 10^{-3} \mathrm{mm}\) during the unloading stage to finish a complete loading-unloading cycle containing totally 32 load steps. For the shear test, during the loading stage, the horizontal displacement \(u_{x}\) increases to \(15.0 \times 10^{-3} \mathrm{mm}\) with the step size \(1.0 \times 10^{-3} \mathrm{mm}\). During the unloading stage, the horizontal displacement decreases to zero with the step size \(-1.0 \times 10^{-3} \mathrm{mm}\) to finish a complete loading-unloading cycle with totally 30 load steps. For the bending test, the vertical displacement \(u_{y}\) first increases to \(1.0 \mathrm{mm}\) with the step size \(2.0 \times 10^{-2} \mathrm{mm}\) during the loading stage. Then, \(u_{y}\) decreases to zero with the step size \(-2.0 \times 10^{-2} \mathrm{mm}\) during the unloading stage with totally 100 load steps for the whole loading-unloading cycle.；For the staggered scheme using the alternate minimization, inside one iteration, the nonlinear displacement sub-problem is firstly solved via the Newton-Raphson method. Then, the obtained displacement field is fixed and the linear phase-field sub-problem is solved subsequently. For the two monolithic schemes, each iteration only involves one solve of the monolithic linear system. For the proposed L-BFGS-B solver, the detailed convergence criteria are listed under Section 3.1. For the staggered scheme and the limited-memory BFGS scheme using the history variable to enforce irreversibility, the following residual-based convergence criteria are adopted:；where \(\pmb{r}_{u}\) and \(\pmb{r}_{d}\) represent the residuals of the displacement and the phase-field sub-problems, and \(\Delta \pmb {u}\) and \(\Delta d\) represent the increments of the displacement and the phase-field solutions during one iteration. To ensure a fair comparison regarding the solver convergence behaviors, all the three solvers adopt the same tolerance value tol \(= 10^{-6}\). Notice that for the proposed L-BFGS-B approach, besides checking the \(l_{2}\)-norms of the residual in the form of the projected gradient and the solution increment, the active set has to remain unchanged between two consecutive iterations.；Similar to the proposed L-BFGS-B solver, the source code and the input files for the staggered approach (A.M.-hist.) and the limited-memory BFGS approach (L-BFGS-hist.) are also hosted on GitHub in support of reproducible research. Fig. 22 reports the number of iterations required for convergence in each load step during a full loading-unloading cycle for the three (tensile, shear, bending) test problems shown in Fig. 21. Both the proposed L-BFGS-B scheme and the BFGS scheme using the history variable for irreversibility (L-BFGS-hist.) require fewer iterations than the staggered scheme (A.M.-hist.) for the shear problem in Fig. 21(b) and the bending problem in Fig. 21(c). For the tensile problem, even though the staggered approach (A.M.-hist.) requires fewer alternate minimization (outer-level) iterations, considering the fact that each alternate minimization iteration involves a nonlinear solve of the displacement sub-problem using 2 to 5 Newton iterations, the total number of required linear solves in the staggered approach is still larger than the counterparts of the two monolithic schemes (L-BFGS-B and L-BFGS-hist.). Between the two monolithic schemes, the proposed L-BFGS-B solver requires slightly more iterations inside each load step than the BFGS solver using the history variable approach (L-BFGS-hist.). This is not surprising because besides satisfying the tolerances of the residuals and the solution increments, the convergence criteria of the former also require the convergence of all the active constraints. Moreover, the iterates have to always stay feasible due to the gradient projection in the proposed L-BFGS-B method. On the upside, the proposed L-BFGS-B solver has the advantage of ensuring that the phase-field is always between 0 and 1, whereas the solvers based on the history variable approach obtain phase-field values larger than 1 or smaller than 0, which is a known drawback for this technique.
- **必要推导过程**：
  1. 根据本式中的等号关系执行代入、求导或约束投影，可从上一式直接得到当前式。
- **与其他公式的内在联系**：
  - 该式在局部推导链中承接前式并提供下一式所需变量定义。

### 公式 109

- **出现位置**：`explain.md` 第 973-975 行；章节：[ ## Appendix ]；文中编号：未显式编号。
- **公式完整表达**：
[
\hat{\mathbf{A}} = \mathbf{A} + \mathbf{U}\mathbf{V}^{\mathrm{T}},
]
- **公式含义**：
  - 该式是低秩逆更新恒等式，解释为何可把大规模逆运算化为小规模逆运算。
  - 本式位于章节 [## Appendix]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：## Appendix；Let \(\mathbf{U}\in \mathbb{R}^{n\times p}\) and \(\mathbf{V}\in \mathbb{R}^{n\times p}\) represent two matrices for some \(p\) between 1 and \(n\). If；according to the Sherman–Morrison–Woodbury formula [45], \(\hat{\mathbf{A}}\) is non-singular if and only if \((\mathbf{I} + \mathbf{V}^{\mathrm{T}}\mathbf{A}^{-1}\mathbf{U})\) is non-singular. In this case,；\hat{\mathbf{A}}^{-1} = \mathbf{A}^{-1} - \mathbf{A}^{-1}\mathbf{U}\left(\mathbf{I} + \mathbf{V}^{\mathrm{T}}\mathbf{A}^{-1}\mathbf{U}\right)^{-1}\mathbf{V}^{\mathrm{T}}\mathbf{A}^{-1}.
- **必要推导过程**：
  1. 把矩阵分解为 [ \hat{\mathbf{A}}=\mathbf{A}+\mathbf{U}\mathbf{V}^T ]。
  2. 应用 Sherman–Morrison–Woodbury 公式把逆转化为低维修正项。
- **与其他公式的内在联系**：
  - 与 L-BFGS 紧凑表示形成“低秩更新 + 快速逆”配套关系。

### 公式 110

- **出现位置**：`explain.md` 第 979-981 行；章节：[ ## Appendix ]；文中编号：未显式编号。
- **公式完整表达**：
[
\hat{\mathbf{A}}^{-1} = \mathbf{A}^{-1} - \mathbf{A}^{-1}\mathbf{U}\left(\mathbf{I} + \mathbf{V}^{\mathrm{T}}\mathbf{A}^{-1}\mathbf{U}\right)^{-1}\mathbf{V}^{\mathrm{T}}\mathbf{A}^{-1}.
]
- **公式含义**：
  - 该式是低秩逆更新恒等式，解释为何可把大规模逆运算化为小规模逆运算。
  - 本式位于章节 [## Appendix]，其符号结构决定它在“建模—离散—优化”链条中的具体角色。
  - 邻近语境提示：\hat{\mathbf{A}} = \mathbf{A} + \mathbf{U}\mathbf{V}^{\mathrm{T}},；according to the Sherman–Morrison–Woodbury formula [45], \(\hat{\mathbf{A}}\) is non-singular if and only if \((\mathbf{I} + \mathbf{V}^{\mathrm{T}}\mathbf{A}^{-1}\mathbf{U})\) is non-singular. In this case,；## Data availability；All the source code and input files are hosted at https://github.com/taojinlln/Phasefield_gradient_projection_monolithic_solver.
- **必要推导过程**：
  1. 由不等式约束定义投影算子 [ \mathrm{Proj}_C ]。
  2. 将任意试探点截断到区间 [ \mathrm{lb}_i,\mathrm{ub}_i ]。
  3. 把矩阵分解为 [ \hat{\mathbf{A}}=\mathbf{A}+\mathbf{U}\mathbf{V}^T ]。
  4. 应用 Sherman–Morrison–Woodbury 公式把逆转化为低维修正项。
- **与其他公式的内在联系**：
  - 与广义 Cauchy 点和子空间最小化共同构成约束处理闭环。
  - 与 L-BFGS 紧凑表示形成“低秩更新 + 快速逆”配套关系。

## 关键编号公式的跨章节对应关系（压缩版索引）

- 建模链：[ (1) ] → [ (2) ] → [ (3) ] → [ (4) ]。
- 变分离散链：[ (7) ] → [ (8) ] → [ (9) ] → [ (10) ] → [ (11) ]。
- 约束优化链：[ (12) ] → [ (13) ] → [ (15) ] → [ (16) ] → [ (17) ] → [ (18) ] → [ (19) ]。
- L-BFGS/GCP 链：[ (20) ]、[ (21) ] → [ (22) ]、[ (23) ]、[ (24) ]。
- 子空间/KKT 链：[ (25) ]、[ (26) ] → [ (27) ]、[ (28) ]、[ (29) ] → [ (30) ]、[ (31) ]、[ (32) ]。
- 工程量链：[ (33) ] 与边界位移公式共同解释数值结果。

---

如需继续细化，我可以在下一版把每条“必要推导过程”再展开成**逐行代数推导**（尤其是 [ (7) ]、[ (10) ]、[ (20) ]、[ (30) ]）。
