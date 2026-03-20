# Gradient projection method for enforcing crack irreversibility as box constraints in a robust monolithic phase-field scheme

Tao Jin

Department of Mechanical Engineering, University of Ottawa, Ottawa Ontario, K1N 6N5, Canada

## ARTICLE INFO

Dataset link: https://github.com/taojinlln/Phasefield_gradient_projection_monolithic_solver

Keywords: Phase-field fracture Monolithic scheme Irreversibility Box constraints Gradient projection method L-BFGS-B method

## ABSTRACT

A phase-field monolithic scheme based on the gradient projection method is developed to model crack propagation in brittle materials under cyclic loading. As a type of active set method, the gradient projection method is particularly attractive to enforce the irreversibility condition imposed on the phase-field variables as bound constraints, or box constraints. This method has the advantages of allowing the rapid change of active constraints during iterations and computing the projected gradient with a negligible cost. The gradient projection method is further combined with the limited-memory BFGS (L-BFGS) method to overcome the convergence difficulties arising from the non-convex energy functional. A compact representation of the BFGS matrix is adopted as the limited-memory feature to avoid the storage of fully dense matrices, making this method practical for large-scale finite element simulations. By locating the generalized Cauchy point on the piecewise linear path formed by the projected gradient, the active set of box constraints can be determined. The variables in the active set, which are at the boundary of the box constraints, are kept fixed to form a subspace minimization problem. A primal approach and a dual approach are presented to solve this subspace minimization problem for the remaining free variables at the generalized Cauchy point. Several two-dimensional (2D) and three-dimensional (3D) examples are provided to demonstrate the capabilities of the proposed monolithic scheme, particularly in enforcing the phase-field irreversibility during crack propagation under cyclic loading. In these numerical examples, the proposed monolithic scheme is combined with an adaptive mesh refinement technique to alleviate the heavy computational cost incurred by the fine mesh resolution required around the crack region. The proposed method is further compared with two other phase-field solving techniques regarding the convergence behavior. To ensure a fair comparison, the same problem settings and implementation techniques are adopted. The proposed monolithic scheme provides a unified framework to overcome the numerical difficulties associated with the non-convex energy functional, effectively enforce the phase-field irreversibility to ensure the thermodynamic consistency, and alleviate the heavy computational cost through adaptive mesh refinement in 2D and 3D phase-field crack simulations.

## 1. Introduction

During the past decade, the phase-field method becomes a popular technique for modeling fracture propagation, particularly in brittle materials, due to its capability of naturally handling complex crack geometry and propagation paths such as branching and merging [1-10]. Comparing with the traditional level-set based tracking strategy [11-13] that attempts to represent the discontinuity using an auxiliary scalar field, the phase-field method relies on a variational approach [14-16] and regularizes the sharp crack geometry in a diffusive manner. According to Francfort and Marigo [14], the energy functional of a fractured quasi-static elastic solid system is expressed as

\[
\Pi (\pmb {u},d) = \int_{\Omega}\psi (\pmb {e}(\pmb {u}),d)\mathrm{d}\Omega + g_{c}\Gamma_{l}(d) - \int_{\Omega}\pmb {b}\cdot \pmb {u}\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot \pmb {u}\mathrm{d}\Gamma , \quad (1)
\]

where \(\pmb {u}(\pmb {x})\) represents the unknown vector displacement field, \(d(\pmb {x})\) represents the unknown scalar phase-field, \(\pmb{b}\) is the body force, \(\pmb{t}\) is the traction load, \(\pmb {e} = \nabla^{(s)}\pmb {u}\) is the small deformation linear strain tensor, \(\psi\) represents the strain energy density function, \(g_{c}\) is the critical energy release rate, and \(\Gamma_{l}\) is an approximation of the crack surface area. Particularly, the approximated crack surface \(\Gamma_{l}\) is defined as [1]

\[
\Gamma_{l}(d) = \int_{\Omega}\gamma (d,\nabla d)\mathrm{d}\Omega = \int_{\Omega}\frac{1}{2l}\left(d^{2} + l^{2}\nabla d\cdot \nabla d\right)\mathrm{d}\Omega , \quad (2)
\]

where \(\gamma (d,\nabla d)\) is considered as the crack surface density function, and \(l\) is the phase-field length-scale parameter.

For a quasi-static problem, let \(t\) represent the pseudo time instead of the real time. The pseudo time \(t\) enters the problem as the load step through time-dependent boundary conditions, for instance, the displacement \(\pmb {u} = \hat{\pmb{u}} (t)\) prescribed as the Dirichlet boundary condition or the pressure \(p = \hat{p} (t)\) prescribed as the Neumann boundary condition. During a typical time (load) step \(\left[t_{n},t_{n + 1}\right]\), let \(\left(\pmb{u}_{n},d_{n}\right)\in \mathbf{V}\times \mathbf{W}\) represent the solution from the previous time step, where \(\mathbf{V} = \mathbf{H}_{0}^{1}(\Omega)\) and \(\mathbf{W} = \mathbf{H}^{1}(\Omega)\). The primary unknown fields \(\left(\pmb{u}_{n + 1},d_{n + 1}\right)\in \mathbf{V}\times \mathbf{W}\) can be obtained from the following minimization [4,17]:

\[
\left(\pmb{u}_{n + 1},d_{n + 1}\right) = \arg \min \Pi (\pmb {u},d) \quad (3)
\]

subject to the inequality constraints

\[
0\leq d_{n}\leq d_{n + 1}\leq 1. \quad (4)
\]

The above inequality constraints represent the following considerations. First, the phase-field cannot decrease in order to ensure the thermodynamic consistency. Second, the phase-field should be between \(d = 0\) (the undamaged state) and \(d = 1\) (the fully damaged state). The initial damage state \(d_{0}\) of the computational domain \(\Omega\) is known. For convenience, initially \((t = 0)\) the computational domain \(\Omega\) is assumed to be in the undamaged state, that is, \(d_{0}(\pmb{x}) = 0\).

Applying the phase-field approach to model fracture propagation encounters the following three challenges. First, it is well known that the energy functional shown in Eq. (1) is non-convex [15,16,18]. As a result, the Newton-based monolithic scheme usually faces convergence difficulties. Therefore, carefully designed numerical techniques are required to solve the nonlinear system. Second, cracks do not self-heal after formation. The phase-field (damage) irreversibility needs to be enforced for the thermodynamic consistency, which renders the phase-field formulation into a constrained minimization with inequality constraints, as shown in Eqs. (3) and (4). Third, the phase-field approach to represent fracture is intrinsically expensive, since the phase-field length-scale needs to be resolved around the crack region. Discretizing the entire computational domain with the same level of high mesh resolution is impractical, particularly for 3D problems. The adoption of adaptive mesh refinement technique is necessary. In the literature, a large body of work has been devoted to overcome the aforementioned three challenges.

To overcome the numerical challenges associated with the non-convex energy functional, several modified Newton methods were proposed [6,19,20]. The common idea of these modified Newton methods is to modify the Jacobian (Hessian) matrix to restore its positive definiteness. However, the modification process involves several heuristic parameters that need to be selected. Other efforts for the development of robust monolithic schemes include a fracture energy based arc-length method [21], a multilevel trust region method [22], and a novel line-search assisted procedure [23]. Particularly, both Wu et al. [8] as well as Kristensen and Martinez-Pameda [24] advocated the adoption of the Broyden-Fletcher-Goldfarb-Shanno (BFGS) method and successfully demonstrated its robustness and efficiency. Moreover, Khalil et al. [25] adopted the BFGS method to investigate the fatigue crack growth in metals considering the material elastic-plastic responses. Recently, Jin et al. [26] proposed a monolithic phase-field approach based on the so-called limited-memory BFGS (L-BFGS) method [27-29]. Integrated with a predictor-corrector adaptive mesh refinement technique [4], the L-BFGS monolithic approach is shown to be extremely robust and computationally efficient both for 2D and 3D problems. On a different front, staggered approaches [2,30-32] based on the alternate minimization process are developed to circumvent the numerical difficulties arising from the non-convexity of the energy functional. However, these approaches require a large amount of iterations between the two convex sub-problems (the displacement sub-problem and the phase-field sub-problem).

The phase-field (damage) irreversibility condition needs to be enforced to ensure the thermodynamic consistency, which subjects the minimization of the energy functional to inequality constraints. Several technical approaches to enforce the phase-field irreversibility are highlighted below. The history variable approach, popularized by Miehe et al. [2], uses the local maximum positive strain energy as the driving force for crack propagation. In this approach, only the maximum positive strain energy (scalar) needs to be stored locally at each quadrature point as a history variable, and the original inequality-constrained minimization is converted to an unconstrained minimization. However, because the positive strain energy is replaced by the history variable (maximum positive strain energy in history), the derived weak form for the finite element formulation is inconsistent with the underlying energy functional, and consequently, destroys the variational structure. Instead of using the history variable, the penalty method provides an alternative approach to enforce the phase-field irreversibility, for instance, see the early work by Miehe et al. [1]. The difficulty associated with the penalty method is to choose the proper value of the penalty parameter such that the irreversibility constraint is effectively enforced without significantly increasing the condition number of the system. For this purpose, Gerasimov and De Lorenzis [33] analytically derived the lower bound of the penalty parameter that is a function of the fracture toughness and the phase-field length-scale. Wheeler et al. [17] and Wick [6] further adopted the augmented Lagrangian method to enforce the irreversibility condition. Departing from the penalty approach, Heister et al. [4] proposed a primal-dual active set method, which is a type of semi-smooth Newton method [34], to enforce the phase-field irreversibility. However, this work relies on an heuristic phase-field extrapolation approach to circumvent the numerical difficulties associated with the non-convexity of the energy functional. Graiser et al. [35] proposed a so-called truncated nonsmooth Newton multigrid method to directly enforce the pointwise phase-field irreversibility constraints. They further proved the convergence of the solver regardless of any load and initial iterate. To directly solve the phase-field crack formulation as an inequality-constrained minimization problem, Wambacq et al. [19] proposed an interior-point method to treat the phase-field irreversibility as a variational inequality.

The phase-field crack formulation is computationally expensive since highly refined meshes around the crack path are needed to represent the sharp crack in a diffusive way. Globally refining the entire domain to achieve the required mesh resolution around the crack region is expensive for 2D problems and simply impractical for 3D problems. Several research efforts are devoted to develop the adaptive mesh refinement technique [4,36-38]. Particularly, Heister et al. [4] proposed a predictor-corrector local mesh adaptivity scheme. This scheme does not rely on the knowledge of the crack path a priori. Rather, it solely uses the phase-field solution as the local refinement criterion and performs the adaptive mesh refinement through multiple prediction-correction cycles. This scheme is appealing due to its computational efficiency, since the phase-field is part of the primary unknowns solved from the nonlinear system and no extra quantities are required.

The gradient projection method [39-42] is a special type of active set methods that can be used to solve the inequality-constrained optimization problem. This method has the advantage of allowing the set of active constraints to change rapidly during each iteration. Fig. 1 illustrates the basic idea of the gradient projection method. Let \(f(x)\) represent the objective function to be minimized. Let \(C\) represent the feasible region that is convex and formed by a series of inequality constraints. The gradient projection method typically involves an operation that projects the descent direction \(\nabla f(x)\) onto the convex feasible region \(C\), that is,

\[
\mathrm{Proj}_C(\mathbf{x}_k = \mathbf{a}_k\nabla f(\mathbf{x}_k)).
\]

For a general convex feasible region, for example, as shown in Fig. 1(a), the projection operation \(\mathrm{Proj}_C(\cdot)\) incurs significant computational cost. Therefore, this method is not competitive comparing with other inequality-constrained minimization techniques such as the interior-point method. However, for inequality constraints that are imposed as bounds on the unknown vector \(\mathbf{x}\) as shown in Fig. 1(b), the gradient projection operation \(\mathrm{Proj}_C(\cdot)\) becomes trivial, making the gradient projection method extremely appealing [41,43,44]. The bound constraints, also known as the box constraints, can be expressed in the following component form,

\[
\mathrm{lb}_i\leq x_i\leq \mathrm{ub}_i,
\]

where \(x_{i}\) represents the \(i\) th component of the unknown vector \(\mathbf{x}\) \(\mathrm{lb}_i\) and \(\mathrm{ub}_i\) represent its corresponding lower and upper bounds, respectively. Comparing the above box constraints with the inequality constraints shown in Eq. (4) imposed on the minimization of the phase-field energy functional, we can easily observe that after the spatial and temporal discretizations, the lower and upper bounds imposed on the discretized phase-field (nodal value) become

\[
\mathrm{lb}_i = d_i^{(n)}\leq d_i^{(n + 1)}\leq 1 = \mathrm{ub}_i,
\]

where \(d_i^{(n)}\) is known at the beginning of the current time (load) step, and \(d_i^{(n + 1)}\) is the unknown phase-field value that needs to be solved.

In this paper, a phase-field monolithic scheme based on the gradient projection method is proposed to enforce the phase-field irreversibility condition as box constraints. Comparing with the history variable approach, the gradient projection method preserves the variational structure of the energy functional and ensures that the obtained phase-field value is between 0 and 1. Unlike the penalty method that requires the proper adjustment of the penalty parameter, the proposed method does not involve extra algorithmic parameters. Furthermore, this method is combined with the limited-memory BFGS (L-BFGS) method to overcome the convergence difficulties arising from the non-convex energy functional. In the field of numerical optimization, the combination of the gradient projection method and the limited-memory BFGS method forms the so-called L-BFGS-B method [41,45], in which “L” stands for limited-memory and “B” stands for box constraints. Comparing with the existing phase-field monolithic schemes in the literature, such as the primal–dual active set method [4], the augmented Lagrangian method [6], and the interior-point method [19], the proposed monolithic scheme based on the L-BFGS-B method can rigorously and efficiently enforce the phase-field irreversibility. Moreover, the BFGS quasi-Newton approach can naturally handle non-convex optimizations without the necessity of introducing any heuristic parameters. The proposed monolithic scheme is further combined with the predictor–corrector adaptive mesh refinement technique [4]. Therefore, this paper presents an integrated monolithic scheme that is able to overcome the convergence difficulties associated with the non-convex energy functional, rigorously and efficiently enforce the phase-field irreversibility, and alleviate the computational cost through the adaptive mesh refinement.

The remaining paper is organized as follows. Section 2 briefly summarizes the phase-field formulation adopted in this work and provides the expressions of the finite element discretized system. Section 3 presents in detail the proposed phase-field monolithic scheme based on the gradient projection, which includes the compact representation of the L-BFGS matrix, the procedure to calculate the generalized Cauchy point and detect the active constraints, and the solving techniques for the subsequently formed subspace minimization. Section 4 provides several numerical examples of 2D and 3D crack propagation problems. Particularly, cyclic loading conditions are applied to demonstrate the effectiveness of the developed gradient projection method on enforcing the phase-field irreversibility condition. In Section 5, the convergence behavior of the proposed L-BFGS-B method is compared with a staggered approach based on the alternate minimization and a BFGS scheme using the history variable for irreversibility. To ensure a fair comparison, the same problem settings and implementation techniques are adopted. Several concluding remarks are provided in Section 6. In support of the efforts of open-source software and reproducible research, all the source codes and input files developed in this paper are hosted on GitHub.¹

## 2. Phase-field formulation and finite element discretization

The phase-field formulation used in this work follows the typical approach widely adopted in the literature. For the sake of completeness, this section still briefly summarizes the formulation. Several important quantities, including the total energy as the objective function to be minimized, the gradient of the total energy as the residual, and the Hessian of the objective function are provided after the finite element discretization. These quantities are needed for the monolithic phase-field scheme based on the gradient projection method discussed in Section 3.

### 2.1. Phase-field formulation

In order to consider the fracture tension–compression asymmetry, the strain energy density function used in the total energy functional in Eq. (1) is additively decomposed into two parts [1],

\[
\psi(\pmb{\epsilon}, d) = [g(d) + k] \psi^{+}(\pmb{\epsilon}) + \psi^{-}(\pmb{\epsilon}), \quad (5)
\]

where \(\psi^{+}\) is the positive strain energy, \(\psi^{-}\) is the negative strain energy, \(k\) is a small non-negative number, and \(g(d)\) represents the phase-field degradation function and adopts the following form,

\[
g(d) = (1 - d)^{2}. \quad (6)
\]

In this work, the small non-negative number \(k\) is set as zero (\(k = 0\)) in all the numerical examples.

The following operators are introduced to describe the constitutive relationship based on the additive decomposition of the strain energy density function,

\[
\langle x\rangle_{+} = \frac{1}{2}(x + |x|), \quad \langle x\rangle_{-} = \frac{1}{2}(x - |x|), \quad H(x) = \begin{cases} 1 & x \geq 0, \\ 0 & x < 0. \end{cases}
\]

The spectrum decomposition of the strain tensor \(\pmb{\epsilon}\) is expressed as

\[
\pmb{\epsilon} = \sum_{\alpha} \epsilon_{\alpha} \mathbf{M}_{\alpha}, \quad \mathbf{M}_{\alpha} = \pmb{n}_{\alpha} \otimes \pmb{n}_{\alpha},
\]

where \(\epsilon_{\alpha}\) and \(\pmb{n}_{\alpha}\) represent a pair of eigenvalue and eigenvector. The positive and negative parts of the strain tensor are defined as,

\[
\pmb{\epsilon}^{+} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{+} \mathbf{M}_{\alpha}, \quad \pmb{\epsilon}^{-} = \sum_{\alpha} \langle \epsilon_{\alpha}\rangle_{-} \mathbf{M}_{\alpha}.
\]

Using the above definitions, the positive and negative parts of the strain energy are expressed as

\[
\psi^{+}(\pmb{\epsilon}) = \frac{1}{2}\lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{+}^{2} + \mu \pmb{\epsilon}^{+} : \pmb{\epsilon}^{+}, \quad \psi^{-}(\pmb{\epsilon}) = \frac{1}{2} \lambda \langle \mathrm{tr}\pmb{\epsilon}\rangle_{-}^{2} + \mu \pmb{\epsilon}^{-} : \pmb{\epsilon}^{-},
\]

where \(\lambda\) and \(\mu\) are the Lame parameters, and \(\mathrm{tr}\pmb{\epsilon}\) is the trace of the strain tensor. The stress tensor \(\pmb{\sigma}\) is derived as

\[
\pmb{\sigma} = \frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial \pmb{\epsilon}} = [g(d) + k]\frac{\partial\psi^{+}(\pmb{\epsilon})}{\partial \pmb{\epsilon}} + \frac{\partial\psi^{-}(\pmb{\epsilon})}{\partial \pmb{\epsilon}} = [g(d) + k]\pmb{\sigma}^{+} + \pmb{\sigma}^{-},
\]

where

\[
\pmb{\sigma}^{+} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{+}\mathbf{I} + 2\mu \pmb{\epsilon}^{+}, \quad \pmb{\sigma}^{-} = \lambda \langle \mathrm{tr}\pmb{\epsilon} \rangle_{-}\mathbf{I} + 2\mu \pmb{\epsilon}^{-}.
\]

Due to the crack tension-compression asymmetry, the phase-field degradation function is exclusively applied to the positive part of the strain energy. Therefore, the stress-strain relationship is nonlinear, and the material tangent modulus is written as

\[
\frac{\partial \pmb{\sigma}}{\partial \pmb{\epsilon}} = [g(d) + k]\frac{\partial \pmb{\sigma}^{+}}{\partial \pmb{\epsilon}} + \frac{\partial \pmb{\sigma}^{-}}{\partial \pmb{\epsilon}} = [g(d) + k]\left[\lambda H(\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes \mathbf{I} + 2\mu \mathbb{P}^{+}\right] + \left[\lambda H(-\mathrm{tr}\pmb{\epsilon})\mathbf{I}\otimes \mathbf{I} + 2\mu \mathbb{P}^{-}\right].
\]

Inside the above tangent modulus, the two fourth-order projection tensors \(\mathbb{P}^{+}\) and \(\mathbb{P}^{-}\) are defined as

\[
\mathbb{P}^{+} = \frac{\partial \pmb{\epsilon}^{+}}{\partial \pmb{\epsilon}}, \quad \mathbb{P}^{-} = \frac{\partial \pmb{\epsilon}^{-}}{\partial \pmb{\epsilon}},
\]

the specific expressions of which can be found in [26].

Using the directional derivative, the first variation of the energy functional is written as

\[
\begin{array}{rl} 
\delta \Pi (\pmb {u},d) &= D_{(\delta \pmb {u},\delta d)}\Pi (\pmb {u},d) = \left.\frac{\mathrm{d}}{\mathrm{d}\epsilon}\right|_{\epsilon = 0}\Pi (\pmb {u} + \epsilon \delta \pmb {u},d + \epsilon \delta d) \\ 
&= \int_{\Omega}\left(\frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial \pmb{\epsilon}}:\pmb{\epsilon}(\delta \pmb {u}) + \frac{\partial\psi(\pmb{\epsilon}(\pmb{u}),d)}{\partial d}\delta d\right)\mathrm{d}\Omega \\ 
&\quad +\int_{\Omega}\frac{g_{c}}{l}\left(d\delta d + l^{2}\nabla d\cdot \nabla \delta d\right)\mathrm{d}\Omega -\int_{\Omega}\pmb {b}\cdot \delta \pmb {u}\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot \delta \pmb {u}\mathrm{d}\Gamma \\ 
&= (\nabla^{(s)}\delta \pmb {u},\pmb {\sigma}) - (\delta \pmb {u},\pmb {b}) - (\delta \pmb {u},\pmb{t})_{\Gamma_{t}} + (\delta d,\frac{g_{c}}{l} d) + (\nabla \delta d,g_{c}l\nabla d) + (\delta d,g^{\prime}(d)\psi^{+}(\pmb{\epsilon})). 
\end{array} \quad (7)
\]

Therefore, the weak form of the phase-field formulation is expressed as

\[
\left\{ \begin{array}{ll} r_{\pmb{u}}(\pmb{u},d) = (\nabla^{(s)}\delta \pmb{u},\pmb{\sigma}) - (\delta \pmb{u},\pmb{b}) - (\delta \pmb{u},\pmb{t})_{\Gamma_{t}} = 0, \\ r_{d}(\pmb{u},d) = (\delta d,\frac{g_{c}}{l} d) + (\nabla \delta d,g_{c}l\nabla d) + (\delta d,g^{\prime}(d)\psi^{+}(\pmb{\epsilon})) = 0, \end{array} \right.
\]

subject to the phase-field irreversibility condition.

### 2.2. Finite element discretization

Let \(\mathbf{N}_{u_{A}}\) represent the vector-valued shape function for the displacement field \(\pmb{u}\) at node \(A\) and \(N_{d_{A}}\) represent the scalar-valued shape function for the phase-field \(d\) at node \(A\). Based on the shape functions \(\{N_{u_{A}},N_{d_{A}}\}\) and the nodal values \(\{\pmb{u}_{A},d_{A}\}\), the displacement field and the phase-field can be expressed as

\[
\pmb {u} = \pmb{N}_{u_{A}}\pmb{u}_{A}, \quad \text{and} \quad d = N_{d_{A}}d_{A}.
\]

Correspondingly, the displacement variation and the phase-field variation can be expressed as

\[
\delta \pmb {u} = \pmb{N}_{u_{A}}\delta \pmb{u}_{A}, \quad \text{and} \quad \delta d = N_{d_{A}}\delta d_{A},
\]

where the Einstein summation is used. Plug the above expressions into the total energy defined in Eq. (1), the objective function for the inequality-constrained minimization after the finite element discretization is expressed as

\[
\begin{array}{rl} 
\Pi (\pmb {u}_{A},d_{A}) = &\int_{\Omega}\psi \left(\pmb{\epsilon}(\pmb {N}_{u_{A}}\pmb {u}_{A}),N_{d_{A}}d_{A}\right)\mathrm{d}\Omega \\ 
&+\int_{\Omega}\frac{g_{c}}{2l}\left((N_{d_{A}}d_{A})^{2} + l^{2}(\nabla N_{d_{A}}d_{A})\cdot (\nabla N_{d_{A}}d_{A})\right)\mathrm{d}\Omega \\ 
&-\int_{\Omega}\pmb {b}\cdot (\pmb {N}_{u_{A}}\pmb {u}_{A})\mathrm{d}\Omega -\int_{\partial \Omega}\pmb {t}\cdot (\pmb {N}_{u_{A}}\pmb {u}_{A})\mathrm{d}\Gamma, 
\end{array} \quad (8)
\]

The gradient of the discretized energy functional is derived as

\[
\begin{array}{rl} 
\pmb{r} = \nabla \Pi = (r_{\pmb{u}},r_{d})^{\mathrm{T}},\\ 
r_{\pmb{u}_{A}} = \frac{\partial\Pi}{\partial \pmb{u}_{A}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},\pmb{\sigma}\right) - \left(\pmb{N}_{u_{A}},\pmb{b}\right) - \left(\pmb{N}_{u_{A}},\pmb{t}\right)_{\Gamma_{t}},\\ 
r_{d_{A}} = \frac{\partial\Pi}{\partial d_{A}} = \left(N_{d_{A}},\frac{g_{c}}{l} d + g^{\prime}(d)\psi^{+}\right) + (\nabla N_{d_{A}},g_{c}l\nabla d). 
\end{array} \quad (9)
\]

The Hessian matrix of the total energy after the finite element discretization is

\[
\mathbf{K} = \nabla^{2}\Pi = \left[ \begin{array}{ll}\mathbf{K}_{uu} & \mathbf{K}_{ud}\\ \mathbf{K}_{du} & \mathbf{K}_{dd} \end{array} \right],
\]

where

\[
\begin{array}{rl} 
\mathbf{K}_{u_{A}u_{B}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},\frac{\partial \pmb{\sigma}}{\partial \pmb{\epsilon}}:\nabla^{(s)}\pmb{N}_{u_{B}}\right),\qquad \mathbf{K}_{u_{A}d_{B}} = \left(\nabla^{(s)}\pmb{N}_{u_{A}},g^{\prime}(d)\pmb{\sigma}^{+}N_{d_{B}}\right),\\ 
\mathbf{K}_{d_{A}u_{B}} = \left(N_{d_{A}},g^{\prime}(d)\pmb{\sigma}^{+}:\nabla^{(s)}\pmb{N}_{u_{B}}\right),\qquad \mathbf{K}_{d_{A}d_{B}} = \left(N_{d_{A}},\left(\frac{g_{c}}{l} +g^{\prime \prime}(d)\psi^{+}\right)N_{d_{B}}\right) + \left(\nabla N_{d_{A}},g_{c}l\nabla N_{d_{B}}\right). 
\end{array} \quad (10)
\]

In the above block matrices, the subscripts \(A\) and \(B\) represent the finite element nodal number, respectively. During a typical time (load) step \([t_n,t_{n + 1}]\), the phase-field irreversibility condition can be expressed as the following box constraints applied at individual finite element node, that is,

\[
d_A^{(n)}\leq d_A\leq 1,
\]

where the subscript \(A\) represents the finite element nodal number, and \(d_A^{(n)}\) is known at the beginning of the current time step. Lastly, a diagonal block matrix is defined as,

\[
\hat{\mathbf{K}} = \left[ \begin{array}{cc}\mathbf{K}_{uu} & \mathbf{0}\\ \mathbf{0} & \mathbf{K}_{dd} \end{array} \right], \quad (11)
\]

which is needed in the monolithic scheme presented in Section 3.

## 3. Gradient projection based monolithic scheme

This section is at the core of this work and presents the theoretical, algorithmic, and implementation details of the proposed phase-field monolithic scheme. This scheme uses the gradient projection method to enforce the box constraints and the limited-memory BFGS method to overcome the convergence difficulties associated with non-convexity of the discretized energy functional. In the field of numerical optimization, the above approach is called the L-BFGS-B method [41,45] and has three major components, including the compact representation of the limited-memory BFGS matrix, the calculation of the generalized Cauchy point to detect the active set of the box constraints, and the formulation and solving techniques of the resulted subspace minimization. Section 3.1 firstly provides an overview of the entire monolithic scheme and how each component is connected with one another. Then, Sections 3.2 to 3.4 provide the corresponding details for each algorithmic component.

### 3.1. Algorithm overview

For the quasi-static phase-field fracture propagation problem, the pseudo time \(t\) is used to represent the actual load step. Inside the time (load) step \([t_n,t_{n + 1}]\), the solution \((\cdot)^{(n)}\) at the beginning is known. The goal is to determine the solution \((\cdot)^{(n + 1)}\) at the end of the current time step by solving the following constrained minimization after the finite element discretization,

\[
\left(\pmb{u}_A,d_A\right) = \arg \min \Pi (\pmb {u}_A,d_A) \quad (12)
\]

subject to the inequality constraints

\[
d_A^{(n)}\leq d_A\leq 1 \quad (13)
\]

and the linear constraints in the following form

\[
\left[\pmb {u}_A,d_A\right]^{\mathrm{T}} = \mathbf{C}\left[\pmb {u}_A,d_A\right]^{\mathrm{T}} + \pmb {k}. \quad (14)
\]

In the above constrained minimization problem, \(\{\pmb {u}_A,d_A\}\) represents the finite element nodal solution of the displacement field and the phase-field. The objective function is shown in Eq. (8), and the gradient of the objective function is shown in Eq. (9). The linear constraints shown in Eq. (14) include the essential (Dirichlet) boundary conditions imposed on the domain boundary and the hanging-node constraints if the adaptive mesh refinement is adopted. The constraint coefficient matrix \(\mathbf{C}\) and the inhomogeneous constraints \(\mathbf{k}\) can be obtained from the essential (Dirichlet) boundary conditions and the adaptive mesh refinement process, for instance, see [26].

To overcome the convergence difficulties arising from the non-convexity of the objective function in Eq. (12), the BFGS method is adopted. As a type of quasi-Newton method, the BFGS method is originally based on the idea from Davidon [46,47] and further developed by Broyden, Fletcher, Goldfarb, and Shanno [48-51]. The basic idea is to construct a series of quadratic models to approximate the original objective function,

\[
m_{k}(\pmb{x}) = \Pi_{k} + \pmb{r}_{k}^{\mathrm{T}}(\pmb{x} - \pmb{x}_{k}) + \frac{1}{2} (\pmb{x} - \pmb{x}_{k})^{\mathrm{T}}\mathbf{B}_{k}(\pmb{x} - \pmb{x}_{k}) \quad (15)
\]

where the subscript \(k\) represents the \(k\) th iteration, \(\mathbf{x}\) represents the unknowns \(\{\pmb{u},d\}\), and \(\pmb{r}_k = \nabla \Pi_k\) represents the gradient of the objective function. The BFGS matrix \(\mathbf{B}_k\) is constructed from two rank-one updates on top of the previous matrix \(\mathbf{B}_{k - 1}\) and is positive definite [45], regardless whether the original objective function is convex or not. However, in the context of the finite element method, the BFGS matrix is fully dense and loses the underlying sparsity pattern associated with the finite element spatial discretization. Therefore, the conventional BFGS method cannot be directly applied to large-scale finite element simulations due to the memory required to store the fully dense matrix. This limitation motivated the development of the “limited-memory” BFGS (L-BFGS) method [27-29], which is further discussed in Section 3.2.

In order to consider the phase-field irreversibility condition, the quadratic model in Eq. (15) needs to satisfy the box constraints shown in Eq. (13). Let \(\mathbf{lb}_A = d_A^{(n)}\) and \(\mathbf{ub}_A = 1.0\) represent the lower and upper bounds of the phase-field at node \(A\) during time step \([t_n,t_{n + 1}]\). The gradient projection method [39-42], which is a type of active set method, can be used to effectively and efficiently enforce the box constraints. Define the projection operator \(\mathrm{Proj}_c(\cdot)\) as

\[
\mathrm{Proj}_c(x_i,\mathrm{lb}_i,\mathrm{ub}_i) = \left\{ \begin{array}{ll}\mathrm{lb}_i & \mathrm{if}\; x_i < \mathrm{lb}_i,\\ x_i & \mathrm{if}\; \mathrm{lb}_i \leq x_i \leq \mathrm{ub}_i,\\ \mathrm{ub}_i & \mathrm{if}\; x_i > \mathrm{ub}_i, \end{array} \right. \quad (16)
\]

where the feasible region \(C\) formed by the phase-field inequality constraints is obviously convex. Based on the above projection operator, at the \(k\) th iteration, the steepest descent direction \(\pmb {r}_k = \nabla \Pi_k\) can be projected onto the feasible region with a negligible computational cost. The projected gradient forms a piecewise linear path, for instance, see Fig. 1(b), and can be written as

\[
\pmb {x}(t) = \mathrm{Proj}_c(\pmb {x}_k - t\pmb {r}_k,\mathbf{lb},\mathbf{ub}),\quad t > 0, \quad (17)
\]

in which the projection operator is applied to the vector \(\pmb{x}_k - t\pmb{r}_k\) component wise.

Following the piecewise linear path formed by the projected gradient, the quadratic model in the BFGS method, as shown in Eq. (15), is transformed into a univariate piecewise quadratic model,

\[
p_k(t) = m_k(\pmb{x}(t)) = \Pi_k + \pmb{r}_k^\mathrm{T}(\pmb{x}(t) - \pmb{x}_k) + \frac{1}{2} (\pmb{x}(t) - \pmb{x}_k)^\mathrm{T}\mathbf{B}_k(\pmb{x}(t) - \pmb{x}_k). \quad (18)
\]

The first local minimizer \(\pmb{x}^c\) of the above model is called the generalized Cauchy point, the calculation of which is presented in Section 3.3.

For an arbitrary solution state \(\pmb{x}\) the active set \(\mathcal{A}(\pmb {x})\) of the inequality constraints shown in Eq. (13) include those components whose values are at either the lower bound or the upper bound. Note that the active set \(\mathcal{A}(\pmb {x})\) actually only contains component index, that is,

\[
\mathcal{A}(\pmb {x}) = \{i:x_i = \mathrm{lb}_i\} \cup \{i:x_i = \mathrm{ub}_i\}.
\]

Based on the generalized Cauchy point \(\pmb{x}^c\), the variables in the active set \(\mathcal{A}_k(\pmb{x}^c)\) (located at the boundary of the box constraints) are held fixed. The solution \(\pmb{x}^*\) is obtained by minimizing the quadratic model shown in Eq. (15) in the subspace of the free variables (located inside the box constraints) subject to the box constraints, that is,

\[
\pmb{x}^* = \arg \min m_k(\pmb {x})
\]

subject to

\[
x_{i} = x_{i}^{c},\;\forall i\in \mathcal{A}_{k}(\pmb{x}^{c})\quad \text{and}\quad \mathrm{lb}_{i}\leq x_{i}\leq \mathrm{ub}_{i},\;\forall i\notin \mathcal{A}_{k}(\pmb{x}^{c}).
\]

Firstly ignoring the inequality constraints on the free variables, the above minimization can be solved using either the primal approach on the subspace of the free variables or the dual approach enforcing the active constraints via the Lagrange multipliers. Then, \(\pmb{x}^*\) is obtained by truncating the above solution to satisfy the box constraints. Section 3.4 discusses in detail the numerical techniques to solve the subspace minimization problem.

Based on the solution \(\pmb{x}^*\) of the above minimization, a search direction \(\pmb{p}_k\) is defined as

\[
\pmb{p}_k = \pmb{x}^* - \pmb{x}_k.
\]

The search direction \(\pmb{p}_k\) is a descent direction of the objective function \(\Pi (\pmb{x})\), see [41] for the proof. The new solution \(\pmb{x}_{k + 1}\) is obtained as

\[
\pmb{x}_{k + 1} = \pmb{x}_k + \alpha_k\pmb{p}_k.
\]

The positive step length \(\alpha_k > 0\) is determined using the line search based on the strong Wolfe conditions, which include the sufficient decrease condition

\[
\Pi_{k + 1} = \Pi (\pmb{x}_k + \alpha_k\pmb {p}_k) \leq \Pi_k + c_1\alpha_k\pmb {p}_k^\mathrm{T}\pmb{r}_k
\]

and the curvature condition

\[
\left|\pmb{r}_{k + 1}^\mathrm{T}\pmb {p}_k\right| = \left|\pmb {r}(\pmb{x}_k + \alpha_k\pmb {p}_k)^\mathrm{T}\pmb {p}_k\right| \leq c_2\left|\pmb{r}_k^\mathrm{T}\pmb {p}_k\right|.
\]

To find a step length \(\alpha_k\) satisfying the strong Wolfe conditions, this work uses the implementation proposed by More and Thuente [52] with the parameters \(c_1 = 10^{-4}\) and \(c_2 = 0.9\).

Once the solution \(\pmb{x}_{k + 1}\) is obtained, a new iteration starts by updating the L-BFGS matrix \(\mathbf{B}_{k + 1}\), calculating the generalized Cauchy point \(\pmb{x}^c\) from the univariate piecewise quadratic model \(p_{k + 1}(t)\), updating the active set \(\mathcal{A}_{k + 1}(\pmb{x}^c)\), solving the subspace minimization to obtain the search direction \(\pmb{p}_{k + 1}\), and performing the line search based on the strong Wolfe conditions. The nonlinear iterations stop when the following convergence criteria are met simultaneously:

1.  The active set remains unchanged between two consecutive iterations:
    \[
    \mathcal{A}_{k + 1}(\pmb{x}^c) = \mathcal{A}_k(\pmb{x}^c).
    \]
2.  The \(l_{2}\)-norm of the projected gradient is smaller than the prescribed tolerance:
    \[
    \| \mathrm{Proj}_C(\mathbf{x}_{k + 1} - \mathbf{r}_{k + 1},\mathbf{lb},\mathbf{ub}) - \mathbf{x}_{k + 1}\|_{2} < \mathrm{tol}.
    \]
3.  The \(l_{2}\)-norm of the solution increment is smaller than the prescribed tolerance:
    \[
    \| \Delta \mathbf{x}_{k + 1}\|_{2} = \| \mathbf{x}_{k + 1} - \mathbf{x}_{k}\|_{2} < \mathrm{tol}.
    \]

Recall that the projection operator is defined in Eq. (16). For the \(i\) th component of the projected gradient,

\[
\mathrm{Proj}_C(x_i - r_i,\mathrm{lb}_i,\mathrm{ub}_i) - x_i = \left\{ \begin{array}{ll}\mathrm{lb}_i - x_i & \mathrm{if}\; x_i - r_i < \mathrm{lb}_i,\\ - r_i & \mathrm{if}\; \mathrm{lb}_i \leq x_i - r_i \leq \mathrm{ub}_i,\\ \mathrm{ub}_i - x_i & \mathrm{if}\; x_i - r_i > \mathrm{ub}_i. \end{array} \right. \quad (19)
\]

The second convergence criterion based on the projected gradient means that if the \(i\) th component of the projected gradient is within the feasible region (the box constraints), the \(l_{2}\)-norm uses the gradient value \(r_i\) itself. If the \(i\) th component of the projected gradient is outside the feasible region, the \(l_{2}\)-norm uses the distance between the solution component \(x_i\) and the boundary of the box constraints \(\mathrm{lb}_i\) or \(\mathrm{ub}_i\). This distance measures how well the inequality constraint is enforced.

**Comment 1.** The linear constraints shown in Eq. (14) include the essential (Dirichlet) boundary conditions imposed on the domain boundary and the hanging-node constraints if the adaptive mesh refinement is applied. Let \(\mathbf{x} = \{\mathbf{u},d\}\) represent the unknowns. Then, the linear constraints can be written in the following form

\[
\mathbf{x} = \mathbf{C}\mathbf{x} + \mathbf{k}.
\]

For a linear system \(\mathbf{A}\mathbf{x} = \mathbf{b}\) with the above set of constraints, we can instead solve the following modified linear system [53]

\[
\left(\mathbf{C}^{\mathrm{T}}\mathbf{A}\mathbf{C} + \mathbf{I}_{d_{c}}\right)\hat{\mathbf{x}} = \mathbf{C}^{\mathrm{T}}(\mathbf{b} - \mathbf{A}\mathbf{k})
\]

and then recover the true solution \(\mathbf{x}\) as

\[
\mathbf{x} = \mathbf{C}\hat{\mathbf{x}} + \mathbf{k}.
\]

In the modified linear system, the matrix \(\mathbf{I}_{d_c}\) is defined as
\[
(\mathbf{I}_{d_c})_{ii} = \begin{cases} 1 & \text{if DoF } i \text{ is constrained}, \\ 0 & \text{otherwise}, \end{cases}
\]
where \(\mathcal{T}\) represents the set of the degrees of freedom at the constrained nodes, including the hanging nodes and the nodes prescribed with essential boundary conditions. For a detailed discussion about how to apply the above linear constraints, see [26].

### 3.2. Compact representation of limited-memory BFGS matrix

The detailed approach of the standard BFGS method can be found in many textbooks in the field of numerical optimization, for instance, see the classical textbook by Nocedal and Wright [45]. Let \(\mathbf{x} = \{\mathbf{u}_A,d_A\}\) represent the solutions of the displacement field and the phase-field after the finite element discretization. Recall that \(\Pi (\mathbf{x})\) represents the objective function (in this case the discretized total energy) to be minimized, and \(\mathbf{r} = \nabla \Pi (\mathbf{x})\) represents its gradient. At the \(k\) th BFGS iteration, the following two vectors are defined:

\[
\mathbf{s}_k = \mathbf{x}_{k + 1} - \mathbf{x}_k,\quad \mathbf{y}_k = \mathbf{r}_{k + 1} - \mathbf{r}_k,
\]

and the BFGS matrix is updated as

\[
\mathbf{B}_{k + 1} = \mathbf{B}_k - \frac{\mathbf{B}_k\mathbf{s}_k\mathbf{s}_k^{\mathrm{T}}\mathbf{B}_k}{\mathbf{s}_k^{\mathrm{T}}\mathbf{B}_k\mathbf{s}_k} + \frac{\mathbf{y}_k\mathbf{y}_k^{\mathrm{T}}}{\mathbf{y}_k^{\mathrm{T}}\mathbf{s}_k}.
\]

If the old BFGS matrix \(\mathbf{B}_k\) is positive definite and the vector pairs \(\mathbf{s}_k,\mathbf{y}_k\) satisfy the curvature condition,

\[
\mathbf{s}_k^{\mathrm{T}}\mathbf{y}_k > 0,
\]

then the updated BFGS matrix \(\mathbf{B}_{k + 1}\) is also positive definite [45]. Since \(\mathbf{s}_k\) and \(\mathbf{y}_k\) are two vectors with \(n\) components, where \(n\) represents the total number of degrees of freedom of the discretized system, each of the two terms \(\mathbf{B}_k\mathbf{s}_k\mathbf{s}_k^{\mathrm{T}}\mathbf{B}_k\) and \(\mathbf{y}_k\mathbf{y}_k^{\mathrm{T}}\) is a rank-one fully dense \(n\times n\) matrix, respectively. Therefore, in the BFGS matrix update, regardless whether \(\mathbf{B}_k\) is a sparse matrix or not, the updated matrix \(\mathbf{B}_{k + 1}\) is a fully dense matrix. Storing fully dense matrices is impractical due to the required memory. As a result, the conventional BFGS method cannot be directly used in large-scale finite element simulations. This limitation is the motivation for the development of the “limited-memory” feature for the BFGS method in the literature of numerical optimization, see [27-29].

In this work, the compact representation of the limited-memory BFGS (L-BFGS) matrix, originally proposed by Byrd et al. [29], is adopted to avoid the storage of a fully dense matrix. At the \(k\) th iteration, the BFGS matrix is essentially composed of an initial matrix \(\mathbf{B}^0_k\) and a series of vector-pairs \(\{\mathbf{s}_i, \mathbf{y}_i\}_{i=1}^{k-1}\). These vector-pairs contain the curvature information of the objective function. According to Nocedal and Wright [45], “the main idea of this method (L-BFGS) is to use curvature information from only the most recent iterations to construct the Hessian approximation. Curvature information from earlier iterations, which is less likely to be relevant to the actual behavior of the Hessian at the current iteration, is discarded in the interest of saving storage”. Instead of storing the entire history of the vector pairs, only the most recent \(m\) vector-pairs \(\{\mathbf{s}_i, \mathbf{y}_i\}_{i=k-m}^{k-1}\) need to be stored, where \(m\) is a small number. Obviously, if \(k < m\), all the vector-pairs are stored. According to the compact representation proposed by Byrd et al. [29], the following two \(n\times m\) correction matrices are defined as

\[
\mathbf{S}_k = [\mathbf{s}_{k-m} \cdots \mathbf{s}_{k-1}]
\]
and
\[
\mathbf{Y}_k = [\mathbf{y}_{k-m} \cdots \mathbf{y}_{k-1}].
\]

Then, the L-BFGS matrix can be represented in the following compact form,

\[
\mathbf{B}_k = \mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k, \quad (20)
\]

where

\[
\mathbf{W}_k = [\mathbf{Y}_k \quad \mathbf{B}^0_k \mathbf{S}_k] \in \mathbb{R}^{n\times 2m}
\]
and

\[
\mathbf{M}_k = \begin{bmatrix} -\mathbf{D}_k & \mathbf{L}^{\mathrm{T}}_k \\ \mathbf{L}_k & \mathbf{S}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{S}_k \end{bmatrix}^{-1} \in \mathbb{R}^{2m\times 2m}.
\]

The matrix \(\mathbf{M}_k\) is a \(2m\times 2m\) matrix with the block matrices defined as

\[
\mathbf{D}_k = \mathrm{diag}\{\mathbf{s}^{\mathrm{T}}_{k-m}\mathbf{y}_{k-m}, \ldots, \mathbf{s}^{\mathrm{T}}_{k-1}\mathbf{y}_{k-1}\} \in \mathbb{R}^{m\times m}
\]
and
\[
\mathbf{L}_k \in \mathbb{R}^{m\times m}, \quad (\mathbf{L}_k)_{i,j} = \begin{cases} \mathbf{s}^{\mathrm{T}}_{i+k-m-1}\mathbf{y}_{j+k-m-1} & \text{if } i > j, \\ 0 & \text{if } i \leq j. \end{cases}
\]

At the \(k\) th iteration of the phase-field monolithic scheme, the initial BFGS matrix \(\mathbf{B}^0_k\) takes the diagonal block matrix defined in Eq. (11), that is,

\[
\mathbf{B}^0_k = \hat{\mathbf{K}}^{(k)} = \begin{bmatrix} \mathbf{K}_{uu} & \mathbf{0} \\ \mathbf{0} & \mathbf{K}_{dd} \end{bmatrix}^{(k)}. \quad (21)
\]

The matrices \(\mathbf{K}_{uu}\) and \(\mathbf{K}_{dd}\) are both assembled from the finite element procedure, and therefore, are sparse. Also, both \(\mathbf{K}_{uu}\) and \(\mathbf{K}_{dd}\) are positive definite. Therefore, \(\mathbf{B}^0_k\) is sparse and positive definite, and can be stored in the memory without any issue. On the other hand, the second term \(\mathbf{W}_k\mathbf{M}_k\mathbf{W}^{\mathrm{T}}_k\) in Eq. (20) results a fully dense \(n\times n\) matrix and is never calculated or stored as a whole. Rather, only the matrices \(\mathbf{W}_k \in \mathbb{R}^{n\times 2m}\) and \(\mathbf{M}_k \in \mathbb{R}^{2m\times 2m}\) are stored individually. Recall that \(m \ll n\) is only a small number and typically ranges from 10 to 40 [45].

In summary, according to the compact representation of the L-BFGS matrix, at the \(k\) th iteration, only the following matrices need to be stored in the memory, including the sparse and positive definite matrix \(\mathbf{B}^0_k = \hat{\mathbf{K}}\), the \(2m\times 2m\) matrix \(\mathbf{M}_k\), and the \(n\times 2m\) matrix \(\mathbf{W}_k\). At the end of the \(k\) th iteration, the oldest vector-pair \(\{\mathbf{s}_{k-m}, \mathbf{y}_{k-m}\}\) is discarded, and the newest pair \(\{\mathbf{s}_k, \mathbf{y}_k\}\) is calculated and added into the list.

**Comment 2.** In the numerical implementation of the compact representation of L-BFGS matrix, the two \(n\times m\) correction matrices \(\mathbf{Y}\) and \(\mathbf{S}\) do not need to be recalculated from scratch during each iteration. Since only the oldest vector-pair is removed and the newly formed vector-pair is added, the rest of the vector-pairs remain unchanged. The best data structure to store \(\mathbf{Y}\) or \(\mathbf{S}\) is not as a matrix, but as a linked list, in which the complexity of adding or removing an element into/from the memory is constant \(\mathcal{O}(1)\).

### 3.3. Generalized Cauchy point

Once the L-BFGS matrix is formed using the compact representation, the next step is to locate the generalized Cauchy point, which is defined as the first local minimizer along the piecewise linear path formed by the projected gradient. The calculation of the generalized Cauchy point mainly follows the steps detailed in [41]. However, several places are modified to consider the impact of the initial BFGS matrix \(\mathbf{B}^0 = \hat{\mathbf{K}}\) on the calculation. At the \(k\) th nonlinear iteration, we use \(\mathbf{x}^0\) to represent the solution at the beginning of the current iteration. For notation convenience, we drop the iteration index \(k\) in all the vectors and matrices such that \(\mathbf{r}= \mathbf{r}_k\) and \(\mathbf{B}= \mathbf{B}_k\). Instead, the subscript index is used to represent the component of a vector. For instance, \(r_i\) represents the \(i\) th component of the gradient vector \(\mathbf{r}\). Also, the concept of “break points” is introduced to represent the intersections between the projected gradient and the box constraints. The schematic in Fig. 2 illustrates the path of the projected gradient and the break points. For the \(i\) th coordinate direction, the break point along the projected gradient can be represented by a positive scalar as

\[
t_i = \begin{cases} (x^0_i - \mathrm{ub}_i)/r_i & \text{if } r_i < 0, \\ (x^0_i - \mathrm{lb}_i)/r_i & \text{if } r_i > 0, \\ +\infty & \text{if } r_i = 0. \end{cases}
\]

Note that if there is no upper or lower bound on the \(i\) th direction, then set \(\mathrm{ub}_i = +\infty\) or \(\mathrm{lb}_i = -\infty\). After all the break points are computed for each direction \(\{t_i\}_{i = 1}^n\), they need to be reordered to form an increasing order set \(\{t^{(j)}: t^{(j)} \leq t^{(j+1)}\}_{j=1}^n\). Note that the subscript \(i\) represents the \(i\) th vector component, while the superscript \(j\) represents the \(j\) th element in the ordered set. Since the break points need to be sorted multiple times during each iteration, it is important to efficiently perform the sort operation using the appropriate data structure and algorithm, for instance, the heap sort algorithm with a complexity of \(\mathcal{O}(n\log n)\). Using the break points, the piecewise linear path formed by the projected gradient \(\mathrm{Proj}_C(\mathbf{x}^0 - t\mathbf{r},\mathbf{lb},\mathbf{ub})\) can be expressed in the following component form

\[
x_i(t) = x_i^{0} - \min \{t,t_i\} r_i.
\]

Inside an interval \([t^{(j-1)}, t^{(j)}]\) along the piecewise linear path formed by the projected gradient, the solution vector can be expressed as a function of the scalar \(t\)

\[
\mathbf{x}(t) = \mathbf{x}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)},
\]

where \(\mathbf{x}^{(j-1)} = \mathbf{x}(t^{(j-1)})\), \(\Delta t = t - t^{(j-1)}\), and \(\mathbf{d}^{(j-1)}\) is the direction of the projected gradient in that interval. Let \(\mathbf{z}^{(j-1)} = \mathbf{x}^{(j-1)} - \mathbf{x}^0\), and the univariate quadratic model on this interval can be expressed as

\[
\begin{array}{rl} 
p(t) &= m(\mathbf{x}(t)) = \Pi + \mathbf{r}^{\mathrm{T}}(\mathbf{x}(t) - \mathbf{x}^0) + \frac{1}{2} (\mathbf{x}(t) - \mathbf{x}^0)^{\mathrm{T}}\mathbf{B}(\mathbf{x}(t) - \mathbf{x}^0) \\
&= \Pi + \mathbf{r}^{\mathrm{T}}(\mathbf{z}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)}) + \frac{1}{2} (\mathbf{z}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)})^{\mathrm{T}}\mathbf{B}(\mathbf{z}^{(j-1)} + \Delta t \mathbf{d}^{(j-1)}) \\
&= \left(\Pi + \mathbf{r}^{\mathrm{T}}\mathbf{z}^{(j-1)} + \frac{1}{2} \mathbf{z}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j-1)}\right) \\
&\quad + \left(\mathbf{r}^{\mathrm{T}}\mathbf{d}^{(j-1)} + \mathbf{d}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j-1)}\right)\Delta t + \frac{1}{2}\left(\mathbf{d}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j-1)}\right)\Delta t^{2} \\
&= f_{j-1} + f_{j-1}^{\prime}\Delta t + \frac{1}{2} f_{j-1}^{\prime \prime}\Delta t^{2} = \hat{p} (\Delta t), 
\end{array} \quad (22)
\]

where \(f_{j-1}, f_{j-1}^{\prime}\) and \(f_{j-1}^{\prime \prime}\) represent the coefficients of the quadratic polynomial \(\hat{p}(\Delta t)\).

As discussed in Section 3.2, the updated L-BFGS matrix \(\mathbf{B}\) is positive definite if the vector-pairs \(\{\mathbf{s}_i,\mathbf{y}_i\}_{i=k-m}^{k-1}\) satisfy the curvature condition shown in Eq. (18). Therefore, the coefficient \(f_{j-1}^{\prime \prime}\) is positive, that is,

\[
f_{j-1}^{\prime \prime} = \mathbf{d}^{(j-1)}{}^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j-1)} > 0.
\]

The quadratic polynomial \(\hat{p} (\Delta t)\) reaches to its minimum value at the critical point \(\Delta t^*\)

\[
\hat{p}'(\Delta t) = 0 \Rightarrow \Delta t^* = -f_{j-1}' / f_{j-1}^{\prime \prime}.
\]

Depending on the value of \(\Delta t^*\), there are the following three scenarios:

1.  If \(\Delta t^* \in [0, t^{(j)} - t^{(j-1)})\), the generalized Cauchy point \(\mathbf{x}^c = \mathbf{x}(t^{(j-1)} + \Delta t^*)\).
2.  If \(\Delta t^* < 0\), the generalized Cauchy point \(\mathbf{x}^c = \mathbf{x}(t^{(j-1)})\).
3.  If \(\Delta t^* \geq t^{(j)} - t^{(j-1)}\), the search of the generalized Cauchy point moves onto the next interval (line segment) of the projected gradient \([t^{(j)}, t^{(j+1)}]\).

For the first two scenarios, the set of active constraints include those components whose break point \(t\)-value is no larger than \(t^{(j-1)}\), that is,
\[
\mathcal{A}(\mathbf{x}^c) = \{i : t_i \leq t^{(j-1)}, i = 1, \ldots, n\}.
\]

During the next interval \([t^{(j)}, t^{(j+1)}]\), it is not necessary to calculate the values of \(f'_j\) and \(f''_j\) from scratch. Rather, they can be obtained from updating the old values of \(f'_{j-1}\) and \(f''_{j-1}\), which could significantly lower the computational cost. Let
\[
\Delta t^{(j-1)} = t^{(j)} - t^{(j-1)}, \quad \mathbf{x}^{(j)} = \mathbf{x}^{(j-1)} + \Delta t^{(j-1)}\mathbf{d}^{(j-1)}.
\]

Moving from the previous interval \([t^{(j-1)}, t^{(j)}]\) to the new interval \([t^{(j)}, t^{(j+1)}]\), assume that the constraint associated with the \(b\) th vector component becomes active, that is, \(t_b = t^{(j)}\). The direction of the projected gradient in the current interval becomes
\[
\mathbf{d}^{(j)} = \mathbf{d}^{(j-1)} + r_b \mathbf{e}_b,
\]
where \(\mathbf{e}_b\) is the \(b\) th unit vector, and \(r_b\) is the \(b\) th component of the gradient vector \(\mathbf{r}\). Similar to Eq. (22), the new values \(f'_j\) and \(f''_j\) are written as

\[
f'_j = \mathbf{r}^{\mathrm{T}}\mathbf{d}^{(j)} + \mathbf{d}^{(j)}{}^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j)} = f'_{j-1} + \Delta t^{(j-1)} f''_{j-1} + r^2_b + r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{z}^{(j)} \quad (23)
\]
and
\[
f''_j = \mathbf{d}^{(j)}{}^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j)} = f''_{j-1} + 2r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{d}^{(j-1)} + r^2_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}\mathbf{e}_b. \quad (24)
\]

Based on the compact representation of the L-BFGS matrix, as shown in Eq. (20),
\[
\mathbf{B} = \mathbf{B}^0 - \mathbf{W}\mathbf{M}\mathbf{W}^{\mathrm{T}}.
\]
Let
\[
\mathbf{w}_b = \mathbf{W}^{\mathrm{T}}\mathbf{e}_b, \quad \mathbf{p}^{(j)} = \mathbf{W}^{\mathrm{T}}\mathbf{d}^{(j)} = \mathbf{W}^{\mathrm{T}}(\mathbf{d}^{(j-1)} + r_b \mathbf{e}_b) = \mathbf{p}^{(j-1)} + r_b \mathbf{w}_b, \quad \mathbf{c}^{(j)} = \mathbf{W}^{\mathrm{T}}\mathbf{z}^{(j)} = \mathbf{W}^{\mathrm{T}}(\mathbf{z}^{(j-1)} + \Delta t^{(j-1)}\mathbf{d}^{(j-1)}) = \mathbf{c}^{(j-1)} + \Delta t^{(j-1)}\mathbf{p}^{(j-1)}.
\]

Plug all the above relationships into Eqs. (23) and (24), the coefficients of the quadratic polynomial \(\hat{p}(\Delta t)\) for the interval \([t^{(j)}, t^{(j+1)}]\) are written as
\[
f'_j = f'_{j-1} + \Delta t^{(j-1)} f''_{j-1} + r^2_b + r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{z}^{(j)} - r_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{c}^{(j)}
\]
and
\[
f''_j = f''_{j-1} + 2r_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{d}^{(j-1)} - 2r_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{p}^{(j-1)} + r^2_b \mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{e}_b - r^2_b \mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{w}_b.
\]

Similar to the steps in the previous interval, the critical point of the quadratic polynomial \(\hat{p}(\Delta t)\) is calculated as \(\Delta t^* = -f'_j / f''_j\). If \(\Delta t^* \in [0, t^{(j+1)} - t^{(j)})\), the generalized Cauchy point \(\mathbf{x}^c = \mathbf{x}(t^{(j)} + \Delta t^*)\). If \(\Delta t^* < 0\), the generalized Cauchy point \(\mathbf{x}^c = \mathbf{x}(t^{(j)})\). In the above two scenarios, the set of active constraints include those components whose break point \(t\)-value is no larger than \(t^{(j)}\), that is,
\[
\mathcal{A}(\mathbf{x}^c) = \{i : t_i \leq t^{(j)}, i = 1, \ldots, n\}.
\]
If \(\Delta t^* \geq t^{(j+1)} - t^{(j)}\), the search of the generalized Cauchy point moves onto the next interval (line segment) of the projected gradient. The above processes stop until the critical point \(\Delta t^*\) is smaller than the interval length formed by two neighboring break points.

**Comment 3.** The calculation of the generalized Cauchy point is extremely efficient. Recall that \(\mathbf{B}^0 = \hat{\mathbf{K}}\) is a sparse matrix assembled from the standard finite element procedure, and \(\mathbf{e}_b\) is the \(b\) th unit vector. Assume the \(b\) th row of the sparse matrix \(\mathbf{B}^0\) has \(n_b\) non-zero entries, where \(n_b\) is way smaller than the degrees of freedom \(n\) of the discretized system. When only the multiplication is counted, the operation counts of \(\mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{d}^{(j-1)}\) and \(\mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{z}^{(j)}\) are both \(\mathcal{O}(n_b)\), and the operation count of \(\mathbf{e}_b^{\mathrm{T}}\mathbf{B}^0\mathbf{e}_b\) is \(\mathcal{O}(1)\). Also, since \(\mathbf{M} \in \mathbb{R}^{2m\times 2m}\), and \(\mathbf{w}_b, \mathbf{c}^{(j)}, \mathbf{p}^{(j-1)} \in \mathbb{R}^{2m}\), the operation counts of \(\mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{c}^{(j)}\), \(\mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{p}^{(j-1)}\), and \(\mathbf{w}^{\mathrm{T}}_b \mathbf{M} \mathbf{w}_b\) are all \(\mathcal{O}(4m^2)\), where \(m\) is the small number of vector-pairs stored for the compact representation of the L-BFGS matrix.

### 3.4. Subspace minimization

Once the generalized Cauchy point \(\mathbf{x}^c\) is located inside the interval \([\mathbf{x}(t^{(j)}), \mathbf{x}(t^{(j+1)})]\) on the piecewise linear path formed by the projected gradient, the active set of the box constraints is known, that is,
\[
\mathcal{A}(\mathbf{x}^c) = \{i : x^c_i = \mathrm{lb}_i\} \cup \{i : x^c_i = \mathrm{ub}_i\} = \{i : t_i \leq t^{(j)}, i = 1, \ldots, n\}.
\]

The components of the solution vector \(\mathbf{x}\) inside the active set are also known, each of which takes either the lower bound or the upper bound value,
\[
x_i = x^c_i = \begin{cases} \mathrm{ub}_i & \text{if } r_i < 0 \\ \mathrm{lb}_i & \text{if } r_i > 0 \end{cases}, \quad \forall i \in \mathcal{A}(\mathbf{x}^c).
\]

Assume that the size of the active set \(\mathcal{A}(\mathbf{x}^c)\) is \(q\), that is, there are totally \(q\) components of the generalized Cauchy point \(\mathbf{x}^c\) located at the boundary of the box constraints. These \(q\) components are fixed in the subsequent process. Assume that the \(i\) th component of the vector \(\mathbf{x}^c\) is fixed (at the boundary of the box constraints), let \(\mathbf{e}_i\) represent the corresponding unit vector. Then, a matrix \(\mathbf{Q} = \{\mathbf{e}_i\} \in \mathbb{R}^{n\times q}\) can be defined such that its columns are these unit vectors spanning the subspace of the fixed variables at \(\mathbf{x}^c\). Similarly, assume that the \(j\) th component of the vector \(\mathbf{x}^c\) is located inside the box constraints (the corresponding constraint is inactive), let \(\mathbf{e}_j\) represent the corresponding unit vector. Then, a matrix \(\mathbf{Z} \in \mathbb{R}^{n\times (n-q)}\) can be defined such that its columns are these unit vectors spanning the subspace of the free variables at \(\mathbf{x}^c\).

During the \(k\) th L-BFGS iteration, since at the generalized Cauchy point \(\mathbf{x}^c\), the variables contained in the active set of the box constraints \(\mathcal{A}(\mathbf{x}^c)\) are known, the goal is to solve the free variables through the minimization while keeping the components inside the active set of the box constraints fixed, that is,

\[
\mathbf{x}^* = \arg \min m_k(\mathbf{x}) \quad (25)
\]

subject to

\[
x_i = x^c_i, \quad \forall i \in \mathcal{A}(\mathbf{x}^c)
\]
and
\[
\mathrm{lb}_i \leq x_i \leq \mathrm{ub}_i, \quad \forall i \notin \mathcal{A}(\mathbf{x}^c). \quad (26)
\]

Recall that the objective function \(m_k(\mathbf{x})\) is the quadratic function defined in Eq. (15). There are generally two approaches to solve the above minimization, the primal approach that directly solves the free variables through the subspace minimization, or the dual approach that enforces the fixed variables via the Lagrange multipliers. Sections 3.4.1 and 3.4.2 present the primal approach via a direct matrix factorization and a conjugate gradient method, respectively. Section 3.4.3 lays out the Schur complement method for the dual approach and discusses the difficulties of using this method to solve the current problem.

#### 3.4.1. Direct matrix factorization for the primal approach

At the \(k\) th L-BFGS iteration, let \(\mathbf{x}^c \in \mathbb{R}^n\) represent the generalized Cauchy point, and \(\mathbf{Z}_k \in \mathbb{R}^{n\times (n-q)}\) represent the matrix whose columns span the subspace of the free variables at \(\mathbf{x}^c\). Assume that \(q\) components of the generalized Cauchy point \(\mathbf{x}^c\) are at the boundary of the box constraints. Then, the unknown vector \(\mathbf{x}\) can be expressed in the following form,

\[
\mathbf{x} = \mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}},
\]

where \(\hat{\mathbf{x}} \in \mathbb{R}^{n-q}\) is the vector of the free variables. The quadratic model in Eq. (15) can be transformed as

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

Since the L-BFGS matrix \(\mathbf{B}_k\) is positive definite, the above quadratic function \(\hat{m}_k(\hat{\mathbf{x}})\) of the free variables \(\hat{\mathbf{x}}\) reaches to its minimum at the critical point,

\[
\hat{m}'_k(\hat{\mathbf{x}}) = 0 \implies \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k \hat{\mathbf{x}} = -\mathbf{Z}^{\mathrm{T}}_k [\mathbf{r}_k + \mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k)].
\]

Let \(\hat{\mathbf{B}}_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k\) and \(\hat{\mathbf{r}}_k = \mathbf{Z}^{\mathrm{T}}_k [\mathbf{r}_k + \mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k)]\), the reduced linear system of the subspace minimization becomes

\[
\hat{\mathbf{B}}_k \hat{\mathbf{x}} = -\hat{\mathbf{r}}_k. \quad (27)
\]

Recall that \(\mathbf{B}_k\) is not directly available, since it is not stored in the memory component wise. Rather, \(\mathbf{B}_k\) is only available in the compact representation form as discussed in Section 3.2, that is,
\[
\mathbf{B}_k = \mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k.
\]
Therefore,
\[
\hat{\mathbf{B}}_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}_k \mathbf{Z}_k = \mathbf{Z}^{\mathrm{T}}_k (\mathbf{B}^0_k - \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k) \mathbf{Z}_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{Z}_k - \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k. \quad (28)
\]

Similar to the L-BFGS matrix \(\mathbf{B}_k\), the reduced matrix \(\hat{\mathbf{B}}_k\) is not stored in the memory component wise. As a result, the direct solvers based on the matrix factorization, such as the LU decomposition and the Cholesky decomposition, cannot be directly used to solve the reduced system in Eq. (27). However, \(\mathbf{B}^0_k\) is sparse (recall that it is assembled from the finite element procedure). Therefore, the reduced matrix \(\hat{\mathbf{B}}^0_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{Z}_k\) is also sparse and can be factorized using these techniques.

Using the Sherman–Morrison–Woodbury formula [45], as stated in the Appendix, the vector of the free variables can be solved as

\[
\hat{\mathbf{x}} = -\hat{\mathbf{B}}^{-1}_k \hat{\mathbf{r}}_k,
\]
\[
\hat{\mathbf{B}}^{-1}_k = (\hat{\mathbf{B}}^0_k)^{-1} + (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \left[\mathbf{I} - \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \right]^{-1} \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1}, \quad (29)
\]

where
\[
\mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k = \mathbf{Z}^{\mathrm{T}}_k [\mathbf{Y}_k \quad \mathbf{B}^0_k \mathbf{S}_k] = [\mathbf{Z}^{\mathrm{T}}_k \mathbf{y}_{k-m} \cdots \mathbf{Z}^{\mathrm{T}}_k \mathbf{y}_{k-1} \quad \mathbf{Z}^{\mathrm{T}}_k(\mathbf{B}^0_k \mathbf{s}_{k-m}) \cdots \mathbf{Z}^{\mathrm{T}}_k(\mathbf{B}^0_k \mathbf{s}_{k-1})].
\]

In order to calculate the inverse of \(\hat{\mathbf{B}}_k\), the inverse of \(\hat{\mathbf{B}}^0_k\) needs to be repeatedly applied to the column vectors in \(\mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k\). This can be achieved through the sparse LU decomposition of \(\hat{\mathbf{B}}^0_k\). Also, notice that in Eq. (29), the following term
\[
\mathbf{I} - \mathbf{M}_k \mathbf{W}^{\mathrm{T}}_k \mathbf{Z}_k (\hat{\mathbf{B}}^0_k)^{-1} \mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k \in \mathbb{R}^{2m\times 2m}
\]
is a \(2m\times 2m\) matrix, therefore, its inverse can be directly calculated with a negligible computational cost.

#### 3.4.2. Conjugate gradient method for the primal approach

As mentioned before, \(\mathbf{B}_k\) is not stored in the memory component wise and is only available in the compact representation form shown in Eq. (20). Therefore, the reduced matrix \(\hat{\mathbf{B}}_k\) is not directly available and only has the compact form shown in Eq. (28). Recall that \(\mathbf{B}_k\) is positive definite as long as the vector-pairs \(\{\mathbf{s}_i, \mathbf{y}_i\}_{i=k-m}^{k-1}\) satisfy the curvature condition Eq. (18). Therefore, \(\hat{\mathbf{B}}_k\) is also positive definite, making the conjugate gradient method an ideal candidate to solve the reduced linear system Eq. (27). Moreover, as a type of iterative solver, the conjugate gradient method does not require the component form of the reduced matrix \(\hat{\mathbf{B}}_k\). Rather, it only needs to know how to perform the matrix–vector multiplication between the matrix \(\hat{\mathbf{B}}_k\) and an arbitrary vector. Indeed, this matrix–vector multiplication is well defined according to Eq. (28), in which the matrix form of \(\hat{\mathbf{B}}^0_k = \mathbf{Z}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{Z}_k\), \(\mathbf{Z}^{\mathrm{T}}_k \mathbf{W}_k\), and \(\mathbf{M}_k\) are all available. Therefore, besides the direct matrix factorization approach provided in Section 3.4.1, the conjugate gradient method provides another option to solve the reduced linear system in order to obtain the free variables \(\hat{\mathbf{x}}\) at the generalized Cauchy point \(\mathbf{x}^c\).

It is well known that an effective preconditioner could significantly improve the performance of an iterative linear solver such as the conjugate gradient method. A good preconditioner should balance between the computational cost incurred for its own construction and the performance boost provided to the iterative solver. Let \(\mathbf{P}^{-1}\) represent the preconditioner. Considering the similarity between \(\hat{\mathbf{B}}^0_k\) and \(\hat{\mathbf{B}}_k\), in this work, the incomplete LU decomposition of \(\hat{\mathbf{B}}^0_k\) is chosen as the preconditioner for the conjugate gradient method to solve the reduced linear system, that is,

\[
\mathbf{P}^{-1} = \mathrm{ILU}(\hat{\mathbf{B}}^0_k) = \mathrm{ILU}(\mathbf{Z}^{\mathrm{T}}_k \mathbf{B}^0_k \mathbf{Z}_k).
\]

**Comment 4.** The reduced matrix \(\hat{\mathbf{B}}^0_k\) is positive definite since the curvature condition in Eq. (18) is satisfied due to the line search based on the strong Wolfe conditions. Therefore, the Cholesky decomposition is more efficient than the LU decomposition in the direct matrix factorization approach. Similarly, in the conjugate gradient method, the preconditioner based on the incomplete Cholesky decomposition is less expensive to construct compared with the counterpart based on the incomplete LU decomposition. However, since the UMFPACK used in the deal.II library is built on the direct sparse LU factorization, and the incomplete LU preconditioner is also readily available, these techniques are chosen here mainly for implementation convenience.

#### 3.4.3. Schur complement for the dual approach

Instead of formulating the minimization problem in Eq. (25) subject to the constraints in Eq. (26) as a subspace minimization, an alternative approach is to formulate the above problem via a dual approach using the Lagrange multipliers. Recall that \(\mathbf{x}_k\) is known at the beginning of the current iteration, and the generalized Cauchy point \(\mathbf{x}^c\) is already located on the piecewise linear path formed by the projected gradient. The unknown vector \(\mathbf{x}\) is expressed as

\[
\mathbf{x} = \mathbf{x}_k + \Delta \mathbf{x}_k.
\]

Plug the above equation into the quadratic model shown in Eq. (15),

\[
m_k(\mathbf{x}) = \Pi_k + \mathbf{r}_k^{\mathrm{T}} \Delta \mathbf{x}_k + \frac{1}{2} \Delta \mathbf{x}_k^{\mathrm{T}} \mathbf{B}_k \Delta \mathbf{x}_k.
\]

Also, at the generalized Cauchy point \(\mathbf{x}^c\), for the variables in the active set of the box constraints,
\[
x_i = x^c_i, \quad \forall i \in \mathcal{A}_k(\mathbf{x}^c),
\]
where \(x^c_i\) takes the upper bound value \(\mathrm{ub}_i\) or the lower bound value \(\mathrm{lb}_i\) of the box constraints. Recall that the matrix \(\mathbf{Q}_k = \{\mathbf{e}_i\} \in \mathbb{R}^{n\times q}\) is defined such that its columns are the unit vectors spanning the subspace of the fixed variables at \(\mathbf{x}^c\). Therefore,
\[
\mathbf{Q}^{\mathrm{T}}_k \Delta \mathbf{x}_k = \mathbf{Q}^{\mathrm{T}}_k (\mathbf{x} - \mathbf{x}_k) = \mathbf{Q}^{\mathrm{T}}_k (\mathbf{x}^c - \mathbf{x}_k).
\]

The minimization problem in Eq. (25) and the corresponding constraints in Eq. (26) can be rewritten as
\[
\min \left\{ \Pi_k + \mathbf{r}_k^{\mathrm{T}} \Delta \mathbf{x}_k + \frac{1}{2} \Delta \mathbf{x}_k^{\mathrm{T}} \mathbf{B}_k \Delta \mathbf{x}_k \right\}
\]
subject to
\[
\mathbf{Q}_k^{\mathrm{T}} \Delta \mathbf{x}_k = \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k)
\]
and the box constraints
\[
\mathbf{lb} - \mathbf{x}_k \leq \Delta \mathbf{x}_k \leq \mathbf{ub} - \mathbf{x}_k.
\]

For now, ignore the box constraints and let \(\lambda\) represent the Lagrange multiplier for the equality constraints. The optimality condition of the above constrained minimization is

\[
\begin{bmatrix} \mathbf{B}_k & \mathbf{Q}_k \\ \mathbf{Q}_k^{\mathrm{T}} & \mathbf{0} \end{bmatrix} \begin{bmatrix} \Delta \mathbf{x}_k \\ \lambda \end{bmatrix} = \begin{bmatrix} -\mathbf{r}_k \\ \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k) \end{bmatrix}. \quad (30)
\]

A popular method to solve the above linear system is to use the Schur complement approach by firstly solving \(\lambda\) from

\[
(\mathbf{Q}_k^{\mathrm{T}} \mathbf{B}_k^{-1} \mathbf{Q}_k) \lambda = -\mathbf{Q}_k^{\mathrm{T}} \mathbf{B}_k^{-1} \mathbf{r}_k - \mathbf{Q}_k^{\mathrm{T}} (\mathbf{x}^c - \mathbf{x}_k). \quad (31)
\]

Then, solve \(\Delta \mathbf{x}_k\) from

\[
\mathbf{B}_k \Delta \mathbf{x}_k = -(\mathbf{r}_k + \mathbf{Q}_k \lambda). \quad (32)
\]

Notice that \(\mathbf{B}_k^{-1}\) is contained in Eq. (31). Recall that \(\mathbf{B}_k\) only has the compact representation form. Therefore, similar to the primal approach, either the direct matrix factorization based on the Sherman–Morrison–Woodbury formula or the conjugate gradient method needs to be applied to obtain \(\mathbf{B}_k^{-1}\). To make the situation more complex, on top of the matrix factorization or the conjugate gradient method required for \(\mathbf{B}_k^{-1}\), \((\mathbf{Q}_k^{\mathrm{T}} \mathbf{B}_k^{-1} \mathbf{Q}_k)^{-1}\) has to be formed either via another direct linear solver or iterative linear solver to solve for the Lagrange multiplier \(\lambda\). Therefore, the Schur complement approach, which generally works well to solve a Karush–Kuhn–Tucker (KKT) system, is not practical for the current problem due to the lack of the component form of the L-BFGS matrix \(\mathbf{B}_k\).

Once the free variables \(\hat{\mathbf{x}}\) at the generalized Cauchy point \(\mathbf{x}^c\) are solved from the reduced linear system shown in Eq. (27), the search direction at the current \(k\) th L-BFGS iteration can be defined as

\[
\mathbf{p}_k = \mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}} - \mathbf{x}_k.
\]

The updated solution is obtained as

\[
\mathbf{x}_{k + 1} = \mathbf{x}_k + \alpha_k \mathbf{p}_k,
\]

where the step length parameter \(\alpha_k\) is determined by simultaneously satisfying the following two conditions:

1.  The updated solution \(\mathbf{x}_{k + 1}\) should stay feasible, that is,
    \[
    \mathbf{lb} \leq \mathbf{x}_k + \alpha_k \mathbf{p}_k \leq \mathbf{ub}.
    \]
2.  The positive step length \(\alpha_k\) should satisfy the strong Wolfe conditions such that the curvature condition in Eq. (18) is satisfied to ensure the positive-definiteness of the updated L-BFGS matrix \(\mathbf{B}_{k + 1}\).

Then, a new vector-pair \(\{\mathbf{s}_k, \mathbf{y}_k\}\) is computed. The oldest vector-pair \(\{\mathbf{s}_{k-m}, \mathbf{y}_{k-m}\}\) is removed from the vector-pair list \(\{\mathbf{s}_i, \mathbf{y}_i\}_{i=k-m+1}^{k-1}\) and the newly computed vector-pair \(\{\mathbf{s}_k, \mathbf{y}_k\}\) is inserted to form the new vector-pair list \(\{\mathbf{s}_i, \mathbf{y}_i\}_{i=k-m+1}^{k}\) for the next L-BFGS iteration. Algorithm 1 summaries the major steps of the gradient projection L-BFGS-B scheme.

---

**Algorithm 1** Gradient projection L-BFGS-B monolithic scheme

1:  Set \(\mathbf{x}_0 = [\mathbf{u}_A, d_A]^{(n)}\). Initial guess uses the solution at \(t_n\).
2:  For the phase-field, set the lower bound \(\mathbf{lb} = [d_A]^{(n)}\) and the upper bound \(\mathbf{ub} = [1.0]\).
3:  **for** \(k = 0 : N_{\mathrm{max}}\) **do** ▷ L-BFGS-B iterations
4:    Calculate the gradient (residual) \(\mathbf{r}_k = [\mathbf{r}_{u_k}, \mathbf{r}_{d_k}]\) according to Eq. (9).
5:    **if** \(\mathcal{A}_k(\mathbf{x}^c) = \mathcal{A}_{k-1}(\mathbf{x}^c)\) and \(\| \Delta \mathbf{x}_k\| _2 < \mathrm{tol}\) and \(\| \mathrm{Proj}_c(\mathbf{x}_k - \mathbf{r}_k,\mathbf{lb},\mathbf{ub}) - \mathbf{x}_k\| _2 < \mathrm{tol}\) **then**
6:      break.
7:    **end if**
8:    Calculate \(\mathbf{B}_k^0 = \hat{\mathbf{K}}\) according to Eq. (11). ▷ \(\hat{\mathbf{K}}\) could be formed every 3 to 5 iterations instead
9:    Calculate \(\mathbf{W}_k\) and \(\mathbf{M}_k\) for the compact representation of \(\mathbf{B}_k = \mathbf{B}_k^0 - \mathbf{W}_k\mathbf{M}_k\mathbf{W}_k^{\mathrm{T}}\) according to Eq. (20).
10:   Compute the generalized Cauchy point \(\mathbf{x}^c\) and the active set \(\mathcal{A}_k(\mathbf{x}^c)\) according to Sect. 3.3.
11:   Form the subspace matrices \(\mathbf{Z}_k\) for the free variables and \(\mathbf{Q}_k\) for the fixed variables at \(\mathbf{x}^c\).
12:   Form the reduced matrix \(\hat{\mathbf{B}}_k = \mathbf{Z}_k^{\mathrm{T}}\mathbf{B}_k\mathbf{Z}_k\) and the reduced vector \(\hat{\mathbf{r}}_k = \mathbf{Z}_k^{\mathrm{T}} [\mathbf{r}_k + \mathbf{B}_k(\mathbf{x}^c - \mathbf{x}_k)]\).
13:   Solve the reduced linear system \(\hat{\mathbf{B}}_k \hat{\mathbf{x}} = -\hat{\mathbf{r}}_k\) for the free variables \(\hat{\mathbf{x}}\) at the generalized Cauchy point \(\mathbf{x}^c\).
14:   Compute the search direction \(\mathbf{p}_k = \mathbf{x}^c + \mathbf{Z}_k \hat{\mathbf{x}} - \mathbf{x}_k\).
15:   Calculate the step length \(\alpha_k\) according to the strong Wolfe conditions and ensure \(\mathbf{x}_k + \alpha_k \mathbf{p}_k\) is still feasible.
16:   Set \(\mathbf{x}_{k + 1} = \mathbf{x}_k + \alpha_k \mathbf{p}_k\) and compute \(\mathbf{r}_{k + 1}\).
17:   Compute \(\mathbf{s}_k = \mathbf{x}_{k + 1} - \mathbf{x}_k\) and \(\mathbf{y}_k = \mathbf{r}_{k + 1} - \mathbf{r}_k\).
18:   **if** \(\mathbf{s}_k^{\mathrm{T}}\mathbf{y}_k > 0\) **then**
19:     **if** \(k > m\) **then** ▷ \(m\) is the prescribed number for “limited-memory”
20:       Remove the oldest vector-pair \(\{\mathbf{s}_{k-m}, \mathbf{y}_{k-m}\}\) from the vector-pair list \(\{\mathbf{s}_i, \mathbf{y}_i\}_{i=k-m}^{k-1}\).
21:     **end if**
22:     Insert \(\{\mathbf{s}_k, \mathbf{y}_k\}\) into the vector-pair list \(\{\mathbf{s}_i, \mathbf{y}_i\}\).
23:   **end if**
24: **end for**

---

## 4. Numerical examples

In this section, several 2D and 3D numerical examples are provided to demonstrate the capabilities of the proposed gradient projection L-BFGS-B phase-field monolithic scheme. To emphasize the scheme's effectiveness in enforcing the phase-field irreversibility condition, displacement-controlled cyclic loading is applied. In order to alleviate the heavy computational cost, the proposed L-BFGS-B monolithic scheme is combined with a predictor-corrector adaptive mesh refinement technique originally proposed by Heister et al. [4]. During the finite element formulation, linear constraints are introduced by the applied Dirichlet boundary conditions and the hanging-node constraints during the adaptive mesh refinement. The algorithmic details regarding how to handle these linear constraints in the L-BFGS method can be found in [26]. The L-BFGS-B monolithic scheme is implemented in deal.II [54], which is an open-source C++ finite element library.

### 4.1. Cyclic tension-compression test

The 2D example shown in Fig. 3(a) is widely adopted in the literature of the phase-field fracture simulation to test the robustness of the numerical method. Under the displacement-controlled load shown in Fig. 3(b) that is vertically applied at the top edge of the unit square, the preexisting crack is expected to propagate from the domain center to the right edge in a unstable fashion. Since the crack propagation is known a priori, a preferred mesh shown in Fig. 3(a) is sufficient and no adaptive mesh refinement is necessary. The material parameters include the Lamé parameters \(\lambda = 121.15 \mathrm{kN/mm^2}\) and \(\mu = 80.77 \mathrm{kN/mm^2}\), the critical energy release rate \(g_c = 2.7 \times 10^{-3} \mathrm{kN/mm}\), and the small positive parameter \(k = 0.0\) in the degradation function \(g(d)\). The phase-field length-scale parameter is chosen as \(l = 0.0075 \mathrm{mm}\), and the mesh size in the preferred region is \(h = 0.00265 \mathrm{mm}\) such that the ratio is \(h/l \approx 1/3\).

Fig. 4 shows the phase-field distribution in the unit square under the vertical cyclic load at various pseudo time steps. When the vertical displacement increases from \(u_{y} = 5.9 \times 10^{-3} \mathrm{mm}\) to \(u_{y} = 6.0 \times 10^{-3} \mathrm{mm}\), the crack abruptly propagates from the preexisting crack tip at the domain center to the right edge of the domain. Moreover, after the crack is fully developed, the phase-field distribution remains unchanged in the subsequent load step, and the phase-field nodal value never decreases during the unloading and compression stages due to the enforced irreversibility condition. Also, the phase-field value is strictly in the range of \([0, 1]\). Fig. 5(a) shows the relationship between the vertical reaction force at the bottom edge and the pseudo time. Fig. 5(b) shows the relationship between the vertical reaction force at the bottom edge and the vertical displacement at the top edge. The sample loses the tensile load bearing capacity after the crack is fully developed at \(t = 6.0 \times 10^{-3}\). Subsequently, the crack never self-heals during the unloading stage due to the enforced phase-field (damage) irreversibility. The sample starts to be compressed after \(t = 14.0 \times 10^{-3}\). Due to the tension-compression asymmetry, the sample still retains the fully elastic behavior under compression. Fig. 5(c) compares the crack dissipation energy calculated as

\[
E_{\mathrm{crack}} = g_{\mathrm{c}}\Gamma_{l}(d) = g_{\mathrm{c}}\int_{\Omega}\gamma (d,\nabla d)\mathrm{d}\Omega = g_{\mathrm{c}}\int_{\Omega}\frac{1}{2l}\left(d^{2} + l^{2}\nabla d\cdot \nabla d\right)\mathrm{d}\Omega \quad (33)
\]

with and without enforcing the phase-field irreversibility condition. When the irreversibility condition is properly enforced by the proposed L-BFGS-B scheme, the crack dissipation energy remains constant after the crack is fully developed. When the irreversibility condition is not enforced, the crack dissipation energy decreases to zero during the unloading stage, which obviously violates the thermodynamic consistency.

Fig. 6 shows the total energy of the system, the \(l_{2}\)-norm of the projected gradient of the displacement field \(\| \mathrm{Proj}_{c}^{\pmb{u}}\|_{2}\) and the phase-field \(\| \mathrm{Proj}_{c}^{d}\|_{2}\) according to Eq. (19), and the number of the active constraints during the L-BFGS-B iterations at the critical load step when the crack unstably propagates. Due to the line search method based on the strong Wolfe conditions, the total energy of the system consistently decreases during the iterations and eventually reaches to a plateau. The convergence criteria are satisfied when the \(l_{2}\)-norms of the projected gradient and the solution increments are reduced to the prescribed tolerance \((10^{-6})\) and the active constraints remain unchanged during two consecutive iterations.

### 4.2. Cyclic shear test

The domain of the unit square shown in Fig. 7(a) has a preexisting crack and undergoes a displacement-controlled cyclic load shown in Fig. 7(b) applied horizontally at the top edge. Since in this case the crack propagation path is not known a priori, the mesh is only preferred around the preexisting crack tip and the adaptive mesh refinement is adopted during the crack propagation.

The material parameters include the Lamé parameters \(\lambda = 121.15 \mathrm{kN/mm^2}\) and \(\mu = 80.77 \mathrm{kN/mm^2}\), the critical energy release rate \(g_c = 2.7 \times 10^{-3} \mathrm{kN/mm}\), and the small positive parameter \(k = 0.0\) in the degradation function \(g(d)\). The phase-field length-scale parameter is chosen as \(l = 0.0075 \mathrm{mm}\). The mesh is adaptively refined during the crack propagation such that the maximum allowed ratio between the mesh size and the phase-field length-scale is \((h/l)_{\mathrm{max}} = 0.5\), which guarantees that there are at least three elements within the phase-field length-scale.

Fig. 8 reports the phase-field distributions at various load steps under the horizontally applied cyclic loading. Due to the effectively enforced phase-field irreversibility condition, the phase-field value is strictly in the range of \([0, 1]\). During the loading stage, the crack fully propagates downward in the lower-right region when the horizontal displacement at the top edge is \(u_{x} = 0.015 \mathrm{mm}\), as shown in Fig. 8(c). Afterwards, the displacement-controlled load is reversed in the opposite direction until another crack fully propagates upward in the upper-right region, as shown in Fig. 8(f). The symmetry of the two cracks formed during the applied cyclic loading shown in Fig. 7(b) indicates the effectiveness of the gradient projection method in enforcing the phase-field irreversibility. Furthermore, Fig. 9 shows the constraint status of the phase-field degree of freedom (DoF) at each finite element node. Let \(d\) represent the phase-field value at a node and \(d^{(n)}\) represent the phase-field value at the same node from the previous time step. There are four different statuses in total. The constraint at the phase-field DoF is inactive \((d^{(n)} < d < 1.0)\), active at the lower bound \((d = d^{(n)})\), active at the upper bound \((d = 1.0)\), or active when the lower and upper bounds overlap \((d^{(n)} = d = 1.0)\) indicating that the phase-field value at the previous time step is already 1.0. The lower bounds of the box constraints prevent the phase-field from decreasing, while the upper bounds of the box constraints prevent the phase-field from going over 1.0.

Fig. 10(a) reports the relationship between the horizontal reaction force at the domain bottom edge and the pseudo time. Fig. 10(b) reports the relationship between the horizontal reaction force at the domain bottom edge and the horizontal displacement at the domain top edge. The sample experiences a full cycle of loading, unloading, reverse loading, and reverse unloading. The effective enforcement of the phase-field irreversibility condition prevents the crack from self-healing during the unloading stages. Moreover, the force-displacement relationship associated with the downward crack (in the lower-right corner) developed during the loading stage is the same as the one associated with the upward crack (in the upper-right corner) developed during the reverse loading stage. Fig. 10(c) exhibits the crack dissipation energy, calculated according to Eq. (33), with and without the enforcement of the phase-field irreversibility condition. It is obvious that the crack dissipation energy should only monotonically increase during the loading stage and remain constant during the unloading stage. Also, the crack energy dissipation occurred during the loading stage is the same as the one occurred during the reverse loading stage. These observations demonstrate the effectiveness of the proposed L-BFGS-B method in enforcing the phase-field irreversibility condition.

Fig. 11 shows the total energy of the system, the \(l_{2}\)-norm of the projected gradient of the displacement field \(\| \mathrm{Proj}_{c}^{\pmb{u}}\|_{2}\) and the phase-field \(\| \mathrm{Proj}_{c}^{d}\|_{2}\) according to Eq. (19), and the numbers of active constraints during the L-BFGS-B iterations inside a typical time step of the cyclic shear test. Due to the line search method based on the strong Wolfe conditions, the total energy of the system consistently decreases during the iterations and eventually reaches to a plateau. Recall that the active lower/upper bound means that the lower bound and the upper bound at a phase-field degree of freedom overlap. This case happens when the phase-field value from the previous time step already reaches to 1.0. The L-BFGS-B iteration converges when the \(l_{2}\)-norms of the projected gradient and the solution increments are reduced to the prescribed tolerance \((10^{-6})\) and the active lower bounds, the active upper bounds, and the active lower/upper bounds remain constant during two consecutive iterations.

Tables 1 and 2 report the convergence history of the L-BFGS-B iterations in a typical time step during the loading and unloading stages in the cyclic shear test, respectively. During the loading stage, as shown in Table 1, only part of the phase-field DoFs are at the boundary of the box constraints. The rest of the phase-field DoFs still need to be solved during the subspace minimization. On the other hand, during the unloading stage, as shown in Table 2, since the phase-field value never decreases due to the enforced irreversibility condition, all of the phase-field DoFs are at the lower bound of the box constraints at the generalized Cauchy point. Therefore, the \(l_{2}\)-norm of the projected gradient of the phase-field \(\| \mathrm{Proj}_{c}^{d}\|_{2}\), calculated according to Eq. (19), becomes zero, and only the displacement field needs to be calculated via the subspace minimization.

**Table 1**
Convergence history of the L-BFGS-B iterations in a typical time step during the loading stage in the cyclic shear test.

| Itr. | L.B.\(^a\) | U.B.\(^b\) | L.U.B.\(^c\) | CG-itr.\(^d\) | Energy | \(\|\mathrm{Proj}_c^{\pmb{u}}\|_2\) | \(\|\mathrm{Proj}_c^{d}\|_2\) | \(\|\Delta \pmb{u}\|_2\) | \(\|\Delta d\|_2\) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 2885 | 0 | 347 | 115 | 5.388e-03 | 3.900e-01 | 1.712e-04 | 2.800e-02 | 4.798e-01 |
| 2 | 1759 | 10 | 347 | 87 | 5.351e-03 | 4.831e-03 | 1.677e-04 | 8.286e-04 | 5.110e-01 |
| ... | | | | | | | | | |
| 73 | 1934 | 28 | 347 | 7 | 5.345e-03 | 1.694e-06 | 4.608e-09 | 3.233e-06 | 4.043e-04 |
| 74 | 1934 | 28 | 347 | 7 | 5.345e-03 | 1.237e-06 | 3.050e-09 | 3.051e-06 | 3.716e-04 |
| 75 | 1934 | 28 | 347 | 0 | 5.345e-03 | 9.287e-07 | 2.978e-09 | 4.297e-09 | 1.378e-11 |

\(^a\) Number of active constraints at lower bound.
\(^b\) Number of active constraints at upper bound.
\(^c\) Number of active constraints at lower/upper bound.
\(^d\) Number of conjugate-gradient iterations needed for the subspace minimization.

**Table 2**
Convergence history of the L-BFGS-B iterations in a typical time step during the unloading stage in the cyclic shear test.

| Itr. | L.B. | U.B. | L.U.B. | CG-itr. | Energy | \(\|\mathrm{Proj}_c^{\pmb{u}}\|_2\) | \(\|\mathrm{Proj}_c^{d}\|_2\) | \(\|\Delta \pmb{u}\|_2\) | \(\|\Delta d\|_2\) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 4813 | 0 | 444 | 122 | 4.292e-03 | 3.902e-01 | 7.317e-05 | 2.981e-02 | 2.226e-01 |
| 2 | 4863 | 0 | 444 | 79 | 4.284e-03 | 1.548e-03 | 8.148e-05 | 3.116e-04 | 2.226e-01 |
| ... | | | | | | | | | |
| 10 | 4863 | 0 | 444 | 6 | 4.284e-03 | 1.287e-05 | 0.000e+00 | 3.725e-07 | 0.000e+00 |
| 11 | 4863 | 0 | 444 | 7 | 4.284e-03 | 2.158e-06 | 0.000e+00 | 4.120e-07 | 0.000e+00 |
| 12 | 4863 | 0 | 444 | 1 | 4.284e-03 | 7.824e-07 | 0.000e+00 | 4.220e-08 | 0.000e+00 |

### 4.3. Cyclic L-shape bending test

In this example, an L-shape sample with the geometrical dimensions shown in Fig. 12(a) is fixed at the bottom edge. A displacement-controlled cyclic load, as shown in Fig. 12(b), is applied upward at a location \(30 \mathrm{mm}\) to the sample's right edge. Again, since the crack propagation path is not known a priori, the adaptive mesh refinement is adopted during the crack propagation. The material parameters include the Lamé parameters \(\lambda = 6.16 \mathrm{kN/mm^2}\) and \(\mu = 10.95 \mathrm{kN/mm^2}\), the critical energy release rate \(g_{c} = 9.5 \times 10^{-5} \mathrm{kN/mm}\), and the small positive parameter \(k = 0.0\) in the degradation function \(g(d)\). The phase-field length-scale parameter is chosen as \(l = 1.0 \mathrm{mm}\). The mesh is adaptively refined during the crack propagation such that the maximum allowed ratio between the mesh size and the phase-field length-scale is \((h/l)_{\mathrm{max}} = 0.5\), which guarantees that there are at least three elements within the phase-field length-scale.

Fig. 13 shows the phase-field distribution and the corresponding adaptively refined mesh in the L-shape sample at various time steps under the displacement-controlled cyclic loading. Under the load applied upward, the crack initiates at the corner of the L-shape sample and propagates along a curved path. Also, due to the enforced phase-field irreversibility condition, the phase-field value never goes over the upper bound 1.0. During the unloading stage of the cyclic loading, the phase-field value never decreases, indicating that there is no self-healing.

Fig. 14 further reports the constraint status of the phase-field DoFs in two consecutive loading steps of the L-shape test. Fig. 14(a) shows the constraint status when the crack is firstly initiated at the L-shape corner. In this step, the upper bound of the phase-field DoFs inside the crack is activated such that the converged phase-field solution does not go over the upper limit 1.0. In the subsequent loading step, as shown in Fig. 14(b), the constraint status of the phase-field DoFs in the already developed crack changes from the “active (U.B.)” to the “active (L.U.B.)”, indicating the lower bound and the upper bound overlap and are both equal to 1.0. Meanwhile, the crack keeps propagating forward, and more phase-field DoFs reach to the upper limit value 1.0.

Fig. 15(a) reports the relationship between the vertical reaction force at the sample bottom edge and the pseudo time. Fig. 15(b) reports the relationship between the same reaction force and the vertically applied displacement-controlled cyclic load. During each of the three unloading stages, the force-displacement curve follows a unloading path with a reduced tangent due to the accumulated damage inside the L-shape sample. As the crack propagates and the damage accumulates, the tangent of the unloading path becomes smaller, indicating that the enforced phase-field irreversibility condition successfully prevents the crack from self-healing. Moreover, due to the tension-compression asymmetry in the adopted phase-field formulation, once the vertically applied displacement load reduces to zero and starts to point toward the reversed direction, the already developed crack experiences compression internally. As a result, the L-shape sample still exhibits fully elastic responses under the reversed loading. Fig. 15(c) shows the crack dissipation energy, calculated according to Eq. (33), with and without enforcing the phase-field irreversibility condition. Due to the effectiveness of the L-BFGS-B scheme, the crack dissipation energy remains constant during the unloading stages and only starts to increase when the crack further propagates inside the L-shape sample.

Fig. 16 shows the total energy of the system, the \(l_{2}\)-norm of the projected gradient of the displacement field \(\| \mathrm{Proj}_{c}^{\pmb{u}}\|_{2}\) and the phase-field \(\| \mathrm{Proj}_{c}^{d}\|_{2}\) according to Eq. (19), and the numbers of active constraints during the L-BFGS-B iterations inside a typical time step of the L-shape sample under the displacement-controlled cyclic load. Due to the line search method based on the strong Wolfe conditions, the total energy of the system consistently decreases during the iterations and eventually reaches to a plateau. Recall that the active lower/upper bound means that the lower bound and the upper bound at a phase-field degree of freedom overlap. This case happens when the phase-field value from the previous time step already reaches to 1.0, which is the lower bound for the current time step. The L-BFGS-B iteration converges when the \(l_{2}\)-norms of the projected gradient and the solution increments are reduced to the prescribed tolerance \((10^{-6})\) and the active lower bounds, the active upper bounds, and the active lower/upper bounds remain constant during two consecutive iterations, respectively.

Tables 3 and 4 report the convergence history of the L-BFGS-B iterations in a typical time step during the loading and unloading stages in the cyclic L-shape bending test, respectively. Similar to the previous example, during the loading stage, as shown in Table 3, only part of the phase-field DoFs are at the boundary of the box constraints detected by the gradient projection method. The rest of the phase-field DoFs still need to be solved during the subspace minimization. On the other hand, during the unloading stage, as shown in Table 4, since the phase-field value never decreases due to the enforced irreversibility condition, all of the phase-field DoFs are at the lower bound of the box constraints at the generalized Cauchy point. Therefore, the \(l_{2}\)-norm of the projected gradient of the phase-field \(\| \mathrm{Proj}_{c}^{d}\|_{2}\), calculated according to Eq. (19), becomes zero, and only the displacement field needs to be calculated via the subspace minimization.

**Table 3**
Convergence history of the L-BFGS-B iterations in a time step during the loading stage in the cyclic L-shape bending test.

| Itr. | L.B.\(^a\) | U.B.\(^b\) | L.U.B.\(^c\) | CG-itr.\(^d\) | Energy | \(\|\mathrm{Proj}_c^{\pmb{u}}\|_2\) | \(\|\mathrm{Proj}_c^{d}\|_2\) | \(\|\Delta \pmb{u}\|_2\) | \(\|\Delta d\|_2\) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 3743 | 0 | 200 | 137 | 2.858e-02 | 3.690e-01 | 1.686e-04 | 3.848e-01 | 3.185e-02 |
| 2 | 2054 | 5 | 200 | 102 | 2.857e-02 | 3.565e-04 | 2.863e-04 | 6.071e-04 | 2.319e-01 |
| ... | | | | | | | | | |
| 73 | 2177 | 15 | 200 | 5 | 2.853e-02 | 8.095e-07 | 6.515e-09 | 3.704e-06 | 5.522e-05 |
| 74 | 2177 | 15 | 200 | 1 | 2.853e-02 | 1.481e-06 | 8.394e-09 | 8.375e-07 | 1.727e-05 |
| 75 | 2177 | 15 | 200 | 0 | 2.853e-02 | 6.913e-07 | 7.250e-09 | 2.705e-08 | 2.837e-10 |

\(^a\) Number of active constraints at lower bound.
\(^b\) Number of active constraints at upper bound.
\(^c\) Number of active constraints at lower/upper bound.
\(^d\) Number of conjugate-gradient iterations needed for the subspace minimization.

**Table 4**
Convergence history of the L-BFGS-B iterations in a time step during the unloading stage in the cyclic L-shape bending test.

| Itr. | L.B. | U.B. | L.U.B. | CG-itr. | Energy | \(\|\mathrm{Proj}_c^{\pmb{u}}\|_2\) | \(\|\mathrm{Proj}_c^{d}\|_2\) | \(\|\Delta \pmb{u}\|_2\) | \(\|\Delta d\|_2\) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 8190 | 0 | 588 | 223 | 2.821e-02 | 3.691e-01 | 9.896e-05 | 4.502e-01 | 1.888e-02 |
| 2 | 8192 | 0 | 588 | 75 | 2.820e-02 | 9.749e-06 | 1.952e-04 | 2.367e-04 | 1.888e-02 |
| ... | | | | | | | | | |
| 88 | 8192 | 0 | 588 | 4 | 2.820e-02 | 1.544e-05 | 0.000e+00 | 1.843e-06 | 0.000e+00 |
| 98 | 8192 | 0 | 588 | 7 | 2.820e-02 | 6.380e-06 | 0.000e+00 | 5.370e-06 | 0.000e+00 |
| 108 | 8192 | 0 | 588 | 0 | 2.820e-02 | 9.044e-07 | 0.000e+00 | 3.092e-08 | 0.000e+00 |

### 4.4. Three-dimensional torsion test

In the last example, a 3D beam shown in Fig. 17 undergoes a torsional load. The beam has a square cross section of \(50 \times 50 \mathrm{mm}\) and a length of \(200 \mathrm{mm}\). In the middle of the beam, there is a preexisting \(45^{\circ}\) notch. The width of the notch opening is \(5 \mathrm{mm}\), and the depth of the notch is \(25 \mathrm{mm}\). The right end of the beam is fixed in all three directions. A displacement controlled rotation is applied at the left end of the beam so that it undergoes torsion. At an arbitrary node \((0,y,z)\) located at the left surface of the beam, the Dirichlet boundary conditions are applied as

\[
u_{y} = z \tan t, \quad u_{z} = -y \tan t,
\]

where \(t\) is the pseudo time representing the load step. The initial mesh is pre-refined near the bottom region of the notch where the crack is expected to initialize. The material parameters include the Lamé parameters \(\lambda = 9.72 \mathrm{kN/mm^2}\) and \(\mu = 14.58 \mathrm{kN/mm^2}\), the critical energy release rate \(g_{c} = 1.1 \times 10^{-4} \mathrm{kN/mm}\), and the small positive parameter \(k = 0.0\) in the degradation function \(g(d)\). The phase-field length-scale parameter is chosen as \(l = 1.0 \mathrm{mm}\).

Fig. 18 shows different views of the crack surface (the phase-field value \(d \geq 0.8\)) and the corresponding adaptively refined mesh when the rotation angle (pseudo time) \(t = 8.0 \times 10^{-3}\). Obviously, the crack surface is non-planar under the torsional load. With the adaptive mesh refinement technique, the total number of DoFs is still more than one million. From this example, we can see that adaptive mesh refinement is absolutely necessary for any practical 3D phase-field fracture simulations to control the computational cost. Furthermore, Fig. 19 shows the constraint status of the phase-field degree of freedom (DoF) at each finite element node. The lower bounds of the box constraints prevent the phase-field from decreasing, while the upper bounds of the box constraints prevent the phase-field from going over 1.0.

Fig. 20 shows the total energy of the system, the \(l_{2}\)-norm of the projected gradient of the displacement field \(\| \mathrm{Proj}_{c}^{\pmb{u}}\|_{2}\) and the phase-field \(\| \mathrm{Proj}_{c}^{d}\|_{2}\) according to Eq. (19), and the numbers of active constraints during the L-BFGS-B iterations inside a typical time step of the 3D beam torsion test. Similar to the previous 2D cases, due to the line search method based on the strong Wolfe conditions, the total energy of the system consistently decreases during the iterations and eventually reaches to a plateau. Recall that the active lower/upper bound means that the lower bound and the upper bound at a phase-field degree of freedom overlap, that is, \(d^{(n)} = d = 1.0\). This case happens when the phase-field value at the previous time step already reaches to 1.0, which forms the lower bound for the current time step. The L-BFGS-B iteration converges when the \(l_{2}\)-norms of the projected gradient and the solution increments are reduced to the prescribed tolerance \((10^{-6})\), while at the same time the active lower bounds, the active upper bounds, and the active lower/upper bounds remain constant during two consecutive iterations, respectively. Table 5 reports the convergence history of the L-BFGS-B iterations in a typical time step of the 3D torsion test. When the iterations converge, the numbers of active lower bounds, active upper bounds, and active lower/upper bounds remain constant.

**Table 5**
Convergence history of the L-BFGS-B iterations in a time step during the loading stage in the 3D torsion test.

| Itr. | L.B.\(^a\) | U.B.\(^b\) | L.U.B.\(^c\) | CG-itr.\(^d\) | Energy | \(\|\mathrm{Proj}_c^{\pmb{u}}\|_2\) | \(\|\mathrm{Proj}_c^{d}\|_2\) | \(\|\Delta \pmb{u}\|_2\) | \(\|\Delta d\|_2\) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| 1 | 5920 | 40 | 5586 | 23 | 5.298e-01 | 5.266e+00 | 1.860e-02 | 8.872e-01 | 2.217e+00 |
| 2 | 1864 | 221 | 3558 | 62 | 5.154e-01 | 1.642e-01 | 1.605e-02 | 2.877e-02 | 2.376e+00 |
| ... | | | | | | | | | |
| 350 | 2643 | 166 | 2558 | 6 | 5.132e-01 | 7.876e-07 | 7.380e-09 | 1.381e-06 | 4.274e-05 |
| 351 | 2643 | 166 | 2558 | 6 | 5.132e-01 | 9.710e-07 | 5.189e-09 | 8.288e-07 | 2.266e-05 |
| 352 | 2643 | 166 | 2558 | 0 | 5.132e-01 | 6.891e-07 | 5.981e-09 | 2.877e-12 | 2.497e-14 |

\(^a\) Number of active constraints at lower bound.
\(^b\) Number of active constraints at upper bound.
\(^c\) Number of active constraints at lower/upper bound.
\(^d\) Number of conjugate-gradient iterations needed for the subspace minimization.

## 5. Wall-clock time and convergence behavior

This section aims to accomplish the following two tasks. First, the wall-clock time of each numerical example provided in Section 4 is reported. Additionally, the percentages (with respect to the total time) of the important steps in the proposed L-BFGS-B solver, including the Cauchy point calculation, the matrix formation and update, the quadrature data update, and the subspace minimization, are presented in detail. Second, the convergence behavior of the proposed L-BFGS-B solver is compared with the counterparts of two existing phase-field solvers in the literature, a staggered scheme based on the alternate minimization [15,16] and a limited-memory BFGS monolithic solver [26] using the history variable of the maximum positive strain energy \(H\) [2] for the enforcement of the phase-field irreversibility. In order to ensure a fair comparison, all the three solvers are implemented under the finite element library deal.II with the same programming techniques and style. For each example problem, the simulations using the three solvers are based on the same pre-refined mesh without mesh adaptivity, adopt the same type of linear solver and the same tolerances of the residual-based convergence criteria. In the proposed method, even though the diagonal block matrix \(\hat{\mathbf{K}}\) in Eq. (11) only needs to be updated every few iterations, it is updated every iteration here to ensure a fair comparison among the three solvers.

### 5.1. Wall-clock time

All the numerical simulations provided in Section 4 are performed on a workstation laptop with 16 CPUs (8 cores, 2 threads per core) of 11th Gen Intel(R) Core(TM) i7-11800H @ 2.30 GHz. The only adopted parallelism in each simulation is based on the Threading Building Blocks (TBB) library. The TBB (multi-threading) feature is only invoked during the formation of the element stiffness matrix and the right-hand side vector as well as the global assembly stage when looping over all the elements. The rest part of the code is purely sequential. As discussed in Section 3.4, once the set of active constraints are detected via the calculation of the generalized Cauchy point, there are three different approaches to solve for the remaining free variables, including the direct matrix factorization, the conjugate gradient method with an appropriate preconditioner, and the Schur complement approach using the Lagrange multipliers. As discussed in Section 3.4.3, since the L-BFGS matrix \(\mathbf{B}_k\) only has the compact representation form, it is quite cumbersome to implement the Schur complement approach. Therefore, only the direct matrix factorization approach and the conjugate gradient approach are implemented and used for the provided numerical examples. The direct matrix factorization is based on the sparse LU decomposition, and the conjugate gradient method uses the incomplete LU decomposition as the preconditioner. The conjugate gradient method is consistently more efficient than the direct matrix factorization in all the provided numerical examples.

Table 6 reports the wall-clock time of the numerical examples provided in Section 4. Particularly, the percentage with respect to the total time spent on each important step in the proposed numerical method is also listed, including the Cauchy point calculation, the matrix formation and update, the quadrature data update, and the subspace minimization. For all the provided numerical examples, the Cauchy point calculations only take a tiny percentage \((< 1\%)\) of the total time. The majority of the wall-clock time is spent on solving the subspace minimization problem and updating the data at the quadrature points.

**Table 6**
Wall-clock time of the provided numerical examples.

| Example | Growth of No. DoFs | Wall-clock time (s) | Cauchy point calculation | Matrix formation and update | Quadrature data update | Subspace minimization\(^c\) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| Tension\(^a\) | 16401–16401 | \(3.121 \times 10^2\) | 0.20% | 6.5% | 3.64% | 84.1% |
| Shear\(^b\) | 3792–15921 | \(7.783 \times 10^2\) | 0.36% | 21.3% | 17.8% | 43.3% |
| L-shape bending\(^b\) | 6258–32973 | \(7.378 \times 10^2\) | 0.25% | 17.1% | 14.4% | 51.8% |
| Torsion (3D)\(^b\) | 155892–137536 | \(6.541 \times 10^4\) | 0.14% | 18.0% | 31.8% | 38.5% |

\(^a\) Pre-refined mesh without mesh adaptivity.
\(^b\) Adaptive mesh refinement.
\(^c\) Conjugate gradient method using the incomplete LU (ILU) decomposition as the preconditioner.

### 5.2. Comparison of convergence behaviors

The convergence behavior of the proposed L-BFGS-B monolithic scheme is compared with the counterparts of a staggered scheme based on the alternate minimization [15,16] and a limited-memory BFGS monolithic scheme [26]. Both the staggered scheme and the BFGS scheme adopt the history variable approach [2] to enforce the phase-field irreversibility, while recall that the proposed L-BFGS-B scheme uses the gradient projection method to treat the phase-field box constraints for the same purpose. In order to ensure a fair comparison, all the three phase-field solvers are implemented under the same finite element library deal.II. During the implementations, the same programming techniques and style as well as the same type of linear solver (the conjugate-gradient method) are adopted. In order to eliminate the potential impact of mesh adaptivity on the solver behavior, all the simulations are performed using the same set of pre-refined meshes as shown in Fig. 21. The boundary conditions of the tensile, shear, and bending tests are exactly the same as the counterparts used in the numerical examples provided in Sections 4.1, 4.2, and 4.3, respectively.

All the three numerical tests go through a complete loading-unloading cycle, and the required number of iterations is reported for each load step. For the tensile test, during the loading stage, the vertical displacement \(u_{y}\) first increases to \(5.0 \times 10^{-3} \mathrm{mm}\) with the step size \(1.0 \times 10^{-3} \mathrm{mm}\). Then, \(u_{y}\) further increases to \(7.0 \times 10^{-3} \mathrm{mm}\) with the step size \(1.0 \times 10^{-4} \mathrm{mm}\). Afterwards, \(u_{y}\) decreases to zero with the step size \(-1.0 \times 10^{-3} \mathrm{mm}\) during the unloading stage to finish a complete loading-unloading cycle containing totally 32 load steps. For the shear test, during the loading stage, the horizontal displacement \(u_{x}\) increases to \(15.0 \times 10^{-3} \mathrm{mm}\) with the step size \(1.0 \times 10^{-3} \mathrm{mm}\). During the unloading stage, the horizontal displacement decreases to zero with the step size \(-1.0 \times 10^{-3} \mathrm{mm}\) to finish a complete loading-unloading cycle with totally 30 load steps. For the bending test, the vertical displacement \(u_{y}\) first increases to \(1.0 \mathrm{mm}\) with the step size \(2.0 \times 10^{-2} \mathrm{mm}\) during the loading stage. Then, \(u_{y}\) decreases to zero with the step size \(-2.0 \times 10^{-2} \mathrm{mm}\) during the unloading stage with totally 100 load steps for the whole loading-unloading cycle.

For the staggered scheme using the alternate minimization, inside one iteration, the nonlinear displacement sub-problem is firstly solved via the Newton-Raphson method. Then, the obtained displacement field is fixed and the linear phase-field sub-problem is solved subsequently. For the two monolithic schemes, each iteration only involves one solve of the monolithic linear system. For the proposed L-BFGS-B solver, the detailed convergence criteria are listed under Section 3.1. For the staggered scheme and the limited-memory BFGS scheme using the history variable to enforce irreversibility, the following residual-based convergence criteria are adopted:

\[
\| \pmb {r}_u\| _2 < \mathrm{tol}, \quad \| \pmb {r}_d\| _2 < \mathrm{tol}, \quad \| \Delta \pmb {u}\| _2 < \mathrm{tol}, \quad \| \Delta d\| _2 < \mathrm{tol},
\]

where \(\pmb{r}_{u}\) and \(\pmb{r}_{d}\) represent the residuals of the displacement and the phase-field sub-problems, and \(\Delta \pmb {u}\) and \(\Delta d\) represent the increments of the displacement and the phase-field solutions during one iteration. To ensure a fair comparison regarding the solver convergence behaviors, all the three solvers adopt the same tolerance value tol \(= 10^{-6}\). Notice that for the proposed L-BFGS-B approach, besides checking the \(l_{2}\)-norms of the residual in the form of the projected gradient and the solution increment, the active set has to remain unchanged between two consecutive iterations.

Similar to the proposed L-BFGS-B solver, the source code and the input files for the staggered approach (A.M.-hist.) and the limited-memory BFGS approach (L-BFGS-hist.) are also hosted on GitHub in support of reproducible research. Fig. 22 reports the number of iterations required for convergence in each load step during a full loading-unloading cycle for the three (tensile, shear, bending) test problems shown in Fig. 21. Both the proposed L-BFGS-B scheme and the BFGS scheme using the history variable for irreversibility (L-BFGS-hist.) require fewer iterations than the staggered scheme (A.M.-hist.) for the shear problem in Fig. 21(b) and the bending problem in Fig. 21(c). For the tensile problem, even though the staggered approach (A.M.-hist.) requires fewer alternate minimization (outer-level) iterations, considering the fact that each alternate minimization iteration involves a nonlinear solve of the displacement sub-problem using 2 to 5 Newton iterations, the total number of required linear solves in the staggered approach is still larger than the counterparts of the two monolithic schemes (L-BFGS-B and L-BFGS-hist.). Between the two monolithic schemes, the proposed L-BFGS-B solver requires slightly more iterations inside each load step than the BFGS solver using the history variable approach (L-BFGS-hist.). This is not surprising because besides satisfying the tolerances of the residuals and the solution increments, the convergence criteria of the former also require the convergence of all the active constraints. Moreover, the iterates have to always stay feasible due to the gradient projection in the proposed L-BFGS-B method. On the upside, the proposed L-BFGS-B solver has the advantage of ensuring that the phase-field is always between 0 and 1, whereas the solvers based on the history variable approach obtain phase-field values larger than 1 or smaller than 0, which is a known drawback for this technique.

It is worth emphasizing that the number of iterations reported in the above three test problems should not be directly compared to the counterparts of similar problems in the literature. This is because there are always differences in problem settings, choices of phase-field models and linear solvers, details of numerical implementations, and so on, which could significantly impact the solver convergence behavior. For instance, all the simulations in this work do not use any regularization or artificial parameter \((k = 0\) in Eq. (5)) in the phase-field model. Experiences show that the introduction of a small viscosity parameter could significantly reduce the required iterations for brittle crack problems. Additionally, the choice of convergence criteria also has a large impact on the required number of iterations for convergence.

## 6. Conclusions

In this paper, a phase-field monolithic scheme based on the gradient projection method and the limited-memory BFGS approach is proposed to simulate the crack propagation in brittle materials. As a type of active set method, the gradient projection method is particularly attractive to handle bound constraints, or box constraints, imposed on the primary variables. This method has the advantage of allowing the rapid change of active constraints between iterations and can compute the projected gradient with a negligible computational cost. Therefore, it is an ideal method to enforce the phase-field irreversibility condition to ensure the thermodynamic consistency required by the phase-field formulation. On a different front, it is well known that the underlying energy functional for the phase-field crack formulation is non-convex, causing convergence issues in the Newton-based method. As a type of quasi-Newton method, the limited-memory BFGS (L-BFGS) method can effectively overcome the numerical difficulties associated with the non-convexity of the energy functional. Furthermore, the limited-memory feature can avoid storing the fully dense BFGS matrix, making it practical for large-scale finite element simulations.

The proposed phase-field monolithic scheme has three major ingredients. First, a compact representation of the updated BFGS matrix is adopted as the limited-memory feature to avoid storing fully dense matrix in the memory. Instead, only a few vector-pairs obtained from the most recent iterations need to be stored, making the method practical for finite element simulations. Second, the gradient projection method is used to locate the generalized Cauchy point on the piecewise linear path formed by the projected gradient and the box constraints. The variables at the generalized Cauchy point that are located on the boundary of the box constraints form the active set. These variables belonging to the active set are fixed in the subsequent subspace minimization problem, which is used to calculate the remaining free variables at the generalized Cauchy point. Lastly, in order to solve the subspace minimization problem, a primal approach based on the direct matrix factorization and the conjugate gradient method as well as a dual approach using the Lagrange multipliers are discussed.

Several 2D and 3D examples are provided to demonstrate the capabilities of the proposed monolithic approach in modeling crack propagation under cyclic loading conditions. In these numerical examples, the proposed monolithic scheme is further combined with a predictor-corrector adaptive mesh refinement technique to alleviate the computational cost. The monolithic scheme presented in this paper provides a unified framework to overcome the numerical difficulties arising from the non-convex energy functional, effectively enforce the phase-field irreversibility to ensure the thermodynamic consistency, and alleviate the heavy computational cost through adaptive mesh refinement in 2D and 3D phase-field crack simulations.

This work also provides a detailed comparison regarding the convergence behavior of the proposed L-BFGS-B method with two other phase-field solvers using the history variable approach for irreversibility. In the author's opinion, this is the fair way to compare solver behaviors, since both the problem settings (finite element mesh, load step, material parameters, etc.) and the numerical implementations (choice of linear solver, convergence criteria and tolerance, software environment) are kept the same. In the phase-field literature, there are multiple phase-field solvers that aim to address the numerical issues related to the non-convex energy functional and the phase-field irreversibility, such as the modified Newton method by Wick [6], the augmented Lagrangian method by Wheeler et al. [17], the interior-point method by Wambacq et al. [19], and the primal-dual active set method by Heister et al. [4], to name a few. However, it is quite challenging to fairly compare the convergence behavior of all these solvers. As a long-term project, the author plans to adopt a similar approach to compare the aforementioned methods, which could be useful to the phase-field community by providing a fair evaluation about various methods.

## CRediT authorship contribution statement

**Tao Jin:** Writing – review & editing, Writing – original draft, Visualization, Software, Resources, Methodology, Investigation, Funding acquisition, Formal analysis, Conceptualization.

## Declaration of competing interest

The authors declare that they have no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

## Acknowledgments

Tao Jin is supported by the Natural Sciences and Engineering Research Council of Canada (NSERC) under the Discovery Grants Program (funding reference number: RGPIN-2021-02561). Their support is greatly appreciated.

## Appendix

Let \(\mathbf{U}\in \mathbb{R}^{n\times p}\) and \(\mathbf{V}\in \mathbb{R}^{n\times p}\) represent two matrices for some \(p\) between 1 and \(n\). If

\[
\hat{\mathbf{A}} = \mathbf{A} + \mathbf{U}\mathbf{V}^{\mathrm{T}},
\]

according to the Sherman–Morrison–Woodbury formula [45], \(\hat{\mathbf{A}}\) is non-singular if and only if \((\mathbf{I} + \mathbf{V}^{\mathrm{T}}\mathbf{A}^{-1}\mathbf{U})\) is non-singular. In this case,

\[
\hat{\mathbf{A}}^{-1} = \mathbf{A}^{-1} - \mathbf{A}^{-1}\mathbf{U}\left(\mathbf{I} + \mathbf{V}^{\mathrm{T}}\mathbf{A}^{-1}\mathbf{U}\right)^{-1}\mathbf{V}^{\mathrm{T}}\mathbf{A}^{-1}.
\]

## Data availability

All the source code and input files are hosted at https://github.com/taojinlln/Phasefield_gradient_projection_monolithic_solver.

## References

[1] C. Miehe, F. Welschinger, M. Hofacker, Thermodynamically consistent phase-field models of fracture: Variational principles and multi-field FE implementations, Internat. J. Numer. Methods Engrg. 83 (10) (2010) 1273-1311, http://dx.doi.org/10.1002/nme.2861.
[2] C. Miehe, M. Hofacker, F. Welschinger, A phase field model for rate-independent crack propagation: Robust algorithmic implementation based on operator splits, Comput. Methods Appl. Mech. Engrg. 199 (45) (2010) 2765-2778, http://dx.doi.org/10.1016/j.cma.2010.04.011.
[3] M. Ambati, T. Gerasimov, L. De Lorenzis, A review on phase-field models of brittle fracture and a new fast hybrid formulation, Comput. Mech. 55 (2) (2015) 383-405, http://dx.doi.org/10.1007/s00466-014-1109-y.
[4] T. Heister, M.F. Wheeler, T. Wick, A primal-dual active set method and predictor-corrector mesh adaptivity for computing fracture propagation using a phase-field approach, Comput. Methods Appl. Mech. Engrg. 290 (2015) 466-495, http://dx.doi.org/10.1016/j.cma.2015.03.009.
[5] J.-Y. Wu, A unified phase-field theory for the mechanics of damage and quasi-brittle failure, J. Mech. Phys. Solids 103 (2017) 72-99, http://dx.doi.org/10.1016/j.jmps.2017.03.015.
[6] T. Wick, Modified Newton methods for solving fully monolithic phase-field quasi-static brittle fracture propagation, Comput. Methods Appl. Mech. Engrg. 325 (2017) 577-611, http://dx.doi.org/10.1016/j.cma.2017.07.026.
[7] E. Martinez-Pañeda, A. Golahmar, C.F. Niordson, A phase field formulation for hydrogen assisted cracking, Comput. Methods Appl. Mech. Engrg. 342 (2018) 742-761, http://dx.doi.org/10.1016/j.cma.2018.07.021.
[8] J.-Y. Wu, Y. Huang, V.P. Nguyen, On the BFGS monolithic algorithm for the unified phase field damage theory, Comput. Methods Appl. Mech. Engrg. 360 (2020) 112704, http://dx.doi.org/10.1016/j.cma.2019.112704.
[9] L. Svolos, C.A. Bronkhorst, H. Waisman, Thermal-conductivity degradation across cracks in coupled thermo-mechanical systems modeled by the phase-field fracture method, J. Mech. Phys. Solids 137 (2020) 103861, http://dx.doi.org/10.1016/j.jmps.2019.103861.
[10] A. Costa, M. Cusini, T. Jin, R. Settgast, J.E. Dolbow, A multi-resolution approach to hydraulic fracture simulation, Int. J. Fract. 237 (1) (2022) 165-188, http://dx.doi.org/10.1007/s10704-022-00662-y.
[11] J. Oliver, A.E. Huespe, E. Samaniego, E.W.V. Chaves, Continuum approach to the numerical simulation of material failure in concrete, Int. J. Numer. Anal. Methods Geomech. 28 (7-8) (2004) 609-632, http://dx.doi.org/10.1002/nag.365.
[12] F. Armero, J. Kim, Three-dimensional finite elements with embedded strong discontinuities to model material failure in the infinitesimal range, Internat. J. Numer. Methods Engrg. 91 (12) (2012) 1291-1330, http://dx.doi.org/10.1002/nme.4314.
[13] T. Jin, H.M. Mourad, C.A. Bronkhorst, A comparative study of shear band tracking strategies in three-dimensional finite elements with embedded weak discontinuities, Finite Elem. Anal. Des. 155 (2019) 11-31, http://dx.doi.org/10.1016/j.finel.2018.11.001.
[14] G. Francfort, J.-J. Marigo, Revisiting brittle fracture as an energy minimization problem, J. Mech. Phys. Solids 46 (8) (1998) 1319-1342, http://dx.doi.org/10.1016/S0022-5096(98)00034-9.
[15] B. Bourdin, G. Francfort, J.-J. Marigo, Numerical experiments in revisited brittle fracture, J. Mech. Phys. Solids 48 (4) (2000) 797-826, http://dx.doi.org/10.1016