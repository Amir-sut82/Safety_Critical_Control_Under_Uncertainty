# ISS Safety-Critical Control under Actuation Uncertainty

A robust control framework for a nonlinear MIMO attitude system with uncertain actuator effectiveness.

This project studies safety-critical stabilization when the control input is affected by unknown multiplicative uncertainties. Two complementary approaches are investigated:

1. **ISS-based modular adaptive control** using parameter estimation and CLF/CBF analysis.
2. **Robust min-max control** that provides safety guarantees under bounded actuator uncertainty.

The objective is to compare accuracy, feasibility, control effort, robustness, and computational complexity between adaptive and robust safety-critical approaches.

---

# 1. Problem Formulation

The nonlinear uncertain system is represented as:

$$
\dot{x}=f(x)+g(x)(I+\Delta)u
$$

where:

- $x$ is the system state
- $u$ is the control input
- $f(x)$ represents nonlinear dynamics
- $g(x)$ represents the input matrix
- $\Delta$ represents unknown actuator uncertainty

The considered system is a six-state attitude model:

$$
x=[\phi,\theta,\psi,p,q,r]^T
$$

where the first three states represent Euler angles and the last three states represent body angular velocities.

The uncertainty is modeled as a bounded multiplicative actuator uncertainty.

---

# 2. Modular Adaptive Control with ISS

The first architecture separates:

- parameter estimation
- controller design

The estimation error is treated as an input to the closed-loop system.

The Input-to-State Stability (ISS) formulation guarantees bounded behavior despite imperfect parameter estimation.

The controller uses an ISS-CLF formulation:

$$
\dot V \leq -cV+\gamma(||\tilde{\theta}||)
$$

where:

- $V$ is the Lyapunov function
- $\tilde{\theta}$ is the parameter estimation error
- $\gamma$ represents the effect of estimation uncertainty

---

# 3. Parameter Estimation Algorithms

A history-stack based learning framework is implemented to improve parameter convergence.

The regression model is written as:

$$
Y_k=\Phi_k\theta+\nu_k
$$

A history stack stores previous measurements and improves excitation.

Implemented estimators:

- Gradient adaptation
- Recursive Least Squares (RLS)
- Constant forgetting RLS
- Variable forgetting RLS

The variable forgetting RLS improves the balance between convergence speed and robustness by adapting the forgetting factor according to prediction error.

---

# 4. Safety-Critical Control using CLF and HOCBF

## Control Lyapunov Function

The CLF condition ensures convergence:

$$
\dot V(x)+cV(x)\leq0
$$

The controller is obtained through quadratic programming optimization.

---

## Higher Order Control Barrier Function

The safety objective is maintaining the pitch angle within a safe region:

$$
C=\{x:h(x)\geq0\}
$$

Because the safety constraint has relative degree two, a Higher Order CBF is used:

$$
\psi_0=h(x)
$$

$$
\psi_1=\dot{\psi}_0+\alpha_1(\psi_0)
$$

$$
\psi_2=\dot{\psi}_1+\alpha_2(\psi_1)
$$

The final safety condition is:

$$
\psi_2(x,u)\geq0
$$

---

# 5. Robust Control under Actuator Uncertainty

The second approach directly handles uncertainty through worst-case optimization.

The robust CLF condition considers the maximum possible uncertainty effect:

$$
\max_{\Delta\in\Omega} \Delta u^T \nabla V
$$

The uncertainty is converted into tractable quadratic programming constraints using support functions and convex approximations.

This provides guaranteed stability and safety without requiring parameter convergence.

---

# 6. Online SMID Uncertainty Reduction

A Set Membership Identification (SMID) method is used to reduce conservatism.

Instead of assuming a fixed uncertainty bound, SMID updates the certified uncertainty set online.

The adaptive uncertainty radius contracts as more information becomes available.

Reported results show that SMID reduces the uncertainty radius while maintaining safety guarantees.

---

# 7. Experimental Comparison

The project compares:

| Method | Main Advantage |
|-|-|
| ISS Adaptive Control | Exploits learnable uncertainty |
| Robust Min-Max Control | Immediate worst-case guarantees |
| SMID + Robust Control | Reduced conservatism with safety guarantees |

Main observations:

- Variable forgetting RLS provides a strong balance between convergence speed and robustness.
- Robust control guarantees safety under actuator uncertainty.
- Online uncertainty reduction decreases conservatism and control effort.

---

# 8. Key Contributions

This project demonstrates:

- Nonlinear uncertain system modeling
- ISS-based adaptive control
- History-stack parameter learning
- CLF-based stabilization
- HOCBF safety constraints
- Robust min-max optimization
- Online uncertainty set reduction

---

# 9. Repository Structure

```
ISS_action_uncertainity/
│
├── Controller implementations
├── Estimation algorithms
├── Simulation files
├── Results
└── README.md
```

---

# References

- Khalil, H. K. Nonlinear Systems.
- Ames, A. D. et al. Control Barrier Function Based Quadratic Programs for Safety Critical Systems.
- Adaptive and Robust Safety-Critical Control literature.
