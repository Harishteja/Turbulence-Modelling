In this project, the Finite Volume Method (FVM) is used to numerically solve the governing equations for a fully developed turbulent channel flow using MATLAB. The turbulence closure is modeled using Eddy Viscosity Models, specifically the k–ε and k–ω models. The results are compared with Direct Numerical Simulation (DNS) data from Kim, Moser, and Mansour (1999). MATLAB codes developed for the project are available in the MATLAB Codes folder, and the results are summarized in the reports. A brief explanation of the governing equations, turbulence models, and the FVM discretization is provided below, with detailed discussions available in [2].

---

## Governing Equations

For a **fully developed turbulent channel flow**:

* \$x\$: streamwise direction
* \$y\$: wall-normal direction
* \$z\$: spanwise direction

**Assumptions:**

* **Statistical stationarity:** \$\dfrac{\partial \overline{\phi}}{\partial t} = 0\$ for any mean quantity \$\overline{\phi}\$
* **Spanwise homogeneity:** \$\dfrac{\partial \overline{\phi}}{\partial z} = 0\$

Applying these, the **Reynolds-Averaged Navier–Stokes (RANS)** equations with the Boussinesq approximation reduce to:

**Continuity:**

```math
\frac{\partial U}{\partial x} + \frac{\partial V}{\partial y} = 0
```

**Momentum (x-direction):**

```math
U \frac{\partial U}{\partial x} + V \frac{\partial U}{\partial y}
= -\frac{1}{\rho} \frac{\partial P}{\partial x}
+ \frac{\partial}{\partial x} \Big[ (\nu + \nu_t) \frac{\partial U}{\partial x} \Big]
+ \frac{\partial}{\partial y} \Big[ (\nu + \nu_t) \frac{\partial U}{\partial y} \Big]
```

**Momentum (y-direction):**

```math
U \frac{\partial V}{\partial x} + V \frac{\partial V}{\partial y}
= -\frac{1}{\rho} \frac{\partial P}{\partial y}
+ \frac{\partial}{\partial x} \Big[ (\nu + \nu_t) \frac{\partial V}{\partial x} \Big]
+ \frac{\partial}{\partial y} \Big[ (\nu + \nu_t) \frac{\partial V}{\partial y} \Big]
```

For a **fully developed flow**:

* \$\dfrac{\partial \overline{U}}{\partial x} = 0\$
* \$V = W = 0\$

This simplifies to a **1D equation**:

```math
0 = -\frac{1}{\rho} \frac{\partial P}{\partial x}
+ \frac{\partial}{\partial y} \Big[ (\nu + \nu_t) \frac{\partial U}{\partial y} \Big]
```

---

## Turbulence Modeling

In **tensor notation**, the RANS equations are written as:

```math
\frac{\partial \overline{u_i}}{\partial t} +
\overline{u_j} \frac{\partial \overline{u_i}}{\partial x_j}
= -\frac{1}{\rho} \frac{\partial \overline{p}}{\partial x_i}
+ \nu \frac{\partial^2 \overline{u_i}}{\partial x_j \partial x_j}
+ \frac{\partial \overline{u_i' u_j'}}{\partial x_j}
```

The additional term \$\overline{u\_i' u\_j'}\$ represents the **Reynolds stresses**, which are unknown.

* To close the equations, these stresses must be modeled.
* This process is called the **turbulence closure problem**.
* Common approaches include the **k–ε** and **k–ω** turbulence models.

---

## Eddy Viscosity Models

One widely used approach to turbulence modeling is the **Boussinesq approximation**, which expresses the Reynolds stresses as:

```math
\overline{u_i' u_j'} =
- \nu_t \left( \frac{\partial \overline{u_i}}{\partial x_j}
+ \frac{\partial \overline{u_j}}{\partial x_i} \right)
+ \frac{2}{3} k \, \delta_{ij}
```

This reduces the **six unknown Reynolds stresses** to just two quantities:

* **Turbulent kinetic energy** (\$k\$)
* **Turbulent (eddy) viscosity** (\$\nu\_t\$)

Substituting the Boussinesq approximation into the RANS equations gives the **general form of eddy viscosity models**.

* One-equation models: e.g., **Prandtl’s mixing length model**
* Two-equation models: e.g., **k–ε**, **k–ω**

In this project, we focus on the **two-equation models (k–ε, k–ω)**.

>  Limitation: By reducing six unknowns to two, these models assume turbulence is isotropic, and therefore cannot fully capture anisotropic effects.

---
### k–ε Model

The **k–ε model** solves two transport equations:

1. A transport equation for turbulent kinetic energy $k$ (derived from the exact $k$ equation with modeling assumptions).  
2. A transport equation for the **dissipation rate** of turbulent kinetic energy, $\epsilon$.

The turbulent viscosity is computed from $k$ and $\epsilon$:

$$
\nu_t = C_\mu \frac{k^2}{\epsilon}
$$

**Governing equations (fully developed turbulent channel flow):**

Momentum (streamwise):

$$
\frac{\partial}{\partial y}\!\left[(\nu+\nu_t)\frac{\partial U}{\partial y}\right] -\frac{1}{\rho}\frac{\partial P}{\partial y}=0
$$

Turbulent kinetic energy $(k)$:

$$
\frac{\partial}{\partial y}\!\left[\left(\nu+\frac{\nu_t}{\sigma_k}\right)\frac{\partial k}{\partial y}\right] + P_k - \epsilon = 0 
$$

Dissipation rate $(\epsilon)$:

$$
\frac{\partial}{\partial y}\!\left[\left(\nu+\frac{\nu_t}{\sigma_\epsilon}\right)\frac{\partial \epsilon}{\partial y}\right] + C_1 \frac{\epsilon}{k} P_k - C_2 \frac{\epsilon^2}{k} = 0
$$

Production term:

$$
P_k = \nu_t \left(\frac{\partial U}{\partial y}\right)^2
$$

---

### k–ω Model

The **k–ω model** uses the **specific dissipation rate** (\$\omega\$) instead of \$ε\$:

```math
\omega = \frac{\epsilon}{\beta^* k}
```

Turbulent viscosity:

```math
\nu_t = \frac{k}{\omega} 
```

**Equations (for fully developed turbulent channel flow):**

Momentum:

```math
\frac{\partial}{\partial y} \Big[ (\nu + \nu_t) \frac{\partial U}{\partial y} \Big]
- \frac{1}{\rho} \frac{\partial P}{\partial y} = 0 
```

Turbulent kinetic energy (\$k\$):

```math
\frac{\partial}{\partial y} \Big[ \Big(\nu + \frac{\nu_t}{\sigma_k}\Big) \frac{\partial k}{\partial y} \Big]
+ P_k - \beta^* k \omega = 0 
```

Specific dissipation rate (\$\omega\$):

```math
\frac{\partial}{\partial y} \Big[ \Big(\nu + \frac{\nu_t}{\sigma_\omega}\Big) \frac{\partial \omega}{\partial y} \Big]
+ \alpha \frac{\omega}{k} P_k - \beta \omega^2 = 0 
```

Production term:

```math
P_k = \nu_t \left( \frac{\partial U}{\partial y} \right)^2 
```

---






