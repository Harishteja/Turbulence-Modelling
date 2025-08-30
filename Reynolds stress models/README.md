
## Reynolds Stress Models (RSM)

Unlike eddy-viscosity models (e.g., k–ε, k–ω) that assume **isotropy** of turbulence, **Reynolds Stress Models (RSM)** directly solve transport equations for each component of the Reynolds stress tensor. This eliminates the isotropy limitation and provides a more detailed description of turbulent flows.

The general transport equation for the Reynolds stresses is:

$$
\frac{\partial}{\partial x_k} \left( \rho U_k \overline{u_i' u_j'} \right) = \mu \frac{\partial^2 \overline{u_i' u_j'}}{\partial x_k \partial x_k} + P_{ij} + \Phi_{ij} + D_{ij} - \rho \varepsilon_{ij}
$$

Where the terms on the RHS represent:

* **Viscous diffusion**
* **Production** of Reynolds stresses
* **Pressure redistribution** (split into slow and fast terms, modeled using Rotta and IP models)
* **Turbulent diffusion**
* **Dissipation**

Since these equations are complex, detailed formulations can be found in standard turbulence modeling references \[2].

Additionally, a separate **transport equation for the dissipation of turbulent kinetic energy** is required to close the RSM.

For near-wall treatment, this project implements **wall functions** within the Reynolds Stress Model (details are available in the attached codes).


