# NMPDE Project - Cabiati, Muscarnera, Scarpelli
Numerical Methods for Partial Differential Equations Project by Cabiati, Muscarnera and Scarpelli. Academic Year 2022-2023

# Description 
The project is an implementation of a Navier-Stokes solver.
## Strong formulation
$$ \begin{cases}
  -\mu \nabla \vec{u} + (\rho \vec{u} \cdot \nabla) \vec{u} + \nabla p = \vec{f} & \text{in }\Omega \\
  \nabla \cdot \vec{u} = 0 & \text{in }\Omega \\
  \mu \nabla \vec{u} \cdot \widehat{\vec{n}} - p \cdot \widehat{\vec{n}} = -p_\text{out} \cdot \widehat{\vec{n}} & \text{on } \Gamma_0 \\
  \vec{u} = \vec{u}_\text{in} & \text{on } \Gamma_1 \\
  \vec{u} \cdot \widehat{\vec{n}} = \vec{0} & \text{on } \Gamma_2 \cup \Gamma_3 \cup \Gamma_4 \cup \Gamma_5 \\
  (\mu \nabla \vec{u} \cdot \widehat{\vec{n}} - p  \cdot \widehat{\vec{n}}) \cdot \widehat{\vec{t}}_i = \vec{0} & \text{on } \Gamma_i,\ i=2,3,4,5
\end{cases} $$

# References
- <https://nora.nerc.ac.uk/id/eprint/519363/1/Richardson.pdf> for the velocity
