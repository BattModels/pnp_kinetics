# Poisson-Nernst-Planck Solve

Solving coupled Poisson-Nernst-Planck equation in an axisymmetric nanopipette crossection using Fenics 2019.1.0.
We model four concentration species (K+, Cl-, Ru3+ and Ru2+), where Ruhex (Ruthenium hexamine) is actively underdoing redox at the bottom surface ($z=0$) and KCl is the supporting electrolyte. 

## Governing equations

1) $$\Delta \phi = \frac{-F}{\epsilon \epsilon_{o}} \sum_{j}z_{j}c_{j}$$

2) $$D_{i}\Delta c_{i} = -\frac{z_{i}eD_{i}}{k_{B}T}\left( c_{i}\nabla(\phi) \right) $$ 
where the four species are indexed by $i=1,2,3,4$.
Here $i=1$ : $K^{+}$
$i=2$ : $Cl^{-}$
$i=3$ : $Ru^{3+}$
$i=4$ : $Ru^{2+}$

Boundary Conditions (BCs), <\br>

1) $\phi = 0$ at $z=l$, $0<r<a$
2) $\phi = \phi_{o}$ at $z=0$, $0<r<a_{s}$

(Dirichlet) Charge neutrality at the top surface $z=l$, $0<r<a$ <\br>
3) $c_{1} = 100$ mM , $c_{2} = 106$ mM, $c_{3} = 2$ mM and $c_{4} = 0$ mM <\br>

Neumann condition <\br>
4) $$\frac{\partial c_{i}}{\partial n} = 0$$ (At the bottom and side surfaces if i=1,2. Only at the side surface if i=3,4)

Coupled Robin condition at the bottom surface (for i=3,4) <\br>
5) $$-D_{o}\frac{\partial c_{3}}{\partial n} = D_{o}\frac{\partial c_{4}}{\partial n} = c_{3}k_{red} - c_{4}k_{ox}$$

where $k_{red}$ and $k_{ox}$ are the rate constants derived from the selected rate theory.  

## Usage
The directory contains the following files,
```
funcs_pnp.py # Formulation of PNP equation in UFL 

run_pnp.py # Loads functions from funcs_pnp.py to iteratively converging to solution

check_sol.py # Checks final concentration (c) and potential (phi) solution on boundary conditions and convergence.

inital_c0.pkl # Loads initial guess vector for conc.

pippete_pnp.ipynb # Benchmark and testing in jupyter notebook on nanopipette mesh

test_pnp.ipynb # Testing in jupyter notebook on simple square mesh
```