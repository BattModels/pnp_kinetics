# Poisson-Nernst-Planck Solve on 2D and 3D Nanopippete Meshes

Solving coupled Poisson-Nernst-Planck equation in an axisymmetric nanopipette crossection using Fenics 2019.1.0. 
For illustration, we model four concentration species (K$`^{+}`$, Cl$`^{-}`$, Ru$`^{3+}`$ and Ru$`^{2+}`$), where Ruhex (Ruthenium hexamine) is actively underdoing redox at the bottom surface ($`z=0`$) with the electrode (twisted bilayer graphene). KCl is the supporting electrolyte. 

## Mesh

3D picture of the nanopippete mesh is shown,

$`a_{s} = 100`$ nm is the bottom radius

$`a_{t} = 3.5\times10^{-6}`$ m is the bottom radius

$`z_{t} = 2\times10^{-5}`$ m is the top height

## Governing equations

The following partial differential equations are solved simultaneously,

1) $$ \Delta \phi = \frac{-F}{\epsilon \epsilon_{o}} \sum_{j}z_{j}c_{j} $$

2) $$D_{i}\Delta c_{i} = -\frac{z_{i}eD_{i}}{k_{B}T}\left( \nabla \left( c_{i}\nabla \phi\right) \right) $$

where the four species are indexed by $`i=1,2,3,4`$.  
Here,  
$`i=1`$ : K$`^{+}`$  
$`i=2`$ : Cl$`^{-}`$  
$`i=3`$ : Ru$`^{3+}`$  
$`i=4`$ : Ru$`^{2+}`$  

Boundary Conditions (BCs), 

1) $` \phi = 0 `$ at $` z=l `$, $` 0 \le r \le a_{s}`$
2) $`\phi = \phi_{o}`$ at $`z=0`$, $` 0 \le r \le a_{s}`$

At the top surface, the concentration of K$`^{+}`$ and Ru$`^{3+}`$ are fixed and there is no redox reaction. Cl$`^{-}`$ concentration is obtained from charge neutrality. 
These constraints are set Dirichlet boundary conditions at the top surface $`z=l`$, $` 0 \le r \le a_{t}`$.

3) $`c_{1} = 100`$ mM , $`c_{2} = 106`$ mM, $`c_{3} = 2`$ mM and $`c_{4} = 0`$ mM  

Neumann condition  
4) $`\frac{\partial c_{i}}{\partial n} = 0`$ (At the bottom and side surfaces if $`i=1,2`$. Only at the side surface if $`i=3,4`$)  

Coupled Robin condition at the bottom surface (for i=3,4)  
5) $`-D_{o}\frac{\partial c_{3}}{\partial n} = D_{r}\frac{\partial c_{4}}{\partial n} = c_{3}k_{red} - c_{4}k_{ox}`$  

where $`k_{red}`$ and $`k_{ox}`$ are the rate constants derived from the selected rate theory, and $`D_{o}`$ and $`D_{r}`$ are diffusivity coefficients of the oxidized and reduced species respectively.


## Dependencies


1. Make new conda environment with legacy version of fenics (2019)
```
conda create -n fenicsproject -c conda-forge fenics
source activate fenicsproject
```
2. Install dependencies; gmsh, meshio, dragonfly, shyaml, shapely, mat73
```
pip install --upgrade gmsh
pip install meshio
pip install dragonfly-opt -v
pip install shyaml
pip install ase
pip install mat73
```
3. [Install Julia](https://julialang.org/downloads/). Clone private fork of [ElectrochemicalKinetic.jl](https://github.com/mbabar09/ElectrochemicalKinetics/) to local system. Make a new environment and add packages.
```
git clone https://github.com/mbabar09/ElectrochemicalKinetics/ 
cd ElectrochemicalKinetics.jl
activate .
add ArgParse, YAML, CSV, MAT, Glob, Interpolations, QuadGK, ClusterManagers, SharedArrays, Distributed, Zygote, SpecialFunctions
add Dierckx, Optim, NLsolve, DataStructures, Plots
```

## Usage

The following files are common in all directories,

1. `funcs_pnp.py` : Helper functions for solve, saving files and post processing. 
2. `pnp_par.py` : Main fenics code which solves the PNP equation. It loads mesh, constants from config.yml, specifies boundary conditions and the weak form. The code works in parallel with MPI.
3. `solve_vmg.sh` : Main bash file that runs pnp_par.py while attaching necessary rate model data or calculating it using ElectrochemicalKinetics.jl. It makes directories and appends necessary paths in the config file.
4. `get_current.py` : Post processing (PP) code for calculating the total current from redox reaction at the bottom surface for the given overpotential `Vapp`
5. `process_vmg.sh` : Post processing bash file for running PP codes like get_current.py. User can add user-specific PP codes. 
6. `run.sh` : Executable file that runs solve_vmg.sh 

These files can be used for different calculations 
```
funcs_pnp.py # Formulation of PNP equation in UFL 

run_pnp.py # Loads functions from funcs_pnp.py to iteratively converging to solution

check_sol.py # Checks final concentration (c) and potential (phi) solution on boundary conditions and convergence.

inital_c0.pkl # Loads initial guess vector for conc.

pippete_pnp.ipynb # Benchmark and testing in jupyter notebook on nanopipette mesh

test_pnp.ipynb # Testing in jupyter notebook on simple square mesh
```

## Notes 
-Use your Julia executable path where .jl script is run. 
-Change the $num_cores variable in .sh files for assigning different number of processes for parallelization.

## Contributing

We are happy to take pull requests for new features, new models, etc. by pull request! It is suggested (though not required) that you first open an issue to discuss.
