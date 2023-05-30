# MHC-LDOS model with PNP

Implementing local DOS into MHC rate theory as BC to PNP equation in a 3D mesh. The local DOS does not obey axisymmetry, for which 2D interpolants for $`k_{ox}`$, $`k_{red}`$, $`V_{dl}`$ are calculated at the bottom surface. The interpolants output the value of these variables at the specified $`(r, \theta)`$. Accordingly, `pnp_par.py` is modified to support the increased dimensionality of the problem. New 2D interpolant functions are added to `funcs_pnp.py`.

Rest all steps are similar to `dos_incorp/`,

1. Set location of mesh (`.xml`), LDOS (`.mat`) and other parameters in `config.yml`. The DOS matlab keys should be the same as the example file `ldos-th_2_nsamps_30_knum_24.mat`. Specify overpotential range in `solve_vmg.sh`.

2. Julia script `script.jl` is used to compute oxidation and reduction rates based on MHC-DOS theory for each $`V_{app}`$, and are updated in the `Sol/config_$$` file. The new config file is loaded to solve the PNP equations with new values to the boundary conditions at the bottom surface of the nanopippete. MHC-DOS and other helper functions are called from private fork ([ElectrochemicalKinetic.jl](https://github.com/mbabar09/ElectrochemicalKinetics/)). 

3. Execute `./run.sh Sol` : Runs PNP equation solve on a range of $`V_{app}`$ specified in `solve_vmg.sh`. Final currents and overpotentials are collected in `Sol/iv_data.txt`. Solution files in `.h5` format are also saved in `Sol/`, and can be used for further analysis.

4. Plot `iv_data.txt`

## Notes

* Make bash files executable by running `chmod +x run.sh solve_vmg.sh process_vmg.sh`