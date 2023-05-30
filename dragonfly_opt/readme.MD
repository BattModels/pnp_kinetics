# Bayesian optimization using Dragonfly 

All steps are same as for 3D mesh,

1. Set parameters in `config.yml` and overpotential range in `solve_vmg.sh`. 

2. Designate an analytical function or an interpolator of experimental data in `optimize.py` as ground truth for fitting. Presently a sigmoid function is defined. Also specify max iterations in `optimize.py`.

3. Specify the variables to be tuned and their ranges in discrete numerics in `optimize.json`. In the directory `A_opt/`, MHC-prefactor is varied from 0 to 100. Moreover, redox rates are derived from a Julia script using `ElectrochemicalKinetics.jl` package. The script used is defined in `config.yml`, and executed in `solve_vmg.sh`. The MHC-DOS rate model employs the DOS file in `dos_incorp/` for 2 deg twisted bilayer graphene.  

4. Execute `python optimize.py` : It guesses the first set of hyperparameters and makes a new directory. Then generates new config files with the given range of $`V_{app}`$, and obtains the IV curve for that set. The loss function in `optimize.py` outputs the RMSE error with the ground truth supplied. Based on the error, the optimizer chooses a new set of hyperparameters and runs the cycle again until max iteration is reached. Loss function per iteration is pickled in `dragonfly_opt.p` or can be individually calculated from `iv_data.txt` using `python get_rmse.py folder_name`, where folder_name contains `iv_data.txt`.

One can specify multiple variables in the `optimize.json` and make subsequent changes in `optimize.py` file for multi-variable optimization. 
 
5. Execute `python get_rmse.py folder_name` to calculate RMSE error with the analytical/experimental data. The directory `folder_name` should have a complete `iv_data.txt` file. 

As an example, the optimizer finds an accurate BV rate constant $`k_{o}`$ which fits the sigmoid function within 10 search iterations,

<img src="../../img/opt.png" alt="iv" width="500">
<img src="../../img/rmse.png" alt="iv" width="500">


## Notes

* Make bash files executable by running `chmod +x run.sh solve_vmg.sh process_vmg.sh`
