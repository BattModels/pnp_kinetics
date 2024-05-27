import numpy as np
import os, sys, re, pickle, glob, yaml, mat73
from scipy.io import loadmat, savemat
from math import log10 , floor

## User functions
def read_yml(ymlfile):
    with open(ymlfile) as file:
        out_dict = yaml.load(file, Loader=yaml.FullLoader)
    return out_dict

def write_yml(dict_file, fname='test.yaml'):
    with open(fname, 'w') as file:
        yaml.dump(dict_file, file)

def write_jobfile(sh_name, entry_id, cores=64):
    with open(sh_name, "r") as f:
        txt = f.read()
    with open(sh_name, "w") as f:
            f.write(txt.format(entry_id, cores))

## User specifics
domain = 'AA'
Vapp_list = np.round(np.linspace(-0.5, 0.3, 41), 3)

dos_files = np.array(glob.glob('/ocean/projects/cts180021p/mbabar/PhD/tBLG_Carr/kp_tblg/MATLAB_vers/TBG_data_nr_48/ldos*.mat'))

## Theta input
#theta_list = [0.77, 1.1, 2, 3, 4.6]
theta_list = [1.1, 2.6]
ids = []
for th in theta_list:
    for idx in range(len(dos_files)):
        if th == np.double(re.split('th_|_nsamps_', dos_files[idx])[1]):
            print("theta = {}, idx = {}".format(th, idx))
            ids.append(idx)
            continue

print(ids)
## Run jobs

for i in range(len(dos_files[ids])):
    ldos_data = loadmat(dos_files[ids[i]])
    theta = np.round(np.squeeze(ldos_data['theta']), 3)
    os.system('mkdir '+str(theta))
    os.chdir(str(theta))
    os.system('cp ../job.sh .')
    os.system('cp ../pnp_par.py .')
    os.system('cp ../funcs_pnp.py .')
    os.system('cp ../get_current.py .')
    os.system('cp ../config.yml .')
    os.system('cp ../solve_vmg.sh .')
    os.system('cp ../run_pnp.py .')
    
    for j in range(len(Vapp_list)):
        print('Running DOS file : {} at Vapp = {}'.format(dos_files[ids[i]], Vapp_list[j]))
        rate_loc = '/ocean/projects/cts180021p/mbabar/PhD/PDE/PNP_solve/3D_vmg/ldos_incorp/k_data_nr_48/{}/k_data_Vapp_{}.mat'.format(str(theta), str(Vapp_list[j]))
        pref_loc = '/ocean/projects/cts180021p/mbabar/PhD/PDE/PNP_solve/3D_vmg/prefactor/pref_data/pref_sc_{}.mat'.format(str(theta))
       
        config_loc = 'config_th_{}_Vapp_{}_{}.yml'.format(str(theta), str(Vapp_list[j]), domain)   
        constants = read_yml(config_loc)
        constants['dos_file'] = str(dos_files[ids[i]])
        constants['rate_loc'] = rate_loc
        constants['MHC_prefactor'] = pref_loc
        constants['domain'] = domain
        #constants['sol_dir'] = os.getcwd() # Specified in solve_vmg.sh
        write_yml(constants, fname=config_loc)
        #write_jobfile('job.sh', str(theta), cores=32)
        #os.system('sbatch job.sh')
        print('Written config file: {}'.format(config_loc))

    os.chdir('../')
