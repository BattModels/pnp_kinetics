import numpy as np
import os, sys, re, pickle, glob, yaml, mat73
from scipy.io import loadmat, savemat
from math import log10 , floor

## User functions
def round_it(x, sig):
    return round(x, sig-int(floor(log10(abs(x))))-1)

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

## Load config data
config_loc = "config.yml" # ARGS[1]

dos_files = glob.glob('/ocean/projects/cts180021p/mbabar/PhD/tBLG_Carr/kp_tblg/MATLAB_vers/TBG_data_nr_48/ldos*.mat')

for i in range(len(dos_files)):
    print('Running DOS file : {}'.format(dos_files[i]))
    ldos_data = loadmat(dos_files[i])
    theta = np.round(np.squeeze(ldos_data['theta']), 3)
    os.system('mkdir '+str(theta))
    os.chdir(str(theta))
    os.system('cp ../job.sh .')
    os.system('cp ../eta_run_script.jl .')
    os.system('cp ../config.yml .')
    constants = read_yml(config_loc)
    constants['dos_file'] = dos_files[i]
    constants['sol_dir'] = os.getcwd()
    write_yml(constants, fname=config_loc)
    write_jobfile('job.sh', str(theta), cores=64)
    os.system('sbatch job.sh')
    print('Created job file')
    os.chdir('../')
