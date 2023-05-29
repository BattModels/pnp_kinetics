from dragonfly import load_config_file, maximise_function
from argparse import Namespace
import numpy as np
import pickle, os, re, sys
import time, math
from pexpect.exceptions import TIMEOUT
from funcs_pnp import *
import yaml

def analytical_model(x, shift=0.01, factor=25): 
    # Using sigmoid voltammogram for simplicity, normalized
    # factor and shifts for better comparison in the Vapp values
    return (1/(1 + np.exp(-(x + shift)*factor))) - 1

def loss_error(run_dir):
    iv_file = run_dir+'/iv_data.txt'
    print(iv_file)
    iv_data = np.loadtxt(iv_file, skiprows=1)
    
    argids = np.argsort(iv_data[:,1])
    Vapp_list = iv_data[:,1][argids]
    I = iv_data[:,0][argids]
    y_predicted = I/np.abs(I[0])
    y_actual = analytical_model(Vapp_list)
    MSE = np.square(np.subtract(y_actual,y_predicted)).mean() 
    RMSE = math.sqrt(MSE)
    return -RMSE

def experiment(x):
    #x is list of variables defined in json file, in order
    #Take in variables and return objective
    #First try optimizing ko
    
    ## Load base configuration file
    loc = ''
    with open(loc+'config.yml', 'r') as file:
	    constants = yaml.safe_load(file)
        
    try:
        ## Change involved parameters in the config file 
        run_dir = 'A_'+str(round_it(x[0], 3))
        constants['MHC_prefactor'] = float(x[0])
        write_yml(constants, fname='config.yml')
        os.system('./run.sh '+run_dir)
        objective=loss_error(run_dir)
    except(TIMEOUT):
	    objective=experiment(x)
            
    print(x)
    print(objective)
    return objective

def run_opt(SAVE_AND_LOAD_PROGRESS=True):
    config=load_config_file('optimize.json')
    options=Namespace(build_new_model_every=3,report_results_every=1,report_model_on_each_build=True)
    if SAVE_AND_LOAD_PROGRESS:
        options.progress_load_from_and_save_to = 'dragonfly_opt.p'
        options.progress_save_every = 1

    ##Optimize
    max_num_evals=10 #Optimization budget
    #opt_method="rand"
    opt_method="bo" #"ea"
    opt_pt, opt_val, history = maximise_function(
	experiment, config.domain, max_num_evals, config=config, options=options,opt_method=opt_method)
    print("Optimum {} at {}".format(opt_val,opt_pt))

if __name__ == '__main__':
    run_opt(SAVE_AND_LOAD_PROGRESS=True)
