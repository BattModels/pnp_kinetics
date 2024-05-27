import numpy as np
import os, sys, re, glob
import pickle, time
from funcs_pnp import *
from numpy.random import rand
from numpy.linalg import norm
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib as mpl
from scipy.io import loadmat, savemat
from mpi4py import MPI
import yaml, mat73

num_cores = 17
config_loc = 'config.yml'

constants = read_yml(config_loc)
globals().update(constants)

with open(iv_file, 'w') as f:
    f.write('idx I Vapp Center(x) Center(y)\n')

ldos_data = loadmat(constants["dos_file"])
rscx = ldos_data["rscx"]
rscy = ldos_data["rscy"]
theta = np.squeeze(ldos_data["theta"])
moir_l, sc_area = tblg_dim(theta)
moir_a = np.pi/3
coords = moir_l*np.array([rscx.flatten(), rscy.flatten()])
coords_data = {'theta': theta, 'coords': coords}
savemat("moire_coords.mat", coords_data)

for j in range(len(coords.T)): #len(coords.T)
    t0 = time.time()
    center = np.array(coords[:, j])
    print('\nCenter = {}'.format(center))
    os.system('mpirun -n '+str(num_cores)+' python pnp_par.py '+config_loc+' '+str(j))
    os.system('python get_current.py '+config_loc+' '+str(j))
    t1 = time.time()
    print('Time for iteration {} : {} sec'.format(j, t1-t0))

