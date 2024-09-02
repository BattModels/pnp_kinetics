from __future__ import print_function
import numpy as np
from funcs_pnp import *
import os, sys, re, glob
from dolfin import *
from fenics import *
from numpy.random import rand
from numpy.linalg import norm
from scipy.io import loadmat, savemat
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib as mpl
import time, yaml, logging
import pandas as pd

# Load configuration 
config_loc = sys.argv[1] # Location of config file
#center_id = int(sys.argv[2]) # Center of nanopipette radius

with open(config_loc, 'r') as file:
    constants = yaml.safe_load(file)
    globals().update(constants)

fname = sol_file+'.h5'
# Load and scale mesh
mesh = Mesh(mesh_file)
m_coords = mesh.coordinates()
scale_factor = np.float64(1e-9)
m_coords[:, :] *= scale_factor
mesh.bounding_box_tree().build(mesh) # development version

a_s, a_t, z_t = get_mesh_geo(mesh)
print(a_s, a_t, z_t)

P1 = FiniteElement('P', tetrahedron, 1)
element = MixedElement([P1, P1, P1, P1, P1])
V = FunctionSpace(mesh, element)
sol = Function(V)
load_sol(fname, sol, mesh)
_phi, _c_1, _c_2, _c_3, _c_4 = split(sol)

rate_data = mat73.loadmat(constants["rate_loc"])
Vapp = np.round(np.double(rate_data["Vapp"]),3)
theta = np.round(np.double(rate_data["theta"]),3)
print('theta={}, Vapp={}'.format(theta, Vapp),'\n')

curr = get_current_3D(_c_3, mesh, D_o, a_s = a_s)

print('I = {}'.format(curr))

print("{}, I={}, Vapp={}".format(domain, curr, Vapp))

#with open(iv_file, 'a') as f:
#    f.write(domain+' '+str(curr)+' '+str(Vapp)+'\n')

iv_df = pd.read_csv(iv_file, sep=' ', header=None, names=['Domain', 'I', 'Vapp'], skiprows=1)


idx = list(iv_df.index[(iv_df['Domain'] == domain) & (iv_df['Vapp'] == Vapp)])

iv_df.at[int(np.squeeze(idx)), 'I'] = curr

print(iv_df.loc[idx])

iv_df.to_csv(iv_file, sep=' ', index=False)
