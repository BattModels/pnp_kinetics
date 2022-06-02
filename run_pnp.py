from __future__ import print_function
from funcs_pnp import *
import numpy as np
import os, sys, re, glob
from dolfin import *
from fenics import *
import pickle
from numpy.random import rand
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib as mpl
font = {'family' : 'Times New Roman',
        'weight' : 'regular',
        'size'   : 16}

mpl.rc('font', **font)
mpl.rcParams['figure.figsize'] = (10.0, 6.0)
mpl.rc('axes', linewidth=2)
mpl.rcParams['lines.linewidth'] = 2
#mpl.rcParams['xtick.major.size'] = 20                                                                                                                                     
mpl.rcParams['xtick.major.width'] = 2
mpl.rcParams['ytick.major.width'] = 2

# Load Mesh
loc = '../'
mesh = Mesh(loc+'mesh.xml');
cd=MeshFunction('size_t',mesh,loc+'mesh_physical_region.xml')
fd=MeshFunction('size_t',mesh,loc+'mesh_facet_region.xml')

eta_list = np.linspace(-0.5, 0.5, 30)
Vapp_list = eta_list + Eo
i_list = []
idx = 0
c_0 = load_obj('initial_c0.pkl')
phi_z0 = -1
Vapp = Vapp_list[idx]

phi, c, trial_f, bcs_t = initial_solve(Vapp, phi_z0, c_0, mesh)
c_0 = c.copy()
phi_0 = phi.copy()

#print(c.vector()[:])
phi, err = Phi_solve(phi, c, trial_f, bcs_t)
#phi_0 = phi.copy()
print(err, '\n')


#for i in range(3):
while err>1e-12:
    c = NP_solve(phi, c, trial_f, bcs_t, Vapp)
    #phi_0 = phi.copy()
    #c_0 = c.copy()
    phi, err = Phi_solve(phi, c, trial_f, bcs_t)
    print(err, '\n')

c = NP_solve(phi, c, trial_f, bcs_t, Vapp)

save_obj(phi.vector()[:], 'phi_sol.pkl')
save_obj(c.vector()[:], 'c_sol.pkl')

_c_1, _c_2, _c_3, _c_4 = c.split()
current = get_current(_c_3, mesh, a_s = 1e-7)
i_list.append(current)
print(Vapp_list[idx], current)
