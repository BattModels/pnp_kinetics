import numpy as np
from scipy.io import loadmat, savemat
import os, re, pickle, sys
from scipy.interpolate import CloughTocher2DInterpolator
from numpy.linalg import norm
import matplotlib.pyplot as plt
from ase import Atoms
import mat73

## Load data
dos_loc = 'ldos-th_2_nsamps_30_knum_24.mat' # Location of ldos file 
rate_loc = 'k_data_η_0.0.mat'
ldos_data = loadmat(dos_loc)
k_data = mat73.loadmat(rate_loc)
# ldos = ldos_data["data"]
rscx = ldos_data["rscx"]
rscy = ldos_data["rscy"]
kox_data = k_data["kox_data"]
kred_data = k_data["kred_data"]
theta = np.squeeze(ldos_data["theta"])
nsamps = ldos_data["nsamps"]

# Calculate constants
alpha = 2.47
moir_a = np.pi/3 # Moire angle 60deg
sc_alpha = alpha/(2*np.sin(np.deg2rad(theta/2)))
sc_area = (sc_alpha**2)*np.sin(moir_a)*1e-2 #area in nm^2
moir_l = sc_alpha # Moire supercell length Ang
#ldos[:, idy, idx]
#ldos[5499, idy, idx]
orig = sc_alpha*np.array([rscx[0,0], rscy[0,0]]) # Origin of the moire supercell
bvec1 = np.array([moir_l*np.cos(moir_a/2), moir_l*np.sin(moir_a/2)]) # basis vector 1
bvec2 = np.array([moir_l*np.cos(moir_a/2), -moir_l*np.sin(moir_a/2)]) # basis vector 2

## User functions
def query_point(point, rscx, rscy, moir_l):
    # Query indices closest to point in rscx and rscy
    coords = moir_l*np.array([rscx.flatten(), rscy.flatten()])
    length = len(rscx)
    diff = coords - np.reshape(np.array(point), (2,1))
    id = np.argmin(norm(diff, axis=0))
    idx, idy = id % length, id // length
    return idx, idy

def wrap_pos(point, bvec1, bvec2, orig=[0,0]):
    atoms = Atoms('CC', positions=[(orig[0], orig[1], 0), (point[0], point[1], 0)])
    cell = np.array([(bvec1[0], bvec1[1], 0), (bvec2[0], bvec2[1], 0), (0, 0, 0)])
    cell[:2,:] = cell[:2,:] + np.array([orig[0], orig[1], 0])
    atoms.set_cell(cell)
    atoms.set_pbc([1,1,0])
    atoms.wrap()
    return atoms.get_positions()[0,:2], atoms.get_positions()[1,:2]

def interpolator_2d(x, y, z):
    # Input x, y and z (all same lengths and flattened)
    interpolator = CloughTocher2DInterpolator(np.array([x,y]).T, z)
    return interpolator

# 2D Interpolation
coords = moir_l*np.array([rscx.flatten(), rscy.flatten()])
x, y = coords[0,:], coords[1,:]
z = kred_data.flatten()

ip = interpolator_2d(x, y, z)

point = [31, -101] # xy coordinates of query in Ang
print("Coordinates:",point)
wrap_point = wrap_pos(point, bvec1, bvec2, orig=orig)[1]
print("Wrap point:",wrap_point)
idx, idy = query_point(wrap_point, rscx, rscy, moir_l)
wrap_test = [moir_l*rscx[idx, idy], moir_l*rscy[idy, idx]]
print("Wrap closest point to rsc", wrap_test)
print("Interpolation", ip(wrap_point[0], wrap_point[1]))
print("Interpolation to rsc", ip(wrap_test[0], wrap_test[1]))
print("Exact value", z)
