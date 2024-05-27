import numpy as np
import os, sys, re, glob
from dolfin import *
from fenics import *
import pickle, mat73
from funcs_pnp import *
from numpy.random import rand
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib as mpl
from scipy.interpolate import CloughTocher2DInterpolator as CT
from scipy.io import loadmat, savemat 
from numpy.linalg import norm

# Load and scale mesh
mesh = Mesh('/Users/mbabar/Desktop/PhD/Analysis/PDE/PNP_solve/3D_vmg/Mesh/Testing/mesh_2.xml')
m_coords = mesh.coordinates()
scal_r = np.float64(1e-9) # scale to nm
m_coords[:, :] *= scal_r
mesh.bounding_box_tree().build(mesh) # development version

a_s, a_t, z_t = get_mesh_geo(mesh)
print('a_s = {}, a_t = {}, z_t = {}'.format(a_s, a_t, z_t))

# Load solution
fname = sys.argv[1]
P1 = FiniteElement('P', tetrahedron, 1)
vec = Function(FunctionSpace(mesh, P1))
load_sol(fname+'.h5', vec, mesh)
vec.set_allow_extrapolation(True)
#_phi, _c_1, _c_2, _c_3, _c_4 = split(vec)

dos_file = '/Users/mbabar/Desktop/PhD/Analysis/PDE/PNP_solve/3D_vmg/ldos_incorp/ldos-th_2_nsamps_30_knum_24.mat'
rate_loc = '/Users/mbabar/Desktop/PhD/Analysis/PDE/PNP_solve/3D_vmg/ldos_incorp/k_data_Vapp_0.0.mat'
ldos_data = loadmat(dos_file)
rate_data = mat73.loadmat(rate_loc)
rscx = ldos_data["rscx"]
rscy = ldos_data["rscy"]
kred_data, kox_data = rate_data["kred_data"], rate_data["kox_data"]
Vdl_data = rate_data["Vdl_data"]
theta = np.squeeze(ldos_data["theta"])
nsamps = ldos_data["nsamps"]
Vapp = rate_data["Vapp"]
print('Vapp={}'.format(Vapp),'\n')

moir_l, sc_area = tblg_dim(theta) # Moire supercell length and area Ang
moir_a = np.pi/3 # Moire angle 60deg
orig = moir_l*np.array([rscx[0,0], rscy[0,0]]) # Origin of the moire supercell
bvec1 = np.array([moir_l*np.cos(moir_a/2), moir_l*np.sin(moir_a/2)]) # basis vector 1
bvec2 = np.array([moir_l*np.cos(moir_a/2), -moir_l*np.sin(moir_a/2)]) # basis vector 2

theta_list = np.linspace( 0 , 2 * np.pi , 150 ) 
radius = 1e10*a_s
a = radius * np.cos(theta_list)
b = radius * np.sin(theta_list)

size = 100
xc, yc = np.linspace(-radius, radius,100), np.linspace(-radius, radius,100)
X, Y = np.meshgrid(xc, yc)

Z2 = np.zeros(np.shape(X.flatten()))
for i in range(len(Z2)):
    # wrap_pt = wrap_pos([X.flatten()[i], Y.flatten()[i]], bvec1, bvec2, orig=orig)[1]
    # p = 1e-10*np.array([wrap_pt[0], wrap_pt[1], 0.0]) # Convert to meter
    p = 1e-10*np.array([X.flatten()[i], Y.flatten()[i], 0.0])
    #print(p)
    Z2[i] = vec(*Point(p)) #ip_red(wrap_pt[0], wrap_pt[1])[1]
        
Z2 = Z2.reshape((size, size))
plt.contourf(X, Y, Z2)
plt.plot(a, b, 'k--')
plt.colorbar()
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.axis('square')
#plt.xlim([-radius, radius])
#plt.ylim([-radius, radius])
plt.xlabel('x (A)')
plt.ylabel('y (A)')
plt.tight_layout()
plt.savefig(fname+'.png')