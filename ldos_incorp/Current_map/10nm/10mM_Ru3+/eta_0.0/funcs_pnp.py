from __future__ import print_function
from calendar import c
import numpy as np
import os, sys, re, glob
from dolfin import *
from fenics import *
from scipy.interpolate import CloughTocher2DInterpolator as CT
from mpi4py import MPI 
from numpy.random import rand
from numpy.linalg import norm
import pickle, yaml, mat73
from scipy.io import loadmat, savemat
from math import log10 , floor
from ase import Atoms

comm = MPI.COMM_WORLD

# Functions

def get_mesh_geo(mesh): 
    # Assuming mesh is circular at the bottom and top surfaces with origin at (0,0,0)
    #mesh = Mesh(mesh_file)
    bmesh = BoundaryMesh(mesh, "exterior", True)
    bm_coords = bmesh.coordinates()
    bot_coords = np.array([bm_coords[i] for i in range(len(bm_coords)) if np.any(np.absolute(bm_coords[i,2]) < DOLFIN_EPS)])
    a_s = np.max(norm(bot_coords[:,:2], axis=1))
    z = np.max(bm_coords[:,2])
    top_coords = np.array([bm_coords[i] for i in range(len(bm_coords)) if np.allclose(bm_coords[i,2],z)])
    a_t = np.max(norm(top_coords[:,:2], axis=1))
    return a_s, a_t, z

def BV_rates(ko, eta=0.0, A=1.0, alpha=0.5):
    kbT = 0.0259
    kred = A*ko*np.exp(-(alpha*(eta))/kbT)
    kox = A*ko*np.exp(((1-alpha)*(eta))/kbT)
    return kox, kred

def get_current_3D(func, mesh, D_o, a_s = 1e-7):
    # func: _c_3 in our case (Ru3+), a_s: 100nm pippete radius
    P1 = FiniteElement('P', tetrahedron, 1)
    F = 96485.332
    r_list = np.linspace(0,a_s,100)
    th_list = np.linspace(0,2*np.pi,100)
    Jo_mat = np.zeros((len(r_list), len(th_list)))
    vec = project(func.dx(2), FunctionSpace(mesh, P1))
    vec.set_allow_extrapolation(True)
    
    def xy_coord(r, th):
        return [r*np.cos(th), r*np.sin(th)]
    
    for j in range(len(r_list)):
        for k in range(len(th_list)):
            xy = xy_coord(r_list[j], th_list[k])
            point = np.array([xy[0], xy[1], 0.0])
            xy_point = Point(point)
            Jo = -D_o*vec(xy_point)
            Jo_mat[j, k] = Jo
    
    y = np.trapz(Jo_mat, x=th_list, axis=1)*r_list
    i = F*np.trapz(y, x=r_list)
    return i

def get_ko(pos, a_s=1e-7): # Testing rate model function
    # pos is at the bottom boundary
    return 0.0002*(1 + 0.001*(np.sin(4*np.pi*pos[0]/a_s) + np.cos(4*np.pi*pos[1]/a_s)))

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
    interpolator = CT(np.array([x,y]).T, z)
    return interpolator

def my_CT(xy, z):
    """CT interpolator + nearest-neighbor extrapolation.

    Parameters
    ----------
    xy : ndarray, shape (ndim, npoints)
        Coordinates of data points
    z : ndarray, shape (npoints)
        Values at data points

    Returns
    -------
    func : callable
        A callable object which mirrors the CT behavior,
        with an additional neareast-neighbor extrapolation
        outside of the data range.
    """
    x = xy[0, :]
    y = xy[1, :]
    f = CT(np.array([x,y]).T, z)

    # this inner function will be returned to a user
    def new_f(xx, yy):
        # evaluate the CT interpolator. Out-of-bounds values are nan.
        xx, yy = np.array(xx), np.array(yy)
        zz = f(xx, yy)
        nans = np.isnan(zz)

        if nans.any():
            # for each nan point, find its nearest neighbor
            inds = np.argmin((x[:, None] - xx[nans])**2 + (y[:, None] - yy[nans])**2, axis=0)
            # ... and use its value
            zz[nans] = z[inds]
        return zz

    return new_f

class Model_rates(UserExpression):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.kred = args[0] # kred value or interpolator
        self.kox = args[1] # kox value or interpolator
        self.Vdl = args[2] # Vdl value or interpolator
        # Moire cell parameters
        self.bvec1 = args[3] # basis vector 1
        self.bvec2 = args[4] # basis vector 2
        self.orig = args[5] # origin of the moire supercell
        self.center = args[6] # Center of nanopipette
        
    def eval(self, values, x):
        #kox, kred = BV_rates(get_ko(x, a_s=self.a_s), eta=self.Vapp) # Insert model
        if x[2] < 1e-9: # Keeping it close to the bottom surface
            pos_Ang = (1e10*np.array([x[0], x[1]])) + np.array(self.center)# Convert x from m to Ang before interpolation
            wrap_pt = wrap_pos(pos_Ang, self.bvec1, self.bvec2, orig=self.orig)[1] # Wrap position inside moire supercell
            values[0] = self.kred(wrap_pt[0], wrap_pt[1])
            values[1] = self.kox(wrap_pt[0], wrap_pt[1])
            values[2] = self.Vdl(wrap_pt[0], wrap_pt[1])
        else:
            values[0] = 0.0
            values[1] = 0.0
            values[2] = 0.0         
        
    def value_shape(self):
        return (3,)

def par_eval(mesh, vec, x, bbt = None, val_size=1):
    if bbt == None:
        bbt = mesh.bounding_box_tree()
    x_point = Point(x)
    bbt_id = bbt.compute_first_entity_collision(x_point)
    if bbt_id < mesh.num_cells():
        val = np.zeros(val_size)
        vec.eval_cell(val, x, Cell(mesh, bbt_id))
        return val
    else:
        return None
    
def tblg_dim(theta):
    alpha = 2.47
    sind = lambda degrees: np.sin(np.deg2rad(degrees))
    cosd = lambda degrees: np.cos(np.deg2rad(degrees))
    sc_alpha = alpha/(2*sind(theta/2))
    sc_area = (sc_alpha**2)*sind(60)*1e-2
    return sc_alpha, sc_area

def save_sol(fname, u, mesh):
    # u is a function                   
    output_file = HDF5File(mesh.mpi_comm(), fname+".h5", "w")
    output_file.write(u, "solution")
    output_file.close()
    
def load_sol(fname, U, mesh):
    # U is a function                                
    input_file = HDF5File(mesh.mpi_comm(), fname, "r")
    input_file.read(U, "solution")
    input_file.close()
    return U

def save_obj(obj, name ):
    with open(name, 'wb') as f:
        pickle.dump(obj, f)

def load_obj(name ):
    with open(name, 'rb') as f:
        return pickle.load(f)

def read_yml(ymlfile):
    with open(ymlfile) as file:
        out_dict = yaml.load(file, Loader=yaml.FullLoader)
    return out_dict

def write_yml(dict_file, fname='test.yaml'):
    with open(fname, 'w') as file:
        yaml.dump(dict_file, file)

def round_it(x, sig):
    return round(x, sig-int(floor(log10(abs(x))))-1)