from __future__ import print_function
from calendar import c
import numpy as np
import os, sys, re, glob
from dolfin import *
from fenics import *
from mpi4py import MPI 
from numpy.random import rand
from numpy.linalg import norm
import pickle
from scipy.io import loadmat, savemat
from math import log10 , floor
import yaml

comm = MPI.COMM_WORLD

# Functions

def get_mesh_geo(mesh): 
    # Assuming mesh is circular at the bottom and top surfaces with origin at (0,0,0)
    #mesh = Mesh(mesh_file)
    bmesh = BoundaryMesh(mesh, "exterior", True)
    bm_coords = bmesh.coordinates()
    bot_coords = np.array([bm_coords[i] for i in range(len(bm_coords)) if np.allclose(bm_coords[i,2],0)])
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

def get_current_2D(func, mesh, D_o, a_s = 1e-7):
    # func: _c_3 in our case (Ru3+), a_s: 100nm pippete radius
    P1 = FiniteElement('P', tetrahedron, 1)
    F = 96485.332
    r = np.linspace(0,a_s,100)
    Jo_list = []
    vec = project(func.dx(2), FunctionSpace(mesh, P1))
    for j in range(len(r)):
        point = np.array([r[j], 0.0])
        x = np.array(point)
        x_point = Point(x)
        Jo = -D_o*vec(x_point)
        Jo_list.append(Jo)
    i = 2*np.pi*F*np.trapz(Jo_list*r, x=r)
    return i

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

def get_ko(pos, a_s=1e-7):
    # pos is at the bottom boundary
    return 0.0002*(1 + 0.001*(np.sin(4*np.pi*pos[0]/a_s) + np.cos(4*np.pi*pos[1]/a_s)))
    
class Model_rates(UserExpression):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.kred = args[0] # kred from config_file
        self.kox = args[1] # kox
        #self.Vapp = args[2] # Vapp
        #self.a_s = args[3] # Nanopipette radius
        
    def eval(self, values, x):
        #kox, kred = BV_rates(get_ko(x, a_s=self.a_s), eta=self.Vapp) # Insert model
        values[0] = self.kred
        values[1] = self.kox
        
    def value_shape(self):
        return (2,)

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

def check_sol1():
    point = [a_s, 0.5e-7]
    x = np.array(point)
    x_point = Point(*x) 
    P1 = FiniteElement('P', triangle, 1)
    vec1 = project(_c_4.dx(0), FunctionSpace(mesh, P1))
    vec2 = project(_c_3.dx(0), FunctionSpace(mesh, P1))
    vec3 = project(_c_2.dx(0), FunctionSpace(mesh, P1))
    vec4 = project(_c_1.dx(0), FunctionSpace(mesh, P1))
    print(norm(np.array([vec1(x_point), vec2(x_point), vec3(x_point), vec4(x_point)])))

def check_sol2():
    point = [0.5e-7, 0]
    x = np.array(point)
    x_point = Point(*x) 
    P1 = FiniteElement('P', triangle, 1)
    vec1 = project(_c_2.dx(1), FunctionSpace(mesh, P1))
    vec2 = project(_c_1.dx(1), FunctionSpace(mesh, P1))
    print(norm(np.array([vec1(x_point), vec2(x_point)])))
    
def tblg_dim(theta):
    alpha = 2.47
    sind = lambda degrees: np.sin(np.deg2rad(degrees))
    cosd = lambda degrees: np.cos(np.deg2rad(degrees))
    sc_alpha = alpha/(2*sind(theta/2))
    sc_area = sc_alpha^2*sind(60)*1e-2
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