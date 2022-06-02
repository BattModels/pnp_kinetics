from __future__ import print_function
from calendar import c
import numpy as np
import os, sys, re, glob
from dolfin import *
from fenics import *
from numpy.random import rand
from numpy.linalg import norm
import pickle

# Load Mesh
loc = '../'
mesh = Mesh(loc+'mesh.xml');
cd=MeshFunction('size_t',mesh,loc+'mesh_physical_region.xml')
fd=MeshFunction('size_t',mesh,loc+'mesh_facet_region.xml')

# Constants Radial
z_1, z_2, z_3, z_4 = 1, -1, 3, 2 # Charge numbers K+, Cl-, Ru3+, Ru2+
eps = 80
eps0 = 8.854187e-12 # Vacuum permittivity F/m
e = 1.60218e-19 # Coulombs 
kbT = 0.0259 # eV, room temp
F = 96485.332 # C/mol
Eo = -0.07
D_o, D_r = 6.5e-10, 6.5e-10 # Diff constants (m^2/s)
a = (1/kbT)
b = a*(F/(eps*eps0))
#a, b = 1, 10e14

# Boundaries
tol = 1E-14
def boundary_1(x, on_boundary):
    return on_boundary and near(x[1], 0, tol)

def boundary_2(x, on_boundary):
    return on_boundary and near(x[0], 1e-7, tol)# and (x[1] < 1e-7))

def boundary_3(x, on_boundary):
    return on_boundary and (1e-7 < x[0] < 2e-7)

def boundary_4(x, on_boundary):
    return on_boundary and (2e-7 < x[0] < 3.5e-6)

def boundary_5(x, on_boundary):
    return on_boundary and near(x[1], 2e-5, tol)

def save_obj(obj, name ):
    with open(name, 'wb') as f:
        pickle.dump(obj, f)

def load_obj(name ):
    with open(name, 'rb') as f:
        return pickle.load(f)


def BV_rates(ko, Eo, A=1.0, Vapp=0.0, alpha=0.5):
    kbT = 0.0259
    kox = A*ko*np.exp(-(alpha*(Vapp - Eo))/kbT)
    kred = A*ko*np.exp(-((1-alpha)*(Vapp - Eo))/kbT)
    return kox, kred


def get_current(func, mesh, a_s = 1e-7):
    # func: _c_3 in our case (Ru3+), a_s: 100nm pippete radius
    P1 = FiniteElement('P', triangle, 1)
    F = 96485.332
    r = np.linspace(0,a_s,1000)
    Jo_list = []
    vec = project(func.dx(1), FunctionSpace(mesh, P1))
    for j in range(len(r)):
        point = np.array([r[j], 0.0])
        x = np.array(point)
        x_point = Point(*x) 
        Jo = -D_o*vec(x_point)
        Jo_list.append(Jo)
    i = 2*np.pi*F*np.trapz(Jo_list*r, x=r)
    return i


def initial_solve(Vapp, phi_z0, c_0, mesh):

    ## initialize phi
    V_phi = FunctionSpace(mesh, 'P', 1)
    P1 = FiniteElement('P', triangle, 1)
    element = MixedElement([P1, P1, P1, P1])
    V_c = FunctionSpace(mesh, element)

    phi_D = Expression(str(0.0), degree=1)
    phi_L = Expression(str(phi_z0), degree=1)

    phi = TrialFunction(V_phi)
    v = TestFunction(V_phi)


    ## initialize conc

    c_ruhex = 2 #mM
    c_KCl = 100 #mM
    # Note 1mM = 1mol/m3 (SI units)
    c_a = Expression(str(c_KCl), degree=1)
    c_b = Expression(str(z_1*c_KCl + z_3*c_ruhex), degree=1)
    c_c = Expression(str(c_ruhex), degree=1)
    c_d = Expression(str(0.0), degree=1)

    q_1, q_2, q_3, q_4 = TestFunctions(V_c)
    c = Function(V_c)
    c_1, c_2, c_3, c_4 = split(c)

    # phi BCs
    bc_1= DirichletBC(V_phi, phi_D, boundary_5)
    bc_2= DirichletBC(V_phi, phi_L, boundary_1)
    bcs_phi = [bc_1, bc_2] # All BCs are Neumann

    # conc BCs
    bc_5a= DirichletBC(V_c.sub(0), c_a, boundary_5)
    bc_5b= DirichletBC(V_c.sub(1), c_b, boundary_5)
    bc_5c= DirichletBC(V_c.sub(2), c_c, boundary_5)
    bc_5d= DirichletBC(V_c.sub(3), c_d, boundary_5)
    bcs_c = [bc_5a, bc_5b, bc_5c, bc_5d]

    ## Different Neumann BC for bottom BC
    # create a mesh function which assigns an unsigned integer (size_t) to each edge
    mf = MeshFunction("size_t", mesh, 1) # 3rd argument is dimension of an edge
    mf.set_all(0) # initialize the function to zero
    class BottomBoundary(SubDomain):
        def inside(self, x, on_boundary):
            return near(x[1], 0.0, tol) and on_boundary

    bottomboundary = BottomBoundary() # instantiate it

    # use this bottomboundary object to set values of the mesh function to 1 in the subdomain
    bottomboundary.mark(mf, 1)

    # define a new measure ds based on this mesh function
    ds = Measure("ds", domain=mesh, subdomain_data=mf)
    #ds = Measure("ds")(subdomain_data=mf)

    ### Compute solution
    ## Define problem for phi
    r = Expression('x[0]', degree=1)
    nabla_phi = (z_1*c_1 + z_2*c_2 + z_3*c_3 + z_4*c_4)
    g = Expression('0.0', degree=1)
    F_phi = (dot(grad(phi),grad(v)))*r*dx() 
    #L_phi = 0
    L_phi = (F/(eps*eps0))*nabla_phi*v*r*dx() + g*v*r*ds()
    phi = Function(V_phi)
    solve(F_phi == L_phi, phi, bcs_phi)

    ## Define problem for conc
    r = Expression('x[0]', degree=1)
    #f = Constant(0.0)
    #g = Expression('-4*x[1]', degree=1)
    g_1 = Expression('0.0', degree=1)
    g_2 = Expression('0.0', degree=1)
    g_3 = Expression('0.0', degree=1)
    g_4 = Expression('0.0', degree=1)
    n = FacetNormal(mesh)
    #m1 = dot(grad(c_3), n)
    #m1 = Dx(c_3,1) does not give correct solution
    kox, kred = BV_rates(1e-3, -0.07, A=1.0, Vapp=Vapp)
    m1 = -(kred*c_3 - kox*c_4)/D_o # Rate theory input
    m2 = dot(grad(c_4), n)
    F_c = ((dot(grad(c_1), grad(q_1))) - ((z_1*a)*div(c_1*grad(phi))*q_1))*r*dx() \
        + ((dot(grad(c_2), grad(q_2))) - ((z_2*a)*div(c_2*grad(phi))*q_2))*r*dx() \
        + ((dot(grad(c_3), grad(q_3))) - ((z_3*a)*div(c_3*grad(phi))*q_3))*r*dx() \
        + ((dot(grad(c_4), grad(q_4))) - ((z_4*a)*div(c_4*grad(phi))*q_4))*r*dx() \
        - g_1*q_1*r*ds() - g_2*q_2*r*ds() - g_3*q_3*r*ds(0) - g_4*q_4*r*ds(0) \
        - m1*q_3*r*ds(1) + (D_o/D_r)*m1*q_4*r*ds(1)

    L_c = 0
    c.vector()[:] = c_0
    solve(F_c == L_c, c, bcs_c)
    trial_f = [v, q_1, q_2, q_3, q_4]
    bcs_t = [bcs_phi, bcs_c]
    #print(phi.vector()[:].copy())
    return phi, c, trial_f, bcs_t

## Using outputs of initial_solve
def Phi_solve(phi, c, trial_f, bcs_t):
    phi_0 = phi.copy()
    phi_vec = phi_0.vector()[:].copy()
    c_1, c_2, c_3, c_4 = split(c)
    bcs_phi, bcs_c = bcs_t[0], bcs_t[1]
    v = trial_f[0]

    nabla_phi = (z_1*c_1 + z_2*c_2 + z_3*c_3 + z_4*c_4)
    r = Expression('x[0]', degree=1)
    g = Expression('0.0', degree=1)
    F_phi = (dot(grad(phi),grad(v)))*r*dx() - (F/(eps*eps0))*nabla_phi*v*r*dx() - g*v*r*ds()
    L_phi = 0
    phi.vector()[:] = phi_0.vector()[:]
    solve(F_phi == L_phi, phi, bcs_phi)
    #print(errornorm(phi_0, phi, 'L2'))
    #print(phi.vector()[:].copy())
    #print(phi_vec)
    err = norm(np.array(phi.vector()[:].copy() - phi_vec))
    return phi, err

def NP_solve(phi, c, trial_f, bcs_t, Vapp):
    c_0 = c.copy()
    c_1, c_2, c_3, c_4 = c.split()
    q_1, q_2, q_3, q_4 = trial_f[1], trial_f[2], trial_f[3], trial_f[4]
    bcs_phi, bcs_c = bcs_t[0], bcs_t[1]
    nabla_phi = (z_1*c_1 + z_2*c_2 + z_3*c_3 + z_4*c_4)

    r = Expression('x[0]', degree=1)
    g_1 = Expression('0.0', degree=1)
    g_2 = Expression('0.0', degree=1)
    g_3 = Expression('0.0', degree=1)
    g_4 = Expression('0.0', degree=1)
    n = FacetNormal(mesh)
    #m1 = dot(grad(c_3), n)
    #m1 = Dx(c_3,1) does not give correct solution
    kox, kred = BV_rates(1e-3, Eo, A=1.0, Vapp=Vapp)
    m1 = -(kred*c_3 - kox*c_4)/D_o # Rate theory input
    #m2 = dot(grad(c_4), n)
    F_c = ((dot(grad(c_1), grad(q_1))) - ((z_1*a)*div(c_1*grad(phi))*q_1))*r*dx() \
        + ((dot(grad(c_2), grad(q_2))) - ((z_2*a)*div(c_2*grad(phi))*q_2))*r*dx() \
        + ((dot(grad(c_3), grad(q_3))) - ((z_3*a)*div(c_3*grad(phi))*q_3))*r*dx() \
        + ((dot(grad(c_4), grad(q_4))) - ((z_4*a)*div(c_4*grad(phi))*q_4))*r*dx() \
        - g_1*q_1*r*ds() - g_2*q_2*r*ds() - g_3*q_3*r*ds(0) - g_4*q_4*r*ds(0) \
        - m1*q_3*r*ds(1) + (D_o/D_r)*m1*q_4*r*ds(1)

    L_c = 0
    c.vector()[:] = c_0.vector()[:]
    solve(F_c == L_c, c, bcs_c)
    return c


def iter_solve(phi, c, trial_f, bcs_t, Vapp):
    phi_0 = phi.copy()
    phi_vec = phi_0.vector()[:].copy()
    c_1, c_2, c_3, c_4 = split(c)
    bcs_phi, bcs_c = bcs_t[0], bcs_t[1]
    v, q_1, q_2, q_3, q_4 = trial_f[0], trial_f[1], trial_f[2], trial_f[3], trial_f[4]

    nabla_phi = (z_1*c_1 + z_2*c_2 + z_3*c_3 + z_4*c_4)
    r = Expression('x[0]', degree=1)
    g = Expression('0.0', degree=1)
    F_phi = (dot(grad(phi),grad(v)))*r*dx() - (F/(eps*eps0))*nabla_phi*v*r*dx() - g*v*r*ds()
    L_phi = 0

    #phi.vector()[:] = phi_0.vector()[:]
    solve(F_phi == L_phi, phi, bcs_phi)
    #print(errornorm(phi_0, phi, 'L2'))
    #print(phi.vector()[:].copy())
    #print(phi_vec)
    err = norm(np.array(phi.vector()[:] - phi_vec))
    print('phi_err', err, '\n')

    # Solve c
    c_0 = c.copy()
    c_vec = c_0.vector()[:].copy()
    g_1 = Expression('0.0', degree=1)
    g_2 = Expression('0.0', degree=1)
    g_3 = Expression('0.0', degree=1)
    g_4 = Expression('0.0', degree=1)
    n = FacetNormal(mesh)
    #m1 = dot(grad(c_3), n)
    #m1 = Dx(c_3,1) does not give correct solution
    kox, kred = BV_rates(1e-3, Eo, A=1.0, Vapp=Vapp)
    m1 = -(kred*c_3 - kox*c_4)/D_o # Rate theory input
    #m2 = dot(grad(c_4), n)
    F_c = ((dot(grad(c_1), grad(q_1))) - ((z_1*a)*div(c_1*grad(phi))*q_1))*r*dx() \
        + ((dot(grad(c_2), grad(q_2))) - ((z_2*a)*div(c_2*grad(phi))*q_2))*r*dx() \
        + ((dot(grad(c_3), grad(q_3))) - ((z_3*a)*div(c_3*grad(phi))*q_3))*r*dx() \
        + ((dot(grad(c_4), grad(q_4))) - ((z_4*a)*div(c_4*grad(phi))*q_4))*r*dx() \
        - g_1*q_1*r*ds() - g_2*q_2*r*ds() - g_3*q_3*r*ds(0) - g_4*q_4*r*ds(0) \
        - m1*q_3*r*ds(1) + (D_o/D_r)*m1*q_4*r*ds(1)

    L_c = 0
    c.vector()[:] = c_0.vector()[:]
    solve(F_c == L_c, c, bcs_c)
    err = norm(np.array(c.vector()[:] - c_vec))
    print('c_err', err, '\n')
    return phi, c