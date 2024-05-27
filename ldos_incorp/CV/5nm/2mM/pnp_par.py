import numpy as np
import os, sys, re, glob
from dolfin import *
from fenics import *
import pickle
from funcs_pnp import *
from numpy.random import rand
from numpy.linalg import norm
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib as mpl
from scipy.io import loadmat, savemat
from mpi4py import MPI
import yaml, mat73

comm = MPI.COMM_WORLD
rank = comm.rank

config_loc = sys.argv[1] # Location of config file
#center_id = int(sys.argv[2]) # Center of nanopipette radius

# Load configuration and rate model data
with open(config_loc, 'r') as file:
    constants = yaml.safe_load(file)
    globals().update(constants)

# Load Mesh and calculate geometry
mesh = Mesh(mesh_file)

# Scale mesh
m_coords = mesh.coordinates()
scale_factor = np.float64(1e-9)
m_coords[:, :] *= scale_factor
mesh.bounding_box_tree().build(mesh) # development version

bmesh = BoundaryMesh(mesh, "exterior", True)
bm_list = []
bm_coords = bmesh.coordinates()

data_bm = np.array(comm.gather(bm_coords, root=0), dtype=object)
if rank==0:
    for j in range(len(data_bm)):
        bm_list += list(data_bm[j])
    bm_coords = np.array(bm_list)
    bot_coords = np.array([bm_coords[i] for i in range(len(bm_coords)) if np.any(np.absolute(bm_coords[i,2]) < DOLFIN_EPS)])
    a_s = np.max(norm(bot_coords[:,:2], axis=1))
    z_t = np.max(bm_coords[:,2])
    top_coords = np.array([bm_coords[i] for i in range(len(bm_coords)) if np.allclose(bm_coords[i,2],z_t)])
    a_t = np.max(norm(top_coords[:,:2], axis=1))
    geo_mesh = {'a_s':round_it(a_s, 3), 'a_t':round_it(a_t, 3), 'z_t':round_it(z_t, 3)}
else:
    geo_mesh = None

geo_mesh = comm.bcast(geo_mesh, root=0)
a_s, a_t, z_t = geo_mesh['a_s'], geo_mesh['a_t'], geo_mesh['z_t']
if rank==0:
    print("Mesh file: ",mesh_file)
    print("Mesh geometry: a_s={}, a_t={}, z_t={}".format(a_s, a_t, z_t))
#a_s, a_t, z_t = 1e-7, 3.5e-6, 2e-5

## initialize phi                                                                                                                       
# V_phi = FunctionSpace(mesh, 'P', 1)

deg = 1
phi_D = Expression(str(0.0), degree=deg)

def boundary(x, on_boundary):
    return on_boundary

# phi = TrialFunction(V_phi)
# v = TestFunction(V_phi)

## initialize conc                                                                                                                      
P1 = FiniteElement('P', tetrahedron, 1)
element = MixedElement([P1, P1, P1, P1, P1])
V = FunctionSpace(mesh, element)

c_a = Expression(str(c_KCl), degree=deg)
c_b = Expression(str(z_1*c_KCl + z_3*c_ruhex), degree=deg)
c_c = Expression(str(c_ruhex), degree=deg)
c_d = Expression(str(0.0), degree=deg)

v, q_1, q_2, q_3, q_4 = TestFunctions(V)
sol = Function(V)
dsol = TrialFunction(V)
phi, c_1, c_2, c_3, c_4 = split(sol)
#c_ini = load_sol('Sol/_c',c,mesh)
#c.vector()[:] = c_ini.vector()[:]
#sol.vector().set_local(rand(sol.vector().size())) # For starting calculations
#sol.vector().apply("")

## BCs
tol = DOLFIN_EPS
def boundary_5(x, on_boundary):
    return on_boundary and near(x[2], z_t, tol)

# phi BCs                                                                                                                               
bc_1= DirichletBC(V.sub(0), phi_D, boundary_5)                                                                                                                      

# conc BCs                                                                                                                              
bc_5a= DirichletBC(V.sub(1), c_a, boundary_5)
bc_5b= DirichletBC(V.sub(2), c_b, boundary_5)
bc_5c= DirichletBC(V.sub(3), c_c, boundary_5)
bc_5d= DirichletBC(V.sub(4), c_d, boundary_5)
bcs = [bc_1, bc_5a, bc_5b, bc_5c, bc_5d]

## Different Neumann BC for bottom BC                                                                                                   
## create a mesh function which assigns an unsigned integer (size_t) to each edge                                                        
mf = MeshFunction("size_t", mesh, mesh.topology().dim() - 1) # 3rd argument is dimension of an edge                                                             
mf.set_all(0) # initialize the function to zero                                                                                         
class BottomBoundary(SubDomain):
    def inside(self, x, on_boundary):
        tol = DOLFIN_EPS
        return near(x[2], 0.0, tol) and on_boundary # Use x[2] for 3D.

bottomboundary = BottomBoundary() # instantiate it 

# use this bottomboundary object to set values of the mesh function to 1 in the subdomain                                               
bottomboundary.mark(mf, 1)

## Define a new measure ds based on this mesh function                                                                                   
ds = Measure("ds", domain=mesh, subdomain_data=mf)
#File("mf.pvd").write(mf) # Save mf for Paraview                                                                                        

### Compute solution                                                                                                                    
## Define problem for phi                                                                                                               
n = FacetNormal(mesh)
#r = Expression('x[0]', degree=deg)
#nabla_phi = (z_1*c_1 + z_2*c_2 + z_3*c_3 + z_4*c_4)                                                                                    
dphi_1 = Expression('0.0', degree=deg)

## Define problem for conc                                                                                                                                                                                                                 
g_1 = Expression('0.0', degree=deg)
g_2 = Expression('0.0', degree=deg)
g_3 = Expression('0.0', degree=deg)
g_4 = Expression('0.0', degree=deg)
#m1 = dot(grad(c_3), n)                                                                                                                 
#m1 = Dx(c_3,1) does not give correct solution  
# kox, kred = BV_rates(k_o, eta=Vapp)
# m1 = -(kred*c_3 - kox*c_4)/D_o # Rate theory input 
# m2 = dot(grad(c_4), n)

## Input kinetic rates
ldos_data = loadmat(constants["dos_file"])
rate_data = mat73.loadmat(constants["rate_loc"])
rscx = ldos_data["rscx"]
rscy = ldos_data["rscy"]
if type(MHC_prefactor) == str:
    pref_data= loadmat(MHC_prefactor)
    A = np.squeeze(pref_data['pref_sc'])
else:
    A = float(MHC_prefactor)

kred_data, kox_data = A*rate_data["kred_data"], A*rate_data["kox_data"]
Vdl_data = rate_data["Vdl_data"]
theta = np.squeeze(ldos_data["theta"])
nsamps = ldos_data["nsamps"]
Vapp = np.round(np.double(rate_data["Vapp"]),3)
if rank==0:
    print('Vapp={}'.format(Vapp),'\n')

moir_l, sc_area = tblg_dim(theta)
moir_a = np.pi/3
orig = moir_l*np.array([rscx[0,0], rscy[0,0]]) # Set origin in Angstrom
# bvec1 = np.array([np.amax(moir_l*rscx)/2, np.amax(moir_l*rscy)])
# bvec2 = np.array([np.amax(moir_l*rscx)/2, -np.amax(moir_l*rscy)])
bvec1 = np.array([moir_l*np.cos(moir_a/2), moir_l*np.sin(moir_a/2)]) # basis vector 1
bvec2 = np.array([moir_l*np.cos(moir_a/2), -moir_l*np.sin(moir_a/2)]) # basis vector 2
coords = moir_l*np.array([rscx.flatten(), rscy.flatten()])
z_red, z_ox = kred_data.flatten(), kox_data.flatten()
ip_red, ip_ox = my_CT(coords, z_red), my_CT(coords, z_ox) # 2D interpolators of kred and kox
ip_Vdl = my_CT(coords, Vdl_data.flatten()) # Interpolators take coords in Angstroms

center_id, center = center_domain(coords, bvec1, domain)
# if rank==0:
#     print('\nCenter = {}'.format(center))
V_k = FunctionSpace(mesh, MixedElement([P1, P1, P1]))
k_model = Function(V_k)
u_bc = Model_rates(ip_red, ip_ox, ip_Vdl, bvec1, bvec2, orig, center)
k_model = interpolate(u_bc, V_k)
kred_f, kox_f, Vdl_f = k_model[0], k_model[1], k_model[2]
# vec_red = project(kred_f, FunctionSpace(mesh, P1))
# vec_ox = project(kox_f, FunctionSpace(mesh, P1))
# vec_Vdl = project(Vdl_f, FunctionSpace(mesh, P1))
# vec_red.set_allow_extrapolation(True)
# vec_ox.set_allow_extrapolation(True)
# vec_Vdl.set_allow_extrapolation(True)
# save_sol('kred_Vapp_'+str(Vapp), vec_red, mesh)
# save_sol('kox_Vapp_'+str(Vapp), vec_ox, mesh)
# save_sol('Vdl_Vapp_'+str(Vapp), vec_Vdl, mesh)
# File("kred_f.pvd").write(vec_red)
# File("kox_f.pvd").write(vec_ox)
# File("Vdl_f.pvd").write(vec_Vdl)
if rank==0:
    print("Interpolated kred, kox, Vdl")

## Solve 3D-PNP
Func = (inner(grad(phi),grad(v)))*dx() - ((Vdl_f - phi)/d_h)*v*ds(1) \
    - (F/(eps*eps0))*(z_1*c_1 + z_2*c_2 + z_3*c_3 + z_4*c_4)*v*dx() + dphi_1*v*ds(0) \
    + ((inner(grad(c_1), grad(q_1))) - ((z_1/kbT)*div(c_1*grad(phi))*q_1))*dx() \
    + ((inner(grad(c_2), grad(q_2))) - ((z_2/kbT)*div(c_2*grad(phi))*q_2))*dx() \
    + ((inner(grad(c_3), grad(q_3))) - ((z_3/kbT)*div(c_3*grad(phi))*q_3))*dx() \
    + ((inner(grad(c_4), grad(q_4))) - ((z_4/kbT)*div(c_4*grad(phi))*q_4))*dx() \
    - g_1*q_1*ds() - g_2*q_2*ds() - g_3*q_3*ds(0) - g_4*q_4*ds(0) \
    - (-(kred_f*c_3 - kox_f*c_4)/D_o)*q_3*ds(1) \
    + (-(kred_f*c_3 - kox_f*c_4)/D_r)*q_4*ds(1) \
    #- m1*q_3*ds(1) + (D_o/D_r)*m1*q_4*ds(1) \

#a = lhs(Func)
#L = rhs(Func)
L = 0

#c.vector()[:] = csol.vector()[:]
#c.vector().set_local(rand(c.vector().size())) # For starting calculations
#c.vector().apply("")

# Solver configs                                                                                                                        
ffc_options = {"optimize": True, "quadrature_degree": 6}

solver_parameters = {"nonlinear_solver": "snes",
                    "snes_solver" : {
                        "report" : True,
                        "error_on_nonconvergence" : False,
                        "line_search" : "bt",
                        "linear_solver" : "bicgstab", # cg, bicgstab, gmres, lu, superlu_dist                                          
                        "preconditioner" : "sor", # sor, amg, hypre_euclid, petsc_amg, ilu, icc                                                                   
                        "absolute_tolerance" : 1E-7,
                        "relative_tolerance" : 1E-6,                                                                                  
                        "maximum_iterations" : 20000,                                                                                 
                        "krylov_solver" : {
                            "maximum_iterations" : 40000,
                            #"monitor_convergence" : True,                                                                             
                            "nonzero_initial_guess" : False
                        }}}

J = derivative(Func, sol, dsol)
problem = NonlinearVariationalProblem(Func, sol, bcs, J, form_compiler_parameters=ffc_options)
solver = NonlinearVariationalSolver(problem)
solver.parameters.update(solver_parameters)
solver.solve()
_phi, _c_1, _c_2, _c_3, _c_4 = sol.split()

## Save solution
sol_0 = sol.copy()
save_sol(sol_file, sol, mesh) # Use solution to get_current.py

## Check solution 
# point = [0.0e-7, 0.0e-7, 0]
# p = np.array(point)

# bbt = mesh.bounding_box_tree()
# d_c4 = project(_c_4.dx(2), FunctionSpace(mesh, P1))
# d_c3 = project(_c_3.dx(2), FunctionSpace(mesh, P1))

# bbt_id = bbt.compute_first_entity_collision(Point(point))
# if bbt_id<mesh.num_cells():
#     val = np.zeros(1)
#     _c_4.eval_cell(val, p, Cell(mesh, bbt_id))
#     print("Process {}: _c_4(x) = {}".format(rank, val[0]))
#     _c_3.eval_cell(val, p, Cell(mesh, bbt_id))
#     print("Process {}: _c_3(x) = {}".format(rank, val[0]))
#     d_c4.eval_cell(val, p, Cell(mesh, bbt_id))
#     print("Process {}: dc4_dy (z=0) = {}".format(rank, val[0]))
#     d_c3.eval_cell(val, p, Cell(mesh, bbt_id))
#     print("Process {}: dc3_dy (z=0) = {}".format(rank, val[0]))
    

## Get current in parallel -- Lacks slight accuracy in result, use post-processing get_current.py                                       
# r_list = np.linspace(0,a_s,100)
# th_list = np.linspace(0,2*np.pi,100)
# Jo_mat = np.zeros((len(r_list), len(th_list)))
# d_c3.set_allow_extrapolation(True)

# def xy_coord(r, th):
#     return [r*np.cos(th), r*np.sin(th)]

# for j in range(len(r_list)):
#     for k in range(len(th_list)):
#         xy = xy_coord(r_list[j], th_list[k])
#         point = np.array([xy[0], xy[1], 0.0])
#         xy_point = Point(point)
#         bbt_id = bbt.compute_first_entity_collision(xy_point)
#         val = np.zeros(1)
#         if bbt_id<mesh.num_cells():
#             d_c3.eval_cell(val, point, Cell(mesh, bbt_id))
#             Jo = -D_o*val[0]
#             Jo_mat[j, k] = Jo
#     else:
#             Jo_mat[j, k] = 0

# data_j = comm.gather(Jo_mat, root=0)
# #print(data_j)                                                                                                                          

# if rank==0:
#     #print(np.shape(data_j))
#     #print(np.shape(np.sum(data_j, axis=0)), '\n')
#     #print(Jo_mat)
#     Jo_mat = np.sum(data_j, axis=0)                                                                                                                
#     y = np.trapz(Jo_mat, x=th_list, axis=1)*r_list
#     curr = F*np.trapz(y, x=r_list)
#     print("I = {} Vapp = {} \n".format(curr, Vapp))
#     comm.bcast(curr, root=0)
#     #print('\n')
#     with open(iv_file, 'a') as f:
#         f.write(str(curr)+' '+str(Vapp)+' '+str(center[0])+' '+str(center[1])+'\n')
