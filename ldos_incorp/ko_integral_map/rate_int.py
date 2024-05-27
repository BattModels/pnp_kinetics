import numpy as np
import os, re, pickle, glob, sys, mat73
from scipy.io import loadmat, savemat
from funcs_pnp import *
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
#mpl.rcParams['xtick.minor.size'] = 10
#mpl.rcParams['xtick.minor.width'] = 2

plt.legend(frameon=False)

## User functions
def get_rate_integral(interp_func, center = [0,0], a_s = 100):
    # 2D interp function (takes x,y as args), and nanopipette radius in Ang
    r_list = np.linspace(0,a_s,200)
    th_list = np.linspace(0,2*np.pi,200)
    Jo_mat = np.zeros((len(r_list), len(th_list)))
    
    def xy_coord(r, th):
        return [r*np.cos(th), r*np.sin(th)]
    
    for j in range(len(r_list)):
        for k in range(len(th_list)):
            xy = xy_coord(r_list[j], th_list[k])
            point = np.array(xy)
            Jo = interp_func(point[0] + center[0], point[1] + center[1])
            Jo_mat[j, k] = Jo
    
    y = np.trapz(Jo_mat, x=th_list, axis=1)*r_list
    int = np.trapz(y, x=r_list)
    return int

#MHC_prefactor = 1
#center = center_domain(bvec1, domain='AB')
#center = np.array([0, 0]) # Center of nanopipette in Angstroms     
#rate_int = get_rate_integral(ip_ox, center = center, a_s = a_s)
#print(bvec1)
domains = ['AA', 'AB']

thetas = [0.77, 1.1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8]
a_s = 50 # Nanopipette radius in Angstroms
AA_int, AB_int = [], []

for theta in thetas:
    print('Theta = {}'.format(theta))
    for dom in domains:
        rate_data = mat73.loadmat('/Users/mbabar/Desktop/PhD/Analysis/PDE/PNP_solve/3D_vmg/ldos_incorp/k_data_nr_48/{}/k_data_Vapp_0.0.mat'.format(str(theta)))
        MHC_prefactor = '/Users/mbabar/Desktop/PhD/Analysis/PDE/PNP_solve/3D_vmg/prefactor/pref_data/pref_sc_{}.mat'.format(str(theta))
        
        if type(MHC_prefactor) == str:
            pref_data= loadmat(MHC_prefactor)
            A = np.squeeze(pref_data['pref_sc'])
        else:
            A = float(MHC_prefactor)
        kred_data, kox_data = A*rate_data["kred_data"], A*rate_data["kox_data"]
        #kred_data, kox_data = A*np.ones(np.shape(kred_data)), A*np.ones(np.shape(kox_data))
        rscx, rscy = rate_data['rscx'], rate_data['rscy']
        Vapp = rate_data['Vapp']

        moir_l, sc_area = tblg_dim(theta)
        moir_a = np.pi/3
        orig = moir_l*np.array([rscx[0,0], rscy[0,0]]) # Set origin in Angstrom                                                                       
        bvec1 = np.array([moir_l*np.cos(moir_a/2), moir_l*np.sin(moir_a/2)]) # basis vector 1                                                      
        bvec2 = np.array([moir_l*np.cos(moir_a/2), -moir_l*np.sin(moir_a/2)]) # basis vector 2                                                                                                                            
        coords = moir_l*np.array([rscx.flatten(), rscy.flatten()])
        z_red, z_ox = kred_data.flatten(), kox_data.flatten()
        ip_red, ip_ox = my_CT(coords, z_red), my_CT(coords, z_ox) # 2D interpolators of kred and kox 

        center = center_domain(bvec1, domain=dom)
        rate_int = get_rate_integral(ip_ox, center = center, a_s = a_s)
        if dom == 'AA':
            AA_int.append(rate_int)
            print('AA int = {}'.format(rate_int))
        elif dom == 'AB':
            AB_int.append(rate_int)
            print('AB int = {}'.format(rate_int))


        #print('Theta = {}, Domain = {}, Rate Integral= {}'.format(theta, dom,rate_int))

AA_int = np.array(AA_int)
AB_int = np.array(AB_int)

plt.figure()
plt.plot(thetas, AA_int/AB_int, 'ro-')
plt.xlabel('theta')
plt.ylabel('AA/AB integral ratio')
plt.show()
