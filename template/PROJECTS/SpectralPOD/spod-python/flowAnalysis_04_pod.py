# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# Import data provided by Towne et al. and perform Spectral-POD                #
#                                                                              #
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

from import_hdf5 import *

import h5py
import time
import numpy as np
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt

from pod import pod_time
from spod import spod

# === Analysis parameters ===
npod = 300
nt_proj = 2500
npod_plot = 6

# === Import fields ===
folder = '/home/davide/PROJECTS/SpectralPOD/spod-towne/' + \
         'spod_matlab-master/jet_data/'
filen = folder + 'jetLES.mat'

f = import_hdf5(filen,True) 

# === plot mean pressure field ===
cmap = plt.get_cmap('PiYG')
fig, ax = plt.subplots()
ax.axis('equal')
ax.autoscale(enable=True, axis='y', tight=True)
cs = ax.contourf(f["x"],f["r"],f["p_mean"])

# === POD: Proper Orthogoanl Decomposition of p ===
pres = f["p"] ; nt = pres.shape[2]   # nt: n. of time intervals 
dt = f['dt'][0,0]
uPOD, s, vh = pod_time(pres[:,:,0:npod],dt)
plt.figure(101, figsize=(4,9))
for ip in range(npod_plot):
    plt.subplot(npod_plot,1,ip+1)
    plt.contourf(f["x"],f["r"],-uPOD[:,:,ip])
    plt.axis('scaled')
    ax.autoscale(enable=True, axis='y', tight=True)
    plt.ylabel('r')
    if ( ip == 0):
        plt.title(' First %d POD modes ' % (npod_plot) )
    if ( ip == npod_plot-1 ):
        plt.xlabel('x')


# === Project flow field onto the POD === 
# - Fields are reshaped as monodimensional arrays (for each time interval or
#   POD mode)
# - Discrete inner product is used for projection (u,v) = \sum_i u_i* v_i:
#   other definitions of inner products may be used, here and in pod computation
#   to obtain a consistent process.
print(' pres.shape: ', pres.shape)
p   = np.reshape(pres, ( np.product(pres.shape[0:2]), nt) )
pod = np.reshape(uPOD, ( np.product(uPOD.shape[0:2]), npod) ) 
print(' p.shape: ', p.shape)
print(' pod.shape: ', pod.shape)

# - Initialize the array of the projection coefficients of the flow
# - remove mean field, relying on broadcasting
p_prime = p - np.reshape( np.mean(p, axis=1) , ( p.shape[0] , 1) )
podAmpl = np.zeros( (npod, nt_proj) )
for it in range(nt_proj):
    print(' it: %d / %d' % ( it+1 , nt_proj ), end='\r')
    podAmpl[:,it] = np.matmul( pod.T , p[:,it])

# === Plots ===
# avoid first POD ~ mean field
# - Plot time evolution of the coefficients
plt.figure(201, figsize=(8,5))
plt.plot(podAmpl[1:npod_plot,:].T) , plt.grid(True)
plt.xlabel('it'), plt.ylabel('a(t)')
plt.title(' POD amplitude ')

# - L2 norm of the projection coefficients
podAmpl_norm = np.sum( podAmpl**2. , axis=1 )
plt.figure(202, figsize=(8,5))
plt.semilogy(podAmpl_norm[1:],'o') , plt.grid(True)
plt.xlabel('i.POD'), plt.ylabel('$|a_i(t)|$')
plt.title('$L_2$-norm of POD ampl')

# - Fourier analysis of the projection coefficients on PODs
T = dt * nt_proj
dOm = 2.*np.pi / T # Frequency Resolution
om = np.arange(0,nt_proj) * dOm
P = np.fft.fft(podAmpl[1:npod_plot+1,0:nt_proj]) / nt_proj
plt.figure(203, figsize=(4,9))
for ip in range(npod_plot):
    plt.subplot(npod_plot,1,ip+1)
    plt.plot(om, np.abs(P[ip,:]))
    plt.ylabel('|A|')
    if ( ip == 0):
        plt.title(' FFT of amplitude of first %d POD modes ' % (npod_plot) )
    if ( ip == npod_plot-1 ):
        plt.xlabel('\omega')



plt.show()
