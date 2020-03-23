# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# Import data provided by Towne et al. and perform Spectral-POD       #
#                                                                     #
# do it in Python and train yourself                                  #
#                                                                     #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# import libraries -------------
#import sys
import h5py
import time
import numpy as np
import matplotlib
matplotlib.use('TkAgg') # do this before importing pylab
import matplotlib.pyplot as plt

from pod import pod_time
from spod import spod

# beginning of the program -----
print('\n Start of the program.\n')

# read data --------------------
filen = '/home/davide/PROJECTS/SpectralPOD/spod-towne/spod_matlab-master/jet_data/jetLES.mat'
f = h5py.File(filen, 'r') # h5py.File acts like a Python dictionary

print('list(f.keys()):', list(f.keys()))
print('f:', f)

for i1 in range(0,len(f)):
  print(' field', i1,'. key: ',list(f.keys())[i1], \
        ' shape: ',f[list(f.keys())[i1]].shape)

# plot mean pressure field -----
cmap = plt.get_cmap('PiYG')
fig, ax = plt.subplots()
ax.axis('equal')
ax.autoscale(enable=True, axis='y', tight=True)
cs = ax.contourf(f["x"],f["r"],f["p_mean"])


# Proper Orthogoanl Decomposition of p
pres = f["p"]
dt = f['dt'][0,0]
print('dt:', dt)
# only on the first 100 snapshots
# pres = pres[:,:,0:100]
uPOD, s, vh = pod_time(pres[:,:,0:300],dt)

# SPOD -------------------------
nSnap1 = 256 ; nOverlap = 128 ; nSnap = 5000  # pres.shape[len(pres.shape)-1]
PPsi, LLambda = spod(pres[:,:,0:nSnap],nSnap1,nOverlap,dt)

T = nSnap1 * dt
dOm = 2*np.pi / T
om = np.arange(0,nSnap1) * dOm
print(' dOmega: ', dOm)

print('LLambda.shape: ', LLambda.shape)

plt.figure(101)
plt.semilogy(om,LLambda[:,:])
# plt.plot(om,LLambda[:,:])

# plot sum of lambda
plt.figure(102)
plt.semilogy(om,np.sum(LLambda,axis=1))

# plot mode contributions at peak frequency
isort = np.argsort(LLambda[0:int(np.ceil(nSnap1/2)),0]) ; isort = isort[::-1]
    
# print(' omega: ', dOm*indmax )
# plt.figure(103)
# plt.semilogy(np.arange(0,LLambda.shape[1]),LLambda[indmax,:],'o')

nplot = 2
for ip in range(0,nplot):
    indMode = isort[ip]
    # Plot mode contribution at the considered frequency
    plt.figure(51+ip)
    plt.semilogy(np.arange(0,LLambda.shape[1]),LLambda[indMode,:],'o')
    plt.title(' Mode contributions at $\omega$ = %8.3f' % (dOm*indMode))

    # Plot optimal mode
    # reshape back to bi-dimensional field
    OptMode = np.reshape(PPsi[:,indMode,1], pres.shape[0:2])
    plt.figure(111+ip, figsize=(9, 4.5))
    plt.subplot(211)
    plt.contourf(f["x"],f["r"],np.real(OptMode))
    plt.axis('scaled')
    ax.autoscale(enable=True, axis='y', tight=True)
    plt.title(' Optimal mode: $\omega$ = %8.3f. Real part' % (dOm*indMode))
    plt.ylabel('r')
    
    plt.subplot(212)
    plt.contourf(f["x"],f["r"],np.imag(OptMode))
    plt.axis('scaled')
    ax.autoscale(enable=True, axis='y', tight=True)
    plt.title(' Optimal mode: $\omega$ = %8.3f. Imag part' % (dOm*indMode))
    plt.xlabel('x')
    plt.ylabel('r')
    

plt.figure(201, figsize=(5, 10))
nPlot = 9
stepPlot = 2
for ip in range(1,nPlot):
    Mode = np.reshape(PPsi[:,ip*stepPlot,1], pres.shape[0:2])
    ax = plt.subplot(str(nPlot-1)+str(1)+str(ip))
    plt.contourf(f["x"],f["r"],np.imag(Mode))
    plt.axis('scaled')
    ax.autoscale(enable=True, axis='y', tight=True)
    plt.title('$\omega$ = %8.3f , $\lambda$ = %8.3f' \
            % (dOm*ip*stepPlot,LLambda[ip*stepPlot,0]))
    ax.title.set_fontsize(10)

plt.show()

# # low-dimensional example ----------------------
# data = np.array( [[[0, 1],
#                    [2, 3],
#                    [4, 5]],
#                   [[6, 7],
#                    [8, 9], 
#                    [10, 11]]] )
# datav_C = np.ravel(data,order='C')
# datav_F = np.ravel(data,order='F')
# print('data.shape: ', data.shape)
# print('datav_C.shape: ', datav_C.shape)
# print('datav_C      : ', datav_C      )
# print('datav_F.shape: ', datav_F.shape)
# print('datav_F      : ', datav_F      )
# pod_time(data)
# # low-dimensional example ----------------------

# # plt.show()
# 
# # plot movie pressure time evolution ----
# dt = f['dt'][0][0]
# nt = f['p'].shape[2]  # f['nt'][0][0]
# print('nt: ',f['p'].shape[2],', dt: ',dt)
# pmin = 4.33963 # np.amin(f['p'])
# pmax = 4.51078 # np.amax(f['p'])
# print('min(p):', pmin)
# print('max(p):', pmax)
# 
# fig_mov = plt.figure()
# nclevs = 30
# #plt.ion()
# def animate():
#     for i in range(0,int(nt),20):
#         print('i: ',i,'. t = ',i*f['dt'][0,0])
# #       im=plt.imshow(f["p"][:,:,i])
#         im=plt.contourf(f["x"],f["r"],f["p"][:,:,i], nclevs, \
#              vmin=pmin, vmax=pmax , \
#              cmap=plt.cm.bone )
#         fig_mov.canvas.draw()
# #       time.sleep(1.0)
# 
# win = fig_mov.canvas.manager.window
# fig_mov.canvas.manager.window.after(1000, animate)
# #plt.axis('equal')
# plt.axis([0, 20, 0, 2.9],'equal')
# plt.gca().set_aspect('equal', adjustable='box')
# plt.show()
# 
# 
# 
# print('\n End of the program. Bye!\n')
# # end of the program -----
