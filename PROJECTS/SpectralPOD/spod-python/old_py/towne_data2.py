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
pres = pres[:,:,0:2000]
# uPOD, s, vh = pod_time(pres,dt)

# # Show some PODs
# for ipod in range(0,2):
#     fig, ax = plt.subplots()
#     ax.axis('equal')
#     ax.autoscale(enable=True, axis='y', tight=True)
#     cs = ax.contourf(f["x"],f["r"],uPOD[:,:,ipod])
# plt.show()

# FFT --------------------------
# fft of a block
# Reshape the input as a 2-dimensional array:
n = len(pres.shape)-1  # n. of 'non-time' dimensions
# 1st index for the unravelled 'non-time' dimensions
# 2nd index for the time variable 
p2D = np.reshape(pres, (np.product(pres.shape[0:n]),pres.shape[n]) )
n2D = p2D.shape[1]
print(' len(p2D.shape): ', len(p2D.shape))
print('     p2D.shape : ',     p2D.shape)

Nf = 512 
q1 = p2D[:,0:Nf]
meanq1 = np.mean(q1,axis=1)

#
# Reshape POD mode back to the original dimensions
meanq13D = np.reshape(meanq1, pres.shape[0:2])
cmap = plt.get_cmap('PiYG')
fig, ax = plt.subplots()
ax.axis('equal')
ax.autoscale(enable=True, axis='y', tight=True)
cs = ax.contourf(f["x"],f["r"],meanq13D)

print(' meanq1.shape : ', meanq1.shape)
for i in range(0,q1.shape[1]):    # remove mean value
    q1[:,i] = q1[:,i] - meanq1

Q1 = np.fft.fft(q1)
print(' q1.shape : ', q1.shape)
print(' Q1.shape : ', Q1.shape)
normQ1 = np.zeros(q1.shape[1])
for i in range(0,q1.shape[1]):
    normQ1[i] = np.sqrt(np.sum(np.conjugate(Q1[:,i]) * Q1[:,i]))
    print(' normQ1[',i,']:',normQ1[i])

T = dt * Nf           # observation period
dOm = 2*np.pi / T     # frequency resolution
om = np.arange(0,Q1.shape[1]) * dOm
plt.figure(101)
plt.plot(om,normQ1)
plt.show()

# find dominant mode
indmax = np.argmax(normQ1)
print(' indmax: ', indmax)
print(' indmax: ', indmax, '. omega: ', indmax*dOm)
# Reshape dominant mode back to the original dimensions
mode = np.reshape(Q1[:,indmax], pres.shape[0:2])
cmap = plt.get_cmap('PiYG')
plt.figure(102)
plt.subplot(211)
plt.contourf(f["x"],f["r"],np.real(mode))
plt.axis('scaled')
plt.subplot(212)
plt.contourf(f["x"],f["r"],np.imag(mode))
plt.axis('scaled')
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
