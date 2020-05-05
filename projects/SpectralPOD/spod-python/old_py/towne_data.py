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

# beginning of the program -----
print('\n Start of the program.\n')

# read data --------------------
filen = '/home/davide/Software/spod-towne/spod_matlab-master/jet_data/jetLES.mat'
f = h5py.File(filen, 'r') # h5py.File acts like a Python dictionary

print('list(f.keys()):', list(f.keys()))
print('f:', f)

for i1 in range(0,len(f)):
  print(' field', i1,'. key: ',list(f.keys())[i1], \
        ' shape: ',f[list(f.keys())[i1]].shape)

# plot mean pressure field -----
#cmap = plt.get_cmap('PiYG')
#fig, ax = plt.subplots()
#ax.axis('equal')
#ax.autoscale(enable=True, axis='y', tight=True)
#cs = ax.contourf(f["x"],f["r"],f["p_mean"])
#plt.show()

# plot movie pressure time evolution ----
dt = f['dt'][0][0]
nt = f['p'].shape[2]  # f['nt'][0][0]
print('nt: ',f['p'].shape[2],', dt: ',dt)
fig_mov = plt.figure()
pmin = 4.33963 # np.amin(f['p'])
pmax = 4.51078 # np.amax(f['p'])
print('min(p):', pmin)
print('max(p):', pmax)

nclevs = 30
#plt.ion()
def animate():
    for i in range(0,int(nt),20):
        print('i: ',i,'. t = ',i*f['dt'][0,0])
#       im=plt.imshow(f["p"][:,:,i])
        im=plt.contourf(f["x"],f["r"],f["p"][:,:,i], nclevs, \
             vmin=pmin, vmax=pmax , \
             cmap=plt.cm.bone )
        fig_mov.canvas.draw()
#       time.sleep(1.0)

win = fig_mov.canvas.manager.window
fig_mov.canvas.manager.window.after(1000, animate)
#plt.axis('equal')
plt.axis([0, 20, 0, 2.9],'equal')
plt.gca().set_aspect('equal', adjustable='box')
plt.show()



print('\n End of the program. Bye!\n')
# end of the program -----
