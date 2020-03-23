# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# Import data provided by Towne et al. and show the temporal evolution of the  #
# system                                                                       #
#                                                                              #
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

from import_hdf5 import *

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

# === Analysis ===
# showMovie = True    # False
# saveMovie = True    # False

# === Import fields ===
folder = '/home/davide/PROJECTS/SpectralPOD/spod-towne/' + \
         'spod_matlab-master/jet_data/'
filen = folder + 'jetLES.mat'

f = import_hdf5(filen,True)

# === Show and save a movie ===
FFMpegWriter = manimation.writers['ffmpeg']
metadata = dict(title='Movie Test', artist='Matplotlib',
                comment='Movie support!')
writer = FFMpegWriter(fps=25, metadata=metadata)

fig = plt.figure()
plt.xlim(np.amin(f['x']), np.amax(f['x']))
plt.ylim(np.amin(f['r']), np.amax(f['r']))
dp = np.amax(f['p']) - np.amin(f['p'])
minp = np.amin(f['p']) + .2 * dp
maxp = np.amax(f['p']) - .2 * dp

with writer.saving(fig, "writer_test.mp4", 100):  # fig, filen, dpi
    for i in range(0,500,5):
        print(' i: ', i)
        plt.contourf(f['x'], f['r'], f['p'][:,:,i], 25)
        plt.clim(minp, maxp)
        plt.axis('scaled')
        plt.autoscale(enable=True, axis='y', tight=True)
        writer.grab_frame()













# nx = f['x'].shape[0]
# nr = f['x'].shape[1]
# print(' x :')
# print(f['x'][0:3,0:2])
# print(f['x'][nx-3:nx,nr-2:nr])
# print(' r :')
# print(f['r'][0:3,0:2])
# print(f['r'][nx-3:nx,nr-2:nr])
# 
# # - Initialise domain and axis 
# fig = plt.figure()
# # plt.axes(xlim=(np.amin(f['x']), np.amax(f['x'])), \
# #          ylim=(np.amin(f['r']), np.amax(f['r'])) )  
# plt.contourf(f['x'], f['r'], f['p'][:,:,0], 25)
# plt.clim(np.amin(f['p']), np.amax(f['p']))
# plt.axis('scaled')
# plt.autoscale(enable=True, axis='y', tight=True)
# plt.xlabel(r'x') , plt.ylabel(r'r')
# 
# # animation function
# def animate(i):
#     print(' i: ', i)
#     cont = plt.contourf(f['x'], f['r'], f['p'][:,:,i], 25)
# 
#     return cont  
# 
# anim = animation.FuncAnimation(fig, animate, frames=10)
# 
# anim.save('sim_Jet.avi')
# 
# plt.show()
# 



