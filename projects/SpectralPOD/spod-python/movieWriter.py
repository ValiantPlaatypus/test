# This example uses a MovieWriter directly to grab individual frames and
# write them to a file. This avoids any event loop integration, but has
# the advantage of working with even the Agg backend. This is not recommended
# for use in an interactive setting.
# -*- noplot -*-


from import_hdf5 import *

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

# === Analysis ===
showMovie = True    # False
saveMovie = True    # False

# === Import fields ===
folder = '/home/davide/PROJECTS/SpectralPOD/spod-towne/' + \
         'spod_matlab-master/jet_data/'
filen = folder + 'jetLES.mat'

f = import_hdf5(filen,True)


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




