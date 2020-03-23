# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# Import data provided by Towne et al. and compute the spectrum of pointwise   #
# measurements                                                                 #
#                                                                              #
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

from import_hdf5 import *
from window_functions import *

import numpy as np
import matplotlib.pyplot as plt

# === Analysis ===
# 1. Simple Fourier analysis of the whole time history
analysis_fourier_simple = False
# 2. Welch's method to get a low-noise estimation of the spectrum
analysis_fourier_welch  = True
blockSize = 512 
nOverlap = 256 
window_type = 'hamming'

# === Desired probes ===
rProbes = np.array([[ 5., 0. ],
                    [ 5.,  .5],
                    [10., 0. ],
                    [10.,  .5]] )
nProbes = rProbes.shape[0]
print('\n nProbes: ', nProbes)

# === Import fields ===
folder = '/home/davide/PROJECTS/SpectralPOD/spod-towne/' + \
         'spod_matlab-master/jet_data/'
filen = folder + 'jetLES.mat'

f = import_hdf5(filen,False)  # verbose=False

# === Actual probes ===
# - coinciding with the closest grid nodes to the desired ones:
#   -> save the ix,ir indices in iProbes array
iProbes = np.zeros(rProbes.shape, dtype=int)
for i in range(nProbes):
    iProbes[i,0] = np.argmin(np.abs(f['x'][:,0]-rProbes[i,0]))
    iProbes[i,1] = np.argmin(np.abs(f['r'][1,:]-rProbes[i,1]))

# === Time and Angular Velocity vectors ===
# - vectors of time intervals where the sampled function is known and
# - vectors of angular velocity induced by the observation period tend
dt = f['dt'][:,:].item(0,0)  # HDF5 dataset ~ NumPy array
                             # .item to get dt as a scalar (not an np.array)
t0 = 0. ; Nt = f['p'].shape[2] ; tend = t0 + dt*(Nt-1)
dOm = 2.*np.pi / tend        #  Frequency resolution 
t = np.linspace(t0,tend+dt,num=Nt, endpoint=False)
om = np.arange(0,t.shape[0]) * dOm

# === Fourier Analysis ===
# 1. Raw analysis on the whole sampled process
if ( analysis_fourier_simple ):

    print('\n === Simple Fourier Analysis === \n')

    plt.figure(1, figsize=(8,8) )
    for i in range(nProbes):
    
        ix = iProbes[i,0] ; ir = iProbes[i,1]
        p = f['p'][ix, ir, :]
        print(' probe n. ', i,' (x,r) = (',f['x'][ix,0],',',f['r'][0,ir],')')
        
        # fft of the pressure measurements
        p = p - np.mean(p)
        P = np.fft.fft(p) / Nt
    
        plt.subplot(2,2,i+1)
        plt.plot(om,np.abs(P))
        plt.title('(x,y)=(%5.2f,%5.2f)' % (f['x'][ix,0], f['r'][0,ir]) )
        plt.xlabel('$\omega$'), plt.grid(True)
    
    print('\n')

# 2. Welch's method
if ( analysis_fourier_welch ):

    print('\n === Welch\'s method === \n')

    # - Compute n. of blocks 
    nBlocks = int( np.floor(Nt-blockSize)/(blockSize-nOverlap) + 1 )
    # - Initialize array for average Fourier transform
    welch_P_ave = np.zeros(blockSize)
    welch_dOm = 2.*np.pi / (dt*blockSize)     # Frequency Resolution (lower)
    welch_om = np.arange(0,blockSize) * welch_dOm 
    print(' Total n. of time snapshots : ', Nt)
    print(' Block n. of time snapshots : ', blockSize)
    print(' N. of overlapping snashots : ', nOverlap)
    print('               ---> nBlocks : ', nBlocks)
    print('               ---> freq.Res: ', welch_dOm)
    
    win = window(blockSize,window_type)

    plt.figure(2, figsize=(8,8) )
    for i in range(nProbes):
    
        ix = iProbes[i,0] ; ir = iProbes[i,1]
        p = f['p'][ix, ir, :]
    
        for ib in range(nBlocks):
        
            # First and last index of the block for slicing ---
            ind0 =  ib*(blockSize-nOverlap)    # Python starts from 0
            ind1 = ind0+blockSize              # Python counts the separators
        
            # Slicing and Windowing ---
            p_block = p[ind0:ind1] 
            p_block = p_block - np.mean(p_block)  # remove mean value
            p_block = p_block * win
            welch_P = np.fft.fft(p_block) / blockSize
        
            welch_P_ave = welch_P_ave + welch_P
    
        welch_P_ave = welch_P_ave / nBlocks
    
        plt.subplot(2,2,i+1)
        plt.stem(welch_om,np.abs(welch_P_ave), markerfmt=' ')
        plt.title('(x,y)=(%5.2f,%5.2f)' % (f['x'][ix,0], f['r'][0,ir]) )
        plt.xlabel('$\omega$'), plt.grid(True)



plt.show()



