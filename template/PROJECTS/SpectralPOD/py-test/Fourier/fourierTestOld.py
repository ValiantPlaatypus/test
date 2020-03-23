# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# Fourier Trasform in Python. Test numpy functions                             #
#                                                                              #
# fft transform and digital signals.                                           #
# - Test for aliasing, spectral leakage and scalloping                         #
# -"Physical" DFT to:                                                          #
#   - get the right physical units                                             #
#   - retrieve the continuous case in the limit dt -> 0                        #
#                                                                              #
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

import numpy as np
import matplotlib.pyplot as plt


# define the signal in time ------------
# - time vector ------
t0 = 0.0 ; tend = 50.*2*np.pi ; T = tend - t0 ; N = 2000 ; dt = tend/N
t = np.arange(t0,tend+dt,dt) ; N = len(t)
print(' T, dt, N: ', T,', ', dt,', ', N)
# - signal -----------
Am = np.array( [ 2.0 ,  1.0 ] ) 
Om = np.array( [ 1.0 ,  2.0 ] ) 
Ph = np.array( [ 0.0 , 45.0 ] ) # * np.pi/180.0 ;
f = np.zeros(t.shape)
for i in range(0,Am.shape[0]):
    f = f + Am[i]*np.sin(Om[i]*t-Ph[i]) 
# print(' f : ',f)

# plot signal
plt.figure(1, figsize=(9, 3))
plt.subplot(121)
plt.plot(t,f)
plt.axis('scaled'), plt.axis([t0, tend, -5, 5])

# fft transform of the time signal -----
F = np.fft.fft(f) / N #  N/2 
print(' f.shape : ', f.shape)
print(' F.shape : ', F.shape)

dOm = 2*np.pi 
om = np.arange(0,F.shape[0]) * dOm
plt.subplot(122)
plt.plot(om,np.real(F))
plt.plot(om,np.imag(F))
plt.plot(om,np.abs(F))
# plt.axis('square') # plt.axis([0, F.shape[0], F.min, F.max])

# Fourier transform as amplitude and phase
if np.mod(N,2) == 0:   # N even
    # ---
    print(' F.shape: ', F.shape)
    nOm = int(np.floor(N/2+1))
    print(' n.om, N: ', nOm, N)
    omReal = om[0:nOm+1]
    Fabs = F[0:nOm+1]
    Fabs[1:nOm] = np.abs(Fabs[1:nOm]) * 2. # + F[N:nOm-2:-1] 
    
elif np.mod(N,2) == 1:   # N odd
    # ---
    print(' F.shape: ', F.shape)
    nOm = int(np.floor(N/2+1))
    print(' n.om: ', nOm)
    omReal = om[0:nOm+1]
    Fabs = F[0:nOm+1]
    print(F[1:nOm])
    print(F[N:nOm-1:-1])
    Fabs[1:nOm] = np.abs(Fabs[1:nOm]) * 2. #  + np.conjugate(F[N:nOm-1:-1])

plt.figure(101)
plt.plot(omReal,np.abs(Fabs))

plt.show()

