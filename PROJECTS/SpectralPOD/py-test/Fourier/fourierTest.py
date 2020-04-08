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

# === Noise amplitude ===
# - Amplitude of a (Psuedo)random noise with
#   - a uniform distribution (np.rand()) 
#   - a normal distribution (np.randn())
AmNoise = 0.1

# === Define the signal in time ===
# - Define the analytical (time continuous) signal as the sum of the sine fncs
# with
Am = np.array( [ 2. ,  1. ,  1. ] )             # amplitudes
Om = np.array( [ 1. ,  2. , 15. ] )             # angular velocities
Ph = np.array( [ 0. , 45. ,  0. ] ) *np.pi/180. # phases

# - Vector of time instants when the signal is sampled
t0 = 0.0 ; tend = 200. ; T = tend - t0 ; dt = .10
# -time observation containing an integer number of periods -> no aliasing
# t0 = 0.0 ; tend = 50.*2*np.pi ; T = tend - t0 ; N = 2000 ; dt = tend/N
t = np.arange(t0,tend,dt) ; N = len(t)

print('\n === Sampled signal === ')
print(' Observation time,   T: %8.3f' % (T) )
print(' Sampling time   ,  dt: ', dt)
print(' N.of samples    ,   N: ', N)
print(' Sampling freq   ,  fs: %8.3f ' % ( 1./dt) )
print(' Sampling ang.vel, Oms: %8.3f ' % (2.*np.pi/dt) )
print(' -> Frequency resolution, dOm: %8.3f ' % ( 1./T ) )
print(' -> Max. non-aliased freq.   : %8.3f ' % ( .5/dt ) )
print(' -> Max. non-aliased ang.vel : %8.3f ' % ( np.pi/dt ) )

# - Sampled (discrete time) signal
f = np.zeros(t.shape)
for i in range(0,Am.shape[0]):
    f = f + Am[i]*np.sin(Om[i]*t-Ph[i]) 

# - Add noise
f = f + 2.*AmNoise*( np.random.rand(t.shape[0]) - .5 )

# - plot signal
plt.figure(1, figsize=(9, 4))
plt.subplot(121)
plt.plot(t,f)
# plt.axis('scaled'), plt.axis([t0, tend, -5, 5])
plt.xlabel('t'), plt.title('x(t)')

# === Fourier transform of the signal ===
# - Remove average value
fAve = np.mean(f) ; f = f - fAve
# - Fast Fourier Transform (fft) is an anlgorithm to efficiently compute the 
# Discrete Fourier Transform (DFT) of a sequence {x_n}_{n=0:N-1}
#  F[k] = sum_{n=0}^{N-1} f[n] e^{i 2*pi*k*n/N} , k = 0:N-1             # (Eq.1)
F = np.fft.fft(f)
# - This fft produces a result that is not consistent with the physical units
# - To obtain a "phyisical-dimensional" correct fft, the multiplication by a
#   time is needed. If (Eq.1) is multiplied by dt, the DFT converges to the 
#   continuous transform for dt -> 0
# F = F * dt                                                            # (Eq.2)
# - In order to retrieve the "correct amplitude" of the Fourier transform, it is
#   needed to divide (Eq.1) by N. This is equal to the average of the Fourier
#   transform of (Eq.2) over the observation period T, since T = N*dt,
#   coinciding with the coefficients of the Fourier series.
F = F / N                                                             # (Eq.3)
# !!! It is importart to remind that the product of the multiplicative
# !!! coefficients of the DFT and the IDFT must be equal to 1/N       


# - Observation time T induces the Frequency resolution dOm
dOm = 2.*np.pi/T ; om = np.arange(0,F.shape[0]) * dOm

# - Plot Fourier transform of the signal: real and imag part and absolute value
plt.subplot(122)
plt.plot(om,np.real(F), om,np.imag(F), om,np.abs(F))
plt.xlabel('$\omega$'), plt.title('$\hat{x}(\omega)$')

# # - Plot Fourier transform of the signal: absolute value and phase
# plt.figure(2, figsize=(5,5))
# plt.subplot(211)
# plt.plot(om,np.abs(F))
# plt.xlabel('$\omega$'), plt.title('$|\hat{x}|(\omega)$')
# plt.subplot(212)
# plt.plot(om,np.angle(F),'o')
# plt.xlabel('$\omega$'), plt.title('$\\angle\hat{x}(\omega)$')


plt.show()

