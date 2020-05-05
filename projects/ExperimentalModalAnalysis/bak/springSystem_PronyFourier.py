
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal


npzfile = np.load('./data/springTime.npz')
print(npzfile.files)
t = npzfile['t'] ; print('t.shape: ') ; print(t.shape) # time
u = npzfile['u'] ; print('u.shape: ') ; print(u.shape) # input
x = npzfile['x'] ; print('x.shape: ') ; print(x.shape) # state (x,v)

N = t.shape[0]

# # === Select 2^n timesteps for fft ===
# # slice for fft (removing the initial transient or the final part with the
# # largest frequencies?)
# n_fft = int(2**(np.floor(np.log(N)/np.log(2))))
# print(' n_fft: ', n_fft)
# t = t[0:n_fft]    #  [N-n_fft:N]
# u = u[0:n_fft]    #  [N-n_fft:N]
# x = x[0:n_fft,:]  #  [N-n_fft:N,:]
# print(' x.shape: ', x.shape)

# from time to frequency domain
dt = t[1]-t[0] ; N = t.shape[0] ; To = (N+1)*dt
print(' dt: %8.3f , N: %d , To: %8.3f' % (dt,N,To))
df = 1. / To    # Frequency Resolution
f = np.arange(0,N) * df
print(' df: %8.3f' % (df))
print(' max(f)/2: %8.3f' % (.5/dt) )

U = np.fft.fft(u) ; print(' U.shape:', U.shape)
X = np.fft.fft(x.T) ; print(' X.shape:', X.shape)

# === Compute FRF and IRF ===
H = X / U  # broadcasting 

# inverse Fourier tranform to obtain the IRF
ut = np.fft.ifft(U) ; print(' ut.shape:', ut.shape)
ht = np.fft.ifft(H) ; print(' xt.shape:', ht.shape)

maximag = np.amax(np.abs(np.imag(ht)))
print('max imag(IRF(t)): ', maximag )
if ( maximag < 1.e-7 ):
    ht = np.real(ht)
else:
    print(' error. maximag > 1.e-7')
    stop

np.savez('./data/springFRF.npz', f=f, U=U, X=X, H=H)
np.savez('./data/springIRF.npz', t=t, ht=ht)

# Plot FRFs ----
plt.figure(1, figsize=(5,9))
plt.subplot(511)
plt.loglog(f,np.abs(H[0,:])), plt.grid(True)
plt.subplot(512)
plt.loglog(f,np.abs(H[1,:])), plt.grid(True)
plt.subplot(513)
plt.loglog(f,np.abs(H[2,:])), plt.grid(True)
plt.subplot(514)
plt.loglog(f,np.abs(H[3,:])), plt.grid(True)
plt.subplot(515)
plt.loglog(f,np.abs(H[4,:])), plt.grid(True), plt.xlabel('f')
# Plot IRFs ----
plt.figure(101, figsize=(5,9))
plt.subplot(511), plt.title(' IRF ')
plt.plot(t,ht[0,:]), plt.grid(True)
plt.subplot(512)                   
plt.plot(t,ht[1,:]), plt.grid(True)
plt.subplot(513)                   
plt.plot(t,ht[2,:]), plt.grid(True)
plt.subplot(514)                   
plt.plot(t,ht[3,:]), plt.grid(True)
plt.subplot(515)                   
plt.plot(t,ht[4,:]), plt.grid(True), plt.xlabel('t')

plt.show()


# plt.figure
# plt.plot(f,np.abs(U))
# 
# plt.figure(figsize=(4,9))
# plt.subplot(511)
# plt.plot(f,np.abs(X[0,:]))
# plt.subplot(512)
# plt.plot(f,np.abs(X[1,:]))
# plt.subplot(513)
# plt.plot(f,np.abs(X[2,:]))
# plt.subplot(514)
# plt.plot(f,np.abs(X[3,:]))
# plt.subplot(515)
# plt.plot(f,np.abs(X[4,:]))


# H3 = X[2,:]/U
# plt.figure(201)
# plt.subplot(211), plt.loglog(om,np.abs(H3)), plt.grid(True)
# plt.subplot(212), plt.semilogx(om,np.unwrap(np.angle(H3)),'-o'), plt.grid(True)
# plt.figure(202)
# plt.subplot(211), plt.plot(om,np.real(H3)), plt.grid(True)
# plt.subplot(212), plt.plot(om,np.imag(H3)), plt.grid(True)
# plt.figure(203)
# plt.plot(np.real(H3),np.imag(H3),'-o'), plt.grid(True)

plt.show()


