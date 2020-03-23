# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# Solve the Lorenz system and represent the solution of the Boussinesq eqn     #
#  for convections.                                                            #
#                                                                              #
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #


import numpy as np
from scipy.integrate import odeint
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# === Lorenz dynamical system ===
# dX/dt = - sig*X + sig*Y
# dX/dt = -     Y + X*(rho-Z)
# dZ/dt = -beta*Z + X*Y
# 
# === Parameters ===
# beta = 4 / ( k^2 + 1 )
# rho = Ra / pi * k^2 / ( k^2 + 1 )^3
# sig = Pr
# Reduced parameters: let sig=10 and beta=8/3.
sig = 10. ; rho = 28. ; beta=8./3.
# Physical parameters
k  = np.sqrt( 4. / beta - 1. )                # wavelength  ( max beta = 4. !!!)
Ra = np.pi**4 * ( k**2 + 1 )**3 / k**2 * rho  # Rayleigh number
Pr = sig                                      # Prandtl number

# Integrate the Lorenz dynamical system
def lorenz( x,t, sig,rho,beta ):
    dxdt = np.array([ sig*(x[1]-x[0]) , \
                      -x[1]+x[0]*(rho-x[2]) , \
                      -beta*x[2] + x[0]*x[1] ] ) 
    return(dxdt)

x0 = np.array([ -10., 10., 5.])line_ani.save('line_animation_3d_funcanimation.mp4', writer='ffmpeg',fps=1000/100)
line_ani.save('line_animation_3d_funcanimation.gif', writer='imagemagick',fps=1000/100)
t0 = 0. ; tend =  30. ; dt = .01  ; tv = np.arange( t0, tend, dt)

x = odeint( lorenz, x0,tv, args=(sig, rho, beta) )
print(' x.shape: ', x.shape)

fig = plt.figure(101)
ax = fig.gca(projection='3d')
ax.plot(x[:,0],x[:,1],x[:,2]) # ,projection='3d')
ax.grid(True)
plt.draw()

plt.show()
stop


# Physical problem (non-dimensional variables)
# Fourier-Galerkin truncated approximation:
#   psi = a * sin(pi*z) * sin(k*pi*x)
#   tau = b * sin(pi*z) * cos(k*pi*x) + c * sin(2*pi*z)
# u = - d psi / dz =
# w =   d psi / dx =
# T = tau + 1 * z

a = np.sqrt(2.) *( k**2. + 1. ) / k * x[:,0]
b = np.sqrt(2.) *( k**2. + 1. ) / k * \
    ( np.pi**3 * ( k**2. + 1. )**2. ) / ( k     * Ra ) * x[:,1]
c = ( np.pi**3 * ( k**2. + 1. )**3. ) / ( k**2. * Ra ) * x[:,2]

fig = plt.figure(102)
plt.plot(tv,x[:,2]) # ,projection='3d')
plt.grid(True)
plt.draw()

# domain (x,z) \in (0,nx*2/k) x (0,1) , nx in \mathbb{N}
nx = 4.
x0 = 0. ; x1 = nx * 2. / k
z0 = 0. ; z1 = 2.

# grid
dx = 0.05
xv = np.arange(x0,x1+0.001,dx)
zv = np.arange(z0,z1+0.001,dx)
X, Z = np.meshgrid(xv, zv)
print(' X.shape: ', X.shape)
print(' Z.shape: ', Z.shape)

# # plot a snapshot
# it = 0 # initial conditions
# tau = b[it] * np.sin( np.pi * Z ) * np.cos( k * np.pi * X ) + \
#       c[it] * np.sin( 2.*np.pi * Z )
# T = tau + 1. * Z
# print(' tau.shape: ', tau.shape)

for it in range(0,tv.shape[0] , 10):
    tau = b[it] * np.sin( np.pi * Z ) * np.cos( k * np.pi * X ) + \
          c[it] * np.sin( 2.*np.pi * Z )
    T = tau + 1. * Z

    plt.figure(201)
    plt.subplot(121), plt.contourf( X, Z, tau, 15), plt.title(' %d ' % (it) )
    plt.clim(-0.2,0.2) # , plt.axis('scaled')
#   plt.autoscale(enable=True, axis='y', tight=True)
#   plt.colorbar()
#   plt.subplot(212), plt.contourf( X, Z, T )
#   plt.clim(-0.1,2.1), plt.axis('scaled')
#   plt.autoscale(enable=True, axis='y', tight=True)
    plt.subplot(122) # ; ax = fig.gca(projection='3d')
    plt.plot(tv[ :],x[ :,2]), plt.grid(True)
    plt.plot(tv[it],x[it,2], 'o')
    plt.draw()
    plt.pause(0.02)

plt.show()
