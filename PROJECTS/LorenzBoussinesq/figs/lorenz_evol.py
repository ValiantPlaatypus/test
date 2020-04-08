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

x01 = np.array([ -10., 10., 5.])
x02 = np.array([ -10., 10., 5.]) + 1.e-9
t0 = 0. ; tend = 100. ; dt = .01  ; tv = np.arange( t0, tend, dt)

x1 = odeint( lorenz, x01,tv, args=(sig, rho, beta) )
x2 = odeint( lorenz, x02,tv, args=(sig, rho, beta) )
print(' x1.shape: ', x1.shape)

fig = plt.figure(101)
ax = fig.gca(projection='3d')
plt.plot(x1[:,0],x1[:,1],x1[:,2], linewidth=.5) # ,projection='3d')
plt.plot(x2[:,0],x2[:,1],x2[:,2], linewidth=.5) # ,projection='3d')
plt.grid(True)
ax.set_xlabel('X'), ax.set_ylabel('Y'), ax.set_zlabel('Z')
ax.set_aspect('equal')

fig = plt.figure(102)
plt.plot(tv,x1[:,0]) # ,projection='3d')
plt.plot(tv,x2[:,0]) # ,projection='3d')
plt.grid(True)
plt.autoscale(enable=True, axis='x', tight=True)
plt.xlabel('t')
plt.ylabel('X(t)')

plt.show()
