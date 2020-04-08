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
import matplotlib.animation as animation


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

x0 = np.array([ -10., 10., 5.])
t0 = 0. ; tend =  30. ; dt = .01  ; tv = np.arange( t0, tend, dt)

x = odeint( lorenz, x0,tv, args=(sig, rho, beta) )
print(' x.shape: ', x.shape)

# fig = plt.figure(101)
# ax = fig.gca(projection='3d')
# ax.plot(x[:,0],x[:,1],x[:,2]) # ,projection='3d')
# ax.grid(True)
# plt.draw()


fig = plt.figure(101)   # subplots()
ax = fig.add_subplot(1,1,1, projection='3d')
line, = ax.plot(x[:,0], x[:,1], x[:,2]) # , x[:,2], color='k')
ax.set_title(' aaa ')

def update(num, x, tv, ax, line):
    line.set_data( x[:num,0:2] ) # , x[:num,2])
    line.set_3d_properties(data[2, :num])
    ax.set_title(' t: %d' % (num) )
#   line.axes.axis([0, 10, 0, 1])
    return line,

ani = animation.FuncAnimation(fig, update, x.shape[0], fargs=[x, tv, ax, line],
                              interval=25, blit=False)
# ani.save('test.gif')
plt.show()
