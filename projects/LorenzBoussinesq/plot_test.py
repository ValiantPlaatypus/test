
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from matplotlib.gridspec import GridSpec

# plt.ion()

# Lines to plot in 3D
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

# === Plot section ===
xv = x[:,0]
yv = x[:,1]
zv = x[:,2]

# fig, ax = plt.subplots(, figsize=plt.figaspect(1./2.))
# ax0 = fig.add_subplot(1,2,1)
# ax1 = fig.add_subplot(1,2,2,projection='3d')

# domain (x,z) \in (0,nx*2/k) x (0,1) , nx in \mathbb{N}
a = np.sqrt(2.) *( k**2. + 1. ) / k * x[:,0]
b = np.sqrt(2.) *( k**2. + 1. ) / k * \
    ( np.pi**3 * ( k**2. + 1. )**2. ) / ( k     * Ra ) * x[:,1]
c = ( np.pi**3 * ( k**2. + 1. )**3. ) / ( k**2. * Ra ) * x[:,2]

nx = 4.
x0 = 0. ; x1 = nx * 2. / k
z0 = 0. ; z1 = 2.

# grid
dx = 0.05
xx = np.arange(x0,x1+0.001,dx)
zz = np.arange(z0,z1+0.001,dx)
X, Z = np.meshgrid(xx, zz)


fig = plt.figure(101, figsize=plt.figaspect(1./2.))

# gs1 = GridSpec(1,4)
# gs1.update(left=0.05, right=0.48, wspace=0.05)
ax0 = fig.add_subplot(1,2,1)
ax1 = fig.add_subplot(1,2,2, projection='3d')
ax0.grid(True), ax0.set_aspect(1.0)
ax1.grid(True), ax1.set_aspect(1.0)
ax0.set_title(' aaa ')
ax1.set_title(' bbb ')
# line0, = ax0.contourf([],[],[],15)
line1, = ax1.plot(x[:,0],x[:,1],x[:,2])

def update_lines(it, xv,yv,zv, X,Z, a,b,c, k, line1):
    num = 5 *it
    print(' num: ', num,' it: ', it)
    tau = b[num] * np.sin( np.pi * Z ) * np.cos( k * np.pi * X ) + \
          c[num] * np.sin( 2.*np.pi * Z )
#   line0.set_data( xv[:num] , yv[:num] )

    ax0.cla() 
    ax0.contourf(X,Z, tau,15)
#   line0.set_data( X, Z, tau )
    line1.set_data( xv[:num],  yv[:num] )
    line1.set_3d_properties(zv[:num])
#     line.axes.axis([-50, 50, -50, 50])
    return line1

line_ani = animation.FuncAnimation(fig, update_lines,  300, \
          fargs=[ xv, yv, zv, X,Z, a,b,c,k, line1], interval=25, blit=False, repeat=False )

line_ani.save('test.mp4')


plt.show()
