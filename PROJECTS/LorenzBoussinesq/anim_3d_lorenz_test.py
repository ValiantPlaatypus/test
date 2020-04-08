# 3d lines animation from
# https://pythonmatplotlibtips.blogspot.com/2017/12/draw-3d-line-animation-using-python-matplotlib-funcanimation.html

# [1]
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from IPython.display import HTML

# [2]
def update_lines(num, dataLines, lines):
    print(' num: ', num)
    for line, data in zip(lines, dataLines):
        # NOTE: there is no .set_data() for 3 dim data...
        line.set_data(data[0:2, :num])
        line.set_3d_properties(data[2, :num])
#       line.set_marker("o")
    return lines

# [3]
# Attaching 3D axis to the figure
fig = plt.figure()
ax = p3.Axes3D(fig)

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
data = np.array([np.transpose(x)])
print(data.shape)
# --- example ----
#t = np.linspace(-2*np.pi,2*np.pi,50)
#x1, y1, z1 = np.cos(t), np.sin(t), t/t.max()
#x2, y2, z2 = t/t.max(), np.cos(t), np.sin(t)
#data = np.array([[x1,y1,z1],[x2,y2,z2]])
#print(data.shape)

# NOTE: Can't pass empty arrays into 3d version of plot()
lines = [ax.plot(dat[0, 0:1], dat[1, 0:1], dat[2, 0:1])[0] for dat in data]

ax.set_xlim(-50.,50.)
ax.set_ylim(-50.,50.)
ax.set_zlim(-50.,50.)
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
plt.rcParams['animation.html'] = 'html5'

line_ani = animation.FuncAnimation(fig, update_lines, 3000, fargs=(data, lines),
                                   interval=25, blit=True, repeat=True)

# The empty figure shown below is unknown.
# Does anyone know how to delete this?

# [4]
line_ani

# [5]
line_ani.save('line_animation_3d_funcanimation.mp4', writer='ffmpeg',fps=25)
# line_ani.save('line_animation_3d_funcanimation.gif', writer='imagemagick',fps=1000/100)

plt.show()



