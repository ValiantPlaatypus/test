
import numpy as np
import matplotlib.pyplot as plt

sig = 10.
beta = 8./3.
rho = np.arange(0.,2.05, .05)

print(' rho: ', rho)

eigs = np.zeros( (rho.shape[0],3), dtype=float )
eigs[:,0] = - beta
eigs[:,1] = -.5*(sig+1.) - .5*np.sqrt((sig+1.)**2.-4.*sig*(1-rho) )
eigs[:,2] = -.5*(sig+1.) + .5*np.sqrt((sig+1.)**2.-4.*sig*(1-rho) )

rhos = np.arange(0.,1.05, .05)
rhou = np.arange(1.,2., .01)
xeqs = np.zeros( (rhos.shape[0],1), dtype=float )
xequ = np.zeros( (rhou.shape[0],3), dtype=float )
xeqs[:,0] =  .0    
xequ[:,0] =  .0    
xequ[:,1] = - np.sqrt( beta * ( rhou - 1 ) )
xequ[:,2] = + np.sqrt( beta * ( rhou - 1 ) )

fig = plt.figure(1, figsize=plt.figaspect(1.5))
ax0 = fig.add_subplot(2,1,1)
ax1 = fig.add_subplot(2,1,2)

ax0.plot(rho,eigs[:,0] , 'C0' , linestyle='--'), ax0.grid(True)
ax0.plot(rho,eigs[:,1] , 'C0' , linestyle='-')
ax0.plot(rho,eigs[:,2] , 'C0' , linestyle='-.')
ax0.set_ylabel(' $\lambda^{E1}_i$ ')
ax0.legend(['$\lambda^{E1}_1$','$\lambda^{E1}_2$','$\lambda^{E1}_3$'], \
        loc='lower right')
ax0.fill_between(rho,0,4, color='k' , alpha=0.2)
ax0.axis([0., 2., -15., 4.])
ax0.autoscale(enable=True, axis='x', tight=True)
ax0.set_axisbelow(True)

ax1.plot(rhou,xequ[:,0] , 'C0' , linestyle='--' ), ax1.grid(True)
ax1.plot(rhou,xequ[:,1] , 'C1' )
ax1.plot(rhou,xequ[:,2] , 'C2' )
ax1.plot(rhos,xeqs      , 'C0' )
ax1.set_xlabel('$\\rho$'), ax1.set_ylabel('$\overline{x}$')
ax1.legend(['$x^{E1}$','$x^{E2}$','$x^{E3}$'])

ax1.fill_between(rhos,0,2, color='C3' , alpha=0.3)
ax1.fill_between(rhou,xequ[:,2],2, color='C3' , alpha=0.3)
ax1.fill_between(rhou,xequ[:,1],0, color='C3' , alpha=0.3)
ax1.fill_between(rhos,0,-2, color='C4' , alpha=0.3)
ax1.fill_between(rhou,xequ[:,1],-2, color='C4' , alpha=0.3)
ax1.fill_between(rhou,xequ[:,2], 0, color='C4' , alpha=0.3)
ax1.autoscale(enable=True, tight=True)
ax1.arrow(.75,-1.5,.0, 1., width=.01, head_length=.1, facecolor='k')
ax1.arrow(.75, 1.5,.0,-1., width=.01, head_length=.1, facecolor='k')
ax1.arrow(1.7,-0.4,.0,-.5, width=.01, head_length=.1, facecolor='k')
ax1.arrow(1.7, 0.4,.0, .5, width=.01, head_length=.1, facecolor='k')
ax1.set_axisbelow(True)

plt.show()
