
import numpy as np
import matplotlib.pyplot as plt

sig = 10.
beta = 8./3.
rho = np.arange(0.,2.05, .05)

epsv = np.arange(-.5,.505,.002)
epsv1 = np.arange(-.5,.002,.002) ; epsv1[-1] = .0
epsv2 = np.arange( .0,.505,.002)

al_r = 0.0302
be_r = 0.003

req1 = np.zeros( (epsv1.shape[0],2), dtype=float )
req2 = np.zeros( (epsv2.shape[0],1), dtype=float )

req1[:,0] = .0
req1[:,1] = np.sqrt( -al_r * epsv1 / be_r )
req2[:,0] = .0

fig = plt.figure(1, figsize=plt.figaspect(0.75))
ax0 = fig.add_subplot(1,1,1)
ax0.plot(epsv1,req1[:,0] , 'C0' ), ax0.grid(True)
ax0.plot(epsv1,req1[:,1] , 'C1' , linestyle='--' , linewidth=2.)
ax0.plot(epsv2,req2[:,0] , 'C0' , linestyle='--' , linewidth=2.)
ax0.set_xlabel('$\\varepsilon$'), ax0.set_ylabel('$\overline{r}$')
ax0.legend(['$\overline{r}^{E2,3}_1$','$\overline{r}^{E2,3}_2$'])
ax0.axis([-.5, .5 , -.01 , 2.5])

ax0.fill_between(epsv1,0.0,req1[:,1], color='C3' , alpha=0.3)
ax0.fill_between(epsv1,2.5,req1[:,1], color='C4' , alpha=0.3)
ax0.fill_between(epsv2,2.5,req2[:,0], color='C4' , alpha=0.3)
# ax0.fill_between(rhou,xequ[:,1],0, color='C3' , alpha=0.3)
# ax0.fill_between(rhos,0,-2, color='C4' , alpha=0.3)
# ax0.fill_between(rhou,xequ[:,1],-2, color='C4' , alpha=0.3)
# ax0.fill_between(rhou,xequ[:,2], 0, color='C4' , alpha=0.3)
ax0.autoscale(enable=True, axis='x', tight=True)
ax0.arrow(-.30, 1.3,.0,-0.6, width=.005, head_length=.05, facecolor='k')
ax0.arrow( .10, 1.2,.0, 0.6, width=.005, head_length=.05, facecolor='k')
ax0.set_axisbelow(True)

plt.show()


