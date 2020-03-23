#!/usr/bin/python

from numpy import *
from scipy.integrate import odeint

xv = arange( 0.0 , 10.0 , 0.01 )

def fBlasius(y,t):
   return [ y[1] , y[2] , -0.5*y[0]*y[2] ]

maxLoop=20

for i1 in arange(0,maxLoop) :
  print i1

tol = 1.0 * 10.0**(-9)
res = 100.0
nIter = 1

aold = 0.0
#a    = aold
#y0   = [ 0.0 , 0.0 , 0.332 ]
#y    = odeint( fBlasius , y0 , xv )
#derb = y[y[:,1].size,[1]]
#while ( res > tol and nIter < maxLoop) :
#  print nIter , maxLoop
#  nIter = nIter + 1
#  res = res / 10.0
  




a0 = 0.332
y0 = [ 0.0 , 0.0 , a0 ]

y = odeint( fBlasius, y0 , xv )

import matplotlib.pyplot as plt
plt.plot(y[:,0],xv,y[:,1],xv,y[:,2],xv)
plt.title('Blasius')
plt.legend(['g', 'g''','g"'])
#plt.axes([0 2 0 2])
plt.show()
