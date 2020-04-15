
import sys
import numpy as np


def naca4digit(M, P, SS, c, n, spacing='cosine'):

  m  =  M / 100.;  p  =  P / 10.;  t  = SS / 100.;
  if ( m == 0 ): p = 1.

  xv = np.linspace(0., c, num=n+1)
  if ( spacing == 'half-cosine' ):
    xv =    c * ( 1. - np.cos(0.5*np.pi*xv/c) )
  elif ( spacing == 'cosine' ):
    xv = .5*c * ( 1. - np.cos(np.pi*xv/c) )
  else:
    sys.exit('So far, only spacing=''cosine'' implemented')

  #> Thickness
  tv = 5. * t * c * ( \
       0.2969 * (xv/c)**.5 - 0.1260 * (xv/c) - 0.3516 * (xv/c)**2 + \
       0.2843 * (xv/c)**3  - 0.1015 * (xv/c)**4 )

  #> Mean line and its derivative
  yv = np.zeros(n+1);  dv = np.zeros(n+1)
  for i in np.arange(n+1):
    if ( xv[i] <= p*c):
      yv[i] = c * ( m/p**2 * (xv[i]/c) * ( 2.*p - (xv[i]/c) ) )
      dv[i] = m/p**2 * 2.*( p - (xv[i]/c) )
    else:
      yv[i] = c * ( m/(1.-p)**2 * ( 1.+ ( 2.*p - xv[i]/c)) * (xv[i]/c) - 2.*p )
      dv[i] = m/(1.-p)**2 * 2.*( p - (xv[i]/c) )

  th = np.arctan2( dv, 1. )

  #> Coordinates of the lower and upper surfaces of the airfoil
  xU = xv - tv*np.sin(th);  xL = xv + tv*np.sin(th)
  yU = yv + tv*np.cos(th);  yL = yv - tv*np.cos(th)

  #> Sort the points on the airfoil surface
  x = np.append( xL[n:0:-1], xU )
  y = np.append( yL[n:0:-1], yU )

  rr = np.array([ x, y ])
  return rr

