#
#
#

import numpy as np
import build_geometry as geo
import matplotlib.pyplot as plt
import mod_linsys as linsys

# debug flags
check_geometry = False

deg2rad = np.pi/180.

#> Set Reynolds number of the analysis, Re1 = rho * U * c1 / mu, 
# for airfoils w/ unitary chord c1 = 1.:
# - infinite value sets inviscid analysis
# - finite values turn on viscous corrections
Re1 = np.inf

#> Initialize freeStream dictionary with input data for inviscid
# computation
freeStream = { 'p': 0., 'rho': 1., 'U': 1., 'al': 2.*deg2rad, \
               'viscous_correction':False }

#> Update freeStream dictionary with useful quantities
freeStream['Uvec'] = np.array([ np.cos(freeStream['al']) , \
                                np.sin(freeStream['al'])]) *freeStream['U']

if ( Re1 != np.inf ):
  freeStream['viscous_corrections'] = True
  freeStream['Re1'] = Re1
  freeStream['dyn_visc'] = freeStream['rho']*freeStream['U']/Re1
  freeStream['kin_visc'] = freeStream['dyn_visc']/freeStream['rho']

# Airfoil input
bodies = { 'Airfoil_1':[] }
bodies['Airfoil_1'] = { 'id':1, 'airfoil':'NACA0012', 'chord':1., \
                      'theta': 0.*deg2rad, 'xc_ref_point':0.25, \
                      'ref_point':np.array([ 0., 0. ]), \
                      'n_chord_pan': 30 }

#> === Build geometry ===========================================
elems, rr, ee, ii_te = geo.build_geometry( bodies['Airfoil_1'] )

#> Check geometry ---
if ( check_geometry ):
  print(' ee   : '); print(ee)
  print(' ii_te: '); print(ii_te)
  print(' elems: ')
  for i in np.arange(len(elems)):
    print(i, ': ', elems[i])
  print(' rr: '); print(rr)
  plt.plot(rr[0,:], rr[1,:], '-o'); plt.axis('equal')

#> === Build linear system ======================================
A, b, Au, Av = linsys.build_linsys( freeStream, elems, ii_te )

#> === Solve linear system ======================================
sing = np.linalg.solve( A, b )
u = np.matmul( Au, sing ) + freeStream['Uvec'][0]
v = np.matmul( Av, sing ) + freeStream['Uvec'][1]
u = u[:,0]; v = v[:,0]

vt = np.zeros( len(u) ); vn = np.zeros( len(u) )
for i in np.arange(len(u)):
  vt[i] = np.dot( np.array([ u[i], v[i] ]), elems[i]['tan'] )
  vn[i] = np.dot( np.array([ u[i], v[i] ]), elems[i]['nor'] )

# if max(abs(vn)) > eps, something went wrong. b.c. not satisfied
vnorm = np.abs(vt)
cp = 1. - vnorm**2 / freeStream['U']

xc = np.zeros( len(u) )
for i in np.arange( len(u) ):
  xc[i] = elems[i]['cen'][0]

plt.plot(xc, cp); plt.axis('equal'); plt.grid(True)


plt.show()
