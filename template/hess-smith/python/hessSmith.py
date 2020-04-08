#
#
#

import numpy as np
import build_geometry as geo
import matplotlib.pyplot as plt

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
freeStream['Uvec']:np.array([ np.cos(freeStream['al']) , \
                              np.sin(freeStream['al'])]) *freeStream['U']
if ( Re1 != np.inf ):
  freeStream['viscous_corrections']:True
  freeStream['Re1']:Re1
  freeStream['dyn_visc']:freeStream['rho']*freeStream['U']/Re1
  freeStream['kin_visc']:freeStream['dyn_visc']/freeStream['rho']

# Airfoil input
bodies = { 'Airfoil_1':[] }
bodies['Airfoil_1'] = { 'id':1, 'airfoil':'NACA0012', 'chord':1., \
                      'theta': 0.*deg2rad, 'xc_ref_point':0.25, \
                      'ref_point':np.array([ 0., 0. ]), \
                      'n_chord_pan':30 }

elems, rr, ee, ee_te = geo.build_geometry( bodies['Airfoil_1'] )

print(' ee   : '); print(ee)
print(' ee_te: '); print(ee_te)
print(' elems: ')
for i in np.arange(len(elems)):
  print(i, ': ', elems[i])
# rr = naca4.naca4digit(0, 0, 12, 1., 30)
print(' rr: '); print(rr)
plt.plot(rr[0,:], rr[1,:], '-o')
plt.axis('equal')

plt.show()
