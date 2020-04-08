# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#                                                                              #
# Hess-Smith method                                                            #
# for 2-dimensional steady incompressible irrotational flows around an airfoil #
#                                                                              #
#                                                                              #
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#                                                                              #
# The implementation exploits python dictionaries (DICT) and lists (LIST)      #
# They make the implementation easier to understand (or they should, at least) #
#                                                                              #
# Some suggestions for newbies:                                                #
# - python counts starting from 0, in an awkward way. Some examples:           #
#  a[0] is the first elem of the array a                                       #
#  a[0:4] are the first 4 (!!!!) elems of the array a                          #
#  if the index of the array exceed the dimension of the array, the count      #
#   restarts from zero, without any error (some "circular" indices)            #
# - lots of math operations come with numpy. To use them, you need to call     #
#  them with numpy (or its alias, e.g. import numpy as np, np is the alias).   #
#  E.g. numpy.cos(theta) for the cos() function, or np.cos() if the alias is   #
#  defined                                                                     #
# - Dictionaries are 'indexed' through keys. If a dictionary is defined as     #
#  dict_example = { 'name':'paolino', 'surname':'paperino', 'age':100, \       #
#                   'height':1.20 }                                            #
#  you can access its contents, using its KEYs. If you want to know the value  #
#  or the height field,                                                        #
#                                                                              #
#  paperino_height = dict_example['height']                                    #
#                                                                              #
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

import numpy as np                  # py library for math
import matplotlib.pyplot as plt     # py library for plot

import build_geometry as geo        # file with the function that builds the geom

deg2rad = np.pi/180.

#> === Free-stream conditions ===
#> Set Reynolds number of the analysis, Re1 = rho * U * c1 / mu, 
# for airfoils w/ unitary chord c1 = 1.:
# -> infinite value sets inviscid analysis
# -> finite values  turn on viscous corrections
Re1 = np.inf

#> Initialize freeStream DICT with input data for an inviscid analysis
freeStream = { 'p': 0., 'rho': 1., 'U': 1., 'al': 2.*deg2rad, \
               'viscous_correction':False }

#> Update freeStream DICT, with some useful quantities:
# free stream velocity vector
freeStream['Uvec']:np.array([ np.cos(freeStream['al']) , \
                              np.sin(freeStream['al'])]) *freeStream['U']

# viscosity, if viscous corrections are required (with finite value of Re1)
if ( Re1 != np.inf ):
  freeStream['viscous_corrections']:True
  freeStream['Re1']:Re1
  freeStream['dyn_visc']:freeStream['rho']*freeStream['U']/Re1
  freeStream['kin_visc']:freeStream['dyn_visc']/freeStream['rho']

#> === Airfoil input ===
#> Initialize bodies DICT, containing the list of the DICTs of the bodies
# of the simulation: in our lab, we will use just one body, whose name
# (KEY) is 'Airfoil_1'. At the initialization its value is void, []
bodies = { 'Airfoil_1':[] }
#> Initialize bodies['Airfoil_1'] DICT with airfoil data
bodies['Airfoil_1'] = { 'id':1, 'airfoil':'NACA0012', 'chord':1., \
                      'theta': 0.*deg2rad, 'xc_ref_point':0.25, \
                      'ref_point':np.array([ 0., 0. ]), \
                      'n_chord_pan':15 }

#> === Define the geometry of the body ( see build_geometry.py ) ===
elems, rr, ee, ii_te = geo.build_geometry( bodies['Airfoil_1'] )

#> Check output of build_geometry()
print(' np.shape(elems): '); print(np.shape(elems))
print(' np.shape(rr)   : '); print(np.shape(rr))
print(' np.shape(ee)   : '); print(np.shape(ee))
print(' np.shape(ee_te): '); print(np.shape(ee_te))
# print(' elems: ')
# for i in np.arange(len(elems)):
#   print(i, ': ', elems[i])
# print(' rr: '); print(rr)
# print(' ee: '); print(ee)

#> Check geometry
plt.plot(rr[0,:], rr[1,:], '-o')
for i in np.arange( np.size(elems,0) ):
  plt.quiver( elems[i]['cen'][0], elems[i]['cen'][1], \
              elems[i]['nor'][0], elems[i]['nor'][1], \
              color='black', width=0.005, headwidth=3 )
  plt.quiver( elems[i]['cen'][0], elems[i]['cen'][1], \
              elems[i]['tan'][0], elems[i]['tan'][1], \
              color='red', width=0.005, headwidth=3 )
plt.axis('equal')



#> === Build the linear system ===
# A, b, Au, Av = linsys.build_linsys( freeStream, elems, ii_te )

#> === Solve linear system ===
# ...

#> === Retrieve physical quantities ===
#> compontents of the velocity
# u = ...
# v = ...

#> pressure, pressure coefficients 
# ...

#> === Integral loads ===
#> aerodynamic coefficients
# ...


#> === Viscous correction ===
# ...













plt.show()
