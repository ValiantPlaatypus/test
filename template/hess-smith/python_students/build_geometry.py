# TODO:
# 1. build the required geometrical quantities
# 2. uncomment last lines to build the elems DICTionary

# Input:
# - body (DICTionary)
#   'id'             (integer): id. number of the airfoil (almost useless for just one airfoil)
#   'airfoil'        (string) : identifying the airfoil, e.g 'NACA0012'
#   'chord'          (real)   : chord
#   'theta'          (real)   : angle of attack [rad]
#   'xc_ref_point'   (real)   : position along the chord (in %) of the reference point, e.g. 0.25
#   'ref_point'[2]   (real)   : numpy array of size [2] containing x0,y0 coord. of the refPoint in 
#                               in the "global" reference system 
#   'n_chord_pan'    (integer): n. of surface elements for the discretisation of the upper and lower
#                               surface of the airfoil: the total n. of surface elements will be 2* 
# Output:
# - rr[2,npoints]  (real)   : numpy.array containing the x,y coords of the points.
#                             rr[0,i-1], rr[1,i-1] contain the x,y coord of the i-th point
# - ee[2,nelems]   (integer): numpy.array for the connectivity matrix.
#                             ee[0,i-1], ee[1,i-1] contain the id of the first and second node of the 
#                             i-th elem
# - ii_te[2]       (integer): vector containing the id of the elems at the trailing edge
# - elems[nelems]     (DICT): structure containing all the relevant informations of the elements.
#                             see last lines of the function.
#                             Uncomment and fill these lines, after having computed all the other
#                             outputs 

import sys
import numpy as np

def build_geometry( body ):

  #> === Some input checks ===
  if ( body['airfoil'][0:4] != 'NACA' ):
    sys.exit('Only NACA airfoil implemented so far.')
  if ( len(body['airfoil']) != 8 ):
    sys.exit('Only 4-digit NACA airfoil implemented so far.')

  #> === Read NACA 4-digit code ===
  M = int( body['airfoil'][4]   )
  P = int( body['airfoil'][5]   )
  SS= int( body['airfoil'][6:8] )
  c = body['chord']
  n = body['n_chord_pan']

  #> === Build the geometry ===
  # ... create rr, ee, ii_te numpy arrays
  # ...
  # ...
  # ...
  # ...
  # ...
  # ...

  #> === Create elems structure ===
  # print(' n_e: ', n_e )
  # elems = [];
  # for ie in np.arange(n_e):
  #   e = {}
  #   e['id'] = ie
  #   e['ver1'] = rr[ :, ee[0,ie] ] ; print(e['ver1'])
  #   e['ver2'] = rr[ :, ee[1,ie] ] ; print(e['ver2'])
  #   e['cen']  = 0.5 * ( e['ver1'] + e['ver2'] )
  #   e['len']  = np.linalg.norm( e['ver1'] - e['ver2'] )
  #   e['tan']  = ( e['ver2'] - e['ver1'] ) / e['len']
  #   e['nor']  = ( [ -e['tan'][1], e['tan'][0] ] )
  #   elems.append(e);


  return elems, rr, ee, ii_te
