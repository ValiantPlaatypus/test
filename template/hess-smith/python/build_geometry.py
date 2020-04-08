
import sys
import numpy as np
import naca4digit as naca4

def build_geometry( body ):

  elems = []; ee_te = []

  if ( body['airfoil'][0:4] != 'NACA' ):
    sys.exit('Only NACA airfoil implemented so far.')
  if ( len(body['airfoil']) != 8 ):
    sys.exit('Only 4-digit NACA airfoil implemented so far.')

  M = int( body['airfoil'][4]   )
  P = int( body['airfoil'][5]   )
  SS= int( body['airfoil'][6:8] )
  c = body['chord']
  n = body['n_chord_pan']

  rr = naca4.naca4digit( M, P, SS, c, n, )

  n_p = np.size(rr, 1)   # n.of points
  ee  = np.array( [ np.arange( 0, n_p-1, dtype=int ), \
                    np.arange( 1, n_p  , dtype=int ) ] )
  n_e = np.size(ee, 1)
  ee_te = np.array( [ 0, n_e-1 ] )

  print(' n_e: ', n_e )
  elems = [];
  for ie in np.arange(n_e):
    e = {}
    e['id'] = ie
    e['ver1'] = rr[ :, ee[0,ie] ] ; print(e['ver1'])
    e['ver2'] = rr[ :, ee[1,ie] ] ; print(e['ver2'])
    e['cen']  = 0.5 * ( e['ver1'] + e['ver2'] )
    e['len']  = np.linalg.norm( e['ver1'] - e['ver2'] )
    e['tan']  = ( e['ver2'] - e['ver1'] ) / e['len']
    e['nor']  = ( [ -e['tan'][1], e['tan'][0] ] )
    elems.append(e);


  return elems, rr, ee, ee_te
