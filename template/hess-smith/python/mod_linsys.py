
import sys
import numpy as np
from numpy import linalg as la

def build_linsys( freeStream, elems, ii_te ):
  
  #> Problem dimension
  nelems = len(elems)
  n_te = 1;   # hardcoded

  #> Initialization
  # matrices of the linsys
  A = np.zeros((nelems+n_te, nelems+n_te))
  b = np.zeros((nelems+n_te, 1))
  # auxiliary matrices
  Au = np.zeros((nelems, nelems+n_te))
  Av = np.zeros((nelems, nelems+n_te))

  for jj in np.arange(nelems):  # active element ( columns )
    #> b.c.
    for ii in np.arange(nelems):  # passive element ( rows )

      vsij = compute_velocity_source( elems[jj], elems[ii] )
      vvij = compute_velocity_vortex( elems[jj], elems[ii] )

      #> Matrix of the linear system
      A[ ii,jj]      = np.dot( elems[ii]['nor'], vsij )
      A[ ii,nelems] += np.dot( elems[ii]['nor'], vvij )

      #> RHS
      if ( jj == 0 ):
        b[ii] = -np.dot( elems[ii]['nor'], freeStream['Uvec'] )

      #> Auxiliary matrices for on-body velocity evaluation
      Au[ii,jj] = vsij[0] ;    Au[ii,nelems] += vvij[0]
      Av[ii,jj] = vsij[1] ;    Av[ii,nelems] += vvij[1]

      #> check numpy.isin
      if ( np.isin( elems[ii]['id'], ii_te ) ):
        #> kutta condition
        A[nelems,    jj] += np.dot( elems[ii]['tan'] , vsij )
        A[nelems,nelems] += np.dot( elems[ii]['tan'] , vvij )
        #> RHS
        if ( jj == 0 ):
          b[nelems] -= np.dot( elems[ii]['tan'], freeStream['Uvec'] )

  return A, b, Au, Av

def compute_velocity_source( el1, el2 ):
  # el1: active elem
  # el2: passive elem

  r1 = el2['cen'] - el1['ver1'] 
  r2 = el2['cen'] - el1['ver2']

  if ( el1['id'] == el2['id'] ):
    sij = 1.0; bij = np.pi
  else:
    sij = la.norm(r2)/la.norm(r1)
    vcross = r1[0]*r2[1] - r1[1]*r2[0]
    sinbij = vcross / ( la.norm(r1) * la.norm(r2) )
    cosbij = np.sum( r1*r2 ) / ( la.norm(r1) * la.norm(r2) )
    bij = np.arctan2(sinbij, cosbij)

  vstar = np.array(( -np.log(sij), bij)) / (2.*np.pi)
  
  #> rotation: from "local" coordinates to global coordinates
  Rot = np.array([ el1['tan'] , el1['nor'] ]).T
  v = np.matmul( Rot, vstar )

  return v

def compute_velocity_vortex( el1, el2 ):
  # el1: active elem
  # el2: passive elem

  r1 = el2['cen'] - el1['ver1'] # ; print(r1)
  r2 = el2['cen'] - el1['ver2'] # ; print(r2)

  if ( el1['id'] == el2['id'] ):
    sij = 1.0; bij = np.pi
  else:
    sij = la.norm(r2)/la.norm(r1)                            # ; print(sij)
    vcross = r1[0]*r2[1] - r1[1]*r2[0]                       # ; print(vcross)
    sinbij = vcross / ( la.norm(r1) * la.norm(r2) )          # ; print(sinbij)
    cosbij = np.sum( r1*r2 ) / ( la.norm(r1) * la.norm(r2) ) # ; print(cosbij)
    bij = np.arctan2(sinbij, cosbij)                         # ; print(bij)
    # stop

  vstar = np.array(( -bij, -np.log(sij) )) / (2.*np.pi)
  
  #> rotation: from "local" coordinates to global coordinates
  Rot = np.array([ el1['tan'] , el1['nor'] ]).T
  v = np.matmul( Rot, vstar )

  return v
