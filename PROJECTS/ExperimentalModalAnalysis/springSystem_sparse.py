# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#                                                                              #
#                                                                              #
#                                                                              #
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

import numpy as np
import scipy.sparse as spm
from scipy.sparse import *
from scipy.sparse.linalg import inv , spsolve, eigs
import matplotlib.pyplot as plt

# === Parameters ===
neigs = 8 
# === Define the srping system ===
# uniform spacing of the nodes of the undeformed system
nSpring = 20
nPoints = nSpring+1
xx = np.arange(0.,nPoints)
yy = np.zeros(nPoints)
rr = np.array([xx,yy])
# uniform mass, damping and stiffness distributions
m1 = 1. ; c1 = .00 ; k1 = 1.
mass = np.ones(nPoints) * m1 
damp = np.ones(nSpring) * c1 
stif = np.ones(nSpring) * k1 

print('\n === Spring system details ===')
print(' rr: ')
print(rr)
print(' mass: ', mass)
print(' damp: ', damp)
print(' stif: ', stif)
print(' -----------------------------')

plt.figure(1, figsize=(7,3) )
plt.plot(rr[0,:],rr[1,:],'-o')
plt.grid(True) , plt.axis('equal')

# === Build mass, damping and stiffness matrices ===
# - in sparse format (take a look at SciPy documentation for sparse mat format)
# -> Strategy 1.
#    - Build a coo_sparse mat (~ Matlab, incremental def by indices and values)
#    - Convert to csr_ or csc_sparse for arithmetic operations and slicing
ke = np.zeros(0) ; ce = np.zeros(0)
i1 = np.zeros(0) ; j1 = np.zeros(0)
for ie in range(nSpring):
    i1 = np.append( i1 , [ie,ie,ie+1,ie+1] )
    j1 = np.append( j1 , [ie,ie+1,ie,ie+1] )
    ke = np.append( ke , [k1,-k1,-k1,k1] )
    ce = np.append( ce , [c1,-c1,-c1,c1] )

print(' i1: ', i1)
print(' j1: ', j1)

K = spm.coo_matrix( (ke,(i1,j1)) )
C = spm.coo_matrix( (ce,(i1,j1)) )

me = m1*np.ones(nPoints) ; im = np.arange(0,nPoints) ; jm = im
M = spm.coo_matrix( (me,(im,jm)) )

# Constrained system
#
#    |o----o----o----o-- ... --o----o
#
#  [ Kuu Kud ] d/ [ u   ] + [ Cuu Cud ] d/ [ u   ] + [ Kuu Kud ][ u   ] = [ fu ]
#  [ Kdu Kdd ]dt^2[ u_d ]   [ Cdu Cdd ] dt [ u_d ]   [ Kdu Kdd ][ u_d ]   [ fd ]
# u  : free dofs
# u_d: constrained dofs (known)
# --> 
# Kuu u'' + Cuu u' + Kuu u = fu - Kud u_d'' + Cud u_d' + Kud u_d

K = spm.csr_matrix(K)
C = spm.csr_matrix(C)
M = spm.csr_matrix(M)
# Constrained dofs
i_d = [ 0 ]   # indices
u_d = [ 0.]   # values
i_n = np.setdiff1d( np.arange(0,nPoints), i_d )
print(' i_d: ', i_d)
print(' i_n: ', i_n)
Kuu = K[:,i_n] ; Kuu = Kuu[i_n,:]
Cuu = C[:,i_n] ; Cuu = Cuu[i_n,:]
Muu = M[:,i_n] ; Muu = Muu[i_n,:]

print(' Kuu.toarray(): ')
print(Kuu)
# print(K)
print(Kuu.toarray())
print(spm.coo_matrix(Kuu))
print(Kuu.tocoo)
print(Kuu.shape)
Asp = bmat([[None, spm.eye(len(i_n))], \
            [-spm.linalg.spsolve(Muu,Kuu) , \
             -spm.linalg.spsolve(Muu,Cuu)]])
print(' Asp: ')
print(Asp.toarray())

# check ----
# print(spm.csr_matrix(K))
print(' C: ')
# print(C)
print(C.toarray())
print(' M: ')
# print(M)
print(M.toarray())
# check ----

# === Eigenproblem of A, extended matrix ===
vals, vecs = eigs(Asp,k=neigs,which='SM')
print(' vals: ')
print(vals)
print(' vecs: ')
print(vecs)

plt.figure(101)
plt.plot(np.real(vals),np.imag(vals),'o')
plt.axis([-3, 3, -3, 3])
plt.grid(True)

plt.figure(102)
plt.plot(rr[0,1:],np.imag(vecs[0:nPoints-1,:]))

print(' vecs.shape[] :' ,vecs.shape)
print('  Muu.shape[] :' , Muu.shape)
m = np.zeros(neigs, dtype=complex)
for ie in range(neigs):
    m[ie] = np.conjugate(vecs[0:nPoints-1,ie].T) @ vecs[0:nPoints-1,ie]

print(' m: ', m)

plt.show()

