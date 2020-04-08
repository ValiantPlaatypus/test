# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#                                                                              #
#                                                                              #
#                                                                              #
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

import numpy as np
import matplotlib.pyplot as plt

# === Parameters ===
neigs = 8 
# === Define the srping system ===
# uniform spacing of the nodes of the undeformed system
nSpring = 5 
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
# - Matrix (full format) initialisation
K = np.zeros((nPoints,nPoints),dtype=float)
C = np.zeros((nPoints,nPoints),dtype=float) # NumPy assigns by references
for ie in range(nSpring):
    K[[ie,ie,ie+1,ie+1],[ie,ie+1,ie,ie+1]] +=  np.array([k1,-k1,-k1,k1])
    C[[ie,ie,ie+1,ie+1],[ie,ie+1,ie,ie+1]] +=  np.array([c1,-c1,-c1,c1])
M = np.diag(np.ones(nPoints))

# check ----
# print(' K: ') ; print(  K   )
# print(' M: ') ; print(  M   )
# check ----

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

# Constrained dofs
i_d = [ 0 ]   # indices
u_d = [ 0.]   # values
i_n = np.setdiff1d( np.arange(0,nPoints), i_d )
print(' i_d: ', i_d)
print(' i_n: ', i_n)
Kuu = K[:,i_n] ; Kuu = Kuu[i_n,:]
Cuu = C[:,i_n] ; Cuu = Cuu[i_n,:]
Muu = M[:,i_n] ; Muu = Muu[i_n,:]

# check ----
print(' Kuu: ') ; print(Kuu)
print(' Cuu: ') ; print(Cuu)
print(' Muu: ') ; print(Muu)
# check ----

# === State-space representation of the dynamical system ===
# --- A: extended matrix
A = np.block([[np.zeros((len(i_n),len(i_n)),dtype=float) , np.eye(len(i_n))], \
              [-np.linalg.solve(Muu,Kuu) , \
               -np.linalg.solve(Muu,Cuu)]])

# --- B: forcing at the free end
B = np.zeros((2*(nPoints-1), 1), dtype=float)
B[nPoints-1] = 1.

# --- C: position of the nodes
C = np.block( [ np.eye(nPoints-1) , \
                np.zeros((nPoints-1, nPoints-1), dtype=float) ] )
# check ----
print(' A: ') ; print(A)
print(' B: ') ; print(B)
print(' C: ') ; print(C)
# check ----

# === Eigenproblem of A, extended matrix ===
# - eigensolution sorted for decreasing module of the eigenvalues
vals, vecs = np.linalg.eig(A)
# - eigenvectors requires normalisation (amplitude and phase) for plots TODO
# check ----
print(' vals: ')
print(vals)
# print(' vecs: ')
# print(vecs)
# check ----

plt.figure(101)
plt.plot(np.real(vals),np.imag(vals),'o')
plt.axis([-3, 3, -3, 3])
plt.grid(True)

plt.figure(102)
plt.plot(rr[0,1:],np.abs(vecs[0:nPoints-1,:]))

m = np.zeros(nPoints-1, dtype=complex)
for ie in range(nPoints-1):
    m[ie] = np.conjugate(vecs[0:nPoints-1,ie].T) @ vecs[0:nPoints-1,ie]

print(' m: ', m)

plt.show()

