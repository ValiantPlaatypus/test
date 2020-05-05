# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#                                                                              #
#                                                                              #
#                                                                              #
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

import numpy as np
# from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from experimental_inputs import *

# === Parameters ===
neigs = 8 

# === Define the srping system ===
# uniform spacing of the nodes of the undeformed system
nSpring = 5 
nPoints = nSpring+1
xx = np.arange(0.,nPoints)
yy = np.zeros(nPoints)
rr = np.array([xx,yy])
# Define input and output
i_r = 4   # output  index
i_s = 5   # forcing index

# uniform mass, damping and stiffness distributions
m1 = 1. ; c1 = 0.1  ; k1 = 100.
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
# SISO state space dynamical system
# A \in R[2*(nPoints-1)]^2
# B \in R[2*(nPoints-1),1]
# C \in R[1,2*(nPoints-1)]
# --- A: extended matrix
A = np.block([[np.zeros((len(i_n),len(i_n)),dtype=float) , np.eye(len(i_n))], \
              [-np.linalg.solve(Muu,Kuu) , \
               -np.linalg.solve(Muu,Cuu)]])

# --- B: forcing at the free end
B = np.zeros((2*(nPoints-1), 1), dtype=float)
B[i_s-len(i_d)+nSpring] = 1.

# --- C: position of the nodes
C = np.block( [ np.eye(nPoints-1) , \
                np.zeros((nPoints-1, nPoints-1), dtype=float) ] )
# C = np.block( [ np.zeros((nPoints-1, nPoints-1), dtype=float) , \
#                 np.eye(nPoints-1) ] )

C = C[i_r-len(i_d),:]

# check ----
print(' A: ') ; print(A)
print(' B: ') ; print(B)
print(' C: ') ; print(C)
# check ----

# === Emulation of an experimental test ===
# define Chirp (sine sweep) signal: linear chirp
f0 = 0.1 ; f1 = 20. ; T = 2000.
# k = ( f1 - f0 ) / T # linear sweep "chirpiness"
k = f1 / f0  # exponential sweep "chirpiness" .5 
# time vector
t0 = 0. ; tend = T ; dt = 0.010
tv = np.arange( t0, tend, dt )

# forcing
Timp = 0.1 ; ampl = 1. # ; ampl_noise = .01 
def u_inp(t, uA=ampl,f0=f0,f1=f1,T=T,k=k):
#   uA = .1 ; f0 = 1. ; f1 = 20. ; T = 1000. ; k = (f1-f0)/T
#   u = uA * np.sin( 2*np.pi * ( f0 + k * t ) * t )
#   u = sweep(t, T,f0,f1,k,uA,'linear')
#   u = impulsive(t, Timp,ampl, 'cosine')
    u = sweep(t, T,f0,f1,k,uA,'exponential')
#   print(' t, k, f0, f1, T, uA, u:', t, k, f0, f1, k, uA, u)
    return(u)
uv = np.zeros(tv.shape)
for it in range(tv.shape[0]):
    uv[it] = u_inp(tv[it], uA=ampl,f0=f0,f1=f1,T=T,k=k )

print('tv: ') ; print(tv)
print('uv: ') ; print(uv)
plt.figure(101)
plt.plot(tv,uv), plt.title('Input')
plt.grid(True), plt.xlabel('t'), plt.ylabel('u')

# Numerical integration
def model(x,t ,A,B,C,u):
#   print('t:',t)
#   print('u(t):',u(t))
#   print('A:') ; print(A)
#   print('B:') ; print(B)
#   print('C:') ; print(C)
#   print('x:') ; print(x)
#   print('np.matmul( A , x ):') ; print(np.matmul( A , x ))
#   print('np.matmul( A , x ).shape:') ; print(np.matmul( A , x ).shape)
    dxdt = np.matmul( A , x ) +  B[:,0] * u(t) 
#   print('np.matmul(A,x).shape: ') ; print(np.matmul(A,x).shape)
#   print('B.shape: ') ; print(B.shape)
#   print('B*u(t)') ; print(B*u(t))
#   print('dxdt.shape: ') ; print(dxdt.shape)
    return(dxdt)

x0 = np.zeros( (A.shape[0]) ,dtype=float).T
tspan = np.array([t0, tend])

print(' x0: ', x0)
print(' tspan: ', tspan)
x = odeint(model,x0,tv, args=(A,B,C, u_inp)) 

print(' x.shape: ', x.shape)

plt.figure(102)
plt.subplot(211)
plt.plot(tv,x[:,0:nSpring]), plt.title('State')
plt.grid(True), plt.xlabel('t'), plt.ylabel('x')
plt.subplot(212)
plt.plot(tv,x[:,nSpring:2*nSpring])
plt.grid(True), plt.xlabel('t'), plt.ylabel('v')

np.savez('./data/springTime.npz', t=tv, u=uv, x=x)



# # === Eigenproblem of A, extended matrix ===
print(' === Eigenproblem === ')
# - eigensolution sorted for decreasing module of the eigenvalues
vals, vecs = np.linalg.eig(A)
# - eigenvectors requires normalisation (amplitude and phase) for plots TODO
# check ----

print(' vals: ')
print(vals)
print('  f_0:' )
print(np.abs(vals)/(2.*np.pi))
print(' om_0:' )
print(np.abs(vals))
print(' csi :' )
print(-np.real(vals)/np.abs(vals))

# sort eigensolutions ----
isort = np.argsort(np.imag(vals))
vals = vals[isort]
vecs = vecs[:,isort]

plt.figure(501, figsize=(3,9))
for ip in range(5):
    plt.subplot(5,1,ip+1)
    plt.plot(np.real(vecs[0:nSpring,ip]))
    plt.plot(np.imag(vecs[0:nSpring,ip]))
    plt.title('s=%8.3f+j%8.3f' % (np.real(vals[ip]),np.imag(vals[ip])) )
    plt.grid(True) 


# # print(' vecs: ')
# # print(vecs)
# # check ----
# 
# plt.figure(101)
# plt.plot(np.real(vals),np.imag(vals),'o')
# plt.axis([-3, 3, -3, 3])
# plt.grid(True)
# 
# plt.figure(102)
# plt.plot(rr[0,1:],np.abs(vecs[0:nPoints-1,:]))
# 
# m = np.zeros(nPoints-1, dtype=complex)
# for ie in range(nPoints-1):
#     m[ie] = np.conjugate(vecs[0:nPoints-1,ie].T) @ vecs[0:nPoints-1,ie]
# 
# print(' m: ', m)

plt.show()
