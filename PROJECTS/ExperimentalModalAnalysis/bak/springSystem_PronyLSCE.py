# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#                                                                              #
# Use IRF and LSCE method to compute an approximation of the eigenvalues       #
# of the system                                                                #
# - IRF : impulsive respose function                                           #
# - LSCE: least-squares complex exponential (~Prony's method)                  #
#                                                                              #
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# 
# TODO:
# 1. compute coherence function to assess the validity of the identified modes
#       

import numpy as np
import matplotlib.pyplot as plt

# === Order of the model ===
# with a system of order = 80, the egienvalues reads
#  eig |   Numerical (Exact,Theoretical)   |     Experimental(Identified)   
# -----+-----------------------------------+-----------------------------------
#   0  |                                   |                                   
#   1  |        -0.00405 +-  2.8463j       |        -0.00405 +-  2.8463j       
#   2  |        -0.0304  +-  8.3082j       |        -0.03045 +-  8.3082j       
#   3  |        -0.0857  +- 13.0969j       |        -0.08579 +- 13.0969j       
#   4  |        -0.1415  +- 16.8245j       |        -0.14183 +- 16.8245j       
#   5  |        -0.1841  +- 19.1889j       |        -0.18472 +- 19.1893j       
#  ... |                ...                |                ...                

# === Parameters of the LSCE method ===
orderInit = 80 # 5  # order of the model 
orderMax  = 80      # order of the model
nt    = 100000      # n. of timesteps used n the ARMA linear system
                    # to compute the coefficients of the determinant
itInit = 1000       # initial time of the analysis 

# === Import files ===
frffile = np.load('./data/springFRF.npz')
print('\n reading FRF file: ',frffile.files)
f = frffile['f'] ; print('f.shape: ', f.shape) # frequency
U = frffile['U'] ; print('U.shape: ', U.shape) # F(input)
X = frffile['X'] ; print('X.shape: ', X.shape) # F(state) (x,v)
H = frffile['H'] ; print('H.shape: ', H.shape) # FRF (x,v)

irffile = np.load('./data/springIRF.npz')
print('\n reading IRF file: ',irffile.files)
t = irffile['t' ] ; print(' t.shape: ',  t.shape) # time
ht= irffile['ht'] ; print('ht.shape: ', ht.shape) # irf  

dt = t[1]-t[0] ; print(' dt: ', dt)

t  =  t[  itInit:]
ht = ht[:,itInit:]


Ndoftot = ht.shape[0]
t = t[itInit:]
Ndof = int( ht.shape[0]/2 )

StabilisationOrder = np.array( (0), dtype=int  )
StabilisationFreq  = np.array( (0), dtype=float)
for order in range(orderInit,orderMax+1,5):
    
    # === Assemble the ARMA overdetermined linear system ===
    print(' Ndof, nt, order: ', Ndof, nt, order)
    A = np.zeros((Ndof*nt,order), dtype=float)
    B = np.zeros((Ndof*nt), dtype=float)
    
    irow = 0
    for idof in range(Ndof):
        for it in range(nt):
            A[irow,:] = np.array(ht[idof,it:it+order])
            B[irow] = - ht[idof,it+order]
            irow += 1
    
    # === Solve the overdetermined linear system in LS sense ===
    print('\n === Solution of the LS linear system ===')
    print(' A.shape, B.shape: ', A.shape, B.shape)
    alpha, residuals, rank, s = np.linalg.lstsq(A, B, rcond=-1)
    
    print(' alpha: ', alpha)
    print(' residuals: ', residuals)
    print(' ')
    
    # Find the zeros of the polynomial p(z) 
    # with p(z) = z^order + alpha[0] z^(order-1) + alpha[1] z^(order-2) + ...
    #                   ... alpha[order-2] z + alpha[order-1]
    # print('np.array( 1 , alpha): ', np.append( 1 , alpha ) )
    print(' order: ', order)
    print(' alpha.size: ', alpha.shape )
    al = np.append( 1 , alpha[::-1] )
    print(' al.size: ', al.shape )
    z = np.roots( al )
    print(' z.shape: ', z.shape )
    print(' z: ', z )
    
    s = np.log(np.abs(z)) / dt + 1j * np.angle(z)/dt
    
    print(' s: ', s)
    
    StabilisationOrder = np.append( StabilisationOrder , \
                                    order * np.ones((order)) )
    StabilisationFreq  = np.append( StabilisationFreq  , \
                                    np.imag(s) )

# === Save the identified eigenvalues ===
np.savez('./data/springIdEig.npy', s=s)

# === Stabilisation diagram ===
plt.figure(101)
plt.plot(StabilisationFreq,StabilisationOrder,'o')

plt.show()



