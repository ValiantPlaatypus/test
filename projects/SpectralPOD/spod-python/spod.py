import numpy as np

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Take a numpy array as an input:
# inp(:,:,...,:):
# the last index represents the time variable
# the first n-1 indices represent the multidimensional field
def spod(inp,nSnap1,nOverlap,dt):
    print('--spod()--------------')

    # Input analysis
    print('\n Input analysis ----- ')
    print('   len(inp.shape): ', len(inp.shape))
    print('       inp.shape : ',     inp.shape)
    n = len(inp.shape)-1  # n. of 'non-time' dimensions
    
    # Reshape the input as a 2-dimensional array:
    # 1st index for the unravelled 'non-time' dimensions
    # 2nd index for the time variable 
    inp2D = np.reshape(inp, (np.product(inp.shape[0:n]),inp.shape[n]) )
    print('     inp2D.shape : ',     inp2D.shape)

    # SPOD input analysis
    nDofs = inp2D.shape[0]
    nSnap = inp2D.shape[1]
    print(' Input: total number non-time dofs     : ', nDofs)
    print('        total number of snapshots      : ', nSnap)
    print('        number of snapshots per block  : ', nSnap1)
    print('        number of overlapping snapshots: ', nOverlap)

    # check input values
    if nSnap < nSnap1:
        print(' ERROR:')
        print(' n. of total snapshots lower than n. of spanshots per block.')
        return(-1)
    # ... other checks ...

    # find number of blocks
    nBlocks = int(np.floor( (nSnap-nSnap1)/(nSnap1-nOverlap) + 1 ))
    print(' -----> nBlocks: ', nBlocks)

    # Compute FFT of all the realisations --------------------------------------
    print(' Computing FFT of all the realisations ...')
    # initialise fft matrix for all the realisations:
    QQ = np.zeros((nDofs,nSnap1,nBlocks), dtype=complex)
    print(' QQ.shape: ', QQ.shape)
    for ib in range(0,nBlocks):
        print(' Realisation n.',ib+1,'/',nBlocks,end='\r')
        ind0 = ib*(nSnap1-nOverlap)        # Python starts from 0
        ind1 = ib*(nSnap1-nOverlap)+nSnap1 # Python counts the separators

        # fast fourier transform
        # remove mean value ( broadcasting )
        meanField = np.mean(inp2D[:,ind0:ind1],axis=1)
        q = np.zeros((nDofs,ind1-ind0))
        for i in range(ind0,ind1):
            q[:,i-ind0] = inp2D[:,i] - meanField 

        QQ[:,:,ib] = np.fft.fft(q)

    print(' ... done.')

    # Compute the SPOD modes for all the frequencies ---------------------------
    print(' Computing SPOD modes for all the frequencies ...')
    LLambda = np.zeros((nSnap1,nBlocks))  # <- must be real , dtype=complex)
    PPsi = np.zeros((nDofs,nSnap1,nBlocks), dtype=complex)
    for k in range(0,nSnap1):
        print('\r Frequency n.',k+1,'/',nSnap1,end='\r')

        # Build M matrix ( multiplication by dt/(s*nBlocks) missing )
        M = np.matmul( np.transpose(np.conjugate(QQ[:,k,:])) , QQ[:,k,:] ) \
                * dt/nBlocks
#       print(' M.shape: ', M.shape)

        # Compute eigenvalue decomposition M = Theta * Lambda * Theta'
        Lambda, Theta = np.linalg.eig(M)

        # Sort eigenvalues and eigenvectors ----
        isort = np.argsort(np.real(Lambda))   # ; print(isort)
        isort = isort[::-1]                   # ; print(isort)
        Lambda = Lambda[isort]
        Theta = Theta[:,isort]

        # Obtain and store SPOD fot frequency k
        Psi = np.matmul( QQ[:,k,:], Theta )
        Psi = np.matmul( Psi , np.diag(np.power(np.real(Lambda),-.5)) )

        PPsi[:,k,:] = Psi
        LLambda[k,:] = np.real(Lambda)

    print(' ... done.')

    print('--end of spod()-------')

    return(PPsi, LLambda)
