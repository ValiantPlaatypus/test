import numpy as np

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
# Take a numpy array as an input:
# inp(:,:,...,:):
# the last index represents the time variable
# the first n-1 indices represent the multidimensional field
def pod_time(inp,dt):
    print('--pod_time()----------')
    
    # Input analysis
    print('\n Input analysis ----- ')
    print('   len(inp.shape): ', len(inp.shape))
    print('       inp.shape : ',     inp.shape)
    n = len(inp.shape)-1  # n. of 'non-time' dimensions
    
    # Reshape the input as a 2-dimensional array:
    # 1st index for the unravelled 'non-time' dimensions
    # 2nd index for the time variable 
    inp2D = np.reshape(inp, (np.product(inp.shape[0:n]),inp.shape[n]) )
    n2D = inp2D.shape[1]
#   print(' len(inp2D.shape): ', len(inp2D.shape))
    print('     inp2D.shape : ',     inp2D.shape)
    
    # Multiply snapshots by integration weights
    inp2D[:,1:n2D-1] = inp2D[:,1:n2D-1]*np.sqrt(dt)
    inp2D[:,0      ] = inp2D[:,0      ]*np.sqrt(dt*0.5)
    inp2D[:,  n2D-1] = inp2D[:,  n2D-1]*np.sqrt(dt*0.5)

    # SVD computation 
    print('\n SVD ---------------- ')
    print('   svd computation... ')
    u, s, vh = np.linalg.svd(inp2D, full_matrices=False)
    print('   ...done.')
#   print('  u.shape:',  u.shape)
#   print('  s.shape:',  s.shape)
#   print(' vh.shape:', vh.shape)

    # Reshape POD mode back to the original dimensions
    uPOD = np.reshape(u, inp.shape)
#   print('  uPOD.shape:',  uPOD.shape)

    print('--end of pod_time()---')
    return(uPOD, s, vh)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
