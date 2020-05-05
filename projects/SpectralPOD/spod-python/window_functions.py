# === Define windowing functions ===

import numpy as np

# https://ccrma.stanford.edu/~jos/sasp/Hamming_Window.html
def window(N,win_name='None'):
    if ( ( win_name == 'None' ) or ( win_name == 'rectangular' ) ):
        window = np.ones(N)
    elif ( win_name == 'hamming' ):
        window = .5 * ( 1. - np.cos( 2*np.pi * np.arange(0,N) / ( N - 1.) ) )
    elif ( win_name == 'hann' ):
        alpha = 25./46. ; beta = 1. - alpha
        window = alpha - beta * np.cos( 2*np.pi * np.arange(0,N) / ( N-1.) ) 
    else:
        print(' WARNING: the available window functions are: \n' + \
              ' - \'None\' \n' + \
              ' - \'rectangular\' \n' + \
              ' - \'hamming\' \n' + \
              ' - \'hann\' \n' + \
              ' No windowing applied.\n' +
              ' You could add your window function in window_functions.py\n')
        window = np.ones(N)

    return(window)

def hamming(N):
    alpha = 25./46. ; beta = 1.-alpha
    window = .5 * ( 1. - np.cos( 2*np.pi * np.arange(0,N) / ( N-1.) ) )
    return(window)

def hann(N):
    window = .5 * ( 1. - np.cos( 2*np.pi * np.arange(0,N) / ( N - 1.) ) )
    return(window)


