import numpy as np

# Signals to be used in Experimental Modal Analysis:
# - Impulsive or quasi-impulsive function
# - Frequency sweep
# - ...


# Quasi-impulsive functions
def impulsive(t, T,ampl, itype):
    if ( itype == 'gaussian' ):
        u = ampl * exp(-(ampl/(2.*T))**2)
    elif ( itype == 'cosine' ):
        if ( t < T ):
            u = ampl * .5 * ( np.cos( np.pi*t/T ) + 1. ) 
        else:
            u = 0.
    else:
        print(' impulsive(t, T,ampl, itype). ' + 
              ' itype must be ''gaussian'' or ''cosine'' ')
    return(u)

# Constant amplitude frequency sweep. Reference time is t=0
def sweep(t,T,f1,f2,S,ampl,stype):
    if ( stype == 'linear' ):
        u = ampl * linear_sweep(t, T,f1,f2,S) 
    elif ( stype == 'exponential'):
        u = ampl * exponential_sweep(t, T,f1,f2,S) 
    else:
        print(' sweep(t,f1,f2,S,ampl,stype). ' + 
              ' stype must be ''linear'' or ''exponential'' ')
    return(u)

def linear_sweep(t,T,f1,f2,S):
    u = np.sin( 2*np.pi * ( .5*S*t**2 + f1*t ) )
    return(u)

# To be checked -----
def exponential_sweep(t,T,f1,f2,S):
    # [S] octave/min
#   u = np.sin( 2*np.pi * ( 60.*f1/(S*np.log(2.)) * (2.**(S*t/60.)-1.) ) )
#   print(' S:', S)
#   print( 60.*f1/(S*np.log(2.)) * (2.**(S*t/60.)-1.) )
#   print( (2.**(S*t/60.)-1.) )
    u = np.sin( 2*np.pi*f1*T/np.log(f2/f1) * ( f2/f1 * np.exp(t/T) - 1. ) )
    return(u)
