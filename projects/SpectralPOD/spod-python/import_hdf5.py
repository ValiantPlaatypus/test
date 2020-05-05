# Some details about h5py from:   docs.h5py.org/en/latest/quick.html
#
# hdf5 file is a container of twi kinds of objects:
# - <datasets>: array-like collection of data --> HDF5 dataset ~ NumPy arrays
# - <groups>  : folder-like containers of datasets and groups --> dictionaries
# <attributes> are metadata/description that can be stored along with the
# corresponding data.
#
# Import a hdf5 file, composed of only one group containing datasets

import h5py

def import_hdf5(filen,verbose):

    print('\n Importing hdf5 file:\n'+'   '+filen+'\n')

    f = h5py.File(filen,'r')

    if ( verbose == True ):

        print(' list(f.keys()):', list(f.keys()))
        
        for i1 in range(0,len(f)):
            print(' field', i1,'. key: ',list(f.keys())[i1], \
                  ' shape: ',f[list(f.keys())[i1]].shape)

    print('\n ...done.\n')

    return(f)    
