import numpy as np
import h5py as h5

def rdh5(h5_filename):
    x = np.array([])
    y = np.array([])
    z = np.array([])
    f = np.array([])

    h5file = h5.File(h5_filename, 'r')
    f = h5file['Data']
    dims = f.shape
    ndims = np.ndim(f)

    #Get the scales if they exist:
    for i in range(0,ndims):
        if i == 0:
            if (len(h5file['Data'].dims[0].keys())!=0):
                x = h5file['Data'].dims[0][0]
        elif i == 1:
            if (len(h5file['Data'].dims[1].keys())!=0):
                y = h5file['Data'].dims[1][0]
        elif i == 2:
            if (len(h5file['Data'].dims[2].keys())!=0):
                z = h5file['Data'].dims[2][0]

    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    f = np.array(f)

    h5file.close()

    return (x,y,z,f)

def rdhdf(hdf_filename):

    x,y,z,f = rdh5(hdf_filename)
    return (x,y,z,f)


def rdhdf_1d(hdf_filename):

    x,y,z,f = rdhdf(hdf_filename)
    return (x,f)

def rdhdf_2d(hdf_filename):

    x,y,z,f = rdhdf(hdf_filename)
    return(x,y,f)

def rdhdf_3d(hdf_filename):

    x,y,z,f = rdhdf(hdf_filename)
    return(x,y,z,f)


def wrh5(h5_filename, x, y, z, f):

    h5file = h5.File(h5_filename, 'w')

    # Create the dataset (Data is the name used by the psi data)).
    h5file.create_dataset("Data", data=f)

    # Make sure the scales are desired by checking x type, which can
    # be None or None converted by np.asarray (have to trap seperately)
    if x is None: 
        x = np.array([], dtype=f.dtype)
        y = np.array([], dtype=f.dtype)
        z = np.array([], dtype=f.dtype)
    if x.any() == None:
        x = np.array([], dtype=f.dtype)
        y = np.array([], dtype=f.dtype)
        z = np.array([], dtype=f.dtype)

    # Make sure scales are the same precision as data.
    x=x.astype(f.dtype)
    y=y.astype(f.dtype)
    z=z.astype(f.dtype)

    #Get number of dimensions:
    ndims = np.ndim(f)

    #Set the scales:
    for i in range(0,ndims):
        if i == 0 and len(x) != 0:
            dim = h5file.create_dataset("dim1", data=x)
#            h5file['Data'].dims.create_scale(dim,'dim1')
            dim.make_scale('dim1')
            h5file['Data'].dims[0].attach_scale(dim)
            h5file['Data'].dims[0].label = 'dim1'
        if i == 1 and len(y) != 0:
            dim = h5file.create_dataset("dim2", data=y)
#            h5file['Data'].dims.create_scale(dim,'dim2')
            dim.make_scale('dim2')
            h5file['Data'].dims[1].attach_scale(dim)
            h5file['Data'].dims[1].label = 'dim2'
        elif i == 2 and len(z) != 0:
            dim = h5file.create_dataset("dim3", data=z)
#            h5file['Data'].dims.create_scale(dim,'dim3')
            dim.make_scale('dim3')
            h5file['Data'].dims[2].attach_scale(dim)
            h5file['Data'].dims[2].label = 'dim3'

    # Close the file:
    h5file.close()

def wrhdf(hdf_filename, x, y, z, f):

    wrh5(hdf_filename, x, y, z, f)

def wrhdf_1d(hdf_filename,x,f):

    x = np.asarray(x)
    y = np.array([])
    z = np.array([])
    f = np.asarray(f)
    wrhdf(hdf_filename,x,y,z,f)


def wrhdf_2d(hdf_filename,x,y,f):

    x = np.asarray(x)
    y = np.asarray(y)
    z = np.array([])
    f = np.asarray(f)
    wrhdf(hdf_filename,x,y,z,f)


def wrhdf_3d(hdf_filename,x,y,z,f):

    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    f = np.asarray(f)
    wrhdf(hdf_filename,x,y,z,f)

