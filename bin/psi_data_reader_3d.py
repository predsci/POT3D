#!/usr/bin/env python3
import numpy as np
import argparse
import psi_io as ps

def argParsing():

    parser = argparse.ArgumentParser(description='Read PSI POT3D 3D hdf5 data.')

    parser.add_argument("psi_3D_hdf5_file_name",
        help='Name of 3D PSI POT3D HDF5 file (e.g. br.h5).')

    args = parser.parse_args()

    return parser.parse_args()

def main():

    ## Get input agruments:
    args = argParsing()

    rvec, tvec, pvec, data = ps.rdhdf_3d(args.psi_3D_hdf5_file_name)

    rmin = np.min(rvec)
    rmax = np.max(rvec)
    tmin = np.min(tvec)
    tmax = np.max(tvec)
    pmin = np.min(pvec)
    pmax = np.max(pvec)

    NR = rvec.size
    NT = tvec.size
    NP = pvec.size

    print('Opened file:', args.psi_3D_hdf5_file_name)
    print('NR:         ', NR)
    print('NT:         ', NT)
    print('NP:         ', NP)
    print('min(r):     ', np.min(rvec))
    print('max(r):     ', np.max(rvec))
    print('min(theta): ', np.min(tvec))
    print('max(theta): ', np.max(tvec))
    print('min(phi):   ', np.min(pvec))
    print('max(phi):   ', np.max(pvec))
    print('min(data):  ', np.min(data))
    print('max(data):  ', np.max(data))
    print('mean(data): ', np.mean(data))
    print('Example data point:')
    print('r[3]', 'theta[4]', 'phi[5]', 'data[5,4,3]', sep='\t')
    print(rvec[3], tvec[4], pvec[5], data[5,4,3], sep='\t')

if __name__ == '__main__':
    main()
