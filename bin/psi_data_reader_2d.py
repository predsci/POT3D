#!/usr/bin/env python3
import numpy as np
import argparse
import psi_io as ps

def argParsing():

    parser = argparse.ArgumentParser(description='Read PSI 2D hdf5 data.')

    parser.add_argument("psi_2d_hdf5_file_name",
        help='Name of 2D PSI HDF5 file (e.g. br_photo.h5).')

    args = parser.parse_args()

    return parser.parse_args()

def main():

    ## Get input agruments:
    args = argParsing()

    xvec, yvec, data = ps.rdhdf_2d(args.psi_2d_hdf5_file_name)

    # If the data is in tp format, transpose to pt:
    # This makes plotting better/easier.
    if (np.max(yvec) > 3.5):
      tmpt = xvec
      tmpp = yvec
      yvec = tmpt
      xvec = tmpp
      data = np.transpose(data)

    # If the data is in sinlat, convert scale to theta:
    # No interpolation!  Rather, data is associated with nonuniform theta coordinates)
    if (np.min(yvec) < -0.5):
      yvec = np.arccos(yvec)

    tvec = yvec
    pvec = xvec

    tmin = np.min(tvec)
    tmax = np.max(tvec)
    pmin = np.min(pvec)
    pmax = np.max(pvec)

    NT = tvec.size
    NP = pvec.size

    print('Opened file:', args.psi_2d_hdf5_file_name)
    print('NT:         ', NT)
    print('NP:         ', NP)
    print('min(theta): ', np.min(tvec))
    print('max(theta): ', np.max(tvec))
    print('min(phi):   ', np.min(pvec))
    print('max(phi):   ', np.max(pvec))
    print('min(data):  ', np.min(data))
    print('max(data):  ', np.max(data))
    print('mean(data): ', np.mean(data))
    print('Example data point at data[5,5]:')
    print('theta', 'phi', 'value', sep='\t')
    print(tvec[5], pvec[5], data[5,5], sep='\t')


if __name__ == '__main__':
    main()
