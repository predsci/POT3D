#!/usr/bin/env python3
import h5py as h5
import numpy as np
import os
import sys
import signal
import argparse
import psi_io as ps

def signal_handler(signal, frame):
  print('You pressed Ctrl+C! Stopping!')
  sys.exit(0)

signal.signal(signal.SIGINT, signal_handler)

def s2c(r,t,p):
    x = r*np.sin(t)*np.cos(p)
    y = r*np.sin(t)*np.sin(p)
    z = r*np.cos(t)
    return (x,y,z)

def interp_3dcart_4pt(x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,f1,f2,f3,f4,xv,yv,zv):

    w1 = np.sqrt((xv-x1)*(xv-x1) + (yv-y1)*(yv-y1) + (zv-z1)*(zv-z1))
    w2 = np.sqrt((xv-x2)*(xv-x2) + (yv-y2)*(yv-y2) + (zv-z2)*(zv-z2))
    w3 = np.sqrt((xv-x3)*(xv-x3) + (yv-y3)*(yv-y3) + (zv-z3)*(zv-z3))
    w4 = np.sqrt((xv-x4)*(xv-x4) + (yv-y4)*(yv-y4) + (zv-z4)*(zv-z4))
 
    if (w1 == 0.0):
        w1=1.0
        w2=0.0
        w3=0.0
        w4=0.0
    elif (w2 == 0.0):
        w1=0.0
        w2=1.0
        w3=0.0
        w4=0.0
    elif (w3 == 0.0):
        w1=0.0
        w2=0.0
        w3=1.0
        w4=0.0
    elif (w4 == 0.0):
        w1=0.0
        w2=0.0
        w3=0.0
        w4=1.0
    else:
        w1=1.0/w1
        w2=1.0/w2
        w3=1.0/w3
        w4=1.0/w4

    w = w1 + w2 + w3 + w4

    value = (w1/w)*f1 + (w2/w)*f2 + (w3/w)*f3 + (w4/w)*f4

    return (value)


def argParsing():

    parser = argparse.ArgumentParser(description='Interpolate 2D spherical surface variables to a set of provided theta,phi points.')

    parser.add_argument("h5file",
        help='Name of 2D spherical surface h5 file.')

    parser.add_argument("-t",
        required=True,
        nargs='+',
        metavar=('theta0', 'theta1'),
        dest='theta',
        help="Space-separated theta (co-latitude) coordinate(s) (radians) to get interpolated values for.", 
        type=float, default=None)

    parser.add_argument("-p",
        nargs='+',
        required=True,
        metavar=('phi0', 'phi1'),
        dest='phi',
        help="Space-separated phi (longitude) coordinate(s) (radians) to get interpolated values for.", 
        type=float, default=None)

    args = parser.parse_args()

    return parser.parse_args()

def main():

    ## Get input agruments:
    args = argParsing()

    if (len(args.phi) != len(args.theta)):
       print("ERROR.  You must supply the same number of theta and phi points.")
       exit(1)

    xvec, yvec, data = ps.rdhdf_2d(args.h5file)

    # If the data is in tp format, transpose to pt:
    if (np.max(yvec) > 3.5):
      tmpt = xvec
      tmpp = yvec
      yvec = tmpt
      xvec = tmpp
      data = np.transpose(data)

    # If the data is in sinlat, convert scale to theta:
    if (np.min(yvec) < -0.5):
      yvec = np.arccos(yvec)

    tvec = yvec
    pvec = xvec

    tmin = np.min(tvec)
    tmax = np.max(tvec)
    pmin = np.min(pvec)
    pmax = np.max(pvec)

    Nt = tvec.size
    Np = pvec.size

    points = np.vstack((np.array(args.theta), np.array(args.phi))).T

    r=1.0
    eps=1e-7

    print('{}\t{}\t{}'.format('theta', 'phi', 'value'))

    for ii, point in enumerate(points):
        t = point[0]
        p = point[1]
        
        if (t>tmax+eps or t<tmin-eps):
            print('ERROR!  t is outside tvec!')
            exit(1)

        if (p>pmax+eps or p<pmin-eps):
            print('ERROR!  p is outside pvec!')
            exit(1)

        # Find indices that surround the point:
        i = np.searchsorted(tvec, t, side='left')
        j = np.searchsorted(pvec, p, side='left')

        im1=i-1
        jm1=j-1

        if (i==0):
            im1=0
        elif (i==Nt):
            i=Nt-1

        t1 = tvec[im1]
        t2 = tvec[i]
        
        if (j==0):
            jm1=0
        elif (j==Np):
            j=Np-1
            
        p1 = pvec[jm1]
        p2 = pvec[j] 

        x1,y1,z1 = s2c (r,t1,p1)
        f1 = data[im1,jm1]
        x2,y2,z2 = s2c (r,t2,p1)
        f2 = data[i,jm1]
        x3,y3,z3 = s2c (r,t1,p2)
        f3 = data[im1,j]
        x4,y4,z4 = s2c (r,t2,p2)
        f4 = data[i,j]
        
        xv,yv,zv = s2c (r,t,p)        
        
        value = interp_3dcart_4pt(x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,f1,f2,f3,f4,xv,yv,zv)

        print('{}\t{}\t{}'.format(t, p, value))


if __name__ == '__main__':
    main()


