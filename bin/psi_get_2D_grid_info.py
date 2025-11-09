#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import psi_io as ps
import argparse

def argParsing():

    parser = argparse.ArgumentParser(description='grid_info: Show some grid info about a 2D file.')

    parser.add_argument('iFile',
                        help='Name of 2D file')

    return parser.parse_args()

def main():
    ## Get input agruments:
    args = argParsing()

    #Read the MAS file:
    pvec,tvec,data = ps.rdhdf_2d(args.iFile)

    Nt=tvec.size
    Np=pvec.size

    #Get stats on grid:
    tmin = np.min(tvec)
    tmax = np.max(tvec)
    pmin = np.min(pvec)
    pmax = np.max(pvec)

    dt_vec = np.diff(tvec)
    dp_vec = np.diff(pvec)

    ddt_vec = np.abs(np.diff(dt_vec))
    ddp_vec = np.abs(np.diff(dp_vec))

    per_dt_change_vec = (ddt_vec/dt_vec[1:Nt])
    per_dp_change_vec = (ddp_vec/dp_vec[1:Np])

    #Surface elements
    dt2d,dp2d = np.meshgrid(dt_vec,dp_vec)
    t2d,p2d = np.meshgrid(tvec,pvec)
    sint2d=np.sin(t2d[1:Np,1:Nt])
    vol = np.abs(sint2d*dt2d*dp2d)

    print()
    print(" File:  "+args.iFile)
    print()
    print(" T domain: [%5.2f, %5.2f]" % (tmin,tmax))
    print(" P domain: [%5.2f, %5.2f]" % (pmin,pmax))
    print()
    print(" Nt: %4d" % Nt)
    print(" Np: %4d" % Np)
    print()
    print("%10s %10s %10s %10s %10s %10s" % ('Quantity','Minimum','Maximum','Max/Min','Mean','STD'))
    print("-----------------------------------------------------------------")
    format_string = "%10s %10.2e %10.2e %10.2e %10.2e %10.2e"
    vec = dt_vec
    print(format_string % ('dt',np.min(vec),np.max(vec),np.max(vec)/np.min(vec),np.mean(vec),np.std(vec)))
    vec = dp_vec
    print(format_string % ('dp',np.min(vec),np.max(vec),np.max(vec)/np.min(vec),np.mean(vec),np.std(vec)))
    print()
    vec = vol
    print(format_string % ('dV',np.min(vec),np.max(vec),np.max(vec)/np.min(vec),np.mean(vec),np.std(vec)))
    print()
    format_string = "%10s %10.2e %10.2e %10s %10.2e %10.2e"
    vec = ddt_vec
    print(format_string % ('|ddt|',np.min(vec),np.max(vec),' ',np.mean(vec),np.std(vec)))
    vec = ddp_vec
    print(format_string % ('|ddp|',np.min(vec),np.max(vec),' ',np.mean(vec),np.std(vec)))
    print()
    format_string = "%10s %10.3f %10.3f %10s %10.3f %10.3f"
    vec = per_dt_change_vec
    print(format_string % ('|ddt|/dt',np.min(vec),np.max(vec),' ',np.mean(vec),np.std(vec)))
    vec = per_dp_change_vec
    print(format_string % ('|ddp|/dp',np.min(vec),np.max(vec),' ',np.mean(vec),np.std(vec)))
    print()


if __name__ == '__main__':
    main()
