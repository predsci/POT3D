#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import psi_io as ps
import argparse

def argParsing():

    parser = argparse.ArgumentParser(description='grid_info: Show some grid info about a 3D file.')

    parser.add_argument('file',
                        help='Name of 3D file')

    parser.add_argument('-plot',
                        help='Make plots of grid.',
                        dest='plot',
                        action='store_true',
                        default=False,
                        required=False)

    parser.add_argument('-cm',
                        help='Show stats in cm units',
                        dest='cm',
                        action='store_true',
                        default=False,
                        required=False)

    parser.add_argument('-km',
                        help='Show stats in km units',
                        dest='km',
                        action='store_true',
                        default=False,
                        required=False)

    return parser.parse_args()

## Get input agruments:
args = argParsing()

#Read the file:
rvec,tvec,pvec,data = ps.rdhdf_3d(args.file)

if (args.cm):
    units = 6.96000000e+10
    unit_str="cm"
elif (args.km):
    units = 6.96000000e+5
    unit_str="km"
else:
    units = 1.0
    unit_str="Rs"

rvec = rvec*units

Nr=rvec.size
Nt=tvec.size
Np=pvec.size

#Get stats on grid:
rmin = np.min(rvec)
rmax = np.max(rvec)
tmin = np.min(tvec)
tmax = np.max(tvec)
pmin = np.min(pvec)
pmax = np.max(pvec)

dr_vec = np.diff(rvec)
dt_vec = np.diff(tvec)
dp_vec = np.diff(pvec)

ddr_vec = np.abs(np.diff(dr_vec))
ddt_vec = np.abs(np.diff(dt_vec))
ddp_vec = np.abs(np.diff(dp_vec))

per_dr_change_vec = (ddr_vec/dr_vec[1:Nr])
per_dt_change_vec = (ddt_vec/dt_vec[1:Nt])
per_dp_change_vec = (ddp_vec/dp_vec[1:Np])

#Volume elements
dr3d,dt3d,dp3d = np.meshgrid(dr_vec,dt_vec,dp_vec)
r3d,t3d,p3d = np.meshgrid(rvec,tvec,pvec)
sint3d=np.sin(t3d[1:Nt,1:Nr,1:Np])
r3d2=r3d[1:Nt,1:Nr,1:Np]*r3d[1:Nt,1:Nr,1:Np]
vol = np.abs(r3d2*sint3d*dr3d*dt3d*dp3d)

print(" ")
print(" MAS file:  %s" % args.file)
print(" ")
print(" R domain: [%5.2f, %5.2f] %s" % (rmin,rmax,unit_str))
print(" T domain: [%5.2f, %5.2f]" % (tmin,tmax))
print(" P domain: [%5.2f, %5.2f]" % (pmin,pmax))
print(" ")
print(" Nr: %4d" % Nr)
print(" Nt: %4d" % Nt)
print(" Np: %4d" % Np)
print(" ")
print("%10s %10s %10s %10s %10s %10s" % ('Quantity','Minimum','Maximum','Max/Min','Mean','STD'))
print("-----------------------------------------------------------------")
vec = dr_vec
print("%10s %10.2e %10.2e %10.2e %10.2e %10.2e" % ('dr',np.min(vec),np.max(vec),np.max(vec)/np.min(vec),np.mean(vec),np.std(vec)))
vec = dt_vec
print("%10s %10.2e %10.2e %10.2e %10.2e %10.2e" % ('dt',np.min(vec),np.max(vec),np.max(vec)/np.min(vec),np.mean(vec),np.std(vec)))
vec = dp_vec
print("%10s %10.2e %10.2e %10.2e %10.2e %10.2e" % ('dp',np.min(vec),np.max(vec),np.max(vec)/np.min(vec),np.mean(vec),np.std(vec)))
print(" ")
vec = vol
print("%10s %10.2e %10.2e %10.2e %10.2e %10.2e" % ('dV ('+unit_str+'^3)',np.min(vec),np.max(vec),np.max(vec)/np.min(vec),np.mean(vec),np.std(vec)))
volsun = (4.0/3.0)*np.pi*units*units*units
vec = vol/volsun
print("%10s %10.2e %10.2e %10.2e %10.2e %10.2e" % ('dV/Vsun',np.min(vec),np.max(vec),np.max(vec)/np.min(vec),np.mean(vec),np.std(vec)))
print(" ")
vec = ddr_vec
print("%10s %10.2e %10.2e %10s %10.2e %10.2e" % ('|ddr|',np.min(vec),np.max(vec),' ',np.mean(vec),np.std(vec)))
vec = ddt_vec
print("%10s %10.2e %10.2e %10s %10.2e %10.2e" % ('|ddt|',np.min(vec),np.max(vec),' ',np.mean(vec),np.std(vec)))
vec = ddp_vec
print("%10s %10.2e %10.2e %10s %10.2e %10.2e" % ('|ddp|',np.min(vec),np.max(vec),' ',np.mean(vec),np.std(vec)))
print(" ")
vec = per_dr_change_vec
print("%10s %10.3f %10.3f %10s %10.3f %10.3f" % ('|ddr|/dr',np.min(vec),np.max(vec),' ',np.mean(vec),np.std(vec)))
vec = per_dt_change_vec
print("%10s %10.3f %10.3f %10s %10.3f %10.3f" % ('|ddt|/dt',np.min(vec),np.max(vec),' ',np.mean(vec),np.std(vec)))
vec = per_dp_change_vec
print("%10s %10.3f %10.3f %10s %10.3f %10.3f" % ('|ddp|/dp',np.min(vec),np.max(vec),' ',np.mean(vec),np.std(vec)))
print(" ")

if args.plot:
    fsize=20
    lw=5.0

    print(" ")
    print("Generating plots...")

    fig=plt.figure(num=None, figsize=[14,7], dpi=96, facecolor='w')
    h=plt.plot(range(1,Nr+1),rvec/units,'r-',linewidth=lw)
    ax=fig.gca()
    ax.set_ylabel("$r$ (Rs)",{'fontsize': fsize})
    ax.set_xlabel("Grid index",{'fontsize': fsize})
    plt.tick_params(axis='both', labelsize=fsize, size=fsize)
    ax.grid()
    plt.xlim(xmin=1,xmax=Nr)
    plt.ylim(ymin=1,ymax=2.5)
    fig.savefig("rgrid_index.eps",bbox_inches='tight',dpi=150)
    plt.close('all')

    fig=plt.figure(num=None, figsize=[14,7], dpi=96, facecolor='w')
    h=plt.plot(range(1,Nt+1),tvec,'g-',linewidth=lw)
    ax=fig.gca()
    ax.set_ylabel("$\\theta$",{'fontsize': fsize})
    ax.set_xlabel("Grid index",{'fontsize': fsize})
    plt.yticks((0,np.pi/4,np.pi/2,3*np.pi/4,np.pi),('0','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'))
    plt.ylabel('Colatitude ($\\theta$)',{'fontsize': fsize})
    plt.tick_params(axis='both', labelsize=fsize, size=fsize)
    ax.grid()
    plt.xlim(xmin=1,xmax=Nt)
    plt.ylim(ymin=0,ymax=np.pi)
    fig.savefig("tgrid_index.eps",bbox_inches='tight',dpi=150)
    plt.close('all')

    fig=plt.figure(num=None, figsize=[14,7], dpi=96, facecolor='w')
    h=plt.plot(range(1,Np+1),pvec,'b-',linewidth=lw)
    ax=fig.gca()
    ax.set_ylabel("$\phi$",{'fontsize': fsize})
    ax.set_xlabel("Grid index",{'fontsize': fsize})
    plt.yticks((0,np.pi/2,np.pi,3*np.pi/2,2*np.pi), ('0','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'))
    plt.ylabel('$\phi$',{'fontsize': fsize})
    plt.tick_params(axis='both', labelsize=fsize, size=fsize)
    ax.grid()
    plt.xlim(xmin=1,xmax=Np)
    plt.ylim(ymin=0,ymax=2*np.pi)
    fig.savefig("pgrid_index.eps",bbox_inches='tight',dpi=150)
    plt.close('all')

    fig=plt.figure(num=None, figsize=[14,7], dpi=96, facecolor='w')
    h=plt.plot(rvec_at_dr/units,dr_vec,'r--',linewidth=lw)
    ax=fig.gca()
    ax.set_ylabel("$\Delta r$ (km)",{'fontsize': fsize})
    ax.set_xlabel("$r$ (Rs)",{'fontsize': fsize})
    plt.tick_params(axis='both', labelsize=fsize, size=fsize)
    ax.grid()
    plt.xlim(xmin=1,xmax=2.5)
    #plt.ylim(ymin=1,ymax=2.5)
    fig.savefig("rgrid_dr.eps",bbox_inches='tight',dpi=150)
    plt.close('all')

    fig=plt.figure(num=None, figsize=[14,7], dpi=96, facecolor='w')
    h=plt.plot(tvec_at_dt,dt_vec,'g--',linewidth=lw)
    ax=fig.gca()
    ax.set_ylabel("$\Delta \\theta$",{'fontsize': fsize})
    ax.set_xlabel("\\theta",{'fontsize': fsize})
    plt.xticks((0,np.pi/4,np.pi/2,3*np.pi/4,np.pi),('0','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'))
    plt.xlabel('Colatitude ($\\theta$)',{'fontsize': fsize})
    plt.tick_params(axis='both', labelsize=fsize, size=fsize)
    ax.grid()
    plt.xlim(xmin=0,xmax=np.pi)
    fig.savefig("tgrid_dt.eps",bbox_inches='tight',dpi=150)
    plt.close('all')

    fig=plt.figure(num=None, figsize=[14,7], dpi=96, facecolor='w')
    h=plt.plot(pvec_at_dp,dp_vec,'b--',linewidth=lw)
    ax=fig.gca()
    ax.set_ylabel("$\Delta \phi$",{'fontsize': fsize})
    ax.set_xlabel("\phi",{'fontsize': fsize})
    plt.xticks((0,np.pi/2,np.pi,3*np.pi/2,2*np.pi), ('0','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'))
    plt.xlabel('$\phi$',{'fontsize': fsize})
    plt.tick_params(axis='both', labelsize=fsize, size=fsize)
    ax.grid()
    plt.xlim(xmin=0,xmax=2*np.pi)
    fig.savefig("pgrid_dp.eps",bbox_inches='tight',dpi=150)
    plt.close('all')

    fig=plt.figure(num=None, figsize=[14,7], dpi=96, facecolor='w')
    plot_h=plt.pcolormesh(pvec,tvec,np.log(np.squeeze(vol[0,:,:])))
    plt.set_cmap("rainbow")
    plt.xticks((0,np.pi/2,np.pi,3*np.pi/2,2*np.pi), ('0','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'))
    plt.xlabel('$\phi$',{'fontsize': fsize})
    plt.yticks((0,np.pi/4,np.pi/2,3*np.pi/4,np.pi),('0','$\pi/4$','$\pi/2$','$3\pi/4$','$\pi$'))
    plt.ylabel('Colatitude ($\\theta$)',{'fontsize': fsize})
    plt.tick_params(axis='both', labelsize=fsize, size=fsize)
    cb=plt.colorbar(plot_h,fraction=0.024, pad=0.02, aspect=20)
    cb.set_label("log10 [dv ($km^3$)]",fontsize=fsize)
    cb.ax.yaxis.set_tick_params(labelsize=fsize)
    plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), size=fsize)
    ax.set_aspect('equal')
    plt.gca().invert_yaxis()
    fig.savefig("theta-phi-grid.eps",bbox_inches='tight',dpi=150)
    plt.close('all')

    fig=plt.figure(num=None, figsize=[14,7], dpi=96, facecolor='w')
    #Create 2D coordinates:
    P,R = np.meshgrid(pvec,rvec/units,indexing='ij')
    # Get cart coords:
    xvec_plot,yvec_plot = p2c(P,R)
    tidx=np.abs(tvec - np.pi/2.0).argmin()
    plot_h=plt.pcolormesh(xvec_plot,yvec_plot,np.transpose(np.log(np.squeeze(vol[:,tidx,:]))))
    ax=fig.gca()
    plt.set_cmap("rainbow")
    plt.tick_params(axis='both', labelsize=fsize, size=fsize)
    cb=plt.colorbar(plot_h,fraction=0.025, pad=0.02, aspect=38)
    cb.set_label("log10 [dv ($km^3$)]",fontsize=fsize)
    plt.xlabel('x (Rs)',{'fontsize': fsize})
    plt.ylabel('y (Rs)',{'fontsize': fsize})
    cb.ax.yaxis.set_tick_params(labelsize=fsize)
    plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), size=fsize)
    ax.set_aspect('equal')
    fig.savefig("r-phi-grid.eps",bbox_inches='tight',dpi=150)
    plt.close('all')

















