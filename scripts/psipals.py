import numpy as np
import os
import matplotlib as mpl
import matplotlib.cm as cm
#from matplotlib.colors import ListedColormap
from matplotlib.colors import LinearSegmentedColormap
import sys

def load(N=1024):
    # Check for tools or corhel directory and set path root accordingly.
    ierr=0
    rootdir=os.environ.get('PS_TOOLS_HOME')
    if (rootdir is not None):
        rgbdir="color_palettes/rgb/"
        ierr=1
    else:
        rootdir=os.environ.get('CORHEL_HOME')
        if (rootdir is not None):
            rgbdir="/tools/ps_rsrc/rgb_color_palettes/rgb/"
            ierr=1
        else:
            rootdir=sys.path[0]+"/"
            rgbdir="psi_color_palettes/"
            ierr=1
    
    if(ierr==0):
        print("Error in psipals!")
    if (N<1):
        print("Error in psipals!  Must have N>=1")
        ierr=0

    if (ierr==1):
# Loop over dat files, extract name and make colormaps.
        for filename in os.listdir(os.path.join(rootdir, rgbdir)):
            if filename.endswith(".dat"):
                pal_dat_file=os.path.join(os.path.join(rootdir, rgbdir), filename)
                cmap_name="psi_"+os.path.splitext(filename)[0]
                pal = np.loadtxt(pal_dat_file)
            # Because of wierd IDL colormap format,
            # need to remove first and last entries:
            pal = np.delete(pal, (0), axis=0)
            pal = np.delete(pal, (len(pal)-1), axis=0)
            pal = pal/255.0
            #cmap_pal = ListedColormap(pal, name=cmap_name, N=None)
            cmap_pal = LinearSegmentedColormap.from_list(cmap_name, pal, N=N)
            mpl_version = mpl.__version__
            major,minor,patch = map(int,mpl_version.split('.')[:3])
            if major > 3 or (major == 3 and (minor >=4)):
                mpl.colormaps.register(cmap=cmap_pal,name=cmap_name)
            else:
                cm.register_cmap(name=cmap_name, cmap=cmap_pal)
    else:
        print("Warning!  PSIPALS did not load!")

