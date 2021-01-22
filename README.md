![POT3D](pot3d_logo.png)

# POT3D: High Performance Potential Field Solver #  
Predictive Science Inc.  
www.predsci.com

## OVERVIEW ##

POT3D is a Fortran code that computes potential field solutions 
to approximate the solar coronal magnetic field using observed 
photospheric magnetic fields as a boundary condition.  It is also 
used for computing Open Field and PFCS models.  
It has been (and continues to be) used for numerous studies of coronal
 structure and dynamics and is the potential field solver for the 
WSA model in the CORHEL software suite publicly hosted at 
NASA's Community Coordinated Modeling Center (CCMC)
(https://ccmc.gsfc.nasa.gov).  
POT3D has previously been released as part of the 
Standard Performance Evaluation Corporation's (SPEC) 
beta version of the SPEChpc(TM) 2021 
benchmark suites (https://www.spec.org/hpc2021).

## HOW TO BUILD ##

Modify "build.sh" to set the proper HDF5 library paths/flags   
and compiler flags for your system.  
Then, simply run "./build.sh".  
See comments inside build.sh for more details.

## HOW TO RUN SMALL TEST VALIDATION ##

After building the code, simply run `./validate.sh`.  
This will perform 2 runs of a small case using 1 and 2 MPI ranks respectively.

The runs are performed in `testsuite/validation/run/` and the second run overwrites the first.

Each result will be checked against a reference solution (in `/runs/validation/validation`) and  
a PASS/FAIL message will be displayed.

## POT3D DOCUMENTATION ##

POT3D uses a namelist in an input file called `pot3d.dat` to set all parameters of a run.  
It also requires an input 2D data set in HDF5 format to use for the lower boundary condition.  
This file in the tests provided are called `br_input_<TEST>.h5` where <TEST> is the specific test name.

### RUNNING POT3D ###

Simply run the code as:  
 > `mpiexec -np <NtotalRanks> ./pot3d`  
(or ibrun, mpirun, etc).

For CPU runs, please ensure `ifprec=2` is set in `pot3d.dat`  
in order to ensure optimal computation speed.

### RUNNING ON GPUs ###

For standard cases, one should launch the code such that  
the number of MPI ranks per node is equal to the number of GPUs per node  
e.g.  
> `mpiexec -np <NtotalRanks> --ntasks-per-node 4 ./pot3d`  
or  
> `mpiexec -np <NtotalRanks> --npersocket 2 ./pot3d`   

To run efficiently, it is critical that `ifprec=1` is set in pot3d.dat.

### RUN OUTPUT ###

At the end of the run, POT3D (as setup in the tests) outputs the following files:

 - `pot3d.out`      An output log showing grid information and derived magnetic energy.  
 - `timing.out`     Results of the MPI timing recorded in the code.

### VALIDATE ###

We have included a validation script (`pot3d_validate.sh`) that takes two `pot3d.out` files and  
compares the magnetic energy values in them.  If they match, the test passes,  
if they do not, their values are displayed for further investigation.  

Example:  
> `./pot3d\_validate.sh ./runs/test/run/pot3d.out ./runs/test/validation/pot3d.out`  

Within each test directory, there is a `validation` folder  
with the pot3d.out file for that run in it.

Information on the computer, compiler, and setup used to generate the  
validation data is in the `validation_run_information.txt` file in  
each validation folder.

### MEMORY FOOTPRINT ###

To estimate how many GB of RAM is needed for a run, compute:  
> `memory-needed = nr*nt*np*8*15/1024/1000/1000`  
