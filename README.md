![POT3D](pot3d_logo.png)

# POT3D: High Performance Potential Field Solver #  
Predictive Science Inc.  
www.predsci.com

## OVERVIEW ##

`POT3D` is a Fortran code that computes potential field solutions 
to approximate the solar coronal magnetic field using observed 
photospheric magnetic fields as a boundary condition.  It is also 
used for computing Open Field and Potential Field Current Sheet (PFCS) models.  It has been (and continues to be) used for numerous studies of coronal
 structure and dynamics.  The code is highly parallelized using MPI and is GPU-accelerated using MPI+[OpenACC] (https://www.openacc.org/).  
 
`POT3D` is the potential field solver for the 
WSA model in the CORHEL software suite publicly hosted at 
the [Community Coordinated Modeling Center (CCMC)] (https://ccmc.gsfc.nasa.gov).  
A version of POT3D that additionaly includes GPU-offload with MPI+OpenMP was released as part of the 
Standard Performance Evaluation Corporation's (SPEC) 
beta version of the [SPEChpc(TM) 2021 
benchmark suites] (https://www.spec.org/hpc2021).  

Details of the `POT3D` code can be found in the following publications/presentations:  
 - *Variations in Finite Difference Potential Fields*.  
 Caplan, R.M., Downs, C., Linker, J.A., and Mikic, Z.  Submitted to Ap.J. (2021)  
 - *From MPI to MPI+OpenACC: Conversion of a legacy FORTRAN PCG solver for the spherical Laplace equation*.  
 Caplan, R.M., Mikic, Z., and Linker, J.L.  [arXiv:1709.01126] (https://arxiv.org/abs/1709.01126) (2017)

## HOW TO BUILD ##

Modify `build.sh` to set the proper `HDF5` library paths/flags 
and compiler flags for your system.  
Then, run `./build.sh`.  
See comments inside `build.sh` for more details.  

## HOW TO RUN SMALL TEST VALIDATION ##

After building the code, simply run `./validate.sh`.  
This will perform 2 runs of a small case using 1 and 2 MPI ranks respectively.

The runs are performed in `testsuite/validation/run/` and the second run overwrites the first.

Each result will be checked against a reference solution (in `/runs/validation/validation`) and  
a PASS/FAIL message will be displayed.

--------------------

## POT3D DOCUMENTATION ##

POT3D uses a namelist in an input file called `pot3d.dat` to set all parameters of a run.  
It also requires an input 2D data set in HDF5 format to use for the lower boundary condition.  
This file in the tests provided are called `br_input_<TEST>.h5` where <TEST> is the specific test name.

### RUNNING POT3D ###

To run `POT3D`, set the desired run parameters in an input  
text file named `pot3d.dat` (see below for parameter details),  
then copy or link the `pot3d` executable into the same directory as `pot3d.dat`  
and run the comand :  
  
    <MPI_LAUNCHER> -np <N> ./pot3d
  
where `<N>` is the total number of MPI ranks to use (typically CPU cores).  
and `<MPI_LAUNCHER>` is your MPI run command (`mpiexec`,`mpirun`, `ibrun`, `srun`, etc).  
For example:  `mpiexec -np 1024 ./pot3d`

**Important!**  
For CPU runs, make sure `ifprec=2` is set in the `pot3d.dat` input file.  
For GPU runs, make sure `ifprec=1` is set in the `pot3d.dat` input file.

### RUNNING ON GPUs ###

For standard cases, one should launch the code such that  
the number of MPI ranks per node is equal to the number of GPUs per node  
e.g.  
`mpiexec -np <N> --ntasks-per-node 4 ./pot3d`  
or  
`mpiexec -np <N> --npersocket 2 ./pot3d`   

To run efficiently, it is critical that `ifprec=1` is set in pot3d.dat.
      
### MEMORY REQUIREMENTS ###

To estimate how much memory (RAM) is needed for a run, compute:  
    `memory-needed = nr*nt*np*8*15/1024/1000/1000 GB`  
where `nr`, `nt`, and `np` are the chosen problem sizes in the `r`, `theta`, and `phi` dimension.

### RUN OUTPUT ###

Depending on the input parameters, `POT3D` can have various outputs.  
Typically, the three componeents of the potential magnetic field is output as `HDF5` files.  
In every run, the following two text files are ouptu as well:  

 - `pot3d.out`      An output log showing grid information and magnetic energy diagnostics.  
 - `timing.out`     Time profile information of the run.
      

### TESTSUITE and EXAMPLES ###

We have included a testsuite and examples.
      
We have included a validation script (`pot3d_validate.sh`) that takes two `pot3d.out` files and  
compares the magnetic energy values in them.  If they match, the test passes,  
if they do not, their values are displayed for further investigation.  

Example:  
`./pot3d\_validate.sh ./runs/test/run/pot3d.out ./runs/test/validation/pot3d.out`  

Within each test directory, there is a `validation` folder  
with the pot3d.out file for that run in it.

Information on the computer, compiler, and setup used to generate the  
validation data is in the `validation_run_information.txt` file in  
each validation folder.


### INPUT PARAMETER OPTIONS FOR POT3D.DAT ###
      
See the file "Documentation.txt: for a list of supported  
input parameter options able to be set in `pot3d.dat`.

      
      
