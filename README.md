![POT3D](pot3d_logo.png)
  
# POT3D: High Performance Potential Field Solver #
Predictive Science Inc.  
www.predsci.com
  
## OVERVIEW ##
  
`POT3D` is a Fortran code that computes potential field solutions to approximate the solar coronal magnetic field using observed photospheric magnetic fields as a boundary condition.  It can be used to generate potential field source surface (PFSS), potential field current sheet (PFCS), and open field (OF) models. It has been (and continues to be) used for numerous studies of coronal structure and dynamics.  The code is highly parallelized using [MPI](https://www.mpi-forum.org) and is GPU-accelerated using MPI+[OpenACC](https://www.openacc.org/), along with an option to use the [NVIDIA cuSparse library](https://developer.nvidia.com/cusparse). The [HDF5](https://www.hdfgroup.org/solutions/hdf5) file format is used for input/output.
  
`POT3D` is the potential field solver for the WSA model in the CORHEL software suite publicly hosted at the [Community Coordinated Modeling Center (CCMC)](https://ccmc.gsfc.nasa.gov/models/modelinfo.php?model=CORHEL/MAS/WSA/ENLIL).  
A version of `POT3D` that includes GPU-acceleration with both MPI+OpenACC and MPI+[OpenMP](https://www.openmp.org//) was released as part of the Standard Performance Evaluation Corporation's (SPEC) beta version of the [SPEChpc(TM) 2021 benchmark suites](https://www.spec.org/hpc2021).  
  
Details of the `POT3D` code can be found in the following publications:  
  
 - *Variations in Finite Difference Potential Fields*.  
 Caplan, R.M., Downs, C., Linker, J.A., and Mikic, Z.  [Ap.J. 915,1 44 (2021)](https://iopscience.iop.org/article/10.3847/1538-4357/abfd2f)
 - *From MPI to MPI+OpenACC: Conversion of a legacy FORTRAN PCG solver for the spherical Laplace equation*.  
 Caplan, R.M., Mikic, Z., and Linker, J.L.  [arXiv:1709.01126](https://arxiv.org/abs/1709.01126) (2017)
  
--------------------------------
  
## HOW TO BUILD POT3D ##
  
Copy the file `build.sh` to `my_build.sh`.  
Modify `my_build.sh` to set the `HDF5` library paths/flags and compiler flags compatible with your system environment.  
Then, run `./my_build.sh`.
  
See comments in `build.sh` for more details.
  
### Validate Installation ###
  
After building the code, you can test it is working by running `./validate.sh`.  
This will perform 2 runs of a small case using 1 and 2 MPI ranks respectively.
  
The runs are performed in `testsuite/validation/run/` and the second run overwrites the first.
  
Each result will be checked against a reference solution (in `/runs/validation/validation`) and a PASS/FAIL message will be displayed.
  
--------------------------------
  
## HOW TO USE POT3D ##
  
### Setting Input Options
  
POT3D uses a namelist in an input text file called `pot3d.dat` to set all parameters of a run.  See the provided `pot3d_input_documentation.txt` file for details on the various parameter options.  For any run, an input 2D data set in HDF5 format is required for the lower radial magnetic field (`Br`) boundary condition.  Examples of this file are contained in the `examples` and `testsuite` folders.
  
### Launching the Code ###
  
To run `POT3D`, set the desired run parameters in a `pot3d.dat` text file, then copy or link the `pot3d` executable into the same directory as `pot3d.dat`
and run the command:  
  `<MPI_LAUNCHER> -np <N> ./pot3d `  
where `<N>` is the total number of MPI ranks to use (typically equal to the number of CPU cores) and `<MPI_LAUNCHER>` is your MPI run command (e.g. `mpiexec`,`mpirun`, `ibrun`, `srun`, etc).  
For example:  `mpiexec -np 1024 ./pot3d`
  
**Important!**  
For CPU runs, set `ifprec=2` in the `pot3d.dat` input file.  
For GPU runs, set `ifprec=1` in the `pot3d.dat` input file, unless you build with the `cuSparse` library option, in which case you should set `ifprec=2`.
  
### Running POT3D on GPUs ###
  
For standard cases, one should launch the code such that the number of MPI ranks per node is equal to the number of GPUs per node  
e.g.  
`mpiexec -np <N> --ntasks-per-node 4 ./pot3d`  
or  
`mpiexec -np <N> --npersocket 2 ./pot3d`  

If the `cuSparse` library option was used to build the code, than set `ifprec=2` in `pot3d.dat`.  
If the `cuSparse` library option was NOT used to build the code, it is critical to set `ifprec=1` for efficient performance.
  
### Memory Requirements ###
  
To estimate how much memory (RAM) is needed for a run, compute:  
  
`memory-needed = nr*nt*np*8*15/1024/1000/1000 GB`  
  
where `nr`, `nt`, and `np` are the chosen problem sizes in the `r`, `theta`, and `phi` dimension.  
Note that this estimate is when using `ifprec=1`.  If using `ifprec=2`, the required memory is ~2x higher on the CPU, and even higher when using `cuSparse` on the GPU.
  
### Solution Output ###
  
Depending on the input parameters, `POT3D` can have various outputs. Typically, the three components of the potential magnetic field is output as `HDF5` files. In every run, the following two text files are output:

 - `pot3d.out`      An output log showing grid information and magnetic energy diagnostics.
 - `timing.out`     Time profile information of the run.
  
### Helpful Scripts ###
  
Some useful python scripts for reading and plotting the POT3D input data, and reading the output data can be found in the  `scripts` folder.  
  
-----------------------------
  
## EXAMPLES and TESTSUITE ##
  
### Examples ###
  
In the `examples` folder, we provide ready-to-run examples of three use cases of `POT3D` in the following folders:
  
1. **`/potential_field_source_surface`**  
A standard PFSS run with a source surface radii of 2.5 Rsun.
2. **`/potential_field_current_sheet`**  
A standard PFCS run using the outer boundary of the PFSS example as its inner boundary condition, with a domain that extends to 30 Rsun. The magnetic field solution produced is unsigned.  
3. **`/open_field`**  
An example of computing the "open field" model from the solar surface out to 30 Rsun using the same input surface Br as the PFSS example. The magnetic field solution produced is unsigned.  

### Testsuite ###

In the `testsuite` folder, we provide test cases of various sizes that can be used to validate and test the performance of `POT3D`.  
Each test case contains an `input` folder with the run input files, a `run` folder used to run the test, and a `validation` folder containing the output diagnotics used to validate the test, as well as a text file named `validation_run_information.txt`  containing information on how the validation run was computed (system, compiler, number of ranks, etc.) with performance details.  Note that all tests are set to use `ifprec=1` only.  An option to use `ifprec=2` will be added later.

To run a test, use the included script `run_test.sh` as:  
`run_test.sh <TEST> <NP>`  
where `<TEST>` is the test folder name and `<NP>` is the number of MPI ranks to use.  The test will run and then use the included script `scripts/pot3d_validate.sh` that takes two `pot3d.out` files and compares their magnetic energy values in order to validate the run results.  

The following is a list of the included tests, and their problem size and memory requirements:

1. **`validation`**  
Grid size:  63x91x225 = 1.28 million cells  
Memory (RAM) needed (using `ifprec=1`):  ~1 GB  
2. **`small`**  
Grid size:  133x361x901 = 43.26 million cells   
Memory (RAM) needed (using `ifprec=1`):  ~6 GB
3. **`medium`**  
Grid size:  267x721x1801 = 346.7 million cells  
Memory (RAM) needed (using `ifprec=1`):  ~41 GB  
4. **`large`**  
Grid size:  535x1441x3601 = 2.78 billion cells  
Memory (RAM) needed (using `ifprec=1`):  ~330 GB 

Note that these tests will *not* output the 3D magnetic field results of the run, so no extra disk space is needed.  
Instead, the validation is done with the magnetic energy diagnostics in the `pot3d.out` file.  




