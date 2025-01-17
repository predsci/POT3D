## Steps to build POT3D on macOS

make sure we've `hdf5` installed as a dependency. You can do so via homebrew as:
```bash
$ brew install hdf5
```

it uses `FC=mpif90` in `build_examples/build_cpu_mpi-only_gcc_macos_homebrew_lf.sh`

```bash
$ git clone https://github.com/gxyd/POT3D
$ cd POT3D
$ git checkout -t origin/build_pot3d_unmodified
$ git checkout b547560acbb5e92de9f4220ce3a5eed5a79ac1b7
$ ./build_examples/build_cpu_mpi-only_gcc_macos_homebrew_lf.sh
$ export PATH=$PWD/bin:$PATH
```

## Running with a sample input

and then run:
```bash
$ ./validate.sh
```

the above command will produce an output like:
```console
Running POT3D with 1 MPI rank...
Done!
Wall clock time:                15.784941 seconds
 
Run has PASSED validation!
 
Running POT3D with 2 MPI ranks...
Done!
Wall clock time:                12.093250 seconds
 
Run has PASSED validation!
```
