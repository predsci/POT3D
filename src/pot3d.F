c#######################################################################
c  _____   ____ _______ ____  _____
c |  __ \ / __ \__   __|___ \|  __ \
c | |__) | |  | | | |    __) | |  | |
c |  ___/| |  | | | |   |__ <| |  | |
c | |    | |__| | | |   ___) | |__| |
c |_|     \____/  |_|  |____/|_____/
c
c ****** POT3D: 3D potential magnetic field outside a sphere.
c
c ****** This program can find the classical potential field, the
c ****** fully open field, the source-surface field, and the
c ****** source-surface plus current-sheet field.
c
c        Authors:  Zoran Mikic
c                  Ronald M. Caplan
c                  Jon A. Linker
c                  Roberto Lionello
c                  Miko Stulajter
c
c        Predictive Science Inc.
c        www.predsci.com
c        San Diego, California, USA 92121
c
c#######################################################################
c Copyright 2021 Predictive Science Inc.
c
c Licensed under the Apache License, Version 2.0 (the "License");
c you may not use this file except in compliance with the License.
c You may obtain a copy of the License at
c
c    http://www.apache.org/licenses/LICENSE-2.0
c
c Unless required by applicable law or agreed to in writing, software
c distributed under the License is distributed on an "AS IS" BASIS,
c WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
c implied.
c See the License for the specific language governing permissions and
c limitations under the License.
c#######################################################################
c
c#######################################################################
c
      module ident
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c ****** Code name.
c-----------------------------------------------------------------------
c
      character(*), parameter :: idcode='POT3D'
      character(*), parameter :: vers  ='r4.0.0'
      character(*), parameter :: update='03/02/2022'
c
      end module
c#######################################################################
      module number_types_pc
c
c-----------------------------------------------------------------------
c
      use number_types
      use iso_fortran_env
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer, parameter :: r_typ_pc=REAL32
c
      end module
c#######################################################################
      module constants
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
      real(r_typ), parameter :: pi=3.1415926535897932_r_typ
c
      end module
c#######################################################################
      module global_dims
c
c-----------------------------------------------------------------------
c ****** Global number of mesh points.
c-----------------------------------------------------------------------
c
      implicit none
c
c ****** Global mesh size.
c
      integer :: nr_g,nrm1_g
      integer :: nt_g,ntm1_g
      integer :: np_g,npm1_g
c
c ****** Flag to indicate an axisymmetric run.
c
      logical :: axisymmetric=.false.
c
      end module
c#######################################################################
      module local_dims_r
c
c-----------------------------------------------------------------------
c ****** Local number of mesh points and indices in the r direction
c ****** for the radial mesh.
c-----------------------------------------------------------------------
c
      implicit none
c
c ****** Local mesh size.
c
      integer :: nr,nrm1
c
c ****** Dimensions of arrays on the "main" mesh.
c
      integer :: nrm
c
c ****** Indices of start and end points in the global mesh
c ****** belonging to this processor.
c
      integer :: i0_g,i1_g
c
c ****** Flags to indicate whether this processor has points
c ****** on the physical boundaries.
c
      logical :: rb0,rb1
c
      end module
c#######################################################################
      module local_dims_tp
c
c-----------------------------------------------------------------------
c ****** Local number of mesh points and indices in the theta
c ****** and phi dimensions.
c-----------------------------------------------------------------------
c
      implicit none
c
c ****** Local mesh size.
c
      integer :: nt,ntm1
      integer :: np,npm1
c
c ****** Dimensions of arrays on the "main" mesh.
c
      integer :: ntm
      integer :: npm
c
c ****** Indices of start and end points in the global mesh
c ****** belonging to this processor.
c
      integer :: j0_g,j1_g
      integer :: k0_g,k1_g
c
c ****** Flags to indicate whether this processor has points
c ****** on the physical boundaries.
c
      logical :: tb0,tb1
c
      end module
c#######################################################################
      module global_mesh
c
c-----------------------------------------------------------------------
c ****** Global mesh.
c-----------------------------------------------------------------------
c
      use number_types
      use constants
c
      implicit none
c
      real(r_typ), dimension(:), allocatable :: r_g,rh_g,dr_g,drh_g
      real(r_typ), dimension(:), allocatable :: t_g,th_g,dt_g,dth_g
      real(r_typ), dimension(:), allocatable :: p_g,ph_g,dp_g,dph_g
c
      real(r_typ), dimension(:), allocatable :: st_g,ct_g,sth_g,cth_g
      real(r_typ), dimension(:), allocatable :: sp_g,cp_g,sph_g,cph_g
c
c ****** Physical mesh size.
c
      real(r_typ), parameter, private :: one=1._r_typ
      real(r_typ), parameter, private :: two=2._r_typ
c
      real(r_typ) :: r0=1._r_typ
      real(r_typ) :: r1=30._r_typ
      real(r_typ), parameter :: t0=0.
      real(r_typ), parameter :: t1=pi
      real(r_typ), parameter :: p0=0.
      real(r_typ), parameter :: p1=two*pi
c
      real(r_typ), parameter :: pl=p1-p0
      real(r_typ), parameter :: pl_i=one/pl
c
      end module
c#######################################################################
      module local_mesh_r
c
c-----------------------------------------------------------------------
c ****** Local mesh for the r dimension.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
      real(r_typ), dimension(:), allocatable :: r,r2,rh,dr,drh
      real(r_typ) :: dr1,drn
c
c ****** Inverse quantities (for efficiency).
c
      real(r_typ), dimension(:), allocatable :: r_i,rh_i
      real(r_typ), dimension(:), allocatable :: dr_i,drh_i
c
      end module
c#######################################################################
      module local_mesh_tp
c
c-----------------------------------------------------------------------
c ****** Local mesh for the theta and phi dimensions.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
      real(r_typ), dimension(:), allocatable :: t,th,dt,dth
      real(r_typ), dimension(:), allocatable :: p,ph,dp,dph
c
      real(r_typ), dimension(:), allocatable :: st,ct,sth,cth
      real(r_typ), dimension(:), allocatable :: sp,cp,sph,cph
c
c ****** Inverse quantities (for efficiency).
c
      real(r_typ), dimension(:), allocatable :: dt_i,dth_i
      real(r_typ), dimension(:), allocatable :: st_i,sth_i
      real(r_typ), dimension(:), allocatable :: dp_i,dph_i
c
      end module
c#######################################################################
      module mpidefs
c
c-----------------------------------------------------------------------
c ****** MPI variables, processor topology, and processor information.
c-----------------------------------------------------------------------
c
      use mpi
c
      implicit none
c
c ****** Total number of processors.
c
      integer :: nproc
c
c ****** Total number of processors per node.
c
      integer :: nprocsh
c
c ****** Processor rank of this process in communicator
c ****** MPI_COMM_WORLD.
c
      integer :: iprocw
c
c ****** Processor rank of this process in communicator
c ****** comm_shared.
c
      integer :: iprocsh
c
c ****** Flag to designate that this is the processor with
c ****** rank 0 in communicator MPI_COMM_WORLD.
c
      logical :: iamp0
c
c ****** Communicator over all processors in the Cartesian topology.
c
      integer :: comm_all
c
c ****** Processor rank of this process in communicator
c ****** COMM_ALL.
c
      integer :: iproc
c
c ****** Processor rank in communicator COMM_ALL for the
c ****** processor that has rank 0 in MPI_COMM_WORLD.
c
      integer :: iproc0
c
c ****** Communicators over all processors in the phi dimension.
c
      integer :: comm_phi
c
c ****** Communicator over all shared processors on the node.
c
      integer :: comm_shared
c
c ****** Communicators over all processors in the r dimension.
c
      integer :: comm_r
c
c ****** Processor coordinate indices of this process
c ****** in the Cartesian topology.
c
      integer :: iproc_r,iproc_t,iproc_p
c
c ****** Processor coordinate indices of the neighboring
c ****** processors in the Cartesian topology.
c
      integer :: iproc_rm,iproc_rp
      integer :: iproc_tm,iproc_tp
      integer :: iproc_pm,iproc_pp
c
c ****** Number of processors along r, theta, and phi.
c
      integer :: nproc_r,nproc_t,nproc_p
c
c ****** Number type for REALs to be used in MPI calls.
c
      integer :: ntype_real
c
c ****** Total number of GPUs/node (DEPRICATED).
c
      integer :: gpn=0
c
      end module
c#######################################################################
      module decomposition_params
c
c-----------------------------------------------------------------------
c ****** Input parameters that define the domain decomposition
c ****** among processors.
c-----------------------------------------------------------------------
c
      implicit none
c
c ****** Number of processors per dimension.
c
      integer, dimension(3) :: nprocs=(/-1,-1,-1/)
c
      end module
c#######################################################################
      module decomposition
c
c-----------------------------------------------------------------------
c ****** Information defining the domain decomposition.
c-----------------------------------------------------------------------
c
      implicit none
c
c ****** Define the structure type for mapping local arrays
c ****** into global arrays.
c
      type :: map_struct
        integer :: n
        integer :: i0
        integer :: i1
        integer :: offset
      end type
c
c ****** Mapping structures for the different mesh types.
c
      type(map_struct), dimension(:), pointer :: map_rh
      type(map_struct), dimension(:), pointer :: map_rm
      type(map_struct), dimension(:), pointer :: map_th
      type(map_struct), dimension(:), pointer :: map_tm
      type(map_struct), dimension(:), pointer :: map_ph
      type(map_struct), dimension(:), pointer :: map_pm
c
      end module
c#######################################################################
      module meshdef
c
c-----------------------------------------------------------------------
c ****** Variables that define the mesh-point distributions.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
      integer, parameter :: nmseg=30
c
      real(r_typ), dimension(nmseg) :: drratio=0.
      real(r_typ), dimension(nmseg) :: dtratio=0.
      real(r_typ), dimension(nmseg) :: dpratio=0.
      real(r_typ), dimension(nmseg) :: rfrac=0.
      real(r_typ), dimension(nmseg) :: tfrac=0.
      real(r_typ), dimension(nmseg) :: pfrac=0.
c
      integer :: nfrmesh=0
      integer :: nftmesh=0
      integer :: nfpmesh=0
c
      real(r_typ) :: phishift=0.
c
      end module
c#######################################################################
      module fields
c
c-----------------------------------------------------------------------
c ****** Local field arrays.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
c ****** Potential solution array and cg temp array.
c
      real(r_typ), dimension(:,:,:), allocatable :: phi
      real(r_typ), dimension(:,:,:), allocatable :: x_ax
c
c ****** Boundary radial magnetic field arrays.
c
      real(r_typ), dimension(:,:), allocatable :: br0
      real(r_typ), dimension(:,:), allocatable :: br1
c
c ****** Arrays used in polar boundary conditions.
c
      real(r_typ), dimension(:), allocatable :: sum0,sum1
c
c ****** Arrays used for final magnetic field.
c
      real(r_typ), dimension(:,:,:), allocatable :: br,bt,bp
c
      end module
c#######################################################################
      module cgcom
c
      use number_types
c
      implicit none
c
c-----------------------------------------------------------------------
c ****** Number of equations to solve in the CG solver.
c-----------------------------------------------------------------------
c
      integer :: ncgeq
c
c-----------------------------------------------------------------------
c ****** CG field solver parameters.
c-----------------------------------------------------------------------
c
      integer :: ifprec=1
      integer :: ncgmax=5000
      integer :: ncghist=0
      integer :: ncgflush=25
      real(r_typ) :: epscg=1.e-9
c
c-----------------------------------------------------------------------
c ****** CG field solver variables.
c-----------------------------------------------------------------------
c
      integer :: ncg
      real(r_typ) :: epsn
c
      end module
c#######################################################################
      module vars
c
c-----------------------------------------------------------------------
c ****** Miscellaneous input variables.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
      character(256) :: outfile='pot3d.out'
      character(256) :: phifile='default'
      character(256) :: br0file='default'
      character(256) :: brfile='default'
      character(256) :: btfile='default'
      character(256) :: bpfile='default'
      character(256) :: br_photo_file='default'
      character(256) :: br_photo_original_file='default'
c
c ****** Type of field solution.
c ****** Select between 'potential', 'open', and 'source-surface'.
c
      character(16) :: option='potential'
c
c ****** Interval at which to dump diagonstics during the
c ****** iteration for the source-surface plus current-sheet
c ****** solution.
c
      integer :: ndump=0
c
c ****** Flag to skip the balancing of the flux (for PFSS and
c ****** OPEN field options only).

      logical :: do_not_balance_flux=.false.
c
c ****** Set format for output files.
c
      character(3) :: fmt='h5'
c
      logical :: hdf32=.true.
c
c ***** Validation run (tilted dipole).
c
      logical :: validation_run=.false.
c
      real(r_typ) :: dipole_angle=0.7853981633974483_r_typ
c
      end module
c#######################################################################
      module solve_params
c
c-----------------------------------------------------------------------
c ****** Parameters used in the solver.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
c ****** Boundary condition switch at r=R1.
c
      real(r_typ) :: pm_r1
c
      end module
c#######################################################################
      module timer
c
c-----------------------------------------------------------------------
c ****** Timer stack.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
      integer, parameter :: nstack=10
      integer :: istack=0
      real(r_typ), dimension(nstack) :: tstart=0.
c
      end module
c#######################################################################
      module timing
c
c-----------------------------------------------------------------------
c ****** Timing variables.
c-----------------------------------------------------------------------
c
      use number_types
c
      implicit none
c
      real(r_typ) :: t_startup=0.
      real(r_typ) :: t_solve=0.
      real(r_typ) :: t_pc_load=0.
      real(r_typ) :: t_pc=0.
      real(r_typ) :: t_ax=0.
      real(r_typ) :: t_io=0.
      real(r_typ) :: c_seam=0.
      real(r_typ) :: c_cgdot=0.
      real(r_typ) :: c_sumphi=0.
      real(r_typ) :: t_wall=0.
c
      end module
c#######################################################################
      module debug
c
c-----------------------------------------------------------------------
c ****** Debugging level.
c-----------------------------------------------------------------------
c
      implicit none
c
      integer :: idebug=0
c
      end module
c#######################################################################
      module assemble_array_interface
      interface
        subroutine assemble_array (map_r,map_t,map_p,a,a_g)
        use number_types
        use decomposition
        use mpidefs
        implicit none
        type(map_struct), dimension(0:nproc-1) :: map_r,map_t,map_p
        real(r_typ), dimension(:,:,:) :: a,a_g
        end subroutine
      end interface
      end module
c#######################################################################
      module cusparse_interface
c
      use, intrinsic :: iso_c_binding
c
#ifdef CUSPARSE
      interface
        subroutine load_lusol_cusparse(CSR_A,CSR_AI,CSR_AJ,N,M)
     &    BIND(C, name="load_lusol_cusparse")
          use, intrinsic :: iso_c_binding
          implicit none
          integer(C_INT), value :: N,M
          type(C_PTR), value :: CSR_A,CSR_AI,CSR_AJ
        end subroutine load_lusol_cusparse
c
        subroutine lusol_cusparse(x)
     &    BIND(C, name="lusol_cusparse")
          use, intrinsic :: iso_c_binding
          implicit none
          type(C_PTR), value :: x
        end subroutine lusol_cusparse
c
        subroutine unload_lusol_cusparse()
     &    BIND(C, name="unload_lusol_cusparse")
        end subroutine unload_lusol_cusparse
      end interface
#endif
c
      integer(c_int) :: cN,cM
c
      end module
c#######################################################################
      module matrix_storage_pot3d_solve
c
c-----------------------------------------------------------------------
c ****** Storage for the matrix/preconditioners of the solve.
c-----------------------------------------------------------------------
c
      use number_types
      use number_types_pc
c
      implicit none
c
      real(r_typ), dimension(:,:,:,:), allocatable :: a
      real(r_typ_pc), dimension(:), allocatable :: a_i
c
      integer, dimension(7) :: a_offsets

      integer :: N,M
      real(r_typ_pc), dimension(:), allocatable :: a_csr
      real(r_typ_pc), dimension(:), allocatable :: lu_csr
      real(r_typ_pc), dimension(:), allocatable :: a_csr_d
      integer, dimension(:), allocatable :: lu_csr_ja
      integer, dimension(:), allocatable :: a_csr_ia
      integer, dimension(:), allocatable :: a_csr_ja
      integer, dimension(:), allocatable :: a_N1
      integer, dimension(:), allocatable :: a_N2
      integer, dimension(:), allocatable :: a_csr_dptr
c
      end module
c#######################################################################
      program POT3D
c
c-----------------------------------------------------------------------
c
      use ident
      use mpidefs
      use vars
      use solve_params
      use timing
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: ierr
c
c-----------------------------------------------------------------------
c
c ****** Initialize MPI.
c
      call init_mpi
c
c ****** Start the wall-clock timer.
c
      call timer_on
c
c ****** Write the code name and version.
c
      if (iamp0) then
        write (*,*)
        write (*,*) idcode,' Version ',vers,', updated on ',update
      end if
c
      call timer_on
c
c ****** Read the input file.
c
      call read_input_file
c
c ****** Create the output file.
c
      if (iamp0) then
        call ffopen (9,outfile,'rw',ierr)
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in POT3D:'
          write (*,*) '### Could not create the output file.'
          write (*,*) 'File name: ',trim(outfile)
        end if
      end if
      call check_error_on_p0 (ierr)
c
c ****** Check the input parameters.
c
      call check_input
c
c ****** Check the processor topology.
c
      call check_proc_topology
c
c ****** Decompose the domain.
c
      call decompose_domain
c
c ****** Allocate global arrays.
c
      call allocate_global_arrays
c
c ****** Set the global meshes.
c
      call set_global_mesh
c
c ****** Decompose the mesh.
c
      call decompose_mesh_r
      call decompose_mesh_tp
c
c ****** Allocate local arrays.
c
      call allocate_local_arrays_r
      call allocate_local_arrays_tp
c
c ****** Set the local meshes.
c
      call set_local_mesh_r
      call set_local_mesh_tp
c
c ****** Print decomposition diagnostics.
c
      call decomp_diags
c
c ****** Initialize the flux and if validating, write analytic solution.
c
      if (validation_run) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### COMMENT from POT3D:'
          write (*,*) '### Validation run requested.'
          write (*,*) '### Ignoring br input file,'
          write (*,*) '### setting HDF32=.FALSE.,'
          write (*,*) '### and overriding output filenames.'
          write (9,*)
          write (9,*) '### COMMENT from POT3D:'
          write (9,*) '### Validation run requested.'
          write (9,*) '### Ignoring br input file'
          write (9,*) '### setting HDF32=.FALSE.,'
          write (9,*) '### and overriding output filenames.'
        end if
        hdf32=.false.
        brfile='br_solved.'//trim(fmt)
        btfile='bt_solved.'//trim(fmt)
        bpfile='bp_solved.'//trim(fmt)
        phifile='phi_solved.'//trim(fmt)
        call set_validation_flux
        call write_validation_solution
      else
        call set_flux
      end if
c
      call timer_off (t_startup)
c
c ****** Find the solution.
c
      if (iamp0) then
        write (*,*)
        write (*,*) '### COMMENT from POT3D:'
        write (*,*) '### Starting PCG solve.'
        write (9,*)
        write (9,*) '### COMMENT from POT3D:'
        write (9,*) '### Starting PCG solve.'
      end if
      call flush_output_file(outfile,9)
c
      call timer_on
c
      call potfld
c
      call timer_off (t_solve)
c
c ****** Compute B.
c
      call getb
c
c ****** Write solution to file.
c
      call write_solution
c
c ****** Magnetic energy diagnostics.
c
      call magnetic_energy
c
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      call timer_off (t_wall)
c
      call write_timing
c
      call endrun (.false.)
c
      end
c#######################################################################
      subroutine flush_output_file (filename,fsid)
c
c-----------------------------------------------------------------------
c
c ****** Flush output file by closing it and re-opening it in append.
c ****** filename is a character(256) which contains the filename
c ****** and fsid is the file's unit number.
c
c-----------------------------------------------------------------------
c
      use mpidefs
      use vars
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: fsid
      character(256) :: filename
c
c-----------------------------------------------------------------------
c
      integer :: ierr
c
c-----------------------------------------------------------------------
c
      if (iamp0) then
        close(fsid)
        call ffopen (fsid,trim(filename),'a',ierr)
      end if
      call check_error_on_p0 (ierr)
c
      end
c#######################################################################
      subroutine read_input_file
c
c-----------------------------------------------------------------------
c
c ****** Read the input file.
c
c-----------------------------------------------------------------------
c
      use global_dims
      use global_mesh
      use mpidefs
      use meshdef
      use cgcom
      use debug
      use vars
      use decomposition_params
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c ****** Values for the global mesh size.
c ****** Since these names conflict with those in LOCAL_DIMS*, it is
c ****** important not to use these modules here.
c
      integer :: nr=0
      integer :: nt=0
      integer :: np=0
c
c-----------------------------------------------------------------------
c
      namelist /topology/ nprocs,nr,nt,np,gpn
c
c-----------------------------------------------------------------------
c
      namelist /inputvars/ r0,r1,fmt,
     &                     drratio,dtratio,dpratio,
     &                     rfrac,tfrac,pfrac,
     &                     nfrmesh,nftmesh,nfpmesh,
     &                     phishift,
     &                     ifprec,ncgmax,ncghist,epscg,
     &                     idebug,br0file,phifile,
     &                     brfile,btfile,bpfile,br_photo_file,
     &                     br_photo_original_file,
     &                     option,
     &                     do_not_balance_flux,hdf32,
     &                     validation_run,dipole_angle,
     &                     ncgflush
c
c-----------------------------------------------------------------------
c
      integer :: ierr
      character(80) :: infile='pot3d.dat'
c
c-----------------------------------------------------------------------
c
c ****** Read the input file.
c
      call ffopen (8,infile,'r',ierr)
c
      if (ierr.ne.0) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in READ_INPUT_FILE:'
          write (*,*) '### Could not open the input file.'
          write (*,*) 'File name: ',trim(infile)
        end if
        call endrun (.true.)
      end if
c
      read (8,topology)
c
      read (8,inputvars)
c
      close (8)
c
      if (trim(fmt).ne.'h5') then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in READ_INPUT_FILE:'
          write (*,*) '### "fmt" must be "h5".'
          write (*,*) 'fmt: ',trim(fmt)
        end if
        call endrun (.true.)
      end if
c
c ****** Check if output names were overwritten.
c ****** If not, set default names with format fmt.
c
      if (trim(phifile).eq.'default') then
        phifile='phi.'//trim(fmt)
      end if
      if (trim(br0file).eq.'default') then
        br0file='br0.'//trim(fmt)
      end if
      if (trim(brfile).eq.'default') then
        brfile='br.'//trim(fmt)
      end if
      if (trim(btfile).eq.'default') then
        btfile='bt.'//trim(fmt)
      end if
      if (trim(bpfile).eq.'default') then
        bpfile='bp.'//trim(fmt)
      end if
      if (trim(br_photo_file).eq.'default') then
        br_photo_file='br_photo.'//trim(fmt)
      end if
      if (trim(br_photo_original_file).eq.'default') then
        br_photo_original_file='br_photo_original.'//trim(fmt)
      end if
c
      nr_g=nr
      nt_g=nt
      np_g=np
c
c ****** Check if the specified mesh dimensions are valid.
c
      call check_mesh_dimensions (nr_g,nt_g,np_g)
c
      nrm1_g=nr_g-1
      ntm1_g=nt_g-1
      npm1_g=np_g-1
c
      end subroutine
c#######################################################################
      subroutine check_error_on_p0 (ierr0)
c
c-----------------------------------------------------------------------
c
c ****** Check if the error flag IERR0 on processor 0 in
c ****** MPI_COMM_WORLD (i.e., processor IPROC0 in COMM_ALL)
c ****** indicates that the code should exit.
c
c ****** If IERR0 is non-zero, all the processors are directed
c ****** to call ENDRUN to terminate the code.
c
c-----------------------------------------------------------------------
c
      use mpidefs
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: ierr0
c
c-----------------------------------------------------------------------
c
c ****** MPI error return.
c
      integer :: ierr
c
c-----------------------------------------------------------------------
c
c ****** Broadcast IERR0 to all processors.
c
      call MPI_Bcast (ierr0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c
c ****** Call ENDRUN if IERR0 is non-zero.
c
      if (ierr0.ne.0) then
        call endrun (.true.)
      end if
c
      return
      end
c#######################################################################
      subroutine endrun (ifstop)
c
c-----------------------------------------------------------------------
c
c ****** End the run and exit the code.
c
c-----------------------------------------------------------------------
c
      use mpidefs
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      logical :: ifstop
c
c-----------------------------------------------------------------------
c
c ****** MPI error return.
c
      integer :: ierr
c
c-----------------------------------------------------------------------
c
c ****** Close the output file.
c
      if (iamp0) then
        close (9)
      end if
c
c ****** Exit MPI gracefully.
c
      call MPI_Finalize (ierr)
c
c ****** Call the STOP statement if requested.
c
      if (ifstop) then
        stop
      end if
c
      return
      end
c#######################################################################
      subroutine init_mpi
c
c-----------------------------------------------------------------------
c
c ****** Initialize MPI.
c
c-----------------------------------------------------------------------
c
      use number_types
      use mpidefs
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c ****** MPI error return.
c
      integer :: ierr,tcheck
c
c-----------------------------------------------------------------------
c
c ****** Real number to determine the KIND of REALs.
c
      real(r_typ) :: def_real
c
c-----------------------------------------------------------------------
c
      call MPI_Init_thread (MPI_THREAD_FUNNELED,tcheck,ierr)
c
c ****** Get the total number of processors.
c
      call MPI_Comm_size (MPI_COMM_WORLD,nproc,ierr)
c
c ****** Get the index (rank) of the local processor in
c ****** communicator MPI_COMM_WORLD in variable IPROCW.
c
      call MPI_Comm_rank (MPI_COMM_WORLD,iprocw,ierr)
c
c ****** Create a shared communicator for all ranks in the node.
c
      call MPI_Comm_split_type (MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,
     &                          MPI_INFO_NULL,comm_shared,ierr)
c
c ****** Get the total number of processors in node.
c
      call MPI_Comm_size (comm_shared,nprocsh,ierr)
c
c ****** Get the index (rank) of the local processor in the local node.
c
      call MPI_Comm_rank (comm_shared,iprocsh,ierr)
c
c ****** Set the flag to designate whether this processor
c ****** has rank 0 in communicator MPI_COMM_WORLD.
c
      if (iprocw.eq.0) then
        iamp0=.true.
      else
        iamp0=.false.
      end if
c
c ****** Set the number type for communicating REAL
c ****** numbers in MPI calls.
c
      if (kind(def_real).eq.KIND_REAL_4) then
        ntype_real=MPI_REAL4
      else if (kind(def_real).eq.KIND_REAL_8) then
        ntype_real=MPI_REAL8
      else
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in INIT_MPI:'
          write (*,*) '### Unrecognized default REAL number KIND:'
          write (*,*) 'KIND(default_real) = ',kind(def_real)
          write (*,*) 'This is a fatal error.'
        end if
        call endrun (.true.)
      end if
c
c ****** Set the GPU device number based on local rank.
c ****** NOTE! This assumes than #GPUs per node = #MPI ranks per node.
c
!$acc set device_num(iprocsh)
c
      end subroutine
c#######################################################################
      subroutine check_input
c
c-----------------------------------------------------------------------
c
c ****** Check the validity of the input parameters.
c
c-----------------------------------------------------------------------
c
      use number_types
      use vars
      use solve_params
      use mpidefs
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
c ****** Check that OPTION is valid and set the boundary condition
c ****** switches accordingly.
c
      if (option.eq.'potential') then
c
c ****** For a potential field, set d(phi)/dr to zero at r=R1
c ****** (i.e., the field is tangential to the boundary).
c
        pm_r1=one
c
      else if (option.eq.'open') then
c
c ****** For an open field, set phi to zero at r=R1
c ****** (i.e., the field is radial there).
c
        pm_r1=-one
c
      else if (option.eq.'ss') then
c
c ****** For a source surface field, set phi to zero at r=R1
c ****** (i.e., the field is radial there).
c
        pm_r1=-one
c
      else
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in CHECK_INPUT:'
          write (*,*) '### Invalid OPTION:'
          write (*,*)
          write (*,*) 'OPTION = ',trim(option)
          write (*,*)
          write (*,*) 'The options allowed are:'
          write (*,*) '''potential'''
          write (*,*) '''open'''
          write (*,*) '''ss'''
        end if
        call endrun (.true.)
      end if
c
      if (iamp0) then
        write (*,*)
        write (*,*) '### COMMENT from CHECK_INPUT:'
        write (*,*) '### Field solve type:'
        write (*,*)
        write (*,*) 'OPTION = ',option
        write (9,*)
        write (9,*) '### COMMENT from CHECK_INPUT:'
        write (9,*) '### Field solve type:'
        write (9,*)
        write (9,*) 'OPTION = ',option
      end if
c
      return
      end
c#######################################################################
      subroutine set_proc_topology
c
c-----------------------------------------------------------------------
c
c ****** Set the optimal values of the MPI rank topology
c ****** in dimensions not set by user.
c
c-----------------------------------------------------------------------
c
      use mpidefs
      use decomposition_params
      use number_types
      use global_dims
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1.0_r_typ
      real(r_typ), parameter :: zero=0.0_r_typ
      real(r_typ), parameter :: bigval=HUGE(1.0_r_typ)
c
c-----------------------------------------------------------------------
c
      integer, dimension(:), allocatable :: factors
      integer, dimension(:,:), allocatable :: rank_factors
      real(r_typ), dimension(:,:), allocatable :: nperrank
      real(r_typ), dimension(:), allocatable :: penalty
c
      integer :: i,j,k,num_fac,num_rank_fac,best_idx
      integer :: fr=1,ft=1,fp=1
      real(r_typ) :: a12,a13,a23
c
c-----------------------------------------------------------------------
c
c ****** Extract nproc values.  A value of -1 indicates the dimension
c ****** should be autoset.
c
      nproc_r=nprocs(1)
      nproc_t=nprocs(2)
      nproc_p=nprocs(3)
c
c ****** If no dimensions are to be autoset, return.
c
      if(nproc_r.ne.-1.and.nproc_t.ne.-1.and.nproc_p.ne.-1) return
c
c ****** Get all factors of nproc and store them in factors array.
c
      i=1
      num_fac=0
      do while(i.le.nproc)
        if (MOD(nproc,i).eq.0) then
          num_fac=num_fac+1
        endif
        i=i+1
      enddo
      allocate (factors(num_fac))
      i=1
      num_fac=0
      do while(i.le.nproc)
        if (MOD(nproc,i).eq.0) then
          num_fac=num_fac+1
          factors(num_fac)=i
        endif
        i=i+1
      enddo
c
c ****** Set penalty function parameters and any fixed dimensions
c ****** based on which dimensions are to be autoset.
c
      a12=one
      a13=one
      a23=one
c
      if (nproc_r.ne.-1) then
        fr=nproc_r
        a12=zero
        a13=zero
      end if
      if (nproc_t.ne.-1) then
        ft=nproc_t
        a12=zero
        a23=zero
      end if
      if (nproc_p.ne.-1) then
        fp=nproc_p
        a13=zero
        a23=zero
      end if
c
c ****** Loop over all combinations of factors and save those that
c ****** yield the correct number of MPI ranks into rank_factors array.
c
      num_rank_fac=0
      do k=1,num_fac
        do j=1,num_fac
          do i=1,num_fac
            if(nproc_r.eq.-1) fr=factors(i)
            if(nproc_t.eq.-1) ft=factors(j)
            if(nproc_p.eq.-1) fp=factors(k)
            if (fr*ft*fp.eq.nproc) then
              num_rank_fac=num_rank_fac+1
            end if
          enddo
        enddo
      enddo
c
      if (num_rank_fac.eq.0) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in SET_PROC_TOPOLOGY:'
          write (*,*) '### Processor topology specification error.'
          write (*,*) 'No valid topologies found for selected options.'
          write (*,*) 'Number of MPI ranks = ',nproc
          write (*,*) 'NPROC_R = ',nproc_r
          write (*,*) 'NPROC_T = ',nproc_t
          write (*,*) 'NPROC_P = ',nproc_p
        end if
        call endrun (.true.)
      end if
c
      allocate(rank_factors(num_rank_fac,3))
      allocate(nperrank(num_rank_fac,3))
      allocate(penalty(num_rank_fac))
c
      rank_factors(:,:)=-1
      penalty(:)=bigval
c
      num_rank_fac=0
      do k=1,num_fac
        do j=1,num_fac
          do i=1,num_fac
            if(nproc_r.eq.-1) fr=factors(i)
            if(nproc_t.eq.-1) ft=factors(j)
            if(nproc_p.eq.-1) fp=factors(k)
            if (fr*ft*fp.eq.nproc) then
              num_rank_fac=num_rank_fac+1
              rank_factors(num_rank_fac,1)=fr
              rank_factors(num_rank_fac,2)=ft
              rank_factors(num_rank_fac,3)=fp
            end if
          enddo
        enddo
      enddo
c
c ****** Get number of grid points per rank for each dimension.
c
      nperrank(:,1)=real(nr_g)/rank_factors(:,1)
      nperrank(:,2)=real(nt_g)/rank_factors(:,2)
      nperrank(:,3)=real(np_g)/rank_factors(:,3)
c
c ****** Compute penalty function.
c
      penalty(:)=a12*(nperrank(:,1)-nperrank(:,2))**2
     &          +a23*(nperrank(:,2)-nperrank(:,3))**2
     &          +a13*(nperrank(:,3)-nperrank(:,1))**2
c
c ****** Eliminate any choices that yield less than a minimum number
c ****** of grid points per rank.
c
      do i=1,num_rank_fac
        if (nperrank(i,1).lt.4) penalty(i)=bigval
        if (nperrank(i,2).lt.4) penalty(i)=bigval
        if (nperrank(i,3).lt.3) penalty(i)=bigval
      enddo
c
c ****** Find optimal topology.
c
      best_idx=MINLOC(penalty,1)
c
      if (penalty(best_idx).eq.bigval) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in SET_PROC_TOPOLOGY:'
          write (*,*) '### Processor topology specification error.'
          write (*,*) 'No valid topologies found for selected options'
          write (*,*) 'with selected grid.  '
          write (*,*) 'It is likely you are using too many MPI ranks.'
          write (*,*) 'Number of MPI ranks = ',nproc
          write (*,*) 'NPROC_R = ',nproc_r
          write (*,*) 'NPROC_T = ',nproc_t
          write (*,*) 'NPROC_P = ',nproc_p
          write (*,*) 'NR = ',nr_g
          write (*,*) 'NT = ',nt_g
          write (*,*) 'NP = ',np_g
        end if
        call endrun (.true.)
      end if
c
c ****** Set optimal topology.
c
      nprocs(1)=rank_factors(best_idx,1)
      nprocs(2)=rank_factors(best_idx,2)
      nprocs(3)=rank_factors(best_idx,3)
c
      deallocate(factors)
      deallocate(rank_factors)
      deallocate(nperrank)
      deallocate(penalty)
c
      end subroutine
c#######################################################################
      subroutine check_proc_topology
c
c-----------------------------------------------------------------------
c
c ****** Check the validity of the requested processor topology.
c
c-----------------------------------------------------------------------
c
      use mpidefs
      use decomposition_params
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: i,nreq
c
c-----------------------------------------------------------------------
c
c ****** Check the processor topology.
c
      do i=1,3
        if (nprocs(i).lt.1.and.nprocs(i).ne.-1) then
          if (iamp0) then
            write (*,*)
            write (*,*) '### ERROR in CHECK_PROC_TOPOLOGY:'
            write (*,*) '### Processor topology specification error.'
            write (*,*) 'Invalid number of processors requested.'
            write (*,*) 'Dimension = ',i
            write (*,*) 'Number of processors requested = ',
     &                  nprocs(i)
          end if
          call endrun (.true.)
        end if
      enddo
c
c ****** Set the optimal values of the topology if requested.
c
      call set_proc_topology
c
c ****** Check that the number of processors available
c ****** matches the number requested.
c
      nreq=nprocs(1)*nprocs(2)*nprocs(3)
c
      if (nreq.ne.nproc) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in CHECK_PROC_TOPOLOGY:'
          write (*,*) '### Processor topology specification error.'
          write (*,*) 'The number of processors requested does not'//
     &                ' equal the number available.'
          write (*,*) 'Number of processors requested = ',nreq
          write (*,*) 'Number of processors available = ',nproc
        end if
        call endrun (.true.)
      end if
c
      return
      end
c#######################################################################
      subroutine decompose_domain
c
c-----------------------------------------------------------------------
c
c ****** Decompose the domain into a Cartesian MPI topology.
c
c-----------------------------------------------------------------------
c
      use mpidefs
      use decomposition_params
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: ierr
c
c-----------------------------------------------------------------------
c
      integer, parameter :: ndim=3
      integer, dimension(ndim) :: coords
      logical, dimension(ndim) :: periodic
      logical :: reorder
      logical, dimension(ndim) :: keep_dim
c
c-----------------------------------------------------------------------
c
c ****** Create a communicator over all processors, COMM_ALL,
c ****** that has a Cartesian topology.
c
c ****** Specify the periodicity of the coordinate system.
c
      periodic(1)=.false.
      periodic(2)=.false.
      periodic(3)=.true.
c
c ****** Allow re-ordering in the Cartesian topology.
c
      reorder=.true.
c
      call MPI_Cart_create (MPI_COMM_WORLD,ndim,nprocs,
     &                      periodic,reorder,comm_all,ierr)
c
c ****** Get the index (rank) of the local processor in
c ****** communicator COMM_ALL in variable IPROC.
c
c ****** IMPORTANT NOTE:
c ****** If re-odering was allowed in the Cartesian topology
c ****** creation (above), then the rank of the local processor
c ****** in communicator COMM_ALL may be different from its rank
c ****** in communicator MPI_COMM_WORLD.
c
      call MPI_Comm_rank (comm_all,iproc,ierr)
c
c ****** Set the processor rank IPROC0 in communicator COMM_ALL
c ****** for the processor that has rank 0 in MPI_COMM_WORLD.
c ****** This value is broadcast to all the processors.
c
      if (iamp0) then
        iproc0=iproc
      end if
      call MPI_Bcast (iproc0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c
c ****** Get the coordinate indices of this processor in the
c ****** Cartesian MPI topology.
c
      call MPI_Cart_coords (comm_all,iproc,ndim,coords,ierr)
c
      iproc_r=coords(1)
      iproc_t=coords(2)
      iproc_p=coords(3)
c
      nproc_r=nprocs(1)
      nproc_t=nprocs(2)
      nproc_p=nprocs(3)
c
c ****** Get the rank of the neighboring processors in the
c ****** Cartesian MPI topology.
c
      call MPI_Cart_shift (comm_all,0,1,iproc_rm,iproc_rp,ierr)
      call MPI_Cart_shift (comm_all,1,1,iproc_tm,iproc_tp,ierr)
      call MPI_Cart_shift (comm_all,2,1,iproc_pm,iproc_pp,ierr)
c
c ****** Create communicators for operations involving all
c ****** processors in the phi dimension.  These communicators
c ****** are stored in COMM_PHI (and generally represent different
c ****** communicators on different processors).
c
      keep_dim(1)=.false.
      keep_dim(2)=.false.
      keep_dim(3)=.true.
c
      call MPI_Cart_sub (comm_all,keep_dim,comm_phi,ierr)
c
c ****** Create communicators for operations involving
c ****** all processors in the r dimension.
c ****** These communicators are stored in COMM_R
c ****** (and generally represent different communicators on
c ****** different processors).
c
      keep_dim(1)=.true.
      keep_dim(2)=.false.
      keep_dim(3)=.false.
c
      call MPI_Cart_sub (comm_all,keep_dim,comm_r,ierr)
c
      return
      end
c#######################################################################
      subroutine decompose_mesh_r
c
c-----------------------------------------------------------------------
c
c ****** Decompose the r mesh between processors.
c
c-----------------------------------------------------------------------
c
      use global_dims
      use local_dims_r
      use decomposition
      use solve_params
      use mpidefs
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: ierr,i,npts
      integer :: i0_h,i1_h,i0_m,i1_m
      integer, dimension(nproc_r) :: mp_r
c
c-----------------------------------------------------------------------
c
c ****** Decompose the r dimension.
c
      npts=nr_g
c
      call decompose_dimension (npts,nproc_r,mp_r,ierr)
      if (ierr.ne.0) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in DECOMPOSE_MESH:'
          write (*,*) '### Anomaly in decomposing the mesh'//
     &                ' between processors.'
          write (*,*) '### Could not decompose the r mesh.'
          write (*,*) 'Number of mesh points in r = ',npts
          write (*,*) 'Number of processors along r = ',nproc_r
        end if
        call endrun (.true.)
      end if
c
c ****** Check that the resulting mesh topology is valid.
c
      call check_mesh_topology (nproc_r,mp_r,1,'r')
c
c ****** Compute the mapping between the processor decomposition
c ****** and the global mesh.
c
c ****** Note that there is a two-point overlap in the mesh
c ****** between adjacent processors in r.
c
      i0_g=1
      do i=1,iproc_r
        i0_g=i0_g+mp_r(i)
      enddo
      nr=mp_r(iproc_r+1)+2
      i1_g=i0_g+nr-1
c
      nrm1=nr-1
c
c ****** Set the flags to indicate whether this processor has
c ****** points on the physical boundaries.
c
      if (iproc_r.eq.0) then
        rb0=.true.
      else
        rb0=.false.
      end if
c
      if (iproc_r.eq.nproc_r-1) then
        rb1=.true.
      else
        rb1=.false.
      end if
c
c ****** Set the dimensions of arrays for the "main" meshes
c ****** (i.e., the "m" mesh) for which normal derivatives are
c ****** needed (e.g., v).  These vary on different processors,
c ****** depending if they are left-boundary, internal, or
c ****** right-boundary processors.
c
      if (rb1) then
        nrm=nrm1
      else
        nrm=nr
      end if
c
c ****** Store the mapping structure (for this processor).
c
      allocate (map_rh(0:nproc-1))
      allocate (map_rm(0:nproc-1))
c
      if (rb0) then
        i0_h=1
      else
        i0_h=2
      end if
      if (rb1) then
        i1_h=nr
      else
        i1_h=nrm1
      end if
c
      if (rb0) then
        i0_m=1
      else
        i0_m=2
      end if
      i1_m=nrm1
c
      map_rh(iproc)%i0=i0_h
      map_rh(iproc)%i1=i1_h
c
      map_rm(iproc)%i0=i0_m
      map_rm(iproc)%i1=i1_m
c
      map_rh(iproc)%offset=i0_g+map_rh(iproc)%i0-1
      map_rh(iproc)%n=map_rh(iproc)%i1-map_rh(iproc)%i0+1
c
      map_rm(iproc)%offset=i0_g+map_rm(iproc)%i0-1
      map_rm(iproc)%n=map_rm(iproc)%i1-map_rm(iproc)%i0+1
c
c ****** Gather the mapping information by communicating among
c ****** all processors.
c
      call gather_mapping_info (map_rh)
      call gather_mapping_info (map_rm)
c
      return
      end
c#######################################################################
      subroutine decompose_mesh_tp
c
c-----------------------------------------------------------------------
c
c ****** Decompose the theta and phi mesh between processors.
c
c-----------------------------------------------------------------------
c
      use global_dims
      use local_dims_tp
      use decomposition
      use mpidefs
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: ierr,j,k
      integer :: j0_h,j1_h,j0_m,j1_m
      integer :: k0_h,k1_h,k0_m,k1_m
      integer, dimension(nproc_t) :: mp_t
      integer, dimension(nproc_p) :: mp_p
c
c-----------------------------------------------------------------------
c
c ****** Decompose the t dimension.
c
      call decompose_dimension (nt_g,nproc_t,mp_t,ierr)
      if (ierr.ne.0) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in DECOMPOSE_MESH_TP:'
          write (*,*) '### Anomaly in decomposing the mesh'//
     &                ' between processors.'
          write (*,*) '### Could not decompose the theta mesh.'
          write (*,*) 'Number of mesh points in theta = ',nt_g
          write (*,*) 'Number of processors along theta = ',nproc_t
        end if
        call endrun (.true.)
      end if
c
c ****** Decompose the p dimension.
c
      call decompose_dimension (np_g,nproc_p,mp_p,ierr)
      if (ierr.ne.0) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in DECOMPOSE_MESH_TP:'
          write (*,*) '### Anomaly in decomposing the mesh'//
     &                ' between processors.'
          write (*,*) '### Could not decompose the phi mesh.'
          write (*,*) 'Number of mesh points in phi = ',np_g
          write (*,*) 'Number of processors along phi = ',nproc_p
        end if
        call endrun (.true.)
      end if
c
c ****** Check that the resulting mesh topology is valid.
c
      call check_mesh_topology (nproc_t,mp_t,1,'theta')
      call check_mesh_topology (nproc_p,mp_p,1,'phi')
c
c ****** Set the flag for an axisymmetric run (requested by
c ****** setting NP_G = 3).
c
      if (np_g.eq.3) then
        axisymmetric=.true.
      else
        axisymmetric=.false.
      end if
c
c ****** Compute the mapping between the processor decomposition
c ****** and the global mesh.
c
c ****** Note that there is a two-point overlap in the mesh
c ****** between adjacent processors in theta and phi.
c
      j0_g=1
      do j=1,iproc_t
        j0_g=j0_g+mp_t(j)
      enddo
      nt=mp_t(iproc_t+1)+2
      j1_g=j0_g+nt-1
c
      k0_g=1
      do k=1,iproc_p
        k0_g=k0_g+mp_p(k)
      enddo
      np=mp_p(iproc_p+1)+2
      k1_g=k0_g+np-1
c
      ntm1=nt-1
      npm1=np-1
c
c ****** Set the flags to indicate whether this processor has
c ****** points on the physical boundaries.
c
      if (iproc_t.eq.0) then
        tb0=.true.
      else
        tb0=.false.
      end if
c
      if (iproc_t.eq.nproc_t-1) then
        tb1=.true.
      else
        tb1=.false.
      end if
c
c ****** Set the dimensions of arrays for the "main" meshes
c ****** (i.e., the "m" mesh) for which normal derivatives are
c ****** needed (e.g., v).  These vary on different processors,
c ****** depending if they are left-boundary, internal, or
c ****** right-boundary processors.
c
      if (tb1) then
        ntm=ntm1
      else
        ntm=nt
      end if
c
c ****** Since the phi dimension is periodic, all processors
c ****** have the same mesh limits.
c
      npm=np
c
c ****** Store the mapping structure (for this processor).
c
      allocate (map_th(0:nproc-1))
      allocate (map_tm(0:nproc-1))
      allocate (map_ph(0:nproc-1))
      allocate (map_pm(0:nproc-1))
c
      if (tb0) then
        j0_h=1
      else
        j0_h=2
      end if
      if (tb1) then
        j1_h=nt
      else
        j1_h=ntm1
      end if
c
      if (tb0) then
        j0_m=1
      else
        j0_m=2
      end if
      j1_m=ntm1
c
      if (iproc_p.eq.0) then
        k0_m=1
      else
        k0_m=2
      end if
      k1_m=npm1
c
      if (iproc_p.eq.0) then
        k0_h=1
      else
        k0_h=2
      end if
      if (iproc_p.eq.nproc_p-1) then
        k1_h=np
      else
        k1_h=npm1
      end if
c
      map_th(iproc)%i0=j0_h
      map_th(iproc)%i1=j1_h
c
      map_tm(iproc)%i0=j0_m
      map_tm(iproc)%i1=j1_m
c
      map_ph(iproc)%i0=k0_h
      map_ph(iproc)%i1=k1_h
c
      map_pm(iproc)%i0=k0_m
      map_pm(iproc)%i1=k1_m
c
      map_th(iproc)%offset=j0_g+map_th(iproc)%i0-1
      map_th(iproc)%n=map_th(iproc)%i1-map_th(iproc)%i0+1
c
      map_tm(iproc)%offset=j0_g+map_tm(iproc)%i0-1
      map_tm(iproc)%n=map_tm(iproc)%i1-map_tm(iproc)%i0+1
c
      map_ph(iproc)%offset=k0_g+map_ph(iproc)%i0-1
      map_ph(iproc)%n=map_ph(iproc)%i1-map_ph(iproc)%i0+1
c
      map_pm(iproc)%offset=k0_g+map_pm(iproc)%i0-1
      map_pm(iproc)%n=map_pm(iproc)%i1-map_pm(iproc)%i0+1
c
c ****** Gather the mapping information by communicating among
c ****** all processors.
c
      call gather_mapping_info (map_th)
      call gather_mapping_info (map_tm)
      call gather_mapping_info (map_ph)
      call gather_mapping_info (map_pm)
c
      return
      end
c#######################################################################
      subroutine check_mesh_dimensions (nr_g,nt_g,np_g)
c
c-----------------------------------------------------------------------
c
c ****** Check that the requested (global) mesh dimensions are
c ****** valid.
c
c-----------------------------------------------------------------------
c
      use mpidefs
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: nr_g,nt_g,np_g
c
c-----------------------------------------------------------------------
c
      if (nr_g.lt.4) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in CHECK_MESH_DIMENSIONS:'
          write (*,*) '### Invalid number of r mesh points'//
     &              ' requested.'
          write (*,*) '### The minimum number of mesh points is 4.'
          write (*,*)
          write (*,*) 'Number of mesh points requested = ',nr_g
        end if
        call endrun (.true.)
      end if
c
      if (nt_g.lt.4) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in CHECK_MESH_DIMENSIONS:'
          write (*,*) '### Invalid number of theta mesh points'//
     &              ' requested.'
          write (*,*) '### The minimum number of mesh points is 4.'
          write (*,*)
          write (*,*) 'Number of mesh points requested = ',nt_g
        end if
        call endrun (.true.)
      end if
c
      if (np_g.lt.3) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in CHECK_MESH_DIMENSIONS:'
          write (*,*) '### Invalid number of phi mesh points'//
     &              ' requested.'
          write (*,*) '### The minimum number of mesh points is 3.'
          write (*,*)
          write (*,*) 'Number of mesh points requested = ',np_g
        end if
        call endrun (.true.)
      end if
c
      return
      end
c#######################################################################
      subroutine decompose_dimension (nx,np,mp,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Decompose the mesh points NX along NP processors.
c
c ****** The decomposed mesh points are returned in array MP.
c
c-----------------------------------------------------------------------
c
c ****** This routine attempts to assign the mesh points as equally
c ****** as possible between the processors.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: nx
      integer :: np
      integer, dimension(np) :: mp
      integer :: ierr
c
c-----------------------------------------------------------------------
c
      integer :: nxm2,mpav,nrem
c
c-----------------------------------------------------------------------
c
      ierr=0
c
      nxm2=nx-2
c
      if (nxm2.le.0) then
        ierr=1
        return
      end if
c
      if (np.le.0) then
        ierr=2
        return
      end if
c
      mpav=nxm2/np
c
      mp(:)=mpav
c
      nrem=nxm2-mpav*np
c
      mp(1:nrem)=mp(1:nrem)+1
c
      return
      end
c#######################################################################
      subroutine check_mesh_topology (np,mp,min_pts,coord)
c
c-----------------------------------------------------------------------
c
c ****** Check the validity of the requested mesh topology.
c
c-----------------------------------------------------------------------
c
      use mpidefs
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: np
      integer, dimension(np) :: mp
      integer :: min_pts
      character(*) :: coord
c
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c
c ****** Check that the number of mesh points on each processor
c ****** is valid.
c
      do i=1,np
        if (mp(i).lt.min_pts) then
          if (iamp0) then
            write (*,*)
            write (*,*) '### ERROR in CHECK_MESH_TOPOLOGY:'
            write (*,*) '### Mesh topology specification error.'
            write (*,*) 'Invalid number of ',coord,
     &                  ' mesh points requested.'
            write (*,*) 'Processor index = ',i
            write (*,*) 'Number of mesh points requested = ',mp(i)
            write (*,*) 'Minimum number of mesh points allowed = ',
     &                  min_pts
          end if
          call endrun (.true.)
        end if
      enddo
c
      return
      end
c#######################################################################
      subroutine gather_mapping_info (map)
c
c-----------------------------------------------------------------------
c
c ****** Gather a mapping information array by communicating
c ****** among all processors.
c
c-----------------------------------------------------------------------
c
      use mpidefs
      use decomposition
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(map_struct), dimension(0:nproc-1) :: map
c
c-----------------------------------------------------------------------
c
c ****** Buffers for packing the data.
c
      integer, parameter :: lbuf=4
      integer, dimension(lbuf) :: sbuf
      integer, dimension(lbuf,0:nproc-1) :: rbuf
c
c-----------------------------------------------------------------------
c
      integer :: ierr,irank
c
c-----------------------------------------------------------------------
c
c ****** Put the local section of the mapping information
c ****** array into a buffer.
c
      sbuf(1)=map(iproc)%n
      sbuf(2)=map(iproc)%i0
      sbuf(3)=map(iproc)%i1
      sbuf(4)=map(iproc)%offset
c
c ****** Communicate among all processors.  After this call, all
c ****** processors have the complete mapping information.
c
      call MPI_Allgather (sbuf,lbuf,MPI_INTEGER,
     &                    rbuf,lbuf,MPI_INTEGER,comm_all,ierr)
c
c ****** Extract the mapping information from the buffer.
c
      do irank=0,nproc-1
        map(irank)%n     =rbuf(1,irank)
        map(irank)%i0    =rbuf(2,irank)
        map(irank)%i1    =rbuf(3,irank)
        map(irank)%offset=rbuf(4,irank)
      enddo
c
      return
      end
c#######################################################################
      subroutine decomp_diags
c
c-----------------------------------------------------------------------
c
c ****** Print diagnostics about the mesh decomposition.
c
c-----------------------------------------------------------------------
c
      use global_dims
      use global_mesh
      use local_dims_r
      use local_mesh_r
      use local_dims_tp
      use local_mesh_tp
      use mpidefs
      use solve_params
      use debug
      use decomposition
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: ierr
      integer :: irank
      real(r_typ) :: n_per_grid_min,n_per_grid_max
c
c-----------------------------------------------------------------------
c
      if (iamp0) then
c
        n_per_grid_min=floor(real(nr_g)/nproc_r)
     &                *floor(real(nt_g)/nproc_t)
     &                *floor(real(np_g)/nproc_p)
c
        n_per_grid_max=ceiling(real(nr_g)/nproc_r)
     &                *ceiling(real(nt_g)/nproc_t)
     &                *ceiling(real(np_g)/nproc_p)
c
        write (*,*)
        write (*,*) 'Total number of MPI ranks = ',nproc
        write (*,*)
        write (*,*) 'Number of MPI ranks in r = ',nproc_r
        write (*,*) 'Number of MPI ranks in t = ',nproc_t
        write (*,*) 'Number of MPI ranks in p = ',nproc_p
        write (*,*)
        write (*,*) 'Global mesh dimension in r = ',nr_g
        write (*,*) 'Global mesh dimension in t = ',nt_g
        write (*,*) 'Global mesh dimension in p = ',np_g
        write (*,*)
        write (*,'(A,F6.1)') ' Average # of r mesh pts per rank = ',
     &               real(nr_g)/nproc_r
        write (*,'(A,F6.1)') ' Average # of t mesh pts per rank = ',
     &               real(nt_g)/nproc_t
        write (*,'(A,F6.1)') ' Average # of p mesh pts per rank = ',
     &               real(np_g)/nproc_p
        write (*,*)
        write (*,'(A,F6.2,A)') ' Estimated load imbalance = ',
     &      100.0*(1.0-real(n_per_grid_min)/real(n_per_grid_max)),' %'
c
        write (9,*)
        write (9,*) 'Total number of MPI ranks = ',nproc
        write (9,*)
        write (9,*) 'Number of MPI ranks in r = ',nproc_r
        write (9,*) 'Number of MPI ranks in t = ',nproc_t
        write (9,*) 'Number of MPI ranks in p = ',nproc_p
        write (9,*)
        write (9,*) 'Global mesh dimension in r = ',nr_g
        write (9,*) 'Global mesh dimension in t = ',nt_g
        write (9,*) 'Global mesh dimension in p = ',np_g
        write (9,*)
        write (9,'(A,F6.1)') ' Average # of r mesh pts per rank = ',
     &               real(nr_g)/nproc_r
        write (9,'(A,F6.1)') ' Average # of t mesh pts per rank = ',
     &               real(nt_g)/nproc_t
        write (9,'(A,F6.1)') ' Average # of p mesh pts per rank = ',
     &               real(np_g)/nproc_p
        write (9,*)
        write (9,'(A,F6.2,A)') ' Estimated load imbalance = ',
     &      100.0*(1.0-real(n_per_grid_min)/real(n_per_grid_max)),' %'
c
      end if
c
      if (idebug.le.1) return
c
      do irank=0,nproc-1
        call MPI_Barrier (comm_all,ierr)
        if (irank.eq.iproc) then
          write (*,*)
          write (*,*)
          write (*,*) 'Rank id = ',iproc
          write (*,*) 'nr = ',nr
          write (*,*) 'nt = ',nt
          write (*,*) 'np = ',np
          write (*,*) 'i0_g = ',i0_g
          write (*,*) 'i1_g = ',i1_g
          write (*,*) 'j0_g = ',j0_g
          write (*,*) 'j1_g = ',j1_g
          write (*,*) 'k0_g = ',k0_g
          write (*,*) 'k1_g = ',k1_g
          write (*,*) 'Rank index in r    = ',iproc_r
          write (*,*) 'Rank index in t    = ',iproc_t
          write (*,*) 'Rank index in p    = ',iproc_p
          write (*,*) 'Rank to left  in r = ',iproc_rm
          write (*,*) 'Rank to right in r = ',iproc_rp
          write (*,*) 'Rank to left  in t = ',iproc_tm
          write (*,*) 'Rank to right in t = ',iproc_tp
          write (*,*) 'Rank to left  in p = ',iproc_pm
          write (*,*) 'Rank to right in p = ',iproc_pp
          write (*,*)
          write (*,*) 'Rank in MPI_COMM_WORLD = ',iprocw
          write (*,*) 'Rank in COMM_ALL       = ',iproc
          if (idebug.gt.2) then
            write (*,*)
            write (*,*) 'r mesh:'
            write (*,*) r
            write (*,*)
            write (*,*) 'theta mesh:'
            write (*,*) t
            write (*,*)
            write (*,*) 'phi mesh:'
            write (*,*) p
          end if
        end if
      enddo
c
      return
      end
c#######################################################################
      subroutine allocate_global_arrays
c
c-----------------------------------------------------------------------
c
c ****** Allocate global arrays.
c
c-----------------------------------------------------------------------
c
      use global_dims
      use global_mesh
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c ****** Allocate global mesh arrays.
c
      allocate (r_g (nrm1_g))
      allocate (dr_g(nrm1_g))
c
      allocate (rh_g (nr_g))
      allocate (drh_g(nr_g))
c
      allocate (t_g (ntm1_g))
      allocate (dt_g(ntm1_g))
c
      allocate (th_g (nt_g))
      allocate (dth_g(nt_g))
c
      allocate (p_g  (np_g))
      allocate (dp_g (np_g))
      allocate (ph_g (np_g))
      allocate (dph_g(np_g))
c
      allocate (st_g(ntm1_g))
      allocate (ct_g(ntm1_g))
c
      allocate (sth_g(nt_g))
      allocate (cth_g(nt_g))
c
      allocate (sp_g (np_g))
      allocate (cp_g (np_g))
      allocate (sph_g(np_g))
      allocate (cph_g(np_g))
c
      return
      end
c#######################################################################
      subroutine allocate_local_arrays_r
c
c-----------------------------------------------------------------------
c
c ****** Allocate local arrays for the r dimension and 3D arrays.
c
c-----------------------------------------------------------------------
c
      use local_dims_r
      use local_mesh_r
      use local_dims_tp
      use fields
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      allocate (r (nrm))
      allocate (r2 (nrm))
      allocate (dr(nrm))
c
      allocate (rh (nr))
      allocate (drh(nr))
c
c ****** Allocate inverse quantities.
c
      allocate (r_i (nrm))
      allocate (dr_i(nrm))
c
      allocate (rh_i (nr))
      allocate (drh_i(nr))
c
c ****** Allocate the potential array and cg scratch array.
c
      allocate (phi(nr,nt,np))
      allocate (x_ax(nr,nt,np))
      phi(:,:,:)=0.
      x_ax(:,:,:)=0.
!$acc enter data copyin(phi,x_ax)
c
c ****** Allocate polar boundary arrays.
c
      allocate (sum0(nr))
      allocate (sum1(nr))
      sum0(:)=0.
      sum1(:)=0.
!$acc enter data copyin(sum0,sum1)
c
c ****** Allocate the local magnetic field arrays.
c
      allocate (br(nrm,nt,np))
      allocate (bt(nr,ntm,np))
      allocate (bp(nr,nt,npm))
      br(:,:,:)=0.
      bt(:,:,:)=0.
      bp(:,:,:)=0.
c
      end subroutine
c#######################################################################
      subroutine allocate_local_arrays_tp
c
c-----------------------------------------------------------------------
c
c ****** Allocate local arrays for the theta and phi dimensions.
c
c-----------------------------------------------------------------------
c
      use local_dims_tp
      use local_mesh_tp
      use fields
      use vars
      use solve_params
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      allocate (t (ntm))
      allocate (dt(ntm))
c
      allocate (th (nt))
      allocate (dth(nt))
c
      allocate (p (np))
      allocate (dp(np))
c
      allocate (ph (np))
      allocate (dph(np))
c
      allocate (st(ntm))
      allocate (ct(ntm))
c
      allocate (sth(nt))
      allocate (cth(nt))
c
      allocate (sp (np))
      allocate (cp (np))
      allocate (sph(np))
      allocate (cph(np))
c
c ****** Allocate inverse quantities.
c
      allocate (dt_i(ntm))
      allocate (st_i(ntm))
c
      allocate (dth_i(nt))
      allocate (sth_i(nt))
c
      allocate (dp_i (np))
      allocate (dph_i(np))
c
c ****** Allocate the boundary radial magnetic field array(s).
c
      allocate (br0(nt,np))
      br0(:,:)=0.
!$acc enter data copyin(br0)
c
      if (validation_run) then
        allocate (br1(nt,np))
        br1(:,:)=0.
!$acc enter data copyin(br1)
      end if
c
      return
      end
c#######################################################################
      subroutine set_global_mesh
c
c-----------------------------------------------------------------------
c
c ****** Define the global mesh arrays.
c
c-----------------------------------------------------------------------
c
      use number_types
      use global_dims
      use global_mesh
      use meshdef
      use mpidefs
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: zero=0._r_typ
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: half=.5_r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,j,k
c
c-----------------------------------------------------------------------
c
c ****** Define the radial mesh.
c
      call genmesh (9,'r',nrm1_g,r0,r1,nmseg,rfrac,drratio,
     &              nfrmesh,.false.,zero,r_g)
c
      do i=2,nrm1_g
        rh_g(i)=.5*(r_g(i)+r_g(i-1))
        drh_g(i)=r_g(i)-r_g(i-1)
      enddo
      rh_g(1)=rh_g(2)-drh_g(2)
      rh_g(nr_g)=rh_g(nrm1_g)+drh_g(nrm1_g)
      drh_g(1)=drh_g(2)
      drh_g(nr_g)=drh_g(nrm1_g)
c
      do i=1,nrm1_g
        dr_g(i)=rh_g(i+1)-rh_g(i)
      enddo
c
c ****** Define the theta mesh.
c
      call genmesh (9,'t',ntm1_g,t0,t1,nmseg,tfrac,dtratio,
     &              nftmesh,.false.,zero,t_g)
c
      do j=2,ntm1_g
        th_g(j)=.5*(t_g(j)+t_g(j-1))
        dth_g(j)=t_g(j)-t_g(j-1)
      enddo
      th_g(1)=th_g(2)-dth_g(2)
      th_g(nt_g)=th_g(ntm1_g)+dth_g(ntm1_g)
      dth_g(1)=dth_g(2)
      dth_g(nt_g)=dth_g(ntm1_g)
c
      do j=1,ntm1_g
        dt_g(j)=th_g(j+1)-th_g(j)
      enddo
c
c ****** Define the periodic phi mesh.
c
      call genmesh (9,'p',npm1_g,p0,p1,nmseg,pfrac,dpratio,
     &              nfpmesh,.true.,phishift,p_g)
      p_g(np_g)=p_g(2)+pl
c
      do k=2,np_g
        ph_g(k)=half*(p_g(k)+p_g(k-1))
        dph_g(k)=p_g(k)-p_g(k-1)
      enddo
      ph_g(1)=ph_g(npm1_g)-pl
      dph_g(1)=dph_g(npm1_g)
c
      do k=1,npm1_g
        dp_g(k)=ph_g(k+1)-ph_g(k)
      enddo
      dp_g(np_g)=dp_g(2)
c
c ****** Enforce exact periodicity to protect symmetry properties
c ****** from round-off errors (especially for axisymmetric cases).
c
      dph_g(np_g)=dph_g(2)
      dp_g(1)=dp_g(npm1_g)
c
c ****** Define global auxiliary mesh-related arrays.
c
      st_g(:)=sin(t_g(:))
      ct_g(:)=cos(t_g(:))
      sth_g(:)=sin(th_g(:))
      cth_g(:)=cos(th_g(:))
c
      sp_g(:)=sin(p_g(:))
      cp_g(:)=cos(p_g(:))
      sph_g(:)=sin(ph_g(:))
      cph_g(:)=cos(ph_g(:))
c
c ****** For an axisymmetric case, set the exact values of
c ****** sin(phi) and cos(phi) to preserve symmetry properties
c ****** in the presence of round-off errors.
c
      if (axisymmetric) then
        sp_g(2)=0.
        cp_g(2)=one
        sph_g(2)=0.
        cph_g(2)=-one
      end if
c
c ****** Enforce exact periodicity to protect symmetry properties
c ****** from round-off errors (especially for axisymmetric cases).
c
      sph_g(1)=sph_g(npm1_g)
      sph_g(np_g)=sph_g(2)
      cph_g(1)=cph_g(npm1_g)
      cph_g(np_g)=cph_g(2)
      sp_g(1)=sp_g(npm1_g)
      sp_g(np_g)=sp_g(2)
      cp_g(1)=cp_g(npm1_g)
      cp_g(np_g)=cp_g(2)
c
      return
      end
c#######################################################################
      subroutine set_local_mesh_r
c
c-----------------------------------------------------------------------
c
c ****** Define the local r mesh arrays.
c
c-----------------------------------------------------------------------
c
      use number_types
      use global_dims
      use global_mesh
      use local_dims_r
      use local_mesh_r
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c
c ****** Define the local meshes.
c
      do i=1,nrm
        r(i)=r_g(i0_g+i-1)
        dr(i)=dr_g(i0_g+i-1)
      enddo
c
      dr1=dr(1)
      drn=dr(nrm)
c
      do i=1,nr
        rh(i)=rh_g(i0_g+i-1)
        drh(i)=drh_g(i0_g+i-1)
      enddo
c
c ****** Define local auxiliary mesh-related arrays.
c
      r2(:)=r(:)**2
      r_i(:)=one/r(:)
      dr_i(:)=one/dr(:)
      rh_i(:)=one/rh(:)
      drh_i(:)=one/drh(:)
c
!$acc enter data copyin(r,r2,r_i,dr,dr_i,dr1,drn,rh,rh_i,drh,drh_i)
      return
      end
c#######################################################################
      subroutine set_local_mesh_tp
c
c-----------------------------------------------------------------------
c
c ****** Define the local theta and phi mesh arrays.
c
c-----------------------------------------------------------------------
c
      use number_types
      use global_dims
      use global_mesh
      use local_dims_tp
      use local_mesh_tp
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: j,k,j0,j1
c
c-----------------------------------------------------------------------
c
c ****** Define the local meshes.
c
      do j=1,ntm
        t(j)=t_g(j0_g+j-1)
        dt(j)=dt_g(j0_g+j-1)
      enddo
c
      do j=1,nt
        th(j)=th_g(j0_g+j-1)
        dth(j)=dth_g(j0_g+j-1)
      enddo
c
      do k=1,npm
        p(k)=p_g(k0_g+k-1)
        dp(k)=dp_g(k0_g+k-1)
      enddo
c
      do k=1,np
        ph(k)=ph_g(k0_g+k-1)
        dph(k)=dph_g(k0_g+k-1)
      enddo
c
c ****** Define local auxiliary mesh-related arrays.
c
      do j=1,ntm
        st(j)=st_g(j0_g+j-1)
        ct(j)=ct_g(j0_g+j-1)
      enddo
c
      do j=1,nt
        sth(j)=sth_g(j0_g+j-1)
        cth(j)=cth_g(j0_g+j-1)
      enddo
c
      do k=1,npm
        sp(k)=sp_g(k0_g+k-1)
        cp(k)=cp_g(k0_g+k-1)
      enddo
c
      do k=1,np
        sph(k)=sph_g(k0_g+k-1)
        cph(k)=cph_g(k0_g+k-1)
      enddo
c
      dt_i(:)=one/dt(:)
      dth_i(:)=one/dth(:)
      sth_i(:)=one/sth(:)
      dp_i(:)=one/dp(:)
      dph_i(:)=one/dph(:)
c
c ****** Prevent division by zero at the poles for sin(theta).
c
      if (tb0) then
        j0=2
      else
        j0=1
      end if
      if (tb1) then
        j1=ntm1-1
      else
        j1=ntm1
      end if
c
      st_i(:)=0.
      do j=j0,j1
        st_i(j)=one/st(j)
      enddo
c
!$acc enter data copyin(t,th,dt,dth,p,ph,dp,dph,st,ct,sth,cth,
!$acc&                  sp,cp,sph,cph,dt_i,dth_i,st_i,sth_i,dp_i,dph_i)
      return
      end
c#######################################################################
      subroutine genmesh (io,label,nc,c0,c1,nseg,frac,dratio,
     &                    nfilt,periodic,c_shift,c)
c
c-----------------------------------------------------------------------
c
c ****** Generate a one-dimensional mesh.
c
c-----------------------------------------------------------------------
c
c ****** Input arguments:
c
c          IO      : [integer]
c                    Fortran file unit number to which to write
c                    mesh diagnostics.  Set IO=0 if diagnostics
c                    are not of interest.  It is assumed that
c                    unit IO has been connected to a file prior
c                    to calling this routine.
c
c          LABEL   : [character(*)]
c                    Name for the mesh coordinate (example: 'x').
c
c          NC      : [integer]
c                    Number of mesh points to load.
c
c          C0      : [real]
c                    The starting location for the coordinate.
c
c          C1      : [real]
c                    The ending location for the coordinate.
c                    It is required that C1.gt.C0.
c
c          NSEG    : [integer]
c                    Maximum number of mesh segments.
c                    The mesh spacing in each segment varies
c                    exponentially with a uniform amplification
c                    factor.  The actual number of mesh segments
c                    used is NSEG or less.  It is obtained from the
c                    information in array FRAC.
c
c          FRAC    : [real array, dimension NSEG]
c                    The normalized positions of the mesh segment
c                    boundaries (as a fraction of the size of the
c                    domain).  For a non-periodic mesh, the first
c                    value of FRAC specified must equal 0. and the
c                    last value must equal 1.  For a periodic mesh,
c                    FRAC must not contain both 0. and 1., since
c                    these represent the same point.
c
c          DRATIO  : [real array, dimension NSEG]
c                    The ratio of the mesh spacing at the end of a
c                    segment to that at the beginning.
c
c          NFILT   : [integer]
c                    The number of times to filter the mesh-point
c                    distribution array.  Set NFILT=0 if filtering
c                    is not desired.  Filtering can reduce
c                    discontinuities in the derivative of the mesh
c                    spacing.
c
c          PERIODIC: [logical]
c                    A flag to indicate whether the mesh to be
c                    generated represents a periodic coordinate.
c                    If the coordinate is specified as periodic,
c                    the range [C0,C1] should be the whole periodic
c                    interval; the first mesh point is set at C0
c                    and the last mesh point, C(NC), is set at C1.
c
c          C_SHIFT : [real]
c                    Amount by which to shift the periodic coordinate.
c                    C_SHIFT is only used when PERIODIC=.true.,
c                    and is ignored otherwise.  A positive C_SHIFT
c                    moves the mesh points to the right.
c
c ****** Output arguments:
c
c          C       : [real array, dimension NC]
c                    The locations of the mesh points.
c
c-----------------------------------------------------------------------
c
c ****** The arrays DRATIO and FRAC define the mesh as follows.
c
c ****** For example, suppose that a (non-periodic) mesh with three
c ****** segments is desired.  Suppose the domain size is c=[0:2].
c ****** In the first segment (with c between 0 and .5) the mesh
c ****** spacing is decreasing with c, such that DC at c=.5 is half
c ****** DC at c=0.  From c=.5 to c=1, the mesh is uniform.  From c=1
c ****** to c=2, the mesh spacing is increasing with c such that DC at
c ****** c=2 is 10 times DC at c=1.  This mesh would be specified by:
c ******
c ******     FRAC=0.,.25,.5,1.
c ******     DRATIO=.5,1.,10.
c ******
c ****** The variable C_SHIFT can be used to shift the mesh point
c ****** distribution for a periodic coordinate.  For example,
c ****** suppose C represents mesh points in the interval [0,2*pi].
c ****** C_SHIFT=.5*pi would move the distribution of mesh points
c ****** so that the original mesh point with C(1)=0. would be
c ****** close to .5*pi in the new mesh.
c
c-----------------------------------------------------------------------
c
      use number_types
      use mpidefs
      use debug
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer, intent(in) :: io
      character(*), intent(in) :: label
      integer, intent(in) :: nc
      real(r_typ), intent(in) :: c0,c1
      integer, intent(in) :: nseg
      real(r_typ), dimension(nseg), intent(in) :: frac,dratio
      integer, intent(in) :: nfilt
      logical, intent(in) :: periodic
      real(r_typ), intent(in) :: c_shift
      real(r_typ), dimension(nc), intent(out) :: c
c
c-----------------------------------------------------------------------
c
c ****** Storage for the coordinate transformation.
c
      integer :: ns
      real(r_typ), dimension(:), allocatable :: xi,cs,a,r,f
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: zero=0._r_typ
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: half=.5_r_typ
      real(r_typ), parameter :: eps=1.e-5_r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,j,nf,nr,ll,j0
      real(r_typ) :: alpha,dr,fac,d,dxi,xiv,cshft,xi_shift
      real(r_typ), dimension(:), allocatable :: dc,rdc
c
c-----------------------------------------------------------------------
c
c ****** Check that the number of mesh points is valid.
c
      if (nc.lt.2) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in GENMESH:'
          write (*,*) '### Invalid number of mesh points specified.'
          write (*,*) '### There must be at least two mesh points.'
          write (*,*) 'Mesh coordinate: ',label
          write (*,*) 'Number of mesh points specified = ',nc
        end if
        call endrun (.true.)
      end if
c
c ****** Check that a positive mesh interval has been specified.
c
      if (c0.ge.c1) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in GENMESH:'
          write (*,*) '### Invalid mesh interval specified.'
          write (*,*) '### C1 must be greater than C0.'
          write (*,*) 'Mesh coordinate: ',label
          write (*,*) 'C0 = ',c0
          write (*,*) 'C1 = ',c1
        end if
        call endrun (.true.)
      end if
c
c ****** Find the number of values of FRAC specified.
c
      do nf=1,nseg-1
        if (frac(nf+1).eq.zero) exit
      enddo
c
c ****** When no values have been specified (NF=1, the default),
c ****** a uniform mesh is produced.
c
      if (nf.eq.1.and.frac(1).eq.zero) then
        ns=1
        allocate (cs(ns+1))
        allocate (r(ns))
        cs(1)=c0
        cs(2)=c1
        r(1)=one
        go to 100
      end if
c
c ****** Check that the specified values of FRAC are monotonically
c ****** increasing.
c
      do i=2,nf
        if (frac(i).lt.frac(i-1)) then
          if (iamp0) then
            write (*,*)
            write (*,*) '### ERROR in GENMESH:'
            write (*,*) '### Invalid mesh specification.'
            write (*,*) 'Mesh coordinate: ',label
            write (*,*) 'The values in FRAC must increase'//
     &                  ' monotonically.'
            write (*,*) 'FRAC = ',frac(1:nf)
          end if
          call endrun (.true.)
        end if
      enddo
c
c ****** Check the specified values of FRAC.
c
      if (periodic) then
c
c ****** A periodic mesh requires the specified values of FRAC
c ****** to be in the range 0. to 1.
c
        if (frac(1).lt.zero.or.frac(nf).gt.one) then
          if (iamp0) then
            write (*,*)
            write (*,*) '### ERROR in GENMESH:'
            write (*,*) '### Invalid mesh specification.'
            write (*,*) 'Mesh coordinate: ',label
            write (*,*) 'For a periodic coordinate, the values in'//
     &                  ' FRAC must be between 0. and 1.'
            write (*,*) 'FRAC = ',frac(1:nf)
          end if
          call endrun (.true.)
        end if
c
c ****** A periodic mesh cannot contain both 0. and 1. in FRAC,
c ****** since these represent the same point.
c
        if (frac(1).eq.zero.and.frac(nf).eq.one) then
          if (iamp0) then
            write (*,*)
            write (*,*) '### ERROR in GENMESH:'
            write (*,*) '### Invalid mesh specification.'
            write (*,*) 'Mesh coordinate: ',label
            write (*,*) 'For a periodic coordinate, FRAC must not'//
     &                  ' contain both 0. and 1.'
            write (*,*) 'FRAC = ',frac(1:nf)
          end if
          call endrun (.true.)
        end if
c
      else
c
c ****** A non-periodic mesh requires the first specified value
c ****** of FRAC to be 0., and the last value to equal 1.
c
        if (frac(1).ne.zero) then
          if (iamp0) then
            write (*,*)
            write (*,*) '### ERROR in GENMESH:'
            write (*,*) '### Invalid mesh specification.'
            write (*,*) 'Mesh coordinate: ',label
            write (*,*) 'For a non-periodic coordinate, the first'//
     &                ' value of FRAC must equal 0.'
            write (*,*) 'FRAC = ',frac(1:nf)
          end if
          call endrun (.true.)
        end if
c
        if (frac(nf).ne.one) then
          if (iamp0) then
            write (*,*)
            write (*,*) '### ERROR in GENMESH:'
            write (*,*) '### Invalid mesh specification.'
            write (*,*) 'Mesh coordinate: ',label
            write (*,*) 'For a non-periodic coordinate, the last'//
     &                  ' value of FRAC must equal 1.'
            write (*,*) 'FRAC = ',frac(1:nf)
          end if
          call endrun (.true.)
        end if
c
      end if
c
c ****** Check that the required values of DRATIO have been set,
c ****** and are positive.
c
      if (periodic) then
        nr=nf
      else
        nr=nf-1
      end if
c
      do i=1,nr
        if (dratio(i).le.zero) then
          if (iamp0) then
            write (*,*)
            write (*,*) '### ERROR in GENMESH:'
            write (*,*) '### Invalid mesh specification.'
            write (*,*) 'Mesh coordinate: ',label
            write (*,*) 'A required value in DRATIO has not been'//
     &                  ' set or is not positive.'
            write (*,*) 'DRATIO = ',dratio(1:nr)
          end if
          call endrun (.true.)
        end if
      enddo
c
c ****** Check that an inherently discontinuous mesh has not been
c ****** specified inadvertently.
c
      if (periodic.and.nr.eq.1.and.dratio(1).ne.one) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### WARNING from GENMESH:'
          write (*,*) '### Discontinuous mesh specification.'
          write (*,*) 'Mesh coordinate: ',label
          write (*,*) 'An inherently discontinuous mesh has been'//
     &                ' specified.'
          write (*,*) 'FRAC = ',frac(1:nf)
          write (*,*) 'DRATIO = ',dratio(1:nr)
        end if
      end if
c
c ****** Set the number of segments.
c
      ns=nf-1
c
c ****** For a periodic coordinate, add points at XI=0. and XI=1.
c ****** if they are not already present.
c
      if (periodic) then
        if (frac(1).ne.zero) ns=ns+1
        if (frac(nf).ne.one) ns=ns+1
      end if
c
      allocate (cs(ns+1))
      allocate (r(ns))
c
c ****** Set up the coordinate limits of the segments.
c
      if (periodic) then
        if (frac(1).ne.zero) then
          cs(1)=c0
          cs(2:nf+1)=c0+(c1-c0)*frac(1:nf)
          if (frac(nf).ne.one) then
            alpha=(one-frac(nf))/(frac(1)+one-frac(nf))
            r(1)=dratio(nr)/(one+alpha*(dratio(nr)-one))
            r(2:nr+1)=dratio(1:nr)
            cs(ns+1)=c1
            r(ns)=one+alpha*(dratio(nr)-one)
          else
            r(1)=dratio(nr)
            r(2:nr)=dratio(1:nr-1)
          end if
        else
          cs(1:nf)=c0+(c1-c0)*frac(1:nf)
          r(1:nr)=dratio(1:nr)
          cs(ns+1)=c1
        end if
      else
        cs(1:nf)=c0+(c1-c0)*frac(1:nf)
        r(1:nr)=dratio(1:nr)
      end if
c
  100 continue
c
      allocate (xi(ns+1))
      allocate (a(ns))
      allocate (f(ns))
c
c ****** Compute the XI values at the segment limits.
c
      do i=1,ns
        dr=r(i)-one
        if (abs(dr).lt.eps) then
          f(i)=(cs(i+1)-cs(i))*(one+half*dr)
        else
          f(i)=(cs(i+1)-cs(i))*log(r(i))/dr
        end if
      enddo
c
      fac=zero
      do i=ns,1,-1
        fac=fac/r(i)+f(i)
      enddo
c
      d=f(1)/fac
      xi(1)=zero
      do i=2,ns
        xi(i)=xi(i-1)+d
        if (i.lt.ns) d=d*f(i)/(f(i-1)*r(i-1))
      enddo
      xi(ns+1)=one
c
c ****** Set the amplification factor for each segment.
c
      do i=1,ns
        a(i)=log(r(i))/(xi(i+1)-xi(i))
      enddo
c
c ****** For a periodic coordinate, find the XI shift corresponding
c ****** to a shift of C_SHIFT in the coordinate.
c ****** Note that a positive value of C_SHIFT moves the mesh
c ****** points to the right.
c
      if (periodic) then
        cshft=-c_shift
        call map_c_to_xi (periodic,ns,xi,cs,a,r,cshft,xi_shift)
      else
        xi_shift=0.
      end if
c
c ****** Compute the location of the mesh points in array C
c ****** by mapping from the XI values.
c
      dxi=one/(nc-one)
c
      c(1)=c0
      do j=2,nc-1
        xiv=(j-1)*dxi
        call map_xi_to_c (periodic,ns,xi,cs,a,r,
     &                    cshft,xi_shift,xiv,c(j))
      enddo
      c(nc)=c1
c
c ****** Filter the mesh if requested.
c
      if (nfilt.gt.0) then
        do i=1,nfilt
          if (periodic) then
            call filter_coord_periodic (c1-c0,nc,c)
          else
            call filter_coord (nc,c)
          end if
        enddo
      end if
c
c ****** Write out the mesh information.
c
      if (io.gt.0.and.iamp0) then
c
        write (io,*)
        write (io,*) '### COMMENT from GENMESH:'
        write (io,*) '### Mesh information for coordinate ',label,':'
c
        if (idebug.gt.0) then
          write (io,*)
          write (io,*) 'Flag to indicate a periodic mesh: ',periodic
          write (io,*) 'Number of mesh points = ',nc
          write (io,*) 'Lower mesh limit = ',c0
          write (io,*) 'Upper mesh limit = ',c1
          write (io,*) 'Number of times to filter the mesh = ',nfilt
          if (periodic) then
            write (io,*) 'Amount to shift the mesh = ',c_shift
          end if
        end if
c
        write (io,*)
        write (io,*) 'Number of mesh segments = ',ns
c
        ll=len_trim(label)
c
        write (io,900) 'Segment      xi0       xi1'//
     &                 repeat (' ',16-ll)//label//'0'//
     &                 repeat (' ',16-ll)//label//'1'//
     &                 '            ratio'
        do i=1,ns
          write (io,910) i,xi(i),xi(i+1),cs(i),cs(i+1),r(i)
        enddo
c
        allocate (dc(nc))
        allocate (rdc(nc))
c
        dc=c-cshift(c,-1)
        if (periodic) dc(1)=dc(nc)
        rdc=dc/cshift(dc,-1)
        if (periodic) rdc(1)=rdc(nc)
c
        write (io,*)
        write (io,920) 'Mesh-point locations:'
        write (io,920) '     i'//
     &                 repeat (' ',18-ll)//label//
     &                 repeat (' ',17-ll)//'d'//label//
     &                 '             ratio'
c
        if (periodic) then
          j0=1
        else
          j0=3
          write (io,930) 1,c(1)
          write (io,930) 2,c(2),dc(2)
        end if
        do j=j0,nc
          write (io,930) j,c(j),dc(j),rdc(j)
        enddo
c
        deallocate (dc)
        deallocate (rdc)
c
      end if
c
  900 format (/,tr1,a)
  910 format (tr1,i4,2x,2f10.6,4f17.8)
  920 format (tr1,a)
  930 format (tr1,i6,3f18.8)
c
      deallocate (cs)
      deallocate (r)
      deallocate (xi)
      deallocate (a)
      deallocate (f)
c
      return
      end
c#######################################################################
      subroutine map_xi_to_c (periodic,ns,xi,cs,a,r,
     &                        cshft,xi_shift,xiv,cv)
c
c-----------------------------------------------------------------------
c
c ****** Get the mesh coordinate value CV for the specified
c ****** xi value XIV.
c
c ****** Set PERIODIC=.true. for a periodic coordinate.
c ****** NS is the number of segments in the mesh definition.
c ****** The arrays XI, CS, A, and R define the mesh mapping.
c
c ****** CSHFT represents the amount to shift a periodic coordinate.
c ****** XI_SHIFT represents the corresponding amount to shift xi.
c
c ****** This is a utility routine for GENMESH.
c
c-----------------------------------------------------------------------
c
      use number_types
      use mpidefs
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      logical, intent(in) :: periodic
      integer, intent(in) :: ns
      real(r_typ), dimension(ns+1), intent(in) :: xi,cs
      real(r_typ), dimension(ns), intent(in) :: a,r
      real(r_typ), intent(in) :: cshft,xi_shift
      real(r_typ), intent(in) :: xiv
      real(r_typ), intent(out) :: cv
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: zero=0._r_typ
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: half=.5_r_typ
      real(r_typ), parameter :: eps=1.e-5_r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i
      real(r_typ) :: xiv_p,d,d1,da,da1,fac
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: fold
c
c-----------------------------------------------------------------------
c
c ****** Find the index of the segment to which XIV belongs.
c
      if (periodic) then
c
c ****** Shift XIV by XI_SHIFT.
c
        xiv_p=xiv+xi_shift
c
c ****** Fold XIV_P into the main interval.
c
        xiv_p=fold(zero,one,xiv_p)
c
      else
c
        xiv_p=xiv
c
      end if
c
      do i=1,ns
        if (xiv_p.ge.xi(i).and.xiv_p.le.xi(i+1)) exit
      enddo
c
      if (i.gt.ns) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_XI_TO_C:'
          write (*,*) '### Error in finding the XI segment.'
          write (*,*) '### Could not find XIV in the XI table.'
          write (*,*) '[Utility routine for GENMESH.]'
          write (*,*) '[This is an internal error.]'
          write (*,*) 'XI = ',xi
          write (*,*) 'XIV = ',xiv
          write (*,*) 'XIV_P = ',xiv_p
        end if
        call endrun (.true.)
      end if
c
      d =xiv_p  -xi(i)
      d1=xi(i+1)-xi(i)
c
      da =a(i)*d
      da1=a(i)*d1
c
c ****** Interpolate the mapping function at XIV_P.
c
      if (abs(da1).lt.eps) then
        fac=(d*(one+half*da))/(d1*(one+half*da1))
      else
        fac=(exp(da)-one)/(r(i)-one)
      end if
c
      cv=cs(i)+(cs(i+1)-cs(i))*fac
c
      if (periodic) then
c
c ****** Shift CV by the amount CSHFT.
c
        cv=cv-cshft
c
c ****** Fold CV into the main interval.
c
        cv=fold(cs(1),cs(ns+1),cv)
c
      end if
c
      return
      end
c#######################################################################
      subroutine map_c_to_xi (periodic,ns,xi,cs,a,r,cv,xiv)
c
c-----------------------------------------------------------------------
c
c ****** Get the xi value XIV for the specified coordinate value CV.
c
c ****** Set PERIODIC=.true. for a periodic coordinate.
c ****** NS is the number of segments in the mesh definition.
c ****** The arrays XI, CS, A, and R define the mesh mapping.
c
c ****** This is a utility routine for GENMESH.
c
c-----------------------------------------------------------------------
c
      use number_types
      use mpidefs
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      logical, intent(in) :: periodic
      integer, intent(in) :: ns
      real(r_typ), dimension(ns+1), intent(in) :: xi,cs
      real(r_typ), dimension(ns), intent(in) :: a,r
      real(r_typ), intent(in) :: cv
      real(r_typ), intent(out) :: xiv
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: eps=1.e-5_r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i
      real(r_typ) :: cv_p,d,da,fac
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: fold
c
c-----------------------------------------------------------------------
c
c ****** Find the index of the segment to which CV belongs.
c
      if (periodic) then
c
c ****** Fold CV_P into the main interval.
c
        cv_p=fold(cs(1),cs(ns+1),cv)
c
      else
c
        cv_p=cv
c
      end if
c
      do i=1,ns
        if (cv_p.ge.cs(i).and.cv_p.le.cs(i+1)) exit
      enddo
c
      if (i.gt.ns) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in MAP_C_TO_XI:'
          write (*,*) '### Error in finding the CS segment.'
          write (*,*) '### Could not find CV in the CS table.'
          write (*,*) '[Utility routine for GENMESH.]'
          write (*,*) '[This is an internal error.]'
          write (*,*) 'CS = ',cs
          write (*,*) 'CV = ',cv
          write (*,*) 'CV_P = ',cv_p
        end if
        call endrun (.true.)
      end if
c
      d=(cv_p-cs(i))/(cs(i+1)-cs(i))
      da=(r(i)-one)*d
c
c ****** Interpolate the mapping function at XIV_P.
c
      if (abs(da).lt.eps) then
        fac=d*(xi(i+1)-xi(i))
      else
        fac=log(da+one)/a(i)
      end if
c
      xiv=xi(i)+fac
c
      return
      end
c#######################################################################
      subroutine filter_coord (n,f)
c
c-----------------------------------------------------------------------
c
c ****** Apply a "(1,2,1)/4" low-pass digital filter to a
c ****** 1D coordinate.
c
c ****** The end-points F(1) and F(N) are not changed.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: n
      real(r_typ), dimension(n) :: f
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: two=2._r_typ
      real(r_typ), parameter :: quarter=.25_r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(n) :: ff
c
c-----------------------------------------------------------------------
c
      integer :: i
c
c-----------------------------------------------------------------------
c
c ****** Make a copy of the function.
c
      ff=f
c
c ****** Apply the filter.
c
      do i=2,n-1
        f(i)=quarter*(ff(i-1)+two*ff(i)+ff(i+1))
      enddo
c
      return
      end
c#######################################################################
      subroutine filter_coord_periodic (xl,n,f)
c
c-----------------------------------------------------------------------
c
c ****** Apply a "(1,2,1)/4" low-pass digital filter to a
c ****** periodic 1D coordinate.
c
c-----------------------------------------------------------------------
c
c ****** XL is the periodic interval for the coordinate.
c
c ****** The filtered coordinate is translated so that F(1)
c ****** is preserved.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: xl
      integer :: n
      real(r_typ), dimension(n) :: f
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: two=2._r_typ
      real(r_typ), parameter :: quarter=.25_r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(0:n+1) :: ff
c
c-----------------------------------------------------------------------
c
      integer :: i
      real(r_typ) :: f1old,f1new
c
c-----------------------------------------------------------------------
c
c ****** Save the value of F(1).
c
      f1old=f(1)
c
c ****** Make a periodic copy of the function.
c
      ff(1:n)=f(:)
c
      ff(0)=f(n-1)-xl
      ff(n+1)=f(2)+xl
c
c ****** Apply the filter.
c
      do i=1,n
        f(i)=quarter*(ff(i-1)+two*ff(i)+ff(i+1))
      enddo
c
c ****** Translate F so that F(1) is preserved.
c
      f1new=f(1)
      do i=1,n
        f(i)=f(i)-f1new+f1old
      enddo
c
      return
      end
c#######################################################################
      function fold (x0,x1,x)
c
c-----------------------------------------------------------------------
c
c ****** "Fold" X into the periodic interval [X0,X1].
c
c ****** On return, X is such that X0.le.X.lt.X1.
c
c-----------------------------------------------------------------------
c
c ****** It is assumed that X0 does not equal X1, as is physically
c ****** necessary.  If X0 and X1 are equal, the routine just
c ****** returns with FOLD=X.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: fold
      real(r_typ) :: x0,x1,x
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: xl
c
c-----------------------------------------------------------------------
c
      fold=x
c
      if (x0.eq.x1) return
c
      xl=x1-x0
c
      fold=mod(x-x0,xl)+x0
c
      if (fold.lt.x0) fold=fold+xl
      if (fold.ge.x1) fold=fold-xl
c
      return
      end
c#######################################################################
      subroutine set_flux
c
c-----------------------------------------------------------------------
c
c ****** Set the radial magnetic field at the photosphere.
c
c-----------------------------------------------------------------------
c
      use number_types
      use global_dims
      use global_mesh
      use local_dims_r
      use local_mesh_r
      use local_dims_tp
      use local_mesh_tp
      use fields
      use vars
      use mpidefs
      use timing
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c ****** Global flux array.
c
      real(r_typ), dimension(:,:), allocatable :: br0_g
c
c-----------------------------------------------------------------------
c
      integer :: j,k,ierr
c
c-----------------------------------------------------------------------
c
      call timer_on
c
      allocate (br0_g(nt_g,np_g))
c
c ****** Define the global flux array.
c
c ****** Read the flux from file BR0FILE (only on processor IPROC0).
c
      if (iamp0) then
        call readbr (br0file,br0_g,ierr)
      end if
      call check_error_on_p0 (ierr)
c
c ****** Broadcast BR0_G to all the processors.
c
      call MPI_Bcast (br0_g,nt_g*np_g,ntype_real,0,comm_all,ierr)
c
c ****** For a fully open field, reverse negative Br
c ****** (i.e., use the monopole trick).
c
      if (option.eq.'open') then
c
c ****** Write the boundary flux (before the sign flip) to a file
c ****** if requested.
c
        if (iamp0) then
          if (br_photo_original_file.ne.'') then
            write (*,*)
            write (*,*) '### COMMENT from SET_FLUX:'
            write (*,*)
            write (*,*) 'Writing BR0 (before sign flip) to file: ',
     &                  trim(br_photo_original_file)
            write (9,*)
            write (9,*) '### COMMENT from SET_FLUX:'
            write (9,*)
            write (9,*) 'Writing BR0 (before sign flip) to file: ',
     &                  trim(br_photo_original_file)
            call wrhdf_2d (br_photo_original_file,
     &                     .true.,nt_g,np_g,
     &                     br0_g,th_g,ph_g,hdf32,ierr)
          end if
        end if
c
c ****** Reverse Br.
c
        br0_g(:,:)=abs(br0_g(:,:))
c
      end if
c
c ****** Write the boundary flux to a file if requested.
c
      if (iamp0) then
        if (br_photo_file.ne.' ') then
          write (*,*)
          write (*,*) '### COMMENT from SET_FLUX:'
          write (*,*)
          write (*,*) 'Writing BR0 to file: ',trim(br_photo_file)
          write (9,*)
          write (9,*) '### COMMENT from SET_FLUX:'
          write (9,*)
          write (9,*) 'Writing BR0 to file: ',trim(br_photo_file)
          call wrhdf_2d (br_photo_file,.true.,nt_g,np_g,
     &                   br0_g,th_g,ph_g,hdf32,ierr)
        end if
      end if
c
      do j=1,nt
        do k=1,np
          br0(j,k)=br0_g(j0_g+j-1,k0_g+k-1)
        enddo
      enddo
!$acc update device(br0)
c
      deallocate(br0_g)
c
      call timer_off (t_io)
c
      end subroutine
c#######################################################################
      subroutine potfld
c
c-----------------------------------------------------------------------
c
c ****** Find the potential field solution.
c
c-----------------------------------------------------------------------
c
      use number_types
      use local_dims_r
      use local_mesh_r
      use local_dims_tp
      use local_mesh_tp
      use fields
      use cgcom
      use solve_params
      use mpidefs
      use debug
      use timing
      use matrix_storage_pot3d_solve
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: ierr,nrm2,ntm2,npm2,i
c
      real(r_typ), dimension(:), allocatable :: rhs_cg,x_cg
c
c-----------------------------------------------------------------------
c
c ****** Load matrix and preconditioner.
c
      nrm2=nrm1-1
      ntm2=ntm1-1
      npm2=npm1-1
c
      a_offsets(1)=-nrm2*ntm2
      a_offsets(2)=-nrm2
      a_offsets(3)=-1
      a_offsets(4)= 0
      a_offsets(5)= 1
      a_offsets(6)= nrm2
      a_offsets(7)= nrm2*ntm2
c
c ****** Allocate cg 1D vectors.
c
      N=nrm2*ntm2*npm2
c
c ****** Prepare the guess, and rhs for the solve.
c
      allocate(rhs_cg(N))
      allocate(x_cg(N))
      rhs_cg(:)=0.
      x_cg(:)=0.
!$acc enter data copyin(rhs_cg,x_cg)
c
      call getM (N,a_offsets,M)
      call alloc_pot3d_matrix_coefs
      call load_matrix_pot3d_solve
!$acc enter data copyin(a)
      call load_preconditioner_pot3d_solve
!$acc enter data copyin(a_i)
c
c ****** Use a trick to accumulate the contribution of the
c ****** boundary conditions (i.e., the inhomogenous part).
c
      call set_boundary_points (x_ax,one)
      call seam (x_ax,nr,nt,np)
      call delsq (x_ax,rhs_cg)
c
c ****** Original rhs is zero so just use negative of boundary
c        trick contributions:
c
      do concurrent (i=1:N)
        rhs_cg(i)=-rhs_cg(i)
      enddo
c
c ****** Solve for the potential.
c
      if (idebug.gt.0.and.iamp0) then
        write (*,*)
        write (*,*) '### COMMENT from POTFLD:'
        write (*,*) '### Doing a solution:'
      end if
c
      call solve (x_cg,rhs_cg,N,ierr)
c
      if (ierr.ne.0) then
        call endrun (.true.)
      end if
c
      call unpack_scalar (phi,x_cg)
c
      call set_boundary_points (phi,one)
      call seam (phi,nr,nt,np)
c
!$acc exit data delete(rhs_cg,x_cg,a,a_i)
      call dealloc_pot3d_matrix_coefs
      deallocate(rhs_cg)
      deallocate(x_cg)
c
      return
      end
c#######################################################################
      subroutine set_validation_flux
c
c-----------------------------------------------------------------------
c
c ****** Set the radial magnetic field at the photosphere.
c ****** This uses the tilted dipole analytic solution for validation.
c
c-----------------------------------------------------------------------
c
      use number_types
      use global_dims
      use global_mesh
      use local_dims_r
      use local_mesh_r
      use local_dims_tp
      use local_mesh_tp
      use fields
      use vars
      use mpidefs
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: two=2.0_r_typ
c
c-----------------------------------------------------------------------
c
c ****** Global Br boundary arrays.
c
      real(r_typ), dimension(:,:), allocatable :: br0_g
      real(r_typ), dimension(:,:), allocatable :: br1_g
c
c-----------------------------------------------------------------------
c
      integer :: j,k
      real(r_typ) :: sum00,sum11
c
c-----------------------------------------------------------------------
c
c ****** Define the global boundary condition arrays.
c
      allocate (br0_g(nt_g,np_g))
      allocate (br1_g(nt_g,np_g))
c
c ****** Set the tilted dipole analytic Br.
c
      do j=1,nt_g
        do k=1,np_g
          br0_g(j,k)=(two/r0**3)*(cth_g(j)*cos(dipole_angle) +
     &                            sth_g(j)*cph_g(k)*sin(dipole_angle))
          br1_g(j,k)=(two/r1**3)*(cth_g(j)*cos(dipole_angle) +
     &                            sth_g(j)*cph_g(k)*sin(dipole_angle))
        enddo
      enddo
c
c ****** Set Br to be periodic (should not be nesessary?).
c
      br0_g(:,1)=br0_g(:,npm1_g)
      br0_g(:,np_g)=br0_g(:,2)
      br1_g(:,1)=br1_g(:,npm1_g)
      br1_g(:,np_g)=br1_g(:,2)
c
c ****** Set BCs at the poles.
c
      sum00=sum(br0_g(     2,2:npm1_g)*dph_g(2:npm1_g))*pl_i
      sum11=sum(br0_g(ntm1_g,2:npm1_g)*dph_g(2:npm1_g))*pl_i
c
      br0_g(1   ,:)=two*sum00-br0_g(     2,:)
      br0_g(nt_g,:)=two*sum11-br0_g(ntm1_g,:)
c
      sum00=sum(br1_g(     2,2:npm1_g)*dph_g(2:npm1_g))*pl_i
      sum11=sum(br1_g(ntm1_g,2:npm1_g)*dph_g(2:npm1_g))*pl_i
c
      br1_g(1   ,:)=two*sum00-br1_g(     2,:)
      br1_g(nt_g,:)=two*sum11-br1_g(ntm1_g,:)
c
      do j=1,nt
        do k=1,np
          br0(j,k)=br0_g(j0_g+j-1,k0_g+k-1)
          br1(j,k)=br1_g(j0_g+j-1,k0_g+k-1)
        enddo
      enddo
!$acc update device(br0,br1)
c
      deallocate(br0_g)
      deallocate(br1_g)
c
      end subroutine
c#######################################################################
      subroutine write_validation_solution
c
c-----------------------------------------------------------------------
c
c ****** Write out the tilted dipole analytic solution for validation.
c ****** This writes out the phi, B, and B computed from phi.
c
c-----------------------------------------------------------------------
c
      use number_types
      use global_dims
      use global_mesh
      use local_dims_r
      use local_mesh_r
      use local_dims_tp
      use local_mesh_tp
      use fields
      use vars
      use mpidefs
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1.0_r_typ
      real(r_typ), parameter :: two=2.0_r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,j,k
      character(256) :: fnamephi,fnamebr,fnamebt,fnamebp
c
c-----------------------------------------------------------------------
c
c ****** Set phi.
c
      do concurrent (k=1:np, j=1:nt, i=1:nr)
        phi(i,j,k)=(-one/rh(i)**2)*(cth(j)*cos(dipole_angle)+
     &                             sth(j)*cph(k)*sin(dipole_angle))
      enddo
c
      call set_boundary_points (phi,one)
      call seam (phi,nr,nt,np)
c
!$acc enter data create(br,bt,bp)
c
c ****** Set Br.
c
      do concurrent (k=1:np, j=1:nt, i=1:nrm1)
        br(i,j,k)=(two/r(i)**3)*(cth(j)*cos(dipole_angle) +
     &                           sth(j)*cph(k)*sin(dipole_angle))
      enddo
c
c ****** Set Bt.
c
      do concurrent (k=1:np, j=1:ntm1, i=1:nr)
        bt(i,j,k)=(one/rh(i)**3)*(st(j)*cos(dipole_angle) -
     &                            ct(j)*cph(k)*sin(dipole_angle))
      enddo
c
c ****** Set Bp.
c
      do concurrent (k=1:npm1, j=1:nt, i=1:nr)
        bp(i,j,k)=(one/rh(i)**3)*(sph(k)*sin(dipole_angle))
      enddo
c
c ****** Due to half mesh trig issues on periodic domain, we need
c ****** to seam these even though seaming is not needed in getb().
c
      call seam (br,nrm,nt,np)
      call seam (bt,nr,ntm,np)
      call seam (bp,nr,nt,npm)
c
c ****** Backup fnames and set new ones so we can reuse write routine.
c
      fnamephi=phifile
      phifile='phi_exact.'//trim(fmt)
      fnamebr=brfile
      brfile='br_exact.'//trim(fmt)
      fnamebt=btfile
      btfile='bt_exact.'//trim(fmt)
      fnamebp=bpfile
      bpfile='bp_exact.'//trim(fmt)
c
      call write_solution
c
c ****** Restore original output names.
c
      phifile=fnamephi
      brfile=fnamebr
      btfile=fnamebt
      bpfile=fnamebp
c
c ****** Reset phi.
c
      phi(:,:,:)=0.
!$acc update device(phi)
!$acc exit data delete(br,bt,bp)
c
      end subroutine
c#######################################################################
      subroutine solve (x,rhs,N,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Solve the implicit equations iteratively.
c
c-----------------------------------------------------------------------
c
c ****** Return IERR=0 if the iteration converges; otherwise,
c ****** IERR is set to a nonzero value.
c
c ****** X is the initial guess at the solution.
c ****** RHS is the right-hand side.
c
c-----------------------------------------------------------------------
c
      use number_types
      use cgcom
      use mpidefs
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: N
      real(r_typ), dimension(N) :: x,rhs
      integer :: ierr
c
c-----------------------------------------------------------------------
c
c ****** Solve the equations using the CG method.
c
      call cgsolve (x,rhs,N,ierr)
c
c ****** Check for convergence.
c
      if (ierr.ne.0) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in SOLVE:'
          write (*,*) '### The field solver did not converge.'
          write (*,*) 'IERR = ',ierr
          write (*,100) ncg,epsn
  100     format (1x,'N = ',i6,' EPSN = ',1pe13.6)
        end if
      else
        if (iamp0) then
          write (*,*)
          write (*,*) '### COMMENT from SOLVE:'
          write (*,*) '### The field solver converged.'
          write (*,*) 'Number of iterations = ',ncg
          write (9,*)
          write (9,*) '### COMMENT from SOLVE:'
          write (9,*) '### The field solver converged.'
          write (9,*) 'Number of iterations = ',ncg
        end if
      end if
c
      return
      end
c#######################################################################
      subroutine cgsolve (x,r,N,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Solve the linear system:
c
c            A * x = b
c
c ****** using the classical Conjugate Gradient method for symmetric
c ****** and positive-definite matrices.
c
c-----------------------------------------------------------------------
c
c ****** On input, X(N) contains a guess at the solution, and
c ****** R(N) contains the right-hand side, b.
c
c ****** On exit, X contains an estimate to the solution, and
c ****** R contains the residual (b-Ax).
c
c ****** IERR=0 indicates that the solution converged to the
c ****** requested accuracy.  Other values indicate that the
c ****** iteration did not converge for the given maximum number
c ****** of iterations.
c
c-----------------------------------------------------------------------
c
      use number_types
      use cgcom
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: N
      real(r_typ), dimension(N) :: x,r
      integer :: ierr,i
c
c-----------------------------------------------------------------------
c
c ****** Scratch space for the CG iteration vectors.
c
      real(r_typ), dimension(N) :: p,ap
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: cgdot
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: bdotb,rdotr,pdotap,alphai,rdotr_old,betai
c
c-----------------------------------------------------------------------
c
      ncg=0
!$acc enter data create(p,ap)
c
c ****** Get the norm of the RHS.
c
      do concurrent (i=1:N)
        p(i)=r(i)
      enddo
c
      call prec_inv (p)
      bdotb=cgdot(r,p,N)
c
c ****** If the RHS is zero, return with a zero solution.
c
      if (bdotb.eq.0.) then
        do concurrent (i=1:N)
          x(i)=0.
        enddo
        epsn=0.
        ierr=0
        return
      end if
c
c-----------------------------------------------------------------------
c ****** Initialization.
c-----------------------------------------------------------------------
c
      call ax (x,ap,N)
c
      do concurrent (i=1:N)
        r(i)=r(i)-ap(i)
        p(i)=r(i)
      enddo
c
c ****** Find the initial error norm.
c
      call prec_inv (p)
      rdotr=cgdot(r,p,N)
c
      call ernorm (bdotb,rdotr,ierr)
      if (ierr.ge.0) return
c
c-----------------------------------------------------------------------
c ****** Main iteration loop.
c-----------------------------------------------------------------------
c
      do
        ncg=ncg+1
c
        call ax (p,ap,N)
c
        pdotap=cgdot(p,ap,N)
        alphai=rdotr/pdotap
c
        do concurrent (i=1:N)
          x(i)=x(i)+alphai*p(i)
          r(i)=r(i)-alphai*ap(i)
          ap(i)=r(i)
        enddo
c
        call prec_inv (ap)
        rdotr_old=rdotr
        rdotr=cgdot(r,ap,N)
c
c ****** Check for convergence.
c
        call ernorm (bdotb,rdotr,ierr)
        if (ierr.ge.0) exit
c
        betai=rdotr/rdotr_old
c
        do concurrent (i=1:N)
          p(i)=betai*p(i)+ap(i)
        enddo
c
      enddo
c
!$acc exit data delete(p,ap)
      end subroutine
c#######################################################################
      subroutine ernorm (bdotb,rdotr,ierr)
c
c-----------------------------------------------------------------------
c
c ****** This subroutine checks if the iterative solver has
c ****** converged or if the maximum allowed number of iterations,
c ****** NCGMAX, has been exceeded.
c
c-----------------------------------------------------------------------
c
c ****** Convergence is deemed to have occurred when:
c ******
c ******     ||R||/||B|| .lt. EPSCG
c ******
c ****** where ||R|| is the norm of the (preconditioned)
c ****** residual, ||B|| is the norm of the (preconditioned)
c ****** RHS, and EPSCG is the specified convergence criterion.
c
c ****** Set IERR=0 if the error is below the error criterion
c ****** (i.e., the solution has converged).
c ****** Set IERR=-1 if the error does not yet meet the error
c ****** criterion and the number of iterations is less than NCGMAX.
c ****** Set IERR=1 if the maximum number of iterations has
c ****** been exceeded without convergence.
c
c-----------------------------------------------------------------------
c
c ****** On input, BDOTB has the dot product of the RHS vector
c ****** with itself, weighted by the preconditioning matrix.
c ****** Similarly, RDOTR has the dot product of the residual vector
c ****** with itself, weighted by the preconditioning matrix.
c ****** This is used to normalize the error estimate.
c
c-----------------------------------------------------------------------
c
      use number_types
      use cgcom
      use mpidefs
      use vars
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: bdotb,rdotr
      integer :: ierr
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: epssq
c
c-----------------------------------------------------------------------
c
      ierr=-1
c
      epssq=rdotr/bdotb
      epsn=sign(one,epssq)*sqrt(abs(epssq))
c
      if (ncghist.gt.0) then
        if (iamp0.and.mod(ncg,50).eq.0) then
          call flush_output_file(outfile,9)
        end if
c
        if (mod(ncg,ncghist).eq.0) then
          if (iamp0) then
            if (ncg.eq.0) then
              write (9,*)
              write (9,*) '### Comment from ERNORM:'
              write (9,*) '### Convergence information:'
            end if
            write (9,100) ncg,epsn
            if (mod(ncg,ncgflush).eq.0) then
              call flush_output_file(outfile,9)
            end if
  100       format (1x,'N = ',i6,' EPSN = ',1pe13.6)
          end if
        end if
      end if
c
c ****** Check for convergence.
c
      if (epsn.lt.epscg) then
        if (ncghist.gt.0) then
          if (iamp0) then
            write (9,*)
            write (9,*) '### Comment from ERNORM:'
            write (9,*) '### The CG solver has converged.'
            write (9,100) ncg,epsn
            call flush_output_file(outfile,9)
          end if
        end if
        ierr=0
      else if (ncg.ge.ncgmax) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in ERNORM:'
          write (*,*) '### Exceeded maximum number of iterations.'
          write (*,*) 'NCGMAX = ',ncgmax
          write (*,*) 'EPSN = ',epsn
          write (9,*)
          write (9,*) '### ERROR in ERNORM:'
          write (9,*) '### Exceeded maximum number of iterations.'
          write (9,*) 'NCGMAX = ',ncgmax
          write (9,*) 'EPSN = ',epsn
          call flush_output_file(outfile,9)
        end if
        ierr=1
      end if
c
      return
      end
c#######################################################################
      subroutine alloc_pot3d_matrix_coefs
c
c-----------------------------------------------------------------------
c
c ****** Allocate the arrays in which the matrix coefficients
c ****** for the pot3d solve are stored.
c
c-----------------------------------------------------------------------
c
      use matrix_storage_pot3d_solve
      use cgcom
      use local_dims_r
      use local_dims_tp
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      allocate (a(2:nrm1,2:ntm1,2:npm1,7))
      a=0.
      allocate (a_i(N))
      a_i=0.
c
      if (ifprec.eq.2) then
        allocate (a_csr(M))
        allocate (lu_csr(M))
        allocate (lu_csr_ja(M))
        allocate (a_csr_ja(M))
        allocate (a_csr_ia(1+N))
        allocate (a_N1(N))
        allocate (a_N2(N))
        allocate (a_csr_d(N))
        allocate (a_csr_dptr(N))
      endif
c
      return
      end
c#######################################################################
      subroutine dealloc_pot3d_matrix_coefs
c
c-----------------------------------------------------------------------
c
c ****** Deallocate the arrays in which the matrix coefficients
c ****** for the pot3d solve are stored.
c
c-----------------------------------------------------------------------
c
      use matrix_storage_pot3d_solve
      use cgcom
      use cusparse_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      deallocate (a)
      deallocate (a_i)
!$acc exit data delete(a,a_i)
c
      if (ifprec.eq.2) then
#ifdef CUSPARSE
        call unload_lusol_cusparse
!$acc exit data delete(a_csr,a_csr_ja,a_csr_ia)
#endif
        deallocate (a_csr)
        deallocate (lu_csr)
        deallocate (lu_csr_ja)
        deallocate (a_csr_ia)
        deallocate (a_csr_ja)
        deallocate (a_csr_d)
        deallocate (a_N1)
        deallocate (a_N2)
        deallocate (a_csr_dptr)
      endif
c
      end subroutine
c#######################################################################
      subroutine load_matrix_pot3d_solve
c
c-----------------------------------------------------------------------
c
c ****** Load the matrix coefficients for the pot3d solve.
c
c-----------------------------------------------------------------------
c
      use number_types
      use matrix_storage_pot3d_solve
      use local_dims_r
      use local_mesh_r
      use local_dims_tp
      use local_mesh_tp
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: i,j,k
c
c-----------------------------------------------------------------------
c
c ****** Set matrix coefs
c
      do k=2,npm1
        do j=2,ntm1
          do i=2,nrm1
c           a*ps(i,j,k-1):
            a(i,j,k,1)=-drh(i)*dth(j)*sth_i(j)*dp_i(k-1)
c           a*ps(i,j-1,k):
            a(i,j,k,2)=-drh(i)*dph(k)*st(j-1)*dt_i(j-1)
c           a*ps(i-1,j,k):
            a(i,j,k,3)=-sth(j)*dth(j)*dph(k)*r2(i-1)*dr_i(i-1)
c           a*ps(i+1,j,k):
            a(i,j,k,5)=-sth(j)*dth(j)*dph(k)*r2(i  )*dr_i(i  )
c           a*ps(i,j+1,k):
            a(i,j,k,6)=-drh(i)*dph(k)*st(j  )*dt_i(j  )
c           a*ps(i,j,k+1):
            a(i,j,k,7)=-drh(i)*dth(j)*sth_i(j)*dp_i(k  )
c
c           a*ps(i,j,k):
            a(i,j,k,4)=-(a(i,j,k,1)+a(i,j,k,2)+a(i,j,k,3)+
     &                   a(i,j,k,5)+a(i,j,k,6)+a(i,j,k,7))
          enddo
        enddo
      enddo
c
      return
      end
c#######################################################################
      subroutine load_preconditioner_pot3d_solve
c
c-----------------------------------------------------------------------
c
c ****** Load the preconditioner for the pot3d solve.
c
c-----------------------------------------------------------------------
c
      use number_types
      use matrix_storage_pot3d_solve
      use cgcom
      use local_dims_r
      use local_dims_tp
      use, intrinsic :: iso_c_binding
      use cusparse_interface
      use timing
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: one=1._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,j,k,icode,ii
c
c-----------------------------------------------------------------------
c
      call timer_on
c
      if (ifprec.eq.0) return
c
      if (ifprec.eq.1) then
c
c ****** Diagonal scaling:
c
        ii=0
        do k=2,npm1
          do j=2,ntm1
            do i=2,nrm1
              ii=ii+1
              a_i(ii)=real(one/a(i,j,k,4),r_typ_pc)
            enddo
          enddo
        enddo
c
      elseif (ifprec.eq.2) then
c
c ****** Convert A matrix into CSR format:
c
        call diacsr (N,M,a,a_offsets,a_csr,a_csr_ja,
     &                     a_csr_ia,a_csr_dptr)
#ifdef CUSPARSE
      cN=N
      cM=M
!$acc enter data copyin(a_csr,a_csr_ja,a_csr_ia)
!$acc host_data use_device(a_csr,a_csr_ja,a_csr_ia)
      call load_lusol_cusparse (C_LOC(a_csr(1)),
     &                          C_LOC(a_csr_ia(1)),
     &                          C_LOC(a_csr_ja(1)),cN,cM)
!$acc end host_data
#else
c
c ****** Overwrite CSR A with preconditioner L and U matrices:
c
c ****** Incomplete LU (ILU)
c
        icode=0
        call ilu0 (N,M,a_csr,a_csr_ja,a_csr_ia,a_csr_dptr,icode)
c
        if (icode.ne.0) then
          print*, '### ERROR IN ILU FORMATION'
        endif
c
c ****** Convert LU stored in A to LU matrix in optimized layout.
c
        call lu2luopt (N,M,lu_csr,a_csr,a_csr_ia,a_csr_ja,lu_csr_ja,
     &                 a_csr_dptr,a_N1,a_N2)
c
c ****** Store inverse of diagonal of LU matrix.
c
        do i=1,N
          a_csr_d(i)=one/a_csr(a_csr_dptr(i))
        enddo
c
#endif
      endif
c
      call timer_off (t_pc_load)
c
      end subroutine
c#######################################################################
      subroutine ilu0 (N,M,A,JA,IA,A_da,icode)
c
c-----------------------------------------------------------
c
c     Set-up routine for ILU(0) preconditioner. This routine
c     computes the L and U factors of the ILU(0) factorization
c     of a general sparse matrix A stored in CSR format with
c     1-based indexing. Since
c     L is unit triangular, the L and U factors can be stored
c     as a single matrix which occupies the same storage as A.
c     New ja and ia arrays are not needed for the LU matrix
c     since the pattern of the LU matrix is identical with
c     that of A.
c
c     Original Author:  Yousef Saad
c            Iterative Methods for Sparse Linear Systems 2nd Ed. pg. 309
c     Modified by R.M. Caplan
c
c-----------------------------------------------------------
c     INPUT:
c     N         : Dimension of matrix
c     A, JA, IA : Sparse matrix in CSR sparse storage format
c     A_da      : Pointers to the diagonal elements in the CSR
c                 data structure luval
c
c     OUTPUT:
c     A     : L/U matrices stored together. On return A,
c             JA, and IA are the combined CSR data structure for
c             the L and U factors.
c     icode : Integer indicating error code on return:
c             (0): Normal return.
c             (k): Encountered a zero pivot at step k.
c------------------------------------------------------------
c
      use number_types
      use number_types_pc
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: JA(M),IA(N+1),A_da(N),iw(N)
      integer :: icode,N,M
      real(r_typ_pc) :: A(M)
c
c-----------------------------------------------------------------------
c
      integer :: i,ik,kj,k,ij,IA_i,IA_ip1m1
      real(r_typ_pc) :: Aik
c
c-----------------------------------------------------------------------
c
      icode=0
c     Initialize scratch index array:
      iw(:)=0
c
      do i=2,N
c       Store index of (i,j) in A in scratch array of iw(j=1:N)
c       This allows lookup given a column index (j) in row (k)
c       to see if the column is in row (i).
        IA_i    =IA(i)
        IA_ip1m1=IA(i+1)-1
c
        do ij=IA_i,IA_ip1m1
          iw(JA(ij))=ij
        enddo
c
c       Loop from first element in row i to 1 less than diagonal elem:
        do ik=IA_i,A_da(i)-1     !IA(i+1) !ik is index of (i,k) in A[]
          k    =JA(ik)           !Actual column index in matrix (k)
          Aik  =A(ik)/A(A_da(k)) !Save Aik for next loop as an optim.
          A(ik)=Aik
c
c         Loop from 1 more than diag elem to last elem in row k:
          do kj=A_da(k)+1,IA(k+1)-1 !kj is index of (k,j) in A[]
c            Get ij location from scratch array (if 0, no ij present)
             ij=iw(JA(kj))
             if (ij .ne. 0) then
               A(ij)=A(ij)-Aik*A(kj)
             endif
          enddo
        enddo
c
        if (A(ik).eq.0) then
          icode=i
          exit
        endif
c
c       Reset scratch index array:
        do ij=IA_i,IA_ip1m1
          iw(JA(ij))=0
        enddo
      enddo
c
      end subroutine
c#######################################################################
      subroutine lu2luopt (N,M,LU,A,IA,JA,LUJA,A_da,N1,N2)
c
c-----------------------------------------------------------------------
c
c ****** Re-order elements of LU matrix in CSR format into custom,
c ****** optimized format for use with lusol().
c ****** (Eventually, this could be merged with the ilu0 and/or diacsr)
c
c-----------------------------------------------------------------------
c
      use number_types
      use number_types_pc
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: N,M
      integer :: JA(M),LUJA(M),IA(N+1),A_da(N)
      integer :: N1(N),N2(N)
      real(r_typ_pc) :: A(M),LU(M)
c
c-----------------------------------------------------------------------
c
      integer :: i,k,ii
c
c-----------------------------------------------------------------------
c
      ii=0
c
      do i=1,N
        do k=IA(i),A_da(i)-1
           ii=ii+1
           LU(ii)=A(k)
           LUJA(ii)=JA(k)
        enddo
c
c       Store k1 and k2 ranges for lusolve:
c
        N1(i)=A_da(i)-1-IA(i)
        N2(i)=IA(i+1)-2-A_da(i)
      enddo
c
      do i=N,1,-1
        do k=A_da(i)+1,IA(i+1)-1
           ii=ii+1
           LU(ii)=A(k)
           LUJA(ii)=JA(k)
        enddo
      enddo
c
      end subroutine
c#######################################################################
      subroutine diacsr (N,M,Adia,ioff,Acsr,JA,IA,Adptr)
c
c-----------------------------------------------------------------------
c
c *** DIACSR_INNER converts a solver matrix in a MAS-style
c     diagonal format to standard compressed sparse row (CSR)
c     including periodic coefficents when nproc_p=1.
c
c     Author of original diacsr: Youcef Saad
c     Modifications for MAS:     RM Caplan
c
c     Input:
c                     N: Size of the matrix (NxN)
c                     M: Number of non-zero entries in matrix
c                        (computed with getM_tc())
c         Adia(IDIAG,N): The matrix in modified "DIA" format
c           ioff(IDIAG): Offsets of the diagonals in A.
c
c     Output:
c            Acsr(M), JA(M), IA(N+1): The matrix A in CSR.
c                           Adptr(N): Pointers to diag elements in A,
c                                     [e.g. A(i,i) == A(Adptr(i))]
c
c-----------------------------------------------------------------------
c
      use number_types
      use number_types_pc
      use local_dims_r
      use local_dims_tp
      use mpidefs
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer, parameter :: IDIAG=7
c
c-----------------------------------------------------------------------
c
      real (r_typ_pc) :: Acsr(M)
      real (r_typ) :: Adia(N,IDIAG)
      integer :: N,M
      integer :: Adptr(N)
      integer :: IA(N+1)
      integer :: JA(M)
      integer :: ioff(IDIAG)
c
c-----------------------------------------------------------------------
c
      integer :: i,j,jj,mi,mj,mk,ko,x
      integer :: ioffok(IDIAG)
c
c-----------------------------------------------------------------------
c
      x=0
c
      IA(1)=1
      ko=1
      i=0
c
      do mk=2,npm1
        do mj=2,ntm1
          do mi=2,nrm1
c ********* Set index of value and column indicies array:
            i=i+1
c
c ********* Do not add coefs that multiply boundaries:
c           For each boundary, there is a sub-set of coefs in the
c           matrix row that should not be added.
c           This makes "local" matrices have no bc info
c
c ********* Reset "i-offset-ok-to-use-coef-jj" array:
c
            ioffok(:)=1
c
            if (mi.eq.2) then
              ioffok(3)=0;
            endif
c
            if (mi.eq.nrm1) then
              ioffok(5)=0;
            endif
c
            if (mj.eq.2) then
              ioffok(2)=0;
            endif
c
            if (mj.eq.ntm1) then
              ioffok(6)=0;
            endif
c
c ********* Eliminate periodic ceofs in the case nproc_p>1
c
            if (nproc_p.gt.1) then
              if (mk.eq.2) then
                ioffok(1)=0
              endif
              if (mk.eq.npm1) then
                ioffok(7)=0
              endif
            endif
c
c ********* To handle periodicity of phi in nproc_p=1 case:
c           We want CSR matrix to be in order so
c           have to sweep three times to avoid sorting:
c
c ********* Add periodic coefs of "right side":
c
            do jj=1,IDIAG
              if (ioffok(jj).eq.1) then
                j=i+ioff(jj)-x
                if (j.gt.N-x) then
                  j=j-N
                  Acsr(ko)=real(Adia(i,jj),r_typ_pc)
                  JA(ko)=j
                  ko=ko+1
                endif
              endif
            enddo
c
c ********* Now do non-periodic coefs:
c
            do jj=1,IDIAG
              if (ioffok(jj).eq.1) then
                j=i+ioff(jj)-x
                if (j.ge.1.and.j.le.N-x) then
c                 Store pointer to diagonal elements in A:
                  if (jj.eq.4) Adptr(i)=ko
                  Acsr(ko)=real(Adia(i,jj),r_typ_pc)
                  JA(ko)=j
                  ko=ko+1
                endif
              endif
            enddo
c
c ********* Now do periodic coefs of "left side":
c
            do jj=1,IDIAG
              if (ioffok(jj).eq.1) then
                j=i+ioff(jj)-x
                if (j.lt.1) then
                  j=N+j
                  Acsr(ko)=real(Adia(i,jj),r_typ_pc)
                  JA(ko)=j
                  ko=ko+1
                endif
              endif
            enddo
c
c ********* Set row offset:
c
            IA(i+1)=ko-x
          enddo
        enddo
      enddo
c
      return
      end
c#######################################################################
      subroutine getM (N, ioff, M)
c
c-----------------------------------------------------------------------
c
c *** This routine computes the number of non-zeros in the
c     solver matrix for use with allocating the matrices.
c     See diacsr() for description of inputs.
c
c     Output:  M  # of nonzeros.
c
c-----------------------------------------------------------------------
c
      use mpidefs
      use local_dims_r
      use local_dims_tp
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer, parameter :: IDIAG=7
      integer :: N,M,i,j,jj,ko,mi,mj,mk,x
      integer :: ioff(IDIAG)
      integer :: ioffok(IDIAG)
c
      x=0
c
      ko=1
      i=0
c
      do mk=2,npm1
        do mj=2,ntm1
          do mi=2,nrm1
c
            ioffok(:)=1
c
            if (mi.eq.2) then
              ioffok(3)=0;
            endif
c
            if (mi.eq.nrm1) then
              ioffok(5)=0;
            endif
c
            if (mj.eq.2) then
              ioffok(2)=0;
            endif
c
            if (mj.eq.ntm1) then
              ioffok(6)=0;
            endif
c
c ********* Eliminate periodic ceofs in the case nproc_p>1
c
            if (nproc_p.gt.1) then
              if (mk.eq.2) then
                ioffok(1)=0
              endif
              if (mk.eq.npm1) then
                ioffok(7)=0
              endif
            endif
c
            do jj=1,IDIAG
              if (ioffok(jj).eq.1) then
                j=i+ioff(jj)-x
                if (j.gt.N-x) then
                  ko=ko+1
                endif
              endif
            enddo
c
            do jj=1,IDIAG
              if (ioffok(jj).eq.1) then
                j=i+ioff(jj)-x
                if (j.ge.1.and.j.le.N-x) then
                  ko=ko+1
                endif
              endif
            enddo
c
            do jj=1,IDIAG
              if (ioffok(jj).eq.1) then
                j=i+ioff(jj)-x
                if (j.lt.1) then
                  ko=ko+1
                endif
              endif
            enddo
          enddo
        enddo
      enddo
c
c *** Save number of non-zeros of matrix:
c
      M=ko-1
c
      return
      end
c#######################################################################
      subroutine ax (x,y,N)
c
c-----------------------------------------------------------------------
c
c ****** Set y = A * x.
c
c-----------------------------------------------------------------------
c
      use number_types
      use local_dims_r
      use local_dims_tp
      use fields, ONLY : x_ax
      use timing
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: zero=0._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: N
      real(r_typ), dimension(N) :: x,y
c
c-----------------------------------------------------------------------
c
c ****** Expand X array to allow for boundary and seam values.
c
      call unpack_scalar (x_ax,x)
c
c ****** Set the boundary values of X.
c
      call set_boundary_points (x_ax,zero)
c
c ****** Seam along edges between processors.
c
      call seam (x_ax,nr,nt,np)
c
c ****** Get the matrix-vector product.
c
      call timer_on
c
      call delsq (x_ax,y)
c
      call timer_off (t_ax)
c
      end subroutine
c#######################################################################
      subroutine prec_inv (x)
c
c-----------------------------------------------------------------------
c
c ****** Apply preconditioner: x := M(inv) * x.
c
c-----------------------------------------------------------------------
c
      use number_types
      use number_types_pc
      use cgcom
      use solve_params
      use matrix_storage_pot3d_solve
      use, intrinsic :: iso_c_binding
      use cusparse_interface
      use timing
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(N) :: x
      real(r_typ_pc), dimension(N) :: x_32
      integer :: i
c
c-----------------------------------------------------------------------
c
      if (ifprec.eq.0) return
c
      call timer_on
c
      if (ifprec.eq.1) then
c
c ****** Point-Jacobi (diagonal scaling):
c
        do concurrent (i=1:N)
          x(i)=a_i(i)*x(i)
        enddo
c
      elseif (ifprec.eq.2) then
c
c ****** ILU0 Partial-Block-Jacobi:
c
#ifdef CUSPARSE
!$acc host_data use_device(x)
        call lusol_cusparse(C_LOC(x(1)))
!$acc end host_data
c
#else
!$acc update self(x)
c
c ****** Convert input array to single precision.
c
        do i=1,N
          x_32(i) = real(x(i),r_typ_pc)
        enddo
c
        call lusol (N,M,x_32,lu_csr,lu_csr_ja,a_N1,a_N2,a_csr_d)
c
c ****** Convert result back to double precision.
c
        do i=1,N
          x(i) = real(x_32(i),r_typ)
        enddo
c
!$acc update device(x)
#endif
      endif
c
      call timer_off (t_pc)
c
      end subroutine
c#######################################################################
      subroutine lusol (N,M,x,LU,LU_ja,N1,N2,LUd_i)
c
c-----------------------------------------------------------
c
c     Performs a forward and a backward solve for the sparse system
c     (LU) x=y  where LU is in an optimized custom CSR format
c                                              (see lu2luopt())
c
c     For use where LU is an ILU or SSOR/SGS factorization.
c
c     Author of original lusol: Yousef Saad
c           Iterative Methods for Sparse Linear Systems 2nd Ed. pg. 299
c
c     Modified by RM Caplan to include optimized memory access
c     as described in
c     B. Smith, H. Zhang  Inter. J. of High Perf. Comp. Appl.
c     Vol. 25 #4 pg. 386-391 (2011)
c
c-----------------------------------------------------------
c     PARAMETERS:
c     N     : Dimension of problem
c     x     : At input, x is rhs (y), at output x is the solution.
c     LU    : Values of the LU matrix. L and U are stored together in
c             order of access in this routine.
c     LU_ja : Column indices of elements in LU.
c     N1    : Row-start indicies in original CSR LU.
c     N2    : Indices of diagonal elements in orig CSR LU
c     LUd_i : Inverse diagonal elements of U
c------------------------------------------------------------
c
      use number_types
      use number_types_pc
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ_pc) :: x(N),LUd_i(N),LU(M)
      integer :: N1(N),N2(N),LU_ja(M)
      integer :: N,M
c
c-----------------------------------------------------------------------
c
      integer :: i,k,k1,k2
c
c-----------------------------------------------------------------------
c
c ****** FORWARD SOLVE: Solve L x'=y
c
      k2=0
      do i=1,N
c       Compute x(i) := x(i) - sum L(i,j) * x(j)
        k1=k2+1
        k2=k1+N1(i)
        do k=k1,k2
          x(i)=x(i)-LU(k)*x(LU_ja(k))
        enddo
c       Diagonal is always 1 for L so no division here is nessesary.
      enddo
c
c ****** BACKWARD SOLVE: Solve U x=x'
c
      do i=N,1,-1
c       Compute x(i) := x(i) - sum U(i,j) * x(j)
        k1=k2+1
        k2=k1+N2(i)
        do k=k1,k2
          x(i)=x(i)-LU(k)*x(LU_ja(k))
        enddo
c       Compute x(i) := x(i) / U(i,i)
        x(i)=x(i)*LUd_i(i)
      enddo
c
      end subroutine
c#######################################################################
      subroutine unpack_scalar (s,x)
c
c-----------------------------------------------------------------------
c
c ****** Unpack the scalar x into
c ****** three-dimensional array s leaving room for boundaries.
c
c-----------------------------------------------------------------------
c
      use number_types
      use local_dims_r
      use local_dims_tp
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(nr,nt,np) :: s
      real(r_typ), dimension(2:nrm1,2:ntm1,2:npm1) :: x
c
c-----------------------------------------------------------------------
c
      integer :: i,j,k
c
c-----------------------------------------------------------------------
c
      do concurrent (k=2:npm1, j=2:ntm1, i=2:nrm1)
        s(i,j,k)=x(i,j,k)
      enddo
c
      return
      end
c#######################################################################
      subroutine delsq (x,y)
c
c-----------------------------------------------------------------------
c
c ****** Set Y = - (dV * del-squared X) at the internal points.
c
c-----------------------------------------------------------------------
c
      use number_types
      use local_dims_r
      use local_dims_tp
      use matrix_storage_pot3d_solve
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(nr,nt,np) :: x
      real(r_typ), dimension(2:nrm1,2:ntm1,2:npm1) :: y
c
c-----------------------------------------------------------------------
c
      integer :: i,j,k
c
c-----------------------------------------------------------------------
c
      do concurrent (k=2:npm1, j=2:ntm1, i=2:nrm1)
        y(i,j,k)=a(i,j,k,1)*x(i  ,j  ,k-1)
     &          +a(i,j,k,2)*x(i  ,j-1,k  )
     &          +a(i,j,k,3)*x(i-1,j  ,k  )
     &          +a(i,j,k,4)*x(i  ,j  ,k  )
     &          +a(i,j,k,5)*x(i+1,j  ,k  )
     &          +a(i,j,k,6)*x(i  ,j+1,k  )
     &          +a(i,j,k,7)*x(i  ,j  ,k+1)
      enddo
c
      return
      end
c#######################################################################
      subroutine set_boundary_points (x,vmask)
c
c-----------------------------------------------------------------------
c
c ****** Set boundary points of X at the physical boundaries.
c
c-----------------------------------------------------------------------
c
      use number_types
      use global_mesh
      use local_dims_r
      use local_mesh_r
      use local_dims_tp
      use local_mesh_tp
      use fields
      use solve_params
      use vars, ONLY : validation_run
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(nr,nt,np) :: x
      real(r_typ) :: vmask
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: two=2._r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,j,k
c
c-----------------------------------------------------------------------
c
c ****** Set X at the radial boundaries.
c
      if (rb0) then
        do concurrent (k=2:npm1, j=2:ntm1)
          x( 1,j,k)=x(2,j,k)-vmask*br0(j,k)*dr1
        enddo
      end if
c
      if (rb1) then
        if (validation_run) then
          do concurrent (k=2:npm1, j=2:ntm1)
            x(nr,j,k)=x(nrm1,j,k)+vmask*br1(j,k)*drn
          enddo
        else
          do concurrent (k=2:npm1, j=2:ntm1)
            x(nr,j,k)=pm_r1*x(nrm1,j,k)
          enddo
        end if
      end if
c
c ****** If this processor does not contain any points at the
c ****** pole, return.
c
      if (.not.(tb0.or.tb1)) return
c
c ****** Get the m=0 component of X at the poles.
c
      if (tb0) then
        do concurrent (i=1:nr)
          sum0(i)=0
        enddo
!$acc parallel loop collapse(2) default(present) async(1)
        do k=2,npm1
          do i=1,nr
!$acc atomic
            sum0(i)=sum0(i)+x(i,2,k)*dph(k)*pl_i
          enddo
        enddo
      end if
c
      if (tb1) then
        do concurrent (i=1:nr)
          sum1(i)=0
        enddo
!$acc parallel loop collapse(2) default(present) async(2)
        do k=2,npm1
          do i=1,nr
!$acc atomic
            sum1(i)=sum1(i)+x(i,ntm1,k)*dph(k)*pl_i
          enddo
        enddo
      end if
!$acc wait(1,2)
c
c ****** Sum over all processors.
c
      call sum_over_phi (nr,sum0,sum1)
c
c ****** Set X to have only an m=0 component at the poles.
c
      if (tb0) then
        do concurrent (k=2:npm1, i=1:nr)
          x(i,1,k)=two*sum0(i)-x(i,2,k)
        enddo
      end if
c
      if (tb1) then
        do concurrent (k=2:npm1, i=1:nr)
          x(i,nt,k)=two*sum1(i)-x(i,ntm1,k)
        enddo
      end if
c
      return
      end subroutine
c#######################################################################
      subroutine sum_over_phi (n,a0,a1)
c
c-----------------------------------------------------------------------
c
c ****** Sum the contribution over all processors in the phi
c ****** dimension (only for processors with points on the poles).
c
c ****** The sum is performed for all N points in the vectors
c ****** SUM0(N) and SUM1(N), at the North and South pole,
c ****** respectively.
c
c-----------------------------------------------------------------------
c
      use number_types
      use local_dims_tp
      use mpidefs
      use timing
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: n
      real(r_typ), dimension(n) :: a0,a1
c
c-----------------------------------------------------------------------
c
c ****** MPI error return.
c
      integer :: ierr
c
c-----------------------------------------------------------------------
c
      call timer_on
c
!$acc host_data use_device(a0,a1)
      if (tb0) then
        call MPI_Allreduce (MPI_IN_PLACE,a0,n,ntype_real,
     &                        MPI_SUM,comm_phi,ierr)
      end if
c
      if (tb1) then
        call MPI_Allreduce (MPI_IN_PLACE,a1,n,ntype_real,
     &                        MPI_SUM,comm_phi,ierr)
      end if
!$acc end host_data
c
      call timer_off (c_sumphi)
c
      return
      end
c#######################################################################
      subroutine zero_boundary_points (x)
c
c-----------------------------------------------------------------------
c
c ****** Set the boundary points at the physical boundaries
c ****** of X to zero.
c
c-----------------------------------------------------------------------
c
      use number_types
      use local_dims_r
      use local_dims_tp
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(nr,nt,np) :: x
c
c-----------------------------------------------------------------------
c
      if (rb0) x( 1,:,:)=0.
      if (rb1) x(nr,:,:)=0.
      if (tb0) x(:, 1,:)=0.
      if (tb1) x(:,nt,:)=0.
c
      return
      end
c#######################################################################
      function cgdot (x,y,N)
c
c-----------------------------------------------------------------------
c
c ****** Get the dot product of the vectors X and Y.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: N,i
      real(r_typ) :: cgdot
      real(r_typ), dimension(N) :: x,y
c
c-----------------------------------------------------------------------
c
      cgdot=0.
c
!$acc parallel loop present(x,y) reduction(+:cgdot)
      do i=1,N
        cgdot=cgdot+x(i)*y(i)
      enddo
c
c ****** Sum over all the processors.
c
      call global_sum (cgdot)
c
      return
      end function
c#######################################################################
      subroutine global_sum (x)
c
c-----------------------------------------------------------------------
c
c ****** Overwrite X by the its sum over all processors.
c
c-----------------------------------------------------------------------
c
      use number_types
      use mpidefs
      use timing
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: x
c
c-----------------------------------------------------------------------
c
c ****** MPI error return.
c
      integer :: ierr
c
c-----------------------------------------------------------------------
c
      call timer_on
c
c ****** Take the sum over all the processors.
c
      call MPI_Allreduce (MPI_IN_PLACE,x,1,ntype_real,
     &                    MPI_SUM,comm_all,ierr)
c
      call timer_off (c_cgdot)
c
      end subroutine
c#######################################################################
      subroutine seam (a,n1,n2,n3)
c
c-----------------------------------------------------------------------
c
c ****** Seam the boundary points of 3D (r,t,p) array A between
c ****** adjacent processors.
c
c ****** This routine assumes that there is a two-point
c ****** overlap between processors in each dimension.
c
c-----------------------------------------------------------------------
c
c ****** This version uses non-blocking MPI sends and receives
c ****** whenever possible in order to overlap communications.
c
c-----------------------------------------------------------------------
c
      use number_types
      use mpidefs
      use timing
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(n1,n2,n3) :: a
c
c-----------------------------------------------------------------------
c
      real(r_typ), dimension(n2,n3) :: sbuf11,rbuf11
      real(r_typ), dimension(n2,n3) :: sbuf12,rbuf12
      real(r_typ), dimension(n1,n3) :: sbuf21,rbuf21
      real(r_typ), dimension(n1,n3) :: sbuf22,rbuf22
c
c-----------------------------------------------------------------------
c
c ****** MPI error return.
c
      integer :: ierr
c
c ****** MPI tag for MPI_ISEND and MPI_IRECV (not tagged).
c
      integer :: tag=0
c
c-----------------------------------------------------------------------
c
      integer :: lbuf,i,j
      integer :: n1,n2,n3
      integer :: reqs(4)
c
c-----------------------------------------------------------------------
c
      call timer_on
c
c ****** Seam the third (periodic) dimension.
c ****** Since halo data is stride-1, no need for buffers.
c
      lbuf=n1*n2
c
!$acc host_data use_device(a)
      call MPI_Isend (a(:,:,n3-1),lbuf,ntype_real,iproc_pp,tag,
     &                comm_all,reqs(1),ierr)
c
      call MPI_Isend (a(:,:,   2),lbuf,ntype_real,iproc_pm,tag,
     &                comm_all,reqs(2),ierr)
c
      call MPI_Irecv (a(:,:, 1),lbuf,ntype_real,iproc_pm,tag,
     &                comm_all,reqs(3),ierr)
c
      call MPI_Irecv (a(:,:,n3),lbuf,ntype_real,iproc_pp,tag,
     &                comm_all,reqs(4),ierr)
c
      call MPI_Waitall (4,reqs,MPI_STATUSES_IGNORE,ierr)
!$acc end host_data
c
c ****** Seam the first dimension.
c
      if (nproc_r.gt.1) then
c
!$acc enter data create(sbuf11,sbuf12,rbuf11,rbuf12)
c
        lbuf=n2*n3
c
        do concurrent (j=1:n3, i=1:n2)
          sbuf11(i,j)=a(n1-1,i,j)
          sbuf12(i,j)=a(   2,i,j)
        enddo
c
!$acc host_data use_device(sbuf11,sbuf12,rbuf11,rbuf12)
        call MPI_Isend (sbuf11,lbuf,ntype_real,iproc_rp,tag,
     &                  comm_all,reqs(1),ierr)
c
        call MPI_Isend (sbuf12,lbuf,ntype_real,iproc_rm,tag,
     &                  comm_all,reqs(2),ierr)
c
        call MPI_Irecv (rbuf11,lbuf,ntype_real,iproc_rm,tag,
     &                  comm_all,reqs(3),ierr)
c
        call MPI_Irecv (rbuf12,lbuf,ntype_real,iproc_rp,tag,
     &                  comm_all,reqs(4),ierr)
c
        call MPI_Waitall (4,reqs,MPI_STATUSES_IGNORE,ierr)
!$acc end host_data
c
        if (iproc_rm.ne.MPI_PROC_NULL) then
          do concurrent (j=1:n3, i=1:n2)
            a( 1,i,j)=rbuf11(i,j)
          enddo
        end if
c
        if (iproc_rp.ne.MPI_PROC_NULL) then
          do concurrent (j=1:n3, i=1:n2)
            a(n1,i,j)=rbuf12(i,j)
          enddo
        end if
c
!$acc exit data delete(sbuf11,sbuf12,rbuf11,rbuf12)
      end if
c
c ****** Seam the second dimension.
c
      if (nproc_t.gt.1) then
c
!$acc enter data create(sbuf21,sbuf22,rbuf21,rbuf22)
c
        lbuf=n1*n3
c
        do concurrent (j=1:n3, i=1:n1)
          sbuf21(i,j)=a(i,n2-1,j)
          sbuf22(i,j)=a(i,   2,j)
        enddo
c
!$acc host_data use_device(sbuf21,sbuf22,rbuf21,rbuf22)
        call MPI_Isend (sbuf21,lbuf,ntype_real,iproc_tp,tag,
     &                  comm_all,reqs(1),ierr)
c
        call MPI_Isend (sbuf22,lbuf,ntype_real,iproc_tm,tag,
     &                  comm_all,reqs(2),ierr)
c
        call MPI_Irecv (rbuf21,lbuf,ntype_real,iproc_tm,tag,
     &                  comm_all,reqs(3),ierr)
c
        call MPI_Irecv (rbuf22,lbuf,ntype_real,iproc_tp,tag,
     &                  comm_all,reqs(4),ierr)
c
        call MPI_Waitall (4,reqs,MPI_STATUSES_IGNORE,ierr)
!$acc end host_data
c
        if (iproc_tm.ne.MPI_PROC_NULL) then
          do concurrent (j=1:n3, i=1:n1)
            a(i, 1,j)=rbuf21(i,j)
          enddo
        end if
c
        if (iproc_tp.ne.MPI_PROC_NULL) then
          do concurrent (j=1:n3, i=1:n1)
            a(i,n2,j)=rbuf22(i,j)
          enddo
        end if
c
!$acc exit data delete(sbuf21,sbuf22,rbuf21,rbuf22)
      end if
c
      call timer_off (c_seam)
c
      end subroutine
c#######################################################################
      subroutine write_solution
c
c-----------------------------------------------------------------------
c
c ****** Write the global solution.
c
c-----------------------------------------------------------------------
c
      use number_types
      use global_dims
      use global_mesh
      use fields
      use vars
      use solve_params
      use mpidefs
      use decomposition
      use assemble_array_interface
      use timing
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the global arrays.
c ****** These arrays are only allocated on processor IPROC0.
c
      real(r_typ), dimension(:,:,:), allocatable :: phi_g
      real(r_typ), dimension(:,:,:), allocatable :: br_g
      real(r_typ), dimension(:,:,:), allocatable :: bt_g
      real(r_typ), dimension(:,:,:), allocatable :: bp_g
c
c-----------------------------------------------------------------------
c
      integer :: ierr
      character(256) :: fname
c
c-----------------------------------------------------------------------
c
      call timer_on
c
c ****** Potential.
c
      if (phifile.ne.'') then
!$acc update self(phi)
c
c ****** Allocate the global array PHI_G (on processor IPROC0).
c
        if (iamp0) then
          allocate (phi_g(nr_g,nt_g,np_g))
        else
          allocate (phi_g(1,1,1))
        end if
c
c ****** Assemble the global PHI array.
c
        call assemble_array (map_rh,map_th,map_ph,phi,phi_g)
c
        fname=phifile
c
c ****** Write out the potential to a file.
c
        if (iamp0) then
            write (*,*)
            write (*,*) '### COMMENT from WRITE_SOLUTION:'
            write (*,*)
            write (*,*) 'Writing the potential to file: ',trim(fname)
            call wrhdf_3d (fname,.true.,nr_g,nt_g,np_g,
     &                     phi_g,rh_g,th_g,ph_g,hdf32,ierr)
        end if
c
        deallocate (phi_g)
c
      end if
c
c ****** Br.
c
      if (brfile.ne.'') then
!$acc update self(br)
c
        fname=brfile
c
        if (iamp0) then
          allocate (br_g(nrm1_g,nt_g,np_g))
        else
          allocate (br_g(1,1,1))
        end if
c
c ****** Assemble the global PHI array.
c
        call assemble_array (map_rm,map_th,map_ph,br,br_g)
c
        if (iamp0) then
          write (*,*)
          write (*,*) '### COMMENT from WRITE_SOLUTION:'
          write (*,*)
          write (*,*) 'Writing Br to file: ',trim(fname)
          call wrhdf_3d (fname,.true.,nrm1_g,nt_g,np_g,
     &                   br_g,r_g,th_g,ph_g,hdf32,ierr)
        end if
c
        deallocate (br_g)
c
      end if
c
c ****** Bt.
c
      if (btfile.ne.'') then
!$acc update self(bt)
c
        fname=btfile
c
        if (iamp0) then
          allocate (bt_g(nr_g,ntm1_g,np_g))
        else
          allocate (bt_g(1,1,1))
        end if
c
c ****** Assemble the global PHI array.
c
        call assemble_array (map_rh,map_tm,map_ph,bt,bt_g)
c
        if (iamp0) then
          write (*,*)
          write (*,*) '### COMMENT from WRITE_SOLUTION:'
          write (*,*)
          write (*,*) 'Writing Bt to file: ',trim(fname)
          call wrhdf_3d (fname,.true.,nr_g,ntm1_g,np_g,
     &                   bt_g,rh_g,t_g,ph_g,hdf32,ierr)
c
        end if
c
        deallocate (bt_g)
c
      end if
c
c ****** Bp.
c
      if (bpfile.ne.'') then
!$acc update self(bp)
c
        fname=bpfile
c
        if (iamp0) then
          allocate (bp_g(nr_g,nt_g,npm1_g))
        else
          allocate (bp_g(1,1,1))
        end if
c
c ****** Assemble the global PHI array.
c
        call assemble_array (map_rh,map_th,map_pm,bp,bp_g)
c
        if (iamp0) then
          write (*,*)
          write (*,*) '### COMMENT from WRITE_SOLUTION:'
          write (*,*)
          write (*,*) 'Writing Bp to file: ',trim(fname)
          call wrhdf_3d (fname,.true.,nr_g,nt_g,npm1_g,
     &                   bp_g,rh_g,th_g,p_g,hdf32,ierr)
c
        end if
c
        deallocate (bp_g)
c
      end if
c
      call timer_off (t_io)
c
      end subroutine
c#######################################################################
      subroutine getb
c
c-----------------------------------------------------------------------
c
c ****** Calculate B from grad-phi.
c
c-----------------------------------------------------------------------
c
      use number_types
      use global_dims
      use global_mesh
      use vars
      use fields
      use local_dims_r
      use local_dims_tp
      use local_mesh_r
      use local_mesh_tp
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: i,j,k
c
c-----------------------------------------------------------------------
c
!$acc enter data create(br,bt,bp)
c
c ****** Get Br.
c
      do concurrent (k=1:np, j=1:nt, i=1:nrm1)
        br(i,j,k)=(phi(i+1,j,k)-phi(i,j,k))/dr(i)
      enddo
c
c ****** Get Bt.
c
      do concurrent (k=1:np, j=1:ntm1, i=1:nr)
        bt(i,j,k)=(phi(i,j+1,k)-phi(i,j,k))/(rh(i)*dt(j))
      enddo
c
c ****** Get Bp.
c
      do concurrent (k=1:npm1, j=1:nt, i=1:nr)
        bp(i,j,k)=(phi(i,j,k+1)-phi(i,j,k))/(rh(i)*sth(j)*dp(k))
      enddo
c
      end subroutine
c#######################################################################
      subroutine magnetic_energy
c
c-----------------------------------------------------------------------
c
c ****** Calculate magnetic energy from B.
c
c-----------------------------------------------------------------------
c
      use number_types
      use global_dims
      use global_mesh
      use vars
      use fields
      use mpidefs
      use local_dims_r
      use local_dims_tp
      use local_mesh_r
      use local_mesh_tp
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: half=.5_r_typ
c
c-----------------------------------------------------------------------
c
      integer :: i,j,k,ierr
      real(r_typ) :: brav,btav,bpav,dv
      real(r_typ) :: wr,wt,wp
      real(r_typ), dimension(3) :: w
      character(32) :: fmtstr
c
c-----------------------------------------------------------------------
c
      if (hdf32) then
        fmtstr="(A, ES14.8)"
      else
        fmtstr="(A,ES22.16)"
      end if
c
      wr=0.
      wt=0.
      wp=0.
!$acc parallel loop default(present) collapse(3) reduction(+:wr,wt,wp)
      do k=2,npm1
        do j=2,ntm1
          do i=2,nrm1
            dv=rh(i)**2*drh(i)*dth(j)*sth(j)*dph(k)
            brav=half*(br(i,j,k)+br(i-1,j,k))
            btav=half*(bt(i,j,k)+bt(i,j-1,k))
            bpav=half*(bp(i,j,k)+bp(i,j,k-1))
            wr=wr+half*brav**2*dv
            wt=wt+half*btav**2*dv
            wp=wp+half*bpav**2*dv
          enddo
        enddo
      enddo
c
c ****** Sum up all processors into final values and print.
c
      w(1)=wr
      w(2)=wt
      w(3)=wp
      call MPI_Allreduce(MPI_IN_PLACE,w,3,ntype_real,
     &                   MPI_SUM,comm_all,ierr)
c
      if (iamp0) then
        write (*,*)
        write (*,*) '### COMMENT from GETB:'
        write (*,*) '### Magnetic energy diagnostic:'
        write (*,*)
        write (*,trim(fmtstr)) 'Energy in Br**2 = ',w(1)
        write (*,trim(fmtstr)) 'Energy in Bt**2 = ',w(2)
        write (*,trim(fmtstr)) 'Energy in Bp**2 = ',w(3)
        write (*,trim(fmtstr)) 'Magnetic energy = ',SUM(w)
        write (9,*)
        write (9,*) '### COMMENT from GETB:'
        write (9,*) '### Magnetic energy diagnostic:'
        write (9,*)
        write (9,trim(fmtstr)) 'Energy in Br**2 = ',w(1)
        write (9,trim(fmtstr)) 'Energy in Bt**2 = ',w(2)
        write (9,trim(fmtstr)) 'Energy in Bp**2 = ',w(3)
        write (9,trim(fmtstr)) 'Magnetic energy = ',SUM(w)
      end if
c
      end subroutine
c#######################################################################
      subroutine assemble_array (map_r,map_t,map_p,a,a_g)
c
c-----------------------------------------------------------------------
c
c ****** Assemble a global array (into A_G) on processor IPROC0 by
c ****** fetching the local sections (A) from all the processors.
c
c-----------------------------------------------------------------------
c
      use number_types
      use decomposition
      use mpidefs
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(map_struct), dimension(0:nproc-1) :: map_r,map_t,map_p
      real(r_typ), dimension(:,:,:) :: a,a_g
c
c-----------------------------------------------------------------------
c
c ****** Storage for the buffers.
c
      integer :: lbuf,lsbuf
      real(r_typ), dimension(:), allocatable :: sbuf
      real(r_typ), dimension(:), allocatable :: rbuf
c
c-----------------------------------------------------------------------
c
      integer :: tag=0
      integer :: irank,l1,l2,l3,i,j,k,ii
      integer :: i0,j0,k0,i1,j1,k1
      integer :: i0g,j0g,k0g
      integer :: ierr
c
c-----------------------------------------------------------------------
c
      l1=map_r(iproc)%n
      l2=map_t(iproc)%n
      l3=map_p(iproc)%n
      lsbuf=l1*l2*l3
c
      i0=map_r(iproc)%i0
      i1=map_r(iproc)%i1
      j0=map_t(iproc)%i0
      j1=map_t(iproc)%i1
      k0=map_p(iproc)%i0
      k1=map_p(iproc)%i1
c
c ***** Extract 1D array of non-overlapping values from local array.
c
      allocate (sbuf(lsbuf))
c
      sbuf=reshape(a(i0:i1,j0:j1,k0:k1),(/lsbuf/))
c
c ****** If proc0, recieve/store local arrays into global array.
c
      if (iamp0) then
        do irank=0,nproc-1
c
          l1=map_r(irank)%n
          l2=map_t(irank)%n
          l3=map_p(irank)%n
          lbuf=l1*l2*l3
c
          i0g=map_r(irank)%offset
          j0g=map_t(irank)%offset
          k0g=map_p(irank)%offset
c
c ****** If proc0 is the current rank in loop, simply copy local array.
          if (iproc==irank) then
            do k=1,l3
              do j=1,l2
                do i=1,l1
                  ii=l2*l1*(k-1)+l1*(j-1)+i
                  a_g(i0g+i-1,j0g+j-1,k0g+k-1)=sbuf(ii)
                enddo
              enddo
            enddo
c ****** Otherwise recieve data:
          else
            allocate (rbuf(lbuf))
            call MPI_Recv (rbuf,lbuf,ntype_real,irank,tag,
     &                     comm_all,MPI_STATUS_IGNORE,ierr)
            do k=1,l3
              do j=1,l2
                do i=1,l1
                  ii=l2*l1*(k-1)+l1*(j-1)+i
                  a_g(i0g+i-1,j0g+j-1,k0g+k-1)=rbuf(ii)
                enddo
              enddo
            enddo
            deallocate(rbuf)
          end if
        enddo
      else
c
c ****** Send local array to iproc0.
c
      call MPI_Ssend (sbuf,lsbuf,ntype_real,iproc0,tag,comm_all,ierr)
c
      end if
      deallocate (sbuf)
c
      return
      end
c#######################################################################
      subroutine timer_on
c
c-----------------------------------------------------------------------
c
c ****** Push an entry onto the timing stack and initialize
c ****** a timing event.
c
c-----------------------------------------------------------------------
c
c ****** This routine can be called in a nested way to measure
c ****** multiple timing events.  Calls to TIMER_ON and TIMER_OFF
c ****** need to be nested like do-loops in FORTRAN.
c
c-----------------------------------------------------------------------
c
      use mpidefs
      use timer
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      if (istack.ge.nstack) then
        write (*,*)
        write (*,*) '### WARNING from TIMER_ON:'
        write (*,*) '### Timing stack depth exceeded.'
        write (*,*) 'This may signal an incorrect nesting of '//
     &              'TIMER_ON/TIMER_OFF calls.'
        write (*,*) 'Timing information will not be valid.'
        return
      else
        istack=istack+1
      end if
c
      tstart(istack)=MPI_Wtime()
c
      return
      end
c#######################################################################
      subroutine timer_off (tused)
c
c-----------------------------------------------------------------------
c
c ****** Increment the CPU time used since the call to TIMER_ON
c ****** in variable TUSED, and pop an entry off the timing
c ****** stack.
c
c-----------------------------------------------------------------------
c
c ****** This routine can be called in a nested way to measure
c ****** multiple timing events.  Calls to TIMER_ON and TIMER_OFF
c ****** need to be nested like do-loops in FORTRAN.
c
c-----------------------------------------------------------------------
c
      use number_types
      use mpidefs
      use timer
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: tused
c
c-----------------------------------------------------------------------
c
      if (istack.le.0) then
        write (*,*)
        write (*,*) '### WARNING from TIMER_OFF:'
        write (*,*) '### Timing stack cannot be popped.'
        write (*,*) 'This may signal an incorrect nesting of '//
     &              'TIMER_ON/TIMER_OFF calls.'
        write (*,*) 'Timing information will not be valid.'
        return
      else
        istack=istack-1
      end if
c
      tused=tused+MPI_Wtime()-tstart(istack+1)
c
      return
      end
c#######################################################################
      subroutine write_timing
c
c-----------------------------------------------------------------------
c
c ****** Write out the timing info.
c
c-----------------------------------------------------------------------
c
      use number_types
      use mpidefs
      use timing
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c ****** Timing buffers.
c
      integer, parameter :: lbuf=10
      real(r_typ), dimension(lbuf) :: sbuf
      real(r_typ), dimension(lbuf,0:nproc-1) :: tbuf
c
c ****** Timing statistics.
c
      real(r_typ), dimension(lbuf) :: tmin,tmax,tavg,tsdev
c
c-----------------------------------------------------------------------
c
      integer :: ierr,irank
      real(r_typ) :: t_tot_avg,c_tot_avg,c_tot
c
      character(80) :: tfile='timing.out'
c
c-----------------------------------------------------------------------
c
c ****** Gather the timing information for all processors into TBUF.
c
      sbuf(1)=t_solve
      sbuf(2)=t_startup
      sbuf(3)=t_io
      sbuf(4)=t_ax
      sbuf(5)=t_pc_load
      sbuf(6)=t_pc
      sbuf(7)=c_seam
      sbuf(8)=c_cgdot
      sbuf(9)=c_sumphi
      sbuf(10)=t_wall
c
      call MPI_Allgather (sbuf,lbuf,ntype_real,
     &                    tbuf,lbuf,ntype_real,comm_all,ierr)
c
c ****** Calculate the timing statistics.
c
      tavg=sum(tbuf,dim=2)/nproc
      tmin=minval(tbuf,dim=2)
      tmax=maxval(tbuf,dim=2)
c
      tsdev(:)=0.
      do irank=0,nproc-1
        tsdev(:)=tsdev(:)+(tbuf(:,irank)-tavg(:))**2
      enddo
      tsdev(:)=sqrt(tsdev(:)/nproc)
c
      t_tot_avg=tavg(10)
      c_tot_avg=tavg(7)+tavg(8)+tavg(9)
c
      if (iamp0) then
c
        call ffopen (1,tfile,'rw',ierr)
c
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### WARNING from WRITE_TIMING:'
          write (*,*) '### Could not create the timing file.'
          write (*,*) 'File name: ',trim(tfile)
        end if
c
        do irank=0,nproc-1
          c_tot=tbuf(7,irank)+tbuf(8,irank)+tbuf(9,irank)
          write (1,*)
          write (1,100)
          write (1,*)
          write (1,*) 'Processor id = ',irank
          write (1,*)
          write (1,200) 'Comm. time in SEAM    = ',tbuf(7,irank)
          write (1,200) 'Comm. time in CGDOT   = ',tbuf(8,irank)
          write (1,200) 'Comm. time in SUMPHI  = ',tbuf(9,irank)
          write (1,*)   '------------------------------------'
          write (1,200) 'Total comm. time      = ',c_tot
          write (1,*)
          write (1,200) 'Time used in start-up = ',tbuf(2,irank)
          write (1,200) 'Time used in i/o      = ',tbuf(3,irank)
          write (1,200) 'Time used in POTFLD   = ',tbuf(1,irank)
          write (1,200) 'Time used in AX       = ',tbuf(4,irank)
          write (1,200) 'Time used in PCLOAD   = ',tbuf(5,irank)
          write (1,200) 'Time used in PC       = ',tbuf(6,irank)
          write (1,*)   '------------------------------------'
          write (1,200) 'Total time used       = ',tbuf(10,irank)
  100     format (80('-'))
  200     format (1x,a,f12.6)
        enddo
        write (1,*)
        write (1,100)
c
        write (1,*)
        write (1,*) 'Average times:'
        write (1,*) '-------------'
        write (1,*)
        write (1,300) 'Avg         Min         Max      S. Dev'
        write (1,300) '---         ---         ---      ------'
        write (1,400) 'Avg comm. time in SEAM    = ',
     &                tavg(7),tmin(7),tmax(7),tsdev(7)
        write (1,400) 'Avg comm. time in CGDOT   = ',
     &                tavg(8),tmin(8),tmax(8),tsdev(8)
        write (1,400) 'Avg comm. time in SUMPHI  = ',
     &                tavg(9),tmin(9),tmax(9),tsdev(9)
        write (1,400) 'Avg time used in start-up = ',
     &                tavg(2),tmin(2),tmax(2),tsdev(2)
        write (1,400) 'Avg time used in i/o      = ',
     &                tavg(3),tmin(3),tmax(3),tsdev(3)
        write (1,400) 'Avg time used in POTFLD   = ',
     &                tavg(1),tmin(1),tmax(1),tsdev(1)
        write (1,400) 'Avg time used in AX       = ',
     &                tavg(4),tmin(4),tmax(4),tsdev(4)
        write (1,400) 'Avg time used in PCLOAD   = ',
     &                tavg(5),tmin(5),tmax(5),tsdev(5)
        write (1,400) 'Avg time used in PC       = ',
     &                tavg(6),tmin(6),tmax(6),tsdev(6)
        write (1,400) 'Avg total time            = ',
     &                tavg(10),tmin(10),tmax(10),tsdev(10)
  300   format (1x,33x,a)
  400   format (1x,a,4f12.3)
c
        write (1,*)
        write (1,200) 'Average time used per proc  = ',t_tot_avg
        write (1,200) 'Average comm. time per proc = ',c_tot_avg
        write (1,*)
        write (1,100)
        write (1,*)
c
        close (1)
c
      end if
c
      end subroutine
c#######################################################################
      subroutine readbr (fname,br0_g,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read in the radial magnetic field at the photosphere
c ****** and interpolate it into array BR0_G.
c
c ****** FNAME is the name of the file to read.
c
c-----------------------------------------------------------------------
c
      use number_types
      use global_dims
      use global_mesh
      use rdhdf_2d_interface
      use vars
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      real(r_typ), dimension(nt_g,np_g) :: br0_g
      integer :: ierr
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: two=2._r_typ
c
c-----------------------------------------------------------------------
c
c ****** Br array read in and its scales.
c
      real(r_typ), dimension(:,:), pointer :: bn
      real(r_typ), dimension(:), pointer :: tn,pn
c
c-----------------------------------------------------------------------
c
      integer :: ntn,npn,j,k
      logical :: scale
      real(r_typ) :: sum0,sum1,area,fluxp,fluxm,da,br00err
c
c-----------------------------------------------------------------------
c
      ierr=0
c
c ****** Read in the normal field.
c
      write (*,*)
      write (*,*) '### COMMENT from READBR:'
      write (*,*) '### Reading Br file: ',trim(fname)
c
      call rdhdf_2d (fname,scale,ntn,npn,bn,tn,pn,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READBR:'
        write (*,*) '### The flux file has the wrong format.'
        write (*,*) 'IERR (from RDHDF_2D) = ',ierr
        write (*,*) 'File name: ',trim(fname)
        ierr=1
        return
      end if
c
c ****** Check that the arrays has scales.
c
      if (.not.scale) then
        write (*,*)
        write (*,*) '### ERROR in READBR:'
        write (*,*) '### The flux file does not have scales.'
        write (*,*) 'File name: ',trim(fname)
        ierr=2
        return
      end if
c
c ****** Interpolate the field to the code mesh (into array BR0_G).
c
      call intrp2d (ntn,npn,tn,pn,bn,
     &              nt_g-2,np_g-2,th_g(2:ntm1_g),ph_g(2:npm1_g),
     &              br0_g(2:ntm1_g,2:npm1_g),ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READBR:'
        write (*,*) '### The scales in the Br file are invalid.'
        write (*,*) 'File name: ',trim(fname)
        ierr=3
        return
      end if
c
c ****** De-allocate the memory for the BN array and its scales.
c
      deallocate (bn)
      deallocate (tn)
      deallocate (pn)
c
c ****** Set Br to be periodic.
c
      br0_g(:,1)=br0_g(:,npm1_g)
      br0_g(:,np_g)=br0_g(:,2)
c
c ****** Set BCs at the poles.
c
c ****** Br has only an m=0 component there.
c
      sum0=sum(br0_g(     2,2:npm1_g)*dph_g(2:npm1_g))*pl_i
      sum1=sum(br0_g(ntm1_g,2:npm1_g)*dph_g(2:npm1_g))*pl_i
c
      br0_g(1   ,:)=two*sum0-br0_g(     2,:)
      br0_g(nt_g,:)=two*sum1-br0_g(ntm1_g,:)
c
c ****** Calculate the total flux.
c
      area=0.
      fluxp=0.
      fluxm=0.
      do j=2,ntm1_g
        do k=2,npm1_g
          da=sth_g(j)*dth_g(j)*dph_g(k)
          if (br0_g(j,k).gt.0.) then
            fluxp=fluxp+br0_g(j,k)*da
          else
            fluxm=fluxm+br0_g(j,k)*da
          end if
          area=area+da
        enddo
      enddo
c
      write (*,*)
      write (*,*) '### COMMENT from READBR:'
      write (*,*) '### Computed flux balance:'
      write (*,*)
      write (*,*) 'Positive flux = ',fluxp
      write (*,*) 'Negative flux = ',fluxm
c
c ****** Fix the magnetic field so that the total flux is zero
c ****** (unless this has not been requested).
c
      if (.not.((option.eq.'ss'.or.option.eq.'open')
     &          .and.do_not_balance_flux)) then
c
        br00err=(fluxp+fluxm)/area
c
        do k=1,np_g
          do j=1,nt_g
            br0_g(j,k)=br0_g(j,k)-br00err
          enddo
        enddo
c
        write (*,*)
        write (*,*) '### COMMENT from READBR:'
        write (*,*) '### Flux balance correction:'
        write (*,*)
        write (*,*) 'BR00 (monopole Br field magnitude) = ',br00err
c
      end if
c
      return
      end
c#######################################################################
      subroutine intrp2d (nxi,nyi,xi,yi,fi,nx,ny,x,y,f,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Interpolate a 2D field from array FI(NXI,NYI), defined
c ****** on the mesh XI(NXI) x YI(NYI), into the array F(NX,NY),
c ****** defined on the mesh X(NX) x Y(NY).
c
c ****** Note that if a data point is outside the bounds of
c ****** the XI x YI mesh, IERR=2 is returned.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: nxi,nyi
      real(r_typ), dimension(nxi) :: xi
      real(r_typ), dimension(nyi) :: yi
      real(r_typ), dimension(nxi,nyi) :: fi
      integer :: nx,ny
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      real(r_typ), dimension(nx,ny) :: f
      integer :: ierr
c
c-----------------------------------------------------------------------
c
      real(r_typ), parameter :: zero=0._r_typ
c
c-----------------------------------------------------------------------
c
      real(r_typ), external :: flint
c
c-----------------------------------------------------------------------
c
      integer :: i,j,ii,jj
      real(r_typ) :: dum(1)=0.
      real(r_typ) :: dummy,xv,yv,ax,ay
c
c-----------------------------------------------------------------------
c
      ierr=0
c
c ****** Check that the scales XI and YI are monotonic.
c
      dummy=flint(zero,nxi,xi,dum,1,ierr)
      if (ierr.ne.0) go to 900
c
      dummy=flint(zero,nyi,yi,dum,1,ierr)
      if (ierr.ne.0) go to 900
c
c ****** Interpolate the data.
c
      do j=1,ny
        yv=y(j)
        if (yv.lt.yi(1).or.yv.gt.yi(nyi)) then
          go to 910
        end if
        call interp (yi,nyi,yv,jj,ay)
        do i=1,nx
          xv=x(i)
          if (xv.lt.xi(1).or.xv.gt.xi(nxi)) then
            go to 910
          end if
          call interp (xi,nxi,xv,ii,ax)
          f(i,j)= (1.-ax)*((1.-ay)*fi(ii  ,jj  )+ay*fi(ii  ,jj+1))
     &           +    ax *((1.-ay)*fi(ii+1,jj  )+ay*fi(ii+1,jj+1))
        enddo
      enddo
c
      return
c
c ****** Error exit: invalid scales.
c
  900 continue
c
      write (*,*)
      write (*,*) '### ERROR in INTRP2D:'
      write (*,*) '### Scales are not monotonically increasing.'
      ierr=1
c
      return
c
c ****** Error exit: data outside range of scales.
c
  910 continue
c
      write (*,*)
      write (*,*) '### ERROR in INTRP2D:'
      write (*,*) '### An interpolation was attempted outside the'//
     &            ' range of the defined scales.'
      ierr=2
c
      return
      end
c#######################################################################
      function flint (x,n,xn,fn,icheck,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Interpolate a function linearly.
c
c-----------------------------------------------------------------------
c
c ****** The funcion is defined at N nodes, with values given by
c ****** FN(N) at positions XN(N).  The function value returned is
c ****** the linear interpolant at X.
c
c ****** Note that if X.lt.XN(1), the function value returned
c ****** is FN(1), and if X.gt.XN(N), the function value returned
c ****** is FN(N).
c
c ****** Call once with ICHECK.ne.0 to check that the values
c ****** in XN(N) are monotonically increasing.  In this mode
c ****** the array XN(N) is checked, and X and FN(N) are not
c ****** accessed.  If the check is passed, IERR=0 is returned.
c ****** Otherwise, IERR=1 is returned.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      real(r_typ) :: flint
      real(r_typ) :: x
      integer :: n
      real(r_typ), dimension(n) :: xn,fn
      integer :: icheck
      integer :: ierr
c
c-----------------------------------------------------------------------
c
      integer :: i,j
      real(r_typ) :: x1,x2,alpha
c
c-----------------------------------------------------------------------
c
      ierr=0
      flint=0.
c
c ****** If ICHECK.ne.0, check the abscissa table.
c
      if (icheck.ne.0) then
        if (n.le.0) then
          write (*,*)
          write (*,*) '### ERROR in FLINT:'
          write (*,*) '### Bad dimension of abscissa table.'
          write (*,*) 'N = ',n
          ierr=1
          return
        end if
        do 100 i=1,n-1
          if (xn(i+1).le.xn(i)) then
            write (*,*)
            write (*,*) '### ERROR in FLINT:'
            write (*,*) '### Bad data in abscissa table.'
            write (*,*) 'N = ',n
            write (*,*) 'XN = '
            write (*,*) (xn(j),j=1,n)
            ierr=1
            return
          end if
  100   continue
        return
      end if
c
c ****** Get the interpolated value.
c
      if (x.le.xn(1)) then
        flint=fn(1)
      else if (x.gt.xn(n)) then
        flint=fn(n)
      else
        do 200 i=1,n-1
          if (x.ge.xn(i).and.x.lt.xn(i+1)) go to 300
  200   continue
  300   continue
        x1=xn(i)
        x2=xn(i+1)
        alpha=(x-x1)/(x2-x1)
        flint=fn(i)*(1.-alpha)+fn(i+1)*alpha
      end if
c
      return
      end
c#######################################################################
      subroutine interp (x,n,xv,i,alpha)
c
c-----------------------------------------------------------------------
c
c ****** Get interpolation factor ALPHA and table index I.
c
c ****** This routine does not do the actual interpolation.  Use the
c ****** returned values of I and ALPHA to get the interpolated value.
c
c-----------------------------------------------------------------------
c
      use number_types
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: n
      real(r_typ), dimension(n) :: x
      real(r_typ) :: xv
      integer :: i
      real(r_typ) :: alpha
c
c-----------------------------------------------------------------------
c
      do 100 i=1,n-1
        if (xv.ge.x(i).and.xv.le.x(i+1)) then
          alpha=(xv-x(i))/(x(i+1)-x(i))
          go to 200
        end if
  100 continue
c
c ****** Value not found --- signal error and stop.
c
      write (*,*)
      write (*,*) '### ERROR in INTERP:'
      write (*,*) '### Value not found in table.'
      write (*,*) 'Value requested = ',xv
      write (*,*) 'Min table value = ',x(1)
      write (*,*) 'Max table value = ',x(n)
      call endrun (.true.)
c
  200 continue
c
      return
      end
c#######################################################################
c
c ****** Revision history:
c
c ### Version 1.00, 03/02/2006, file pot3d.f, modified by ZM:
c
c       - Cleaned up the previous version of POT3D.
c
c ### Version 1.01, 03/06/2006, file pot3d.f, modified by ZM:
c
c       - Added the ability to do a "source-surface plus
c         current-sheet" solution.  Select this by setting
c         OPTION='ss+cs'.
c
c ### Version 1.02, 06/18/2007, file pot3d.f, modified by ZM:
c
c       - Fixed a bug that caused the code to hang when an error
c         was encountered (when running a parallel job).
c
c ### Version 1.03, 03/17/2009, file pot3d.f, modified by ZM:
c
c       - Added the ability to write the boundary flux before the
c         sign flip for current-sheet solutions (i.e., OPTION='open').
c         Set the variable BR_PHOTO_ORIGINAL_FILE to the desired
c         file name to request this.
c
c ### Version 1.50, 01/25/2016, file pot3d.f, modified by RC:
c
c       - Added new (much faster) BILU0 preconditioner to CG solver.
c         To activate, set ifprec=2 in pot3d.dat file.
c       - Modified CG solve to use 1D arrays
c         for SAXPY and DOT operations.
c
c ### Version 2.00, 06/06/2017, file pot3d.f, modified by RC:
c
c       - Added OpenACC support to the code.
c         - OpenACC support is currently ONLY on 'standard'
c           pot3d runs (not inner-outer-iteratative mode)
c           and is only efficient on GPUs when using ifprec=1.
c         - OpenACC adds support for running the code on
c           Nvidia GPUs/Intel KNL/x86-multicore/OpenPower.
c         - To use OpenACC, simply compile the code with a compiler
c           that supports OpenACC with the correct flags activated.
c         - Multi-gpu support is included by using the new
c           ngpus_per_node input parameter.  This should be set
c           to the number of GPUs available per node.
c           The number of MPI ranks per node should match the
c           number of gpus per node.  This can be launched with
c           "mpiexec -np <np> -ntasks-per-node <ngpus_per_node>".
c         - The GPU features of the code are fully portable, i.e.
c           the code can be compiled/used as before on CPUs with no
c           changes in compilation or run-time.
c       - Modified some routines to be "nicer" for OpenACC
c         and optimized some MPI communications.
c       - Added wall-clock timer and corrected placement of
c         MPI_Finalize().  The wall clock timer now reflects the
c         true runtime.
c
c ### Version 2.01, 10/02/2017, file pot3d.f, modified by RC:
c
c       - Optimized OpenACC.
c       - Renamed cgsolv() to cgsolve().
c
c ### Version 2.10, 01/15/2018, file pot3d.f, modified by ZM+RC:
c
c       - Added the ability to skip the balancing of the flux
c         when doing a PFSS or OPEN field.  To invoke this, set
c         DO_NOT_BALNCE_FLUX=.true..
c       - Changed some pointers to allocatables for better
c         vectorization.
c
c ### Version 2.11, 03/19/2018, file pot3d.f, modified by RC:
c
c       - Added 'fmt' input parameter to set output file type.
c         Set fmt='h5' to write out HDF5 (default is 'hdf').
c
c ### Version 2.12, 10/08/2018, file pot3d.f, modified by RC:
c
c       - COMPATIBILITY CHANGE! Renamed gpus_per_node to gpn.
c         gpn is default 0 which will set gpn to the number of
c         MPI ranks in the local node.
c         Setting gpn manually is not recommended and only used for
c         oversubscribing the GPUs.
c       - Added MPI shared communicator for automatically setting gpn.
c         This requires an MPI-3 capable MPI library.
c       - Changed layout of matrix coefficient arrays to be more
c         vector-friendly instead of cache-friendly.
c
c ### Version 2.13, 11/19/2018, file pot3d.f, modified by RC:
c
c       - Small modifications to polar boundary condition calculations.
c       - Updated array ops and ACC clauses to be F2003 optimized.
c
c ### Version 2.20, 01/09/2019, file pot3d.f, modified by RC:
c
c       - Added double precision output option.
c         Set hdf32=.false. to activate 64-bit output.
c       - Updated magnetic field computation.  B is now computed
c         in parallel.  3D output fields now gathered to rank 0
c         using Sends and Receives instead of Gatherv in order
c         to allow very large resolutions.
c       - Added automatic topology.  Now, nprocs is optional.
c         One can specify one or more topology dimensions and
c         use the flag value "-1" to indicate dimensions to auto-set.
c         It is recommended to simply delete nprocs from input files.
c       - Added output file flushing so CG iterations can be monitored.
c       - Added new MPI rank diagnostics including
c         estimated load imbalance.
c       - Processor topology and magnetic energy output now written to
c         pot3d.out as well as the terminal.
c
c ### Version 2.21, 01/11/2019, file pot3d.f, modified by RC:
c
c       - Small updates to magnetic_energy routine.
c
c ### Version 2.22, 11/27/2019, file pot3d.f, modified by RC:
c
c       - Optimized some OpenACC directives.  Expanded some
c         array-syntax lines to full loops.
c
c ### Version 2.23, 08/11/2020, file pot3d.f, modified by RC:
c
c       - Small bug fix for default output file names and
c         format option fmt.
c
c ### Version 3.0.0, 02/10/2021, file pot3d.f, modified by RC:
c
c       - Changed version number scheme to semantic versioning.
c
c ### Version 3.1.0, 08/11/2020, file pot3d.f, modified by RC:
c
c       - Removed ss+cs inner-outer iterative mode.
c         It had not been used in many years and never worked
c         quite right.  Algorithm is saved in previous releases
c         on the github and svn.
c       - Added a validation run mode.  To activate, set:
c         VALIDATION_RUN=.TRUE.
c         This will set up a tilted dipole exact solution, write it
c         out, and then solve for it.  It will ignore input br files.
c         To change the angle of the dipole, set DIPOLE_ANGLE
c         (default is pi/4).
c       - Added NCGFLUSH parameter to flush output file during solver
c         iterations if NCGHIST>0 every NCGFLUSH iteration (default 25).
c         This is useful to track progress of the solver.
c
c ### Version 3.1.1, 11/02/2021, modified by RC:
c
c       - Fixed default output name mechanics.  The default names
c         were ignoring "fmt", making "fmt" almost useless.
c
c ### Version 3.2.0, 12/21/2021, modified by RC:
c
c       - GPN is now depricated.  Now, # of GPUs per node is assumed
c         to be equal to # of MPI ranks per node.
c         This removes the ability to oversubscribe GPUs, but this
c         was not being used ever, so no great loss.
c         This change was done to allow setting the gpu device
c         number before any allocations so that it works correctly
c         when using NVIDIA unified managed memory.
c
c ### Version 3.3.0, 02/14/2022, modified by RC:
c
c       - Added new timers for AX, PC, PCLOAD, and modified some timers
c         to be more useful.
c       - Updated ILU0 PC2 preconditioner to use single precision.
c         This speeds up the solve both on CPUs and GPUs.
c
c ### Version 4.0.0, 03/02/2022, modified by RC/MS:
c
c       - Added cuSparse library calls for ifprec=2 on NVIDIA GPUs.
c       - Added #ifdef
c         preprocessor directives which will either compile a cusparse
c         GPU version of the code or a CPU compatible version
c         without cusparse calls and no c code references
c
c#######################################################################
