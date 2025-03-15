!#######################################################################
!  _____   ____ _______ ____  _____
! |  __ \ / __ \__   __|___ \|  __ \
! | |__) | |  | | | |    __) | |  | |
! |  ___/| |  | | | |   |__ <| |  | |
! | |    | |__| | | |   ___) | |__| |
! |_|     \____/  |_|  |____/|_____/
!
! ****** POT3D: 3D potential magnetic field outside a sphere.
!
! ****** This program can find the classical potential field, the
! ****** fully open field, the source-surface field, and the
! ****** source-surface plus current-sheet field.
!
!        Authors:  Zoran Mikic
!                  Ronald M. Caplan
!                  Jon A. Linker
!                  Roberto Lionello
!                  Miko Stulajter
!
!        Predictive Science Inc.
!        www.predsci.com
!        San Diego, California, USA 92121
!
!#######################################################################
! Copyright 2021 Predictive Science Inc.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!    http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
! implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!#######################################################################
!
!#######################################################################
!
module ident
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
! ****** Code name.
!-----------------------------------------------------------------------
!
      character(*), parameter :: idcode='POT3D'
      character(*), parameter :: vers  ='4.3.1'
      character(*), parameter :: update='12/05/2024'
!
end module
!#######################################################################
module number_types
!
!-----------------------------------------------------------------------
! ****** Basic number types.
! ****** This module is used to set the default precision for REALs.
!-----------------------------------------------------------------------
!
      use iso_fortran_env
!
!-----------------------------------------------------------------------
!
      implicit none
!
      integer, parameter :: KIND_REAL_4=REAL32
      integer, parameter :: KIND_REAL_8=REAL64
      integer, parameter :: KIND_REAL_16=max(REAL128,REAL64)
!
      integer, parameter :: r_typ=KIND_REAL_8
!
end module
!#######################################################################
module number_types_pc
!
!-----------------------------------------------------------------------
!
      use number_types
      use iso_fortran_env
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer, parameter :: r_typ_pc=REAL32
!
end module
!#######################################################################
module constants
!
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
      real(r_typ), parameter :: pi=3.1415926535897932_r_typ
!
end module
!#######################################################################
module global_dims
!
!-----------------------------------------------------------------------
! ****** Global number of mesh points.
!-----------------------------------------------------------------------
!
      implicit none
!
! ****** Global mesh size.
!
      integer :: nr_g,nrm1_g
      integer :: nt_g,ntm1_g
      integer :: np_g,npm1_g
!
! ****** Flag to indicate an axisymmetric run.
!
      logical :: axisymmetric=.false.
!
end module
!#######################################################################
module local_dims
!
!-----------------------------------------------------------------------
! ****** Local number of mesh points.
!-----------------------------------------------------------------------
!
      implicit none
!
! ****** Flags to indicate whether this processor has points
! ****** on the physical boundaries.
!
      logical :: rb0,rb1,tb0,tb1
!
! ****** Local mesh size.
!
      integer :: nr,nrm1
      integer :: nt,ntm1
      integer :: np,npm1
!
! ****** Dimensions of arrays on the "main" mesh.
!
      integer :: nrm
      integer :: ntm
      integer :: npm
!
! ****** Indices of start and end points in the global mesh
! ****** belonging to this processor.
!
      integer :: i0_g,i1_g
      integer :: j0_g,j1_g
      integer :: k0_g,k1_g
!
end module
!#######################################################################
module global_mesh
!
!-----------------------------------------------------------------------
! ****** Global mesh.
!-----------------------------------------------------------------------
!
      use number_types
      use constants
!
      implicit none
!
      real(r_typ), dimension(:), allocatable :: r_g,rh_g,dr_g,drh_g
      real(r_typ), dimension(:), allocatable :: t_g,th_g,dt_g,dth_g
      real(r_typ), dimension(:), allocatable :: p_g,ph_g,dp_g,dph_g
!
      real(r_typ), dimension(:), allocatable :: st_g,ct_g,sth_g,cth_g
      real(r_typ), dimension(:), allocatable :: sp_g,cp_g,sph_g,cph_g
!
! ****** Physical mesh size.
!
      real(r_typ), parameter, private :: one=1._r_typ
      real(r_typ), parameter, private :: two=2._r_typ
!
      real(r_typ) :: r0=1._r_typ
      real(r_typ) :: r1=30._r_typ
      real(r_typ), parameter :: t0=0.
      real(r_typ), parameter :: t1=pi
      real(r_typ), parameter :: p0=0.
      real(r_typ), parameter :: p1=two*pi
!
      real(r_typ), parameter :: pl=p1-p0
      real(r_typ), parameter :: pl_i=one/pl
!
end module
!#######################################################################
module local_mesh
!
!-----------------------------------------------------------------------
! ****** Local mesh.
!-----------------------------------------------------------------------
!
      use number_types
!
      implicit none
!
      real(r_typ), dimension(:), allocatable :: r,r2,rh,dr,drh
      real(r_typ) :: dr1,drn
!
      real(r_typ), dimension(:), allocatable :: t,th,dt,dth
      real(r_typ), dimension(:), allocatable :: p,ph,dp,dph
!
      real(r_typ), dimension(:), allocatable :: st,ct,sth,cth
      real(r_typ), dimension(:), allocatable :: sp,cp,sph,cph
!
! ****** Inverse quantities (for efficiency).
!
      real(r_typ), dimension(:), allocatable :: r_i,rh_i
      real(r_typ), dimension(:), allocatable :: dr_i,drh_i
      real(r_typ), dimension(:), allocatable :: dt_i,dth_i
      real(r_typ), dimension(:), allocatable :: st_i,sth_i
      real(r_typ), dimension(:), allocatable :: dp_i,dph_i
!
end module
!#######################################################################
module mpidefs
!
!-----------------------------------------------------------------------
! ****** MPI variables, processor topology, and processor information.
!-----------------------------------------------------------------------
!
      use mpi
!
      implicit none
!
! ****** Total number of processors.
!
      integer :: nproc
!
! ****** Total number of processors per node.
!
      integer :: nprocsh
!
! ****** Processor rank of this process in communicator
! ****** MPI_COMM_WORLD.
!
      integer :: iprocw
!
! ****** Processor rank of this process in communicator
! ****** comm_shared.
!
      integer :: iprocsh
!
! ****** Flag to designate that this is the processor with
! ****** rank 0 in communicator MPI_COMM_WORLD.
!
      logical :: iamp0
!
! ****** Communicator over all processors in the Cartesian topology.
!
      integer :: comm_all
!
! ****** Processor rank of this process in communicator
! ****** COMM_ALL.
!
      integer :: iproc
!
! ****** Processor rank in communicator COMM_ALL for the
! ****** processor that has rank 0 in MPI_COMM_WORLD.
!
      integer :: iproc0
!
! ****** Communicators over all processors in the phi dimension.
!
      integer :: comm_phi
!
! ****** Communicator over all shared processors on the node.
!
      integer :: comm_shared
!
! ****** Communicators over all processors in the r dimension.
!
      integer :: comm_r
!
! ****** Processor coordinate indices of this process
! ****** in the Cartesian topology.
!
      integer :: iproc_r,iproc_t,iproc_p
!
! ****** Processor coordinate indices of the neighboring
! ****** processors in the Cartesian topology.
!
      integer :: iproc_rm,iproc_rp
      integer :: iproc_tm,iproc_tp
      integer :: iproc_pm,iproc_pp
!
! ****** Number of processors along r, theta, and phi.
!
      integer :: nproc_r,nproc_t,nproc_p
!
! ****** Number type for REALs to be used in MPI calls.
!
      integer :: ntype_real
!
end module
!#######################################################################
module decomposition_params
!
!-----------------------------------------------------------------------
! ****** Input parameters that define the domain decomposition
! ****** among processors.
!-----------------------------------------------------------------------
!
      implicit none
!
! ****** Number of processors per dimension.
!
      integer, dimension(3) :: nprocs=(/-1,-1,-1/)
!
end module
!#######################################################################
module decomposition
!
!-----------------------------------------------------------------------
! ****** Information defining the domain decomposition.
!-----------------------------------------------------------------------
!
      implicit none
!
! ****** Define the structure type for mapping local arrays
! ****** into global arrays.
!
      type :: map_struct
        integer :: n
        integer :: i0
        integer :: i1
        integer :: offset
      end type
!
! ****** Mapping structures for the different mesh types.
!
      type(map_struct), dimension(:), pointer :: map_rh
      type(map_struct), dimension(:), pointer :: map_rm
      type(map_struct), dimension(:), pointer :: map_th
      type(map_struct), dimension(:), pointer :: map_tm
      type(map_struct), dimension(:), pointer :: map_ph
      type(map_struct), dimension(:), pointer :: map_pm
!
end module
!#######################################################################
module meshdef
!
!-----------------------------------------------------------------------
! ****** Variables that define the mesh-point distributions.
!-----------------------------------------------------------------------
!
      use number_types
!
      implicit none
!
      integer, parameter :: nmseg=30
!
      real(r_typ), dimension(nmseg) :: drratio=0.
      real(r_typ), dimension(nmseg) :: dtratio=0.
      real(r_typ), dimension(nmseg) :: dpratio=0.
      real(r_typ), dimension(nmseg) :: rfrac=0.
      real(r_typ), dimension(nmseg) :: tfrac=0.
      real(r_typ), dimension(nmseg) :: pfrac=0.
!
      integer :: nfrmesh=0
      integer :: nftmesh=0
      integer :: nfpmesh=0
!
      real(r_typ) :: phishift=0.
!
end module
!#######################################################################
module fields
!
!-----------------------------------------------------------------------
! ****** Local field arrays.
!-----------------------------------------------------------------------
!
      use number_types
!
      implicit none
!
! ****** Potential solution array and cg temp array.
!
      real(r_typ), dimension(:,:,:), allocatable :: phi
      real(r_typ), dimension(:,:,:), allocatable :: x_ax
!
! ****** Boundary radial magnetic field arrays.
!
      real(r_typ), dimension(:,:), allocatable :: br0
      real(r_typ), dimension(:,:), allocatable :: br1
!
! ****** Arrays used in polar boundary conditions.
!
      real(r_typ), dimension(:), allocatable :: sum0,sum1
!
! ****** Arrays used for final magnetic field.
!
      real(r_typ), dimension(:,:,:), allocatable :: br,bt,bp
!
end module
!#######################################################################
module cgcom
!
      use number_types
!
      implicit none
!
!-----------------------------------------------------------------------
! ****** Number of equations to solve in the CG solver.
!-----------------------------------------------------------------------
!
      integer :: ncgeq
!
!-----------------------------------------------------------------------
! ****** CG field solver parameters.
!-----------------------------------------------------------------------
!
      integer :: ifprec=1
      integer :: ncgmax=1000000
      integer :: ncghist=100
      real(r_typ) :: epscg=1.e-9
!
!-----------------------------------------------------------------------
! ****** CG field solver variables.
!-----------------------------------------------------------------------
!
      integer :: ncg
      real(r_typ) :: epsn
!
! ****** Seam buffers.
!
      real(r_typ), allocatable, dimension(:,:) :: sbuf_rt1,sbuf_rt2
      real(r_typ), allocatable, dimension(:,:) :: sbuf_tp1,sbuf_tp2
      real(r_typ), allocatable, dimension(:,:) :: sbuf_rp1,sbuf_rp2
!
      real(r_typ), allocatable, dimension(:,:) :: rbuf_rt1,rbuf_rt2
      real(r_typ), allocatable, dimension(:,:) :: rbuf_tp1,rbuf_tp2
      real(r_typ), allocatable, dimension(:,:) :: rbuf_rp1,rbuf_rp2
!
end module
!#######################################################################
module vars
!
!-----------------------------------------------------------------------
! ****** Miscellaneous input variables.
!-----------------------------------------------------------------------
!
      use number_types
!
      implicit none
!
      character(256) :: outfile='pot3d.out'
      character(256) :: phifile='default'
      character(256) :: br0file='default'
      character(256) :: brfile='default'
      character(256) :: btfile='default'
      character(256) :: bpfile='default'
      character(256) :: br_photo_file='default'
      character(256) :: br_photo_original_file='default'
!
! ****** Type of field solution.
! ****** Select between 'potential', 'open', and 'source-surface'.
!
      character(16) :: option='potential'
!
! ****** Interval at which to dump diagonstics during the
! ****** iteration for the source-surface plus current-sheet
! ****** solution.
!
      integer :: ndump=0
!
! ****** Flag to skip the balancing of the flux (for PFSS and
! ****** OPEN field options only).

      logical :: do_not_balance_flux=.false.
!
! ****** Set format for output files.
!
      character(3) :: fmt='h5'
!
      logical :: hdf32=.true.
!
! ***** Validation run (tilted dipole).
!
      logical :: validation_run=.false.
!
      real(r_typ) :: dipole_angle=0.7853981633974483_r_typ
!
end module
!#######################################################################
module solve_params
!
!-----------------------------------------------------------------------
! ****** Parameters used in the solver.
!-----------------------------------------------------------------------
!
      use number_types
!
      implicit none
!
! ****** Boundary condition switch at r=R1.
!
      real(r_typ) :: pm_r1
!
end module
!#######################################################################
module timer
!
!-----------------------------------------------------------------------
! ****** Timer stack.
!-----------------------------------------------------------------------
!
      use number_types
!
      implicit none
!
      integer, parameter :: nstack=10
      integer :: istack=0
      real(r_typ), dimension(nstack) :: tstart=0.
!
end module
!#######################################################################
module timing
!
!-----------------------------------------------------------------------
! ****** Timing variables.
!-----------------------------------------------------------------------
!
      use number_types
!
      implicit none
!
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
!
end module
!#######################################################################
module debug
!
!-----------------------------------------------------------------------
! ****** Debugging level.
!-----------------------------------------------------------------------
!
      implicit none
!
      integer :: idebug=0
!
end module
!#######################################################################
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
!#######################################################################
module cusparse_interface
!
      use, intrinsic :: iso_c_binding
!
#ifdef CUSPARSE
      interface
  subroutine load_lusol_cusparse(CSR_A,CSR_AI,CSR_AJ,N,M) &
          BIND(C, name="load_lusol_cusparse")
          use, intrinsic :: iso_c_binding
          implicit none
          integer(C_INT), value :: N,M
          type(C_PTR), value :: CSR_A,CSR_AI,CSR_AJ
  end subroutine load_lusol_cusparse
!
  subroutine lusol_cusparse(x) &
          BIND(C, name="lusol_cusparse")
          use, intrinsic :: iso_c_binding
          implicit none
          type(C_PTR), value :: x
  end subroutine lusol_cusparse
!
  subroutine unload_lusol_cusparse() &
          BIND(C, name="unload_lusol_cusparse")
  end subroutine unload_lusol_cusparse
      end interface
#endif
!
      integer(c_int) :: cN,cM
!
end module
!#######################################################################
module matrix_storage_pot3d_solve
!
!-----------------------------------------------------------------------
! ****** Storage for the matrix/preconditioners of the solve.
!-----------------------------------------------------------------------
!
      use number_types
      use number_types_pc
!
      implicit none
!
      real(r_typ), dimension(:,:,:,:), allocatable :: a
      real(r_typ_pc), dimension(:), allocatable :: a_i
!
      integer, dimension(7) :: a_offsets

      integer :: N,M
      real(r_typ_pc), dimension(:), allocatable, target :: a_csr
      real(r_typ_pc), dimension(:), allocatable :: lu_csr
      real(r_typ_pc), dimension(:), allocatable :: a_csr_d
      integer, dimension(:), allocatable :: lu_csr_ja
      integer, dimension(:), allocatable, target :: a_csr_ia
      integer, dimension(:), allocatable, target :: a_csr_ja
      integer, dimension(:), allocatable :: a_N1
      integer, dimension(:), allocatable :: a_N2
      integer, dimension(:), allocatable :: a_csr_dptr
!
end module
!#######################################################################
program POT3D
!
!-----------------------------------------------------------------------
!
      use ident
      use mpidefs
      use vars
      use solve_params
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: ierr
!
!-----------------------------------------------------------------------
!
! ****** Initialize MPI.
!
      call init_mpi
!
! ****** Start the wall-clock timer.
!
      call timer_on
!
! ****** Write the code name and version.
!
      if (iamp0) then
        write (*,*)
        write (*,*) idcode,' Version ',vers,', updated on ',update
      end if
!
      call timer_on
!
! ****** Read the input file.
!
      call read_input_file
!
! ****** Create the output file.
!
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
!
! ****** Check the input parameters.
!
      call check_input
!
! ****** Check the processor topology.
!
      call check_proc_topology
!
! ****** Decompose the domain.
!
      call decompose_domain
!
! ****** Allocate global arrays.
!
      call allocate_global_arrays
!
! ****** Set the global meshes.
!
      call set_global_mesh
!
! ****** Decompose the mesh.
!
      call decompose_mesh_r
      call decompose_mesh_tp
!
! ****** Allocate local arrays.
!
      call allocate_local_arrays_r
      call allocate_local_arrays_tp
!
! ****** Set the local meshes.
!
      call set_local_mesh_r
      call set_local_mesh_tp
!
! ****** Print decomposition diagnostics.
!
      call decomp_diags
!
! ****** Set up seam for solver.
!
      call seam_setup
!
! ****** Initialize the flux and if validating, write analytic solution.
!
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
!
      call timer_off (t_startup)
!
! ****** Find the solution.
!
      if (iamp0) then
        write (*,*)
        write (*,*) '### COMMENT from POT3D:'
        write (*,*) '### Starting PCG solve.'
        call FLUSH(OUTPUT_UNIT)
        write (9,*)
        write (9,*) '### COMMENT from POT3D:'
        write (9,*) '### Starting PCG solve.'
        call FLUSH(9)
      end if
!
      call timer_on
!
      call potfld
!
      call timer_off (t_solve)
!
! ****** Compute B.
!
      call getb
!
! ****** Write solution to file.
!
      call write_solution
!
! ****** Magnetic energy diagnostics.
!
      call magnetic_energy
!
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      call timer_off (t_wall)
!
      call write_timing
!
      call endrun (.false.)
!
end program pot3d
!#######################################################################
subroutine read_input_file
!
!-----------------------------------------------------------------------
!
! ****** Read the input file.
!
!-----------------------------------------------------------------------
!
      use global_dims
      use global_mesh
      use mpidefs
      use meshdef
      use cgcom
      use debug
      use vars
      use decomposition_params
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
! ****** Values for the global mesh size.
! ****** Since these names conflict with those in LOCAL_DIMS*, it is
! ****** important not to use these modules here.
!
      integer :: nr=0
      integer :: nt=0
      integer :: np=0
!
!-----------------------------------------------------------------------
!
      namelist /topology/ &
        nprocs,           &! MPI topology triplet.  Default -1,-1,-1
                           !   automatically sets "best" topology.
        nr,               &! Grid resolution in the `r` direction.
        nt,               &! Grid resolution in the `theta` direction.
        np                 ! Grid resolution in the `phi` direction.
!
!-----------------------------------------------------------------------
!
      namelist /inputvars/       &
        r0,                      &! Lower radial boundary.
        r1,                      &! Upper radial boundary.
        drratio,                 &! Ratio of grid spacing at the end
                                  ! of each segment to that at the
                                  ! beginning for the radial grid
                                  ! [ length(drratio) = length(rfrac)-1 ].
        dtratio,                 &! Ratio of grid spacing (theta)
        dpratio,                 &! Ratio of grid spacing (phi)
        rfrac,                   &! Normalized positions of grid segment
                                  ! boundaries (frac of domain size)
                                  ! for radial grid.
        tfrac,                   &! Normalized positions of grid (theta)
        pfrac,                   &! Normalized positions of grid (phi)
        nfrmesh,                 &! Number of times to filter/smooth
                                  ! the radial grid spacing.
        nftmesh,                 &! Number of times to filter (theta)
        nfpmesh,                 &! Number of times to filter (phi)
        phishift,                &! Apply an optional phi shift
                                  ! (radians) to the input Br at r0.
        ifprec,                  &! Preconditioner method:
                                  !  1: Diagonal (use for GPU runs)
                                  !  2: ILU0     (use for CPU runs or
                                  !   GPU runs when built with cusparse)
        ncgmax,                  &! Maximum alowed solver iterations.
        ncghist,                 &! Iteration information.
                                  !  0: Only write # total iterations.
                                  ! >0: Write every NCGHIST iteration.
        epscg,                   &! Solver convergence tolerance,
                                  !|residual|/|right-hand-side|.
        idebug,                  &! Output debugging info during run.
        br0file,                 &! Filename of input 2D (tp) Br
        phifile,                 &! Filename to write 3D PHI potential.
        brfile,                  &! Filename to write 3D Br field.
        btfile,                  &! Filename to write 3D Bt field.
        bpfile,                  &! Filename to write 3D Bp field.
        br_photo_file,           &! Filename to write 2D Br@r=r0
                                  !(after interp/flux-balance).
        br_photo_original_file,  &! Filename to write 2D Br@r=r0
                                  ! (option='open' only, writes Br
                                  !  before sign change).
        option,                  &!  'ss'        PFSS
                                  !  'potential' PF with closed-wall
                                  !              (requires flux balance).
                                  !  'open'      Open field. Used for
                                  !              current sheet and fully
                                  !              open field runs.
                                  !              B field will be unsigned.
        do_not_balance_flux,     &! Do not balance flux of input Br@r0.
        hdf32,                   &! Output precision:
                                  !  .true.   Single (32-bit) output.
                                  !  .false.  Double (64-bit) output.
        validation_run,          &!  Set this to run a validation test
                                  ! with an analytic tilted dipole
                                  ! solution. Overrides other inputs.
        dipole_angle              ! Tilt angle for validation run dipole.
!
!-----------------------------------------------------------------------
!
      integer :: ierr
      character(80) :: infile='pot3d.dat'
!
!-----------------------------------------------------------------------
!
! ****** Read the input file.
!
      call ffopen (8,infile,'r',ierr)
!
      if (ierr.ne.0) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in READ_INPUT_FILE:'
          write (*,*) '### Could not open the input file.'
          write (*,*) 'File name: ',trim(infile)
        end if
        call endrun (.true.)
      end if
!
      read (8,topology)
!
      read (8,inputvars)
!
      close (8)
!
      if (trim(fmt).ne.'h5') then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in READ_INPUT_FILE:'
          write (*,*) '### "fmt" must be "h5".'
          write (*,*) 'fmt: ',trim(fmt)
        end if
        call endrun (.true.)
      end if
!
! ****** Check if output names were overwritten.
! ****** If not, set default names with format fmt.
!
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
!
      nr_g=nr
      nt_g=nt
      np_g=np
!
! ****** Check if the specified mesh dimensions are valid.
!
      call check_mesh_dimensions (nr_g,nt_g,np_g)
!
      nrm1_g=nr_g-1
      ntm1_g=nt_g-1
      npm1_g=np_g-1
!
end subroutine
!#######################################################################
subroutine check_error_on_p0 (ierr0)
!
!-----------------------------------------------------------------------
!
! ****** Check if the error flag IERR0 on processor 0 in
! ****** MPI_COMM_WORLD (i.e., processor IPROC0 in COMM_ALL)
! ****** indicates that the code should exit.
!
! ****** If IERR0 is non-zero, all the processors are directed
! ****** to call ENDRUN to terminate the code.
!
!-----------------------------------------------------------------------
!
      use mpidefs
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: ierr0
!
!-----------------------------------------------------------------------
!
! ****** MPI error return.
!
      integer :: ierr
!
!-----------------------------------------------------------------------
!
! ****** Broadcast IERR0 to all processors.
!
      call MPI_Bcast (ierr0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!
! ****** Call ENDRUN if IERR0 is non-zero.
!
      if (ierr0.ne.0) then
        call endrun (.true.)
      end if
!
end subroutine
!#######################################################################
subroutine endrun (ifstop)
!
!-----------------------------------------------------------------------
!
! ****** End the run and exit the code.
!
!-----------------------------------------------------------------------
!
      use mpidefs
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      logical :: ifstop
!
!-----------------------------------------------------------------------
!
! ****** MPI error return.
!
      integer :: ierr
!
!-----------------------------------------------------------------------
!
! ****** Close the output file.
!
      if (iamp0) then
        close (9)
      end if
!
! ****** Exit MPI gracefully.
!
      call MPI_Finalize (ierr)
!
! ****** Call the STOP statement if requested.
!
      if (ifstop) then
        stop
      end if
!
end subroutine
!#######################################################################
subroutine init_mpi
!
!-----------------------------------------------------------------------
!
! ****** Initialize MPI.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
! ****** MPI error return.
!
      integer :: ierr,tcheck
!
!-----------------------------------------------------------------------
!
! ****** Real number to determine the KIND of REALs.
!
      real(r_typ) :: def_real
!
!-----------------------------------------------------------------------
!
      call MPI_Init_thread (MPI_THREAD_FUNNELED,tcheck,ierr)
!
! ****** Get the total number of processors.
!
      call MPI_Comm_size (MPI_COMM_WORLD,nproc,ierr)
!
! ****** Get the index (rank) of the local processor in
! ****** communicator MPI_COMM_WORLD in variable IPROCW.
!
      call MPI_Comm_rank (MPI_COMM_WORLD,iprocw,ierr)
!
! ****** Create a shared communicator for all ranks in the node.
!
      call MPI_Comm_split_type (MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0, &
                                MPI_INFO_NULL,comm_shared,ierr)
!
! ****** Get the total number of processors in node.
!
      call MPI_Comm_size (comm_shared,nprocsh,ierr)
!
! ****** Get the index (rank) of the local processor in the local node.
!
      call MPI_Comm_rank (comm_shared,iprocsh,ierr)
!
! ****** Set the flag to designate whether this processor
! ****** has rank 0 in communicator MPI_COMM_WORLD.
!
      if (iprocw.eq.0) then
        iamp0=.true.
      else
        iamp0=.false.
      end if
!
! ****** Set the number type for communicating REAL
! ****** numbers in MPI calls.
!
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
!
! ****** Set the GPU device number based on local rank.
! ****** NOTE! This assumes than #GPUs per node = #MPI ranks per node.
!
!$    call omp_set_default_device (iprocsh)
!$acc set device_num(iprocsh)
!
end subroutine
!#######################################################################
subroutine check_input
!
!-----------------------------------------------------------------------
!
! ****** Check the validity of the input parameters.
!
!-----------------------------------------------------------------------
!
      use number_types
      use vars
      use solve_params
      use mpidefs
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
!
!-----------------------------------------------------------------------
!
! ****** Check that OPTION is valid and set the boundary condition
! ****** switches accordingly.
!
      if (option.eq.'potential') then
!
! ****** For a potential field, set d(phi)/dr to zero at r=R1
! ****** (i.e., the field is tangential to the boundary).
!
        pm_r1=one
!
      else if (option.eq.'open') then
!
! ****** For an open field, set phi to zero at r=R1
! ****** (i.e., the field is radial there).
!
        pm_r1=-one
!
      else if (option.eq.'ss') then
!
! ****** For a source surface field, set phi to zero at r=R1
! ****** (i.e., the field is radial there).
!
        pm_r1=-one
!
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
!
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
!
end subroutine
!#######################################################################
subroutine set_proc_topology
!
!-----------------------------------------------------------------------
!
! ****** Set the optimal values of the MPI rank topology
! ****** in dimensions not set by user.
!
!-----------------------------------------------------------------------
!
      use mpidefs
      use decomposition_params
      use number_types
      use global_dims
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1.0_r_typ
      real(r_typ), parameter :: zero=0.0_r_typ
      real(r_typ), parameter :: bigval=HUGE(1.0_r_typ)
!
!-----------------------------------------------------------------------
!
      integer, dimension(:), allocatable :: factors
      integer, dimension(:,:), allocatable :: rank_factors
      real(r_typ), dimension(:,:), allocatable :: nperrank
      real(r_typ), dimension(:), allocatable :: penalty
!
      integer :: i,j,k,num_fac,num_rank_fac,best_idx
      integer :: fr=1,ft=1,fp=1
      real(r_typ) :: a12,a13,a23
!
!-----------------------------------------------------------------------
!
! ****** Extract nproc values.  A value of -1 indicates the dimension
! ****** should be autoset.
!
      nproc_r=nprocs(1)
      nproc_t=nprocs(2)
      nproc_p=nprocs(3)
!
! ****** If no dimensions are to be autoset, return.
!
      if(nproc_r.ne.-1.and.nproc_t.ne.-1.and.nproc_p.ne.-1) return
!
! ****** Get all factors of nproc and store them in factors array.
!
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
!
! ****** Set penalty function parameters and any fixed dimensions
! ****** based on which dimensions are to be autoset.
!
      a12=one
      a13=one
      a23=one
!
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
!
! ****** Loop over all combinations of factors and save those that
! ****** yield the correct number of MPI ranks into rank_factors array.
!
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
!
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
!
      allocate(rank_factors(num_rank_fac,3))
      allocate(nperrank(num_rank_fac,3))
      allocate(penalty(num_rank_fac))
!
      rank_factors(:,:)=-1
      penalty(:)=bigval
!
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
!
! ****** Get number of grid points per rank for each dimension.
!
      nperrank(:,1)=real(nr_g)/rank_factors(:,1)
      nperrank(:,2)=real(nt_g)/rank_factors(:,2)
      nperrank(:,3)=real(np_g)/rank_factors(:,3)
!
! ****** Compute penalty function.
!
      penalty(:)=a12*(nperrank(:,1)-nperrank(:,2))**2 &
                +a23*(nperrank(:,2)-nperrank(:,3))**2 &
                +a13*(nperrank(:,3)-nperrank(:,1))**2
!
! ****** Eliminate any choices that yield less than a minimum number
! ****** of grid points per rank.
!
      do i=1,num_rank_fac
        if (nperrank(i,1).lt.4) penalty(i)=bigval
        if (nperrank(i,2).lt.4) penalty(i)=bigval
        if (nperrank(i,3).lt.3) penalty(i)=bigval
      enddo
!
! ****** Find optimal topology.
!
      best_idx=MINLOC(penalty,1)
!
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
!
! ****** Set optimal topology.
!
      nprocs(1)=rank_factors(best_idx,1)
      nprocs(2)=rank_factors(best_idx,2)
      nprocs(3)=rank_factors(best_idx,3)
!
      deallocate(factors)
      deallocate(rank_factors)
      deallocate(nperrank)
      deallocate(penalty)
!
end subroutine
!#######################################################################
subroutine check_proc_topology
!
!-----------------------------------------------------------------------
!
! ****** Check the validity of the requested processor topology.
!
!-----------------------------------------------------------------------
!
      use mpidefs
      use decomposition_params
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: i,nreq
!
!-----------------------------------------------------------------------
!
! ****** Check the processor topology.
!
      do i=1,3
        if (nprocs(i).lt.1.and.nprocs(i).ne.-1) then
          if (iamp0) then
            write (*,*)
            write (*,*) '### ERROR in CHECK_PROC_TOPOLOGY:'
            write (*,*) '### Processor topology specification error.'
            write (*,*) 'Invalid number of processors requested.'
            write (*,*) 'Dimension = ',i
            write (*,*) 'Number of processors requested = ', &
                        nprocs(i)
          end if
          call endrun (.true.)
        end if
      enddo
!
! ****** Set the optimal values of the topology if requested.
!
      call set_proc_topology
!
! ****** Check that the number of processors available
! ****** matches the number requested.
!
      nreq=nprocs(1)*nprocs(2)*nprocs(3)
!
      if (nreq.ne.nproc) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in CHECK_PROC_TOPOLOGY:'
          write (*,*) '### Processor topology specification error.'
          write (*,*) 'The number of processors requested does not'// &
                     ' equal the number available.'
          write (*,*) 'Number of processors requested = ',nreq
          write (*,*) 'Number of processors available = ',nproc
        end if
        call endrun (.true.)
      end if
!
end subroutine
!#######################################################################
subroutine decompose_domain
!
!-----------------------------------------------------------------------
!
! ****** Decompose the domain into a Cartesian MPI topology.
!
!-----------------------------------------------------------------------
!
      use mpidefs
      use decomposition_params
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      integer, parameter :: ndim=3
      integer, dimension(ndim) :: coords
      logical, dimension(ndim) :: periodic
      logical :: reorder
      logical, dimension(ndim) :: keep_dim
!
!-----------------------------------------------------------------------
!
! ****** Create a communicator over all processors, COMM_ALL,
! ****** that has a Cartesian topology.
!
! ****** Specify the periodicity of the coordinate system.
!
      periodic(1)=.false.
      periodic(2)=.false.
      periodic(3)=.true.
!
! ****** Allow re-ordering in the Cartesian topology.
!
      reorder=.true.
!
      call MPI_Cart_create (MPI_COMM_WORLD,ndim,nprocs, &
                            periodic,reorder,comm_all,ierr)
!
! ****** Get the index (rank) of the local processor in
! ****** communicator COMM_ALL in variable IPROC.
!
! ****** IMPORTANT NOTE:
! ****** If re-odering was allowed in the Cartesian topology
! ****** creation (above), then the rank of the local processor
! ****** in communicator COMM_ALL may be different from its rank
! ****** in communicator MPI_COMM_WORLD.
!
      call MPI_Comm_rank (comm_all,iproc,ierr)
!
! ****** Set the processor rank IPROC0 in communicator COMM_ALL
! ****** for the processor that has rank 0 in MPI_COMM_WORLD.
! ****** This value is broadcast to all the processors.
!
      if (iamp0) then
        iproc0=iproc
      end if
      call MPI_Bcast (iproc0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!
! ****** Get the coordinate indices of this processor in the
! ****** Cartesian MPI topology.
!
      call MPI_Cart_coords (comm_all,iproc,ndim,coords,ierr)
!
      iproc_r=coords(1)
      iproc_t=coords(2)
      iproc_p=coords(3)
!
      nproc_r=nprocs(1)
      nproc_t=nprocs(2)
      nproc_p=nprocs(3)
!
! ****** Get the rank of the neighboring processors in the
! ****** Cartesian MPI topology.
!
      call MPI_Cart_shift (comm_all,0,1,iproc_rm,iproc_rp,ierr)
      call MPI_Cart_shift (comm_all,1,1,iproc_tm,iproc_tp,ierr)
      call MPI_Cart_shift (comm_all,2,1,iproc_pm,iproc_pp,ierr)
!
! ****** Create communicators for operations involving all
! ****** processors in the phi dimension.  These communicators
! ****** are stored in COMM_PHI (and generally represent different
! ****** communicators on different processors).
!
      keep_dim(1)=.false.
      keep_dim(2)=.false.
      keep_dim(3)=.true.
!
      call MPI_Cart_sub (comm_all,keep_dim,comm_phi,ierr)
!
! ****** Create communicators for operations involving
! ****** all processors in the r dimension.
! ****** These communicators are stored in COMM_R
! ****** (and generally represent different communicators on
! ****** different processors).
!
      keep_dim(1)=.true.
      keep_dim(2)=.false.
      keep_dim(3)=.false.
!
      call MPI_Cart_sub (comm_all,keep_dim,comm_r,ierr)
!
end subroutine
!#######################################################################
subroutine decompose_mesh_r
!
!-----------------------------------------------------------------------
!
! ****** Decompose the r mesh between processors.
!
!-----------------------------------------------------------------------
!
      use global_dims
      use local_dims
      use decomposition
      use solve_params
      use mpidefs
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: ierr,i,npts
      integer :: i0_h,i1_h,i0_m,i1_m
      integer, dimension(nproc_r) :: mp_r
!
!-----------------------------------------------------------------------
!
! ****** Decompose the r dimension.
!
      npts=nr_g
!
      call decompose_dimension (npts,nproc_r,mp_r,ierr)
      if (ierr.ne.0) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in DECOMPOSE_MESH:'
          write (*,*) '### Anomaly in decomposing the mesh'// &
                      ' between processors.'
          write (*,*) '### Could not decompose the r mesh.'
          write (*,*) 'Number of mesh points in r = ',npts
          write (*,*) 'Number of processors along r = ',nproc_r
        end if
        call endrun (.true.)
      end if
!
! ****** Check that the resulting mesh topology is valid.
!
      call check_mesh_topology (nproc_r,mp_r,1,'r')
!
! ****** Compute the mapping between the processor decomposition
! ****** and the global mesh.
!
! ****** Note that there is a two-point overlap in the mesh
! ****** between adjacent processors in r.
!
      i0_g=1
      do i=1,iproc_r
        i0_g=i0_g+mp_r(i)
      enddo
      nr=mp_r(iproc_r+1)+2
      i1_g=i0_g+nr-1
!
      nrm1=nr-1
!
! ****** Set the flags to indicate whether this processor has
! ****** points on the physical boundaries.
!
      if (iproc_r.eq.0) then
        rb0=.true.
      else
        rb0=.false.
      end if
!
      if (iproc_r.eq.nproc_r-1) then
        rb1=.true.
      else
        rb1=.false.
      end if
!
! ****** Set the dimensions of arrays for the "main" meshes
! ****** (i.e., the "m" mesh) for which normal derivatives are
! ****** needed (e.g., v).  These vary on different processors,
! ****** depending if they are left-boundary, internal, or
! ****** right-boundary processors.
!
      if (rb1) then
        nrm=nrm1
      else
        nrm=nr
      end if
!
! ****** Store the mapping structure (for this processor).
!
      allocate (map_rh(0:nproc-1))
      allocate (map_rm(0:nproc-1))
!
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
!
      if (rb0) then
        i0_m=1
      else
        i0_m=2
      end if
      i1_m=nrm1
!
      map_rh(iproc)%i0=i0_h
      map_rh(iproc)%i1=i1_h
!
      map_rm(iproc)%i0=i0_m
      map_rm(iproc)%i1=i1_m
!
      map_rh(iproc)%offset=i0_g+map_rh(iproc)%i0-1
      map_rh(iproc)%n=map_rh(iproc)%i1-map_rh(iproc)%i0+1
!
      map_rm(iproc)%offset=i0_g+map_rm(iproc)%i0-1
      map_rm(iproc)%n=map_rm(iproc)%i1-map_rm(iproc)%i0+1
!
! ****** Gather the mapping information by communicating among
! ****** all processors.
!
      call gather_mapping_info (map_rh)
      call gather_mapping_info (map_rm)
!
end subroutine
!#######################################################################
subroutine decompose_mesh_tp
!
!-----------------------------------------------------------------------
!
! ****** Decompose the theta and phi mesh between processors.
!
!-----------------------------------------------------------------------
!
      use global_dims
      use local_dims
      use decomposition
      use mpidefs
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: ierr,j,k
      integer :: j0_h,j1_h,j0_m,j1_m
      integer :: k0_h,k1_h,k0_m,k1_m
      integer, dimension(nproc_t) :: mp_t
      integer, dimension(nproc_p) :: mp_p
!
!-----------------------------------------------------------------------
!
! ****** Decompose the t dimension.
!
      call decompose_dimension (nt_g,nproc_t,mp_t,ierr)
      if (ierr.ne.0) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in DECOMPOSE_MESH_TP:'
          write (*,*) '### Anomaly in decomposing the mesh'// &
                      ' between processors.'
          write (*,*) '### Could not decompose the theta mesh.'
          write (*,*) 'Number of mesh points in theta = ',nt_g
          write (*,*) 'Number of processors along theta = ',nproc_t
        end if
        call endrun (.true.)
      end if
!
! ****** Decompose the p dimension.
!
      call decompose_dimension (np_g,nproc_p,mp_p,ierr)
      if (ierr.ne.0) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in DECOMPOSE_MESH_TP:'
          write (*,*) '### Anomaly in decomposing the mesh'// &
                      ' between processors.'
          write (*,*) '### Could not decompose the phi mesh.'
          write (*,*) 'Number of mesh points in phi = ',np_g
          write (*,*) 'Number of processors along phi = ',nproc_p
        end if
        call endrun (.true.)
      end if
!
! ****** Check that the resulting mesh topology is valid.
!
      call check_mesh_topology (nproc_t,mp_t,1,'theta')
      call check_mesh_topology (nproc_p,mp_p,1,'phi')
!
! ****** Set the flag for an axisymmetric run (requested by
! ****** setting NP_G = 3).
!
      if (np_g.eq.3) then
        axisymmetric=.true.
      else
        axisymmetric=.false.
      end if
!
! ****** Compute the mapping between the processor decomposition
! ****** and the global mesh.
!
! ****** Note that there is a two-point overlap in the mesh
! ****** between adjacent processors in theta and phi.
!
      j0_g=1
      do j=1,iproc_t
        j0_g=j0_g+mp_t(j)
      enddo
      nt=mp_t(iproc_t+1)+2
      j1_g=j0_g+nt-1
!
      k0_g=1
      do k=1,iproc_p
        k0_g=k0_g+mp_p(k)
      enddo
      np=mp_p(iproc_p+1)+2
      k1_g=k0_g+np-1
!
      ntm1=nt-1
      npm1=np-1
!
! ****** Set the flags to indicate whether this processor has
! ****** points on the physical boundaries.
!
      if (iproc_t.eq.0) then
        tb0=.true.
      else
        tb0=.false.
      end if
!
      if (iproc_t.eq.nproc_t-1) then
        tb1=.true.
      else
        tb1=.false.
      end if
!
! ****** Set the dimensions of arrays for the "main" meshes
! ****** (i.e., the "m" mesh) for which normal derivatives are
! ****** needed (e.g., v).  These vary on different processors,
! ****** depending if they are left-boundary, internal, or
! ****** right-boundary processors.
!
      if (tb1) then
        ntm=ntm1
      else
        ntm=nt
      end if
!
! ****** Since the phi dimension is periodic, all processors
! ****** have the same mesh limits.
!
      npm=np
!
! ****** Store the mapping structure (for this processor).
!
      allocate (map_th(0:nproc-1))
      allocate (map_tm(0:nproc-1))
      allocate (map_ph(0:nproc-1))
      allocate (map_pm(0:nproc-1))
!
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
!
      if (tb0) then
        j0_m=1
      else
        j0_m=2
      end if
      j1_m=ntm1
!
      if (iproc_p.eq.0) then
        k0_m=1
      else
        k0_m=2
      end if
      k1_m=npm1
!
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
!
      map_th(iproc)%i0=j0_h
      map_th(iproc)%i1=j1_h
!
      map_tm(iproc)%i0=j0_m
      map_tm(iproc)%i1=j1_m
!
      map_ph(iproc)%i0=k0_h
      map_ph(iproc)%i1=k1_h
!
      map_pm(iproc)%i0=k0_m
      map_pm(iproc)%i1=k1_m
!
      map_th(iproc)%offset=j0_g+map_th(iproc)%i0-1
      map_th(iproc)%n=map_th(iproc)%i1-map_th(iproc)%i0+1
!
      map_tm(iproc)%offset=j0_g+map_tm(iproc)%i0-1
      map_tm(iproc)%n=map_tm(iproc)%i1-map_tm(iproc)%i0+1
!
      map_ph(iproc)%offset=k0_g+map_ph(iproc)%i0-1
      map_ph(iproc)%n=map_ph(iproc)%i1-map_ph(iproc)%i0+1
!
      map_pm(iproc)%offset=k0_g+map_pm(iproc)%i0-1
      map_pm(iproc)%n=map_pm(iproc)%i1-map_pm(iproc)%i0+1
!
! ****** Gather the mapping information by communicating among
! ****** all processors.
!
      call gather_mapping_info (map_th)
      call gather_mapping_info (map_tm)
      call gather_mapping_info (map_ph)
      call gather_mapping_info (map_pm)
!
end subroutine
!#######################################################################
subroutine check_mesh_dimensions (nr_g,nt_g,np_g)
!
!-----------------------------------------------------------------------
!
! ****** Check that the requested (global) mesh dimensions are
! ****** valid.
!
!-----------------------------------------------------------------------
!
      use mpidefs
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: nr_g,nt_g,np_g
!
!-----------------------------------------------------------------------
!
      if (nr_g.lt.4) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in CHECK_MESH_DIMENSIONS:'
          write (*,*) '### Invalid number of r mesh points'// &
                    ' requested.'
          write (*,*) '### The minimum number of mesh points is 4.'
          write (*,*)
          write (*,*) 'Number of mesh points requested = ',nr_g
        end if
        call endrun (.true.)
      end if
!
      if (nt_g.lt.4) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in CHECK_MESH_DIMENSIONS:'
          write (*,*) '### Invalid number of theta mesh points'// &
                    ' requested.'
          write (*,*) '### The minimum number of mesh points is 4.'
          write (*,*)
          write (*,*) 'Number of mesh points requested = ',nt_g
        end if
        call endrun (.true.)
      end if
!
      if (np_g.lt.3) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### ERROR in CHECK_MESH_DIMENSIONS:'
          write (*,*) '### Invalid number of phi mesh points'// &
                    ' requested.'
          write (*,*) '### The minimum number of mesh points is 3.'
          write (*,*)
          write (*,*) 'Number of mesh points requested = ',np_g
        end if
        call endrun (.true.)
      end if
!
end subroutine
!#######################################################################
subroutine decompose_dimension (nx,np,mp,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Decompose the mesh points NX along NP processors.
!
! ****** The decomposed mesh points are returned in array MP.
!
!-----------------------------------------------------------------------
!
! ****** This routine attempts to assign the mesh points as equally
! ****** as possible between the processors.
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: nx
      integer :: np
      integer, dimension(np) :: mp
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      integer :: nxm2,mpav,nrem
!
!-----------------------------------------------------------------------
!
      ierr=0
!
      nxm2=nx-2
!
      if (nxm2.le.0) then
        ierr=1
        return
      end if
!
      if (np.le.0) then
        ierr=2
        return
      end if
!
      mpav=nxm2/np
!
      mp(:)=mpav
!
      nrem=nxm2-mpav*np
!
      mp(1:nrem)=mp(1:nrem)+1
!
end subroutine
!#######################################################################
subroutine check_mesh_topology (np,mp,min_pts,coord)
!
!-----------------------------------------------------------------------
!
! ****** Check the validity of the requested mesh topology.
!
!-----------------------------------------------------------------------
!
      use mpidefs
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: np
      integer, dimension(np) :: mp
      integer :: min_pts
      character(*) :: coord
!
!-----------------------------------------------------------------------
!
      integer :: i
!
!-----------------------------------------------------------------------
!
! ****** Check that the number of mesh points on each processor
! ****** is valid.
!
      do i=1,np
        if (mp(i).lt.min_pts) then
          if (iamp0) then
            write (*,*)
            write (*,*) '### ERROR in CHECK_MESH_TOPOLOGY:'
            write (*,*) '### Mesh topology specification error.'
            write (*,*) 'Invalid number of ',coord, &
                        ' mesh points requested.'
            write (*,*) 'Processor index = ',i
            write (*,*) 'Number of mesh points requested = ',mp(i)
            write (*,*) 'Minimum number of mesh points allowed = ', &
                        min_pts
          end if
          call endrun (.true.)
        end if
      enddo
!
end subroutine
!#######################################################################
subroutine gather_mapping_info (map)
!
!-----------------------------------------------------------------------
!
! ****** Gather a mapping information array by communicating
! ****** among all processors.
!
!-----------------------------------------------------------------------
!
      use mpidefs
      use decomposition
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(map_struct), dimension(0:nproc-1) :: map
!
!-----------------------------------------------------------------------
!
! ****** Buffers for packing the data.
!
      integer, parameter :: lbuf=4
      integer, dimension(lbuf) :: sbuf
      integer, dimension(lbuf,0:nproc-1) :: rbuf
!
!-----------------------------------------------------------------------
!
      integer :: ierr,irank
!
!-----------------------------------------------------------------------
!
! ****** Put the local section of the mapping information
! ****** array into a buffer.
!
      sbuf(1)=map(iproc)%n
      sbuf(2)=map(iproc)%i0
      sbuf(3)=map(iproc)%i1
      sbuf(4)=map(iproc)%offset
!
! ****** Communicate among all processors.  After this call, all
! ****** processors have the complete mapping information.
!
      call MPI_Allgather (sbuf,lbuf,MPI_INTEGER, &
                          rbuf,lbuf,MPI_INTEGER,comm_all,ierr)
!
! ****** Extract the mapping information from the buffer.
!
      do irank=0,nproc-1
        map(irank)%n     =rbuf(1,irank)
        map(irank)%i0    =rbuf(2,irank)
        map(irank)%i1    =rbuf(3,irank)
        map(irank)%offset=rbuf(4,irank)
      enddo
!
end subroutine
!#######################################################################
subroutine decomp_diags
!
!-----------------------------------------------------------------------
!
! ****** Print diagnostics about the mesh decomposition.
!
!-----------------------------------------------------------------------
!
      use global_dims
      use global_mesh
      use local_dims
      use local_mesh
      use mpidefs
      use solve_params
      use debug
      use decomposition
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: ierr
      integer :: irank
      real(r_typ) :: n_per_grid_min,n_per_grid_max
!
!-----------------------------------------------------------------------
!
      if (iamp0) then
!
        n_per_grid_min=floor(real(nr_g)/nproc_r) &
                      *floor(real(nt_g)/nproc_t) &
                      *floor(real(np_g)/nproc_p)
!
        n_per_grid_max=ceiling(real(nr_g)/nproc_r) &
                      *ceiling(real(nt_g)/nproc_t) &
                      *ceiling(real(np_g)/nproc_p)
!
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
        write (*,'(A,F6.1)') ' Average # of r mesh pts per rank = ', &
                     real(nr_g)/nproc_r
        write (*,'(A,F6.1)') ' Average # of t mesh pts per rank = ', &
                     real(nt_g)/nproc_t
        write (*,'(A,F6.1)') ' Average # of p mesh pts per rank = ', &
                     real(np_g)/nproc_p
        write (*,*)
        write (*,'(A,F6.2,A)') ' Estimated load imbalance = ',       &
            100.0*(1.0-real(n_per_grid_min)/real(n_per_grid_max)),' %'
!
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
        write (9,'(A,F6.1)') ' Average # of r mesh pts per rank = ', &
                     real(nr_g)/nproc_r
        write (9,'(A,F6.1)') ' Average # of t mesh pts per rank = ', &
                     real(nt_g)/nproc_t
        write (9,'(A,F6.1)') ' Average # of p mesh pts per rank = ', &
                     real(np_g)/nproc_p
        write (9,*)
        write (9,'(A,F6.2,A)') ' Estimated load imbalance = ',       &
            100.0*(1.0-real(n_per_grid_min)/real(n_per_grid_max)),' %'
!
      end if
!
      if (idebug.le.1) return
!
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
!
end subroutine
!#######################################################################
subroutine allocate_global_arrays
!
!-----------------------------------------------------------------------
!
! ****** Allocate global arrays.
!
!-----------------------------------------------------------------------
!
      use global_dims
      use global_mesh
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
! ****** Allocate global mesh arrays.
!
      allocate (r_g (nrm1_g))
      allocate (dr_g(nrm1_g))
!
      allocate (rh_g (nr_g))
      allocate (drh_g(nr_g))
!
      allocate (t_g (ntm1_g))
      allocate (dt_g(ntm1_g))
!
      allocate (th_g (nt_g))
      allocate (dth_g(nt_g))
!
      allocate (p_g  (np_g))
      allocate (dp_g (np_g))
      allocate (ph_g (np_g))
      allocate (dph_g(np_g))
!
      allocate (st_g(ntm1_g))
      allocate (ct_g(ntm1_g))
!
      allocate (sth_g(nt_g))
      allocate (cth_g(nt_g))
!
      allocate (sp_g (np_g))
      allocate (cp_g (np_g))
      allocate (sph_g(np_g))
      allocate (cph_g(np_g))
!
end subroutine
!#######################################################################
subroutine allocate_local_arrays_r
!
!-----------------------------------------------------------------------
!
! ****** Allocate local arrays for the r dimension and 3D arrays.
!
!-----------------------------------------------------------------------
!
      use local_dims
      use local_mesh
      use fields
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      allocate (r (nrm))
      allocate (r2 (nrm))
      allocate (dr(nrm))
!
      allocate (rh (nr))
      allocate (drh(nr))
!
! ****** Allocate inverse quantities.
!
      allocate (r_i (nrm))
      allocate (dr_i(nrm))
!
      allocate (rh_i (nr))
      allocate (drh_i(nr))
!
! ****** Allocate the potential array and cg scratch array.
!
      allocate (phi(nr,nt,np))
      allocate (x_ax(nr,nt,np))
      phi(:,:,:)=0.
      x_ax(:,:,:)=0.
!$omp target enter data map(to:phi,x_ax)
!
! ****** Allocate polar boundary arrays.
!
      allocate (sum0(nr))
      allocate (sum1(nr))
      sum0(:)=0.
      sum1(:)=0.
!$omp target enter data map(to:sum0,sum1)
!
! ****** Allocate the local magnetic field arrays.
!
      allocate (br(nrm,nt,np))
      allocate (bt(nr,ntm,np))
      allocate (bp(nr,nt,npm))
      br(:,:,:)=0.
      bt(:,:,:)=0.
      bp(:,:,:)=0.
!
end subroutine
!#######################################################################
subroutine allocate_local_arrays_tp
!
!-----------------------------------------------------------------------
!
! ****** Allocate local arrays for the theta and phi dimensions.
!
!-----------------------------------------------------------------------
!
      use local_dims
      use local_mesh
      use fields
      use vars
      use solve_params
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      allocate (t (ntm))
      allocate (dt(ntm))
!
      allocate (th (nt))
      allocate (dth(nt))
!
      allocate (p (np))
      allocate (dp(np))
!
      allocate (ph (np))
      allocate (dph(np))
!
      allocate (st(ntm))
      allocate (ct(ntm))
!
      allocate (sth(nt))
      allocate (cth(nt))
!
      allocate (sp (np))
      allocate (cp (np))
      allocate (sph(np))
      allocate (cph(np))
!
! ****** Allocate inverse quantities.
!
      allocate (dt_i(ntm))
      allocate (st_i(ntm))
!
      allocate (dth_i(nt))
      allocate (sth_i(nt))
!
      allocate (dp_i (np))
      allocate (dph_i(np))
!
! ****** Allocate the boundary radial magnetic field array(s).
!
      allocate (br0(nt,np))
      br0(:,:)=0.
!$omp target enter data map(to:br0)
!
      if (validation_run) then
        allocate (br1(nt,np))
        br1(:,:)=0.
!$omp target enter data map(to:br1)
      end if
!
end subroutine
!#######################################################################
subroutine set_global_mesh
!
!-----------------------------------------------------------------------
!
! ****** Define the global mesh arrays.
!
!-----------------------------------------------------------------------
!
      use number_types
      use global_dims
      use global_mesh
      use meshdef
      use mpidefs
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: zero=0._r_typ
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: half=.5_r_typ
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k
!
!-----------------------------------------------------------------------
!
! ****** Define the radial mesh.
!
      call genmesh (9,'r',nrm1_g,r0,r1,nmseg,rfrac,drratio,nfrmesh, &
                    .false.,zero,r_g)
!
      do i=2,nrm1_g
        rh_g(i)=.5*(r_g(i)+r_g(i-1))
        drh_g(i)=r_g(i)-r_g(i-1)
      enddo
      rh_g(1)=rh_g(2)-drh_g(2)
      rh_g(nr_g)=rh_g(nrm1_g)+drh_g(nrm1_g)
      drh_g(1)=drh_g(2)
      drh_g(nr_g)=drh_g(nrm1_g)
!
      do i=1,nrm1_g
        dr_g(i)=rh_g(i+1)-rh_g(i)
      enddo
!
! ****** Define the theta mesh.
!
      call genmesh (9,'t',ntm1_g,t0,t1,nmseg,tfrac,dtratio, &
                    nftmesh,.false.,zero,t_g)
!
      do j=2,ntm1_g
        th_g(j)=.5*(t_g(j)+t_g(j-1))
        dth_g(j)=t_g(j)-t_g(j-1)
      enddo
      th_g(1)=th_g(2)-dth_g(2)
      th_g(nt_g)=th_g(ntm1_g)+dth_g(ntm1_g)
      dth_g(1)=dth_g(2)
      dth_g(nt_g)=dth_g(ntm1_g)
!
      do j=1,ntm1_g
        dt_g(j)=th_g(j+1)-th_g(j)
      enddo
!
! ****** Define the periodic phi mesh.
!
      call genmesh (9,'p',npm1_g,p0,p1,nmseg,pfrac,dpratio, &
                    nfpmesh,.true.,phishift,p_g)
      p_g(np_g)=p_g(2)+pl
!
      do k=2,np_g
        ph_g(k)=half*(p_g(k)+p_g(k-1))
        dph_g(k)=p_g(k)-p_g(k-1)
      enddo
      ph_g(1)=ph_g(npm1_g)-pl
      dph_g(1)=dph_g(npm1_g)
!
      do k=1,npm1_g
        dp_g(k)=ph_g(k+1)-ph_g(k)
      enddo
      dp_g(np_g)=dp_g(2)
!
! ****** Enforce exact periodicity to protect symmetry properties
! ****** from round-off errors (especially for axisymmetric cases).
!
      dph_g(np_g)=dph_g(2)
      dp_g(1)=dp_g(npm1_g)
!
! ****** Define global auxiliary mesh-related arrays.
!
      st_g(:)=sin(t_g(:))
      ct_g(:)=cos(t_g(:))
      sth_g(:)=sin(th_g(:))
      cth_g(:)=cos(th_g(:))
!
      sp_g(:)=sin(p_g(:))
      cp_g(:)=cos(p_g(:))
      sph_g(:)=sin(ph_g(:))
      cph_g(:)=cos(ph_g(:))
!
! ****** For an axisymmetric case, set the exact values of
! ****** sin(phi) and cos(phi) to preserve symmetry properties
! ****** in the presence of round-off errors.
!
      if (axisymmetric) then
        sp_g(2)=0.
        cp_g(2)=one
        sph_g(2)=0.
        cph_g(2)=-one
      end if
!
! ****** Enforce exact periodicity to protect symmetry properties
! ****** from round-off errors (especially for axisymmetric cases).
!
      sph_g(1)=sph_g(npm1_g)
      sph_g(np_g)=sph_g(2)
      cph_g(1)=cph_g(npm1_g)
      cph_g(np_g)=cph_g(2)
      sp_g(1)=sp_g(npm1_g)
      sp_g(np_g)=sp_g(2)
      cp_g(1)=cp_g(npm1_g)
      cp_g(np_g)=cp_g(2)
!
end subroutine
!#######################################################################
subroutine set_local_mesh_r
!
!-----------------------------------------------------------------------
!
! ****** Define the local r mesh arrays.
!
!-----------------------------------------------------------------------
!
      use number_types
      use global_dims
      use global_mesh
      use local_dims
      use local_mesh
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
!
!-----------------------------------------------------------------------
!
      integer :: i
!
!-----------------------------------------------------------------------
!
! ****** Define the local meshes.
!
      do i=1,nrm
        r(i)=r_g(i0_g+i-1)
        dr(i)=dr_g(i0_g+i-1)
      enddo
!
      dr1=dr(1)
      drn=dr(nrm)
!
      do i=1,nr
        rh(i)=rh_g(i0_g+i-1)
        drh(i)=drh_g(i0_g+i-1)
      enddo
!
! ****** Define local auxiliary mesh-related arrays.
!
      r2(:)=r(:)**2
      r_i(:)=one/r(:)
      dr_i(:)=one/dr(:)
      rh_i(:)=one/rh(:)
      drh_i(:)=one/drh(:)
!
!$omp target enter data map(to:r,r2,r_i,dr,dr_i,dr1,drn,rh,rh_i,drh,drh_i)
end subroutine
!#######################################################################
subroutine set_local_mesh_tp
!
!-----------------------------------------------------------------------
!
! ****** Define the local theta and phi mesh arrays.
!
!-----------------------------------------------------------------------
!
      use number_types
      use global_dims
      use global_mesh
      use local_dims
      use local_mesh
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
!
!-----------------------------------------------------------------------
!
      integer :: j,k,j0,j1
!
!-----------------------------------------------------------------------
!
! ****** Define the local meshes.
!
      do j=1,ntm
        t(j)=t_g(j0_g+j-1)
        dt(j)=dt_g(j0_g+j-1)
      enddo
!
      do j=1,nt
        th(j)=th_g(j0_g+j-1)
        dth(j)=dth_g(j0_g+j-1)
      enddo
!
      do k=1,npm
        p(k)=p_g(k0_g+k-1)
        dp(k)=dp_g(k0_g+k-1)
      enddo
!
      do k=1,np
        ph(k)=ph_g(k0_g+k-1)
        dph(k)=dph_g(k0_g+k-1)
      enddo
!
! ****** Define local auxiliary mesh-related arrays.
!
      do j=1,ntm
        st(j)=st_g(j0_g+j-1)
        ct(j)=ct_g(j0_g+j-1)
      enddo
!
      do j=1,nt
        sth(j)=sth_g(j0_g+j-1)
        cth(j)=cth_g(j0_g+j-1)
      enddo
!
      do k=1,npm
        sp(k)=sp_g(k0_g+k-1)
        cp(k)=cp_g(k0_g+k-1)
      enddo
!
      do k=1,np
        sph(k)=sph_g(k0_g+k-1)
        cph(k)=cph_g(k0_g+k-1)
      enddo
!
      dt_i(:)=one/dt(:)
      dth_i(:)=one/dth(:)
      sth_i(:)=one/sth(:)
      dp_i(:)=one/dp(:)
      dph_i(:)=one/dph(:)
!
! ****** Prevent division by zero at the poles for sin(theta).
!
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
!
      st_i(:)=0.
      do j=j0,j1
        st_i(j)=one/st(j)
      enddo
!
!$omp target enter data map(to:t,th,dt,dth,p,ph,dp,dph,st,ct,sth,cth, &
!$omp                   sp,cp,sph,cph,dt_i,dth_i,st_i,sth_i,dp_i,dph_i)
end subroutine
!#######################################################################
subroutine genmesh (io,label,nc,c0,c1,nseg,frac,dratio,nfilt,periodic, &
                    c_shift,c)
!
!-----------------------------------------------------------------------
!
! ****** Generate a one-dimensional mesh.
!
!-----------------------------------------------------------------------
!
! ****** Input arguments:
!
!          IO      : [integer]
!                    Fortran file unit number to which to write
!                    mesh diagnostics.  Set IO=0 if diagnostics
!                    are not of interest.  It is assumed that
!                    unit IO has been connected to a file prior
!                    to calling this routine.
!
!          LABEL   : [character(*)]
!                    Name for the mesh coordinate (example: 'x').
!
!          NC      : [integer]
!                    Number of mesh points to load.
!
!          C0      : [real]
!                    The starting location for the coordinate.
!
!          C1      : [real]
!                    The ending location for the coordinate.
!                    It is required that C1.gt.C0.
!
!          NSEG    : [integer]
!                    Maximum number of mesh segments.
!                    The mesh spacing in each segment varies
!                    exponentially with a uniform amplification
!                    factor.  The actual number of mesh segments
!                    used is NSEG or less.  It is obtained from the
!                    information in array FRAC.
!
!          FRAC    : [real array, dimension NSEG]
!                    The normalized positions of the mesh segment
!                    boundaries (as a fraction of the size of the
!                    domain).  For a non-periodic mesh, the first
!                    value of FRAC specified must equal 0. and the
!                    last value must equal 1.  For a periodic mesh,
!                    FRAC must not contain both 0. and 1., since
!                    these represent the same point.
!
!          DRATIO  : [real array, dimension NSEG]
!                    The ratio of the mesh spacing at the end of a
!                    segment to that at the beginning.
!
!          NFILT   : [integer]
!                    The number of times to filter the mesh-point
!                    distribution array.  Set NFILT=0 if filtering
!                    is not desired.  Filtering can reduce
!                    discontinuities in the derivative of the mesh
!                    spacing.
!
!          PERIODIC: [logical]
!                    A flag to indicate whether the mesh to be
!                    generated represents a periodic coordinate.
!                    If the coordinate is specified as periodic,
!                    the range [C0,C1] should be the whole periodic
!                    interval; the first mesh point is set at C0
!                    and the last mesh point, C(NC), is set at C1.
!
!          C_SHIFT : [real]
!                    Amount by which to shift the periodic coordinate.
!                    C_SHIFT is only used when PERIODIC=.true.,
!                    and is ignored otherwise.  A positive C_SHIFT
!                    moves the mesh points to the right.
!
! ****** Output arguments:
!
!          C       : [real array, dimension NC]
!                    The locations of the mesh points.
!
!-----------------------------------------------------------------------
!
! ****** The arrays DRATIO and FRAC define the mesh as follows.
!
! ****** For example, suppose that a (non-periodic) mesh with three
! ****** segments is desired.  Suppose the domain size is c=[0:2].
! ****** In the first segment (with c between 0 and .5) the mesh
! ****** spacing is decreasing with c, such that DC at c=.5 is half
! ****** DC at c=0.  From c=.5 to c=1, the mesh is uniform.  From c=1
! ****** to c=2, the mesh spacing is increasing with c such that DC at
! ****** c=2 is 10 times DC at c=1.  This mesh would be specified by:
! ******
! ******     FRAC=0.,.25,.5,1.
! ******     DRATIO=.5,1.,10.
! ******
! ****** The variable C_SHIFT can be used to shift the mesh point
! ****** distribution for a periodic coordinate.  For example,
! ****** suppose C represents mesh points in the interval [0,2*pi].
! ****** C_SHIFT=.5*pi would move the distribution of mesh points
! ****** so that the original mesh point with C(1)=0. would be
! ****** close to .5*pi in the new mesh.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
      use debug
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
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
!
!-----------------------------------------------------------------------
!
! ****** Storage for the coordinate transformation.
!
      integer :: ns
      real(r_typ), dimension(:), allocatable :: xi,cs,a,r,f
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: zero=0._r_typ
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: half=.5_r_typ
      real(r_typ), parameter :: eps=1.e-5_r_typ
!
!-----------------------------------------------------------------------
!
      integer :: i,j,nf,nr,ll,j0
      real(r_typ) :: alpha,dr,fac,d,dxi,xiv,cshft,xi_shift
      real(r_typ), dimension(:), allocatable :: dc,rdc
!
!-----------------------------------------------------------------------
!
! ****** Check that the number of mesh points is valid.
!
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
!
! ****** Check that a positive mesh interval has been specified.
!
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
!
! ****** Find the number of values of FRAC specified.
!
      do nf=1,nseg-1
        if (frac(nf+1).eq.zero) exit
      enddo
!
! ****** When no values have been specified (NF=1, the default),
! ****** a uniform mesh is produced.
!
      if (nf.eq.1.and.frac(1).eq.zero) then
        ns=1
        allocate (cs(ns+1))
        allocate (r(ns))
        cs(1)=c0
        cs(2)=c1
        r(1)=one
        go to 100
      end if
!
! ****** Check that the specified values of FRAC are monotonically
! ****** increasing.
!
      do i=2,nf
        if (frac(i).lt.frac(i-1)) then
          if (iamp0) then
            write (*,*)
            write (*,*) '### ERROR in GENMESH:'
            write (*,*) '### Invalid mesh specification.'
            write (*,*) 'Mesh coordinate: ',label
            write (*,*) 'The values in FRAC must increase'// &
                       ' monotonically.'
            write (*,*) 'FRAC = ',frac(1:nf)
          end if
          call endrun (.true.)
        end if
      enddo
!
! ****** Check the specified values of FRAC.
!
      if (periodic) then
!
! ****** A periodic mesh requires the specified values of FRAC
! ****** to be in the range 0. to 1.
!
        if (frac(1).lt.zero.or.frac(nf).gt.one) then
          if (iamp0) then
            write (*,*)
            write (*,*) '### ERROR in GENMESH:'
            write (*,*) '### Invalid mesh specification.'
            write (*,*) 'Mesh coordinate: ',label
            write (*,*) 'For a periodic coordinate, the values in'// &
                       ' FRAC must be between 0. and 1.'
            write (*,*) 'FRAC = ',frac(1:nf)
          end if
          call endrun (.true.)
        end if
!
! ****** A periodic mesh cannot contain both 0. and 1. in FRAC,
! ****** since these represent the same point.
!
        if (frac(1).eq.zero.and.frac(nf).eq.one) then
          if (iamp0) then
            write (*,*)
            write (*,*) '### ERROR in GENMESH:'
            write (*,*) '### Invalid mesh specification.'
            write (*,*) 'Mesh coordinate: ',label
            write (*,*) 'For a periodic coordinate, FRAC must not'// &
                        ' contain both 0. and 1.'
            write (*,*) 'FRAC = ',frac(1:nf)
          end if
          call endrun (.true.)
        end if
!
      else
!
! ****** A non-periodic mesh requires the first specified value
! ****** of FRAC to be 0., and the last value to equal 1.
!
        if (frac(1).ne.zero) then
          if (iamp0) then
            write (*,*)
            write (*,*) '### ERROR in GENMESH:'
            write (*,*) '### Invalid mesh specification.'
            write (*,*) 'Mesh coordinate: ',label
            write (*,*) 'For a non-periodic coordinate, the first'// &
                      ' value of FRAC must equal 0.'
            write (*,*) 'FRAC = ',frac(1:nf)
          end if
          call endrun (.true.)
        end if
!
        if (frac(nf).ne.one) then
          if (iamp0) then
            write (*,*)
            write (*,*) '### ERROR in GENMESH:'
            write (*,*) '### Invalid mesh specification.'
            write (*,*) 'Mesh coordinate: ',label
            write (*,*) 'For a non-periodic coordinate, the last'// &
                        ' value of FRAC must equal 1.'
            write (*,*) 'FRAC = ',frac(1:nf)
          end if
          call endrun (.true.)
        end if
!
      end if
!
! ****** Check that the required values of DRATIO have been set,
! ****** and are positive.
!
      if (periodic) then
        nr=nf
      else
        nr=nf-1
      end if
!
      do i=1,nr
        if (dratio(i).le.zero) then
          if (iamp0) then
            write (*,*)
            write (*,*) '### ERROR in GENMESH:'
            write (*,*) '### Invalid mesh specification.'
            write (*,*) 'Mesh coordinate: ',label
            write (*,*) 'A required value in DRATIO has not been'// &
                        ' set or is not positive.'
            write (*,*) 'DRATIO = ',dratio(1:nr)
          end if
          call endrun (.true.)
        end if
      enddo
!
! ****** Check that an inherently discontinuous mesh has not been
! ****** specified inadvertently.
!
      if (periodic.and.nr.eq.1.and.dratio(1).ne.one) then
        if (iamp0) then
          write (*,*)
          write (*,*) '### WARNING from GENMESH:'
          write (*,*) '### Discontinuous mesh specification.'
          write (*,*) 'Mesh coordinate: ',label
          write (*,*) 'An inherently discontinuous mesh has been'// &
                      ' specified.'
          write (*,*) 'FRAC = ',frac(1:nf)
          write (*,*) 'DRATIO = ',dratio(1:nr)
        end if
      end if
!
! ****** Set the number of segments.
!
      ns=nf-1
!
! ****** For a periodic coordinate, add points at XI=0. and XI=1.
! ****** if they are not already present.
!
      if (periodic) then
        if (frac(1).ne.zero) ns=ns+1
        if (frac(nf).ne.one) ns=ns+1
      end if
!
      allocate (cs(ns+1))
      allocate (r(ns))
!
! ****** Set up the coordinate limits of the segments.
!
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
!
  100 continue
!
      allocate (xi(ns+1))
      allocate (a(ns))
      allocate (f(ns))
!
! ****** Compute the XI values at the segment limits.
!
      do i=1,ns
        dr=r(i)-one
        if (abs(dr).lt.eps) then
          f(i)=(cs(i+1)-cs(i))*(one+half*dr)
        else
          f(i)=(cs(i+1)-cs(i))*log(r(i))/dr
        end if
      enddo
!
      fac=zero
      do i=ns,1,-1
        fac=fac/r(i)+f(i)
      enddo
!
      d=f(1)/fac
      xi(1)=zero
      do i=2,ns
        xi(i)=xi(i-1)+d
        if (i.lt.ns) d=d*f(i)/(f(i-1)*r(i-1))
      enddo
      xi(ns+1)=one
!
! ****** Set the amplification factor for each segment.
!
      do i=1,ns
        a(i)=log(r(i))/(xi(i+1)-xi(i))
      enddo
!
! ****** For a periodic coordinate, find the XI shift corresponding
! ****** to a shift of C_SHIFT in the coordinate.
! ****** Note that a positive value of C_SHIFT moves the mesh
! ****** points to the right.
!
      if (periodic) then
        cshft=-c_shift
        call map_c_to_xi (periodic,ns,xi,cs,a,r,cshft,xi_shift)
      else
        xi_shift=0.
      end if
!
! ****** Compute the location of the mesh points in array C
! ****** by mapping from the XI values.
!
      dxi=one/(nc-one)
!
      c(1)=c0
      do j=2,nc-1
        xiv=(j-1)*dxi
        call map_xi_to_c (periodic,ns,xi,cs,a,r,cshft,xi_shift,xiv,c(j))
      enddo
      c(nc)=c1
!
! ****** Filter the mesh if requested.
!
      if (nfilt.gt.0) then
        do i=1,nfilt
          if (periodic) then
            call filter_coord_periodic (c1-c0,nc,c)
          else
            call filter_coord (nc,c)
          end if
        enddo
      end if
!
! ****** Write out the mesh information.
!
      if (io.gt.0.and.iamp0) then
!
        write (io,*)
        write (io,*) '### COMMENT from GENMESH:'
        write (io,*) '### Mesh information for coordinate ',label,':'
!
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
!
        write (io,*)
        write (io,*) 'Number of mesh segments = ',ns
!
        ll=len_trim(label)
!
        write (io,900) 'Segment      xi0       xi1'//   &
                       repeat (' ',16-ll)//label//'0'// &
                       repeat (' ',16-ll)//label//'1'// &
                       '            ratio'
        do i=1,ns
          write (io,910) i,xi(i),xi(i+1),cs(i),cs(i+1),r(i)
        enddo
!
        allocate (dc(nc))
        allocate (rdc(nc))
!
        dc=c-cshift(c,-1)
        if (periodic) dc(1)=dc(nc)
        rdc=dc/cshift(dc,-1)
        if (periodic) rdc(1)=rdc(nc)
!
        write (io,*)
        write (io,920) 'Mesh-point locations:'
        write (io,920) '     i'//                         &
                       repeat (' ',18-ll)//label//        &
                       repeat (' ',17-ll)//'d'//label//   &
                       '             ratio'
!
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
!
        deallocate (dc)
        deallocate (rdc)
!
      end if
!
  900 format (/,tr1,a)
  910 format (tr1,i4,2x,2f10.6,4f17.8)
  920 format (tr1,a)
  930 format (tr1,i6,3f18.8)
!
      deallocate (cs)
      deallocate (r)
      deallocate (xi)
      deallocate (a)
      deallocate (f)
!
end subroutine
!#######################################################################
subroutine map_xi_to_c (periodic,ns,xi,cs,a,r,cshft,xi_shift,xiv,cv)
!
!-----------------------------------------------------------------------
!
! ****** Get the mesh coordinate value CV for the specified
! ****** xi value XIV.
!
! ****** Set PERIODIC=.true. for a periodic coordinate.
! ****** NS is the number of segments in the mesh definition.
! ****** The arrays XI, CS, A, and R define the mesh mapping.
!
! ****** CSHFT represents the amount to shift a periodic coordinate.
! ****** XI_SHIFT represents the corresponding amount to shift xi.
!
! ****** This is a utility routine for GENMESH.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      logical, intent(in) :: periodic
      integer, intent(in) :: ns
      real(r_typ), dimension(ns+1), intent(in) :: xi,cs
      real(r_typ), dimension(ns), intent(in) :: a,r
      real(r_typ), intent(in) :: cshft,xi_shift
      real(r_typ), intent(in) :: xiv
      real(r_typ), intent(out) :: cv
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: zero=0._r_typ
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: half=.5_r_typ
      real(r_typ), parameter :: eps=1.e-5_r_typ
!
!-----------------------------------------------------------------------
!
      integer :: i
      real(r_typ) :: xiv_p,d,d1,da,da1,fac
!
!-----------------------------------------------------------------------
!
      real(r_typ), external :: fold
!
!-----------------------------------------------------------------------
!
! ****** Find the index of the segment to which XIV belongs.
!
      if (periodic) then
!
! ****** Shift XIV by XI_SHIFT.
!
        xiv_p=xiv+xi_shift
!
! ****** Fold XIV_P into the main interval.
!
        xiv_p=fold(zero,one,xiv_p)
!
      else
!
        xiv_p=xiv
!
      end if
!
      do i=1,ns
        if (xiv_p.ge.xi(i).and.xiv_p.le.xi(i+1)) exit
      enddo
!
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
!
      d =xiv_p  -xi(i)
      d1=xi(i+1)-xi(i)
!
      da =a(i)*d
      da1=a(i)*d1
!
! ****** Interpolate the mapping function at XIV_P.
!
      if (abs(da1).lt.eps) then
        fac=(d*(one+half*da))/(d1*(one+half*da1))
      else
        fac=(exp(da)-one)/(r(i)-one)
      end if
!
      cv=cs(i)+(cs(i+1)-cs(i))*fac
!
      if (periodic) then
!
! ****** Shift CV by the amount CSHFT.
!
        cv=cv-cshft
!
! ****** Fold CV into the main interval.
!
        cv=fold(cs(1),cs(ns+1),cv)
!
      end if
!
end subroutine
!#######################################################################
subroutine map_c_to_xi (periodic,ns,xi,cs,a,r,cv,xiv)
!
!-----------------------------------------------------------------------
!
! ****** Get the xi value XIV for the specified coordinate value CV.
!
! ****** Set PERIODIC=.true. for a periodic coordinate.
! ****** NS is the number of segments in the mesh definition.
! ****** The arrays XI, CS, A, and R define the mesh mapping.
!
! ****** This is a utility routine for GENMESH.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      logical, intent(in) :: periodic
      integer, intent(in) :: ns
      real(r_typ), dimension(ns+1), intent(in) :: xi,cs
      real(r_typ), dimension(ns), intent(in) :: a,r
      real(r_typ), intent(in) :: cv
      real(r_typ), intent(out) :: xiv
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
      real(r_typ), parameter :: eps=1.e-5_r_typ
!
!-----------------------------------------------------------------------
!
      integer :: i
      real(r_typ) :: cv_p,d,da,fac
!
!-----------------------------------------------------------------------
!
      real(r_typ), external :: fold
!
!-----------------------------------------------------------------------
!
! ****** Find the index of the segment to which CV belongs.
!
      if (periodic) then
!
! ****** Fold CV_P into the main interval.
!
        cv_p=fold(cs(1),cs(ns+1),cv)
!
      else
!
        cv_p=cv
!
      end if
!
      do i=1,ns
        if (cv_p.ge.cs(i).and.cv_p.le.cs(i+1)) exit
      enddo
!
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
!
      d=(cv_p-cs(i))/(cs(i+1)-cs(i))
      da=(r(i)-one)*d
!
! ****** Interpolate the mapping function at XIV_P.
!
      if (abs(da).lt.eps) then
        fac=d*(xi(i+1)-xi(i))
      else
        fac=log(da+one)/a(i)
      end if
!
      xiv=xi(i)+fac
!
end subroutine
!#######################################################################
subroutine filter_coord (n,f)
!
!-----------------------------------------------------------------------
!
! ****** Apply a "(1,2,1)/4" low-pass digital filter to a
! ****** 1D coordinate.
!
! ****** The end-points F(1) and F(N) are not changed.
!
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: n
      real(r_typ), dimension(n) :: f
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: two=2._r_typ
      real(r_typ), parameter :: quarter=.25_r_typ
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(n) :: ff
!
!-----------------------------------------------------------------------
!
      integer :: i
!
!-----------------------------------------------------------------------
!
! ****** Make a copy of the function.
!
      ff=f
!
! ****** Apply the filter.
!
      do i=2,n-1
        f(i)=quarter*(ff(i-1)+two*ff(i)+ff(i+1))
      enddo
!
end subroutine
!#######################################################################
subroutine filter_coord_periodic (xl,n,f)
!
!-----------------------------------------------------------------------
!
! ****** Apply a "(1,2,1)/4" low-pass digital filter to a
! ****** periodic 1D coordinate.
!
!-----------------------------------------------------------------------
!
! ****** XL is the periodic interval for the coordinate.
!
! ****** The filtered coordinate is translated so that F(1)
! ****** is preserved.
!
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: xl
      integer :: n
      real(r_typ), dimension(n) :: f
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: two=2._r_typ
      real(r_typ), parameter :: quarter=.25_r_typ
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(0:n+1) :: ff
!
!-----------------------------------------------------------------------
!
      integer :: i
      real(r_typ) :: f1old,f1new
!
!-----------------------------------------------------------------------
!
! ****** Save the value of F(1).
!
      f1old=f(1)
!
! ****** Make a periodic copy of the function.
!
      ff(1:n)=f(:)
!
      ff(0)=f(n-1)-xl
      ff(n+1)=f(2)+xl
!
! ****** Apply the filter.
!
      do i=1,n
        f(i)=quarter*(ff(i-1)+two*ff(i)+ff(i+1))
      enddo
!
! ****** Translate F so that F(1) is preserved.
!
      f1new=f(1)
      do i=1,n
        f(i)=f(i)-f1new+f1old
      enddo
!
end subroutine
!#######################################################################
function fold (x0,x1,x)
!
!-----------------------------------------------------------------------
!
! ****** "Fold" X into the periodic interval [X0,X1].
!
! ****** On return, X is such that X0.le.X.lt.X1.
!
!-----------------------------------------------------------------------
!
! ****** It is assumed that X0 does not equal X1, as is physically
! ****** necessary.  If X0 and X1 are equal, the routine just
! ****** returns with FOLD=X.
!
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: fold
      real(r_typ) :: x0,x1,x
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: xl
!
!-----------------------------------------------------------------------
!
      fold=x
!
      if (x0.eq.x1) return
!
      xl=x1-x0
!
      fold=mod(x-x0,xl)+x0
!
      if (fold.lt.x0) fold=fold+xl
      if (fold.ge.x1) fold=fold-xl
!
      return
end function
!#######################################################################
subroutine set_flux
!
!-----------------------------------------------------------------------
!
! ****** Set the radial magnetic field at the photosphere.
!
!-----------------------------------------------------------------------
!
      use number_types
      use global_dims
      use global_mesh
      use local_dims
      use local_mesh
      use fields
      use vars
      use mpidefs
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
! ****** Global flux array.
!
      real(r_typ), dimension(:,:), allocatable :: br0_g
!
!-----------------------------------------------------------------------
!
      integer :: j,k,ierr
!
!-----------------------------------------------------------------------
!
      call timer_on
!
      allocate (br0_g(nt_g,np_g))
!
! ****** Define the global flux array.
!
! ****** Read the flux from file BR0FILE (only on processor IPROC0).
!
      if (iamp0) then
        call readbr (br0file,br0_g,ierr)
      end if
      call check_error_on_p0 (ierr)
!
! ****** Broadcast BR0_G to all the processors.
!
      call MPI_Bcast (br0_g,nt_g*np_g,ntype_real,0,comm_all,ierr)
!
! ****** For a fully open field, reverse negative Br
! ****** (i.e., use the monopole trick).
!
      if (option.eq.'open') then
!
! ****** Write the boundary flux (before the sign flip) to a file
! ****** if requested.
!
        if (iamp0) then
          if (br_photo_original_file.ne.'') then
            write (*,*)
            write (*,*) '### COMMENT from SET_FLUX:'
            write (*,*)
            write (*,*) 'Writing BR0 (before sign flip) to file: ', &
                        trim(br_photo_original_file)
            write (9,*)
            write (9,*) '### COMMENT from SET_FLUX:'
            write (9,*)
            write (9,*) 'Writing BR0 (before sign flip) to file: ', &
                        trim(br_photo_original_file)
            ! [XX: HDF5 file write]
            ! call wrhdf_2d (br_photo_original_file,.true.,nt_g,np_g, &
            !                br0_g,th_g,ph_g,hdf32,ierr)
          end if
        end if
!
! ****** Reverse Br.
!
        br0_g(:,:)=abs(br0_g(:,:))
!
      end if
!
! ****** Write the boundary flux to a file if requested.
!
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
      ! [XX: HDF5 file write]
      !     call wrhdf_2d (br_photo_file,.true.,nt_g,np_g, &
      !                    br0_g,th_g,ph_g,hdf32,ierr)
        end if
      end if
!
      do j=1,nt
        do k=1,np
          br0(j,k)=br0_g(j0_g+j-1,k0_g+k-1)
        enddo
      enddo
!$omp target update to(br0)
!
      deallocate(br0_g)
!
      call timer_off (t_io)
!
end subroutine
!#######################################################################
subroutine potfld
!
!-----------------------------------------------------------------------
!
! ****** Find the potential field solution.
!
!-----------------------------------------------------------------------
!
      use number_types
      use local_dims
      use local_mesh
      use fields
      use cgcom
      use solve_params
      use mpidefs
      use debug
      use timing
      use matrix_storage_pot3d_solve
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
!
!-----------------------------------------------------------------------
!
      integer :: ierr,nrm2,ntm2,npm2,i
!
      real(r_typ), dimension(:), allocatable :: rhs_cg,x_cg
!
!-----------------------------------------------------------------------
!
! ****** Load matrix and preconditioner.
!
      nrm2=nrm1-1
      ntm2=ntm1-1
      npm2=npm1-1
!
      a_offsets(1)=-nrm2*ntm2
      a_offsets(2)=-nrm2
      a_offsets(3)=-1
      a_offsets(4)= 0
      a_offsets(5)= 1
      a_offsets(6)= nrm2
      a_offsets(7)= nrm2*ntm2
!
! ****** Allocate cg 1D vectors.
!
      N=nrm2*ntm2*npm2
!
! ****** Prepare the guess, and rhs for the solve.
!
      allocate(rhs_cg(N))
      allocate(x_cg(N))
      rhs_cg(:)=0.
      x_cg(:)=0.
!$omp target enter data map(to:rhs_cg,x_cg)
!
      call getM (N,a_offsets,M)
      call alloc_pot3d_matrix_coefs
      call load_matrix_pot3d_solve
!$omp target enter data map(to:a)
      call load_preconditioner_pot3d_solve
!$omp target enter data map(to:a_i)
!
! ****** Use a trick to accumulate the contribution of the
! ****** boundary conditions (i.e., the inhomogenous part).
!
      call set_boundary_points (x_ax,one)
      call seam_hhh (x_ax)
      call delsq (x_ax,rhs_cg)
!
! ****** Original rhs is zero so just use negative of boundary
!        trick contributions:
!
      do concurrent (i=1:N)
        rhs_cg(i)=-rhs_cg(i)
      enddo
!
! ****** Solve for the potential.
!
      if (idebug.gt.0.and.iamp0) then
        write (*,*)
        write (*,*) '### COMMENT from POTFLD:'
        write (*,*) '### Doing a solution:'
      end if
!
      call solve (x_cg,rhs_cg,N,ierr)
!
      if (ierr.ne.0) then
        call endrun (.true.)
      end if
!
      call unpack_scalar (phi,x_cg)
!
      call set_boundary_points (phi,one)
      call seam_hhh (phi)
!
!$omp target exit data map(delete:rhs_cg,x_cg,a,a_i)
      call dealloc_pot3d_matrix_coefs
      deallocate(rhs_cg)
      deallocate(x_cg)
!
end subroutine
!#######################################################################
subroutine set_validation_flux
!
!-----------------------------------------------------------------------
!
! ****** Set the radial magnetic field at the photosphere.
! ****** This uses the tilted dipole analytic solution for validation.
!
!-----------------------------------------------------------------------
!
      use number_types
      use global_dims
      use global_mesh
      use local_dims
      use local_mesh
      use fields
      use vars
      use mpidefs
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: two=2.0_r_typ
!
!-----------------------------------------------------------------------
!
! ****** Global Br boundary arrays.
!
      real(r_typ), dimension(:,:), allocatable :: br0_g
      real(r_typ), dimension(:,:), allocatable :: br1_g
!
!-----------------------------------------------------------------------
!
      integer :: j,k
      real(r_typ) :: sum00,sum11
!
!-----------------------------------------------------------------------
!
! ****** Define the global boundary condition arrays.
!
      allocate (br0_g(nt_g,np_g))
      allocate (br1_g(nt_g,np_g))
!
! ****** Set the tilted dipole analytic Br.
!
      do j=1,nt_g
        do k=1,np_g
          br0_g(j,k)=(two/r0**3)*(cth_g(j)*cos(dipole_angle) + &
                                  sth_g(j)*cph_g(k)*sin(dipole_angle))
          br1_g(j,k)=(two/r1**3)*(cth_g(j)*cos(dipole_angle) + &
                                  sth_g(j)*cph_g(k)*sin(dipole_angle))
        enddo
      enddo
!
! ****** Set Br to be periodic (should not be nesessary?).
!
      br0_g(:,1)=br0_g(:,npm1_g)
      br0_g(:,np_g)=br0_g(:,2)
      br1_g(:,1)=br1_g(:,npm1_g)
      br1_g(:,np_g)=br1_g(:,2)
!
! ****** Set BCs at the poles.
!
      sum00=sum(br0_g(     2,2:npm1_g)*dph_g(2:npm1_g))*pl_i
      sum11=sum(br0_g(ntm1_g,2:npm1_g)*dph_g(2:npm1_g))*pl_i
!
      br0_g(1   ,:)=two*sum00-br0_g(     2,:)
      br0_g(nt_g,:)=two*sum11-br0_g(ntm1_g,:)
!
      sum00=sum(br1_g(     2,2:npm1_g)*dph_g(2:npm1_g))*pl_i
      sum11=sum(br1_g(ntm1_g,2:npm1_g)*dph_g(2:npm1_g))*pl_i
!
      br1_g(1   ,:)=two*sum00-br1_g(     2,:)
      br1_g(nt_g,:)=two*sum11-br1_g(ntm1_g,:)
!
      do j=1,nt
        do k=1,np
          br0(j,k)=br0_g(j0_g+j-1,k0_g+k-1)
          br1(j,k)=br1_g(j0_g+j-1,k0_g+k-1)
        enddo
      enddo
!$omp target update to(br0,br1)
!
      deallocate(br0_g)
      deallocate(br1_g)
!
end subroutine
!#######################################################################
subroutine write_validation_solution
!
!-----------------------------------------------------------------------
!
! ****** Write out the tilted dipole analytic solution for validation.
! ****** This writes out the phi, B, and B computed from phi.
!
!-----------------------------------------------------------------------
!
      use number_types
      use global_dims
      use global_mesh
      use local_dims
      use local_mesh
      use fields
      use vars
      use mpidefs
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1.0_r_typ
      real(r_typ), parameter :: two=2.0_r_typ
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k
      character(256) :: fnamephi,fnamebr,fnamebt,fnamebp
!
!-----------------------------------------------------------------------
!
! ****** Set phi.
!
      do concurrent (k=1:np, j=1:nt, i=1:nr)
        phi(i,j,k)=(-one/rh(i)**2)*(cth(j)*cos(dipole_angle)+ &
                                   sth(j)*cph(k)*sin(dipole_angle))
      enddo
!
      call set_boundary_points (phi,one)
      call seam_hhh (phi)
!
!$omp target enter data map(alloc:br,bt,bp)
!
! ****** Set Br.
!
      do concurrent (k=1:np, j=1:nt, i=1:nrm1)
        br(i,j,k)=(two/r(i)**3)*(cth(j)*cos(dipole_angle) + &
                                 sth(j)*cph(k)*sin(dipole_angle))
      enddo
!
! ****** Set Bt.
!
      do concurrent (k=1:np, j=1:ntm1, i=1:nr)
        bt(i,j,k)=(one/rh(i)**3)*(st(j)*cos(dipole_angle) - &
                                  ct(j)*cph(k)*sin(dipole_angle))
      enddo
!
! ****** Set Bp.
!
      do concurrent (k=1:npm1, j=1:nt, i=1:nr)
        bp(i,j,k)=(one/rh(i)**3)*(sph(k)*sin(dipole_angle))
      enddo
!
! ****** Due to half mesh trig issues on periodic domain, we need
! ****** to seam these even though seaming is not needed in getb().
!
      call seam_gen (br,nrm,nt,np)
      call seam_gen (bt,nr,ntm,np)
      call seam_gen (bp,nr,nt,npm)
!
! ****** Backup fnames and set new ones so we can reuse write routine.
!
      fnamephi=phifile
      phifile='phi_exact.'//trim(fmt)
      fnamebr=brfile
      brfile='br_exact.'//trim(fmt)
      fnamebt=btfile
      btfile='bt_exact.'//trim(fmt)
      fnamebp=bpfile
      bpfile='bp_exact.'//trim(fmt)
!
      call write_solution
!
! ****** Restore original output names.
!
      phifile=fnamephi
      brfile=fnamebr
      btfile=fnamebt
      bpfile=fnamebp
!
! ****** Reset phi.
!
      phi(:,:,:)=0.
!$omp target update to(phi)
!$omp target exit data map(delete:br,bt,bp)
!
end subroutine
!#######################################################################
subroutine solve (x,rhs,N,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Solve the implicit equations iteratively.
!
!-----------------------------------------------------------------------
!
! ****** Return IERR=0 if the iteration converges; otherwise,
! ****** IERR is set to a nonzero value.
!
! ****** X is the initial guess at the solution.
! ****** RHS is the right-hand side.
!
!-----------------------------------------------------------------------
!
      use number_types
      use cgcom
      use mpidefs
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: N
      real(r_typ), dimension(N) :: x,rhs
      integer :: ierr
!
!-----------------------------------------------------------------------
!
! ****** Solve the equations using the CG method.
!
      call cgsolve (x,rhs,N,ierr)
!
! ****** Check for convergence.
!
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
!
end subroutine
!#######################################################################
subroutine cgsolve (x,r,N,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Solve the linear system:
!
!            A * x = b
!
! ****** using the classical Conjugate Gradient method for symmetric
! ****** and positive-definite matrices.
!
!-----------------------------------------------------------------------
!
! ****** On input, X(N) contains a guess at the solution, and
! ****** R(N) contains the right-hand side, b.
!
! ****** On exit, X contains an estimate to the solution, and
! ****** R contains the residual (b-Ax).
!
! ****** IERR=0 indicates that the solution converged to the
! ****** requested accuracy.  Other values indicate that the
! ****** iteration did not converge for the given maximum number
! ****** of iterations.
!
!-----------------------------------------------------------------------
!
      use number_types
      use cgcom
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: N
      real(r_typ), dimension(N) :: x,r
      integer :: ierr,i
!
!-----------------------------------------------------------------------
!
! ****** Scratch space for the CG iteration vectors.
!
      real(r_typ), dimension(N), target :: p,ap
!
!-----------------------------------------------------------------------
!
      real(r_typ), external :: cgdot
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: bdotb,rdotr,pdotap,alphai,rdotr_old,betai
!
!-----------------------------------------------------------------------
!
      ncg=0
!$omp target enter data map(alloc:p,ap)
!
! ****** Get the norm of the RHS.
!
      do concurrent (i=1:N)
        p(i)=r(i)
      enddo
!
      call prec_inv (p)
      bdotb=cgdot(r,p,N)
!
! ****** If the RHS is zero, return with a zero solution.
!
      if (bdotb.eq.0.) then
        do concurrent (i=1:N)
          x(i)=0.
        enddo
        epsn=0.
        ierr=0
        return
      end if
!
!-----------------------------------------------------------------------
! ****** Initialization.
!-----------------------------------------------------------------------
!
      call ax (x,ap,N)
!
      do concurrent (i=1:N)
        r(i)=r(i)-ap(i)
        p(i)=r(i)
      enddo
!
! ****** Find the initial error norm.
!
      call prec_inv (p)
      rdotr=cgdot(r,p,N)
!
      call ernorm (bdotb,rdotr,ierr)
      if (ierr.ge.0) return
!
!-----------------------------------------------------------------------
! ****** Main iteration loop.
!-----------------------------------------------------------------------
!
      do
        ncg=ncg+1
!
        call ax (p,ap,N)
!
        pdotap=cgdot(p,ap,N)
        alphai=rdotr/pdotap
!
        do concurrent (i=1:N)
          x(i)=x(i)+alphai*p(i)
          r(i)=r(i)-alphai*ap(i)
          ap(i)=r(i)
        enddo
!
        call prec_inv (ap)
        rdotr_old=rdotr
        rdotr=cgdot(r,ap,N)
!
! ****** Check for convergence.
!
        call ernorm (bdotb,rdotr,ierr)
        if (ierr.ge.0) exit
!
        betai=rdotr/rdotr_old
!
        do concurrent (i=1:N)
          p(i)=betai*p(i)+ap(i)
        enddo
!
      enddo
!
!$omp target exit data map(delete:p,ap)
end subroutine
!#######################################################################
subroutine ernorm (bdotb,rdotr,ierr)
!
!-----------------------------------------------------------------------
!
! ****** This subroutine checks if the iterative solver has
! ****** converged or if the maximum allowed number of iterations,
! ****** NCGMAX, has been exceeded.
!
!-----------------------------------------------------------------------
!
! ****** Convergence is deemed to have occurred when:
! ******
! ******     ||R||/||B|| .lt. EPSCG
! ******
! ****** where ||R|| is the norm of the (preconditioned)
! ****** residual, ||B|| is the norm of the (preconditioned)
! ****** RHS, and EPSCG is the specified convergence criterion.
!
! ****** Set IERR=0 if the error is below the error criterion
! ****** (i.e., the solution has converged).
! ****** Set IERR=-1 if the error does not yet meet the error
! ****** criterion and the number of iterations is less than NCGMAX.
! ****** Set IERR=1 if the maximum number of iterations has
! ****** been exceeded without convergence.
!
!-----------------------------------------------------------------------
!
! ****** On input, BDOTB has the dot product of the RHS vector
! ****** with itself, weighted by the preconditioning matrix.
! ****** Similarly, RDOTR has the dot product of the residual vector
! ****** with itself, weighted by the preconditioning matrix.
! ****** This is used to normalize the error estimate.
!
!-----------------------------------------------------------------------
!
      use number_types
      use cgcom
      use mpidefs
      use vars
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: bdotb,rdotr
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: epssq
!
!-----------------------------------------------------------------------
!
      ierr=-1
!
      epssq=rdotr/bdotb
      epsn=sign(one,epssq)*sqrt(abs(epssq))
!
      if (ncghist.gt.0) then
!
        if (mod(ncg,ncghist).eq.0) then
          if (iamp0) then
            if (ncg.eq.0) then
              write (9,*)
              write (9,*) '### Convergence information:'
              write (*,*)
            end if
            write (*,100) ncg,epsn
            FLUSH(OUTPUT_UNIT)
            write (9,100) ncg,epsn
            FLUSH(9)
  100       format (1x,'Iteration: ',i8,'   Residual: ',1pe23.15)
          end if
        end if
      end if
!
! ****** Check for convergence.
!
      if (epsn.lt.epscg) then
        if (ncghist.gt.0) then
          if (iamp0) then
            write (*,100) ncg,epsn
            FLUSH(OUTPUT_UNIT)
            write (9,*)
            write (9,*) '### The CG solver has converged.'
            write (9,100) ncg,epsn
            FLUSH(9)
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
          FLUSH(OUTPUT_UNIT)
          write (9,*)
          write (9,*) '### ERROR in ERNORM:'
          write (9,*) '### Exceeded maximum number of iterations.'
          write (9,*) 'NCGMAX = ',ncgmax
          write (9,*) 'EPSN = ',epsn
          FLUSH(9)
        end if
        ierr=1
      end if
!
end subroutine
!#######################################################################
subroutine alloc_pot3d_matrix_coefs
!
!-----------------------------------------------------------------------
!
! ****** Allocate the arrays in which the matrix coefficients
! ****** for the pot3d solve are stored.
!
!-----------------------------------------------------------------------
!
      use matrix_storage_pot3d_solve
      use cgcom
      use local_dims
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      allocate (a(2:nrm1,2:ntm1,2:npm1,7))
      a=0.
      allocate (a_i(N))
      a_i=0.
!
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
!
end subroutine
!#######################################################################
subroutine dealloc_pot3d_matrix_coefs
!
!-----------------------------------------------------------------------
!
! ****** Deallocate the arrays in which the matrix coefficients
! ****** for the pot3d solve are stored.
!
!-----------------------------------------------------------------------
!
      use matrix_storage_pot3d_solve
      use cgcom
      use cusparse_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      deallocate (a)
      deallocate (a_i)
!$omp target exit data map(delete:a,a_i)
!
      if (ifprec.eq.2) then
#ifdef CUSPARSE
        call unload_lusol_cusparse
!$omp target exit data map(delete:a_csr,a_csr_ja,a_csr_ia)
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
!
end subroutine
!#######################################################################
subroutine load_matrix_pot3d_solve
!
!-----------------------------------------------------------------------
!
! ****** Load the matrix coefficients for the pot3d solve.
!
!-----------------------------------------------------------------------
!
      use number_types
      use matrix_storage_pot3d_solve
      use local_dims
      use local_mesh
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k
!
!-----------------------------------------------------------------------
!
! ****** Set matrix coefs
!
      do k=2,npm1
        do j=2,ntm1
          do i=2,nrm1
!           a*ps(i,j,k-1):
            a(i,j,k,1)=-drh(i)*dth(j)*sth_i(j)*dp_i(k-1)
!           a*ps(i,j-1,k):
            a(i,j,k,2)=-drh(i)*dph(k)*st(j-1)*dt_i(j-1)
!           a*ps(i-1,j,k):
            a(i,j,k,3)=-sth(j)*dth(j)*dph(k)*r2(i-1)*dr_i(i-1)
!           a*ps(i+1,j,k):
            a(i,j,k,5)=-sth(j)*dth(j)*dph(k)*r2(i  )*dr_i(i  )
!           a*ps(i,j+1,k):
            a(i,j,k,6)=-drh(i)*dph(k)*st(j  )*dt_i(j  )
!           a*ps(i,j,k+1):
            a(i,j,k,7)=-drh(i)*dth(j)*sth_i(j)*dp_i(k  )
!
!           a*ps(i,j,k):
            a(i,j,k,4)=-(a(i,j,k,1)+a(i,j,k,2)+a(i,j,k,3)+ &
                         a(i,j,k,5)+a(i,j,k,6)+a(i,j,k,7))
          enddo
        enddo
      enddo
!
end subroutine
!#######################################################################
subroutine load_preconditioner_pot3d_solve
!
!-----------------------------------------------------------------------
!
! ****** Load the preconditioner for the pot3d solve.
!
!-----------------------------------------------------------------------
!
      use number_types
      use matrix_storage_pot3d_solve
      use cgcom
      use local_dims
      use, intrinsic :: iso_c_binding
      use cusparse_interface
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: one=1._r_typ
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k,icode,ii
!
!-----------------------------------------------------------------------
!
      call timer_on
!
      if (ifprec.eq.0) return
!
      if (ifprec.eq.1) then
!
! ****** Diagonal scaling:
!
        ii=0
        do k=2,npm1
          do j=2,ntm1
            do i=2,nrm1
              ii=ii+1
              a_i(ii)=real(one/a(i,j,k,4),r_typ_pc)
            enddo
          enddo
        enddo
!
      elseif (ifprec.eq.2) then
!
! ****** Convert A matrix into CSR format:
!
        call diacsr (N,M,a,a_offsets,a_csr,a_csr_ja,a_csr_ia,a_csr_dptr)
#ifdef CUSPARSE
      cN=N
      cM=M
!$omp target enter data map(to:a_csr,a_csr_ja,a_csr_ia)
!$omp target data use_device_ptr(a_csr,a_csr_ja,a_csr_ia)
      call load_lusol_cusparse (C_LOC(a_csr(1)),          &
                                C_LOC(a_csr_ia(1)),       &
                                C_LOC(a_csr_ja(1)),cN,cM)
!$omp end target data
#else
!
! ****** Overwrite CSR A with preconditioner L and U matrices:
!
! ****** Incomplete LU (ILU)
!
        icode=0
        call ilu0 (N,M,a_csr,a_csr_ja,a_csr_ia,a_csr_dptr,icode)
!
        if (icode.ne.0) then
          print*, '### ERROR IN ILU FORMATION'
        endif
!
! ****** Convert LU stored in A to LU matrix in optimized layout.
!
        call lu2luopt (N,M,lu_csr,a_csr,a_csr_ia,a_csr_ja,lu_csr_ja, &
                       a_csr_dptr,a_N1,a_N2)
!
! ****** Store inverse of diagonal of LU matrix.
!
        do i=1,N
          a_csr_d(i)=one/a_csr(a_csr_dptr(i))
        enddo
!
#endif
      endif
!
      call timer_off (t_pc_load)
!
end subroutine
!#######################################################################
subroutine ilu0 (N,M,A,JA,IA,A_da,icode)
!
!-----------------------------------------------------------
!
!     Set-up routine for ILU(0) preconditioner. This routine
!     computes the L and U factors of the ILU(0) factorization
!     of a general sparse matrix A stored in CSR format with
!     1-based indexing. Since
!     L is unit triangular, the L and U factors can be stored
!     as a single matrix which occupies the same storage as A.
!     New ja and ia arrays are not needed for the LU matrix
!     since the pattern of the LU matrix is identical with
!     that of A.
!
!     Original Author:  Yousef Saad
!            Iterative Methods for Sparse Linear Systems 2nd Ed. pg. 309
!     Modified by R.M. Caplan
!
!-----------------------------------------------------------
!     INPUT:
!     N         : Dimension of matrix
!     A, JA, IA : Sparse matrix in CSR sparse storage format
!     A_da      : Pointers to the diagonal elements in the CSR
!                 data structure luval
!
!     OUTPUT:
!     A     : L/U matrices stored together. On return A,
!             JA, and IA are the combined CSR data structure for
!             the L and U factors.
!     icode : Integer indicating error code on return:
!             (0): Normal return.
!             (k): Encountered a zero pivot at step k.
!------------------------------------------------------------
!
      use number_types
      use number_types_pc
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: JA(M),IA(N+1),A_da(N),iw(N)
      integer :: icode,N,M
      real(r_typ_pc) :: A(M)
!
!-----------------------------------------------------------------------
!
      integer :: i,ik,kj,k,ij,IA_i,IA_ip1m1
      real(r_typ_pc) :: Aik
!
!-----------------------------------------------------------------------
!
      icode=0
!     Initialize scratch index array:
      iw(:)=0
!
      do i=2,N
!       Store index of (i,j) in A in scratch array of iw(j=1:N)
!       This allows lookup given a column index (j) in row (k)
!       to see if the column is in row (i).
        IA_i    =IA(i)
        IA_ip1m1=IA(i+1)-1
!
        do ij=IA_i,IA_ip1m1
          iw(JA(ij))=ij
        enddo
!
!       Loop from first element in row i to 1 less than diagonal elem:
        do ik=IA_i,A_da(i)-1     !IA(i+1) !ik is index of (i,k) in A[]
          k    =JA(ik)           !Actual column index in matrix (k)
          Aik  =A(ik)/A(A_da(k)) !Save Aik for next loop as an optim.
          A(ik)=Aik
!
!         Loop from 1 more than diag elem to last elem in row k:
          do kj=A_da(k)+1,IA(k+1)-1 !kj is index of (k,j) in A[]
!            Get ij location from scratch array (if 0, no ij present)
             ij=iw(JA(kj))
             if (ij .ne. 0) then
               A(ij)=A(ij)-Aik*A(kj)
             endif
          enddo
        enddo
!
        if (A(ik).eq.0) then
          icode=i
          exit
        endif
!
!       Reset scratch index array:
        do ij=IA_i,IA_ip1m1
          iw(JA(ij))=0
        enddo
      enddo
!
end subroutine
!#######################################################################
subroutine lu2luopt (N,M,LU,A,IA,JA,LUJA,A_da,N1,N2)
!
!-----------------------------------------------------------------------
!
! ****** Re-order elements of LU matrix in CSR format into custom,
! ****** optimized format for use with lusol().
! ****** (Eventually, this could be merged with the ilu0 and/or diacsr)
!
!-----------------------------------------------------------------------
!
      use number_types
      use number_types_pc
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: N,M
      integer :: JA(M),LUJA(M),IA(N+1),A_da(N)
      integer :: N1(N),N2(N)
      real(r_typ_pc) :: A(M),LU(M)
!
!-----------------------------------------------------------------------
!
      integer :: i,k,ii
!
!-----------------------------------------------------------------------
!
      ii=0
!
      do i=1,N
        do k=IA(i),A_da(i)-1
           ii=ii+1
           LU(ii)=A(k)
           LUJA(ii)=JA(k)
        enddo
!
!       Store k1 and k2 ranges for lusolve:
!
        N1(i)=A_da(i)-1-IA(i)
        N2(i)=IA(i+1)-2-A_da(i)
      enddo
!
      do i=N,1,-1
        do k=A_da(i)+1,IA(i+1)-1
           ii=ii+1
           LU(ii)=A(k)
           LUJA(ii)=JA(k)
        enddo
      enddo
!
end subroutine
!#######################################################################
subroutine diacsr (N,M,Adia,ioff,Acsr,JA,IA,Adptr)
!
!-----------------------------------------------------------------------
!
! *** DIACSR_INNER converts a solver matrix in a MAS-style
!     diagonal format to standard compressed sparse row (CSR)
!     including periodic coefficents when nproc_p=1.
!
!     Author of original diacsr: Youcef Saad
!     Modifications for MAS:     RM Caplan
!
!     Input:
!                     N: Size of the matrix (NxN)
!                     M: Number of non-zero entries in matrix
!                        (computed with getM_tc())
!         Adia(IDIAG,N): The matrix in modified "DIA" format
!           ioff(IDIAG): Offsets of the diagonals in A.
!
!     Output:
!            Acsr(M), JA(M), IA(N+1): The matrix A in CSR.
!                           Adptr(N): Pointers to diag elements in A,
!                                     [e.g. A(i,i) == A(Adptr(i))]
!
!-----------------------------------------------------------------------
!
      use number_types
      use number_types_pc
      use local_dims
      use mpidefs
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer, parameter :: IDIAG=7
!
!-----------------------------------------------------------------------
!
      real (r_typ_pc) :: Acsr(M)
      real (r_typ) :: Adia(N,IDIAG)
      integer :: N,M
      integer :: Adptr(N)
      integer :: IA(N+1)
      integer :: JA(M)
      integer :: ioff(IDIAG)
!
!-----------------------------------------------------------------------
!
      integer :: i,j,jj,mi,mj,mk,ko,x
      integer :: ioffok(IDIAG)
!
!-----------------------------------------------------------------------
!
      x=0
!
      IA(1)=1
      ko=1
      i=0
!
      do mk=2,npm1
        do mj=2,ntm1
          do mi=2,nrm1
! ********* Set index of value and column indicies array:
            i=i+1
!
! ********* Do not add coefs that multiply boundaries:
!           For each boundary, there is a sub-set of coefs in the
!           matrix row that should not be added.
!           This makes "local" matrices have no bc info
!
! ********* Reset "i-offset-ok-to-use-coef-jj" array:
!
            ioffok(:)=1
!
            if (mi.eq.2) then
              ioffok(3)=0;
            endif
!
            if (mi.eq.nrm1) then
              ioffok(5)=0;
            endif
!
            if (mj.eq.2) then
              ioffok(2)=0;
            endif
!
            if (mj.eq.ntm1) then
              ioffok(6)=0;
            endif
!
! ********* Eliminate periodic ceofs in the case nproc_p>1
!
            if (nproc_p.gt.1) then
              if (mk.eq.2) then
                ioffok(1)=0
              endif
              if (mk.eq.npm1) then
                ioffok(7)=0
              endif
            endif
!
! ********* To handle periodicity of phi in nproc_p=1 case:
!           We want CSR matrix to be in order so
!           have to sweep three times to avoid sorting:
!
! ********* Add periodic coefs of "right side":
!
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
!
! ********* Now do non-periodic coefs:
!
            do jj=1,IDIAG
              if (ioffok(jj).eq.1) then
                j=i+ioff(jj)-x
                if (j.ge.1.and.j.le.N-x) then
!                 Store pointer to diagonal elements in A:
                  if (jj.eq.4) Adptr(i)=ko
                  Acsr(ko)=real(Adia(i,jj),r_typ_pc)
                  JA(ko)=j
                  ko=ko+1
                endif
              endif
            enddo
!
! ********* Now do periodic coefs of "left side":
!
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
!
! ********* Set row offset:
!
            IA(i+1)=ko-x
          enddo
        enddo
      enddo
!
end subroutine
!#######################################################################
subroutine getM (N, ioff, M)
!
!-----------------------------------------------------------------------
!
! *** This routine computes the number of non-zeros in the
!     solver matrix for use with allocating the matrices.
!     See diacsr() for description of inputs.
!
!     Output:  M  # of nonzeros.
!
!-----------------------------------------------------------------------
!
      use mpidefs
      use local_dims
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer, parameter :: IDIAG=7
      integer :: N,M,i,j,jj,ko,mi,mj,mk,x
      integer :: ioff(IDIAG)
      integer :: ioffok(IDIAG)
!
      x=0
!
      ko=1
      i=0
!
      do mk=2,npm1
        do mj=2,ntm1
          do mi=2,nrm1
!
            ioffok(:)=1
!
            if (mi.eq.2) then
              ioffok(3)=0;
            endif
!
            if (mi.eq.nrm1) then
              ioffok(5)=0;
            endif
!
            if (mj.eq.2) then
              ioffok(2)=0;
            endif
!
            if (mj.eq.ntm1) then
              ioffok(6)=0;
            endif
!
! ********* Eliminate periodic ceofs in the case nproc_p>1
!
            if (nproc_p.gt.1) then
              if (mk.eq.2) then
                ioffok(1)=0
              endif
              if (mk.eq.npm1) then
                ioffok(7)=0
              endif
            endif
!
            do jj=1,IDIAG
              if (ioffok(jj).eq.1) then
                j=i+ioff(jj)-x
                if (j.gt.N-x) then
                  ko=ko+1
                endif
              endif
            enddo
!
            do jj=1,IDIAG
              if (ioffok(jj).eq.1) then
                j=i+ioff(jj)-x
                if (j.ge.1.and.j.le.N-x) then
                  ko=ko+1
                endif
              endif
            enddo
!
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
!
! *** Save number of non-zeros of matrix:
!
      M=ko-1
!
end subroutine
!#######################################################################
subroutine ax (x,y,N)
!
!-----------------------------------------------------------------------
!
! ****** Set y = A * x.
!
!-----------------------------------------------------------------------
!
      use number_types
      use local_dims
      use fields, ONLY : x_ax
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: zero=0._r_typ
!
!-----------------------------------------------------------------------
!
      integer :: N
      real(r_typ), dimension(N) :: x,y
!
!-----------------------------------------------------------------------
!
! ****** Expand X array to allow for boundary and seam values.
!
      call unpack_scalar (x_ax,x)
!
! ****** Set the boundary values of X.
!
      call set_boundary_points (x_ax,zero)
!
! ****** Seam along edges between processors.
!
      call seam_hhh (x_ax)
!
! ****** Get the matrix-vector product.
!
      call timer_on
!
      call delsq (x_ax,y)
!
      call timer_off (t_ax)
!
end subroutine
!#######################################################################
subroutine prec_inv (x)
!
!-----------------------------------------------------------------------
!
! ****** Apply preconditioner: x := M(inv) * x.
!
!-----------------------------------------------------------------------
!
      use number_types
      use number_types_pc
      use cgcom
      use solve_params
      use matrix_storage_pot3d_solve
      use, intrinsic :: iso_c_binding
      use cusparse_interface
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(N) :: x
      real(r_typ_pc), dimension(N) :: x_32
      integer :: i
!
!-----------------------------------------------------------------------
!
      if (ifprec.eq.0) return
!
      call timer_on
!
      if (ifprec.eq.1) then
!
! ****** Point-Jacobi (diagonal scaling):
!
        do concurrent (i=1:N)
          x(i)=a_i(i)*x(i)
        enddo
!
      elseif (ifprec.eq.2) then
!
! ****** ILU0 Partial-Block-Jacobi:
!
#ifdef CUSPARSE
!$omp target data use_device_ptr(x)
        call lusol_cusparse(C_LOC(x(1)))
!$omp end target data
!
#else
!$omp target update from(x)
!
! ****** Convert input array to single precision.
!
        do i=1,N
          x_32(i) = real(x(i),r_typ_pc)
        enddo
!
        call lusol (N,M,x_32,lu_csr,lu_csr_ja,a_N1,a_N2,a_csr_d)
!
! ****** Convert result back to double precision.
!
        do i=1,N
          x(i) = real(x_32(i),r_typ)
        enddo
!
!$omp target update to(x)
#endif
      endif
!
      call timer_off (t_pc)
!
end subroutine
!#######################################################################
subroutine lusol (N,M,x,LU,LU_ja,N1,N2,LUd_i)
!
!-----------------------------------------------------------
!
!     Performs a forward and a backward solve for the sparse system
!     (LU) x=y  where LU is in an optimized custom CSR format
!                                              (see lu2luopt())
!
!     For use where LU is an ILU or SSOR/SGS factorization.
!
!     Author of original lusol: Yousef Saad
!           Iterative Methods for Sparse Linear Systems 2nd Ed. pg. 299
!
!     Modified by RM Caplan to include optimized memory access
!     as described in
!     B. Smith, H. Zhang  Inter. J. of High Perf. Comp. Appl.
!     Vol. 25 #4 pg. 386-391 (2011)
!
!-----------------------------------------------------------
!     PARAMETERS:
!     N     : Dimension of problem
!     x     : At input, x is rhs (y), at output x is the solution.
!     LU    : Values of the LU matrix. L and U are stored together in
!             order of access in this routine.
!     LU_ja : Column indices of elements in LU.
!     N1    : Row-start indicies in original CSR LU.
!     N2    : Indices of diagonal elements in orig CSR LU
!     LUd_i : Inverse diagonal elements of U
!------------------------------------------------------------
!
      use number_types
      use number_types_pc
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ_pc) :: x(N),LUd_i(N),LU(M)
      integer :: N1(N),N2(N),LU_ja(M)
      integer :: N,M
!
!-----------------------------------------------------------------------
!
      integer :: i,k,k1,k2
!
!-----------------------------------------------------------------------
!
! ****** FORWARD SOLVE: Solve L x'=y
!
      k2=0
      do i=1,N
!       Compute x(i) := x(i) - sum L(i,j) * x(j)
        k1=k2+1
        k2=k1+N1(i)
        do k=k1,k2
          x(i)=x(i)-LU(k)*x(LU_ja(k))
        enddo
!       Diagonal is always 1 for L so no division here is nessesary.
      enddo
!
! ****** BACKWARD SOLVE: Solve U x=x'
!
      do i=N,1,-1
!       Compute x(i) := x(i) - sum U(i,j) * x(j)
        k1=k2+1
        k2=k1+N2(i)
        do k=k1,k2
          x(i)=x(i)-LU(k)*x(LU_ja(k))
        enddo
!       Compute x(i) := x(i) / U(i,i)
        x(i)=x(i)*LUd_i(i)
      enddo
!
end subroutine
!#######################################################################
subroutine unpack_scalar (s,x)
!
!-----------------------------------------------------------------------
!
! ****** Unpack the scalar x into
! ****** three-dimensional array s leaving room for boundaries.
!
!-----------------------------------------------------------------------
!
      use number_types
      use local_dims
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(nr,nt,np) :: s
      real(r_typ), dimension(2:nrm1,2:ntm1,2:npm1) :: x
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k
!
!-----------------------------------------------------------------------
!
      do concurrent (k=2:npm1, j=2:ntm1, i=2:nrm1)
        s(i,j,k)=x(i,j,k)
      enddo
!
end subroutine
!#######################################################################
subroutine delsq (x,y)
!
!-----------------------------------------------------------------------
!
! ****** Set Y = - (dV * del-squared X) at the internal points.
!
!-----------------------------------------------------------------------
!
      use number_types
      use local_dims
      use matrix_storage_pot3d_solve
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(nr,nt,np) :: x
      real(r_typ), dimension(2:nrm1,2:ntm1,2:npm1) :: y
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k
!
!-----------------------------------------------------------------------
!
      do concurrent (k=2:npm1, j=2:ntm1, i=2:nrm1)
        y(i,j,k)=a(i,j,k,1)*x(i  ,j  ,k-1) &
                +a(i,j,k,2)*x(i  ,j-1,k  ) &
                +a(i,j,k,3)*x(i-1,j  ,k  ) &
                +a(i,j,k,4)*x(i  ,j  ,k  ) &
                +a(i,j,k,5)*x(i+1,j  ,k  ) &
                +a(i,j,k,6)*x(i  ,j+1,k  ) &
                +a(i,j,k,7)*x(i  ,j  ,k+1)
      enddo
!
end subroutine
!#######################################################################
subroutine set_boundary_points (x,vmask)
!
!-----------------------------------------------------------------------
!
! ****** Set boundary points of X at the physical boundaries.
!
!-----------------------------------------------------------------------
!
      use number_types
      use global_mesh
      use local_dims
      use local_mesh
      use fields
      use solve_params
      use vars, ONLY : validation_run
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(nr,nt,np) :: x
      real(r_typ) :: vmask,temp_sum0
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: two=2._r_typ
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k
!
!-----------------------------------------------------------------------
!
! ****** Set X at the radial boundaries.
!
      if (rb0) then
        do concurrent (k=2:npm1, j=2:ntm1)
          x( 1,j,k)=x(2,j,k)-vmask*br0(j,k)*dr1
        enddo
      end if
!
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
!
! ****** If this processor does not contain any points at the
! ****** pole, return.
!
      if (.not.(tb0.or.tb1)) return
!
! ****** Get the m=0 component of X at the poles.
!
      if (tb0) then
        do concurrent (i=1:nr)
          sum0(i)=0
        enddo
!$omp target teams loop
        do i=1,nr
          temp_sum0=0.
!$omp loop reduction(+:temp_sum0)
          do k=2,npm1
            temp_sum0=temp_sum0+x(i,2,k)*dph(k)*pl_i
          enddo
          sum0(i)=temp_sum0
        enddo
      end if
!
      if (tb1) then
        do concurrent (i=1:nr)
          sum1(i)=0
        enddo
!$omp target teams loop
        do i=1,nr
          temp_sum0=0.
!$omp loop reduction(+:temp_sum0)
          do k=2,npm1
            temp_sum0=temp_sum0+x(i,ntm1,k)*dph(k)*pl_i
          enddo
          sum1(i)=temp_sum0
        enddo
      end if
!
! ****** Sum over all processors.
!
      call sum_over_phi (nr,sum0,sum1)
!
! ****** Set X to have only an m=0 component at the poles.
!
      if (tb0) then
        do concurrent (k=2:npm1, i=1:nr)
          x(i,1,k)=two*sum0(i)-x(i,2,k)
        enddo
      end if
!
      if (tb1) then
        do concurrent (k=2:npm1, i=1:nr)
          x(i,nt,k)=two*sum1(i)-x(i,ntm1,k)
        enddo
      end if
!
      return
end subroutine
!#######################################################################
subroutine sum_over_phi (n,a0,a1)
!
!-----------------------------------------------------------------------
!
! ****** Sum the contribution over all processors in the phi
! ****** dimension (only for processors with points on the poles).
!
! ****** The sum is performed for all N points in the vectors
! ****** SUM0(N) and SUM1(N), at the North and South pole,
! ****** respectively.
!
!-----------------------------------------------------------------------
!
      use number_types
      use local_dims
      use mpidefs
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: n
      real(r_typ), dimension(n) :: a0,a1
!
!-----------------------------------------------------------------------
!
! ****** MPI error return.
!
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      call timer_on
!
!$omp target data use_device_ptr(a0,a1)
      if (tb0) then
        call MPI_Allreduce (MPI_IN_PLACE,a0,n,ntype_real, &
                            MPI_SUM,comm_phi,ierr)
      end if
!
      if (tb1) then
        call MPI_Allreduce (MPI_IN_PLACE,a1,n,ntype_real, &
                            MPI_SUM,comm_phi,ierr)
      end if
!$omp end target data
!
      call timer_off (c_sumphi)
!
end subroutine
!#######################################################################
subroutine zero_boundary_points (x)
!
!-----------------------------------------------------------------------
!
! ****** Set the boundary points at the physical boundaries
! ****** of X to zero.
!
!-----------------------------------------------------------------------
!
      use number_types
      use local_dims
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(nr,nt,np) :: x
!
!-----------------------------------------------------------------------
!
      if (rb0) x( 1,:,:)=0.
      if (rb1) x(nr,:,:)=0.
      if (tb0) x(:, 1,:)=0.
      if (tb1) x(:,nt,:)=0.
!
end subroutine
!#######################################################################
function cgdot (x,y,N)
!
!-----------------------------------------------------------------------
!
! ****** Get the dot product of the vectors X and Y.
!
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: N,i
      real(r_typ) :: cgdot
      real(r_typ), dimension(N) :: x,y
!
!-----------------------------------------------------------------------
!
      cgdot=0.
!
!$omp target teams loop reduction(+:cgdot)
      do i=1,N
        cgdot=cgdot+x(i)*y(i)
      enddo
!
! ****** Sum over all the processors.
!
      call global_sum (cgdot)
!
      return
end function
!#######################################################################
subroutine global_sum (x)
!
!-----------------------------------------------------------------------
!
! ****** Overwrite X by the its sum over all processors.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: x
!
!-----------------------------------------------------------------------
!
! ****** MPI error return.
!
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      call timer_on
!
! ****** Take the sum over all the processors.
!
      call MPI_Allreduce (MPI_IN_PLACE,x,1,ntype_real, &
                          MPI_SUM,comm_all,ierr)
!
      call timer_off (c_cgdot)
!
end subroutine
!#######################################################################
subroutine seam_setup
!
!-----------------------------------------------------------------------
!
! ****** Allocate seam buffers.
! ****** Doing this here avoids GPU allocate/deallocates.
!
!-----------------------------------------------------------------------
!
      use number_types
      use cgcom
      use local_mesh
      use local_dims
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      allocate(sbuf_rt1(nr,nt))
      allocate(sbuf_rt2(nr,nt))
      allocate(sbuf_tp1(nt,np))
      allocate(sbuf_tp2(nt,np))
      allocate(sbuf_rp1(nr,np))
      allocate(sbuf_rp2(nr,np))
!
      allocate(rbuf_rt1(nr,nt))
      allocate(rbuf_rt2(nr,nt))
      allocate(rbuf_tp1(nt,np))
      allocate(rbuf_tp2(nt,np))
      allocate(rbuf_rp1(nr,np))
      allocate(rbuf_rp2(nr,np))
!
!$omp target enter data map(alloc:sbuf_rt1,sbuf_rt2)
!$omp target enter data map(alloc:sbuf_tp1,sbuf_tp2)
!$omp target enter data map(alloc:sbuf_rp1,sbuf_rp2)
!$omp target enter data map(alloc:rbuf_rt1,rbuf_rt2)
!$omp target enter data map(alloc:rbuf_tp1,rbuf_tp2)
!$omp target enter data map(alloc:rbuf_rp1,rbuf_rp2)
!
end subroutine
!#######################################################################
subroutine seam_hhh (a)
!
!-----------------------------------------------------------------------
!
! ****** Seam the boundary points of 3D (r,t,p) array A between
! ****** adjacent processors.
!
! ****** This routine assumes that there is a two-point
! ****** overlap between processors in each dimension.
!
!-----------------------------------------------------------------------
!
! ****** This version uses non-blocking MPI sends and receives
! ****** whenever possible in order to overlap communications.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
      use timing
      use local_mesh
      use local_dims
      use cgcom
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(nr,nt,np) :: a
!
!-----------------------------------------------------------------------
!
! ****** MPI error return.
!
      integer :: ierr
!
! ****** MPI tag for MPI_ISEND and MPI_IRECV (not tagged).
!
      integer :: tag=0
!
!-----------------------------------------------------------------------
!
      integer :: lbuf,i,j
      integer :: reqs(4)
!
!-----------------------------------------------------------------------
!
      call timer_on
!
! ****** Seam the third (periodic) dimension.
! ****** Since halo data is stride-1, no need for buffers.
!
      lbuf=nr*nt
!
!$omp target data use_device_ptr(a)
      call MPI_Isend (a(:,:,np-1),lbuf,ntype_real,iproc_pp,tag, &
                      comm_all,reqs(1),ierr)
!
      call MPI_Isend (a(:,:,   2),lbuf,ntype_real,iproc_pm,tag, &
                      comm_all,reqs(2),ierr)
!
      call MPI_Irecv (a(:,:, 1),lbuf,ntype_real,iproc_pm,tag,   &
                      comm_all,reqs(3),ierr)
!
      call MPI_Irecv (a(:,:,np),lbuf,ntype_real,iproc_pp,tag,   &
                      comm_all,reqs(4),ierr)
!
      call MPI_Waitall (4,reqs,MPI_STATUSES_IGNORE,ierr)
!$omp end target data
!
! ****** Seam the first dimension.
!
      if (nproc_r.gt.1) then
!
        lbuf=nt*np
!
        do concurrent (j=1:np, i=1:nt)
          sbuf_tp1(i,j)=a(nr-1,i,j)
          sbuf_tp2(i,j)=a(   2,i,j)
        enddo
!
!$omp target data use_device_ptr(sbuf_tp1,sbuf_tp2,rbuf_tp1,rbuf_tp2)
        call MPI_Isend (sbuf_tp1,lbuf,ntype_real,iproc_rp,tag, &
                        comm_all,reqs(1),ierr)
!
        call MPI_Isend (sbuf_tp2,lbuf,ntype_real,iproc_rm,tag, &
                        comm_all,reqs(2),ierr)
!
        call MPI_Irecv (rbuf_tp1,lbuf,ntype_real,iproc_rm,tag, &
                        comm_all,reqs(3),ierr)
!
        call MPI_Irecv (rbuf_tp2,lbuf,ntype_real,iproc_rp,tag, &
                        comm_all,reqs(4),ierr)
!
        call MPI_Waitall (4,reqs,MPI_STATUSES_IGNORE,ierr)
!$omp end target data
!
        if (iproc_rm.ne.MPI_PROC_NULL) then
          do concurrent (j=1:np, i=1:nt)
            a( 1,i,j)=rbuf_tp1(i,j)
          enddo
        end if
!
        if (iproc_rp.ne.MPI_PROC_NULL) then
          do concurrent (j=1:np, i=1:nt)
            a(nr,i,j)=rbuf_tp2(i,j)
          enddo
        end if
      end if
!
! ****** Seam the second dimension.
!
      if (nproc_t.gt.1) then
!
        lbuf=nr*np
!
        do concurrent (j=1:np, i=1:nr)
          sbuf_rp1(i,j)=a(i,nt-1,j)
          sbuf_rp2(i,j)=a(i,   2,j)
        enddo
!
!$omp target data use_device_ptr(sbuf_rp1,sbuf_rp2,rbuf_rp1,rbuf_rp2)
        call MPI_Isend (sbuf_rp1,lbuf,ntype_real,iproc_tp,tag, &
                        comm_all,reqs(1),ierr)
!
        call MPI_Isend (sbuf_rp2,lbuf,ntype_real,iproc_tm,tag, &
                        comm_all,reqs(2),ierr)
!
        call MPI_Irecv (rbuf_rp1,lbuf,ntype_real,iproc_tm,tag, &
                        comm_all,reqs(3),ierr)
!
        call MPI_Irecv (rbuf_rp2,lbuf,ntype_real,iproc_tp,tag, &
                        comm_all,reqs(4),ierr)
!
        call MPI_Waitall (4,reqs,MPI_STATUSES_IGNORE,ierr)
!$omp end target data
!
        if (iproc_tm.ne.MPI_PROC_NULL) then
          do concurrent (j=1:np, i=1:nr)
            a(i, 1,j)=rbuf_rp1(i,j)
          enddo
        end if
!
        if (iproc_tp.ne.MPI_PROC_NULL) then
          do concurrent (j=1:np, i=1:nr)
            a(i,nt,j)=rbuf_rp2(i,j)
          enddo
        end if
!
      end if
!
      call timer_off (c_seam)
!
end subroutine
!#######################################################################
subroutine seam_gen (a,n1,n2,n3)
!
!-----------------------------------------------------------------------
!
! ****** Seam the boundary points of 3D (r,t,p) array A between
! ****** adjacent processors.
!
! ****** This routine assumes that there is a two-point
! ****** overlap between processors in each dimension.
!
!-----------------------------------------------------------------------
!
! ****** This version uses non-blocking MPI sends and receives
! ****** whenever possible in order to overlap communications.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(n1,n2,n3) :: a
!
!-----------------------------------------------------------------------
!
      real(r_typ), dimension(n2,n3) :: sbuf11,rbuf11
      real(r_typ), dimension(n2,n3) :: sbuf12,rbuf12
      real(r_typ), dimension(n1,n3) :: sbuf21,rbuf21
      real(r_typ), dimension(n1,n3) :: sbuf22,rbuf22
!
!-----------------------------------------------------------------------
!
! ****** MPI error return.
!
      integer :: ierr
!
! ****** MPI tag for MPI_ISEND and MPI_IRECV (not tagged).
!
      integer :: tag=0
!
!-----------------------------------------------------------------------
!
      integer :: lbuf,i,j
      integer :: n1,n2,n3
      integer :: reqs(4)
!
!-----------------------------------------------------------------------
!
      call timer_on
!
! ****** Seam the third (periodic) dimension.
! ****** Since halo data is stride-1, no need for buffers.
!
      lbuf=n1*n2
!
!$omp target data use_device_ptr(a)
      call MPI_Isend (a(:,:,n3-1),lbuf,ntype_real,iproc_pp,tag, &
                      comm_all,reqs(1),ierr)
!
      call MPI_Isend (a(:,:,   2),lbuf,ntype_real,iproc_pm,tag, &
                      comm_all,reqs(2),ierr)
!
      call MPI_Irecv (a(:,:, 1),lbuf,ntype_real,iproc_pm,tag,   &
                      comm_all,reqs(3),ierr)
!
      call MPI_Irecv (a(:,:,n3),lbuf,ntype_real,iproc_pp,tag,   &
                      comm_all,reqs(4),ierr)
!
      call MPI_Waitall (4,reqs,MPI_STATUSES_IGNORE,ierr)
!$omp end target data
!
! ****** Seam the first dimension.
!
      if (nproc_r.gt.1) then
!
!$omp target enter data map(alloc:sbuf11,sbuf12,rbuf11,rbuf12)
!
        lbuf=n2*n3
!
        do concurrent (j=1:n3, i=1:n2)
          sbuf11(i,j)=a(n1-1,i,j)
          sbuf12(i,j)=a(   2,i,j)
        enddo
!
!$omp target data use_device_ptr(sbuf11,sbuf12,rbuf11,rbuf12)
        call MPI_Isend (sbuf11,lbuf,ntype_real,iproc_rp,tag, &
                        comm_all,reqs(1),ierr)
!
        call MPI_Isend (sbuf12,lbuf,ntype_real,iproc_rm,tag, &
                        comm_all,reqs(2),ierr)
!
        call MPI_Irecv (rbuf11,lbuf,ntype_real,iproc_rm,tag, &
                        comm_all,reqs(3),ierr)
!
        call MPI_Irecv (rbuf12,lbuf,ntype_real,iproc_rp,tag, &
                        comm_all,reqs(4),ierr)
!
        call MPI_Waitall (4,reqs,MPI_STATUSES_IGNORE,ierr)
!$omp end target data
!
        if (iproc_rm.ne.MPI_PROC_NULL) then
          do concurrent (j=1:n3, i=1:n2)
            a( 1,i,j)=rbuf11(i,j)
          enddo
        end if
!
        if (iproc_rp.ne.MPI_PROC_NULL) then
          do concurrent (j=1:n3, i=1:n2)
            a(n1,i,j)=rbuf12(i,j)
          enddo
        end if
!
!$omp target exit data map(delete:sbuf11,sbuf12,rbuf11,rbuf12)
      end if
!
! ****** Seam the second dimension.
!
      if (nproc_t.gt.1) then
!
!$omp target enter data map(alloc:sbuf21,sbuf22,rbuf21,rbuf22)
!
        lbuf=n1*n3
!
        do concurrent (j=1:n3, i=1:n1)
          sbuf21(i,j)=a(i,n2-1,j)
          sbuf22(i,j)=a(i,   2,j)
        enddo
!
!$omp target data use_device_ptr(sbuf21,sbuf22,rbuf21,rbuf22)
        call MPI_Isend (sbuf21,lbuf,ntype_real,iproc_tp,tag, &
                        comm_all,reqs(1),ierr)
!
        call MPI_Isend (sbuf22,lbuf,ntype_real,iproc_tm,tag, &
                        comm_all,reqs(2),ierr)
!
        call MPI_Irecv (rbuf21,lbuf,ntype_real,iproc_tm,tag, &
                        comm_all,reqs(3),ierr)
!
        call MPI_Irecv (rbuf22,lbuf,ntype_real,iproc_tp,tag, &
                        comm_all,reqs(4),ierr)
!
        call MPI_Waitall (4,reqs,MPI_STATUSES_IGNORE,ierr)
!$omp end target data
!
        if (iproc_tm.ne.MPI_PROC_NULL) then
          do concurrent (j=1:n3, i=1:n1)
            a(i, 1,j)=rbuf21(i,j)
          enddo
        end if
!
        if (iproc_tp.ne.MPI_PROC_NULL) then
          do concurrent (j=1:n3, i=1:n1)
            a(i,n2,j)=rbuf22(i,j)
          enddo
        end if
!
!$omp target exit data map(delete:sbuf21,sbuf22,rbuf21,rbuf22)
      end if
!
      call timer_off (c_seam)
!
end subroutine
!#######################################################################
subroutine write_solution
!
!-----------------------------------------------------------------------
!
! ****** Write the global solution.
!
!-----------------------------------------------------------------------
!
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
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
! ****** Declaration for the global arrays.
! ****** These arrays are only allocated on processor IPROC0.
!
      real(r_typ), dimension(:,:,:), allocatable :: phi_g
      real(r_typ), dimension(:,:,:), allocatable :: br_g
      real(r_typ), dimension(:,:,:), allocatable :: bt_g
      real(r_typ), dimension(:,:,:), allocatable :: bp_g
!
!-----------------------------------------------------------------------
!
      integer :: ierr
      character(256) :: fname
!
!-----------------------------------------------------------------------
!
      call timer_on
!
! ****** Potential.
!
      if (phifile.ne.'') then
!$omp target update from(phi)
!
! ****** Allocate the global array PHI_G (on processor IPROC0).
!
        if (iamp0) then
          allocate (phi_g(nr_g,nt_g,np_g))
        else
          allocate (phi_g(1,1,1))
        end if
!
! ****** Assemble the global PHI array.
!
        call assemble_array (map_rh,map_th,map_ph,phi,phi_g)
!
        fname=phifile
!
! ****** Write out the potential to a file.
!
        if (iamp0) then
            write (*,*)
            write (*,*) '### COMMENT from WRITE_SOLUTION:'
            write (*,*)
            write (*,*) 'Writing the potential to file: ',trim(fname)
            ! [XX: HDF5 file write]
            ! call wrhdf_3d (fname,.true.,nr_g,nt_g,np_g, &
            !                phi_g,rh_g,th_g,ph_g,hdf32,ierr)
        end if
!
        deallocate (phi_g)
!
      end if
!
! ****** Br.
!
      if (brfile.ne.'') then
!$omp target update from(br)
!
        fname=brfile
!
        if (iamp0) then
          allocate (br_g(nrm1_g,nt_g,np_g))
        else
          allocate (br_g(1,1,1))
        end if
!
! ****** Assemble the global PHI array.
!
        call assemble_array (map_rm,map_th,map_ph,br,br_g)
!
        if (iamp0) then
          write (*,*)
          write (*,*) '### COMMENT from WRITE_SOLUTION:'
          write (*,*)
          write (*,*) 'Writing Br to file: ',trim(fname)
      ! [XX: HDF5 file write]
      !     call wrhdf_3d (fname,.true.,nrm1_g,nt_g,np_g, &
      !                    br_g,r_g,th_g,ph_g,hdf32,ierr)
        end if
!
        deallocate (br_g)
!
      end if
!
! ****** Bt.
!
      if (btfile.ne.'') then
!$omp target update from(bt)
!
        fname=btfile
!
        if (iamp0) then
          allocate (bt_g(nr_g,ntm1_g,np_g))
        else
          allocate (bt_g(1,1,1))
        end if
!
! ****** Assemble the global PHI array.
!
        call assemble_array (map_rh,map_tm,map_ph,bt,bt_g)
!
        if (iamp0) then
          write (*,*)
          write (*,*) '### COMMENT from WRITE_SOLUTION:'
          write (*,*)
          write (*,*) 'Writing Bt to file: ',trim(fname)
      !     call wrhdf_3d (fname,.true.,nr_g,ntm1_g,np_g, &
      !                    bt_g,rh_g,t_g,ph_g,hdf32,ierr)
!
        end if
!
        deallocate (bt_g)
!
      end if
!
! ****** Bp.
!
      if (bpfile.ne.'') then
!$omp target update from(bp)
!
        fname=bpfile
!
        if (iamp0) then
          allocate (bp_g(nr_g,nt_g,npm1_g))
        else
          allocate (bp_g(1,1,1))
        end if
!
! ****** Assemble the global PHI array.
!
        call assemble_array (map_rh,map_th,map_pm,bp,bp_g)
!
        if (iamp0) then
          write (*,*)
          write (*,*) '### COMMENT from WRITE_SOLUTION:'
          write (*,*)
          write (*,*) 'Writing Bp to file: ',trim(fname)
      !     call wrhdf_3d (fname,.true.,nr_g,nt_g,npm1_g, &
      !                    bp_g,rh_g,th_g,p_g,hdf32,ierr)
!
        end if
!
        deallocate (bp_g)
!
      end if
!
      call timer_off (t_io)
!
end subroutine
!#######################################################################
subroutine getb
!
!-----------------------------------------------------------------------
!
! ****** Calculate B from grad-phi.
!
!-----------------------------------------------------------------------
!
      use number_types
      use global_dims
      use global_mesh
      use vars
      use fields
      use local_dims
      use local_mesh
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k
!
!-----------------------------------------------------------------------
!
!$omp target enter data map(alloc:br,bt,bp)
!
! ****** Get Br.
!
      do concurrent (k=1:np, j=1:nt, i=1:nrm1)
        br(i,j,k)=(phi(i+1,j,k)-phi(i,j,k))/dr(i)
      enddo
!
! ****** Get Bt.
!
      do concurrent (k=1:np, j=1:ntm1, i=1:nr)
        bt(i,j,k)=(phi(i,j+1,k)-phi(i,j,k))/(rh(i)*dt(j))
      enddo
!
! ****** Get Bp.
!
      do concurrent (k=1:npm1, j=1:nt, i=1:nr)
        bp(i,j,k)=(phi(i,j,k+1)-phi(i,j,k))/(rh(i)*sth(j)*dp(k))
      enddo
!
end subroutine
!#######################################################################
subroutine magnetic_energy
!
!-----------------------------------------------------------------------
!
! ****** Calculate magnetic energy from B.
!
!-----------------------------------------------------------------------
!
      use number_types
      use global_dims
      use global_mesh
      use vars
      use fields
      use mpidefs
      use local_dims
      use local_mesh
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: half=.5_r_typ
!
!-----------------------------------------------------------------------
!
      integer :: i,j,k,ierr
      real(r_typ) :: brav,btav,bpav,dv
      real(r_typ) :: wr,wt,wp
      real(r_typ), dimension(3) :: w
      character(32) :: fmtstr
!
!-----------------------------------------------------------------------
!
      if (hdf32) then
        fmtstr="(A, ES14.8)"
      else
        fmtstr="(A,ES22.16)"
      end if
!
      wr=0.
      wt=0.
      wp=0.
!$omp target teams loop collapse(3) reduction(+:wr,wt,wp)
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
!
! ****** Sum up all processors into final values and print.
!
      w(1)=wr
      w(2)=wt
      w(3)=wp
      call MPI_Allreduce(MPI_IN_PLACE,w,3,ntype_real, &
                         MPI_SUM,comm_all,ierr)
!
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
!
end subroutine
!#######################################################################
subroutine assemble_array (map_r,map_t,map_p,a,a_g)
!
!-----------------------------------------------------------------------
!
! ****** Assemble a global array (into A_G) on processor IPROC0 by
! ****** fetching the local sections (A) from all the processors.
!
!-----------------------------------------------------------------------
!
      use number_types
      use decomposition
      use mpidefs
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      type(map_struct), dimension(0:nproc-1) :: map_r,map_t,map_p
      real(r_typ), dimension(:,:,:) :: a,a_g
!
!-----------------------------------------------------------------------
!
! ****** Storage for the buffers.
!
      integer :: lbuf,lsbuf
      real(r_typ), dimension(:), allocatable :: sbuf
      real(r_typ), dimension(:), allocatable :: rbuf
!
!-----------------------------------------------------------------------
!
      integer :: tag=0
      integer :: irank,l1,l2,l3,i,j,k,ii
      integer :: i0,j0,k0,i1,j1,k1
      integer :: i0g,j0g,k0g
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      l1=map_r(iproc)%n
      l2=map_t(iproc)%n
      l3=map_p(iproc)%n
      lsbuf=l1*l2*l3
!
      i0=map_r(iproc)%i0
      i1=map_r(iproc)%i1
      j0=map_t(iproc)%i0
      j1=map_t(iproc)%i1
      k0=map_p(iproc)%i0
      k1=map_p(iproc)%i1
!
! ***** Extract 1D array of non-overlapping values from local array.
!
      allocate (sbuf(lsbuf))
!
      sbuf=reshape(a(i0:i1,j0:j1,k0:k1),(/lsbuf/))
!
! ****** If proc0, recieve/store local arrays into global array.
!
      if (iamp0) then
        do irank=0,nproc-1
!
          l1=map_r(irank)%n
          l2=map_t(irank)%n
          l3=map_p(irank)%n
          lbuf=l1*l2*l3
!
          i0g=map_r(irank)%offset
          j0g=map_t(irank)%offset
          k0g=map_p(irank)%offset
!
! ****** If proc0 is the current rank in loop, simply copy local array.
          if (iproc==irank) then
            do k=1,l3
              do j=1,l2
                do i=1,l1
                  ii=l2*l1*(k-1)+l1*(j-1)+i
                  a_g(i0g+i-1,j0g+j-1,k0g+k-1)=sbuf(ii)
                enddo
              enddo
            enddo
! ****** Otherwise recieve data:
          else
            allocate (rbuf(lbuf))
            call MPI_Recv (rbuf,lbuf,ntype_real,irank,tag, &
                           comm_all,MPI_STATUS_IGNORE,ierr)
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
!
! ****** Send local array to iproc0.
!
      call MPI_Ssend (sbuf,lsbuf,ntype_real,iproc0,tag,comm_all,ierr)
!
      end if
      deallocate (sbuf)
!
end subroutine
!#######################################################################
subroutine timer_on
!
!-----------------------------------------------------------------------
!
! ****** Push an entry onto the timing stack and initialize
! ****** a timing event.
!
!-----------------------------------------------------------------------
!
! ****** This routine can be called in a nested way to measure
! ****** multiple timing events.  Calls to TIMER_ON and TIMER_OFF
! ****** need to be nested like do-loops in FORTRAN.
!
!-----------------------------------------------------------------------
!
      use mpidefs
      use timer
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      if (istack.ge.nstack) then
        write (*,*)
        write (*,*) '### WARNING from TIMER_ON:'
        write (*,*) '### Timing stack depth exceeded.'
        write (*,*) 'This may signal an incorrect nesting of '// &
                    'TIMER_ON/TIMER_OFF calls.'
        write (*,*) 'Timing information will not be valid.'
        return
      else
        istack=istack+1
      end if
!
      tstart(istack)=MPI_Wtime()
!
end subroutine
!#######################################################################
subroutine timer_off (tused)
!
!-----------------------------------------------------------------------
!
! ****** Increment the CPU time used since the call to TIMER_ON
! ****** in variable TUSED, and pop an entry off the timing
! ****** stack.
!
!-----------------------------------------------------------------------
!
! ****** This routine can be called in a nested way to measure
! ****** multiple timing events.  Calls to TIMER_ON and TIMER_OFF
! ****** need to be nested like do-loops in FORTRAN.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
      use timer
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: tused
!
!-----------------------------------------------------------------------
!
      if (istack.le.0) then
        write (*,*)
        write (*,*) '### WARNING from TIMER_OFF:'
        write (*,*) '### Timing stack cannot be popped.'
        write (*,*) 'This may signal an incorrect nesting of '// &
                    'TIMER_ON/TIMER_OFF calls.'
        write (*,*) 'Timing information will not be valid.'
        return
      else
        istack=istack-1
      end if
!
      tused=tused+MPI_Wtime()-tstart(istack+1)
!
end subroutine
!#######################################################################
subroutine write_timing
!
!-----------------------------------------------------------------------
!
! ****** Write out the timing info.
!
!-----------------------------------------------------------------------
!
      use number_types
      use mpidefs
      use timing
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
! ****** Timing buffers.
!
      integer, parameter :: lbuf=10
      real(r_typ), dimension(lbuf) :: sbuf
      real(r_typ), dimension(lbuf,0:nproc-1) :: tbuf
!
! ****** Timing statistics.
!
      real(r_typ), dimension(lbuf) :: tmin,tmax,tavg,tsdev
!
!-----------------------------------------------------------------------
!
      integer :: ierr,irank
      real(r_typ) :: t_tot_avg,c_tot_avg,c_tot
!
      character(80) :: tfile='timing.out'
!
!-----------------------------------------------------------------------
!
! ****** Gather the timing information for all processors into TBUF.
!
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
!
      call MPI_Allgather (sbuf,lbuf,ntype_real, &
                          tbuf,lbuf,ntype_real,comm_all,ierr)
!
! ****** Calculate the timing statistics.
!
      tavg=sum(tbuf,dim=2)/nproc
      tmin=minval(tbuf,dim=2)
      tmax=maxval(tbuf,dim=2)
!
      tsdev(:)=0.
      do irank=0,nproc-1
        tsdev(:)=tsdev(:)+(tbuf(:,irank)-tavg(:))**2
      enddo
      tsdev(:)=sqrt(tsdev(:)/nproc)
!
      t_tot_avg=tavg(10)
      c_tot_avg=tavg(7)+tavg(8)+tavg(9)
!
      if (iamp0) then
!
        call ffopen (1,tfile,'rw',ierr)
!
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### WARNING from WRITE_TIMING:'
          write (*,*) '### Could not create the timing file.'
          write (*,*) 'File name: ',trim(tfile)
        end if
!
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
!
        write (1,*)
        write (1,*) 'Average times:'
        write (1,*) '-------------'
        write (1,*)
        write (1,300) '    Avg         Min         Max      S. Dev'
        write (1,300) '    ---         ---         ---      ------'
        write (1,400) 'Avg comm. time in SEAM    = ', &
                      tavg(7),tmin(7),tmax(7),tsdev(7)
        write (1,400) 'Avg comm. time in CGDOT   = ', &
                      tavg(8),tmin(8),tmax(8),tsdev(8)
        write (1,400) 'Avg comm. time in SUMPHI  = ', &
                      tavg(9),tmin(9),tmax(9),tsdev(9)
        write (1,400) 'Avg time used in start-up = ', &
                      tavg(2),tmin(2),tmax(2),tsdev(2)
        write (1,400) 'Avg time used in i/o      = ', &
                      tavg(3),tmin(3),tmax(3),tsdev(3)
        write (1,400) 'Avg time used in POTFLD   = ', &
                      tavg(1),tmin(1),tmax(1),tsdev(1)
        write (1,400) 'Avg time used in AX       = ', &
                      tavg(4),tmin(4),tmax(4),tsdev(4)
        write (1,400) 'Avg time used in PCLOAD   = ', &
                      tavg(5),tmin(5),tmax(5),tsdev(5)
        write (1,400) 'Avg time used in PC       = ', &
                      tavg(6),tmin(6),tmax(6),tsdev(6)
        write (1,400) 'Avg total time            = ', &
                      tavg(10),tmin(10),tmax(10),tsdev(10)
  300   format (1x,33x,a)
  400   format (1x,a,4f12.3)
!
        write (1,*)
        write (1,200) 'Average time used per proc  = ',t_tot_avg
        write (1,200) 'Average comm. time per proc = ',c_tot_avg
        write (1,*)
        write (1,100)
        write (1,*)
!
        close (1)
!
        write(*,*)
        write(*,200) 'Total time:               ',tmax(10)
        write(*,200) 'Max Comm/Imbalance time:  ', &
                      tmax(7)+tmax(8)+tmax(9)
        write(*,*)
!
      end if
!
end subroutine
!#######################################################################
subroutine readbr (fname,br0_g,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Read in the radial magnetic field at the photosphere
! ****** and interpolate it into array BR0_G.
!
! ****** FNAME is the name of the file to read.
!
!-----------------------------------------------------------------------
!
      use number_types
      use global_dims
      use global_mesh
      use vars
      use rdhdf_2d_interface
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      character(*) :: fname
      real(r_typ), dimension(nt_g,np_g) :: br0_g
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: two=2._r_typ
!
!-----------------------------------------------------------------------
!
! ****** Br array read in and its scales.
!
      real(r_typ), dimension(:,:), pointer :: bn
      real(r_typ), dimension(:), pointer :: tn,pn
!
!-----------------------------------------------------------------------
!
      integer :: ntn,npn,j,k
      logical :: scale
      real(r_typ) :: sum0,sum1,area,fluxp,fluxm,da,br00err
!
!-----------------------------------------------------------------------
!
      ierr=0
!
! ****** Read in the normal field.
!
      write (*,*)
      write (*,*) '### COMMENT from READBR:'
      write (*,*) '### Reading Br file: ',trim(fname)
!
      call rdhdf_2d (fname,scale,ntn,npn,bn,tn,pn,ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READBR:'
        write (*,*) '### The flux file has the wrong format.'
        write (*,*) 'IERR (from RDHDF_2D) = ',ierr
        write (*,*) 'File name: ',trim(fname)
        ierr=1
        return
      end if
!
! ****** Check that the arrays has scales.
!
      if (.not.scale) then
        write (*,*)
        write (*,*) '### ERROR in READBR:'
        write (*,*) '### The flux file does not have scales.'
        write (*,*) 'File name: ',trim(fname)
        ierr=2
        return
      end if
!
! ****** Interpolate the field to the code mesh (into array BR0_G).
!
      call intrp2d (ntn,npn,tn,pn,bn, &
                    nt_g-2,np_g-2,th_g(2:ntm1_g),ph_g(2:npm1_g), &
                    br0_g(2:ntm1_g,2:npm1_g),ierr)
!
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in READBR:'
        write (*,*) '### The scales in the Br file are invalid.'
        write (*,*) 'File name: ',trim(fname)
        ierr=3
        return
      end if
!
! ****** De-allocate the memory for the BN array and its scales.
!
      deallocate (bn)
      deallocate (tn)
      deallocate (pn)
!
! ****** Set Br to be periodic.
!
      br0_g(:,1)=br0_g(:,npm1_g)
      br0_g(:,np_g)=br0_g(:,2)
!
! ****** Set BCs at the poles.
!
! ****** Br has only an m=0 component there.
!
      sum0=sum(br0_g(     2,2:npm1_g)*dph_g(2:npm1_g))*pl_i
      sum1=sum(br0_g(ntm1_g,2:npm1_g)*dph_g(2:npm1_g))*pl_i
!
      br0_g(1   ,:)=two*sum0-br0_g(     2,:)
      br0_g(nt_g,:)=two*sum1-br0_g(ntm1_g,:)
!
! ****** Calculate the total flux.
!
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
!
      write (*,*)
      write (*,*) '### COMMENT from READBR:'
      write (*,*) '### Computed flux balance:'
      write (*,*)
      write (*,*) 'Positive flux = ',fluxp
      write (*,*) 'Negative flux = ',fluxm
!
! ****** Fix the magnetic field so that the total flux is zero
! ****** (unless this has not been requested).
!
      if (.not.((option.eq.'ss'.or.option.eq.'open') &
                .and.do_not_balance_flux)) then
!
        br00err=(fluxp+fluxm)/area
!
        do k=1,np_g
          do j=1,nt_g
            br0_g(j,k)=br0_g(j,k)-br00err
          enddo
        enddo
!
        write (*,*)
        write (*,*) '### COMMENT from READBR:'
        write (*,*) '### Flux balance correction:'
        write (*,*)
        write (*,*) 'BR00 (monopole Br field magnitude) = ',br00err
!
      end if
!
end subroutine
!#######################################################################
subroutine intrp2d (nxi,nyi,xi,yi,fi,nx,ny,x,y,f,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Interpolate a 2D field from array FI(NXI,NYI), defined
! ****** on the mesh XI(NXI) x YI(NYI), into the array F(NX,NY),
! ****** defined on the mesh X(NX) x Y(NY).
!
! ****** Note that if a data point is outside the bounds of
! ****** the XI x YI mesh, IERR=2 is returned.
!
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: nxi,nyi
      real(r_typ), dimension(nxi) :: xi
      real(r_typ), dimension(nyi) :: yi
      real(r_typ), dimension(nxi,nyi) :: fi
      integer :: nx,ny
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      real(r_typ), dimension(nx,ny) :: f
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      real(r_typ), parameter :: zero=0._r_typ
!
!-----------------------------------------------------------------------
!
      real(r_typ), external :: flint
!
!-----------------------------------------------------------------------
!
      integer :: i,j,ii,jj
      real(r_typ) :: dum(1)=0.
      real(r_typ) :: dummy,xv,yv,ax,ay
!
!-----------------------------------------------------------------------
!
      ierr=0
!
! ****** Check that the scales XI and YI are monotonic.
!
      dummy=flint(zero,nxi,xi,dum,1,ierr)
      if (ierr.ne.0) go to 900
!
      dummy=flint(zero,nyi,yi,dum,1,ierr)
      if (ierr.ne.0) go to 900
!
! ****** Interpolate the data.
!
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
          f(i,j)= (1.-ax)*((1.-ay)*fi(ii  ,jj  )+ay*fi(ii  ,jj+1)) &
                 +    ax *((1.-ay)*fi(ii+1,jj  )+ay*fi(ii+1,jj+1))
        enddo
      enddo
!
      return
!
! ****** Error exit: invalid scales.
!
  900 continue
!
      write (*,*)
      write (*,*) '### ERROR in INTRP2D:'
      write (*,*) '### Scales are not monotonically increasing.'
      ierr=1
!
      return
!
! ****** Error exit: data outside range of scales.
!
  910 continue
!
      write (*,*)
      write (*,*) '### ERROR in INTRP2D:'
      write (*,*) '### An interpolation was attempted outside the'// &
                  ' range of the defined scales.'
      ierr=2
!
end subroutine
!#######################################################################
function flint (x,n,xn,fn,icheck,ierr)
!
!-----------------------------------------------------------------------
!
! ****** Interpolate a function linearly.
!
!-----------------------------------------------------------------------
!
! ****** The funcion is defined at N nodes, with values given by
! ****** FN(N) at positions XN(N).  The function value returned is
! ****** the linear interpolant at X.
!
! ****** Note that if X.lt.XN(1), the function value returned
! ****** is FN(1), and if X.gt.XN(N), the function value returned
! ****** is FN(N).
!
! ****** Call once with ICHECK.ne.0 to check that the values
! ****** in XN(N) are monotonically increasing.  In this mode
! ****** the array XN(N) is checked, and X and FN(N) are not
! ****** accessed.  If the check is passed, IERR=0 is returned.
! ****** Otherwise, IERR=1 is returned.
!
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      real(r_typ) :: flint
      real(r_typ) :: x
      integer :: n
      real(r_typ), dimension(n) :: xn,fn
      integer :: icheck
      integer :: ierr
!
!-----------------------------------------------------------------------
!
      integer :: i,j
      real(r_typ) :: x1,x2,alpha
!
!-----------------------------------------------------------------------
!
      ierr=0
      flint=0.
!
! ****** If ICHECK.ne.0, check the abscissa table.
!
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
!
! ****** Get the interpolated value.
!
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
!
      return
end function
!#######################################################################
subroutine interp (x,n,xv,i,alpha)
!
!-----------------------------------------------------------------------
!
! ****** Get interpolation factor ALPHA and table index I.
!
! ****** This routine does not do the actual interpolation.  Use the
! ****** returned values of I and ALPHA to get the interpolated value.
!
!-----------------------------------------------------------------------
!
      use number_types
!
!-----------------------------------------------------------------------
!
      implicit none
!
!-----------------------------------------------------------------------
!
      integer :: n
      real(r_typ), dimension(n) :: x
      real(r_typ) :: xv
      integer :: i
      real(r_typ) :: alpha
!
!-----------------------------------------------------------------------
!
      do 100 i=1,n-1
        if (xv.ge.x(i).and.xv.le.x(i+1)) then
          alpha=(xv-x(i))/(x(i+1)-x(i))
          go to 200
        end if
  100 continue
!
! ****** Value not found --- signal error and stop.
!
      write (*,*)
      write (*,*) '### ERROR in INTERP:'
      write (*,*) '### Value not found in table.'
      write (*,*) 'Value requested = ',xv
      write (*,*) 'Min table value = ',x(1)
      write (*,*) 'Max table value = ',x(n)
      call endrun (.true.)
!
  200 continue
!
end subroutine
!#######################################################################
!
! ****** Revision history:
!
! ### Version 1.00, 03/02/2006, file pot3d.f, modified by ZM:
!
!       - Cleaned up the previous version of POT3D.
!
! ### Version 1.01, 03/06/2006, file pot3d.f, modified by ZM:
!
!       - Added the ability to do a "source-surface plus
!         current-sheet" solution.  Select this by setting
!         OPTION='ss+cs'.
!
! ### Version 1.02, 06/18/2007, file pot3d.f, modified by ZM:
!
!       - Fixed a bug that caused the code to hang when an error
!         was encountered (when running a parallel job).
!
! ### Version 1.03, 03/17/2009, file pot3d.f, modified by ZM:
!
!       - Added the ability to write the boundary flux before the
!         sign flip for current-sheet solutions (i.e., OPTION='open').
!         Set the variable BR_PHOTO_ORIGINAL_FILE to the desired
!         file name to request this.
!
! ### Version 1.50, 01/25/2016, file pot3d.f, modified by RC:
!
!       - Added new (much faster) BILU0 preconditioner to CG solver.
!         To activate, set ifprec=2 in pot3d.dat file.
!       - Modified CG solve to use 1D arrays
!         for SAXPY and DOT operations.
!
! ### Version 2.00, 06/06/2017, file pot3d.f, modified by RC:
!
!       - Added OpenACC support to the code.
!         - OpenACC support is currently ONLY on 'standard'
!           pot3d runs (not inner-outer-iteratative mode)
!           and is only efficient on GPUs when using ifprec=1.
!         - OpenACC adds support for running the code on
!           Nvidia GPUs/Intel KNL/x86-multicore/OpenPower.
!         - To use OpenACC, simply compile the code with a compiler
!           that supports OpenACC with the correct flags activated.
!         - Multi-gpu support is included by using the new
!           ngpus_per_node input parameter.  This should be set
!           to the number of GPUs available per node.
!           The number of MPI ranks per node should match the
!           number of gpus per node.  This can be launched with
!           "mpiexec -np <np> -ntasks-per-node <ngpus_per_node>".
!         - The GPU features of the code are fully portable, i.e.
!           the code can be compiled/used as before on CPUs with no
!           changes in compilation or run-time.
!       - Modified some routines to be "nicer" for OpenACC
!         and optimized some MPI communications.
!       - Added wall-clock timer and corrected placement of
!         MPI_Finalize().  The wall clock timer now reflects the
!         true runtime.
!
! ### Version 2.01, 10/02/2017, file pot3d.f, modified by RC:
!
!       - Optimized OpenACC.
!       - Renamed cgsolv() to cgsolve().
!
! ### Version 2.10, 01/15/2018, file pot3d.f, modified by ZM+RC:
!
!       - Added the ability to skip the balancing of the flux
!         when doing a PFSS or OPEN field.  To invoke this, set
!         DO_NOT_BALNCE_FLUX=.true..
!       - Changed some pointers to allocatables for better
!         vectorization.
!
! ### Version 2.11, 03/19/2018, file pot3d.f, modified by RC:
!
!       - Added 'fmt' input parameter to set output file type.
!         Set fmt='h5' to write out HDF5 (default is 'hdf').
!
! ### Version 2.12, 10/08/2018, file pot3d.f, modified by RC:
!
!       - COMPATIBILITY CHANGE! Renamed gpus_per_node to gpn.
!         gpn is default 0 which will set gpn to the number of
!         MPI ranks in the local node.
!         Setting gpn manually is not recommended and only used for
!         oversubscribing the GPUs.
!       - Added MPI shared communicator for automatically setting gpn.
!         This requires an MPI-3 capable MPI library.
!       - Changed layout of matrix coefficient arrays to be more
!         vector-friendly instead of cache-friendly.
!
! ### Version 2.13, 11/19/2018, file pot3d.f, modified by RC:
!
!       - Small modifications to polar boundary condition calculations.
!       - Updated array ops and ACC clauses to be F2003 optimized.
!
! ### Version 2.20, 01/09/2019, file pot3d.f, modified by RC:
!
!       - Added double precision output option.
!         Set hdf32=.false. to activate 64-bit output.
!       - Updated magnetic field computation.  B is now computed
!         in parallel.  3D output fields now gathered to rank 0
!         using Sends and Receives instead of Gatherv in order
!         to allow very large resolutions.
!       - Added automatic topology.  Now, nprocs is optional.
!         One can specify one or more topology dimensions and
!         use the flag value "-1" to indicate dimensions to auto-set.
!         It is recommended to simply delete nprocs from input files.
!       - Added output file flushing so CG iterations can be monitored.
!       - Added new MPI rank diagnostics including
!         estimated load imbalance.
!       - Processor topology and magnetic energy output now written to
!         pot3d.out as well as the terminal.
!
! ### Version 2.21, 01/11/2019, file pot3d.f, modified by RC:
!
!       - Small updates to magnetic_energy routine.
!
! ### Version 2.22, 11/27/2019, file pot3d.f, modified by RC:
!
!       - Optimized some OpenACC directives.  Expanded some
!         array-syntax lines to full loops.
!
! ### Version 2.23, 08/11/2020, file pot3d.f, modified by RC:
!
!       - Small bug fix for default output file names and
!         format option fmt.
!
! ### Version 3.0.0, 02/10/2021, file pot3d.f, modified by RC:
!
!       - Changed version number scheme to semantic versioning.
!
! ### Version 3.1.0, 08/11/2020, file pot3d.f, modified by RC:
!
!       - Removed ss+cs inner-outer iterative mode.
!         It had not been used in many years and never worked
!         quite right.  Algorithm is saved in previous releases
!         on the github and svn.
!       - Added a validation run mode.  To activate, set:
!         VALIDATION_RUN=.TRUE.
!         This will set up a tilted dipole exact solution, write it
!         out, and then solve for it.  It will ignore input br files.
!         To change the angle of the dipole, set DIPOLE_ANGLE
!         (default is pi/4).
!       - Added NCGFLUSH parameter to flush output file during solver
!         iterations if NCGHIST>0 every NCGFLUSH iteration (default 25).
!         This is useful to track progress of the solver.
!
! ### Version 3.1.1, 11/02/2021, modified by RC:
!
!       - Fixed default output name mechanics.  The default names
!         were ignoring "fmt", making "fmt" almost useless.
!
! ### Version 3.2.0, 12/21/2021, modified by RC:
!
!       - GPN is now depricated.  Now, # of GPUs per node is assumed
!         to be equal to # of MPI ranks per node.
!         This removes the ability to oversubscribe GPUs, but this
!         was not being used ever, so no great loss.
!         This change was done to allow setting the gpu device
!         number before any allocations so that it works correctly
!         when using NVIDIA unified managed memory.
!
! ### Version 3.3.0, 02/14/2022, modified by RC:
!
!       - Added new timers for AX, PC, PCLOAD, and modified some timers
!         to be more useful.
!       - Updated ILU0 PC2 preconditioner to use single precision.
!         This speeds up the solve both on CPUs and GPUs.
!
! ### Version 4.0.0, 03/02/2022, modified by RC/MS:
!
!       - Added cuSparse library calls for ifprec=2 on NVIDIA GPUs.
!       - Added #ifdef
!         preprocessor directives which will either compile a cusparse
!         GPU version of the code or a CPU compatible version
!         without cusparse calls and no c code references
!
! ### Version 4.0.1, 03/31/2022, modified by RC:
!
!       - Added 'target' attribute to arrays used with C_LOC.
!         This fixes a compiler error with GFORTRAN when
!         using cuSparse.
!
! ### Version 4.1.0, 10/24/2022, modified by RC:
!
!       - Removed NCGFLUSH input parameter.  The output is
!         flushed at every write now.
!       - Changed behavior of solver residual reporting
!         to also output to terminal.
!       - Removed flush_output_file() routine in favor of Fortran 2003
!         intrinsic FLUSH.
!
! ### Version 4.2.0, 08/21/2024, modified by MS/RC:
!
!       - Replaced OpenACC with OpenMP Target directives.
!         There remains a single OpenACC directive for device selection
!         in order to work with nvfortran.
!
! ### Version 4.2.1, 09/19/2024, modified by RC:
!
!       - Changed array reduction loops to use
!         OpenMP reductions instead of atomics.
!       - Updated OpenMP target directives to use the
!         more streamlined "loop" syntax.
!
! ### Version 4.2.2, 12/05/2024, modified by RC:
!
!       - Fixed bug with OpenMP target loops.
!         Was using 'target loop' instead of 'target teams loop'
!         nvfortran seemed to allow this, but ifx does not.
!
! ### Version 4.3.0, 12/05/2024, modified by RC:
!
!       - Converted code to free-format and renamed file from
!         pot3d.F to pot3d.F90.
!
! ### Version 4.3.1, 12/05/2024, modified by RC:
!
!       - Code cleanup of modules.
!       - Made solver seam buffers global to avoid too many
!         allocates and deallocates for GPU.
!
!#######################################################################
