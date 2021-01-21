c
c-----------------------------------------------------------------------
c
c ****** Source to build the SDS library.
c ****** These routines are used by Zoran Mikic's tools.
c
c-----------------------------------------------------------------------
c
c **********************************************************************
c
c Copyright 2018 Predictive Science Inc.
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
c
c **********************************************************************
c
c        07/29/2003, ZM, Version 1.00:
c
c         - Original version of the SDS library.
c           This library was put together to facilitate the
c           development of ZM's tools.
c           It includes the new read/write routines for
c           scientific data sets (both text and HDF format).
c           The code was cleaned up to use standard FORTRAN90.
c
c        02/20/2004, ZM, Version 1.01:
c
c         - Added the ability to specify the format in writing
c           floating point numbers in routine WRFP.  This is used
c           in writing SDS text files using routine WRTXT.
c           For 32-bit data, WRTXT specifies 7 digits of precision,
c           whereas for 64-bit data, 15 digits of precision are
c           used.
c
c        04/02/2005, ZM, Version 1.02:
c
c         - Added a call to DFSDclear() in routine WRHDF.  This
c           initializes the SDS interface to the default state
c           for each new file.  This is needed to prevent settings
c           from previous files from interfering with each other.
c
c        06/16/2006, ZM, Version 1.03:
c
c         - Fixed some pointer allocation and manipulation issues in
c           the HDF read and write routines to make them obey
c           the FORTRAN 90 standard.  They were misbehaving on
c           the IFORT compiler.
c
c        02/24/2009, ZM, Version 1.04:
c
c         - Made a small change to the way an SDS is deallocated.
c         - Added a routine to initialize an SDS.  This is useful
c           when deallocating SDS structures.
c
c        08/30/2016, RC, Version 2.00:
c
c         - Added ability to read and write hdf5 files.
c           Now, rdhdf and wrhdf will read or write an hdf5 file
c           if the given fname ends in ".h5".
c           The library now needs to be linked to the hdf5 libraries.
c         - Modified rdhdf to be compatible with hdf4 files made using
c           the SD API instead of the DFSD API.
c
c        09/05/2016, RC, Version 2.01:
c
c         - Fixed problem with hdf5 writes when using 32-bit data.
c
c        05/22/2017, RC, Version 2.02:
c
c         - Fixed problem with 1D and 2D hdf5 writes when
c           the s%dims() are not set to 1 for the unused dimensions.
c
c        05/03/2019, RC, Version 2.03:
c
c         - Bug fix for HDF5 reads.
c
c-----------------------------------------------------------------------
c
c#######################################################################
      module sdslib_ident
c
      character(*), parameter :: cname='SDSLIB'
      character(*), parameter :: cvers='2.03'
      character(*), parameter :: cdate='05/03/2019'
c
      end module
c#######################################################################
      module assign_ptr_1d_interface
      interface
        subroutine assign_ptr_1d (from,to)
        use number_types
        implicit none
        real(r_typ), dimension(:), target :: from
        real(r_typ), dimension(:), pointer :: to
        end subroutine
      end interface
      end module
c#######################################################################
      module assign_ptr_3d_interface
      interface
        subroutine assign_ptr_3d (from,to)
        use number_types
        implicit none
        real(r_typ), dimension(:,:,:), target :: from
        real(r_typ), dimension(:,:,:), pointer :: to
        end subroutine
      end interface
      end module
c#######################################################################
      subroutine assign_ptr_1d (from,to)
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
      real(r_typ), dimension(:), target :: from
      real(r_typ), dimension(:), pointer :: to
c
c-----------------------------------------------------------------------
c
      to=>from
c
      return
      end subroutine
c#######################################################################
      subroutine assign_ptr_3d (from,to)
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
      real(r_typ), dimension(:,:,:), target :: from
      real(r_typ), dimension(:,:,:), pointer :: to
c
c-----------------------------------------------------------------------
c
      to=>from
c
      return
      end subroutine
c#######################################################################
      subroutine init_sds_pointer_status (s)
c
c-----------------------------------------------------------------------
c
c ****** Disassociate all the pointers in the SDS in structure S.
c
c-----------------------------------------------------------------------
c
c ****** This is useful when subsequently querying the association
c ****** status of these pointers (e.g., when deallocating storage).
c
c-----------------------------------------------------------------------
c
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
      nullify (s%f)
c
      nullify (s%scales(1)%f)
      nullify (s%scales(2)%f)
      nullify (s%scales(3)%f)
c
      return
      end subroutine
c#######################################################################
      subroutine deallocate_sds (s)
c
c-----------------------------------------------------------------------
c
c ****** Deallocate the memory used by the SDS in structure S.
c
c-----------------------------------------------------------------------
c
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
      if (associated(s%f)) deallocate (s%f)
c
      if (associated(s%scales(1)%f)) deallocate (s%scales(1)%f)
      if (associated(s%scales(2)%f)) deallocate (s%scales(2)%f)
      if (associated(s%scales(3)%f)) deallocate (s%scales(3)%f)
c
      return
      end subroutine
c#######################################################################
      subroutine rdhdf_1d (fname,scale,nx,f,x,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read a 1D scientific data set from an HDF file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine RDHDF to read the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx
      real(r_typ), dimension(:), pointer :: f
      real(r_typ), dimension(:), pointer :: x
      integer :: ierr
      intent(in) :: fname
      intent(out) :: scale,nx,ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Read the data set.
c
      call rdhdf (fname,s,ierr)
c
      if (ierr.ne.0) return
c
c ****** Check that this is a 1D data set.
c
      if (s%ndim.ne.1) then
        write (*,*)
        write (*,*) '### ERROR in RDHDF_1D:'
        write (*,*) '### The HDF file does not contain a 1D data set.'
        write (*,*) 'File name: ',trim(fname)
        ierr=3
        return
      end if
c
c ****** Set the output arguments.
c
      nx=s%dims(1)
      scale=s%scale
      x=>s%scales(1)%f
c
      allocate (f(nx))
      f=s%f(:,1,1)
      deallocate (s%f)
c
      return
      end
c#######################################################################
      subroutine rdhdf_2d (fname,scale,nx,ny,f,x,y,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read a 2D scientific data set from an HDF file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine RDHDF to read the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx,ny
      real(r_typ), dimension(:,:), pointer :: f
      real(r_typ), dimension(:), pointer :: x,y
      integer :: ierr
      intent(in) :: fname
      intent(out) :: scale,nx,ny,ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Read the data set.
c
      call rdhdf (fname,s,ierr)
c
      if (ierr.ne.0) return
c
c ****** Check that this is a 2D data set.
c
      if (s%ndim.ne.2) then
        write (*,*)
        write (*,*) '### ERROR in RDHDF_2D:'
        write (*,*) '### The HDF file does not contain a 2D data set.'
        write (*,*) 'File name: ',trim(fname)
        ierr=3
        return
      end if
c
c ****** Set the output arguments.
c
      nx=s%dims(1)
      ny=s%dims(2)
      scale=s%scale
      x=>s%scales(1)%f
      y=>s%scales(2)%f
c
      allocate (f(nx,ny))
      f(:,:)=s%f(:,:,1)
      deallocate (s%f)
c
      return
      end
c#######################################################################
      subroutine rdhdf_3d (fname,scale,nx,ny,nz,f,x,y,z,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read a 3D scientific data set from an HDF file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine RDHDF to read the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx,ny,nz
      real(r_typ), dimension(:,:,:), pointer :: f
      real(r_typ), dimension(:), pointer :: x,y,z
      integer :: ierr
      intent(in) :: fname
      intent(out) :: scale,nx,ny,nz,ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Read the data set.
c
      call rdhdf (fname,s,ierr)
c
      if (ierr.ne.0) return
c
c ****** Check that this is a 3D data set.
c
      if (s%ndim.ne.3) then
        write (*,*)
        write (*,*) '### ERROR in RDHDF_3D:'
        write (*,*) '### The HDF file does not contain a 3D data set.'
        write (*,*) 'File name: ',trim(fname)
        ierr=3
        return
      end if
c
c ****** Set the output arguments.
c
      nx=s%dims(1)
      ny=s%dims(2)
      nz=s%dims(3)
      scale=s%scale
      x=>s%scales(1)%f
      y=>s%scales(2)%f
      z=>s%scales(3)%f
      f=>s%f
c
      return
      end
c#######################################################################
      subroutine wrhdf_1d (fname,scale,nx,f,x,hdf32,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Write a 1D scientific data set to an HDF file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine WRHDF to write the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
      use assign_ptr_1d_interface
      use assign_ptr_3d_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx
      real(r_typ), dimension(nx,1,1) :: f
      real(r_typ), dimension(nx) :: x
      logical :: hdf32
      integer :: ierr
      intent(in) :: fname,scale,nx,f,x,hdf32
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Set the structure components.
c
      s%ndim=1
      s%dims(1)=nx
      s%dims(2)=1
      s%dims(3)=1
      s%scale=scale
      s%hdf32=hdf32
      if (scale) then
        call assign_ptr_1d (x,s%scales(1)%f)
      else
        nullify (s%scales(1)%f)
      end if
      nullify (s%scales(2)%f)
      nullify (s%scales(3)%f)
      call assign_ptr_3d (f,s%f)
c
c ****** Write the data set.
c
      call wrhdf (fname,s,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRHDF_1D:'
        write (*,*) '### Could not write the 1D data set.'
        write (*,*) 'File name: ',trim(fname)
        return
      end if
c
      return
      end
c#######################################################################
      subroutine wrhdf_2d (fname,scale,nx,ny,f,x,y,hdf32,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Write a 2D scientific data set to an HDF file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine WRHDF to write the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
      use assign_ptr_1d_interface
      use assign_ptr_3d_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx,ny
      real(r_typ), dimension(nx,ny,1) :: f
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      logical :: hdf32
      integer :: ierr
      intent(in) :: fname,scale,nx,ny,f,x,y,hdf32
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Set the structure components.
c
      s%ndim=2
      s%dims(1)=nx
      s%dims(2)=ny
      s%dims(3)=1
      s%scale=scale
      s%hdf32=hdf32
      if (scale) then
        call assign_ptr_1d (x,s%scales(1)%f)
        call assign_ptr_1d (y,s%scales(2)%f)
      else
        nullify (s%scales(1)%f)
        nullify (s%scales(2)%f)
      end if
      nullify (s%scales(3)%f)
      call assign_ptr_3d (f,s%f)
c
c ****** Write the data set.
c
      call wrhdf (fname,s,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRHDF_2D:'
        write (*,*) '### Could not write the 2D data set.'
        write (*,*) 'File name: ',trim(fname)
        return
      end if
c
      return
      end
c#######################################################################
      subroutine wrhdf_3d (fname,scale,nx,ny,nz,f,x,y,z,hdf32,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Write a 3D scientific data set to an HDF file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine WRHDF to write the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
      use assign_ptr_1d_interface
      use assign_ptr_3d_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx,ny,nz
      real(r_typ), dimension(nx,ny,nz) :: f
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      real(r_typ), dimension(nz) :: z
      logical :: hdf32
      integer :: ierr
      intent(in) :: fname,scale,nx,ny,nz,f,x,y,z,hdf32
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Set the structure components.
c
      s%ndim=3
      s%dims(1)=nx
      s%dims(2)=ny
      s%dims(3)=nz
      s%scale=scale
      s%hdf32=hdf32
      if (scale) then
        call assign_ptr_1d (x,s%scales(1)%f)
        call assign_ptr_1d (y,s%scales(2)%f)
        call assign_ptr_1d (z,s%scales(3)%f)
      else
        nullify (s%scales(1)%f)
        nullify (s%scales(2)%f)
        nullify (s%scales(3)%f)
      end if
      call assign_ptr_3d (f,s%f)
c
c ****** Write the data set.
c
      call wrhdf (fname,s,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRHDF_3D:'
        write (*,*) '### Could not write the 3D data set.'
        write (*,*) 'File name: ',trim(fname)
        return
      end if
c
      return
      end
c#######################################################################
      subroutine rdsds (fmt,fname,s,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read a scientific data set from file FNAME into
c ****** SDS structure S.
c
c ****** Use routine RDTXT or RDHDF, depending on the format
c ****** specified by FMT.
c
c-----------------------------------------------------------------------
c
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fmt
      character(*) :: fname
      type(sds) :: s
      integer :: ierr
      intent(in) :: fmt,fname
      intent(out) :: s,ierr
c
c-----------------------------------------------------------------------
c
      if (fmt.eq.'text') then
        call rdtxt (fname,s,ierr)
      else if (fmt.eq.'hdf'.or.fmt.eq.'h5') then
        call rdhdf (fname,s,ierr)
      else
        write (*,*)
        write (*,*) '### ERROR in RDSDS:'
        write (*,*) '### Invalid file format specified.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) 'Format = ',fmt
        ierr=5
        return
      end if
c
      return
      end
c#######################################################################
      subroutine wrsds (fmt,fname,s,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Write a scientific data set from SDS structure S to
c ****** file FNAME.
c
c ****** Use routine WRTXT or WRHDF, depending on the format
c ****** specified by FMT.
c
c-----------------------------------------------------------------------
c
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fmt
      character(*) :: fname
      type(sds) :: s
      integer :: ierr
      intent(in) :: fmt,fname,s
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
      if (fmt.eq.'text') then
        call wrtxt (fname,s,ierr)
      else if (fmt.eq.'hdf'.or.fmt.eq.'h5') then
        call wrhdf (fname,s,ierr)
      else
        write (*,*)
        write (*,*) '### ERROR in WRSDS:'
        write (*,*) '### Invalid file format specified.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) 'Format = ',fmt
        ierr=5
        return
      end if
c
      return
      end
c#######################################################################
      subroutine rdhdf (fname,s,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read a 1D, 2D, or 3D scientific data set from an HDF file.
c ****** This routine uses the new SD API instead of the
c ****** outdated DFSD API.
c
c-----------------------------------------------------------------------
c
c ****** This routine allocates the required memory and returns
c ****** pointers to the data and scale arrays.
c
c-----------------------------------------------------------------------
c
c ****** Input arguments:
c
c          FNAME   : [character(*)]
c                    HDF data file name to read from.
c
c ****** Output arguments:
c
c          S       : [structure of type SDS]
c                    A structure that holds the field, its
c                    dimensions, and the scales, with the
c                    components described below.
c
c          IERR    : [integer]
c                    IERR=0 is returned if the data set was read
c                    successfully.  Otherwise, IERR is set to a
c                    nonzero value.
c
c ****** Components of structure S:
c
c          NDIM    : [integer]
c                    Number of dimensions found in the data set.
c
c          DIMS    : [integer, dimension(3)]
c                    Number of points in the data set dimensions.
c                    For a 1D data set, DIMS(2)=DIMS(3)=1.
c                    For a 2D data set, DIMS(3)=1.
c
c          SCALE   : [logical]
c                    Flag to indicate the presence of scales (axes)
c                    in the data set.  SCALE=.false. means that scales
c                    were not found; SCALE=.true. means that scales
c                    were found.
c
c          HDF32   : [logical]
c                    Flag to indicate the precision of the data set
c                    read in.  HDF32=.true. means that the data is
c                    32-bit; HDF32=.false. means that the data is
c                    64-bit.
c
c          SCALES  : [structure of type RP1D, dimension(3)]
c                    This array holds the pointers to the scales
c                    when SCALE=.true., and is undefined otherwise.
c
c          F       : [real, pointer to a rank-3 array]
c                    This array holds the data set values.
c
c ****** The storage for the arrays pointed to by F, and the
c ****** scales (if present) in structure SCALES, is allocated by
c ****** this routine.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      character, dimension(256) :: sds_name, dim_name
      type(sds) :: s
      integer :: ierr,i
      intent(in) :: fname
      intent(out) :: s,ierr
c
c-----------------------------------------------------------------------
c
      ierr=0
c
c-----------------------------------------------------------------------
c
c ****** Read hdf5 file if fname ends in '.h5'.
c
      i=index(fname,'.h');
      if (fname(i+1:i+2).eq.'h5') then
        call rdh5 (fname,s,ierr)
        return
      else
        print*,"HDF4 has been disabled"
        ierr=-1
      end if
c
      return
      end
c#######################################################################
      subroutine rdh5 (fname,s,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read a 1D, 2D, or 3D scientific data set from an HDF5 file.
c ****** The HDF5 file is currently assumed to contain only one
c ****** dataset (1D,2d,3D), with or without scales, in group "/",
c ****** and has no other data members.
c
c-----------------------------------------------------------------------
c
c ****** This routine allocates the required memory and returns
c ****** pointers to the data and scale arrays.
c
c-----------------------------------------------------------------------
c
c ****** Input arguments:
c
c          FNAME   : [character(*)]
c                    HDF5 data file name to read from.
c
c ****** Output arguments:
c
c          S       : [structure of type SDS]
c                    A structure that holds the field, its
c                    dimensions, and the scales, with the
c                    components described below.
c
c          IERR    : [integer]
c                    IERR=0 is returned if the data set was read
c                    successfully.  Otherwise, IERR is set to a
c                    nonzero value.
c
c ****** Components of structure S:
c
c          NDIM    : [integer]
c                    Number of dimensions found in the data set.
c
c          DIMS    : [integer, dimension(3)]
c                    Number of points in the data set dimensions.
c                    For a 1D data set, DIMS(2)=DIMS(3)=1.
c                    For a 2D data set, DIMS(3)=1.
c
c          SCALE   : [logical]
c                    Flag to indicate the presence of scales (axes)
c                    in the data set.  SCALE=.false. means that scales
c                    were not found; SCALE=.true. means that scales
c                    were found.
c
c          HDF32   : [logical]
c                    Flag to indicate the precision of the data set
c                    read in.  HDF32=.true. means that the data is
c                    32-bit; HDF32=.false. means that the data is
c                    64-bit.
c
c          SCALES  : [structure of type RP1D, dimension(3)]
c                    This array holds the pointers to the scales
c                    when SCALE=.true., and is undefined otherwise.
c
c          F       : [real, pointer to a rank-3 array]
c                    This array holds the data set values.
c
c ****** The storage for the arrays pointed to by F, and the
c ****** scales (if present) in structure SCALES, is allocated by
c ****** this routine.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
      use hdf5
      use h5ds
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      type(sds) :: s
      character(*) :: fname
c
c-----------------------------------------------------------------------
c
      integer :: ierr
c
c-----------------------------------------------------------------------
c
      integer :: i,obj_type,n_members
c
      integer(HID_T) :: file_id       ! File identifier
      integer(HID_T) :: dset_id       ! Dataset identifier
      integer(HID_T) :: dspace_id     ! Dataspace identifier
      integer(HID_T) :: dim_id        ! Dimension identifiers
      integer(HID_T) :: datatype_id   ! Datatype identifiers
c
      integer(SIZE_T) :: prec
c
      integer(HSIZE_T),dimension(:), allocatable :: s_dims,maxpts
      integer(HSIZE_T),dimension(1) :: s_dims_i
c
      real(KIND_REAL_4), dimension(:,:,:), allocatable :: f4
      real(KIND_REAL_4), dimension(:),     allocatable :: f4dim
      real(KIND_REAL_8), dimension(:,:,:), allocatable :: f8
      real(KIND_REAL_8), dimension(:),     allocatable :: f8dim
c
      character(512) :: obj_name
      character(4), parameter :: cname='RDH5'
c
      logical :: is_scale
c
c-----------------------------------------------------------------------
c
c ****** Initialize dimension count and arrays.
c
      s%ndim=0
      s%dims(:)=1
c
c ****** Initialize hdf5 interface.
c
      call h5open_f (ierr)
c
c ****** Open hdf5 file.
c
      call h5Fopen_f (trim(fname),H5F_ACC_RDONLY_F,file_id,ierr)
c
c ****** Get information about the hdf5 file.
c
      call h5Gn_members_f (file_id,"/",n_members,ierr)
c
c ****** Make sure there is (at maximum) one 3D dataset with scales.
c
      if (n_members.eq.0.or.n_members.gt.4) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Input file contains too few/many datasets.'
        write (*,*) 'File name: ',trim(fname)
        return
      endif
c
c ****** Assume the Dataset is in index 0 and get its name.
c
      call h5Gget_obj_info_idx_f (file_id,"/",0,obj_name,obj_type,ierr)
c
c ****** Open Dataset.
c
      call h5Dopen_f (file_id,trim(obj_name),dset_id,ierr)
c
c ****** Make sure the Dataset is not a scale.
c
      call h5DSis_scale_f(dset_id,is_scale,ierr)
      if (is_scale) then
        write (*,*)
        write (*,*) '### ERROR in ',cname,':'
        write (*,*) '### Input file Dataset at index 0 is a scale.'
        write (*,*) 'File name: ',trim(fname)
        return
      endif
c
c ****** Get dimensions (need s_dims array for format requirements).
c
      call h5Dget_space_f (dset_id,dspace_id,ierr)
      call h5Sget_simple_extent_ndims_f (dspace_id,s%ndim,ierr)
c
      allocate(s_dims(s%ndim))
c
      allocate(maxpts(s%ndim))
      call h5Sget_simple_extent_dims_f (dspace_id,s_dims,maxpts,ierr)
      deallocate(maxpts)
c
      s%dims(1:s%ndim)=s_dims(1:s%ndim)
c
c ****** Get the floating-point precision of the data and set flag.
c
      call h5Dget_type_f (dset_id,datatype_id,ierr)
      call h5Tget_precision_f (datatype_id,prec,ierr)
c
      if (prec.eq.32) then
        s%hdf32=.true.
      elseif (prec.eq.64) then
        s%hdf32=.false.
      end if
c
c ****** Allocate the memory for the Dataset array in s.
c
      allocate (s%f(s%dims(1),s%dims(2),s%dims(3)))
c
c ****** Need to read the file in its own datatype, and then convert
c ****** to datatype of s%f.
c
      if (s%hdf32) then
        allocate (f4(s%dims(1),s%dims(2),s%dims(3)))
        call h5Dread_f (dset_id,datatype_id,f4,s_dims,ierr)
        s%f(:,:,:)=f4(:,:,:)
        deallocate (f4)
      else
        allocate (f8(s%dims(1),s%dims(2),s%dims(3)))
        call h5Dread_f (dset_id,datatype_id,f8,s_dims,ierr)
        s%f(:,:,:)=f8(:,:,:)
        deallocate (f8)
      end if
c
      deallocate(s_dims)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in RDH5:'
        write (*,*) '### Error while reading the dataset.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) '[Error return (from H5DREAD_F) = ',ierr,']'
        ierr=4
        return
      end if
c
c ****** Close the hdf5 type descriptor.
c
      call h5Tclose_f (datatype_id,ierr)
c
c ****** Check if there might be scales present, if so, read them.
c
      if (n_members.gt.1) then
c
c ***** First check that the number of scale datasets match the # dim.
c
        if (n_members-1.ne.s%ndim) then
          write (*,*)
          write (*,*) '### ERROR in RDH5:'
          write (*,*) '### # scales does not match # dims.'
          write (*,*) 'File name: ',trim(fname)
          return
        end if
c
        s%scale=.true.
c
c ****** Loop through scales, make sure each is a scale, and read them.
c
        do i=1,n_members-1
c
c ****** Get the name of scale dataset.
c
          call h5Gget_obj_info_idx_f (file_id,"/",i,
     &                               obj_name,obj_type,ierr)
c
c ****** Open scale dataset.
c
          call h5Dopen_f (file_id,trim(obj_name),dim_id,ierr)
c
c ****** Make sure the scale is a scale.
c
          call h5DSis_scale_f (dim_id,is_scale,ierr)
          if (.not.is_scale) then
            write (*,*)
            write (*,*) '### ERROR in RDH5:'
            write (*,*) '### Scale is not a scale.'
            write (*,*) 'File name: ',trim(fname)
            return
          end if
c
c ****** Get dimension of scale.
c
          s_dims_i=s%dims(i)
c
c ****** Allocate scale.
c
          allocate (s%scales(i)%f(s_dims_i(1)))
c
c ****** Get the floating-point precision of the scale.
c
          call h5Dget_type_f (dim_id,datatype_id,ierr)
          call h5Tget_precision_f (datatype_id,prec,ierr)
c
c ****** Read in the scale data.
c
          if (s%hdf32) then
            allocate (f4dim(s_dims_i(1)))
            call h5Dread_f (dim_id,datatype_id,f4dim,s_dims_i,ierr)
            s%scales(i)%f(:)=f4dim(:)
            deallocate (f4dim)
          else
            allocate (f8dim(s_dims_i(1)))
            call h5Dread_f (dim_id,datatype_id,f8dim,s_dims_i,ierr)
            s%scales(i)%f(:)=f8dim(:)
            deallocate (f8dim)
          end if
c
c ****** Close the scale dataset.
c
          call h5Dclose_f (dim_id,ierr)
c
        enddo
c
c ****** Allocate dummy scales (of length 1) for empty dimensions.
c
        do i=s%ndim+1,3
          allocate (s%scales(i)%f(1))
        enddo
      else
c
c ****** If scales are not present, allocate dummy
c ****** scales (of length 1) so that the pointers to the scales
c ****** are valid.
c
        s%scale=.false.
c
        allocate (s%scales(1)%f(1))
        allocate (s%scales(2)%f(1))
        allocate (s%scales(3)%f(1))
      end if
c
c ****** Close the dataset.
c
      call h5Dclose_f (dset_id,ierr)
c
c ****** Close the file.
c
      call h5Fclose_f (file_id,ierr)
c
c ****** Close FORTRAN interface.
c
      call h5close_f (ierr)
c
      return
      end subroutine
c#######################################################################
      subroutine wrhdf (fname,s,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Write a 1D, 2D, or 3D scientific data set to an HDF file.
c
c-----------------------------------------------------------------------
c
c ****** Input arguments:
c
c          FNAME   : [character(*)]
c                    HDF data file name to write to.
c
c          S       : [structure of type SDS]
c                    A structure that holds the field, its
c                    dimensions, and the scales, with the
c                    components described below.
c
c ****** Output arguments:
c
c          IERR    : [integer]
c                    IERR=0 is returned if the data set was written
c                    successfully.  Otherwise, IERR is set to a
c                    nonzero value.
c
c ****** Components of structure S:
c
c          NDIM    : [integer]
c                    Number of dimensions in the data set.
c
c          DIMS    : [integer, dimension(3)]
c                    Number of points in the data set dimensions.
c                    Only DIMS(1 .. NDIM) are referenced.
c
c          SCALE   : [logical]
c                    Flag to indicate the presence of scales (axes)
c                    in the data set.  SCALE=.false. means that scales
c                    are not being supplied; SCALE=.true. means that
c                    scales are being supplied.
c
c          HDF32   : [logical]
c                    Flag to specify the precision of the data to
c                    be written to the file.  Set HDF32=.true. to
c                    write 32-bit data, and HDF32=.false. to write
c                    64-bit data.
c
c          SCALES  : [structure of type RP1D, dimension(3)]
c                    This array holds the pointers to the scales
c                    when SCALE=.true., and is not referenced
c                    otherwise.
c
c          F       : [real, pointer to a rank-3 array]
c                    This array holds the data set values.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      type(sds) :: s
      integer :: ierr
      intent(in) :: fname,s
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
      integer :: iret,i,n
c
c-----------------------------------------------------------------------
c
      ierr=0
c
c ****** Check the number of dimensions.
c
      if (s%ndim.le.0.or.s%ndim.gt.3) then
        write (*,*)
        write (*,*) '### ERROR in WRHDF:'
        write (*,*) '### Could not write the SDS data.'
        write (*,*) 'Invalid number of dimensions.'
        write (*,*) 'Number of dimensions = ',s%ndim
        write (*,*) 'File name: ',trim(fname)
        ierr=1
        return
      end if
c
c ****** Write hdf5 file if fname ends in '.h5'.
c
      i=index(fname,'.h')
      if (fname(i+1:i+2).eq.'h5') then
        call wrh5 (fname,s,ierr)
      else
        print*,"HDF4 disabled."
        ierr=-1
      end if
c
      return
      end
c#######################################################################
      subroutine wrh5 (fname,s,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Write a 1D, 2D, or 3D scientific data set to an HDF5 file.
c
c-----------------------------------------------------------------------
c
c ****** Input arguments:
c
c          FNAME   : [character(*)]
c                    HDF data file name to write to.
c
c          S       : [structure of type SDS]
c                    A structure that holds the field, its
c                    dimensions, and the scales, with the
c                    components described below.
c
c ****** Output arguments:
c
c          IERR    : [integer]
c                    IERR=0 is returned if the data set was written
c                    successfully.  Otherwise, IERR is set to a
c                    nonzero value.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
      use hdf5
      use h5ds
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      type(sds) :: s
      integer :: ierr
      intent(in) :: fname,s
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
      character(8) ::   dimname
      integer :: i
      integer(HID_T) :: file_id       ! File identifier
      integer(HID_T) :: dset_id       ! Dataset identifier
      integer(HID_T) :: dspace_id,dspacedim_id   ! Dataspace identifiers
      integer(HID_T) :: dim_id        ! Dimension identifiers
      integer(HSIZE_T),dimension(3) :: s_dims
      integer(HSIZE_T),dimension(1) :: s_dims_i
c
      real(KIND_REAL_4), dimension(:,:,:), allocatable :: f4
      real(KIND_REAL_4), dimension(:),     allocatable :: f4dim
      real(KIND_REAL_8), dimension(:,:,:), allocatable :: f8
      real(KIND_REAL_8), dimension(:),     allocatable :: f8dim
c
c-----------------------------------------------------------------------
c
c ****** HDF5 calls are picky about the integer format for the dims
c ****** so the s%dims need to be converted to HSIZE_T integers.
c
c ****** Also, sometimes calls to wrhdf() for 1D and 2D datasets
c ****** do not have the unused dims(i) set to 1 (unset).
c ****** To avoid needing a function call to implicitly reshape
c ****** f(n), set the dims here.
c
      do i=1,3
         if (i.le.s%ndim) then
           s_dims(i)=s%dims(i)
         else
           s_dims(i)=1
         endif
      end do
c
c ****** Initialize hdf5 interface.
c
      call h5open_f (ierr)
c
c ****** Create the file.
c
      call h5Fcreate_f (trim(fname),H5F_ACC_TRUNC_F,file_id,ierr)
c
c ****** Create the dataspace.
c
      call h5Screate_simple_f (s%ndim,s_dims,dspace_id,ierr)
c
c ****** Create and write the dataset (convert s%f to proper type).
c
      if (s%hdf32) then
        allocate (f4(s_dims(1),s_dims(2),s_dims(3)))
        f4(:,:,:)=s%f(:,:,:)
        call h5Dcreate_f (file_id,'Data',H5T_NATIVE_REAL,
     &                    dspace_id,dset_id,ierr)
        call h5Dwrite_f (dset_id,H5T_NATIVE_REAL,f4,s_dims,ierr)
        deallocate (f4)
      else
        allocate (f8(s_dims(1),s_dims(2),s_dims(3)))
        f8(:,:,:)=s%f(:,:,:)
        call h5Dcreate_f (file_id,'Data',H5T_NATIVE_DOUBLE,
     &                    dspace_id,dset_id,ierr)
        call h5Dwrite_f (dset_id,H5T_NATIVE_DOUBLE,f8,s_dims,ierr)
        deallocate (f8)
      endif
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRH5:'
        write (*,*) '### Could not write the dataset.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) '[Error return (from h5Dwrite_f) = ',ierr,']'
        ierr=4
        return
      end if
c
c ****** Check for scales.  If present, add them to the hdf5 dataset.
c
      if (s%scale) then
        do i=1,s%ndim
          if (i.eq.1) then
            dimname='dim1'
          elseif (i.eq.2) then
            dimname='dim2'
          elseif (i.eq.3) then
            dimname='dim3'
          endif
          s_dims_i=s_dims(i)
          call h5Screate_simple_f(1,s_dims_i,dspacedim_id,ierr)
          if (s%hdf32) then
            allocate (f4dim(s_dims_i(1)))
            f4dim(:)=s%scales(i)%f(:)
            call h5Dcreate_f (file_id,dimname,H5T_NATIVE_REAL,
     &                        dspacedim_id,dim_id,ierr)
            call h5Dwrite_f (dim_id,H5T_NATIVE_REAL,
     &                       f4dim,s_dims_i,ierr)
            deallocate (f4dim)
          else
            allocate (f8dim(s_dims_i(1)))
            f8dim(:)=s%scales(i)%f(:)
            call h5Dcreate_f (file_id,dimname,H5T_NATIVE_DOUBLE,
     &                        dspacedim_id,dim_id,ierr)
            call h5Dwrite_f (dim_id,H5T_NATIVE_DOUBLE,
     &                       f8dim,s_dims_i,ierr)
            deallocate (f8dim)
          endif
          call h5DSset_scale_f (dim_id,ierr,dimname)
          call h5DSattach_scale_f (dset_id,dim_id,i,ierr)
          call h5DSset_label_f(dset_id, i, dimname, ierr)
          call h5Dclose_f (dim_id,ierr)
          call h5Sclose_f (dspacedim_id,ierr)
        end do
        if (ierr.ne.0) then
          write (*,*)
          write (*,*) '### ERROR in WRH5:'
          write (*,*) '### Could not write the scales.'
          write (*,*) 'File name: ',trim(fname)
          ierr=5
          return
        endif
      endif
c
c ****** Close the dataset.
c
      call h5Dclose_f (dset_id,ierr)
c
c ****** Close the dataspace.
c
      call h5Sclose_f (dspace_id,ierr)
c
c ****** Close the file.
c
      call h5Fclose_f (file_id,ierr)
c
c ****** Close the hdf5 interface.
c
      call h5close_f (ierr)
c
      end subroutine
c#######################################################################
      subroutine rdtxt_1d (fname,scale,nx,f,x,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read a 1D scientific data set from a text file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine RDTXT to read the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx
      real(r_typ), dimension(:), pointer :: f
      real(r_typ), dimension(:), pointer :: x
      integer :: ierr
      intent(in) :: fname
      intent(out) :: scale,nx,ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Read the data set.
c
      call rdtxt (fname,s,ierr)
c
      if (ierr.ne.0) return
c
c ****** Check that this is a 1D data set.
c
      if (s%ndim.ne.1) then
        write (*,*)
        write (*,*) '### ERROR in RDTXT_1D:'
        write (*,*) '### The test file does not contain a 1D data set.'
        write (*,*) 'File name: ',trim(fname)
        ierr=3
        return
      end if
c
c ****** Set the output arguments.
c
      nx=s%dims(1)
      scale=s%scale
      x=>s%scales(1)%f
c
      allocate (f(nx))
      f=s%f(:,1,1)
      deallocate (s%f)
c
      return
      end
c#######################################################################
      subroutine rdtxt_2d (fname,scale,nx,ny,f,x,y,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read a 2D scientific data set from a text file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine RDTXT to read the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx,ny
      real(r_typ), dimension(:,:), pointer :: f
      real(r_typ), dimension(:), pointer :: x,y
      integer :: ierr
      intent(in) :: fname
      intent(out) :: scale,nx,ny,ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Read the data set.
c
      call rdtxt (fname,s,ierr)
c
      if (ierr.ne.0) return
c
c ****** Check that this is a 2D data set.
c
      if (s%ndim.ne.2) then
        write (*,*)
        write (*,*) '### ERROR in RDTXT_2D:'
        write (*,*) '### The text file does not contain a 2D data set.'
        write (*,*) 'File name: ',trim(fname)
        ierr=3
        return
      end if
c
c ****** Set the output arguments.
c
      nx=s%dims(1)
      ny=s%dims(2)
      scale=s%scale
      x=>s%scales(1)%f
      y=>s%scales(2)%f
c
      allocate (f(nx,ny))
      f=s%f(:,:,1)
      deallocate (s%f)
c
      return
      end
c#######################################################################
      subroutine rdtxt_3d (fname,scale,nx,ny,nz,f,x,y,z,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read a 3D scientific data set from a text file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine RDTXT to read the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx,ny,nz
      real(r_typ), dimension(:,:,:), pointer :: f
      real(r_typ), dimension(:), pointer :: x,y,z
      integer :: ierr
      intent(in) :: fname
      intent(out) :: scale,nx,ny,nz,ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Read the data set.
c
      call rdtxt (fname,s,ierr)
c
      if (ierr.ne.0) return
c
c ****** Check that this is a 3D data set.
c
      if (s%ndim.ne.3) then
        write (*,*)
        write (*,*) '### ERROR in RDTXT_3D:'
        write (*,*) '### The text file does not contain a 3D data set.'
        write (*,*) 'File name: ',trim(fname)
        ierr=3
        return
      end if
c
c ****** Set the output arguments.
c
      nx=s%dims(1)
      ny=s%dims(2)
      nz=s%dims(3)
      scale=s%scale
      x=>s%scales(1)%f
      y=>s%scales(2)%f
      z=>s%scales(3)%f
      f=>s%f
c
      return
      end
c#######################################################################
      subroutine wrtxt_1d (fname,scale,nx,f,x,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Write a 1D scientific data set to a text file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine WRTXT to write the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
      use assign_ptr_1d_interface
      use assign_ptr_3d_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx
      real(r_typ), dimension(nx,1,1) :: f
      real(r_typ), dimension(nx) :: x
      integer :: ierr
      intent(in) :: fname,scale,nx,f,x
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Set the structure components.
c
      s%ndim=1
      s%dims(1)=nx
      s%dims(2)=1
      s%dims(3)=1
      s%scale=scale
      if (scale) then
        call assign_ptr_1d (x,s%scales(1)%f)
      else
        nullify (s%scales(1)%f)
      end if
      nullify (s%scales(2)%f)
      nullify (s%scales(3)%f)
      call assign_ptr_3d (f,s%f)
c
c ****** Write the data set.
c
      call wrtxt (fname,s,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRTXT_1D:'
        write (*,*) '### Could not write the 1D data set.'
        write (*,*) 'File name: ',trim(fname)
        return
      end if
c
      return
      end
c#######################################################################
      subroutine wrtxt_2d (fname,scale,nx,ny,f,x,y,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Write a 2D scientific data set to a text file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine WRTXT to write the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
      use assign_ptr_1d_interface
      use assign_ptr_3d_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx,ny
      real(r_typ), dimension(nx,ny,1) :: f
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      integer :: ierr
      intent(in) :: fname,scale,nx,ny,f,x,y
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Set the structure components.
c
      s%ndim=2
      s%dims(1)=nx
      s%dims(2)=ny
      s%dims(3)=1
      s%scale=scale
      if (scale) then
        call assign_ptr_1d (x,s%scales(1)%f)
        call assign_ptr_1d (y,s%scales(2)%f)
      else
        nullify (s%scales(1)%f)
        nullify (s%scales(2)%f)
      end if
      nullify (s%scales(3)%f)
      call assign_ptr_3d (f,s%f)
c
c ****** Write the data set.
c
      call wrtxt (fname,s,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRTXT_2D:'
        write (*,*) '### Could not write the 2D data set.'
        write (*,*) 'File name: ',trim(fname)
        return
      end if
c
      return
      end
c#######################################################################
      subroutine wrtxt_3d (fname,scale,nx,ny,nz,f,x,y,z,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Write a 3D scientific data set to a text file.
c
c-----------------------------------------------------------------------
c
c ****** This routine calls routine WRTXT to write the file.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
      use assign_ptr_1d_interface
      use assign_ptr_3d_interface
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      logical :: scale
      integer :: nx,ny,nz
      real(r_typ), dimension(nx,ny,nz) :: f
      real(r_typ), dimension(nx) :: x
      real(r_typ), dimension(ny) :: y
      real(r_typ), dimension(nz) :: z
      integer :: ierr
      intent(in) :: fname,scale,nx,ny,nz,f,x,y,z
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
c ****** Declaration for the SDS structure.
c
      type(sds) :: s
c
c-----------------------------------------------------------------------
c
c ****** Set the structure components.
c
      s%ndim=3
      s%dims(1)=nx
      s%dims(2)=ny
      s%dims(3)=nz
      s%scale=scale
      if (scale) then
        call assign_ptr_1d (x,s%scales(1)%f)
        call assign_ptr_1d (y,s%scales(2)%f)
        call assign_ptr_1d (z,s%scales(3)%f)
      else
        nullify (s%scales(1)%f)
        nullify (s%scales(2)%f)
        nullify (s%scales(3)%f)
      end if
      call assign_ptr_3d (f,s%f)
c
c ****** Write the data set.
c
      call wrtxt (fname,s,ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRTXT_3D:'
        write (*,*) '### Could not write the 3D data set.'
        write (*,*) 'File name: ',trim(fname)
        return
      end if
c
      return
      end
c#######################################################################
      subroutine rdtxt (fname,s,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read a 1D, 2D, or 3D scientific data set from a text file.
c
c-----------------------------------------------------------------------
c
c ****** This routine allocates the required memory and returns
c ****** pointers to the data and scale arrays.
c
c-----------------------------------------------------------------------
c
c ****** Input arguments:
c
c          FNAME   : [character(*)]
c                    Text data file name to read from.
c
c ****** Output arguments:
c
c          S       : [structure of type SDS]
c                    A structure that holds the field, its
c                    dimensions, and the scales, with the
c                    components described below.
c
c          IERR    : [integer]
c                    IERR=0 is returned if the data set was read
c                    successfully.  Otherwise, IERR is set to a
c                    nonzero value.
c
c ****** Components of structure S:
c
c          NDIM    : [integer]
c                    Number of dimensions found in the data set.
c
c          DIMS    : [integer, dimension(3)]
c                    Number of points in the data set dimensions.
c                    For a 1D data set, DIMS(2)=DIMS(3)=1.
c                    For a 2D data set, DIMS(3)=1.
c
c          SCALE   : [logical]
c                    Flag to indicate the presence of scales (axes)
c                    in the data set.  SCALE=.false. means that scales
c                    were not found; SCALE=.true. means that scales
c                    were found.
c
c          HDF32   : [logical]
c                    This flag is is not relevant to text files;
c                    it is used for HDF data.
c                    It is arbitrarily set to HDF32=.false..
c
c          SCALES  : [structure of type RP1D, dimension(3)]
c                    This array holds the pointers to the scales
c                    when SCALE=.true., and is undefined otherwise.
c
c          F       : [real, pointer to a rank-3 array]
c                    This array holds the data set values.
c
c ****** The storage for the arrays pointed to by F, and the
c ****** scales (if present) in structure SCALES, is allocated by
c ****** this routine.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      type(sds) :: s
      integer :: ierr
      intent(in) :: fname
      intent(out) :: s,ierr
c
c-----------------------------------------------------------------------
c
      integer, parameter :: mxdim=3
      integer :: ifscale,i,n
      integer :: int_tmp(1)
c
c-----------------------------------------------------------------------
c
      ierr=0
c
c ****** Open the file for reading.
c
      call ffopen (1,fname,'r',ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in RDTXT:'
        write (*,*) '### Could not open the text file.'
        write (*,*) 'File name: ',trim(fname)
        ierr=1
        return
      end if
c
c ****** Get the number of dimensions.
c
      int_tmp(1)=s%ndim
      call rdint (1,1,int_tmp,ierr)
      if (ierr.ne.0) go to 910
c
      if (s%ndim.lt.1.or.s%ndim.gt.mxdim) then
        write (*,*)
        write (*,*) '### ERROR in RDTXT:'
        write (*,*) '### Invalid number of dimensions in file.'
        write (*,*) 'File name: ',trim(fname)
        write (*,*) 'Number of dimensions = ',s%ndim
        write (*,*) 'Maximum number of dimensions = ',mxdim
        ierr=2
        return
      end if
c
c ****** Read the dimensions.
c
      s%dims(:)=1
c
      do i=1,s%ndim
        call rdint (1,1,s%dims(i),ierr)
        if (ierr.ne.0) go to 910
        if (s%dims(i).le.0) go to 920
      enddo
c
c ****** Check if the scales are present.
c
      int_tmp(1)=ifscale
      call rdint (1,1,int_tmp,ierr)
      if (ierr.ne.0) go to 910
c
      s%scale=ifscale.ne.0
c
c ****** Allocate memory and read the scales (if present).
c
c ****** If scales are not present,  allocate dummy scales
c ****** (of length 1) so that the pointers to the scales
c ****** are valid.
c
      if (s%scale) then
        do i=1,s%ndim
          allocate (s%scales(i)%f(s%dims(i)))
          call rdfp (1,s%dims(i),s%scales(i)%f,ierr)
          if (ierr.ne.0) go to 910
        enddo
        do i=s%ndim+1,3
          allocate (s%scales(i)%f(1))
        enddo
      else
        allocate (s%scales(1)%f(1))
        allocate (s%scales(2)%f(1))
        allocate (s%scales(3)%f(1))
      end if
c
c ****** Allocate memory for the array.
c
      allocate (s%f(s%dims(1),s%dims(2),s%dims(3)))
c
c ****** Read the data array.
c
      n=product(s%dims(1:s%ndim))
      call rdfp (1,n,s%f,ierr)
      if (ierr.ne.0) go to 910
c
      s%hdf32=.false.
c
      close (1)
c
      return
c
  910 continue
c
      write (*,*)
      write (*,*) '### ERROR in RDTXT:'
      write (*,*) '### Error while reading text data.'
      write (*,*) 'File name: ',trim(fname)
      ierr=3
      return
c
  920 continue
c
      write (*,*)
      write (*,*) '### ERROR in RDTXT:'
      write (*,*) '### Invalid value for dimension.'
      write (*,*) 'Dimension number = ',i
      write (*,*) 'Dimension value = ',s%dims(i)
      ierr=4
c
      return
      end
c#######################################################################
      subroutine rdint (iun,n,i,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read N words of INTEGER data into array I from unit IUN
c ****** using a free format text read.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: iun
      integer :: n
      integer, dimension(n) :: i
      integer :: ierr
      intent(in) :: iun,n
      intent(out) :: i,ierr
c
c-----------------------------------------------------------------------
c
      ierr=0
c
      read (iun,*,err=100,end=100) i
c
      return
c
  100 continue
c
c ****** Error in reading the data.
c
      ierr=1
c
      return
      end
c#######################################################################
      subroutine rdfp (iun,n,f,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Read N words of REAL data into array F from unit IUN
c ****** using a free format text read.
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
      integer :: iun
      integer :: n
      real(r_typ), dimension(n) :: f
      integer :: ierr
      intent(in) :: iun,n
      intent(out) :: f,ierr
c
c-----------------------------------------------------------------------
c
      ierr=0
c
      read (iun,*,err=100,end=100) f
c
      return
c
  100 continue
c
c ****** Error in reading the data.
c
      ierr=1
c
      return
      end
c#######################################################################
      subroutine wrtxt (fname,s,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Write a 1D, 2D, or 3D scientific data set to a text file.
c
c-----------------------------------------------------------------------
c
c ****** Input arguments:
c
c          FNAME   : [character(*)]
c                    Text data file name to write to.
c
c          S       : [structure of type SDS]
c                    A structure that holds the field, its
c                    dimensions, and the scales, with the
c                    components described below.
c
c ****** Output arguments:
c
c          IERR    : [integer]
c                    IERR=0 is returned if the data set was written
c                    successfully.  Otherwise, IERR is set to a
c                    nonzero value.
c
c ****** Components of structure S:
c
c          NDIM    : [integer]
c                    Number of dimensions in the data set.
c
c          DIMS    : [integer, dimension(3)]
c                    Number of points in the data set dimensions.
c                    Only DIMS(1 .. NDIM) are referenced.
c
c          SCALE   : [logical]
c                    Flag to indicate the presence of scales (axes)
c                    in the data set.  SCALE=.false. means that scales
c                    are not being supplied; SCALE=.true. means that
c                    scales are being supplied.
c
c          HDF32   : [logical]
c                    Flag that indicates the precision of the data.
c                    This flag is used to determine the format for data
c                    written to the text file.  When HDF32=.TRUE., the
c                    data is assumed to originate from a 32-bit HDF data
c                    file, and is written with 7 digits to the text file.
c                    Otherwise, the data is assumed to originate from a
c                    64-bit HDF data file, and is written with 14 digits
c                    to the text file.
c
c          SCALES  : [structure of type RP1D, dimension(3)]
c                    This array holds the pointers to the scales
c                    when SCALE=.true., and is not referenced
c                    otherwise.
c
c          F       : [real, pointer to a rank-3 array]
c                    This array holds the data set values.
c
c-----------------------------------------------------------------------
c
      use number_types
      use sds_def
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character(*) :: fname
      type(sds) :: s
      integer :: ierr
      integer :: int_tmp(1)
      intent(in) :: fname,s
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
c ****** Declarations for temporary variables.
c
      integer :: i
      integer :: n
      character(32) :: fmt
c
c-----------------------------------------------------------------------
c
      ierr=0
c
c ****** Open the file for writing.
c
      call ffopen (1,fname,'rw',ierr)
c
      if (ierr.ne.0) then
        write (*,*)
        write (*,*) '### ERROR in WRTXT:'
        write (*,*) '### Could not open the text file for writing.'
        write (*,*) 'File name: ',trim(fname)
        ierr=1
        return
      end if
c
c ****** Check the number of dimensions.
c
      if (s%ndim.le.0.or.s%ndim.gt.3) then
        write (*,*)
        write (*,*) '### ERROR in WRTXT:'
        write (*,*) '### Could not write the SDS data.'
        write (*,*) 'Invalid number of dimensions.'
        write (*,*) 'NDIM = ',s%ndim
        write (*,*) 'File name: ',trim(fname)
        ierr=1
        return
      end if
c
c ****** Construct the format string for writing floating point
c ****** numbers to the output file.
c
      if (s%hdf32) then
        fmt='(5(1x,1pe13.6))'
      else
        fmt='(3(1x,1pe21.14))'
      end if
c
c ****** Write the number of dimensions.
c
      int_tmp=s%ndim
      call wrint (1,1,int_tmp,ierr)
      if (ierr.ne.0) go to 900
c
c ****** Write the dimensions.
c
      do i=1,s%ndim
        call wrint (1,1,s%dims(i),ierr)
        if (ierr.ne.0) go to 900
      enddo
c
c ****** Write the scales.
c
      if (s%scale) then
        int_tmp=1
        call wrint (1,1,int_tmp,ierr)
        if (ierr.ne.0) go to 900
        do i=1,s%ndim
          call wrfp (1,s%dims(i),s%scales(i)%f,fmt,ierr)
          if (ierr.ne.0) go to 900
        enddo
      else
        int_tmp=0
        call wrint (1,1,int_tmp,ierr)
        if (ierr.ne.0) go to 900
      end if
c
c ****** Write the array.
c
      n=product(s%dims(1:s%ndim))
      call wrfp (1,n,s%f,fmt,ierr)
      if (ierr.ne.0) go to 900
c
      close (1)
c
      return
c
  900 continue
c
      write (*,*)
      write (*,*) '### ERROR in WRTXT:'
      write (*,*) '### Error in writing data to the text file.'
      write (*,*) 'File name: ',trim(fname)
      ierr=2
c
      return
      end
c#######################################################################
      subroutine wrint (iun,n,i,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Write N words of INTEGER data from array I to the file
c ****** connected to unit IUN using a free format text write.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer :: iun
      integer :: n
      integer, dimension(n) :: i
      integer :: ierr
      intent(in) :: iun,n,i
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
      ierr=0
c
      write (iun,*,err=100) i
c
      return
c
  100 continue
c
c ****** Error in writing the data.
c
      ierr=1
c
      return
      end
c#######################################################################
      subroutine wrfp (iun,n,f,fmt,ierr)
c
c-----------------------------------------------------------------------
c
c ****** Write N words of REAL data from array F to the file
c ****** connected to unit IUN.
c
c ****** FMT specifies the format string to use.
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
      integer :: iun
      integer :: n
      real(r_typ), dimension(n) :: f
      character(*) :: fmt
      integer :: ierr
      intent(in) :: iun,n,f,fmt
      intent(out) :: ierr
c
c-----------------------------------------------------------------------
c
      ierr=0
c
      write (iun,fmt=fmt,err=100) f
c
      return
c
  100 continue
c
c ****** Error in writing the data.
c
      ierr=1
c
      return
      end
