c#######################################################################
      module number_types
c
c-----------------------------------------------------------------------
c ****** Basic number types.
c ****** This module is used to set the default precision for REALs.
c-----------------------------------------------------------------------
c
      use iso_fortran_env
c
c-----------------------------------------------------------------------
c
      implicit none
c
      integer, parameter :: KIND_REAL_4=REAL32
      integer, parameter :: KIND_REAL_8=REAL64
      integer, parameter :: KIND_REAL_16=max(REAL128,REAL64)
c
      integer, parameter :: r_typ=KIND_REAL_8
c
      end module
c#######################################################################
