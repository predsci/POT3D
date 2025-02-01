module mpi
    implicit none
    integer, parameter :: MPI_THREAD_FUNNELED = 1
    ! not sure if this is correct really
    integer, parameter :: MPI_INTEGER = 0
    integer, parameter :: MPI_REAL4 = 0
    integer, parameter :: MPI_REAL8 = 1
    integer, parameter :: MPI_COMM_TYPE_SHARED = 1
    integer, parameter :: MPI_PROC_NULL = -1

    integer, parameter :: MPI_COMM_WORLD = 0
    real(8), parameter :: MPI_IN_PLACE = 1
    integer, parameter :: MPI_SUM = 1
    integer, parameter :: MPI_INFO_NULL = 0
    integer :: MPI_STATUS_IGNORE = 0
    ! NOTE: I've no idea for how to implement this, refer
    ! see section 2.5.4 page 21 of mpi40-report.pdf
    ! this is probably not correct right now
    integer, allocatable :: MPI_STATUSES_IGNORE(:)

    ! not used in pot3d.F90
    interface MPI_Init
        module procedure MPI_Init_proc
    end interface MPI_Init

    interface MPI_Init_thread
        module procedure MPI_Init_thread_proc
    end interface MPI_Init_thread

    interface MPI_Finalize
        module procedure MPI_Finalize_proc
    end interface MPI_Finalize

    interface MPI_Comm_size
        module procedure MPI_Comm_size_proc
    end interface MPI_Comm_size

    interface MPI_Bcast
        module procedure MPI_Bcast_int
        module procedure MPI_Bcast_real
    end interface MPI_Bcast

    interface MPI_Allgather
        module procedure MPI_Allgather_int
        module procedure MPI_Allgather_real
    end interface MPI_Allgather

    interface MPI_Isend
        module procedure MPI_Isend_2d
        module procedure MPI_Isend_3d
    end interface

    interface MPI_IRecv
        module procedure MPI_IRecv_proc
    end interface

    interface MPI_Allreduce
        module procedure MPI_Allreduce_scalar
        module procedure MPI_Allreduce_1d
    end interface

    interface MPI_Wtime
        module procedure MPI_Wtime_proc
    end interface

    interface MPI_Barrier
        module procedure MPI_Barrier_proc
    end interface

    interface MPI_Comm_rank
        module procedure MPI_Comm_rank_proc
    end interface

    interface MPI_Comm_split_type
        module procedure MPI_Comm_split_type_proc
    end interface

    interface MPI_Recv
        module procedure MPI_Recv_proc
    end interface

    interface MPI_Waitall
        module procedure MPI_Waitall_proc
    end interface

    interface MPI_Ssend
        module procedure MPI_Ssend_proc
    end interface

    interface MPI_Cart_create
        module procedure MPI_Cart_create_proc
    end interface

    interface MPI_Cart_sub
        module procedure MPI_Cart_sub_proc
    end interface

    interface MPI_Cart_shift
        module procedure MPI_Cart_shift_proc
    end interface

    interface MPI_Dims_create
        module procedure MPI_Dims_create_proc
    end interface

    interface MPI_Cart_coords
        module procedure MPI_Cart_coords_proc
    end interface

   contains

    subroutine MPI_Init_proc(ierr)
        use mpi_c_bindings, only: c_mpi_init
        use iso_c_binding, only : c_int
        integer, optional, intent(out) :: ierr
        integer :: local_ierr
        if (present(ierr)) then
            call c_mpi_init(ierr)
        else
            call c_mpi_init(local_ierr)
            if (local_ierr /= 0) then
                print *, "MPI_Init failed with error code: ", local_ierr
            end if
        end if
    end subroutine

    subroutine MPI_Init_thread_proc(required, provided, ierr)
        use mpi_c_bindings, only : c_mpi_init_thread
        use iso_c_binding, only: c_int
        integer, intent(in) :: required
        integer, intent(in) :: provided
        integer, optional, intent(out) :: ierr
        integer :: local_ierr
        if (present(ierr)) then
            call c_mpi_init_thread(required, provided, ierr)
        else
            call c_mpi_init_thread(required, provided, local_ierr)
            if (local_ierr /= 0) then
                print *, "MPI_Init_thread failed with error code: ", local_ierr
            end if
        end if
    end subroutine

    subroutine MPI_Finalize_proc(ierr)
        use mpi_c_bindings, only: c_mpi_finalize
        use iso_c_binding, only: c_int
        integer, optional, intent(out) :: ierr
        integer :: local_ierr
        if (present(ierr)) then
            call c_mpi_finalize(ierr)
        else
            call c_mpi_finalize(local_ierr)
            if (local_ierr /= 0) then
                print *, "MPI_Finalize failed with error code: ", local_ierr
            end if
        end if
    end subroutine

    subroutine MPI_Comm_size_proc(comm, size, ierr)
        use mpi_c_bindings, only: c_mpi_comm_size
        integer, intent(in) :: comm
        integer, intent(out) :: size
        integer, optional, intent(out) :: ierr
        integer :: local_ierr
        if (present(ierr)) then
            call c_mpi_comm_size(comm, size, ierr)
        else
            call c_mpi_comm_size(comm, size, local_ierr)
            if (local_ierr /= 0) then
                print *, "MPI_Comm_size failed with error code: ", local_ierr
            end if
        end if
    end subroutine

    subroutine MPI_Bcast_int(buffer, count, datatype, root, comm, ierror)
        use mpi_c_bindings, only: c_mpi_bcast_int
        integer :: buffer
        integer, intent(in) :: count, root
        integer, intent(in) :: datatype
        integer, intent(in) :: comm
        integer, optional, intent(out) :: ierror
        call c_mpi_bcast_int(buffer, count, datatype, root, comm, ierror)
    end subroutine

    subroutine MPI_Bcast_real(buffer, count, datatype, root, comm, ierror)
        use mpi_c_bindings, only: c_mpi_bcast_real
        real(8), dimension(:, :) :: buffer
        integer, intent(in) :: count, root
        integer, intent(in) :: datatype
        integer, intent(in) :: comm
        integer, optional, intent(out) :: ierror
        call c_mpi_bcast_real(buffer, count, datatype, root, comm, ierror)
    end subroutine

    subroutine MPI_Allgather_int(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierror)
        use mpi_c_bindings, only: c_mpi_allgather_int
        integer, dimension(:), intent(in) :: sendbuf
        integer, dimension(:, :) :: recvbuf
        integer, intent(in) :: sendcount, recvcount
        integer, intent(in) :: sendtype, recvtype
        integer, intent(in) :: comm
        integer, optional, intent(out) :: ierror
        call c_mpi_allgather_int(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierror)
    end subroutine

    subroutine MPI_Allgather_real(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierror)
        use mpi_c_bindings, only: c_mpi_allgather_real
        real(8), dimension(:), intent(in) :: sendbuf
        real(8), dimension(:, :) :: recvbuf
        integer, intent(in) :: sendcount, recvcount
        integer, intent(in) :: sendtype, recvtype
        integer, intent(in) :: comm
        integer, optional, intent(out) :: ierror
        call c_mpi_allgather_real(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierror)
    end subroutine

    subroutine MPI_Isend_2d(buf, count, datatype, dest, tag, comm, request, ierror)
        use mpi_c_bindings, only: c_mpi_isend
        real(8), dimension(:, :), intent(in) :: buf
        integer, intent(in) :: count, dest, tag
        integer, intent(in) :: datatype
        integer, intent(in) :: comm
        integer, intent(out) :: request
        integer, optional, intent(out) :: ierror
        call c_mpi_isend(buf, count, datatype, dest, tag, comm, request, ierror)
    end subroutine

    subroutine MPI_Isend_3d(buf, count, datatype, dest, tag, comm, request, ierror)
        use mpi_c_bindings, only: c_mpi_isend
        real(8), dimension(:, :, :), intent(in) :: buf
        integer, intent(in) :: count, dest, tag
        integer, intent(in) :: datatype
        integer, intent(in) :: comm
        integer, intent(out) :: request
        integer, optional, intent(out) :: ierror
        call c_mpi_isend(buf, count, datatype, dest, tag, comm, request, ierror)
    end subroutine

    subroutine MPI_Irecv_proc(buf, count, datatype, source, tag, comm, request, ierror)
        use mpi_c_bindings, only: c_mpi_irecv
        real(8), dimension(:,:) :: buf
        integer, intent(in) :: count, source, tag
        integer, intent(in) :: datatype
        integer, intent(in) :: comm
        integer, intent(out) :: request
        integer, optional, intent(out) :: ierror
        call c_mpi_irecv(buf, count, datatype, source, tag, comm, request, ierror)
    end subroutine

    subroutine MPI_Allreduce_scalar(sendbuf, recvbuf, count, datatype, op, comm, ierror)
        use mpi_c_bindings, only: c_mpi_allreduce_scalar
        real(8), intent(in) :: sendbuf
        real(8), intent(out) :: recvbuf
        integer, intent(in) :: count, datatype, op, comm
        integer, intent(out), optional :: ierror
        call c_mpi_allreduce_scalar(sendbuf, recvbuf, count, datatype, op, comm, ierror)
    end subroutine

    subroutine MPI_Allreduce_1d(sendbuf, recvbuf, count, datatype, op, comm, ierror)
        use mpi_c_bindings, only: c_mpi_allreduce_1d
        real(8), intent(in) :: sendbuf
        real(8), dimension(:), intent(out) :: recvbuf
        integer, intent(in) :: count, datatype, op, comm
        integer, intent(out), optional :: ierror
        call c_mpi_allreduce_1d(sendbuf, recvbuf, count, datatype, op, comm, ierror)
    end subroutine

    function MPI_Wtime_proc() result(time)
        use mpi_c_bindings, only: c_mpi_wtime
        real(8) :: time
        time = c_mpi_wtime()
    end function

    subroutine MPI_Barrier_proc(comm, ierror)
        use mpi_c_bindings, only: c_mpi_barrier
        integer, intent(in) :: comm
        integer, intent(out), optional :: ierror
        call c_mpi_barrier(comm, ierror)
    end subroutine

    subroutine MPI_Comm_rank_proc(comm, rank, ierror)
        use mpi_c_bindings, only: c_mpi_comm_rank
        integer, intent(in) :: comm
        integer, intent(out) :: rank
        integer, optional, intent(out) :: ierror
        call c_mpi_comm_rank(comm, rank, ierror)
    end subroutine

    subroutine MPI_Comm_split_type_proc(comm, split_type, key, info, newcomm, ierror)
        use mpi_c_bindings, only: c_mpi_comm_split_type
        integer :: comm
        integer, intent(in) :: split_type, key
        integer, intent(in) :: info
        integer, intent(out) :: newcomm
        integer, optional, intent(out) :: ierror
        call c_mpi_comm_split_type(comm, split_type, key, info, newcomm, ierror)
    end subroutine

    subroutine MPI_Recv_proc(buf, count, datatype, source, tag, comm, status, ierror)
        use mpi_c_bindings, only: c_mpi_recv
        real(8), dimension(:) :: buf
        integer, intent(in) :: count, source, tag
        integer, intent(in) :: datatype
        integer, intent(in) :: comm
        integer, intent(out) :: status
        integer, optional, intent(out) :: ierror
        call c_mpi_recv(buf, count, datatype, source, tag, comm, status, ierror)
    end subroutine

    subroutine MPI_Waitall_proc(count, array_of_requests, array_of_statuses, ierror)
        use mpi_c_bindings, only: c_mpi_waitall
        integer, intent(in) :: count
        integer, intent(inout) :: array_of_requests(count)
        integer :: array_of_statuses(*)
        integer, optional, intent(out) :: ierror
        call c_mpi_waitall(count, array_of_requests, array_of_requests, ierror)
    end subroutine

    subroutine MPI_Ssend_proc(buf, count, datatype, dest, tag, comm, ierror)
        use mpi_c_bindings, only: c_mpi_ssend
        real(8), dimension(*), intent(in) :: buf
        integer, intent(in) :: count, dest, tag
        integer, intent(in) :: datatype
        integer, intent(in) :: comm
        integer, optional, intent(out) :: ierror
        call c_mpi_ssend(buf, count, datatype, dest, tag, comm, ierror)
    end subroutine

    subroutine MPI_Cart_create_proc(comm, ndims, dims, periods, reorder, newcomm, ierror)
        use mpi_c_bindings, only: c_mpi_cart_create
        use iso_c_binding, only: c_int
        integer, intent(in) :: ndims, dims(ndims)
        logical, intent(in) :: periods(ndims), reorder
        integer, intent(in) :: comm
        integer, intent(out) :: newcomm
        integer, optional, intent(out) :: ierror
        integer(c_int) :: ndims_c, reorder_c, dims_c(ndims), periods_c(ndims)
            ndims_c = ndims
            if (reorder) then
                reorder_c = 1
            else
                reorder_c = 0
            end if
            dims_c = dims
            where (periods)
                periods_c = 1
            elsewhere
                periods_c = 0
            end where
        call c_mpi_cart_create(comm, ndims, dims_c, periods_c, reorder_c, newcomm, ierror)
    end subroutine

    subroutine MPI_Cart_coords_proc(comm, rank, maxdims, coords, ierror)
        use mpi_c_bindings, only: c_mpi_cart_coords
        integer, intent(in) :: comm
        integer, intent(in) :: rank, maxdims
        integer, intent(out) :: coords(maxdims)
        integer, optional, intent(out) :: ierror
        call c_mpi_cart_coords(comm, rank, maxdims, coords, ierror)
    end subroutine

    subroutine MPI_Cart_shift_proc(comm, direction, disp, rank_source, rank_dest, ierror)
        use mpi_c_bindings, only: c_mpi_cart_shift
        integer, intent(in) :: comm
        integer, intent(in) :: direction, disp
        integer, intent(out) :: rank_source, rank_dest
        integer, optional, intent(out) :: ierror
        call c_mpi_cart_shift(comm, direction, disp, rank_source, rank_dest, ierror)
    end subroutine

    subroutine MPI_Dims_create_proc(nnodes, ndims, dims, ierror)
        use mpi_c_bindings, only: c_mpi_dims_create
        integer, intent(in) :: nnodes, ndims
        integer, intent(out) :: dims(ndims)
        integer, optional, intent(out) :: ierror
        call c_mpi_dims_create(nnodes, ndims, dims, ierror)
    end subroutine

    subroutine MPI_Cart_sub_proc (comm, remain_dims, newcomm, ierror)
        use mpi_c_bindings, only: c_mpi_cart_sub
        integer, intent(in) :: comm
        logical, intent(in) :: remain_dims(:)
        integer, intent(out) :: newcomm
        integer, optional, intent(out) :: ierror
        integer :: remain_dims_i(size(remain_dims))
        where (remain_dims)
            remain_dims_i = 1
        elsewhere
            remain_dims_i = 0
        end where
        call c_mpi_cart_sub(comm, remain_dims_i, newcomm, ierror)
    end subroutine
end module mpi

! program main
!     use mpi
!     implicit none
!     integer :: ierr
!     integer :: required
!     integer :: provided
!     integer :: tcheck
!     integer :: ierr0
!     integer :: nproc1
!     integer :: rank = 0
!     integer, parameter :: nt_g = 2
!     integer, parameter :: np_g = 3
!     real(8), dimension(:,:), allocatable :: br0_g
!     integer, parameter :: lbuf=4
!     integer, dimension(lbuf) :: sbuf
!     integer, parameter :: nproc = 2
!     integer, dimension(lbuf,0:nproc-1) :: rbuf
!     integer :: ntype_real = MPI_REAL8
!     integer :: comm_all, newcomm_all, newcomm_after_sub
!     integer, parameter :: nr = 2
!     integer, parameter :: nt = 3
!     integer, parameter :: np = 2
!     real(8), dimension(nr,nt,np) :: a
!     integer, parameter :: n1 = 2
!     real(8), dimension(n1) :: a0
!     integer :: tag=0
!     real(8), dimension(:), allocatable :: rbuf4
!     integer, parameter :: lbuf4 = 4
!     integer :: irank4 = 2
    
!     integer, parameter :: maxdim = 2
!     integer :: coords(maxdim)
!     integer :: tag4 = 0

!     integer, parameter :: lbuf2=10
!     real(8), dimension(lbuf2) :: sbuf2
!     real(8), dimension(lbuf2,0:nproc-1) :: tbuf
!     integer, parameter :: lbuf3 = 10
!     integer :: comm_phi
!     integer :: iproc_pp
!     integer :: reqs(4)
!     integer, parameter :: nstack=10
!     real(8), dimension(nstack) :: tstart=0.
!     integer :: istack=1
!     integer :: iprocw
!     integer :: comm_shared
!     real(8), dimension(:), allocatable :: sbuf5
!     integer, parameter :: lsbuf5 = 24
!     real(8), dimension(2,3,4) :: a5

!     integer :: iproc05 = 10
!     integer, parameter :: ndims = 2
!     integer :: dims(ndims) = 1
!     logical :: periods(ndims), remain_dims(ndims) = .false.
!     integer :: direction, displ, source, dest

!     allocate (br0_g(nt_g,np_g))
!     allocate (rbuf4(lbuf4))


!     allocate (sbuf5(lsbuf5))
!     sbuf5=reshape(a5(1:2,1:3,1:4),(/lsbuf5/))

!     ! NOTE: called in pot3d.F90 as:
!     call MPI_Init_thread (MPI_THREAD_FUNNELED,tcheck,ierr)
!     ! call MPI_Init_thread(required, provided, ierr)
!     if (ierr /= 0) error stop "MPI_Init_thread failed"

!     ierr = -1
!     ! NOTE: called in pot3d.F90 as:
!     call MPI_Comm_size (MPI_COMM_WORLD,nproc1,ierr)
!     if (ierr /= 0) error stop
!     print *, "Number of processes:", nproc1

!     ierr = -1
!     call MPI_Bcast (ierr0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!     if (ierr /= 0) error stop

!     ierr = -1
!     call MPI_Bcast(br0_g,nt_g*np_g,ntype_real,0,comm_all,ierr)
!     if (ierr /= 0) error stop


!     ierr = -1
!     call MPI_Allgather (sbuf,lbuf,MPI_INTEGER, &
!         rbuf,lbuf,MPI_INTEGER,comm_all,ierr)
    
!     if (ierr /= 0) error stop

!     ierr = -1
!     call MPI_Allgather (sbuf2,lbuf2,ntype_real, &
!         tbuf,lbuf2,ntype_real,comm_all,ierr)

!     if (ierr /= 0) error stop

!     ierr = -1
!     call MPI_Isend (a(:,:,np-1),lbuf3,ntype_real,iproc_pp,tag, &
!         comm_all,reqs(1),ierr)
!     if (ierr /= 0) error stop

!     ierr = -1
!     call MPI_Irecv (a(:,:, 1),lbuf3,ntype_real,iproc_pp,tag,   &
!         comm_all,reqs(3),ierr)
!     if (ierr /= 0) error stop

!     ierr = -1
!     call MPI_Allreduce (MPI_IN_PLACE,a0,n1,ntype_real, &
!         MPI_SUM,comm_phi,ierr)
!     if (ierr /= 0) error stop

!     tstart(istack)=MPI_Wtime()

!     ierr = -1
!     call MPI_Barrier(MPI_COMM_WORLD,ierr)
!     if (ierr /= 0) error stop

!     ierr = -1
!     call MPI_Barrier (comm_all,ierr)
!     if (ierr /= 0) error stop

!     ierr = -1
!     call MPI_Comm_rank (MPI_COMM_WORLD,iprocw,ierr)
!     if (ierr /= 0) error stop
!     print *, iprocw

!     ierr = -1
!     call MPI_Comm_split_type (MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0, &
!         MPI_INFO_NULL,comm_shared,ierr)
!     if (ierr /= 0) error stop

!     ! I'm not sure of the exact conditions to test "MPI_Recv",
!     ! currently it fails on program execution
!     ! ierr = -1
!     ! call MPI_Recv (rbuf4,lbuf4,ntype_real,irank4,tag4, &
!     !                comm_all,MPI_STATUS_IGNORE,ierr)
!     ! if (ierr /= 0) error stop

!     ! maybe this required a very specific "rank" for the arguments
!     ! ierr = -1
!     ! call MPI_Ssend (sbuf5,lsbuf5,ntype_real,iproc05,tag,comm_all,ierr)
!     ! if (ierr /= 0) error stop

!     ! ierr = -1
!     ! I've no idea for why the below works, I don't understand
!     ! things here, cause I've no idea for why declaring
!     ! MPI_STATUSES_IGNORE as an integer allocatable makes it work here
!     ! call MPI_Waitall (4,reqs,MPI_STATUSES_IGNORE,ierr)
!     ! if (ierr /= 0) error stop

!     ierr = -1
!     call MPI_Dims_create (nproc1, ndims, dims, ierr)
!     if (ierr /= 0) error stop
!     print *, "Computed dimensions:", dims
!     ! Here Ideally we would need MPI_Dims_create(size,2,dims) as dims() value can't be zero it have to be initialized
!     ierr = -1
!     ! Just experimental
!     ! periods(1) = .TRUE.
!     call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, .FALSE., newcomm_all, ierr)
!     if (ierr /= 0) error stop
!     print *, "Cartesian communicator created"

!     ierr = -1
!     call MPI_Comm_rank (newcomm_all,iprocw,ierr)
!     if (ierr /= 0) error stop
!     print *, "Cartesian rank:", iprocw

!     ierr = -1
!     call MPI_Cart_coords(newcomm_all, iprocw, 2, coords, ierr)
!     if (ierr /= 0) error stop
!     print *, "Coordinates:", coords

!     ierr = -1
!     call MPI_Cart_shift(newcomm_all, direction, displ, source, dest, ierr)
!     if (ierr /= 0) error stop
!     print *, "Shift results - Source:", source, "Dest:", dest

!     ierr = -1
!     call MPI_Cart_sub(newcomm_all, remain_dims, newcomm_after_sub, ierr)
!     if (ierr /= 0) error stop

!     print *, "Initial Communicator created:", newcomm_all
!     print *, "Subcommunicator created:", newcomm_after_sub

!     ! called in pot3d.F90 as
!     ierr = -1
!     call MPI_Finalize(ierr)
!     if (ierr /= 0) error stop "MPI_Finalize failed"

! end program main
