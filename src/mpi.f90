module mpi
    implicit none
    integer, parameter :: MPI_THREAD_FUNNELED = 1
    ! not sure if this is correct really
    integer, parameter :: MPI_INTEGER = 0
    integer, parameter :: MPI_REAL4 = 0
    integer, parameter :: MPI_REAL8 = 1
    integer, parameter :: MPI_COMM_TYPE_SHARED = 1
    integer, parameter :: MPI_PROC_NULL = -1
    integer, parameter :: MPI_SUCCESS = 0

    integer, parameter :: MPI_COMM_WORLD = 0
    real(8), parameter :: MPI_IN_PLACE = -1
    integer, parameter :: MPI_SUM = 1
    integer, parameter :: MPI_INFO_NULL = 0
    integer :: MPI_STATUS_IGNORE = 0
    ! NOTE: I've no idea for how to implement this, refer
    ! see section 2.5.4 page 21 of mpi40-report.pdf
    ! this is probably not correct right now
    integer :: MPI_STATUSES_IGNORE(1024)

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
        module procedure MPI_Allreduce_array
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
        use iso_c_binding, only: c_int, c_ptr, c_null_ptr
        integer, optional, intent(out) :: ierr
        integer(c_int) :: local_ierr
        integer(c_int) :: argc
        type(c_ptr) :: argv = c_null_ptr
        argc = 0
        ! Call C MPI_Init directly with argc=0, argv=NULL
        local_ierr = c_mpi_init(argc, argv)

        if (present(ierr)) then
            ierr = int(local_ierr)
        else if (local_ierr /= 0) then
            print *, "MPI_Init failed with error code: ", local_ierr
        end if
    end subroutine MPI_Init_proc

    subroutine MPI_Init_thread_proc(required, provided, ierr)
        use mpi_c_bindings, only : c_mpi_init_thread
        use iso_c_binding, only: c_int
        integer, intent(in) :: required
        integer, intent(out) :: provided
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

    subroutine MPI_Allreduce_array(sendbuf, recvbuf, count, datatype, op, comm, ierror)
        use mpi_c_bindings, only: c_mpi_allreduce_array
        ! Declare both send and recv as arrays:
        real(8), dimension(:), intent(in)  :: sendbuf
        real(8), dimension(:), intent(out) :: recvbuf
        integer, intent(in) :: count, datatype, op, comm
        integer, intent(out), optional :: ierror
        call c_mpi_allreduce_array(sendbuf, recvbuf, count, datatype, op, comm, ierror)
    end subroutine MPI_Allreduce_array

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
        call c_mpi_waitall(count, array_of_requests, array_of_statuses, ierror)
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