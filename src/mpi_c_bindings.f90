module mpi_c_bindings
    implicit none

    interface
        subroutine c_mpi_init(ierr) bind(C, name="mpi_init_wrapper")
            use iso_c_binding, only: c_int
            integer(c_int), intent(out) :: ierr
        end subroutine c_mpi_init

        subroutine c_mpi_init_thread(required, provided, ierr) bind(C, name="mpi_init_thread_wrapper")
            use iso_c_binding, only: c_int
            integer(c_int), intent(in) :: required
            integer(c_int), intent(in) :: provided
            integer(c_int), intent(out) :: ierr
        end subroutine c_mpi_init_thread

        subroutine c_mpi_finalize(ierr) bind(C, name="mpi_finalize_wrapper")
            use iso_c_binding, only : c_int
            integer(c_int), intent(out) :: ierr
        end subroutine c_mpi_finalize

        subroutine c_mpi_comm_size(comm, size, ierr) bind(C, name="mpi_comm_size_wrapper")
            use iso_c_binding, only: c_int
            integer(c_int), intent(in) :: comm
            integer(c_int), intent(out) :: size
            integer(c_int), intent(out) :: ierr
        end subroutine c_mpi_comm_size

        subroutine c_mpi_bcast_int(buffer, count, datatype, root, comm, ierror) bind(C, name="mpi_bcast_int_wrapper")
            use iso_c_binding, only: c_int
            integer(c_int) :: buffer
            integer(c_int), intent(in) :: count, root
            integer(c_int), intent(in) :: datatype
            integer(c_int), intent(in) :: comm
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_bcast_real(buffer, count, datatype, root, comm, ierror) bind(C, name="mpi_bcast_real_wrapper")
            use iso_c_binding, only : c_int, c_double
            real(c_double), dimension(*) :: buffer
            integer(c_int), intent(in) :: count, root
            integer(c_int), intent(in) :: datatype
            integer(c_int), intent(in) :: comm
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_allgather_int(sendbuf, sendcount, sendtype, recvbuf, &
            recvcount, recvtype, comm, ierror) bind(C, name="mpi_allgather_int_wrapper")
            use iso_c_binding, only: c_int
            integer(c_int), dimension(*), intent(in) :: sendbuf
            integer(c_int), dimension(*) :: recvbuf
            integer(c_int), intent(in) :: sendcount, recvcount
            integer(c_int), intent(in) :: sendtype, recvtype
            integer(c_int), intent(in) :: comm
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_allgather_real(sendbuf, sendcount, sendtype, recvbuf, &
                    recvcount, recvtype, comm, ierror) bind(C, name="mpi_allgather_real_wrapper")
            use iso_c_binding, only: c_int, c_double
            real(c_double), dimension(*), intent(in) :: sendbuf
            real(c_double), dimension(*) :: recvbuf
            integer(c_int), intent(in) :: sendcount, recvcount
            integer(c_int), intent(in) :: sendtype, recvtype
            integer(c_int), intent(in) :: comm
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_isend(buf, count, datatype, dest, tag, comm, request, ierror) bind(C, name="mpi_isend_wrapper")
            use iso_c_binding, only: c_int, c_double
            real(c_double), dimension(*), intent(in) :: buf
            integer(c_int), intent(in) :: count, dest, tag
            integer(c_int), intent(in) :: datatype
            integer(c_int), intent(in) :: comm
            integer(c_int), intent(out) :: request
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_irecv(buf, count, datatype, source, tag, comm, request, ierror) bind(C, name="mpi_irecv_wrapper")
            use iso_c_binding, only: c_int, c_double
            real(c_double), dimension(*) :: buf
            integer(c_int), intent(in) :: count, source, tag
            integer(c_int), intent(in) :: datatype
            integer(c_int), intent(in) :: comm
            integer(c_int), intent(out) :: request
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_allreduce_scalar(sendbuf, recvbuf, count, datatype, op, comm, ierror) bind(C, name="mpi_allreduce_wrapper")
            use iso_c_binding, only: c_int, c_double
            real(c_double), intent(in) :: sendbuf
            real(c_double), intent(out) :: recvbuf
            integer(c_int), intent(in) :: count, datatype, op, comm
            integer(c_int), intent(out), optional :: ierror
        end subroutine

        subroutine c_mpi_allreduce_1d(sendbuf, recvbuf, count, datatype, op, comm, ierror) bind(C, name="mpi_allreduce_wrapper")
            use iso_c_binding, only: c_int, c_double
            real(c_double), intent(in) :: sendbuf
            real(c_double), dimension(*), intent(out) :: recvbuf
            integer(c_int), intent(in) :: count, datatype, op, comm
            integer(c_int), intent(out), optional :: ierror
        end subroutine

        function c_mpi_wtime() result(time) bind(C, name="mpi_wtime_wrapper")
            use iso_c_binding, only: c_double
            real(c_double) :: time
        end function

        subroutine c_mpi_barrier(comm, ierror) bind(C, name="mpi_barrier_wrapper")
            use iso_c_binding, only: c_int
            integer(c_int), intent(in) :: comm
            integer(c_int), intent(out), optional :: ierror
        end subroutine

        subroutine c_mpi_comm_rank(comm, rank, ierror) bind(C, name="mpi_comm_rank_wrapper")
            use iso_c_binding, only: c_int
            integer(c_int), intent(in) :: comm
            integer(c_int), intent(out) :: rank
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_comm_split_type(comm, split_type, key, info, newcomm, ierror) bind(C, name="mpi_comm_split_type_wrapper")
            use iso_c_binding, only: c_int
            integer(c_int) :: comm
            integer(c_int), intent(in) :: split_type, key
            integer(c_int), intent(in) :: info
            integer(c_int), intent(out) :: newcomm
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_recv(buf, count, datatype, source, tag, comm, status, ierror) bind(C, name="mpi_recv_wrapper")
            use iso_c_binding, only: c_int, c_double
            real(c_double), dimension(*) :: buf
            integer(c_int), intent(in) :: count, source, tag
            integer(c_int), intent(in) :: datatype
            integer(c_int), intent(in) :: comm
            integer(c_int), intent(out) :: status
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_waitall(count, array_of_requests, array_of_statuses, ierror) bind(C, name="mpi_waitall_wrapper")
            use iso_c_binding, only: c_int
            integer(c_int), intent(in) :: count
            integer(c_int), intent(inout) :: array_of_requests(count)
            integer(c_int) :: array_of_statuses(*)
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_ssend(buf, count, datatype, dest, tag, comm, ierror) bind(C, name="mpi_ssend_wrapper")
            use iso_c_binding, only: c_int, c_double
            real(c_double), dimension(*), intent(in) :: buf
            integer(c_int), intent(in) :: count, dest, tag
            integer(c_int), intent(in) :: datatype
            integer(c_int), intent(in) :: comm
            integer(c_int), optional, intent(out) :: ierror
        end subroutine

        subroutine c_mpi_cart_create(comm, ndims, dims, periods, reorder, newcomm, ierror) bind(C, name="mpi_cart_create_wrapper")
            use iso_c_binding, only: c_int
            integer(c_int), intent(in) :: comm, ndims, reorder
            integer(c_int), intent(in) :: dims(*), periods(*)
            integer(c_int), intent(out) :: newcomm, ierror
        end subroutine

        subroutine c_mpi_cart_coords(comm, rank, maxdims, coords, ierror) bind(C, name="mpi_cart_coords_wrapper")
            use iso_c_binding,  only: c_int
            integer(c_int), intent(in) :: comm, rank, maxdims
            integer(c_int), intent(out) :: coords(*), ierror
        end subroutine 

        subroutine c_mpi_cart_shift(comm, direction, disp, rank_source, rank_dest, ierror) bind(C, name="mpi_cart_shift_wrapper")
            use iso_c_binding, only: c_int
            integer(c_int), intent(in) :: comm, direction, disp
            integer(c_int), intent(out) :: rank_source, rank_dest, ierror
        end subroutine

        subroutine c_mpi_dims_create(nnodes, ndims, dims, ierror) bind(C, name="mpi_dims_create_wrapper")
            use iso_c_binding, only: c_int
            integer(c_int), intent(in) :: nnodes, ndims
            integer(c_int), intent(inout) :: dims(*), ierror
        end subroutine

        subroutine c_mpi_cart_sub(comm, remain_dims, newcomm, ierror) bind(C, name ="mpi_cart_sub_wrapper")
            use iso_c_binding, only: c_int
            integer(c_int), intent(in) :: comm
            integer(c_int), intent(in) :: remain_dims(*)
            integer(c_int), intent(out) :: newcomm, ierror
        end subroutine
    end interface
end module mpi_c_bindings
