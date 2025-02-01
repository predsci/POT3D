#include <mpi.h>
#include <stdlib.h>

#define MPI_STATUS_SIZE 5

void mpi_init_wrapper(int *ierr) {
    int argc = 0;
    char **argv = NULL;
    *ierr = MPI_Init(&argc, &argv);
}

void mpi_init_thread_wrapper(int *required, int *provided, int *ierr) {
    int argc = 0;
    char **argv = NULL;
    *ierr = MPI_Init_thread(&argc, &argv, *required, provided);
}

void mpi_finalize_wrapper(int *ierr) {
    *ierr = MPI_Finalize();
}

void mpi_comm_size_wrapper(int *comm_f, int *size, int *ierr) {
    MPI_Comm comm = MPI_Comm_f2c(*comm_f);
    *ierr = MPI_Comm_size(comm, size);
}

void mpi_bcast_int_wrapper(int *buffer, int *count, int *datatype_f, int *root, int *comm_f, int *ierror) {
    MPI_Comm comm = MPI_Comm_f2c(*comm_f);
    MPI_Datatype datatype;
    switch (*datatype_f) {
        case 0:
            datatype = MPI_INT;
            break;
        case 1:
            datatype = MPI_FLOAT;
            break;
        case 2:
            datatype = MPI_DOUBLE;
            break;
        default:
            *ierror = -1;
            return;
    }
    *ierror = MPI_Bcast(buffer, *count, datatype, *root, comm);
}

void mpi_bcast_real_wrapper(double *buffer, int *count, int *datatype_f, int *root, int *comm_f, int *ierror) {
    MPI_Comm comm = MPI_Comm_f2c(*comm_f);
    MPI_Datatype datatype;
    switch (*datatype_f) {
        case 0:
            datatype = MPI_FLOAT;
            break;
        case 1:
            datatype = MPI_DOUBLE;
            break;
        default:
            *ierror = -1;
            return;
    }
    *ierror = MPI_Bcast(buffer, *count, datatype, *root, comm);
}

void mpi_allgather_int_wrapper(const int *sendbuf, int *sendcount, int *sendtype_f,
                               int *recvbuf, int *recvcount, int *recvtype_f, 
                               int *comm_f, int *ierror) {
    MPI_Comm comm = MPI_Comm_f2c(*comm_f);

    MPI_Datatype sendtype, recvtype;
    switch (*sendtype_f) {
        case 0:
            sendtype = MPI_INT;
            break;
        default:
            *ierror = -1;
            return;
    }

    switch (*recvtype_f) {
        case 0:
            recvtype = MPI_INT;
            break;
        default:
            *ierror = -1;
            return;
    }

    *ierror = MPI_Allgather(sendbuf, *sendcount, sendtype,
                            recvbuf, *recvcount, recvtype, comm);
}

void mpi_allgather_real_wrapper(const double *sendbuf, int *sendcount, int *sendtype_f,
                               double *recvbuf, int *recvcount, int *recvtype_f, 
                               int *comm_f, int *ierror) {
    MPI_Comm comm = MPI_Comm_f2c(*comm_f);

    MPI_Datatype sendtype, recvtype;
    switch (*sendtype_f) {
        case 0:
            sendtype = MPI_FLOAT;
            break;
        case 1:
            sendtype = MPI_DOUBLE;
            break;
        default:
            *ierror = -1;
            return;
    }

    switch (*recvtype_f) {
        case 0:
            recvtype = MPI_FLOAT;
            break;
        case 1:
            recvtype = MPI_DOUBLE;
            break;
        default:
            *ierror = -1;
            return;
    }

    *ierror = MPI_Allgather(sendbuf, *sendcount, sendtype,
                            recvbuf, *recvcount, recvtype, comm);
}

void mpi_isend_wrapper(const double *buf, int *count, int *datatype_f,
                        int *dest, int *tag, int *comm_f, int *request_f,
                        int *ierror) {
    MPI_Comm comm = MPI_Comm_f2c(*comm_f);
    MPI_Datatype datatype;
    switch (*datatype_f) {
        case 0:
            datatype = MPI_FLOAT;
            break;
        case 1:
            datatype = MPI_DOUBLE;
            break;
        default:
            *request_f = -1;
            return;
    }

    MPI_Request request;
    *ierror = MPI_Isend(buf, *count, datatype, *dest, *tag, comm, &request);
    *request_f = MPI_Request_c2f(request);
}

void mpi_irecv_wrapper(double *buf, int *count, int *datatype_f,
                        int *source, int *tag, int *comm_f, int *request_f,
                        int *ierror) {
    MPI_Comm comm = MPI_Comm_f2c(*comm_f);
    MPI_Datatype datatype;
    switch (*datatype_f) {
        case 0:
            datatype = MPI_FLOAT;
            break;
        case 1:
            datatype = MPI_DOUBLE;
            break;
        default:
            *request_f = -1;
            return;
    }

    MPI_Request request;
    *ierror = MPI_Irecv(buf, *count, datatype, *source, *tag, comm, &request);
    *request_f = MPI_Request_c2f(request);
}

void mpi_allreduce_wrapper(const double *sendbuf, double *recvbuf, int *count,
                            int *datatype_f, int *op_f, int *comm_f, int *ierror) {
    MPI_Comm comm = MPI_Comm_f2c(*comm_f);
    MPI_Datatype datatype;
    switch (*datatype_f) {
        case 0:
            datatype = MPI_FLOAT;
            break;
        case 1:
            datatype = MPI_DOUBLE;
            break;
        default:
            *ierror = -1;
            return;
    }

    // I'm a little doubtful: as how would it identify that this part
    // is supposed to be for MPI_SUM?
    MPI_Op op = MPI_Op_f2c(*op_f);
    // I've hard-coded as "MPI_SUM" for now, as in POT3D codebase, it's always
    // used as MPI_SUM
    *ierror = MPI_Allreduce(sendbuf, recvbuf, *count, datatype, MPI_SUM, comm);
}

double mpi_wtime_wrapper() {
    return MPI_Wtime();
}

void mpi_barrier_wrapper(int *comm_f, int *ierror) {
    MPI_Comm comm = MPI_Comm_f2c(*comm_f);
    *ierror = MPI_Barrier(comm);
}

void mpi_comm_rank_wrapper(int *comm_f, int *rank, int *ierror) {
    MPI_Comm comm = MPI_Comm_f2c(*comm_f);
    *ierror = MPI_Comm_rank(comm, rank);
}

void mpi_comm_split_type_wrapper(int *comm_f, int *split_type, int *key,
                                int *info_f, int *newcomm_f, int *ierror) {
    MPI_Comm comm = MPI_Comm_f2c(*comm_f);
    MPI_Info info = MPI_Info_f2c(*info_f);
    MPI_Comm newcomm = MPI_Comm_f2c(*newcomm_f);
    *ierror = MPI_Comm_split_type( comm, *split_type, *key , info, &newcomm);
}

void mpi_recv_wrapper(double *buf, int *count, int *datatype_f, int *source,
                    int *tag, int *comm_f, int *status_f, int *ierror) {
    MPI_Datatype datatype;
    switch (*datatype_f) {
        case 0:
            datatype = MPI_FLOAT;
            break;
        case 1:
            datatype = MPI_DOUBLE;
            break;
        default:
            *ierror = -1;
            return;
    }

    MPI_Comm comm = MPI_Comm_f2c(*comm_f);
    MPI_Status status;
    *ierror = MPI_Recv(buf, *count, datatype, *source, *tag, comm, &status);
    if (*ierror == MPI_SUCCESS) {
        MPI_Status_c2f(&status, status_f);
    }
}

void mpi_waitall_wrapper(int *count, int *array_of_requests_f,
                        int *array_of_statuses_f, int *ierror) {
    MPI_Request *array_of_requests;
    MPI_Status *array_of_statuses;
    array_of_requests = (MPI_Request *)malloc((*count) * sizeof(MPI_Request));
    array_of_statuses = (MPI_Status *)malloc((*count) * sizeof(MPI_Status));
    if (array_of_requests == NULL || array_of_statuses == NULL) {
        *ierror = MPI_ERR_NO_MEM;
        return;
    }
    for (int i = 0; i < *count; i++) {
        array_of_requests[i] = MPI_Request_f2c(array_of_requests_f[i]);
    }

    *ierror = MPI_Waitall(*count, array_of_requests, array_of_statuses);
    for (int i = 0; i < *count; i++) {
        array_of_requests_f[i] = MPI_Request_c2f(array_of_requests[i]);
    }

    for (int i = 0; i < *count; i++) {
        MPI_Status_c2f(&array_of_statuses[i], &array_of_statuses_f[i * MPI_STATUS_SIZE]);
    }

    free(array_of_requests);
    free(array_of_statuses);
}

void mpi_ssend_wrapper(double *buf, int *count, int *datatype_f, int *dest,
                       int *tag, int *comm_f, int *ierror) {
    MPI_Datatype datatype;
    switch (*datatype_f) {
        case 0:
            datatype = MPI_FLOAT;
            break;
        case 1:
            datatype = MPI_DOUBLE;
            break;
        default:
            *ierror = -1;
            return;
    }

    MPI_Comm comm = MPI_Comm_f2c(*comm_f);
    *ierror = MPI_Ssend(buf, *count, datatype, *dest, *tag, comm);
}

void mpi_cart_create_wrapper(int * comm_f, int * ndims, int * dims, int * periods, int * reorder, int * newcomm_f, int * ierror){
    MPI_Comm newcomm = MPI_COMM_NULL;
    MPI_Comm comm = MPI_Comm_f2c(*comm_f);
    *ierror = MPI_Cart_create(comm, *ndims, dims, periods, *reorder, &newcomm);
    *newcomm_f = MPI_Comm_c2f(newcomm);
}

void mpi_cart_coords_wrapper(int * comm_f, int * rank, int * maxdims, int * coords, int * ierror)
{
    MPI_Comm comm = MPI_Comm_f2c(*comm_f);
    *ierror = MPI_Cart_coords(comm, *rank, *maxdims, coords);
}

void mpi_cart_shift_wrapper(int * comm_f, int * dir, int * disp, int * rank_source, int * rank_dest, int * ierror)
{
    MPI_Comm comm = MPI_Comm_f2c(*comm_f);
    *ierror = MPI_Cart_shift(comm, *dir, *disp, rank_source, rank_dest);
}

void mpi_dims_create_wrapper(int * nnodes, int * ndims, int * dims, int * ierror)
{
    *ierror = MPI_Dims_create(*nnodes, *ndims, dims);
}

void mpi_cart_sub_wrapper(int *  comm_f, int * rmains_dims, int * newcomm_f, int * ierror) {
    MPI_Comm comm = MPI_Comm_f2c(*comm_f);
    MPI_Comm newcomm = MPI_COMM_NULL;
    *ierror = MPI_Cart_sub(comm, rmains_dims, &newcomm);
    *newcomm_f = MPI_Comm_c2f(newcomm);
}
