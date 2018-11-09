#include<stdio.h>
#include<mpi.h>

int main( int argc, char **argv)
{

    int world_size;
    int rank;
    int namelen;
    char procname[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv); // Initialize MPI
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(procname, &namelen);
    printf("MPI Hello World from process %d of %d, which is running on machine: %s\n", rank, world_size, procname);

    MPI_Finalize();
}
