#include<stdio.h>
#include<mpi.h>

int main( int argc, char **argv)
{

    int size;
    int rank;
    int namelen;
    char procname[MPI_MAX_PROCESSOR_NAME];

    int numberOfIterations = argc;

    int counter = 0;
    for( counter = 0 ; counter < argc ; counter++)
    {
        printf("Arg number %d: %s\n", counter, argv[counter]);

    }

    printf("Num of arguments is: %d\n", numberOfIterations);

    printf("N is: %d\n", numberOfIterations);

    MPI_Init(&argc, &argv); // Initialize MPI
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(procname, &namelen);

    printf("MPI Hello World from process %d of %d, which is running on machine: %s\n", rank, size, procname);


    MPI_Finalize();

}
