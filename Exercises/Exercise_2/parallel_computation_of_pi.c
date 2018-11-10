#include<stdio.h>
#include<stdlib.h>
#include <math.h>
#include <mpi.h>

int main( int argc, char **argv)
{

    int size;
    int rank;
    int namelen;
    char procname[MPI_MAX_PROCESSOR_NAME];

    int numberOfIterations = atoi(argv[1]);

    printf("N is: %d\n", numberOfIterations);

    int i = 0;
    double piApprox;
    double sumFactor;

    for(i = 1 ; i <= numberOfIterations ; i++)
    {
      sumFactor += 1/(1 + pow( ((i - 0.5)/numberOfIterations) , 2.0));
    }

    printf("Pi factor = %lf\n", sumFactor);

    piApprox = (sumFactor/numberOfIterations)*4;

    printf("Pi approximation with N = %d is: %lf\n", numberOfIterations, piApprox);

    MPI_Init(&argc, &argv); // Initialize MPI
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(procname, &namelen);

    // printf("MPI Hello World from process %d of %d, which is running on machine: %s\n", rank, size, procname);


    MPI_Finalize();

}
