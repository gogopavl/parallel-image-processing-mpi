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

    MPI_Init(&argc, &argv); // Initialize MPI
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Get_processor_name(procname, &namelen);

    printf("MPI Hello World from process %d of %d, which is running on machine: %s\n", rank, size, procname);

    int localNumber = rank+10;
    int personalCounter = 0;
    int buf;

    MPI_Request request;
    MPI_Status status;
    MPI_Comm comm;
    comm = MPI_COMM_WORLD;

    int i;

    for(i = 0 ; i < 4 ; i++)
    {
        if(rank == 0)
        {
            printf("\n\nIteration %d\n\n",i);
            printf("Sum should be: %d\n", 10+11+12+13);
            // receive from 3
            // send to 1
            printf("%d num = %d\n", rank, localNumber);

            MPI_Issend(&localNumber, 1, MPI_INT, 1, 0, comm, &request);
            MPI_Recv(&buf, 1, MPI_INT, 3, 0, comm, &status);
            localNumber = buf;
            personalCounter += buf;
            printf("I am process %d and I received number %d and my current sum is: %d\n", rank, buf, personalCounter);
            MPI_Wait(&request, &status);
        }
        else if (rank == 1)
        {
            // receive from 0
            // send to 2
            printf("%d num = %d\n", rank, localNumber);
            MPI_Issend(&localNumber, 1, MPI_INT, 2, 0, comm, &request);
            MPI_Recv(&buf, 1, MPI_INT, 0, 0, comm, &status);
            localNumber = buf;
            personalCounter += buf;
            printf("I am process %d and I received number %d and my current sum is: %d\n", rank, buf, personalCounter);
            MPI_Wait(&request, &status);
        }
        else if (rank == 2)
        {
            // receive from 1
            // send to 3
            printf("%d num = %d\n", rank, localNumber);
            MPI_Issend(&localNumber, 1, MPI_INT, 3, 0, comm, &request);
            MPI_Recv(&buf, 1, MPI_INT, 1, 0, comm, &status);
            localNumber = buf;
            personalCounter += buf;
            printf("I am process %d and I received number %d and my current sum is: %d\n", rank, buf, personalCounter);
            MPI_Wait(&request, &status);
        }
        else if (rank == 3)
        {
            // receive from 2
            // send to 0
            printf("%d num = %d\n", rank, localNumber);
            MPI_Issend(&localNumber, 1, MPI_INT, 0, 0, comm, &request);
            MPI_Recv(&buf, 1, MPI_INT, 2, 0, comm, &status);
            localNumber = buf;
            personalCounter += buf;
            printf("I am process %d and I received number %d and my current sum is: %d\n", rank, buf, personalCounter);
            MPI_Wait(&request, &status);
        }
    }


    MPI_Finalize();

}
