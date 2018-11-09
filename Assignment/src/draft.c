#include<stdio.h>
#include<mpi.h>

void printError(char *argument);

int main( int argc, char **argv)
{

    MPI_Init(NULL, NULL); // Initialize MPI

    int this_rank; // get rank id
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);

    int world_size; // get world size
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);


    if(world_size < 2){
        printError(argv[0]);
    }

    int number;
    printf("in");
    if (this_rank == 0){
        number = -1;
        MPI_Send(&number, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    } else if(this_rank == 1){
        MPI_Recv(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Process 1 received number %d from process 0\n", number);
    }
    MPI_Finalize();
}

void printError(char *argument){
    fprintf(stderr, "World size must be greater than 1 for %s\n", argument);
    MPI_Abort(MPI_COMM_WORLD, 1);
}
