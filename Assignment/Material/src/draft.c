#include <stdlib.h>
#include<stdio.h>
#include<mpi.h>

void checkNumberOfArgs(char *argument, int world_size);

int main( int argc, char **argv)
{
    MPI_Init(NULL, NULL); // Initialize MPI

    int world_size; // get world size
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    checkNumberOfArgs(argv[0], world_size);

    int this_rank; // get rank id
    MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);
    

    int token;

     if (this_rank != 0) {
        MPI_Recv(&token, 1, MPI_INT, this_rank -1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Process %d received token %d from process %d\n", this_rank, token, this_rank-1);
    } else {
        token = 55;
    } 
    MPI_Send(&token, 1, MPI_INT, (this_rank + 1) % world_size, 0, MPI_COMM_WORLD);

    if (this_rank == 0){
        MPI_Recv(&token, 1, MPI_INT, world_size - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Process %d received token %d from process %d\n", this_rank, token, world_size - 1);
    }


    MPI_Finalize();
}

void checkNumberOfArgs(char *argument, int world_size){
    if(world_size < 2){
        fprintf(stderr, "World size must be greater than 1 for %s\n", argument);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}
