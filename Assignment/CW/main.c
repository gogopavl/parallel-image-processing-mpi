
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

#include "pgmio.h"

#define WIDTH 192
#define HEIGHT 128
#define PROCS 4
#define PWIDTH WIDTH/PROCS
#define PHEIGHT HEIGHT
#define TRUE 1
#define FALSE 0
#define MAX_ITERS 1500

void checkNumberOfArgs(char *argument, int world_size);

int main(int argc, char *argv[])
{  
  clock_t start, end;
  start = clock();
  double cpu_time_used;

  MPI_Comm comm2d;
  int periodic[2], reorder;
  periodic[0] = TRUE;
  periodic[1] = FALSE;
  reorder = FALSE;
  int coord[2], id;
  int rank_up, rank_down, rank_left, rank_right;

  double master_image[WIDTH][HEIGHT]; // Array for master process
  double image[WIDTH][HEIGHT]; // Array for each other process
  double edge[PWIDTH + 2][PHEIGHT + 2];
  double old[PWIDTH + 2][PHEIGHT + 2];
  double new[PWIDTH + 2][PHEIGHT + 2];
  int width, height;

  MPI_Init(NULL, NULL); // Initialize MPI
  
  int world_size; // get world size
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  MPI_Comm comm;
  comm = MPI_COMM_WORLD;

  //checkNumberOfArgs(argv[0], world_size); // checking the number of arguments given

  int this_rank; // get rank id
  MPI_Comm_rank(comm, &this_rank);

  int decomp_params[2]; // used to store the decomposition parameters
  MPI_Dims_create(world_size, 2, &decomp_params[0]); // funtion to decide how to split image among processes
  
  if(this_rank == 0){
    printf("~~~~~~~~~\nSuggested M split = %d and N split = %d\n", decomp_params[0], decomp_params[1]);
    char *filename;
    filename = "img/edge192x128.pgm";    
    pgmsize(filename, &width, &height);
    printf("~~~~~~~~~\n\"%s\" width = %d and height = %d\n", filename, width, height);
    // Import file to image array
    printf("~~~~~~~~~\nReading \"%s\"\n", filename);
    pgmread(filename, master_image, WIDTH, HEIGHT);
  }
  MPI_Cart_create(comm, 2, decomp_params, periodic, reorder, &comm2d);
  
  if(this_rank == 0){
    printf("0     1\n2     3\n");
    MPI_Cart_coords(comm2d, this_rank, 2, coord);
    printf("I am %d and my coords are: %d and %d\n", this_rank, coord[0], coord[1]);
    MPI_Cart_shift(comm2d, 0, 1, &rank_up, &rank_down);
    MPI_Cart_shift(comm2d, 1, 1, &rank_left, &rank_right);
    printf("Neighbors: left = %d right = %d up = %d down = %d\n",rank_left, rank_right, rank_up, rank_down);
  }

  if(this_rank == 0){
    // Export image to output file
    char *outputfile;
    outputfile = "out.pgm";
    printf("~~~~~~~~~\n");
    pgmwrite(outputfile, &master_image[0][0], WIDTH, HEIGHT);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("~~~~~~~~~\nFinished %d iterations in %f seconds\n", MAX_ITERS, cpu_time_used);
  }
  MPI_Finalize();
 
}

void checkNumberOfArgs(char *argument, int world_size){
    if(world_size != PROCS){
        fprintf(stderr, "World size must be 4 for %s\n", argument);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}
