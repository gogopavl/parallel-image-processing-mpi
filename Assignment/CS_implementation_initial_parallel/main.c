
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

#define MAX_ITERS 500000


void checkNumberOfArgs(char *argument, int world_size);

int main(int argc, char *argv[])
{  
  clock_t start, end;
  start = clock();
  double cpu_time_used;

  double master_image[WIDTH][HEIGHT]; // Array for master process
  double image[WIDTH][HEIGHT]; // Array for each other process
  double edge[PWIDTH + 2][PHEIGHT + 2];
  double old[PWIDTH + 2][PHEIGHT + 2];
  double new[PWIDTH + 2][PHEIGHT + 2];

  MPI_Init(NULL, NULL); // Initialize MPI

  int world_size; // get world size
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  checkNumberOfArgs(argv[0], world_size);

  int this_rank; // get rank id
  MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);

  int next_rank, previous_rank;

  next_rank = this_rank + 1;
  previous_rank = this_rank - 1;
  if(next_rank >= world_size){
    next_rank = MPI_PROC_NULL; // Far right and far left processes do not need to communicate sideways
  }
  if(previous_rank < 0){
    previous_rank = MPI_PROC_NULL;
  }
  
  if(this_rank == 0){
    // Declaring input file
    char *filename;
    filename = "edge192x128.pgm";

    // Import file to image array
    printf("\nReading <%s>\n", filename);
    pgmread(filename, master_image, WIDTH, HEIGHT);
    printf("\n");
  }

  MPI_Scatter(&master_image[0][0], PWIDTH*PHEIGHT, MPI_DOUBLE, &image[0][0], PWIDTH*PHEIGHT, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  int i, j = 0;

  // Copying image to edge array
  for(i = 1 ; i < PWIDTH + 1; ++i){
    for(j = 1 ; j < PHEIGHT + 1; ++j) {
      edge[i][j] = image[i-1][j-1];
    }
  }
  
  // Initializing old with white pixels
  for(i = 0 ; i < PWIDTH + 2; ++i){
    for(j = 0; j < PHEIGHT + 2; ++j){
      old[i][j] = 255.0;
    }
  }

  int iter = 0;
  // Calculating new array
  for(iter = 0; iter < MAX_ITERS; ++iter){

    MPI_Sendrecv(&old[PWIDTH][1], PHEIGHT, MPI_DOUBLE, next_rank, 1, &old[0][1], PHEIGHT, MPI_DOUBLE, previous_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&old[1][1], PHEIGHT, MPI_DOUBLE, previous_rank, 2, &old[PWIDTH +1 ][1], PHEIGHT, MPI_DOUBLE, next_rank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for(i = 1; i < PWIDTH + 1; ++i){
      for(j = 1; j < PHEIGHT + 1; ++j){
        new[i][j] = 0.25*(old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - edge[i][j]);
      }
    }
    for(i = 1 ; i < PWIDTH + 1; ++i){
      for(j = 1; j < PHEIGHT + 1; ++j){
        old[i][j] = new[i][j];
      }
    }
  }

  // Copying old/new to image - minus the halo cells
  for(i = 1 ; i < PWIDTH + 1; ++i){
    for(j = 1; j < PHEIGHT + 1; ++j){
      image[i-1][j-1] = old[i][j];
    }
  }

  // MPI GATHER
  MPI_Gather(&image[0][0], PWIDTH*PHEIGHT, MPI_DOUBLE, &master_image[0][0], PWIDTH*PHEIGHT, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  if(this_rank == 0){
    // Export image to output file
    char *outputfile;
    outputfile = "out.pgm";
    printf("\nWriting %s\n", outputfile);
    pgmwrite(outputfile, &master_image[0][0], WIDTH, HEIGHT);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Finished %d iterations in %f seconds\n", MAX_ITERS, cpu_time_used);
  }
  MPI_Finalize();
 
}

void checkNumberOfArgs(char *argument, int world_size){
    if(world_size != PROCS){
        fprintf(stderr, "World size must be 4 for %s\n", argument);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}