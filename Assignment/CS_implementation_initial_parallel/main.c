
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "pgmio.h"

#define WIDTH 192
#define HEIGHT 128
#define PROCS 4
#define PWIDTH WIDTH/PROCS
#define PHEIGHT HEIGHT

#define MAX_ITERS 1000


void checkNumberOfArgs(char *argument, int world_size);

int main(int argc, char const *argv[])
{  
  double master_image[WIDTH][HEIGHT]; // Array for master process
  double image[PWIDTH][PHEIGHT]; // Array for each other process
  double edge[PWIDTH + 2][PHEIGHT + 2];
  double old[PWIDTH + 2][PHEIGHT + 2];
  double new[PWIDTH + 2][PHEIGHT + 2];

  MPI_Init(NULL, NULL); // Initialize MPI

  int world_size; // get world size
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  checkNumberOfArgs(argv[0], world_size);

  int this_rank; // get rank id
  MPI_Comm_rank(MPI_COMM_WORLD, &this_rank);
  
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
  
  /*
  
  
  scattered the data, pick up from here!!!
  
  
  */
  
  

  int i, j = 0;

  // Copying image to edge array
  for(i = 1 ; i < WIDTH + 1; ++i){
    for(j = 1 ; j < PHEIGHT + 1; ++j) {
      edge[i][j] = image[i-1][j-1];
    }
  }
  
  // Initializing old with white pixels
  for(i = 0 ; i < WIDTH + 2; ++i){
    for(j = 0; j < PHEIGHT + 2; ++j){
      old[i][j] = 255.0;
    }
  }

  int iter = 0;
  // Calculating new array
  for(iter = 0; iter < MAX_ITERS; ++iter){
    
    for(i = 1; i < WIDTH + 1; ++i){
      for(j = 1; j < HEIGHT + 1; ++j){
        new[i][j] = 0.25*(old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - edge[i][j]);
      }
    }
    for(i = 1 ; i < WIDTH + 1; ++i){
      for(j = 1; j < HEIGHT + 1; ++j){
        old[i][j] = new[i][j];
      }
    }
  }

  // Copying old/new to image - minus the halo cells
  for(i = 1 ; i < WIDTH + 1; ++i){
    for(j = 1; j < HEIGHT + 1; ++j){
      image[i-1][j-1] = old[i][j];
    }
  }
  
  // Export image to output file
  char *outputfile;
  outputfile = "out.pgm";
  printf("\nWriting %s\n", outputfile);
  pgmwrite(outputfile, image, WIDTH, HEIGHT);
  
  return 0;
}

void checkNumberOfArgs(char *argument, int world_size){
    if(world_size != PROCS){
        fprintf(stderr, "World size must be 4 for %s\n", argument);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}