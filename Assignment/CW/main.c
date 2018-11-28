
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>

#include "pgmio.h"

#define WIDTH 192
#define HEIGHT 128
#define PWIDTH WIDTH/2
#define PHEIGHT HEIGHT/2
#define PROCS 4
#define TRUE 1
#define FALSE 0
#define MAX_ITERS 100000

void checkNumberOfArgs(char *argument, int world_size);
double boundaryval(int i, int m);

int main(int argc, char *argv[])
{
  clock_t start, end;
  start = clock();
  double cpu_time_used;

  double master_image[WIDTH][HEIGHT]; // Array for master process to store initial edge image
  double image[PWIDTH][PHEIGHT]; // Array for each process to store edge images locally
  double old[PWIDTH + 2][PHEIGHT + 2]; // Array for each process used for the calculation
  double new[PWIDTH + 2][PHEIGHT + 2]; // Array for each process used for the calculation
  int width, height; // Image width and height - pixels
  double val;

  MPI_Comm comm2d; // Topology Comm world
  int periodic[2], reorder; // Variables used for the topology creation
  periodic[0] = FALSE; //  Horizontally periodic - set to TRUE LATER
  periodic[1] = FALSE; // Vertically periodic 
  reorder = FALSE; // No reordering of processes in topology
  int this_rank_coords[2];
  int rank_up, rank_down, rank_left, rank_right;

  MPI_Init(NULL, NULL); // Initialize MPI

  int world_size; // Get world size
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  MPI_Comm comm;
  comm = MPI_COMM_WORLD;

  MPI_Status status;
  MPI_Request request;

  //checkNumberOfArgs(argv[0], world_size); // checking the number of arguments given

  int this_rank; // Get rank id
  MPI_Comm_rank(comm, &this_rank);

  int decomp_params[2]; // Used to store the decomposition parameters
  MPI_Dims_create(world_size, 2, &decomp_params[0]); // Funtion to decide how to split image among processes
  MPI_Cart_create(comm, 2, decomp_params, periodic, reorder, &comm2d); // Create Cartesian Topology

  /* Derived Data Type Decleration */
  // Subarrays
  MPI_Datatype array_block;
  MPI_Type_vector(WIDTH/decomp_params[0], HEIGHT/decomp_params[1], HEIGHT, MPI_DOUBLE, &array_block);
  MPI_Type_commit(&array_block);
  // Elements constituting a row in the subarray
  MPI_Datatype row_halo;
  MPI_Type_contiguous(PHEIGHT, MPI_DOUBLE, &row_halo);
  MPI_Type_commit(&row_halo);
  // Elements constituting a column in the subarray
  MPI_Datatype column_halo;
  MPI_Type_vector(PWIDTH, 1, PHEIGHT+2, MPI_DOUBLE, &column_halo);
  MPI_Type_commit(&column_halo);

  if(this_rank == 0){
    printf("~~~~~~~~~\nSuggested M split = %d and N split = %d\n", decomp_params[0], decomp_params[1]);
    char *filename;
    filename = "img/edge192x128.pgm";
    pgmsize(filename, &width, &height);
    printf("~~~~~~~~~\n\"%s\" width = %d and height = %d\n", filename, width, height);
    // Import file to image array
    printf("~~~~~~~~~\nReading \"%s\"\n", filename);
    pgmread(filename, master_image, WIDTH, HEIGHT);

    // Send subarrays to processes
    int w_index, h_index;
    int current_process = 1;
    for(w_index = 0; w_index < decomp_params[0]; ++w_index){
      for(h_index = 0; h_index < decomp_params[1]; ++h_index){
        if((w_index == 0) && (h_index == 0)){
          continue; // Allocate 0,0 locally
        }
        // Change parameters!!! - should be generic
        MPI_Issend(&master_image[w_index*PWIDTH][h_index*PHEIGHT], 1, array_block, current_process++, 0, comm2d, &request);
      }
    }
    // Copy subarray to local buffer (image array)
    int i, j;
    for(i = 0 ; i < WIDTH/2 ; ++i){
      for(j = 0 ; j < HEIGHT/2 ; ++j){
        image[i][j] = master_image[i][j];
      }
    }
  }
  else {
    // Process receives data from 0 and stores it in its local buffer (image array)
    // MPI_Irecv is better performance-wise than MPI_Recv!
    MPI_Irecv(&image[0][0], PWIDTH*PHEIGHT, MPI_DOUBLE, 0, 0, comm2d, &request);
    MPI_Wait(&request, &status);
  }

  // Initialize old subarray with edge image values
  int i, j;
  for(i = 0 ; i < PWIDTH+2 ; ++i){
    for(j = 0 ; j < PHEIGHT+2 ; ++j){
      if((i == 0)||(j == 0)||(i == PWIDTH+1)||(j == PHEIGHT+1)){
        old[i][j] = 255.0; // Boundary elements are set to white
      }
      else{
        old[i][j] = image[i-1][j-1];
      }
    }
  }

  // Rank position in topology & neighbours' IDs
  MPI_Cart_coords(comm2d, this_rank, 2, this_rank_coords);
  MPI_Cart_shift(comm2d, 0, 1, &rank_up, &rank_down);
  MPI_Cart_shift(comm2d, 1, 1, &rank_left, &rank_right);

  printf("I am %d, Coords:  %d and %d, Neighbors: left = %d right = %d up = %d down = %d\n", this_rank, this_rank_coords[0], this_rank_coords[1], rank_left, rank_right, rank_up, rank_down);

  /* SAWTOOTH VALUES - they ruin the image output :P */
  if(rank_up < 0){
    rank_up = MPI_PROC_NULL;

    for (j=1; j < PHEIGHT+1; ++j){
      val = boundaryval(j, PHEIGHT);
      old[0][j] = (int)(255.0*(1.0-val));
    }
  }
  if(rank_down < 0){
    rank_down = MPI_PROC_NULL;

     for (j=1; j < PHEIGHT+1; ++j){
      val = boundaryval(j, PHEIGHT);
      old[PWIDTH+1][j] = (int)(255.0*val);
    }    
  }

  /*Foo print for testing */
  // if(this_rank == 1){
  //   char *procout;
  //   procout = "out_process1.pgm";
  //   printf("~~~~~~~~~\n");
  //   pgmwrite(procout, &old[0][0], PWIDTH+2, PHEIGHT+2);
  // }

  int iter;
  for(iter = 0; iter < MAX_ITERS; ++iter){
    
    /* Halo Swapping */
    MPI_Request request_array[8];
    // periodic
    MPI_Isend(&old[1][1], PHEIGHT, MPI_DOUBLE, rank_up, 1, comm2d, &request_array[0]); // Send to up - tag 1
    MPI_Isend(&old[PWIDTH][1], PHEIGHT, MPI_DOUBLE, rank_down, 2, comm2d, &request_array[1]); // Send to down - tag 2
    // non-periodic
    MPI_Isend(&old[1][1], 1, column_halo, rank_left, 3, comm2d, &request_array[2]); // Send to left - tag 3
    MPI_Isend(&old[1][PHEIGHT], 1, column_halo, rank_right, 4, comm2d, &request_array[3]); // Send to right - tag 4

    // periodic
    MPI_Irecv(&old[0][1], PHEIGHT, MPI_DOUBLE, rank_up, 2, comm2d, &request_array[4]); // Receive from up
    MPI_Irecv(&old[PWIDTH+1][1], PHEIGHT, MPI_DOUBLE, rank_down, 1, comm2d, &request_array[5]); // Receive from down
    // non-periodic
    MPI_Irecv(&old[1][0], 1, column_halo, rank_left, 4, comm2d, &request_array[6]); // Receive from left
    MPI_Irecv(&old[1][PHEIGHT+1], 1, column_halo, rank_right, 3, comm2d, &request_array[7]); // Receiver from right
    
    // Non-boundary element depandant calculations
    for(i = 2; i < PWIDTH; ++i){
        for(j = 2; j < PHEIGHT; ++j){
          new[i][j] = 0.25*(old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - image[i-1][j-1]);
        }
    }

    MPI_Status status_array[8];
    MPI_Waitall(8, request_array, status_array);

    // Boundary element depandant calculations
    for(i = 1; i <= PWIDTH; i += PWIDTH-1){
        for(j = 1; j <= PHEIGHT; ++j){
          new[i][j] = 0.25*(old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - image[i-1][j-1]);
        }
    }

    for(i = 2; i < PWIDTH; ++i){
        for(j = 1; j <= PHEIGHT; j += PHEIGHT-1){
          new[i][j] = 0.25*(old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - image[i-1][j-1]);
        }
    }

    for(i = 1 ; i <= PWIDTH; ++i){
      for(j = 1; j <= PHEIGHT; ++j){
        old[i][j] = new[i][j];
      }
    }
  }

  // Copy old to image buffer in order to send to master process
  for (i=1; i <= PWIDTH; ++i){
    for (j=1; j <= PHEIGHT; ++j){
      image[i-1][j-1] = old[i][j];
    }
  }
  /* Gather results and print output */
  if(this_rank == 0){
    MPI_Request request_array[PROCS-1];
    MPI_Status status_array[PROCS-1];

    int w_index, h_index;
    int current_process = 1;
    for(w_index = 0; w_index < decomp_params[0]; w_index++){
      for(h_index = 0; h_index < decomp_params[1]; h_index++){
        if((w_index == 0) && (h_index == 0)){
          continue;
        }
        // MPI_Recv(&master_image[w_index*PWIDTH][h_index*PHEIGHT], 1, array_block, current_process++, 0, comm2d, &status_array[current_process-1]);
        MPI_Irecv(&master_image[w_index*PWIDTH][h_index*PHEIGHT], 1, array_block, current_process++, 0, comm2d, &request_array[current_process-1]);
      }
    }
    MPI_Waitall(PROCS-1, request_array, status_array);
    // Copy own subarray
    for (i = 0; i < PWIDTH; ++i){
      for (j = 0; j < PHEIGHT; ++j){
        master_image[i][j] = image[i][j];
      }
    }
    // Write final image to file
    char *outputfile;
    outputfile = "out.pgm";
    printf("~~~~~~~~~\n");
    pgmwrite(outputfile, &master_image[0][0], WIDTH, HEIGHT);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("~~~~~~~~~\nFinished %d iterations in %f seconds\n", MAX_ITERS, cpu_time_used);
  }
  else {
    // Send subarray to master process
    MPI_Issend(&image[0][0], PWIDTH*PHEIGHT, MPI_DOUBLE, 0, 0, comm2d, &request);
  }
  MPI_Finalize();
}

void checkNumberOfArgs(char *argument, int world_size){
    if(world_size != PROCS){
        fprintf(stderr, "World size must be 4 for %s\n", argument);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

double boundaryval(int i, int m){
  double val;

  val = 2.0*((double)(i-1))/((double)(m-1));
  if (i >= m/2+1) val = 2.0-val;
  
  return val;
}