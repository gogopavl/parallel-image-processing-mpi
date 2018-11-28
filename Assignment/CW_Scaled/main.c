
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include "pgmio.h"

#define PROCS 4
#define TRUE 1
#define FALSE 0
#define MAX_ITERS 200000

void checkNumberOfArgs(char *argument, int world_size);
double boundaryval(int i, int m);

int main(int argc, char *argv[])
{
  clock_t start, end;
  start = clock();
  double cpu_time_used;

  char *filename;
  filename = "img/edgenew192x128.pgm";
  int dimensions_array[4];

  MPI_Init(NULL, NULL); // Initialize MPI

  MPI_Comm comm;
  comm = MPI_COMM_WORLD;

  int world_size; // Get world size
  MPI_Comm_size(comm, &world_size);

  MPI_Status status;
  MPI_Request request;

  //checkNumberOfArgs(argv[0], world_size); // checking the number of arguments given

  int this_rank; // Get rank id
  MPI_Comm_rank(comm, &this_rank);

  int width, height; // Image width and height - pixels
  int pwidth, pheight; // Process width and height dimensions  
  int decomp_params[2]; // Used to store the decomposition parameters

  if(this_rank == 0){

    MPI_Dims_create(world_size, 2, &decomp_params[0]); // Funtion to decide how to split image among processes

    printf("~~~~~~~~~\nSuggested M (width) split = %d and N (height) split = %d\n", decomp_params[0], decomp_params[1]);
    pgmsize(filename, &width, &height);
    
    pwidth = width/decomp_params[0];
    pheight = height/decomp_params[1];    

    printf("~~~~~~~~~\n\"%s\" width = %d and height = %d\n", filename, width, height);
    printf("pwidth = %d and pheight = %d\n", pwidth, pheight);
  }

  double master_image[width][height]; // Array for master process to store initial edge image

  if(this_rank == 0){
    printf("~~~~~~~~~\nReading \"%s\"\n", filename);
    pgmread(filename, master_image, width, height);
    
    dimensions_array[0] = width;
    dimensions_array[1] = height;
    dimensions_array[2] = pwidth;
    dimensions_array[3] = pheight;
    int process_iterator;
    for(process_iterator = 1; process_iterator < world_size; ++process_iterator){
      MPI_Issend(&dimensions_array, 4, MPI_INT, process_iterator, 0, comm, &request);
    }
  } 
  else {
    MPI_Irecv(&dimensions_array[0], 4, MPI_INT, 0, 0, comm, &request);
    MPI_Wait(&request, &status);
    width = dimensions_array[0];
    height = dimensions_array[1];
    pwidth = dimensions_array[2];
    pheight = dimensions_array[3];
  }

  double image[pwidth][pheight]; // Array for each process to store edge images locally
  double old[pwidth + 2][pheight + 2]; // Array for each process used for the calculation
  double new[pwidth + 2][pheight + 2]; // Array for each process used for the calculation
  double val;

  MPI_Comm comm2d; // Topology Comm world
  int periodic[2], reorder; // Variables used for the topology creation
  periodic[0] = FALSE; //  Horizontally periodic - set to TRUE LATER
  periodic[1] = TRUE; // Vertically periodic
  reorder = FALSE; // No reordering of processes in topology
  int this_rank_coords[2];
  int rank_up, rank_down, rank_left, rank_right;

  decomp_params[0] = width/pwidth;
  decomp_params[1] = height/pheight;

  MPI_Cart_create(comm, 2, decomp_params, periodic, reorder, &comm2d); // Create Cartesian Topology
  

  /* Derived Data Type Decleration */
  // Subarrays
  MPI_Datatype array_block;
  MPI_Type_vector(pwidth, pheight, height, MPI_DOUBLE, &array_block);
  MPI_Type_commit(&array_block);
  // Elements constituting a row in the subarray
  MPI_Datatype row_halo;
  MPI_Type_contiguous(pheight, MPI_DOUBLE, &row_halo);
  MPI_Type_commit(&row_halo);
  // Elements constituting a column in the subarray
  MPI_Datatype column_halo;
  MPI_Type_vector(pwidth, 1, pheight+2, MPI_DOUBLE, &column_halo);
  MPI_Type_commit(&column_halo);

  printf("I am %d %d %d\n", this_rank, pwidth, pheight);

  if(this_rank == 0){
    // Send subarrays to processes
    int w_index, h_index;
    int current_process = 1;

    for(w_index = 0; w_index < decomp_params[0]; ++w_index){
      for(h_index = 0; h_index < decomp_params[1]; ++h_index){
        if((w_index == 0) && (h_index == 0)){
          continue; // Allocate 0,0 locally
        }
        // Change parameters!!! - should be generic
        MPI_Issend(&master_image[w_index*pwidth][h_index*pheight], 1, array_block, current_process++, 0, comm2d, &request);
      }
    }
    // Copy subarray to local buffer (image array)
    int i, j;
    for(i = 0 ; i < pwidth ; ++i){
      for(j = 0 ; j < pheight ; ++j){
        image[i][j] = master_image[i][j];
      }
    }
  }
  else {
    // Process receives data from 0 and stores it in its local buffer (image array)
    // MPI_Irecv is better performance-wise than MPI_Recv!
    MPI_Irecv(&image[0][0], pwidth*pheight, MPI_DOUBLE, 0, 0, comm2d, &request);
    MPI_Wait(&request, &status);
  }

  // Initialize old subarray with edge image values
  int i, j;
  for(i = 0 ; i < pwidth+2 ; ++i){
    for(j = 0 ; j < pheight+2 ; ++j){
        old[i][j] = 255.0; // Elements are set to white
    }
  }

  // Rank position in topology & neighbours' IDs
  MPI_Cart_coords(comm2d, this_rank, 2, this_rank_coords);
  MPI_Cart_shift(comm2d, 0, 1, &rank_up, &rank_down);
  MPI_Cart_shift(comm2d, 1, 1, &rank_left, &rank_right);

  printf("I am %d, Coords:  %d and %d, Neighbors: left = %d right = %d up = %d down = %d\n", this_rank, this_rank_coords[0], this_rank_coords[1], rank_left, rank_right, rank_up, rank_down);

  /* SAWTOOTH VALUES - they ruin the image output :P */
  // My index should be global!!!
  if(rank_up < 0){
    rank_up = MPI_PROC_NULL;

    for (j=1; j < pheight+1; ++j){
      val = boundaryval(j, height);
      old[0][j] = (int)(255.0*(1.0-val));
    }
  }
  if(rank_down < 0){
    rank_down = MPI_PROC_NULL;

     for (j=1; j < pheight+1; ++j){
      val = boundaryval(j, height);
      old[pwidth+1][j] = (int)(255.0*val);
    }
  }

  int iter;
  for(iter = 0; iter < MAX_ITERS; ++iter){

    /* Halo Swapping */
    MPI_Request request_array[8];
    // periodic
    MPI_Isend(&old[1][1], pheight, MPI_DOUBLE, rank_up, 1, comm2d, &request_array[0]); // Send to up - tag 1
    MPI_Isend(&old[pwidth][1], pheight, MPI_DOUBLE, rank_down, 2, comm2d, &request_array[1]); // Send to down - tag 2
    // non-periodic
    MPI_Isend(&old[1][1], 1, column_halo, rank_left, 3, comm2d, &request_array[2]); // Send to left - tag 3
    MPI_Isend(&old[1][pheight], 1, column_halo, rank_right, 4, comm2d, &request_array[3]); // Send to right - tag 4

    // periodic
    MPI_Irecv(&old[0][1], pheight, MPI_DOUBLE, rank_up, 2, comm2d, &request_array[4]); // Receive from up
    MPI_Irecv(&old[pwidth+1][1], pheight, MPI_DOUBLE, rank_down, 1, comm2d, &request_array[5]); // Receive from down
    // non-periodic
    MPI_Irecv(&old[1][0], 1, column_halo, rank_left, 4, comm2d, &request_array[6]); // Receive from left
    MPI_Irecv(&old[1][pheight+1], 1, column_halo, rank_right, 3, comm2d, &request_array[7]); // Receiver from right

    // Non-boundary element depandant calculations
    for(i = 2; i < pwidth; ++i){
        for(j = 2; j < pheight; ++j){
          new[i][j] = 0.25*(old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - image[i-1][j-1]);
        }
    }

    MPI_Status status_array[8];
    MPI_Waitall(8, request_array, status_array);

    // Boundary element depandant calculations
    for(i = 1; i <= pwidth; i += pwidth-1){
        for(j = 1; j <= pheight; ++j){
          new[i][j] = 0.25*(old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - image[i-1][j-1]);
        }
    }

    for(i = 2; i < pwidth; ++i){
        for(j = 1; j <= pheight; j += pheight-1){
          new[i][j] = 0.25*(old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - image[i-1][j-1]);
        }
    }

    for(i = 1 ; i <= pwidth; ++i){
      for(j = 1; j <= pheight; ++j){
        old[i][j] = new[i][j];
      }
    }
  }

  // Copy old to image buffer in order to send to master process
  for (i=1; i <= pwidth; ++i){
    for (j=1; j <= pheight; ++j){
      image[i-1][j-1] = old[i][j];
    }
  }
  /* Gather results and print output */
  if(this_rank == 0){
    MPI_Request request_array[world_size-1];
    MPI_Status status_array[world_size-1];

    int w_index, h_index;
    int current_process = 1;
    for(w_index = 0; w_index < decomp_params[0]; w_index++){
      for(h_index = 0; h_index < decomp_params[1]; h_index++){
        if((w_index == 0) && (h_index == 0)){
          continue;
        }
        // MPI_Recv(&master_image[w_index*pwidth][h_index*pheight], 1, array_block, current_process++, 0, comm2d, &status_array[current_process-1]);
        MPI_Irecv(&master_image[w_index*pwidth][h_index*pheight], 1, array_block, current_process++, 0, comm2d, &request_array[current_process-1]);
      }
    }
    MPI_Waitall(world_size-1, request_array, status_array);
    // Copy own subarray
    for (i = 0; i < pwidth; ++i){
      for (j = 0; j < pheight; ++j){
        master_image[i][j] = image[i][j];
      }
    }
    // Write final image to file
    char *outputfile;
    outputfile = "out.pgm";
    printf("~~~~~~~~~\n");
    pgmwrite(outputfile, &master_image[0][0], width, height);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("~~~~~~~~~\nFinished %d iterations in %f seconds\n", MAX_ITERS, cpu_time_used);
  }
  else {
    // Send subarray to master process
    MPI_Issend(&image[0][0], pwidth*pheight, MPI_DOUBLE, 0, 0, comm2d, &request);
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
