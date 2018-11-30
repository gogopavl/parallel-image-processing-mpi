
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include "pgmio.h"

#define TRUE 1
#define FALSE 0
#define MAX_ITERS 100000
#define CHECK_FREQ 500
#define THRESHOLD 0.1

void check_num_of_args(int argc);
double calculate_delta(double new, double old);
double boundaryval(int i, int m);

int main(int argc, char *argv[])
{
  clock_t start, iotstamp, iterstart, iterend, niterstart, niterend, writeend, end;
  start = clock();
  double sum_of_niters; // Variable used to calculate the mean runtime of the nested iteration

  MPI_Init(&argc, &argv); // Initialize MPI

  int world_size; // Get world size
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  MPI_Comm comm;
  comm = MPI_COMM_WORLD; // Shorten the name of MPI_COMM_WORLD

  MPI_Status status;
  MPI_Request request;
  
  int this_rank; // Get the current rank's id
  MPI_Comm_rank(comm, &this_rank);

  if(this_rank == 0){
    check_num_of_args(argc); // Checking the number of arguments given
  }

  char *filename;
  filename = argv[1]; // Given argument is the path leading to the file

  int width, height; // Image width and height - pixels
  int pwidth, pheight; // Process width and height - pixels  
  int decomp_params[2]; // Used to store the decomposition parameters
  int dimensions_array[4]; // Used to broadcast image dimensions to processes

  // Initialised to 0 to suppress Dims_Create null dimension error
  decomp_params[0] = 0;
  decomp_params[1] = 0; 
  width, height = 0;
  MPI_Dims_create(world_size, 2, decomp_params); // Funtion to decide how to split image
  
  if(this_rank == 0){ // Get image dimensions
    pgmsize(filename, &width, &height);
    pwidth = width/decomp_params[0];
    pheight = height/decomp_params[1];    
  }

  double master_image[width][height]; // Array for master process to store initial edge image

  if(this_rank == 0){ // Pack values for Broadcast
    pgmread(filename, master_image, width, height);
    dimensions_array[0] = width;
    dimensions_array[1] = height;
    dimensions_array[2] = pwidth;
    dimensions_array[3] = pheight;
  }
  
  MPI_Bcast(dimensions_array, 4, MPI_INT, 0, comm); // Broadcasting image and process dimension - pixels

  if(this_rank != 0) { // Unpacking values from Broadcast
    width = dimensions_array[0];
    height = dimensions_array[1];
    pwidth = dimensions_array[2];
    pheight = dimensions_array[3];
  }
  
  double image[pwidth][pheight]; // Array for each process to store edge images locally - buffer
  double old[pwidth + 2][pheight + 2]; // Array for each process used for the calculation - + 2 on each dimension for halo swapping
  double new[pwidth + 2][pheight + 2]; // Array for each process used for the calculation
  double val, delta; // Val used for the sawtooth values calculation and delta for the pixel difference

  MPI_Comm comm2d; // Topology communicator
  int periodic[2], reorder; // Variables used for the topology creation
  periodic[0] = FALSE; //  Fixed values on the image's sides
  periodic[1] = TRUE; // Periodic across long side
  reorder = FALSE; // No reordering of processes in topology
  int this_rank_coords[2]; // To unpack rank's position values - 
  int this_row, this_col, rank_up, rank_down, rank_left, rank_right; // Used for process position in topology and neighbours - Cart coords and shift


  MPI_Cart_create(comm, 2, decomp_params, periodic, reorder, &comm2d); // Create Cartesian Topology

  /* Derived Data Type Decleration */
  // Subarrays
  MPI_Datatype array_block;
  MPI_Type_vector(pwidth, pheight, height, MPI_DOUBLE, &array_block);
  MPI_Type_commit(&array_block);
  // Elements constituting a column in the subarray
  MPI_Datatype column_halo;
  MPI_Type_vector(pwidth, 1, pheight+2, MPI_DOUBLE, &column_halo);
  MPI_Type_commit(&column_halo);

  // Rank position in topology & neighbours' IDs
  MPI_Cart_coords(comm2d, this_rank, 2, this_rank_coords);
  MPI_Cart_shift(comm2d, 0, 1, &rank_up, &rank_down);
  MPI_Cart_shift(comm2d, 1, 1, &rank_left, &rank_right);
  this_row = this_rank_coords[0];
  this_col = this_rank_coords[1];

  /* Image distribution across processes */
  if(this_rank == 0){
    // Send subarrays to processes
    int w_index, h_index;
    int current_process = 1;

    for(w_index = 0; w_index < decomp_params[0]; ++w_index){
      for(h_index = 0; h_index < decomp_params[1]; ++h_index){
        if((w_index == 0) && (h_index == 0)){
          continue; // Allocate 0,0 locally
        }
        MPI_Isend(&master_image[w_index*pwidth][h_index*pheight], 1, array_block, current_process++, 0, comm2d, &request); // Sending array block custom type
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
    MPI_Irecv(&image[0][0], pwidth*pheight, MPI_DOUBLE, 0, 0, comm2d, &request);
    MPI_Wait(&request, &status); // Each process has to definitely receive the sub-array before progressing to the next step
  }
  iotstamp = clock(); // Timestamp for I/O completion
  
  /* Buffer initializations */
  // Initialize old subarray - all processes
  int i, j;
  for(i = 0 ; i < pwidth+2 ; ++i){
    for(j = 0 ; j < pheight+2 ; ++j){
        old[i][j] = 255.0; // Elements are set to white
    }
  }

  // Computing sawtooth values
  if(rank_up < 0){ // MPI_PROC_NULL
    for (j=1; j < pheight+1; ++j){
      val = boundaryval(j + this_col*pheight, height); // Global position of column element
      old[0][j] = (int)(255.0*(1.0-val));
    }
  }
  if(rank_down < 0){ // MPI_PROC_NULL
     for (j=1; j < pheight+1; ++j){
      val = boundaryval(j + this_col*pheight, height); // Global position of column element
      old[pwidth+1][j] = (int)(255.0*val);
    }
  }

  /* Computation start */
  iterstart = clock(); // Timestamp for iteration start - all iterations
  int iter, calc_delta;
  calc_delta = FALSE; // Flag to calculate deltas - difference in values
  for(iter = 1; iter <= MAX_ITERS; ++iter){
    niterstart = clock(); // Timestamp for nested iteration start - single iteration
    if(iter % CHECK_FREQ == 0){ // Trigger flag
      calc_delta = TRUE;
    }
    else {
      calc_delta = FALSE;
    }

    /* Halo Swapping */
    MPI_Request request_array[8]; // Array to store requests from halo swapping send and receives
    // Sending
    MPI_Isend(&old[1][1], pheight, MPI_DOUBLE, rank_up, 1, comm2d, &request_array[0]); // Send to up - tag 1
    MPI_Isend(&old[pwidth][1], pheight, MPI_DOUBLE, rank_down, 2, comm2d, &request_array[1]); // Send to down - tag 2
    MPI_Isend(&old[1][1], 1, column_halo, rank_left, 3, comm2d, &request_array[2]); // Send to left - tag 3
    MPI_Isend(&old[1][pheight], 1, column_halo, rank_right, 4, comm2d, &request_array[3]); // Send to right - tag 4
    // Receiving
    MPI_Irecv(&old[0][1], pheight, MPI_DOUBLE, rank_up, 2, comm2d, &request_array[4]); // Receive from up
    MPI_Irecv(&old[pwidth+1][1], pheight, MPI_DOUBLE, rank_down, 1, comm2d, &request_array[5]); // Receive from down
    MPI_Irecv(&old[1][0], 1, column_halo, rank_left, 4, comm2d, &request_array[6]); // Receive from left
    MPI_Irecv(&old[1][pheight+1], 1, column_halo, rank_right, 3, comm2d, &request_array[7]); // Receiver from right

    
    /* Boundary independent element calculations */
    for(i = 2; i < pwidth; ++i){
        for(j = 2; j < pheight; ++j){
          new[i][j] = 0.25*(old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - image[i-1][j-1]);
          if(calc_delta == TRUE){
            delta += calculate_delta(new[i][j], old[i][j]);
          }
        }
    }

    // Make sure all messages have been delivered and received before moving on
    MPI_Status status_array[8];
    MPI_Waitall(8, request_array, status_array);

    /* Boundary dependent element calculations */
    for(i = 1; i <= pwidth; i += pwidth-1){
        for(j = 1; j <= pheight; ++j){
          new[i][j] = 0.25*(old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - image[i-1][j-1]);
          if(calc_delta == TRUE){
            delta += calculate_delta(new[i][j], old[i][j]);
          }
        }
    }

    for(i = 2; i < pwidth; ++i){
        for(j = 1; j <= pheight; j += pheight-1){
          new[i][j] = 0.25*(old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - image[i-1][j-1]);
          if(calc_delta == TRUE){
            delta += calculate_delta(new[i][j], old[i][j]);
          }
        }
    }
    // Copying values for next iteration
    for(i = 1 ; i <= pwidth; ++i){
      for(j = 1; j <= pheight; ++j){
        old[i][j] = new[i][j];
      }
    }

    /* Calculation of average pixel values and deltas */
    double pixel_val, avg_pixel_val; // Used for average pixel calculation
    double global_delta; // Maximum delta across all processes
    if(calc_delta == TRUE){
      pixel_val = 0.0;
      for(i = 1; i <= pwidth; ++i){
        for(j = 1; j <= pheight; ++j){
          pixel_val += new[i][j];
        }
      }
      pixel_val /= (pwidth*pheight); // Average pixel value in this process
      MPI_Reduce(&pixel_val, &avg_pixel_val, 1, MPI_DOUBLE, MPI_SUM, 0, comm2d); // Process 0 issues a Reduce to get all average pixel values
            
      delta /= (pwidth*pheight); // Average delta value in this process
      MPI_Allreduce(&delta, &global_delta, 1, MPI_DOUBLE, MPI_MAX, comm2d); // Every process gets the maximum delta in the communicator - break condition is the same for all
      
      if(this_rank == 0){
        avg_pixel_val /= (double)world_size; // Process 0 calculates the average value from all averages
        printf("Average pixel value = %f on iteration %d\n", avg_pixel_val, iter);
        printf("Max delta value = %f on iteration %d\n", global_delta, iter);    
      }

      if(global_delta < THRESHOLD){
        break;
      }
    }
    niterend = clock(); // Timestamp for end of nested iteration
    sum_of_niters += (double) (niterend - niterstart); // Accummulate times to calculate average
  }
  iterend = clock(); // Timestamp  for end of all iterations

  /* Copy final result without halo values to buffer */
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
        MPI_Irecv(&master_image[w_index*pwidth][h_index*pheight], 1, array_block, current_process++, 0, comm2d, &request_array[current_process-1]);
      }
    }
    MPI_Waitall(world_size-1, request_array, status_array); // Has to receive all sub-arrays before continuing
    // Copy own sub-array to big buffer
    for (i = 0; i < pwidth; ++i){
      for (j = 0; j < pheight; ++j){
        master_image[i][j] = image[i][j];
      }
    }
    /* Render output */
    char *outputfile;
    // outputfile = "result.pgm";
    outputfile = argv[2];
    pgmwrite(outputfile, &master_image[0][0], width, height);
    end = clock(); // Timestamp for end

    /* Print to console */

    double total_time, beginio, total_iters, mean_niters, write_out;
    beginio = ((double) (iotstamp - start)) / CLOCKS_PER_SEC;
    total_iters = ((double) (iterend - iterstart)) / CLOCKS_PER_SEC;
    mean_niters = ((double) (sum_of_niters/iter)) / CLOCKS_PER_SEC;
    write_out = ((double) (end - iterend)) / CLOCKS_PER_SEC;
    total_time = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Number of Processes = %d\nInput I/O = %f, Total iterations runtime = %f, Single iteration mean runtime = %f, Output I/O = %f, Total runtime = %f\n",\
    world_size, beginio, total_iters, mean_niters, write_out, total_time);
    printf("\n");
  }
  else {
    // Send sub-array to master process
    MPI_Issend(&image[0][0], pwidth*pheight, MPI_DOUBLE, 0, 0, comm2d, &request);
  }
  MPI_Finalize();
  return 0;
}
/* Void function that checks the number of arguments given when executing runnable.
   If the number is incorrect it gives an explanatory error message and exits with status 1.
   
   Arguments
   ---------
   int argc : the number of arguments
*/
void check_num_of_args(int argc){
    if(argc != 3){
      fprintf(stderr, "Incorrect number of arguments.\nExpected 2 arguments but received %d. \
      \nShould receive one argument which is the path leading to the file and one to the output file. \nExample: \
      mpirun -np 4 ./image_exec img/edgenew192x128.pgm output.pgm\n", argc-1);

      exit(1);
    }
}

/* Double function that calculates the absolute difference between two values.
   
   Arguments
   ---------
   double new : the new pixel value
   double old : the old pixel value

   Returns
   -------
   double difference : the absolute value of the difference
*/
double calculate_delta(double new, double old){
  double difference = 0;
  difference = abs(new -old);  

  return difference;
}

/* Double function that calculates the boundary value for a pixel.
   
   Arguments
   ---------
   int i : the i index of the array
   int j : the j index of the array

   Returns
   -------
   double val : the corresponding value of the pixel
*/
double boundaryval(int i, int m){
  double val;

  val = 2.0*((double)(i-1))/((double)(m-1));
  if (i >= m/2+1) val = 2.0-val;

  return val;
}
