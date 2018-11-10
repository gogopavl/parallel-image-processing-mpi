
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pgmio.h"

#define WIDTH 512
#define HEIGHT 384
#define MAX_ITERS 1000

int main(int argc, char const *argv[])
{
  printf("bkpt\n");
  double image[WIDTH][HEIGHT]; // Array containing the original pgm file
  double edge[WIDTH + 2][HEIGHT + 2];
  double old[WIDTH + 2][HEIGHT + 2];
  double new[WIDTH + 2][HEIGHT + 2];

  // Declaring input file
  char *filename;
  filename = "edge512x384.pgm";

  // Import file to image array
  printf("\nReading <%s>\n", filename);
  pgmread(filename, image, WIDTH, HEIGHT);
  printf("\n");

  int i, j = 0;


  // Copying image to edge array
  for(i = 1 ; i < WIDTH + 1; ++i){
    for(j = 1 ; j < HEIGHT + 1; ++j) {
      edge[i][j] = image[i-1][j-1];
    }
  }
  
  // Initializing old with white pixels
  for(i = 0 ; i < WIDTH + 2; ++i){
    for(j = 0; j < HEIGHT + 2; ++j){
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