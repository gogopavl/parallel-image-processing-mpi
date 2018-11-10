
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pgmio.h"

#define WIDTH 192
#define HEIGHT 128

int main(int argc, char const *argv[])
{
  double image[WIDTH][HEIGHT]; // Array containing the original pgm file

  // Declaring input file
  char *filename;  
  filename = "edge192x128.pgm";

  // Import file to image array
  printf("\nReading <%s>\n", filename);
  pgmread(filename, image, WIDTH, HEIGHT);
  printf("\n");
  
  // Export image to output file
  char *outputfile;
  outputfile = "out.pgm";
  printf("\nWriting %s\n", outputfile);
  pgmwrite(outputfile, image, WIDTH, HEIGHT);
  
  return 0;
}
