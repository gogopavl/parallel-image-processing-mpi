
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pgmio.h"

#define M 192
#define N 128

int main(int argc, char const *argv[])
{
  double image[M][N];

  char *filename;  
  filename = "edge192x128.pgm";
  printf("\nReading <%s>\n", filename);
  pgmread(filename, image, M, N);
  printf("\n");

  char *outputfile;
  outputfile = "out.pgm";
  printf("\nWriting %s\n", outputfile);
  pgmwrite(outputfile, image, M, N);
  
  return 0;
}
