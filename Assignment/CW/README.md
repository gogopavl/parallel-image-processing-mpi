# Message Passing Programming Coursework

## Folder contents

- img/ : Folder containing the provided image files
- main.c : C Program with the implementation
- pgmio.c : C Program with I/O functions
- pgmio.h : Header file of pgmio.c
- Makefile : Project Makefile
- image_exec : Program executable
- logfile.txt : Console output of execution
- logfile_loop.txt : Console output of loop execution
- output.pgm : Exemplary output of the prgram
- subimage.pbs : Exemplary Cirrus job submission file
- autorun.sh : Bash script that contains all necessary commands for execution
- autorun_loop.sh : Bash script that contains all necessary commands for loop execution


## Running the program

### First option

From a shell run the autorun.sh script (".\autorun.sh")

### Second option

From a shell run the autorun_loop.sh script (".\autorun_loop.sh")

### Third option

From a shell run:

- make
- mpirun -np num_of_procs ./image_exec img/input_file output_file

Where:

- num_of_procs: the desired number of processes (e.g. 16)
- input_file: path leading to the input file (e.g. img/edgenew512x384.pgm)
- output_file: output file name (e.g. output.pgm) - Caution! Has to be filename, not path and extension should be *.pgm

### Notes

- The autorun.sh script runs the program for all a specified number of processes. This can be changed within the file.

- The autorun_loop.sh script runs the program for all number of processes within the specified iteration loop. This can also be changed within the file.


