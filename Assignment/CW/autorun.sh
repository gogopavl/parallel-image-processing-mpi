#!/bin/bash
ulimit -s unlimited

procs=4
filename=img/edgenew512x384.pgm
output=output.pgm
logfile=logfile.txt

rm $logfile
touch $logfile



mpirun -np $procs ./image_exec $filename $output  >> $logfile


