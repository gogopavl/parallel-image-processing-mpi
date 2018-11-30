#!/bin/bash
ulimit -s unlimited

max_procs=4
filename=img/edgenew512x384.pgm
output=output.pgm
logfile=logfile_loop.txt

rm $logfile
touch $logfile



for j in `seq 1 $max_procs`;
do
	mpirun -np $j ./image_exec $filename $output  >> $logfile
done 

