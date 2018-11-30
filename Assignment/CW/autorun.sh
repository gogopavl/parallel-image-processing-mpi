#!/bin/bash
ulimit -s unlimited

logfile=logfile.txt

rm $logfile
touch $logfile

filename=img/edgenew512x384.pgm
output=output.pgm

for j in `seq 1 4`;
do
	mpirun -np $j ./image_exec $filename $output  >> $logfile
done 

