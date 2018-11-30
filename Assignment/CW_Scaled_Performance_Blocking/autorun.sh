#!/bin/bash
ulimit -s unlimited
rm out.txt
touch out.txt

filename=img/edgenew768x768.pgm

for j in `seq 1 4`;
do
	mpirun -np $j ./image_exec $filename  >> out.txt
done 

