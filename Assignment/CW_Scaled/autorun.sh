#!/bin/bash
ulimit -s unlimited
rm out.txt
touch out.txt

filename=img/edgenew768x768.pgm

for j in `seq 1 4`;
do
	echo "---------------------------------------" >> out.txt
	echo "Number of processes = $j" >> out.txt
	echo "---------------------------------------" >> out.txt
	mpirun -np $j ./image_exec $filename  >> out.txt
	echo "---------------------------------------" >> out.txt
done 

