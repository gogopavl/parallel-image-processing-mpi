#!/bin/bash
rm out.txt
touch out.txt
for j in `seq 1 4`;
do
	mpirun -np $j ./image_exec img/edgenew192x128.pgm  >> out.txt
	echo "---------------------------------------" >> out.txt
	echo "---------------------------------------" >> out.txt
	echo "---------------------------------------" >> out.txt
done 

