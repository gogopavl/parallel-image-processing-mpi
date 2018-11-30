#!/bin/bash
rm out.txt
touch out.txt
for j in `seq 1 10`;
do
	mpirun -np 4 ./image_exec >> out.txt
	echo "---------------------------------------" >> out.txt
	echo "---------------------------------------" >> out.txt
	echo "---------------------------------------" >> out.txt
done 

