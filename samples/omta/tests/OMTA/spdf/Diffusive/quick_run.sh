#!/bin/bash

for d1 in L10
do

cd $d1

for d0 in 4x4 #5x5 6x6
do

cd $d0

for d2 in cf-1
do

cd $d2

ulimit -s unlimited
date "+Began at: %H:%M:%S"

~/Work/trunk/bin/htrans.x

#mpirun -np 32 ~/work_trunk/bin/htrans.x 
date "+Finished at: %H:%M:%S"


cd ../
done
cd ../
done
cd ../
done
