#!/bin/sh

rm dynamicBoundary.8*

cd tmp

rm vtk* 

cd ..

cd build

make -j 32 

cd ..

mpirun -np 8 ./dynamicBoundary in.lmp4cell 1000 10
