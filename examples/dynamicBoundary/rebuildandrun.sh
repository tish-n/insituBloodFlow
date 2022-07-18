#!/bin/sh

cd tmp

rm -rf * 

cd ..

cd build

make -j8 

cd ..

mpirun -np 4 ./dynamicBoundary in.embolism 1000 10