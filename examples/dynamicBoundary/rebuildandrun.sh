#!/bin/sh

cd tmp

rm -rf * 

cd ..

cd build

make -j20 

cd ..

mpirun -np 20 ./dynamicBoundary in.embolism 1000 10
