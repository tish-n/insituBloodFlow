#!/bin/sh

cd tmp

rm -rf * 

cd ..

cd build

make -j20 

cd ..

mpirun -np 8 ./dynamicBoundary in.lmp4cell 1000 10
