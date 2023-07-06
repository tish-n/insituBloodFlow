#!/bin/sh

mpirun -n 8 /home/ntishchenko/test/insituBloodFlow/src/lammps/src/lmp_mpi -in in.lmp4cell #> out.lmp4cell