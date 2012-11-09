LU Factorization - MPI
=============================
The project is part of the Assignment 2 CCS 524 - Parallel Computing Architectures & Algorithms

COMPILE			
=======
make 

manually compile:
mpic++ LU-mpi.cpp -o lu-mpi -DMPICH_IGNORE_CXX_SEEK 

RUN THE PROGRAM		
===============
mpirun -np <num-of-process> lu-mpi <size-of-linear-system>

e.g: mpirun -np 16 lu-mpi 900

