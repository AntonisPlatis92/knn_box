# knn_box

An implementation of MPI (Message Parsing Interface) in C.

It implements the knn algorithm used to find the nearest point in a 3-D Box from a set of Points.

It divides the problem in a power of 2 (e.g. 8,16,128) sub-problems (smaller boxes) making it possible to ran them parallely in different machines.

Tha machines are communicating when needed through MPI Messages.

The code was written in C (.c) using the mpi.h header file.
