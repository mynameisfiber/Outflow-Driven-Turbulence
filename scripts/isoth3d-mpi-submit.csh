#!/bin/csh
#PBS -l nodes=8:ppn=8
#PBS -q workq
#PBS -r n
#PBS -l walltime=48:00:00
#PBS -N isoth3d_mpi_output

cd $PBS_O_WORKDIR
echo "compiling"
make clean && make isoth3d-mpi
echo "starting lam"
lamboot
lamnodes
echo "setting omp threads"
set OMP_NUM_THREADS=16
echo "running"
(mpirun n0,1,2,3,4,5,6,7 ./isoth3d-mpi >& log) && echo "Success" || echo "Fail"
echo "done"
lamhalt
