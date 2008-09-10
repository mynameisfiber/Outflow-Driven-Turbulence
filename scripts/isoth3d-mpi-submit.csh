#!/bin/csh
#PBS -l nodes=4:ppn=8
#PBS -q workq
#PBS -r n
#PBS -l walltime=06:00:00
#PBS -N isoth3d_cts
# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE nodes,ppn,walltime and my_job_name VALUES

cd $PBS_O_WORKDIR
make clean && make isoth3d-mpi
lamboot
set OMP_NUM_THREADS=8
mpirun ./isoth3d-mpi
lamhalt
