#!/bin/csh
#PBS -l nodes=8:ppn=8
#PBS -q workq
#PBS -r n
#PBS -l walltime=48:00:00
#PBS -N col0.05-400dx

cd $PBS_O_WORKDIR

echo "starting lam"
lamboot

echo "creating directories"
lamnodes | cut -f 2 | cut -d "." -f 1 > nodes
foreach node (`lamnodes | cut -f 2 | cut -d "." -f 1`)
  echo "  $node"
  mkdir -p /mnt/scratch/$node/mgorelick/
  rm -rf /mnt/scratch/$node/mgorelick/*
end

echo "setting omp threads and compiling"
set OMP_NUM_THREADS=8
make clean && make

echo "Running"
mpirun n0,1,2,3,4,5,6,7 isoth3d-mpi && echo "Success" || echo "Fail"

echo "Recovering output"
mkdir "$PBS_O_WORKDIR/output/"
foreach node (`lamnodes | cut -f 2 | cut -d "." -f 1`)
  echo "  $node"
  mv /mnt/scratch/$node/mgorelick/* $PBS_O_WORKDIR/output/
end

echo "Done"
lamhalt
