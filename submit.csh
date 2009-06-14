#!/bin/csh
#PBS -l nodes=8:ppn=8
#PBS -q workq
#PBS -r n
#PBS -l walltime=48:00:00
#PBS -N runA

cd $PBS_O_WORKDIR

echo "starting lam"
lamboot

echo "creating directories"
lamnodes | cut -f 2 | cut -d "." -f 1 > nodes
foreach node (`lamnodes | cut -f 2 | cut -d "." -f 1`)
  echo "  $node"
  ssh $node rm -rf /mnt/node_scratch/mgorelick/*
  ssh $node mkdir -p /mnt/node_scratch/mgorelick/output/
end

echo "setting omp threads and compiling"
set OMP_NUM_THREADS=10
make clean && make

echo "Running"
mpirun n0,1,2,3,4,5,6,7 isoth3d-mpi && echo "Success" || echo "Fail"

echo "Recovering output"
mkdir "$PBS_O_WORKDIR/output/"
foreach node (`lamnodes | cut -f 2 | cut -d "." -f 1`)
  echo "  $node"
  scp -r mgorelick@${node}:/mnt/node_scratch/mgorelick/output/ ${PBS_O_WORKDIR}/output/ >& ${PBS_O_WORKDIR}/${node}_transfer.log
  ssh $node "find /mnt/node_scratch/mgorelick/output/ -type f -exec rm -rf {} ;\"
end

echo "Done"
lamhalt
