#!/bin/csh
#PBS -l nodes=1:ppn=8:rack2
#PBS -q workq
#PBS -r n
#PBS -l walltime=06:00:00
#PBS -N isoth3d_cts
# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE nodes,ppn,walltime and my_job_name VALUES
cd $PBS_O_WORKDIR
/bin/bash ./cts_run.bash