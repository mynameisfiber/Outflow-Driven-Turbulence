CC=mpif77.lam
RUN=mpirun.lam

isoth3d-mpi:
	${CC} -O3 -o isoth3d-mpi isoth3d-mpi.f90 analysis.f90

isoth3d:
	gfortran -O3 -g -pg -fopenmp -o isoth3d isoth3d.f90 

runmpi: isoth3d-mpi
	${RUN} -np `egrep -o "procs([ ]?)=([ 0-9]*)" isoth3d-mpi.f90 | tr -d " " | cut -d"=" -f2` isoth3d-mpi

clean:
	rm -rf output-* isoth3d isoth3d2 gmon.out isoth3d-mpi

