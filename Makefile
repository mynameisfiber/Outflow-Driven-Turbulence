CC=mpif77
RUN=mpirun

isoth3d-mpi:
	#-fdefault-real-8  -fopenmp
	${CC} -Wall -fopenmp -ffast-math -O3 -lm -o isoth3d-mpi params.f90 analysis.f90 isoth3d-mpi.f90 

isoth3d:
	gfortran -O3 -g -pg -fopenmp -o isoth3d isoth3d.f90 

runmpi: isoth3d-mpi
	${RUN} -ger -np `egrep -o "procs([ ]?)=([ 0-9]*)" params.f90 | tr -d " " | cut -d"=" -f2` isoth3d-mpi

runmpitv: isoth3d-mpi
	${RUN} -tv isoth3d-mpi

clean:
	rm -rf nodelog-* output-* isoth3d isoth3d2 gmon.out isoth3d-mpi analysis.mod params.mod
