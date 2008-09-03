isoth3d:
	gfortran -O3 -g -pg -fopenmp -o isoth3d isoth3d.f90 

clean:
	rm -rf output-* isoth3d isoth3d2 gmon.out

