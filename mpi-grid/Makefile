gridcontrol:
	mpif90 -o gridcontrol gridcontrol.f90

clean: 
	rm -rf gridcontrol output-test

run: gridcontrol
	mpirun -np 8 gridcontrol
