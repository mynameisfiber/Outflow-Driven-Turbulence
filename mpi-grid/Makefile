gridcontrol:
	mpif90.openmpi -o gridcontrol gridcontrol.f90

clean: 
	rm -rf gridcontrol output-test

run: gridcontrol
	mpirun -np 27 gridcontrol
