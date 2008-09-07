ifeq ("`hostname`","mercury")
	CC=mpif77.lam
	RUN=mpirun.lam
else
	CC=mpif90
	RUN=mpirun
endif

gridcontrol:
	${CC} -o gridcontrol gridcontrol.f90

clean: 
	rm -rf gridcontrol output*-*

run: gridcontrol
	${RUN} -np 8 gridcontrol
