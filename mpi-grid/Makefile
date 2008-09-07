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
	rm -rf gridcontrol output-*-*

run: gridcontrol
	${RUN} -np `egrep -o "procs([ ]?)=([ 0-9]*)" gridcontrol.f90 | tr -d " " | cut -d"=" -f2` gridcontrol
