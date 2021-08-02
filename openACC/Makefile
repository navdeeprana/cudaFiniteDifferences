FLAGS     =
ACCFLAGS = -fast -ta=tesla,host -Minfo=accel -cpp $(FLAGS)

all : openacc
openacc :
	pgfortran $(ACCFLAGS) fd3d_openacc.f90 -o exe_openacc

.phony : clean cleandata cleanall

clean:
	@rm -f exe* a.out

cleanobj:
	@rm -rf obj ../src/*.mod *.mod *.o

cleandata:
	@rm -rf data obj

cleanall: clean cleanobj cleandata

