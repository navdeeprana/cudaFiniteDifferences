FLAGS     =
ACCFLAGS = -fast -ta=tesla,host -Minfo=accel -cpp $(FLAGS)

all : openacc
openacc :
	pgfortran $(ACCFLAGS) fd3d_openacc.f90 -o exe_openacc

mpi:
	mpif90 $(ACCFLAGS) openacc_mpi_simple.f90 -o exe_pgi_simp



.phony : clean cleandata cleanall

clean:
	@rm -f exe* a.out

cleanobj:
	@rm -rf obj ../src/*.mod *.mod *.o

cleandata:
	@rm -rf data obj

cleanall: clean cleanobj cleandata

