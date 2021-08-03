FLAGS     =
ACCFLAGS = -Mpreprocess -Minform=warn -Minfo=accel -Mcuda=cc70,cuda11.0 -fast -ta=tesla,host $(FLAGS)

all : openacc
openacc :
	nvfortran $(ACCFLAGS) fd3d_openacc.f90 -o exe_openacc

mpi:
	mpif90 $(ACCFLAGS) openacc_mpi_simple.f90 -o exe_pgi_simp


fd_mpi:
	mpif90 $(ACCFLAGS) fd3d_openacc_mpi.f90 -o exe_openacc_mpi


.phony : clean cleandata cleanall

clean:
	@rm -f exe* a.out

cleanobj:
	@rm -rf obj ../src/*.mod *.mod *.o

cleandata:
	@rm -rf data obj

cleanall: clean cleanobj cleandata

