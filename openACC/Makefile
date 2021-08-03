FLAGS     =
ACCFLAGS = -Mpreprocess -Minform=warn -Minfo=accel -Mcuda=cc70,cuda11.0 -fast -ta=tesla,host $(FLAGS)


NSYS=/opt/nvidia/hpc_sdk/Linux_x86_64/20.9/compilers/bin/nsys profile
NSYSFLAGS=-t nvtx,cuda --force-overwrite=true

all : openacc

nvtx.o:	 nvtx.f
		nvfortran -c -Mfree nvtx.f

openacc :
	nvfortran $(ACCFLAGS) fd3d_openacc.f90 -o exe_openacc

mpi:
	mpif90 $(ACCFLAGS) openacc_mpi_simple.f90 -o exe_pgi_simp


fd_mpi: nvtx.o
	mpif90 $(ACCFLAGS) nvtx.o fd3d_openacc_mpi.f90 -o exe_openacc_mpi -lnvToolsExt

profile: exe_openacc_mpi
	$(NSYS) $(NSYSFLAGS) -o openacc_mpi_report mpirun -np 1  ./exe_openacc_mpi output_profile 

.phony : clean cleandata cleanall

clean:
	@rm -f exe* a.out

cleanobj:
	@rm -rf obj ../src/*.mod *.mod *.o

cleandata:
	@rm -rf data obj

cleanall: clean cleanobj cleandata

