FLAGS     =
CUDAFLAGS = -Mcuda=cc70,cuda11.0 -Mpreprocess -O3 $(FLAGS)
all : cpu mod halo sharedxy
cpu :
	gfortran $(FLAGS) -O3 fd3d_cpu.f90 -o exe_cpu
mod :
	nvfortran $(CUDAFLAGS) finiteDifference3d.f90 fd3d_mod.f90 -o exe_mod
halo :
	nvfortran $(CUDAFLAGS) finiteDifference3d.f90 fd3d_halo.f90 -o exe_halo
sharedxy :
	nvfortran $(CUDAFLAGS) finiteDifference3d.f90 fd3d_sharedxy.f90 -o exe_sharedxy

.phony : clean cleandata cleanall

clean:
	@rm -f exe* a.out

cleanobj:
	@rm -rf obj ../src/*.mod *.mod *.o

cleandata:
	@rm -rf data obj

cleanall: clean cleanobj cleandata
