FLAGS     =
CUDAFLAGS = -Mcuda=cc70,cuda10.1 -Mpreprocess -O3 $(FLAGS)
all : cpu mod sharedxy
cpu :
	gfortran $(FLAGS) -O3 fd3d_cpu.f90 -o exe_cpu
mod :
	nvfortran $(CUDAFLAGS) finiteDifference3d.f90 fd3d_mod.f90 -o exe_mod
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
