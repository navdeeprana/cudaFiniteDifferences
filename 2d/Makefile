FLAGS = 
CUDAFLAGS = -Mcuda=cc70,cuda10.1 -Mpreprocess -O3 $(FLAGS) 

all : mod padded shared sharedx_registery

mod : fd2d_mod.f90
	pgf90 $(CUDAFLAGS) finiteDifference2d.f90 $^ -o exe_$@
padded: fd2d_padded.f90
	pgf90 $(CUDAFLAGS) finiteDifference2d.f90 $^ -o exe_$@
shared : fd2d_shared.f90
	pgf90 $(CUDAFLAGS) finiteDifference2d.f90 $^ -o exe_$@
sharedx_registery : fd2d_sharedx_registery.f90
	pgf90 $(CUDAFLAGS) finiteDifference2d.f90 $^ -o exe_$@

.phony : clean cleanall

clean:
	@rm -f exe* a.out *.mod *.o

cleanall: clean
