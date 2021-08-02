module mod_mpi
  implicit none
  save
  
  integer:: ierr,myrank,nprocs
  include 'mpif.h'
  
contains
  
  subroutine mpi_initialize()
    
    call mpi_init(ierr)
    !! Get the rank
    call mpi_comm_rank(mpi_comm_world, myrank,ierr)
    !! Get the number of processes
    call mpi_comm_size(mpi_comm_world, nprocs, ierr)
    
  end subroutine mpi_initialize

  
end module mod_mpi



module mod_param_array
  implicit none
  save
  real*8,allocatable,dimension(:)::A
  integer::N
  character(100)::fname
  
  namelist /box_nml/ N
  
end module mod_param_array




program test_mpi_openacc
  use openacc
  use mod_mpi
  use mod_param_array
  implicit none
  integer::i, j
  
  
  call mpi_initialize
! #ifdef OPENACC
  call acc_set_device_num(myrank,acc_device_nvidia)
! #endif

  
  open(unit=10,file="input.dat",status="old")
  read(10,nml=box_nml)
  close(10)
  
  allocate(A(N/nprocs))

  A=0.d0
 !$acc data copy(A,N,nprocs,myrank)
  do j = 1, N/nprocs
  !$acc kernels present(A,myrank,nprocs,N)
  do i=1,N/nprocs
     A(i)=i+myrank*N/nprocs
  enddo
  !$acc end kernels
  end do

  !$acc end data


  write(fname,'(g8.0)') myrank+1
  open(unit=10,file="out"//trim(adjustl(fname)),status="unknown")
  do i=1,N/nprocs
     write(10,*) i,A(i)
  enddo
  close(10)
  
  
  call mpi_finalize(ierr)
  
end program test_mpi_openacc
