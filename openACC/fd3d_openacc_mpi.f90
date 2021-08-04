program fd3d
    use nvtx    
    use openacc
    use mpi
    use derivative_mod_gpu

    implicit none
    integer, parameter  :: dp = kind(1.d0)
    real(dp), parameter :: two_pi = 8.d0 * atan(1.d0)

    integer  :: n(3), num_iters, block_size(3), itime, local_n3
    real(dp) :: length, dl, dl_inv
    real(dp) :: factors(5)
    real(dp) :: timer(2)
    real(dp), dimension(:, :, :), allocatable :: u, du, du_exact
    integer :: ierr, mpi_rank, mpi_procs
 
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world, mpi_rank,ierr)
    call mpi_comm_size(mpi_comm_world, mpi_procs, ierr)


    call box_init('der3d.nml')
    call allocate_arrays()
    call initial_condition()
    call set_factors()

   ! #ifdef OPENACC
    call acc_set_device_num(mpi_rank,acc_device_nvidia)
   ! #endif

    !$acc data copy(u,du,n,factors)

    write(*,*) "Loop begins"
    call cpu_time(timer(1))
    do itime = 0, num_iters
        call nvtxStartRange("MPI")
        call synchronize_gpu_state()
        call nvtxEndRange
        call nvtxStartRange("Derivative")
!!        call derivative()
        call derivative_halo_gpu_dev(u,du,n(1),n(2),local_n3,factors)
        call nvtxEndRange
        call mpi_barrier(mpi_comm_world,ierr)
    end do
    call cpu_time(timer(2))
    write(*,*) "Loop ends"
    !$acc end data

    write (*, *) 'time taken :', timer(2) - timer(1), 'error :', maxval(du_exact - du(1:n(1), 1:n(2), 1:local_n3))

contains

    subroutine box_init(filename)
        character(*), intent(in) :: filename
        namelist /box_nml/ n, length, num_iters

        open (unit=10, file=filename, status='old')
        read (10, nml=box_nml)
        close (10)

        local_n3 = n(3)/mpi_procs; 
        length = length * two_pi
        dl = length / n(1)
        dl_inv = 1.d0 / dl
    end subroutine box_init

    subroutine set_factors()
        factors = dl_inv * (1.d0 / 12.d0)*[1, -8, 0, 8, -1]
    end subroutine set_factors

    subroutine allocate_arrays()
        allocate (u(-1:n(1) + 2, -1:n(2) + 2, -1:local_n3 + 2))
        allocate (du(n(1), n(2), local_n3))
        allocate (du_exact(n(1), n(2), local_n3))

        u = 0.d0
        du = 0.d0
        du_exact = 0.d0
    end subroutine allocate_arrays

    subroutine initial_condition()
        integer  :: i, j, k
        real(dp) :: x, y, z

        u(:, :, :) = 0.d0
        du_exact(:, :, :) = 0.d0

        do k = 1, local_n3
            z = dl * (k - 1 + local_n3*mpi_rank)
            do j = 1, n(2)
                y = dl * (j - 1)
                do i = 1, n(1)
                    x = dl * (i - 1)
                    u(i, j, k) = sin(x) + sin(y) + sin(z)
                end do
            end do
        end do

        do k = 1, local_n3
            z = dl * (k - 1 + local_n3*mpi_rank)
            do j = 1, n(2)
                y = dl * (j - 1)
                do i = 1, n(1)
                    x = dl * (i - 1)
                    du_exact(i, j, k) = cos(x) + cos(y) + cos(z)
                end do
            end do
        end do
    end subroutine initial_condition

    subroutine derivative()
        integer  :: i, j, k
        !$acc parallel loop collapse(3) present(u,du,n,factors)
        do k = 1, local_n3
            do j = 1, n(2)
                do i = 1, n(1)
                    du(i, j, k) = (u(i - 2, j, k) + u(i, j - 2, k) + u(i, j, k - 2)) * factors(1) + &
                        (u(i - 1, j, k) + u(i, j - 1, k) + u(i, j, k - 1)) * factors(2) + &
                        (u(i, j, k) + u(i, j, k) + u(i, j, k)) * factors(3) + &
                        (u(i + 1, j, k) + u(i, j + 1, k) + u(i, j, k + 1)) * factors(4) + &
                        (u(i + 2, j, k) + u(i, j + 2, k) + u(i, j, k + 2)) * factors(5)
                end do
            end do
        end do
        !$acc end parallel
    end subroutine derivative

    subroutine synchronize_gpu_state()
        integer :: i,j,k
        integer :: data_size, right, left, mpi_error_flag
        integer :: mpi_status_flag(mpi_status_size)

        !! Periodic Boundary Setup
        !! X-direction
        
        !$acc parallel loop collapse(3)  present(u)  
        do k=1,local_n3
        do j=1,n(2)
!!        !$acc loop independent
        do i=-1,0
        u(i,j,k)      = u(n(1)+i,j,k)
        u(n(1)+i+2,j,k) = u(i+2,j,k)
        enddo
        enddo
        enddo
        !$acc end parallel 

        !! Y-direction
        !$acc parallel loop collapse(3) present(u)
        do k=1,local_n3
        do j=-1,0 
!!        !$acc loop independent
        do i=1,n(1)
        u(i,j,k)      = u(i,n(2)+j,k)
        u(i,n(2)+j+2,k) = u(i,j+2,k)
         enddo
         enddo
         enddo
        !$acc end parallel 

        !! Z-Direction Synchronisation, If only 1 process, no point using mpi.

        data_size = 2*((n(1)+4)*(n(2)+4))

        !! if (mpi_procs ==1) then
        !!     u(:,:,0)          = u(:,:,local_n3)
        !!     u(:,:,-1)         = u(:,:,local_n3-1)

        !!     u(:,:,local_n3+1) = u(:,:,1)
        !!     u(:,:,local_n3+2) = u(:,:,2)

        !! else

        if (mpi_rank == mpi_procs - 1) then
              right = 0
        else
             right = mpi_rank+1
        end if
        if (mpi_rank == 0) then
              left = mpi_procs - 1
        else
              left = mpi_rank - 1
        end if

        !$acc host_data use_device( u ) 

        !! Send right boundary elements to next processor, and recieve the same from previous.
        call MPI_SENDRECV(u(-1,-1,local_n3-1),data_size,MPI_DOUBLE_PRECISION,right,0,&
            u(-1,-1,-1),data_size,MPI_DOUBLE_PRECISION,left,0,MPI_COMM_WORLD,mpi_status_flag,mpi_error_flag)

        !! Send left boundary elements to previous processor, which stores it at left.
        call MPI_SENDRECV(u(-1,-1,1),data_size,MPI_DOUBLE_PRECISION,left,1,&
            u(-1,-1,local_n3+1),data_size,MPI_DOUBLE_PRECISION,right,1,MPI_COMM_WORLD,mpi_status_flag,mpi_error_flag)
        !! endif
        !$acc end host_data
    end subroutine synchronize_gpu_state


end program fd3d
