program fd3d

    integer,  parameter :: dp     = kind(1.d0)
    real(dp), parameter :: two_pi = 8.d0*atan(1.d0)
    
    integer  :: n(3), num_iters, block_size(3), itime
    real(dp) :: length, dl, dl_inv
    real(dp) :: factors(5)
    real(4)  :: cpuStart, cpuFinish
    real(dp), dimension(:, :, :), allocatable :: u, du, du_exact

    call box_init('der3d.nml')
    call allocate_arrays()
    call initial_condition()
    call set_factors()

    call cpu_time(cpuStart)
    do itime = 0, num_iters
        call divergence_cpu()
    end do
    call cpu_time(cpuFinish)
    write(*,*) 'cpu            :', cpuStart-cpuFinish, maxval(du_exact-du)

contains

    subroutine box_init(filename)
        character(*), intent(in) :: filename
        namelist /box_nml/ n, length, num_iters

        open(unit=10,file=filename,status='old')
        read(10,nml=box_nml)
        close(10)

        length = length*two_pi
        dl     = length/n(1)
        dl_inv = 1.d0/dl
    end subroutine box_init

    subroutine set_factors()
        factors   = dl_inv*(1.d0/12.d0)*[ 1, -8, 0, 8, -1 ]
    end subroutine set_factors

    subroutine allocate_arrays()
        allocate(u (n(1),n(2),n(3)))
        allocate(du(n(1),n(2),n(3)))
        allocate(du_exact(n(1),n(2),n(3)))

        u        = 0.d0
        du       = 0.d0
        du_exact = 0.d0
    end subroutine allocate_arrays

    subroutine initial_condition()
        integer  :: i,j,k
        real(dp) :: x,y,z

        u(:,:,:)        = 0.d0
        du_exact(:,:,:) = 0.d0

        !$OMP PARALLEL DO PRIVATE(i,j,x,y,z)
        do k = 1,n(2)
            z = dl*(k-1)
            do j = 1,n(2)
                y = dl*(j-1)
                do i = 1,n(1)
                    x = dl*(i-1)
                    u(i,j,k)        = sin(x)+sin(y)+sin(z)
                    du_exact(i,j,k) = cos(x)+cos(y)+cos(z)
                end do
            end do
        end do
    end subroutine initial_condition

    subroutine divergence_cpu()
        integer  :: i, j, k, l, xn(5),yn(5),zn(5)
        real(dp) :: ux(5), uy(5), uz(5)

        !$OMP PARALLEL DO PRIVATE(i,j,l,xn,yn,zn,ux,uy,uz)
        do k = 1, n(3)
            do j = 1, n(2)
                do i = 1, n(1)
                    do l = -2, 2
                        xn(l+3) = mod(i+l+n(1)-1,n(1)) + 1
                        yn(l+3) = mod(j+l+n(2)-1,n(2)) + 1
                        zn(l+3) = mod(k+l+n(3)-1,n(3)) + 1
                    end do

                    ux = u(xn,j,k)
                    uy = u(i,yn,k)
                    uz = u(i,j,zn)
                    du(i,j,k) = (ux(1)+uy(1)+uz(1))*factors(1) &
                              + (ux(2)+uy(2)+uz(2))*factors(2) &
                              + (ux(3)+uy(3)+uz(3))*factors(3) &
                              + (ux(4)+uy(4)+uz(4))*factors(4) &
                              + (ux(5)+uy(5)+uz(5))*factors(5)          
                end do
            end do
        end do
    end subroutine divergence_cpu

end program fd3d
