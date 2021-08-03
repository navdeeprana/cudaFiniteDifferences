program fd3d
    implicit none
    integer, parameter  :: dp = kind(1.d0)
    real(dp), parameter :: two_pi = 8.d0 * atan(1.d0)

    integer  :: n(3), num_iters, block_size(3), itime
    real(dp) :: length, dl, dl_inv
    real(dp) :: factors(5)
    real(dp) :: timer(2)
    real(dp), dimension(:, :, :), allocatable :: u, du, du_exact

    call box_init('der3d.nml')
    call allocate_arrays()
    call initial_condition()
    call set_factors()

    !$acc data copy(u,du,n,factors)
    write(*,*) "Profiling start"
    call cpu_time(timer(1))
    do itime = 0, num_iters
        call derivative_halo()
        ! call derivative_mod()
        ! call derivative_mod_ver2()
    end do
    call cpu_time(timer(2))
    !$acc end data
    write (*, *) 'time taken :', timer(2) - timer(1), 'error :', maxval(du_exact - du)
    write(*,*) "Profiling end"

contains

    subroutine box_init(filename)
        character(*), intent(in) :: filename
        namelist /box_nml/ n, length, num_iters

        open (unit=10, file=filename, status='old')
        read (10, nml=box_nml)
        close (10)

        length = length * two_pi
        dl = length / n(1)
        dl_inv = 1.d0 / dl
    end subroutine box_init

    subroutine set_factors()
        factors = dl_inv * (1.d0 / 12.d0)*[1, -8, 0, 8, -1]
    end subroutine set_factors

    subroutine allocate_arrays()
        allocate (u(-1:n(1) + 2, -1:n(2) + 2, -1:n(3) + 2))
        allocate (du(n(1), n(2), n(3)))
        allocate (du_exact(n(1), n(2), n(3)))

        u = 0.d0
        du = 0.d0
        du_exact = 0.d0
    end subroutine allocate_arrays

    subroutine initial_condition()
        integer  :: i, j, k
        real(dp) :: x, y, z

        u(:, :, :) = 0.d0
        du_exact(:, :, :) = 0.d0

        do k = 1, n(3)
            z = dl * (k - 1)
            do j = 1, n(2)
                y = dl * (j - 1)
                do i = 1, n(1)
                    x = dl * (i - 1)
                    u(i, j, k)        = sin(x) + sin(y) + sin(z)
                    du_exact(i, j, k) = cos(x) + cos(y) + cos(z)
                end do
            end do
        end do

        ! Periodic boundary conditions
        u(0,:,:)      = u(n(1),:,:)
        u(-1,:,:)     = u(n(1)-1,:,:)
        u(n(1)+1,:,:) = u(1,:,:)
        u(n(1)+2,:,:) = u(2,:,:)

        u(:,0,:)      = u(:,n(2),:)
        u(:,-1,:)     = u(:,n(2)-1,:)
        u(:,n(2)+1,:) = u(:,1,:)
        u(:,n(2)+2,:) = u(:,2,:)

        u(:,:,0)      = u(:,:,n(3))
        u(:,:,-1)     = u(:,:,n(3)-1)
        u(:,:,n(3)+1) = u(:,:,1)
        u(:,:,n(3)+2) = u(:,:,2)
    end subroutine initial_condition

    subroutine derivative_halo()
        integer  :: i, j, k

        !$acc parallel loop collapse(3) present(u,du,n,factors)
        do k = 1, n(3)
            do j = 1, n(2)
                do i = 1, n(1)
                    du(i, j, k) = (u(i - 2, j, k) + u(i, j - 2, k) + u(i, j, k - 2)) * factors(1) + &
                                  (u(i - 1, j, k) + u(i, j - 1, k) + u(i, j, k - 1)) * factors(2) + &
                                  (u(i, j, k)     + u(i, j, k)     + u(i, j, k)    ) * factors(3) + &
                                  (u(i + 1, j, k) + u(i, j + 1, k) + u(i, j, k + 1)) * factors(4) + &
                                  (u(i + 2, j, k) + u(i, j + 2, k) + u(i, j, k + 2)) * factors(5)
                end do
            end do
        end do
        !$acc end parallel
    end subroutine derivative_halo

    subroutine derivative_mod()
        integer :: i, j, k
        integer :: l, xn(5), yn(5), zn(5)

        !$acc parallel loop collapse(3) present(u,du,n,factors) private(xn,yn,zn)
        do k = 1, n(3)
            do j = 1, n(2)
                do i = 1, n(1)
                    do l = -2, 2
                        xn(l + 3) = mod(i + l + n(1) - 1, n(1)) + 1
                        yn(l + 3) = mod(j + l + n(2) - 1, n(2)) + 1
                        zn(l + 3) = mod(k + l + n(3) - 1, n(3)) + 1
                    end do
                    du(i, j, k) = (u(xn(1), j, k) + u(i, yn(1), k) + u(i, j, zn(1))) * factors(1) + &
                                  (u(xn(2), j, k) + u(i, yn(2), k) + u(i, j, zn(2))) * factors(2) + &
                                  (u(xn(3), j, k) + u(i, yn(3), k) + u(i, j, zn(3))) * factors(3) + &
                                  (u(xn(4), j, k) + u(i, yn(4), k) + u(i, j, zn(4))) * factors(4) + &
                                  (u(xn(5), j, k) + u(i, yn(5), k) + u(i, j, zn(5))) * factors(5)
                end do
            end do
        end do
        !$acc end parallel
    end subroutine derivative_mod

    subroutine derivative_mod_ver2()
        integer :: i, j, k
        integer :: l, xn, yn, zn
        real(dp) :: duijk

        !$acc parallel loop collapse(3) present(u,du,n,factors)
        do k = 1, n(3)
            do j = 1, n(2)
                do i = 1, n(1)
                    duijk = 0.0_dp
                    do l = -2, 2
                        xn    = mod(i + l + n(1) - 1, n(1)) + 1
                        yn    = mod(j + l + n(2) - 1, n(2)) + 1
                        zn    = mod(k + l + n(3) - 1, n(3)) + 1
                        duijk = duijk + (u(xn, j, k) + u(i, yn, k) + u(i, j, zn)) * factors(l+3)
                    end do
                    du(i,j,k) =  duijk
                end do
            end do
        end do
        !$acc end parallel
    end subroutine derivative_mod_ver2
end program fd3d
