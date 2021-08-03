#define XWORK (32)
#define YWORK (4)
#define ZWORK (16)
module finiteDifference3d

    use cudafor

    implicit none

    integer, parameter  :: dp = kind(1.d0)
    real(dp), parameter :: two_pi = 8.d0 * atan(1.d0)
    !! factors

    integer  :: n(3), num_iters, block_size(3), itime
    real(dp) :: length, dl, dl_inv
    !! Input Parameters

    real(dp)  :: factors(5)
    real(dp), constant :: factors_d(5)

    real(dp), dimension(:, :, :), allocatable         :: u, du, du_exact
    real(dp), dimension(:, :, :), device, allocatable :: u_d, du_d

contains

    subroutine get_block_size()
        character(len=128) :: cmdArg
        call get_command_argument(1, cmdArg)
        read (cmdArg, *) block_size(1)
        call get_command_argument(2, cmdArg)
        read (cmdArg, *) block_size(2)
        call get_command_argument(3, cmdArg)
        read (cmdArg, *) block_size(3)
    end subroutine get_block_size

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
        factors_d = factors
    end subroutine set_factors

    subroutine allocate_arrays()
        allocate (u(-1:n(1) + 2, -1:n(2) + 2, -1:n(3) + 2))
        allocate (du(n(1), n(2), n(3)))
        allocate (du_exact(n(1), n(2), n(3)))

        allocate (u_d(-1:n(1) + 2, -1:n(2) + 2, -1:n(3) + 2))
        allocate (du_d(n(1), n(2), n(3)))
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
                    u(i, j, k) = sin(x) + sin(y) + sin(z)
                    du_exact(i, j, k) = cos(x) + cos(y) + cos(z)
                end do
            end do
        end do

        ! Periodic boundary conditions
        u(0, :, :) = u(n(1), :, :)
        u(-1, :, :) = u(n(1) - 1, :, :)
        u(n(1) + 1, :, :) = u(1, :, :)
        u(n(1) + 2, :, :) = u(2, :, :)

        u(:, 0, :) = u(:, n(2), :)
        u(:, -1, :) = u(:, n(2) - 1, :)
        u(:, n(2) + 1, :) = u(:, 1, :)
        u(:, n(2) + 2, :) = u(:, 2, :)

        u(:, :, 0) = u(:, :, n(3))
        u(:, :, -1) = u(:, :, n(3) - 1)
        u(:, :, n(3) + 1) = u(:, :, 1)
        u(:, :, n(3) + 2) = u(:, :, 2)

        !! Copy to device
        u_d(:, :, :) = u(:, :, :)
        !! Set the device output to 0
        du_d(:, :, :) = 0.d0
    end subroutine initial_condition

    attributes(global) subroutine divergence_sharedxy(n1, n2, n3, u, du)
        !! divergence_shared :
        !! Compute the divergence using mod function to access neighbours on boundaries.
        !! The memory access on boundaries in not-contigous.
        integer, value, intent(in)                               :: n1, n2, n3
        real(dp), dimension(-1:n1 + 2, -1:n2, -1:n3), intent(in) :: u
        real(dp), dimension(n1, n2, n3), intent(out)             :: du

        integer          :: i, j, k, iG, jG, kW, kB
        real(dp)         :: front(2), center, behind(2)
        !! Put z-stuff on registers
        real(dp), shared :: u_ds(-1:blockDim%x + 2, -1:blockDim%y + 2)
        !! Shared memory for xy planes

        i = threadIdx%x
        j = threadIdx%y
        k = threadIdx%z
        !! Local thread index

        iG = (blockIdx%x - 1) * blockDim%x + threadIdx%x
        jG = (blockIdx%y - 1) * blockDim%y + threadIdx%y
        !! Global index, will be used to access u_d/write du_d

        kB = ZWORK * ((blockIdx%z - 1) * blockDim%z + threadIdx%z - 1) + 1
        !! Global Index on which threads are going to work

        if (blockIdx%z == 1) then
            behind = [u_d(iG, jG, n3 - 1), u_d(iG, jG, n3)]
        else
            behind = [u_d(iG, jG, kB - 2), u_d(iG, jG, kB - 1)]
        end if
        center = u_d(iG, jG, kB)
        front = [u_d(iG, jG, kB + 1), u_d(iG, jG, kB + 2)]

        do kW = 1, ZWORK
            u_ds(i, j) = center
            if (i <= 2) then
                !! check if block is not on x edge, if yes gotta fetch from edges.
                if (blockIdx%x == 1) then
                    u_ds(i - 2, j) = u_d(n1 + i - 2, jG, kB)        !! Left padding
                    u_ds(i + blockDim%x, j) = u_d(iG + blockDim%x, jG, kB)  !! Right padding
                else if (blockIdx%x == gridDim%x) then
                    u_ds(i - 2, j) = u_d(iG - 2, jG, kB)       !! Left padding
                    u_ds(i + blockDim%x, j) = u_d(i, jG, kB)  !! Right padding
                else
                    u_ds(i - 2, j) = u_d(iG - 2, jG, kB)       !! Left padding
                    u_ds(i + blockDim%x, j) = u_d(iG + blockDim%x, jG, kB)  !! Right padding
                end if
            end if
            if (j <= 2) then
                if (blockIdx%y == 1) then
                    u_ds(i, j - 2) = u_d(iG, n2 + j - 2, kB)                !! Up padding
                    u_ds(i, j + blockDim%y) = u_d(iG, jG + blockDim%y, kB)  !! Down padding
                else if (blockIdx%y == gridDim%y) then
                    u_ds(i, j - 2) = u_d(iG, jG - 2, kB)                  !! Up padding
                    u_ds(i, j + blockDim%y) = u_d(iG, j, kB)            !! Down padding
                else
                    u_ds(i, j - 2) = u_d(iG, jG - 2, kB)                  !! Up padding
                    u_ds(i, j + blockDim%y) = u_d(iG, jG + blockDim%y, kB)  !! Down padding
                end if
            end if
            call syncthreads()
            du_d(iG, jG, kB) = (u_ds(i - 2, j) + u_ds(i, j - 2) + behind(1)) * factors_d(1) &
                               + (u_ds(i - 1, j) + u_ds(i, j - 1) + behind(2)) * factors_d(2) &
                               + (u_ds(i, j) + u_ds(i, j) + center) * factors_d(3) &
                               + (u_ds(i + 1, j) + u_ds(i, j + 1) + front(1)) * factors_d(4) &
                               + (u_ds(i + 2, j) + u_ds(i, j + 2) + front(2)) * factors_d(5)

            behind(1) = behind(2)
            behind(2) = center
            center = front(1)
            front(1) = front(2)

            kB = kB + 1
            if (blockIdx%z == gridDim%z) then
                front(2) = u_d(iG, jG, mod(kB + n3 + 1, n3) + 1)
            else
                front(2) = u_d(iG, jG, kB + 2)
            end if
            call syncthreads()
        end do
    end subroutine divergence_sharedxy

    attributes(global) subroutine divergence_mod(n1, n2, n3, u, du)
        !! divergence_mod :
        !! Compute the divergence using mod function to access neighbours on boundaries.
        !! The memory access on boundaries in not-contigous.
        integer, value, intent(in)                               :: n1, n2, n3
        real(dp), dimension(-1:n1 + 2, -1:n2, -1:n3), intent(in) :: u
        real(dp), dimension(n1, n2, n3), intent(out)             :: du

        integer  :: i, j, k, l, xn(5), yn(5), zn(5)
        real(dp) :: ux(5), uy(5), uz(5)

        i = (blockIdx%x - 1) * blockDim%x + threadIdx%x
        j = (blockIdx%y - 1) * blockDim%y + threadIdx%y
        k = (blockIdx%z - 1) * blockDim%z + threadIdx%z

        do l = -2, 2
            xn(l + 3) = mod(i + l + n1 - 1, n1) + 1
            yn(l + 3) = mod(j + l + n2 - 1, n2) + 1
            zn(l + 3) = mod(k + l + n3 - 1, n3) + 1
        end do

        ux = u_d(xn, j, k)
        uy = u_d(i, yn, k)
        uz = u_d(i, j, zn)

        du_d(i, j, k) = (ux(1) + uy(1) + uz(1)) * factors_d(1) &
                        + (ux(2) + uy(2) + uz(2)) * factors_d(2) &
                        + (ux(3) + uy(3) + uz(3)) * factors_d(3) &
                        + (ux(4) + uy(4) + uz(4)) * factors_d(4) &
                        + (ux(5) + uy(5) + uz(5)) * factors_d(5)
    end subroutine divergence_mod

    attributes(global) subroutine divergence_halo(n1, n2, n3, u, du)
        !! divergence_mod :
        !! Compute the divergence using mod function to access neighbours on boundaries.
        !! The memory access on boundaries in not-contigous.
        integer, value, intent(in)                               :: n1, n2, n3
        real(dp), dimension(-1:n1 + 2, -1:n2, -1:n3), intent(in) :: u
        real(dp), dimension(n1, n2, n3), intent(out)             :: du

        integer  :: i, j, k

        i = (blockIdx%x - 1) * blockDim%x + threadIdx%x
        j = (blockIdx%y - 1) * blockDim%y + threadIdx%y
        k = (blockIdx%z - 1) * blockDim%z + threadIdx%z
        du(i, j, k) = (u(i - 2, j, k) + u(i, j - 2, k) + u(i, j, k - 2)) * factors_d(1) + &
                      (u(i - 1, j, k) + u(i, j - 1, k) + u(i, j, k - 1)) * factors_d(2) + &
                      (u(i, j, k) + u(i, j, k) + u(i, j, k)) * factors_d(3) + &
                      (u(i + 1, j, k) + u(i, j + 1, k) + u(i, j, k + 1)) * factors_d(4) + &
                      (u(i + 2, j, k) + u(i, j + 2, k) + u(i, j, k + 2)) * factors_d(5)
    end subroutine divergence_halo

    subroutine cuda_check_stat(istat)
        integer, intent(in) :: istat
        if (istat /= cudaSuccess) then
            write (*, "('Error Code:',i0,': ')") istat
            write (*, *) cudaGetErrorString(istat)
            stop
        end if
    end subroutine cuda_check_stat

end module finiteDifference3d
