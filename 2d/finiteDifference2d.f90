#define XWORK (256)
#define YWORK (64)

module finiteDifference2d

    use cudafor

    implicit none

    integer,  parameter :: dp     = kind(1.d0)
    real(dp), parameter :: two_pi = 8.d0*atan(1.d0)
    !! factors

    integer  :: nGrid, numIterations, blockSize(2), itime
    real(dp) :: length, dl, dl_inv
    !! Input Parameters

    real(dp)           :: factors(5)
    real(dp), constant :: factors_d(5)

    real(dp), dimension(:, :), allocatable         :: u, du, du_exact
    real(dp), dimension(:, :), allocatable, device :: u_d, du_d

contains

    subroutine get_block_size()
        character(len=128) :: cmdArg
        call get_command_argument(1,cmdArg)
        read(cmdArg,*) blockSize(1)
        call get_command_argument(2,cmdArg)
        read(cmdArg,*) blockSize(2)
    end subroutine get_block_size

    subroutine box_init(filename)
        character(*), intent(in) :: filename
        namelist /box_nml/ nGrid, length, numIterations

        open(unit=10,file=filename,status='old')
        read(10,nml=box_nml)
        close(10)

        length = length*two_pi
        dl     = length/nGrid
        dl_inv = 1.d0/dl
    end subroutine box_init

    subroutine set_factors()
        factors   = dl_inv*(1.d0/12.d0)*[ 1, -8, 0, 8, -1 ]
        factors_d = factors
    end subroutine set_factors

    subroutine allocate_arrays()
        allocate(u(nGrid,nGrid))
        allocate(du(nGrid,nGrid))
        allocate(du_exact(nGrid,nGrid))

        allocate(u_d(nGrid,nGrid))
        allocate(du_d(nGrid,nGrid))
        u        = 0.d0
        du       = 0.d0
        du_exact = 0.d0
    end subroutine allocate_arrays

    subroutine allocate_arrays_padded()
        !! Padding around 2d arrays required for one of the algorithm.
        allocate(u(-1:nGrid+2,-1:nGrid+2))
        allocate(du(-1:nGrid+2,-1:nGrid+2))
        allocate(du_exact(-1:nGrid+2,-1:nGrid+2))

        allocate(u_d(-1:nGrid+2,-1:nGrid+2))
        allocate(du_d(-1:nGrid+2,-1:nGrid+2))
        u        = 0.d0
        du       = 0.d0
        du_exact = 0.d0
    end subroutine allocate_arrays_padded

    subroutine write_device_array()
        integer :: i,j
        u = u_d
        do j = -1, nGrid+2
            do i = -1, nGrid+2
                write(*,20,advance='no') u(i,j)
            enddo
            write(*,*)
        enddo
        20      format (1f10.4)
    end subroutine write_device_array

    subroutine cuda_check_stat(istat)
        !! Subroutine to check cuda errors.
        integer, intent(in) :: istat
        if (istat /= cudaSuccess) then
            write(*, "('Error Code:',i0,': ')") istat
            write(*,*) cudaGetErrorString(istat)
            stop
        end if
    end subroutine cuda_check_stat

    subroutine initial_condition()
        integer  :: i,j
        real(dp) :: x,y

        u(:,:)        = 0.d0
        du_exact(:,:) = 0.d0

        do j = 1,nGrid
            y = dl*(j-1)
            do i = 1,nGrid
                x = dl*(i-1)
                u(i,j)        = sin(x)+sin(y)
                du_exact(i,j) = cos(x)+cos(y)
            end do
        end do

        !! Copy to device
        u_d(:,:) = u(:,:)
        !! Set the device output to 0
        du_d(:,:) = 0.d0
    end subroutine initial_condition

    !! Implementations begin here.

    !! Implementation 00
    attributes(global) subroutine divergence_mod(nGrid)
        integer, value, intent(in) :: nGrid

        integer  :: i, j, l, xn(5),yn(5)
        real(dp) :: ux(5), uy(5)

        i = (blockIdx%x-1) * blockDim%x + threadIdx%x
        j = (blockIdx%y-1) * blockDim%y + threadIdx%y

        do l = -2, 2
            xn(l+3) = mod(i+l+nGrid-1,nGrid) + 1
            yn(l+3) = mod(j+l+nGrid-1,nGrid) + 1
        end do

        ux = u_d(xn,j)
        uy = u_d(i,yn)

        du_d(i,j) = (ux(1)+uy(1))*factors_d(1) + (ux(2)+uy(2))*factors_d(2) &
                  + (ux(3)+uy(3))*factors_d(3) + (ux(4)+uy(4))*factors_d(4) &
                  + (ux(5)+uy(5))*factors_d(5)          
    end subroutine divergence_mod

    !! Implementation 01
    attributes(global) subroutine divergence_padding()
        integer  :: i, j
        real(dp) :: ux(5), uy(5)

        i = (blockIdx%x-1) * blockDim%x + threadIdx%x
        j = (blockIdx%y-1) * blockDim%y + threadIdx%y

        ux = u_d(i-2:i+2,j)
        uy = u_d(i,j-2:j+2)

        du_d(i,j) = (ux(1)+uy(1))*factors_d(1) + (ux(2)+uy(2))*factors_d(2) &
                  + (ux(3)+uy(3))*factors_d(3) + (ux(4)+uy(4))*factors_d(4) &
                  + (ux(5)+uy(5))*factors_d(5)          
    end subroutine divergence_padding

    !! Helper functions for divergence_padding.
    subroutine copy_to_padding()
        !! Copy data to padding.
        type(dim3) :: blockSp, gridSp
        !! Copy first the contiguous x direction, cuda support mpi like copying,
        !! so we can just specify the start of read and write locations and total
        !! number of elements to use.
        call cuda_check_stat(cudaMemCpy(u_d(-1,-1),u_d(-1,nGrid-1), 2*nGrid + 8, &
            cudaMemcpyDeviceToDevice))
        call cuda_check_stat(cudaMemCpy(u_d(-1,nGrid+1),u_d(-1,1) , 2*nGrid + 8, &
            cudaMemcpyDeviceToDevice)) 

        !! Copy the non contiguous y direction
        blockSp = dim3(256,1,1)
        gridSp  = dim3(nGrid/blockSp%x,1,1)
        call copy_y_padding <<< gridSp,blockSp >>> (nGrid)
    end subroutine copy_to_padding

    attributes(global) subroutine copy_y_padding(nGrid)
        !! Need a kernel to copy non-contiguous y-data.
        integer, value, intent(in) :: nGrid
        integer :: i
        i = (blockIdx%x-1) * blockDim%x + threadIdx%x
        u_d(-1,i)   = u_d(nGrid-1,i)
        u_d(0,i)    = u_d(nGrid,i)
        u_d(nGrid+1,i) = u_d(1,i)
        u_d(nGrid+2,i) = u_d(2,i)
    end subroutine copy_y_padding

    !! Implementation 02
    attributes(global) subroutine divergence_shared(nGrid)
        !! Compute derivatives by bringing values to shared memory.

        integer, value, intent(in) :: nGrid

        integer          :: iG, jG, i, j
        real(dp), shared :: u_ds(-1:blockDim%x+2,-1:blockDim%y+2)
        !! Shared memory

        i  = threadIdx%x
        j  = threadIdx%y
        iG = (blockIdx%x-1) * blockDim%x + threadIdx%x
        jG = (blockIdx%y-1) * blockDim%y + threadIdx%y

        u_ds(i,j) = u_d(iG,jG)
        !! now we copy paddings properly
        if (i <= 2) then
            !! check if block is not on x edge, if yes gotta fetch from edges.
            if (blockIdx%x == 1) then
                u_ds(i-2, j)  = u_d(nGrid+i-2,jG)            !! Left padding
                u_ds(i+blockDim%x,j) = u_d(iG+blockDim%x,jG) !! Right padding
            else if (blockIdx%x == gridDim%x) then
                u_ds(i-2, j)  = u_d(iG-2,jG)                 !! Left padding
                u_ds(i+blockDim%x,j) = u_d(i,jG)             !! Right padding
            else
                u_ds(i-2, j)  = u_d(iG-2,jG)                 !! Left padding
                u_ds(i+blockDim%x,j) = u_d(iG+blockDim%x,jG) !! Right padding
            end if
        end if

        if (j <= 2) then
            if (blockIdx%y == 1) then
                u_ds(i, j-2)  = u_d(iG,nGrid+j-2)            !! Up padding
                u_ds(i,j+blockDim%y) = u_d(iG,jG+blockDim%y) !! Down padding
            else if (blockIdx%y == gridDim%y) then
                u_ds(i, j-2)  = u_d(iG,jG-2)                 !! Up padding
                u_ds(i,j + blockDim%y) = u_d(iG,j)           !! Down padding
            else
                u_ds(i, j-2)  = u_d(iG,jG-2)                 !! Up padding
                u_ds(i,j+blockDim%y) = u_d(iG,jG+blockDim%y) !! Down padding
            end if
        endif

        call syncthreads()

        du_d(iG,jG) = ( u_ds(i-2,j) + u_ds(i,j-2) )*factors_d(1) &
                    + ( u_ds(i-1,j) + u_ds(i,j-1) )*factors_d(2) &
                    + ( u_ds(i  ,j) + u_ds(i,j  ) )*factors_d(3) &
                    + ( u_ds(i+1,j) + u_ds(i,j+1) )*factors_d(4) &
                    + ( u_ds(i+2,j) + u_ds(i,j+2) )*factors_d(5)
    end subroutine divergence_shared

    !! Implementation 03
    attributes(global) subroutine divergence_sharedx_registery(nGrid)
        !! Compute using shared memory and registers.
        !! Load x-data to shared memory, load y-data in registers.

        integer, value, intent(in) :: nGrid

        integer          :: i, j, iG, jW, jB

        real(dp)         :: front(2), center, behind(2)
        !! Put y-stuff on registers
        real(dp), shared :: u_ds(-1:blockDim%x+2)
        !! Shared memory for x direction

        i = threadIdx%x
        j = threadIdx%y
        !! Local thread index
        iG = (blockIdx%x-1) * blockDim%x + threadIdx%x
        !! Global index, will be used to access u_d/write du_d
        jB = YWORK*((blockIdx%y-1) * blockDim%y + threadIdx%y-1) + 1
        !! Global Index on which threads are going to work

        if (blockIdx%y == 1) then
            behind = [u_d(iG,nGrid-1),u_d(iG,nGrid)]
        else
            behind = [u_d(iG,jB-2),u_d(iG,jB-1)]
        end if
        center = u_d(iG,jB)
        front  = [u_d(iG,jB+1),u_d(iG,jB+2)]
        !! Preload the y-stencil for first point.
        !! Need to take care of the behind elements for lower most block properly.
        !! center and front can be loaded because threads always start at lowest index.

        !! now we copy x-paddings properly
        !! check if block is not on x edge, if yes gotta fetch from edges.
        !! Having if conditions outside YWORK loop reduces computational time,
        !! but increases code size.
        if (blockIdx%x == 1) then
            do jW = 1, YWORK
                u_ds(i) = center
                if (i <= 2) then
                    u_ds(i-2)  = u_d(nGrid+i-2,jB)             !! Left padding
                    u_ds(i+blockDim%x) = u_d(iG+blockDim%x,jB) !! Right padding
                end if
                call syncthreads()
                du_d(iG,jB) = ( u_ds(i-2) + behind(1))*factors_d(1) &
                                  + ( u_ds(i-1) + behind(2))*factors_d(2) &
                                  + ( u_ds(i  ) + center   )*factors_d(3) &
                                  + ( u_ds(i+1) + front(1) )*factors_d(4) &
                                  + ( u_ds(i+2) + front(2) )*factors_d(5)

                behind(1) = behind(2)
                behind(2) = center
                center    = front(1)
                front(1)  = front(2)
                jB     = jB + 1
                if (blockIdx%y == gridDim%y) then
                    front(2)  = u_d(iG,mod(jB+nGrid+1,nGrid)+1)
                else
                    front(2)  = u_d(iG,jB+2)
                end if
                call syncthreads()
            end do
        else if (blockIdx%x == gridDim%x) then
            do jW = 1, YWORK
                u_ds(i) = center
                if (i <= 2) then
                    u_ds(i-2)  = u_d(iG-2,jB)       !! Left padding
                    u_ds(i+blockDim%x) = u_d(i,jB)  !! Right padding
                end if
                call syncthreads()
                du_d(iG,jB) = ( u_ds(i-2) + behind(1))*factors_d(1) &
                                  + ( u_ds(i-1) + behind(2))*factors_d(2) &
                                  + ( u_ds(i  ) + center   )*factors_d(3) &
                                  + ( u_ds(i+1) + front(1) )*factors_d(4) &
                                  + ( u_ds(i+2) + front(2) )*factors_d(5)

                behind(1) = behind(2)
                behind(2) = center
                center    = front(1)
                front(1)  = front(2)
                jB     = jB + 1
                if (blockIdx%y == gridDim%y) then
                    front(2)  = u_d(iG,mod(jB+nGrid+1,nGrid)+1)
                else
                    front(2)  = u_d(iG,jB+2)
                end if
                call syncthreads()
            end do
        else
            do jW = 1, YWORK
                u_ds(i) = center
                if (i <= 2) then
                    u_ds(i-2)  = u_d(iG-2,jB)                  !! Left padding
                    u_ds(i+blockDim%x) = u_d(iG+blockDim%x,jB) !! Right padding
                end if
                call syncthreads()
                du_d(iG,jB) = ( u_ds(i-2) + behind(1))*factors_d(1) &
                                  + ( u_ds(i-1) + behind(2))*factors_d(2) &
                                  + ( u_ds(i  ) + center   )*factors_d(3) &
                                  + ( u_ds(i+1) + front(1) )*factors_d(4) &
                                  + ( u_ds(i+2) + front(2) )*factors_d(5)

                behind(1) = behind(2)
                behind(2) = center
                center    = front(1)
                front(1)  = front(2)
                jB     = jB + 1
                if (blockIdx%y == gridDim%y) then
                    front(2)  = u_d(iG,mod(jB+nGrid+1,nGrid)+1)
                else
                    front(2)  = u_d(iG,jB+2)
                end if
                call syncthreads()
            end do
        end if
    end subroutine divergence_sharedx_registery

end module finiteDifference2d

