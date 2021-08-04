module derivative_mod_gpu
    use cudafor
    implicit none
    contains

        attributes(global) subroutine derivative_halo_gpu_kernel(u,du,n1, n2, n3,factors)
        !! divergence_mod :
        !! Compute the divergence using mod function to access neighbours on boundaries.
        !! The memory access on boundaries in not-contigous.
        real*8::u(-1:n1+2,-1:n2+2,-1:n3+2)
        real*8::du(n1,n2,n3),factors(5)
        integer, value, intent(in) :: n1, n2, n3

        integer  :: i, j,k, l, xn(5),yn(5),zn(5)
        real*8 :: ux(5), uy(5), uz(5)

        i = (blockIdx%x-1) * blockDim%x + threadIdx%x
        j = (blockIdx%y-1) * blockDim%y + threadIdx%y
        k = (blockIdx%z-1) * blockDim%z + threadIdx%z
  
!!        print *,i,j,k

          du(i, j, k) = (u(i - 2, j, k) + u(i, j - 2, k) + u(i, j, k - 2)) * factors(1) + &
                                  (u(i - 1, j, k) + u(i, j - 1, k) + u(i, j, k - 1)) * factors(2) + &
                                  (u(i, j, k)     + u(i, j, k)     + u(i, j, k)    ) * factors(3) + &
                                  (u(i + 1, j, k) + u(i, j + 1, k) + u(i, j, k + 1)) * factors(4) + &
                                  (u(i + 2, j, k) + u(i, j + 2, k) + u(i, j, k + 2)) * factors(5)



    end subroutine derivative_halo_gpu_kernel 

        subroutine derivative_halo_gpu_dev(u,du,n1,n2,n3,factors)
                integer,value,intent(in)::n1,n2,n3
                real*8::u(-1:n1+2,-1:n2+2,-1:n3+2)
                real*8::du(n1,n2,n3),factors(5)
                integer::block_size(3)
                type(dim3) :: grid_sp,block_sp
                block_size(1)=32;block_size(2)=4;block_size(3)=4
                block_sp   = dim3(block_size(1),block_size(2),block_size(3))
                grid_sp    = dim3(n1/block_sp%x,n2/block_sp%y,n3/block_sp%z)
                !$acc host_data use_device(u,du,factors)
                !!$acc data present(u,du,factors)
                call derivative_halo_gpu_kernel<<<grid_sp,block_sp>>>(u,du,n1,n2,n3,factors)
!! CHECK THIS::                call cuda_get_last_error()
                !$acc wait
                !!$acc end data
                !$acc end host_data
        end subroutine derivative_halo_gpu_dev







end module derivative_mod_gpu 
