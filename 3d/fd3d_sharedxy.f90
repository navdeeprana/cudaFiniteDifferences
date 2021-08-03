#define XWORK (32)
#define YWORK (4)
#define ZWORK (16)

program fd3d

    use finiteDifference3d

    implicit none

    integer         :: total_shared_mem
    real(4)         :: simTime, cpuStart, cpuFinish
    type(dim3)      :: grid_sp, block_sp
    type(cudaEvent) :: startEvent, stopEvent
    !! cuda events to measure timings

    call box_init('der3d.nml')
    call get_block_size()

    call cuda_check_stat(cudaEventCreate(startEvent))
    call cuda_check_stat(cudaEventCreate(stopEvent))
    !! create cuda_events

    call allocate_arrays()
    call initial_condition()
    call set_factors()

    block_sp = dim3(XWORK, YWORK, 1)
    grid_sp = dim3(n(1) / block_sp%x, n(2) / block_sp%y, n(3) / ZWORK)
    total_shared_mem = 8 * (block_sp%x + 4) * (block_sp%y + 4)

    call cuda_check_stat(cudaEventRecord(startEvent, 0))
    do itime = 0, num_iters
        call divergence_sharedxy <<< grid_sp, block_sp, total_shared_mem >>> (n(1), n(2), n(3), u_d, du_d)
    end do
    call cuda_check_stat(cudaEventRecord(stopEvent, 0))
    call cuda_check_stat(cudaEventSynchronize(stopEvent))
    call cuda_check_stat(cudaEventElapsedTime(simTime, startEvent, stopEvent))
    du = du_d
    write (*, *) 'Sharedxy :', XWORK, YWORK, ZWORK, simTime, maxval(du_exact - du)

end program fd3d
