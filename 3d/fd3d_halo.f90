program fd3d

    use finiteDifference3d

    implicit none

    real(4)         :: simTime, cpuStart, cpuFinish
    type(dim3)      :: grid_sp, block_sp
    type(cudaEvent) :: startEvent, stopEvent
    !! cuda events to measure timings

    call box_init('der3d.nml')
    ! call get_block_size()

    call cuda_check_stat(cudaEventCreate(startEvent))
    call cuda_check_stat(cudaEventCreate(stopEvent))
    !! create cuda_events

    call allocate_arrays()
    call initial_condition()
    call set_factors()

    block_size = [32, 4, 4]
    block_sp   = dim3(block_size(1),block_size(2),block_size(3))
    grid_sp    = dim3(n(1)/block_sp%x,n(2)/block_sp%y,n(3)/block_sp%z)

    call cuda_check_stat(cudaEventRecord(startEvent,0))
    do itime = 0, num_iters
        call divergence_halo<<<grid_sp, block_sp>>>(n(1),n(2),n(3),u_d,du_d)
    end do
    call cuda_check_stat(cudaEventRecord(stopEvent,0))
    call cuda_check_stat(cudaEventSynchronize(stopEvent))
    call cuda_check_stat(cudaEventElapsedTime(simTime,startEvent,stopEvent))
    du = du_d
    write(*,*) 'Timings : ', block_sp, simTime, maxval(du_exact-du)

end program fd3d
