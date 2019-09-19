program fd2d

    use finiteDifference2d

    implicit none

    real(4)         :: simTime, cpuStart, cpuFinish
    type(dim3)      :: gridSp, blockSp
    type(cudaEvent) :: startEvent, stopEvent
    !! cuda events to measure timings

    call box_init('der2d.nml')
    call get_block_size()

    call cuda_check_stat(cudaEventCreate(startEvent))
    call cuda_check_stat(cudaEventCreate(stopEvent))
    !! create cuda_events

    call allocate_arrays()
    call initial_condition()
    call set_factors()

    blockSp   = dim3(blockSize(1),blockSize(2),1)
    gridSp    = dim3(nGrid/blockSp%x,nGrid/blockSp%y,1)

    call cuda_check_stat(cudaEventRecord(startEvent,0))
    do itime = 0, numIterations
        call divergence_mod<<<gridSp, blockSp>>>(nGrid)
    end do
    call cuda_check_stat(cudaEventRecord(stopEvent,0))
    call cuda_check_stat(cudaEventSynchronize(stopEvent))
    call cuda_check_stat(cudaEventElapsedTime(simTime,startEvent,stopEvent))
    du = du_d
    write(*,*) 'mod            :', blockSp, simTime, maxval(du_exact-du)

end program fd2d
