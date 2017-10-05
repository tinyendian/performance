!
! This program tests cache performance.
!
! A work loop is repeated a number of times to simulate data reusage,
! which amplifies the effect of cache use on overall floating-point
! performance. Work load is an "axpy"-type computation that can be
! performed as a single FMA operation by the processor.
!
! Loop blocking can be used to test performance optimisation by fitting
! data into a cache level.
!
program cache

  implicit none

  !
  ! Program parameters
  !

  ! CPU parameters: clock rate in GHz, vector width
  real, parameter :: clock_rate = 3.1, vector_width = 8

  ! Floating-point precision
  integer, parameter :: realkind = 4

  ! Reuse factor (number of times work loop repeats)
  integer, parameter :: nreuse = 200

  ! Switch for loop blocking and block size
  ! Try L1, L2, L3 cache sizes divided by sizeof(real)
  logical, parameter :: blocking = .true.
  integer, parameter :: block_size = 32768/realkind

  ! Number of array size steps
  integer, parameter :: size_steps = 54

  !
  ! Program variables
  !

  integer :: k, nelements, nstats, size_step, decade
  integer(kind=8) :: nbyte
  real(kind=realkind), dimension(:), allocatable :: workarray
  real(kind=8) :: maxperformance, performance, peakperformance
  logical :: tick_warning

  !
  ! Run test
  !

  print '(A)', '***********************************************************************'
  print '(A)', 'Cache performance test                                                 '
  print '(A)', '***********************************************************************'
  print *
  print '(A)', '  Arr Size *  Memory  * Actual GFLOPS *  Peak GFLOPS  *   Efficiency  *' 

  ! Compute peak performance
  ! Clock rate * vector width * 2 (FMA = 2 FLOPs)
  peakperformance = clock_rate*vector_width*2.0d0

  ! Iterate array size, using 9 points per decade
  do size_step = 1, size_steps

    ! Compute decade, start with 10
    decade = 10**((size_step+8)/9)

    ! Compute array size
    nelements = (mod(size_step+8, 9)+1)*decade

    ! Compute number of repetitions for statistics, repetitions
    ! decrease with every decade. Need a lot of repetitions for
    ! the shortest array lengths.
    nstats = max(10, 1000000/decade)

    ! Compute number of bytes in array
    nbyte = realkind * nelements

    ! Limit work array size to 1 GB
    if (nbyte > 2**30) then
      write(*,*) "Array to large! Exit..."
      stop
    end if
	
    ! Get memory
    allocate(workarray(nelements))
	
	! Reset warning for short tick counts
	tick_warning = .false.
  
    ! Run worker function in a loop to gather statistics
    maxperformance = 0.d0
    do k = 1, nstats
      call compute(nelements, nreuse, blocking, block_size, workarray, &
	  & performance, tick_warning)
      maxperformance = max(maxperformance, performance)
    end do

    if (tick_warning) then
      print *, 'WARNING: at least one computation was shorter than 10 ticks!'
      print *, 'Performance measurement may not be reliable.'
    end if
    
    ! Print performance results
    print '(2(1X, I10), 1X, 3(7X, F8.4, 1X))', nelements, nbyte, maxperformance, peakperformance, &
         & maxperformance/peakperformance

    deallocate(workarray)

  end do

  !
  ! End of main program
  !

contains

  subroutine compute(nx, nrepeat, use_blocking, blocksize, work, performance, too_few_ticks)
    implicit none
	
    !
    ! Arguments
    !

    integer, intent(in) :: nx, nrepeat
    logical, intent(in) :: use_blocking
    integer, intent(in) :: blocksize
    real(kind=realkind), dimension(nx), intent(out) :: work
    real(kind=8), intent(out) :: performance
    logical, intent(inout) :: too_few_ticks

    !
    ! Local variables
    !

    integer :: i, j, k, elems_per_block, nblocks
    integer, dimension(2) :: irange
    integer(kind=8) :: startclock, stopclock, clockrate
    real(kind=8) :: timing

    ! Set block size to full array width if loop blocking is not used
    if (use_blocking) then
      elems_per_block = blocksize
    else
      elems_per_block = nx
    end if

    ! Fill array with random numbers
    call random_number(work)
 
    ! Get timer
    call system_clock(startclock, clockrate)
  
    ! Compute number of blocks
    nblocks = (nx-1)/elems_per_block + 1

    ! Block loop
    do k = 1, nblocks

      ! Compute index range for this block
      irange(1) = (k-1)*elems_per_block + 1
      irange(2)  = min(k*elems_per_block, nx)

      ! Simulate cache reuse: repeat computation
      do j = 1, nrepeat

        ! Work loop - vectorises
        do i = irange(1), irange(2)

          ! This should be run as one FMA instruction (2 FLOPs)
          work(i) = 2.3*work(i) + 7.13
        end do
      end do
    end do

    ! Get timer
    call system_clock(stopclock, clockrate)

    ! Check if number of ticks is large enough; return no result and
    ! raise flag if not.
    if ((stopclock - startclock) < 10) then
      too_few_ticks = .true.
      performance = 0.
    else

      ! Compute time spent in seconds
      timing = dble(stopclock - startclock)/dble(clockrate)
  
      ! Compute performance in GFLOPS
      performance = 2.d0*dble(nx)*dble(nrepeat)/(timing*1.0d9)

    end if

  end subroutine compute

end program cache
