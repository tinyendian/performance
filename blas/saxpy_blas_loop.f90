program saxpy_blas_loop

  ! Test program that computes a large number of ax+b expressions
  ! in single precision using BLAS

  implicit none

  ! Work array size - needs to be large
  integer(kind=4), parameter :: n = 10**8

  ! Repeat computation to allow for runtime variance
  integer(kind=4), parameter :: nrepeat = 10

  real(kind=4) :: datax(n), datay(n), deltat
  integer(kind=8) :: tstart, tstop, trate
  integer(kind=4) :: repeat
  real(kind=4), parameter :: a = 0.7

  interface
     subroutine saxpy(N, SA, SX, INCX, SY, INCY)
       integer(kind=4), intent(in) :: N, INCX, INCY
       real(kind=4), intent(in) :: SA, SX(1 + (N-1)*abs(INCX))
       real(kind=4), intent(inout) :: SY(1 + (N-1)*abs(INCY))
     end subroutine saxpy
  end interface

  deltat = HUGE(deltat)

  do repeat = 1, nrepeat
     call random_number(datax)
     call random_number(datay)
     call system_clock(tstart, trate)
     call saxpy(n, a, datax, 1_4, datay, 1_4)
     call system_clock(tstop)

     ! Capture best result
     deltat = min(deltat, real(tstop-tstart)/real(trate))
  end do

  write (*, '(A)') "saxpy_loop:"
  write (*, '(A,X,I0,X,A)') "Repeated computation", nrepeat, "times to capture best performance"
  write (*, '(A,X,E10.3,X,A)') "Computed", real(2*n), "FLOPs"
  write (*, '(A,X,E10.3)') "Shortest time spent on computing saxpy [s]:", deltat
  write (*, '(A,X,E10.3,X,A)') "Best performance:", real(2*n)/deltat, "FLOPS"

end program saxpy_blas_loop
