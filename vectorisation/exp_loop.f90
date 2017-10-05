program exp_loop

  ! Test program that computes a large number of exponential functions
  ! using a vectorisable loop

  implicit none

  ! Work array size - needs to be large
  integer(kind=4), parameter :: n = 10**7

  ! Repeat computation to allow for runtime variance
  integer(kind=4), parameter :: nrepeat = 10

  real(kind=4) :: data(n), deltat
  integer(kind=8) :: tstart, tstop, trate
  integer(kind=4) :: repeat

  deltat = HUGE(deltat)

  do repeat = 1, nrepeat
     call random_number(data)
     call system_clock(tstart, trate)
     data(:) = exp(data(:))
     call system_clock(tstop)

     ! Capture best result
     deltat = min(deltat, real(tstop-tstart)/real(trate))
  end do

  write (*, '(A)') "exp_loop:"
  write (*, '(A,X,I0,X,A)') "Repeated computation", nrepeat, "times to capture best performance"
  write (*, '(A,X,E10.3,X,A)') "Computed", real(n), "exponentials"
  write (*, '(A,X,E10.3)') "Shortest time spent on computing exp function [s]:", deltat
  write (*, '(A,X,E10.3,X,A)') "Best performance:", real(n)/deltat, "exp per sec"

end program exp_loop
