program sgemm_loop

  ! Test program that computes a large number of matrix multiplications
  ! in single precision using BLAS

  implicit none

  ! Matrix rank
  integer(kind=4), parameter :: n = 1000

  ! Repeat computation to allow for runtime variance
  integer(kind=4), parameter :: nrepeat = 10

  real(kind=4) :: a(n,n), b(n,n), c(n,n), deltat
  integer(kind=8) :: tstart, tstop, trate
  integer(kind=4) :: repeat
  real(kind=4), parameter :: alpha = 1.0, beta = 1.0

  interface
     subroutine dgemm(TRANSA, TRANSB, M, N, K, ALPHA ,A , LDA, B, LDB, BETA, C, LDC)
       CHARACTER(1), intent(in) :: TRANSA, TRANSB
       integer(kind=4), intent(in) :: M, N, K
       real(kind=8), intent(in) :: ALPHA
       real(kind=8), intent(in) :: A(LDA, *)
       integer(kind=4), intent(in) :: LDA
       real(kind=8), intent(in) :: B(LDB, *)
       integer(kind=4), intent(in) :: LDB
       real(kind=8), intent(in) :: BETA
       real(kind=8), intent(inout) :: C(LDC, N)
       integer(kind=4), intent(in) :: LDC
     end subroutine dgemm
  end interface

  deltat = HUGE(deltat)

  do repeat = 1, nrepeat
     call random_number(a)
     call random_number(b)
     call random_number(c)
     call system_clock(tstart, trate)
     call sgemm('N', 'N', n, n, n, alpha, a, n, b, n, beta, c, n)
     call system_clock(tstop)

     ! Capture best result
     deltat = min(deltat, real(tstop-tstart)/real(trate))
  end do

  write (*, '(A)') "sgemm_loop:"
  write (*, '(A,X,I0,X,A)') "Repeated computation", nrepeat, "times to capture best performance"
  write (*, '(A,X,E10.3,X,A)') "Computed", real(2*n**3), "FLOPs"
  write (*, '(A,X,E10.3)') "Shortest time spent on computing sgemm [s]:", deltat
  write (*, '(A,X,E10.3,X,A)') "Best performance:", real(2*n**3)/deltat, "FLOPS"

end program sgemm_loop
