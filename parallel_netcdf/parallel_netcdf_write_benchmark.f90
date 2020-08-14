program parallel_netcdf_write_benchmark

  use netcdf, only: nf90_independent

  implicit none

  !
  ! This benchmark writes out a 2D-decomposed 3D array in parallel in a single
  ! or multiple (one per rank) netCDF files, mimicking the behaviour of many
  ! 2D+1D simulation codes. Performance is measured using wall clock timings
  ! for netCDF and an explicit MPI barrier to catch load imbalance caused by
  ! IO resource contention.
  !
  ! Note that you will need to run this benchmark on nprocx*nprocy processors,
  ! domain decomposition is not automatic.
  !

  type config
     ! Seed for random number generator
     integer :: rng_seed = 123456

     ! Subdomain dimensions - adjusted to 1GB of double precision data per rank
     integer :: nx = 1024, ny = 1024, nz = 128

     ! Number of processes in xy direction
     integer :: nprocx = 8, nprocy = 5

     ! Chunking flag and sizes
     logical :: use_chunking = .false.
     integer :: chunk_nx = 1024, chunk_ny = 512, chunk_nz = 64
     !  integer, parameter :: chunk_nx = 1300, chunk_ny = 500, chunk_nz = 30

     ! Parallel access method (MPI-IO independent or collective)
     integer :: paraccess = nf90_independent

     ! Use one file per rank
     logical :: use_multiple_files = .false.

     character(len=255) :: filename = 'nc_par_test.nc'
  end type config

  call run_benchmark()

contains

  subroutine handle_err(status)
    use netcdf, only: nf90_strerror
    implicit none
    integer, intent (in) :: status
    print '(A)', trim(nf90_strerror(status))
    stop
  end subroutine handle_err

  ! --------------------------------------------------------------------------

  subroutine get_cmd_arguments(myrank, bm_setup)
    use netcdf, only: nf90_collective
    implicit none
    integer, intent(in) :: myrank
    type(config), intent(out) :: bm_setup
    integer :: i
    character(len=255) :: cmd_argument

    ! Defaults
    bm_setup = config()

    ! Read command line arguments
    if (command_argument_count() > 0) then
      do i = 1,command_argument_count()
        call get_command_argument(number=i, value=cmd_argument)
        if (i == command_argument_count()) then
          bm_setup%filename = trim(cmd_argument)
        else if (trim(cmd_argument) == "--collective") then
          bm_setup%paraccess = nf90_collective
        else if (trim(cmd_argument) == "--chunking") then
          bm_setup%use_chunking = .true.
        else if (trim(cmd_argument) == "--multiple") then
          bm_setup%use_multiple_files = .true.
        end if
      end do
    else
      if (myrank == 0) print '(A)', &
        'Usage: parallel_netcdf_write_benchmark.x [--collective] [--chunking] [--multiple] FILE'
      stop
    end if
  end subroutine get_cmd_arguments

  ! --------------------------------------------------------------------------

  subroutine create_nc_file(bm_setup, myrank, comm, ncfileid)
    use mpi, only: MPI_INFO_NULL
    use netcdf, only: nf90_create, nf90_clobber, nf90_netcdf4, nf90_mpiio, nf90_noerr
    implicit none
    type(config), intent(in) :: bm_setup
    integer, intent(in) :: myrank, comm
    integer, intent(out) :: ncfileid
    character(len=255) :: rank_suffix
    integer :: stat

    if (bm_setup%use_multiple_files) then
      write (rank_suffix,'(I0.3)') myrank
      stat = nf90_create(trim(bm_setup%filename) // '_'// trim(rank_suffix), &
                         IOR(nf90_clobber, nf90_netcdf4), ncfileid)
      if (stat /= nf90_noerr) call handle_err(stat)
    else
      stat = nf90_create(trim(bm_setup%filename), IOR(nf90_clobber, IOR(nf90_netcdf4, nf90_mpiio)), &
                         ncfileid, comm=comm, info=MPI_INFO_NULL)
      if (stat /= nf90_noerr) call handle_err(stat)
    end if

  end subroutine create_nc_file

  ! --------------------------------------------------------------------------

  subroutine setup_nc_file(bm_setup, ncfileid, datavarId)
    use netcdf, only: nf90_def_dim, nf90_put_att, nf90_def_var, nf90_var_par_access, &
                      nf90_double, nf90_global, nf90_noerr
    implicit none
    type(config), intent(in) :: bm_setup
    integer, intent(in) :: ncfileid
    integer, intent(out) :: datavarId
    integer :: stat, dims(3), dimids(3)

    ! Set variable size to subdomain or domain size
    if (bm_setup%use_multiple_files) then
      dims = (/bm_setup%nx, bm_setup%ny, bm_setup%nz/)
    else
      dims = (/bm_setup%nprocx*bm_setup%nx, bm_setup%nprocy*bm_setup%ny, bm_setup%nz/)
    end if

    stat = nf90_def_dim(ncfileid, "dim1", dims(1), dimids(1))
    if (stat /= nf90_noerr) call handle_err(stat)

    stat = nf90_def_dim(ncfileid, "dim2", dims(2), dimids(2))
    if (stat /= nf90_noerr) call handle_err(stat)

    stat = nf90_def_dim(ncfileid, "dim3", dims(3), dimids(3))
    if (stat /= nf90_noerr) call handle_err(stat)

    stat = nf90_put_att(ncfileid, nf90_global, "nx", bm_setup%nx)
    if (stat /= nf90_noerr) call handle_err(stat)

    stat = nf90_put_att(ncfileid, nf90_global, "ny", bm_setup%ny)
    if (stat /= nf90_noerr) call handle_err(stat)

    stat = nf90_put_att(ncfileid, nf90_global, "nz", bm_setup%nz)
    if (stat /= nf90_noerr) call handle_err(stat)

    stat = nf90_put_att(ncfileid, nf90_global, "nprocx", bm_setup%nprocx)
    if (stat /= nf90_noerr) call handle_err(stat)

    stat = nf90_put_att(ncfileid, nf90_global, "nprocy", bm_setup%nprocy)
    if (stat /= nf90_noerr) call handle_err(stat)

    stat = nf90_put_att(ncfileid, nf90_global, "seed", bm_setup%rng_seed)
    if (stat /= nf90_noerr) call handle_err(stat)

    if (bm_setup%use_chunking) then
      stat = nf90_def_var(ncfileid, "v1", nf90_double, dimids, datavarId, &
                          chunksizes=(/bm_setup%chunk_nx,bm_setup%chunk_ny,bm_setup%chunk_nz/))
    else
      stat = nf90_def_var(ncfileid, "v1", nf90_double, dimids, datavarId)
    end if
    if (stat /= nf90_noerr) call handle_err(stat)

    ! Parallel access mode only makes sense when writing to a single file
    if (.not. bm_setup%use_multiple_files) then
      stat = nf90_var_par_access(ncfileid, datavarId, bm_setup%paraccess)
      if (stat /= nf90_noerr) call handle_err(stat)
    end if

  end subroutine setup_nc_file

  ! --------------------------------------------------------------------------

  subroutine generate_data(bm_setup, arr)
  implicit none
  type(config), intent(in) :: bm_setup
  real(kind=8), intent(out) :: arr(:,:,:)
  integer :: seed_dim
  integer, allocatable :: seed(:)

    call random_seed(size=seed_dim)
    allocate(seed(seed_dim))
    seed = bm_setup%rng_seed
    call random_seed(put=seed)
    deallocate(seed)
    call random_number(arr)

  end subroutine generate_data

  ! --------------------------------------------------------------------------

  subroutine run_benchmark

  use mpi, only: MPI_Init, MPI_Comm_size, MPI_Comm_rank, MPI_Barrier, MPI_Finalize, &
                 MPI_COMM_WORLD
  use netcdf, only: nf90_put_var, nf90_enddef, nf90_close, nf90_strerror, &
                    nf90_independent, nf90_collective, nf90_double, nf90_noerr

  implicit none

  type(config) :: bm_setup
  integer :: ierr, status, ncid
  integer :: mpi_size, mpi_rank, cart_comm, coord3D(3)
  integer(kind=4) :: varId, start(3), count(3)
  real(kind=8), allocatable :: data(:,:,:)
  integer(kind=8) :: nbytes_written, starttime(7), count_rate, count_max

  !
  ! Initialise MPI and set up Cartesian topology
  !

  call MPI_Init(ierr)

  call MPI_Comm_size(MPI_COMM_WORLD, mpi_size, ierr);
  call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, ierr);

  if (mpi_rank == 0) then
    print *
    print '(A)', '*** Parallel netCDF write benchmark ***'
    print *
  end if

  ! MPI domain decomposition is fixed for simplicity
  if (mpi_size /= bm_setup%nprocx*bm_setup%nprocy) then
    print '(2(A,X,I3,X))', 'ERROR - requested', bm_setup%nprocx*bm_setup%nprocy, &
          'processes but running on', mpi_size
    stop
  end if

  call MPI_Cart_create(MPI_COMM_WORLD, 3, (/bm_setup%nprocx, bm_setup%nprocy, 1/), &
                       (/0,0,0/), 0, cart_comm, ierr)
  call MPI_Cart_coords(cart_comm, mpi_rank, 3, coord3D, ierr)

  !
  ! Configuration
  !

  call get_cmd_arguments(mpi_rank, bm_setup)

  if (mpi_rank == 0) then
    print '(A30,X,I4)', 'Number of MPI ranks:', bm_setup%nprocx*bm_setup%nprocy
    print '(A30,X,2(I4,X))', 'Domain decomposition:', bm_setup%nprocx, bm_setup%nprocy
    if (bm_setup%use_multiple_files) then
      print '(A30,X,A)', 'Single file output:', 'NO'
      print '(A30,X,A)', 'Filename:', trim(bm_setup%filename) // '_<rank>'
    else
      print '(A30,X,A)', 'Single file output:', 'YES'
      print '(A30,X,A)', 'Filename:', trim(bm_setup%filename)
      if (bm_setup%paraccess == nf90_collective) then
        print '(A30,X,A)','MPI-IO mode:', 'COLLECTIVE'
      else
        print '(A30,X,A)','MPI-IO mode:', 'INDEPENDENT'
      end if
    end if
    print '(A30,X,3(I5,X))', 'Domain dimensions:', bm_setup%nprocx*bm_setup%nx, &
                                                   bm_setup%nprocy*bm_setup%ny, bm_setup%nz
    print '(A30,X,3(I5,X))', 'Subdomain dimensions:', bm_setup%nx, bm_setup%ny, bm_setup%nz
    if (bm_setup%use_chunking) then
      print '(A30,X,A)', 'Chunking:', 'YES'
      print '(A30,X,3(I5,X))', 'Chunksizes:', bm_setup%chunk_nx, bm_setup%chunk_ny, bm_setup%chunk_nz
    else
      print '(A30,X,A)', 'Chunking:', 'NO'
    end if
  end if

  !
  ! Set up output netCDF file
  !

  call create_nc_file(bm_setup, mpi_rank, cart_comm, ncid)

  call setup_nc_file(bm_setup, ncid, varId)

  if (mpi_rank == 0) then
    print *
    print '(A)', 'Writing netCDF metadata...'
  end if

  ! It is not normally necessary to call nf90_enddef for netCDF4/HDF5 files, but
  ! we'll do it anyway for measuring metadata write times
  call system_clock(starttime(1), count_rate, count_max)
  status = nf90_enddef(ncid)
  if (status /= nf90_noerr) call handle_err(status)
  call system_clock(starttime(2), count_rate, count_max)
  if (mpi_rank == 0) then

  end if

  !
  ! Create output data
  !

  if (mpi_rank == 0) print '(A)', 'Generating data...'

  allocate(data(bm_setup%nx,bm_setup%ny,bm_setup%nz))
  call generate_data(bm_setup, data)

  !
  ! Write output data
  !

  if (mpi_rank == 0) print '(A)', 'Writing variable data...'

  ! Set up subdomain location and size in netCDF variable
  if (bm_setup%use_multiple_files) then
    start = (/1,1,1/)
  else
    start = (/coord3D(1)*bm_setup%nx+1,coord3D(2)*bm_setup%ny+1,1/)
  end if
  count = (/bm_setup%nx,bm_setup%ny,bm_setup%nz/)

  call system_clock(starttime(3), count_rate, count_max)
  status = nf90_put_var(ncid, varId, data, start=start, count=count)
  if (status /= nf90_noerr) call handle_err(status)
  call system_clock(starttime(4), count_rate, count_max)

  ! Make sure that all ranks have finished and measure potential load imbalance
  call MPI_Barrier(cart_comm, ierr)
  call system_clock(starttime(5), count_rate, count_max)

  ! Measure potential wait times
  status = nf90_close(ncid)
  if (status /= nf90_noerr) call handle_err(status)
  call system_clock(starttime(6), count_rate, count_max)

  ! Measure potential wait times
  call MPI_Finalize(ierr)
  call system_clock(starttime(7), count_rate, count_max)

  if (mpi_rank == 0) then
    print *
    print '(A)', 'TIMING RESULTS:'
    print '(A30,X,E12.6)', 'nf90_enddef timing [s]:',     dble(starttime(2)-starttime(1))/dble(count_rate)
    print '(A30,X,E12.6)', 'nf90_put_var timing [s]:',    dble(starttime(4)-starttime(3))/dble(count_rate)
    print '(A30,X,E12.6)', 'mpi_barrier timing [s]:',     dble(starttime(5)-starttime(4))/dble(count_rate)
    print '(A30,X,E12.6)', 'nf90_close timing [s]:',      dble(starttime(6)-starttime(5))/dble(count_rate)
    print '(A30,X,E12.6)', 'mpi_finalize timing [s]:',    dble(starttime(7)-starttime(6))/dble(count_rate)
    print '(A30,X,E12.6)', 'TOTAL WRITE timing [s]:', dble(starttime(7)-starttime(3))/dble(count_rate)
    print *
    nbytes_written = size(data, kind=8) * storage_size(data, kind=8)*int(mpi_size, kind=8)/8_8
    print '(A30,X,I16,X,A,F5.1,X,A)', 'Total bytes written:', nbytes_written, '(', &
          dble(nbytes_written)/1024.**3, 'GiB)'
    print '(A30,X,F10.3)', 'Effective data rate [MiB/s]:', &
          dble(nbytes_written)/1024./1024./dble(starttime(7)-starttime(3))*dble(count_rate)
  end if

  deallocate(data)

  end subroutine run_benchmark

end program parallel_netcdf_write_benchmark
