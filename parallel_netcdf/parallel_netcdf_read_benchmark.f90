program parallel_netcdf_read_benchmark

  use netcdf, only: nf90_independent

  implicit none

  !
  ! IMPORTANT: Running this benchmark requires running the "write" benchmark
  !            first to create netCDF file(s) with the correct setup! You
  !            will also need to use the same number of processors with
  !            which the file(s) have been written, see netCDF attributes
  !            "nprocx" and "nprocy" in the file(s).
  !
  ! This benchmark reads a 2D-decomposed 3D array in parallel from a single
  ! or multiple (one per rank) netCDF files, mimicking the behaviour of many
  ! 2D+1D simulation codes. Performance is measured using wall clock timings
  ! for netCDF and an explicit MPI barrier to catch load imbalance caused by
  ! IO resource contention.
  !


  ! See parallel_netcdf_write_benchmark for description of parameters
  type config
     integer :: rng_seed = -1
     integer :: nx=-1, ny=-1, nz=-1
     integer :: nprocx=-1, nprocy=-1
     logical :: use_chunking=.false.
     integer :: chunk_nx=-1, chunk_ny=-1, chunk_nz=-1
     integer :: paraccess = nf90_independent
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
        else if (trim(cmd_argument) == "--multiple") then
          bm_setup%use_multiple_files = .true.
        end if
      end do
    else
      if (myrank == 0) print '(A)', &
      'Usage: parallel_netcdf_read_benchmark.x [--collective] [--multiple] FILE'
      stop
    end if
  end subroutine get_cmd_arguments

  ! --------------------------------------------------------------------------

  subroutine get_input_file_specs(myrank, comm, bm_setup)
    use mpi
    use netcdf, only: nf90_open, nf90_get_att, nf90_close, nf90_inq_varid, &
                      nf90_inquire_variable, nf90_nowrite, nf90_noerr, &
                      nf90_global
    implicit none
    integer, intent(in) :: myrank, comm
    type(config), intent(inout) :: bm_setup
    integer :: status, ncid, ierr, chunksizes(3), varId
    character(len=255) :: filename
    logical :: contiguous_var

    ! Need to determine domain decomposition used by write benchmark
    ! Root rank reads metadata and broadcasts parameters to the other ranks
    if (myrank == 0) then
      if (bm_setup%use_multiple_files) then
        filename = trim(bm_setup%filename) // '_000'
      else
        filename = bm_setup%filename
      end if
      status = nf90_open(filename, nf90_nowrite, ncid)
      if (status /= nf90_noerr) call handle_err(status)

      status = nf90_get_att(ncid, nf90_global, "nprocx", bm_setup%nprocx)
      if (status /= nf90_noerr) call handle_err(status)

      status = nf90_get_att(ncid, nf90_global, "nprocy", bm_setup%nprocy)
      if (status /= nf90_noerr) call handle_err(status)

      status = nf90_get_att(ncid, nf90_global, "nx", bm_setup%nx)
      if (status /= nf90_noerr) call handle_err(status)

      status = nf90_get_att(ncid, nf90_global, "ny", bm_setup%ny)
      if (status /= nf90_noerr) call handle_err(status)

      status = nf90_get_att(ncid, nf90_global, "nz", bm_setup%nz)
      if (status /= nf90_noerr) call handle_err(status)

      status = nf90_get_att(ncid, nf90_global, "seed", bm_setup%rng_seed)
      if (status /= nf90_noerr) call handle_err(status)

      status = nf90_inq_varid(ncid, "v1", varId)
      if (status /= nf90_noerr) call handle_err(status)

      status = nf90_inquire_variable(ncid, varId, contiguous = contiguous_var)
      if(status /= nf90_noerr) call handle_err(status)

      bm_setup%use_chunking = .not. contiguous_var

      if (.not. contiguous_var) then
        status = nf90_inquire_variable(ncid, varId, chunksizes = chunksizes)
        if(status /= nf90_noerr) call handle_err(status)
      else
        chunksizes = 0
      end if

      bm_setup%chunk_nx = chunksizes(1)
      bm_setup%chunk_ny = chunksizes(2)
      bm_setup%chunk_nz = chunksizes(3)

      status = nf90_close(ncid)
      if (status /= nf90_noerr) call handle_err(status)

    end if

    call MPI_Bcast(bm_setup%nprocx, 1, MPI_INT, 0, comm, ierr)
    call MPI_Bcast(bm_setup%nprocy, 1, MPI_INT, 0, comm, ierr)
    call MPI_Bcast(bm_setup%nx, 1, MPI_INT, 0, comm, ierr)
    call MPI_Bcast(bm_setup%ny, 1, MPI_INT, 0, comm, ierr)
    call MPI_Bcast(bm_setup%nz, 1, MPI_INT, 0, comm, ierr)
    call MPI_Bcast(bm_setup%chunk_nx, 1, MPI_INT, 0, comm, ierr)
    call MPI_Bcast(bm_setup%chunk_ny, 1, MPI_INT, 0, comm, ierr)
    call MPI_Bcast(bm_setup%chunk_nz, 1, MPI_INT, 0, comm, ierr)
    call MPI_Bcast(bm_setup%rng_seed, 1, MPI_INT, 0, comm, ierr)

    if (bm_setup%chunk_nx > 0) then
      bm_setup%use_chunking = .true.
    else
      bm_setup%use_chunking = .false.
    end if

  end subroutine get_input_file_specs

  ! --------------------------------------------------------------------------

  subroutine open_nc_file(bm_setup, myrank, comm, ncid)
    use mpi, only: MPI_INFO_NULL
    use netcdf, only: nf90_open, nf90_nowrite, nf90_noerr, nf90_mpiio
    implicit none
    type(config), intent(in) :: bm_setup
    integer, intent(in) :: myrank, comm
    integer, intent(out) :: ncid
    character(len=255) :: rank_suffix
    integer :: status

    if (bm_setup%use_multiple_files) then
      write(rank_suffix,'(I0.3)') myrank
      status = nf90_open(trim(bm_setup%filename) // '_' // trim(rank_suffix), nf90_nowrite, ncid)
      if (status /= nf90_noerr) call handle_err(status)
    else
      status = nf90_open(trim(bm_setup%filename), IOR(nf90_nowrite, nf90_mpiio), ncid, &
                       comm=comm, info=MPI_INFO_NULL)
      if (status /= nf90_noerr) call handle_err(status)
    end if

  end subroutine open_nc_file

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

  subroutine run_benchmark()

    use mpi, only: MPI_Init, MPI_Comm_size, MPI_Comm_rank, MPI_COMM_WORLD
    use netcdf, only: nf90_get_var, nf90_inq_varid, nf90_var_par_access, &
                      nf90_close, nf90_collective, nf90_noerr

    implicit none

    type(config) :: bm_setup
    integer :: ierr, status, ncid
    integer :: mpi_size, mpi_rank, cart_comm, coord3D(3)
    integer(kind=4) :: varId, start(3), count(3)
    real(kind=8), allocatable :: data(:,:,:)
    real(kind=8) :: checksum
    integer(kind=8) :: nbytes_read, starttime(6), count_rate, count_max

    !
    ! Initialise MPI
    !

    call MPI_Init(ierr)

    call MPI_Comm_size(MPI_COMM_WORLD, mpi_size, ierr);
    call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, ierr);

    if (mpi_rank == 0) then
      print *
      print '(A)', '*** Parallel netCDF read benchmark ***'
      print *
    end if

    !
    ! Configuration
    !

    call get_cmd_arguments(mpi_rank, bm_setup)
    call get_input_file_specs(mpi_rank, MPI_COMM_WORLD, bm_setup)

    if (mpi_rank == 0) then
      print '(A30,X,I4)', 'Number of MPI ranks:', bm_setup%nprocx*bm_setup%nprocy
      print '(A30,X,2(I4,X))', 'Domain decomposition:', bm_setup%nprocx, bm_setup%nprocy
      if (bm_setup%use_multiple_files) then
        print '(A30,X,A)', 'Single file input:', 'NO'
        print '(A30,X,A)', 'Filename:', trim(bm_setup%filename) // '_<rank>'
      else
        print '(A30,X,A)', 'Single file input:', 'YES'
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
    ! Set up Cartesian topology
    !

    if (mpi_size /= bm_setup%nprocx*bm_setup%nprocy) then
      print '(2(A,X,I3,X))', 'ERROR - requested', bm_setup%nprocx*bm_setup%nprocy, &
           'processes but running on', mpi_size
      stop
    end if

    call MPI_Cart_create(MPI_COMM_WORLD, 3, (/bm_setup%nprocx, bm_setup%nprocy, 1/), &
                         (/0,0,0/), 0, cart_comm, ierr)
    call MPI_Cart_coords(cart_comm, mpi_rank, 3, coord3D, ierr)

    !
    ! Open netCDF file and set parallel access mode (if using single file IO)
    !

    call open_nc_file(bm_setup, mpi_rank, cart_comm, ncid)

    status = nf90_inq_varid(ncid, "v1", varId)
    if (status /= nf90_noerr) call handle_err(status)

    if (.not. bm_setup%use_multiple_files) then
      status = nf90_var_par_access(ncid, varId, bm_setup%paraccess)
      if (status /= nf90_noerr) call handle_err(status)
    end if

    if (mpi_rank == 0) then
      print *
      print '(A)', 'Generating checksum for variable data...'
    end if

    allocate(data(bm_setup%nx,bm_setup%ny,bm_setup%nz))
    call generate_data(bm_setup, data)
    checksum = sum(data)
    data = -1.0d0

    if (mpi_rank == 0) print '(A)', 'Reading variable data...'

    ! Set up 2D+1D subdomain location and size in netCDF variable
    if (bm_setup%use_multiple_files) then
      start = (/1,1,1/)
    else
      start = (/coord3D(1)*bm_setup%nx+1,coord3D(2)*bm_setup%ny+1,1/)
    end if
    count = (/bm_setup%nx,bm_setup%ny,bm_setup%nz/)

    call system_clock(starttime(1), count_rate, count_max)
    status = nf90_get_var(ncid, varId, data, start=start, count=count)
    if (status /= nf90_noerr) call handle_err(status)
    call system_clock(starttime(2), count_rate, count_max)

    ! Make sure that all ranks have finished
    call MPI_Barrier(cart_comm, ierr)
    call system_clock(starttime(3), count_rate, count_max)

    status = nf90_close(ncid)
    if (status /= nf90_noerr) call handle_err(status)
    call system_clock(starttime(4), count_rate, count_max)

    call MPI_Finalize(ierr)
    call system_clock(starttime(5), count_rate, count_max)

    ! Compare input data with checksum
    if (sum(data) /= checksum) then
      print '(A4,X,I3,X,A20)', 'Rank', mpi_rank, ': Input data INVALID'
    else
      if (mpi_rank == 0) then
        print *
        print '(A)', 'Input data VALID'
      end if
    end if

    if (mpi_rank == 0) then
      print *
      print '(A)', 'TIMING RESULTS:'
      print '(A30,X,E12.6)', 'nf90_get_var timing [s]:  ', dble(starttime(2)-starttime(1))/dble(count_rate)
      print '(A30,X,E12.6)', 'mpi_barrier timing [s]:   ', dble(starttime(3)-starttime(2))/dble(count_rate)
      print '(A30,X,E12.6)', 'nf90_close timing [s]:    ', dble(starttime(4)-starttime(3))/dble(count_rate)
      print '(A30,X,E12.6)', 'mpi_finalize timing [s]:  ', dble(starttime(5)-starttime(4))/dble(count_rate)
      print '(A30,X,E12.6)', 'TOTAL VAR READ timing [s]:', dble(starttime(5)-starttime(1))/dble(count_rate)
      print *
      nbytes_read = size(data, kind=8) * storage_size(data, kind=8)*int(mpi_size, kind=8)/8_8
      print '(A30,X,I16,X,A,F5.1,X,A)', 'Total bytes read:', nbytes_read, '(', &
           dble(nbytes_read)/1024.**3, 'GiB)'
      print '(A30,X,F10.3)', 'Effective data rate [MiB/s]:', &
           dble(nbytes_read)/1024./1024./dble(starttime(5)-starttime(1))*dble(count_rate)
    end if

    deallocate(data)

  end subroutine run_benchmark

end program parallel_netcdf_read_benchmark
