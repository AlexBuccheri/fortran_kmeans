! Note, I previously had this defined as src/run_kmeans.f90, and for some reason, use kmeans_m was failing
! Assume from some name mangling + shadowing?

! ./cmake-build-debug/run_kmeans
program run
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib

    use mpi_m, only: mpi_t, mpi_init, MPI_Finalize, distribute_elements_start_stop
    use grids_m,  only: linspace, linspace_to_grid
    use kmeans_m, only: weighted_kmeans
    implicit none

    ! MPI
    type(mpi_t), allocatable :: comm
    integer, allocatable :: start_end_indices(:, :), ir_indices(:) !< Indices for data distributed on the grid

    ! Files
    ! Root
    character(len=200) :: default_root = "/Users/alexanderbuccheri/Codes/kmeans/"
    character(len=200) :: file, root
    integer :: fid

    ! Grid
    integer,  parameter  :: n_dim = 2        !< Hard-coded to 2D
    integer              :: sampling(n_dim)
    real(dp)             :: ranges(n_dim, 2)
    integer              :: nr               !< Total number of grid points
    real(dp), allocatable :: grid(:, :)

    integer               :: n_centroids
    real(dp), allocatable :: weight(:)
    real(dp), allocatable :: centroids(:, :)

    integer :: ir, ic, niter, ierr, i
    logical :: verbose

    ! Initialise MPI env
    call mpi_init(ierr)
    comm = mpi_t() 

    write(*, *) 'comm got initialised'

    ! For convenience, have each process do its own IO

    ! Check if a command-line argument is provided
    if (command_argument_count() > 0) then
        ! Retrieve the first command-line argument
        call get_command_argument(1, root)

        ! Ensure trailing /
        if (root(len_trim(root):len_trim(root)) /= '/') then
            root = trim(root) // '/'
        endif

        write(*, *) 'Root directory:', trim(root)
    else
        write(*, *) 'No directory path provided. Defaulting to ', default_root
        root = trim(default_root)
    end if
    
    write(*, *) 'Parse command line'

    ! Parse grid settings
    file = trim(adjustl(root)) // "jupyter/grid_settings.dat"
    open(newunit=fid, file=trim(adjustl(file)))
    ! Header
    read(fid, *)
    read(fid, *) ranges(:, 1), sampling(1)
    read(fid, *) ranges(:, 2), sampling(2)
    close(fid)

    write(*, *) 'Parse grid settings'

    ! Generate grid from this
    nr = product(sampling)
    allocate(grid(n_dim, nr))
    call construct_2d_grid(sampling, ranges, grid)

    ! Print grid. Looks consistent with python setup
    if (comm%is_root()) then
        write(*, *) 'Outputting grid'
        do ir = 1, nr
            write(201, *) grid(:, ir)
        enddo
    endif

    ! Parse weights
    allocate(weight(nr))
    file = trim(adjustl(root)) // "jupyter/weights.dat"
    open(newunit=fid, file=trim(adjustl(file)))
    do ir = 1, nr
        read(fid, *) weight(ir)
    enddo
    close(fid)

    ! Distribute grid and weight - quantities on the grid
    allocate(start_end_indices(2, comm%np))
    start_end_indices = distribute_elements_start_stop(comm%np, nr)
    ir_indices = [(i, i= start_end_indices(1, comm%rank+1), start_end_indices(2, comm%rank+1))]

    ! Parse initial centroids
    file = trim(adjustl(root)) // "jupyter/init_centroids.dat"
    open(newunit=fid, file=trim(adjustl(file)))
    read(fid, *) n_centroids
    allocate(centroids(n_dim, n_centroids))
    do ic = 1, n_centroids
        read(fid, *) centroids(:, ic)
    enddo
    close(fid)

    niter = 50
    verbose = .true.

    call weighted_kmeans(comm, grid(:, ir_indices), weight(ir_indices), centroids, n_iter=niter, verbose=verbose)

    ! Output final centroids from root
    if (comm%is_root()) then
        file = trim(adjustl(root)) // "jupyter/final_centroids.dat"
        open(newunit=fid, file=trim(adjustl(file)), action='write')
        write(fid, *) n_centroids
        do ic = 1, n_centroids
            write(fid, *) centroids(:, ic)
        enddo
        close(fid)
    endif

    call mpi_finalize(ierr)

contains
    subroutine construct_2d_grid(sampling, ranges, grid)
        integer,  intent(in)               :: sampling(:)
        real(dp), intent(in)               :: ranges(:, :)
        real(dp), allocatable, intent(out) :: grid(:, :)
        real(dp), allocatable  :: x(:), y(:)

        integer :: nx, ny
        
        nx = sampling(1)
        ny = sampling(1)
        allocate(x(nx), y(ny), grid(2, nx * ny))
        call linspace(ranges(1, 1), ranges(2, 1), nx, x)
        call linspace(ranges(1, 2), ranges(2, 2), ny, y)
        call linspace_to_grid(x, y, grid)

    end subroutine construct_2d_grid

end program run