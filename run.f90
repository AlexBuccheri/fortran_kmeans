! Note, I previously had this defined as src/run_kmeans.f90, and for some reason, use kmeans_m was failing
! Assume from some name mangling + shadowing?

! ./cmake-build-debug/run_kmeans
program run
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib

    use grids_m,  only: linspace, linspace_to_grid
    use kmeans_m, only: weighted_kmeans
    implicit none

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

    integer :: ir, ic, niter
    logical :: verbose


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
    
    ! Parse grid settings
    file = trim(adjustl(root)) // "jupyter/grid_settings.dat"
    open(newunit=fid, file=trim(adjustl(file)))
    ! Header
    read(fid, *)
    read(fid, *) ranges(:, 1), sampling(1)
    read(fid, *) ranges(:, 2), sampling(2)
    close(fid)

    ! Generate grid from this
    nr = product(sampling)
    allocate(grid(n_dim, nr))
    call construct_2d_grid(sampling, ranges, grid)

    ! Print grid. Looks consistent with python setup
    ! do ir = 1, nr
    !     write(*, *) grid(:, ir)
    ! enddo

    ! Parse weights
    allocate(weight(nr))
    file = trim(adjustl(root)) // "jupyter/weights.dat"
    open(newunit=fid, file=trim(adjustl(file)))
    do ir = 1, nr
        read(fid, *) weight(ir)
    enddo
    close(fid)

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
    call weighted_kmeans(grid, weight, centroids, n_iter=niter, verbose=verbose)

    ! Output final centroids
    file = trim(adjustl(root)) // "jupyter/final_centroids.dat"
    open(newunit=fid, file=trim(adjustl(file)), action='write')
    write(fid, *) n_centroids
    do ic = 1, n_centroids
        write(fid, *) centroids(:, ic)
    enddo
    close(fid)

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