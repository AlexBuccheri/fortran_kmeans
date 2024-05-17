!> Grid routines
module grids_m
    use, intrinsic :: iso_fortran_env
    use omp_lib
    implicit none
    private

    ! Exposed routines
    public :: generate_real_space_grid, linspace, linspace_to_grid, generate_gaussian, &
            &  discretise_values_to_grid

    interface linspace_to_grid
        module procedure :: linspace_to_grid2d, linspace_to_grid3d, linspace_to_grid_suboptimal
    end interface linspace_to_grid

contains


    ! TODO(Alex) Does not make the most sense to provide both sampling and spacing => One has to work
    ! out what their box limits are
    !> @brief Generate a real-space grid.
    !! 
    !! Produce 1D, 2D and 3D grids.
    !! 
    !!  !One could also implement as integer grid construction, then convert to real spacing:
    !!  ir = 0
    !!  do k = 1, sampling(3)
    !!      do j = 1, sampling(2)
    !!          do i = 1, sampling(1)
    !!              ir = ir + 1
    !!              integer_grid(:, ir) = [i, j, k]
    !!          enddo   
    !!      enddo
    !!  enddo 
    !!  call dgemm('n', 'n', n_dims, n_points, n_dims, &
    !!     1._real64, spacings, n_dims, integer_grid, n_dims, 0._real64, grid, n_dims)
    subroutine generate_real_space_grid(sampling, spacings, grid, origin)
        integer,                intent(in)  :: sampling(:)
        real(real64),           intent(in)  :: spacings(:, :)
        real(real64), optional, intent(in)  :: origin(:)
        real(real64),           intent(out) :: grid(:, :)  !< real-space grid

        integer :: i, j, k, jk, ir, n_points, n_dims
        integer, allocatable :: int_point(:)
        real(real64), allocatable :: org(:)

        n_dims = size(sampling)
        n_points = product(sampling)
        allocate(int_point(n_dims))

        if (.not. present(origin)) then
            allocate(org(n_dims), source=0._real64)
        else
            allocate(org(n_dims), source=origin)
        endif

        if (n_dims == 1) then
            !$omp parallel do simd default(shared) private(i)
            do i = 1, sampling(1)
                grid(:, i) = real(i, kind=real64) * spacings(1, 1) - org
            enddo
            !$omp end parallel do simd

        else if(n_dims == 2) then
            !$omp parallel do simd collapse(2) default(shared) private(i, j, ir)
            do j = 1, sampling(2)
                do i = 1, sampling(1)
                    ir = i + (j - 1) * sampling(1)
                    grid(:, ir) = matmul(spacings, [i, j]) - org
                enddo   
            enddo
            !$omp end parallel do simd

        else if(n_dims == 3) then
            !$omp parallel do simd collapse(3) default(shared) private(i, j, k, ir)
            do k = 1, sampling(3)
                do j = 1, sampling(2)
                    do i = 1, sampling(1)
                        ! Have to evaluate it here to allow loop collapsing
                        jk = j + (k - 1) * sampling(2)
                        ir = i + (jk - 1) * sampling(1)
                        grid(:, ir) = matmul(spacings, [i, j, k]) - org
                    enddo   
                enddo
            enddo 
            !$omp end parallel do simd

        else
            write(*, *) 'Grid construction with dimensions', n_dims, "is not implemented"
            stop
        endif

    end subroutine generate_real_space_grid


    !> @ brief Return evenly spaced numbers over a specified interval.
    !! Returns num evenly spaced samples, calculated over the interval [start, stop].
    !!
    !! Spacing not returned, but is trivially grid(2) - grid(1)
    !! 
    !! Based on https://numpy.org/doc/stable/reference/generated/numpy.linspace.html
    subroutine linspace(start, stop, n_points, grid, end_point)
        real(real64),      intent(in)  :: start
        real(real64),      intent(in)  :: stop
        integer,           intent(in)  :: n_points
        logical, optional, intent(in)  :: end_point      !< The endpoint of the interval can optionally be excluded.
        real(real64),      intent(out) :: grid(n_points) !< real-space grid

        integer :: i
        real(real64) :: dx
        logical :: include_end_point

        include_end_point = .true.
        if (present(end_point)) include_end_point = end_point

        if (include_end_point) then
            dx = (stop - start) / real(n_points - 1, real64)
        else
            dx = (stop - start) / real(n_points, real64)
        endif

        !$omp parallel do simd default(shared) private(i)
        do i = 1, n_points
            grid(i) = start + (i - 1) * dx
        enddo
        !$omp end parallel do simd

    end subroutine linspace


    !> @brief Create an ND grid from linear-spaced data
    !!
    !! Data is passed as a contiguous 1D array, so one also
    !! needs to pass the number of x, y, z, ... points for
    !! the routine to distinguish dimensions. 
    !! THIS IS BAD - requires discontiguous access of the `ranges` or `grid` array, 
    !! which for large grids will create numerous cache misses.
    !! Unless I return grid(n_points, n_dim)... which is usually not desirable.
    !!
    !! ```fortran
    !! call linspace(1._dp, 10.0_dp, 10, x)
    !! call linspace(1._dp, 20.0_dp, 10, y)
    !! call linspace_to_grid3D([x, y], grid_2d)
    !! ```
    subroutine linspace_to_grid_suboptimal(ranges, n_points, grid)
        real(real64), contiguous, intent(in)  :: ranges(:)
        integer,                  intent(in)  :: n_points(:)
        real(real64),             intent(out) :: grid(:, :)  !< real-space grid

        integer :: idim, ix, ir

        ! assert size(grid, 1) == product(n_points)
        ir = 0
        do idim = 1, size(n_points)
            do ix = 1, n_points(idim)
                ir = ir + 1
                grid(ir, idim) = ranges(ix)
            enddo
        enddo

    end subroutine linspace_to_grid_suboptimal


    !> @brief Create a 2D grid from linear-spaced data
    subroutine linspace_to_grid2d(x, y, grid)
        real(real64), contiguous, intent(in)  :: x(:), y(:)
        real(real64),             intent(out) :: grid(:, :)  !< real-space grid
        integer :: nx, ny, i, j, ir

        nx = size(x)
        ny = size(y)
        ! assert size(grid, 2) == nx * ny

        !$omp parallel do simd collapse(2) default(shared) private(i, j, ir)
        do j = 1, ny
            do i = 1, nx
                ir = i + (j - 1) * nx
                ! y(j) accessed per inner loop so one can collapse the loops with OMP
                grid(:, ir) = [x(i), y(j)]
            enddo
        enddo
        !$omp end parallel do simd

    end subroutine linspace_to_grid2d


    !> @brief Create a 3D grid from linear-spaced data
    subroutine linspace_to_grid3d(x, y, z, grid)
        real(real64), contiguous, intent(in)  :: x(:), y(:), z(:)
        real(real64),             intent(out) :: grid(:, :)  !< real-space grid

        integer :: nx, ny, nz, i, j, k, ir, jk

        nx = size(x)
        ny = size(y)
        nz = size(z)
        ! assert size(grid, 2) == nx * ny * nz

        !$omp parallel do simd collapse(3) default(shared) private(i, j, k, ir)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    ! jk, y(j), z(k) accessed per inner loop so one can collapse the loops with OMP
                    jk = j + (k - 1) * ny
                    ir = i + (jk - 1) * nx
                    grid(:, ir) = [x(i), y(j), z(k)]
                enddo
            enddo
        enddo
        !$omp end parallel do simd

    end subroutine linspace_to_grid3d


    !> @brief Generate a Gaussian on a grid
    !!
    !! TODO Give maths
    subroutine generate_gaussian(mean, sigma, grid, gaussian)
        real(real64), intent(in)  :: mean(:)        !< Centre of the Gaussian
        real(real64), intent(in)  :: sigma(:)       !< Width of the Gaussian
        real(real64), intent(in)  :: grid(:, :)     !< Grid on which to evaluate the Gaussian
        real(real64), intent(out) :: gaussian(:)  !< Gaussian

        real(real64), parameter :: pi = 3.14159265359_real64
        integer :: ir
        real(real64) :: n, n_dim, exponent, norm

        n_dim = real(size(mean), real64)
        n =  0.5 * n_dim
        norm = 1._real64  / (product(sigma) * (2._real64 * pi)**n)

        !$omp parallel do simd default(shared) private(ir, exponent)
        do ir = 1, size(grid, 2)
            ! (r -r0) / sigma
            exponent = sum((grid(:, ir) - mean(:))**2._real64 / sigma(:)**2._real64)
            gaussian(ir) =  norm * exp(-0.5_real64 * exponent)
        enddo
       !$omp end parallel do simd

    end subroutine generate_gaussian


    ! Serial reference
    subroutine discretise_values_to_grid(values, grid)
        real(real64), intent(inout) :: values(:, :)  !< In: Continuous values  (n_dim, M) 
        !                                               Out: Discretised values (n_dim, M) 
        real(real64), intent(in) :: grid(:, :)        !< Discrete grid (n_dim, Nr)
        
        integer :: ir, iv, nr, iopt, n_dim
        real(real64) :: norm, new_norm
        real(real64), allocatable :: val(:)

        n_dim = size(grid, 1)
        nr = size(grid, 2)
        allocate(val(n_dim))

        do iv = 1, size(values, 2)
            val = values(:, iv)
            ! Initialise 
            norm = sum((grid(:, 1) - val(:))**2)
            iopt = -1
            do ir = 2, nr
                new_norm = sum((grid(:, 1) - val(:))**2)
                if (new_norm < norm) then
                    norm = new_norm
                    iopt = ir
                endif
            enddo
            values(:, iopt) = grid(:, iopt)
        enddo

    end subroutine discretise_values_to_grid


    ! For parallel, one should work with the integer representation of the grid, not the real.
    ! Better than suggestion below. Simpler than k-d tree search.

    ! !> @brief Given a set of values that span the range of a grid, 
    ! !! assign each value to the nearest, discrete grid point.

    ! !! Proper approach: Parallelised k-d tree search, but I'm not going to partition
    ! !! the grid up according to this. Particularly as this ultimately needs to work with metis
    ! !! (not clear that my idea will)
    ! !!
    ! !! Outline of my idea:
    ! !! a) Approximate or exactly compute min_norm per process
    ! !!   If approximate, one does so on a random sampling of points
    ! !!   Store mni_norm(1:Nv)
    ! !! b) Return process of min_norm for each iv
    ! !!   So what one wants is a subset of the global iv indices associated with each process
    ! !! c) Compute the min_norm for that process and take it to be the global minimum
    ! !!    - Store associated grid points
    ! !! d) allgatherv the grid points for discretised values

    ! subroutine discretise_values_to_grid_mpi(comm, values, grid)
    !     real(real64), intent(inout), :: values(:, :)  !< In: Continuous values  (n_dim, M) 
    !     !                                               Out: Discretised values (n_dim, M) 
    !     real(real64), intent(in) :: grid(:, :)        !< Discrete grid (n_dim, Nr)

    !     integer :: n_dim, nr, ir, iv
    !     real(real64) :: tmp_norm
    !     real(real64), allocatable :: val(:)

    !     n_dim = size(grid, 1)
    !     nr = size(grid, 2)
    !     allocate(val(n_dim))

    !     ! a) Compute min_norm for each value, using either all of the subgrid, 
    !     ! or random sampling of the subgrid 
    !     ! Could colapse this with OMP
    !     do iv = 1, size(values, 2)
    !         ! val = values(:, iv)
    !         do ir = 1, nr
    !             min_norm(i) = norm2(grid(:, ir) - values(:, i))
    !         enddo
    !     enddo

    !     ! b) For each process, store the iv of min_norm(iv) that are minimised over all processes
    ! NOT FINISHED

    ! end subroutine discretise_values_to_grid_mpi

end module grids_m
