!> Utility routines
module utils_m
    use, intrinsic :: iso_fortran_env

    implicit none
    private

    !> @brief Place-holder for MPI type
    type, public :: mpi_t
        integer :: comm  !< Communicator
        integer :: rank  !< Process/rank
        integer :: np    !< number of processes for this communicator
    end type mpi_t

    public :: generate_real_space_grid, linspace !, generate_gaussian

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
        integer, allocatable :: int_point(:), integer_grid(:, :)
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


    ! !> @brief Create a 3D grid from linear-spaced data
    ! subroutine linspace_to_grid3D(x, y, z, n_points, grid, origin)
    !     real(real64), contiguous, intent(in)  :: x(:)
    !     real(real64), contiguous, intent(in)  :: y(:)
    !     real(real64), contiguous, intent(in)  :: z(:)

    !     real(real64),           intent(out) :: grid(:, :)  !< real-space grid

    !     integer                    :: n_dims, i
    !     real(real64), allocatable  :: spacings(:, :)
    !     real(real64), allocatable  :: origin(:)

    !     n_dims = size(limits, 2)
        
    !     do k = 1, 
    !     do j = 1, 
    !     do i = 1, 

         


    ! end subroutine linspace_to_grid3D

    ! subroutine generate_gaussian()
    ! end subroutine generate_gaussian


end module utils_m