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

    public :: generate_real_space_grid !, generate_gaussian

contains

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


    ! subroutine generate_gaussian()
    ! end subroutine generate_gaussian


end module utils_m