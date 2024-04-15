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

    !> @brief Generate a real-space grid
    !! 
    ! Produce 1D, 2D and 3D grids
    subroutine generate_real_space_grid(sampling, spacings, grid, origin)
        integer,                intent(in)  :: sampling(:)
        real(real64),           intent(in)  :: spacings(:, :)
        real(real64), optional, intent(in)  :: origin(:)
        real(real64),           intent(out) :: grid(:, :)  !< real-space grid

        integer :: i, j, k, ir, n_points, n_dims
        integer, allocatable :: integer_grid(:, :)

        n_dims = size(sampling)
        n_points = product(sampling)
        allocate(integer_grid(n_dims, n_points))

        if (n_dims == 1) then
            do i = 1, sampling(1)
                integer_grid(:, i) = [i]
            enddo
            ! Reduces to scaling of a vector by a constant
            grid = spacings(1, 1) * real(integer_grid, kind=real64)

        else if(n_dims == 2) then
            ir = 0
            do j = 1, sampling(2)
                do i = 1, sampling(1)
                    ir = ir + 1
                    integer_grid(:, ir) = [i, j]
                enddo   
            enddo            
             call dgemm('n', 'n', n_dims, n_points, n_dims, &
                       1._real64, spacings, n_dims, integer_grid, n_dims, 0._real64, grid, n_dims)
    
        else if(n_dims == 3) then
            ir = 0
            do k = 1, sampling(3)
                do j = 1, sampling(2)
                    do i = 1, sampling(1)
                        ir = ir + 1
                        integer_grid(:, ir) = [i, j, k]
                    enddo   
                enddo
            enddo 
            call dgemm('n', 'n', n_dims, n_points, n_dims, &
            1._real64, spacings, n_dims, integer_grid, n_dims, 0._real64, grid, n_dims)
    
        else
            write(*, *) 'Grid construction with dimensions', n_dims, "is not implemented"
            stop
        endif


        ! Decide best way to implement - probs an OMP loop
        ! if(present(origin)) then
        !     call daxpy (n, da, dx, incx, dy, incy)
        !     grid = grid - origin
        ! endif

    end subroutine generate_real_space_grid


    ! subroutine generate_gaussian()
    ! end subroutine generate_gaussian


end module utils_m