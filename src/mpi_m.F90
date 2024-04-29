!> MPI wrappers
module mpi_m
#ifdef USE_MPI    
    use mpi
#endif    
    implicit none
    ! Let things propagate
    public


    !> @brief Place-holder for MPI type
    type, public :: mpi_t
        integer :: comm  !< Communicator
        integer :: rank  !< Process/rank
        integer :: np    !< number of processes for this communicator
        integer :: ierr  !< MPI error
        integer :: root  !< Root process id
    contains
        procedure :: is_root 
    end type mpi_t

   ! Overload default constructor
    interface mpi_t
        procedure mpi_t_constructor
    end interface mpi_t

    ! interface mpi_bcast
    !     procedure mpi_bcast_char
    ! end interface mpi_bcast


! Expose for MPI or mock for serial
#ifdef USE_MPI 
    public :: MPI_COMM_WORLD
#else   
    integer, parameter :: mpi_comm_world = 100
#endif   


contains

! Serial overloads
#ifndef USE_MPI    
    !> @brief Serial overload for init
    subroutine mpi_init(ierr)
        integer, intent(inout) :: ierr
        ierr = 0
    end subroutine mpi_init

    !> @brief Serial overload for init
    subroutine MPI_Finalize(ierr)
        integer, intent(inout) :: ierr
        ierr = 0
    end subroutine MPI_Finalize
#endif   

    !> @brief Constructor for mpi_t
    function mpi_t_constructor(comm, root) result(mpi_instance)
        integer, intent(in), optional :: comm  !< Communicator
        integer, intent(in), optional :: root  !< Root process id
        class(mpi_t), allocatable :: mpi_instance
        
        allocate(mpi_instance)
        
        if (present(comm)) then
            mpi_instance%comm = comm
        else
            mpi_instance%comm = mpi_comm_world
        endif

#ifdef USE_MPI
        call MPI_Comm_rank(mpi_instance%comm, mpi_instance%rank, mpi_instance%ierr)
        call MPI_Comm_size(mpi_instance%comm, mpi_instance%np, mpi_instance%ierr)
#else
        mpi_instance%rank = 0
        mpi_instance%np = 1
        mpi_instance%ierr = 0
#endif      

        if (present(root)) then
            mpi_instance%root = root
        else
            mpi_instance%root = 0
        endif

    end function mpi_t_constructor


    !> @brief Get whether this process is root
    logical function is_root(this)
        class(mpi_t), intent(in) :: this
        is_root = this%rank == this%root
    end function is_root  


!     subroutine mpi_bcast_char(comm, x)
!         type(mpi_t), intent(in) :: comm  !< Communicator
!         character(len=*), intent(in) :: x

! #ifdef USE_MPI
!         ! Should check it's size, and not len()
!         call mpi_bcast(x, size(x), MPI_CHAR, comm%root, comm%comm)
! #else
!         continue
! #endif    

!     end subroutine mpi_bcast_char


    !> @brief Number of elements on process 'rank' for an array of n elements,
    !! distributed over np processes.
    integer function distribute_elements(np, n, rank)
        integer, intent(in) :: n               !< Number of elements to distribute
        integer, intent(in) :: np              !< Number of processes to distribute across
        integer, intent(in) :: rank            !< Process to compute local number of elements
        integer, allocatable :: indices(:, :)  !< Start and stop indices for each each process
        indices = distribute_elements_start_stop(np, n)
        distribute_elements = indices(2, rank + 1) - indices(1, rank + 1) + 1
    end function distribute_elements


    !> @brief Distribute N elements over NP processes.
    !!
    !! Example 1: n = 12 and np = 3
    !! [1, 2, 3, 4 | 5, 6, 7, 8 | 9, 10, 11, 12]
    !!
    !! Example 2: n = 14 and np = 4
    !! If n is not integer-divisible by np, assign the remainders equally between the first n_remainder processes
    !! [1, 2, 3, 4 | 5, 6, 7, 8 |9, 10, 11 |12, 13, 14]
    !!
    function distribute_elements_start_stop(np, n) result(indices)
        integer, intent(in) :: n       !< Number of elements to distribute
        integer, intent(in) :: np      !< Number of processes to distribute across
        integer :: elements_per_process, n_remainder
        integer :: n_assigned_processes, n_assigned_elements
        integer :: iprocess, istart, istop
        integer :: indices(2, np)

        ! When load-balanced
        elements_per_process = int(n / np)
        n_remainder = mod(n, np)

        ! Assign processes with extra elements first
        iprocess = 0
        do while (n_remainder > 0)
            iprocess = iprocess + 1
            istart = 1 + (iprocess - 1)*(elements_per_process + 1)
            istop = iprocess * (elements_per_process + 1)
            indices(:, iprocess) = [istart, istop]
            n_remainder = n_remainder - 1
        end do

        n_assigned_processes = iprocess
        n_assigned_elements = (elements_per_process + 1) * n_assigned_processes

        do iprocess = n_assigned_processes + 1, np
            istart = 1 + ((iprocess - n_assigned_processes - 1) * elements_per_process) + n_assigned_elements
            istop = ((iprocess - n_assigned_processes) * elements_per_process) + n_assigned_elements
            indices(:, iprocess) = [istart, istop]
        end do

    end function distribute_elements_start_stop

end module mpi_m
