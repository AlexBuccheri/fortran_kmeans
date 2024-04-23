!> MPI wrappers
module mpi_m
    use mpi
    private

    !> @brief Place-holder for MPI type
    type, public :: mpi_t
        integer :: comm  !< Communicator
        integer :: rank  !< Process/rank
        integer :: np    !< number of processes for this communicator
    end type mpi_t

end module mpi_m
