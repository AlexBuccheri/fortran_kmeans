program test_kmeans_m
    use, intrinsic :: iso_fortran_env
    ! Test framework
    use fortuno_serial, only : execute_serial_cmd_app, is_equal, & 
       test => serial_case_item,&
       check => serial_check
    ! Required modules
    use utils_m, only: mpi_t, generate_real_space_grid
    ! Module under test
    use kmeans_m, only : assign_points_to_centroids
    implicit none

    ! Register tests by providing name and subroutine to run for each test.
    ! Note: this routine does not return but stops the program with the right exit code.
    call execute_serial_cmd_app(testitems=[&
            test("Assign points to centroids", test_assign_points_to_centroids) &
        ])
  
contains

    subroutine test_assign_points_to_centroids
        type(mpi_t) :: serial_comm

        ! Grid
        integer, parameter :: n_dims = 2
        integer :: n_total

        integer :: sampling(n_dims)
        real(real64) :: spacings(n_dims, n_dims)
        real(real64), allocatable :: grid(:, :), centroids(:, :)

        ! Clusters, centred on centroids
        integer :: n_centroids
        integer,     allocatable  :: clusters(:, :), cluster_sizes(:)

        ! Inputs
        serial_comm = mpi_t(1, 1, 1)
        ! 2D Grid
        sampling = [5, 5]
        spacings = reshape([0.1_real64, 0.00_real64, &
                            0.00_real64, 0.1_real64], [n_dims, n_dims])
        n_total = product(sampling)
        allocate(grid(n_dims, n_total))            
        call generate_real_space_grid(sampling, spacings, grid)

        ! TODO(Alex) Reproduce the python tests I wrote
        n_centroids = 2
        allocate(centroids(n_dims, n_centroids))
        call assign_points_to_centroids(serial_comm, grid, centroids, clusters, cluster_sizes)

    end subroutine test_assign_points_to_centroids


end program test_kmeans_m
