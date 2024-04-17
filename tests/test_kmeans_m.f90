program test_kmeans_m
    use, intrinsic :: iso_fortran_env, only: dp => real64
    ! Test framework
    use fortuno_serial, only : execute_serial_cmd_app, is_equal, & 
       test => serial_case_item,&
       check => serial_check

       use fortuno_mpi, only : as_char, global_comm  !, is_equal, 
    !    test => mpi_case_item,&
    !    & check => mpi_check, test_item, this_rank

    ! Required modules
    use utils_m, only: mpi_t, linspace, linspace_to_grid

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
        integer, parameter    :: n_dims = 2
        integer               :: nx, ny
        real(dp), allocatable :: x(:), y(:), grid(:, :)
        ! Centroids
        real(dp), allocatable :: centroids(:, :)

        ! Clusters, centred on centroids
        integer                :: n_centroids
        integer,  allocatable  :: clusters(:, :), cluster_sizes(:)

        integer :: ic, ir, ig

        ! Inputs
        serial_comm = mpi_t(1, 1, 1)

        ! 2D Grid
        nx = 5
        ny = 5
        allocate(x(nx), y(ny), grid(n_dims, nx * ny))
        call linspace(1._dp, 5.0_dp, nx, x)
        call linspace(1._dp, 5.0_dp, ny, y)
        call linspace_to_grid(x, y, grid)

        n_centroids = 2
        allocate(centroids(n_dims, n_centroids))
        centroids = reshape([2._dp, 3._dp, &  ! centroid 1
                             4._dp, 3._dp],&  ! centroid 2
                            [n_dims, n_centroids])
        call assign_points_to_centroids(serial_comm, grid, centroids, clusters, cluster_sizes)

        ! Central 5 points are equidistant from the centroids. In the case that >=2 values are the same, 
        ! minloc appears to return the first instance
        call check(cluster_sizes(1) == 15, msg='Expect all points on the left, plus the line along the &
                                                middle to assign to first cluster')
        call check(cluster_sizes(2) == 10, msg='Expect all points on the right to assign to first cluster')

        call check(all(clusters(1, 1:cluster_sizes(1)) == [1, 2, 3, 6, 7, 8, 11, 12, 13, 16, 17, 18, 21, 22, 23]), &
            msg='indices of points assigned to cluster 1')

        call check(all(clusters(2, 1:cluster_sizes(2)) == [4, 5, 9, 10, 14, 15, 19, 20, 24, 25]), &
            msg='indices of points assigned to cluster 2')

        ! Inspection
        ! do ic = 1, n_centroids
        !     do ig = 1, cluster_sizes(ic)
        !         ir = clusters(ic, ig)
        !         ! If ic ==1, then |c1 - ri| should be smaller
        !         ! If ic ==2, then |c2 - ri| should be smaller
        !         write(*, *) 'cluster, grid index, |c1 - ri|, |c2 - ri|', ic, ir, &
        !         norm2(centroids(:, 1) - grid(:, ir)), norm2(centroids(:, 2) - grid(:, ir))
        !     enddo
        ! enddo

    end subroutine test_assign_points_to_centroids


end program test_kmeans_m
