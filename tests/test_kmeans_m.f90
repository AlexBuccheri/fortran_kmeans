program test_kmeans_m
    use, intrinsic :: iso_fortran_env, only: dp => real64
    
    ! TODO Overload serial and mpi versions
    ! Test framework
    use fortuno_serial, only : execute_serial_cmd_app, is_equal, & 
       test => serial_case_item,&
       check => serial_check
    use fortuno_mpi, only : as_char, global_comm  !, is_equal, 
    !    test => mpi_case_item,&
    !    & check => mpi_check, test_item, this_rank

    ! Required modules
    use mpi_m,   only: mpi_t   
    use grids_m, only: linspace, linspace_to_grid
    use maths_m, only: all_close

    ! Module under test
    use kmeans_m, only : assign_points_to_centroids, update_centroids, compute_grid_difference, & 
                         report_differences_in_grids
    implicit none

    ! Register tests
    call execute_serial_cmd_app(testitems=[&
            test("Assign points to centroids", test_assign_points_to_centroids),     &
            test("Update centroids - no change", test_update_centroids_no_movement), &
            test("Difference in two sets of points", test_points_are_converged)      &
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

        call check(all(clusters(1:cluster_sizes(1), 1) == [1, 2, 3, 6, 7, 8, 11, 12, 13, 16, 17, 18, 21, 22, 23]), &
            msg='indices of points assigned to cluster 1')

        call check(all(clusters(1:cluster_sizes(2), 2) == [4, 5, 9, 10, 14, 15, 19, 20, 24, 25]), &
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


    subroutine test_update_centroids_no_movement()
        type(mpi_t) :: serial_comm
        real(dp), allocatable :: grid(:, :), centroids(:, :), uniform_weighting(:)
        integer,  allocatable :: clusters(:, :), cluster_sizes(:)

        integer :: i, Nr

        ! Mocked comm
        serial_comm = mpi_t(1, 1, 1)

        Nr = 8
        grid = reshape([[0, 0], [0, 1], [1, 0], [1, 1], &
                        [2, 0], [3, 0], [2, 1], [3, 1]], [2, Nr])

        ! Place centroids such that the grid points are evenly distributed 
        ! between the two of them
        centroids = reshape([[0.5, 0.5], &
                             [2.5, 0.5]], [2, 2])

        call assign_points_to_centroids(serial_comm, grid, centroids, clusters, cluster_sizes)

        call check(all(cluster_sizes == [4, 4]), msg='Each cluster gets half of the grid points')
        call check(all(clusters(:, 1) == [1, 2, 3, 4]), msg='Each cluster gets half of the grid points')
        call check(all(clusters(:, 2) == [5, 6, 7, 8]), msg='Each cluster gets half of the grid points')

        ! With a uniform weighting, the centroids should not change
        allocate(uniform_weighting(Nr))
        do i = 1, Nr
            uniform_weighting(i) = 1._dp / real(Nr, dp)
        enddo
        
        call update_centroids(grid, uniform_weighting, clusters, cluster_sizes, centroids)
        call check(all_close(centroids(:, 1), [0.5_dp, 0.5_dp]))
        call check(all_close(centroids(:, 2), [2.5_dp, 0.5_dp]))

    end subroutine test_update_centroids_no_movement


    ! TODO(Alex) Add 1-2 more tests for updating centroid positions


    subroutine test_points_are_converged()
        real(dp) :: x(10), y(10)
        real(dp), allocatable :: grid(:, :), second_grid(:, :)
        logical, allocatable :: points_differ(:)
        integer  :: ir, nr, n_diff
        real(dp) :: tol, noise, random_num
        logical  :: no_points_differ

        nr = 100
        tol = 1.e-6_dp

        allocate(grid(2, nr), second_grid(2, nr))
        call linspace(1._dp, 10.0_dp, size(x), x)
        call linspace(1._dp, 10.0_dp, size(y), y)
        call linspace_to_grid(x, y, grid)
        allocate(points_differ(nr))

        ! All points within tolerance
        do ir = 1, nr
            ! random_num in [0, 1]
            call random_number(random_num)
            noise = tol * random_num
            second_grid(:, ir) = grid(:, ir) + noise
        enddo

        call compute_grid_difference(grid, second_grid, tol, points_differ)
        n_diff = count(points_differ .eqv. .true.)
        call check(n_diff == 0, "No points should differ")
        no_points_differ = all(points_differ .eqv. .false.)
        call check(no_points_differ, "All points should be equivalent")

        ! One point not within tolerance
        second_grid = grid
        second_grid(:, 1) = grid(:, 1) + 1.1_dp * tol
        call compute_grid_difference(grid, second_grid, tol, points_differ)
        n_diff = count(points_differ .eqv. .true.)
        call check(n_diff == 1, "Only one point should differ")
        
        if (any(points_differ)) then
            call report_differences_in_grids(grid, second_grid, tol, points_differ)
        endif

    end subroutine test_points_are_converged


end program test_kmeans_m
