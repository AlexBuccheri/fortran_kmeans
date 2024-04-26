program test_kmeans_m
    use, intrinsic :: iso_fortran_env, only: dp => real64
    
    ! Get via wrapper
    ! use fortuno_serial, only : execute_serial_cmd_app, is_equal, & 
    !    test => serial_case_item,&
    !    check => serial_check

    use fortuno_interface_m, only: execute_cmd_app, test, check, is_equal
    use grids_m,             only: linspace, linspace_to_grid, generate_gaussian
    use maths_m,             only: all_close

    ! Module under test
    use kmeans_m, only : assign_points_to_centroids, update_centroids, compute_grid_difference, & 
                         report_differences_in_grids, weighted_kmeans
    implicit none

    ! Register tests
    call execute_cmd_app(testitems=[&
            test("Assign points to centroids", test_assign_points_to_centroids),     &
            test("Update centroids - no change", test_update_centroids_no_movement), &
            test("Difference in two sets of points", test_points_are_converged),     &  
            test("Run weighted k-means", test_weighted_kmeans)     & 
        ])

contains

    subroutine test_assign_points_to_centroids
        ! Grid
        integer, parameter    :: n_dims = 2
        integer               :: nx, ny
        real(dp), allocatable :: x(:), y(:), grid(:, :)
        ! Centroids
        real(dp), allocatable :: centroids(:, :)

        ! Clusters, centred on centroids
        integer                :: n_centroids
        integer,  allocatable  :: clusters(:, :), cluster_sizes(:)

        integer :: ic, ir, ig, ierr

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
        call assign_points_to_centroids(grid, centroids, clusters, cluster_sizes)

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
        real(dp), allocatable :: grid(:, :), centroids(:, :), uniform_weighting(:)
        integer,  allocatable :: clusters(:, :), cluster_sizes(:)

        integer :: i, Nr

        Nr = 8
        grid = reshape([[0, 0], [0, 1], [1, 0], [1, 1], &
                        [2, 0], [3, 0], [2, 1], [3, 1]], [2, Nr])

        ! Place centroids such that the grid points are evenly distributed 
        ! between the two of them
        centroids = reshape([[0.5, 0.5], &
                             [2.5, 0.5]], [2, 2])

        call assign_points_to_centroids(grid, centroids, clusters, cluster_sizes)

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


    ! Setup routine
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


    ! Setup routine
    subroutine construct_gaussians_on_grid(mean, sigma, grid, total_gaussian)
        real(dp), intent(in)   :: mean(:, :)  !< mean, shape (n_dim, n_gauss)
        real(dp), intent(in)   :: sigma(:)    !< Same broadening for all Gaussians, shape (n_dim)
        real(dp), intent(in)   :: grid(:, :)  !< shape(n_dim, nr)

        integer                :: n_dim, nr, n_gauss, i
        real(dp), allocatable  :: gaussian(:)
        real(dp), intent(out)  :: total_gaussian(:)

        n_dim = size(grid, 1)
        nr = size(grid, 2)
        n_gauss = size(mean, 2)
        total_gaussian = 0._dp

        allocate(gaussian(nr))

        do i = 1, n_gauss
            call generate_gaussian(mean(:, i), sigma, grid, gaussian)
            total_gaussian = total_gaussian + gaussian
        enddo

    end subroutine construct_gaussians_on_grid


    ! TODO(Alex) Make this work as required
    !> Generate `n` random integers from the range [1, m]
    ! subroutine random_points(m, n, points, fixed_seed)
    !     integer, intent(in) :: m          !< Integer range of values [1: m]
    !     integer, intent(in) :: n          !< Number of random integers to generate
    !     integer, intent(out) :: points(n) !< Random integers
    !     integer, optional, intent(in) :: fixed_seed  !< Fix seed for reproducibility

    !     integer              :: i, rand_index, my_size
    !     real(dp)             :: rand_val
    !     integer, allocatable :: state(:)

    !     if (present(fixed_seed)) then
    !         call random_seed(size=my_size)
    !         allocate(state(my_size))
    !         call random_seed(put=state)
    !     endif

    !     ! TODO Should probably do some seed init for randomness
    !     ! call random_init(.true., .false.)

    !     do i = 1, n
    !         ! Generate a random real number [0, 1]
    !         call random_number(rand_val)
    !         ! Scale rand_val to range [1, m]
    !         rand_index = 1 + floor(rand_val * m)
    !         points(i) = rand_index
    !     end do

    ! end subroutine random_points


    subroutine test_weighted_kmeans()
        ! Grid
        integer,  parameter  :: n_dim = 2
        integer              :: nr
        integer              :: sampling(n_dim)
        real(dp)             :: ranges(n_dim, 2)
        real(dp), allocatable :: grid(:, :)
        ! Gaussian
        integer,  parameter  :: n_gauss = 4
        real(dp)  :: mean(n_dim, n_gauss), sigma(n_dim)
        real(dp), allocatable :: gaussian(:), total_gaussian(:)
        ! Centroids
        integer,  allocatable :: indices(:)
        real(dp), allocatable :: centroids(:, :)

        integer :: ix, iy, ir, n_iter, i, n_centroid, ierr
        logical :: verbose
        integer, parameter :: seed = 20180815

        ! Grid
        sampling = [50, 50]
        nr = product(sampling)
        ranges(:, 1) = [1._dp, 5.0_dp]
        ranges(:, 2) = [1._dp, 5.0_dp]
        call construct_2d_grid(sampling, ranges, grid)

        ! Gaussians on a grid
        sigma = [0.4_dp, 0.4_dp]
        mean = reshape([2._dp, 2._dp, &
                        2._dp, 4._dp, &
                        4._dp, 2._dp, &
                        4._dp, 4._dp], [n_dim, n_gauss])
        allocate(total_gaussian(nr))
        call construct_gaussians_on_grid(mean, sigma, grid, total_gaussian)

        ! Output Gaussians for gnuplot
        ir = 0
        do iy = 1 , sampling(2)
            do ix = 1, sampling(1)
                ir = ir + 1
                write(100, *) grid(:, ir), total_gaussian(ir)
            enddo
            write(100, *)
        enddo

        ! Optimal centroids
        n_centroid = 25
        n_iter = 50
        verbose = .true.
        allocate(centroids(n_dim, n_centroid))

        ! Pick random points, with a fixed seed such that it's reproducible... or just hard-code for now
        allocate(indices(n_centroid))
        indices = [1178, 294, 1637, 282, 2353, 986, 496, 998, 1563, 1522, 1013, 1566, 1534, 2339, &
                   679, 43, 559, 257, 2218, 321, 790, 2021, 738, 283, 798]

        ! Assign initial centroids
        do i = 1, n_centroid
            ir = indices(i)
            centroids(:, i) = grid(:, ir)
            ! Output initial centroids for plotting over Gaussians
            ! 3rd column is a dummy value for plotting
            write(101, *) grid(:, ir), 0._dp
        enddo
    
        call weighted_kmeans(grid, total_gaussian, centroids, n_iter, verbose=verbose)

        ! Output final centroids for plotting over Gaussians
        do i = 1, n_centroid
            write(102, *) centroids(:, i), 0._dp
        enddo

        ! GNUPLOTing:
        ! set term x11
        ! splot 'fort.100' u 1:2:3 w l title'Gaussian weighting'
        ! replot 'fort.101' u 1:2:3 w p lc 'red' title'Initial centroids'
        ! replot 'fort.102' u 1:2:3 w p lc blue title'Optimised centroids'

        ! TODOs:
        ! Assert on the number of iterations used
        ! Test https://docs.lfortran.org/en/installation/ for notebook 
        ! Assert final centroid positions -> Regression test rather than unit
        ! Check if they fall on grid points -> If not, pin them to the closest grid points
        ! Move on to a) Parallelisation and b) greedy k-means

    end subroutine test_weighted_kmeans


end program test_kmeans_m
