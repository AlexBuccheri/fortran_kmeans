!> @brief Module implementing k-means clustering algorithms
module kmeans_m
    use, intrinsic :: iso_fortran_env
    use omp_lib

    use mpi_m, only: mpi_t
    implicit none

    private
    public :: assign_points_to_centroids, update_centroids, compute_grid_difference, &
              report_differences_in_grids, weighted_kmeans

contains

    !> @brief Assign each grid point to the closest centroid. 
    !! A centroid and its set of nearest grid points defines a cluster.
    !!
    !! TODOs
    !!  * Add maths
    !!  * Look at replacing norm2 call - needs to be efficient to operate on n_points
    !!  * Extend to OMP? NOTE: Inner loop DOES NOT lend itself to OMP if the region is instantiated before outer loop
    !!
    subroutine assign_points_to_centroids(comm, grid_points, centroids, clusters, cluster_sizes)
        type(mpi_t),  intent(in) :: comm                            !< MPI instance
        real(real64), intent(in) :: grid_points(:, :)               !< Real-space grid (n_dims, N)
        real(real64), intent(in) :: centroids(:, :)                 !< Centroid positions (n_dims, Ncentroids)
        integer,      intent(out), allocatable :: clusters(:, :)    !< Cluster assignment for each grid point (Ncentroids, max_cpoints)
        integer,      intent(out), allocatable :: cluster_sizes(:)  !< Keep track of the number of points assigned to each cluster

        integer :: ir, ic, icen
        integer :: n_dims                                  !< System dimensions
        integer :: n_points                                !< Number of grid points
        integer :: n_centroids                             !< Number of  centroids
        integer :: max_cpoints                             !< Upper bound for number of points per centroid

        real(real64), allocatable :: displacements(:, :)   !< Displacements between centroids and grid point i
        real(real64), allocatable :: dmatrix(:)            !< Distance matrix
        integer,      allocatable :: work_clusters(:, :)   !< Work array for cluster assignment for each grid point (Ncentroids, upper_bound)

        n_dims = size(grid_points, 1)
        n_centroids = size(centroids, 2)
        n_points = size(grid_points, 2)

        ! We don't know how many points will be assigned per cluster
        ! Assuming a uniform distribution, one expects ~ N / Ncentroids
        ! Allocate a conservative upper bound (could replace with linked list)
        max_cpoints = int(3 * n_points / n_centroids)
        allocate(cluster_sizes(n_centroids), source=0)

        ! Work arrays
        allocate(work_clusters(max_cpoints, n_centroids))
        allocate(displacements(n_dims, n_centroids))
        allocate(dmatrix(n_centroids))

        do ir = 1, n_points
            do ic = 1, n_centroids
                displacements(:, ic) = centroids(:, ic) - grid_points(:, ir)
            enddo
            dmatrix = norm2(displacements, dim=1)
            icen = minloc(dmatrix, dim=1)
            ! Track increasing size of each cluster (also acts as an index)
            cluster_sizes(icen) = cluster_sizes(icen) + 1
            ! Assign point ir to the closest centroid
            work_clusters(cluster_sizes(icen), icen) = ir
        enddo 
        
        ! Resize. Note, there will still be redundant memory, hence we require cluster_sizes
        max_cpoints = maxval(cluster_sizes)
        allocate(clusters(max_cpoints, n_centroids), source=work_clusters(1: max_cpoints, 1: n_centroids))

    end subroutine assign_points_to_centroids
    

    !> @bried Compute a new set of centroids.
    subroutine update_centroids(grid, weight, clusters, cluster_sizes, centroids)
        real(real64), intent(in)    :: grid(:, :)               !< Real-space grid (n_dims, N)
        real(real64), intent(in)    :: weight(:)                !< Weights (N)
        integer,      intent(in)    :: clusters(:, :)           !< Cluster assignment for each grid point (max_cpoints, Ncentroids)
        integer,      intent(in)    :: cluster_sizes(:)         !< Keep track of the number of points assigned to each cluster
        
        real(real64), intent(inout) :: centroids(:, :)          !< In: centroid positions (n_dims, Ncentroids)
        !                                                         Out: Updated centroid positions
        
        integer                   :: n_dims, n_centroids, icen, j, ir
        real(real64), allocatable :: numerator(:)
        real(real64)              :: denominator

        n_dims = size(grid, 1)
        n_centroids = size(cluster_sizes)
        allocate(numerator(n_dims))

        ! TODO Make OMP
        do icen = 1, n_centroids
            numerator = 0._real64
            denominator = 0._real64
            ! Iterate over all points in current centroid
            do j = 1, cluster_sizes(icen)
                ir = clusters(j, icen)
                numerator(:) = numerator(:) + (grid(:, ir) * weight(ir))
                denominator = denominator + weight(ir)
            enddo
            centroids(:, icen) = numerator(:) / denominator
        enddo

    end subroutine update_centroids


    ! TODO(Alex) Document exactly what this is computing/returning
    subroutine compute_grid_difference(points, updated_points, tol, points_differ)
        real(real64), intent(in)  :: points(:, :)      !< Real-space grids (n_dims, N)
        real(real64), intent(in)  :: updated_points(:, :)
        real(real64), intent(in)  :: tol
        logical,      intent(out) :: points_differ(:)     !< |a_i - b_i|
       
        integer      ::  i, n_dim
        real(real64), allocatable :: diff(:)

        n_dim = size(points, 1)
        allocate(diff(n_dim))

        !$omp parallel do simd default(shared) private(i, diff)
        do i = 1, size(points, 2)
            points_differ(i) = .false.
            diff(:) = abs(updated_points(:, i) - points(:, i))
            points_differ(i) = any(diff > tol)
        enddo
        !$omp end parallel do simd

    end subroutine compute_grid_difference
  

    ! This routine can be less efficient as it should not be called in production runs
    ! or perhaps more accurately, when things are working
    subroutine report_differences_in_grids(points, updated_points, tol, points_differ)
        real(real64), intent(in)  :: points(:, :)      !< Real-space grids (n_dims, N)
        real(real64), intent(in)  :: updated_points(:, :)
        real(real64), intent(in)  :: tol
        logical,      intent(in)  :: points_differ(:)   !< If any element of point i > tol
        
        integer, allocatable :: indices(:)
        integer :: i, j, n_unconverged
        real(real64), allocatable :: diff(:)
        character(len=100) :: fmt

        indices = pack([(i, i=1,size(points_differ))], points_differ)
        n_unconverged = size(indices)
        allocate(diff(size(points, 1)))

        ! TODO Generalise for the pretty print... although it's absolutely not worth the hassle
        write(*, *) "# Current Point  ,  Prior Point  ,  |ri - r_{i-1}|  ,  tol"
        do j = 1, n_unconverged
            i = indices(j)
            diff(:) = abs(updated_points(:, i) - points(:, i))
            write(*, '(2(F16.10, X), 2(F16.10, X), 2(F16.10, X), F16.10)') updated_points(:, i), points(:, i), diff, tol
        enddo
        write(*, *) "Summary:", n_unconverged, "of out", size(points, 2), "are not converged"
   
    end subroutine report_differences_in_grids


    !> @brief Weighted K-means clustering.
    !!
    !! TODO(Alex) Add maths
    subroutine weighted_kmeans(comm, grid, weight, centroids, n_iter, centroid_tol, verbose)
        type(mpi_t),  intent(in) :: comm                         !< MPI instance
        real(real64), intent(in) :: grid(:, :)                   !< Grid    (n_dim, n_points)
        real(real64), intent(in) :: weight(:)                    !< Weights (n_points)
        real(real64), intent(inout) :: centroids(:, :)           !< In: Initial centroids (n_dim, n_centroid)
        !                                                           Out: Final centroids
        integer,      optional, intent(in) :: n_iter             !< Optional max number of iterations
        real(real64), optional, intent(in) :: centroid_tol       !< Optional convergence criterion
        logical,      optional, intent(in) :: verbose            !< Verbosity
  
        logical          :: print_out
        integer          :: n_iterations, nr, n_dim, n_centroid, i
        real(real64)     :: tol

        integer,      allocatable :: clusters(:, :), cluster_sizes(:)
        real(real64), allocatable :: prior_centroids(:, :)
        logical,      allocatable :: points_differ(:)

        ! Optional arg assignment
        n_iterations = 200
        if (present(n_iter)) n_iterations = n_iter

        tol = 1.e-6_real64
        if (present(centroid_tol)) tol = centroid_tol

        print_out = .false.
        if (present(verbose)) print_out = verbose

        nr = size(grid, 2)
        n_dim = size(grid, 1)
        n_centroid = size(centroids, 2)

        ! Consistency checks
        if (n_iterations < 1) then
            write(*, *) 'n_iter must be > 0'
            error stop 101
        endif
        if (size(weight) /= nr) then
            write(*, *) 'Number of weights inconsistent with number of grid points'
            error stop 102
        endif
        if (size(centroids, 1) /= n_dim) then
            write(*, *) 'Centroids 1st dim inconsistent grid points 1st dim'
            error stop 103
        endif

        ! Work arrays
        allocate(prior_centroids(n_dim, size(centroids, 2)), source=centroids)
        allocate(cluster_sizes(n_centroid))
        allocate(clusters(n_dim, n_centroid))
        allocate(points_differ(nr))

        do i = 1, n_iterations
            if (print_out) write(*, *) 'Iteration ', i
            call assign_points_to_centroids(comm, grid, centroids, clusters, cluster_sizes)
            call update_centroids(grid, weight, clusters, cluster_sizes, centroids)
            call compute_grid_difference(prior_centroids, centroids, tol, points_differ)
            if (any(points_differ)) then
                if (print_out) call report_differences_in_grids(prior_centroids, centroids, tol, points_differ)
                prior_centroids = centroids
            else
                ! Break loop - effectively exit routine
                exit 
            endif
        enddo

        ! TODO Return number of iterations used

    end subroutine weighted_kmeans

end module kmeans_m
