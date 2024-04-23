!> @brief Module implementing k-means clustering algorithms
module kmeans_m
    use, intrinsic :: iso_fortran_env
    use omp_lib

    use mpi_m, only: mpi_t
    implicit none

    private
    public :: assign_points_to_centroids, update_centroids, points_are_converged

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
            do j = 1, cluster_sizes(icen)
                ir = clusters(j, icen)
                numerator(:) = numerator(:) + (grid(:, ir) * weight(ir))
                denominator = denominator + weight(ir)
            enddo
            centroids(:, icen) = numerator(:) / denominator
        enddo

    end subroutine update_centroids


    !> @brief Given the difference in two sets of points, determine whether the updated
    !!        points are sufficiently close to the prior points.
    function points_are_converged(updated_points, points, tol, verbose) result(converged)
        real(real64), intent(in)    :: updated_points(:, :), points(:, :) !< Real-space grids (n_dims, N)
        real(real64), intent(in)    :: tol
        logical, intent(in), optional :: verbose
        logical :: converged

        logical :: print_out
        integer :: n_dims, n_points, n_unconverged, i
        real(real64), allocatable :: norm(:)   !< norm of each vector (a_i - b_i)

        n_dims = size(updated_points, 1)
        n_points = size(updated_points, 2)
        allocate(norm(n_points))

        print_out = .false.
        if (present(verbose)) then
            if (verbose) print_out = .true.
        endif

        !$omp parallel do simd default(shared) private(i)
        do i = 1, n_points
            norm(i) = norm2(updated_points(:, i) - points(:, i))
        enddo
        !$omp end parallel do simd

        if(print_out) then
            write(*, *) "# Current Point  ,  Prior Point  ,  |ri - r_{i-1}|  ,  tol"
            n_unconverged = 0
            do i = 1, n_points
                if (norm(i) > tol) then
                    write(*, *) updated_points(:, i), points(:, i), norm(i), tol
                    n_unconverged = n_unconverged + 1
                endif
            enddo
            write(*, *) "Summary:", n_unconverged, "of out", n_points, "are not converged"
            converged = n_unconverged == 0
        else
            converged = .not. any(norm > tol)
        endif

    end function points_are_converged

    ! !> @brief Weighted K-means clustering.
    ! subroutine weighted_kmeans()
    ! end subroutine weighted_kmeans

end module kmeans_m
