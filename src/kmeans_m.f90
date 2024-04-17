!> @brief Module implementing k-means clustering algorithms
module kmeans_m
    use, intrinsic :: iso_fortran_env
    use utils_m, only: mpi_t
    implicit none

    private
    public :: assign_points_to_centroids

contains

    !> @brief Assign each grid point to the closest centroid. 
    !! A centroid and its set of nearest grid points defines a cluster.
    !!
    !! TODO Add maths
    !!      Distribute routine
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
        allocate(work_clusters(n_centroids, max_cpoints))
        allocate(displacements(n_dims, n_centroids))
        allocate(dmatrix(n_centroids))

        ! TODO Distribute outer loop
        ! NOTE: Inner loop DOES NOT lend itself to OMP if the region is instantiated here
        do ir = 1, n_points
            do ic = 1, n_centroids
                displacements(:, ic) = centroids(:, ic) - grid_points(:, ir)
            enddo
            ! TODO Look at replacing - needs to be efficient to operate on n_points
            dmatrix = norm2(displacements, dim=1)
            icen = minloc(dmatrix, dim=1)
            ! Track increasing size of each cluster (also acts as an index)
            cluster_sizes(icen) = cluster_sizes(icen) + 1
            ! Assign point ir to the closest centroid
            work_clusters(icen, cluster_sizes(icen)) = ir
        enddo 
        
        ! Resize. Note, there will still be redundant memory, hence we require cluster_sizes
        max_cpoints = maxval(cluster_sizes)
        allocate(clusters(n_centroids, max_cpoints), source=work_clusters(1: n_centroids, 1: max_cpoints))

    end subroutine assign_points_to_centroids
    
    ! subroutine update_centroids()
    ! end subroutine update_centroids

    ! function points_are_converged() result(converged)
    !     logical :: converged
    ! end function points_are_converged

    ! !> @brief Weighted K-means clustering.
    ! subroutine weighted_kmeans()
    ! end subroutine weighted_kmeans

end module kmeans_m
