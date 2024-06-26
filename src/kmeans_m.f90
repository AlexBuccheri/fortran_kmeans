!> @brief Module implementing k-means clustering algorithms
module kmeans_m
    use, intrinsic :: iso_fortran_env
    use omp_lib
#ifdef USE_MPI          
    use mpi, only: MPI_DOUBLE_PRECISION, MPI_SUM, MPI_ALLREDUCE, MPI_IN_PLACE
#endif 
    use mpi_m, only: mpi_t

    implicit none

    private
    public :: assign_points_to_centroids, update_centroids, compute_grid_difference, &
              report_differences_in_grids, weighted_kmeans

contains

    !> @brief Assign each grid point to the closest centroid. 
    !! A centroid and its set of nearest grid points defines a cluster.
    !!
    !! This can be mathematically expressed as:
    !!\f[ 
    !!  C_\mu=\left\{Z\left(\mathbf{r}_i\right) \mid \operatorname{dist}\left(Z\left(\mathbf{r}_i\right), Z\left(\mathbf{r}_\mu\right)\right) 
    !!      \leq \operatorname{dist}\left(Z\left(\mathbf{r}_i\right), Z\left(\mathbf{r}_m\right)\right) \text { for all } m \neq \mu\right\}
    !!\f]
    !! where \f$ Z(\mathbf{r}_i) \f$ are data points and \f$ Z(\mathbf{r}_\mu) \f$ is the centroid. The distance metric, \f$\operatorname{dist}\f$,
    !! is defined as the dot product in this implementation.    
    !!
    !! See eqs 10-13 in [Complex-valued K-means clustering of interpolative separable density fitting algorithm for large-scale hybrid functional 
    !! enabled ab initio  molecular dynamics simulations within plane waves](https://doi.org/10.48550/arXiv.2208.07731)
    subroutine assign_points_to_centroids(grid_points, centroids, ir_to_ic)
        real(real64), intent(in)  :: grid_points(:, :)     !< Real-space grid (n_dims, N)
        real(real64), intent(in)  :: centroids(:, :)       !< Centroid positions (n_dims, Ncentroids)
        integer,      intent(out) :: ir_to_ic(:)           !< Index array that maps grid indices to centroid indices

        integer :: ir, ic, icen, n_unassigned
        integer :: n_dims                                  !< System dimensions
        integer :: n_points                                !< Number of grid points
        integer :: n_centroids                             !< Number of centroids
        integer, allocatable :: npoints_per_cluster(:)
        real(real64) :: dist, min_dist

        n_dims = size(grid_points, 1)
        n_centroids = size(centroids, 2)
        n_points = size(grid_points, 2)

        if (size(centroids, 1) /= n_dims) then
            write(*, *) 'Size of centroid vector inconsistent with size of grid point vector'
            error stop 101
        endif

        if (size(ir_to_ic) /= n_points) then
            write(*, *) 'Size of ir_to_ic map inconsistent with number of grid points'
            error stop 102
        endif

        ! For sanity checking only: Check if any centroids have not been assigned
        ! If this is the case, it suggests that one or more input centroids
        ! are duplicates - this can happen if one does not take care with random
        ! number generation used to do the initial sampling of centroids
        allocate(npoints_per_cluster(n_centroids), source=0)

        do ir = 1, n_points
            icen = 1
            min_dist = sum((centroids(:, 1) - grid_points(:, ir))**2)
            do ic = 2, n_centroids
                dist = sum((centroids(:, ic) - grid_points(:, ir))**2)
                if (dist < min_dist) then
                    min_dist = dist
                    icen = ic
                endif
            enddo
            ir_to_ic(ir) = icen
            npoints_per_cluster(icen) = npoints_per_cluster(icen) + 1
        enddo 
        
        ! Sanity checking
        n_unassigned = 0
        do ic = 1, n_centroids
            if (npoints_per_cluster(icen) == 0) then
                n_unassigned = n_unassigned + 1
                write(*, *) 'Centroid', ic, 'was not assigned any grid points'
            endif
        enddo
        if (n_unassigned > 0) then
            write(*, *) n_unassigned, 'centroids were not assigned grid points'
            write(*, *) 'This can happen when there are two or more centroids with the same point,'
            write(*, *) 'which can occur if the same random number is generated multiple times (for example)'
            write(*, *) 'in the initial sampling of centroid points'
            error stop 102
        endif
        deallocate(npoints_per_cluster)

    end subroutine assign_points_to_centroids


    !> @brief Compute a new set of centroids.
    !!
    !! A centroid is defined as:
    !! \f[
    !!    \mathbf{r}_\mu = \frac{\sum_{\mathbf{r}_j \in C_\mu} \mathbf{r}_j w(\mathbf{r}_j)}{\sum_{\mathbf{r}_j \in C_\mu} w(\mathbf{r}_j)}
    !! \f]
    !! where \f$\mathbf{r}_j\f$ and \f$w(\mathbf{r}_j\f$ are grid points and weights restricted to the cluster \f$ C_\mu \f$, respectively.
    subroutine update_centroids(comm, grid, weight, ir_to_ic, centroids)
        type(mpi_t),  intent(inout) :: comm                     !< MPI instance
        real(real64), intent(in)    :: grid(:, :)               !< Real-space grid (n_dims, N)
        real(real64), intent(in)    :: weight(:)                !< Weights (N)
        integer,      intent(in)    :: ir_to_ic(:)              !< Map grid indices to centroid indices (N)
        
        real(real64), intent(inout) :: centroids(:, :)          !< In: centroid positions (n_dims, Ncentroids)
        !                                                         Out: Updated centroid positions
        
        integer                   :: n_dims, nr, n_centroids, ir, ic
        real(real64)              :: tmp
        real(real64), parameter :: zero_tol = 1.e-10_real64     !< Empirically defined
        real(real64), allocatable :: denominator(:)

        n_dims = size(grid, 1)
        nr = size(grid, 2)

        n_centroids = size(centroids, 2)
        allocate(denominator(n_centroids))

        ! The indexing of weight and grid must be consistent => must have the same distribution
        ! or there must be a local-to-global index map (which clearly I am not using)
        if (size(weight) /= nr) then
            write(*, *) "The number of grid points /= weight. This might indicate that their distributions &
            & differ, hence the indexing will not be consistent"
            error stop 101
        endif

        !$omp parallel default(shared)

        !$omp do simd
        do ic = 1, n_centroids
            centroids(:, ic) = 0._real64
            denominator(ic) = 0._real64
        enddo

        ! Iterate over all grid points
        !$omp do simd private(ic) reduction(+ : centroids, denominator)
        do ir = 1, nr
            ic = ir_to_ic(ir)
            ! Initially accumulate the numerator in `centroids`
            centroids(:, ic) = centroids(:, ic) + (grid(:, ir) * weight(ir))
            denominator(ic) = denominator(ic) + weight(ir)
        enddo

        !$omp end parallel

#ifdef USE_MPI         
        call MPI_ALLREDUCE(MPI_IN_PLACE, centroids, size(centroids), MPI_DOUBLE_PRECISION, MPI_SUM, comm%comm, comm%ierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, denominator, size(denominator), MPI_DOUBLE_PRECISION, MPI_SUM, comm%comm, comm%ierr)
#endif         

        ! Finish defining centroids as numerator / denominator
        ! If the code crashes here due to division by zero, it implies that the sum(weight) = 0 for all grid points
        ! in cluster ic. This can occur if the initial centroids are poorly chosen i.e.
        ! do not choice centroids in regions of no weight (such as the vacuum of a crystal cell)
        !$omp parallel do simd private(tmp) reduction(* : centroids)
        do ic = 1, n_centroids
            tmp = 1._real64 / denominator(ic)
            centroids(:, ic) = centroids(:, ic) * tmp
        enddo
        !$omp end parallel do simd

    end subroutine update_centroids

    
    !> @brief Compute the difference in two grids as abs(\mathbf{g}_1 - \mathbf{g}_2).
    !!
    !! If an component of the difference vector for grid point i is greater than
    !! the tolerance, the element points_differ(i) is set to true.
    subroutine compute_grid_difference(points, updated_points, tol, points_differ)
        real(real64), intent(in)    :: points(:, :)      !< Real-space grids (n_dims, N)
        real(real64), intent(in)    :: updated_points(:, :)
        real(real64), intent(in)    :: tol
        logical,      intent(out)   :: points_differ(:)     !< |a_i - b_i|
       
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
    !> @brief Report differences returned from `compute_grid_difference`.
    subroutine report_differences_in_grids(points, updated_points, tol, points_differ)
        real(real64), intent(in)  :: points(:, :)      !< Real-space grids (n_dims, N)
        real(real64), intent(in)  :: updated_points(:, :)
        real(real64), intent(in)  :: tol
        logical,      intent(in)  :: points_differ(:)   !< If any element of point i > tol
        
        integer, allocatable :: indices(:)
        integer :: i, j, n_unconverged
        real(real64), allocatable :: diff(:)

        indices = pack([(i, i=1,size(points_differ))], points_differ)
        n_unconverged = size(indices)
        allocate(diff(size(points, 1)))

        write(*, *) "# Current Point  ,  Prior Point  ,  |ri - r_{i-1}|  ,  tol"
        do j = 1, n_unconverged
            i = indices(j)
            diff(:) = abs(updated_points(:, i) - points(:, i))
            write(*, '(2(F16.10, X), 2(F16.10, X), 2(F16.10, X), F16.10)') updated_points(:, i), points(:, i), diff, tol
        enddo
        write(*, *) "Summary:", n_unconverged, "of out", size(points, 2), "are not converged"
   
    end subroutine report_differences_in_grids


    !> @brief Check that the sum of the weight for a given centroid/cluster is finite.
    !!
    !! If a centroid is chosen such that the total weight associated with the 
    !! cluster''s grid points is zero, updating its position will fail due to a division by zero.
    !! This can occur if an initial centroid is chosen in a region of a cell, far away
    !! from any density (for example) i.e. in the vacuum of a cell.
    subroutine check_weight_at_centroids(comm, weight, ir_to_ic, centroids, null_centroid_indices)
        type(mpi_t),  intent(inout) :: comm                     !< MPI instance
        real(real64), intent(in)    :: weight(:)                !< Weights (np)
        integer,      intent(in)    :: ir_to_ic(:)              !< Map grid indices to centroid indices (np)
        real(real64), intent(in)    :: centroids(:, :)          !< Centroids (n_dims, Ncentroids)
        integer,      allocatable, intent(out)   :: null_centroid_indices(:) !< Centroids defined where weight < tol

        integer                   :: np, n_centroids, ir, ic, n_empty
        integer,      allocatable :: tmp_indices(:)
        real(real64), allocatable :: summed_weight(:)           !< Summed weight function, for each cluster
        real(real64), parameter   :: tol = 1.e-8_real64

        np = size(weight)
        n_centroids = size(centroids, 2)
        allocate(summed_weight(n_centroids))

        !$omp parallel default(shared)

        !$omp do simd
        do ic = 1, n_centroids
            summed_weight(ic) = 0._real64
        enddo

        ! Iterate over all grid points
        !$omp do simd private(ic) reduction(+ : summed_weight)
        do ir = 1, np
            ic = ir_to_ic(ir)
            summed_weight(ic) = summed_weight(ic) + weight(ir)
        enddo

        !$omp end parallel

#ifdef USE_MPI         
        call MPI_ALLREDUCE(MPI_IN_PLACE, summed_weight, np, MPI_DOUBLE_PRECISION, MPI_SUM, comm%comm, comm%ierr)
#endif         

        allocate(tmp_indices(n_centroids))
        n_empty = 0
        do ic = 1, n_centroids
            if (abs(summed_weight(ic)) < tol) then
                n_empty = n_empty + 1
                tmp_indices(n_empty) = ic
            endif
        enddo
        if (n_empty > 0) then
            allocate(null_centroid_indices(n_empty), source=tmp_indices(1:n_empty))
        endif
        deallocate(tmp_indices)

    end subroutine check_weight_at_centroids


    !> @brief Weighted K-means clustering.
    !!
    !! The K-means algorithm divides a set of \f$N_r\f$ samples (in this case, grid points) into \f$N_\mu\f$ disjoin clusters \f$C\f$.
    !! The mean of each cluster defines the cluster centroid. Note that centroids are not, in general, points from the discrte `grid`
    !! - they can take any contiuous values that span the grid limits.
    !!
    !! ## Theory
    !! Given a grid, and some initial guess at centroids, the K-means algorithm aims to choose centroids that minimise the inertia, 
    !! or within-cluster sum-of-squares criterion:
    !! \f[
    !!    argmin \sum^{N_\mu}_{\mu=1} \sum_{\mathbf{r}_\mathbf{k} \in C_\mu} || Z(\mathbf{r}_k) - Z(\mathbf{r}_\mu)||^2
    !! \f]
    !! This implementation is based on Algorithm 1. given in [Complex-valued K-means clustering of interpolative separable density fitting algorithm 
    !! for large-scale hybrid functional enabled ab initio  molecular dynamics simulations within plane waves](https://doi.org/10.48550/arXiv.2208.07731),
    !! however it is equivalent to implementations found in packages such as [scikit-learn](https://scikit-learn.org/stable/modules/clustering.html#k-means).
    !!
    !! ## Algorithm Description
    !! The K-means algorithm consists of looping between two steps. 
    !!  1. The first step assigns each sample to its nearest centroid. See `assign_points_to_centroids`
    !!  2. The second step creates new centroids by taking the mean value of all of the samples assigned to each previous centroid. See `update_centroids`
    !! The difference between the old and the new centroids are computed, and the algorithm repeats these last two steps until this value is less than a threshold.
    subroutine weighted_kmeans(comm, grid, weight, centroids, n_iter, centroid_tol, verbose)
        type(mpi_t),  intent(inout) :: comm                      !< MPI instance
        real(real64), intent(in)    :: grid(:, :)                !< Grid    (n_dim, n_points)
        real(real64), intent(in)    :: weight(:)                 !< Weights (n_points)
        real(real64), intent(inout) :: centroids(:, :)           !< In: Initial centroids (n_dim, n_centroid)
        !                                                           Out: Final centroids
        integer,      optional, intent(in) :: n_iter             !< Optional max number of iterations
        real(real64), optional, intent(in) :: centroid_tol       !< Optional convergence criterion
        logical,      optional, intent(in) :: verbose            !< Verbosity
  
        logical          :: print_out
        integer          :: n_iterations, nr, n_dim, n_centroid, i, j, ic
        real(real64)     :: tol

        integer,      allocatable :: ir_to_ic(:), null_centroid_indices(:)
        real(real64), allocatable :: prior_centroids(:, :)
        logical,      allocatable :: points_differ(:)

        ! Optional arg assignment
        n_iterations = 200
        if (present(n_iter)) n_iterations = n_iter

        tol = 1.e-6_real64
        if (present(centroid_tol)) tol = centroid_tol

        print_out = .false.
        if (present(verbose)) then
            if (comm%is_root()) print_out = verbose
        endif

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
        allocate(ir_to_ic(nr))
        allocate(points_differ(n_centroid))

        do i = 1, n_iterations
            if (print_out) write(*, *) 'Iteration ', i

            call assign_points_to_centroids(grid, centroids, ir_to_ic)
            ! Sanity check for input centroids
            if (i == 1) then
                call check_weight_at_centroids(comm, weight, ir_to_ic,  centroids, null_centroid_indices)
                if (allocated(null_centroid_indices)) then
                    do j = 1, size(null_centroid_indices)
                        ic = null_centroid_indices(j)
                        write(*, *) 'Centroid ' , ic , 'at position' , centroids(:, ic) , ' has no associated weight'
                    enddo
                    error stop 104
                endif
            endif

            call update_centroids(comm, grid, weight, ir_to_ic, centroids)
            call compute_grid_difference(prior_centroids, centroids, tol, points_differ)

            if (any(points_differ)) then
                if (print_out) call report_differences_in_grids(prior_centroids, centroids, tol, points_differ)
                prior_centroids = centroids
            else
                ! TODO Return number of iterations used
                if (print_out) write(*, *) 'All points converged'
                return 
            endif

        enddo

        ! TODO Return number of iterations used

    end subroutine weighted_kmeans



end module kmeans_m
