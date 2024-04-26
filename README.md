# Test Program for Weighted K-Means Algorithm, plus the Fortuno Test Framework (in Fortran)

CMake configure, build and test commands:

```shell
cmake --fresh -B cmake-build-debug   
cmake --build cmake-build-debug
ctest --test-dir ./cmake-build-debug --output-on-failure --verbose
```

One can also run each test binary manually:

```shell
./cmake-build-debug/test_kmeans
```

Fortuno is available [here](https://github.com/fortuno-repos/fortuno)

### Centroid Algorithm

#### MPI Implementation

In the MPI implementation, `assign_points_to_centroids` does not change. All centroids are defined on all processes.
Instead, `assign_points_to_centroids` expects a grid or sub-grid, local to a process. It returns `cluster_sizes` and `clusters`,
which assign all points of the sub-grid to clusters.

`update_centroids` also remains the same. It returns `centroids` of shape `(n_dim, n_centroids)`, which have only been computed
using a local set of grid points. The MPI overload of this routine `allreduce`s the centroids(n_dim, n_centroids), such that the
final centroids contain contributions from all grid points, over all processes. `centroids` are now complete on all processes.

#### TODOs

* Implement "kmeans++" to find initial weights - see my python reference v
  * Extend to "greedy kmeans++"

### QR Decomposition

* Implement QR decomposition with pivoting - see my python reference implementation
