# Test Program for Weighted K-Means Algorithm, plus the Fortuno Test Framework (in Fortran)

## Serial (Well, threaded, no MPI)

CMake configure, build and test commands:

```shell
cmake --fresh -B serial-cmake-build-debug -DCMAKE_BUILD_TYPE=Debug -DMPI=Off
cmake --build serial-cmake-build-debug
ctest --test-dir ./serial-cmake-build-debug --output-on-failure --verbose
```

Run the application test with:

```shell
./serial-cmake-build-debug/run_kmeans
```

## MPI

```shell
cmake --fresh -B cmake-build-debug -DCMAKE_BUILD_TYPE=Debug  
cmake --build cmake-build-debug
ctest --test-dir ./cmake-build-debug --output-on-failure
```

ctest will fail for MPI build because one needs to use the MPI routines provided by the framework:

```shell
# Executing test items
.*** The MPI_Comm_f2c() function was called before MPI_INIT was invoked.
*** This is disallowed by the MPI standard.
*** Your MPI job will now abort.
```

One can also run each test binary manually:

```shell
./cmake-build-debug/test_kmeans
```

Run the application test with:

```shell
mpirun -np 4 cmake-build-debug/run_kmeans
```

Fortuno is available [here](https://github.com/fortuno-repos/fortuno)

### Build Issues

* Unit tests for debug build all crash when building with `-ffpe-trap=invalid,zero,overflow,underflow` and/or
  `-finit-real=nan`. This needs investigating. One assumes the issue is in the grid routines.

### Centroid Algorithm

#### MPI Implementation

In the MPI implementation, `assign_points_to_centroids` does not change. All centroids are defined on all processes.
Instead, `assign_points_to_centroids` expects a grid or sub-grid, local to a process. It returns `cluster_sizes` and `clusters`,
which assign all points of the sub-grid to clusters.

`update_centroids` also remains the same. It returns `centroids` of shape `(n_dim, n_centroids)`, which have only been computed
using a local set of grid points. The MPI overload of this routine `allreduce`s the centroids(n_dim, n_centroids), such that the
final centroids contain contributions from all grid points, over all processes. `centroids` are now complete on all processes.

#### TODOs

* Implement "kmeans++" to find initial weights - see my python reference
  * Extend to "greedy kmeans++"

### QR Decomposition

* Implement QR decomposition with pivoting - see my python reference implementation
