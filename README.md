# Test Program for Weighted K-Means Algorithm, plus the Fortuno Test Framework (in Fortran)

This project implements an MPI-OMP parallelised version of the weighted K-Means clustering using Lloyd’s algorithm.
In practice, the k-means algorithm is one of the fastest clustering algorithms available, but it falls in local minima. 
That’s why it can be useful to restart it several times, or utilise (greedy) kmeans++ to improve the robustness of the
initial seeding.

Additional information can be found at the scikit learn [docs](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html)

## Serial (Well, no MPI but threaded)

```shell
# Define the project root
export KMEANS_ROOT=/Users/alexanderbuccheri/Codes/kmeans/
```

CMake configure, build, install and test commands:

```shell
cmake --fresh -B serial-cmake-build-debug -DCMAKE_BUILD_TYPE=Debug -DMPI=Off -DCMAKE_INSTALL_PREFIX=${KMEANS_ROOT}/serial_install
cmake --build serial-cmake-build-debug
ctest --test-dir ./serial-cmake-build-debug --output-on-failure
cmake --install serial-cmake-build-debug
```

Run the application test with:

```shell
./serial-cmake-build-debug/run_kmeans
```

One can also run each unit test binary manually. For example:

```shell
./serial-cmake-build-debug/test_grids "Real-space 1D grid"
```

## MPI

```shell
cmake --fresh -B cmake-build-debug -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=${KMEANS_ROOT}/mpi_install 
cmake --build cmake-build-debug
ctest --test-dir ./cmake-build-debug --output-on-failure
cmake --install cmake-build-debug
```

One can run the application test with:

```shell
mpirun -np 4 cmake-build-debug/run_kmeans
```

or have the notebook execute it.

### References

Additional `ctest` args:

* `--show-only`   List all tests
* `--verbose`     Always write to stdout

Fortuno is available [here](https://github.com/fortuno-repos/fortuno)

The MPI version is available [here](https://github.com/fortuno-repos/fortuno-mpi/)

## Centroid Algorithm MPI Implementation Description

In the MPI implementation, `assign_points_to_centroids` does not change. All centroids are defined on all processes.
Instead, `assign_points_to_centroids` expects a grid or sub-grid, local to a process. It returns `cluster_sizes` and `clusters`,
which assign all points of the sub-grid to clusters.

`update_centroids` also remains the same. It returns `centroids` of shape `(n_dim, n_centroids)`, which have only been computed
using a local set of grid points. The MPI overload of this routine `allreduce`s the centroids(n_dim, n_centroids), such that the
final centroids contain contributions from all grid points, over all processes. `centroids` are now complete on all processes.

## TODOs

* Unit tests for debug build all crash when building with `-ffpe-trap=invalid,zero,overflow,underflow` and/or
  `-finit-real=nan`. This needs investigating. One assumes the issue is in the grid routines.

* Implement "kmeans++" to find initial weights - see my python reference
  * Extend to "greedy kmeans++"

* Implement QR decomposition with pivoting - see my python reference implementation
