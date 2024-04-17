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

## TODOs

### Centroid algorithm

* Add MPI distribution routine and test `assign_points_to_centroids` with
  distribution of outer loop
  * Add routine to distribute loops
  * I'm not happy with the first pass implementation
  * `norm2` should *probably* be replaced
  * Evaluate the use of OMP
* Implement routine `update_centroids`
* Implement function `points_are_converged`
* Implement routine `weighted_kmeans`, which calls those above
* Test on the Gaussian distribution demonstrated by `test_generate_gaussian_2d` in `test_utils.f90`
* Implement "kmeans++" to find initial weights - see my python reference v
  * Extend to "greedy kmeans++"

### QR Decomposition

* Implement QR decomposition with pivoting - see my python reference implementation
