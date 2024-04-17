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
