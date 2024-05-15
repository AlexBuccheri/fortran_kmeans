# Source code
target_sources(lib_kmeans PUBLIC
       kmeans_m.F90
       grids_m.f90
       maths_m.f90
       mpi_m.F90
       )

target_sources(fortuno_wrapper PRIVATE
       fortuno_interface_m.F90
       )
