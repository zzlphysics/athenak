# 3060 GPU
cmake -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_AMPERE86=On -D CMAKE_CXX_COMPILER=/home/zzl/Project/athenak/kokkos/bin/nvcc_wrapper -D PROBLEM=gr_torus ../

# AMD 5800H CPU
cmake -D PROBLEM=gr_torus -D Athena_ENABLE_MPI=ON -D Kokkos_ARCH_ZEN3=On ../