# 3060 GPU
cmake -D Kokkos_ENABLE_CUDA=On -D Kokkos_ARCH_AMPERE86=On -D CMAKE_CXX_COMPILER=/home/zzl/Project/athenak/kokkos/bin/nvcc_wrapper -D PROBLEM=gr_torus ../

# AMD 5800H CPU
cmake -D PROBLEM=gr_torus -D Athena_ENABLE_MPI=ON -D Kokkos_ARCH_ZEN3=On ../

# 服务器 Mi50 GPU HIP
module load cmake/3.29.2 rocm/6.1.3
cmake -D Kokkos_ENABLE_HIP=On -D Kokkos_ARCH_AMD_GFX906=On -D CMAKE_CXX_COMPILER=clang++ -D CMAKE_C_COMPILER=clang -D PROBLEM=gr_torus ..

# 服务器CPU
module load cmake/3.29.2 openmpi/5.0.3
cmake -D PROBLEM=gr_torus -D Athena_ENABLE_MPI=ON -D Kokkos_ARCH_NATIVE=On ../