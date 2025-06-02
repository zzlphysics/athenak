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

# 服务器Mi50 GPU HIP with OpenMPI-GPU
module load cmake/3.29.2 ucc/1.3-rocm ucx/1.19-rocm openmpi/5.0.5-rocm rocm/6.1.3 hdf5/1.14.3-openmpi-5.0.5-rocm
# OpenMPI-GPU支持
export UCX_MEMTYPE_CACHE=n
export UCX_IB_GPU_DIRECT_RDMA=n  # 添加这行
export UCX_TLS=sm,self,rocm_copy,rocm_ipc
export UCX_RNDV_SCHEME=put_zcopy
export UCX_RNDV_THRESH=16384
export UCX_MAX_RNDV_RAILS=1
export UCX_MEMTYPE_REG_WHOLE_ALLOC_TYPES=rocm  # 添加这行
# UCX设置 - 单机优化
export UCX_TLS=sm,self,rocm_copy  # 保留这个，因为需要GPU间通信
export UCX_MEMTYPE_CACHE=n

# 单机优化的MPI设置
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=^tcp,openib  # 禁用不需要的传输层
export OMPI_MCA_osc=ucx
cmake -D Athena_ENABLE_MPI=ON -D Kokkos_ENABLE_HIP=ON -D Kokkos_ARCH_AMD_GFX906=ON -D CMAKE_CXX_COMPILER=clang++ -D CMAKE_C_COMPILER=clang -D PROBLEM=gr_torus ..