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


cat cmake.sh 
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


cmake -D CMAKE_BUILD_TYPE=Release -D Kokkos_ENABLE_CUDA=ON -D Kokkos_ARCH_VOLTA70=ON -D CMAKE_CXX_COMPILER=/home/vipuser/Projects/athenak/kokkos/bin/nvcc_wrapper -D Athena_SINGLE_PRECISION=OFF -D Athena_ENABLE_MPI=ON -D PROBLEM=gr_torus ..

步骤 1: 检查内核是否支持 CMA
vipuser@ubuntu22:~/Projects/athenak/output$ ompi_info --param btl vader | grep -A 5 -B 5 cma
vipuser@ubuntu22:~/Projects/athenak/output$ grep CONFIG_CROSS_MEMORY_ATTACH /boot/config-$(uname -r)
CONFIG_CROSS_MEMORY_ATTACH=y
vipuser@ubuntu22:~/Projects/athenak/output$ cat /proc/sys/kernel/yama/ptrace_scope
1
如果 ptrace_scope 不为 0，则 CMA 被阻止。
修复：
*.bash
Shell
echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
vipuser@ubuntu22:~/Projects/athenak/output$ echo 0 | sudo tee /proc/sys/kernel/yama/ptrace_scope
[sudo] password for vipuser: 
0

mpirun \
  --mca pml ucx \          # 强制使用 UCX PML
  --mca btl ^vader,smcuda \ # 禁用旧的 BTL
  -x UCX_TLS=sm,cuda_copy,cuda_ipc,gdr_copy \ # 启用 CUDA 相关传输
  -x UCX_CUDA_IPC_CACHE=n \ # 避免 IPC 句柄缓存问题（可选）
  your_program


  安装带 UCX 支持的 Open MPI：
*.bash
Shell
spack install ucx +cuda cuda_arch=70
spack install openmpi +cuda fabrics=ucx cuda_arch=70
验证安装：
*.bash
Shell
spack load openmpi
ompi_info --parsable --all | grep -i ucx
使用 UCX 运行：
*.bash
Shell
mpirun --mca pml ucx --mca btl ^vader,smcuda -x UCX_TLS=sm,cuda_copy,cuda_ipc your_program

编译运行时
# 1. 启用 CUDA 相关传输层
export UCX_TLS=sm,cuda_copy,cuda_ipc,mm

# 2. 禁用可能冲突的传输层（如 IB，除非你有 InfiniBand）
export UCX_TLS_ALLOW_LIST=sm,cuda_copy,cuda_ipc,mm

# 3. 避免 CUDA IPC 句柄缓存问题（云环境常见）
export UCX_CUDA_IPC_CACHE=n

# 4. 设置 CUDA 设备（强烈建议在代码中显式调用 cudaSetDevice）
#    但也可通过环境变量辅助（非必需）
export CUDA_VISIBLE_DEVICES=0,1

mpirun -n 2 --mca pml ucx --mca btl ^vader,smcuda -x UCX_TLS=sm,cuda_copy,cuda_ipc athena.sm70.ucx -i gr_fm_torus_mad_192_7_1024.athinput 





zhangzelin@ubuntux99:~/Projects/athenak/build$ cmake -D CMAKE_BUILD_TYPE=Release -D Kokkos_ENABLE_CUDA=ON -D Kokkos_ARCH_VOLTA70=ON -D CMAKE_CXX_COMPILER=/home/zhangzelin/Projects/athenak/kokkos/bin/nvcc_wrapper -D Athena_SINGLE_PRECISION=OFF -D Athena_ENABLE_MPI=ON -D PROBLEM=gr_torus ..

