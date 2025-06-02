#!/bin/bash

#SBATCH --job-name=ath-hip
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=2
#SBATCH --mem=32000M
#SBATCH --time=30-00:00:00
#SBATCH --partition=hip
#SBATCH --nodelist=node3
#SBATCH --output=res_kharma_%j.out
#SBATCH --error=res_kharma_%j.err

# 加载模块
module load rocm/6.1.3
module load openmpi/5.0.5-rocm
module load ucx/1.19-rocm
module load ucc/1.3-rocm
module load cmake/3.29.2
module load hdf5/1.14.3-openmpi-5.0.5-rocm

# 基础设置
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=vader,self
export OMPI_MCA_btl_base_warn_component_unused=0

# 共享内存传输设置
export OMPI_MCA_btl_vader_single_copy_mechanism=none
export OMPI_MCA_btl_sm_eager_limit=32768
export OMPI_MCA_btl_vader_eager_limit=32768

# GPU设置和亲和性
export OMPI_MCA_rmaps_base_mapping_policy=dist:span    # 分散进程映射
export OMPI_MCA_hwloc_base_binding_policy=numa         # NUMA级别绑定
export OMPI_MCA_rmaps_base_stride=1
export OMPI_MCA_rmaps_base_device_arc_distance=1       # 考虑设备架构距离

# GPU可见性设置 - 每个进程只看到最近的GPU
export ROCR_VISIBLE_DEVICES="0,1"
export HIP_VISIBLE_DEVICES="0,1"
# 设置GPU亲和性映射
export HSA_GPU_AFFINITY="0:0,1:1"                      # 进程0使用GPU0，进程1使用GPU1

# ROCm性能设置
export HSA_ENABLE_SDMA=0
export GPU_MAX_HW_QUEUES=8
export GPU_USE_SYNC_OBJECTS=1
export HSA_ENABLE_INTERRUPT=0
export HSA_NUM_GPU_QUEUE=8                             # 限制每个GPU的队列数
export HSA_ENABLE_VM=0                                 # 禁用虚拟内存管理
export HSA_LOCAL_MEMORY_SIZE=8589934592               # 设置本地内存大小(8GB)

# 性能优化
export OMPI_MCA_mpi_yield_when_idle=0
export OMPI_MCA_mpi_event_tick_rate=-1
export OMPI_MCA_mpi_show_handle_leaks=1
export OMPI_MCA_mpi_leave_pinned=1

# 禁用CUDA支持
export OMPI_MCA_mpi_cuda_support=0
export OMPI_MCA_opal_cuda_support=0

# 运行命令
mpirun --map-by numa:PE=8 \
       --bind-to numa \
       --report-bindings \
       --display-map \
       --display-allocation \
       -x ROCR_VISIBLE_DEVICES \
       -x HIP_VISIBLE_DEVICES \
       -x HSA_GPU_AFFINITY \
       -x HSA_ENABLE_SDMA \
       -x GPU_MAX_HW_QUEUES \
       -x GPU_USE_SYNC_OBJECTS \
       -x HSA_ENABLE_INTERRUPT \
       -x HSA_NUM_GPU_QUEUE \
       -x HSA_ENABLE_VM \
       -x HSA_LOCAL_MEMORY_SIZE \
       -x OMPI_MCA_btl \
       -x OMPI_MCA_btl_vader_single_copy_mechanism \
       -x OMPI_MCA_btl_sm_eager_limit \
       -x OMPI_MCA_btl_vader_eager_limit \
       -x OMPI_MCA_mpi_leave_pinned \
       -n 2 /home/zhangzelin/projects/athenak/output/athena.hip -i /home/zhangzelin/projects/athenak/output/20241112_hip_sane_a+09735_192-192-192-7_gamma13/gr_fm_torus_sane_8_4.athinput