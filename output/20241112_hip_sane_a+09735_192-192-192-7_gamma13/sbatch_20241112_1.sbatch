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
export OMPI_MCA_pml=ob1                 # 使用ob1点对点消息层
export OMPI_MCA_btl=vader,self          # 使用vader和self点对点消息层
export OMPI_MCA_btl_base_warn_component_unused=0     # 禁用未使用的组件警告

# 共享内存传输设置
export OMPI_MCA_btl_vader_single_copy_mechanism=none  # 禁用vader的单拷贝机制
export OMPI_MCA_btl_sm_eager_limit=32768             # 设置共享内存eager传输限制
export OMPI_MCA_btl_vader_eager_limit=32768            # 设置vader传输限制

# GPU设置
export HIP_VISIBLE_DEVICES=0,1           # 设置可见的GPU设备
export ROCR_VISIBLE_DEVICES=0,1          # 设置可见的GPU设备
export HSA_ENABLE_SDMA=0                  # 禁用SDMA
export GPU_MAX_HW_QUEUES=8                # 设置GPU硬件队列数
export GPU_USE_SYNC_OBJECTS=1               # 使用同步对象

# MPI进程绑定和亲和性
export OMPI_MCA_rmaps_base_mapping_policy=slot  # 设置映射策略
export OMPI_MCA_hwloc_base_binding_policy=core   # 设置绑定策略
export OMPI_MCA_rmaps_base_stride=1             # 设置步长
export OMPI_MCA_mpi_warn_on_fork=0              # 禁用fork警告

# 性能优化
export OMPI_MCA_mpi_yield_when_idle=0        # 禁用空闲时让出CPU
export OMPI_MCA_mpi_event_tick_rate=-1       # 禁用事件tick
export OMPI_MCA_mpi_show_handle_leaks=1      # 显示句柄泄露
export OMPI_MCA_mpi_leave_pinned=1           # 禁用pinned内存

# 禁用一些可能导致问题的功能
export HSA_ENABLE_INTERRUPT=0                # 禁用中断
export OMPI_MCA_mpi_cuda_support=0           # 禁用CUDA支持
export OMPI_MCA_opal_cuda_support=0          # 禁用CUDA支持

# 运行命令
mpirun --map-by slot:PE=8 \
       --bind-to core \
       --report-bindings \
       -x HIP_VISIBLE_DEVICES \
       -x ROCR_VISIBLE_DEVICES \
       -x HSA_ENABLE_SDMA \
       -x GPU_MAX_HW_QUEUES \
       -x GPU_USE_SYNC_OBJECTS \
       -x HSA_ENABLE_INTERRUPT \
       -x OMPI_MCA_btl \
       -x OMPI_MCA_btl_vader_single_copy_mechanism \
       -x OMPI_MCA_btl_sm_eager_limit \
       -x OMPI_MCA_btl_vader_eager_limit \
       -x OMPI_MCA_mpi_leave_pinned \
       -x OMPI_MCA_mpi_cuda_support \
       -x OMPI_MCA_opal_cuda_support \
       --display-map \
       --display-allocation \
       -n 2 /home/zhangzelin/projects/athenak/output/athena.hip -i /home/zhangzelin/projects/athenak/output/20241112_hip_sane_a+09735_192-192-192-7_gamma13/gr_fm_torus_sane_8_4.athinput