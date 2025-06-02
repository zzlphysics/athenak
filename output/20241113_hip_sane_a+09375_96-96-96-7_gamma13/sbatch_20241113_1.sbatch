#!/bin/bash

#SBATCH --job-name=ath-hip
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=1
#SBATCH --mem=32000M
#SBATCH --time=30-00:00:00
#SBATCH --partition=hip
#SBATCH --nodelist=node2
#SBATCH --output=res_kharma_%j.out
#SBATCH --error=res_kharma_%j.err

module load cmake/3.29.2 ucc/1.3-rocm ucx/1.19-rocm openmpi/5.0.5-rocm rocm/6.1.3 hdf5/1.14.3-openmpi-5.0.5-rocm

# export OMP_PROC_BIND=spread
# export OMP_PLACES=threads
export OMP_NUM_THREADS=4
export UCX_MEMTYPE_CACHE=n
export UCX_IB_GPU_DIRECT_RDMA=n  # 添加这行
export UCX_TLS=sm,self,rocm_copy,rocm_ipc
export UCX_RNDV_SCHEME=put_zcopy
export UCX_RNDV_THRESH=16384
export UCX_MAX_RNDV_RAILS=1
export UCX_MEMTYPE_REG_WHOLE_ALLOC_TYPES=rocm  # 添加这行

srun -N 1 --gpu-bind=closest /home/zhangzelin/projects/athenak/output/athena.hip -i /home/zhangzelin/projects/athenak/output/20241113_hip_sane_a+09375_96-96-96-7_gamma13/gr_fm_torus_sane_8_4.athinput