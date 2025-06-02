#!/bin/bash

#SBATCH --job-name=ath-hip
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=1
#SBATCH --mem=32000M
#SBATCH --time=30-00:00:00
#SBATCH --partition=hip
#SBATCH --nodelist=node2,node3
#SBATCH --output=res_kharma_%j.out
#SBATCH --error=res_kharma_%j.err

module load cmake/3.29.2 rocm/6.1.3

# export OMP_PROC_BIND=spread
# export OMP_PLACES=threads
export OMP_NUM_THREADS=4

srun -N 1 --gpu-bind=closest /home/zhangzelin/projects/athenak/output/athena.hip -i /home/zhangzelin/projects/athenak/output/20240929_hip_sane_a+09375_128-128-5_gamma13/gr_fm_torus_sane_8_4.athinput