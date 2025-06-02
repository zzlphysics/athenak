#!/bin/bash

#SBATCH --job-name=ath-k
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=18
#SBATCH --mem=410000M
#SBATCH --time=60-00:00:00
#SBATCH --partition=hip
#SBATCH --nodelist=node3
#SBATCH --output=res_kharma_%j.out
#SBATCH --error=res_kharma_%j.err

export OMP_NUM_THREADS=18

# 加载模块
module load cmake/3.29.2 openmpi/5.0.3
echo "Number of threads: $OMP_NUM_THREADS"
srun -N 1 -n 4 /home/zhangzelin/projects/athenak/output/athena.cpu \
       -i /home/zhangzelin/projects/athenak/output/20241114_cpu_sane_a+09375_2048-2048-2048-3_gamma13/gr_fm_torus_sane_8_4.athinput