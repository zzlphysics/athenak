#!/bin/bash

#SBATCH --job-name=ath-debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --cpus-per-task=2
#SBATCH --mem=128000M
#SBATCH --time=30-00:00:00
#SBATCH --partition=hip
#SBATCH --nodelist=node2
#SBATCH --output=res_kharma_%j.out
#SBATCH --error=res_kharma_%j.err

module load cmake/3.29.2 openmpi/5.0.3

# export OMP_PROC_BIND=spread
# export OMP_PLACES=threads
# export OMP_NUM_THREADS=4

srun -N 1 /home/zhangzelin/projects/athenak/output/athena.cpu -i /home/zhangzelin/projects/athenak/output/test/gr_fm_torus_sane_8_4.athinput time/nlim=100