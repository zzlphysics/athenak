#!/bin/bash

#SBATCH --job-name=ath-hip
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=65536M
#SBATCH --time=30-00:00:00
#SBATCH --partition=hip
#SBATCH --nodelist=node2
#SBATCH --output=res_kharma_%j.out
#SBATCH --error=res_kharma_%j.err

module load cmake/3.29.2 rocm/6.1.3

# export OMP_PROC_BIND=spread
# export OMP_PLACES=threads
# export OMP_NUM_THREADS=4
export OMP_NUM_THREADS=4
export KOKKOS_NUM_DEVICES=$SLURM_NTASKS_PER_NODE

srun -N 1 --partition=hip /home/zhangzelin/projects/athenak/output/athena.hip -r /home/zhangzelin/projects/athenak/output/mad/rst/torus.00070.rst