#!/bin/bash
#SBATCH --job-name=sccv3_array
#SBATCH --chdir=/dartfs/rc/nosnapshots/V/VaickusL-nb/EDIT_Students/users/alos/colon_work/adata
#SBATCH --partition=v100_vaickus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=04:00:00
#SBATCH --array=0-40
#SBATCH --output=logs/sccv3_%A_%a.out
#SBATCH --error=logs/sccv3_%A_%a.err

cd /dartfs/rc/nosnapshots/V/VaickusL-nb/EDIT_Students/users/alos/colon_work/adata
mkdir -p logs

module purge
module load R/4.5.1
module load intel-compilers/23.2

# Required: newer libstdc++ for Rfast/SpatialCellChat
export LD_LIBRARY_PATH=/dartfs-hpc/rc/home/k/f006fpk/anaconda3/envs/cellchat_env/lib:$LD_LIBRARY_PATH

echo "=============================="
echo "ARRAY TASK: $SLURM_ARRAY_TASK_ID"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start: $(date)"
echo "=============================="

Rscript 04_spatialcellchat_v3_array.R $SLURM_ARRAY_TASK_ID

echo "Finished: $(date)"
