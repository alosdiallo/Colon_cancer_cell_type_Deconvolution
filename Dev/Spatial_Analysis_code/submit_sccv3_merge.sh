#!/bin/bash
#SBATCH --job-name=sccv3_merge
#SBATCH --chdir=/dartfs/rc/nosnapshots/V/VaickusL-nb/EDIT_Students/users/alos/colon_work/adata
#SBATCH --partition=v100_vaickus
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=logs/sccv3_merge_%j.out
#SBATCH --error=logs/sccv3_merge_%j.err

cd /dartfs/rc/nosnapshots/V/VaickusL-nb/EDIT_Students/users/alos/colon_work/adata
mkdir -p logs

module purge
module load R/4.5.1
module load intel-compilers/23.2

export LD_LIBRARY_PATH=/dartfs-hpc/rc/home/k/f006fpk/anaconda3/envs/cellchat_env/lib:$LD_LIBRARY_PATH

echo "=============================="
echo "MERGE AND COMPARE"
echo "Start: $(date)"
echo "=============================="

Rscript 05_merge_and_compare.R

echo "Finished: $(date)"
