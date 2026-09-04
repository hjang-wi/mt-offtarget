#!/bin/bash
#SBATCH --job-name=reditool
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --array=<XX>
#SBATCH -n 7
#SBATCH -N 1
#SBATCH -p 20
#SBATCH --mem=40G
#SBATCH --time=6-12:00:00
#SBATCH --mail-user=<YOUR_EMAIL>
#SBATCH --mail-type=ALL

# Initialize conda
source ~/miniforge3/etc/profile.d/conda.sh
conda activate reditools_py2

# Confirm paths
echo "Python: $(which python)"
echo "BWA: $(which bwa)"
echo "Samtools: $(which samtools)"

# Get sample name from list
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" samples.txt)
echo "Processing sample: $SAMPLE"

python mito_offtarget_single_sample_v2.py $SAMPLE
