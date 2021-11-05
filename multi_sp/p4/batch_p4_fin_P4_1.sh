#!/bin/bash
#SBATCH --account=def-leithen
#SBATCH --job-name=occu
#SBATCH --array=1-270
#SBATCH --mem=30G
#SBATCH --time=00-24:00

echo "Working directory is $(pwd)"
echo "Starting R at $(date)."

module load gcc r
export R_LIBS=~/Rlibs

echo "Established connection with R-libraries. Running model..."

Rscript multi_sp/p4/02_run_model_P4_1.R $SLURM_ARRAY_TASK_ID
echo "Finished R at $(date)."