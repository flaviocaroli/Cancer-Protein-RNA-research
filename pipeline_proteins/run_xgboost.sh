#!/bin/bash
#SBATCH --job-name=pipeline_enet          # Job name
#SBATCH --output=enet_output_%a.txt       # Output file for each subset (%a for array index)
#SBATCH --error=enet_error_%a.txt         # Error file for each subset
#SBATCH --time=04:00:00                   # Max runtime
#SBATCH --nodes=1                         # Number of nodes
#SBATCH --gpus-per-node=1                 # Number of GPUs per node
#SBATCH --partition=gpu                   # Use the GPU partition
#SBATCH --mem=128GB                       # Memory required
#SBATCH --array=1-10                      # Array of jobs for subsets 1 to 10
#SBATCH --mail-type=END                   # Email notifications
#SBATCH --mail-user=edoardocalder@gmail.com

module load miniconda3
eval "$(conda shell.bash hook)"

# Activate the Conda environment
conda activate compbio_env

# Run the Python script with the subset number passed as SLURM_ARRAY_TASK_ID
python -u pipeline_proteins/src/pipeline_enet.py --subset ${SLURM_ARRAY_TASK_ID}
