#!/bin/bash
#SBATCH --job-name=pipeline_job_%a       # Job name (%a for array index)
#SBATCH --output=slurm_output_%a_%j.txt  # Output file (%j = Job ID)
#SBATCH --error=slurm_error_%a_%j.txt    # Error file
#SBATCH --nodes=1                        # Number of nodes
#SBATCH --ntasks=1                       # Number of tasks (processes)
#SBATCH --cpus-per-task=4                # Number of CPU cores per task
#SBATCH --time=04:00:00                  # Max runtime
#SBATCH --mem=32G                        # Memory required
#SBATCH --array=1-10                     # Array of jobs for subsets 1 to 10
#SBATCH --mail-type=END                  # Email on completion
#SBATCH --mail-user=flavio.caroli@unibocconi.it  # Replace with your email

# Load necessary modules
module load miniconda3  # Adjust Python version if necessary

# conda create --name test_env python=3.8
# conda activate test_env
# pip install -r requirements.txt

# Run the Python script with the subset number
python ~/pipeline/src/pipeline.py --subset ${SLURM_ARRAY_TASK_ID}
