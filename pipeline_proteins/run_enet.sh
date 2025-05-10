#!/bin/bash

CHUNK_SIZE=500
TOTAL_PROTEINS=5113 # 3000
DATA_DIR=~/pipeline_proteins/data
SRC_DIR=~/pipeline_proteins/src
PIPELINE_TYPE="fixed_percentile"  # Change this to 'fixed_percentile', 'variable_p', or 'fixed_p' as needed
MODEL_TYPE="elasticnet"           # Change this to 'elasticnet' or 'xgboost' as needed

# Ensure the protein list file exists
if [ ! -f "$DATA_DIR/protein_list.txt" ]; then
    echo "Error: Protein list file not found at $DATA_DIR/protein_list.txt"
    exit 1
fi

# Loop through chunks of proteins
for START in $(seq 1 $CHUNK_SIZE $TOTAL_PROTEINS); do
    END=$((START + CHUNK_SIZE - 1))
    if [ $END -gt $TOTAL_PROTEINS ]; then
        END=$TOTAL_PROTEINS
    fi

    # Submit a job array for the current chunk
    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=pipeline_job_%a
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --array=${START}-${END}
#SBATCH --output=${DATA_DIR}/logs/pipeline_job_%a_%A.out
#SBATCH --error=${DATA_DIR}/logs/pipeline_job_%a_%A.err

# Load necessary modules
module load miniconda3
source activate base
conda activate test_env

# Navigate to the source directory
cd $SRC_DIR

# Get the protein name for the current array task
PROTEIN_NAME=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $DATA_DIR/protein_list.txt)

# Run the pipeline for the current protein
echo "Running pipeline for protein: \$PROTEIN_NAME"
python pipeline.py --protein "\$PROTEIN_NAME" --pipeline $PIPELINE_TYPE --model $MODEL_TYPE

# Check if the pipeline ran successfully
if [ \$? -eq 0 ]; then
    echo "Pipeline completed successfully for protein: \$PROTEIN_NAME"
else
    echo "Pipeline failed for protein: \$PROTEIN_NAME"
fi
EOF

done