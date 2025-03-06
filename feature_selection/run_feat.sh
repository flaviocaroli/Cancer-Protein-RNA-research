#!/bin/bash

# Adjust these as needed
CHUNK_SIZE=500
TOTAL_PROTEINS=5113
DATA_DIR="$(pwd)/data"          # Points to the 'data' subfolder in your current directory
LOG_DIR="${DATA_DIR}/logs"      # Where SLURM logs will be written
SCRIPT_DIR="$(pwd)"             # The directory containing your .py scripts
PROTEIN_LIST="${DATA_DIR}/protein_list.txt"
PERCENTILE=10                   # Fixed percentile for feature selection

# Create logs folder if it doesn't exist
mkdir -p "${LOG_DIR}"

# Make sure the protein list file exists
if [ ! -f "$PROTEIN_LIST" ]; then
    echo "Error: Protein list file not found at $PROTEIN_LIST"
    exit 1
fi

# Loop through the protein list in chunks
for START in $(seq 1 $CHUNK_SIZE $TOTAL_PROTEINS); do
    END=$((START + CHUNK_SIZE - 1))
    if [ $END -gt $TOTAL_PROTEINS ]; then
        END=$TOTAL_PROTEINS
    fi

    # Submit a SLURM job array covering [START..END]
    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=feat_sel_%a
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --array=${START}-${END}
#SBATCH --output=${LOG_DIR}/feat_sel_%a_%A.out
#SBATCH --error=${LOG_DIR}/feat_sel_%a_%A.err

# Load your conda or other modules as needed
module load miniconda3
conda activate test_env

# Move to the directory with your Python scripts
cd $SCRIPT_DIR

# Figure out which protein to process for this task
PROTEIN_NAME=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $PROTEIN_LIST)

echo "Running feature selection for protein: \$PROTEIN_NAME"

# Run the Python script that performs feature selection
python selection_script.py --protein "\$PROTEIN_NAME" --percentile $PERCENTILE

if [ \$? -eq 0 ]; then
    echo "Feature selection completed successfully for protein: \$PROTEIN_NAME"
else
    echo "Feature selection failed for protein: \$PROTEIN_NAME"
fi
EOF

done
