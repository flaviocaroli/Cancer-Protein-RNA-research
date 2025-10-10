#!/bin/bash

CHUNK_SIZE=500
TOTAL_PROTEINS=5113
DATA_DIR="$(pwd)/data"
SRC_DIR="$(pwd)"

# Ensure the protein list file exists.
if [ ! -f "$DATA_DIR/protein_list.txt" ]; then
    echo "Error: Protein list file not found at $DATA_DIR/protein_list.txt"
    exit 1
fi

# Loop through chunks of proteins.
for START in $(seq 1 $CHUNK_SIZE $TOTAL_PROTEINS); do
    END=$((START + CHUNK_SIZE - 1))
    if [ $END -gt $TOTAL_PROTEINS ]; then
        END=$TOTAL_PROTEINS
    fi

    # Submit a job array for the current chunk.
    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=corr_job_%a
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --array=${START}-${END}
#SBATCH --output=${DATA_DIR}/logs/corr_job_%a_%A.out
#SBATCH --error=${DATA_DIR}/logs/corr_job_%a_%A.err

module load miniconda3
conda activate test_env

cd $SRC_DIR

# Get the protein name for the current array task.
PROTEIN_NAME=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $DATA_DIR/protein_list.txt)

echo "Computing correlation for protein: \$PROTEIN_NAME"
python corr_script.py --protein "\$PROTEIN_NAME"

if [ \$? -eq 0 ]; then
    echo "Successfully computed correlation for protein: \$PROTEIN_NAME"
else
    echo "Failed to compute correlation for protein: \$PROTEIN_NAME"
fi
EOF

done
