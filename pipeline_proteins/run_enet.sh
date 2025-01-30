#!/bin/bash

CHUNK_SIZE=500
TOTAL_PROTEINS=5113
DATA_DIR=~/pipeline_proteins/data
SRC_DIR=~/pipeline_proteins/src
PIPELINE_TYPE="single"  # Change this to 'single', 'sample', or 'protein' as needed

for START in $(seq 1 $CHUNK_SIZE $TOTAL_PROTEINS); do
    END=$((START + CHUNK_SIZE - 1))
    if [ $END -gt $TOTAL_PROTEINS ]; then
        END=$TOTAL_PROTEINS
    fi

    sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=pipeline_job_%a
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --array=${START}-${END}

module load miniconda3
conda activate test_env

cd $SRC_DIR

PROTEIN_NAME=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $DATA_DIR/protein_list.txt)

python pipeline_enet.py --protein "\$PROTEIN_NAME" --pipeline $PIPELINE_TYPE
EOF

done