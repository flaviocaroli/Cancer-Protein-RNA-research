#!/bin/bash

CHUNK_SIZE=500
TOTAL_PROTEINS=5113
DATA_DIR=~/pipeline_proteins/data
SRC_DIR=~/pipeline_proteins/src
PIPELINE_TYPE="single"  # Change this to 'single', 'sample', or 'protein' as needed
MODEL_TYPE="xgboost"

for START in $(seq 1 $CHUNK_SIZE $TOTAL_PROTEINS); do
    END=$((START + CHUNK_SIZE - 1))
    if [ $END -gt $TOTAL_PROTEINS ]; then
        END=$TOTAL_PROTEINS
    fi

    sbatch <<EOF
#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=03:00:00
#SBATCH --partition=medium_gpu
#SBATCH --gpus=1
#SBATCH --mem=32G
#SBATCH --array=${START}-${END}


module load miniconda3
conda activate test_env

# Set CUDA environment variables
export CUDA_VISIBLE_DEVICES=0
export XGB_USE_CUDA=1

cd $SRC_DIR

PROTEIN_NAME=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" $DATA_DIR/protein_list.txt)

python pipeline.py --protein "\$PROTEIN_NAME" --pipeline $PIPELINE_TYPE --model $MODEL_TYPE
EOF

done