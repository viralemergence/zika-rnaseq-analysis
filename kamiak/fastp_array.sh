#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=fastp_array
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=7-00:00:00
#SBATCH --array=0-16:1%10
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G

ZIKA_IMAGE="/home/alexander.brown/zika-rnaseq-analysis/src/singularity_image/zika.sif"

HOST_APP_DIR=/home/alexander.brown/zika-rnaseq-analysis/src/app

FASTQ_DIR="/scratch/user/alexander.brown/20241022_231236/fastq"
PROCESSED_FASTQ_DIR="/scratch/user/alexander.brown/20241022_231236/processed_fastq"

module load singularity
singularity exec --pwd /src \
    --no-home \
    --bind $HOST_APP_DIR:/src/app \
    --bind $FASTQ_DIR:/src/data/fastq \
    --bind $PROCESSED_FASTQ_DIR:/src/data/processed_fastq \
    $ZIKA_IMAGE \
    python3 -u /src/app/fastp_orchestration.py \
    -i /src/data/fastq -o /src/data/processed_fastq \
    -j $SLURM_ARRAY_TASK_ID