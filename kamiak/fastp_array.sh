#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=fastp_array
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=7-00:00:00
#SBATCH --array=0-94:1%14
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=30G

ZIKA_IMAGE="/home/alexander.brown/zika-rnaseq-analysis/src/singularity_image/zika.sif"

HOST_APP_DIR="/home/alexander.brown/zika-rnaseq-analysis/src/app"

SCRATCH_DIR="/scratch/user/alexander.brown/20241028_213454"
FASTQ_DIR=$SCRATCH_DIR"/fastq"
PROCESSED_FASTQ_DIR=$SCRATCH_DIR"/processed_fastq"

R1_ADAPTER_SEQ="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
R2_ADAPTER_SEQ="AGATCGGAAGAGCGTCGTGTAGGGAAAGA"

module load singularity
singularity exec --pwd /src \
    --no-home \
    --bind $HOST_APP_DIR:/src/app \
    --bind $FASTQ_DIR:/src/data/fastq \
    --bind $PROCESSED_FASTQ_DIR:/src/data/processed_fastq \
    $ZIKA_IMAGE \
    python3 -u /src/app/fastp_orchestration.py \
    -i /src/data/fastq -o /src/data/processed_fastq \
    -j $SLURM_ARRAY_TASK_ID \
    -r1a $R1_ADAPTER_SEQ -r2a $R2_ADAPTER_SEQ