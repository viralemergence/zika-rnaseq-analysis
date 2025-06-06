#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=star_array_new_world
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=7-00:00:00
#SBATCH --array=0-35:1%10
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=30G

SINGULARITY_IMAGE="/home/alexander.brown/zika-rnaseq-analysis/src/singularity_image/zika.sif"

HOST_APP_DIR="/home/alexander.brown/zika-rnaseq-analysis/src/app"
HOST_GENOME_DIR="/home/alexander.brown/zika-rnaseq-analysis/src/genomes"

SCRATCH_DIR="/scratch/user/alexander.brown/20241028_213454"
PROCESSED_FASTQ_DIR=$SCRATCH_DIR"/processed_fastq"
STAR_RESULTS_DIR=$SCRATCH_DIR"/star"

HOST_DATA_STAR_DIR="/home/alexander.brown/zika-rnaseq-analysis/src/data/star"

module load singularity
singularity exec --pwd /src \
    --no-home \
    --bind $HOST_APP_DIR:/src/app \
    --bind $HOST_GENOME_DIR:/src/genomes \
    --bind $PROCESSED_FASTQ_DIR:/src/data/processed_fastq \
    --bind $STAR_RESULTS_DIR:/src/data/star_results \
    --bind $HOST_DATA_STAR_DIR:/src/data/star \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/star_orchestration.py \
    -i /src/data/processed_fastq -f /src/data/star/sample_file_names_new_world.txt \
    -g /src/genomes/artibeus_star_index -o /src/data/star_results -j $SLURM_ARRAY_TASK_ID