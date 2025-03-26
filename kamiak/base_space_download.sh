#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=base_space_download
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G

SINGULARITY_IMAGE="/home/alexander.brown/zika-rnaseq-analysis/src/singularity_image/zika.sif"

SCRATCH_DIR="/scratch/user/alexander.brown/20250324_141326"
FASTQ_DIR=$SCRATCH_DIR"/fastq"

ENV_FILE="/home/alexander.brown/zika-rnaseq-analysis/.env"

module load singularity
singularity exec --pwd /src \
    --no-home \
    --bind $FASTQ_DIR:/src/data/fastq \
    --bind /etc:/etc \
    --env-file $ENV_FILE \
    $SINGULARITY_IMAGE \
    bs download project -i 423121844 -o /src/data/fastq --extension=fastq.gz
# 432190937