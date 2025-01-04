#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=goatools_array
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=1-00:00:00
#SBATCH --array=0,1,5,7%10
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=10G

SINGULARITY_IMAGE="/home/alexander.brown/zika-rnaseq-analysis/src/singularity_image/zika.sif"

HOST_APP_DIR="/home/alexander.brown/zika-rnaseq-analysis/src/app"

HOST_DATA_DESEQ_DIR="/home/alexander.brown/zika-rnaseq-analysis/src/data/pydeseq2"
HOST_DATA_GO_DIR="/home/alexander.brown/zika-rnaseq-analysis/src/data/go_analysis"

module load singularity
singularity exec --pwd /src \
    --no-home \
    --bind $HOST_APP_DIR:/src/app \
    --bind $HOST_DATA_DESEQ_DIR:/src/data/pydeseq2 \
    --bind $HOST_DATA_GO_DIR:/src/data/go_analysis \
    $SINGULARITY_IMAGE \
    python3 -u /src/app/goatools_analysis.py \
    -t 9407 -g $SLURM_ARRAY_TASK_ID \
    -b "/src/data/go_analysis/ncbi_gene_results_9407.txt" \
    -s "/src/data/pydeseq2/degpatterns/R06E_MR_vs_PRV_gene_clusters.csv" \
    -o "/src/data/go_analysis/R06E_MR_vs_PRV"