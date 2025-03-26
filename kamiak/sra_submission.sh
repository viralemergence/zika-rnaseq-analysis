#!/bin/bash
#SBATCH --partition=remi
#SBATCH --job-name=sra_submission
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alexander.brown@wsu.edu
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G

ENV_FILE=$1
. $ENV_FILE

cd $HOST_SRA_FASTQ_DIR

ftp -inv $SRA_FTP_ADDRESS <<EOF
user $SRA_FTP_USER $SRA_FTP_PASSWORD
cd $SRA_REMOTE_DIR/zika_rnaseq
mput *.fastq.gz
exit
EOF