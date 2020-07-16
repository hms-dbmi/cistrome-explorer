#!/bin/bash
#SBATCH -c 1 # number of cores
#SBATCH -N 1 # number of nodes
#SBATCH -t 1-00:00 # runtime in D-HH:MM format
#SBATCH -p medium # partition
#SBATCH --mem=2000 # memory in MB (for all cores)
#SBATCH -o hostname_%j.out # file for STDOUT with job ID
#SBATCH -e hostname_%j.err # file for STDERR with job ID

source ~/.bashrc_mark
conda activate cistrome-to-multivec-pipeline
snakemake --profile cistrome-explorer --config filetype=mv5 user=mk596
