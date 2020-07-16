#!/bin/bash
#SBATCH -c 1 # number of cores
#SBATCH -N 1 # number of nodes
#SBATCH -t 1-00:00 # runtime in D-HH:MM format
#SBATCH -p medium # partition
#SBATCH --mem=2000 # memory in MB (for all cores)
#SBATCH -o snakemake-%j.out # file for STDOUT with job ID
#SBATCH -e snakemake-%j.err # file for STDERR with job ID

source ~/.bashrc_mark # to load the `conda activate` command

conda activate cistrome-to-multivec-pipeline
snakemake --profile cistrome-explorer --local-cores 4 --keep-going --config filetype=$1 user=$2
