#!/bin/bash

#SBATCH --job-name snake

#SBATCH -N 1

#SBATCH -n 1 

#SBATCH --mem=10G
#SBATCH --gpus-per-node=1
#SBATCH --ntasks-per-gpu=1
#SBATCH --nodelist=ressrv4ai8111,ressrv6ai8111

#SBATCH --mail-type=ALL

#SBATCH --mail-user=go2432@wayne.edu

#SBATCH -o output_%j.out

#SBATCH -e errors_%j.err

#SBATCH -t 7-0:0:0
eval "$(conda shell.bash hook)"
source /home6/go2432/Wolf_GOMC/workflow/mambaforge/etc/profile.d/conda.sh
conda activate snakemake-minimal
snakemake --unlock
#cd /home/scratch/Wolf_GOMC/workflow
#snakemake -c 1 --use-conda
snakemake -c 1 --use-conda --rerun-triggers mtime  --rerun-incomplete
echo $HOSTNAME
