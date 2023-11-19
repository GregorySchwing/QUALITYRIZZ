#!/bin/bash

#SBATCH --job-name nxfl

#SBATCH -N 1

#SBATCH -n 1 
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=go2432@wayne.edu
#SBATCH -o output_%j.out
#SBATCH -e errors_%j.err
#SBATCH -t 7-0:0:0
source "${HOME}/mambaforge/etc/profile.d/mamba.sh"
# Contains nextflow and singularity
source activate nextflow
export NXF_EXECUTOR=slurm
export NXF_OPTS="-Xms6G -Xmx8G" 
mkdir -p ${HOME}/singularity_cache
export NXF_SINGULARITY_CACHEDIR=${HOME}/singularity_cache
mkdir -p ${HOME}/xdr
export XDG_RUNTIME_DIR=${HOME}/xdr
mkdir -p ${HOME}/tmpnf
export SINGULARITY_TMPDIR=${HOME}/tmpnf
nextflow run -profile potoff . --param_name nextflow.config -resume
