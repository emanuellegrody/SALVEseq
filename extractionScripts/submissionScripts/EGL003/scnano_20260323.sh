#!/bin/bash

## IO
CONFIG='/home/egy2296/SALVEseq/extractionScripts/nf-core/scnanoseq/scnanoseq_20260323.config'
PROFILE='singularity,cluster'     # always include singularity
VERSION='1.2.2'                   # pipeline version (for reproducibility)


echo 'loading conda...'          # load required software
module purge
module load anaconda3
source ~/.bashrc                  # init conda for terminal (optional)
echo 'loading singularity...'
module load singularity/3.8.1    # load singularity for slurm execution
echo 'loading environment...'
conda activate nextflow           # activate nextflow environment for pipeline


echo 'running pipeline...'       # run the pipeline
nextflow run nf-core/scnanoseq -r $VERSION -c $CONFIG -profile $PROFILE
