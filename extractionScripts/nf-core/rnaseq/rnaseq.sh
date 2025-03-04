#!/bin/bash


## IO
CONFIG='rnaseq.config' 		  # specify your .config file here
PROFILE='singularity,cluster'	  # always include singularity;
  					  # remove "cluster" to run on your local node          				  	  # instead of SLURM
VERSION='3.17.0'            	  # pipeline version (for reproducibility)	 


echo 'loading conda...' 		  # load required software
module purge
module load anaconda3               
source ~/.bashrc                    # init conda for terminal (optional)
echo 'loading singularity...'
module load singularity/3.8.1       # load singularity for slurm execution
echo 'loading environment...'
conda activate nextflow             # activate nextflow environment for pipeline


echo 'running pipeline...'          # run the pipeline
nextflow run nf-core/rnaseq -r $VERSION -c $CONFIG -profile $PROFILE