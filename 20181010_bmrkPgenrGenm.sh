#!/bin/bash
## Job Name
#SBATCH --job-name=BuildPgenrGenome
## Allocation Definition 
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes 
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=0-23:30:00
## Memory per node
#SBATCH --mem=500G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=strigg@uw.edu
## Specify the working directory for this job
#SBATCH --workdir=/gscratch/srlab/strigg/data/Pgenr

#build genome
%%bash

/gscratch/srlab/programs/Bismark-0.19.0/bismark_genome_preparation \
--path_to_bowtie /gscratch/srlab/programs/bowtie2-2.3.4.1-linux-x86_64 \
2> /gscratch/srlab/strigg/data/Pgenr/genomeprep.err

