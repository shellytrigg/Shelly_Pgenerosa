#!/bin/bash
## Job Name
#SBATCH --job-name=BismarkCoverage2Cytosine
## Allocation Definition 
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes 
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=10-15:30:00
## Memory per node
#SBATCH --mem=100G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=strigg@uw.edu
## Specify the working directory for this job
#SBATCH --workdir=/gscratch/srlab/strigg/analyses/20190415_10K

#Run coverage2cytosine command to generate cytosine coverage files

%%bash

find /gscratch/srlab/strigg/analyses/20190415_10K/*.cov.gz \
| xargs basename -s _R1_001_bismark_bt2_pe.deduplicated.bismark.cov.gz | xargs -I{} /gscratch/srlab/programs/Bismark-0.19.0/coverage2cytosine --gzip --genome_folder /gscratch/srlab/strigg/data/Pgenr/Pgenr_070_10K/ -o /gscratch/srlab/strigg/analyses/20190415_10K/{}_cytosine_CpG_cov_report /gscratch/srlab/strigg/analyses/20190415_10K/{}_R1_001_bismark_bt2_pe.deduplicated.bismark.cov.gz

