#!/bin/bash
## Job Name
#SBATCH --job-name=Sort_dedup_bams
## Allocation Definition 
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes 
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=10-15:30:00
## Memory per node
#SBATCH --mem=500G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=strigg@uw.edu
## Specify the working directory for this job
#SBATCH --workdir=/gscratch/srlab/strigg/analyses/20190415_10K


#compile and sort bams
find /gscratch/srlab/strigg/analyses/20190415_10K/*deduplicated.bam| \
xargs basename -s _R1_001_bismark_bt2_pe.deduplicated.bam | xargs -I{} /gscratch/srlab/programs/samtools-1.9/samtools \
sort /gscratch/srlab/strigg/analyses/20190415_10K/{}_R1_001_bismark_bt2_pe.deduplicated.bam \
-o /gscratch/srlab/strigg/analyses/20190415_10K/{}_dedup.sorted.bam
