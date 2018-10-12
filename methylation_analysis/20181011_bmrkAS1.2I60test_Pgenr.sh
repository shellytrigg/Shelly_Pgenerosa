#!/bin/bash
## Job Name
#SBATCH --job-name=BismarkAlignAS1.2I60test_Pgenr
## Allocation Definition 
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes 
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=04-15:30:00
## Memory per node
#SBATCH --mem=500G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=strigg@uw.edu
## Specify the working directory for this job
#SBATCH --workdir=/gscratch/srlab/strigg/analyses/20181011


#align with bismark
%%bash

find /gscratch/srlab/strigg/data/Pgenr/FASTQS/EPI-*R1* \
| xargs basename -s _001_val_1.fq.gz | xargs -I{} /gscratch/srlab/programs/Bismark-0.19.0/bismark \
--path_to_bowtie /gscratch/srlab/programs/bowtie2-2.3.4.1-linux-x86_64 \
--samtools_path /gscratch/srlab/programs/samtools-1.9 \
--score_min L,0,-1.2 \
-I 60 \
-u 10000 \
--non_directional \
--genome /gscratch/srlab/strigg/data/Pgenr/ \
-1 /gscratch/srlab/strigg/data/Pgenr/FASTQS/{}_001_val_1.fq.gz \
-2 /gscratch/srlab/strigg/data/Pgenr/FASTQS/{}_001_val_2.fq.gz \
-o /gscratch/srlab/strigg/analyses/20181011/ 


#run deduplicaiton
%%bash
/gscratch/srlab/programs/Bismark-0.19.0/deduplicate_bismark \
--bam -p \
/gscratch/srlab/strigg/analyses/20181011/*.bam \
-o /gscratch/srlab/strigg/analyses/20181011/ \
2> /gscratch/srlab/strigg/analyses/20181011/dedup.err \
--samtools_path /gscratch/srlab/programs/samtools-1.9/

#create summary report
cat /gscratch/srlab/strigg/analyses/20181011/*PE_report.txt | \
grep 'Mapping\ efficiency\:' | \
cat - /gscratch/srlab/strigg/analyses/20181011/*.deduplication_report.txt > /gscratch/srlab/strigg/analyses/20181011/EPI_mapping_dedup_summary.txt
#clean up summary report
sed 's/Mapping\ efficiency\://g' /gscratch/srlab/strigg/analyses/20181011/EPI_mapping_dedup_summary.txt | \
sed 's/Total\ number\ duplicated\ alignments\ removed\://g' | \
sed 's/ //g' | awk '{print $1}' > /gscratch/srlab/strigg/analyses/20181011/EPI_mapping_dedup_summary_clean.txt


