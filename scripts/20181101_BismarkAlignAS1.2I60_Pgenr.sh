#!/bin/bash
## Job Name
#SBATCH --job-name=BismarkAlignAS1.2I60_Pgenr
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
#SBATCH --workdir=/gscratch/srlab/strigg/analyses/20181101


#align with bismark
%%bash

find /gscratch/srlab/strigg/data/Pgenr/FASTQS/EPI-*R1* \
| xargs basename -s _R1_001_val_1.fq.gz | xargs -I{} /gscratch/srlab/programs/Bismark-0.19.0/bismark \
--path_to_bowtie /gscratch/srlab/programs/bowtie2-2.3.4.1-linux-x86_64 \
--samtools_path /gscratch/srlab/programs/samtools-1.9 \
--score_min L,0,-1.2 \
-I 60 \
-p 28 \
--genome /gscratch/srlab/strigg/data/Pgenr/ \
-1 /gscratch/srlab/strigg/data/Pgenr/FASTQS/{}_R1_001_val_1.fq.gz \
-2 /gscratch/srlab/strigg/data/Pgenr/FASTQS/{}_R2_001_val_2.fq.gz \
-o /gscratch/srlab/strigg/analyses/20181101/ 


#run deduplicaiton
%%bash
/gscratch/srlab/programs/Bismark-0.19.0/deduplicate_bismark \
--bam -p \
/gscratch/srlab/strigg/analyses/20181101/*.bam \
-o /gscratch/srlab/strigg/analyses/20181101/ \
2> /gscratch/srlab/strigg/analyses/20181101/dedup.err \
--samtools_path /gscratch/srlab/programs/samtools-1.9/

#create summary report
cat /gscratch/srlab/strigg/analyses/20181101/*PE_report.txt | \
grep 'Mapping\ efficiency\:' | \
cat - /gscratch/srlab/strigg/analyses/20181101/*.deduplication_report.txt | \
grep 'Mapping\ efficiency\:\|removed' \
> /gscratch/srlab/strigg/analyses/20181101/EPI_mapping_dedup_summary.txt

#compile and sort bams
find /gscratch/srlab/strigg/analyses/20181101/*deduplicated.bam| \
xargs basename -s _R1_001_val_1_bismark_bt2_pe.deduplicated.bam | xargs -I{} /gscratch/srlab/programs/samtools-1.9/samtools \
sort /gscratch/srlab/strigg/analyses/20181101/{}_R1_001_val_1_bismark_bt2_pe.deduplicated.bam \
-o /gscratch/srlab/strigg/analyses/20181101/{}_dedup.sorted.bam

#run methylation extractor
/gscratch/srlab/programs/Bismark-0.19.0/bismark_methylation_extractor \
--paired-end --bedGraph --counts --scaffolds \
--multicore 28 \
/gscratch/srlab/strigg/analyses/20181101/*deduplicated.bam \
-o /gscratch/srlab/strigg/analyses/20181101/ \
--samtools /gscratch/srlab/programs/samtools-1.9/samtools \
2> /gscratch/srlab/strigg/analyses/20181101/bme.err

#create bismark reports for individual samlpes
/gscratch/srlab/programs/Bismark-0.19.0/bismark2report

#create bismark summary report for all samples
/gscratch/srlab/programs/Bismark-0.19.0/bismark2summary
