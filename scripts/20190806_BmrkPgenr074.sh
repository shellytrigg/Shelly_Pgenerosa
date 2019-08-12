#!/bin/bash
## Job Name
#SBATCH --job-name=BismarkAlignAS1.2I60_Pgenr074
## Allocation Definition 
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes 
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=6-00:00:00
## Memory per node
#SBATCH --mem=500G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=strigg@uw.edu
## Specify the working directory for this job
#SBATCH --workdir=/gscratch/scrubbed/strigg/analyses/20190806_v074

#align with bismark
%%bash

find /gscratch/scrubbed/strigg/Pgen_FQs/*_L004_R1_001_val_1.fq.gz \
| xargs basename -s _L004_R1_001_val_1.fq.gz | xargs -I{} /gscratch/srlab/programs/Bismark-0.19.0/bismark \
--path_to_bowtie /gscratch/srlab/programs/bowtie2-2.3.4.1-linux-x86_64 \
--samtools_path /gscratch/srlab/programs/samtools-1.9 \
--score_min L,0,-1.2 \
-I 60 \
-p 4 \
--genome /gscratch/srlab/sr320/data/geoduck/v074 \
-1 /gscratch/scrubbed/strigg/Pgen_FQs/{}_L004_R1_001_val_1.fq.gz \
-2 /gscratch/scrubbed/strigg/Pgen_FQs/{}_L004_R2_001_val_2.fq.gz \
-o /gscratch/scrubbed/strigg/analyses/20190806_v074


#run deduplicaiton
%%bash
/gscratch/srlab/programs/Bismark-0.19.0/deduplicate_bismark \
--bam -p \
/gscratch/scrubbed/strigg/analyses/20190806_v074/*.bam \
-o /gscratch/scrubbed/strigg/analyses/20190806_v074/ \
2> /gscratch/scrubbed/strigg/analyses/20190806_v074/dedup.err \
--samtools_path /gscratch/srlab/programs/samtools-1.9/

#create summary report
cat /gscratch/scrubbed/strigg/analyses/20190806_v074/*PE_report.txt | \
grep 'Mapping\ efficiency\:' | \
cat - /gscratch/scrubbed/strigg/analyses/20190806_v074/*.deduplication_report.txt | \
grep 'Mapping\ efficiency\:\|removed' \
> /gscratch/scrubbed/strigg/analyses/20190806_v074/mapping_dedup_summary.txt

#run methylation extractor
/gscratch/srlab/programs/Bismark-0.19.0/bismark_methylation_extractor \
--paired-end --bedGraph --counts --scaffolds \
--multicore 28 \
/gscratch/scrubbed/strigg/analyses/20190806_v074/*deduplicated.bam \
-o /gscratch/scrubbed/strigg/analyses/20190806_v074/ \
--samtools /gscratch/srlab/programs/samtools-1.9/samtools \
2> /gscratch/scrubbed/strigg/analyses/20190806_v074/bme.err

#create bismark reports for individual samlpes
/gscratch/srlab/programs/Bismark-0.19.0/bismark2report

#create bismark summary report for all samples
/gscratch/srlab/programs/Bismark-0.19.0/bismark2summary

#Run coverage2cytosine command to generate cytosine coverage files

find /gscratch/scrubbed/strigg/analyses/20190806_v074/*.cov.gz \
| xargs basename -s _R1_001_bismark_bt2_pe.deduplicated.bismark.cov.gz \
| xargs -I{} /gscratch/srlab/programs/Bismark-0.19.0/coverage2cytosine --gzip \
--genome_folder /gscratch/srlab/sr320/data/geoduck/v074 \
-o /gscratch/scrubbed/strigg/analyses/20190806_v074/{}_cytosine_CpG_cov_report \
/gscratch/scrubbed/strigg/analyses/20190806_v074/{}_R1_001_bismark_bt2_pe.deduplicated.bismark.cov.gz

#compile and sort bams for methylkit
find /gscratch/scrubbed/strigg/analyses/20190806_v074/*deduplicated.bam| \
xargs basename -s _R1_001_bismark_bt2_pe.deduplicated.bam | xargs -I{} /gscratch/srlab/programs/samtools-1.9/samtools \
sort /gscratch/scrubbed/strigg/analyses/20190806_v074/{}_R1_001_bismark_bt2_pe.deduplicated.bam \
-o /gscratch/scrubbed/strigg/analyses/20190806_v074/{}_dedup.sorted.bam
