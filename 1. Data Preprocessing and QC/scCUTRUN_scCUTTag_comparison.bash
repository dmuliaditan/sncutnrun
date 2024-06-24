#!/usr/bin/env bash

#SBATCH -J snCUTRUN_scCUTTag
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -p normal
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 20
#SBATCH --mem-per-cpu 3500
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err


## Below is the analysis for scCUTTag GBM (primary and recurrent) data

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#Set variables and directories
longLine="--------------------"
THREADS=13

mamba activate py3
module load samtools
module load gatk

cd /scratch/users/astar/gis/muliaditand/sncutrun/public_data/sccuttag/K27me3_GBM_primary/

msg="Map with bowtie2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
bowtie2 -p "$THREADS" --end-to-end --very-sensitive --no-mixed --no-discordant \
--met-file 'K27me3_GBM_primary.txt' \
-q --phred33 -I 10 -X 700 -x /scratch/users/astar/gis/muliaditand/programmes/bowtie2/GRCh38_noalt_decoy_as/GRCh38_noalt_decoy_as \
-1 K27me3_GBM_primary_S1_L001_R1_001.fastq \
-2 K27me3_GBM_primary_S1_L001_R2_001.fastq > K27me3_GBM_primary.sam

msg="Convert sam to bam"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools view -Sb K27me3_GBM_primary.sam -@ "$THREADS" > K27me3_GBM_primary_unsort.bam

msg="Sort bam with "$SAMTOOLS""; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools sort K27me3_GBM_primary_unsort.bam -@ "$THREADS" > K27me3_GBM_primary_sorted.bam

msg="Mark and remove duplicates with Picard and sort with "$SAMTOOLS""; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
gatk MarkDuplicates \
-I K27me3_GBM_primary_sorted.bam \
-O K27me3_GBM_primary_rmdup.bam \
-M K27me3_GBM_primary_sorted_metrics.txt \
--ASSUME_SORTED true --REMOVE_DUPLICATES true

samtools sort K27me3_GBM_primary_rmdup.bam -@ "$THREADS" > K27me3_GBM_primary.bam
samtools index K27me3_GBM_primary.bam -@ "$THREADS"
  
msg="Remove old and unsorted bam and sam, and preprocessing beds"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
rm K27me3_GBM_primary_unsort.bam
rm K27me3_GBM_primary.sam
rm K27me3_GBM_primary_sorted.bam
rm K27me3_GBM_primary_rmdup.bam
 
msg="Gzip fastq"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
gzip K27me3_GBM_primary_S1_L001_R1_001.fastq 
gzip K27me3_GBM_primary_S1_L001_R2_001.fastq

msg="Call peaks with MACS2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
macs2 callpeak -t K27me3_GBM_primary.bam -f AUTO -n K27me3_GBM_primary \
--nomodel -p 5e-2 --min-length 500 --max-gap 400 --SPMR --call-summits -B

cd /scratch/users/astar/gis/muliaditand/sncutrun/public_data/sccuttag/K27me3_GBM_recurrent_rep1/

msg="Map with bowtie2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
bowtie2 -p "$THREADS" --end-to-end --very-sensitive --no-mixed --no-discordant \
--met-file 'K27me3_GBM_recurrent_rep1.txt' \
-q --phred33 -I 10 -X 700 -x /scratch/users/astar/gis/muliaditand/programmes/bowtie2/GRCh38_noalt_decoy_as/GRCh38_noalt_decoy_as \
-1 K27me3_GBM_recurrent_rep1_S1_L001_R1_001.fastq \
-2 K27me3_GBM_recurrent_rep1_S1_L001_R2_001.fastq > K27me3_GBM_recurrent_rep1.sam

msg="Convert sam to bam"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools view -Sb K27me3_GBM_recurrent_rep1.sam -@ "$THREADS" > K27me3_GBM_recurrent_rep1_unsort.bam

msg="Sort bam with "$SAMTOOLS""; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools sort K27me3_GBM_recurrent_rep1_unsort.bam -@ "$THREADS" > K27me3_GBM_recurrent_rep1_sorted.bam

msg="Mark and remove duplicates with Picard and sort with "$SAMTOOLS""; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
gatk MarkDuplicates \
-I K27me3_GBM_recurrent_rep1_sorted.bam \
-O K27me3_GBM_recurrent_rep1_rmdup.bam \
-M K27me3_GBM_recurrent_rep1_sorted_metrics.txt \
--ASSUME_SORTED true --REMOVE_DUPLICATES true

samtools sort K27me3_GBM_recurrent_rep1_rmdup.bam -@ "$THREADS" > K27me3_GBM_recurrent_rep1.bam
samtools index K27me3_GBM_recurrent_rep1.bam -@ "$THREADS"
  
msg="Remove old and unsorted bam and sam, and preprocessing beds"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
rm K27me3_GBM_recurrent_rep1_unsort.bam
rm K27me3_GBM_recurrent_rep1.sam
rm K27me3_GBM_recurrent_rep1_sorted.bam
rm K27me3_GBM_recurrent_rep1_rmdup.bam
 
msg="Gzip fastq"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
gzip K27me3_GBM_recurrent_rep1_S1_L001_R1_001.fastq 
gzip K27me3_GBM_recurrent_rep1_S1_L001_R2_001.fastq

msg="Call peaks with MACS2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
macs2 callpeak -t K27me3_GBM_recurrent_rep1.bam -f AUTO -n K27me3_GBM_recurrent_rep1 \
--nomodel -p 5e-2 --min-length 500 --max-gap 400 --SPMR --call-summits -B

cd /scratch/users/astar/gis/muliaditand/sncutrun/public_data/sccuttag/K27me3_GBM_recurrent_rep2/

msg="Map with bowtie2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
bowtie2 -p "$THREADS" --end-to-end --very-sensitive --no-mixed --no-discordant \
--met-file 'K27me3_GBM_recurrent_rep2.txt' \
-q --phred33 -I 10 -X 700 -x /scratch/users/astar/gis/muliaditand/programmes/bowtie2/GRCh38_noalt_decoy_as/GRCh38_noalt_decoy_as \
-1 K27me3_GBM_recurrent_rep2_S1_L001_R1_001.fastq \
-2 K27me3_GBM_recurrent_rep2_S1_L001_R2_001.fastq > K27me3_GBM_recurrent_rep2.sam

msg="Convert sam to bam"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools view -Sb K27me3_GBM_recurrent_rep2.sam -@ "$THREADS" > K27me3_GBM_recurrent_rep2_unsort.bam

msg="Sort bam with "$SAMTOOLS""; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools sort K27me3_GBM_recurrent_rep2_unsort.bam -@ "$THREADS" > K27me3_GBM_recurrent_rep2_sorted.bam

msg="Mark and remove duplicates with Picard and sort with "$SAMTOOLS""; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
gatk MarkDuplicates \
-I K27me3_GBM_recurrent_rep2_sorted.bam \
-O K27me3_GBM_recurrent_rep2_rmdup.bam \
-M K27me3_GBM_recurrent_rep2_sorted_metrics.txt \
--ASSUME_SORTED true --REMOVE_DUPLICATES true

samtools sort K27me3_GBM_recurrent_rep2_rmdup.bam -@ "$THREADS" > K27me3_GBM_recurrent_rep2.bam
samtools index K27me3_GBM_recurrent_rep2.bam -@ "$THREADS"
  
msg="Remove old and unsorted bam and sam, and preprocessing beds"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
rm K27me3_GBM_recurrent_rep2_unsort.bam
rm K27me3_GBM_recurrent_rep2.sam
rm K27me3_GBM_recurrent_rep2_sorted.bam
rm K27me3_GBM_recurrent_rep2_rmdup.bam
 
msg="Gzip fastq"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
gzip K27me3_GBM_recurrent_rep2_S1_L001_R1_001.fastq 
gzip K27me3_GBM_recurrent_rep2_S1_L001_R2_001.fastq

msg="Call peaks with MACS2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
macs2 callpeak -t K27me3_GBM_recurrent_rep2.bam -f AUTO -n K27me3_GBM_recurrent_rep2 \
--nomodel -p 5e-2 --min-length 500 --max-gap 400 --SPMR --call-summits -B
