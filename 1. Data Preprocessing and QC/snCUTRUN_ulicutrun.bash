#!/usr/bin/env bash

#SBATCH -J snCUTRUN_ulicutrun
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

cd /scratch/users/astar/gis/muliaditand/sncutrun/public_data/ulicut_run

msg="Map with bowtie2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
bowtie2 -p "$THREADS" --end-to-end --very-sensitive --no-mixed --no-discordant \
--met-file 'Uli_Cut_Run_Mouse_embryo_H3K4me3_rep1.txt' \
-q --phred33 -I 10 -X 700 -x /scratch/users/astar/gis/muliaditand/programmes/bowtie2/GRCm38/GRCm38 \
-1 Uli_Cut_Run_Mouse_embryo_H3K4me3_rep1_1.fastq \
-2 Uli_Cut_Run_Mouse_embryo_H3K4me3_rep1_2.fastq > Uli_Cut_Run_Mouse_embryo_H3K4me3_rep1.sam

msg="Convert sam to bam"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools view -Sb Uli_Cut_Run_Mouse_embryo_H3K4me3_rep1.sam -@ "$THREADS" > Uli_Cut_Run_Mouse_embryo_H3K4me3_rep1.bam

msg="Sort bam with "$SAMTOOLS""; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools sort Uli_Cut_Run_Mouse_embryo_H3K4me3_rep1.bam -@ "$THREADS" > Uli_Cut_Run_Mouse_embryo_H3K4me3_rep1_sorted.bam

msg="Mark and remove duplicates with Picard and sort with "$SAMTOOLS""; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
gatk MarkDuplicates \
-I Uli_Cut_Run_Mouse_embryo_H3K4me3_rep1_sorted.bam \
-O Uli_Cut_Run_Mouse_embryo_H3K4me3_rep1_rmdup.bam \
-M Uli_Cut_Run_Mouse_embryo_H3K4me3_rep1_sorted_metrics.txt \
--ASSUME_SORTED true --REMOVE_DUPLICATES true

samtools sort Uli_Cut_Run_Mouse_embryo_H3K4me3_rep1_rmdup.bam -@ "$THREADS" > Uli_Cut_Run_Mouse_embryo_H3K4me3_rep1.bam
samtools index Uli_Cut_Run_Mouse_embryo_H3K4me3_rep1.bam -@ "$THREADS"

msg="Remove old and unsorted bam and sam, and preprocessing beds"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
rm Uli_Cut_Run_Mouse_embryo_H3K4me3_rep1_unsort.bam
rm Uli_Cut_Run_Mouse_embryo_H3K4me3_rep1.sam
rm Uli_Cut_Run_Mouse_embryo_H3K4me3_rep1_sorted.bam
rm Uli_Cut_Run_Mouse_embryo_H3K4me3_rep1_rmdup.bam

msg="Gzip fastq"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
gzip Uli_Cut_Run_Mouse_embryo_H3K4me3_rep1_1.fastq
gzip Uli_Cut_Run_Mouse_embryo_H3K4me3_rep1_2.fastq

msg="Call peaks with MACS2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
macs2 callpeak -t Uli_Cut_Run_Mouse_embryo_H3K4me3_rep1.bam -f AUTO -n Uli_Cut_Run_Mouse_embryo_H3K4me3_rep1 \
--nomodel -p 5e-2 --min-length 500 --max-gap 400 --SPMR --call-summits -B

msg="Map with bowtie2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
bowtie2 -p "$THREADS" --end-to-end --very-sensitive --no-mixed --no-discordant \
--met-file 'Uli_Cut_Run_Mouse_embryo_H3K4me3_rep2.txt' \
-q --phred33 -I 10 -X 700 -x /scratch/users/astar/gis/muliaditand/programmes/bowtie2/GRCm38/GRCm38 \
-1 Uli_Cut_Run_Mouse_embryo_H3K4me3_rep2_1.fastq \
-2 Uli_Cut_Run_Mouse_embryo_H3K4me3_rep2_2.fastq > Uli_Cut_Run_Mouse_embryo_H3K4me3_rep2.sam

msg="Convert sam to bam"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools view -Sb Uli_Cut_Run_Mouse_embryo_H3K4me3_rep2.sam -@ "$THREADS" > Uli_Cut_Run_Mouse_embryo_H3K4me3_rep2.bam

msg="Sort bam with "$SAMTOOLS""; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools sort Uli_Cut_Run_Mouse_embryo_H3K4me3_rep2.bam -@ "$THREADS" > Uli_Cut_Run_Mouse_embryo_H3K4me3_rep2_sorted.bam

msg="Mark and remove duplicates with Picard and sort with "$SAMTOOLS""; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
gatk MarkDuplicates \
-I Uli_Cut_Run_Mouse_embryo_H3K4me3_rep2_sorted.bam \
-O Uli_Cut_Run_Mouse_embryo_H3K4me3_rep2_rmdup.bam \
-M Uli_Cut_Run_Mouse_embryo_H3K4me3_rep2_sorted_metrics.txt \
--ASSUME_SORTED true --REMOVE_DUPLICATES true

samtools sort Uli_Cut_Run_Mouse_embryo_H3K4me3_rep2_rmdup.bam -@ "$THREADS" > Uli_Cut_Run_Mouse_embryo_H3K4me3_rep2.bam
samtools index Uli_Cut_Run_Mouse_embryo_H3K4me3_rep2.bam -@ "$THREADS"

msg="Remove old and unsorted bam and sam, and preprocessing beds"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
rm Uli_Cut_Run_Mouse_embryo_H3K4me3_rep2_unsort.bam
rm Uli_Cut_Run_Mouse_embryo_H3K4me3_rep2.sam
rm Uli_Cut_Run_Mouse_embryo_H3K4me3_rep2_sorted.bam
rm Uli_Cut_Run_Mouse_embryo_H3K4me3_rep2_rmdup.bam

msg="Gzip fastq"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
gzip Uli_Cut_Run_Mouse_embryo_H3K4me3_rep2_1.fastq
gzip Uli_Cut_Run_Mouse_embryo_H3K4me3_rep2_2.fastq

msg="Call peaks with MACS2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
macs2 callpeak -t Uli_Cut_Run_Mouse_embryo_H3K4me3_rep2.bam -f AUTO -n Uli_Cut_Run_Mouse_embryo_H3K4me3_rep2 \
--nomodel -p 5e-2 --min-length 500 --max-gap 400 --SPMR --call-summits -B

