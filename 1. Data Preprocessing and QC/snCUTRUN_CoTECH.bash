#!/usr/bin/env bash

#SBATCH -J snCUTRUN_CoTECH_bulk
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -p normal
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 20
#SBATCH --mem-per-cpu 3500
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

#Set variables and directories
longLine="--------------------"
THREADS=13

mamba activate py3
module load samtools
module load gatk

cd /scratch/users/astar/gis/muliaditand/sncutrun/public_data/cotech/

for i in cotech_mouse_NIH3T3_H3K4me3_rep1 cotech_mouse_NIH3T3_H3K4me3_rep2

do 

msg="Map with bowtie2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
bowtie2 -p "$THREADS" --end-to-end --very-sensitive --no-mixed --no-discordant \
--met-file "$i"'.txt' \
-q --phred33 -I 10 -X 700 -x /scratch/users/astar/gis/muliaditand/programmes/bowtie2/GRCm38/GRCm38 \
-1 "$i"'_1.fastq' \
-2 "$i"'_2.fastq' > "$i"'.sam'

msg="Convert sam to bam"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools view -Sb  "$i"'.sam' -@ "$THREADS" >  "$i"'_unsort.bam'

msg="Sort bam with "$SAMTOOLS""; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools sort "$i"'_unsort.bam' -@ "$THREADS" > "$i"'_sorted.bam'

msg="Mark and remove duplicates with Picard and sort with "$SAMTOOLS""; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
gatk MarkDuplicates \
-I "$i"'_sorted.bam' \
-O "$i"'_rmdup.bam' \
-M "$i"'_sorted_metrics.txt' \
--ASSUME_SORTED true --REMOVE_DUPLICATES true

samtools sort "$i"'_rmdup.bam' -@ "$THREADS" > "$i"'.bam'
samtools index "$i"'.bam' -@ "$THREADS"

msg="Remove old and unsorted bam and sam, and preprocessing beds"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
rm "$i"'.sam'
rm "$i"'_unsort.bam'
rm "$i"'_sorted.bam'
rm "$i"'_rmdup.bam'

msg="Gzip fastq"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
gzip "$i"'_1.fastq'
gzip "$i"'_2.fastq'

msg="Call peaks with MACS2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
macs2 callpeak -t"$i"'.bam' -f AUTO -n "$i" \
--nomodel -p 5e-2 --min-length 500 --max-gap 400 --SPMR --call-summits -B

done

for i in cotech_human_k562_H3K4me3_nopreamplification_rep1 cotech_human_k562_H3K4me3_nopreamplification_rep2 cotech_human_k562_H3K4me3_nopreamplification_rep3 cotech_human_k562_H3K4me3_nopreamplification_rep4 cotech_human_k562_H3K4me3_nopreamplification_rep5 cotech_human_k562_H3K4me3_nopreamplification_rep6 cotech_human_k562_H3K4me3_nopreamplification_rep7 cotech_human_k562_H3K4me3_nopreamplification_rep8 cotech_human_k562_H3K4me3_nopreamplification_rep9 cotech_human_k562_H3K4me3_nopreamplification_rep10 cotech_human_k562_H3K4me3_nopreamplification_rep11 cotech_human_k562_H3K4me3_nopreamplification_rep2cotech_human_k562_H3K4me3_nopreamplification_rep12 cotech_human_k562_H3K4me3_nopreamplification_rep13 cotech_human_k562_H3K4me3_nopreamplification_rep14 cotech_human_k562_H3K4me3_nopreamplification_rep15 cotech_human_k562_H3K4me3_preamplification_rep1 cotech_human_k562_H3K4me3_preamplification_rep2 cotech_human_k562_H3K4me3_preamplification_rep3 cotech_human_k562_H3K4me3_preamplification_rep4 cotech_human_k562_H3K4me3_preamplification_rep5 cotech_human_k562_H3K4me3_preamplification_rep6 cotech_human_k562_H3K4me3_preamplification_rep7 cotech_human_k562_H3K4me3_preamplification_rep8 cotech_human_k562_H3K4me3_preamplification_rep9 cotech_human_k562_H3K4me3_preamplification_rep10 cotech_human_k562_H3K4me3_preamplification_rep11 cotech_human_k562_H3K4me3_preamplification_rep12 cotech_human_k562_H3K4me3_preamplification_rep12 cotech_human_k562_H3K4me3_preamplification_rep13 cotech_human_k562_H3K4me3_preamplification_rep14 cotech_human_k562_H3K4me3_preamplification_rep15  

do 

msg="Map with bowtie2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
bowtie2 -p "$THREADS" --end-to-end --very-sensitive --no-mixed --no-discordant \
--met-file "$i"'.txt' \
-q --phred33 -I 10 -X 700 -x /scratch/users/astar/gis/muliaditand/programmes/bowtie2/GRCh38_noalt_decoy_as/GRCh38_noalt_decoy_as \
-1 "$i"'_1.fastq' \
-2 "$i"'_2.fastq' > "$i"'.sam'

msg="Convert sam to bam"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools view -Sb  "$i"'.sam' -@ "$THREADS" >  "$i"'_unsort.bam'

msg="Sort bam with "$SAMTOOLS""; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools sort "$i"'_unsort.bam' -@ "$THREADS" > "$i"'_sorted.bam'

msg="Mark and remove duplicates with Picard and sort with "$SAMTOOLS""; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
gatk MarkDuplicates \
-I "$i"'_sorted.bam' \
-O "$i"'_rmdup.bam' \
-M "$i"'_sorted_metrics.txt' \
--ASSUME_SORTED true --REMOVE_DUPLICATES true

samtools sort "$i"'_rmdup.bam' -@ "$THREADS" > "$i"'.bam'
samtools index "$i"'.bam' -@ "$THREADS"

msg="Remove old and unsorted bam and sam, and preprocessing beds"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
rm "$i"'.sam'
rm "$i"'_unsort.bam'
rm "$i"'_sorted.bam'
rm "$i"'_rmdup.bam'

msg="Gzip fastq"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
gzip "$i"'_1.fastq'
gzip "$i"'_2.fastq'

msg="Call peaks with MACS2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
macs2 callpeak -t"$i"'.bam' -f AUTO -n "$i" \
--nomodel -p 5e-2 --min-length 500 --max-gap 400 --SPMR --call-summits -B

done

