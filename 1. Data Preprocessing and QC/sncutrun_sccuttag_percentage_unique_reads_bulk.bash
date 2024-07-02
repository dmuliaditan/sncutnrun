#!/usr/bin/env bash

#SBATCH -J snCUTRUN_sccuttag_calculate_mapping_rate_bulk
#SBATCH -t 72:00:00
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

mamba init
mamba activate py3
module load samtools
module load gatk
module load bedtools2

cd /scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/bulk

touch 'snCUTRUN_sccuttag_mapping_metrics_bulk.txt'

for i in `cat 'filenames.txt'`

do 
msg="Sample: $i"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"

gzip -d "$i"'_R1_001.fastq.gz'
gzip -d "$i"'_R2_001.fastq.gz'

msg="Count number of fastq reads"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
echo $(cat "$i"'_R1_001.fastq'|wc -l)/4|bc > "$i"'_fastq_read1_count.txt'
echo $(cat "$i"'_R2_001.fastq'|wc -l)/4|bc > "$i"'_fastq_read2_count.txt'

msg="Map with bowtie2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
(bowtie2 -p "$THREADS" --end-to-end --very-sensitive --no-mixed --no-discordant \
--met-file "$i"'.txt' \
-q --phred33 -I 10 -X 700 -x /scratch/users/astar/gis/muliaditand/programmes/bowtie2/GRCh38_noalt_decoy_as/GRCh38_noalt_decoy_as \
-1 "$i"'_R1_001.fastq' \
-2 "$i"'_R2_001.fastq' -S "$i"'.sam') 2>>alignment_metrics.txt 

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
gzip "$i"'_R1_001.fastq'
gzip "$i"'_R2_001.fastq'

#Compile barcode file per cell and append to master barcode file
msg="Compile barcode file per cell and append to master barcode file"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
echo "$i" > "$i"'_barcode.txt'

#1. Calculate Unique Mapped Reads (UMRs)
samtools view -c -f 1 -F 12 "$i"'.bam' > "$i"'.txt'

paste "$i"'_barcode.txt' "$i"'.txt' "$i"'_fastq_read1_count.txt' "$i"'_fastq_read2_count.txt' > "$i"'_barcodes_precount.txt'
cat "$i"'_barcodes_precount.txt' >> 'snCUTRUN_sccuttag_mapping_metrics_bulk.txt'

rm "$i"'.txt'
rm "$i"'_fastq_read1_count.txt'
rm "$i"'_fastq_read2_count.txt'
rm "$i"'_barcodes_precount.txt'
rm "$i"'_sorted_metrics.txt'
rm "$i"'_barcode.txt'

done

cd /scratch/users/astar/gis/muliaditand/sncutrun/public_data/sccuttag_kaya_okur_2019

for i in `cat 'filenames.txt'`

do 
msg="Sample: $i"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"

gzip -d "$i"'_1.fastq.gz'
gzip -d "$i"'_2.fastq.gz'

msg="Count number of fastq reads"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
echo $(cat "$i"'_1.fastq'|wc -l)/4|bc > "$i"'_fastq_read1_count.txt'
echo $(cat "$i"'_2.fastq'|wc -l)/4|bc > "$i"'_fastq_read2_count.txt'

#Compile barcode file per cell and append to master barcode file
msg="Compile barcode file per cell and append to master barcode file"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
echo "$i" > "$i"'_barcode.txt'

#1. Calculate Unique Mapped Reads (UMRs)
samtools view -c -f 1 -F 12 "$i"'.bam' > "$i"'.txt'

paste "$i"'_barcode.txt' "$i"'.txt' "$i"'_fastq_read1_count.txt' "$i"'_fastq_read2_count.txt' > "$i"'_barcodes_precount.txt'
cat "$i"'_barcodes_precount.txt' >> 'snCUTRUN_sccuttag_mapping_metrics_bulk.txt'

rm "$i"'.txt'
rm "$i"'_fastq_read1_count.txt'
rm "$i"'_fastq_read2_count.txt'
rm "$i"'_barcodes_precount.txt'
rm "$i"'_sorted_metrics.txt'
rm "$i"'_barcode.txt'

done
