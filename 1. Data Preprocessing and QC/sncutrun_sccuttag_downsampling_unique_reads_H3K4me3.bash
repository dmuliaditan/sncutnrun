#!/usr/bin/env bash

#SBATCH -J after_downsampling_snCUTRUN_calculate_mapping_rate_H3K4me3
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

mamba init
mamba activate py3
module load samtools
module load gatk
module load bedtools2

cd /scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/H3K4me3

touch 'after_downsampling_snCUTRUN_H3K4me3_mapping_metrics.txt'

for i in `cat 'filenames.txt'`

do 
msg="Cell: $i"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"

gzip -d "$i"'_L001_R1_001.fastq.gz'
gzip -d "$i"'_L001_R2_001.fastq.gz'

msg="Downsample cell: $i"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
/scratch/users/astar/gis/muliaditand/programmes/seqtk/seqtk/seqtk sample -s100 "$i"'_L001_R1_001.fastq' 0.75 >  "$i"'_L001_R1_001_downsampled.fastq' 
/scratch/users/astar/gis/muliaditand/programmes/seqtk/seqtk/seqtk sample -s100 "$i"'_L001_R2_001.fastq' 0.75 >  "$i"'_L001_R2_001_downsampled.fastq' 

msg="Map with bowtie2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
(bowtie2 -p "$THREADS" --end-to-end --very-sensitive --no-mixed --no-discordant \
--met-file "$i"'_downsampled.txt' \
-q --phred33 -I 10 -X 700 -x /scratch/users/astar/gis/muliaditand/programmes/bowtie2/GRCh38_noalt_decoy_as/GRCh38_noalt_decoy_as \
-1 "$i"'_L001_R1_001_downsampled.fastq' \
-2 "$i"'_L001_R2_001_downsampled.fastq' -S "$i"'_downsampled.sam') 2>>alignment_metrics_downsampled.txt 

msg="Convert sam to bam"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools view -Sb  "$i"'_downsampled.sam' -@ "$THREADS" >  "$i"'_downsampled_unsort.bam'

msg="Sort bam with "$SAMTOOLS""; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools sort "$i"'_downsampled_unsort.bam' -@ "$THREADS" > "$i"'_downsampled_sorted.bam'

msg="Mark and remove duplicates with Picard and sort with "$SAMTOOLS""; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
gatk MarkDuplicates \
-I "$i"'_downsampled_sorted.bam' \
-O "$i"'_downsampled_rmdup.bam' \
-M "$i"'_downsampled_sorted_metrics.txt' \
--ASSUME_SORTED true --REMOVE_DUPLICATES true

samtools sort "$i"'_downsampled_rmdup.bam' -@ "$THREADS" > "$i"'_downsampled.bam'
samtools index "$i"'_downsampled.bam' -@ "$THREADS"

msg="Remove old and unsorted bam and sam, and preprocessing beds"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
rm "$i"'_downsampled.sam'
rm "$i"'_downsampled_unsort.bam'
rm "$i"'_downsampled_sorted.bam'
rm "$i"'_downsampled_rmdup.bam'

msg="Count number of fastq reads"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
echo $(cat "$i"'_L001_R1_001_downsampled.fastq'|wc -l)/4|bc > "$i"'_fastq_read1_count_downsampled.txt'
echo $(cat "$i"'_L001_R2_001_downsampled.fastq'|wc -l)/4|bc > "$i"'_fastq_read2_count_downsampled.txt'

msg="Gzip fastq"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
gzip "$i"'_L001_R1_001.fastq'
gzip "$i"'_L001_R2_001.fastq'
gzip "$i"'_L001_R1_001_downsampled.fastq'
gzip "$i"'_L001_R2_001_downsampled.fastq'

#Compile barcode file per cell and append to master barcode file
msg="Compile barcode file per cell and append to master barcode file"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
echo "$i" > "$i"'_barcode.txt'

#1. Calculate Unique Mapped Reads (UMRs)
samtools view -c -f 1 -F 12 "$i"'_downsampled.bam' > "$i"'_downsampled.txt'

paste "$i"'_barcode.txt' "$i"'_downsampled.txt' "$i"'_fastq_read1_count_downsampled.txt' "$i"'_fastq_read2_count_downsampled.txt' > "$i"'_barcodes_precount_downsampled.txt'
cat "$i"'_barcodes_precount_downsampled.txt' >> 'after_downsampling_snCUTRUN_H3K4me3_mapping_metrics.txt'

rm "$i"'_downsampled.txt'
rm "$i"'_fastq_read1_count_downsampled.txt'
rm "$i"'_fastq_read2_count_downsampled.txt'
rm "$i"'_barcodes_precount_downsampled.txt'
rm "$i"'_sorted_metrics_downsampled.txt'
rm "$i"'_barcode.txt'

done

