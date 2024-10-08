#!/usr/bin/env bash

#SBATCH -J 20241008_benchmarking_90_cells_percentage_H3K4me3_in_promoters
#SBATCH -t 02:00:00
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
touch 'snCUTRUN_90_cells_reads_in_promoters.txt'
for i in `cat 'filenames.txt'`

do
#Compile barcode file per cell and append to master barcode file
msg="Compile barcode file per cell and append to master barcode file"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
echo "$i" > "$i"'_barcode.txt'

#1. Calculate Unique Mapped Reads (UMRs)
msg="Calculate Unique Mapped Reads (UMRs)"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools view -c -f 1 -F 12 "$i"'.bam' > "$i"'_raw_reads.txt'

#Convert single-cell bams to beds
bedtools bamtobed -i "$i"'.bam' -bedpe \
| awk '{if($2!="-1") print}' | sort -k1,1 -k2,2n | cut -f1,2,6,7 > "$i"'_preproc.bed'
bedtools sort -i "$i"'_preproc.bed' -faidx /scratch/users/astar/gis/muliaditand/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai > "$i"'_postsort.bed'
nrow=$(wc -l < "$i"'_postsort.bed' )
for ((c=1; c<=nrow; c++)) ; do echo "$i" ; done > "$i"'.txt'
awk 'FNR==NR{a[FNR]=$1;next};{$NF=a[FNR]};1' "$i"'.txt' "$i"'_postsort.bed' > "$i"'.bed'
perl -p -i -e 's/ /\t/g' "$i"'.bed'

#Intersect bed files with bedtools
bedtools intersect -u -a "$i"'.bed' -b promoters_hg38.bed | wc -l > "$i"'_reads_in_promoters.txt'

paste "$i"'_barcode.txt' "$i"'_raw_reads.txt' "$i"'_reads_in_promoters.txt' > "$i"'raw_reads_in_promoter.txt'
cat "$i"'raw_reads_in_promoter.txt' >> 'snCUTRUN_90_cells_reads_in_promoters.txt'
done
