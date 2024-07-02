#!/usr/bin/env bash

#SBATCH -J sccuttag_calculate_mapping_rate
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

cd /scratch/users/astar/gis/muliaditand/sncutrun/public_data/sccuttag_kaya_okur_2019/K562_H3K4me2_single_cells/

touch 'sccuttag_K562_H3K4me2_mapping_metrics.txt'

for i in `cat 'filenames.txt'`

do 
msg="Cell: $i"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"

gzip -d "$i"'_1.fastq.gz'
gzip -d "$i"'_2.fastq.gz'

msg="Count number of fastq reads"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
echo $(cat "$i"'_1.fastq'|wc -l)/4|bc > "$i"'_fastq_read1_count.txt'
echo $(cat "$i"'_2.fastq'|wc -l)/4|bc > "$i"'_fastq_read2_count.txt'

msg="Gzip fastq"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
gzip "$i"'_1.fastq'
gzip "$i"'_2.fastq'

#Compile barcode file per cell and append to master barcode file
msg="Compile barcode file per cell and append to master barcode file"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
echo "$i" > "$i"'_barcode.txt'

#1. Calculate Unique Mapped Reads (UMRs)
samtools view -c -f 1 -F 12 "$i"'.bam' > "$i"'.txt'

paste "$i"'_barcode.txt' "$i"'.txt' "$i"'_fastq_read1_count.txt' "$i"'_fastq_read2_count.txt' > "$i"'_barcodes_precount.txt'
cat "$i"'_barcodes_precount.txt' >> 'sccuttag_K562_H3K4me2_mapping_metrics.txt'

rm "$i"'.txt'
rm "$i"'.bed'
rm "$i"'_fastq_read1_count.txt'
rm "$i"'_fastq_read2_count.txt'
rm "$i"'_barcodes_precount.txt'
rm "$i"'_sorted_metrics.txt'
rm "$i"'_barcode.txt'

done

cd /scratch/users/astar/gis/muliaditand/sncutrun/public_data/sccuttag_kaya_okur_2019/K562_H3K27me3_single_cells/

touch 'sccuttag_K562_H3K27me3_mapping_metrics.txt'

for i in `cat 'filenames.txt'`

do 
msg="Cell: $i"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"

gzip -d "$i"'_1.fastq.gz'
gzip -d "$i"'_2.fastq.gz'

msg="Count number of fastq reads"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
echo $(cat "$i"'_1.fastq'|wc -l)/4|bc > "$i"'_fastq_read1_count.txt'
echo $(cat "$i"'_2.fastq'|wc -l)/4|bc > "$i"'_fastq_read2_count.txt'

msg="Gzip fastq"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
gzip "$i"'_1.fastq'
gzip "$i"'_2.fastq'

#Compile barcode file per cell and append to master barcode file
msg="Compile barcode file per cell and append to master barcode file"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
echo "$i" > "$i"'_barcode.txt'

#1. Calculate Unique Mapped Reads (UMRs)
samtools view -c -f 1 -F 12 "$i"'.bam' > "$i"'.txt'

paste "$i"'_barcode.txt' "$i"'.txt' "$i"'_fastq_read1_count.txt' "$i"'_fastq_read2_count.txt' > "$i"'_barcodes_precount.txt'
cat "$i"'_barcodes_precount.txt' >> 'sccuttag_K562_H3K27me3_mapping_metrics.txt'

rm "$i"'.txt'
rm "$i"'.bed'
rm "$i"'_fastq_read1_count.txt'
rm "$i"'_fastq_read2_count.txt'
rm "$i"'_barcodes_precount.txt'
rm "$i"'_sorted_metrics.txt'
rm "$i"'_barcode.txt'

done

cd /scratch/users/astar/gis/muliaditand/sncutrun/public_data/sccuttag_kaya_okur_2019/H1_K562_H3K27me3_pool2_single_cells/

touch 'sccuttag_H1_K562_H3K27me3_pool2_mapping_metrics.txt'

for i in `cat 'filenames.txt'`

do 
msg="Cell: $i"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"

gzip -d "$i"'_1.fastq.gz'
gzip -d "$i"'_2.fastq.gz'

msg="Count number of fastq reads"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
echo $(cat "$i"'_1.fastq'|wc -l)/4|bc > "$i"'_fastq_read1_count.txt'
echo $(cat "$i"'_2.fastq'|wc -l)/4|bc > "$i"'_fastq_read2_count.txt'

msg="Gzip fastq"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
gzip "$i"'_1.fastq'
gzip "$i"'_2.fastq'

#Compile barcode file per cell and append to master barcode file
msg="Compile barcode file per cell and append to master barcode file"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
echo "$i" > "$i"'_barcode.txt'

#1. Calculate Unique Mapped Reads (UMRs)
samtools view -c -f 1 -F 12 "$i"'.bam' > "$i"'.txt'

paste "$i"'_barcode.txt' "$i"'.txt' "$i"'_fastq_read1_count.txt' "$i"'_fastq_read2_count.txt' > "$i"'_barcodes_precount.txt'
cat "$i"'_barcodes_precount.txt' >> 'sccuttag_H1_K562_H3K27me3_pool2_mapping_metrics.txt'

rm "$i"'.txt'
rm "$i"'.bed'
rm "$i"'_fastq_read1_count.txt'
rm "$i"'_fastq_read2_count.txt'
rm "$i"'_barcodes_precount.txt'
rm "$i"'_sorted_metrics.txt'
rm "$i"'_barcode.txt'

done
