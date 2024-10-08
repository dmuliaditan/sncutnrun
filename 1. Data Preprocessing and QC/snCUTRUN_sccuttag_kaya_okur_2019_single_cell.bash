#!/usr/bin/env bash

#SBATCH -J snCUTRUN_sccuttag_kaya_okur_2019_single_cell
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

cd /scratch/users/astar/gis/muliaditand/sncutrun/public_data/sccuttag_kaya_okur_2019/K562_H3K4me2_single_cells/

for i in `cat 'filenames.txt'`

do 
msg="Cell: $i"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
msg="Map with bowtie2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
(bowtie2 -p "$THREADS" --end-to-end --very-sensitive --no-mixed --no-discordant \
--met-file "$i"'.txt' \
-q --phred33 -I 10 -X 700 -x /scratch/users/astar/gis/muliaditand/programmes/bowtie2/GRCh38_noalt_decoy_as/GRCh38_noalt_decoy_as \
-1 "$i"'_1.fastq' \
-2 "$i"'_2.fastq' -S "$i"'.sam') 2>>alignment_metrics.txt 

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

done

#Aggregate .bam filenames 
msg="Aggregate individual single-cell .bam names to cell line-antibody .bam list"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
ls -d *.bam > 'bamslist.txt'

msg="Merge all .bams to form pseudobulk, sort and index"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools merge -f -b 'bamslist.txt' 'K562_H3K4me2_merged_unsort.bam' --threads "$THREADS"
samtools sort 'K562_H3K4me2_merged_unsort.bam' > 'K562_H3K4me2_merged_sorted.bam'
samtools index 'K562_H3K4me2_merged_sorted.bam'

#Call peaks on each merged bam using MACS2
msg="Call peaks on each merged bam using MACS2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
macs2 callpeak -t 'K562_H3K4me2_merged_sorted.bam' -f BAMPE \
-n K562_H3K4me2 --seed 1234 -g hs -f AUTO --nomodel -B -p 5e-2 --min-length 500 --max-gap 400 --SPMR --call-summits

msg="Finished"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"

cd /scratch/users/astar/gis/muliaditand/sncutrun/public_data/sccuttag_kaya_okur_2019/K562_H3K27me3_single_cells/

for i in `cat 'filenames.txt'`

do 
msg="Cell: $i"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
msg="Map with bowtie2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
(bowtie2 -p "$THREADS" --end-to-end --very-sensitive --no-mixed --no-discordant \
--met-file "$i"'.txt' \
-q --phred33 -I 10 -X 700 -x /scratch/users/astar/gis/muliaditand/programmes/bowtie2/GRCh38_noalt_decoy_as/GRCh38_noalt_decoy_as \
-1 "$i"'_1.fastq' \
-2 "$i"'_2.fastq' -S "$i"'.sam') 2>>alignment_metrics.txt 

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

done

#Aggregate .bam filenames 
msg="Aggregate individual single-cell .bam names to cell line-antibody .bam list"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
ls -d *.bam > 'bamslist.txt'

msg="Merge all .bams to form pseudobulk, sort and index"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools merge -f -b 'bamslist.txt' 'K562_H3K27me3_merged_unsort.bam' --threads "$THREADS"
samtools sort 'K562_H3K27me3_merged_unsort.bam' > 'K562_H3K27me3_merged_sorted.bam'
samtools index 'K562_H3K27me3_merged_sorted.bam'

#Call peaks on each merged bam using MACS2
msg="Call peaks on each merged bam using MACS2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
macs2 callpeak -t 'K562_H3K27me3_merged_sorted.bam' -f BAMPE \
-n K562_H3K27me3 --seed 1234 -g hs -f AUTO --nomodel -B -p 5e-2 --min-length 500 --max-gap 400 --SPMR --call-summits

msg="Finished"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"

cd /scratch/users/astar/gis/muliaditand/sncutrun/public_data/sccuttag_kaya_okur_2019/H1_K562_H3K27me3_pool2_single_cells/

for i in `cat 'filenames.txt'`

do 
msg="Cell: $i"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
msg="Map with bowtie2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
(bowtie2 -p "$THREADS" --end-to-end --very-sensitive --no-mixed --no-discordant \
--met-file "$i"'.txt' \
-q --phred33 -I 10 -X 700 -x /scratch/users/astar/gis/muliaditand/programmes/bowtie2/GRCh38_noalt_decoy_as/GRCh38_noalt_decoy_as \
-1 "$i"'_1.fastq' \
-2 "$i"'_2.fastq' -S "$i"'.sam') 2>>alignment_metrics.txt 

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

done

#Aggregate .bam filenames 
msg="Aggregate individual single-cell .bam names to cell line-antibody .bam list"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
ls -d *.bam > 'bamslist.txt'

msg="Merge all .bams to form pseudobulk, sort and index"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools merge -f -b 'bamslist.txt' 'H1_K562_H3K27me3_pool2_merged_unsort.bam' --threads "$THREADS"
samtools sort 'H1_K562_H3K27me3_pool2_merged_unsort.bam' > 'H1_K562_H3K27me3_pool2_merged_sorted.bam'
samtools index 'H1_K562_H3K27me3_pool2_merged_sorted.bam'

#Call peaks on each merged bam using MACS2
msg="Call peaks on each merged bam using MACS2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
macs2 callpeak -t 'H1_K562_H3K27me3_merged_sorted.bam' -f BAMPE \
-n H1_K562_H3K27me3_pool2 --seed 1234 -g hs -f AUTO --nomodel -B -p 5e-2 --min-length 500 --max-gap 400 --SPMR --call-summits

