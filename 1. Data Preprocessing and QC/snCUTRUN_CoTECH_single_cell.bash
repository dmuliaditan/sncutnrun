#!/usr/bin/env bash

#SBATCH -J snCUTRUN_sccuttag_cotech_single_cell
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

cd /scratch/users/astar/gis/muliaditand/sncutrun/public_data/cotech/mESC_H3K4me3_H3K27me3_single_cells/

for i in `cat '/filenames.txt'`

do 
msg="Cell: $i"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
msg="Map with bowtie2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
(bowtie2 -p "$THREADS" --end-to-end --very-sensitive --no-mixed --no-discordant \
--met-file "$i"'.txt' \
-q --phred33 -I 10 -X 700 -x /scratch/users/astar/gis/muliaditand/programmes/bowtie2/GRCm38/GRCm38 \
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

#Aggregate .bam filenames 
msg="Aggregate individual single-cell .bam names to cell line-antibody .bam list"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
ls -d *.bam > 'bamslist.txt'

msg="Merge all .bams to form pseudobulk, sort and index"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools merge -f -b 'bamslist.txt' 'mESC_H3K4me3_H3K27me3_merged_unsort.bam' --threads "$THREADS"
samtools sort 'mESC_H3K4me3_H3K27me3_merged_unsort.bam' > 'mESC_H3K4me3_H3K27me3_merged_sorted.bam'
samtools index 'mESC_H3K4me3_H3K27me3_merged_sorted.bam'

#Call peaks on each merged bam using MACS2
msg="Call peaks on each merged bam using MACS2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
macs2 callpeak -t 'mESC_H3K4me3_H3K27me3_merged_sorted.bam' -f BAMPE \
-n mESC_H3K4me3_H3K27me3 --seed 1234 -g mm -f AUTO --nomodel -B -p 5e-2 --min-length 500 --max-gap 400 --SPMR --call-summits

msg="Finished"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"

done

cd /scratch/users/astar/gis/muliaditand/sncutrun/public_data/cotech/AGM_H3K27ac_single_cells/

for i in `cat '/filenames.txt'`

do 
msg="Cell: $i"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
msg="Map with bowtie2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
(bowtie2 -p "$THREADS" --end-to-end --very-sensitive --no-mixed --no-discordant \
--met-file "$i"'.txt' \
-q --phred33 -I 10 -X 700 -x /scratch/users/astar/gis/muliaditand/programmes/bowtie2/GRCm38/GRCm38 \
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

#Aggregate .bam filenames 
msg="Aggregate individual single-cell .bam names to cell line-antibody .bam list"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
ls -d *.bam > 'bamslist.txt'

msg="Merge all .bams to form pseudobulk, sort and index"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools merge -f -b 'bamslist.txt' 'AGM_H3K27ac_merged_unsort.bam' --threads "$THREADS"
samtools sort 'AGM_H3K27ac_merged_unsort.bam' > 'AGM_H3K27ac_merged_sorted.bam'
samtools index 'AGM_H3K27ac_merged_sorted.bam'

#Call peaks on each merged bam using MACS2
msg="Call peaks on each merged bam using MACS2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
macs2 callpeak -t 'AGM_H3K27ac_merged_sorted.bam' -f BAMPE \
-n AGM_H3K27ac --seed 1234 -g mm -f AUTO --nomodel -B -p 5e-2 --min-length 500 --max-gap 400 --SPMR --call-summits

done

cd /scratch/users/astar/gis/muliaditand/sncutrun/public_data/cotech/YS_H3K27ac_single_cells/

for i in `cat '/filenames.txt'`

do 
msg="Cell: $i"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
msg="Map with bowtie2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
(bowtie2 -p "$THREADS" --end-to-end --very-sensitive --no-mixed --no-discordant \
--met-file "$i"'.txt' \
-q --phred33 -I 10 -X 700 -x /scratch/users/astar/gis/muliaditand/programmes/bowtie2/GRCm38/GRCm38 \
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

#Aggregate .bam filenames 
msg="Aggregate individual single-cell .bam names to cell line-antibody .bam list"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
ls -d *.bam > 'bamslist.txt'

msg="Merge all .bams to form pseudobulk, sort and index"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
samtools merge -f -b 'bamslist.txt' 'YS_H3K27ac_merged_unsort.bam' --threads "$THREADS"
samtools sort 'YS_H3K27ac_merged_unsort.bam' > 'YS_H3K27ac_merged_sorted.bam'
samtools index 'YS_H3K27ac_merged_sorted.bam'

#Call peaks on each merged bam using MACS2
msg="Call peaks on each merged bam using MACS2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
macs2 callpeak -t 'YS_H3K27ac_merged_sorted.bam' -f BAMPE \
-n YS_H3K27ac --seed 1234 -g mm -f AUTO --nomodel -B -p 5e-2 --min-length 500 --max-gap 400 --SPMR --call-summits

done
