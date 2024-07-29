#!/usr/bin/env bash

#SBATCH -J snCUTRUN_chromHMM_reanalysis
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

mamba activate py3
module load bedtools2

####Gather necessary data
#Call H3K27me3 peaks
cd /scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/bulk

#for i in HN120Pri HN120Met HN137Pri HN137Met
#do

#msg="Call peaks with MACS2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
#macs2 callpeak -t "$i"'_K27me3_rep1.bam' "$i"'_K27me3_rep2.bam' -c "$i"'_IgG.bam' -f AUTO -n "$i"'_K27me3' \
#--nomodel -p 5e-2 --min-length 500 --max-gap 400 --SPMR --call-summits -B

#done

#for i in HN120PCR HN137PCR
#do

#msg="Call peaks with MACS2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
#macs2 callpeak -t "$i"'_K27me3.bam' -c "$i"'_IgG.bam' -f AUTO -n "$i"'_K27me3' \
#--nomodel -p 5e-2 --min-length 500 --max-gap 400 --SPMR --call-summits -B

#done

#Intersect beds and sort
#bedtools intersect -a HN120Pri_K27me3_peaks.narrowPeak -b HN120Met_K27me3_peaks.narrowPeak HN120PCR_K27me3_peaks.narrowPeak \
#HN137Pri_K27me3_peaks.narrowPeak HN137Met_K27me3_peaks.narrowPeak HN137PCR_K27me3_peaks.narrowPeak -sorted > HN120_HN137_H3K27me3_shared_peaks.bed
#sort -k1,1V -k2,2n HN120_HN137_H3K27me3_shared_peaks.bed > HN120_HN137_H3K27me3_shared_peaks_sorted.bed

#Binarize the normalised IgG control .bams and normalised sample .bams with BinarizeBam to create the input of LearnModel
#java -Xms8g -Xmx12g -jar /scratch/users/astar/gis/muliaditand/programmes/chromHMM/ChromHMM/ChromHMM.jar BinarizeBed -center \
#/scratch/users/astar/gis/muliaditand/programmes/chromHMM/ChromHMM/CHROMSIZES/hg38.txt /scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/chromHMM/beds \
#/scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/chromHMM/chromHMM_cell_mark_table_shared_peaks.txt /scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/chromHMM 

#Run ChromHMM LearnModel on the binarized .beds
#Rerun at different number of states and compare biological significance of the results
#mkdir -p /scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/chromHMM/5_states
#java -Xms8g -Xmx12g -jar /scratch/users/astar/gis/muliaditand/programmes/chromHMM/ChromHMM/ChromHMM.jar LearnModel \
#-noautoopen -p 8 -init random /scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/chromHMM/binary_files \
#/scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/chromHMM/5_states 5 hg38

#Do separately on HN120 
#Binarize the normalised IgG control .bams and normalised sample .bams with BinarizeBam to create the input of LearnModel
#java -Xms8g -Xmx12g -jar /scratch/users/astar/gis/muliaditand/programmes/chromHMM/ChromHMM/ChromHMM.jar BinarizeBam \
#-paired -gzip -c /scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/chromHMM/HN120/controls \
#-o /scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/chromHMM/HN120/control_output \
#/scratch/users/astar/gis/muliaditand/programmes/chromHMM/ChromHMM/CHROMSIZES/hg38.txt \
#/scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/chromHMM/HN120/input_bams \
#/scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/chromHMM/HN120/chromHMM_cell_mark_table_HN120.txt \
#/scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/chromHMM/HN120/binary_files 

#Run ChromHMM LearnModel on the binarized .bams
#Rerun at different number of states and compare biological significance of the results
mkdir -p /scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/chromHMM/HN120/5_states
java -Xms8g -Xmx12g -jar /scratch/users/astar/gis/muliaditand/programmes/chromHMM/ChromHMM/ChromHMM.jar LearnModel \
-noautoopen -p 8 /scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/chromHMM/HN120/binary_files  \
/scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/chromHMM/HN120/5_states 5 hg38

#Do separately on HN137
#Binarize the normalised IgG control .bams and normalised sample .bams with BinarizeBam to create the input of LearnModel
#java -Xms8g -Xmx12g -jar /scratch/users/astar/gis/muliaditand/programmes/chromHMM/ChromHMM/ChromHMM.jar BinarizeBam \
#-paired -gzip \
#-c /scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/chromHMM/HN137/controls \
#-o /scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/chromHMM/HN137/control_output \
#/scratch/users/astar/gis/muliaditand/programmes/chromHMM/ChromHMM/CHROMSIZES/hg38.txt \
#/scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/chromHMM/HN137/input_bams \
#/scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/chromHMM/HN137/chromHMM_cell_mark_table_HN137.txt \
#/scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/chromHMM/HN137/binary_files 

#Run ChromHMM LearnModel on the binarized .bams
#Rerun at different number of states and compare biological significance of the results
mkdir -p /scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/chromHMM/HN137/5_states
java -Xms8g -Xmx12g -jar /scratch/users/astar/gis/muliaditand/programmes/chromHMM/ChromHMM/ChromHMM.jar LearnModel \
-noautoopen -p 8 /scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/chromHMM/HN137/binary_files  \
/scratch/users/astar/gis/muliaditand/sncutrun/sequence_runs/chromHMM/HN137/5_states 5 hg38
