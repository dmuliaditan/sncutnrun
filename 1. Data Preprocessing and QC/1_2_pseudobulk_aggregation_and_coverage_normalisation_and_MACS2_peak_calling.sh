#!/usr/bin/env bash

#Set variables and directories
#Here, after pre-processing all single-cell data to .bam format, we will merge the single cell .bams to form a pseudobulk to call peaks

#Set variables
longLine="--------------------"
SAMTOOLS='samtools'
MACS2='macs2'
AGGR_RESULTS_DIR='/mnt/d/snCUT_RUN/data'
THREADS=8

mkdir -p "$AGGR_RESULTS_DIR"'/tmp_dir'
TMP_DIR="$AGGR_RESULTS_DIR"'/tmp_dir'

for p in 120pri_K4me3 120pri_k27ac etc
do
#First, we merge the single cell bams to pseudobulk
msg="Merge all .bams to form pseudobulk, sort and index"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
"$SAMTOOLS" merge -f -b "$AGGR_RESULTS_DIR"'/bamlists/'"$p"'_bamslist.txt' "$AGGR_RESULTS_DIR"'/final_bams/pseudobulk/'"$p"'_unsort.bam' --threads "$THREADS"
"$SAMTOOLS" sort "$AGGR_RESULTS_DIR"'/final_bams/pseudobulk/'"$p"'_unsort.bam' > "$AGGR_RESULTS_DIR"'/final_bams/pseudobulk/'"$p"'_merged_sorted.bam'
"$SAMTOOLS" index "$AGGR_RESULTS_DIR"'/final_bams/pseudobulk/'"$p"'_merged_sorted.bam'

#Next we normalised based on coverage to the sample with lowest coverage
#Reference: https://doi.org/10.1016/j.ymeth.2020.03.005	Nakato (2020)
#"Simple total read normalization is commonly used, which scales the sample read number to be the same. 
#The underlying assumption is that the difference in mapped reads among samples is sufficiently smaller than the total read number."

#First count number of reads in the pseudobulk
echo "$p"
"$SAMTOOLS" view -c -f 1 -F 12 "$AGGR_RESULTS_DIR"'/final_bams/pseudobulk/'"$p"'_merged_sorted.bam'
  
done

#Example: if highest .bam has read number of 100 and lowest 75
#Next manually normalise to the .bam with lowest coverage e.g.:
"$SAMTOOLS" view -bs 45.75 "$AGGR_RESULTS_DIR"'/final_bams/pseudobulk/'"$p"'_merged_sorted.bam' \
> "$AGGR_RESULTS_DIR"'/final_bams/pseudobulk/'"$p"'_norm.bam'
"$SAMTOOLS" sort "$AGGR_RESULTS_DIR"'/final_bams/pseudobulk/'"$p"'_norm.bam' > "$AGGR_RESULTS_DIR"'/final_bams/pseudobulk/'"$p"'_norm_sort.bam'
"$SAMTOOLS" index "$AGGR_RESULTS_DIR"'/final_bams/pseudobulk/'"$p"'_norm_sort.bam'

#Next call peaks with MACS2:

for p in 120pri_K4me3 120pri_k27ac etc
do

#Call peaks on each merged bam using MACS2
msg="Call peaks on each merged bam using "$MACS2""; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
"$MACS2" callpeak -t "$AGGR_RESULTS_DIR"'/final_bams/pseudobulk/'"$p"'_norm_sort.bam' --outdir "$AGGR_RESULTS_DIR"'/macs2' \
-n "$p"'_norm_sort' --seed 1234 -g hs -f BAMPE --nomodel -B -p 5e-2 --min-length 500 --max-gap 400 --SPMR --call-summits
  
done

