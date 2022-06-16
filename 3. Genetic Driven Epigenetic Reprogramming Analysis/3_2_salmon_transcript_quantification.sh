#!/usr/bin/bash

#3.2. Pipeline to analyze snCUT&RUN bulk RNA-seq data
#Tutorial RNA-seq workflow: gene-level exploratory analysis and differential expression, Michael I. Love, Simon Anders, Vladislav Kim and Wolfgang Huber, 16 October, 2019
#https://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

#3.2.1. This section describes the script used to quantify transcripts from raw RNAseq-reads using salmon
#Version 16/06/2022
#Daniel Muliaditan

#Required input:
#-Raw RNAseq input reads
#-Gencode annotation file (for this script, gencode v38 was used

#Required tools: salmon, for install see: https://combine-lab.github.io/salmon/
#Prequisities:
#Install salmon: https://combine-lab.github.io/salmon/getting_started/ using conda
#Activate conda environment: conda activate salmon

#Output: 
#-Salmon index file
#-Salmon .quant files: these are the transcript quantification files used for downstream analysis

#Set working directories and variables
reference_dir=/mnt/raid5/cutnrun/reference/hg38
gencode_version=v38

#Salmon transcript quantification
#Index gencode transcripts to use for salmon
salmon index --gencode -t "$reference_dir"'/gencode.'"$gencode_version"'.transcripts.fa.gz' -i "$reference_dir"'/gencode.'"$gencode_version"'_salmon_1.5.2'
salmon_index="$reference_dir"'/gencode.'"$gencode_version"'_salmon_1.5.2'

for i in HN120Pri_a HN120Pri_b HN120Pri_c HN120Met_a HN120Met_b HN120Met_c HN120PCR_a HN120PCR_b HN120PCR_c HN137Pri_a HN137Met_a HN137Met_b HN137Met_c HN137PCR_a HN137PCR_b HN137PCR_c 
do

#Quantify with salmon
R1='/mnt/raid5/cutnrun/rnaseq/fastqdir/'"$i"'_1.fq.gz'
R2='/mnt/raid5/cutnrun/rnaseq/fastqdir/'"$i"'_2.fq.gz'
mkdir '/mnt/raid5/cutnrun/rnaseq/salmon_output/'"$i"
outdir=/'/mnt/raid5/cutnrun/rnaseq/salmon_output/'"$i"

salmon quant -i "$salmon_index" -l A -p 8 \
         --gcBias --seqBias --numGibbsSamples 20 -o "$outdir" \
         -1 "$R1" -2 "$R2"
         
done
