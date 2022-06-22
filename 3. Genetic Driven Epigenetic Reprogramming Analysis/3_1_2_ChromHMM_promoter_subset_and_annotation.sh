#!/usr/bin/env bash

#3.1.2. This section describes the script used to subset binned data to include only promoter bins and to annotate the bins with the nearest genes
#Version 22/06/2022
#Daniel Muliaditan

#Required input:
#-Output .bed files for each transition as created in section 3.1.1.
#-Reference .bed files with promoter regions. This file was attained using the promoters() function in R.
#-Reference GENCODE annotations with coordinates of the genes

#Required tools: bedtools

#Output: .bed files of each chromatin state transition, annotated with gene names and subset to promoter only regions.

for j in HN120Pri HN120Met HN120PCR HN137Pri HN137Met HN137PCR
do

echo $j

#Subset regions to include only promoter regions
bedtools intersect -a '/mnt/d/snCUT_RUN/results/ChromHMM/outputdir/5_states/'"$j"'_5_segments_binned_sorted.bed' \
-b /mnt/d/genome_references/hg38_reference/promoter_regions_hg38_sorted.bed -sorted \
-g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai > \
'/mnt/d/snCUT_RUN/results/ChromHMM/'"$j"'_binned_promoters.bed'

#Annotate with the closest genes in promoters
bedtools closest -a '/mnt/d/snCUT_RUN/results/ChromHMM/'"$j"'_binned_promoters.bed' \
-b /mnt/d/genome_references/hg38_reference/gencode.v32.basic.annotation.bed \
-g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai -d \
> '/mnt/d/snCUT_RUN/results/ChromHMM/'"$j"'_binned_promoters_nearest_gene.bed'

done
