#Bash code to run ChromHMM and subsetting after R processing

#Intersect ChromHMM regions
#Intersect E3-E4 regions
#Intersect regions that are in promoters
bedtools intersect -a /mnt/e/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E3_E4_regions.bed \
-b /mnt/e/genome_references/hg38_reference/promoter_regions_hg38_sorted.bed -sorted -g /mnt/e/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai >\
/mnt/e/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E3_E4_promoters.bed

#Closest genes in promoters
bedtools closest -a /mnt/e/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E3_E4_promoters.bed \
-b /mnt/e/genome_references/hg38_reference/gencode.v32.basic.annotation.bed -g /mnt/e/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai -d \
> /mnt/e/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E3_E4_regions_nearest_gene.bed

#Process to final bed
awk '{ a[$8]++ } END { for (b in a) { print b } }' /mnt/e/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E3_E4_regions_nearest_gene.bed > \
/mnt/e/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E3_E4_nearest_genes.bed
