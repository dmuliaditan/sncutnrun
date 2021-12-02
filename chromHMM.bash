#Bash code to run ChromHMM and subsetting after R processing

#Intersect ChromHMM regions
#Intersect E3-E4 regions
#Intersect regions that are in promoters
bedtools intersect -a /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E3_E4_regions.bed \
-b /mnt/d/genome_references/hg38_reference/promoter_regions_hg38_sorted.bed -sorted -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai >\
/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E3_E4_promoters.bed

#Closest genes in promoters
bedtools closest -a /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E3_E4_promoters.bed \
-b /mnt/d/genome_references/hg38_reference/gencode.v32.basic.annotation.bed -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai -d \
> /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E3_E4_regions_nearest_gene.bed

#Process to final bed
awk '{ a[$8]++ } END { for (b in a) { print b } }' /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E3_E4_regions_nearest_gene.bed > \
/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E3_E4_nearest_genes.bed

#Intersect E1-E3 regions
#Intersect regions that are in promoters
bedtools intersect -a /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E1_E3_regions.bed \
-b /mnt/d/genome_references/hg38_reference/promoter_regions_hg38_sorted.bed -sorted -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai >\
/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E1_E3_promoters.bed

#Closest genes in promoters
bedtools closest -a /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E1_E3_promoters.bed \
-b /mnt/d/genome_references/hg38_reference/gencode.v32.basic.annotation.bed -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai -d \
> /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E1_E3_regions_nearest_gene.bed

#Process to final bed
awk '{ a[$8]++ } END { for (b in a) { print b } }' /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E1_E3_regions_nearest_gene.bed > \
/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E1_E3_nearest_genes.bed


#Intersect E1_E4 regions
#Intersect regions that are in promoters
bedtools intersect -a /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E1_E4_regions.bed \
-b /mnt/d/genome_references/hg38_reference/promoter_regions_hg38_sorted.bed -sorted -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai >\
/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E1_E4_promoters.bed

#Closest genes in promoters
bedtools closest -a /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E1_E4_promoters.bed \
-b /mnt/d/genome_references/hg38_reference/gencode.v32.basic.annotation.bed -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai -d \
> /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E1_E4_regions_nearest_gene.bed

#Process to final bed
awk '{ a[$8]++ } END { for (b in a) { print b } }' /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E1_E4_regions_nearest_gene.bed > \
/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E1_E4_nearest_genes.bed


#Intersect E1-E5 regions
#Intersect regions that are in promoters
bedtools intersect -a /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E1_E5_regions.bed \
-b /mnt/d/genome_references/hg38_reference/promoter_regions_hg38_sorted.bed -sorted -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai >\
/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E1_E5_promoters.bed

#Closest genes in promoters
bedtools closest -a /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E1_E5_promoters.bed \
-b /mnt/d/genome_references/hg38_reference/gencode.v32.basic.annotation.bed -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai -d \
> /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E1_E5_regions_nearest_gene.bed

#Process to final bed
awk '{ a[$8]++ } END { for (b in a) { print b } }' /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E1_E5_regions_nearest_gene.bed > \
/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E1_E5_nearest_genes.bed

#Intersect E2-E1 regions
#Intersect regions that are in promoters
bedtools intersect -a /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E2_E1_regions.bed \
-b /mnt/d/genome_references/hg38_reference/promoter_regions_hg38_sorted.bed -sorted -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai >\
/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E2_E1_promoters.bed

#Closest genes in promoters
bedtools closest -a /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E2_E1_promoters.bed \
-b /mnt/d/genome_references/hg38_reference/gencode.v32.basic.annotation.bed -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai -d \
> /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E2_E1_regions_nearest_gene.bed

#Process to final bed
awk '{ a[$8]++ } END { for (b in a) { print b } }' /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E2_E1_regions_nearest_gene.bed > \
/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E2_E1_nearest_genes.bed

#Intersect E2-E3 regions
#Intersect regions that are in promoters
bedtools intersect -a /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E2_E3_regions.bed \
-b /mnt/d/genome_references/hg38_reference/promoter_regions_hg38_sorted.bed -sorted -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai >\
/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E2_E3_promoters.bed

#Closest genes in promoters
bedtools closest -a /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E2_E3_promoters.bed \
-b /mnt/d/genome_references/hg38_reference/gencode.v32.basic.annotation.bed -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai -d \
> /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E2_E3_regions_nearest_gene.bed

#Process to final bed
awk '{ a[$8]++ } END { for (b in a) { print b } }' /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E2_E3_regions_nearest_gene.bed > \
/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E2_E3_nearest_genes.bed

#Intersect E2-E4 regions
#Intersect regions that are in promoters
bedtools intersect -a /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E2_E4_regions.bed \
-b /mnt/d/genome_references/hg38_reference/promoter_regions_hg38_sorted.bed -sorted -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai >\
/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E2_E4_promoters.bed

#Closest genes in promoters
bedtools closest -a /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E2_E4_promoters.bed \
-b /mnt/d/genome_references/hg38_reference/gencode.v32.basic.annotation.bed -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai -d \
> /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E2_E4_regions_nearest_gene.bed

#Process to final bed
awk '{ a[$8]++ } END { for (b in a) { print b } }' /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E2_E4_regions_nearest_gene.bed > \
/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E2_E4_nearest_genes.bed

#Intersect E4-E2 regions
#Intersect regions that are in promoters
bedtools intersect -a /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E4_E2_regions.bed \
-b /mnt/d/genome_references/hg38_reference/promoter_regions_hg38_sorted.bed -sorted -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai >\
/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E4_E2_promoters.bed

#Closest genes in promoters
bedtools closest -a /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E4_E2_promoters.bed \
-b /mnt/d/genome_references/hg38_reference/gencode.v32.basic.annotation.bed -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai -d \
> /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E4_E2_regions_nearest_gene.bed

#Process to final bed
awk '{ a[$8]++ } END { for (b in a) { print b } }' /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E4_E2_regions_nearest_gene.bed > \
/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E4_E2_nearest_genes.bed

#Intersect E4-E5 regions
#Intersect regions that are in promoters
bedtools intersect -a /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E4_E5_regions.bed \
-b /mnt/d/genome_references/hg38_reference/promoter_regions_hg38_sorted.bed -sorted -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai >\
/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E4_E5_promoters.bed

#Closest genes in promoters
bedtools closest -a /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E4_E5_promoters.bed \
-b /mnt/d/genome_references/hg38_reference/gencode.v32.basic.annotation.bed -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai -d \
> /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E4_E5_regions_nearest_gene.bed

#Process to final bed
awk '{ a[$8]++ } END { for (b in a) { print b } }' /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E4_E5_regions_nearest_gene.bed > \
/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E4_E5_nearest_genes.bed


#Intersect E5-E4 regions
#Intersect regions that are in promoters
bedtools intersect -a /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E5_E4_regions.bed \
-b /mnt/d/genome_references/hg38_reference/promoter_regions_hg38_sorted.bed -sorted -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai >\
/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E5_E4_promoters.bed

#Closest genes in promoters
bedtools closest -a /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E5_E4_promoters.bed \
-b /mnt/d/genome_references/hg38_reference/gencode.v32.basic.annotation.bed -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai -d \
> /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E5_E4_regions_nearest_gene.bed

#Process to final bed
awk '{ a[$8]++ } END { for (b in a) { print b } }' /mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E5_E4_regions_nearest_gene.bed > \
/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_E5_E4_nearest_genes.bed




