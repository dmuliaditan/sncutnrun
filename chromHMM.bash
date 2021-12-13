#Bash code to run ChromHMM and subsetting after R processing

#Intersect ChromHMM regions


for j in E1 E2 E3 E4 E5
do

  for k in E1 E2 E3 E4 E5
  
  do
  
  echo "Intersect regions that are in promoters:" "$j" ">" "$k"
  bedtools intersect -a '/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_'"$j"'_'"$k"'_regions.bed' \
  -b /mnt/d/genome_references/hg38_reference/promoter_regions_hg38_sorted.bed -sorted -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai >\
  '/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_'"$j"'_'"$k"'_promoters.bed'
  
  #Closest genes in promoters
  bedtools closest -a '/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_'"$j"'_'"$k"'_promoters.bed' \
  -b /mnt/d/genome_references/hg38_reference/gencode.v32.basic.annotation.bed -g /mnt/d/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai -d \
  > '/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_'"$j"'_'"$k"'_regions_nearest_gene.bed'
  
  #Process to final bed
  awk '{ a[$8]++ } END { for (b in a) { print b } }' '/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_'"$j"'_'"$k"'_regions_nearest_gene.bed' > \
  '/mnt/d/snCUT_RUN/results/ChromHMM/HN137Pri_HN137Met_'"$j"'_'"$k"'_nearest_genes.bed'
  
  done
 
done





