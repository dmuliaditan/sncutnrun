#CNVkit v0.9.8

cnvkit.py access -o /home/daniel/genome_references/hg38_reference/hg38_accessible_regions_5kb.bed /home/daniel/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta

cnvkit.py batch -m hybrid --drop-low-coverage \
-p 8 -n -f /home/daniel/genome_references/hg38_reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta \
-t /home/daniel/genome_references/hg38_reference/S07604514_Padded_Agilent.bed --annotate /home/daniel/genome_references/hg38_reference/refFlat.txt \
-g /home/daniel/genome_references/hg38_reference/hg38_accessible_regions_5kb.bed \
--output-reference /home/daniel/snCUT_RUN/wes/CNV/cnvkit/hg38_WES_reference.cnn -d /home/daniel/snCUT_RUN/wes/CNV/cnvkit

for j in HN120PRI_subsampled_sorted HN120MET_subsampled_sorted HN120PCR_subsampled_sorted HN137PRI_subsampled_sorted HN137MET_subsampled_sorted HN137PCR_subsampled_sorted
do

#CNVkit run on Docker
#Create reference

cnvkit.py batch --drop-low-coverage -p 8 \
-r /home/daniel/snCUT_RUN/wes/CNV/cnvkit/hg38_WES_reference.cnn \
-d /home/daniel/snCUT_RUN/wes/CNV/cnvkit '/home/daniel/snCUT_RUN/wes/results_dir/final_bams/'"$j"'.bam' --scatter --diagram

done 
