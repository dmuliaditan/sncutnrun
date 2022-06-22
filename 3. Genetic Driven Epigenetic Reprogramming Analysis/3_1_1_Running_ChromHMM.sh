#!/usr/bin/env bash

#These series of code describes the analysis to correlate gene expression, gene copy number and gene chromatin state.
#The following workflow was used for this analysis:
#3.1. Genome-wide chromatin state annotation with ChromHMM
#3.2. Gene expression data analysis using RNAseq analysis
#3.3. Gene copy number analysis with CNVKit
#3.4. Static correlation between HN137Met gene CN and gene expression
#3.5. Static correlation between HN137Met chromatin state and gene expression
#3.6. Static correlation between HN137Met gene CN, chromatin state and gene expression
#3.7. Alluvial plot detailing the global changes of chromatin state occurring between HN137Pri and HN137Met
#3.8. Looking at enrichment of chromatin state comparing up- vs. downregulated genes in the HN137Pri > HN137Met transition
#3.9. Looking at enrichment of chromatin state in regions affected by CN change during the HN137Pri > HN137Met transition

#3.1.1. This section describes the script used to run ChromHMM to annotate genomic regions with chromatin state.
#Version 22/06/2022
#Daniel Muliaditan

#Required input:
#-Reference file with chromosome names and chromosomes sizes (in this analysis hg38 chromosome annotations were used)
#-Normalised .bam files of the samples and sample controls (in this case, rabbit IgG bulk CUT&RUN data was used)
#-Metadata cell mark table detailing the following columns: sample, histone mark, path to sample .bam, path to matched control .bam of that sample

#Install ChromHMM: http://compbio.mit.edu/ChromHMM/

#Set reference directory
REFERENCE_DIR='/mnt/d/genome_references/hg38_reference'
mkdir /mnt/d/snCUT_RUN/results/ChromHMM/control

#Binarize the normalised IgG control .bams and normalised sample .bams with BinarizeBam to create the input of LearnModel
java -Xms8g -Xmx12g -jar /mnt/d/programmes/ChromHMM/ChromHMM.jar BinarizeBam \
-paired -gzip \
-c /mnt/d/snCUT_RUN/data/final_bams/normalised_bams/control \ #Put the control .bams in this directory
-o /mnt/d/snCUT_RUN/results/ChromHMM/control \ #Output directory for the control .bams
"$REFERENCE_DIR"'/hg38.chrom.sizes.txt' \ #Reference chromosome size file
/mnt/d/snCUT_RUN/data/final_bams/normalised_bams \ #Input directory of the sample .bams
/mnt/d/snCUT_RUN/scripts/chromHMM_cell_mark_table.txt \ #Metadata detailing the input as described above
/mnt/d/snCUT_RUN/results/ChromHMM 

#Run ChromHMM LearnModel on the binarized .bams
#Rerun at different number of states and compare biological significance of the results
mkdir -p /mnt/d/snCUT_RUN/results/ChromHMM/outputdir/5_states
java -Xms8g -Xmx12g -jar /mnt/d/programmes/ChromHMM/ChromHMM.jar LearnModel \
-noautoopen -p 8 /mnt/d/snCUT_RUN/results/ChromHMM/inputdir \
/mnt/d/snCUT_RUN/results/ChromHMM/outputdir/5_states 5 hg38

for k in HN120Pri HN120Met HN120PCR HN137Pri HN137Met HN137PCR
do
echo "$k"
#Used bedtools makewindows to bin the ChromHMM output to 200bp
bedtools makewindows -b '/mnt/d/snCUT_RUN/results/ChromHMM/outputdir/5_states/'"$k"'_5_segments.bed' \
-i src -w 200 > '/mnt/d/snCUT_RUN/results/ChromHMM/outputdir/5_states/'"$k"'_5_segments_binned.bed'

#Remove random, unknown and chrM and sort and index to make final bam
cat '/mnt/d/snCUT_RUN/results/ChromHMM/outputdir/5_states/'"$k"'_5_segments_binned.bed' | grep -v random | grep -v Un | grep -v alt | grep -v chrM | grep -v HLA | grep -v EBV \
> '/mnt/d/snCUT_RUN/results/ChromHMM/outputdir/5_states/'"$k"'_5_segments_binned_presort.bed'
bedtools sort -i '/mnt/d/snCUT_RUN/results/ChromHMM/outputdir/5_states/'"$k"'_5_segments_binned_presort.bed' \
-faidx "$REFERENCE_DIR"'/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai' > \
'/mnt/d/snCUT_RUN/results/ChromHMM/outputdir/5_states/'"$k"'_5_segments_binned_sorted.bed'

done
