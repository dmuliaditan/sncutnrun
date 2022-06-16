#!/usr/bin/env bash

REFERENCE_DIR='/mnt/d/genome_references/hg38_reference'
mkdir /mnt/d/snCUT_RUN/results/ChromHMM/control

#Binarize the normalised IgG control .bams and normalised sample .bams with BinarizeBam to create the input of LearnModel
java -Xms8g -Xmx12g -jar /mnt/d/programmes/ChromHMM/ChromHMM.jar BinarizeBam \
-paired -gzip \
-c /mnt/d/snCUT_RUN/data/final_bams/normalised_bams/control \
-o /mnt/d/snCUT_RUN/results/ChromHMM/control \
"$REFERENCE_DIR"'/hg38.chrom.sizes.txt' \
/mnt/d/snCUT_RUN/data/final_bams/normalised_bams \
/mnt/d/snCUT_RUN/scripts/chromHMM_cell_mark_table.txt \
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
