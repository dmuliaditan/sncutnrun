#!/usr/bin/env bash

#Set variables and directories
longLine="--------------------"
SAMTOOLS='samtools'
GATK_DIR='/mnt/d/programmes/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar'
JAVA_DIR='java'
REFERENCE_DIR='/mnt/d/genome_references/hg38_reference'
BOWTIE="bowtie2"
BOWTIE_INDEX='/mnt/d/genome_references/bowtie2_indexes/hg38'
BEDTOOLS='bedtools'
MACS2='macs2'
SCRIPT_DIR='/mnt/d/snCUT_RUN/scripts'
SEQUENCE_DIR='/mnt/d/snCUT_RUN/sequence_runs'
AGGR_RESULTS_DIR='/mnt/d/snCUT_RUN/data'
THREADS=8

mkdir -p "$AGGR_RESULTS_DIR"'/tmp_dir'
TMP_DIR="$AGGR_RESULTS_DIR"'/tmp_dir'

for p in 120met_k27ac
do
  echo "$p"
  
  for o in `cat "$AGGR_RESULTS_DIR"'/bamlists/'"$p"'_bamslist.txt'`

  do

    echo "$o"

    z=${o#"$SEQUENCE_DIR"'/'}
    q=${z%_.bam}
    r=${q/"results_dir/final_bams/"}
    s=${r#*"/"}_${r%"/"*}
    t=${r%/*}
    u=$SEQUENCE_DIR'/'$t'/results_dir/final_bams/'
    v=${r#*"/"}
    
    #Convert single-cell bams to beds
    "$BEDTOOLS" bamtobed -i "$o" -bedpe \
    | awk '{if($2!="-1") print}' | sort -k1,1 -k2,2n | cut -f1,2,6,7 > "$TMP_DIR"'/'"$s"'_preproc.bed'
    "$BEDTOOLS" sort -i "$TMP_DIR"'/'"$s"'_preproc.bed' \
    -faidx "$REFERENCE_DIR"'/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai' > \
    "$TMP_DIR"'/'"$s"'_postsort.bed'
    nrow=$(wc -l < "$TMP_DIR"'/'"$s"'_postsort.bed')
    for ((c=1; c<=nrow; c++)) ; do echo "$s" ; done > "$TMP_DIR"'/'"$s"'.txt'
    awk 'FNR==NR{a[FNR]=$1;next};{$NF=a[FNR]};1' "$TMP_DIR"'/'"$s"'.txt' "$TMP_DIR"'/'"$s"'_postsort.bed' > "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'.bed'
    perl -p -i -e 's/ /\t/g' "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'.bed' 
    
    #1. Calculate Unique Mapped Reads (UMRs)
    # To calculate the number of Unique Mapped Reads for each paired ended and deduplicated single cell bam file, the following script was used (as per following reference).
    # http://qnot.org/2012/04/14/counting-the-number-of-reads-in-a-bam-file/
    # -c = count, -f 1 = only reads which are paired in sequencing, -F 12 means to include all reads where neither flag 0x0004 or 0x0008 is set, where 0x0004 is not unmapped reads and 0x0008 is where the mate is not unmapped (only include reads where it maps and its mate also map). 
    samtools view -c -f 1 -F 12 "$o" > "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'.txt'
    
    #2. Calculate Number of Reads in Peaks (FRiP)
    #Convert single-cell bams to beds
    #For each paired ended and deduplicated single cell bam file, a tagAlign bed file was created, intersected with reference macs2 narrowPeak file, then then the number of intersections were counted
    bedtools bamtobed -i "$o" | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' > "$u""$v"'.tagAlign'
    bedtools sort -i "$AGGR_RESULTS_DIR"'/macs2/'"$p"'_norm_sort_peaks.narrowPeak' | bedtools merge -i stdin | bedtools intersect -u -a "$u""$v"'.tagAlign' -b stdin | wc -l > "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'_overlapping_peaks.txt'

    #2. Calculate Number Reads in Blacklist (FRiB)
    #For each paired ended and deduplicated single cell bam file, a tagAlign bed file was created, intersected with reference Blacklist file, then then the number of intersections were counted
    bedtools sort -i "$REFERENCE_DIR"'/hg38.blacklist.bed' | bedtools merge -i stdin | bedtools intersect -u -a "$u""$v"'.tagAlign' -b stdin | wc -l > "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'_overlapping_blacklist.txt'
    
    #Compile barcode file per cell and append to master barcode file
    msg="Compile barcode file per cell and append to master barcode file"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
    echo "$s" > "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'_barcode.txt'
    echo "$t" > "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'_experiment.txt'

    paste "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'_barcode.txt' \
    "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'.txt' \
    "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'_experiment.txt' \
    "$AGGR_RESULTS_DIR"'/signac_fragments/'"$p"'_experiment.txt' \
    "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'_overlapping_peaks.txt' \
    "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'_overlapping_blacklist.txt' \
    > "$AGGR_RESULTS_DIR"'/signac_fragments/'"$s"'_barcodes_precount.txt'
    cat "$AGGR_RESULTS_DIR"'/signac_fragments/'"$s"'_barcodes_precount.txt' >> "$AGGR_RESULTS_DIR"'/signac_fragments/'"$p"'_barcodes_precount.txt'
    
    rm "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'_barcode.txt'
    rm "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'.txt'
    rm "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'_experiment.txt'
    rm "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'_overlapping_peaks.txt'
    rm "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'_overlapping_blacklist.txt'
    rm "$AGGR_RESULTS_DIR"'/signac_fragments/'"$s"'_barcodes_precount.txt'
   
   
  done
  
  #Combine beds into large fragment file
  msg="Combine beds into large bed file"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  cat "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'*'.bed' > "$AGGR_RESULTS_DIR"'/signac_fragments/'"$p"'_combined_preproc.bed'
    
  #Calculate % Reads in Peaks and % Reads in Blacklist
  msg="Calculate FRiP and % in blacklist"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  awk '{$7 = $5 / $2 * 100}1' "$AGGR_RESULTS_DIR"'/signac_fragments/'"$p"'_barcodes_precount.txt' | sort -n -k 1 | column -t > "$AGGR_RESULTS_DIR"'/signac_fragments/'"$p"'_barcodes_postcount.txt'
  awk '{$8 = $6 / $2 * 100}1' "$AGGR_RESULTS_DIR"'/signac_fragments/'"$p"'_barcodes_postcount.txt' | sort -n -k 1 | column -t > "$AGGR_RESULTS_DIR"'/signac_fragments/'"$p"'_barcodes.txt'

  #Remove processing FRiP files
  msg="Remove processing FRiP files"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  rm "$AGGR_RESULTS_DIR"'/signac_fragments/'"$p"'_barcodes_precount.txt'
  rm "$AGGR_RESULTS_DIR"'/signac_fragments/'"$p"'_barcodes_postcount.txt'
  rm "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'*'.bed'
 
  msg="Next AB/Cell line!"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  
done
  
#Done with the barcode file
msg="Done with the barcode file!"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"

#Finally, concatenate all the .bed files to make the Signac fragment
msg="Concatenate bedfiles and process to fragment file"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"

msg="H3K4me3"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
cat "$AGGR_RESULTS_DIR"'/signac_fragments/120pri_k4me3_combined_preproc.bed' \
"$AGGR_RESULTS_DIR"'/signac_fragments/120met_k4me3_combined_preproc.bed' \
"$AGGR_RESULTS_DIR"'/signac_fragments/120pcr_k4me3_combined_preproc.bed' \
"$AGGR_RESULTS_DIR"'/signac_fragments/137pcr_k4me3_combined_preproc.bed' \
"$AGGR_RESULTS_DIR"'/signac_fragments/137pri_k4me3_combined_preproc.bed' \
"$AGGR_RESULTS_DIR"'/signac_fragments/137met_k4me3_combined_preproc.bed' \
> "$AGGR_RESULTS_DIR"'/signac_fragments/H3K4me3_combined_preproc.bed'
cat "$AGGR_RESULTS_DIR"'/signac_fragments/H3K4me3_combined_preproc.bed' | awk '{print $0 "\t" "1"}' > "$AGGR_RESULTS_DIR"'/signac_fragments/H3K4me3_combined_presort.bed'
sort -k1,1 -k2,2n "$AGGR_RESULTS_DIR"'/signac_fragments/H3K4me3_combined_presort.bed' > "$AGGR_RESULTS_DIR"'/signac_fragments/H3K4me3_combined_sorted.bed'
bgzip -@ 8 "$AGGR_RESULTS_DIR"'/signac_fragments/H3K4me3_combined_sorted.bed'
tabix -p bed "$AGGR_RESULTS_DIR"'/signac_fragments/H3K4me3_combined_sorted.bed.gz'

#Remove fragment preprocessing files
msg="Remove fragment preprocessing files"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
rm "$AGGR_RESULTS_DIR"'/signac_fragments/H3K4me3_combined_preproc.bed'
rm "$AGGR_RESULTS_DIR"'/signac_fragments/H3K4me3_combined_presort.bed'

#Combine barcode file
msg="Combine barcode files"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
cat "$AGGR_RESULTS_DIR"'/signac_fragments/120pri_k4me3_barcodes.txt' \
"$AGGR_RESULTS_DIR"'/signac_fragments/120met_k4me3_barcodes.txt' \
"$AGGR_RESULTS_DIR"'/signac_fragments/120pcr_k4me3_barcodes.txt' \
"$AGGR_RESULTS_DIR"'/signac_fragments/137pcr_k4me3_barcodes.txt' \
"$AGGR_RESULTS_DIR"'/signac_fragments/137pri_k4me3_barcodes.txt' \
"$AGGR_RESULTS_DIR"'/signac_fragments/137met_k4me3_barcodes.txt' \
> "$AGGR_RESULTS_DIR"'/signac_fragments/H3K4me3_barcodes.txt'

msg="H3K27ac"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
cat "$AGGR_RESULTS_DIR"'/signac_fragments/120pri_k27ac_combined_preproc.bed' \
"$AGGR_RESULTS_DIR"'/signac_fragments/120met_k27ac_combined_preproc.bed' \
"$AGGR_RESULTS_DIR"'/signac_fragments/120pcr_k27ac_combined_preproc.bed' \
"$AGGR_RESULTS_DIR"'/signac_fragments/137pcr_k27ac_combined_preproc.bed' \
"$AGGR_RESULTS_DIR"'/signac_fragments/137pri_k27ac_combined_preproc.bed' \
"$AGGR_RESULTS_DIR"'/signac_fragments/137met_k27ac_combined_preproc.bed' \
> "$AGGR_RESULTS_DIR"'/signac_fragments/H3K27ac_combined_preproc.bed'
cat "$AGGR_RESULTS_DIR"'/signac_fragments/H3K27ac_combined_preproc.bed' | awk '{print $0 "\t" "1"}' > "$AGGR_RESULTS_DIR"'/signac_fragments/H3K27ac_combined_presort.bed'
sort -k1,1 -k2,2n "$AGGR_RESULTS_DIR"'/signac_fragments/H3K27ac_combined_presort.bed' > "$AGGR_RESULTS_DIR"'/signac_fragments/H3K27ac_combined_sorted.bed'
bgzip -@ 8 "$AGGR_RESULTS_DIR"'/signac_fragments/H3K27ac_combined_sorted.bed'
tabix -p bed "$AGGR_RESULTS_DIR"'/signac_fragments/H3K27ac_combined_sorted.bed.gz'

#Remove fragment preprocessing files
msg="Remove fragment preprocessing files"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
rm "$AGGR_RESULTS_DIR"'/signac_fragments/H3K27ac_combined_preproc.bed'
rm "$AGGR_RESULTS_DIR"'/signac_fragments/H3K27ac_combined_presort.bed'

#Combine barcode file
msg="Combine barcode files"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
cat "$AGGR_RESULTS_DIR"'/signac_fragments/120pri_k27ac_barcodes.txt' \
"$AGGR_RESULTS_DIR"'/signac_fragments/120met_k27ac_barcodes.txt' \
"$AGGR_RESULTS_DIR"'/signac_fragments/120pcr_k27ac_barcodes.txt' \
"$AGGR_RESULTS_DIR"'/signac_fragments/137pcr_k27ac_barcodes.txt' \
"$AGGR_RESULTS_DIR"'/signac_fragments/137pri_k27ac_barcodes.txt' \
"$AGGR_RESULTS_DIR"'/signac_fragments/137met_k27ac_barcodes.txt' \
> "$AGGR_RESULTS_DIR"'/signac_fragments/H3K27ac_barcodes.txt'

#At the end of the pipeline, the resulting outputs are:
#1. Fragment files to be used as Signac input file
#2. Barcode files to be used as metadata for Signac analysis
