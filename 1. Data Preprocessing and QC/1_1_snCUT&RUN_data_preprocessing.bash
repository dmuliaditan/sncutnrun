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

mkdir -p "$AGGR_RESULTS_DIR"'/bamlists'
mkdir -p "$AGGR_RESULTS_DIR"'/aggregate_bams'
mkdir -p "$AGGR_RESULTS_DIR"'/macs2'
mkdir -p "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs'
mkdir -p "$AGGR_RESULTS_DIR"'/signac_fragments'

for p in 120pri_K4me3 120pri_k27ac etc
do

#Set sample metadata
#Each *_experiment.txt file has two columns with the following outline
#Experiment_date  Experiment_name e.g.
#DM_20200909_1  HN120Met_K4me3
#DM_20201007  HN120Met_K4me3

table="$SCRIPT_DIR"'/'"$p"'_experiment.txt'
dos2unix "$table"

mkdir -p "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp'

#Create list of bam files
touch "$AGGR_RESULTS_DIR"'/bamlists/'"$p"'_bamslist.txt'
touch "$AGGR_RESULTS_DIR"'/signac_fragments/'"$p"'_barcodes_precount.txt'
echo "$p" > "$AGGR_RESULTS_DIR"'/signac_fragments/'"$p"'_experiment.txt'

#For each cell line - antibody combination, loop through each experiment/batch
cut -f1 "$table" |
while read k; do


  exp=`awk -v j=$k '{
    if($1 ==j)
    {
      print $2
    }
  }' "$table"`


  longLine="--------------------"

  msg="Experiment: $exp"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  msg="Date: $k"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"

  msg="Set up directories and filenames"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  fastqdir="$SEQUENCE_DIR"'/'"$k"'/fastqdir'
  results_dir="$SEQUENCE_DIR"'/'"$k"'/results_dir'
  mkdir -p "$results_dir"'/processed_bams'
  mkdir -p "$results_dir"'/processed_bams/sam'
  mkdir -p "$results_dir"'/processed_bams/bam'
  mkdir -p "$results_dir"'/final_bams'
  mkdir -p "$results_dir"'/beds_and_bedgraphs'
  mkdir -p "$results_dir"'/tmp_dir'

  cd "$fastqdir"
  msg="Unzip fastq"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  gzip -d *.fastq.gz

  ls *R1_001.fastq | sed -e 's/R1_001.fastq//' > "$fastqdir"'/filenames.txt'

  for i in `cat "$fastqdir"'/filenames.txt'`

  do

  msg="Cell: $i"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  msg="Map with bowtie2"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  "$BOWTIE" -p 8 --end-to-end --very-sensitive --no-mixed --no-discordant \
  --met-file "$results_dir"'/processed_bams/sam/'"$i"'.txt' \
  -q --phred33 -I 10 -X 700 -x "$BOWTIE_INDEX"'/hg38' \
  -1 "$fastqdir"'/'"$i"'R1_001.fastq' \
  -2 "$fastqdir"'/'"$i"'R2_001.fastq' > "$results_dir"'/processed_bams/sam/'"$i"'.sam'

  msg="Convert sam to bam"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  "$SAMTOOLS" view -Sb "$results_dir"'/processed_bams/sam/'"$i"'.sam' -@ 8 > "$results_dir"'/processed_bams/bam/'"$i"'_unsort.bam'

  msg="Sort bam with "$SAMTOOLS""; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  "$SAMTOOLS" sort "$results_dir"'/processed_bams/bam/'"$i"'_unsort.bam' -@ 8 > "$results_dir"'/processed_bams/bam/'"$i"'_sorted.bam'

  msg="Mark and remove duplicates with Picard and sort with "$SAMTOOLS""; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  "$JAVA_DIR" -Xms8g -Xmx12g -jar "$GATK_DIR" MarkDuplicates \
  -I "$results_dir"'/processed_bams/bam/'"$i"'_sorted.bam' \
  -O "$results_dir"'/processed_bams/bam/'"$i"'_rmdup.bam' \
  -M "$results_dir"'/processed_bams/bam/'"$i"'_sorted_metrics.txt' \
  --ASSUME_SORTED true --REMOVE_DUPLICATES true
  "$SAMTOOLS" sort "$results_dir"'/processed_bams/bam/'"$i"'_rmdup.bam' -@ 8 > "$results_dir"'/final_bams/'"$i"'.bam'
  "$SAMTOOLS" index "$results_dir"'/final_bams/'"$i"'.bam' -@ 8
  
  # To calculate the number of Unique Mapped Reads for each paired ended and deduplicated single cell bam file, the following script was used (as per following reference).
  # http://qnot.org/2012/04/14/counting-the-number-of-reads-in-a-bam-file/
  # -c = count, -f 1 = only reads which are paired in sequencing, -F 12 means to include all reads where neither flag 0x0004 or 0x0008 is set, where 0x0004 is not unmapped reads and 0x0008 is where the mate is not unmapped (only include reads where it maps and its mate also map). 
  samtools view -c -f 1 -F 12 "$results_dir"'/final_bams/'"$i"'.bam' >
  

  msg="Remove old and unsorted bam and sam, and preprocessing beds"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  rm "$results_dir"'/processed_bams/bam/'"$i"'_unsort.bam'
  rm "$results_dir"'/processed_bams/sam/'"$i"'.sam'
  rm "$results_dir"'/processed_bams/bam/'"$i"'_sorted.bam'
  rm "$results_dir"'/processed_bams/bam/'"$i"'_rmdup.bam'
 
  msg="Gzip fastq"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  gzip "$fastqdir"'/'"$i"'R1_001.fastq'
  gzip "$fastqdir"'/'"$i"'R2_001.fastq'

  done

  
  
  msg="Aggregate individual single-cell .bam names to cell line-antibody .bam list"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  ls -d "$SEQUENCE_DIR"'/'"$k"'/results_dir/final_bams/'*'_.bam' > \
  "$AGGR_RESULTS_DIR"'/bamlists/'"$k"'_bamslist.txt'
  cat "$AGGR_RESULTS_DIR"'/bamlists/'"$k"'_bamslist.txt' >> \
  "$AGGR_RESULTS_DIR"'/bamlists/'"$p"'_bamslist.txt'
  
  msg="Finished"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"

done

  msg="Merge all .bams to form pseudobulk, sort and index"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  "$SAMTOOLS" merge -f -b "$AGGR_RESULTS_DIR"'/bamlists/'"$p"'_bamslist.txt' "$AGGR_RESULTS_DIR"'/aggregate_bams'"$p"'_unsort.bam' --threads "$THREADS"
  "$SAMTOOLS" sort "$AGGR_RESULTS_DIR"'/aggregate_bams'"$p"'_unsort.bam' > "$AGGR_RESULTS_DIR"'/aggregate_bams'"$p"'_sorted.bam'
  "$SAMTOOLS" index "$AGGR_RESULTS_DIR"'/aggregate_bams'"$p"'_sorted.bam'

  msg="Convert aggregated .bam to .bedgraph and .bed"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  "$SAMTOOLS" sort -n -o "$AGGR_RESULTS_DIR"'/aggregate_bams'"$p"'_namesorted.bam' -@ "$THREADS" "$AGGR_RESULTS_DIR"'/aggregate_bams'"$p"'_sorted.bam'
  "$BEDTOOLS" bamtobed -i "$AGGR_RESULTS_DIR"'/aggregate_bams'"$p"'_namesorted.bam' -bedpe \
  | awk '{if($2!="-1") print}' | sort -k1,1 -k2,2n | cut -f1,2,6,7 > "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/'"$p"'_preproc.bed'
  cat "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/'"$p"'_preproc.bed' | grep -v random | grep -v Un | grep -v alt | grep -v chrM | grep -v HLA | grep -v EBV \
  > "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/'"$p"'_presort.bed'
  "$BEDTOOLS" sort -i "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/'"$p"'_presort.bed' \
  -faidx "$REFERENCE_DIR"'/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai' > \
  "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/'"$p"'_pseudobulk.bed'
  perl -p -i -e 's/ /\t/g' "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/'"$p"'_pseudobulk.bed' #Remove empty lines etc, this is the final bed

  rm "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/'"$p"'_preproc.bed'
  rm "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/'"$p"'_presort.bed'

  "$BEDTOOLS" genomecov -i "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/'"$p"'_pseudobulk.bed' \
  -g "$REFERENCE_DIR"'/hg38.chrom.sizes.txt' -bg > "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/'"$p"'_pseudobulk.bedGraph'

  #Call peaks on each merged bam using MACS2
  msg="Call peaks on each merged bam using "$MACS2""; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  "$MACS2" callpeak -t "$AGGR_RESULTS_DIR"'/aggregate_bams'"$p"'_sorted.bam' --outdir "$AGGR_RESULTS_DIR"'/macs2' -n "$p"'_q0_1' \
  -f AUTO -g hs --keep-dup auto -B --verbose 2 --nomodel --shift 0 --ext 200 --d-min 20 --qval 0.1 --SPMR --broad --broad-cutoff 0.1
  
  bedtools bamtobed -i /new/howard_new/sequence_runs/pathToBamFiles/S${i}_PE_rmdup1.bam | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' > /new/howard_new/sequence_runs/pathToBamFiles/S${i}_PE_rmdup1.tagAlign
bedtools sort -i /mnt/raid5/cutnrun/userFiles/howard/macs2normalised/120pri_k27ac_norm_sort_peaks.narrowPeak | bedtools merge -i stdin | bedtools intersect -u -a /new/howard_new/sequence_runs/pathToBamFiles/S${i}_PE_rmdup1.tagAlign -b stdin | wc -l
done

  #Convert peaks to beds
  msg="Convert peaks to beds"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  cat "$AGGR_RESULTS_DIR"'/macs2/'"$p"'_q0_1_peaks.broadPeak' | awk '{ print $1, $2, $3 }' > \
  "$AGGR_RESULTS_DIR"'/macs2/'"$p"'_q0_1_peaks.bed'
  sort -k1V -k2,2n "$AGGR_RESULTS_DIR"'/macs2/'"$p"'_q0_1_peaks.bed' > "$AGGR_RESULTS_DIR"'/macs2/'"$p"'_q0_1_sorted.bed'
  cat "$AGGR_RESULTS_DIR"'/macs2/'"$p"'_q0_1_sorted.bed' | grep -v random | grep -v Un | grep -v alt | grep -v chrM | grep -v HLA | grep -v EBV > "$AGGR_RESULTS_DIR"'/macs2/'"$p"'_peaks.bed'
  perl -p -i -e 's/ /\t/g' "$AGGR_RESULTS_DIR"'/macs2/'"$p"'_peaks.bed' #Final peak file

  #Remove unsorted beds and beds with unknown contigs
  msg="Remove unsorted beds and beds with unknown contigs"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  rm "$AGGR_RESULTS_DIR"'/macs2/'"$p"'_q0_1_peaks.bed'
  rm "$AGGR_RESULTS_DIR"'/macs2/'"$p"'_q0_1_sorted.bed'

  msg="Done with aggregate peak calling!"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"

  for o in `cat "$AGGR_RESULTS_DIR"'/bamlists/'"$p"'_bamslist_namesorted.txt'`

  do

    echo "$o"

    p=${o#"$SEQUENCE_DIR"}
    q=${p%_namesorted.bam}
    r=${q/"results_dir/final_bams/"}
    s=${r#*"/"}_${r%"/"*}
    t=${r%/*}

    #Convert single-cell bams to beds
    "$BEDTOOLS" bamtobed -i "$o" -bedpe \
    | awk '{if($2!="-1") print}' | sort -k1,1 -k2,2n | cut -f1,2,6,7 > "$TMP_DIR"'/'"$s"'_preproc.bed'
    cat "$TMP_DIR"'/'"$s"'_preproc.bed' | grep -v random | grep -v Un | grep -v alt | grep -v chrM | grep -v HLA | grep -v EBV \
    > "$TMP_DIR"'/'"$s"'_presort.bed'
    "$BEDTOOLS" sort -i "$TMP_DIR"'/'"$s"'_presort.bed' \
    -faidx "$REFERENCE_DIR"'/resources_broad_hg38_v0_Homo_sapiens_assembly38.fasta.fai' > \
    "$TMP_DIR"'/'"$s"'_postsort.bed'
    nrow=$(wc -l < "$TMP_DIR"'/'"$s"'_postsort.bed')
    for ((c=1; c<=nrow; c++)) ; do echo "$s" ; done > "$TMP_DIR"'/'"$s"'.txt'
    awk 'FNR==NR{a[FNR]=$1;next};{$NF=a[FNR]};1' "$TMP_DIR"'/'"$s"'.txt' "$TMP_DIR"'/'"$s"'_postsort.bed' > "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/'"$s"'.bed'
    perl -p -i -e 's/ /\t/g' "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'.bed'  #Remove empty lines etc, this is the final bed

    #Remove unsorted and beds with unknown contigs
    msg="Remove unsorted and beds with unknown contigs"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"

    rm "$TMP_DIR"'/'"$s"'_preproc.bed'
    rm "$TMP_DIR"'/'"$s"'_presort.bed'
    rm "$TMP_DIR"'/'"$s"'_postsort.bed'
    rm "$TMP_DIR"'/'"$s"'.txt'

    #Compile barcode file per cell and append to master barcode file
    msg="Compile barcode file per cell and append to master barcode file"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
    echo "$s" > "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'_barcode.txt'
    echo "$t" > "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'_experiment.txt'

    paste "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'_barcode.txt' "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'.txt' \
    "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'_experiment.txt' \
    "$AGGR_RESULTS_DIR"'/signac_fragments/'"$p"'_experiment.txt' \
    "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'overlapping_peaks.txt'\
    "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'overlapping_blacklist.txt' \
    > "$AGGR_RESULTS_DIR"'/signac_fragments/'"$s"'_barcodes_precount.txt'
    cat "$AGGR_RESULTS_DIR"'/signac_fragments/'"$s"'_barcodes_precount.txt' >> "$AGGR_RESULTS_DIR"'/signac_fragments/'"$p"'_barcodes_precount.txt'

    #Remove pre-processing txt files
    msg="Remove pre-processing txt files"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
    rm "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'overlapping_peaks.txt'
    rm "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'overlapping_blacklist.txt'
    rm "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'.txt'
    rm "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'_barcode.txt'
    rm "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'"$s"'_experiment.txt'
    rm "$AGGR_RESULTS_DIR"'/signac_fragments/'"$s"'_barcodes_precount.txt'
    rm "$AGGR_RESULTS_DIR"'/signac_fragments/'"$p"'_experiment.txt'

  done

  #Combine beds into large fragment file
  msg="Combine beds into large bed file"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  cat "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp/'*'.bed' > "$AGGR_RESULTS_DIR"'/signac_fragments/'"$p"'_combined_preproc.bed'

  #Calculate FRiP and % in blacklist
  msg="Calculate FRiP and % in blacklist"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  awk '{$7 = $5 / $2 * 100}1' "$AGGR_RESULTS_DIR"'/signac_fragments/'"$p"'_barcodes_precount.txt' | sort -n -k 1 | column -t > "$AGGR_RESULTS_DIR"'/signac_fragments/'"$p"'_barcodes_postcount.txt'
  awk '{$8 = $6 / $2 * 100}1' "$AGGR_RESULTS_DIR"'/signac_fragments/'"$p"'_barcodes_postcount.txt' | sort -n -k 1 | column -t > "$AGGR_RESULTS_DIR"'/signac_fragments/'"$p"'_barcodes.txt'

  #Remove processing FRiP files
  msg="Remove processing FRiP files"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  rm "$AGGR_RESULTS_DIR"'/signac_fragments/'"$p"'_barcodes_precount.txt'
  rm "$AGGR_RESULTS_DIR"'/signac_fragments/'"$p"'_barcodes_postcount.txt'
  rm -r "$AGGR_RESULTS_DIR"'/aggregate_beds_and_bedgraphs/tmp'

  #Done with the barcode file
  msg="Done with the barcode file!"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"

  msg="Next AB/Cell line!"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"

done

msg="Concatenate bedfiles and process to fragment file"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"

msg="K4me3"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
cat "$results_dir"'/final_fragments/120pri_k4me3/120pri_k4me3_combined_preproc.bed' "$results_dir"'/final_fragments/120met_k4me3/120met_k4me3_combined_preproc.bed' \
"$results_dir"'/final_fragments/120pcr_k4me3/120pcr_k4me3_combined_preproc.bed' "$results_dir"'/final_fragments/137pcr_k4me3/137pcr_k4me3_combined_preproc.bed' \
"$results_dir"'/final_fragments/137pri_k4me3/137pri_k4me3_combined_preproc.bed' "$results_dir"'/final_fragments/137met_k4me3/137met_k4me3_combined_preproc.bed' \
> "$results_dir"'/final_fragments/K4me3_combined_preproc.bed'
cat "$results_dir"'/final_fragments/K4me3_combined_preproc.bed' | awk '{print $0 "\t" "1"}' > "$results_dir"'/final_fragments/K4me3_combined_presort.bed'
sort -k1,1 -k2,2n "$results_dir"'/final_fragments/K4me3_combined_presort.bed' > "$results_dir"'/final_fragments/K4me3_combined_sorted.bed'
bgzip -@ 8 "$results_dir"'/final_fragments/K4me3_combined_sorted.bed'
tabix -p bed "$results_dir"'/final_fragments/K4me3_combined_sorted.bed.gz'

#Remove fragment preprocessing files
msg="Remove fragment preprocessing files"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
rm "$results_dir"'/final_fragments/K4me3_combined_preproc.bed'
rm "$results_dir"'/final_fragments/K4me3_combined_presort.bed'

#Combine barcode file
msg="Combine barcode files"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
cat "$results_dir"'/final_fragments/120pri_k4me3/120pri_k4me3_barcodes.txt' "$results_dir"'/final_fragments/120met_k4me3/120met_k4me3_barcodes.txt' \
"$results_dir"'/final_fragments/120pcr_k4me3/120pcr_k4me3_barcodes.txt' "$results_dir"'/final_fragments/137pcr_k4me3/137pcr_k4me3_barcodes.txt' \
"$results_dir"'/final_fragments/137pri_k4me3/137pri_k4me3_barcodes.txt' "$results_dir"'/final_fragments/137met_k4me3/137met_k4me3_barcodes.txt' \
> "$results_dir"'/final_fragments/K4me3_barcodes.txt'

msg="K27ac"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
cat "$results_dir"'/final_fragments/120pri_k27ac/120pri_k27ac_combined_preproc.bed' "$results_dir"'/final_fragments/120met_k27ac/120met_k27ac_combined_preproc.bed' \
"$results_dir"'/final_fragments/120pcr_k27ac/120pcr_k27ac_combined_preproc.bed' "$results_dir"'/final_fragments/137pcr_k27ac/137pcr_k27ac_combined_preproc.bed' \
"$results_dir"'/final_fragments/137pri_k27ac/137pri_k27ac_combined_preproc.bed' "$results_dir"'/final_fragments/137met_k27ac/137met_k27ac_combined_preproc.bed' \
> "$results_dir"'/final_fragments/K27ac_combined_preproc.bed'
cat "$results_dir"'/final_fragments/K27ac_combined_preproc.bed' | awk '{print $0 "\t" "1"}' > "$results_dir"'/final_fragments/K27ac_combined_presort.bed'
sort -k1,1 -k2,2n "$results_dir"'/final_fragments/K27ac_combined_presort.bed' > "$results_dir"'/final_fragments/K27ac_combined_sorted.bed'
bgzip -@ 8 "$results_dir"'/final_fragments/K27ac_combined_sorted.bed'
tabix -p bed "$results_dir"'/final_fragments/K27ac_combined_sorted.bed.gz'

#Remove fragment preprocessing files
msg="Remove fragment preprocessing files"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
rm "$results_dir"'/final_fragments/K27ac_combined_preproc.bed'
rm "$results_dir"'/final_fragments/K27ac_combined_presort.bed'

#Combine barcode file
msg="Combine barcode files"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
cat "$results_dir"'/final_fragments/120pri_k27ac/120pri_k27ac_barcodes.txt' "$results_dir"'/final_fragments/120met_k27ac/120met_k27ac_barcodes.txt' \
"$results_dir"'/final_fragments/120pcr_k27ac/120pcr_k27ac_barcodes.txt' "$results_dir"'/final_fragments/137pcr_k27ac/137pcr_k27ac_barcodes.txt' \
"$results_dir"'/final_fragments/137pri_k27ac/137pri_k27ac_barcodes.txt' "$results_dir"'/final_fragments/137met_k27ac/137met_k27ac_barcodes.txt' \
> "$results_dir"'/final_fragments/K27ac_barcodes.txt'

#At the end of the pipeline, the resulting outputs are:
#1. Fragment files to be used as Signac input file
#2. Barcode files to be used as metadata for Signac analysis
