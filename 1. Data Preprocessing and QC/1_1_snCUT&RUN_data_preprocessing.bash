#!/usr/bin/env bash

#Set variables and directories
longLine="--------------------"
SAMTOOLS='samtools'
GATK_DIR='/mnt/d/programmes/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar'
JAVA_DIR='java'
REFERENCE_DIR='/mnt/d/genome_references/hg38_reference'
BOWTIE="bowtie2"
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

for p in 120met_k27ac
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
  "$BOWTIE" -p "$THREADS" --end-to-end --very-sensitive --no-mixed --no-discordant \
  --met-file "$results_dir"'/processed_bams/sam/'"$i"'.txt' \
  -q --phred33 -I 10 -X 700 -x "$BOWTIE_INDEX"'/hg38' \
  -1 "$fastqdir"'/'"$i"'R1_001.fastq' \
  -2 "$fastqdir"'/'"$i"'R2_001.fastq' > "$results_dir"'/processed_bams/sam/'"$i"'.sam'

  msg="Convert sam to bam"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  "$SAMTOOLS" view -Sb "$results_dir"'/processed_bams/sam/'"$i"'.sam' -@ "$THREADS" > "$results_dir"'/processed_bams/bam/'"$i"'_unsort.bam'

  msg="Sort bam with "$SAMTOOLS""; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  "$SAMTOOLS" sort "$results_dir"'/processed_bams/bam/'"$i"'_unsort.bam' -@ "$THREADS" > "$results_dir"'/processed_bams/bam/'"$i"'_sorted.bam'

  msg="Mark and remove duplicates with Picard and sort with "$SAMTOOLS""; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  "$JAVA_DIR" -Xms8g -Xmx12g -jar "$GATK_DIR" MarkDuplicates \
  -I "$results_dir"'/processed_bams/bam/'"$i"'_sorted.bam' \
  -O "$results_dir"'/processed_bams/bam/'"$i"'_rmdup.bam' \
  -M "$results_dir"'/processed_bams/bam/'"$i"'_sorted_metrics.txt' \
  --ASSUME_SORTED true --REMOVE_DUPLICATES true
  "$SAMTOOLS" sort "$results_dir"'/processed_bams/bam/'"$i"'_rmdup.bam' -@ "$THREADS" > "$results_dir"'/final_bams/'"$i"'.bam'
  "$SAMTOOLS" index "$results_dir"'/final_bams/'"$i"'.bam' -@ "$THREADS"
  
  msg="Remove old and unsorted bam and sam, and preprocessing beds"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  rm "$results_dir"'/processed_bams/bam/'"$i"'_unsort.bam'
  rm "$results_dir"'/processed_bams/sam/'"$i"'.sam'
  rm "$results_dir"'/processed_bams/bam/'"$i"'_sorted.bam'
  rm "$results_dir"'/processed_bams/bam/'"$i"'_rmdup.bam'
 
  msg="Gzip fastq"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  gzip "$fastqdir"'/'"$i"'R1_001.fastq'
  gzip "$fastqdir"'/'"$i"'R2_001.fastq'

  done
  
  #Aggregate .bam filenames to antibody/cell line specific files
  msg="Aggregate individual single-cell .bam names to cell line-antibody .bam list"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  ls -d "$SEQUENCE_DIR"'/'"$k"'/results_dir/final_bams/'*'_.bam' > \
  "$AGGR_RESULTS_DIR"'/bamlists/'"$k"'_bamslist.txt'
  cat "$AGGR_RESULTS_DIR"'/bamlists/'"$k"'_bamslist.txt' >> \
  "$AGGR_RESULTS_DIR"'/bamlists/'"$p"'_bamslist.txt'
  
  msg="Finished"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
  

done

done

