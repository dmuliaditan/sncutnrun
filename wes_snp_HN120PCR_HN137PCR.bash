#!/usr/bin/env bash

#GATK pipeline for preprocessing
#Convert fastq to unmapped bam

JAVA_DIR=/usr/lib/jvm/java-1.8.0/bin/java
GATK_DIR=/mnt/projects/muliaditand/muliaditand/programmes/gatk-4.2.0.0/gatk-package-4.2.0.0-local.jar
BWA=/mnt/projects/muliaditand/muliaditand/programmes/bwa/bwa
SAMTOOLS=/mnt/projects/muliaditand/muliaditand/programmes/samtools-1.11/samtools
FASTA=/mnt/projects/muliaditand/muliaditand/literature/gatk/Homo_sapiens_assembly38.fasta
DICT=/mnt/projects/muliaditand/muliaditand/literature/gatk/Homo_sapiens_assembly38.dict
REFERENCE_DIR=/mnt/projects/muliaditand/muliaditand/literature/gatk
FASTQ_DIR=/mnt/projects/muliaditand/muliaditand/dglab/HN_WES_data/fastqdir
RESULTS_DIR=/mnt/projects/muliaditand/muliaditand/dglab/snCUTRUN/wes/results_dir
TMP_DIR=/mnt/projects/muliaditand/muliaditand/dglab/HN_WES_data/tmp_dir

mkdir -p "$RESULTS_DIR"'/processed_bam'
mkdir -p "$RESULTS_DIR"'/final_bams'
mkdir -p "$RESULTS_DIR"'/SNVs'
mkdir -p "$RESULTS_DIR"'/SNVs/processing_snvs'
mkdir -p "$RESULTS_DIR"'/SNVs/final_SNVs'

for j in HN120PRI HN120PCR HN137PRI HN137PCR
do

#Convert FASTQ to SAM
"$JAVA_DIR" -Xms12g -jar "$GATK_DIR" FastqToSam \
F1="$FASTQ_DIR"'/'"$j"'_1.fastq.gz' \
F2="$FASTQ_DIR"'/'"$j"'_2.fastq.gz' \
O="$RESULTS_DIR"'/processed_bam/'"$j"'_unaligned.bam' \
SM="$j" \
LB="Illumina_DNA_Truseq_Nano" \
PL="ILLUMINA" \
R="$FASTA" \
RG="$j" \
TMP_DIR="$TMP_DIR"

#Mark Illumina Adapters
"$JAVA_DIR" -Xms16g -jar "$GATK_DIR" MarkIlluminaAdapters \
I="$RESULTS_DIR"'/processed_bam/'"$j"'_unaligned.bam' \
O="$RESULTS_DIR"'/processed_bam/'"$j"'_unaligned_adapter_marked.bam' \
M="$RESULTS_DIR"'/processed_bam/'"$j"'_unaligned_adapter_marked_metrics.txt' \
TMP_DIR="$TMP_DIR"

set -o pipefail

#Convert back to FASTQ, pipe into BWA, and pipe into MergeBamAlignment
"$JAVA_DIR" -Xms16g -jar "$GATK_DIR" SamToFastq \
I="$RESULTS_DIR"'/processed_bam/'"$j"'_unaligned_adapter_marked.bam' \
FASTQ="$RESULTS_DIR"'/processed_bam/'"$j"'.fastq' \
CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
TMP_DIR="$TMP_DIR"

"$BWA" mem -M -t 16 -p -o \
"$RESULTS_DIR"'/processed_bam/'"$j"'.sam' \
"$FASTA" "$RESULTS_DIR"'/processed_bam/'"$j"'.fastq'

"$SAMTOOLS" view -b -@ 13 -q 37 -T "$FASTA" \
"$RESULTS_DIR"'/processed_bam/'"$j"'.sam' \
-o "$RESULTS_DIR"'/processed_bam/'"$j"'_unsort.bam'

#Sort bam
"$SAMTOOLS" sort "$RESULTS_DIR"'/processed_bam/'"$j"'_unsort.bam' \
-o "$RESULTS_DIR"'/processed_bam/'"$j"'_sorted.bam' -@ 13

#Mark and remove duplicates
"$JAVA_DIR" -Xms16g -jar "$GATK_DIR" MarkDuplicates \
I="$RESULTS_DIR"'/processed_bam/'"$j"'_sorted.bam' \
O="$RESULTS_DIR"'/processed_bam/'"$j"'_rmdup.bam' \
M="$RESULTS_DIR"'/processed_bam/'"$j"'_metrics.txt' \
REMOVE_DUPLICATES=true R="$FASTA" \
TMP_DIR="$TMP_DIR"

"$SAMTOOLS" sort "$RESULTS_DIR"'/processed_bam/'"$j"'_rmdup.bam' \
-o "$RESULTS_DIR"'/final_bams/'"$j"'.bam' -@ 13

#Index bam
"$SAMTOOLS" index "$RESULTS_DIR"'/final_bams/'"$j"'.bam' -@ 13

#Remove sam and adapter marked BAM
rm "$RESULTS_DIR"'/processed_bam/'"$j"'_unaligned.bam'
rm "$RESULTS_DIR"'/processed_bam/'"$j"'_unaligned_adapter_marked.bam'
rm "$RESULTS_DIR"'/processed_bam/'"$j"'.fastq'
rm "$RESULTS_DIR"'/processed_bam/'"$j"'.sam'
rm "$RESULTS_DIR"'/processed_bam/'"$j"'_unsort.bam'
rm "$RESULTS_DIR"'/processed_bam/'"$j"'_sorted.bam'
rm "$RESULTS_DIR"'/processed_bam/'"$j"'_rmdup.bam'

done

msg="Preprocess WES intervals with PreprocessIntervals: $cell_line"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
"$JAVA_DIR" -Xms16g -jar "$GATK_DIR" PreprocessIntervals \
-O "$RESULTS_DIR"'/SNVs/WES.preprocessed.interval_list' \
-L /mnt/projects/muliaditand/muliaditand/literature/gatk/S07604514_Regions.bed \
-R "$FASTA" \
-sequence-dictionary "$DICT" \
--bin-length 0 --interval-merging-rule OVERLAPPING_ONLY \
--tmp-dir "$TMP_DIR"

for j in HN120 HN137
do

#Add read group
"$JAVA_DIR" -Xms16g -jar "$GATK_DIR" AddOrReplaceReadGroups \
I="$RESULTS_DIR"'/final_bams/'"$j"'PRI.bam' \
O="$RESULTS_DIR"'/final_bams/'"$j"'PRI_RG_added.bam' \
RGID=1 \
RGLB="$j"'PRI' \
RGPL=ILLUMINA \
RGPU=UNIT1 \
RGSM="$j"'PRI'

"$SAMTOOLS" index "$RESULTS_DIR"'/final_bams/'"$j"'PRI_RG_added.bam' -@ 13

#Call candidate variants with Mutect2
"$JAVA_DIR" -Xms16g -jar "$GATK_DIR" Mutect2 \
-I "$RESULTS_DIR"'/final_bams/'"$j"'PCR.bam' \
-I "$RESULTS_DIR"'/final_bams/'"$j"'PRI_RG_added.bam' \
-tumor "$j"'PCR' \
-normal "$j"'PRI' \
-O "$RESULTS_DIR"'/SNVs/processing_snvs/'"$j"'PCR_unfiltered.vcf.gz' \
-R "$FASTA" \
-L "$RESULTS_DIR"'/SNVs/WES.preprocessed.interval_list' \
--germline-resource "$REFERENCE_DIR"'/af-only-gnomad.hg38.vcf' \
-pon "$REFERENCE_DIR"'/somatic-hg38_1000g_pon.hg38.vcf' \
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter -default-af 0.0000025 \
-bamout "$RESULTS_DIR"'/SNVs/processing_snvs/'"$j"'PCR_m2.bam' \
--f1r2-tar-gz "$RESULTS_DIR"'/SNVs/processing_snvs/'"$j"'PCR_f1r2.tar.gz' \
--tmp-dir "$TMP_DIR"

#msg="LearnReadOrientationModel: $cell_line"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
"$JAVA_DIR" -Xms16g -jar "$GATK_DIR" LearnReadOrientationModel \
-I "$RESULTS_DIR"'/SNVs/processing_snvs/'"$j"'PCR_f1r2.tar.gz' \
-O "$RESULTS_DIR"'/SNVs/processing_snvs/'"$j"'PCR_read_orientation_model.tar.gz'

#msg="Run Getpileupsummaries: $cell_line"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
"$JAVA_DIR" -Xms16g -jar "$GATK_DIR" GetPileupSummaries \
-I "$RESULTS_DIR"'/final_bams/'"$j"'PCR.bam' \
-O "$RESULTS_DIR"'/SNVs/processing_snvs/'"$j"'PCR_getpileupsummaries.table' \
-L "$RESULTS_DIR"'/SNVs/WES.preprocessed.interval_list' \
-V "$REFERENCE_DIR"'/af-only-gnomad.hg38.vcf' \
-R "$FASTA" \
--tmp-dir "$TMP_DIR"

#msg="Estimate contamination with CalculateContamination: $cell_line"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
"$JAVA_DIR" -Xms16g -jar "$GATK_DIR" CalculateContamination \
-I "$RESULTS_DIR"'/SNVs/processing_snvs/'"$j"'PCR_getpileupsummaries.table' \
-O "$RESULTS_DIR"'/SNVs/processing_snvs/'"$j"'PCR_calculatecontamination.table' \
--tmp-dir "$TMP_DIR"

#msg="Filter Mutect2 calls: $cell_line"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
"$JAVA_DIR" -Xms16g -jar "$GATK_DIR" FilterMutectCalls \
-V "$RESULTS_DIR"'/SNVs/processing_snvs/'"$j"'PCR_unfiltered.vcf.gz' \
-O "$RESULTS_DIR"'/SNVs/processing_snvs/'"$j"'PCR_filtered.vcf.gz' \
-R "$FASTA" \
--contamination-table "$RESULTS_DIR"'/SNVs/processing_snvs/'"$j"'PCR_calculatecontamination.table'  \
--ob-priors "$RESULTS_DIR"'/SNVs/processing_snvs/'"$j"'PCR_read_orientation_model.tar.gz' \
-L "$RESULTS_DIR"'/SNVs/WES.preprocessed.interval_list' \
--tmp-dir "$TMP_DIR"

#msg="Functionally annotate variants with Funcotator: $cell_line"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
"$JAVA_DIR" -Xms16g -jar "$GATK_DIR" Funcotator \
--data-sources-path "$REFERENCE_DIR"'/funcotator_dataSources.v1.7.20200521s' \
-O "$RESULTS_DIR"'/SNVs/final_SNVs/'"$j"'PCR_SNV_annotated.vcf' \
--output-file-format VCF --ref-version hg38 --remove-filtered-variants true \
-R "$FASTA" \
-V "$RESULTS_DIR"'/SNVs/processing_snvs/'"$j"'PCR_filtered.vcf.gz' \
-L "$RESULTS_DIR"'/SNVs/WES.preprocessed.interval_list' \
--tmp-dir "$TMP_DIR"

#msg="Convert annotated vcf to table with VariantsToTable: $cell_line"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
"$JAVA_DIR" -Xms16g -jar "$GATK_DIR" VariantsToTable \
-V "$RESULTS_DIR"'/SNVs/final_SNVs/'"$j"'PCR_SNV_annotated.vcf' \
-O "$RESULTS_DIR"'/SNVs/final_SNVs/'"$j"'PCR_SNV.table' \
-F CHROM -F POS -F TYPE -F FUNCOTATION

msg="Extract missense, nonsense, nonstop and frameshift variants: $cell_line"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
grep -w 'NONSENSE' "$RESULTS_DIR"'/SNVs/final_SNVs/'"$j"'PCR_SNV.table' > "$RESULTS_DIR"'/SNVs/final_SNVs/'"$j"'PCR_nonsense.txt'
grep -w 'MISSENSE' "$RESULTS_DIR"'/SNVs/final_SNVs/'"$j"'PCR_SNV.table' > "$RESULTS_DIR"'/SNVs/final_SNVs/'"$j"'PCR_missense.txt'
grep -w 'NONSTOP' "$RESULTS_DIR"'/SNVs/final_SNVs/'"$j"'PCR_SNV.table' > "$RESULTS_DIR"'/SNVs/final_SNVs/'"$j"'PCR_nonstop.txt'
grep -w 'FRAME_SHIFT_INS'  "$RESULTS_DIR"'/SNVs/final_SNVs/'"$j"'PCR_SNV.table' > "$RESULTS_DIR"'/SNVs/final_SNVs/'"$j"'PCR_frameshift_ins.txt'
grep -w 'FRAME_SHIFT_DEL'  "$RESULTS_DIR"'/SNVs/final_SNVs/'"$j"'PCR_SNV.table' > "$RESULTS_DIR"'/SNVs/final_SNVs/'"$j"'PCR_frameshift_del.txt'

done
