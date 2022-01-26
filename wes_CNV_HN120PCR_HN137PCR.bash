#!/usr/bin/env bash


longLine="--------------------"
cd /mnt/projects/muliaditand/muliaditand/dglab/jay


msg="Set up directories and filenames"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
JAVA_DIR=/usr/lib/jvm/java-1.8.0/bin/java
results_dir=/mnt/projects/muliaditand/muliaditand/dglab/snCUTRUN/wes/results_dir
pon_dir=/mnt/projects/muliaditand/muliaditand/literature/gatk/pon_1000g_exome_hg38
gatk=/mnt/projects/muliaditand/muliaditand/programmes/gatk-4.1.4.1/gatk-package-4.1.4.1-local.jar
fasta=/mnt/projects/muliaditand/muliaditand/literature/gatk/Homo_sapiens_assembly38.fasta
dict=/mnt/projects/muliaditand/muliaditand/literature/gatk/Homo_sapiens_assembly38.dict
exac=/mnt/projects/muliaditand/muliaditand/literature/gatk/small_exac_common_3.hg38.vcf
tmp_dir=/mnt/projects/muliaditand/muliaditand/dglab/HN_WES_data/tmp_dir

msg="Set variables"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
mkdir -p "$results_dir"'/CNV'
mkdir -p "$results_dir"'/CNV/pre_processing'

#msg="Preprocess WES intervals with PreprocessIntervals: $cell_line"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
#"$JAVA_DIR" -Xms16g -jar "$gatk" PreprocessIntervals \
#-O "$results_dir"'/CNV/pre_processing/WES.preprocessed.interval_list' \
#-L /mnt/projects/muliaditand/muliaditand/literature/gatk/S07604514_Regions.bed \
#-R "$fasta" \
#-sequence-dictionary "$dict" \
#--bin-length 0 --interval-merging-rule OVERLAPPING_ONLY \
#--tmp-dir "$tmp_dir"

for i in HN120PRI HN120PCR HN137PRI HN137PCR

do

msg="Cell_line: $i"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"

#msg="Collect raw counts data with CollectReadCounts: $cell_line"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
"$JAVA_DIR" -Xms16g -jar "$gatk" CollectReadCounts \
-I "$results_dir"'/final_bams/'"$i"'_RG_added.bam' \
-L "$results_dir"'/CNV/pre_processing/WES.preprocessed.interval_list' \
-O "$results_dir"'/CNV/pre_processing/'"$i"'_readcount.tsv' \
--format TSV --interval-merging-rule OVERLAPPING_ONLY \
-sequence-dictionary "$dict" \
--tmp-dir "$tmp_dir"

#msg="Annotate Exome Target intervals"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
#"$JAVA_DIR" -Xms16g -jar "$gatk" AnnotateIntervals \
#-R "$fasta" \
#-L "$results_dir"'/CNV/pre_processing/WES.preprocessed.interval_list' \
#--interval-merging-rule OVERLAPPING_ONLY \
#-O "$pon_dir"'/cnv_pon/annotated_intervals.tsv'

msg="Standardize and denoise case read counts against the PoN with DenoiseReadCounts: $cell_line"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
"$JAVA_DIR" -Xms16g -jar "$gatk" DenoiseReadCounts \
--denoised-copy-ratios "$results_dir"'/CNV/pre_processing/'"$i"'_clean.denoisedCR.tsv' \
--standardized-copy-ratios "$results_dir"'/CNV/pre_processing/'"$i"'_clean.standardizedCR.tsv' \
-I "$results_dir"'/CNV/pre_processing/'"$i"'_readcount.tsv' \
--annotated-intervals "$pon_dir"'/cnv_pon/annotated_intervals.tsv' \
--count-panel-of-normals "$pon_dir"'/cnv_pon/DM_exome.pon.hdf5' \
--tmp-dir "$tmp_dir"

msg="Plot standardized and denoised copy ratios with PlotDenoisedCopyRatios: $cell_line"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
"$JAVA_DIR" -Xms16g -jar "$gatk" PlotDenoisedCopyRatios \
--denoised-copy-ratios "$results_dir"'/CNV/pre_processing/'"$i"'_clean.denoisedCR.tsv' \
-O "$results_dir"'/CNV/pre_processing/' \
--output-prefix "$i" \
--sequence-dictionary "$dict" \
--standardized-copy-ratios "$results_dir"'/CNV/pre_processing/'"$i"'_clean.standardizedCR.tsv' \
--tmp-dir "$tmp_dir"

msg="Count ref and alt alleles at common germline variant sites using CollectAllelicCounts: $cell_line"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
"$JAVA_DIR" -Xms16g -jar "$gatk" CollectAllelicCounts \
-I "$results_dir"'/final_bams/'"$i"'_RD_added.bam' \
-L "$exac" \
-O "$results_dir"'/CNV/pre_processing/'"$i"'_clean.allelicCounts.tsv' \
-R "$fasta" \
--sequence-dictionary "$dict" \
--tmp-dir "$tmp_dir"

msg="Group contiguous copy ratios into segments with ModelSegments: $cell_line"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
"$JAVA_DIR" -Xms16g -jar "$gatk" ModelSegments \
-O "$results_dir"'/CNV' \
--output-prefix "$i" \
--allelic-counts "$results_dir"'/CNV/pre_processing/'"$i"'_clean.allelicCounts.tsv' \
--denoised-copy-ratios "$results_dir"'/CNV/pre_processing/'"$i"'_clean.denoisedCR.tsv' \
--tmp-dir "$tmp_dir"

msg="Call copy-neutral, amplified and deleted segments with CallCopyRatioSegments: $cell_line"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
"$JAVA_DIR" -Xms16g -jar "$gatk" CallCopyRatioSegments \
-I "$results_dir"'/CNV/'"$i"'.cr.seg' \
-O "$results_dir"'/CNV/'"$i"'.called.seg' \
--neutral-segment-copy-ratio-lower-bound 0.7 \
--neutral-segment-copy-ratio-upper-bound 1.3 \
--tmp-dir "$tmp_dir"

msg="Plot modeled copy ratio and allelic fraction segments with PlotModeledSegments: $cell_line"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"
"$JAVA_DIR" -Xms16g -jar "$gatk" PlotModeledSegments \
-O "$results_dir"'/CNV' \
--output-prefix "$i" \
--segments "$results_dir"'/CNV/'"$i"'.modelFinal.seg' \
--sequence-dictionary "$dict" \
--allelic-counts "$results_dir"'/CNV/'"$i"'.hets.tsv' \
--denoised-copy-ratios "$results_dir"'/CNV/pre_processing/'"$i"'_clean.denoisedCR.tsv' \
--minimum-contig-length 46709983 \
--tmp-dir "$tmp_dir"

msg="finished"; echo "-- $msg $longLine"; >&2 echo "-- $msg $longLine"

done
