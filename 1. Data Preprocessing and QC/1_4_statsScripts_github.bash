
# To calculate the number of Unique Mapped Reads for each paired ended and deduplicated single cell bam file, the following script was used (as per following reference).
# http://qnot.org/2012/04/14/counting-the-number-of-reads-in-a-bam-file/
# -c = count, -f 1 = only reads which are paired in sequencing, -F 12 means to include all reads where neither flag 0x0004 or 0x0008 is set, where 0x0004 is not unmapped reads and 0x0008 is where the mate is not unmapped (only include reads where it maps and its mate also map). 
for i in {1..192}
do
samtools view -c -f 1 -F 12 /new/howard_new/sequence_runs/pathToBamFiles/S${i}_PE_rmdup1.bam
done


#For Fraction of Reads in Peaks, macs2 was used to define the peak regions with the following parameters
-g hs -f BAMPE --nomodel -B -p 5e-2 --min-length 500 --max-gap 400 --SPMR --call-summits


#For each paired ended and deduplicated single cell bam file, a tagAlign bed file was created, intersected with reference macs2 narrowPeak file, then then the number of intersections were counted
for i in {1..192}
do
bedtools bamtobed -i /new/howard_new/sequence_runs/pathToBamFiles/S${i}_PE_rmdup1.bam | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' > /new/howard_new/sequence_runs/pathToBamFiles/S${i}_PE_rmdup1.tagAlign
bedtools sort -i /mnt/raid5/cutnrun/userFiles/howard/macs2normalised/120pri_k27ac_norm_sort_peaks.narrowPeak | bedtools merge -i stdin | bedtools intersect -u -a /new/howard_new/sequence_runs/pathToBamFiles/S${i}_PE_rmdup1.tagAlign -b stdin | wc -l
done


#For each single cell tagAlign bed file, count the number of paired end reads
for i in {1..192}
do
grep -c "." /new/howard_new/sequence_runs/pathToBamFiles/S${i}_PE_rmdup1.tagAlign
done


#For each paired ended and deduplicated single cell bam file, a tagAlign bed file was created, intersected with reference Blacklist file, then then the number of intersections were counted
for i in {1..192}
do
bedtools bamtobed -i /new/howard_new/sequence_runs/pathToBamFiles/S${i}_PE_rmdup1.bam | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' > /new/howard_new/sequence_runs/pathToBamFiles/S${i}_PE_rmdup1.tagAlign
bedtools sort -i /mnt/raid5/cutnrun/reference/regulatoryFeatures/blacklist_20211120/hg38_sorted.blacklist.bed | bedtools merge -i stdin | bedtools intersect -u -a /new/howard_new/sequence_runs/pathToBamFiles/S${i}_PE_rmdup1.tagAlign -b stdin | wc -l
done



