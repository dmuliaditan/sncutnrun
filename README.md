# Single-nucleus CUT&RUN deciphers the function of intrinsic and genomics-driven epigenetic heterogeneity in head and neck cancer progression
Repository for all the code used in the snCUT&amp;RUN project

Daniel Muliaditan, 03 October 2024

Table of Content:
1. Data Preprocessing<br/>
	1.1. Preprocessing individual single-cell .fastq to single-cell .bam<br/>
	1.2. Aggregating single-cell .bam to pseudobulk .bam, pseudobulk .bam normalization, peak calling on pseudobulk .bam with MACS2<br/>
	1.3. Calculating Uniquely Mapped Reads and Fraction Reads in Peaks (FRiP) and Fraction Reads in Blacklist (FRiB) for QC purposes, and getting analysis ready fragment file for Signac<br/>
  	1.4. Comparison between snCUT&RUN and other methods (scCUT&Tag, Droplet Paired-Tag, CoTECH and uliCUT&RUN)

3. Signac analysis<br/>
	2.1. Loading necessary packages and creating a chromatin assay<br/>
	2.2. UMAP dimension reduction and clustering<br/>
	2.3. Peak calling and coverage plots<br/>
	2.4. Transcription Factor Motif analysis<br/>
	2.5. Differential peak calling and TF motif analysis

4. Genetic Driven Epigenetic Reprogramming<br/>
	3.1. Chromatin state annotation with ChromHMM<br/>
	3.2. RNAseq analysis with salmon, tximport and DESeq2<br/>
	3.3. CNV analysis from WES data with CNVkit <br/>
  	3.4. Static correlation between gene CN and gene expression<br/>
  	3.5. Static correlation between chromatin state and gene expression<br/>
  	3.6. Static correlation between gene CN, chromatin state, and gene expression<br/>
 	3.7. Plotting an alluvial plot to visualize global changes in chromatin state from a primary to metastatic tumour<br/>
 	3.8. Analysis of up- and downregulated genes between a primary and metastatic tumour, stratifying genes by chromatin state change<br/>
  	3.9. Analyzing chromatin state changes in genes affected by a change in gene CN

5. Epigenetic Heterogeneity Driven HNSC Progression<br/>
	4.1. Module Score and HN120preMET/HN137prePCR analysis

6. Code used for Supplementary Figures

For questions, please write an e-mail to: Daniel_Muliaditan@gis.a-star.edu.sg

