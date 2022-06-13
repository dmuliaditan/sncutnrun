# Single nucleus CUT&RUN to assess epigenetic heterogeneity in head and neck cancer progression
Repository for all the code used in the snCUT&amp;RUN project

Daniel Muliaditan, 13 June 2022

This repository contains the following:
- Bash scripts used to preprocess raw single-cell fastq files
- R scripts used for Signac/Seurat analysis

Table of Content:
1. Data Preprocessing<br/>
	1.1. Preprocessing individual single-cell .fastq to single-cell .bam, aggregating single-cell .bam to pseudobulk .bam, 
	peak calling on pseduobulk .bam with MACS2, calculating Fraction Reads in Peaks (FRiP) and Fraction Reads in Blacklist (FRiB) for QC purposes,
	and getting analysis ready fragment file for Signac.

2. Signac analysis<br/>
	2.1. Loading necessary packages and creating a chromatin assay<br/>
	2.2. UMAP dimension reduction and clustering<br/>
	2.3. Peak calling and coverage plots<br/>
	2.4. Transcription Factor Motif analysis<br/>
	2.5. Differential peak calling and TF motif analysis

3. Genetic Driven Epigenetic Reprogramming<br/>
  	3.1. Static correlation between gene CN and gene expression<br/>
  	3.2. Static correlation between chromatin state and gene expression<br/>
  	3.3. Static correlation between gene CN, chromatin state, and gene expression<br/>
 	3.4. Plotting an alluvial plot to visualize global changes in chromatin state from a primary to metastatic tumour<br/>
 	3.5. Analysis of up- and downregulated genes between a primary and metastatic tumour, tratifying genes by chromatin state change<br/>
  	3.6. Analyzing chromatin state changes in genes affected by a change in gene CN

4. Epigenetic Heterogeneity Driven HNSC Progression
	4.1. Module Score and HN120preMET/HN137prePCR analysis

For questions, please write an e-mail to: Daniel_Muliaditan@gis.a-star.edu.sg

