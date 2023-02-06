# Transcription_and_genome_features
Exploring genome features (binding/cleavage sites) of *E. coli W3110* genome in the context of transcription.

This repository contains a set of bash and python scripts which were used for RNA-Seq data analysis and visualization. 
Raw sequencing data and some processed files can be retrieved from GEO datasets with accession GSE181687.

If you find this code useful and would like to use it in your own research, please, cite:
Sutormin D, Galivondzhyan A, Musharova O, Travin D, Rusanova A, Obraztsova K, Borukhov S, Severinov K. Interaction between transcribing RNA polymerase and topoisomerase I prevents R-loop formation in E. coli. Nat Commun. 2022 Aug 4;13(1):4524. doi: 10.1038/s41467-022-32106-5. PMID: 35927234; PMCID: PMC9352719.


######################

## RNA-Seq_raw_data_analysis.sh

Shell script that makes QC of the reads before and after the trimming procedure. 
Than script maps trimmed reads to the reference genome, prepares sorted and 
indexed BAM-files suitable for visualization with IGV. Removes PCR duplicates and runs RSeQC for 
FPKM and fragments count calculation.

**Requirements:** factqc, trimmomatic, bwa mem, samtools, shell, RSeQC

**Input:** Raw reads files (FASTQ), Genome file (FASTA), annotation of genomic intervals (genes/operons/transcripts) in BED format

**Output:** FastQC reports, SAM files, sorted and indexed BAM files, xls tables with expression data


######################

## Compare_transcriptomes.py

Takes RPKM data of different RNA-Seq experiments, 
combines them in a dataframe and performes corralation analysis.
Returns xlsx tables with merged data.

**Requirements:** python 3

**Input:** xls tables (output of RSeQC), list of rRNA and tRNA genes/operons/transcripts, other RNA-Seq data (optional)

**Output:** correlation heat maps, xlsx files with merged data, averaged between biological replicas expression data (tab files)


######################

## Genome_intervals_analysis_growth_phases.py

The script analyzes sets of genome intervals (transcription units - TUs, BIME-1s, BIME-2s, IHF sites, Fis sites, H-NS sites, MatP sites, chromosomal macrodomains, etc.)
for the enrichment of GCSs (binomial test), compares their N3E and score with mean GCSs N3E and score (t-test), 
compares intervals mean score with genome mean score (t-test).

**Note: ** Returns some warning messages due to the ommitting statistics (t-test) for too short sets of values.

**Requirements:** python 3

**Input:** TAB files with trusted GCSs and score info, WIG genome score file, TAB transcription units data files, TAB intervals data files

**Output:** TAB file with numbers of GCSs are associated with TUs compartments (USUS, USGB, GBDS, DSDS), 
TAB file with GCSs-TUs association analysis (GCSs number, GCSs N3E, GCSs score, TUs compartments score), TAB file with normalized numbers of GCSs are associated with TUs,
TAB file with the number of GCSs are associated with particular intervals set and statistics (GCSs number, GCSs N3E, GCSs score),
TAB file with the number of GCSs are associated with particular intervals (BIMEs-1, BIMEs-2), TAB file with intervals score statistics


## GCSs_and_transcription_correlation_of_growth_phases.py

The script takes sets of trusted GCSs and correlates their density over the genome with transcription signal (binning approach).
Alternatively or additionally the script takes files with GCSs numbers at DS regions (or any other regions) of TUs and correlates 
them with transcription level of the TUs.

**Requirements:** python 3

**Input:** TAB files with trusted GCSs, TAB file with transcription data, TAB files with GCSs numbers calculated for each TU

**Output:** correlation heatmapsx


