#!bin/bash

##############
##Dmitry Sutormin, 2019##
##RNA-Seq analysis##

#Shell script performs QC of the reads before and after the trimming procedure.
#Than script maps trimmed reads to the reference genome, prepares
#sorted and indexed BAM-files suitable for visualization with IGV.
#Removes PCR duplicates and runs RSeQC for FPKM and fragments count calculation.

#Requirements: factqc, trimmomatic, bwa mem, samtools, RSeQC.
#This variables should be in the path (or replace them with the path to the particular program)
##############



#######
#Variables to be defined.
#######

#Path to the working directory, contains /Raw_data folder with raw reads files.
PWD='/home/cls01/Data_Hi-C/E_coli_RNA-Seq_while_Sutor_is_on_vacation'
echo $PWD
cd $PWD

#Path to the file containing sequencing adapters sequences for trimmomatic uses. Typically in the Trimmomatic-0.36/adapters/XXX.fa
Adapters='/home/cls01/Prog/Trimmomatic-0.38/adapters/All_TruSeq.fa'
trimmomatic='/home/cls01/Prog/Trimmomatic-0.38/trimmomatic-0.38.jar'
#Path to the reference genome.
Ref_genome=$PWD/Genome/E_coli_w3110_G_Mu.fasta
Ref_genome_idx=$PWD/Genome/E_coli_w3110_G_Mu.idx



#######
#Quality control and sequencing data preparation.
#######
echo '
#######################
#Initial quality control is in progress...
#######################
'

#Initial quality control
mkdir $PWD/Fastqc_analysis/
mkdir $PWD/Fastqc_analysis/Initial
fastqc -t 48 -o $PWD/Fastqc_analysis/Initial $PWD/RAW_data/*


#######
#Reads trimming
#######
echo '
#######################
#Reads trimming...
#######################
'

mkdir $PWD/Trimmed/
for i in `ls -a $PWD/RAW_data/ | grep 'fastq' | sed -r "s/(.+)\.fastq/\1/g"`; do
echo ${i}.fastq
java -jar $trimmomatic SE -threads 48 -phred33 $PWD/RAW_data/${i}.fastq $PWD/Trimmed/${i}_trimmed.fastq ILLUMINACLIP:$Adapters:2:30:10 LEADING:0 TRAILING:0 SLIDINGWINDOW:4:0 MINLEN:30 ; done


#######
#Quality control after the trimming procedure
#######
echo '
#######################
#Quality control after trimming...
#######################
'

mkdir $PWD/Fastqc_analysis/Trimmed/
fastqc -t 48 -o $PWD/Fastqc_analysis/Trimmed/ $PWD/Trimmed/*



#######
#Prepare index for reference genome.
#######
echo '
#######################
Reference genome indexing...
#######################
'

bwa index $Ref_genome


#######
#Reads mapping, alignment conversion to IGV-compatible format (sorted indexed BAM).
#######
echo '
#######################
Reads mapping, SAM files generation...
#######################
'

#Reads mapping to the reference genome: make SAM-files
mkdir $PWD/SAM/
mkdir $PWD/SAM_sgm/
for i in `ls -a $PWD/Trimmed/ | grep '_trimmed' | sed -r "s/(.+)_trimmed\.fastq/\1/g" | uniq | sort -d`; do
bwa mem -t 48 $Ref_genome $PWD/Trimmed/${i}_trimmed.fastq > $PWD/SAM/$i.sam; done

#Optionally:
#/home/cls01/Prog/segemehl-0.3.4/segemehl.x -t 48 -A 95 -i $Ref_genome_idx -d $Ref_genome -q $PWD/Trimmed/${i}_trimmed.fastq > $PWD/SAM_sgm/$i.sam; done



#######
#Prepares tracks for IGV: makes BAM-files, sorts them, makes index-files
#######
echo '
#######################
BAM files preparation...
#######################
'

mkdir $PWD/BAM_sorted_sgm/
mkdir $PWD/BAM_sgm/
#Makes BAM-files
for i in `ls -a $PWD/SAM/ | grep '.sam' | sed -r "s/(.+).sam/\1/g"`; do
samtools view -S -b $PWD/SAM/${i}.sam > $PWD/BAM/${i}.bam ; done

#Sorts BAM-files
echo '
########################
#BAM files sorting (by position)...
#######################
'

for i in `ls -a $PWD/BAM/ | grep '.bam' | sed -r "s/(.+).bam/\1/g"`; do
samtools sort -@ 40 $PWD/BAM/${i}.bam -o $PWD/BAM_sorted/${i}_sorted.bam ; done


#Converts bam to bed.
echo '
#######################
BAM to bed conversion...
#######################
'

mkdir $PWD/Cov_depth
for i in `ls -a $PWD/BAM_sorted/ | grep '.bam' | sed -r "s/(.+)_sorted.bam/\1/g"`; do
samtools depth -a $PWD/BAM_sorted/${i}_sorted.bam -d 0 > $PWD/Cov_depth/${i}.bed; done


##Removing PCR-duplicates.

#Sorts BAM-files
echo '
#######################
BAM files sorting (by name)...
#######################
'

mkdir BAM_name_sorted_sgm
for i in `ls -a $PWD/BAM/ | grep '.bam' | sed -r "s/(.+).bam/\1/g"`; do
samtools sort -@ 40 -n $PWD/BAM/${i}.bam -o $PWD/BAM_name_sorted/${i}_name_sorted.bam ; done

#Add ms and MC tags for markdup to use later.
echo '
#######################
Add ms and MC tags for markdup to use later...
#######################
'

mkdir BAM_ns_fixmated_sgm
for i in `ls -a $PWD/BAM_name_sorted/ | grep '.bam' | sed -r "s/(.+).bam/\1/g"`; do
samtools fixmate -m $PWD/BAM_name_sorted/${i}.bam $PWD/BAM_ns_fixmated/${i}_fm.bam; done

#Sorts BAM-files
echo '
#######################
BAM files sorting (by position)...
#######################
'

mkdir BAM_fm_position_sorted_sgm
for i in `ls -a $PWD/BAM_ns_fixmated/ | grep '.bam' | sed -r "s/(.+).bam/\1/g"`; do
samtools sort -@ 40 -o $PWD/BAM_fm_position_sorted/${i}_ps.bam $PWD/BAM_ns_fixmated/${i}.bam ; done

#Marks and removes duplicates.
echo '
#######################
Removes duplicates...
#######################
'

mkdir BAM_nodup_sgm
for i in `ls -a $PWD/BAM_fm_position_sorted/ | grep '.bam' | sed -r "s/(.+).bam/\1/g"`; do
samtools markdup $PWD/BAM_fm_position_sorted/${i}.bam $PWD/BAM_nodup/${i}_nd.bam -r ; done


#Converts bam to bed.
echo '
#######################
BAM to bed conversion...
#######################
'

mkdir $PWD/Cov_depth_nodup
for i in `ls -a $PWD/BAM_nodup/ | grep '.bam' | sed -r "s/(.+).bam/\1/g"`; do
samtools depth -a $PWD/BAM_nodup/${i}.bam > $PWD/Cov_depth_nodup/${i}.bed; done



#Makes index files
echo '
#######################
BAM files indexing...
#######################
'

for i in `ls -a $PWD/BAM_sorted/`; do
samtools index $PWD/BAM_sorted/${i} ; done


#Makes index files
echo '
#######################
BAM files indexing...
#######################
'

for i in `ls -a $PWD/BAM_nodup/`; do
samtools index $PWD/BAM_nodup/${i} ; done


echo '
#######################
FPKM calculation...
#######################
'

for i in `ls -a $PWD/BAM_sorted/ | grep '.bam' | sed -r "s/(.+).bam/\1/g"`; do
FPKM_count.py -i $PWD/BAM_sorted/${i}.bam -o $PWD/FPKM/${i}_genes_del_cor -r $PWD/Genome/DOOR_Mu_del_cor_genes_expression.bed -q 25
FPKM_count.py -i $PWD/BAM_sorted/${i}.bam -o $PWD/FPKM/${i}_operons_del_cor -r $PWD/Genome/DOOR_Mu_del_cor_operons_expression.bed -q 25
FPKM_count.py -i $PWD/BAM_sorted/${i}.bam -o $PWD/FPKM/${i}_UTRs_del_cor -r $PWD/Genome/UTRs_w3110_del_cor.bed -q 25; done

for i in `ls -a $PWD/BAM_nodup/ | grep '.bam' | sed -r "s/(.+).bam/\1/g"`; do
FPKM_count.py -i $PWD/BAM_nodup/${i}.bam -o $PWD/FPKM_nodup/${i}_UTRs_nodup -r Genome/UTRs_w3110.bed -q 25 ; done
FPKM_count.py -i $PWD/BAM_nodup/${i}.bam -o $PWD/FPKM/${i}_operons_del_cor_nodup -r $PWD/Genome/DOOR_Mu_del_cor_operons_expression.bed -q 25
FPKM_count.py -i $PWD/BAM_nodup/${i}.bam -o $PWD/FPKM/${i}_UTRs_del_cor_nodup -r $PWD/Genome/UTRs_w3110_del_cor.bed -q 25; done


#Prepare tar.gz archives
echo '
#######################
tar.gz archive preparation...
#######################
'
tar -czvf Fastqc_analysis.tar.gz $PWD/Fastqc_analysis
tar -czvf RNA-Seq_cov_depth.tar.gz $PWD/Cov_depth
tar -czvf RNA-Seq_cov_depth_nodup.tar.gz $PWD/Cov_depth_nodup
tar -czvf FPKM.tar.gz $PWD/FPKM



echo '
Script ended its work succesfully!
'
