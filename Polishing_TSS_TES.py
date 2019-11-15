###############################################
##Dmitry Sutormin, 2019##
##RNA-Seq analysis##

#Polish data on E. coli TUs using UTR data from regulonDB (MG1655 U00096.3) [downloaded on 16.09.2019].
###############################################

#######
#Packages to be imported.
#######

import os
from os import listdir
import numpy as np
import scipy
import pandas as pd
from pandas import DataFrame
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as sch
from matplotlib import cm as cm
from Bio import SeqIO, SeqUtils
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastnCommandline

#PWD
PWD="F:\TU_starts_ends\\"

#RegulonDB UTR data
Regulon_UTR="F:\RegulonDB_E_coli\\UTR_5_3_sequence_processed.txt"

#MG1655 U00096.3 sequence
MG1655="C:\Sutor\science\DNA-gyrase\Genomes\\E_coli_K12_MG1655_U00096.3.fasta"

#W3110 MuSGS genome.
W3110=PWD+"E_coli_w3110_MuSGS_blast_db\\E_coli_w3110_G_Mu.fasta"

#Path to TUs sequences.
UTRs_sequences_path=PWD+"UTRs_full_sequences.fasta"

#Blast results.
Blast_results=PWD+"UTRs_full_sequences_W3110.blastResult"

#Deletions.
Deletions="C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\\Deletions_w3110_G_Mu_SGS.broadPeak"


#######
#Checks if directory exists and if not - creates.
#######

def Dir_check_create(some_path):
    if not os.path.exists(some_path):
        os.makedirs(some_path)    
    return


#######
#Deletions parser (BED).
#######

def Deletions_parser(Deletions_inpath):
    Deletions_in=open(Deletions_inpath, 'r')
    Deletions_coord=[]
    for line in Deletions_in:
        line=line.rstrip().split('\t')
        Deletions_coord.append([int(line[1])+1, int(line[2])]) #[TAD start, TAD end]
    Deletions_in.close()
    return Deletions_coord

Deletions_data=Deletions_parser(Deletions)

#######
#Read table with UTR data.
#######

def read_UTR(utr_path):
    #UTR data is +1-based (according to RegulonDB description).
    UTRs=pd.read_csv(utr_path, sep='\t')
    UTRs_coordinates_list=UTRs['UTR_coordinates'].tolist()
    TSS=[]
    TES=[]
    for UTR in UTRs_coordinates_list:
        UTR=UTR.split('-')
        TSS.append(UTR[0])
        TES.append(UTR[1])
    UTRs['TSS']=TSS
    UTRs['TES']=TES
    #print(UTRs)           
    return UTRs

UTRs_info=read_UTR(Regulon_UTR)


#######
#Genome sequence parsing.
#######

def genome_seq(genome_path):
    genome=open(genome_path, 'r')
    for record in SeqIO.parse(genome, "fasta"):
        genome_sequence=str(record.seq)
    genome.close()
    print('Whole genome average GC: ' + str(SeqUtils.GC(genome_sequence)))
    print('Whole genome length: ' + str(len(genome_sequence)))        
    return genome_sequence


#######
#Return TUs sequences.
#######

def return_TU_seqs(MG1655_genome_path, UTRs_seqs_path, UTRs):
    MG_genome=genome_seq(MG1655_genome_path)
    UTR_seq_file=open(UTRs_seqs_path, 'w')
    TSS=UTRs['TSS'].tolist()
    TES=UTRs['TES'].tolist()
    Strands=UTRs['Strand'].tolist()
    Names=UTRs['TU_name'].tolist()
    UTRs_seqs=[]
    for i in range(len(TSS)):
        seq=MG_genome[int(TSS[i])-1:int(TES[i])]
        UTRs_seqs.append(seq)
        UTR_seq_file.write(f'>{int(TSS[i])}_{TES[i]}_{Strands[i]}_{Names[i]}\n{seq}\n')
    #print(UTRs_seqs[2])
    UTRs['UTR_sequence']=UTRs_seqs
    return UTRs, UTRs_seqs

UTRs_info_seq, UTRs_sequences=return_TU_seqs(MG1655, UTRs_sequences_path, UTRs_info)


#######
#Blast of UTRs performd via command line, because python-implemented blast is too slow.
#blastn -query F:\TU_starts_ends\UTRs_full_sequences.fasta -subject F:\TU_starts_ends\E_coli_w3110_MuSGS_blast_db\E_coli_w3110_G_Mu.fasta -outfmt=6 -out F:\TU_starts_ends\UTRs_full_sequences_W3110.blastResult -max_target_seqs 1 -max_hsps 1
#Missed yffOP TU was manually added to blast results.
#######

def parse_blast_results(Blast_results_path, UTRs):
    blast_file=open(Blast_results_path, 'r')
    UTR_w3110_starts=[]
    UTR_w3110_ends=[]
    Reversal=[]
    for line in blast_file:
        line=line.rstrip().split('\t')
        query_name=line[0]
        if line[8]!='NaN':
            if int(line[8])<int(line[9]):
                UTR_w3110_start=int(line[8])
                UTR_w3110_end=int(line[9])
                Reversal.append(0)
            elif int(line[8])>int(line[9]):
                UTR_w3110_start=int(line[9])
                UTR_w3110_end=int(line[8])   
                Reversal.append(1)
            UTR_w3110_starts.append(UTR_w3110_start)
            UTR_w3110_ends.append(UTR_w3110_end)
        else:
            UTR_w3110_starts.append('NaN')
            UTR_w3110_ends.append('NaN')
            Reversal.append(2)
    UTRs['W3110_start']=UTR_w3110_starts
    UTRs['W3110_ends']=UTR_w3110_ends
    UTRs_strand=UTRs['Strand'].tolist()
    UTRs_w3110_strand=[]
    for i in range(len(UTRs_strand)):
        if Reversal[i]==0:
            if UTRs_strand[i]=='forward':
                UTRs_w3110_strand.append('+')
            elif UTRs_strand[i]=='reverse':
                UTRs_w3110_strand.append('-')
        elif Reversal[i]==1:
            if UTRs_strand[i]=='forward':
                UTRs_w3110_strand.append('-')
            elif UTRs_strand[i]=='reverse':
                UTRs_w3110_strand.append('+')
        elif Reversal[i]==2:
            if UTRs_strand[i]=='forward':
                UTRs_w3110_strand.append('NaN')
            elif UTRs_strand[i]=='reverse':
                UTRs_w3110_strand.append('NaN')
    UTRs['W3110_strand']=UTRs_w3110_strand
    #print(UTRs)
    return UTRs

UTRs_w3110_info=parse_blast_results(Blast_results, UTRs_info_seq)


#######
#Write BED file for RSeQC.
#######

def write_bed_and_broadPeak(UTRs, deletions_data, pwd):
    print(UTRs.shape)
    UTRs=UTRs[UTRs['W3110_start']!='NaN']
    print(UTRs.shape)
    
    UTRs_bed=pd.DataFrame()
    UTRs_bed['chrom']=['NC_007779.1_w3110_Mu']*UTRs.shape[0] #1
    UTRs_bed['chromStart']=list(np.asarray(UTRs['W3110_start'].tolist())-1) #2
    UTRs_bed['chromEnd']=list(np.asarray(UTRs['W3110_ends'].tolist())) #3
    UTRs_bed['name']=UTRs['TU_name'].tolist() #4
    UTRs_bed['score']=[100]*UTRs.shape[0] #5
    UTRs_bed['strand']=UTRs['W3110_strand'].tolist() #6
    UTRs_bed['thickStart']=list(np.asarray(UTRs['W3110_start'].tolist())-1) #7
    UTRs_bed['thickEnd']=list(np.asarray(UTRs['W3110_ends'].tolist())) #8
    UTRs_bed['itemRgb']=['255,0,0']*UTRs.shape[0] #9
    UTRs_bed['blockCount']=['1']*UTRs.shape[0] #10
    UTRs_bed['blockSizes']=list(np.asarray(UTRs['W3110_ends'].tolist())-np.asarray(UTRs['W3110_start'].tolist())+1) #11
    UTRs_bed['blockStarts']=['0,']*UTRs.shape[0] #12
    UTRs_bed_del_cor=UTRs_bed
    for deletion in deletions_data:
        UTRs_bed_del_cor=UTRs_bed_del_cor[~((deletion[1]>=UTRs_bed_del_cor['chromStart']) & (UTRs_bed_del_cor['chromStart']>=deletion[0])) & ~((deletion[1]>=UTRs_bed_del_cor['chromEnd']) & (UTRs_bed_del_cor['chromEnd']>=deletion[0]))]
    UTRs_bed_del_cor.to_csv(pwd+"UTRs_w3110_del_cor.bed", sep='\t', header=False, index=False)
    UTRs_bed.to_csv(pwd+"UTRs_w3110.bed", sep='\t', header=False, index=False)
    
    UTRs_broadpeak=pd.DataFrame()
    UTRs_broadpeak['chrom']=['NC_007779.1_w3110_Mu']*UTRs.shape[0] #1
    UTRs_broadpeak['chromStart']=list(np.asarray(UTRs['W3110_start'].tolist())-1) #2
    UTRs_broadpeak['chromEnd']=list(np.asarray(UTRs['W3110_ends'].tolist())) #3
    UTRs_broadpeak['name']=UTRs['TU_name'].tolist() #4
    UTRs_broadpeak['score']=[100]*UTRs.shape[0] #5
    UTRs_broadpeak['strand']=UTRs['W3110_strand'].tolist() #6   
    UTRs_broadpeak['signalValue']=[100]*UTRs.shape[0] #5
    UTRs_broadpeak['pValue']=['-1']*UTRs.shape[0] #5
    UTRs_broadpeak['qValue']=['-1']*UTRs.shape[0] #5
    UTRs_broadpeak_del_cor=UTRs_broadpeak
    for deletion in deletions_data:
        UTRs_broadpeak_del_cor=UTRs_broadpeak_del_cor[~((deletion[1]>=UTRs_broadpeak_del_cor['chromStart']) & (UTRs_broadpeak_del_cor['chromStart']>=deletion[0])) & ~((deletion[1]>=UTRs_broadpeak_del_cor['chromEnd']) & (UTRs_broadpeak_del_cor['chromEnd']>=deletion[0]))]    
    UTRs_broadpeak_del_cor.to_csv(pwd+"UTRs_w3110_del_cor.broadPeak", sep='\t', header=False, index=False)
    UTRs_broadpeak.to_csv(pwd+"UTRs_w3110.broadPeak", sep='\t', header=False, index=False)
    return

write_bed_and_broadPeak(UTRs_w3110_info, Deletions_data, PWD)