###############################################
##Dmitry Sutormin, 2020##
##RNA-Seq analysis##

#Takes bed file with expression level assigned to genomic intervals.
#Calculate a histogram representing a range of transcription levels.
###############################################

#######
#Packages to be imported.
#######

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm as cm


#PWD
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_RNA-Seq\Expression_data\DY330_transcripts\Representative_transcripts\\"

#Input data.
Expression_data=PWD + 'DY330_RNA-Seq_transcripts_representative_EP_del_cor.txt'

#Output folder.
Output_folder=PWD + 'Figures\\'



#######
#Reads annotation of particular set of genes .tab BroadPeak-like (determined on a basis of expression level).
#######

def parse_expression_annotation(annot_inpath):
    genes_length=[]
    genes_expression=[]
    filein=open(annot_inpath, 'r')
    for line in filein:
        line=line.rstrip().split('\t')
        if line[0] not in ['GeneID', 'OperonID', 'TU_ID']:
            TU_name=line[1].lstrip('"').rstrip(';"')
            TU_start=int(line[2])
            TU_end=int(line[3])
            TU_strand=line[4]
            TU_expression=float(line[5].replace(',','.'))
            genes_length.append(TU_end-TU_start)
            genes_expression.append(TU_expression)
    filein.close()            
    return genes_length, genes_expression


#######
#Plot distribution of TUs length and expression level.
#######

def plot_distribution(annot_inpath, output_path):
    #Read expression data.
    genes_length, genes_expression = parse_expression_annotation(annot_inpath)
    
    #Plot data.
    fig=plt.figure(figsize=(5,7), dpi=100)
    #Plot distribution of TUs length.
    plt1=fig.add_subplot(3,1,1)     
    plt1.hist(genes_length, bins=100, edgecolor='black', linewidth=0.5, color='#a4c8ff')
    plt1.set_yscale('log')
    plt1.set_xlabel('Length of TUs', size=15)
    plt1.set_ylabel('Number of TUs, log', size=12)
    
    #Plot distribution of TUs expression.
    plt2=fig.add_subplot(3,1,2)     
    plt2.hist(genes_expression, bins=100, edgecolor='black', linewidth=0.5, color='#ffa6d1')
    plt2.set_yscale('log')
    plt2.set_xlabel('Expression level of TUs', size=15)
    plt2.set_ylabel('Number of TUs, log', size=12)
    
    #Plot distribution of TUs expression, zoom-in.
    plt2=fig.add_subplot(3,1,3)     
    plt2.hist(genes_expression, bins=10000, edgecolor='black', linewidth=0.5, color='#ffa6d1')
    plt2.set_yscale('log')
    plt2.set_xlim([0, 1000])
    plt2.set_xlabel('Expression level of TUs', size=15)
    plt2.set_ylabel('Number of TUs, log', size=12)    
    
    plt.tight_layout()
    plt.show()
    plt.savefig(output_path+'Distribution_of_TUs_length_and_exression_level.png', dpi=300)    
    
    return


plot_distribution(Expression_data, Output_folder)