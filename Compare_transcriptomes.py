###############################################
##Dmitry Sutormin, 2019##
##RNA-Seq analysis##

#Takes RPKM data of different RNA-Seq experiments, 
#combines them in a dataframe and performes corralation analysis.
#Returns xlsx tables with merged data.
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


#PWD
PWD="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_RNA-Seq\\"
#Old expression data.
Old_expression_input="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_RNA-Seq\Expression_data\Incarnato_data\DOOR_Mu_del_cor_genes_expression.bed"
#tRNA list.
tRNA_list_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_RNA-Seq\Expression_data\DY330_genes\\tRNA_genes_list.txt"
#rRNA list.
rRNA_list_path="C:\\Users\sutor\OneDrive\ThinkPad_working\Sutor\Science\E_coli_RNA-Seq\Expression_data\DY330_genes\\rRNA_genes_list.txt"


#######
#Iterate RNA-Seq data files and merge them.
#######

def Open_merge_data(input_path, output_path):
    files=listdir(input_path)
    datasets_merged=pd.read_csv(input_path+files[0], sep='\t')
    for file in files[1:]:
        if file[-3:]=="xls":
            print(file)
            dataset=pd.read_csv(input_path+file, sep='\t')
            datasets_merged=pd.merge(datasets_merged, dataset, how='right', on=['#chrom', 'st', 'end', 'accession', 'mRNA_size', 'gene_strand']) #on=['chrom', 'st', 'end', 'accession', 'mRNA_size', 'gene_strand']
            print(datasets_merged.shape) 
    datasets_merged['end']=np.array(datasets_merged['end'].tolist())-1  #Be carefull! Should be removed! Due to bug in annotation!
    datasets_merged.to_excel(output_path+"DY330_RNA-Seq_data_merged.xlsx")
    print(datasets_merged.head())
    datasets_merged.to_excel(output_path+"DY330_RNA-Seq_data_merged_fragment_counts_columns.xlsx", columns=['#chrom', 'st', 'end', 'accession', 'mRNA_size', 'gene_strand', 
                                                                                                          'Frag_count_EP1', 'Frag_count_EP2', 'Frag_count_EP3',                                                                                          
                                                                                                          'Frag_count_ESP1', 'Frag_count_ESP2', 'Frag_count_ESP3',
                                                                                                          'Frag_count_SP1', 'Frag_count_SP2', 'Frag_count_SP3'], index=False)
                                                                                              
    datasets_merged.to_excel(output_path+"DY330_RNA-Seq_data_merged_fragment_counts_columns_for_edgeR.xlsx", columns=['accession', 'Frag_count_EP1', 'Frag_count_EP2', 'Frag_count_EP3',                                                                                          
                                                                                                                    'Frag_count_ESP1', 'Frag_count_ESP2', 'Frag_count_ESP3',
                                                                                                                    'Frag_count_SP1', 'Frag_count_SP2', 'Frag_count_SP3'], index=False)                                                                                         
    datasets_merged.to_excel(output_path+"DY330_RNA-Seq_data_merged_FPKM_columns.xlsx", columns=['#chrom', 'st', 'end', 'accession', 'mRNA_size', 'gene_strand',
                                                                                               'FPKM_EP1', 'FPKM_EP2', 'FPKM_EP3',                                                                                          
                                                                                               'FPKM_ESP1', 'FPKM_ESP2', 'FPKM_ESP3',
                                                                                               'FPKM_SP1', 'FPKM_SP2', 'FPKM_SP3'], index=False)      
    return datasets_merged

#######
#Add old expression data.
#######

def add_old_expression(datasets_merged, old_expression_data_path, output_path):
    print("DY330 RNA-Seq data", datasets_merged.head())
    Old_RNA_seq_data=pd.read_csv(old_expression_data_path, header=None, sep='\t')
    #print(Old_RNA_seq_data)
    Old_RNA_seq_data_named=pd.DataFrame()
    Old_RNA_seq_data_named['#chrom']=Old_RNA_seq_data[0]
    Old_RNA_seq_data_named['st']=Old_RNA_seq_data[1]
    Old_RNA_seq_data_named['end']=Old_RNA_seq_data[2]
    Old_RNA_seq_data_named['accession']=Old_RNA_seq_data[3]
    Old_RNA_seq_data_named['FPKM_EP_old']=Old_RNA_seq_data[4]
    Old_RNA_seq_data_named['gene_strand']=Old_RNA_seq_data[5]
    print("W3110 RNA-Seq data", Old_RNA_seq_data_named.head())
    
    datasets_new_and_old=pd.merge(datasets_merged, Old_RNA_seq_data_named, how='right', on=['#chrom', 'st', 'end', 'accession', 'gene_strand'])
    print("DY330+W3110 RNA-Seq data", datasets_new_and_old.head())
    datasets_new_and_old.to_excel(output_path+"All_RNA-Seq_data_merged_FPKM_columns.xlsx", columns=['#chrom', 'st', 'end', 'accession', 'mRNA_size', 'gene_strand',
                                                                                                    'FPKM_EP1', 'FPKM_EP2', 'FPKM_EP3',                                                                                          
                                                                                                    'FPKM_ESP1', 'FPKM_ESP2', 'FPKM_ESP3',
                                                                                                    'FPKM_SP1', 'FPKM_SP2', 'FPKM_SP3', 'FPKM_EP_old'], index=False)    
    return


#########
##Compute correlation matrix and draw heatmaps.
#########

#Plot diagonal correlation matrix.
def correlation_matrix(df, cor_method, title, outpath):
    fig=plt.figure(figsize=(6,6), dpi=100)
    ax1=fig.add_subplot(111)
    cmap=cm.get_cmap('rainbow', 30)
    #print(df)
    df_cor_matrix=df.corr(method=cor_method)
    cax=ax1.imshow(df_cor_matrix, interpolation="nearest", cmap=cmap, norm=None, vmin=0.9, vmax=1)
    ax1.grid(True, which='minor', linestyle="--", linewidth=0.5, color="black")
    plt.title(title)
    labels=list(df)
    ax1.set_xticks(np.arange(len(labels)))
    ax1.set_yticks(np.arange(len(labels)))    
    ax1.set_xticklabels(labels, fontsize=12, rotation=90)
    ax1.set_yticklabels(labels, fontsize=12)
    #Create text annotation for heatmap pixels.
    for i in range(len(labels)):
        for j in range(len(labels)):
            text = ax1.text(i, j, round(df_cor_matrix[labels[i]][labels[j]], 3), ha="center", va="center", color="black")        
    #Add colorbar, make sure to specify tick locations to match desired ticklabels.
    #Full scale:[-1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55, -0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00])
    fig.colorbar(cax, ticks=[-1.00, -0.90, -0.80, -0.70, -0.60, -0.50, -0.40, -0.30, -0.20, -0.10, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00], shrink=0.7)
    plt.tight_layout()
    plt.savefig(outpath, dpi=400, figsize=(6, 6))
    plt.show()
    plt.close()
    return df_cor_matrix

#######
#Identify clusters in a corralation matrix (hierarchy clustering).
#Code stolen from https://github.com/TheLoneNut/CorrelationMatrixClustering/blob/master/CorrelationMatrixClustering.ipynb
#######

def Clustering(df):
    X = df.corr().values
    d = sch.distance.pdist(X)   # vector of pairwise distances
    L = sch.linkage(d, method='complete')
    ind = sch.fcluster(L, 0.5*d.max(), 'distance')
    columns = [df.columns.tolist()[i] for i in list((np.argsort(ind)))]
    df = df.reindex(columns, axis=1)
    return df


#######
#Read tRNA, rRNA lists.
#######

def Read_RNA_list(list_path):
    RNA_file=open(list_path, 'r')
    RNA_list=[]
    for line in RNA_file:
        line=line.rstrip()
        RNA_list.append(line)
    print(f'Number of genes detected: {len(RNA_list)}')    
    return RNA_list


#######
#Remove tRNA, rRNA TUs.
#######

def Remove_tRNA_rRNA(input_path, tRNA_path, rRNA_path, output_path_no_tRNA, output_path_no_rRNA, output_path_no_tRNA_rRNA):
    datasets_merged=pd.read_excel(input_path)
    print(f'Initial dataframe shape: {datasets_merged.shape}')
    #Obtain tRNA and rRNA genes names.
    tRNA_list=Read_RNA_list(tRNA_path)
    rRNA_list=Read_RNA_list(rRNA_path)
    #Remove tRNA, rRNA and both.
    datasets_merged_no_tRNA=datasets_merged.loc[~datasets_merged['accession'].isin(tRNA_list)]
    print(f'Dataframe shape after tRNA is discarded: {datasets_merged_no_tRNA.shape}')
    datasets_merged_no_rRNA=datasets_merged.loc[~datasets_merged['accession'].isin(rRNA_list)]
    print(f'Dataframe shape after rRNA is discarded: {datasets_merged_no_rRNA.shape}')    
    datasets_merged_no_rRNA_tRNA=datasets_merged_no_rRNA.loc[~datasets_merged_no_rRNA['accession'].isin(tRNA_list)]
    print(f'Dataframe shape after rRNA and tRNA is discarded: {datasets_merged_no_rRNA_tRNA.shape}')      
    datasets_merged_no_tRNA.to_excel(output_path_no_tRNA, index=False)
    datasets_merged_no_rRNA.to_excel(output_path_no_rRNA, index=False)
    datasets_merged_no_rRNA_tRNA.to_excel(output_path_no_tRNA_rRNA, index=False)
    return datasets_merged_no_tRNA, datasets_merged_no_rRNA, datasets_merged_no_rRNA_tRNA


#######
#Read replicas data in DF and compare.
#######

def Read_DF_compare_transcriptomes(input_path, columns_range, output_path, output_path_clusterized):
    datasets_merged=pd.read_excel(input_path, sheet_name='Sheet1', usecols=columns_range)
    #print(datasets_merged)
    correlation_matrix(datasets_merged, 'pearson', 'Correlation of RNA-Seq data', output_path)
    datasets_merged_clusterized=Clustering(datasets_merged)
    correlation_matrix(datasets_merged_clusterized, 'pearson', 'Correlation of RNA-Seq data clusterized', output_path_clusterized)    
    return


#######
#Average biological replicas, return bed and broadPeak files.
#######

def average_replicas_write_door_like_tab(RNA_seq_data_full, Name_ID, data_output_path_EP, data_output_path_ESP, data_output_path_SP):
    E_dataframe=RNA_seq_data_full[['FPKM_EP1', 'FPKM_EP2', 'FPKM_EP3']]
    E_mean=E_dataframe.mean(axis=1).tolist()
    ES_dataframe=RNA_seq_data_full[['FPKM_ESP1', 'FPKM_ESP2', 'FPKM_ESP3']]
    ES_mean=ES_dataframe.mean(axis=1).tolist()
    S_dataframe=RNA_seq_data_full[['FPKM_SP1', 'FPKM_SP2', 'FPKM_SP3']]
    S_mean=S_dataframe.mean(axis=1).tolist()    
    
    Chromosome_ID='NC_007779.1_w3110_Mu'
    
    SES_data=pd.DataFrame()
    SES_data[Name_ID]=range(0, len(S_mean), 1)
    SES_data['Gene_name']=RNA_seq_data_full['accession'].tolist()
    SES_data['Start']=np.array(RNA_seq_data_full['st'].tolist())+1 #Conversion from 0-based bed format to 1-based gff-like
    SES_data['End']=np.array(RNA_seq_data_full['end'].tolist())
    SES_data['Strand']=RNA_seq_data_full['gene_strand'].tolist()
    SES_data['Expression_E']=E_mean
    SES_data['Expression_ES']=ES_mean
    SES_data['Expression_S']=S_mean
    SES_data['Description']=RNA_seq_data_full['accession'].tolist()
    
    E_data=SES_data.drop(['Expression_ES', 'Expression_S'], axis=1)
    ES_data=SES_data.drop(['Expression_E', 'Expression_S'], axis=1)
    S_data=SES_data.drop(['Expression_E', 'Expression_ES'], axis=1)
    
    E_data.to_csv(data_output_path_EP, index=False, header=True, sep='\t')
    ES_data.to_csv(data_output_path_ESP, index=False, header=True, sep='\t')
    S_data.to_csv(data_output_path_SP, index=False, header=True, sep='\t')
    return


#######
#Wrapper function.
#######

def wrapper(data_inpath, data_output_path, old_expression_input, tRNA_inpath, rRNA_inpath, Name_ID):
    #DY330_dataset_merged=Open_merge_data(data_inpath, data_output_path)
    #add_old_expression(DY330_dataset_merged, old_expression_input, data_output_path)
    #
    #datasets_merged_no_tRNA, datasets_merged_no_rRNA, datasets_merged_no_rRNA_tRNA=Remove_tRNA_rRNA(data_output_path+"DY330_RNA-Seq_data_merged_FPKM_columns.xlsx",
    #                                                                                                tRNA_inpath, rRNA_inpath, 
    #                                                                                                data_output_path+"DY330_RNA-Seq_data_merged_FPKM_columns_no_tRNA.xlsx",
    #                                                                                                data_output_path+"DY330_RNA-Seq_data_merged_FPKM_columns_no_rRNA.xlsx",
    #                                                                                                data_output_path+"DY330_RNA-Seq_data_merged_FPKM_columns_no_tRNA_rRNA.xlsx")   
    #
    #Remove_tRNA_rRNA(data_output_path+"All_RNA-Seq_data_merged_FPKM_columns.xlsx",
    #                 tRNA_inpath, rRNA_inpath, 
    #                 data_output_path+"All_RNA-Seq_data_merged_FPKM_columns_no_tRNA.xlsx",
    #                 data_output_path+"All_RNA-Seq_data_merged_FPKM_columns_no_rRNA.xlsx",
    #                 data_output_path+"All_RNA-Seq_data_merged_FPKM_columns_no_tRNA_rRNA.xlsx")    
    #
    #Read_DF_compare_transcriptomes(data_output_path+"DY330_RNA-Seq_data_merged_FPKM_columns.xlsx", 'G:O', 
    #                               data_output_path+"Figures\DY330_RNA-Seq_data_merged_FPKM_correlation.png",
    #                               data_output_path+"Figures\DY330_RNA-Seq_data_merged_counts_correlation_clusterized.png")
    #
    #Read_DF_compare_transcriptomes(data_output_path+"All_RNA-Seq_data_merged_FPKM_columns.xlsx", 'G:P', 
    #                               data_output_path+"Figures\All_RNA-Seq_data_merged_FPKM_correlation.png",
    #                               data_output_path+"Figures\All_RNA-Seq_data_merged_counts_correlation_clusterized.png")
    #
    #Read_DF_compare_transcriptomes(data_output_path+"DY330_RNA-Seq_data_merged_FPKM_columns_no_tRNA_rRNA.xlsx", 'G:O', 
    #                               data_output_path+"Figures\DY330_RNA-Seq_data_merged_FPKM_correlation_no_tRNA_rRNA.png",
    #                               data_output_path+"Figures\DY330_RNA-Seq_data_merged_counts_correlation_clusterized_no_tRNA_rRNA.png")
    #
    #Read_DF_compare_transcriptomes(data_output_path+"All_RNA-Seq_data_merged_FPKM_columns_no_tRNA_rRNA.xlsx", 'G:P', 
    #                               data_output_path+"Figures\All_RNA-Seq_data_merged_FPKM_correlation_no_tRNA_rRNA.png",
    #                               data_output_path+"Figures\All_RNA-Seq_data_merged_counts_correlation_clusterized_no_tRNA_rRNA.png") 
    
    Read_DF_compare_transcriptomes(data_output_path+"DY330_RNA-Seq_data_merged_FPKM_columns_no_tRNA_rRNA.xlsx", 'G:I', 
                                   data_output_path+"Figures\EP_RNA-Seq_data_merged_FPKM_correlation_no_tRNA_rRNA.png",
                                   data_output_path+"Figures\EP_RNA-Seq_data_merged_counts_correlation_clusterized_no_tRNA_rRNA.png")     
    
    #average_replicas_write_door_like_tab(DY330_dataset_merged, Name_ID, data_output_path+"DY330_RNA-Seq_genes_EP_del_cor.txt", 
    #                                     data_output_path+"DY330_RNA-Seq_genes_ESP_del_cor.txt", 
    #                                     data_output_path+"DY330_RNA-Seq_genes_SP_del_cor.txt")
    #
    #average_replicas_write_door_like_tab(datasets_merged_no_rRNA_tRNA, Name_ID, data_output_path+"DY330_RNA-Seq_genes_no_tRNA_rRNA_EP_del_cor.txt", 
    #                                     data_output_path+"DY330_RNA-Seq_genes_no_tRNA_rRNA_ESP_del_cor.txt", 
    #                                     data_output_path+"DY330_RNA-Seq_genes_no_tRNA_rRNA_SP_del_cor.txt")    
    
    return

wrapper(PWD+"FPKM\\Genes_del_cor\\", PWD+"FPKM_analysis_correct\\UTRs_del_cor\\", Old_expression_input, tRNA_list_path, rRNA_list_path, "GeneID")

