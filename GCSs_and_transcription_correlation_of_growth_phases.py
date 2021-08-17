###############################################
##Dmitry Sutormin, 2019##
##RNA-Seq analysis##

#The script takes sets of trusted GCSs and correlates their density over the genome with transcription signal (binning approach).
#Alternatively or additionally the script takes files with GCSs numbers at DS regions (or any other regions) of TUs and correlates 
#them with transcription level of the TUs.
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib.pyplot as plt
from matplotlib import cm as cm
from matplotlib import gridspec
import numpy as np
from scipy.stats import binom
from scipy.stats import pearsonr
import pandas as pd


#######
#Variables to be defined.
#######

print('Variables to be defined:')

#Input data - GCSs, TAB.
path_to_GCSs_files={'Gyrase exponential RifCfx': "F:\E_coli_Topo-Seqs\GCSs_data\Gyrase_EP_RifCfx_GCSs_trusted_h_s_0.01.txt",
                    'Gyrase exponential Cfx': "F:\E_coli_Topo-Seqs\GCSs_data\Gyrase_EP_Cfx_10mkM_GCSs_trusted_h_s_0.01.txt",
                    'Gyrase early stationary Cfx': "F:\E_coli_Topo-Seqs\GCSs_data\Gyrase_Time-course_-3_min_GCSs_trusted_h_s_0.01.txt",
                    'Gyrase stationary Cfx': "F:\E_coli_Topo-Seqs\GCSs_data\Gyrase_Cfx_stat_trusted_GCSs_h_s_0.01.txt", 
                    'Gyrase exponential Micro': "F:\E_coli_Topo-Seqs\GCSs_data\Gyrase_EP_Micro_GCSs_trusted_h_s_0.01.txt",
                    'Gyrase exponential Oxo': "F:\E_coli_Topo-Seqs\GCSs_data\Gyrase_EP_Oxo_GCSs_trusted_h_s_0.01.txt"}   

#Input data - transcription, TAB.
Transcription_data_files={'EP Expression' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\Expression_data\DY330_genes\DY330_RNA-Seq_genes_EP_del_cor.txt",
                          'ESP Expression' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\Expression_data\DY330_genes\DY330_RNA-Seq_genes_ESP_del_cor.txt",
                          'SP Expression' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\Expression_data\DY330_genes\DY330_RNA-Seq_genes_SP_del_cor.txt",
                          'EP Incarnato Expression' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\Expression_data\Incarnato_data\DOOR_Mu_del_cor_genes_expression.txt"}

GCSs_in_DS_regions_pwd="F:\Topo_data_new_expression\TU_based_analysis\Gyrase\\Correction_on_GCSs_calling_thr\\Genes_del_cor_no_rRNA_rpsl_operon_5000bp\\"
#Input data - GCSs numbers in TUs DS regions.
GCSs_in_DS_regions_files={'Exponential phase GCSs': GCSs_in_DS_regions_pwd + 'Exponential phase_numbers_of_associated_GCSs_DS_only.txt',
                          'Early stationary phase GCSs': GCSs_in_DS_regions_pwd + 'Early stationary phase_numbers_of_associated_GCSs_DS_only.txt',
                          'Stationary phase': GCSs_in_DS_regions_pwd + 'Stationary phase_numbers_of_associated_GCSs_DS_only.txt',
                          'Exponential phase Incarnato' : GCSs_in_DS_regions_pwd + "Exponential phase old_numbers_of_associated_GCSs_DS_only.txt"}

#Path out.
Path_out=GCSs_in_DS_regions_pwd + "Figures\\"
if not os.path.exists(Path_out):
    os.makedirs(Path_out)


#########
#########
#Correlation of binned GCSs number with transcription also binned.
#########
#########

#######
#Trusted GCSs data parsing.
#######

def trusted_GCSs_parsing(input_dict):
    GCSs_sets_dict={}
    for k, v in input_dict.items():
        ar=[]
        filein=open(v, 'r')
        for line in filein:
            line=line.rstrip().split('\t')
            if line[0] not in ['GCSs_coordinate']:
                ar.append(int(line[0]))
            else:
                continue
        GCSs_sets_dict[k]=ar
        print(f'Number of trusted GCSs for {k} : {len(ar)}')
    return GCSs_sets_dict

#######
#Convert score/transcription/GC data to histogram.
#######

def bar_convert(ar, bins):
    bar_ar=[]
    for i in range(len(bins)-1):
        if bins[i]>=0:
            bar_ar.append(sum(ar[int(bins[i]):int(bins[i+1])])/(bins[i+1]-bins[i]))
        else:
            bar_ar.append(sum(ar[int(bins[i]):-1]+ar[0:int(bins[i+1])])/(bins[i+1]-bins[i]))
    return bar_ar

#######
#Parsing transcription - TAB file (optional step).
#######

def transcription_data_parser(transcription_path_dict, bins):
    transcription_dict={}
    for condition, transcription_path in transcription_path_dict.items():
        transcription_file=open(transcription_path, 'r')
        transcription=[]
        for i in range(4647454):
            transcription.append(0)
        for line in transcription_file:
            line=line.rstrip().split('\t')
            if line[0] not in ['GeneID', 'OperonID']:
                for j in range(int(line[3])-int(line[2])):
                    transcription[int(line[2])+j]=float(line[5].replace(',','.'))
        transcription_file.close()
        print(f'Whole genome average transcription {condition}: {sum(transcription)/len(transcription)}')     
        transcription_dict[condition]=bar_convert(transcription, bins)  
    return transcription_dict



#######
#Plotting: distribution throughout the genome, bins are evenly distributed.
#######

def Plot_the_distribution(GSCs_data, bins):
    plot_names=['plot1', 'plot2', 'plot3', 'plot4', 'plot5', 'plot6', 'plot7']    
    #GCSs data plotting.
    fig, plot_names=plt.subplots(len(GSCs_data),1, figsize=(11,(2*len(GSCs_data))+1), dpi=100)
    i=0
    Histo_comp_dict={} #Will contain computed histogramm data (bins and values)
    for key, value in GSCs_data.items():
        Histo_comp_dict[key]=plot_names[i].hist(value, bins) #Plot histo and save computed histogramm data (bins and values)
        i+=1
    plt.close()
    return Histo_comp_dict


#########
##Compute correlation matrix and draw heatmaps.
#########

#Plot diagonal correlation matrix.
def correlation_matrix(df, title, scale_min, scale_max, outpath):
    fig=plt.figure(figsize=(8,8), dpi=100)
    ax1=fig.add_subplot(111)
    cmap=cm.get_cmap('rainbow', 30)
    cax=ax1.imshow(df, interpolation="nearest", cmap=cmap, norm=None, vmin=scale_min, vmax=scale_max)
    ax1.grid(True, which='minor', linestyle="--", linewidth=0.5, color="black")
    plt.title(title)
    labels_x=list(df)
    labels_y=df.index.values.tolist()
    print(labels_x, labels_y)
    ax1.set_xticks(np.arange(len(labels_x)))
    ax1.set_yticks(np.arange(len(labels_y)))    
    ax1.set_xticklabels(labels_x, fontsize=12, rotation=90)
    ax1.set_yticklabels(labels_y, fontsize=12)
    #Add colorbar, make sure to specify tick locations to match desired ticklabels.
    #Full scale:[-1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55, -0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00])
    fig.colorbar(cax, ticks=[-1.00, -0.90, -0.80, -0.70, -0.60, -0.50, -0.40, -0.30, -0.20, -0.10, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00], shrink=0.7)
    plt.tight_layout()
    plt.savefig(outpath, dpi=400, figsize=(8, 8))
    plt.show()
    plt.close()
    return


#######
#For evenly distributed bars.
#Computes following correlation: (Cfx GCSs data vs transcription)
#######

def track_corr(GCSs_histo_comp_dict, Transcription_dict):
    cor_df=pd.DataFrame()
    for condition_transcription, transcription_data in Transcription_dict.items():
        row_dict={}
        for condition_gyrase, gyrase_data in GCSs_histo_comp_dict.items():
            gyrase_data_stat=np.array(gyrase_data[0]).tolist()
            print(f'Paerson correlation ({condition_gyrase}, {condition_transcription}) for even bins: ' + str(pearsonr(transcription_data, gyrase_data_stat))) 
            row_dict[condition_gyrase]=pearsonr(transcription_data, gyrase_data_stat)[0]
        row_series=pd.Series(row_dict, name=condition_transcription)
        cor_df=cor_df.append(row_series)    
    #print(cor_df)
    cor_df=cor_df[['Gyrase exponential Micro','Gyrase exponential Oxo','Gyrase exponential RifCfx','Gyrase exponential Cfx','Gyrase early stationary Cfx','Gyrase stationary Cfx']]
    #print(cor_df)
    return cor_df

#######
#Wrapps all the functions together.
#######

def wrap_the_binning_functions(GCSs_input_dict, transcription_path, path_out, genome_len):
    #Even distribution of bins across the genome.
    bins_number=501
    bins=np.linspace(0, genome_len, bins_number).tolist()
    GSCs_data=trusted_GCSs_parsing(GCSs_input_dict)
    GCSs_histo_comp_even_dict=Plot_the_distribution(GSCs_data, bins)
    Trancription_data_dict=transcription_data_parser(transcription_path, bins)
    cor_df=track_corr(GCSs_histo_comp_even_dict, Trancription_data_dict)
    correlation_matrix(cor_df, 'Correlation between GCSs data and transcription\nfor different phases of growth', -0.05, 0.4, f'{path_out}GCSs_density_transcription_corr_{bins_number-1}_bins_all_conditions.png')
    return

#wrap_the_binning_functions(path_to_GCSs_files, Transcription_data_files, Path_out, 4647454)


#########
#########
#Correlation of GCSs number with genes (TU-based calculation).
#########
#########

#######
#Reads files with GCSs numbers.
#######

def read_merge_GCSs_in_DS(GCSs_in_DS_input, additional_condition, set_name):
    New_expression_data=pd.DataFrame()
    for condition_name, condition_path in GCSs_in_DS_input.items():
        condition=pd.read_csv(condition_path, header=0, sep='\t', index_col=False)
        print(condition_name, '\n')
        condition=condition.drop(labels=[f'{set_name}_ID'], axis=1)     
        #print(condition.head())
        if New_expression_data.shape==(0, 0):
            New_expression_data=condition
        else:
            New_expression_data=pd.merge(New_expression_data, condition, how='left', on=[f'{set_name}_name', 'Start', 'End', 'Strand', 'Gyrase exponential Micro DSDS', 'Gyrase exponential Oxo DSDS', 'Gyrase exponential RifCfx DSDS', 'Gyrase exponential Cfx DSDS', 'Gyrase early stationary Cfx DSDS', 'Gyrase stationary Cfx DSDS'])
    print('Initial dataframe after merging', New_expression_data.shape)
    New_expression_data=New_expression_data.dropna()
    print('After NA dropped', New_expression_data.shape)
    if additional_condition=="no_zeroes":
        New_expression_data=New_expression_data[~((New_expression_data['Gyrase exponential Cfx DSDS']==0) & (New_expression_data['Gyrase early stationary Cfx DSDS']==0) & (New_expression_data['Gyrase stationary Cfx DSDS']==0))]
    elif type(additional_condition)==int:
        New_expression_data=New_expression_data[(New_expression_data['Gyrase exponential Cfx DSDS'] + New_expression_data['Gyrase early stationary Cfx DSDS'] + New_expression_data['Gyrase stationary Cfx DSDS'])>additional_condition]
    print(f'After additional condition {additional_condition}', New_expression_data.shape)
    print(list(New_expression_data.columns.values))
    Expression_df=New_expression_data[['EP Expression','ESP Expression','SP Expression', 'EP Incarnato Expression']]
    GCSs_df=New_expression_data[['Gyrase exponential Micro DSDS','Gyrase exponential Oxo DSDS','Gyrase exponential RifCfx DSDS','Gyrase exponential Cfx DSDS','Gyrase early stationary Cfx DSDS','Gyrase stationary Cfx DSDS']]
    
    
    cor_df=pd.DataFrame()
    for expression_condition in list(Expression_df.columns.values):
        print(expression_condition)
        expression_condition_data=Expression_df[expression_condition].tolist()
        row_dict={}
        for GCSs_condition in list(GCSs_df.columns.values):
            #print(GCSs_condition)
            GCSs_condition_data=GCSs_df[GCSs_condition].tolist()
            print(f'Paerson correlation ({expression_condition}, {GCSs_condition}) for even bins: ' + str(pearsonr(expression_condition_data, GCSs_condition_data))) 
            row_dict[GCSs_condition]=pearsonr(expression_condition_data, GCSs_condition_data)[0]
        #print(row_dict)
        row_series=pd.Series(row_dict, name=expression_condition)
        #print(row_series)
        cor_df=cor_df.append(row_series)
    #print(cor_df)
    cor_df=cor_df[['Gyrase exponential Micro DSDS','Gyrase exponential Oxo DSDS','Gyrase exponential RifCfx DSDS','Gyrase exponential Cfx DSDS','Gyrase early stationary Cfx DSDS','Gyrase stationary Cfx DSDS']]
    print('\n', cor_df)
    
    correlation_matrix(cor_df, f'Correlation between GCSs data (GSCs number) and transcription\nfor different phases of growth \n({set_name}s del cor no rRNA, rpsl operon, >{additional_condition}, 5000bp DS)', -0.1, 0.4, f'{Path_out}GCSs_number_transcription_corr_5000bp_DS_all_conditions_thr_{additional_condition}.png')
    return cor_df

set_type="Gene"

read_merge_GCSs_in_DS(GCSs_in_DS_regions_files, -1, set_type)
read_merge_GCSs_in_DS(GCSs_in_DS_regions_files, 3, set_type)
read_merge_GCSs_in_DS(GCSs_in_DS_regions_files, 6, set_type)
read_merge_GCSs_in_DS(GCSs_in_DS_regions_files, 9, set_type)
read_merge_GCSs_in_DS(GCSs_in_DS_regions_files, 12, set_type)
read_merge_GCSs_in_DS(GCSs_in_DS_regions_files, 15, set_type)
read_merge_GCSs_in_DS(GCSs_in_DS_regions_files, 18, set_type)
read_merge_GCSs_in_DS(GCSs_in_DS_regions_files, 20, set_type)
read_merge_GCSs_in_DS(GCSs_in_DS_regions_files, 21, set_type)
read_merge_GCSs_in_DS(GCSs_in_DS_regions_files, 25, set_type)
read_merge_GCSs_in_DS(GCSs_in_DS_regions_files, 30, set_type)
"""
read_merge_GCSs_in_DS(GCSs_in_DS_regions_files, -1, set_type)
read_merge_GCSs_in_DS(GCSs_in_DS_regions_files, 20, set_type)
read_merge_GCSs_in_DS(GCSs_in_DS_regions_files, 40, set_type)
read_merge_GCSs_in_DS(GCSs_in_DS_regions_files, 60, set_type)
read_merge_GCSs_in_DS(GCSs_in_DS_regions_files, 80, set_type)
read_merge_GCSs_in_DS(GCSs_in_DS_regions_files, 100, set_type)
read_merge_GCSs_in_DS(GCSs_in_DS_regions_files, 120, set_type)
read_merge_GCSs_in_DS(GCSs_in_DS_regions_files, 140, set_type)
read_merge_GCSs_in_DS(GCSs_in_DS_regions_files, 300, set_type)
read_merge_GCSs_in_DS(GCSs_in_DS_regions_files, 400, set_type)
"""