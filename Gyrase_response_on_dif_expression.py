###############################################
##Dmitry Sutormin, 2019##
##RNA-Seq analysis##

#The script takes file with number of GCSs calculated for US, GB, DS regions of each TU.
#Takes transcription level of each TU, calculates ratio between transcription level of diffrent growth phases,
#and groups TUs by ratio value.
#Than GCSs number is normalized and is compared for different groups of genes and conditions.
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
import math

#######
#Variables to be defined.
#######

print('Variables to be defined:')

#Path to TUs sets with transcription level.
path_to_TUs_sets_no_rRNA_tRNA={'EP' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\Expression_data\DY330_genes\DY330_RNA-Seq_genes_no_tRNA_rRNA_EP_del_cor.txt",
                               'ESP' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\Expression_data\DY330_genes\DY330_RNA-Seq_genes_no_tRNA_rRNA_ESP_del_cor.txt",
                               'SP' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\Expression_data\DY330_genes\DY330_RNA-Seq_genes_no_tRNA_rRNA_SP_del_cor.txt"}

#Path to GCSs number in US, GB, DS of TUs.
path_to_GCSs_number="F:\Topo_data_new_expression\TU_based_analysis\Gyrase\Correction_on_GCSs_calling_thr\Genes_del_cor_no_rRNA_tRNA_5000bp\GCSs_US_GB_DS\Early stationary phase_numbers_of_associated_GCSs_US_GB_DS.txt"

#Path to differential expression data (edgeR table).
path_dif_exp="F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\FPKM_analysis_correct\Genes_del_cor\EdgeR_analysis_dif_expression\E_vs_ES_DE_no_tRNA_rRNA.txt"

#Genome length, bp
Genome_len=4647454
#Genome length corrected on deletions, bp (deletions are specified in the equiation.
Genome_len_dc=Genome_len-((372148-274500)+(807500-793800)+(1214000-1199000))
#Length of the US and DS window used for GCSs picking, bp.
Window_length=5000
#Total number of GCSs observed for different growth phases.
GCSs_num_dict={'EP' : 1673,
               'ESP': 5407,
               'SP': 3508}

#Outpath.
GCSs_and_transcription_outpath="F:\Topo_data_new_expression\TU_based_analysis\Gyrase\Correction_on_GCSs_calling_thr\Genes_del_cor_no_rRNA_tRNA_5000bp\GCSs_US_GB_DS\Figures\Increasing_SP_vs_EP\\"
if not os.path.exists(GCSs_and_transcription_outpath):
    os.makedirs(GCSs_and_transcription_outpath)
    
 
#######
#Plot distribution of transcription levels ratios.
#######   
    
def plot_transcription_ratio(array, name_1, name_2, outpath):
    #Remove nan and inf values.
    #Stolen from here: https://stackoverflow.com/questions/44606799/how-to-clean-nan-and-inf-in-list-type-data-in-python
    array=[v for v in array if not math.isnan(v) and not math.isinf(v)]
    fig=plt.figure(figsize=(7,3), dpi=100)
    plot0=plt.subplot2grid((1,1),(0,0), rowspan=1, colspan=1)   
    print(array.sort())
    print(min(array), max(array))
    bins=[0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000]
    plot0.hist(array, bins, color='#7FCE79', edgecolor='black', alpha=0.8)
    plot0.set_xlabel(f'Transcription {name_1}/{name_2}', size=22)
    plot0.set_xlim(0.0001, 10000) 
    plot0.set_xticklabels(bins, size=20)
    #plot0.set_yticklabels(range(0, 2500, 500), size=20)
    plot0.set_ylabel('Number of TUs', size=22)
    plot0.tick_params(axis='both', which='major', labelsize=20)
    plot0.set_xscale('log')
    #plot0.legend(loc='center right', frameon=False)
    #plot0.set_title(f'{set_name_1} (1) vs\n{set_name_2} (2) N3E', size=16) 
    plt.tight_layout()
    plt.savefig(f'{outpath}Distribution_of_transcription_ratio_{name_1}_on_{name_2}.png', dpi=400, figsize=(7, 3)) 
    plt.close()
    return

#######
#Read TUs transcription data, merge together.
#######

def read_tr_merge(dict_of_pathes, outpath):
    Merged_transcription_data=pd.DataFrame()
    for name, path in dict_of_pathes.items():
        transcr_data=pd.read_csv(path, header=0, sep='\t', index_col=False)
        if Merged_transcription_data.shape==(0,0):
            Merged_transcription_data=transcr_data
        else:
            Merged_transcription_data=pd.merge(Merged_transcription_data, transcr_data, on=['GeneID', 'Gene_name', 'Start', 'End', 'Strand', 'Description'])
        print(Merged_transcription_data.shape)
    
    Merged_transcription_data['ESP/EP']=Merged_transcription_data['Expression_ES']/Merged_transcription_data['Expression_E'] 
    Merged_transcription_data['EP/ESP']=Merged_transcription_data['Expression_E']/Merged_transcription_data['Expression_ES']
    Merged_transcription_data['SP/ESP']=Merged_transcription_data['Expression_S']/Merged_transcription_data['Expression_ES']
    Merged_transcription_data['ESP/SP']=Merged_transcription_data['Expression_ES']/Merged_transcription_data['Expression_S']
    Merged_transcription_data['SP/EP']=Merged_transcription_data['Expression_S']/Merged_transcription_data['Expression_E']
    
    Merged_transcription_data.replace([np.inf, -np.inf], np.nan)
    Merged_transcription_data=Merged_transcription_data.dropna(axis=0)
    
    plot_transcription_ratio(Merged_transcription_data['ESP/EP'].tolist(), 'ESP', 'EP', outpath)
    plot_transcription_ratio(Merged_transcription_data['SP/ESP'].tolist(), 'SP', 'ESP', outpath)
    plot_transcription_ratio(Merged_transcription_data['SP/EP'].tolist(), 'SP', 'EP', outpath)
    #print(Merged_transcription_data.head(10))
    return Merged_transcription_data


#######
#Read table with GCSs numbers data, merge it with transcription dataframe.
#######

def read_gcss_data_merge_transcription(path_to_GCSs, Merged_transcription_data):
    GCSs_data=pd.read_csv(path_to_GCSs, header=0, sep='\t', index_col=False, usecols=['GeneID', 'Gene_name', 'Start', 'End', 'Strand',
                                                                                      'EP GCSs in US', 'EP GCSs in GB', 'EP GCSs in DS',
                                                                                      'ESP GCSs in US', 'ESP GCSs in GB', 'ESP GCSs in DS',
                                                                                      'SP GCSs in US', 'SP GCSs in GB', 'SP GCSs in DS'])
    GCSs_transcription_data=pd.merge(Merged_transcription_data, GCSs_data, on=['GeneID', 'Gene_name', 'Start', 'End', 'Strand'])
    print(GCSs_transcription_data.shape)
    #print(GCSs_transcription_data.head(10))
    return GCSs_transcription_data


#######
#Calculate expected number of GCSs for US/DS and GB. Calculate enrichment in the number of GSCs.
#######

def calc_GCSs_expected_number_calc_enrichment(transcription_GCSs_data, USDS_window, genome_dc_len, GCSs_numbers_dict):
    for condition, GCSs_num in GCSs_numbers_dict.items():
        GCSs_expected=(USDS_window/float(genome_dc_len))*GCSs_num
        print(f'Number of GCSs expected in US/DS for {condition}: {GCSs_expected}')
        transcription_GCSs_data[f'{condition} DS enr']=transcription_GCSs_data[f'{condition} GCSs in DS']/GCSs_expected
        transcription_GCSs_data[f'{condition} US enr']=transcription_GCSs_data[f'{condition} GCSs in US']/GCSs_expected
        transcription_GCSs_data['Gene_length']=transcription_GCSs_data['End']-transcription_GCSs_data['Start']
        transcription_GCSs_data[f'{condition} GCSs in GB expected']=(transcription_GCSs_data['Gene_length']/float(genome_dc_len))*GCSs_num
        transcription_GCSs_data[f'{condition} GB enr']=transcription_GCSs_data[f'{condition} GCSs in GB']/transcription_GCSs_data[f'{condition} GCSs in GB expected']
    
    #Remove rows containing NA values (division by 0 happend).
    transcription_GCSs_data=transcription_GCSs_data.dropna(axis=0)
    #Remove genes (rows) having zero level of transcription.
    transcription_GCSs_data=transcription_GCSs_data[(transcription_GCSs_data['Expression_E']!=0) & (transcription_GCSs_data['Expression_ES']!=0) & (transcription_GCSs_data['Expression_S']!=0)]
    print(transcription_GCSs_data.shape)
    #print(transcription_GCSs_data.head(10))
    return transcription_GCSs_data


#######
#Plot the distribution of GCSs enrichment.
#######

def set_axis_style(ax, labels):
    ax.get_xaxis().set_tick_params(direction='out')
    ax.xaxis.set_ticks_position('bottom')
    ax.set_xticks(np.arange(1, len(labels)+1))
    ax.set_xticklabels(labels, size=25, rotation=90)
    ax.set_xlim(0.25, len(labels)+0.75)
    return

def GCSs_enrichment_and_TUs(transcription_GCSs_data, transcription_GCSs_data_names, max_enrichment, outpath):
    pos=range(1, len(transcription_GCSs_data)+1, 1)
    print(len(pos))
    #Violin plots.
    fig=plt.figure(figsize=(7,10), dpi=100)
    plt1=fig.add_subplot(1,1,1) 
    violins=plt1.violinplot(transcription_GCSs_data, positions=pos, widths=0.9, showmeans=True, showmedians=True, points=200)
    print(violins)
    for vio in violins['bodies']:
        vio.set_facecolor('#ff7762')
        vio.set_edgecolor('black')
        vio.set_alpha(1)
    vmin=violins['cmins']
    vmin.set_linewidth(1)
    vmin.set_color('black')
    vmin.set_alpha(0.7)
    vmean=violins['cmeans']
    vmean.set_linewidth(1)
    vmean.set_color('black')
    vmean.set_alpha(0.7)
    vmax=violins['cmaxes']
    vmax.set_linewidth(1)
    vmax.set_color('black')
    vmax.set_alpha(0.7)
    vbars=violins['cbars']
    vbars.set_linewidth(2)
    vbars.set_color('black')
    vbars.set_alpha(0.7)
    labels=transcription_GCSs_data_names
    set_axis_style(plt1, labels)
    yticknames1=np.arange(0, 70, 2.5)
    plt1.set_yticks(yticknames1, minor=False)
    plt1.set_yticklabels(yticknames1)
    plt1.set_ylabel('GCSs observed/GCSs expected', size=30)
    plt1.set_ylim(-0.4, 7.5)
    plt.setp(plt1.set_yticklabels(yticknames1), rotation=0, fontsize=25)   
    plt1.annotate(f'{len(transcription_GCSs_data[0])}', xy=(0.85, 0.95), xycoords='axes fraction', size=20, rotation=0)
    #plt1.annotate('Mean RNApol \nsignal='+str(round(np.mean(RNApol_level),2)), xy=(0.6, 8), xycoords='data', size=18, rotation=90)
    #plt1.annotate('Mean RNApol \nsignal='+str(round(np.mean(intervals_param_dict['BIMEs1']),2)), xy=(1.6, 8), xycoords='data', size=18, rotation=90)
    #plt1.annotate('Mean RNApol \nsignal='+str(round(np.mean(intervals_param_dict['BIMEs1 with GCSs']),2)), xy=(2.6, 8), xycoords='data', size=18, rotation=90)
    #plt1.annotate('Mean RNApol \nsignal='+str(round(np.mean(intervals_param_dict['BIMEs1 no GCSs']),2)), xy=(3.6, 8), xycoords='data', size=18, rotation=90)
    #plt1.annotate('Mean RNApol \nsignal='+str(round(np.mean(intervals_param_dict['BIMEs2']),2)), xy=(4.6, 8), xycoords='data', size=18, rotation=90)
    #plt1.annotate('Mean RNApol \nsignal='+str(round(np.mean(intervals_param_dict['BIMEs2 with GCSs']),2)), xy=(5.6, 8), xycoords='data', size=18, rotation=90)
    #plt1.annotate('Mean RNApol \nsignal='+str(round(np.mean(intervals_param_dict['BIMEs2 no GCSs']),2)), xy=(6.6, 8), xycoords='data', size=18, rotation=90)
    plt.tight_layout()
    plt.show()
    plt.savefig(outpath, dpi=400, figsize=(7, 10)) 
    plt.close()    
    return


#######
#Prepare data for violinplot: from dataframe to array of vectors.
#######

def convert_df_to_array_of_vectors(input_dataframe):
    dataframe_columns_names=input_dataframe.columns.values.tolist()
    #print(dataframe_columns_names)
    array_of_lists=[]
    list_of_max=[]
    for name in dataframe_columns_names:
        #print(name)
        column_data=input_dataframe[name].tolist()
        array_of_lists.append(column_data)
        list_of_max.append(max(column_data))
    total_max=max(list_of_max)
    return array_of_lists, dataframe_columns_names, total_max


#######
#Select data by changes in transcription.
#######

def select_data_by_transcr(transcription_GCSs_data, cond2, cond1, threshold):
    ratio_type=f'{cond2}/{cond1}'
    #Sort dataframe by ratio between conditions.
    transcription_GCSs_data_sorted=transcription_GCSs_data.sort_values(by=ratio_type, ascending=True)
    #Dataframe with increasing transcription.
    Group1=transcription_GCSs_data_sorted[transcription_GCSs_data_sorted[ratio_type]>threshold]
    #Dataframe with decreasing transcription of equial length.
    Inc_len=Group1.shape[0]
    Group2=transcription_GCSs_data_sorted.head(Inc_len)
    #Check the data.
    print(Group1.shape, Group2.shape)
    print(Group1.head(5), '\n')
    print(Group2.head(5), '\n')
    cond1='EP'
    cond2='SP'
    #Specify columns with GCSs enrichment for conditions under the scope.
    columns_names=[f'{cond1} US enr', f'{cond2} US enr', f'{cond1} GB enr', f'{cond2} GB enr', f'{cond1} DS enr', f'{cond2} DS enr']
    #Enrichment for increasing transcription.
    Group1_enrichment_data=Group1[columns_names]
    #Enrichment for decreasing transcription.
    Group2_enrichment_data=Group2[columns_names]
    #Convert increasing to array.
    Inc_GCSs_enrichment_array, Inc_GCSs_enrichment_names, Inc_Max_enrichment=convert_df_to_array_of_vectors(Group1_enrichment_data)
    #Convert decreasing to array.
    Dec_GCSs_enrichment_array, Dec_GCSs_enrichment_names, Dec_Max_enrichment=convert_df_to_array_of_vectors(Group2_enrichment_data)
    #Merge data.
    Groups_arrays_together=Dec_GCSs_enrichment_array+Inc_GCSs_enrichment_array
    Groups_names_together=Dec_GCSs_enrichment_names+Inc_GCSs_enrichment_names
    Total_max=max([Dec_Max_enrichment, Inc_Max_enrichment])
    return Inc_GCSs_enrichment_array, Inc_GCSs_enrichment_names, Inc_Max_enrichment


#######
#Plot distribution of GCSs enrichments at US, GB, DS.
#######

def plot_GCSs_enr_distr(cond1, cond2, dataframe, sign, thr_pvalue, thr_logFC, outpath):
    fig=plt.figure(figsize=(9, 3))
    #Plot distribution of GCSs enrichment in US.
    plot0=plt.subplot(131)
    bins0=np.histogram(np.hstack((dataframe[f'{cond1} US enr'].tolist(), dataframe[f'{cond2} US enr'].tolist())), bins=15)[1]
    plot0.hist(dataframe[f'{cond1} US enr'].tolist(), bins0, edgecolor='black', linewidth=0.5, density=True, alpha=0.4, label=cond1)
    plot0.hist(dataframe[f'{cond2} US enr'].tolist(), bins0, edgecolor='black', linewidth=0.5, density=True, alpha=0.4, label=cond2)
    Cond1_med=np.median(dataframe[f'{cond1} US enr'].tolist())
    Cond2_med=np.median(dataframe[f'{cond2} US enr'].tolist())
    plot0.annotate(f'{cond1} median={round(Cond1_med,2)}', xy=(0.2, 0.9), xycoords='axes fraction', size=15)
    plot0.annotate(f'{cond2} median={round(Cond2_med,2)}', xy=(0.2, 0.8), xycoords='axes fraction', size=15)
    Cond1_num=len(dataframe[f'{cond1} US enr'].tolist())
    Cond2_num=len(dataframe[f'{cond2} US enr'].tolist())
    plot0.annotate(f'TUs={Cond1_num}', xy=(0.2, 0.7), xycoords='axes fraction', size=15)
    plot0.set_xlabel('Enrichment of GCSs', size=14)
    plot0.set_ylabel('Number of TUs', size=14)
    plot0.legend(loc='center right', frameon=False, fontsize=15)
    plot0.set_title('US', size=20)
    #Plot distribution of GCSs enrichment in GB.
    plot1=plt.subplot(132)
    bins1=np.histogram(np.hstack((dataframe[f'{cond1} GB enr'].tolist(), dataframe[f'{cond2} GB enr'].tolist())), bins=15)[1]
    plot1.hist(dataframe[f'{cond1} GB enr'].tolist(), bins1, edgecolor='black', linewidth=0.5, density=True, alpha=0.4, label=cond1)
    plot1.hist(dataframe[f'{cond2} GB enr'].tolist(), bins1, edgecolor='black', linewidth=0.5, density=True, alpha=0.4, label=cond2)
    Cond1_med=np.median(dataframe[f'{cond1} GB enr'].tolist())
    Cond2_med=np.median(dataframe[f'{cond2} GB enr'].tolist())
    plot1.annotate(f'{cond1} median={round(Cond1_med,2)}', xy=(0.2, 0.9), xycoords='axes fraction', size=15)
    plot1.annotate(f'{cond2} median={round(Cond2_med,2)}', xy=(0.2, 0.8), xycoords='axes fraction', size=15)    
    plot1.set_xlabel('Enrichment of GCSs', size=14)
    plot1.set_ylabel('Number of TUs', size=14)
    plot1.legend(loc='center right', frameon=False, fontsize=15)
    plot1.set_title('GB', size=20)
    #Plot distribution of GCSs enrichment in DS.
    plot2=plt.subplot(133)
    bins2=np.histogram(np.hstack((dataframe[f'{cond1} DS enr'].tolist(), dataframe[f'{cond2} DS enr'].tolist())), bins=15)[1]
    plot2.hist(dataframe[f'{cond1} DS enr'].tolist(), bins2, edgecolor='black', linewidth=0.5, density=True, alpha=0.4, label=cond1)
    plot2.hist(dataframe[f'{cond2} DS enr'].tolist(), bins2, edgecolor='black', linewidth=0.5, density=True, alpha=0.4, label=cond2)
    Cond1_med=np.median(dataframe[f'{cond1} DS enr'].tolist())
    Cond2_med=np.median(dataframe[f'{cond2} DS enr'].tolist())
    plot2.annotate(f'{cond1} median={round(Cond1_med,2)}', xy=(0.2, 0.9), xycoords='axes fraction', size=15)
    plot2.annotate(f'{cond2} median={round(Cond2_med,2)}', xy=(0.2, 0.8), xycoords='axes fraction', size=15)        
    plot2.set_xlabel('Enrichment of GCSs', size=14)
    plot2.set_ylabel('Number of TUs', size=14) 
    plot2.legend(loc='center right', frameon=False, fontsize=15)
    plot2.set_title('DS', size=20)
    #fig.suptitle(f'Transcription level {cond1}{sign}{cond2}\n\n\n\n  1', fontsize=16)
    plt.tight_layout()
    plt.show() 
    
    translate=0
    if sign==">":
        translate="more"
    elif sign=="<":
        translate="less"
    plt.savefig(f'{outpath}Distribution_of_GCSs_enr_US_GB_DS_{cond1}_{translate}_{cond2}_p_value_{thr_pvalue}_logFC_{thr_logFC}.png', dpi=400, figsize=(9, 3))
    return


#######
#Read table with differentially expressed genes (edgeR output).
#Merge it with datframe containing all expression and GCSs counts data.
#######

def read_difexp_data_merge_transcription_gcss(path_to_diffexpr, Merged_transcription_GCSs_data, thr_pvalue, thr_logFC, cond1, cond2, outpath):
    #Read dif exp data and merge.
    Difexp_data=pd.read_csv(path_to_diffexpr, header=0, sep='\t', index_col=False)
    #print(Difexp_data.head(5))
    GCSs_transcription_difexp_data=pd.merge(Merged_transcription_GCSs_data, Difexp_data, on=['Gene_name'])
    print(GCSs_transcription_difexp_data.shape)
    print(GCSs_transcription_difexp_data.head(10))
    #Plot p-value vs logFC.
    plt.figure(figsize=(9, 3))
    plot0=plt.subplot(131)
    plot0.plot(GCSs_transcription_difexp_data['PValue'], abs(GCSs_transcription_difexp_data['logFC']), 'bo', markersize=2)
    plot0.set_xlabel('p-value')
    plot0.set_ylabel('abs(logFC)')
    #Distibution of logFC
    plot1=plt.subplot(132)
    y, x, _ = plot1.hist(GCSs_transcription_difexp_data['logFC'].tolist(), edgecolor='black', linewidth=0.5)
    plot1.plot([-thr_logFC, -thr_logFC], [0, y.max()], '--', linewidth=2)
    plot1.plot([thr_logFC, thr_logFC], [0, y.max()], '--', linewidth=2)
    plot1.set_xlabel('logFC')
    plot1.set_ylabel('Number of TUs')  
    #Distribution of p-value; 
    plot2=plt.subplot(133)
    y, x, _ = plot2.hist(GCSs_transcription_difexp_data['PValue'].tolist(), edgecolor='black', linewidth=0.5)
    plot2.plot([thr_pvalue, thr_pvalue], [0, y.max()], '--', linewidth=2)
    plot2.set_xlabel('p-value')
    plot2.set_ylabel('Number of TUs')   
    plt.tight_layout()
    plt.show()    
    plt.savefig(f'{outpath}Distribution_of_p_value_logFC_{cond1}_{cond2}.png', dpi=400, figsize=(9, 3))
    
    #Select highly dif exp genes.
    Genes_increasing=GCSs_transcription_difexp_data[(GCSs_transcription_difexp_data['PValue']<thr_pvalue) & (GCSs_transcription_difexp_data['logFC']>thr_logFC)]
    Genes_decreasing=GCSs_transcription_difexp_data[(GCSs_transcription_difexp_data['PValue']<thr_pvalue) & (GCSs_transcription_difexp_data['logFC']<-thr_logFC)] 
    print(Genes_decreasing.sort_values(by='logFC', ascending=True).head(10))
    
    plot_GCSs_enr_distr(cond1, cond2, Genes_increasing, '<', thr_pvalue, thr_logFC, outpath)
    plot_GCSs_enr_distr(cond1, cond2, Genes_decreasing, '>', thr_pvalue, thr_logFC, outpath)
    
    return GCSs_transcription_difexp_data



#######
#Wrapper function.
#######

def wrap_functions(TUs_pathin, GCSs_pathin, window_length, genome_len_dc, GCSs_numbers_dict, threshold, path_to_diffexpr, thr_pvalue, thr_logFC, outpath):
    #Read and merge transcription data.
    Transcription_data=read_tr_merge(TUs_pathin, outpath)
    #Read GCSs data and merge with transcription data.
    Transcription_GCSs_data=read_gcss_data_merge_transcription(GCSs_pathin, Transcription_data)
    #Calculate enrichment of GCSs numbers.
    Transcription_GCSs_data_with_enrichment=calc_GCSs_expected_number_calc_enrichment(Transcription_GCSs_data, window_length, genome_len_dc, GCSs_numbers_dict)
    #Calculate GCSs enrichments distribution for all conditions and regions.
    GCSs_enrichment_data=Transcription_GCSs_data_with_enrichment[['EP US enr', 'EP GB enr', 'EP DS enr', 'ESP US enr', 'ESP GB enr', 'ESP DS enr', 'SP US enr', 'SP GB enr', 'SP DS enr']]
    GCSs_enrichment_array, GCSs_enrichment_names, Max_enrichment=convert_df_to_array_of_vectors(GCSs_enrichment_data)
    GCSs_enrichment_and_TUs(GCSs_enrichment_array, GCSs_enrichment_names, Max_enrichment, outpath + "Enrichment_of_GCSs_all_genes.png")    
    
    #Select subset of genes by expression level or other conditions.
    GCSs_enrichment_array, GCSs_enrichment_names, Max_enrichment=select_data_by_transcr(Transcription_GCSs_data_with_enrichment, 'SP', 'EP', threshold)
    GCSs_enrichment_and_TUs(GCSs_enrichment_array, GCSs_enrichment_names, Max_enrichment, f'{outpath}Enrichment_of_GCSs_SP_vs_EP_thr_{threshold}_eq_len_inc.png') 
    
    #Read dif exp data, select highly expressed genes, plot distributions of GCSs enrichment.
    read_difexp_data_merge_transcription_gcss(path_to_diffexpr, Transcription_GCSs_data_with_enrichment, thr_pvalue, thr_logFC, 'EP', 'SP', outpath)
    return

wrap_functions(path_to_TUs_sets_no_rRNA_tRNA, path_to_GCSs_number, Window_length, Genome_len_dc, GCSs_num_dict, 5, path_dif_exp, 0.001, 2, GCSs_and_transcription_outpath)


#wrap_functions(path_to_TUs_sets_no_rRNA_tRNA, path_to_GCSs_number, Window_length, Genome_len_dc, GCSs_num_dict, 5, GCSs_and_transcription_outpath)
#wrap_functions(path_to_TUs_sets_no_rRNA_tRNA, path_to_GCSs_number, Window_length, Genome_len_dc, GCSs_num_dict, 10, GCSs_and_transcription_outpath)
#wrap_functions(path_to_TUs_sets_no_rRNA_tRNA, path_to_GCSs_number, Window_length, Genome_len_dc, GCSs_num_dict, 20, GCSs_and_transcription_outpath)
#wrap_functions(path_to_TUs_sets_no_rRNA_tRNA, path_to_GCSs_number, Window_length, Genome_len_dc, GCSs_num_dict, 50, GCSs_and_transcription_outpath)
#wrap_functions(path_to_TUs_sets_no_rRNA_tRNA, path_to_GCSs_number, Window_length, Genome_len_dc, GCSs_num_dict, 100, GCSs_and_transcription_outpath)
#wrap_functions(path_to_TUs_sets_no_rRNA_tRNA, path_to_GCSs_number, Window_length, Genome_len_dc, GCSs_num_dict, 200, GCSs_and_transcription_outpath)
#wrap_functions(path_to_TUs_sets_no_rRNA_tRNA, path_to_GCSs_number, Window_length, Genome_len_dc, GCSs_num_dict, 400, GCSs_and_transcription_outpath)
#wrap_functions(path_to_TUs_sets_no_rRNA_tRNA, path_to_GCSs_number, Window_length, Genome_len_dc, GCSs_num_dict, 500, GCSs_and_transcription_outpath)