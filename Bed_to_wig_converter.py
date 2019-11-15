###############################################
##Dmitry Sutormin, 2019##
##ChIP-Seq analysis##

####
#The only purpose - convert bed-like file to wig format and 
#replace Linux line ends (\n) with Windows ones (\r\n)
####

###############################################


#Path to the input file
filein_path_dict={'1' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\WIG\RNA-Seq_cov_depth\Early_Start_1R_S60_R1_001.bed",
                  '2' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\WIG\RNA-Seq_cov_depth\Early_Start_2R_S61_R1_001.bed",
                  '3' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\WIG\RNA-Seq_cov_depth\Early_Start_3R_S62_R1_001.bed",
                  '4' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\WIG\RNA-Seq_cov_depth\Ecoli_Stat_phase_1_prec_S63_R1_001.bed",
                  '5' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\WIG\RNA-Seq_cov_depth\Ecoli_Stat_phase_2_prec_S64_R1_001.bed",
                  '6' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\WIG\RNA-Seq_cov_depth\Ecoli_Stat_phase_3_prec_S65_R1_001.bed",
                  '7' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\WIG\RNA-Seq_cov_depth\Expon_phase_1R_S57_R1_001.bed",
                  '8' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\WIG\RNA-Seq_cov_depth\Expon_phase_2R_S58_R1_001.bed",
                  '9' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\WIG\RNA-Seq_cov_depth\Expon_phase_3R_S59_R1_001.bed",
                  }

#Path to the output file.
fileout_path_dict={'1' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\WIG\RNA-Seq_cov_depth\Early_Start_1R_S60_R1_001.wig",
                   '2' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\WIG\RNA-Seq_cov_depth\Early_Start_2R_S61_R1_001.wig",
                   '3' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\WIG\RNA-Seq_cov_depth\Early_Start_3R_S62_R1_001.wig",
                   '4' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\WIG\RNA-Seq_cov_depth\Ecoli_Stat_phase_1_prec_S63_R1_001.wig",
                   '5' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\WIG\RNA-Seq_cov_depth\Ecoli_Stat_phase_2_prec_S64_R1_001.wig",
                   '6' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\WIG\RNA-Seq_cov_depth\Ecoli_Stat_phase_3_prec_S65_R1_001.wig",
                   '7' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\WIG\RNA-Seq_cov_depth\Expon_phase_1R_S57_R1_001.wig",
                   '8' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\WIG\RNA-Seq_cov_depth\Expon_phase_2R_S58_R1_001.wig",
                   '9' : "F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\WIG\RNA-Seq_cov_depth\Expon_phase_3R_S59_R1_001.wig",}

#ID or short description of the track (will be the name of a track in IGV).
name_dict={'1' : "ES_1",
           '2' : "ES_2",
           '3' : "ES_3",
           '4' : "S_1",
           '5' : "S_2",
           '6' : "S_3",
           '7' : "E_1",
           '8' : "E_2",
           '9' : "E_3",}

#ID of chromosome (for w3110_Mu_SGS: NC_007779.1_w3110_Mu)
Chromosome_name=''
#Mode for Chromosome name writing: 0 - auto detection from bed file provided, 1 - manualy provided by user in Chromosome_name variable.
Auto_or_manual=int(0)


def read_and_convert(filein_path_dict, fileout_path_dict, name_dict, Chromosome_name, Auto_or_manual):
    for sample_name, sample_path in filein_path_dict.items():
        print(f'Now is processing: {sample_path}')
        print(f'Progress: {sample_name}/{len(filein_path_dict)}')
        
        filein=open(filein_path_dict[sample_name], 'r')
        fileout=open(fileout_path_dict[sample_name], 'w')
        
        Ar_of_Cromosome_names=[]
        for line in filein:
            line=line.rstrip().split('\t')
            if line[0] not in Ar_of_Cromosome_names:
                if Auto_or_manual==0:
                    fileout.write('track type=wiggle_0 name="'+name_dict[sample_name]+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+str(line[0])+' start=1 step=1\n')
                elif Auto_or_manual==1:
                    fileout.write('track type=wiggle_0 name="'+name_dict[sample_name]+'" autoScale=off viewLimits=0.0:25.0\nfixedStep chrom='+Chromosome_name+' start=1 step=1\n')
                Ar_of_Cromosome_names.append(line[0])
            else:
                fileout.write(line[2]+'\n')
            
        filein.close()
        fileout.close()    
    return


read_and_convert(filein_path_dict, fileout_path_dict, name_dict, Chromosome_name, Auto_or_manual)