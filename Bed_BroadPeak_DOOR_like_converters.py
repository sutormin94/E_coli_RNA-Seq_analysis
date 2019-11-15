####################
#Parsers converting:
#1. tab to bed; 2. tab to BroadPeak; 3. BroadPeak to bed

##Dmitry Sutormin, 2019.
####################

#######
#Memo on interval writing formats:
#1. GFF: +1-based, start position is included, end position is included. []
#2. BroadPeak: +1-based, start position is not included, end position is included. (] <=> 0-based, start position is included, end position is not included. [)
#3. Bed: +1-based, start position is not included, end position is included. (] <=> 0-based, start position is included, end position is not included. [)

#Example (all intervals are vizualized in a same manner in IGV):
#GFF 190\t250
#BroadPeak 189\t250
#Bed 189\t250 blockSize=250-189=61
#Blast results 190\t250 (not shure about end)
#######



def from_door_like_tab_to_BroadPeak():
    #Door_like_tab format is 1+based format, interval end is closed! GFF-like.
    #BroadPeak is 0based format, interval end is opened! (according to IGV vizualisation and specification).
    tab_in=open("C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\\UTRs_w3110_door_like_del_cor.txt", 'r')
    broadPeak_out=open("C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\\UTRs_w3110_door_like_del_cor.BroadPeak", 'w')
    Chromosome_ID="NC_007779.1_w3110_Mu"
    Name="EP_expression_Incarnato_2018"
    Description="Exponential phase expression Incarnato 2018 data"
    
    Broad_peak_header=f'track type=broadPeak visibility=3 db={Chromosome_ID} name="{Name}" description="{Description}"\nbrowser position {Chromosome_ID}:1-100\n'
    broadPeak_out.write(Broad_peak_header)
    
    type_of_TUs=0
    for line in tab_in:
        line=line.rstrip('\n')
        line=line.split('\t')
        if line[0]=='OperonID':
            type_of_TUs=1 #We are itarating operons
        elif line[0]=='GeneID':  
            type_of_TUs=2 #We are itarating genes
        elif line[0]=='TU_ID':  
            type_of_TUs=3 #We are itarating transcripts        
        if line[0] not in ['OperonID', 'GeneID', 'TU_ID']:
            line[5]=line[5].replace(',', '.')
            if type_of_TUs==1:
                broadPeak_out.write(f'{Chromosome_ID}\t{int(line[2])-1}\t{int(line[3])}\t{line[1][1:-1]}\t{line[5]}\t{line[4]}\t{line[5]}\t-1\t-1\n')
            elif type_of_TUs!=1:
                broadPeak_out.write(f'{Chromosome_ID}\t{int(line[2])-1}\t{int(line[3])}\t{line[1]}\t{line[5]}\t{line[4]}\t{line[5]}\t-1\t-1\n')            
    tab_in.close()
    broadPeak_out.close()
    return

#from_door_like_tab_to_BroadPeak()

def from_door_like_tab_to_bed():
    #Door_like_tab format is 1+based format, interval end is closed! GFF-like.
    #Bed is 0based format, interval end is opened! (according to IGV vizualisation and specification).
    tab_in=open("C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\\UTRs_w3110_door_like_del_cor.txt", 'r')
    bed_out=open("C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\\UTRs_w3110_door_like_del_cor.bed", 'w')
    Chromosome_ID="NC_007779.1_w3110_Mu"
    
    type_of_TUs=0
    for line in tab_in:
        line=line.rstrip('\n')
        line=line.split('\t')
        if line[0]=='OperonID':
            type_of_TUs=1 #We are itarating operons
        elif line[0]=='GeneID':  
            type_of_TUs=2 #We are itarating genes
        elif line[0]=='TU_ID':  
            type_of_TUs=3 #We are itarating transcripts.       
        if line[0] not in ['OperonID', 'GeneID', 'TU_ID']:
            line[5]=line[5].replace(',', '.')
            if type_of_TUs==1:
                bed_out.write(f'{Chromosome_ID}\t{int(line[2])-1}\t{int(line[3])}\t{line[1][1:-1]}\t{line[5]}\t{line[4]}\t{int(line[2])-1}\t{int(line[3])}\t255,0,0\t1\t{int(line[3])-int(line[2])+1}\t0,\n')
            if type_of_TUs!=1:
                bed_out.write(f'{Chromosome_ID}\t{int(line[2])-1}\t{int(line[3])}\t{line[1]}\t{line[5]}\t{line[4]}\t{int(line[2])-1}\t{int(line[3])}\t255,0,0\t1\t{int(line[3])-int(line[2])+1}\t0,\n')            
        
    tab_in.close()
    bed_out.close()
    return

#from_door_like_tab_to_bed()


def from_bed_to_door_like_tab():
    #Bed is 0based format, interval end is opened! (according to IGV vizualisation and specification).
    #Door_like_tab format is 1+based format, interval end is closed! GFF-like.
    bed_in=open("C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\\UTRs_w3110.bed", 'r')
    tab_out=open("C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\\UTRs_w3110_door_like.txt", 'w')
    Chromosome_ID="NC_007779.1_w3110_Mu"
    
    type_of_TUs="TU_ID"
    tab_out.write(f"{type_of_TUs}\tGene_name\tStart\tEnd\tStrand\tExpression\tGene_description\n")
    
    i=0
    for line in bed_in:
        line=line.rstrip('\n')
        line=line.split('\t')
        tab_out.write(f'{i}\t{line[3]}\t{int(line[1])+1}\t{int(line[2])}\t{line[5]}\t{line[4]}\t{line[3]}\n')
        i+=1
    
    bed_in.close()
    tab_out.close()    
    return

from_bed_to_door_like_tab()

def from_broadPeak_to_door_like_tab():
    #BroadPeak is 0based format, interval end is opened! (according to IGV vizualisation and specification).
    #Door_like_tab format is 1+based format, interval end is closed! GFF-like.
    broadPeak_in=open("F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\TUs_maps\\UTRs_w3110.broadPeak", 'r')
    tab_out=open("C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\\UTRs_w3110_door_like_del_cor.txt", 'w')
    Chromosome_ID="NC_007779.1_w3110_Mu"
    
    type_of_TUs="TU_ID"
    tab_out.write(f"{type_of_TUs}\tGene_name\tStart\tEnd\tStrand\tExpression\tGene_description\n")
    
    i=0
    for line in broadPeak_in:
        line=line.rstrip('\n')
        line=line.split('\t')
        tab_out.write(f'{i}\t{line[3]}\t{int(line[1])+1}\t{int(line[2])}\t{line[5]}\t{line[4]}\t{line[3]}\n')
        i+=1
    
    broadPeak_in.close()
    tab_out.close()    
    return

#from_broadPeak_to_door_like_tab()