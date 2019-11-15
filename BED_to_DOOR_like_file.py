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

#PWD.
PWD="F:\E_coli_RNA-Seq\E_coli_DY330_RNA-Seq\TUs\FPKM_analysis\\"


#######
#Read broadPeak file and make DOOR-like file.
#######

def BroadPeak_to_DOOR(inpath, outpath):
    BED_data=pd.read_csv(inpath, sep='\t', skiprows=2, header=None)
    print(BED_data)
    door_like=pd.DataFrame()
    door_like['GeneID']=range(1, BED_data.shape[0]+1) #1
    door_like['Gene_name']=BED_data[3].tolist() #2
    door_like['Start']=BED_data[1].tolist() #3
    door_like['End']=BED_data[2].tolist() #4   
    door_like['Strand']=BED_data[5].tolist() #5
    door_like['Expression']=BED_data[4].tolist() #6
    door_like['Gene_description']=BED_data[3].tolist() #7
    door_like['Operon_ID']=['-']*BED_data.shape[0] #8
    door_like.to_csv(outpath, sep='\t', header=True, index=False)
    return

BroadPeak_to_DOOR(PWD+"Stationary_phase_expression_TU.BroadPeak", PWD+"Stationary_phase_expression_TU_door_like.txt")