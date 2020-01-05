"""
Created on Mon May 20 18:04:57 2019

Plot frequency of variants between paired technical (RT) reps
Can also compute histogram of R^2 values between reps across all samples

Used to create Figure 2 in paper.

@author: david
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import matplotlib.patches as patches
import os
import seaborn as sns
from sklearn.linear_model import LinearRegression
#from Bio import SeqIO
#from Bio.Seq import Seq

def plot_variants(sample,fig_name,var_type='SNV',plot=True):
    
    tsv_file = sample + '_refTF2_pairRepFiltered.tsv'
    tsv_rep1 = sample + '_r1_refTF2.tsv'
    tsv_rep2 = sample + '_r2_refTF2.tsv'
    
    #tsv_file = sample + '_corrAltFreqs_pairRepFiltered.tsv'
    #tsv_rep1 = sample + '_r1_corrAltFreqs.tsv'
    #tsv_rep2 = sample + '_r2_corrAltFreqs.tsv'
    
    filtered = pd.read_table(tsv_file, sep="\t")
    
    "Get average frequency across replicates"
    alt_freq_rep1 = 'ALT_FREQ_' + tsv_rep1
    alt_freq_rep2 = 'ALT_FREQ_' + tsv_rep2
    
    filtered['INDEL'] = filtered.ALT.str.contains('^\+|-') # regexp means first character is '+' or '-' char    
    if (var_type == 'SNV'):
       filtered = filtered.loc[filtered['INDEL'] == False]
    elif (var_type == 'AA'):
        print("AA calling not implemented here yet!!")
    elif (var_type == 'INDEL'):
        filtered = filtered.loc[filtered['INDEL'] == True]
    
    #df_segS = filtered.loc[filtered['REGION'] == "TSWV_segS"] # get rows for segment
    #df_segM = filtered.loc[filtered['REGION'] == "TSWV_segM"] # get rows for segment
    #df_segL = filtered.loc[filtered['REGION'] == "TSWV_segL"] # get rows for segment
    
    "Compute R2"
    regressor = LinearRegression()
    x=filtered[alt_freq_rep1].values.reshape(-1,1)
    y=filtered[alt_freq_rep2].values.reshape(-1,1)
    regressor.fit(x,y) #training the algorithm
    r2_score = regressor.score(x,y)    
    print("R2:",r2_score)
    
    if plot:
    
        #sns.set(style="whitegrid")
        sns.set(color_codes=True)
        fig, ax = plt.subplots(figsize=(4,4))
    
        "Plot correlation between replicates"
        #filtered.plot(x=alt_freq_rep1, y=alt_freq_rep2, ax=ax, kind="scatter", color="black", alpha = 0.3)
        sns.regplot(x=alt_freq_rep1, y=alt_freq_rep2, data=filtered)
        ax.set_xlabel('Frequency Rep 1')
        ax.set_ylabel('Frequency Rep 2')
    
        r2_str = '%.2f' % r2_score
        ax.text(0.05, ax.get_ylim()[1]*0.9, r"$R^2 = $" + r2_str, fontsize=12)
        
        #fig.set_size_inches(10, 8)
        fig.tight_layout()
        plt.show()
        fig.savefig(fig_name, dpi=200)
        
    return r2_score

def plot_sample_hist(vals,fig_name):
    
    #sns.set(style="whitegrid")
    fig, ax = plt.subplots(figsize=(4, 4))
    
    sns.distplot(vals, kde=False, rug=False, ax=ax)
    ax.set_xlabel(r"$R^2$", fontsize=12)
    ax.set_ylabel('Samples', fontsize=12)

    plt.show()
    fig.savefig(fig_name, dpi=200,bbox_inches='tight')
    
if __name__ == '__main__':
    
    batch = False # processing multiple samples?
    
    "For multiple fastq files in a directory"
    if batch:
        
        main_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/samples/'
        results_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/DeepSeqScripts/'
        
        "Get all directories in main_path"
        dir_list = os.listdir(main_path)
        sample_r2 = []
        
        for dir in dir_list:
            
            print(dir)
            
            if (dir[0] == '.'):
                continue
                
            os.chdir(main_path + dir) # change path to sample_path            
            fig_name = results_path + dir + "_pairedRTRep_varFreqs.png"
            if dir == 'P0':           
                r2 = plot_variants(dir,fig_name,var_type='SNV',plot=True)
            else:
                r2 = plot_variants(dir,fig_name,var_type='SNV',plot=False)
            sample_r2.append(r2)
            os.chdir(main_path)

        "Plot freq hist of R2 values across samples"
        fig_name = results_path + 'tswv' + "_pairedRTRep_r2hist.png"
        plot_sample_hist(sample_r2,fig_name)    
            
    else:
        
        main_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/samples/'
        results_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/DeepSeqScripts/'
        dir = 'P0'
        
        fig_name = results_path + dir + "_pairedRTRep_varFreqs.png"
        os.chdir(main_path + dir) # change path to sample_path
        plot_variants(dir,fig_name,var_type='SNV')

#samples = ['WF-P1','P1-L1-7','WF-P2-L1','P2-L1-14','P2-L1-21','WF-P3-L1']
#for s in samples:
    #plot_variants(s)