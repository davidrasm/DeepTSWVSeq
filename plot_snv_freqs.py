"""
Created on Mon May 20 18:04:57 2019

Plot frequency of all variants along each genome segment in a single sample

Used to create Figure 3 and Supp Fig. 2 in paper for TF2/TL2 variants

@author: david
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import os
from call_aa_variants import call_aa_variants

def label_point(x, y, val, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    for i, point in a.iterrows():
        ax.text(point['x'], point['y'], str(point['val']))
        
def add_gene(ax, label, start, end):
    
    "Add gene or feature annotation to plot ax from start to end"
    left = start
    bottom = -0.4 * ax.get_ylim()[1]
    width = end - start
    height = 0.1 * ax.get_ylim()[1]
    p = patches.Rectangle(
        (left, bottom), width, height,
        fill=True, color='gray', transform=ax.transData, clip_on=False
        )
    ax.add_patch(p)
    ax.text(left+(width/2), bottom+(height/2), label,
        horizontalalignment='center',
        verticalalignment='center',
        transform=ax.transData)

def add_var_types(variants):
    
    types = []
    for index, var in variants.iterrows(): 
        if var['INDEL']:
            types.append('INDEL')
        elif var['AA_VAR']:
            types.append('AAV')
        else:
            types.append('SNV')
    variants['VAR_TYPE'] = np.array(types)       
    return variants
    
def plot_variants(tsv_file,tsv_rep1,tsv_rep2,fig_name,var_type='SNV',label=False):
    
    filtered = pd.read_table(tsv_file, sep="\t")
    
    "Get average frequency across replicates"
    alt_freq_rep1 = 'ALT_FREQ_' + tsv_rep1
    alt_freq_rep2 = 'ALT_FREQ_' + tsv_rep2
    filtered['AVG_ALT_FREQ'] = (filtered[alt_freq_rep1] + filtered[alt_freq_rep2]) / 2
    filtered = filtered.loc[filtered['AVG_ALT_FREQ'] < 0.97] # remove "fixed" variants (with 3% margin of error)
    
    aa_variant,aa_subst,cds_annotation = call_aa_variants(filtered)
    filtered['AA_VAR'] = np.array(aa_variant)
    filtered['AA_SUBST'] = np.array(aa_subst)
    filtered['CDS_annotation'] = np.array(cds_annotation)
    
    "Call indels"
    filtered['INDEL'] = filtered.ALT.str.contains('^\+|-') # regexp means first character is '+' or '-' char
    
    #"Compute fraction of SNVs"
    #snv_vars = filtered[(filtered['CDS_annotation'] != '') & (filtered['INDEL'] == False)]
    #aa_vars = snv_vars[(snv_vars['AA_VAR'] == True)]
   
    if (var_type == 'SNV'):
       filtered = filtered.loc[filtered['INDEL'] == False]
    elif (var_type == 'AA'):
        filtered = filtered.loc[filtered['AA_VAR'] == True]
    elif (var_type == 'INDEL'):
        filtered = filtered.loc[filtered['INDEL'] == True]
    elif (var_type == 'ANY'):
        print()
    filtered = add_var_types(filtered) # assign variant types
    
    df_segS = filtered.loc[filtered['REGION'] == "TSWV_segS"] # get rows for segment
    df_segM = filtered.loc[filtered['REGION'] == "TSWV_segM"] # get rows for segment
    df_segL = filtered.loc[filtered['REGION'] == "TSWV_segL"] # get rows for segment
    
    sns.set()
    #sns.set_context("talk")
    fig, axs = plt.subplots(3, 1)
    
    "Plot average freq across replicates"
    type_order = ['SNV','AAV','INDEL']
    sns.scatterplot(x="POS", y="AVG_ALT_FREQ", data=df_segS, ax=axs[0], hue="VAR_TYPE", hue_order=type_order, s=40, alpha=.5)
    #df_segS.plot(x="POS", y='AVG_ALT_FREQ', ax=axs[0], kind="scatter", color="blue", alpha = 0.3)
    axs[0].set_xlabel('')
    axs[0].set_ylabel('Frequency')
    axs[0].grid(True)
    axs[0].set_title('Segment S')
    axs[0].set_xlim([0,2915])
    axs[0].set_ylim([0,1.0])
    if label:
        if (var_type == 'AA'):
            label_point(df_segS.POS, df_segS.AVG_ALT_FREQ, df_segS.AA_SUBST, axs[0])
        else:
            label_point(df_segS.POS, df_segS.AVG_ALT_FREQ, df_segS.ALT, axs[0])
    add_gene(axs[0], 'N', 89, 1483) # Add  N 89..1483
    add_gene(axs[0], 'NSs', 1987, 2763) # Add NSs 1987..2763
    
    "Plot average freq across replicates"
    sns.scatterplot(x="POS", y="AVG_ALT_FREQ", data=df_segM, ax=axs[1], hue="VAR_TYPE", hue_order=type_order, legend=False, s=40, alpha=.5)
    #df_segM.plot(x="POS", y='AVG_ALT_FREQ', ax=axs[1], kind="scatter", color="blue", alpha = 0.3)
    axs[1].set_xlabel('')
    axs[1].set_ylabel('Frequency')
    axs[1].grid(True)
    axs[1].set_title('Segment M')
    axs[1].set_xlim([0,4820])
    axs[1].set_ylim([0,1.0])
    if label:
        if (var_type == 'AA'):
            label_point(df_segM.POS, df_segM.AVG_ALT_FREQ, df_segM.AA_SUBST, axs[1])
        else:
            label_point(df_segM.POS, df_segM.AVG_ALT_FREQ, df_segM.ALT, axs[1])
    add_gene(axs[1], 'NSm', 101, 1009) # Add  NSm 101..1009
    
    "If treating GP as one protein"
    #add_gene(axs[1], 'GP', 1330, 4737) # Add GP 1330..4737
    
    "If dividing GP into Gn/Gc"
    add_gene(axs[1], r'$G_C$', 1330, 3285)
    add_gene(axs[1], r'$G_N$', 3286, 4737)
    
    "Plot average freq across replicates"
    sns.scatterplot(x="POS", y="AVG_ALT_FREQ", data=df_segL, ax=axs[2], hue="VAR_TYPE", hue_order=type_order, legend=False, s=40, alpha=.5)
    #df_segL.plot(x="POS", y='AVG_ALT_FREQ', ax=axs[2], kind="scatter", color="blue", alpha = 0.3)
    axs[2].set_xlabel('')
    axs[2].set_ylabel('Frequency')
    axs[2].grid(True)
    axs[2].set_title('Segment L')
    axs[2].set_xlim([0,8896])
    axs[2].set_ylim([0,1.0])
    if label:
        if (var_type == 'AA'):
            label_point(df_segL.POS, df_segL.AVG_ALT_FREQ, df_segL.AA_SUBST, axs[2])
        else:
            label_point(df_segL.POS, df_segL.AVG_ALT_FREQ, df_segL.ALT, axs[2])
    add_gene(axs[2], 'RdRp', 34, 8661) # Add RdRp 34..8661
    
    fig.set_size_inches(8, 8)
    fig.tight_layout()
    plt.show()
    
    fig.savefig(fig_name, dpi=200)
    
if __name__ == '__main__':
    
    batch = False # processing multiple samples?
    
    "For multiple fastq files in a directory"
    if batch:
        
        main_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/samples/'
        results_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/DeepSeqResults/varFreqs/'
        
        "Get all directories in main_path"
        dir_list = os.listdir(main_path)
        
        for dir in dir_list:
            
            print(dir)
            
            if (dir[0] == '.'):
                continue
                
            os.chdir(main_path + dir) # change path to sample_path            
            fig_name = results_path + dir + "_refTF2_allVarFreqs.png" 
            
            "Set tsv file names"
            tsv_file = dir + '_refTF2_pairRepFiltered.tsv'
            tsv_rep1 = dir + '_r1_refTF2.tsv'
            tsv_rep2 = dir + '_r2_refTF2.tsv'
            
            plot_variants(tsv_file,tsv_rep1,tsv_rep2,fig_name,var_type='ANY',label=False)
            os.chdir(main_path)            
            
    else:
        
        main_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/samples/'
        results_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/samples/TF2/'
        dir = 'TF2' #'P1_L1_7'
        os.chdir(main_path + dir) # change path to sample_path
        
        fig_name = results_path + dir + "_allVarFreqs_noLabels_GnGc.png"
        
        tsv_file = dir + '_refTF2_pairRepFiltered.tsv'
        tsv_rep1 = dir + '_r1_refTF2.tsv'
        tsv_rep2 = dir + '_r2_refTF2.tsv'
        
        plot_variants(tsv_file,tsv_rep1,tsv_rep2,fig_name,var_type='ANY',label=False)
