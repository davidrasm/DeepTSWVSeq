
"""
Created on Tue Oct  8 09:02:22 2019

Computes avg. pairwise distances/divergence among sequences in each sample
Then plot time series of avg. pairwise distance across samples

Used to create diversity/divergence plots in paper (Figure 4)

@author: david
"""
import os
from Bio import SeqIO
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

"Old method with single sample"
def average_dist_defunct(vdf):
    
    site_dists = []
    for rec in records:
        seq = rec.seq
        ref_id = rec.id
        for k in range(len(seq)):
            variants = vdf[(vdf['REGION'] == ref_id) & (vdf['POS'] == k+1)]
            freqs = []
            for index, var in variants.iterrows():
                alt = var['ALT']
                if not ("+" in alt or "-" in alt): # ingore indels
                    freqs.append(var['AVG_ALT_FREQ'])    
            
            "Compute avg. distance/diversity based on freqs"
            if freqs:
                ref_freq = 1 - sum(freqs) # assume freq of ref is 1 - sum(freqs)
                freqs.append(ref_freq)
                freqs = np.array(freqs) # convert to numpy array
                #freqs = freqs / sum(freqs) # normalize
                avg_dist = np.dot(freqs,1-freqs) # compute avg dist as dot product
            else:
                avg_dist = 0.0
            site_dists.append(avg_dist)
    
    sample_avg_dist = sum(site_dists) # / len(site_dists)
    return sample_avg_dist

"New method allowing dist to be computed within or between samples"
def average_dist(vdf,refdf):
    
    "Compute average distance/divergence between reference and sample"
    
    site_dists = []
    for rec in records:
        seq = rec.seq
        ref_id = rec.id
        for k in range(len(seq)):
            
            ref_base = seq[k]
            
            "Get var freqs in sample"
            variants = vdf[(vdf['REGION'] == ref_id) & (vdf['POS'] == k+1)]
            freqs = []
            bases = []
            for index, var in variants.iterrows():
                alt = var['ALT']
                if not ("+" in alt or "-" in alt): # ingore indels
                    freqs.append(var['AVG_ALT_FREQ'])
                    bases.append(alt)
            freqs.append(1-sum(freqs))
            bases.append(ref_base)
            
            if vdf is refdf:
                
                freqs = np.array(freqs)
                not_freqs = 1 - freqs
                
            else:
                    
                "Get var freqs in reference"
                variants = refdf[(refdf['REGION'] == ref_id) & (refdf['POS'] == k+1)]
                ref_freqs = []
                ref_bases = []
                for index, var in variants.iterrows():
                    alt = var['ALT']
                    if not ("+" in alt or "-" in alt): # ingore indels
                        ref_freqs.append(var['AVG_ALT_FREQ'])
                        ref_bases.append(alt)
                ref_freqs.append(1-sum(ref_freqs))
                ref_bases.append(ref_base)
                
                not_freqs = []
                for i,fx in enumerate(freqs):
                    if bases[i] in ref_bases:
                        base_index = ref_bases.index(bases[i])
                        not_freqs.append(1 - ref_freqs[base_index])
                    else:
                        not_freqs.append(1.0)
                        
                freqs = np.array(freqs)
                not_freqs = np.array(not_freqs)
            
            "Compute avg. distance/diversity based on freqs"
            avg_dist = np.dot(freqs,not_freqs) # compute avg dist as dot product
            
            site_dists.append(avg_dist)
    
    sample_avg_dist = sum(site_dists) # / len(site_dists)
    return sample_avg_dist

def get_variant_df(sample):
    
    "Get variants in sample"
    tsv_file = sample + '_refTF2_pairRepFiltered.tsv'
    tsv_rep1 = sample + '_r1_refTF2.tsv'
    tsv_rep2 = sample + '_r2_refTF2.tsv'
    vdf = pd.read_table(tsv_file, sep="\t")
    
    "Get average frequency across replicates"
    alt_freq_rep1 = 'ALT_FREQ_' + tsv_rep1
    alt_freq_rep2 = 'ALT_FREQ_' + tsv_rep2
    vdf['AVG_ALT_FREQ'] = (vdf[alt_freq_rep1] + vdf[alt_freq_rep2]) / 2
    
    return vdf

"Compute divergences between samples and a reference"
def compute_divs(lines):
    
    "Get ref variants"
    ref_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/samples/P0/'
    os.chdir(ref_path)
    ref = 'P0'
    ref_df = get_variant_df(ref)
    
    for l,samples in enumerate(lines):
        avg_divs = []
        data_file = 'tswv_line' + str(l+1) + '_avgPairDivergence_toP0.csv'
        for s in samples:
            print(s)
            os.chdir(main_path + s) # change path to sample_path
            sample_df = get_variant_df(s)
            avg_divs.append(average_dist(sample_df,ref_df))
    
        os.chdir(results_path)
        short_names = [s.replace('_L' + str(l+1),'') for s in samples]
        data = {'Samples':short_names, 'Distances':avg_divs}    
        df = pd.DataFrame(data)
        df.to_csv(data_file, sep=',')

"Compute pairwise distances within samples"        
def compute_dists(lines):

    for l,samples in enumerate(lines):
        avg_dists = []
        data_file = 'tswv_line' + str(l+1) + '_avgPairDists.csv'
        for s in samples:
            print(s)
            os.chdir(main_path + s) # change path to sample_path
            sample_df = get_variant_df(s)
            avg_dists.append(average_dist(sample_df,sample_df))
    
        os.chdir(results_path)
        short_names = [s.replace('_L' + str(l+1),'') for s in samples]
        data = {'Samples':short_names, 'Distances':avg_dists}    
        df = pd.DataFrame(data)
        df.to_csv(data_file, sep=',')  

def test_divs():
    
    "Make sure avg_dist and avg_div return the sample results for TF2 sample"
    
    "Get ref variants"
    ref_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/samples/TF2/'
    os.chdir(ref_path)
    ref = 'TF2'
    ref_df = get_variant_df(ref)
    
    "Should return same results but did not finish implementing test"
    avg_dist = average_dist(ref_df,ref_df)
    
    print("Average distance = " + str(avg_dist))
    print("")

if __name__ == '__main__':
    
    "Set reference genome for indexing purposes"  
    ref_fasta = '../tswv_ref/TF2_consensus.fasta'
    records = list(SeqIO.parse(ref_fasta, "fasta"))
    
    #test_divs() #run test
    
    "Need to iter through lines too"
    samples_line1 = ['P0','WF_P1','P1_L1_7','P1_L1_14','P1_L1_18','P1_L1_28','WF_P2_L1','P2_L1_7','P2_L1_14','P2_L1_21','P2_L1_29','WF_P3_L1','P3_L1_7','P3_L1_14','P3_L1_21','P3_L1_28','WF_P4_L1','P4_L1_7','P4_L1_14','P4_L1_21','P4_L1_28','WF_P5_L1','P5_L1_14','P5_L1_21','P5_L1_28'] # all line 1
    samples_line2 = ['P0','WF_P1','P1_L2_7','P1_L2_14','P1_L2_18','P1_L2_28','WF_P2_L2','P2_L2_7','P2_L2_14','P2_L2_21','P2_L2_29','WF_P3_L2','P3_L2_7','P3_L2_14','P3_L2_21','P3_L2_28','WF_P4_L2','P4_L2_14','P4_L2_21','P4_L2_28','WF_P5_L2','P5_L2_7','P5_L2_14','P5_L2_21','P5_L2_28']
    samples_line3 = ['P0','WF_P1','P1_L3_7','P1_L3_14','P1_L3_18','P1_L3_28','WF_P2_L3','P2_L3_7','P2_L3_14','P2_L3_21','P2_L3_29','WF_P3_L3','P3_L3_7','P3_L3_14','P3_L3_21','P3_L3_28','WF_P4_L3','P4_L3_7','P4_L3_14','P4_L3_21','P4_L3_28','WF_P5_L3','P5_L3_7','P5_L3_14','P5_L3_21','P5_L3_28'] # all line 3
    lines = [samples_line1,samples_line2,samples_line3]
    
    main_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/samples/'
    results_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/DeepSeqResults/diversity/'
    
    "Only need to run compute_dists if we don't already have these saved"
    #compute_divs(lines)
    #compute_dists(lines)
    
    sns.set()
    fig, axs = plt.subplots(3, 1)
    fig.set_size_inches(6,9)

    
    df = pd.read_csv(results_path + 'tswv_line1_avgPairDists.csv', sep=',')
    axs[0].plot(df['Samples'],df['Distances'])
    df = pd.read_csv(results_path + 'tswv_line1_avgPairDivergence_toP0.csv', sep=',')
    axs[0].plot(df['Samples'],df['Distances'],color='orange')
    labels = [r'$P_0$', r'$T_1$',r'$P_1^7$',r'$P_1^{14}$',r'$P_1^{18}$',r'$P_1^{28}$',r'$T_2$',r'$P_2^7$',r'$P_2^{14}$',r'$P_2^{21}$',r'$P_2^{29}$',r'$T_3$',r'$P_3^7$',r'$P_3^{14}$',r'$P_3^{21}$',r'$P_3^{28}$',r'$T_4$',r'$P_4^7$',r'$P_4^{14}$',r'$P_4^{21}$',r'$P_4^{28}$',r'$T_5$',r'$P_5^{14}$',r'$P_5^{21}$',r'$P_5^{28}$']
    axs[0].set_xticklabels(labels,rotation='vertical')    
    my_colors = ['g', 'k', 'b','b','b','b','k', 'g','g','g','g','k','b','b','b','b','k','g','g','g','g','k','b','b','b'] # line 1
    "Color tick labels to emphasize host switches"
    for ticklabel, tickcolor in zip(axs[0].get_xticklabels(), my_colors):
        ticklabel.set_color(tickcolor)
    axs[0].set_title('Line 1: Alternating')
    axs[0].legend(['Diversity','Divergence'])
    #axs[0].set_ylim([15,55])
    
    df = pd.read_csv(results_path + 'tswv_line2_avgPairDists.csv', sep=',')
    axs[1].plot(df['Samples'],df['Distances'])
    df = pd.read_csv(results_path + 'tswv_line2_avgPairDivergence_toP0.csv', sep=',')
    axs[1].plot(df['Samples'],df['Distances'],color='orange')
    #sns.relplot(x="Samples", y="Distances", kind="line", data=df, ax=axs[1])
    labels = [r'$P_0$', r'$T_1$',r'$P_1^7$',r'$P_1^{14}$',r'$P_1^{18}$',r'$P_1^{28}$',r'$T_2$',r'$P_2^7$',r'$P_2^{14}$',r'$P_2^{21}$',r'$P_2^{29}$',r'$T_3$',r'$P_3^7$',r'$P_3^{14}$',r'$P_3^{21}$',r'$P_3^{28}$',r'$T_4$',r'$P_4^{14}$',r'$P_4^{21}$',r'$P_4^{28}$',r'$T_5$',r'$P_5^{7}$',r'$P_5^{14}$',r'$P_5^{21}$',r'$P_5^{28}$']
    axs[1].set_xticklabels(labels,rotation='vertical')   
    my_colors = ['g', 'k', 'g','g','g','g','k', 'g','g','g','g','k','g','g','g','g','k','g','g','g','k','b','b','b','b'] # line 2
    "Color tick labels to emphasize host switches"
    for ticklabel, tickcolor in zip(axs[1].get_xticklabels(), my_colors):
        ticklabel.set_color(tickcolor)
    axs[1].set_ylabel('Average Pairwise Distance/Divergence')
    axs[1].set_title('Line 2: Emilia')
    #axs[1].set_ylim([15,55])
    
    df = pd.read_csv(results_path + 'tswv_line3_avgPairDists.csv', sep=',')
    axs[2].plot(df['Samples'],df['Distances'])
    df = pd.read_csv(results_path + 'tswv_line3_avgPairDivergence_toP0.csv', sep=',')
    axs[2].plot(df['Samples'],df['Distances'],color='orange')
    #sns.relplot(x="Samples", y="Distances", kind="line", data=df, ax=axs[2])
    labels = [r'$P_0$', r'$T_1$',r'$P_1^7$',r'$P_1^{14}$',r'$P_1^{18}$',r'$P_1^{28}$',r'$T_2$',r'$P_2^7$',r'$P_2^{14}$',r'$P_2^{21}$',r'$P_2^{29}$',r'$T_3$',r'$P_3^7$',r'$P_3^{14}$',r'$P_3^{21}$',r'$P_3^{28}$',r'$T_4$',r'$P_4^7$',r'$P_4^{14}$',r'$P_4^{21}$',r'$P_4^{28}$',r'$T_5$',r'$P_5^{7}$',r'$P_5^{14}$',r'$P_5^{21}$',r'$P_5^{28}$']
    axs[2].set_xticklabels(labels,rotation='vertical')  
    my_colors = ['b', 'k', 'b','b','b','b','k', 'b','b','b','b','k','b','b','b','b','k','b','b','b','b','k','g','g','g','g'] # line 2
    "Color tick labels to emphasize host switches"
    for ticklabel, tickcolor in zip(axs[2].get_xticklabels(), my_colors):
        ticklabel.set_color(tickcolor)
    axs[2].set_xlabel('Sample')
    axs[2].set_title('Line 3: Datura')
    #axs[2].set_ylim([15,55])
    
    fig.tight_layout()
    plt.show()
    fig.savefig('tswv_allLines_avgPairDist_divP0.png', dpi=200)
    
    #"To concat dataframes for plotting both div and dist the seaborn way"
    #df_dist = pd.read_csv('tswv_line1_avgPairDists.csv', sep=',')
    #df_dist['MEASURE'] = 'Dist'
    #df_div = pd.read_csv('tswv_line1_avgPairDivergence_toTF2.csv', sep=',')
    #df_div['MEASURE'] = 'Div'
    #df = pd.concat([df_dist,df_div],ignore_index=True)
    #sns.relplot(x="Samples", y="Distances", kind="line", data=df, ax=axs[0])
        