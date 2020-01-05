"""
Created on Mon Jun 10 13:12:06 2019

Build a master list of variants from individual samples and record their freqs across all samples

Updates 10.15.19:
Updated to use new variants called with TF2 reference
Sites are now indexed from one to be consistant with how sites are numbered in variant files
We were not doing this before but this should not have created a large problem: we would have just been ignoring the last site

@author: david
"""
from Bio import SeqIO
import pandas as pd

out_file = 'tswv_L1-3_masterVarList_corrAltFreqs_withDepths.csv' # final master list used for results in manuscript
main_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/samples/'

depths = True # include depth information in addition to allele freqs

"Create list of sample datasets"
#samples = ['P0','WF_P1','P1_L1_7','P1_L1_14','P1_L1_18','P1_L1_28','WF_P2_L1','P2_L1_7','P2_L1_14','P2_L1_21','P2_L1_29','WF_P3_L1','P3_L1_7','P3_L1_14','P3_L1_21','P3_L1_28','WF_P4_L1','P4_L1_7','P4_L1_14','P4_L1_21','P4_L1_28','WF_P5_L1','P5_L1_14','P5_L1_21','P5_L1_28'] # all line 1
#samples = ['P0','WF_P1','P1_L2_7','P1_L2_14','P1_L2_18','P1_L2_28','WF_P2_L2','P2_L2_7','P2_L2_14','P2_L2_21','P2_L2_29','WF_P3_L2','P3_L2_7','P3_L2_14','P3_L2_21','P3_L2_28','WF_P4_L2','P4_L2_7','P4_L2_14','P4_L2_21','P4_L2_28','WF_P5_L2','P5_L2_7','P5_L2_14','P5_L2_21','P5_L2_28']
#samples = ['P0','WF_P1','P1_L3_7','P1_L3_14','P1_L3_18','P1_L3_28','WF_P2_L3','P2_L3_7','P2_L3_14','P2_L3_21','P2_L3_29','WF_P3_L3','P3_L3_7','P3_L3_14','P3_L3_21','P3_L3_28','WF_P4_L3','P4_L3_7','P4_L3_14','P4_L3_21','P4_L3_28','WF_P5_L3','P5_L3_7','P5_L3_14','P5_L3_21','P5_L3_28'] # all line 3

"Concat all samples together - removed P4_L2_7 b/c of low coverage"
samples_L1 = ['P0','WF_P1','P1_L1_7','P1_L1_14','P1_L1_18','P1_L1_28','WF_P2_L1','P2_L1_7','P2_L1_14','P2_L1_21','P2_L1_29','WF_P3_L1','P3_L1_7','P3_L1_14','P3_L1_21','P3_L1_28','WF_P4_L1','P4_L1_7','P4_L1_14','P4_L1_21','P4_L1_28','WF_P5_L1','P5_L1_14','P5_L1_21','P5_L1_28'] # all line 1
samples_L2 = ['P1_L2_7','P1_L2_14','P1_L2_18','P1_L2_28','WF_P2_L2','P2_L2_7','P2_L2_14','P2_L2_21','P2_L2_29','WF_P3_L2','P3_L2_7','P3_L2_14','P3_L2_21','P3_L2_28','WF_P4_L2','P4_L2_14','P4_L2_21','P4_L2_28','WF_P5_L2','P5_L2_7','P5_L2_14','P5_L2_21','P5_L2_28']
samples_L3 = ['P1_L3_7','P1_L3_14','P1_L3_18','P1_L3_28','WF_P2_L3','P2_L3_7','P2_L3_14','P2_L3_21','P2_L3_29','WF_P3_L3','P3_L3_7','P3_L3_14','P3_L3_21','P3_L3_28','WF_P4_L3','P4_L3_7','P4_L3_14','P4_L3_21','P4_L3_28','WF_P5_L3','P5_L3_7','P5_L3_14','P5_L3_21','P5_L3_28'] # all line 3
samples = samples_L1 + samples_L2 + samples_L3

tsv_rep1_str = '_r1_corrAltFreqs.tsv'
tsv_rep2_str = '_r2_corrAltFreqs.tsv'

sample_dfs = {} # sample dataframes
for s in samples:
    file_path = main_path + s + '/'
    tsv_file = file_path + s + '_corrAltFreqs_pairRepFiltered.tsv'
    sample_dfs[s] = pd.read_table(tsv_file, sep="\t")

"Set reference genome -- just for indexing sites"
ref_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/tswv_ref/'  
ref_fasta = ref_path + 'tswv_ref_full.fasta'
records = list(SeqIO.parse(ref_fasta, "fasta"))

"Create master df"
master_columns = ['REGION','POS','REF','ALT']
for s in samples:
    master_columns.append('ALT_FREQ_' + s)
    if depths:
        master_columns.append('REF_DP_' + s)
        master_columns.append('ALT_DP_' + s) 
        
master = pd.DataFrame(columns=master_columns)

for rec in records:
    
    print(rec.name)
    
    seq = rec.seq
    ref_id = rec.id
    for i in range(len(seq)):
        
        print(i)
        
        for s in samples:
            
            print(s)
            
            "Get str names for reps"
            alt_freq_rep1_str = 'ALT_FREQ_' + s + tsv_rep1_str
            alt_freq_rep2_str = 'ALT_FREQ_' + s + tsv_rep2_str
            
            "Get str names for depth variables"
            alt_dp_rep1_str = 'ALT_DP_' + s + tsv_rep1_str
            alt_dp_rep2_str = 'ALT_DP_' + s + tsv_rep2_str
            ref_dp_rep1_str = 'REF_DP_' + s + tsv_rep1_str
            ref_dp_rep2_str = 'REF_DP_' + s + tsv_rep2_str
            
            vdf = sample_dfs[s]
            variants = vdf[(vdf['REGION'] == ref_id) & (vdf['POS'] == i+1)]
            for index, var in variants.iterrows():
                
                ref = var['REF']
                alt = var['ALT']
                alt_rep1 = var[alt_freq_rep1_str]
                alt_rep2 = var[alt_freq_rep2_str]
                avg = (alt_rep1 + alt_rep2) / 2 # take average of reps
                
                if depths:
                    alt_dp_rep1 = var[alt_dp_rep1_str]
                    alt_dp_rep2 = var[alt_dp_rep2_str]
                    alt_dp_total = alt_dp_rep1 + alt_dp_rep2 # get total
                    ref_dp_rep1 = var[ref_dp_rep1_str]
                    ref_dp_rep2 = var[ref_dp_rep2_str]
                    ref_dp_total = ref_dp_rep1 + ref_dp_rep2 # get total

                idx = master.index[(master['REGION'] == ref_id) & (master['POS'] == i+1) & (master['ALT'] == alt)].tolist()
                print(idx)
                if idx: # add data
                    master.at[idx,'ALT_FREQ_' + s] = avg
                    if depths:
                        master.at[idx,'ALT_DP_' + s] = alt_dp_total
                        master.at[idx,'REF_DP_' + s] = ref_dp_total 
                else: # add row to master
                    master = master.append({'REGION' : ref_id, 'POS' : i+1, 'REF' : ref, 'ALT': alt} , ignore_index=True)
                    master.at[master.index.max(),'ALT_FREQ_' + s] = avg
                    if depths:
                        master.at[master.index.max(),'ALT_DP_' + s] = alt_dp_total # was formely idx instead of index.max()
                        master.at[master.index.max(),'REF_DP_' + s] = ref_dp_total 
                    
                    
master = master.fillna(0)
master.to_csv(out_file, sep='\t', index=False)               
                           