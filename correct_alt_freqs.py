"""
Created on Mon Nov 25 17:14:36 2019

Recompute ALT variant frequencies based on depths of all variants and indels

@author: david
"""
import pandas as pd

file_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/samples/P0/'
sample = 'P0_r1_refTF2.tsv'
out_file = file_path + 'P0_r1_correctAltFreqs.tsv'

tsv_file = file_path + sample
vdf = pd.read_table(tsv_file, sep="\t")

for index, var in vdf.iterrows():
    
    "Find other variants at this site"
    region = var['REGION']
    site = var['POS']
    site_vars = vdf[(vdf['REGION'] == region) & (vdf['POS'] == site)]
    site_vars = site_vars.reset_index(drop=True) # reset index so site_vars are indexes from 0
    
    alt_depth = var['ALT_DP']
    ref_depth = var['REF_DP']    
    sum_alt_depths = site_vars['ALT_DP'].sum()
    true_total_depth = ref_depth + sum_alt_depths
    
    vdf.at[index,'ALT_FREQ'] = alt_depth / true_total_depth

vdf.to_csv(out_file, sep='\t', index=False) 