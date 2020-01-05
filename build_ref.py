"""
Created on Thu May 30 13:53:06 2019

Build new reference from variant calls in tsv file
Need an "old" reference genome for indexing
This currently ignores all indels!

Used to create TF2 reference in paper.

@author: david
"""
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
    
def process_variants(variants,ref):
    
    var_bases = []
    var_freqs = []
    for index, var in variants.iterrows():
        alt = var['ALT']
        if ((alt[0] != '+') & (alt[0] != '-')): # if not an indel
            freq_rep1 = var[alt_freq_rep1]
            freq_rep2 = var[alt_freq_rep2]
            var_freqs.append((freq_rep1 + freq_rep2) / 2)
            var_bases.append(alt)
            
    "Add freq and base of reference -- was not doing this before"
    var_freqs.append(1-sum(var_freqs))
    var_bases.append(ref)       
            
    if len(var_freqs) > 0:
        idx = var_freqs.index(max(var_freqs))
        base = var_bases[idx]
    else:
        base = ref
    return base 


"Set up paths to data files"
ref_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/tswv_ref/'
ref_fasta = ref_path + 'tswv_ref_full.fasta'

file_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/samples/P2_L2_14/'
tsv_file = file_path + 'P2_L2_14_refTF2_pairRepFiltered.tsv'

tsv_rep1 = 'P2_L2_14_r1_refTF2.tsv'
tsv_rep2 = 'P2_L2_14_r2_refTF2.tsv'
alt_freq_rep1 = 'ALT_FREQ_' + tsv_rep1
alt_freq_rep2 = 'ALT_FREQ_' + tsv_rep2

out_file = 'P2_L2_14_consensus.fasta'

vdf = pd.read_table(tsv_file, sep="\t")

records = list(SeqIO.parse(ref_fasta, "fasta"))
for record in records:
    ref_id = record.id
    seq = record.seq
    new_seq = ''
    for i in range(len(seq)):
        variants = vdf[(vdf['REGION'] == ref_id) & (vdf['POS'] == i+1)] # had a bug here where POS == i
        base = process_variants(variants,seq[i])
        new_seq = new_seq + base
    record.seq = Seq(new_seq)

SeqIO.write(records, out_file, "fasta")
                
            
            
        
    
    