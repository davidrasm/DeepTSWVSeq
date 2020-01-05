"""
Created on Tue May 28 16:24:56 2019

Call variants in samtools and then filter paired seq replicates using iVar

This is script used to call all variants for main results in paper

@author: david
"""
import os
import sys
import subprocess
import pandas as pd

def call_variants(ref_genome,bam_rep1,bam_rep2,out_rep1,out_rep2,out_filtered):

    """
    ivar options:
    -q sets min quality score in iVar (default is also 20)
    -t sets minimum frequency threshold to call variant (default is 0.03)
    -p is prefix for output .tsv file
    """
    
    q_val = str(20)
    t_val = str(0.03)
    
    "Call variants for rep 1"
    cmd = 'samtools mpileup --reference ' + ref_genome + ' -A -d 1000000 -F 0 -B -Q 0 ' + bam_rep1 + ' | ivar variants -p ' + out_rep1 + ' -q ' + q_val + ' -t ' + t_val
    try:
        output = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
        sys.stdout.write(output)
    except subprocess.CalledProcessError:
        print('Execution of "%s" failed!\n' % cmd)
        sys.exit(1)
    
    "Call variants for rep 2"
    cmd = 'samtools mpileup --reference ' + ref_genome + ' -A -d 1000000 -F 0 -B -Q 0 ' + bam_rep2 + ' | ivar variants -p ' + out_rep2 + ' -q ' + q_val + ' -t ' + t_val
    try:
        output = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
        sys.stdout.write(output)
    except subprocess.CalledProcessError:
        print('Execution of "%s" failed!\n' % cmd)
        sys.exit(1)
        
    "Filter variants across technical replicants in iVar"
    tsv_rep1 = out_rep1 + '.tsv'
    tsv_rep2 = out_rep2 + '.tsv'
    cmd = 'ivar filtervariants -p ' + out_filtered + ' ' + tsv_rep1 + ' ' + tsv_rep2
    try:
        output = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
        sys.stdout.write(output)
    except subprocess.CalledProcessError:
        print('Execution of "%s" failed!\n' % cmd)
        sys.exit(1)
        
        
def correct_alt_freqs(tsv_file,out_file):
    
    "Recompute ALT variant frequencies based on depths of all variants and indels"
    
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
        
if __name__ == '__main__':
    
    batch = True # processing multiple samples?
    correct_freqs = False # correct variants freqs to take into account indel freqs
    
    "Set up paths to data files"
    ref_path = '/Volumes/GoogleDrive/My\ Drive/DeepTSWV/tswv_ref/'
    ref_fasta = 'TF2_consensus.fasta'
    #ref_fasta = 'tswv_ref_full.fasta' # Old reference genome
    ref_genome = ref_path + ref_fasta # Fasta file with GenBank reference genome
    
    "For multiple fastq files in a directory"
    if batch:
        
        main_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/samples/'
        
        "Get all directories in main_path"
        dir_list = os.listdir(main_path)
        
        for dir in dir_list:
            
            print(dir)
            
            if (dir[0] == '.'):
                continue
                
            os.chdir(main_path + dir) # change path to sample_path
            
            bam_rep1 = dir + '_r1.sorted.bam'
            bam_rep2 = dir + '_r2.sorted.bam'
            
            "Prefix for output files"
            out_rep1 = dir + '_r1_refTF2'
            out_rep2 = dir + '_r2_refTF2'
            out_filtered = dir + '_corrAltFreqs_pairRepFiltered'
            
            "If correcting ALT FREQS"
            if correct_freqs:
                tsv_file_rep1 = out_rep1 + '.tsv'
                tsv_file_rep2 = out_rep2 + '.tsv'
                out_file_rep1 = dir + '_r1_corrAltFreqs.tsv'
                out_file_rep2 = dir + '_r2_corrAltFreqs.tsv'
                correct_alt_freqs(tsv_file_rep1, out_file_rep1)
                correct_alt_freqs(tsv_file_rep2, out_file_rep2)
                out_rep1 = dir + '_r1_corrAltFreqs'
                out_rep2 = dir + '_r2_corrAltFreqs'
            
            call_variants(ref_genome,bam_rep1,bam_rep2,out_rep1,out_rep2,out_filtered)
            
            os.chdir(main_path)
            
            
    else:
        
        sample_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/samples/P2_L2_14'
        os.chdir(sample_path) # need to do this within directory if working from Google Drive otherwise paths break
    
        bam_rep1 = 'P2_L2_14_r1_refTF2.sorted.bam'
        bam_rep2 = 'P2_L2_14_r2_refTF2.sorted.bam'
    
        "Prefix for output files"
        tsv_file_rep1 = 'P2_L2_14_r1_refTF2.tsv'
        tsv_file_rep2 = 'P2_L2_14_r2_refTF2.tsv'
        out_filtered = 'P2_L2_14_corrAltFreqs_pairRepFiltered'
        
        "If correcting ALT FREQS"
        if correct_freqs:
            out_file_rep1 = 'P2_L2_14_r1_corrAltFreqs.tsv'
            out_file_rep2 = 'P2_L2_14_r2_corrAltFreqs.tsv'
            correct_alt_freqs(tsv_file_rep1, out_file_rep1)
            correct_alt_freqs(tsv_file_rep2, out_file_rep2)
            out_rep1 = 'P2_L2_14_r1_corrAltFreqs'
            out_rep2 = 'P2_L2_14_r2_corrAltFreqs'
    
        call_variants(ref_genome,bam_rep1,bam_rep2,out_rep1,out_rep2,out_filtered)