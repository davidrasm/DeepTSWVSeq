"""
Created on Fri Nov 15 08:23:09 2019

Embed primers in bed file and use iVar to trim primers from sequences
Note: Never got ivar trim function to work properly for TSWV reads

System dependencies:
    bowtie2
    samtools
    bedtools

@author: david
"""
import sys
import os
import subprocess

def embed_primers():
    
    "TSWV primers do not align to unique positions in reference genome "
    "Need to manually edit resulting bed file to get SF and MF primer sites correct"

    "Set path to ref files"
    ref_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/tswv_ref/'
    os.chdir(ref_path)
    #ref_str = 'tswv_ref'
    ref_file = 'tswv_TF2_consensus'
    
    "Set path to primer fasta file"
    primer_fastq = 'tswv_primers.fastq'
    out_str = 'tswv_primers'
    
    "Align primer sequences against reference genome"
    "Have to reset seed length using -L arg in bowtie to get all primers to align"
    sam_file = out_str + '.sam'
    cmd = 'bowtie2 -D 15 -R 2 -N 0 -L 5 -i S,1,1.15 -x ' + ref_file + ' -U ' + primer_fastq + ' -S ' + sam_file 
    try:
        output = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
        sys.stdout.write(output)
    except subprocess.CalledProcessError:
        print('Execution of "%s" failed!\n' % cmd)
        sys.exit(1)
    
    "Convert sam file to bam files in samtools"
    bam_file = out_str + '.bam'
    cmd = 'samtools view -bS ' + sam_file + ' > ' + bam_file
    try:
        output = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
        sys.stdout.write(output)
    except subprocess.CalledProcessError:
        print('Execution of "%s" failed!\n' % cmd)
        sys.exit(1)
    
    bed_file = out_str + '.bed'
    cmd = 'bedtools bamtobed -i ' + bam_file + ' > ' + bed_file
    try:
        output = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
        sys.stdout.write(output)
    except subprocess.CalledProcessError:
        print('Execution of "%s" failed!\n' % cmd)
        sys.exit(1)
        
def trim_reads(sample_path,bed_file,bam_file,out_str):
    
    os.chdir(sample_path)
    
    new_bam_file = out_str + '.trimmed.bam'
    cmd = 'ivar trim -b ' + bed_file + ' -p ' + new_bam_file + ' -i ' + bam_file
    try:
        output = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
        sys.stdout.write(output)
    except subprocess.CalledProcessError:
        print('Execution of "%s" failed!\n' % cmd)
        sys.exit(1)
    
    "Sort trimmed BAM file"
    sorted_bam_file =  out_str + '.trimmed.sorted.bam'    
    cmd = 'samtools sort ' + new_bam_file + ' -o ' + sorted_bam_file
    try:
        output = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
        sys.stdout.write(output)
    except subprocess.CalledProcessError:
        print('Execution of "%s" failed!\n' % cmd)
        sys.exit(1)
    
    "Do we need to re-index bam file?"
    #samtools index Z52_a.trimmed.sorted.bam
    
if __name__ == '__main__':
    
    #embed_primers()
    
    bed_path = '/Volumes/GoogleDrive/My\ Drive/DeepTSWV/tswv_ref/'
    bed_file = bed_path + 'tswv_primers_edited.bed' #manually edited bed file
    
    main_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/samples/'
    sample_path = main_path + 'P1_L1_7'
    bam_file = 'P1_L1_7_r1.sorted.bam' # file to trim
    out_str = 'P1_L1_7_r1' # output file
    trim_reads(sample_path,bed_file,bam_file,out_str)