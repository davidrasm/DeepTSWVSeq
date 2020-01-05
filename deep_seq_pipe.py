"""
Created on Mon May 27 13:20:31 2019

Run pipeline to map deep seq reads to ref genome. Takes input fastq files and turns them into a SAM/BAM alignment.

System dependencies:
    bowtie2
    samtools
    bwa (if used instead of bowtie2 for read mapping)

@author: david
"""
import os
import sys
import subprocess
import glob
import shutil

def build_ref(ref_path,ref_str):
    
    "Set current working directory to sample data"
    #os.chdir(ref_path)
    
    "Build reference to align against in bowtie2"
    #ref_genome = ref_path + 'tswv_ref_full.fasta' # Fasta file with GenBank reference genome
    ref_genome = ref_path + ref_str + '.fasta' # Fasta file with GenBank reference genome
    cmd = 'bowtie2-build ' + ref_genome + ' ' + ref_path + ref_str
    try:
        output = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
        sys.stdout.write(output)
    except subprocess.CalledProcessError:
        print('Execution of "%s" failed!\n' % cmd)
        sys.exit(1)
        

    
def index_ref_bwa(ref_path,ref_str):
    
    ref_genome = ref_path + 'TF2_consensus.fasta' # Fasta file with GenBank reference genome
    ref_seqs_db =  ref_path + ref_str + '_bwa'
    
    "Build reference in bwa"
    #bwa index -p ref-seqs-db ref-seqs.fasta
    cmd = 'bwa index -p ' + ref_seqs_db + ' ' + ref_genome
    try:
        output = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
        sys.stdout.write(output)
    except subprocess.CalledProcessError:
        print('Execution of "%s" failed!\n' % cmd)
        sys.exit(1)

def run_pipe(ref_file,sample_path,fastq_R1,fastq_R2,out_str):

    "Set current working directory to sample data"
    os.chdir(sample_path)
        
    "Align paired end reads against ref in bowtie2"
    sam_file = out_str + '.sam'
    
    "Different preset alignment options in bowtie2 - we used sensitive-local for all results in paper" 
    #cmd = 'bowtie2 -x ' + ref_file + ' -1 ' + fastq_R1 + ' -2 ' + fastq_R2 + ' -S ' + sam_file 
    #cmd = 'bowtie2 --very-sensitive -x ' + ref_file + ' -1 ' + fastq_R1 + ' -2 ' + fastq_R2 + ' -S ' + sam_file
    cmd = 'bowtie2 --sensitive-local -x ' + ref_file + ' -1 ' + fastq_R1 + ' -2 ' + fastq_R2 + ' -S ' + sam_file
    
    "If using bwa instead of bowtie2"
    #bwa_ref_db = ref_file + '_bwa'
    #cmd = 'bwa mem ' + bwa_ref_db + ' ' + fastq_R1 + ' ' + fastq_R2 + ' > ' + sam_file
    #Align reads in fastq file against ref sequences in database:
    #bwa mem db.prefix reads.fq > out.sam
    
    try:
        output = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
        sys.stdout.write(output)
    except subprocess.CalledProcessError:
        print('Execution of "%s" failed!\n' % cmd)
        sys.exit(1)
    
    "Convert sam file to bam files in samtools"
    bam_file = out_str + '.bam'
    sorted_bam_file =  out_str + '.sorted.bam'
    cmd = 'samtools view -bS ' + sam_file + ' > ' + bam_file
    try:
        output = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
        sys.stdout.write(output)
    except subprocess.CalledProcessError:
        print('Execution of "%s" failed!\n' % cmd)
        sys.exit(1)
        
    cmd = 'samtools sort ' + bam_file + ' -o ' + sorted_bam_file
    try:
        output = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
        sys.stdout.write(output)
    except subprocess.CalledProcessError:
        print('Execution of "%s" failed!\n' % cmd)
        sys.exit(1)
        
    "Get depth of coverage using samtools"
    depth_log_file = out_str + '_depth.log'
    cmd = 'samtools depth ' + sorted_bam_file + ' > ' + depth_log_file
    try:
        output = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
        sys.stdout.write(output)
    except subprocess.CalledProcessError:
        print('Execution of "%s" failed!\n' % cmd)
        sys.exit(1)
        
    "Get consensus sequence"
    #consensus_file = out_str + '_consensus.fasta'
    #cmd = 'samtools mpileup -A -d 1000000 -Q 0 -F 0 ' + sorted_bam_file + ' | ivar consensus -p ' + consensus_file
    #try:
    #    output = subprocess.check_output(cmd, shell=True,stderr=subprocess.STDOUT)
    #    sys.stdout.write(output)
    #except subprocess.CalledProcessError:
    #    print('Execution of "%s" failed!\n' % cmd)
    #    sys.exit(1)


def process_sample():
    
    "Need to do these for both RT replicates"
    sample_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/samples/P2_L2_14' # omit '\' in My Drive directory -- will confuse os command in python
    #os.chdir(sample_path)
    
    fastq_R1 = 'P2_L2_14_r1_S43_L001_R1_001.fastq'
    fastq_R2 = 'P2_L2_14_r1_S43_L001_R2_001.fastq'

    "Prefix for output files"
    out_str = 'P2_L2_14_r1_self-concensus'
    
    run_pipe(ref_file,sample_path,fastq_R1,fastq_R2,out_str)
    
    "Process other RT rep"
    fastq_R1 = 'P2_L2_14_r2_S44_L001_R1_001.fastq'
    fastq_R2 = 'P2_L2_14_r2_S44_L001_R2_001.fastq'

    "Prefix for output files"
    out_str = 'P2_L2_14_r2_self-concensus'
    
    run_pipe(ref_file,sample_path,fastq_R1,fastq_R2,out_str)
    
def process_batch():
    
    main_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/samples/'
    
    "Get all directories in main_path"
    dir_list = os.listdir(main_path)
    
    for dir in dir_list:
        
        print(dir)
            
        if (dir[0] == '.'):
            continue
                
        os.chdir(main_path + dir) # change path to sample_path
        dir_name = main_path + dir
        
        fastq_list = glob.glob('*.fastq') # get all fastq files in dir
    
        "This should call run_pipe on both RT replicates"
        for file in fastq_list:
        
            if ('r1' in file) and ('R1' in file): 
                fastq_r1_R1 = file
            elif ('r1' in file) and ('R2' in file):
                fastq_r1_R2 = file
            elif ('r2' in file) and ('R1' in file):
                fastq_r2_R1 = file
            elif ('r2' in file) and ('R2' in file):
                fastq_r2_R2 = file
            else:
                print("Fastq file not recognized")
                    
        "For RT rep1"
        out_str = dir + '_r1'
        run_pipe(ref_file,dir_name,fastq_r1_R1,fastq_r1_R2,out_str)
        
        "For RT rep2"
        out_str = dir + '_r2'
        run_pipe(ref_file,dir_name,fastq_r2_R1,fastq_r2_R2,out_str)
        
        os.chdir(main_path) # change path back to main_path
    
def process_unsorted_batch():
    
    "Get all .fastq files in directory"
    main_path = '/Users/david/Desktop/fastqfiles_MiSeq2_0129_Rasmussen'
    os.chdir(main_path) # change path to sample_path
    fastq_list = glob.glob('*.fastq') # get all fastq files in dir
    
    "This should call run_pipe on both RT replicates"
    for file in fastq_list:
        
        if not os.path.exists(file): # if we've already moved it
            continue
        
        dir_name = file.split('_r')[0] # parse short name from .fastq file
        
        if not os.path.isdir(dir_name):
            os.mkdir(dir_name)
        
        if 'R1' in file: # if forward reads
            fastq_R1 = file
            fastq_R2 = file.replace('R1', 'R2')
        else:
            fastq_R1 = file.replace('R2', 'R1')
            fastq_R2 = file
            
        # Move .fastq files to new dir
        shutil.move(fastq_R1, dir_name + "/" + fastq_R1)
        shutil.move(fastq_R2, dir_name + "/" + fastq_R2)
        
        if 'r1' in file:
            out_str = dir_name + '_r1' # if RT rep 1
        else:
            out_str = dir_name + '_r2' # if RT rep 2
        run_pipe(ref_file,dir_name,fastq_R1,fastq_R2,out_str)
        os.chdir(main_path) # change path back to main_path
    
    

if __name__ == '__main__':
    
    batch = False # processing multiple samples?
    
    "Set up path to ref files"
    #path = os.path.join(os.environ['HOME'], 'Desktop/DeepTSWV/tswv_ref') # relative path
    #ref_path = '/Volumes/GoogleDrive/My\ Drive/DeepTSWV/tswv_ref/'
    ref_path = '/Volumes/GoogleDrive/My\ Drive/DeepTSWV/samples/P2_L2_14/'
    
    #ref_str = 'tswv_ref'
    #ref_str = 'tswv_TF2_consensus'
    ref_str = 'P2_L2_14_consensus'
    
    "If using bowtie2"
    build_ref(ref_path,ref_str)
    
    "If using bwa"
    #index_ref_bwa(ref_path,ref_str)
    
    ref_file = ref_path + ref_str
    
    #process_batch()
    process_sample()
    
    

    
    

    
    
    