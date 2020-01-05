# DeepTSWVSeq
 Python scipts used to analyze TSWV deep sequence data in:
 
 Ruark-Seward, Kennedy & Rasmussen. (2020). Evolutionary dynamics of *Tomato spotted wilt virus* within and between alternate plant hosts and thrips. *in prep*
 
 Main scripts used for data processing/analysis include:
 
 ***call_rep_variants.py:*** calls SNV/indel variants from paired sequence replicates using iVar
 
 ***deep_seq_pipe.py:*** maps short sequence reads to reference genome using bowtie2 and SAMtools to create BAM/SAM files
 
 ***estimate_fitness.py:*** estimates maximum likelihood fitness effects of variants across hosts based on variant growth rates
 
 ***make_master_table.py:*** creates a master table of variant frequencies across all sampled time points in order to plot time series of variant frequencies.
