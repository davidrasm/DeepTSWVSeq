"""
Created on Sat Dec  7 16:45:42 2019

Compare within-host nucleotide diversity to species level diversity
Diversity is measured in terms of shannon entropy

@author: david
"""

from Bio import AlignIO
import sys
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
from Bio import SeqIO
import numpy as np
import pandas as pd

verbose = True

def shannon_entropy(list_input):
    """Calculate Shannon's Entropy per column of the alignment (H=-\sum_{i=1}^{M} P_i\,log_2\,P_i)"""

    unique_base = set(list_input)
    M   =  len(list_input)
    entropy_list = []
    # Number of residues in column
    for base in unique_base:
        n_i = list_input.count(base) # Number of residues of type i
        P_i = n_i/float(M) # n_i(Number of residues of type i) / M(Number of residues in column)
        entropy_i = P_i*(math.log(P_i,2))
        entropy_list.append(entropy_i)

    sh_entropy = -(sum(entropy_list))

    return sh_entropy


def shannon_entropy_list_msa(alignment):
    """Calculate Shannon Entropy across the whole MSA"""

    shannon_entropy_list = []
    for col_no in range(len(list(alignment[0]))):
        list_input = list(alignment[:, col_no])
        shannon_entropy_list.append(shannon_entropy(list_input))

    return shannon_entropy_list

def running_mean(l, N):
    sum = 0
    result = list(0 for x in l)

    for i in range( 0, N ):
        sum = sum + l[i]
        result[i] = sum / (i+1)

    for i in range( N, len(l) ):
        sum = sum - l[i-N] + l[i]
        result[i] = sum / N

    return result

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

def compute_species_entropy(msa,runningmean):
    
    alignment = AlignIO.read(msa, 'fasta')
    
    # Do a little sanity checking:
    seq_lengths_list = []
    for record in alignment:
       seq_lengths_list.append(len(record))
    
    seq_lengths = set(seq_lengths_list)
    
    if verbose: print("Alignment length is:" + str(list(seq_lengths)))
    
    if len(seq_lengths) != 1:
        sys.stderr.write("Your alignment lengths aren't equal. Check your alignment file.")
        sys.exit(1)
    
    index = range(1, list(seq_lengths)[0]+1)
    
    sel = shannon_entropy_list_msa(alignment)
    
    sel = running_mean(sel, runningmean)
    
    return index, sel

def get_variant_df(sample):
    
    "Get variants in sample"
    main_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/samples/'
    tsv_file = main_path + sample + '/' + sample + '_refTF2_pairRepFiltered.tsv'
    tsv_rep1 = sample + '_r1_refTF2.tsv'
    tsv_rep2 = sample + '_r2_refTF2.tsv'
    vdf = pd.read_table(tsv_file, sep="\t")
    
    "Get average frequency across replicates"
    alt_freq_rep1 = 'ALT_FREQ_' + tsv_rep1
    alt_freq_rep2 = 'ALT_FREQ_' + tsv_rep2
    vdf['AVG_ALT_FREQ'] = (vdf[alt_freq_rep1] + vdf[alt_freq_rep2]) / 2
    
    return vdf

def compute_within_entropy(samples,ref):
    
    seq = ref.seq
    ref_id = ref.id
    
    esa = np.zeros((1,len(seq)))
    for s in samples:
        vdf = get_variant_df(s)
        sel = []
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
            
            entropy_list = []
            for freq in freqs:
                P_i = freq
                if P_i > 0:
                    entropy_i = P_i*(math.log(P_i,2))
                else:
                    entropy_i = 0
                entropy_list.append(entropy_i)
    
            sh_entropy = -(sum(entropy_list))
            sel.append(sh_entropy)
        
        sel = running_mean(sel, runningmean)
        sel_array = np.array(sel,ndmin=2)
        esa = np.concatenate((esa,sel_array))
    
    sel = np.mean(esa,axis=0)
    return sel
    
    
"Set reference genome for indexing purposes"  
ref_fasta = '../tswv_ref/TF2_consensus.fasta'
records = list(SeqIO.parse(ref_fasta, "fasta"))
samples = ['P0','WF_P1','P1_L1_7','P1_L1_14','P1_L1_18','P1_L1_28','WF_P2_L1','P2_L1_7','P2_L1_14','P2_L1_21','P2_L1_29','WF_P3_L1','P3_L1_7','P3_L1_14','P3_L1_21','P3_L1_28','WF_P4_L1','P4_L1_7','P4_L1_14','P4_L1_21','P4_L1_28','WF_P5_L1','P5_L1_14','P5_L1_21','P5_L1_28'] # all line 1

data_path = '../DeepSeqResults/speciesDiversity/'
runningmean = 100 # value used in paper

"Compute and save within-host entropies since these are expensive to compute for many samples"
#host_entropy_segS = compute_within_entropy(samples,records[0])
#df = pd.DataFrame({'host_entropy_segS': host_entropy_segS})
#df.to_csv('tswv_within-host_segS_entropies.csv', sep=',') 
df = pd.read_csv(data_path + 'tswv_within-host_segS_entropies.csv')
host_entropy_segS = df['host_entropy_segS'].values

#host_entropy_segM = compute_within_entropy(samples,records[1])
#host_entropy_segM = host_entropy_segM[:-1] # not sure why this isn't the same length as species_entropy_M
#df = pd.DataFrame({'host_entropy_segM': host_entropy_segM})
#df.to_csv('tswv_within-host_segM_entropies.csv', sep=',') 
df = pd.read_csv(data_path + 'tswv_within-host_segM_entropies.csv')
host_entropy_segM = df['host_entropy_segM'].values

#host_entropy_segL = compute_within_entropy(samples,records[2])
#host_entropy_segL = host_entropy_segL[:-1]
#df = pd.DataFrame({'host_entropy_segL': host_entropy_segL})
#df.to_csv('tswv_within-host_segL_entropies.csv', sep=',')
df = pd.read_csv(data_path + 'tswv_within-host_segL_entropies.csv')
host_entropy_segL = df['host_entropy_segL'].values

"Compute within-host and species-level site entropies for seg S"
msa = data_path + 'TSWV_segS_globalSeqAccessions.fasta'
sites_segS, species_entropy_segS = compute_species_entropy(msa,runningmean)

"Compute within-host and species-level site entropies for seg M"
msa = data_path + 'TSWV_segM_globalSeqAccessions.fasta'
sites_segM, species_entropy_segM = compute_species_entropy(msa,runningmean)


"Compute within-host and species-level site entropies for seg L"
msa = data_path + 'TSWV_RdRp_globalSeqAccessions.fasta'
sites_segL, species_entropy_segL = compute_species_entropy(msa,runningmean)

"Plotting"
sns.set()
fig, axs = plt.subplots(3, 1)

#sites = list(range(df_rep1_segL["Depth"].size))
axs[0].plot(sites_segS, species_entropy_segS,label = "Species-level")
axs[0].plot(sites_segS, host_entropy_segS,label = "Within-host")
axs[0].set_ylabel('Diversity')
axs[0].set_xlim([0,2915])
axs[0].legend()
axs[0].set_title('Segment S')
add_gene(axs[0], 'N', 89, 1483) # Add  N 89..1483
add_gene(axs[0], 'NSs', 1987, 2763) # Add NSs 1987..2763

axs[1].plot(sites_segM, species_entropy_segM,label = "Species-level")
axs[1].plot(sites_segM, host_entropy_segM,label = "Within-host")
axs[1].set_ylabel('Diversity')
axs[1].set_xlim([0,4820])
axs[1].set_title('Segment M')
add_gene(axs[1], 'NSm', 101, 1009) # Add  NSm 101..1009
"If treating GP as one protein"
#add_gene(axs[1], 'GP', 1330, 4737) # Add GP 1330..4737
"If dividing GP into Gn/Gc"
add_gene(axs[1], r'$G_C$', 1330, 3285)
add_gene(axs[1], r'$G_N$', 3286, 4737)

axs[2].plot(sites_segL, species_entropy_segL,label = "Species-level")
axs[2].plot(sites_segL, host_entropy_segL,label = "Within-host")
#axs[2].set_xlabel('Position')
axs[2].set_ylabel('Diversity')
axs[2].set_xlim([0,8896])
axs[2].set_title('Segment L')
add_gene(axs[2], 'RdRp', 34, 8661) # Add RdRp 34..8661

fig.set_size_inches(8, 8)
fig.tight_layout()
plt.show()

out_file = 'tswv_withinhost_vs_species_diversity_GnGc.png'
fig.savefig(out_file, dpi=200)
