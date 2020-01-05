"""
Created on Mon Jun 10 17:31:24 2019

Plot frequency time series for specific variants for multiple lines

Used to create Figure 7 in paper showing freq time series of all three lines

@author: david
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from call_aa_variants import call_aa_variants
import seaborn as sns


def process_variants(vdf):
    
    "Filter by variant type"
    vdf['INDEL'] = vdf.ALT.str.contains('^\+|-') # regexp means first character is '+' or '-' char
    if (var_type == 'SNV'):
       vdf = vdf.loc[vdf['INDEL'] == False]
    elif (var_type == 'AA'):
        aa_variant,aa_subst,cds_annotation = call_aa_variants(vdf)
        vdf['AA_VAR'] = np.array(aa_variant)
        vdf['AA_SUBST'] = np.array(aa_subst)
        vdf['CDS'] = np.array(cds_annotation)
        vdf = vdf.loc[vdf['AA_VAR'] == True]
    elif (var_type == 'INDEL'):
        vdf = vdf.loc[vdf['INDEL'] == True]
    
    if seg:
        variants = vdf[(vdf['REGION'] == seg)]
    else:
        variants = vdf   
    
    return variants

"Specify mutation/region to track"
seg = None #'TSWV_segL'
var_type = 'AA' # 'SNV', 'INDEL' or 'AA"

"Plot by enrichment"
#enrich_by_host = True
#enrichment_delta = 0.05 # enrichment frequency threshold
#enrichment_host = 'PLANT' # 'PLANT', 'THRIP' or 'ANY'

lines = [1,2,3]
main_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/DeepSeqScripts/'

#out_file = 'L1-3_common_enrichedAAVs_delta0.05_corrAltFreqs_freqSeries.png'
out_file = 'L1-3_common_enrichedAAVs_NSm_corrAltFreqs_freqSeries_GnGc.png'

tsv_file = main_path + 'tswv_L1-3_masterVarList_corrAltFreqs_withDepths.csv'
vdf = pd.read_table(tsv_file, sep="\t")
variants_L1 = process_variants(vdf)
samples_L1 = ['P0','WF_P1','P1_L1_7','P1_L1_14','P1_L1_18','P1_L1_28','WF_P2_L1','P2_L1_7','P2_L1_14','P2_L1_21','P2_L1_29','WF_P3_L1','P3_L1_7','P3_L1_14','P3_L1_21','P3_L1_28','WF_P4_L1','P4_L1_7','P4_L1_14','P4_L1_21','P4_L1_28','WF_P5_L1','P5_L1_14','P5_L1_21','P5_L1_28'] # all line 1

tsv_file = main_path + 'tswv_L2_masterVarList_refTF2_withDepths.csv'
vdf = pd.read_table(tsv_file, sep="\t")
variants_L2 = process_variants(vdf)
samples_L2 = ['P0','WF_P1','P1_L2_7','P1_L2_14','P1_L2_18','P1_L2_28','WF_P2_L2','P2_L2_7','P2_L2_14','P2_L2_21','P2_L2_29','WF_P3_L2','P3_L2_7','P3_L2_14','P3_L2_21','P3_L2_28','WF_P4_L2','P4_L2_14','P4_L2_21','P4_L2_28','WF_P5_L2','P5_L2_7','P5_L2_14','P5_L2_21','P5_L2_28']

tsv_file = main_path + 'tswv_L3_masterVarList_refTF2_withDepths.csv'
vdf = pd.read_table(tsv_file, sep="\t")
variants_L3 = process_variants(vdf)
samples_L3 = ['P0','WF_P1','P1_L3_7','P1_L3_14','P1_L3_18','P1_L3_28','WF_P2_L3','P2_L3_7','P2_L3_14','P2_L3_21','P2_L3_29','WF_P3_L3','P3_L3_7','P3_L3_14','P3_L3_21','P3_L3_28','WF_P4_L3','P4_L3_7','P4_L3_14','P4_L3_21','P4_L3_28','WF_P5_L3','P5_L3_7','P5_L3_14','P5_L3_21','P5_L3_28']  

line_variants = [variants_L1,variants_L2,variants_L3]
line_samples = [samples_L1,samples_L2,samples_L3]

"Set variants to plot"
cds = ['NSs','Gc','RdRp','RdRp','RdRp']
muts = ['P181T','E362D','K566R','I1373M','I1602V']

#cds = ['NSm','NSm','NSm','NSm','NSm','NSm']
#muts = ['V17G','Q21L','N22S','N22I','L25F','S30C']

"Set color palette in seaborn"
sns.set(style="darkgrid")
sns.set_palette(sns.hls_palette(len(muts), l=.4, s=.8))
fig, axs = plt.subplots(3, 1)

# Loop through each line
for line in lines:

    variants = line_variants[line-1]
    samples = line_samples[line-1]
    
    for m in range(len(muts)):
        vseries = variants[(variants['CDS'] == cds[m]) & (variants['AA_SUBST'] == muts[m])]
        freqs = [vseries['ALT_FREQ_' + s].values[0] for s in samples]        
        if (var_type == 'AA'):
            axs[line-1].plot(samples, freqs, label = vseries['CDS'].values[0] + ' ' + vseries['AA_SUBST'].values[0])
        else:
            axs[line-1].plot(samples, freqs, label = vseries['REF'] + str(vseries['POS']) + vseries['ALT'])
        
    "Line specific plottign args"
    if line == 1:
        my_colors = ['g', 'k', 'b','b','b','b','k', 'g','g','g','g','k','b','b','b','b','k','g','g','g','g','k','b','b','b'] # line 1
        labels = [r'$P_0$', r'$T_1$',r'$P_1^7$',r'$P_1^{14}$',r'$P_1^{18}$',r'$P_1^{28}$',r'$T_2$',r'$P_2^7$',r'$P_2^{14}$',r'$P_2^{21}$',r'$P_2^{29}$',r'$T_3$',r'$P_3^7$',r'$P_3^{14}$',r'$P_3^{21}$',r'$P_3^{28}$',r'$T_4$',r'$P_4^7$',r'$P_4^{14}$',r'$P_4^{21}$',r'$P_4^{28}$',r'$T_5$',r'$P_5^{14}$',r'$P_5^{21}$',r'$P_5^{28}$']
    elif line == 2:
        my_colors = ['g', 'k', 'g','g','g','g','k', 'g','g','g','g','k','g','g','g','g','k','g','g','g','k','b','b','b','b'] # line 2
        labels = [r'$P_0$', r'$T_1$',r'$P_1^7$',r'$P_1^{14}$',r'$P_1^{18}$',r'$P_1^{28}$',r'$T_2$',r'$P_2^7$',r'$P_2^{14}$',r'$P_2^{21}$',r'$P_2^{29}$',r'$T_3$',r'$P_3^7$',r'$P_3^{14}$',r'$P_3^{21}$',r'$P_3^{28}$',r'$T_4$',r'$P_4^{14}$',r'$P_4^{21}$',r'$P_4^{28}$',r'$T_5$',r'$P_5^{7}$',r'$P_5^{14}$',r'$P_5^{21}$',r'$P_5^{28}$']
    elif line == 3:
        my_colors =  ['b', 'k', 'b','b','b','b','k', 'b','b','b','b','k','b','b','b','b','k','b','b','b','b','k','g','g','g','g'] # line 2
        labels = [r'$P_0$', r'$T_1$',r'$P_1^7$',r'$P_1^{14}$',r'$P_1^{18}$',r'$P_1^{28}$',r'$T_2$',r'$P_2^7$',r'$P_2^{14}$',r'$P_2^{21}$',r'$P_2^{29}$',r'$T_3$',r'$P_3^7$',r'$P_3^{14}$',r'$P_3^{21}$',r'$P_3^{28}$',r'$T_4$',r'$P_4^7$',r'$P_4^{14}$',r'$P_4^{21}$',r'$P_4^{28}$',r'$T_5$',r'$P_5^{7}$',r'$P_5^{14}$',r'$P_5^{21}$',r'$P_5^{28}$']
    
    "Color tick labels to emphasize host switches"
    for ticklabel, tickcolor in zip(axs[line-1].get_xticklabels(), my_colors):
        ticklabel.set_color(tickcolor)

    "Better labeling for xticks"
    axs[line-1].set_xticklabels(labels,rotation='vertical')   

    #axs[line-1].set_xlabel('Sample')
    axs[line-1].set_ylabel('Frequency')
    axs[line-1].grid(True)

axs[2].legend(loc='upper center', bbox_to_anchor=(0.5, -0.24),
          fancybox=True, shadow=False, ncol=5)

#axs[2].set_title('Emilia enriched amino acid variants')

fig.set_size_inches(9, 8)
fig.tight_layout()
plt.show()

fig.savefig(out_file, dpi=200,bbox_inches='tight')