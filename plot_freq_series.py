"""
Created on Mon Jun 10 17:31:24 2019

Plot frequency time series for particular variants across samples

Used to create Figures 6 in paper

@author: david
"""

import pandas as pd
import matplotlib.pyplot as plt


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

main_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/DeepSeqScripts/'

#out_file = 'L1_emiliaDatura_enrichedAAVs_corrAltFreqs_freqSeries.png'
out_file = 'L1_plantThrip_enrichedAAVs_corrAltFreqs_freqSeries_GnGc.png'

tsv_file = main_path + 'tswv_L1-3_masterVarList_corrAltFreqs_withDepths.csv' # final master list used for results in manuscript
line = 1

vdf = pd.read_table(tsv_file, sep="\t")
variants = process_variants(vdf)

if line == 1:
    samples = ['P0','WF_P1','P1_L1_7','P1_L1_14','P1_L1_18','P1_L1_28','WF_P2_L1','P2_L1_7','P2_L1_14','P2_L1_21','P2_L1_29','WF_P3_L1','P3_L1_7','P3_L1_14','P3_L1_21','P3_L1_28','WF_P4_L1','P4_L1_7','P4_L1_14','P4_L1_21','P4_L1_28','WF_P5_L1','P5_L1_14','P5_L1_21','P5_L1_28'] # all line 1
elif line == 2:
    samples = ['P0','WF_P1','P1_L2_7','P1_L2_14','P1_L2_18','P1_L2_28','WF_P2_L2','P2_L2_7','P2_L2_14','P2_L2_21','P2_L2_29','WF_P3_L2','P3_L2_7','P3_L2_14','P3_L2_21','P3_L2_28','WF_P4_L2','P4_L2_14','P4_L2_21','P4_L2_28','WF_P5_L2','P5_L2_7','P5_L2_14','P5_L2_21','P5_L2_28']
elif line == 3:
    samples = ['P0','WF_P1','P1_L3_7','P1_L3_14','P1_L3_18','P1_L3_28','WF_P2_L3','P2_L3_7','P2_L3_14','P2_L3_21','P2_L3_29','WF_P3_L3','P3_L3_7','P3_L3_14','P3_L3_21','P3_L3_28','WF_P4_L3','P4_L3_7','P4_L3_14','P4_L3_21','P4_L3_28','WF_P5_L3','P5_L3_7','P5_L3_14','P5_L3_21','P5_L3_28']

"Set variants to plot for Emilia/Datura comparisons"
#cds = ['NSs','NSm','NSm','RdRp','RdRp','RdRp','RdRp']
#muts = ['I79M','V17G','N22S','K566R','I1373M','I1602V','C2148Y']

"Set variants to plot for Plant/Thrip comparisons"
cds = ['NSs','NSs','Gc','Gn','RdRp','RdRp','RdRp','RdRp','RdRp','RdRp','RdRp','RdRp','RdRp','RdRp','RdRp']
muts = ['P181T','I79M','E362D','V141F','I226S','R290S','N389S','R495Q','K566R','N813D','K863R','I1373M','I1602V','A1799T','C2148Y']

"Set color palette in seaborn"
sns.set(style="darkgrid")
sns.set_palette(sns.hls_palette(len(muts), l=.4, s=.8))
fig, ax = plt.subplots(1, 1)
    
for m in range(len(muts)):
    vseries = variants[(variants['CDS'] == cds[m]) & (variants['AA_SUBST'] == muts[m])]
    freqs = [vseries['ALT_FREQ_' + s].values[0] for s in samples]        
    if (var_type == 'AA'):
        ax.plot(samples, freqs, label = vseries['CDS'].values[0] + ' ' + vseries['AA_SUBST'].values[0])
    else:
        ax.plot(samples, freqs, label = vseries['REF'] + str(vseries['POS']) + vseries['ALT'])
    
"Line specific plottign args"
if line == 1:
    my_colors = ['g', 'k', 'b','b','b','b','k', 'g','g','g','g','k','b','b','b','b','k','g','g','g','g','k','b','b','b'] # line 1
    labels = [r'$P_0$', r'$T_1$',r'$P_1^7$',r'$P_1^{14}$',r'$P_1^{18}$',r'$P_1^{28}$',r'$T_2$',r'$P_2^7$',r'$P_2^{14}$',r'$P_2^{21}$',r'$P_2^{29}$',r'$T_3$',r'$P_3^7$',r'$P_3^{14}$',r'$P_3^{21}$',r'$P_3^{28}$',r'$T_4$',r'$P_4^7$',r'$P_4^{14}$',r'$P_4^{21}$',r'$P_4^{28}$',r'$T_5$',r'$P_5^{14}$',r'$P_5^{21}$',r'$P_5^{28}$']
elif line == 2:
    my_colors = ['g', 'k', 'g','g','g','g','k', 'g','g','g','g','k','g','g','g','g','k','g','g','g','k','g','g','g','g'] # line 2
    labels = [r'$P_0$', r'$T_1$',r'$P_1^7$',r'$P_1^{14}$',r'$P_1^{18}$',r'$P_1^{28}$',r'$T_2$',r'$P_2^7$',r'$P_2^{14}$',r'$P_2^{21}$',r'$P_2^{29}$',r'$T_3$',r'$P_3^7$',r'$P_3^{14}$',r'$P_3^{21}$',r'$P_3^{28}$',r'$T_4$',r'$P_4^{14}$',r'$P_4^{21}$',r'$P_4^{28}$',r'$T_5$',r'$P_5^{7}$',r'$P_5^{14}$',r'$P_5^{21}$',r'$P_5^{28}$']
elif line == 3:
    my_colors =  ['b', 'k', 'b','b','b','b','k', 'b','b','b','b','k','b','b','b','b','k','b','b','b','b','k','b','b','b','b'] # line 2
    labels = [r'$P_0$', r'$T_1$',r'$P_1^7$',r'$P_1^{14}$',r'$P_1^{18}$',r'$P_1^{28}$',r'$T_2$',r'$P_2^7$',r'$P_2^{14}$',r'$P_2^{21}$',r'$P_2^{29}$',r'$T_3$',r'$P_3^7$',r'$P_3^{14}$',r'$P_3^{21}$',r'$P_3^{28}$',r'$T_4$',r'$P_4^7$',r'$P_4^{14}$',r'$P_4^{21}$',r'$P_4^{28}$',r'$T_5$',r'$P_5^{7}$',r'$P_5^{14}$',r'$P_5^{21}$',r'$P_5^{28}$']

"Color tick labels to emphasize host switches"
for ticklabel, tickcolor in zip(plt.gca().get_xticklabels(), my_colors):
    ticklabel.set_color(tickcolor)

"Better labeling for xticks"
ax.set_xticklabels(labels,rotation='vertical')   

#axs[line-1].set_xlabel('Sample')
ax.set_ylabel('Frequency')
ax.grid(True)
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.24),
          fancybox=True, shadow=False, ncol=5)

#axs[2].set_title('Emilia enriched amino acid variants')

fig.set_size_inches(9, 4)
fig.tight_layout()
plt.show()

fig.savefig(out_file, dpi=200,bbox_inches='tight')