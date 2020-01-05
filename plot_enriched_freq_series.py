"""
Created on Mon Jun 10 17:31:24 2019

Plot frequency time series for variants enriched in plants or thrips

@author: david
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from call_aa_variants import call_aa_variants
import seaborn as sns

"Specify mutation/region to track"
seg = None #'TSWV_segL'
gene = None
var_type = 'AA' # 'SNV', 'INDEL' or 'AA"

"Plot by enrichment"
enrich_by_host = True
enrichment_delta = 0.05 # enrichment frequency threshold
enrichment_host = 'ANY' # 'PLANT', 'THRIP' or 'ANY'

line = 1
main_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/DeepSeqScripts/'
#tsv_file = main_path + 'tswv_L1_masterVarList_refTF2_withDepths.csv'
tsv_file = main_path + 'tswv_L1-3_masterVarList_corrAltFreqs_withDepths.csv' # final master list used for results in manuscript
vdf = pd.read_table(tsv_file, sep="\t")

out_file = 'L1_plantThripEnrichedAAVs_corrAltFreqs_freqSeries.png'

"Line 1 samples (alternating)"
samples = ['P0','WF_P1','P1_L1_7','P1_L1_14','P1_L1_18','P1_L1_28','WF_P2_L1','P2_L1_7','P2_L1_14','P2_L1_21','P2_L1_29','WF_P3_L1','P3_L1_7','P3_L1_14','P3_L1_21','P3_L1_28','WF_P4_L1','P4_L1_7','P4_L1_14','P4_L1_21','P4_L1_28','WF_P5_L1','P5_L1_14','P5_L1_21','P5_L1_28'] # all line 1
short_names = [s.replace('_L1','') for s in samples]

thrip_samples = ['WF_P1','WF_P2_L1','WF_P3_L1','WF_P4_L1','WF_P5_L1']
plant_samples = ['P0','P1_L1_7','P1_L1_14','P1_L1_18','P1_L1_28','P2_L1_7','P2_L1_14','P2_L1_21','P2_L1_29','P3_L1_7','P3_L1_14','P3_L1_21','P3_L1_28','P4_L1_7','P4_L1_14','P4_L1_21','P4_L1_28','P5_L1_14','P5_L1_21','P5_L1_28']

"Actually Emilia samples"
#thrip_samples = ['P0','P2_L1_7','P2_L1_14','P2_L1_21','P2_L1_29','P4_L1_7','P4_L1_14','P4_L1_21','P4_L1_28']
"Actually Datura samples"
#plant_samples = ['P1_L1_7','P1_L1_14','P1_L1_18','P1_L1_28','P3_L1_7','P3_L1_14','P3_L1_21','P3_L1_28','P5_L1_14','P5_L1_21','P5_L1_28']

"Line 2 samples (Emilia) -- exluding P4_L2_7 b/c of low coverage"
#samples = ['P0','WF_P1','P1_L2_7','P1_L2_14','P1_L2_18','P1_L2_28','WF_P2_L2','P2_L2_7','P2_L2_14','P2_L2_21','P2_L2_29','WF_P3_L2','P3_L2_7','P3_L2_14','P3_L2_21','P3_L2_28','WF_P4_L2','P4_L2_14','P4_L2_21','P4_L2_28','WF_P5_L2','P5_L2_7','P5_L2_14','P5_L2_21','P5_L2_28']
#short_names = [s.replace('_L2','') for s in samples]
#thrip_samples = ['WF_P1','WF_P2_L2','WF_P3_L2','WF_P4_L2','WF_P5_L2']
#plant_samples = ['P0','P1_L2_7','P1_L2_14','P1_L2_18','P1_L2_28','P2_L2_7','P2_L2_14','P2_L2_21','P2_L2_29','P3_L2_7','P3_L2_14','P3_L2_21','P3_L2_28','P4_L2_14','P4_L2_21','P4_L2_28','P5_L2_7','P5_L2_14','P5_L2_21','P5_L2_28']

"Line 3 samples (Datura)"
#samples = ['P0','WF_P1','P1_L3_7','P1_L3_14','P1_L3_18','P1_L3_28','WF_P2_L3','P2_L3_7','P2_L3_14','P2_L3_21','P2_L3_29','WF_P3_L3','P3_L3_7','P3_L3_14','P3_L3_21','P3_L3_28','WF_P4_L3','P4_L3_7','P4_L3_14','P4_L3_21','P4_L3_28','WF_P5_L3','P5_L3_7','P5_L3_14','P5_L3_21','P5_L3_28']
#short_names = [s.replace('_L3','') for s in samples]
#thrip_samples = ['WF_P1','WF_P2_L3','WF_P3_L3','WF_P4_L3','WF_P5_L3']
#plant_samples = ['P0','P1_L3_7','P1_L3_14','P1_L3_18','P1_L3_28','P2_L3_7','P2_L3_14','P2_L3_21','P2_L3_29','P3_L3_7','P3_L3_14','P3_L3_21','P3_L3_28','P4_L3_7','P4_L3_14','P4_L3_21','P4_L3_28','P5_L3_7','P5_L3_14','P5_L3_21','P5_L3_28']

if enrich_by_host: 
    
    "Compute avg variant freqs in thrips versus plants"
    thrip_freqs = vdf['ALT_FREQ_' + thrip_samples[0]].copy(deep=True)
    for s in range(1,len(thrip_samples)):
        thrip_freqs += vdf['ALT_FREQ_' + thrip_samples[s]]
    vdf['THRIP_FREQS'] = thrip_freqs / len(thrip_samples)
    
    plant_freqs = vdf['ALT_FREQ_' + plant_samples[0]].copy(deep=True)
    for s in range(1,len(plant_samples)):
        plant_freqs += vdf['ALT_FREQ_' + plant_samples[s]]
    vdf['PLANT_FREQS'] = plant_freqs / len(plant_samples)
    
    if (enrichment_host == 'PLANT'):
        vdf['DELTA'] = (vdf['PLANT_FREQS'] - vdf['THRIP_FREQS'])
        vdf = vdf.loc[vdf['DELTA'] >= enrichment_delta]
    elif (enrichment_host == 'THRIP'):
        vdf['DELTA'] = (vdf['THRIP_FREQS'] - vdf['PLANT_FREQS'])
        vdf = vdf.loc[vdf['DELTA'] >= enrichment_delta]
    elif (enrichment_host == 'ANY'):
        vdf['DELTA'] = (abs(vdf['THRIP_FREQS'] - vdf['PLANT_FREQS']))
        vdf = vdf.loc[vdf['DELTA'] >= enrichment_delta]
    
else:
    
    "Compute overall delta between min and max"
    col_list = ['ALT_FREQ_' + s for s in samples]   
    maxValues = vdf[col_list].max(axis=1)
    minValues = vdf[col_list].min(axis=1)
    vdf['DELTA'] = maxValues - minValues
    vdf = vdf.loc[vdf['DELTA'] >= enrichment_delta]

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
    
if gene:
    # var_type must be 'AA'
    variants = variants[(variants['CDS'] == gene)] 

"Set color palette in seaborn"
sns.set(style="darkgrid")
sns.set_palette(sns.hls_palette(len(variants), l=.4, s=.8))

fig, ax = plt.subplots(1, 1)
for index, var in variants.iterrows():
    freqs = [var['ALT_FREQ_' + s] for s in samples] # get freq for each sample in vdf
    if (var_type == 'AA'):
        #ax.plot(short_names,freqs, label = var['AA_SUBST'])
        ax.plot(short_names,freqs, label = var['CDS'] + ' ' + var['AA_SUBST'])
    else:
        ax.plot(short_names,freqs, label = var['REF'] + str(var['POS']) + var['ALT'])

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

ax.set_xlabel('Sample')
ax.set_ylabel('Frequency')
ax.grid(True)

ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.12),
          fancybox=True, shadow=False, ncol=5)
#ax.set_title('Plant-thrip enriched amino acid variants')

fig.set_size_inches(10, 4)
plt.show()

fig.savefig(out_file, dpi=200,bbox_inches='tight')