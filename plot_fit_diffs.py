"""
Created on Tue Sep 24 14:59:52 2019

Plot fitness differences between two hosts mapped onto TSWV genome

Used to plot fitness differences in Figure 9 in paper

@author: david
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.patches as patches


def label_point(x, y, val, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    for i, point in a.iterrows():
        ax.text(point['x'], point['y'], str(point['val']))
        
def add_gene(ax, label, start, end):
    
    "Add gene or feature annotation to plot ax from start to end"
    left = start
    bottom = -1.5 * ax.get_ylim()[1]
    width = end - start
    height = 0.2 * ax.get_ylim()[1]
    p = patches.Rectangle(
        (left, bottom), width, height,
        fill=True, color='gray', transform=ax.transData, clip_on=False
        )
    ax.add_patch(p)
    ax.text(left+(width/2), bottom+(height/2), label,
        horizontalalignment='center',
        verticalalignment='center',
        transform=ax.transData)

fig_name = "L1-3_emiliaDatura_fitnessDiffs_acrossGenome_corrAltFreqs_combined_GnGc.png"
var_type = 'AA'

main_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/DeepSeqResults/dfes/'
#var_file = "tswv_L1_masterVarList_corrAltFreqs_withEstFitEffects.csv"
var_file = main_path + "tswv_L1-3_masterVarList_corrAltFreqs_combined_minDepth100_withEstFitEffects.csv"
variants = pd.read_table(var_file, sep="\t")

if (var_type == 'SNV'):
   variants = variants.loc[variants['INDEL'] == False]
elif (var_type == 'AA'):
    variants = variants.loc[variants['AA_VAR'] == True]
elif (var_type == 'INDEL'):
    variants = variants.loc[variants['INDEL'] == True]
elif (var_type == 'ANY'):
    variants = variants
variants['DELTA_EMILIA_DATURA'] = variants['emilia_FE'] - variants['datura_FE']
variants['DELTA_PLANTS_THRIPS'] = variants['plant_FE'] - variants['thrip_FE']
#variants = variants[variants['DELTA_PLANTS_THRIPS'].notnull()]

df_segS = variants.loc[variants['REGION'] == "TSWV_segS"] # get rows for segment
df_segM = variants.loc[variants['REGION'] == "TSWV_segM"] # get rows for segment
df_segL = variants.loc[variants['REGION'] == "TSWV_segL"] # get rows for segment

sns.set()     
fig, axs = plt.subplots(3, 1)
    
"Plot fitness differences 'delta'"
#df_segS.plot(x="POS", y='DELTA', ax=axs[0], kind="scatter", c=df_segS['DELTA'], cmap=cm.inferno, alpha = 0.5)
df_segS.plot(x="POS", y='DELTA_EMILIA_DATURA', ax=axs[0], kind="scatter", color="green", alpha = 0.3,  label = 'Emilia - Datura')
df_segS.plot(x="POS", y='DELTA_PLANTS_THRIPS', ax=axs[0], kind="scatter", color="purple", alpha = 0.3, label = 'Plants - Thrips')

axs[0].legend()
axs[0].legend(loc='upper right') # Move the legend to an empty part of the plot
axs[0].set_xlabel('')
axs[0].set_ylabel('Fitness difference')
#axs[0].set_yticks([0.])
#axs[0].yaxis.grid(True)
axs[0].grid(True)
axs[0].set_title('TSWV segment S')
axs[0].set_xlim([0,2915])
ymin = -0.5
ymax = 0.5
axs[0].set_ylim([ymin,ymax])
#if (var_type == 'AA'):
    #label_point(df_segS.POS, df_segS.DELTA, df_segS.AA_SUBST, axs[0])
#else:
    #label_point(df_segS.POS, df_segS.DELTA, df_segS.ALT, axs[0])
add_gene(axs[0], 'N', 89, 1483) # Add  N 89..1483
add_gene(axs[0], 'NSs', 1987, 2763) # Add NSs 1987..2763

"Add zero zero line"
x = np.linspace(0, 2915, 2)
axs[0].plot(x, np.zeros(len(x)), 'k', alpha=.25) # dashdot black
    
"Plot fitness differences 'delta'"
#df_segM.plot(x="POS", y='DELTA', ax=axs[1], kind="scatter", c=df_segM['DELTA'], cmap=cm.inferno, alpha = 0.5)
df_segM.plot(x="POS", y='DELTA_EMILIA_DATURA', ax=axs[1], kind="scatter", color="green", alpha = 0.3)
df_segM.plot(x="POS", y='DELTA_PLANTS_THRIPS', ax=axs[1], kind="scatter", color="purple", alpha = 0.3)

axs[1].set_xlabel('')
axs[1].set_ylabel('Fitness difference')
axs[1].grid(True)
axs[1].set_title('TSWV segment M')
axs[1].set_xlim([0,4820])
axs[1].set_ylim([ymin,ymax])
#if (var_type == 'AA'):
    #label_point(df_segM.POS, df_segM.DELTA_EMILIA_DATURA, df_segM.AA_SUBST, axs[1])
#else:
    #label_point(df_segM.POS, df_segM.DELTA, df_segM.ALT, axs[1])
add_gene(axs[1], 'NSm', 101, 1009) # Add  NSm 101..1009
#add_gene(axs[1], 'GP', 1330, 4737) # Add GP 1330..4737
"If dividing GP into Gn/Gc"
add_gene(axs[1], r'$G_C$', 1330, 3285)
add_gene(axs[1], r'$G_N$', 3286, 4737)


"Add zero zero line"
x = np.linspace(0, 4820, 2)
axs[1].plot(x, np.zeros(len(x)), 'k', alpha=.25) # dashdot black
    
"Plot fitness differences 'delta'"
#df_segL.plot(x="POS", y='DELTA', ax=axs[2], kind="scatter", c=df_segL['DELTA'], cmap=cm.inferno, alpha = 0.5)
df_segL.plot(x="POS", y='DELTA_EMILIA_DATURA', ax=axs[2], kind="scatter", color="green", alpha = 0.3)
df_segL.plot(x="POS", y='DELTA_PLANTS_THRIPS', ax=axs[2], kind="scatter", color="purple", alpha = 0.3)

axs[2].set_xlabel('')
axs[2].set_ylabel('Fitness difference')
axs[2].grid(True)
axs[2].set_title('TSWV segment L')
axs[2].set_xlim([0,8896])
axs[2].set_ylim([ymin,ymax])
#if (var_type == 'AA'):
    #label_point(df_segL.POS, df_segL.DELTA, df_segL.AA_SUBST, axs[2])
#else:
    #label_point(df_segL.POS, df_segL.DELTA, df_segL.ALT, axs[2])
add_gene(axs[2], 'RdRp', 34, 8661) # Add RdRp 34..8661

"Add zero zero line"
x = np.linspace(0, 8896, 2)
axs[2].plot(x, np.zeros(len(x)), 'k', alpha=.25) # dashdot black

fig.set_size_inches(8, 8)
fig.tight_layout()
plt.show()

fig.savefig(fig_name, dpi=200)


