"""
Created on Fri Aug 30 10:17:06 2019

Plot joint distribution of fitness effects across host/vectors as scatter plots

Used to create joint DFEs in Figure 8 in paper

@author: david
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

var_file = "tswv_L1-3_masterVarList_corrAltFreqs_combined_minDepth100_withEstFitEffects.csv" # master list including fitness effects used for Figure 8 in paper
variants = pd.read_table(var_file, sep="\t")

sns.set()
#sns.set_context("talk")
fig, axs = plt.subplots(1, 2)
#scatter_kws={"s": 20}

ax = plt.subplot(1,2,1)
sns.scatterplot(x="emilia_FE", y="datura_FE", data=variants, hue="VAR_TYPE", s=20, alpha=.5)
ax.set_xlabel('Fitness in Emilia', fontsize=14)
ax.set_ylabel('Fitness in Datura', fontsize=14)
xmin = -0.5
lim = 1.0
ax.set_xlim((xmin, lim))
ax.set_ylim((xmin, lim))

"Add zero zero line"
x = np.linspace(xmin, lim, 2)
ax.plot(x, np.zeros(len(x)), 'k', alpha=.25) # dashdot black
ax.plot(np.zeros(len(x)), x, 'k', alpha=.25) # dashdot black

ax = plt.subplot(1,2,2)
sns.scatterplot(x="plant_FE", y="thrip_FE", data=variants, hue="VAR_TYPE", s=20, alpha=.5) #,size="col_name_3")
ax.set_xlabel('Fitness in Plants', fontsize=14)
ax.set_ylabel('Fitness in Thrips', fontsize=14)
ax.set_xlim((xmin, lim))
ax.set_ylim((xmin, lim))

"Add zero zero line"
x = np.linspace(xmin, lim, 2)
ax.plot(x, np.zeros(len(x)), 'k', alpha=.25) # dashdot black
ax.plot(np.zeros(len(x)), x, 'k', alpha=.25) # dashdot black

# Move the legend to an empty part of the plot
#ax.legend(loc='lower left')
 
fig.set_size_inches(10, 4)
plt.show()

out_file = 'L1-3_emiliaDaturaPlantThrip_jointDFE_corrAltFreqs_combined_minDepth100.png'
fig.savefig(out_file, dpi=200,bbox_inches='tight')
