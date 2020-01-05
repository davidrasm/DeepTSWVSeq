"""
Created on Tue Aug 27 15:48:53 2019

Estimate maximum likelihood fitness effects for allelic variants based on growth rates across hosts

10.16.19:
Added method to estimate fitness for multi-allelic sites
Wrote test method to verify estimates on simulated data

11.26.19
Added trans_prob_flag to exclude paired samples where at least one variant has an impossible transition e.g. a variant at zero freq has a positve count at next time step
Added min_paired_samples arg so we only estimate fitness effects for variants with a min of at least a certain number of paired samples
Also added min_depth_coverage arg so we only estimate fitness effects for variants at sites with a min depth of coverage

@author: david
"""

import numpy as np
import pandas as pd
from scipy.stats import multinomial
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from call_aa_variants import call_aa_variants

def mn_like(x,q_t0,counts_t1,dt):
    
    sigma = np.array(x)
    like = 0.0
    for n in range(len(dt)):
        q_t1 = q_t0[n] * np.exp(sigma*dt[n])
        q_t1 = q_t1 / np.sum(q_t1) # expected allele freqs at time t1
        
        #print("sigma = ", sigma)
        #print("q_t1 = ", q_t1)
        #print("counts_t1 = ", counts_t1[n])
        
        like += np.log(multinomial.pmf(counts_t1[n], n=np.sum(counts_t1[n]), p=q_t1))
    
    return -like # return negative log likelihood because we are minimizing rather than maximizing function

def add_var_types(variants):
    
    types = []
    for index, var in variants.iterrows(): 
        if var['INDEL']:
            types.append('INDEL')
        elif var['AA_VAR']:
            types.append('AAV')
        else:
            types.append('SNV')
    variants['VAR_TYPE'] = np.array(types)       
    return variants
        
"Biallelic version only considers ALT and REF at each site"
def get_estimates_biallelic(variants,paired_samples,delta_t):
    
    "Need to modify this to handle multiple alleles at a single position"
    estimates = []
    for index, var in variants.iterrows():
        
        q_t0 = [] # allele freqs at time t0
        counts_t1 = [] # counts at time t1
        dt = [] # dt times
        for s in range(len(paired_samples)):
            
            # Freqs @ time t0
            freq_t0 = var['ALT_FREQ_' + paired_samples[s][0]]
            if (freq_t0 == 0.0) | (1.0-freq_t0==0.0): # if alt freq or ref_seq = 0
                continue
            
            # Counts @ time t1
            ref_dp_t1 = var['REF_DP_' + paired_samples[s][1]]
            alt_dp_t1 = var['ALT_DP_' + paired_samples[s][1]]
            if (ref_dp_t1 + alt_dp_t1 == 0): # if sum(counts) = 0
                continue
            
            # Elapsed time dt
            q_t0.append([1-freq_t0,freq_t0]) # [ref_freq, alt_freq]
            counts_t1.append([ref_dp_t1,alt_dp_t1])
            dt.append(delta_t[s])
        
        if q_t0: # if q_t0 is not empty
            q_t0 = np.array(q_t0)
            bnds = ((-10, 10), (-10, 10)) # bounds on sigma
            sigma_0 = (0.0, 0.0) # starting values, as tuple
            res = minimize(mn_like,sigma_0,args=(q_t0,counts_t1,dt),method='SLSQP',bounds=bnds)
            estimates.append(res.x[1])
            print(res.x)
        else:
            estimates.append(float('NaN'))
           
    return estimates

"New multi-allelic version"
def get_estimates(variants,paired_samples,delta_t):
    
    "Need to modify this to handle multiple alleles at a single position"
    estimates = []
    for index, var in variants.iterrows():
        
        "Find other variants at this site"
        region = var['REGION']
        site = var['POS']
        site_vars = variants[(variants['REGION'] == region) & (variants['POS'] == site)]
        site_vars = site_vars.reset_index(drop=True) # reset index so site_vars are indexes from 0
        var_count = len(site_vars.index) + 1
        
        "If site is multi-allelic"
        if var_count > 2:        
            "Find row index of current var"
            idx = site_vars.index[site_vars['ALT'] == var['ALT']].tolist()[0]
            "Move current var to top of site_vars"
            site_vars = pd.concat([site_vars.iloc[[idx],:], site_vars.drop(idx, axis=0)], axis=0)
        
        q_t0 = [] # allele freqs at time t0
        counts_t1 = [] # counts at time t1
        dt = [] # dt times
        for s in range(len(paired_samples)):
            
            "Order of variants in lists is ref,var,other_var,..."
            
            "Does not account for ref strain being at zero freq"
            
            # Freqs @ time t0
            freq_t0 = []
            c_t1 = []
            trans_prob_flag = False
            for idx, svar in site_vars.iterrows():
                var_freq = svar['ALT_FREQ_' + paired_samples[s][0]]
                var_count_t1 = svar['ALT_DP_' + paired_samples[s][1]]
                if (var_freq == 0.0) & (var_count_t1 > 0): # impossible transition
                    trans_prob_flag = True # raise flag so paired sample is not included
                    #print("Impossible variant freq transition")
                freq_t0.append(var_freq)
                c_t1.append(var_count_t1)   

            c_t1.insert(0,var['REF_DP_' + paired_samples[s][1]])
            freq_t0.insert(0,1-sum(freq_t0)) # [ref_freq, alt_freq]
            
            "Conditions for ignoring paired sample"
            if trans_prob_flag: # impossible transition from t0 to t1
                continue
            
            if sum(c_t1) < min_depth_coverage:
                print("Depth of coverage too low")
                continue

            if (freq_t0[0] == 0.0) | (freq_t0[0] == 1.0): # REF FREQ is 0 or 1.0
                continue
            
            if (freq_t0[1] == 0.0) | (freq_t0[1] == 1.0): # ALT FREQ is 0 or 1.0
                continue
            
            if (sum(c_t1) == 0): # no seq reads at c_t1
                continue
            
            # Added paired sample data to lists
            q_t0.append(freq_t0) 
            counts_t1.append(c_t1)
            dt.append(delta_t[s])
        
        if len(q_t0) >= min_paired_samples: # if # of paired samples is greater than min
            q_t0 = np.array(q_t0)
            bnds = tuple([(-10,10) for x in range(var_count)])
            sigma_0 = tuple([0.0 for x in range(var_count)]) # starting values, as tuple
            res = minimize(mn_like,sigma_0,args=(q_t0,counts_t1,dt),method='SLSQP',bounds=bnds)
            
            "New way estimating fitness relative to the reference"
            estimates.append(res.x[1] - res.x[0])
            
            #if np.isnan(res.x).any():
                #print(res.x)
            
            print(res.x)
        else:
            estimates.append(float('NaN'))
           
    return estimates

def test():
    
    "Simulate some data: two alleles at three time points"
    var_count = 3
    freq_t0 = [[0.33,0.33,0.34],[0.33,0.33,0.34],[0.33,0.33,0.34]] # mock initial allele freqs
    sigma = np.array([0.0,0.2,-0.1]) # true sigmas -- need to be in a nparray
    dt = [7.0,7.0,7.0]
    q_t0 = np.array(freq_t0)
    counts_t1 = []
    sample_size = 1000
    
    for n in range(len(dt)):
        q_t1 = q_t0[n] * np.exp(sigma*dt[n])
        q_t1 = q_t1 / np.sum(q_t1) # expected allele freqs at time t1
        counts_t1.append(np.round(sample_size * q_t1))
    
    "Test likelihood optimization"    
    bnds = tuple([(-10,10) for x in range(var_count)])
    sigma_0 = tuple([0.0 for x in range(var_count)]) # starting values, as tuple
    res = minimize(mn_like,sigma_0,args=(q_t0,counts_t1,dt),method='SLSQP',bounds=bnds)
    #estimates.append(res.x[1])
    print(res.x)

    
"Specify mutation/region to track"
seg = 'ANY' #'TSWV_segL'
var_type = 'ANY' # 'SNV', 'INDEL' or 'AA"

"Set minimum number of paired samples we need to estimate fitness effect of variant"
min_paired_samples = 1 # used min_paired_samples = 1 for paper
min_depth_coverage = 100 # used min_depth = 100 for paper

"Get allele freqs from masterVarList"
main_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/DeepSeqScripts/'
#tsv_file = main_path + 'tswv_L1_masterVarList_refTF2_withDepths.csv'
tsv_file = main_path + 'tswv_L1-3_masterVarList_corrAltFreqs_withDepths.csv'
vdf = pd.read_table(tsv_file, sep="\t")

out_file = "tswv_L1-3_masterVarList_corrAltFreqs_combined_minDepth100_withEstFitEffects.csv"

"Filter by variant type"
vdf['INDEL'] = vdf.ALT.str.contains('^\+|-') # regexp means first character is '+' or '-' char
if (var_type == 'SNV'):
   vdf = vdf.loc[vdf['INDEL'] == False]
elif (var_type == 'AA'):
    aa_variant,aa_subst,cds_annotation = call_aa_variants(vdf)
    vdf['AA_VAR'] = np.array(aa_variant)
    vdf['AA_SUBST'] = np.array(aa_subst)
    vdf = vdf.loc[vdf['AA_VAR'] == True]
elif (var_type == 'INDEL'):
    vdf = vdf.loc[vdf['INDEL'] == True]
elif (var_type == 'ANY'):
    aa_variant,aa_subst,cds_annotation = call_aa_variants(vdf)
    vdf['AA_VAR'] = np.array(aa_variant)
    vdf['AA_SUBST'] = np.array(aa_subst)
vdf = add_var_types(vdf) # assign variant types
    
"Filter by REGION"
if seg != 'ANY':
    variants = vdf[(vdf['REGION'] == seg)]
else:
    variants = vdf

"L1 Emilia paired samples"
emilia_paired_samples = [['WF_P2_L1','P2_L1_7'],
                        ['P2_L1_7','P2_L1_14'],
                        ['P2_L1_14','P2_L1_21'],
                        ['P2_L1_21','P2_L1_29'],
                        ['WF_P4_L1','P4_L1_7'],
                        ['P4_L1_7','P4_L1_14'],
                        ['P4_L1_14','P4_L1_21'],
                        ['P4_L1_21','P4_L1_28']]

emilia_delta_t = [7.0,7.0,7.0,8.0,7.0,7.0,7.0,7.0] # elapsed time between sample time points t_k, t_k+1  

"L1 Datura paired samples"
datura_paired_samples = [['WF_P1','P1_L1_7'],
                        ['P1_L1_7','P1_L1_14'],
                        ['P1_L1_14','P1_L1_18'],
                        ['P1_L1_18','P1_L1_28'],
                        ['WF_P3_L1','P3_L1_7'],
                        ['P3_L1_7','P3_L1_14'],
                        ['P3_L1_14','P3_L1_21'],
                        ['P3_L1_21','P3_L1_28'],
                        ['WF_P5_L1','P5_L1_14'],
                        ['P5_L1_14','P5_L1_21'],
                        ['P5_L1_21','P5_L1_28']]

datura_delta_t = [7.0,7.0,4.0,10.0,7.0,7.0,7.0,7.0,14.0,7.0,7.0] # elapsed time between sample time points t_k, t_k+1   

"Thrip paired samples"
thrip_paired_samples_L1 = [['P0','WF_P1'],
                        ['P1_L1_28','WF_P2_L1'],
                        ['P2_L1_29','WF_P3_L1'],
                        ['P3_L1_28','WF_P4_L1'],
                        ['P4_L1_28','WF_P5_L1']]

thrip_delta_t_L1 = [7.0,7.0,7.0,7.0,7.0]

"Can combine Emilia and Datura samples together"
plant_paired_samples_L1 = emilia_paired_samples + datura_paired_samples
plant_delta_t_L1 = emilia_delta_t + datura_delta_t 

"Line 2 samples"

"L2 plant samples"
plant_paired_samples_L2 = [['WF_P1','P1_L2_7'],
                        ['P1_L2_7','P1_L2_14'],
                        ['P1_L2_14','P1_L2_18'],
                        ['P1_L2_18','P1_L2_28'],
                        ['WF_P2_L2','P2_L2_7'],
                        ['P2_L2_7','P2_L2_14'],
                        ['P2_L2_14','P2_L2_21'],
                        ['P2_L2_21','P2_L2_29'],
                        ['WF_P3_L2','P2_L2_7'],
                        ['P3_L2_7','P3_L2_14'],
                        ['P3_L2_14','P3_L2_21'],
                        ['P3_L2_21','P3_L2_28'],
                        ['WF_P4_L2','P4_L2_14'],
                        ['P4_L2_14','P4_L2_21'],
                        ['P4_L2_21','P4_L2_28'],
                        ['WF_P5_L2','P5_L2_7'],
                        ['P5_L2_7','P5_L2_14'],
                        ['P5_L2_14','P5_L2_21'],
                        ['P5_L2_21','P5_L2_28']]

plant_delta_t_L2 = [7.0,7.0,4.0,10.0,7.0,7.0,7.0,8.0,7.0,7.0,7.0,7.0,14.0,7.0,7.0,7.0,7.0,7.0,7.0] # elapsed time between sample time points t_k, t_k+1

"L2 thrip paired samples"
thrip_paired_samples_L2 = [['P0','WF_P1'],
                        ['P1_L2_28','WF_P2_L2'],
                        ['P2_L2_29','WF_P3_L2'],
                        ['P3_L2_28','WF_P4_L2'],
                        ['P4_L2_28','WF_P5_L2']]

thrip_delta_t_L2 = [7.0,7.0,7.0,7.0,7.0]

"Line 3 samples"

"L3 plant samples"
plant_paired_samples_L3 = [['WF_P1','P1_L3_7'],
                        ['P1_L3_7','P1_L3_14'],
                        ['P1_L3_14','P1_L3_18'],
                        ['P1_L3_18','P1_L3_28'],
                        ['WF_P2_L3','P2_L3_7'],
                        ['P2_L3_7','P2_L3_14'],
                        ['P2_L3_14','P2_L3_21'],
                        ['P2_L3_21','P2_L3_29'],
                        ['WF_P3_L3','P3_L3_7'],
                        ['P3_L3_7','P3_L3_14'],
                        ['P3_L3_14','P3_L3_21'],
                        ['P3_L3_21','P3_L3_28'],
                        ['WF_P4_L3','P4_L3_7'],
                        ['P4_L3_7','P4_L3_14'],
                        ['P4_L3_14','P4_L3_21'],
                        ['P4_L3_21','P4_L3_28'],
                        ['WF_P5_L3','P5_L3_7'],
                        ['P5_L3_7','P5_L3_14'],
                        ['P5_L3_14','P5_L3_21'],
                        ['P5_L3_21','P5_L3_28']]

plant_delta_t_L3 = [7.0,7.0,4.0,10.0,7.0,7.0,7.0,8.0,7.0,7.0,7.0,7.0,7.0,7.0,7.0,7.0,7.0,7.0,7.0,7.0] # elapsed time between sample time points t_k, t_k+1

"L3 thrip paired samples"
thrip_paired_samples_L3 = [['P0','WF_P1'],
                        ['P1_L3_28','WF_P2_L3'],
                        ['P2_L3_29','WF_P3_L3'],
                        ['P3_L3_28','WF_P4_L3'],
                        ['P4_L3_28','WF_P5_L3']]

thrip_delta_t_L3 = [7.0,7.0,7.0,7.0,7.0]

"Combine samples across lines"
emilia_paired_samples = emilia_paired_samples + plant_paired_samples_L2
emilia_delta_t = emilia_delta_t + plant_delta_t_L2

datura_paired_samples = datura_paired_samples + plant_paired_samples_L3
datura_delta_t = datura_delta_t + plant_delta_t_L3

plant_paired_samples = plant_paired_samples_L1 + plant_paired_samples_L2 + plant_paired_samples_L3
plant_delta_t = plant_delta_t_L1 + plant_delta_t_L2 + plant_delta_t_L3

thrip_paired_samples = thrip_paired_samples_L1 + thrip_paired_samples_L2[1:] + thrip_paired_samples_L3[1:]
thrip_delta_t = thrip_delta_t_L1 + thrip_delta_t_L2[1:] + thrip_delta_t_L3[1:]


"Paried samples for Emilia"
emilia_fit_ests = get_estimates(variants, emilia_paired_samples, emilia_delta_t)
variants['emilia_FE'] = np.array(emilia_fit_ests)

"Paired sample time points t_k, t_k+1 for L1 in Datura"
datura_fit_ests = get_estimates(variants, datura_paired_samples, datura_delta_t)
variants['datura_FE'] = np.array(datura_fit_ests)

"Paried samples for plants->thrip transfers"
thrip_fit_ests = get_estimates(variants, thrip_paired_samples, thrip_delta_t)
variants['thrip_FE'] = np.array(thrip_fit_ests)

"Paired samples for all plants"
plant_fit_ests = get_estimates(variants, plant_paired_samples, plant_delta_t)
variants['plant_FE'] = np.array(plant_fit_ests)

"Save file of variants"
variants.to_csv(out_file, sep='\t')

fig, ax = plt.subplots()
ax.scatter(variants['plant_FE'],variants['thrip_FE'])
ax.set_xlabel('Fitness in plants', fontsize=14)
ax.set_ylabel('Fitness in thrips', fontsize=14)
lim = 1.0
ax.set_xlim((-lim, lim))
ax.set_ylim((-lim, lim))

fig.set_size_inches(6, 6)
plt.show()

out_file = 'L1_plantThrip_jointDFE_AAVs_L1.png'
fig.savefig(out_file, dpi=200,bbox_inches='tight')

"Need to get this for each time point pair"    
#q_t0 = np.array([0.9,0.1]) # allele freqs at time t0
#counts_t1 = [9,1]
#dt = 1.0

"Test for two time point pairs"
#q_t0 = np.array([[0.5,0.5],[0.5,0.5]]) # allele freqs at time t0
#counts_t1 = [[9,1],[8,2]]
#dt = [1.0,1.0]

#sigma = sigma_0
#like = mn_like(sigma,q_t0,counts_t1,dt) # just for testing
    

    


