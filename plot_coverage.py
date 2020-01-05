"""
Created on Mon May 20 17:32:55 2019

Plot depth of seq coverage across TSWV genome

@author: david
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
import os

def add_gene(ax, label, start, end):
    
    "Add gene or feature annotation to plot ax from start to end"
    left = start
    bottom = -0.3 * ax.get_ylim()[1]
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

def plot_coverage(depth_log_rep1,depth_log_rep2,out_file):
    
    #sns.set()
    #sns.set_context("talk")

    #file_path = "../samples/WF-P1/"
    #depth_log_rep1 = file_path + "WF-P1-r1_depth.log"
    #depth_log_rep2 = file_path + "WF-P1-r2_depth.log"
    
    df_rep1 = pd.read_table(depth_log_rep1, sep = "\t", names = ["Ref", "Pos", "Depth"])
    df_rep2 = pd.read_table(depth_log_rep2, sep = "\t", names = ["Ref", "Pos", "Depth"])
    
    df_rep1_segS = df_rep1.loc[df_rep1['Ref'] == "TSWV_segS"]; # get rows for segment S
    df_rep1_segM = df_rep1.loc[df_rep1['Ref'] == "TSWV_segM"]; # get rows for segment M
    df_rep1_segL = df_rep1.loc[df_rep1['Ref'] == "TSWV_segL"]; # get rows for segment L
    
    df_rep2_segS = df_rep2.loc[df_rep2['Ref'] == "TSWV_segS"]; # get rows for segment S
    df_rep2_segM = df_rep2.loc[df_rep2['Ref'] == "TSWV_segM"]; # get rows for segment M
    df_rep2_segL = df_rep2.loc[df_rep2['Ref'] == "TSWV_segL"]; # get rows for segment L
    
    fig, axs = plt.subplots(3, 1)
    axs[0].plot(df_rep1_segS["Depth"], label = "RT Rep 1")
    axs[0].plot(df_rep2_segS["Depth"], label = "RT Rep 2")
    axs[0].set_xlabel('')
    axs[0].set_ylabel('Depth of coverage')
    axs[0].grid(True)
    axs[0].set_xlim([0,2915])
    axs[0].legend()
    axs[0].set_title('TSWV segment S')
    add_gene(axs[0], 'N', 89, 1483) # Add  N 89..1483
    add_gene(axs[0], 'NSs', 1987, 2763) # Add NSs 1987..2763
    
    sites = list(range(df_rep1_segM["Depth"].size))
    axs[1].plot(sites,df_rep1_segM["Depth"], label = "RT Rep 1")
    sites = list(range(df_rep2_segM["Depth"].size))
    axs[1].plot(sites,df_rep2_segM["Depth"], label = "RT Rep 2")
    axs[1].set_xlabel('')
    axs[1].set_ylabel('Depth of coverage')
    axs[1].grid(True)
    axs[1].set_xlim([0,4820])
    #axs[0].legend()
    axs[1].set_title('TSWV segment M')
    add_gene(axs[1], 'NSm', 101, 1009) # Add  NSm 101..1009
    add_gene(axs[1], 'GP', 1330, 4737) # Add GP 1330..4737
    
    sites = list(range(df_rep1_segL["Depth"].size))
    axs[2].plot(sites,df_rep1_segL["Depth"], label = "RT Rep 1")
    sites = list(range(df_rep2_segL["Depth"].size))
    axs[2].plot(sites,df_rep2_segL["Depth"], label = "RT Rep 2")
    axs[2].set_xlabel('Position')
    axs[2].set_ylabel('Depth of coverage')
    axs[2].grid(True)
    axs[2].set_xlim([0,8896])
    #axs[0].legend()
    axs[2].set_title('TSWV segment L')
    add_gene(axs[2], 'RdRp', 34, 8661) # Add RdRp 34..8661
    
    fig.set_size_inches(8, 6)
    fig.tight_layout()
    plt.show()
    
    fig.savefig(out_file, dpi=200)

if __name__ == '__main__':
    
    batch = False # processing multiple samples?
    
    "For multiple fastq files in a directory"
    if batch:
        
        main_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/samples/'
        
        "Get all directories in main_path"
        dir_list = os.listdir(main_path)
        
        for dir in dir_list:
            
            print(dir)
            
            if (dir[0] == '.'):
                continue
                
            os.chdir(main_path + dir) # change path to sample_path
            depth_log_rep1 = dir + '_r1_depth.log'
            depth_log_rep2 = dir + '_r2_depth.log'
            out_file = dir + '_coverage.png'
            
            plot_coverage(depth_log_rep1,depth_log_rep2,out_file)
            
    else:
        
        main_path = '/Volumes/GoogleDrive/My Drive/DeepTSWV/samples/'
        dir = 'P2_L2_14'
        
        os.chdir(main_path + dir) # change path to sample_path
        depth_log_rep1 = dir + '_r1_self-concensus_depth.log'
        depth_log_rep2 = dir + '_r2_self-concensus_depth.log'
        out_file = dir + '_self-concensus_coverage.png'
        
        plot_coverage(depth_log_rep1,depth_log_rep2,out_file)
