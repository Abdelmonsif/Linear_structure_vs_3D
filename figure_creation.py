
import os,sys,argparse
import pandas as pd
import time
from os import listdir
import numpy as np
import matplotlib.pyplot as plt




def figure_creation(one, three, title, label, sora):

    file1 = pd.read_csv(one,delim_whitespace=True)
    file2 = pd.read_csv(three,delim_whitespace=True)
    dataset_interactions = file1.Median_negative_log_enhancer_score,file2.Median_negative_log_enhancer_score
    fig, ax = plt.subplots(figsize =(12, 14))
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rc('font', family='serif', serif='Times')
    sns.violinplot(data=dataset_interactions, palette=['#77F5BE','#F7A4A6'], ax=ax)
    sns.boxplot(data=dataset_interactions, width=0.33)
    # add x-tick labels
    plt.xticks([0,1],['1D', '3D'], fontsize=24)
    ax.set_ylim(1.5, 6.3)
    ax.set_yticks([2, 3, 4, 5, 6])
    plt.yticks(fontsize=22)
    plt.rcParams.update({'font.size': 20})
    ax.set_title(title,fontsize=30)
    ax.set_ylabel(label, fontsize=28)
    plt.savefig(sora,dpi=600)
    plt.show()
    
    
 def getArgs():
    parser = argparse.ArgumentParser('python')
    parser.add_argument('-one', required=True)
    parser.add_argument('-three', required=True)
    parser.add_argument('-title', required=True)
    parser.add_argument('-label', required=True)
    parser.add_argument('-sora', required=True)
    return parser.parse_args()

if __name__ == "__main__":
    args = getArgs()
    start = time.time()
    filename = nt_comp(args.inputFile, args.outputFile, args.res)
    end = time.time()
    print ('time elapsed:' + str(end - start))  



