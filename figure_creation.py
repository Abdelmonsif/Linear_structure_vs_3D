
import os,sys,argparse
import pandas as pd
import time
from os import listdir
import numpy as np
import matplotlib.pyplot as plt




def figure_creation():

    file1 = pd.read_csv('enh_breast_1D',delim_whitespace=True)
    file2 = pd.read_csv('enh_breast_3D',delim_whitespace=True)
    dataset_interactions = file1.Median_negative_log_enhancer_score,file2.Median_negative_log_enhancer_score
    fig, ax = plt.subplots(figsize =(12, 14))
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rc('font', family='serif', serif='Times')
    sns.violinplot(data=dataset_interactions, palette=['#77F5BE','#F7A4A6'], ax=ax)
    sns.boxplot(data=dataset_interactions, width=0.2)
    # add x-tick labels
    plt.xticks([0,1],['1D', '3D'], fontsize=24)
    ax.set_ylim(1.5, 6.3)
    ax.set_yticks([2, 3, 4, 5, 6])
    plt.yticks(fontsize=22)
    plt.rcParams.update({'font.size': 20})
    ax.set_title('5kb resolution breast cancer\n',fontsize=30)
    ax.set_ylabel('Median DA score for enhancers', fontsize=28)
    plt.savefig('figure2A',dpi=600)
    plt.show()
    
    
    
    file1 = pd.read_csv('enh_prostate_1D',delim_whitespace=True)
    file2 = pd.read_csv('enh_prostate_3D',delim_whitespace=True)
    dataset_interactions = file1.Median_negative_log_enhancer_score,file2.Median_negative_log_enhancer_score
    fig, ax = plt.subplots(figsize =(12, 14))
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rc('font', family='serif', serif='Times')
    sns.violinplot(data=dataset_interactions, palette=['#77F5BE','#F7A4A6'], ax=ax)
    sns.boxplot(data=dataset_interactions, width=0.2)
    # add x-tick labels
    plt.xticks([0,1],['1D', '3D'], fontsize=24)
    ax.set_ylim(1, 6.3)
    ax.set_yticks([1, 2, 3, 4, 5, 6])
    plt.yticks(fontsize=22)
    plt.rcParams.update({'font.size': 20})
    ax.set_title('10kb resolution prostate cancer\n',fontsize=30)
    ax.set_ylabel('Median DA score for enhancers', fontsize=28)
    plt.savefig('figure2B',dpi=600)
    plt.show()
    
    
    
    file1 = pd.read_csv('enh_tg_breast_1D',delim_whitespace=True)
    file2 = pd.read_csv('enh_tg_breast_3D',delim_whitespace=True)
    dataset_interactions = file1.Median_disease_score,file2.Median_disease_score
    fig, ax = plt.subplots(figsize =(12, 14))
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rc('font', family='serif', serif='Times')
    sns.violinplot(data=dataset_interactions, palette=['#77F5BE','#F7A4A6'], ax=ax)
    sns.boxplot(data=dataset_interactions, width=0.2)
    # add x-tick labels
    plt.xticks([0,1],['1D', '3D'], fontsize=24)
    ax.set_ylim(1, 6.3)
    ax.set_yticks([1, 2, 3, 4, 5, 6])
    plt.yticks(fontsize=22)
    plt.rcParams.update({'font.size': 20})
    ax.set_title('5kb resolution breast cancer\n',fontsize=30)
    ax.set_ylabel('Median DA score for enhancer target genes', fontsize=28)
    plt.savefig('figure3A',dpi=600)
    plt.show()
    
    
    file1 = pd.read_csv('enh_tg_prostate_1D',delim_whitespace=True)
    file2 = pd.read_csv('enh_tg_prostate_3D',delim_whitespace=True)
    dataset_interactions = file1.Median_disease_score,file2.Median_disease_score
    fig, ax = plt.subplots(figsize =(12, 14))
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rc('font', family='serif', serif='Times')
    sns.violinplot(data=dataset_interactions, palette=['#77F5BE','#F7A4A6'], ax=ax)
    sns.boxplot(data=dataset_interactions, width=0.2)
    # add x-tick labels
    plt.xticks([0,1],['1D', '3D'], fontsize=24)
    ax.set_ylim(1, 8.4)
    ax.set_yticks([0, 2, 4, 6, 8])
    plt.yticks(fontsize=22)
    plt.rcParams.update({'font.size': 20})
    ax.set_title('10kb resolution prostate cancer\n',fontsize=30)
    ax.set_ylabel('Median DA score for enhancer target genes', fontsize=28)
    plt.savefig('figure3B',dpi=600)
    plt.show()
    
    
    file1 = pd.read_csv('TF_breast_1D',delim_whitespace=True)
    file2 = pd.read_csv('TF_breast_3D',delim_whitespace=True)
    dataset_interactions = file1.Median_TF_score,file2.Median_TF_score
    fig, ax = plt.subplots(figsize =(12, 14))
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rc('font', family='serif', serif='Times')
    sns.violinplot(data=dataset_interactions, palette=['#77F5BE','#F7A4A6'], ax=ax)
    sns.boxplot(data=dataset_interactions, width=0.2)
    # add x-tick labels
    plt.xticks([0,1],['1D', '3D'], fontsize=24)
    ax.set_ylim(0, 23)
    ax.set_yticks([0, 5, 10, 15, 20])
    plt.yticks(fontsize=22)
    plt.rcParams.update({'font.size': 20})
    ax.set_title('5kb resolution breast cancer\n',fontsize=30)
    ax.set_ylabel('Median DA score for TFs', fontsize=28)
    plt.savefig('figure4A',dpi=600)
    plt.show()
    
    
    file1 = pd.read_csv('TF_prostate_1D',delim_whitespace=True)
    file2 = pd.read_csv('TF_prostate_3D',delim_whitespace=True)
    dataset_interactions = file1.TF_frequency,file2.TF_frequency
    fig, ax = plt.subplots(figsize =(12, 14))
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rc('font', family='serif', serif='Times')
    sns.violinplot(data=dataset_interactions, palette=['#77F5BE','#F7A4A6'], ax=ax)
    sns.boxplot(data=dataset_interactions, width=0.2)
    # add x-tick labels
    plt.xticks([0,1],['1D', '3D'], fontsize=24)
    ax.set_ylim(0.2, 0.8)
    ax.set_yticks([0.2, 0.4, 0.6, 0.8])
    plt.yticks(fontsize=22)
    plt.rcParams.update({'font.size': 20})
    ax.set_title('10kb resolution prostate cancer\n',fontsize=30)
    ax.set_ylabel('Median DA score for TFs', fontsize=28)
    plt.savefig('figure4B',dpi=600)
    plt.show()
    
    file1 = pd.read_csv('TF_tg_breast_1D',delim_whitespace=True)
    file2 = pd.read_csv('TF_tg_breast_3D',delim_whitespace=True)
    dataset_interactions = file1.Median_score_disease,file2.Median_score_disease
    fig, ax = plt.subplots(figsize =(12, 14))
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rc('font', family='serif', serif='Times')
    sns.violinplot(data=dataset_interactions, palette=['#77F5BE','#F7A4A6'], ax=ax)
    sns.boxplot(data=dataset_interactions, width=0.2)
    # add x-tick labels
    plt.xticks([0,1],['1D', '3D'], fontsize=24)
    ax.set_ylim(1, 6.6)
    ax.set_yticks([1, 2, 3, 4, 5, 6])
    plt.yticks(fontsize=22)
    plt.rcParams.update({'font.size': 20})
    ax.set_title('5kb resolution breast cancer\n',fontsize=30)
    ax.set_ylabel('Median DA score for TF target genes', fontsize=28)
    plt.savefig('figure5A',dpi=600)
    plt.show()    

    file1 = pd.read_csv('TF_tg_prostate_1D',delim_whitespace=True)
    file2 = pd.read_csv('TF_tg_prostate_3D',delim_whitespace=True)
    dataset_interactions = file1.Median_score_disease,file2.Median_score_disease
    fig, ax = plt.subplots(figsize =(12, 14))
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rc('font', family='serif', serif='Times')
    sns.violinplot(data=dataset_interactions, palette=['#77F5BE','#F7A4A6'], ax=ax)
    sns.boxplot(data=dataset_interactions, width=0.2)
    # add x-tick labels
    plt.xticks([0,1],['1D', '3D'], fontsize=24)
    ax.set_ylim(0, 7)
    ax.set_yticks([0, 2, 4, 6])
    plt.yticks(fontsize=22)
    plt.rcParams.update({'font.size': 20})
    ax.set_title('10kb resolution prostate cancer\n',fontsize=30)
    ax.set_ylabel('Median DA score for TF target genes', fontsize=28)
    plt.savefig('figure5B',dpi=600)
    plt.show()
    
    
    
    file1 = pd.read_csv('enh_breast_noSNPs_TADs',delim_whitespace=True)
    file2 = pd.read_csv('enh_breast_SNPs_TADs',delim_whitespace=True)
    dataset_interactions = file1.Median_negative_log_score,file2.Median_negative_log_score
    fig, ax = plt.subplots(figsize =(12, 14))
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rc('font', family='serif', serif='Times')
    sns.violinplot(data=dataset_interactions, palette=['#77F5BE','#F7A4A6'], ax=ax)
    sns.boxplot(data=dataset_interactions, width=0.2)
    # add x-tick labels
    plt.xticks([0,1],['Control', 'SNP-rich'], fontsize=24)
    ax.set_ylim(1, 6.6)
    ax.set_yticks([2, 3, 4, 5, 6])
    plt.yticks(fontsize=22)
    plt.rcParams.update({'font.size': 20})
    ax.set_title('5kb resolution breast cancer\n',fontsize=30)
    ax.set_ylabel('Median DA score for enhancers', fontsize=28)
    plt.savefig('figure9A',dpi=600)
    plt.show()
    

    file1 = pd.read_csv('enh_prostate_noSNPs_TADs',delim_whitespace=True)
    file2 = pd.read_csv('enh_prostate_SNPs_TADs',delim_whitespace=True)
    dataset_interactions = file1.Median_negative_log_score,file2.Median_negative_log_score
    fig, ax = plt.subplots(figsize =(12, 14))
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rc('font', family='serif', serif='Times')
    sns.violinplot(data=dataset_interactions, palette=['#77F5BE','#F7A4A6'], ax=ax)
    sns.boxplot(data=dataset_interactions, width=0.2)
    plt.xticks([0,1],['Control', 'SNP-rich'], fontsize=24)
    ax.set_ylim(1.5, 6.2)
    ax.set_yticks([2, 3, 4, 5, 6])
    plt.yticks(fontsize=22)
    plt.rcParams.update({'font.size': 20})
    ax.set_title('10kb resolution prostate cancer\n',fontsize=30)
    ax.set_ylabel('Median DA score for enhancers', fontsize=28)
    plt.savefig('figure9B',dpi=600)
    plt.show()
    
    
    file1 = pd.read_csv('TF_breast_noSNPs_TADs',delim_whitespace=True)
    file2 = pd.read_csv('TF_breast_SNPs_TADs',delim_whitespace=True)
    dataset_interactions = file1.median_tf_score,file2.median_tf_score
    fig, ax = plt.subplots(figsize =(12, 14))
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rc('font', family='serif', serif='Times')
    sns.violinplot(data=dataset_interactions, palette=['#77F5BE','#F7A4A6'], ax=ax)
    sns.boxplot(data=dataset_interactions, width=0.002)
    plt.xticks([0,1],['Control', 'SNP-rich'], fontsize=24)
    ax.set_ylim(8.2, 21)
    ax.set_yticks([10, 12, 14, 16, 18, 20])
    plt.yticks(fontsize=22)
    plt.rcParams.update({'font.size': 20})
    ax.set_title('5kb resolution breast cancer\n',fontsize=30)
    ax.set_ylabel('Median DA score for TFs', fontsize=28)
    plt.savefig('figure10A',dpi=600)
    plt.show()    
    
    
    file1 = pd.read_csv('TF_prostate_noSNPs_TADs',delim_whitespace=True)
    file2 = pd.read_csv('TF_prostate_SNPs_TADs',delim_whitespace=True)
    dataset_interactions = file1.TF_frequency,file2.TF_frequency
    fig, ax = plt.subplots(figsize =(12, 14))
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rc('font', family='serif', serif='Times')
    sns.violinplot(data=dataset_interactions, palette=['#77F5BE','#F7A4A6'], ax=ax)
    sns.boxplot(data=dataset_interactions, width=0.2)
    # add x-tick labels
    plt.xticks([0,1],['Control', 'SNP-rich'], fontsize=24)
    ax.set_ylim(0.2, 0.75)
    ax.set_yticks([0.3, 0.4, 0.5, 0.6, 0.7])
    plt.yticks(fontsize=22)
    plt.rcParams.update({'font.size': 20})
    ax.set_title('10kb resolution prostate cancer\n',fontsize=30)
    ax.set_ylabel('Median DA score for TFs', fontsize=28)
    plt.savefig('figure10B',dpi=600)
    plt.show()   

    

if __name__ == "__main__":
    start = time.time()
    filename = figure_creation()
    end = time.time()
    print ('time elapsed:' + str(end - start))  



