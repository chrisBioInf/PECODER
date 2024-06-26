#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 12:45:33 2024

@author: christopher
"""


import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


arch_columns = ['a', 'b', 'a+b', 'a/b']
order = ['a', 'b', 'a+b', 'a/b', 'None']
xlabels = ['alpha interactions', 'beta interactions', 'alpha+beta interactions', 'alpha/beta interactions']


def plot_architecture_comparison(name: str, domain_table: pd.DataFrame):
    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    axes = [x for xs in axes for x in xs]
    
    domain_table = domain_table[domain_table['anchor_domain_arch'] != 'None']
    
    for i in range(0, len(arch_columns)):
        arch = arch_columns[i]
        ax = axes[i]
        
        df_ = domain_table[(domain_table['domain_1_arch'] == arch) & (domain_table['domain_2_arch'].isin(order))]
        
        data = {
            'Domain 2 architecture': 2*list(df_['domain_2_arch']),
            'Residues': list(df_['domain_1']) + list(df_['domain_2']) ,
            'Architecture': len(df_)*[arch] + len(df_)*['Other']
            }
        sns.violinplot(data=data, x='Domain 2 architecture', y='Residues', ax=ax, 
                       split=True, edgecolor='k', hue='Architecture',
                       gap=.15, inner="quart", order=order)
        ax.set_xlabel('Domain 2 architecture')
        ax.set_ylabel('Residues')
    
    if len(name) > 0:
        plt.savefig(name.split('.')[0] + '_alpha_beta.svg', dpi=300)

    plt.show()
    
    