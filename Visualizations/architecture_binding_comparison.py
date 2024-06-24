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
xlabels = ['alpha interactions', 'beta interactions', 'alpha+beta interactions', 'alpha/beta interactions']


def plot_architecture_comparison(name: str, architecture_table: pd.DataFrame):
    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    axes = [x for xs in axes for x in xs]
    
    architecture_table.to_csv('architecture_table.tsv', sep='\t')
    
    for i in range(0, len(arch_columns)):
        arch = arch_columns[i]
        ax = axes[i]
        key = '\{%s}' % arch
        architecture_table[key] = architecture_table['a'] + architecture_table['b'] + architecture_table['a+b'] + architecture_table['a/b'] + architecture_table['other'] - architecture_table[arch]
    
        df_ = architecture_table[architecture_table['anchor domain'] == arch]
        sns.scatterplot(data=df_, x=arch, y=key, ax=ax, edgecolor='k', hue='# of domains')
        ax.set_xlabel(xlabels[i])
        ax.set_ylabel('Other')
        ax.set_ylim((-1, 50))
    
    if len(name) > 0:
        plt.savefig(name.split('.')[0] + '_alpha_beta.svg', dpi=300)

    plt.show()
    
    