#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 14:04:30 2024

@author: christopher
"""


import seaborn as sns
import pandas as pd

from diversity_index import inverse_simpson_index


architecture_types = ['a', 'b', 'a+b', 'a/b']


def draw_domain_scatterplot(fig, ax, domains: pd.DataFrame):
    domains = domains[domains['anchor_domain_arch'] != 'None']
    gl_limit = max(domains['domain_1'].max(), domains['domain_2'].max()) +5
    
    left, bottom, width, height = (0.17, 0.32, 0.12, 0.12)
    ax_inner = fig.add_axes([left, bottom, width, height])
    
    sns.scatterplot(data=domains, x='domain_1', y='domain_2', edgecolor='k', ax=ax[1][0])
    ax[1][0].errorbar(x=domains['domain_1'], y=domains['domain_2'], fmt='none', 
                      xerr=domains['domain_1_err'], yerr=domains['domain_2_err'],
                      color='k')
    ax[1][0].plot([0, gl_limit], [0, gl_limit], 'k--', alpha=0.8)
    ax[1][0].set_xlabel('Domain 1')
    ax[1][0].set_ylabel('Domain 2')
    ax[1][0].set_xlim((-1, gl_limit))
    ax[1][0].set_ylim((-1, gl_limit))
    
    sns.scatterplot(data=domains, x='domain_1', y='domain_2', 
                    marker='+', ax=ax_inner)
    ax_inner.plot([0, gl_limit], [0, gl_limit], 'k--', alpha=0.8)
    ax_inner.set_xlabel('')
    ax_inner.set_ylabel('')
    ax_inner.set_xlim((-1, gl_limit))
    ax_inner.set_ylim((-1, gl_limit))
    

def draw_ab_barplot(ax, df_xgroup: pd.DataFrame):
    counts_architecture = df_xgroup.value_counts("architecture", normalize=True, dropna=False).to_dict()
    
    data = {
        'architectures': architecture_types,
        'fraction': [counts_architecture.get(a, 0) for a in architecture_types],
        }
    sns.barplot(data=data, x='architectures', y='fraction',
                edgecolor='k', color='cyan', width=0.5, ax=ax[1][1])
    ax[1][1].set_xlabel('Domain class')
    ax[1][1].set_ylabel('Fraction of all domains')
    
    
def draw_distance_plot(ax, df_atoms: pd.DataFrame):
    sns.histplot(data=df_atoms, x='distance', kde=True, color='purple',
                 stat='probability', ax=ax[0][0], edgecolor='k')
    ax[0][0].plot([3.5, 3.5], [0, 0.05], 'k--')
    ax[0][0].set_xlabel('Interaction distance [Å]')
    ax[0][0].set_ylabel('Probability')
    
    
def draw_diversity_index(ax, df: pd.DataFrame):
    n_domains = []
    diversity_indices = []
    
    for ligand in df['ligand_uuid'].unique():
        df_ = df[df['ligand_uuid'] == ligand]
        simpson_index = inverse_simpson_index(list(df_['xgroup_proportion']))
        diversity_indices.append(simpson_index)
        n_domains.append(len(df_))
        
    data = {
        'number of domains': n_domains,
        'inverse simpson index': diversity_indices,
        }
    sns.boxplot(data=data, x='number of domains', y='inverse simpson index',
                width=0.5, color='purple', ax=ax[0][1], meanprops={'linecolor': 'k'})
    ax[0][1].set_ylim((0.9, max(data['number of domains'])))
    ax[0][1].set_xlabel('# of interacting domains')
    ax[0][1].set_ylabel('Inverse simpson index (q=2)')

