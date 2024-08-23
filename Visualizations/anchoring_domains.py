#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 10:30:05 2024

@author: christopher
"""


import sys
import os
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from diversity_index import inverse_simpson_index


def generate_fraction_plot():
    xs_d1 = []
    xs_d2 = []
    ys_d1 = []
    ys_d2 = []
    
    for n in range(0, 101):
        fraction = (100 -n) / 100
        xs_d1.append(fraction)
        xs_d2.append(abs(1-fraction))
        ys_d1.append(inverse_simpson_index([fraction, 1-fraction]))
        ys_d2.append(inverse_simpson_index([fraction, 1-fraction]))
        
    ds = ['Domain 1']*len(xs_d1)
    return xs_d1, ys_d1, ds


def aggregate_by_configurations(df, filename):
    df_agg = df.groupby(by=['domain configuration',]).agg(domain_configuration=pd.NamedAgg(column="domain configuration", aggfunc="first"),
                                                        binding_mode=pd.NamedAgg(column="binding mode", aggfunc="first"),
                                                        n_domains=pd.NamedAgg(column="binding domains", aggfunc="first"),
                                                        inverse_simpson=pd.NamedAgg(column="inverse simpson", aggfunc="mean"),
                                                        domain_fraction=pd.NamedAgg(column="domain proportions", aggfunc="max"),
                                                        architecture=pd.NamedAgg(column="architecture", aggfunc="first"),
                                                        )
    fig, axes = plt.subplots(1, 2, figsize=(10,4))
    
    binding_df = df_agg[df_agg['binding_mode'].isin(['A', 'A-A', 'A|A', 'A-B', 'A|B'])]
    binding_mode_counts = {}
    all_configs = len(binding_df)
    for bm in ['A', 'A-A', 'A|A', 'A-B', 'A|B']:
        binding_mode_counts[bm] = round(len(binding_df[binding_df['binding_mode'] == bm]) / all_configs, 4) *100
        
    print(binding_mode_counts)
    F_xs, F_ys, ds = generate_fraction_plot()
    F_vs_simpson = {
        'xs': F_xs,
        'ys': F_ys,
        'Domain': ds
        }
    sns.lineplot(data=F_vs_simpson, x='xs', y='ys', hue='Domain', 
                 palette=['teal'], #marker='D', 
                 ax=axes[0], dashes=False)
    axes[0].plot(F_vs_simpson.get('xs')[::10], F_vs_simpson.get('ys')[::10], 'D', color='teal')
    axes[0].set_xlabel('Fraction of interactions')
    axes[0].set_ylabel('Inverse simpson index (q=2)')
    
    df_agg['n_domains'] = df_agg['n_domains'].astype(int)
    
    sns.boxplot(data=df_agg[df_agg['n_domains'] < 6], x='n_domains', y='inverse_simpson',
                width=0.5, color='purple', ax=axes[1])
    axes[1].set_ylim((0.9, 4.1))
    axes[1].set_xlabel('# Binding domains')
    axes[1].set_ylabel('Inverse simpson index (q=2)')
    
    plt.savefig('%s_anchoring_domains.svg' % filename, dpi=300, bbox_inches='tight')
    plt.savefig('%s_anchoring_domains.png' % filename, dpi=300, bbox_inches='tight')
    plt.show()
    
    
    sns.histplot(data=df_agg[df_agg['n_domains'] == 2], x='inverse_simpson',
                 stat='count', edgecolor='k', bins=20)
    plt.plot([1.1, 1.1], [0, 15], 'k--')
    plt.plot([1.3, 1.3], [0, 15], 'k--')
    plt.plot([1.6, 1.6], [0, 15], 'k--')
    plt.xlabel('Inverse simpson index')
    plt.ylabel('X-group configurations')
    plt.savefig('inverse_simpson_counts.svg', dpi=300, bbox_inches='tight')
    plt.savefig('inverse_simpson_counts.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    fig, axes = plt.subplots(1, 2, figsize=(10,4))
    
    df_bm = df_agg[df_agg['binding_mode'].isin(['A-A', 'A|A', 'A-B', 'A|B'])]
    sns.boxplot(data=df_bm, y='binding_mode', x='inverse_simpson', width=0.5,
                color='teal', ax=axes[0])
    axes[0].set_xlim((0.9, 2.1))
    axes[0].set_ylabel('Binding mode')
    axes[0].set_xlabel('Inverse simpson (q=2)')
    # axes[1][0].set_yticks(['A-A', 'A|A', 'A-B', 'A|B'])
    
    df_arch = df_agg[(df_agg['binding_mode'].isin(['A-A', 'A|A', 'A-B', 'A|B'])) & (df_agg['architecture'].isin(['a', 'b', 'a+b', 'a/b']))]
    
    for i in range(0, len(df_arch)):
        if df_arch['domain_fraction'].iloc[i] >= 0.75:
            continue
        df_arch['architecture'].iloc[i] = 'None'
    
    sns.countplot(data=df_arch, x='architecture', edgecolor='k', width=0.5,
                  order=['None', 'a', 'b', 'a+b', 'a/b'], color='coral',
                  stat='probability', ax=axes[1])
    axes[1].set_xlabel('Anchor domain architecture')
    axes[1].set_ylabel('Fraction')
    axes[1].set_ylim((0, 1.1))
    
    plt.savefig('binding_mode_anchoring_domains.png', dpi=300, bbox_inches='tight')
    plt.show()

    

def main():
    files = []
    
    for f in sys.argv[1:]:
        files += glob.glob(os.path.join(f, '*.tsv'))
    
    filename = 'all'
    dfs = []
    
    for f in files:
        if len(f.split('_')) > 1:
            continue
        dfs.append(pd.read_csv(f, sep='\t'))
    
    df = pd.concat(dfs)
    df = df[~df['domain configuration'].astype(str).str.contains('nan')]
    # df = df[df['binding domains'] == 2]
    aggregate_by_configurations(df, filename)


main()
