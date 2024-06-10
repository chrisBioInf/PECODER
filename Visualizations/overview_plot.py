#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 12:02:52 2024

@author: christopher
"""

"""
For now, the following script serves only for illustration purposes and
is very much NOT cleaned/optimized.

Use at your own considerations.
        
"""


import os
import sys
import glob
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


tertiary_structure_values = {
    'a': (1, 0),
    'b': (0, 1),
    'a+b': (1, 1),
    'a/b': (1, 1),
    'a+b,a/b': (1, 1),
    '-': (0, 0),
    }

architecture_types = ['a', 'b', 'a+b', 'a/b']

working_dir = sys.argv[1]
title = sys.argv[2]
filelist = glob.glob(os.path.join(working_dir, '*.tsv'))

a_s = []
b_s = []
alpha_sum = []
beta_sum = []
n_res = []
n_chains = []
n_domains = []
min_distances = []
max_distances = []

architectures = []

domain_ligands = []
domain_names = []
domain_interactions = []
domain_distances = []

count = 0

for file in filelist:
    df = pd.read_csv(file, sep='\t')
    
    for ligand in df['ligand_id'].unique():
        count += 1
        df_ = df[df['ligand_id'] == ligand]
        
        df_agg = df_.groupby(by=['chain', 'number'], group_keys=False).agg({'architecture':'first',
                                                                            'residue': 'first',
                                                                            'distance': 'mean'})
        alpha = 0
        beta = 0
        
        for i, row in df_agg.iterrows():
            ab_values = tertiary_structure_values.get(row['architecture'], (0, 0))
            alpha += ab_values[0]
            beta += ab_values[1]
            architectures.append(row['architecture'])
        
        alpha_sum.append(alpha)
        beta_sum.append(beta)
        a_s.append(alpha / len(df_agg))
        b_s.append(beta / len(df_agg))
        n_res.append(len(df_agg))
        n_chains.append(len(df_['chain'].unique()))
        n_domains.append(len(df_['ecod_domain_id'].dropna().unique()))
        min_distances.append(df_agg['distance'].min())
        max_distances.append(df_agg['distance'].max())
        
        df_domains = df_.groupby(by=['ecod_domain_id'], group_keys=False).agg({'ecod_t_name': 'first',
                                                                     'distance': 'mean',
                                                                     'atom_residue': 'count'})
        for i, row in df_domains.iterrows():
            domain_ligands.append('ligand_%s' % count)
            domain_names.append(row['ecod_t_name'])
            domain_distances.append(row['distance'])
            domain_interactions.append(row['atom_residue'])
        
    
data = {
        'alpha': a_s,
        'beta': b_s,
        'alpha [sum]': alpha_sum,
        'beta [sum]': beta_sum,
        '# Residues': n_res,
        '# Interacting chains': n_chains,
        '# Interacting domains': n_domains,
        'Minimal distance [Å]': min_distances,
        'Maximal distance [Å]': max_distances,
        }

df_plot = pd.DataFrame(data=data)

fig, ax = plt.subplots(2, 2, figsize=(10, 9))
fig.suptitle("%s domain interaction" % title )

sns.scatterplot(data=df_plot, x='beta [sum]', y='alpha [sum]', size='# Interacting domains', ax=ax[0][0])
ax[0][0].plot([0, max(df_plot['beta [sum]'])], [0, max(df_plot['alpha [sum]'])], 'k--', alpha=0.5)
ax[0][0].set_xlabel('beta [# residues]')
ax[0][0].set_ylabel('alpha [# residues]')


architecture_data = {
    'architectures': architectures,
    }

df_architecture = pd.DataFrame(data=architecture_data)
counts_architecture = df_architecture.value_counts("architectures", normalize=True, dropna=False).to_dict()

architecture_data = {
    'architectures': architecture_types,
    'fraction': [counts_architecture.get(a, 0) for a in architecture_types],
    }

sns.barplot(data=architecture_data, x='architectures', y='fraction', color='lightblue', edgecolor='k', 
              width=0.6, ax=ax[0][1])
ax[0][1].set_xlabel('Tertiary structure')
ax[0][1].set_ylabel('Fraction of interactions')


distance_data = {
    'Distance [Å]': min_distances + max_distances,
    '# Interacting domains': n_domains*2,
    'Agg': ['Minimal']*len(min_distances) + ['Maximal']*len(max_distances),  
    }

sns.boxplot(data=distance_data, x='# Interacting domains', y='Distance [Å]', hue='Agg',
            ax=ax[1][0],)
ax[1][0].set_xlabel('# of interacting domains')
ax[1][0].set_ylabel('Interaction distance [Å]')


distance_data = {
    'Ligand idx': domain_ligands,
    'Domain': domain_names,
    'Distance [Å]': domain_distances,
    'Interactions': domain_interactions,
    }
df_domains = pd.DataFrame(data=distance_data)


interactions_per_ligand = []
diversity_index = []

for idx in df_domains['Ligand idx'].unique():
    df_ = df_domains[df_domains['Ligand idx'] == idx]
    sum_of_interactions = df_.Interactions.sum()
    
    interaction_fraction = []
    
    if len(df_) == 1:
        interactions_per_ligand.append(1)
        diversity_index.append(1)
        continue
    
    for n in range(0, len(df_)):
        interaction_fraction.append(round(df_['Interactions'].iloc[n] / sum_of_interactions, 2))
            
    diversity_set = [abs(a-b) for a, b in itertools.combinations(interaction_fraction, 2)]
    interactions_per_ligand.append(len(df_))
    diversity_index.append(max(diversity_set))
    # print(interaction_fraction, diversity_set, 1 - max(diversity_set))



"""
diversity_data = {
    'Polarity index': diversity_index,
    '# of domains': interactions_per_ligand,
    }

sns.boxplot(data=diversity_data, x='# of domains', y='Polarity index', ax=ax[1][1], width=0.5)
ax[1][1].set_ylim((-0.1, 1.1))
ax[1][1].set_xlabel('# of domains')
ax[1][1].set_ylabel('Polarity index')
"""

"""
sns.scatterplot(data=domain_data, x='Domain 1', y='Domain 2', hue='Domain 3', ax=ax[1][1])
ax[1][1].set_xlim((-0.1, 1.1))
ax[1][1].set_ylim((-0.1, 1.1))
ax[1][1].plot([1,0], [0, 1], 'k--', alpha=0.5) 
"""

plt.savefig('%s_overview_plot.svg' % title, bbox_inches='tight', dpi=300)

plt.show()
    
