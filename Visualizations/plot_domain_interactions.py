#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 14:54:32 2024

@author: christopher
"""

import os
import sys
import glob
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from optparse import OptionParser

from diversity_index import inverse_simpson_index


__version__ = 0.5

tertiary_structure_values = {
    'a': (1, 0),
    'b': (0, 1),
    'a+b': (1, 1),
    'a/b': (1, 1),
    'a+b,a/b': (1, 1),
    '-': (0, 0),
    }

architecture_types = ['a', 'b', 'a+b', 'a/b']


def aggregate_all_tables(directory: str) -> pd.DataFrame:
    filelist = glob.glob(os.path.join(directory, '*.tsv'))
    dfs = []
    
    for file in filelist:
        df = pd.read_csv(file, sep='\t')
        df['ligand_uuid'] = df['pdb'] + '_' + df['ligand_id']
        df['residue_uuid'] =  df['pdb'] + '_' + df['chain'] + '_' + str(df['number'])
        dfs.append(df)
        
    if (len(dfs) == 0):
        print('No valid files found.')
        sys.exit()
        
    return pd.concat(dfs)


def get_ab_values(architectures: list) -> tuple:
    alpha = 0
    beta = 0
    
    for arch in architectures:
        a, b = tertiary_structure_values.get(arch, (0, 0))
        alpha += a
        beta += b
        
    return (alpha, beta)


def aggregate_by_residue(df: pd.DataFrame) -> pd.DataFrame:   
    df_residue = df.groupby(by=['residue_uuid'], group_keys=False).agg({'architecture':'first',
                                                                        'ecod_x_name': 'first',
                                                                        'residue': 'first',
                                                                        'distance': 'min',
                                                                        'ligand_uuid': 'first'})
    return df_residue
 

def aggregate_by_ligand_domain_interactions(df: pd.DataFrame) -> (pd.DataFrame, dict):
    dfs = []
    architecture_dict = {}
    
    for ligand in df['ligand_uuid'].unique():
        df_ = df[df['ligand_uuid'] == ligand]
        
        df_ligand_xgroups = df_.groupby(by=['ecod_x_name'], group_keys=False).agg({'architecture':'first',
                                                                                   'residue': 'count',
                                                                                   'distance': 'mean',
                                                                                   'ligand_uuid': 'first',
                                                                                   'ecod_x_name': 'first'})
        
        df_ligand_xgroups['xgroup_proportion'] = df_ligand_xgroups['residue'] / df_ligand_xgroups['residue'].sum()
        df_ligand_xgroups.sort_index(inplace=True)
        dfs.append(df_ligand_xgroups)
        architecture_dict[';'.join(list(df_ligand_xgroups['ecod_x_name']))] = get_ab_values(list(df_ligand_xgroups['architecture']))
    
    print(architecture_dict)
    return pd.concat(dfs), architecture_dict


def aggregate_by_X_group(df: pd.DataFrame) -> pd.DataFrame:
    df_xgroup = df.groupby(by=['ecod_x_name'], group_keys=False).agg({'architecture':'first',
                                                                      'atom_residue': 'count',
                                                                      'distance': 'mean'})
    return df_xgroup
    
            
def draw_ab_scatterplot(ax, architecture_dict: dict):
    alpha_proportions = []
    beta_proportions = []
    n_domains = []
      
    for key, values in architecture_dict.items():
        totality = sum(values)
        
        if (totality > 0):
            alpha_proportions.append(values[0] / totality)
            beta_proportions.append(values[1] / totality)
            n_domains.append(key.count(';') +1)
    
    data = {
        'alpha': alpha_proportions,
        'beta': beta_proportions,
        '# domains': n_domains,
        }
    sns.scatterplot(data=data, x='alpha', y='beta', size='# domains', edgecolor='k', ax=ax[0][0])
    ax[0][0].plot([0, 1], [0, 1], 'k--', alpha=0.5)
    ax[0][0].set_xlabel('beta [# residues]')
    ax[0][0].set_ylabel('alpha [# residues]')
    ax[0][0].set_title('Non-redundant architecture compositions: %s' % len(data))


def draw_ab_barplot(ax, df_xgroup: pd.DataFrame):
    counts_architecture = df_xgroup.value_counts("architecture", normalize=True, dropna=False).to_dict()
    
    data = {
        'architectures': architecture_types,
        'fraction': [counts_architecture.get(a, 0) for a in architecture_types],
        }
    
    sns.barplot(data=data, x='architectures', y='fraction', 
                edgecolor='k', color='skyblue', width=0.5, ax=ax[0][1])
    ax[0][1].set_xlabel('Tertiary structure')
    ax[0][1].set_ylabel('Fraction')
    
    
def draw_distance_plot(ax, df_atoms: pd.DataFrame):
    sns.histplot(data=df_atoms, x='distance', kde=True, 
                 stat='probability', ax=ax[1][0])
    ax[1][0].set_xlabel('Interaction distance [Å]')
    ax[1][0].set_ylabel('Probability')
    
    
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
                width=0.5, color='darkorange', ax=ax[1][1])
    ax[1][1].set_xlabel('# of interacting domains')
    ax[1][1].set_ylabel('Inverse simpson index (q=2)')
    

def main():
    usage = "\plot_domain_interactions.py -i  "
    parser = OptionParser(usage=usage, version="__version__")
    args = sys.argv
    
    parser.add_option("-o", "--output", action="store", default="", type="string", dest="output", help="Output title prefix for created files, if any.")
    parser.add_option("-i", "--input", action="store", default="", type="string", dest="input", help="Input directory (Required).")
    options, args = parser.parse_args()
    
    if (len(options.input) == 0) or not (os.path.isdir(options.input)):
        print("--input must be a directory.")
        sys.exit()
    
    table = aggregate_all_tables(options.input)
    df_xgroup = aggregate_by_X_group(table)
    df_residues = aggregate_by_residue(table)
    df_ligand_xgroups, architecture_dict = aggregate_by_ligand_domain_interactions(df_residues)
    
    fig, ax = plt.subplots(2, 2, figsize=(10, 9))
    
    if (len(options.output) > 0):
        fig.suptitle("%s domain interaction" % options.output)
    
    draw_ab_scatterplot(ax, architecture_dict)
    draw_ab_barplot(ax, df_xgroup)
    draw_distance_plot(ax, table)
    draw_diversity_index(ax, df_ligand_xgroups)
    
    plt.show()
    
    
if __name__ == "__main__":
    main()

