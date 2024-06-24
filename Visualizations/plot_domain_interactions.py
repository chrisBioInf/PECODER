#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 14:54:32 2024

@author: christopher
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from optparse import OptionParser

from diversity_index import inverse_simpson_index
from architecture_binding_comparison import plot_architecture_comparison
from domain_pie import draw_anchor_domains


__version__ = 0.5

tertiary_structure_values = {
    'a': (1, 0, 0, 0, 0),
    'b': (0, 1, 0, 0, 0),
    'a+b': (0, 0, 1, 0, 0),
    'a/b': (0, 0, 0, 1, 0),
    'other': (0, 0, 0, 0, 1),
    }

architecture_types = ['a', 'b', 'a+b', 'a/b']

two_domain_anchor_cutoff = 0.75
three_domain_anchor_cutoff = 0.67
four_domain_anchor_cutoff = 0.5


def aggregate_all_tables(directory: str) -> pd.DataFrame:
    filelist = glob.glob(os.path.join(directory, '*.tsv'))
    dfs = []
    
    for file in filelist:
        df = pd.read_csv(file, sep='\t')
        
        if (len(df) == 0):
            continue
        
        try:
            df['ligand_uuid'] = df['pdb'] + '_' + df['ligand_id']
            df['residue_uuid'] =  df['pdb'] + '_' + df['chain'] + '_' + str(df['number'])
            dfs.append(df)
        except Exception:
            print("Encountered exception at: %s" % file)
            print(df)
        
    if (len(dfs) == 0):
        print('No valid files found.')
        sys.exit()
    
    return pd.concat(dfs)


def get_ab_values(architectures: list, residues: list) -> tuple:
    a = 0
    b = 0
    a_plus_b = 0
    a_b = 0
    other = 0
    
    for i in range(0, len(architectures)):
        if (architectures[i] == 'a'):
            a += residues[i]
        elif (architectures[i] == 'b'):
            b += residues[i]
        elif (architectures[i] == 'a+b'):
            a_plus_b += residues[i]
        elif (architectures[i] == 'a/b'):
            a_b += residues[i]
        else:
            other += residues[i]
    
    return (a, b, a_plus_b, a_b, other)


def aggregate_by_residue(df: pd.DataFrame) -> pd.DataFrame:   
    df_residue = df.groupby(by=['residue_uuid'], group_keys=False).agg({'architecture':'first',
                                                                        'ecod_x_name': 'first',
                                                                        'residue': 'first',
                                                                        'distance': 'min',
                                                                        'ligand_uuid': 'first'})
    return df_residue
 

def aggregate_by_ligand_domain_interactions(df: pd.DataFrame) -> (pd.DataFrame, pd.DataFrame, pd.DataFrame, dict):
    dfs = []
    
    architecture_dict = {
        'architecture_string': [],
        'anchor domain': [],
        '# of domains': [], 
        'a': [],
        'b': [],
        'a+b': [],
        'a/b': [],
        'other': [],
        }
    two_domain_dict = {
        'Domain 1': [],
        'Domain 2': [],
        }
    anchor_domain_dict = {
        'a': 0,
        'b': 0,
        'a+b': 0,
        'a/b': 0,
        'a+b,a/b': 0,
        'No anchor': 0
        }
    
    interactions_per_ligand = []
    
    for ligand in df['ligand_uuid'].unique():
        df_ = df[df['ligand_uuid'] == ligand]
        interactions_per_ligand.append(len(df_))
        
        # Aggregate by involved domains
        df_ligand_xgroups = df_.groupby(by=['ecod_domain_id'], group_keys=False).agg({'architecture':'first',
                                                                                      'residue': 'count',
                                                                                      'distance': 'mean',
                                                                                      'ligand_uuid': 'first',
                                                                                      'ecod_x_name': 'first'})
        
        df_ligand_xgroups['xgroup_proportion'] = df_ligand_xgroups['residue'] / df_ligand_xgroups['residue'].sum()
        df_ligand_xgroups.sort_index(inplace=True)
        dfs.append(df_ligand_xgroups)
        
        architecture_dict['architecture_string'].append(';'.join(list(df_ligand_xgroups['ecod_x_name'])))
        ab_tuple = get_ab_values(list(df_ligand_xgroups['architecture']), list(df_ligand_xgroups['residue']))
        
        for arch, val in zip(list(tertiary_structure_values.keys()), ab_tuple):
            architecture_dict[arch].append(val)
        
        # Figure out if there is an actual anchor domain
        df_ligand_xgroups.sort_values(by='residue', inplace=True, ascending=False)
        
        if len(df_ligand_xgroups) == 1:
            anchor_domain_dict[df_ligand_xgroups['architecture'].iloc[0]] += 1
            architecture_dict['anchor domain'].append(df_ligand_xgroups['architecture'].iloc[0])
            architecture_dict['# of domains'].append(len(df_ligand_xgroups))
            
        elif (len(df_ligand_xgroups) == 2) and (df_ligand_xgroups['xgroup_proportion'].max() > two_domain_anchor_cutoff):
            anchor_domain_dict[df_ligand_xgroups['architecture'].iloc[0]] += 1
            architecture_dict['anchor domain'].append(df_ligand_xgroups['architecture'].iloc[0])
            architecture_dict['# of domains'].append(len(df_ligand_xgroups))
        
        elif (len(df_ligand_xgroups) == 3) and (df_ligand_xgroups['xgroup_proportion'].max() > three_domain_anchor_cutoff):
            anchor_domain_dict[df_ligand_xgroups['architecture'].iloc[0]] += 1
            architecture_dict['anchor domain'].append(df_ligand_xgroups['architecture'].iloc[0])
            architecture_dict['# of domains'].append(len(df_ligand_xgroups))
        
        elif (len(df_ligand_xgroups) > 3) and (df_ligand_xgroups['xgroup_proportion'].max() > four_domain_anchor_cutoff):
            anchor_domain_dict[df_ligand_xgroups['architecture'].iloc[0]] += 1
            architecture_dict['anchor domain'].append(df_ligand_xgroups['architecture'].iloc[0])
            architecture_dict['# of domains'].append(len(df_ligand_xgroups))
        
        else:
            anchor_domain_dict['No anchor'] += 1
            architecture_dict['anchor domain'].append('None')
            architecture_dict['# of domains'].append(len(df_ligand_xgroups))

        
        # Assemble domain1/domain2 interaction counts
        if (len(df_ligand_xgroups) < 2):
            continue
        two_domain_dict['Domain 1'].append(df_ligand_xgroups['residue'].iloc[0])
        two_domain_dict['Domain 2'].append(df_ligand_xgroups['residue'].iloc[1])
        
    
    architecture_table = pd.DataFrame(data=architecture_dict)
    domain_table = pd.DataFrame(data=two_domain_dict)
    print(architecture_table)
    print(domain_table)
    print("Total number of ligand instances:", len(df['ligand_uuid'].unique()))
    print("Average interactions per ligand:", np.mean(interactions_per_ligand))
    
    return pd.concat(dfs), architecture_table, domain_table, anchor_domain_dict


def aggregate_by_X_group(df: pd.DataFrame) -> pd.DataFrame:
    df_xgroup = df.groupby(by=['ecod_x_name'], group_keys=False).agg({'architecture':'first',
                                                                      'atom_residue': 'count',
                                                                      'distance': 'mean'})
    return df_xgroup
    

def draw_domain_scatterplot(fig, ax, domains: pd.DataFrame):
    gl_limit = max(domains['Domain 1'].max(), domains['Domain 2'].max())
    
    left, bottom, width, height = (0.17, 0.32, 0.12, 0.12)
    ax_inner = fig.add_axes([left, bottom, width, height])
    
    sns.scatterplot(data=domains, x='Domain 1', y='Domain 2', color='cyan', edgecolor='k', ax=ax[1][0])
    ax[1][0].plot([0, gl_limit], [0, gl_limit], 'k--', alpha=0.8)
    ax[1][0].set_xlabel('Domain 1')
    ax[1][0].set_ylabel('Domain 2')
    ax[1][0].set_xlim((-1, gl_limit))
    ax[1][0].set_ylim((-1, gl_limit))
    
    sns.scatterplot(data=domains[domains['Domain 1'] <= 100], x='Domain 1', y='Domain 2', 
                    color='cyan', marker='+', ax=ax_inner)
    ax_inner.plot([0, 100], [0, 100], 'k--', alpha=0.8)
    ax_inner.set_xlabel('')
    ax_inner.set_ylabel('')
    ax_inner.set_xlim((-5, 100))
    ax_inner.set_ylim((-5, 100))
    

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
    ax[0][1].set_ylim((0.9, 4.1))
    ax[0][1].set_xlabel('# of interacting domains')
    ax[0][1].set_ylabel('Inverse simpson index (q=2)')
    

def main():
    usage = "\plot_domain_interactions.py -i  "
    parser = OptionParser(usage=usage, version="__version__")
    args = sys.argv
    
    parser.add_option("-o", "--output", action="store", default="", type="string", dest="output", help="Output title prefix for created files, if any.")
    parser.add_option("-i", "--input", action="store", default="", type="string", dest="input", help="Input directory (Required).")
    parser.add_option("--name", action="store", default="", type="string", dest="name", help="Name of output plot for saving (if any).")
    options, args = parser.parse_args()
    
    if (len(options.input) == 0) or not (os.path.isdir(options.input)):
        print("--input must be a directory.")
        sys.exit()
    
    table = aggregate_all_tables(options.input)
    print("Non-redundant domain types:", len(table['ecod_x_name'].unique()))
    df_xgroup = aggregate_by_X_group(table)
    df_ligand_xgroups, architecture_table, two_domain_table, anchor_domain_dict = aggregate_by_ligand_domain_interactions(table)
    print("Average interactions per domain:", df_ligand_xgroups['residue'].mean())
    
    params = {"ytick.color" : "k",
              "xtick.color" : "k",
              "axes.labelcolor" : "k",
              "axes.edgecolor" : "k"}
    plt.rcParams.update(params)
    
    fig, ax = plt.subplots(2, 2, figsize=(10, 9))
    fig.set_facecolor("w")
    
    if (len(options.output) > 0):
        fig.suptitle("%s domain interaction" % options.output)
    
    draw_domain_scatterplot(fig, ax, two_domain_table)
    draw_ab_barplot(ax, df_xgroup)
    draw_distance_plot(ax, table)
    draw_diversity_index(ax, df_ligand_xgroups)
    
    if len(options.name) > 0:
        plt.savefig(options.name, dpi=300, bbox_inches='tight', transparent=True)
        
    plt.show()
    draw_anchor_domains(options.name, anchor_domain_dict)
    plot_architecture_comparison(options.name, architecture_table)
    
    
if __name__ == "__main__":
    main()

