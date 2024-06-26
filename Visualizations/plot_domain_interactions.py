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
from optparse import OptionParser

from architecture_binding_comparison import plot_architecture_comparison
from domain_pie import draw_anchor_domains
from overview_plots import draw_ab_barplot, draw_distance_plot, draw_diversity_index, draw_domain_scatterplot


__version__ = 0.5

tertiary_structure_values = {
    'a': (1, 0, 0, 0, 0),
    'b': (0, 1, 0, 0, 0),
    'a+b': (0, 0, 1, 0, 0),
    'a/b': (0, 0, 0, 1, 0),
    'other': (0, 0, 0, 0, 1),
    }

two_domain_anchor_cutoff = 0.75
three_domain_anchor_cutoff = 0.67
four_domain_anchor_cutoff = 0.5


def aggregate_all_tables(directory: str) -> pd.DataFrame:
    filelist = glob.glob(os.path.join(directory, '*.tsv'))
    dfs = []
    
    for file in filelist:
        df = pd.read_csv(file, sep='\t', dtype={'pdb': str, 'ligand_id': str, 
                                                'distance': np.float64, 'number': np.int32})
        
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
        'x_group_string': [],
        'domain_ids': [],
        'anchor domain': [], 
        '# of domains': [], 
        'domain 1': [],
        'domain 2': [],
        'domain_1_arch': [],
        'domain_2_arch': [],
        }
    
    for ligand in df['ligand_uuid'].unique():
        df_ = df[df['ligand_uuid'] == ligand]
        
        # Aggregate down to residue - cofactor interactions
        df_residues = df_.groupby(by=['residue', 'number'], group_keys=False).agg({'architecture':'first',
                                                                                   'residue': 'first',
                                                                                   'distance': 'min',
                                                                                   'ligand_uuid': 'first',
                                                                                   'ecod_x_name': 'first',
                                                                                   'ecod_domain_id': 'first',
                                                                                   't_id': 'first'})
        
        # Aggregate by involved domains
        df_ligand_xgroups = df_residues.groupby(by=['ecod_domain_id'], group_keys=False).agg({'architecture':'first',
                                                                                      'residue': 'count',
                                                                                      'distance': 'mean',
                                                                                      'ligand_uuid': 'first',
                                                                                      'ecod_x_name': 'first',
                                                                                      't_id': 'first'})
        
        df_ligand_xgroups['xgroup_proportion'] = df_ligand_xgroups['residue'] / df_ligand_xgroups['residue'].sum()
        df_ligand_xgroups.sort_index(inplace=True)
        dfs.append(df_ligand_xgroups)
        
        # Figure out if there is an actual anchor domain
        df_ligand_xgroups.sort_values(by='residue', inplace=True, ascending=False)
        
        if (len(df_ligand_xgroups) == 0):
            continue
        
        architecture_dict['x_group_string'].append(';'.join(list(df_ligand_xgroups['ecod_x_name'])))
        architecture_dict['domain_ids'].append('_'.join(list(df_ligand_xgroups['t_id'])))
        # ab_tuple = get_ab_values(list(df_ligand_xgroups['architecture']), list(df_ligand_xgroups['residue']))
        
        if len(df_ligand_xgroups) == 1:
            architecture_dict['anchor domain'].append(df_ligand_xgroups['architecture'].iloc[0])
            architecture_dict['# of domains'].append(len(df_ligand_xgroups))
            architecture_dict['domain 1'].append(df_ligand_xgroups['residue'].iloc[0])
            architecture_dict['domain 2'].append(0)
            architecture_dict['domain_1_arch'].append(df_ligand_xgroups['architecture'].iloc[0])
            architecture_dict['domain_2_arch'].append('None')
            
        elif (len(df_ligand_xgroups) == 2) and (df_ligand_xgroups['xgroup_proportion'].max() > two_domain_anchor_cutoff):
            architecture_dict['anchor domain'].append(df_ligand_xgroups['architecture'].iloc[0])
            architecture_dict['# of domains'].append(len(df_ligand_xgroups))
            architecture_dict['domain 1'].append(df_ligand_xgroups['residue'].iloc[0])
            architecture_dict['domain 2'].append(df_ligand_xgroups['residue'].iloc[1])
            architecture_dict['domain_1_arch'].append(df_ligand_xgroups['architecture'].iloc[0])
            architecture_dict['domain_2_arch'].append(df_ligand_xgroups['architecture'].iloc[1])
        
        elif (len(df_ligand_xgroups) == 3) and (df_ligand_xgroups['xgroup_proportion'].max() > three_domain_anchor_cutoff):
            architecture_dict['anchor domain'].append(df_ligand_xgroups['architecture'].iloc[0])
            architecture_dict['# of domains'].append(len(df_ligand_xgroups))
            architecture_dict['domain 1'].append(df_ligand_xgroups['residue'].iloc[0])
            architecture_dict['domain 2'].append(df_ligand_xgroups['residue'].iloc[1])
            architecture_dict['domain_1_arch'].append(df_ligand_xgroups['architecture'].iloc[0])
            architecture_dict['domain_2_arch'].append(df_ligand_xgroups['architecture'].iloc[1])
        
        elif (len(df_ligand_xgroups) > 3) and (df_ligand_xgroups['xgroup_proportion'].max() > four_domain_anchor_cutoff):
            architecture_dict['anchor domain'].append(df_ligand_xgroups['architecture'].iloc[0])
            architecture_dict['# of domains'].append(len(df_ligand_xgroups))
            architecture_dict['domain 1'].append(df_ligand_xgroups['residue'].iloc[0])
            architecture_dict['domain 2'].append(df_ligand_xgroups['residue'].iloc[1])
            architecture_dict['domain_1_arch'].append(df_ligand_xgroups['architecture'].iloc[0])
            architecture_dict['domain_2_arch'].append(df_ligand_xgroups['architecture'].iloc[1])
        
        else:
            architecture_dict['anchor domain'].append('None')
            architecture_dict['# of domains'].append(len(df_ligand_xgroups))
            architecture_dict['domain 1'].append(df_ligand_xgroups['residue'].iloc[0])
            architecture_dict['domain 2'].append(df_ligand_xgroups['residue'].iloc[1])
            architecture_dict['domain_1_arch'].append(df_ligand_xgroups['architecture'].iloc[0])
            architecture_dict['domain_2_arch'].append(df_ligand_xgroups['architecture'].iloc[1])

        
    architecture_table = pd.DataFrame(data=architecture_dict)
    print("Total number of ligand instances:", len(df['ligand_uuid'].unique()))
    
    return pd.concat(dfs), architecture_table


def aggregate_domain_interaction_motifs(architecture_table: pd.DataFrame):
    domain_motif_dict = {
        'anchor_domain_arch': [],
        'domain_motif': [],
        'domain_1': [],
        'domain_2': [],
        'domain_1_err': [],
        'domain_2_err': [],
        'domain_1_arch': [],
        'domain_2_arch': [],
        'sample_size': [],
        }
    for domain_motif in architecture_table['domain_ids'].unique():
        df_ = architecture_table[architecture_table['domain_ids'] == domain_motif]
        domain_motif_dict['domain_motif'].append(df_['domain_ids'].iloc[0])
        domain_motif_dict['anchor_domain_arch'].append(df_['anchor domain'].iloc[0])
        domain_motif_dict['domain_1'].append(df_['domain 1'].mean())
        domain_motif_dict['domain_2'].append(df_['domain 2'].mean())
        domain_motif_dict['domain_1_err'].append((df_['domain 1'].max() - df_['domain 1'].min()) / 2)
        domain_motif_dict['domain_2_err'].append((df_['domain 2'].max() - df_['domain 2'].min()) / 2)
        domain_motif_dict['domain_1_arch'].append(df_['domain_1_arch'].iloc[0])
        domain_motif_dict['domain_2_arch'].append(df_['domain_2_arch'].iloc[0])
        domain_motif_dict['sample_size'].append(len(df_))
    
    domain_motif_table = pd.DataFrame(data=domain_motif_dict)
    domain_motif_table.to_csv('domain_motif_table.tsv', sep='\t')
    print('Non-redundant domain binding motifs:', len(domain_motif_table))
    return domain_motif_table


def aggregate_by_domain_type(df: pd.DataFrame) -> pd.DataFrame:
    df_xgroup = df.groupby(by=['t_id'], group_keys=False).agg({'architecture':'first',
                                                               'atom_residue': 'count',
                                                               'distance': 'mean'})
    return df_xgroup
    

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
    df_xgroup = aggregate_by_domain_type(table)
    df_ligand_xgroups, architecture_table = aggregate_by_ligand_domain_interactions(table)
    df_motif_table = aggregate_domain_interaction_motifs(architecture_table)
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
    
    draw_domain_scatterplot(fig, ax, df_motif_table)
    draw_ab_barplot(ax, df_xgroup)
    draw_distance_plot(ax, table)
    draw_diversity_index(ax, df_ligand_xgroups)
    
    if len(options.name) > 0:
        plt.savefig(options.name, dpi=300, bbox_inches='tight')
        
    plt.show()
    draw_anchor_domains(options.name, df_motif_table)
    plot_architecture_comparison(options.name, df_motif_table)
    
    
if __name__ == "__main__":
    main()

