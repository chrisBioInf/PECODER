#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 14:12:20 2024

@author: christopher
"""

import sys
import os
import glob 
import json
from networkx.readwrite import json_graph
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from diversity_index import inverse_simpson_index


distance_cutoff = 3.6


def read_json_file(filename):
    with open(filename) as f:
        js_graph = json.load(f)
    return json_graph.node_link_graph(js_graph)


def configuration_diversity_index(df: pd.DataFrame, moieties: set, filename='all'):
    xgroup_diversities = {
        'X-group configuration': [],
        'Inverse Simpson': [],
        }
    xconfig_counts = {
        'X-group configuration': [],
        'Count': [],
        }
    heat_data = {'X-group configuration': []}
    
    for group in moieties:
        heat_data[group] = []
    
    # Calculate X-group diversity indices
    for xgroup in df['Domain configuration'].unique():
        df_ = df[(df['Domain configuration'] == xgroup) & (df['Atom group'].isin(moieties))]
        total_interactions = len(df_)
        xconfig_counts['X-group configuration'].append(xgroup)
        xconfig_counts['Count'].append(total_interactions)
        if total_interactions == 0:
            continue
        group_proportions = [len(df_[df_['Atom group'] == atom_group]) / total_interactions 
                             for atom_group in moieties]
        x_diversity = inverse_simpson_index(group_proportions)
        xgroup_diversities['X-group configuration'].append(xgroup)
        xgroup_diversities['Inverse Simpson'].append(x_diversity)
        heat_data['X-group configuration'].append(xgroup)
        
        for atom_group in moieties:
            if atom_group in list(df_['Atom group']):
                heat_data[atom_group].append(len(df_[df_['Atom group'] == atom_group]) / total_interactions)
            else:
                heat_data[atom_group].append(0)
                
    sns.countplot(data=df, y='Domain configuration', edgecolor='k')
    plt.savefig('%s_xconfig_count.pdf' % filename, dpi=300, bbox_inches='tight')
    plt.show()
    
    fig, ax = plt.subplots(1, 2, figsize=(12, 25))
    
    xgroup_diversities = pd.DataFrame(data=xgroup_diversities)
    xgroup_diversities.sort_values(by='Inverse Simpson', ascending=False, inplace=True)
    xgroup_diversities.to_csv('%s_xgroup_diversities.tsv' % filename, index=False, sep='\t')
    order = list(xgroup_diversities['X-group configuration'])
    sns.barplot(data=xgroup_diversities, x='Inverse Simpson', y='X-group configuration', 
                edgecolor='k', width=0.5, order=order, ax=ax[0], color='darkcyan')
    ax[0].set_xlabel('Inverse Simpson index (q=2)')
    ax[0].set_ylabel('X-group configuration')
    
    heat_df = pd.DataFrame(data=heat_data)
    heat_df.set_index('X-group configuration', inplace=True, drop=True)
    heat_df = heat_df.reindex(order)
    heat_df.to_csv('%s_heat_map.tsv' % filename, index=False, sep='\t')
    
    sns.heatmap(data=heat_df, cmap="Blues",
                xticklabels=True, yticklabels=True)
    ax[1].set_xlabel('Moiety')
    ax[1].set_ylabel('')
    ax[1].set_xticklabels(ax[1].get_xticklabels(), rotation=45, ha='right')
    plt.savefig('%s_xconfig_diversity.svg' % filename, bbox_inches='tight', dpi=300)
    plt.savefig('%s_xconfig_diversity.pdf' % filename, bbox_inches='tight', dpi=300)
    plt.show()
    
    # pd.DataFrame(heat_data).to_csv('ECOD/NAD_heatmap.tsv', sep='\t')
    sns.clustermap(data=heat_df, cmap="Blues",
                   xticklabels=True, yticklabels=True, )
    plt.savefig('%s_xconfig_clustermap.svg' % filename, dpi=300, bbox_inches='tight')
    plt.savefig('%s_xconfig_clustermap.pdf' % filename, dpi=300, bbox_inches='tight')
    plt.show()


def draw_xgroup_diversity(df: pd.DataFrame, moieties: list,):
    xgroup_diversities = {
        'X-group': [],
        'Inverse Simpson': [],
        'Anchor': [],
        '# Cofactors': [],
        '# Configurations': [],
        '# Moieties':  [],
        }
    heat_data = {'X-group': []}
    
    for group in moieties:
        heat_data[group] = []
    
    # Calculate X-group diversity indices
    for xgroup in df['X-group'].unique():
        df_ = df[(df['X-group'] == xgroup) & (df['Atom group'].isin(moieties))]
        total_interactions = len(df_)
        if total_interactions == 0:
            continue
        group_proportions = [len(df_[df_['Atom group'] == atom_group]) / total_interactions 
                             for atom_group in moieties]
        x_diversity = inverse_simpson_index(group_proportions)
        xgroup_diversities['X-group'].append(xgroup)
        xgroup_diversities['Inverse Simpson'].append(x_diversity)
        xgroup_diversities['Anchor'].append(df_['Anchor domain'].mean())
        xgroup_diversities['# Cofactors'].append(len(df_['Cofactor'].unique()))
        xgroup_diversities['# Configurations'].append(len(df_['Domain configuration'].unique()))
        xgroup_diversities['# Moieties'].append(len(df_['Moiety'].unique()))
        
        heat_data['X-group'].append(xgroup)
        
        for atom_group in moieties:
            if atom_group in list(df_['Atom group']):
                heat_data[atom_group].append(1)
            else:
                heat_data[atom_group].append(0)
    
    fig, ax = plt.subplots(5, 1, figsize=(25, 18))
    x_df = pd.DataFrame(data=xgroup_diversities)
    x_df.sort_values(by='Inverse Simpson', ascending=False, inplace=True)
    order = list(x_df['X-group'])
    
    sns.barplot(data=x_df, x='X-group', y='Inverse Simpson', 
                edgecolor='k', width=0.7, color='coral', ax=ax[0], order=order)
    ax[0].set_ylabel('Moiety diversity')
    
    sns.barplot(data=x_df, x='X-group', y='Anchor', 
                edgecolor='k', width=0.7, ax=ax[1], order=order)
    ax[1].set_ylabel('Anchor domain fraction')
    
    sns.barplot(data=x_df, x='X-group', y='# Configurations', 
                edgecolor='k', width=0.7, color='teal', ax=ax[2], order=order)
    ax[2].set_ylabel('X-group configurations')
    
    sns.barplot(data=x_df, x='X-group', y='# Cofactors', 
                edgecolor='k', width=0.7, color='plum', ax=ax[3], order=order)
    ax[3].set_ylabel('Interacting cofactors')
    
    sns.barplot(data=x_df, x='X-group', y='# Moieties', 
                edgecolor='k', width=0.7, color='wheat', ax=ax[4], order=order)
    ax[4].set_ylabel('Interacting moieties')
    ax[4].tick_params(axis='x', labelrotation=90)
    
    for a in ax[:-1]:
        a.set_xlabel('')
        a.set(xticklabels=[])
    
    x_df.to_csv('xgroup_diversity.tsv', sep='\t', index=False)
    plt.savefig('xgroup_diversity_overview.svg', dpi=300, bbox_inches='tight')
    plt.savefig('xgroup_diversity_overview.png', dpi=500, bbox_inches='tight')
    plt.show()
    
    
    '''
    xgroup_diversities = pd.DataFrame(data=xgroup_diversities)
    xgroup_diversities.sort_values(by='Inverse Simpson', ascending=False, inplace=True)
    order = list(xgroup_diversities['X-group'])
    sns.barplot(data=xgroup_diversities, x='Inverse Simpson', y='X-group', 
                edgecolor='k', width=0.5, order=order, ax=ax[0], color='darkcyan')
    ax[0].set_xlabel('Inverse Simpson index (q=2)')
    ax[0].set_ylabel('X-group ID')
    
    heat_df = pd.DataFrame(data=heat_data)
    heat_df.set_index('X-group', inplace=True, drop=True)
    heat_df = heat_df.reindex(order)
    print(heat_df)
    
    sns.heatmap(data=heat_df, cmap="Blues",
                xticklabels=True, yticklabels=True)
    ax[1].set_xlabel('Atom group')
    ax[1].set_ylabel('')
    ax[1].set_xticklabels(ax[1].get_xticklabels(), rotation=45, ha='right')
    plt.savefig('%s_xgroup_diversity.svg' % filename, bbox_inches='tight', dpi=300)
    plt.show()
    
    # pd.DataFrame(heat_data).to_csv('ECOD/NAD_heatmap.tsv', sep='\t')
    sns.clustermap(data=heat_df, cmap="Blues",
                   xticklabels=True, yticklabels=True)
    plt.savefig('%s_xgroup_clustermap.svg' % filename, dpi=300, bbox_inches='tight')
    plt.show()
    '''
    

def draw_atom_group_diversity(diversities: list, atom_groups: list):
    data = {
        'Atom group': [],
        'Inverse Simpson': [],
        }
    for i in range(len(diversities)):
        for j in range(0, len(atom_groups)):
            data['Atom group'].append(atom_groups[j])
            data['Inverse Simpson'].append(diversities[i][j])

    sns.boxplot(data=pd.DataFrame(data=data), x='Inverse Simpson', y='Atom group', color='coral')
    plt.xlabel('Inverse Simpson index (q=2)')
    plt.savefig('atom_group_diversity.svg', dpi=300, bbox_inches='tight')
    plt.show()
    
    
def draw_moiety_diversity(df, moieties):
    data = {
        'Moiety': [],
        'Count': [],
        'Value': [],
        }
    for m in moieties:
        df_ = df[df['Atom group'] == m]
        data['Moiety'].append(m)
        data['Moiety'].append(m)
        data['Count'].append(len(df_['X-group'].unique()))
        data['Count'].append(len(df_['Domain configuration'].unique()))
        data['Value'].append('X-group')
        data['Value'].append('Domain config')

    data = pd.DataFrame(data=data)
    order = sorted(data.sort_values(by='Count', ascending=False)['Moiety'].unique())
    sns.barplot(data=data, y='Moiety', x='Count', hue='Value', edgecolor='k', order=order)
    plt.savefig('moiety_partner_count.png', dpi=300, bbox_inches='tight')
    plt.savefig('moiety_partner_count.svg', dpi=300, bbox_inches='tight')
    plt.show()


def calc_atom_group_diversity(df: pd.DataFrame, moiety_atoms: dict) -> list: 
    diversities = []
    
    for atom_group in moiety_atoms.keys():
        df_ = df[df['Moiety'] == atom_group]
        if len(df_) == 0:
            diversities.append(0)
            continue
        total_atoms = len(moiety_atoms.get(atom_group))
        atom_proportions = [len(df_[df_['atom_cofactor'] == atom]) for atom in moiety_atoms.get(atom_group)]
        diversities.append(inverse_simpson_index(atom_proportions))
        
    return diversities
    

def group_by_domains(df: pd.DataFrame, G, cofactor) -> pd.DataFrame:
    df = df[df['distance'] <= distance_cutoff]
    
    if len(df) == 0:
        return []
    
    df = df[(df['binding domains'] == 2) & (~df['domain configuration'].astype(str).str.contains('nan'))]
    group_dict = {
        'Atom': [],
        'Atom group': [],
        'X-group': [],
        'Domain configuration': [],
        
        'Anchor domain': [],
        'Cofactor': [],
        'Moiety': [],
        }
    atom_group_diversities = []
    
    moiety_map = {n[1]['name'] : n[1]['moiety'] for n in G.nodes(data=True)}
    moiety_atoms = {}
    
    for n in G.nodes(data=True):
        m = n[1]['moiety']
        if m not in moiety_atoms:
            moiety_atoms[m] = [n[1]['name']]
        else:
            moiety_atoms[m].append(n[1]['name'])
    
    for ligand in df['ligand_uuid'].unique():
        df_ = df[df['ligand_uuid'] == ligand]
        df_ = df_.dropna()

        # Aggregate by involved domains
        df_atom_domains = df_.groupby(by=['atom_cofactor', 'ecod_domain_id'], group_keys=False).agg({'atom_cofactor':'first',   
                                                                                                     'moiety': 'first',
                                                                                                     'ecod_domain_id': 'first',
                                                                                                     'distance': 'min',
                                                                                                     'ligand_uuid': 'first',
                                                                                                     'ecod_t_name': 'first',
                                                                                                     'x-group': 'first',
                                                                                                     'domain proportions': 'max',
                                                                                                     'domain configuration': 'first',
                                                                                                     })
        # Diversity of atom group
        df_atom_domains['Moiety'] = [moiety_map.get(x) for x in df_atom_domains['atom_cofactor']]
        atom_group_diversities.append(calc_atom_group_diversity(df_atom_domains, moiety_atoms))
        
        for idx, row in df_atom_domains.iterrows():
            atom_group = moiety_map.get(row['atom_cofactor'])
            x_group = row['x-group']
            group_dict['Atom'].append(row['atom_cofactor'])
            group_dict['Atom group'].append(atom_group)
            group_dict['X-group'].append(str(int(x_group)))
            group_dict['Domain configuration'].append(str(row['domain configuration']))
            group_dict['Cofactor'].append(cofactor)
            group_dict['Moiety'].append(row['moiety'])
            
            if row['domain proportions'] > 0.75:
                group_dict['Anchor domain'].append(1)
            else:
                group_dict['Anchor domain'].append(0)
    
    # draw_atom_group_diversity(atom_group_diversities, atom_groups)
    subgroup_df = pd.DataFrame(data=group_dict)
    return subgroup_df
        
        
def main():
    moieties = set()
    dfs = []
    
    for directory in sys.argv[1:]:
        cofactor = directory.split('/')[-1]
        graph_files = glob.glob(os.path.join(directory, '*.json'))
    
        for f in graph_files:
            G = read_json_file(f)
            df = pd.read_csv('%s.tsv' % f.split('_')[0], sep='\t')
            df = group_by_domains(df, G, cofactor)
            if len(df) > 0:
                dfs.append(df)
            else:
                continue
            atom_groups = set([n[1]['moiety'] for n in G.nodes(data=True)])
            
            for moiety in atom_groups:
                moieties.add(moiety)
        
    df = pd.concat(dfs)
    
    draw_moiety_diversity(df, moieties)
    # draw_xgroup_diversity(df, moieties,)
    #configuration_diversity_index(df, moieties,)
        

if __name__ == '__main__':
    main()
