#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 10:51:53 2024

@author: christopher
"""

import sys
import glob
import os
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from xgroup_interactions import group_by_domains
from xgroup_interactions import read_json_file
from binding_distance import jaccard_distance


def match_moiety_atoms():
    moieties = set()
    dfs = []
    
    for directory in sys.argv[1:]:
        cofactor = directory.split('/')[-1]
        graph_files = glob.glob(os.path.join(directory, '*.json'))
        moiety_reference_df = pd.read_csv(os.path.join('Compounds', '%s_map.tsv' % cofactor), sep='\t')
    
        for f in graph_files:
            G = read_json_file(f)
            df = pd.read_csv('%s.tsv' % f.split('_')[0], sep='\t')
            df = group_by_domains(df, G, cofactor)

            if len(df) > 0:
                df = df.merge(moiety_reference_df, on=['Moiety', 'Atom'])
                dfs.append(df)
            else:
                continue
            atom_groups = set([n[1]['moiety'] for n in G.nodes(data=True)])
            
            for moiety in atom_groups:
                moieties.add(moiety)
        
    df = pd.concat(dfs)
    df.to_csv('xgroup_with_moiety_atoms.tsv', sep='\t', index=False)


def main():
    match_moiety_atoms()
    df = pd.read_csv('xgroup_with_moiety_atoms.tsv', sep='\t')
    df.sort_values(by='Reference atom', inplace=True, ascending=True)
    
    #sns.countplot(data=df, y='Reference atom', hue='Moiety', edgecolor='k')
    #plt.show()
    
    # Binding vectors
    atom_list = sorted(df['Reference atom'].unique())
    data = {
        'X-group': [],
        #'Configuration': [],
        'Cofactor': [],
        }    
    for atom in atom_list:
        data[atom] = []
    
    df['agg_key'] = ['%s_%s' % (row['X-group'], row['Cofactor']) for i,row in df.iterrows()]
    
    for key in df['agg_key'].unique():
        xgroup,  cofactor = key.split('_')
        df_ = df[(df['agg_key'] == key)]
        data['X-group'].append(str(xgroup))
        # data['Configuration'].append(str(xconfig))
        data['Cofactor'].append(cofactor)
        
        for atom in atom_list:
            moiety = atom.split('_')[0]
            if moiety not in list(df_['Moiety']):
                data[atom].append(-1)
            elif atom in list(df_['Reference atom']):
                data[atom].append(1)
            else: 
                data[atom].append(0)
    
    df = pd.DataFrame(data=data)
    df['X-group'] = [str(x) for x in df['X-group']]
    df.sort_values(by=['Cofactor', 'X-group'], inplace=True)
    df.to_csv('binding_mode_alignment.tsv', sep='\t', index=False)
    
    print(jaccard_distance(v1=list(df[(df['X-group'] == '2002') & (df['Cofactor'] == 'NAD')].iloc[0][atom_list]), 
                           v2=list(df[(df['X-group'] == '2002') & (df['Cofactor'] == 'COA')].iloc[0][atom_list])))
    
    sns.clustermap(data=df[sorted(atom_list)], row_cluster=True, col_cluster=False,
                   xticklabels=sorted(atom_list), yticklabels=list(df['X-group']), )
    plt.show()
    
    #
    #fig, ax = plt.subplots(1, 1)
    #ax.spy(matrix)
    #plt.show()
    
    '''
    data = {
        'X-group': [],
        'Jaccard distance': [],
        }
    for idx, row in df.iterrows():
        for xgroup in xgroups:
            data['X-group'].append(str(xgroup))
            data['Jaccard distance'].append(jaccard_distance(row.values.tolist(), df.loc[str(xgroup)]))
    
    sns.boxplot(data=data, y='X-group', x='Jaccard distance')
    
    plt.show()
    '''
    
    # Atom Bias
    '''
    data = {
        'X-group': [],
        'Atom bias magnitude': [],
        }
    for xgroup in df['X-group'].unique():
        df_ = df[df['X-group'] == xgroup]
        total_interactions = len(df_)
        atom_proportions = []
        interacting_atoms = len(df_['Atom'].unique())
        null_probability = 1 / interacting_atoms
        
        for atom in df_['Atom'].unique():
            df_atom = df_[df_['Atom'] == atom]
            atom_proportions.append(abs((len(df_atom) / total_interactions) - null_probability))
        
        data['X-group'].append(str(xgroup))
        data['Atom bias magnitude'].append(sum(atom_proportions) / interacting_atoms)
        
    # data = pd.DataFrame(data=data)   
    sns.barplot(data=data, y='X-group', x='Atom bias magnitude', edgecolor='k')
    plt.show()
    '''
    
    
main()
    