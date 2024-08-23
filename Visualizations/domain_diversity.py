#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 11:57:26 2024

@author: christopher
"""


import os
import sys
import glob
import itertools
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn2 , venn2_circles


basepath = 'Results'
modes = ['A', 'A-A', 'A|A', 'A-B', 'A|B']
config_modes = ['-', '|']

cofactor_markers = {
    'COA': 'd',
    'NAD': '^',
    'ATP': 'p',
    'SAH': '<',
    'FAD': '>',
    }


def get_config_seperator(config):
    for sep in config_modes:
        if sep in config:
            return sep
    return ' '


def update_data_sets(args):
    dfs = []
    
    for c in args:
        pathname = os.path.join(basepath, c, '*_1.tsv')
        fs = glob.glob(pathname)
        
        for f in fs:
            reffile = f.split('_')[0] + '.tsv' 
            ligand = f.split('/')[-1].replace('.tsv', '')
            refdf = pd.read_csv(reffile, sep='\t')
            ligand_df = refdf[refdf['ligand_uuid'] == ligand]
            
            if len(ligand_df) == 0:
                continue
            
            pdb_id = ligand_df['pdb'].iloc[0]
            binding_mode = ligand_df['binding mode'].iloc[0]
            domain_configuration = str(ligand_df['domain configuration'].iloc[0])

            df = pd.read_csv(f, sep='\t')
            df['Cofactor'] = c
            df['Binding mode'] = binding_mode
            df['Domain configuration'] = domain_configuration
            df['PDB'] = pdb_id
            dfs.append(df)
            
    df = pd.concat(dfs)
    df.index = list(n for n in range(len(df)))
    df.to_csv('Results/domain_diversities.tsv', index=False, sep='\t')
    
    
def plot_binding_mode_vs_xgroup(args):
    data= {
        'Cofactor': [],
        'Binding mode': [],
        'X-groups': [],
        }
    for c in args:
        data['Cofactor'] += [c]*len(modes)
        data['Binding mode'] += modes
        dfs = []
        pathname = os.path.join(basepath, c, '*.tsv')
        fs = glob.glob(pathname)
        
        for f in fs:
            if len(f.split('_')) > 1:
                continue
            df = pd.read_csv(f, sep='\t')
            df['Cofactor'] = c
            dfs.append(df)
            
        df = pd.concat(dfs)
        total_xgroups = len(df['domain configuration'].unique())
        
        for mode in modes:
            df_ = df[df['binding mode'] == mode]
            
            if len(df_) == 0:
                data['X-groups'].append(0)
            else:
                data['X-groups'].append(len(df_['domain configuration'].unique()))
    
    sns.barplot(data=data, y='Binding mode', x='X-groups',
                edgecolor='k', hue='Cofactor', width=0.7, order=modes)
    plt.ylabel('Binding mode')
    plt.xlabel('X-group configurations')
    plt.savefig('binding_mode_vs_xgroup_config_cofactor.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    sns.barplot(data=data, y='Binding mode', x='X-groups',
                edgecolor='k', width=0.7, order=modes)
    plt.ylabel('Binding mode')
    plt.xlabel('X-group configurations')
    plt.savefig('binding_mode_vs_xgroup_config.png', dpi=300, bbox_inches='tight')
    plt.show()


def plot_domain_diversity(df):
    df = df[(df['Binding mode'].isin(modes)) & (~df['Domain configuration'].astype(str).str.contains('nan'))]
    sns.scatterplot(data=df, x='Inverse Simpson', y='Overlap index', hue='Binding mode',
                  palette=['teal', 'coral', 'purple', 'navy'], markers=cofactor_markers,
                  style='Cofactor', alpha=0.7)
    # g.plot_joint(sns.kdeplot, color="r", zorder=0, levels=6)
    # g.plot_marginals(sns.rugplot, color="r", height=-.15, clip_on=False)
    plt.plot([0, 2], [0.5, 0.5], 'k--')
    plt.plot([1.5, 1.5], [0, 1], 'k--')
    plt.ylim((-0.1, 1.1))
    plt.xlim((0.9, 2.1))
    plt.xlabel('Inverse Simpson (q=2)')
    plt.ylabel('Overlap')
    plt.savefig('domain_diversities_cutoff_2.svg', dpi=300, bbox_inches='tight')
    plt.show()

    
    df_agg = df.groupby(by=['Domain configuration', 'domain 1 atoms', 'domain 2 atoms']).agg({
        'Cofactor': 'first',
        'Domain configuration': 'first',
        'Domain 1': 'first',
        'Domain 2': 'first',
        'Inverse Simpson': 'first',
        'Overlap index': 'first',
        })
    for config in df_agg['Domain configuration'].unique():
        if 'nan' in config:
            continue
        df_ = df_agg[df_agg['Domain configuration'] == config]
        sep = get_config_seperator(config)
        df_['F-config'] = ['%s%s%s' % (dom_1, sep, dom_2) for dom_1,dom_2 in zip(df_['Domain 1'], df_['Domain 2'])]
        fig, ax = plt.subplots(1, 1)
        sns.scatterplot(data=df_, x='Inverse Simpson', y='Overlap index', 
                        hue='F-config', ax=ax, markers=cofactor_markers,
                        style='Cofactor', palette='colorblind', s=75)
        ax.plot([0, 2], [0.5, 0.5], 'k--')
        ax.plot([1.5, 1.5], [0, 1], 'k--')
        ax.set_xlim((0.9, 2.1))
        ax.set_ylim((-0.1, 1.1))
        plt.savefig('%s.svg' % (config), dpi=300, bbox_inches='tight')
    
    
def possible_domain_sets(args):
    existing_configs = []
    existing_xgroups = []
    all_configs = []
    
    for cofactor in args:
        pathname = os.path.join(basepath, cofactor, '*.tsv')
        fs = glob.glob(pathname)
        
        for f in fs:
            if len(f.split('_')) > 1:
                continue
            df = pd.read_csv(f, sep='\t')
            
            if len(df) == 0:
                continue
            
            df = df[df['binding domains'] <= 2]
            df['Cofactor'] = cofactor
            df['x-group'] = pd.Series([str(x) for x in df['x-group']]).astype(str)
            df['domain configuration'] = pd.Series([str(x) for x in df['domain configuration']]).astype(str)
            df = df[(~df['domain configuration'].astype(str).str.contains('nan')) & (~df['x-group'].astype(str).str.contains('nan'))]
            existing_configs += list(df['domain configuration'].unique())
            xgroups = list(df['x-group'].unique())
            existing_xgroups += xgroups
            all_configs += [str().join(config) for config in itertools.product(xgroups, config_modes, xgroups)]
    
    xconfigs = list(set(existing_configs))
    xgroups = list(set(existing_xgroups))
    possible_configs = list(set(all_configs))
    print('Total X-groups:', len(xgroups))
    print('Total X configs:', len(xconfigs))
    print(xgroups)
    print(xconfigs)
    print(possible_configs)
    venn2(subsets=(len(possible_configs), 0, len(xconfigs)),
          set_labels=('Combinatoric \nX-group configurations', 'Observed \nX-group \nconfigurations'),
          set_colors=('coral', 'teal', 'teal'), alpha=0.7)
    venn2_circles(subsets=(len(possible_configs), 0, len(xconfigs)))
    plt.savefig('combined_x_combinations.png' , dpi=300)
    plt.show()


def main():
    possible_domain_sets(sys.argv[1:])
    update_data_sets(sys.argv[1:])
    df = pd.read_csv('Results/domain_diversities.tsv', sep='\t')
    plot_binding_mode_vs_xgroup(sys.argv[1:])
    # df = df[df['Binding mode'].isin(modes[1:])]
    # plot_domain_diversity(df)
    
    
if __name__ == '__main__':
    main()

