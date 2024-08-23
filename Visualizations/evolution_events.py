#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 13:23:53 2024

@author: christopher
"""


import json
import re
import itertools
import glob
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import networkx as nx
from networkx.readwrite import json_graph


mutation_types = ['C-terminus insertion/deletion', 'N-terminus insertion/deletion', 
                  'Homo-oligomerization', 'Hetero-oligomerization', 'Fusion/Fission', 
                  'Duplication', 'Permutation']


def load_data(name):
    dfs = []
    fs = glob.glob('%s/*.tsv' % name)
    
    for f in fs:
        if len(f.split('_')) > 1:
            continue
        dfs.append(pd.read_csv(f, sep='\t'))
    
    return pd.concat(dfs)


def build_mutation_graph(df, configurations, filename):
    G = nx.Graph()
    
    for config in configurations:
        G.add_node(config, label=config)
    
    node_dict = {
        n[1]['label'] : n[0] for n in G.nodes(data=True)
        }
    for idx, row in df.iterrows():
        G.add_edge(node_dict.get(row['Configuration 1']), 
                   node_dict.get(row['Configuration 2']),
                   mutation=row['Operation'],
                   )
    nx.draw(G, with_labels=True, pos=nx.kamada_kawai_layout(G),
            edgecolors=['k']*len(G.nodes))
    plt.savefig('%s_mutation_graph.svg' % filename, dpi=1000)
    plt.savefig('%s_mutation_graph.png' % filename, dpi=1000)
    plt.show()
    
    fig, axes = plt.subplots(4, 1, figsize=(16, 14))
    
    data = {'X-group configuration': [],
            'Degree': []}
    
    for node in G.nodes(data=True):
        data['X-group configuration'].append(node[1]['label'])
        data['Degree'].append(G.degree(node[0]))
        
    sns.barplot(data=data, x='X-group configuration', y='Degree', edgecolor='k',
                ax=axes[0])
    axes[0].set_xlabel('')
    axes[0].set(xticklabels=[])
    
    xgroups = set()
    for config in configurations:
        domains = config.replace('|', '-').split('-')
        for domain in domains:
            xgroups.add(domain)
    
    data = {'X-group configuration': [],
            'Betweenness centrality': [],
            'Degree centrality': [],
            'Eigenvector centrality': []}  
    
    e_centralities = nx.eigenvector_centrality(G)
    d_centralities = nx.degree_centrality(G)
    b_centralities = nx.betweenness_centrality(G)
    
    for node in G.nodes(data=True):
        xgroup = node[1]['label']
        node[1]['eigenvector_centrality'] = e_centralities.get(node[0])
        node[1]['degree_centrality'] = d_centralities.get(node[0])
        node[1]['betweenness_centrality'] = b_centralities.get(node[0])
        centrality = e_centralities.get(node[0])
        data['X-group configuration'].append(xgroup)
        data['Eigenvector centrality'].append(centrality)
        data['Degree centrality'].append(d_centralities.get(node[0]))
        data['Betweenness centrality'].append(b_centralities.get(node[0]))
    
    sns.barplot(data=data, x='X-group configuration', y='Degree centrality', edgecolor='k', 
                color='pink', ax=axes[1])
    axes[1].set_xlabel('')
    axes[1].set(xticklabels=[])
    
    sns.barplot(data=data, x='X-group configuration', y='Betweenness centrality', edgecolor='k', 
                color='teal', ax=axes[2])
    axes[2].set_xlabel('')
    axes[2].set(xticklabels=[])
    
    sns.barplot(data=data, x='X-group configuration', y='Eigenvector centrality', edgecolor='k', 
                color='coral', ax=axes[3])
    axes[3].set_xlabel('')
    axes[3].set_xticklabels(axes[3].get_xticklabels(), rotation=90, ha='right')
    
    plt.savefig('%s_centrality.svg' % filename, dpi=300, bbox_inches='tight')
    plt.show()
    
    # Write data to json
    json_representation = json_graph.node_link_data(G)
    json.dump(obj=json_representation, fp=open(filename + '_configuration_graph.json', 'w'))


def map_domain_changes(domain_configurations):
    config_pairs = itertools.product(domain_configurations, repeat=2)
    configuration_left = []
    configuration_right = []
    mutations = []
    
    for (config1, config2) in config_pairs:
        if (config1 == config2):
            continue
        if len(config1.replace('|', '-').split('-')) == 0 and len(config2.replace('|', '-').split('-')) == 0:
            continue
        
        domains_1 = config1.replace('|', '-').split('-')
        domains_2 = config2.replace('|', '-').split('-')
        
        if len(set(domains_1).intersection(set(domains_2))) == 0:
            continue
        
        seperator_1 = re.sub(r'[0-9]+', '', config1)
        seperator_2 = re.sub(r'[0-9]+', '', config2)
        
        # CASE 1: 2 domains on each side
        if len(domains_1) == len(domains_2):
            if (seperator_1 != seperator_2) and (domains_1 == domains_2):
                mutations.append('Fusion/Fission')
            elif domains_1 == [d for d in reversed(domains_2)]: 
                mutations.append('Permutation')
            else:
                continue
        
        # CASE 2: Left side has only 1 domain
        elif len(domains_1) < len(domains_2):
            if (domains_2[0] == domains_2[1]) and (seperator_2 == '-'):
                mutations.append('Duplication')
            elif (domains_2[0] == domains_2[1]) and (seperator_2 == '|'):
                mutations.append('Homo-oligomerization')
            elif domains_1[0] == domains_2[0] and (seperator_2 == '-'):
                mutations.append('C-terminus insertion/deletion')
            elif domains_1[0] == domains_2[1] and (seperator_2 == '-'):
                mutations.append('N-terminus insertion/deletion')
            else:
                mutations.append('Hetero-oligomerization')
        
        # CASE 3: Right side has only 1 domain
        elif len(domains_1) > len(domains_2):
            if (domains_1[0] == domains_1[1]) and (seperator_1 == '-'):
                mutations.append('Duplication')
            elif (domains_1[0] == domains_1[1]) and (seperator_1 == '|'):
                mutations.append('Homo-oligomerization')
            elif domains_1[0] == domains_2[0] and (seperator_1 == '-'):
                mutations.append('C-terminus insertion/deletion')
            elif domains_1[1] == domains_2[0] and (seperator_1 == '-'):
                mutations.append('N-terminus insertion/deletion')
            else:
                mutations.append('Hetero-oligomerization')
            
        if len(mutations) > len(configuration_left):
            configuration_left.append(config1)
            configuration_right.append(config2)
        
        for i in range(len(mutations)):
            print('%s <-> %s   [%s]' % (configuration_left[i], configuration_right[i], mutations[i]))
        
    df = pd.DataFrame(data={
            'Configuration 1': configuration_left,
            'Configuration 2': configuration_right,
            'Operation': mutations,
            })
    return df
            

def main():
    dfs = []
    
    for arg in sys.argv[1:]:
        df = load_data(arg)
        df['domain configuration'] = [str(x).split('.')[0] for x in df['domain configuration']]
        df = df[(df['binding domains'] < 3) & (~df['domain configuration'].str.contains('nan'))]
        df_agg = df.groupby(by='domain configuration').agg({
            'domain configuration': 'first',
            'binding mode': 'first',
            })
        dfs.append(map_domain_changes(df_agg['domain configuration']))
        dfs[-1]['Cofactor'] = arg.split('/')[-1]
    
    df = pd.concat(dfs)
    filename = arg.split('/')[-1]
    build_mutation_graph(df, df_agg['domain configuration'], filename)
    df.to_csv('Results/mutations.tsv', sep='\t', index=False)
    sns.countplot(data=df, x='Operation', edgecolor='k', hue='Cofactor', width=0.5,
                  order=mutation_types)
    plt.xticks(rotation=45)
    plt.savefig('%s_mutation_counts.svg' % filename, dpi=300, bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    main()

