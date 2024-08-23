#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 12:36:44 2024

@author: christopher
"""


import os
import sys
import glob
import json
import numpy
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import networkx as nx
from networkx.readwrite import json_graph
from xgroup_interactions import group_by_domains
from evolution_events import load_data, map_domain_changes
from diversity_index import inverse_simpson_index


distance_cutoff = 3.6


def read_json_file(filename):
    with open(filename) as f:
        js_graph = json.load(f)
    return json_graph.node_link_graph(js_graph)


def main():
    data = {
        'Cofactor': [],
        'X-group configuration': [],
        'Inverse simpson': [],
        'Betweenness centrality': [],
        'Degree centrality': [],
        'Eigenvector centrality': [],
        }
    for cofactor in sys.argv[1:]:
        G = read_json_file(cofactor + '_configuration_graph.json')
        df = pd.read_csv('%s_xgroup_diversities.tsv' % cofactor, sep='\t')
        
        for node in G.nodes(data=True):
            xconfig = node[1]['label']
            df_ = df[df['X-group configuration'] == xconfig]
            if len(df_) == 0:
                node[1]['inverse_simpson'] = 0.0
                continue
            node[1]['inverse_simpson'] = df_['Inverse Simpson'].iloc[0]
            data['Cofactor'].append(cofactor)
            data['X-group configuration'].append(xconfig)
            data['Inverse simpson'].append(node[1]['inverse_simpson'])
            data['Betweenness centrality'].append(node[1]['betweenness_centrality'])
            data['Degree centrality'].append(node[1]['degree_centrality'])
            data['Eigenvector centrality'].append(node[1]['eigenvector_centrality'])
            
        nx.draw(G, with_labels=True, pos=nx.kamada_kawai_layout(G),
                edgecolors=['k']*len(G.nodes),
                cmap=plt.get_cmap('magma'), 
                node_color=[n[1]['inverse_simpson'] for n in G.nodes(data=True)])
        plt.show()
        
        data = pd.DataFrame(data=data)
        data.sort_values(by='Inverse simpson', inplace=True)
        
        fig, axes = plt.subplots(2, 2, figsize=(10, 10))
        sns.scatterplot(data=data, x='Inverse simpson', y='Degree centrality', 
                        hue='Cofactor', edgecolor='k', ax=axes[0][0])
        sns.scatterplot(data=data, x='Inverse simpson', y='Betweenness centrality', 
                        hue='Cofactor', edgecolor='k', ax=axes[0][1])
        sns.scatterplot(data=data, x='Inverse simpson', y='Eigenvector centrality', 
                        hue='Cofactor', edgecolor='k', ax=axes[1][0])
        
        plt.show()

main()

