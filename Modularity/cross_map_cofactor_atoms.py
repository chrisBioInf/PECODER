#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 10:06:09 2024

@author: christopher
"""


import sys
import glob
import json
import networkx as nx
from random import shuffle
import pandas as pd

import matplotlib.pyplot as plt


cofactors = ['COA', 'NAD', 'ATP', 'SAH', 'FAD']

color_dict = {
    'H': 'blue',
    'O': 'red',
    'C': 'lightblue',
    'S': 'yellow',
    'N': 'plum',
    'P': 'coral',
    }

def read_json_file(filename):
    with open(filename) as f:
        js_graph = json.load(f)
    return nx.readwrite.json_graph.node_link_graph(js_graph)


def get_subgraph_from_nodes(G, nodes, moiety):
    H = G.__class__()
    number = 0
    
    for node in sorted(nodes):
        number += 1
        atom = node[1]['name'][0]
        H.add_node(node[0], label='%s_%s%s' % (moiety, atom, number), atom_type=atom, color=color_dict.get(atom))
        
    H.add_edges_from(
            (n, nbr, d)
            for n, nbrs in G.adj.items()
            if n in [na[0] for na in nodes]
            for nbr, d in nbrs.items()
            if nbr in [na[0] for na in nodes]
    )
    H.graph.update(G.graph)
    nx.draw(H, )
    return H


def build_all_moiety_graphs():
    mapped_moieties = set()
    
    for f in glob.glob('Compounds/*.json'):
        if len(f.split('_')) > 1 or 'atommap' in f:
            continue
        G = read_json_file(f)
        cofactor = f.split('/')[-1].strip('.json')
        reference_map = json.load(open('Compounds/%s_atommap.json' % cofactor, 'r'))
        moiety_atoms = {}
        
        for n in G.nodes(data=True):
            n[1]['moiety'] = reference_map['Moieties'].get(n[1]['name'])
            if n[1]['moiety'] in mapped_moieties:
                continue
            if n[1]['moiety'] not in moiety_atoms:
                moiety_atoms[n[1]['moiety']] = [n]
            else:
                moiety_atoms[n[1]['moiety']].append(n)
                
        for moiety, nodes in moiety_atoms.items():
            H = get_subgraph_from_nodes(G, nodes, moiety)
            nx.draw(H, with_labels=True, pos=nx.kamada_kawai_layout(H), 
                labels={n[0]: n[1]['label'].split('_')[-1] for n in H.nodes(data=True)},
                node_color=[n[1]['color'] for n in H.nodes(data=True)])
            plt.show()
            
            mapped_moieties.add(moiety)
            json_representation = nx.readwrite.json_graph.node_link_data(H)
            json.dump(obj=json_representation, fp=open('Compounds/Moieties/%s.json' % moiety, 'w'))


def align_moiety_to_reference(G, nodes, H, moiety):
    def check_neighborhood_viability(neighbors1, neighbors2):   
        neighbor_atoms = [n['name'][0] for n in neighbors1 if n['moiety'] == moiety]
        for n2 in neighbors2:
            if n2['atom_type'] not in neighbor_atoms:
                return False
        return True
    
    solved = False
    
    while not solved:
        shuffle(nodes)
        M = []
        already_mapped = set()
        # current match
        
        for n in nodes:
            for node_ref in H.nodes(data=True):
                if (not n[1]['name'][0] == node_ref[1]['atom_type']) or (node_ref[0] in already_mapped):
                    continue
                if not check_neighborhood_viability([G.nodes[n_] for n_ in G.neighbors(n[0])], [H.nodes[n_] for n_ in H.neighbors(node_ref[0])]):
                    continue
                M.append((n, node_ref))
                already_mapped.add(node_ref[0])
                break
        
        # CHECK EDGE VIABILITY
        solved = True
        ordered_nodes_G = [n1[0] for n1,n2 in M]
        ordered_nodes_H = [n2[0] for n1,n2 in M]

        for i in range(0, len(ordered_nodes_G)-1):
            for j in range(i+1, len(ordered_nodes_G)):
                if (G.has_edge(ordered_nodes_G[i], ordered_nodes_G[j])) and (not H.has_edge(ordered_nodes_H[i], ordered_nodes_H[j])):
                    solved = False      
    
    for n1,n2 in M:
        print('%s -> %s' % (n1[1]['name'],n2[1]['label']))
        
    return M


def main():
    # BUILD
    # build_all_moiety_graphs()
    # sys.exit()
    
    # TEST
    for cofactor in cofactors:
        G = read_json_file('Compounds/%s.json' % cofactor)
        reference_map = json.load(open('Compounds/%s_atommap.json' % cofactor, 'r'))
        moiety_atoms = {}
        
        for n in G.nodes(data=True):
            n[1]['moiety'] = reference_map['Moieties'].get(n[1]['name'])
            if n[1]['moiety'] not in moiety_atoms:
                moiety_atoms[n[1]['moiety']] = [n]
            else:
                moiety_atoms[n[1]['moiety']].append(n)
                
        atom_map = {
            'Moiety': [],
            'Atom': [],
            'Reference atom': [],
            }
                
        for moiety, nodes in moiety_atoms.items():
            H = read_json_file('Compounds/Moieties/%s.json' % moiety)
            
            if moiety == 'Riboflavin':
                M = []
                for n in nodes:
                    for node_ref in H.nodes(data=True):
                        if not n[0] == node_ref[0]:
                            continue
                        M.append((n, node_ref))
                for n1,n2 in M:
                    print('%s -> %s' % (n1[1]['name'],n2[1]['label']))
            else:
                M = align_moiety_to_reference(G, nodes, H, moiety)
            
            for n1,n2 in M:
                atom_map['Moiety'].append(moiety),
                atom_map['Atom'].append(n1[1]['name']), 
                atom_map['Reference atom'].append(n2[1]['label'])
        
        df = pd.DataFrame(data=atom_map)
        df.to_csv('Compounds/%s_map.tsv' % cofactor, sep='\t', index=False)

main()
