#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 13:40:11 2024

@author: christopher
"""


import matplotlib.pyplot as plt
import json

import networkx as nx
from networkx.readwrite import json_graph
import pubchempy as pcp


color_dict = {
    'H': 'blue',
    'O': 'red',
    'C': 'lightblue',
    'S': 'yellow',
    'N': 'plum',
    'P': 'coral',
    }

def generate_reference(cid: int, atommap_file: str, outfile: str):
    compound = pcp.Compound.from_cid(cid)
    atommap = json.load(open(atommap_file, 'r'))
    G = nx.Graph()
    
    for atom in compound.atoms:
        if (atom.element == 'H'):
            continue
        atom_name = str().join([atom.element, str(atom.aid)])
        G.add_node(atom.aid, name=atommap['Atoms'].get(atom_name), 
                   moiety=atommap['Moieties'].get(atommap['Atoms'].get(atom_name)),
                   color=color_dict.get(atom.element, 'white'))

    for bond in compound.bonds:
        if (bond.aid1 not in G.nodes) or (bond.aid2 not in G.nodes):
            continue
        G.add_edge(bond.aid1, bond.aid2)
        
    json_representation = json_graph.node_link_data(G)
    json.dump(obj=json_representation, fp=open(outfile, 'w'))
    
    nx.draw(G, with_labels=True, pos=nx.kamada_kawai_layout(G), 
            labels={n[0]: n[1]['name'] for n in G.nodes(data=True)},
            node_color=[n[1]['color'] for n in G.nodes(data=True)])
    plt.savefig(outfile.replace('.json', '.svg'), dpi=300)
    plt.show()

