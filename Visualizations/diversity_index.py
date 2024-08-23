#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 16:18:22 2024

@author: christopher
"""


import networkx as nx
from itertools import product


q = 2


def inverse_simpson_index(proportions: list) -> float:
    if (len(proportions) == 0) or set(proportions) == set([0]):
        return 0.0
    
    simpson_index = sum([p_i**q for p_i in proportions])**(1/(1-q))
    return round(simpson_index, 2)


def overlap_index(G, fuzzy_cutoff: int, domain_1_atoms: list, domain_2_atoms: list) -> float:
    if (len(domain_1_atoms) == 0) or (len(domain_2_atoms) == 0):
        return 0.0

    if (fuzzy_cutoff == 0):
        shared_atoms = set(domain_1_atoms).intersection(set(domain_2_atoms))
            
    else:
        node_dict = {a[1]['name'] : a[0] for a in G.nodes(data=True)}
        shared_atoms = []
        for atom1, atom2 in product(domain_1_atoms, domain_2_atoms):
            if len(nx.shortest_path(G, source=node_dict.get(atom1), target=node_dict.get(atom2))) <= fuzzy_cutoff:
                shared_atoms.append(atom1)
                shared_atoms.append(atom2)
                
    return round(len(shared_atoms) / (len(domain_1_atoms) + len(domain_2_atoms)), 2), shared_atoms

