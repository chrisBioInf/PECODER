#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 12:27:47 2024

@author: christopher
"""


import os
import json
import itertools
import pandas as pd
import networkx as nx
from networkx.readwrite import json_graph
import matplotlib.pyplot as plt

from diversity_index import inverse_simpson_index, overlap_index


distance_cutoff = 4.0


def min_zero(iterable) -> float:
    if len(iterable) < 1:
        return 0.0
    else:
        return min(iterable)
    

def max_zero(iterable) -> float:
    if len(iterable) < 1:
        return 0.0
    else:
        return max(iterable)
    
    
def draw_graph_with_interaction_atoms(G: nx.graph, domain_names: list, domain_1_atoms: list, domain_2_atoms: list, shared_atoms: list, 
                                      pdb_id: str, path_profile_str: str, filename: str):
    domain_1_color = 'skyblue'
    domain_2_color = 'darkcyan'
    both_domain_color = 'coral'
    colors = []
    edgecolors = []
    sizes = []
    widths = []
    
    for n in G.nodes(data=True):
        if n[1]['distance'] == 0:
            edgecolors.append('k')
            sizes.append(200)
            widths.append(0.8)
        elif n[1]['distance'] <= 3.6:
            edgecolors.append('k')
            sizes.append(250)
            widths.append(2)
        else:
            edgecolors.append('k')
            sizes.append(200)
            widths.append(1)
        
        if (n[0] in shared_atoms):
            colors.append(both_domain_color)
        elif (n[0] in domain_1_atoms):
            colors.append(domain_1_color)
        elif (n[0] in domain_2_atoms):
            colors.append(domain_2_color)
        else:
            colors.append('white')
    
    nx.draw(G, with_labels=True, pos=nx.kamada_kawai_layout(G),
            node_color=colors, edgecolors=edgecolors, node_size=sizes, 
            linewidths=widths)
    plt.text(x=0.5, y=0.25, s=path_profile_str, fontsize=12)
    plt.savefig(filename +'.svg', dpi=300, bbox_inches='tight')
    plt.clf()
    json_representation = json_graph.node_link_data(G)
    json.dump(obj=json_representation, fp=open(filename + '.json', 'w'))


def calculate_shortest_path_profile(G: nx.graph, graph_file: str, atoms_dict: dict, domain_names: list, shared_atoms: list, domain1: str, 
                                    domain2: str,  pdb_id: str, inverse_simpson: float, overlap: float, binding_mode: str) -> (float, float, float, float):
    node_dict = {a[1]['name'] : a[0] for a in G.nodes(data=True)}
    domain_1_atoms = [node_dict.get(a) for a in atoms_dict[domain1] if a in node_dict.keys()]
    domain_2_atoms = [node_dict.get(a) for a in atoms_dict[domain2] if a in node_dict.keys()]
    shared_atoms = [node_dict.get(a) for a in shared_atoms]
    domain_1_atompairs = list(itertools.combinations(domain_1_atoms, r=2))
    domain_2_atompairs = list(itertools.combinations(domain_2_atoms, r=2))

    domain_1_inter_paths = [len(nx.shortest_path(G, source=a, target=b))-1 for (a, b) in domain_1_atompairs]
    A = max_zero(domain_1_inter_paths) - min_zero(domain_1_inter_paths)
    
    domain_2_inter_paths = [len(nx.shortest_path(G, source=a, target=b))-1 for (a, b) in domain_2_atompairs]
    B = max_zero(domain_2_inter_paths) - min_zero(domain_2_inter_paths)
    
    domain_connecting_paths = [len(nx.shortest_path(G, source=a, target=b))-1 
                               for (a, b) in itertools.product(domain_1_atoms, domain_2_atoms)]
    C_1 = min_zero(domain_connecting_paths)
    C_2 = max_zero(domain_connecting_paths)
    
    path_profile_str = 'Domain 1: %s \nDomain 2: %s \nBinding mode: %s \nIntra domain range A:  %s \nIntra domain range B:  %s \nMin. inter domain distance: %s \nMax. inter domain distance: %s \nAdjusted Simpson: %s \nOverlap: %s' % (domain_names[0], domain_names[1], binding_mode, A, B, C_1, C_2, inverse_simpson, overlap)
    
    draw_graph_with_interaction_atoms(G, domain_names, domain_1_atoms, domain_2_atoms, shared_atoms, 
                                      pdb_id, path_profile_str, graph_file)
    return A, B, C_1, C_2


def calculate_diversity_scores(G: nx.graph, domain_1_atoms: list, domain_2_atoms: list, fuzzy_cutoff: int) -> (float, float, list):
    n_atoms_1, n_atoms_2 = len(domain_1_atoms), len(domain_2_atoms)
    domain_proportions = [n_atoms_1 /(n_atoms_1 + n_atoms_2), n_atoms_2 /(n_atoms_1 + n_atoms_2)]
    inverse_simpson = inverse_simpson_index(domain_proportions)
    overlap, shared_atoms = overlap_index(G, fuzzy_cutoff, domain_1_atoms, domain_2_atoms)
    return inverse_simpson, overlap, shared_atoms


def make_domain_interaction_graph(df: pd.DataFrame, G: nx.graph, pdb_structure, fuzzy_cutoff: int, outfile: str):
    df = df[df['distance'] <= distance_cutoff] 
    atom_name_dict = {
        n[0] : n[1]['name'] for n in G.nodes(data=True)
        }
    
    results_dict = {
        'pdb_id': [],
        'ligand_uuid': [],
        'Domain 1': [],
        'Domain 2': [],
        'domain 1 atoms': [],
        'domain 2 atoms': [],
        'Inverse Simpson': [],
        'Overlap index': [],
        'Intra domain range A': [],
        'Intra domain range B': [],
        'Minimum inter domain distance': [],
        'Maximum inter domain distance': [],
        }
    for ligand in df['ligand_uuid'].unique():
        df_ = df[df['ligand_uuid'] == ligand]
        df_ = df_.dropna()
        
        # Only considering the 2-domain case for now
        if len(df_['ecod_domain_id'].unique()) != 2:
            continue
        pdb_id = df_['pdb'].iloc[0]
        
        if len(df_) == 0:
            continue
        cofactor = df_['cofactor'].iloc[0]
        binding_mode = df['binding mode'].iloc[0]

        # Aggregate by involved domains
        df_atom_domains = df_.groupby(by=['atom_cofactor', 'ecod_domain_id'], group_keys=False).agg({'atom_cofactor':'first',   
                                                                                                     'ecod_domain_id': 'first',
                                                                                                     'distance': 'min',
                                                                                                     'ligand_uuid': 'first',
                                                                                                     'ecod_t_name': 'first',
                                                                                                     't_id': 'first'})
        # Assign distance values to each node
        for node in G.nodes():
            atom_name = atom_name_dict.get(node)
            df_atom = df_[df_['atom_cofactor'] == atom_name]
            
            if len(df_atom) == 0:
                distance = 0.0
            else:
                distance = float(df_atom['distance'].min())
            
            G.nodes[node]['distance'] = distance
        
        # Set domain names - domain 1 has >= interactions than domain 2
        domain_1_name, domain_2_name = df_atom_domains.ecod_domain_id.value_counts().index
        
        atoms_dict = {
            domain_1_name: [],
            domain_2_name: [],
            'domain_1_partners': [],
            'domain_2_partners': [],
            } 
        # Group interacting atoms by domain1/domain2
        for atom in df_atom_domains['atom_cofactor'].unique():
            df_atom = df_atom_domains[df_atom_domains['atom_cofactor'] == atom]
            
            if len(df_atom) == 1:
                local_domain = df_atom['ecod_domain_id'].iloc[0]
                atoms_dict[local_domain].append(atom)
            else:
                df_atom.sort_values(by='distance', inplace=True, ascending=True)
                closer_domain = df_atom['ecod_domain_id'].iloc[0]
                atoms_dict[closer_domain].append(atom)
                secondary_domain = df_atom['ecod_domain_id'].iloc[1]
                atoms_dict[secondary_domain].append(atom)
        
        for chain in pdb_structure[0].get_list():
            for res in chain.get_list():
                if not cofactor in res.get_id()[0]:
                    continue
                for atom in res.get_list():
                    if atom.get_id() in atoms_dict[domain_1_name]:
                        atoms_dict['domain_1_partners'].append(atom)
                    elif atom.get_id() in atoms_dict[domain_2_name]:
                        atoms_dict['domain_2_partners'].append(atom)
        
        # Name for Graph to write as file
        graph_file = os.path.join(outfile, '%s' % (ligand))
        
        # Pen down the results for this ligand
        results_dict['pdb_id'].append(pdb_id)
        results_dict['ligand_uuid'].append(ligand)
        results_dict['Domain 1'].append(df_atom_domains[df_atom_domains['ecod_domain_id'] == domain_1_name]['t_id'].iloc[0])
        results_dict['Domain 2'].append(df_atom_domains[df_atom_domains['ecod_domain_id'] == domain_2_name]['t_id'].iloc[0])
        results_dict['domain 1 atoms'].append(len(atoms_dict[domain_1_name]))
        results_dict['domain 2 atoms'].append(len(atoms_dict[domain_2_name]))
        
        domain_names = [results_dict['Domain 1'][-1], results_dict['Domain 2'][-1]]
        motif_id = '%s_%s_%s' % (pdb_id, domain_names[0], domain_names[1])
        inverse_simpson, overlap, shared_atoms = calculate_diversity_scores(G, atoms_dict[domain_1_name], atoms_dict[domain_2_name], fuzzy_cutoff)
        A, B, C_1, C_2 = calculate_shortest_path_profile(G, graph_file, atoms_dict, domain_names, shared_atoms, domain_1_name, domain_2_name, motif_id, inverse_simpson, overlap, binding_mode)
        
        for (key, val) in zip(('Inverse Simpson', 'Overlap index', 'Intra domain range A', 'Intra domain range B', 'Minimum inter domain distance', 'Maximum inter domain distance'), (inverse_simpson, overlap, A, B, C_1, C_2)):
            results_dict[key].append(val)
            
        data = pd.DataFrame(data=results_dict)
        data.to_csv(os.path.join(outfile, '%s.tsv' % (ligand)), sep='\t')
            
       