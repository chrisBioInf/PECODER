#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 14:24:50 2024

@author: christopher
"""


from Bio.PDB import PDBParser 

import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import itertools

from plot_domain_interactions import aggregate_all_tables


cofactor_family = 'COA'


def min_zero(iterable):
    if len(iterable) < 1:
        return 0
    else:
        return min(iterable)


def cofactor_domain_type_count(df: pd.DataFrame):
    data = {
        'Atom': [],
        '# Domain classes': [],
        }
    
    for atom in df['atom_cofactor'].unique():
        df_ = df[df['atom_cofactor'] == atom]
        data['Atom'].append(atom)
        
        df_atom_domains = df_.groupby(by=['t_id']).agg({'t_id': 'first',
                                                        'number': 'count',
                                                        'distance': 'mean',
                                                        'architecture': 'first',
                                                        'ecod_t_name': 'first',
                                                        })
        data['# Domain classes'].append(len(df_atom_domains))
        print('Atom group', atom)
        if len(df_atom_domains) > 0:
           print(df_atom_domains)
    
    sns.barplot(data=data, y='Atom', x='# Domain classes', edgecolor='k', color='cornflowerblue',)
    plt.show()
    
    
def calculate_distance_profile(atoms_dict: dict):
    # Calculate distance measures for this ligand-domain binding situation
    # The case of |A| = 1 is treated as A=0 right now...
    A = min_zero([abs(a1-a2) for (a1, a2) in itertools.combinations(atoms_dict['domain_1_partners'], r=2)])
    B = min_zero([abs(b1-b2) for (b1, b2) in itertools.combinations(atoms_dict['domain_2_partners'], r=2)])
    C = min_zero([abs(a-b) for (a, b) in itertools.product(atoms_dict['domain_1_partners'], atoms_dict['domain_2_partners'])])
    
    return A, B, C
    
    
    
def group_by_domains(df: pd.DataFrame):
    
    results_dict = {
        'pdb_id': [],
        'ligand_uuid': [],
        'Domain 1': [],
        'Domain 2': [],
        'A': [],
        'B': [],
        'C': [],
        }
    
    for ligand in df['ligand_uuid'].unique():
        df_ = df[df['ligand_uuid'] == ligand]
        df_ = df_.dropna()
        
        # Only considering the 2-domain case for now
        if len(df_['ecod_domain_id'].unique()) != 2:
            continue
        pdb_id = df_['pdb'].iloc[0]

        # Aggregate by involved domains
        df_atom_domains = df_.groupby(by=['atom_cofactor', 'ecod_domain_id'], group_keys=False).agg({'atom_cofactor':'first',   
                                                                                                     'ecod_domain_id': 'first',
                                                                                                     'distance': 'min',
                                                                                                     'ligand_uuid': 'first',
                                                                                                     'ecod_t_name': 'first',
                                                                                                     't_id': 'first'})
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
                local_domain = df_atom['ecod_domain_id'].iloc[0]
                atoms_dict[local_domain].append(atom)
        
        # filter out binding atoms from the Cofactor structure
        parser = PDBParser()
        print("Loading structure", pdb_id)
        try:
            structure = parser.get_structure(pdb_id, 'PDB/CoA/%s.pdb' % pdb_id)
        except Exception:
            print("Failed. Skipping...")
            continue
        
        for chain in structure[0].get_list():
            for res in chain.get_list():
                if not 'COA' in res.get_id()[0]:
                    continue
                
                for atom in res.get_list():
                    if atom.get_id() in atoms_dict[domain_1_name]:
                        atoms_dict['domain_1_partners'].append(atom)
                    elif atom.get_id() in atoms_dict[domain_2_name]:
                        atoms_dict['domain_2_partners'].append(atom)
        
        A, B, C = calculate_distance_profile(atoms_dict)
        
        # Pen down the results for this ligand
        results_dict['pdb_id'].append(pdb_id)
        results_dict['ligand_uuid'].append(ligand)
        results_dict['Domain 1'].append(df_atom_domains[df_atom_domains['ecod_domain_id'] == domain_1_name]['t_id'].iloc[0])
        results_dict['Domain 2'].append(df_atom_domains[df_atom_domains['ecod_domain_id'] == domain_2_name]['t_id'].iloc[0])
        
        for (key, val) in zip(('A', 'B', 'C'), (A, B, C)):
            results_dict[key].append(val)
        
    result_table = pd.DataFrame(data=results_dict)
    result_table.to_csv('Atoms/CoA/Distance_metrics.tsv', sep='\t', index=False)
        
        
def main():
    directory = sys.argv[1]
    table = aggregate_all_tables(directory)
    group_by_domains(table)
    

if __name__ == '__main__':
    main()
