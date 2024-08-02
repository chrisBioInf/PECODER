#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 12:31:54 2024

@author: christopher
"""


"""
Current main functionality script, to be developed further. It basically implements
the following steps: 
    1.) Get a list of all corresponding cofactors for a cofactor familiy (i.e. CoA, etc.))
    2.) Parse a list of PDB crystal structures and locate the given cofactor.
    3.) Calculate, based on spatial distance, a list of all possibly interacting 
        atoms and the corresponding AA residues they belong to.
    4.) Commit these predicted interaction partners to a table and match them 
        with the corresponding ECOD protein doman annotation(s).
    5.) Write the resulting pairings as a table documenting all interacting atoms,
        residues and their domain classifications. 
"""


from Bio.PDB import PDBParser 
from networkx.readwrite import json_graph
import json

import matplotlib.pyplot as plt
import networkx as nx

import os
import sys
import numpy as np
import pandas as pd
import glob

from lookup_tables import aa_chain_backbone, architecture_type, IUPAC_alphabet
from graph_interface import make_domain_interaction_graph
from diversity_index import inverse_simpson_index

THIS_FILE = os.path.abspath(__file__)
THIS_DIR = os.path.dirname(__file__)
ATOM_MAP_DIR  = os.path.join(THIS_DIR, 'Compounds')

if ATOM_MAP_DIR not in sys.path:
    sys.path.append(ATOM_MAP_DIR)



# Hacky workaround to deal with certain numpy versions
np.int = np.int32
np.float = np.float64
np.bool = np.bool_

# Maximum Angstrom distance for interaction
distance_threshold = 4

# Defining of Ecod table columns
ecod_columns = ['pdb', 'chain', 'number', 'pdb_range', 't_id', 'ecod_domain_id',
                'ecod_a_name', 'ecod_x_name', 'ecod_h_name', 'ecod_t_name',
                'ecod_f_name', 'ecod_ligand']

alphabet = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']


def explode_pdb_range(df: pd.DataFrame, connection) -> pd.DataFrame:
    # In the ECOD data base, domains are locally defined as ranges. For easier 
    # db access and matching later, we explode ranges into rows with one residue 
    # position each.
    
    expanded_blocks = []
    
    for idx, row in df.iterrows():
        pdb_range = row['pdb_range']
        ls = pdb_range.split(',')
        rows = []
        
        for entry in ls:
            chain_, range_str = entry.split(':')
            
            try:
                if range_str.startswith('-'):
                    range_str = "0-%s" % range_str.split('-')[-1]

                range_ = range(int(range_str.split('-')[0]), int(range_str.split('-')[1]))
                
            except Exception:
                print('No pdb range given. Skipping... (%s)' % range_str)
                continue
            
            for n in range_:
                row_data = {
                    'pdb': [row['pdb']],
                    'chain': [chain_],
                    'number': [n],
                    'pdb_range': [row['pdb_range']],
                    't_id': [row['t_id']],
                    'ecod_domain_id': [row['ecod_domain_id']],
                    'ecod_a_name': [row['arch_name']],
                    'ecod_x_name': [row['x_name']],
                    'ecod_h_name': [row['h_name']],
                    'ecod_t_name': [row['t_name']],
                    'ecod_f_name': [row['f_name']],
                    'ecod_ligand': [row['ligand']],
                    }
                df_row = pd.DataFrame(data=row_data)
                rows.append(df_row)
            
            if len(rows) > 0:
                expanded_blocks.append(pd.concat(rows))
    
    if len(expanded_blocks) > 0:
        df_ = pd.concat(expanded_blocks)
        df_.index = [n for n in range(0, len(df_))]
    
    else:
        df_ = pd.DataFrame(
            data = {
                k:[] for k in ecod_columns
                }
            )
    df_['pdb'] = df_['pdb'].astype(str)
    df_['chain'] = df_['chain'].astype(str)
    df_['number'] = df_['number'].astype(int)
    return df_


def fetch_cofactor_family(name: str) -> pd.DataFrame:
    ligand_file = os.path.join(THIS_DIR, "Ligands", "ligand_list.tsv")
    df_ligand = pd.read_csv(ligand_file, sep='\t', )
    df_family = df_ligand[df_ligand['Cofactor_Family'] == name]
    return df_family


def spatial_distance(atom1, atom2) -> float:
    return abs(atom1 - atom2)


def get_coordinate_set(structure) -> (np.array, list):
    coord_vectors = []
    res_ids = []
    
    for chain in structure[0].get_list():
        for res in chain.get_list():
            for atom in res.get_list():
                v = atom.get_coord()
                id_res = res.get_id()[1]
                id_atom = atom.get_id()
                coord_vectors.append(v)
                res_ids.append((id_res, res.get_resname(), id_atom))
    
    a = np.asarray(coord_vectors, dtype=np.float64)
    
    return a, res_ids


def calculate_aa_distances(residue1, residue2) -> (list, list, list):
    distances = []
    atoms_res1 = []
    atoms_res2 = []
    
    for atomA in residue1.get_list():
        for atomB in residue2.get_list():
            if (atomB.get_id() in aa_chain_backbone):
                continue
            distances.append(spatial_distance(atomA, atomB))
            atoms_res1.append(atomA)
            atoms_res2.append(atomB)
        
    return distances, atoms_res1, atoms_res2


def get_target_residues(structure, target_list: list, cofactor: str, cofactor_atoms: list) -> (list, list):
    targets = []
    target_ids = []
    
    # Identify residues belonging to target cofactor
    for chain in structure[0].get_list():
        target_instance = 1
        for res in chain.get_list():
            if (res.get_id()[0] in target_list):
                target_id = '%s_%s_%s' % (chain.get_id(), res.get_id()[0], target_instance)
                
                while (target_id in target_ids):
                    target_instance += 1
                    target_id = '%s_%s_%s' % (chain.get_id(), res.get_id()[0], target_instance)
                
                targets.append(res)
                target_ids.append(target_id)
    
    complete_targets = []
    complete_target_ids = []
    
    # Only add chemically -complete- cofactors to analysis
    for i in range(0, len(targets)):
        cofactor_is_complete = True
        detected_atoms = [a.get_id() for a in targets[i].get_list()]
        print(sorted(detected_atoms))
        print(sorted(cofactor_atoms))
        for atom in cofactor_atoms:
            if not atom in detected_atoms:
                cofactor_is_complete = False
                print("Cofactor is incomplete!")
                break
        if cofactor_is_complete:
            complete_targets.append(targets[i])
            complete_target_ids.append(target_ids[i])

    return complete_targets, complete_target_ids


def calculate_domain_proportions(df):
    df['ecod_domain_id'] = df['ecod_domain_id'].fillna(value='NoName')
    count_dict = {id_: df['ecod_domain_id'].value_counts()[id_] for id_ in df['ecod_domain_id'].unique()}
    sum_of_interactions = sum([v for k,v in count_dict.items()])
    proportions_dict = {id_: round(count_dict.get(id_) / sum_of_interactions, 2) for id_ in df['ecod_domain_id'].unique()}
    return [proportions_dict.get(id_) for id_ in df['ecod_domain_id']], [val for val in proportions_dict.values()]


def annote_binding_mode(df, binding_mode_cutoff):
    df['x-group'] = [str(x).split('.')[0] for x in df['t_id']]
    dfs = []
    
    for ligand in df['ligand_uuid'].unique():
        df_ = df[df['ligand_uuid'] == ligand]
        df_['binding domains'] = len(df_['ecod_domain_id'].unique())
        df_['domain proportions'], absolute_proportions = calculate_domain_proportions(df_)
        df_['inverse simpson'] = inverse_simpson_index(absolute_proportions)
        df_['absolute position'] = ['%s_%s' % (row['chain'], row['number']) for idx, row in df_.iterrows()]
        
        if (df_['binding domains'].iloc[0] < 2) or (df_['inverse simpson'].iloc[0] <= binding_mode_cutoff):
            df_['binding mode'] = 'A'
            df_['domain configuration'] = df_.sort_values(by='domain proportions', ascending=False)['x-group'].iloc[0]
            dfs.append(df_)
            continue
        
        binding_mode_str = []
        domain_configuration_str = []
        last_xgroup = ''
        last_chain = ''
        domain_alphabet = {}
        domain_type_N = 0
        df_x = df_.groupby(by=['ecod_domain_id']).agg({'chain': 'first',
                                                       'absolute position': 'first',
                                                       'x-group': 'first',
                                                       'ecod_domain_id': 'first'})
        df_x.sort_values(by='absolute position', inplace=True)
        
        for i, row in df_x.iterrows():
            if last_xgroup == '':
                binding_mode_str.append('A')
                domain_configuration_str.append(row['x-group'])
                domain_alphabet[row['x-group']] = 'A'
                domain_type_N += 1
            elif last_xgroup == row['x-group']:
                if last_chain == row['chain']:
                    binding_mode_str.append('-')
                    domain_configuration_str.append('-')
                else:
                    binding_mode_str.append('|')
                    domain_configuration_str.append('|')
                binding_mode_str.append(domain_alphabet.get(row['x-group']))
                domain_configuration_str.append(row['x-group'])
            else:
                if last_chain == row['chain']:
                    binding_mode_str.append('-')
                    domain_configuration_str.append('-')
                else:
                    binding_mode_str.append('|')
                    domain_configuration_str.append('|')
                domain_alphabet[row['x-group']] = alphabet[domain_type_N]
                binding_mode_str.append(domain_alphabet.get(row['x-group']))
                domain_configuration_str.append(row['x-group'])
                domain_type_N += 1
                
            last_xgroup = row['x-group']
            last_chain = row['chain']
        
        df_['binding mode'] = str().join(binding_mode_str)
        df_['domain configuration'] = str().join(domain_configuration_str)
        dfs.append(df_)
    
    return pd.concat(dfs)


def export_dataframe(filename: str, data: dict) -> pd.DataFrame:
    accepted_amino_acids = set((aa for aa in IUPAC_alphabet.keys()))
    df = pd.DataFrame(data=data)
    df = df[df["residue"].isin(accepted_amino_acids)]
    df['pdb'] = df['pdb'].astype(str)
    df['chain'] = df['chain'].astype(str)
    df['number'] = df['number'].astype(int)
    df.sort_values(by=['chain', 'number'], inplace=True)
    return df 


def assign_architecture_labels(df: pd.DataFrame) -> pd.DataFrame:
    df['architecture'] = [architecture_type.get(a, '-') for a in df['ecod_a_name']]
    return df


def annotate_cofactor_interactions(ecod_file: str, reference_graph: str, pdb_dir: str, cofactor_family: str, 
                                   fuzzy_cutoff: int, binding_mode_cutoff: float, resolution_cutoff: float,
                                   outdir: str):
    parser = PDBParser()
    
    # Load ECOD reference annotation
    df_ecod = pd.read_csv(ecod_file, sep='\t', skiprows=4)
    
    # Load cofactor map
    cofactor_dict = json.load(open(os.path.join(ATOM_MAP_DIR, '%s_atommap.json' % cofactor_family), 'r'))
    cofactor_moiety_map = cofactor_dict.get('Moieties')
    cofactor_atoms = list(cofactor_moiety_map.keys())
    
    # Load compound reference graph
    with open(reference_graph, 'r') as f:
        js_graph = json.load(f)
    
    G = json_graph.node_link_graph(js_graph)
    df_cofactor_family = fetch_cofactor_family(cofactor_family)
    print(df_cofactor_family)
    cofactor_list = list(df_cofactor_family['Cofactor'])
    cofactor_list += ['H_%s' % lig for lig in cofactor_list]

    for i, filename in enumerate(glob.glob(os.path.join(pdb_dir, '*.pdb'))):
        pdb_id = os.path.basename(filename).strip('.pdb').lower()
        ecod_ = explode_pdb_range(df_ecod[df_ecod["pdb"] == pdb_id], connection=None)
        print("Loading structure: %s" % filename)
        
        try:
            structure = parser.get_structure(filename.split('/')[1].split('.')[0], filename)
            targets, target_ids = get_target_residues(structure, target_list=cofactor_list, 
                                                      cofactor=cofactor_family, cofactor_atoms=cofactor_atoms)
        except Exception:
            print("Loading of structure failed for: %s. Skipping..." % filename)
            continue 
        
        # Check for minimal resolution
        if structure.header["resolution"] == None:
            print('No annotated resolution. Skipping...')
            continue
        
        if structure.header["resolution"] > resolution_cutoff:
            print('Insufficient resolution. Skipping...')
            continue
        
        data = {
            'pdb': [],
            'chain': [],
            'number': [],
            'residue': [],
            'distance': [],
            'cofactor': [],
            'ligand_id': [],
            'ligand_uuid': [],
            'atom_cofactor': [],
            'atom_residue': [],
            'moiety': [],
            }
        for target, idx in zip(targets, target_ids):
            id_ = target.get_id()
            print("\nCofactor: %s (res: %s)" % (id_[0], id_[1]))
            print("%s %s" % (target, idx))
            print("Interacting residues:")
            
            for chain in structure[0].get_list():
                for res in chain.get_list():
                    if (res.get_id()[0] in cofactor_list):
                        continue
                    
                    distances, atoms_target, atoms_aa = calculate_aa_distances(target, res)
                    
                    if (len(distances) == 0) or (min(distances) > distance_threshold):
                        continue
                    
                    for distance, atomA, atomB in zip(distances, atoms_target, atoms_aa):
                        if (distance > distance_threshold):
                            continue
                        
                        data['pdb'].append(pdb_id)
                        data['chain'].append(chain.get_id())
                        data['residue'].append(res.get_resname())
                        data['number'].append(res.get_id()[1])
                        data['distance'].append(distance)
                        data['cofactor'].append(id_[0])
                        data['ligand_id'].append(idx)
                        data['ligand_uuid'].append('%s_%s' % (pdb_id, idx))
                        data['atom_cofactor'].append(atomA.get_id())
                        data['atom_residue'].append(atomB.get_id())
                        data['moiety'].append(cofactor_moiety_map.get(atomA.get_id()))
                        print("%s (res: %s): \t%s" % ((res.get_resname(), res.get_id()[1], min(distances))))
    
        df = export_dataframe(filename, data)
        df_merged = df.merge(ecod_, on=["pdb", "chain", "number"], how='left')
        df_merged = assign_architecture_labels(df_merged)
        
        if len(df_merged) != 0:
            df_merged = annote_binding_mode(df_merged, binding_mode_cutoff)
        
        outfile = os.path.join(outdir, '%s.tsv' % pdb_id)
        df_merged.to_csv(outfile, sep='\t', index=None)
        
        # Calculate the interaction graph:
        make_domain_interaction_graph(df_merged, G, structure, fuzzy_cutoff, outdir)

