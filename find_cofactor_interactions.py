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
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.kdtrees import KDTree
from Bio.PDB.NeighborSearch import NeighborSearch

import os
import numpy as np
import pandas as pd
import glob

import sqlite3
from sqlite_interface import create_table_if_not_existing, add_rows_to_table
from lookup_tables import aa_chain_backbone, architecture_type, IUPAC_alphabet


THIS_FILE = os.path.abspath(__file__)
THIS_DIR = os.path.dirname(__file__)


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
    # cursor  = add_rows_to_table(connection, 'ecod', df_)
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


def get_target_residues(structure, target_list: list) -> (list, list):
    targets = []
    target_ids = []
    
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
    
    return targets, target_ids


def export_dataframe(filename: str, chains: list, res_name: list, residue_n: list, 
                     distances: list, interaction_partner: list, atoms_cofactor: list, 
                     atoms_residue: list, family: str, target_id: str) -> pd.DataFrame:
    pdb_id = os.path.basename(filename).strip('.pdb').lower()
    data = {
        'pdb': pdb_id,
        'chain': chains,
        'number': residue_n,
        'residue': res_name,
        'distance': distances,
        'cofactor': interaction_partner,
        'family': family,
        'ligand_id': target_id,
        'atom_cofactor': atoms_cofactor,
        'atom_residue': atoms_residue,
        }
    accepted_amino_acids = set((aa for aa in IUPAC_alphabet.keys()))
    
    df = pd.DataFrame(data=data)
    df = df[df["residue"].isin(accepted_amino_acids)]
    df['pdb'] = df['pdb'].astype(str)
    df['chain'] = df['chain'].astype(str)
    df['number'] = df['number'].astype(int)
    df.sort_values(by=['chain', 'number'], inplace=True)
    return df 


def assign_analysis_labels(df: pd.DataFrame) -> pd.DataFrame:
    df['architecture'] = [architecture_type.get(a, '-') for a in df['ecod_a_name']]
    return df


def annotate_cofactor_interactions(ecod_file: str, pdb_dir: str, cofactor_family: str, outdir: str):
    parser = PDBParser()
    # connection = sqlite3.connect('Cofactor_domains.db')
    # last_row_id = create_table_if_not_existing(connection)
    
    df_ecod = pd.read_csv(ecod_file, sep='\t', skiprows=4)
    
    df_cofactor_family = fetch_cofactor_family(cofactor_family)
    cofactor_list = list(df_cofactor_family['Cofactor'])

    for i, filename in enumerate(glob.glob(os.path.join(pdb_dir, '*.pdb'))):
        pdb_id = os.path.basename(filename).strip('.pdb').lower()
        ecod_ = explode_pdb_range(df_ecod[df_ecod["pdb"] == pdb_id], connection=None)

        print("Loading structure: %s" % filename)
        
        try:
            structure = parser.get_structure(filename.split('/')[1].split('.')[0], filename)
            targets, target_ids = get_target_residues(structure, target_list=cofactor_list)
        except Exception:
            print("Loading of structure failed for: %s. Skipping..." % filename)
            continue
        
        chains = []
        res_name = []
        residue_n = []
        atom_distances = []
        interaction_partner = []
        cofactor_id = []
        atoms_cofactor = []
        atoms_residue = []
        
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
                        
                        chains.append(chain.get_id())
                        res_name.append(res.get_resname())
                        residue_n.append(res.get_id()[1])
                        atom_distances.append(distance)
                        interaction_partner.append(id_[0])
                        cofactor_id.append(idx)
                        atoms_cofactor.append(atomA.get_id())
                        atoms_residue.append(atomB.get_id())
                        
                    print("%s (res: %s): \t%s" % (res_name[-1], residue_n[-1], min(distances)))
    
    
        df = export_dataframe(filename, chains, res_name, residue_n, atom_distances, 
                              interaction_partner, atoms_cofactor, atoms_residue,
                              df_cofactor_family['Cofactor_Family'].iloc[0], cofactor_id)
        df_merged = df.merge(ecod_, on=["pdb", "chain", "number"], how='left')
        df_merged = assign_analysis_labels(df_merged)
        outfile = os.path.join(outdir, '%s.tsv' % pdb_id)
        df_merged.to_csv(outfile, sep='\t', index=None)
        # cursor = add_rows_to_table(connection, table='interactions', df=df_merged)
        # print(df_merged)
    
    # connection.close()


