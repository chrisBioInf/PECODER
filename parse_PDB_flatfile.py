#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 17:12:09 2024

@author: christopher
"""


"""
Parses a flatfile containing a Cofactor -> PDB id mapping, as obtained
from PDBsum, i.e.

> curl https://www.ebi.ac.uk/thornton-srv/databases/pdbsum/data/lig2pdb.lst > lig2pdb.txt

Returns a list with all matching pdb id instances.
"""


def get_ids_for_cofactor(name: str, flatfile: str) -> list:
    pdb_ids = []
    
    with open(flatfile, 'r') as f:
        found_cofactor = False
        
        for line in f.readlines():
            if (line.startswith(' ') and found_cofactor):
                pdb_ids.append(line.strip().strip('\n'))
            
            elif (line.strip('\n') == name):
                found_cofactor = True
                
            else:
                found_cofactor = False
    
    return pdb_ids

