#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 11:16:25 2024

@author: christopher
"""


architecture_type = {
       'alpha arrays': 'a', 
       'alpha bundles': 'a',
       'alpha duplicates or obligate multimers': 'a', 
       'alpha superhelices': 'a',
       'alpha complex topology': 'a', 
       'beta barrels': 'b',
       'beta meanders': 'b', 
       'beta sandwiches': 'b',
       'beta duplicates or obligate multimers': 'b', 
       'beta complex topology': 'b',
       'a+b two layers': 'a+b', 
       'a+b three layers': 'a+b',
       'a+b four layers': 'a+b', 
       'a+b complex topology': 'a+b',
       'a+b duplicates or obligate multimers': 'a+b', 
       'a/b barrels': 'a/b',
       'a/b three-layered sandwiches': 'a/b', 
       'mixed a+b and a/b': 'a+b,a/b',
       'few secondary structure elements': '-', 
       'extended segments': '-',
       }


IUPAC_alphabet = {
    'ALA': 'A',
    'ARG': 'R',
    'ASN': 'N',
    'ASP': 'D',
    'CYS': 'C',
    'GLN': 'Q',
    'GLU': 'E',
    'GLY': 'G',
    'HIS': 'H',
    'ILE': 'I',
    'LEU': 'L',
    'LYS': 'K',
    'MET': 'M',
    'PHE': 'F',
    'PRO': 'P',
    'SER': 'S',
    'THR': 'T',
    'TRP': 'T',
    'TYR': 'Y',
    'VAL': 'V',
    'SEC': 'U',
    'PYL': 'O',
    }


