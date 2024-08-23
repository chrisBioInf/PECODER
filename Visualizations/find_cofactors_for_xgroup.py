#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 14:36:47 2024

@author: christopher
"""


import sys
import os
import glob
import pandas as pd


xgroup = sys.argv[1]

dfs = []
cofactors = ['NAD', 'COA', 'FAD', 'SAH', 'ATP']

for c in cofactors:
    directory = os.path.join('Results', c)
    cofactor = directory.split('/')[-1]
    fs = glob.glob(os.path.join(directory, '*.tsv'))
    
    for f in fs:
        if len(f.split('_')) > 1:
            continue
        df = pd.read_csv(f, sep='\t')
        
        if 'x-group' not in df.columns or str(xgroup) not in [str(x) for x in df['x-group']]:
            continue
        
        df['Cofactor'] = cofactor
        dfs.append(df)
        
df = pd.concat(dfs)
df.to_csv('%s_data.tsv' % xgroup, sep='\t', index=False)

