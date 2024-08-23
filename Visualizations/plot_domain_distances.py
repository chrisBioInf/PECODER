#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 13:32:27 2024

@author: christopher
"""


import os
import sys
import glob
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


dfs = []
cofactors = []

for directory in sys.argv[1:]:
    cofactor = directory.split('/')[-1]
    cofactors.append(cofactor)
    fs = glob.glob(os.path.join(directory, '*_1.tsv'))
    
    for f in fs:
        df = pd.read_csv(f, sep='\t')
        df['Cofactor'] = cofactor
        dfs.append(df)
        
df = pd.concat(dfs)
df.index = list([i for i in range(len(df))])
df['X-group'] = [str(x).split('.')[0] for x in df['Domain 1']]

max_ranges = {}
min_distances = {}
max_distances = {}

for cofactor in cofactors:
    df_ = df[df['Cofactor'] == cofactor]
    max_ranges[cofactor] = df_['Intra domain range A'].max()
    min_distances[cofactor] = df_['Minimum inter domain distance'].min()
    max_distances[cofactor] = df_['Maximum inter domain distance'].max()
    
print(max_ranges)
print(min_distances)
print(max_distances)
     
df['Intra domain range A'] = [x / max_ranges.get(c) for x,c in zip(df['Intra domain range A'], df['Cofactor'])]
df['Minimum inter domain distance'] = [x / max_distances.get(c) for x,c in zip(df['Minimum inter domain distance'], df['Cofactor'])]
df['Maximum inter domain distance'] = [x / max_distances.get(c) for x,c in zip(df['Maximum inter domain distance'], df['Cofactor'])]

fig, ax = plt.subplots(3, 1, figsize=(20, 15))

sns.boxplot(data=df, x='X-group', y='Intra domain range A',
            width=0.5, ax=ax[0])
ax[0].set_xlabel('')
ax[0].set(xticklabels=[])

sns.boxplot(data=df, x='X-group', y='Minimum inter domain distance',
             width=0.5, ax=ax[1])
ax[1].set_xlabel('')
ax[1].set(xticklabels=[])

sns.boxplot(data=df, x='X-group', y='Maximum inter domain distance',
            width=0.5, ax=ax[2])
ax[2].set_xticklabels(ax[2].get_xticklabels(), rotation=90, ha='right')

plt.savefig('domain_distances.png', dpi=600, bbox_inches='tight')
plt.savefig('domain_distances.svg', dpi=600, bbox_inches='tight')
plt.show()

