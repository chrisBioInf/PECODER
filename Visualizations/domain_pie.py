#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 11:43:01 2024

@author: christopher
"""


import pandas as pd
import matplotlib.pyplot as plt


def draw_anchor_domains(name: str, anchor_domains: pd.DataFrame):
    fig, ax = plt.subplots()
    labels, sizes = [], []
    
    anchor_domains_agg = anchor_domains.groupby(by=['anchor_domain_arch']).agg({'anchor_domain_arch': 'first',
                                                                                'domain_1': 'count'})
    
    for idx, row in anchor_domains_agg.iterrows():
        labels.append(row['anchor_domain_arch'])
        sizes.append(row['domain_1'])

    ax.pie(sizes, labels=labels, autopct='%1.1f%%',
           pctdistance=1.25, labeldistance=.6,
           colors=['white', 'cyan', 'mediumorchid', 'lightcoral', 'olivedrab', 'plum'],
           wedgeprops={'linewidth': 1, 'linestyle': 'solid', 'antialiased': True,
                       'edgecolor':'k'})
    
    if len(name) > 0:
        plt.savefig(name.split('.')[0] + '_anchor_pie.svg', dpi=300)
    
    plt.show()
