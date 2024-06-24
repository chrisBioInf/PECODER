#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 11:43:01 2024

@author: christopher
"""


import pandas as pd
import matplotlib.pyplot as plt


def draw_anchor_domains(name: str, anchor_domains: dict):
    fig, ax = plt.subplots()
    labels, sizes = [], []
    
    for key, val in anchor_domains.items():
        labels.append(key)
        sizes.append(val)

    ax.pie(sizes, labels=labels, autopct='%1.1f%%',
           pctdistance=1.25, labeldistance=.6,
           colors=['white', 'cyan', 'mediumorchid', 'lightcoral', 'olivedrab', 'plum'],
           wedgeprops={'linewidth': 1, 'linestyle': 'solid', 'antialiased': True,
                       'edgecolor':'k'})
    
    if len(name) > 0:
        plt.savefig(name.split('.')[0] + '_anchor_pie.svg', dpi=300)
    
