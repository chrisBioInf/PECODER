#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 16:18:22 2024

@author: christopher
"""


import itertools
import numpy as np

fractions = [
    (1.0,),
    (0.75, 0.25),
    (0.25, 0.25, 0.25, 0.25),
    (0.5, 0.167, 0.167, 0.167,)
    ]


for interaction_fraction in fractions:
    if len(interaction_fraction) == 1:
        diversity_set = [1.0]
    else:
        diversity_set = [abs(a-b) for a, b in itertools.combinations(interaction_fraction, 2)]
    print(interaction_fraction, diversity_set, max(diversity_set))