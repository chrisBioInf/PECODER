#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 16:18:22 2024

@author: christopher
"""


q = 2


def inverse_simpson_index(proportions: list) -> float:
    if (len(proportions) == 0) or set(proportions) == set([0]):
        return 0.0
    
    simpson_index = sum([p_i**q for p_i in proportions])**(1/(1-q))
    return simpson_index

