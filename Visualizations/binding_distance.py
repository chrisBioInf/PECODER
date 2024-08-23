#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 17:12:59 2024

@author: christopher
"""


def jaccard_distance(v1: list, v2: list):
    intersection = 0
    union = 0

    for n1,n2 in zip(v1, v2):
        if (n1 == -1 or n2 == -1):
            continue
        if (n1 == 0 and n2 == 0):
            continue
        if (n1 == n2):
            intersection += 1
        union += 1
    
    return intersection / union

    