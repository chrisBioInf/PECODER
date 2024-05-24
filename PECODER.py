#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 13:39:16 2024

@author: christopher
"""

"""
Main Script. Calls the annotation pipeline located in 
find_cofactor_interactions.py

Arguments:
    --ecod  ECOD database as downloaded table (tab-delimited)
    --pdb   Directory with .pdb crystal structures
    -c      Name of the target Cofactor family
    -o      Ouput directory

"""

__version__ = 0.5


import os
import sys
from optparse import OptionParser

from find_cofactor_interactions import annotate_cofactor_interactions


def main():
    usage = "\PECODER.py -ecod [ecod db] -pdb [directory] -c [Cofactor] "
    parser = OptionParser(usage=usage, version="__version__")
    args = sys.argv
    
    parser.add_option("--ecod", action="store", default="", type="string", dest="ecod", help="Downloaded ECOD tab-delimited file (Required).")
    parser.add_option("--pdb", action="store", default="", type="string", dest="pdb_dir", help="Directory with .pdb crystal structure files (Required).")
    parser.add_option("-o", "--outdir", action="store", default="", type="string", dest="outdir", help="Output directory to download files to. Default is cwd.")
    parser.add_option("-c", "--cofactor", action="store", default="", type="string", dest="cofactor", help="Cofactor family to screen for.")
    # parser.add_option("-v", "--verbose", action="store_false", default=True, dest="verbose", help="Display download messages. Default: True.")
    
    options, args = parser.parse_args()
    
    if len(options.ecod) == 0:
        print("ECOD table (--ecod) is required. Exiting...")
        sys.exit()
        
    ECOD = options.ecod
        
    if len(options.pdb_dir) == 0:
        print("Directory with .pdb crystal structures (--pdb) required. Exiting...")
        sys.exit()
    
    elif not os.path.isdir(options.pdb_dir):
        print("Argument provided for --pdb is no directory. Exiting...")
        sys.exit()
        
    PDB_DIR = options.pdb_dir
    
    if len(options.outdir) == 0:
        OUTDIR = str(os.getcwd())
    else:
        OUTDIR = options.outdir
        if not os.path.isdir(OUTDIR):
            print("Error: Output directory does not exist.")
            sys.exit()
            
    if len(options.cofactor) == 0:
        print("Cofactor family name is required. Exiting...")
        sys.exit()
    
    # Call main script
    annotate_cofactor_interactions(ecod_file=ECOD, pdb_dir=PDB_DIR, 
                                   cofactor_family=options.cofactor, 
                                   outdir=OUTDIR)


if __name__ == "__main__":
    main()
