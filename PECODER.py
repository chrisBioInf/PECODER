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
from compound_graph import generate_reference


def calculate_interactions(parser):
    parser.add_option("--ecod", action="store", default="", type="string", dest="ecod", help="Downloaded ECOD tab-delimited file (Required).")
    parser.add_option("--pdb", action="store", default="", type="string", dest="pdb_dir", help="Directory with .pdb crystal structure files (Required).")
    parser.add_option("-g", "--graph", action="store", default="", type="string", dest="graph", help="Path to chemical reference graph (Required).")
    parser.add_option("-o", "--outdir", action="store", default="", type="string", dest="outdir", help="Output directory to download files to. Default is cwd.")
    parser.add_option("-c", "--cofactor", action="store", default="", type="string", dest="cofactor", help="Cofactor family to screen for.")
    parser.add_option("-z", "--fuzzy-cutoff", action="store", default=1, type="int", dest="fuzzy_cutoff", help="Distance cutoff for calculation of domain overlap. If 0, only actually identical atoms are considered shatred between atoms. Default: 1.")
    parser.add_option("-b", "--bindingmode-cutoff", action="store", default=1.0, type="float", dest="binding_mode_cutoff", help="Inverse simpson index minimum to consider a shared ligand to be actually in a multi-domain binding configuration. Default: 1 (i.e. no limit).") 
    parser.add_option("-r", "--resolution-cutoff", action="store", default=2.5, type="float", dest="resolution_cutoff", help="Minimal PDB structure resolution to consider for analysis (Default: 2.5).")    
    options, args = parser.parse_args()
    
    if len(options.ecod) == 0:
        print("ECOD table (--ecod) is required. Exiting...")
        sys.exit()
        
    ECOD = options.ecod
    
    if len(options.graph) == 0:
        print("Reference graph (--graph) is required. Exiting...")
        sys.exit()
        
    GRAPH = options.graph
        
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
    annotate_cofactor_interactions(ecod_file=ECOD, 
                                   reference_graph=GRAPH,
                                   pdb_dir=PDB_DIR, 
                                   cofactor_family=options.cofactor, 
                                   fuzzy_cutoff=options.fuzzy_cutoff,
                                   binding_mode_cutoff=options.binding_mode_cutoff,
                                   resolution_cutoff=options.resolution_cutoff,
                                   outdir=OUTDIR)


def make_reference_graph(parser):
    parser.add_option("-o", "--outfile", action="store", default="", type="string", dest="outfile", help="Name of output file.")
    parser.add_option("-i", "--cid", action="store", type="int", dest="cid", help="Pubchem cid to download from (Required).")
    parser.add_option("-a", "--atommap", action="store", default="", type="string", dest="atommap", help="Atom map for naming convention (Required).")
    options, args = parser.parse_args()
    
    if len(options.outfile) == 0:
        OUTFILE = 'graph.json'
    else:
        OUTFILE = options.outfile
        
    if not options.cid:
        print("A punchem CID is required.")
        sys.exit()
    
    CID = options.cid
        
    if len(options.atommap) == 0:
        print("Path to a JSON file containing a map of atom residues (Pubchem -> PDB) is required.")
        sys.exit()
    
    MAP = options.atommap
    generate_reference(CID, MAP, OUTFILE)


def main():
    usage = "PECODER.py [command]" #-ecod [ecod db] -pdb [directory] -c [Cofactor] "
    parser = OptionParser(usage=usage, version="__version__")
    args = sys.argv
    
    if len(args) < 2:
        print("Command (calculate, refgraph...) must be provided.")
        sys.exit()

    if (args[1] == 'calculate'):
        calculate_interactions(parser)
        
    elif (args[1] == 'refgraph'):
        make_reference_graph(parser)
        
    else:
        print("Function not recognized. Exiting...")
        sys.exit()
        

if __name__ == "__main__":
    main()
