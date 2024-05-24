#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 11:29:26 2024

@author: christopher
"""


"""
Helper script that downloads all PDB structures associated with one of
a list of cofactors. 

Arguments:
    -f Flatfile from PDBsum categorizing PDB ids by associated ligands
    -o Output directory (must exist)
    
The required mapping flatfile can be obtained via

> curl https://www.ebi.ac.uk/thornton-srv/databases/pdbsum/data/lig2pdb.lst > lig2pdb.txt

"""

__version__ = 1.0


from parse_PDB_flatfile import get_ids_for_cofactor

from optparse import OptionParser
import sys
import os
import urllib3
import shutil


http = urllib3.PoolManager()
BASE_URL="https://files.rcsb.org/download"


def download_id(pdb_id: str, outdir: str, filetype: str, verbose: bool) -> int:                                                                                                                                                                                                         
    url = "{base_url}/{id}.pdb.gz".format(
            base_url=BASE_URL,
            id=pdb_id,
            )
    outfile = os.path.join(outdir, "%s.%s.gz" % (pdb_id, filetype))
    
    if verbose:
        print(url)
        
    try:
        with http.request('GET', url, preload_content=False) as r, open(outfile, 'wb') as out_file:        
            shutil.copyfileobj(r, out_file)
        return 0
    
    except Exception:
        print("Failed to download: %s" % url)
        return 1


def download_structures_with_cofactor(cofactor: str, flatfile: str, outdir: str, filetype: str, verbose: bool):
    pdb_id_list = get_ids_for_cofactor(cofactor, flatfile)
    n_ids = len(pdb_id_list)
    return_codes = []
    
    for i, idx in enumerate(pdb_id_list):
        if verbose:
            print("Now downloading (%s/%s):" % (i, n_ids))
        
        ret_ = download_id(pdb_id=idx, outdir=outdir, filetype=filetype, verbose=verbose)
        return_codes.append(ret_)
    
    if verbose:
        print("Successfully downloaded: %s/%s" % (n_ids - sum(return_codes), n_ids))


def main():
    usage = "\download_from_PDB.py -f [ligand flatfile] [Cofactor1] [Cofactor2] ..."
    parser = OptionParser(usage=usage, version="__version__")
    args = sys.argv
    
    parser.add_option("-f", "--flatfile", action="store", default="", type="string", dest="flatfile", help="Ligand -> PDB mapping flatfile, as obtained from PDBsum (Required).")
    parser.add_option("-o", "--outdir", action="store", default="", type="string", dest="outdir", help="Output directory to download files to. Default is cwd.")
    parser.add_option("-t", "--filetype", action="store", default="pdb", type="string", dest="filetype", help="File type to download from PDB (i.e. pdb, cif, etc.). Default: pdb.")
    parser.add_option("-v", "--verbose", action="store_false", default=True, dest="verbose", help="Display download messages. Default: True.")
    
    options, args = parser.parse_args()
    
    if len(options.flatfile) == 0:
        print("Flatfile with Ligand->PDB mapping is required. Exiting...")
        sys.exit()
    
    FLATFILE = options.flatfile
    
    if len(options.outdir) == 0:
        OUTDIR = str(os.getcwd())
    else:
        OUTDIR = options.outdir
        if not os.path.isdir(OUTDIR):
            print("Error: Output directory does not exist.")
            sys.exit()
            
    for argument in args:
        if options.verbose:
            print("Cofactor is: %s" % argument)
        
        download_structures_with_cofactor(cofactor=argument, flatfile=FLATFILE, 
                                          outdir=OUTDIR, filetype=options.filetype,
                                          verbose=options.verbose,
                                          )


if __name__ == "__main__":
    main()

