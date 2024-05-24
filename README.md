# PECODER
Pipeline for analyzing Domain-Cofactor interaction patterns using PDB and ECOD


# Basic Usage

> git clone https://github.com/chrisBioInf/PECODER/
> cd PECODER

(You will need, at a minimum, a functioning Pandas and Biopython installation in your Python PATH)


1.) Get the ECOD reference data base:

> cd ECOD
> sh download_ecod_tsv.sh 291
> cd ..

(you can replace the 291 with any version number you want)

2.) Get a Cofactor -> PDB ID mapping from PDBsum:

> curl https://www.ebi.ac.uk/thornton-srv/databases/pdbsum/data/lig2pdb.lst > lig2pdb.txt

3.) Download PDB crystal structures corresponding to some Cofactor (here: COA):

> python download_from_PDB.py -f lig2pdb.txt -o PDB COA

4.) Scan and annotate atom-level interactions with PECODER:

mkdir -p Results
> python PECODER.py --ecod ECOD/ecod_domains291.tsv --pdb PDB -o Results -c COA
