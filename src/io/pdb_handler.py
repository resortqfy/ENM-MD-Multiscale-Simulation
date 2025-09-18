import MDAnalysis as mda
from Bio.PDB import PDBParser
from typing import Dict

def parse_pdb(pdb_file: str) -> Dict:
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)
    # 进一步处理
    return {'structure': structure}

def mda_load(pdb_file: str) -> mda.Universe:
    return mda.Universe(pdb_file)

