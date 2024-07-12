from Bio import PDB

def extract_residues(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('antibody', pdb_file)
    residue_coords = []
    residue_ids = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() in ['TYR', 'TRP']:
                    for atom in residue:
                        residue_coords.append(atom.get_coord())
                        residue_ids.append(f"{residue.get_resname()}_{residue.id[1]}")

    return residue_coords, residue_ids
