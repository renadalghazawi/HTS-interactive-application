from rdkit import Chem
from rdkit.Chem import Draw

def smiles_to_png(smiles, output_file):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        img = Draw.MolToImage(mol)
        img.save(output_file)
    else:
        raise ValueError("Invalid SMILES string")

