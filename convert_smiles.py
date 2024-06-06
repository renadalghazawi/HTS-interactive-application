from rdkit import Chem
from rdkit.Chem import Draw

def smiles_to_png(smiles, output_file):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    img = Draw.MolToImage(mol)
    img.save(output_file)

