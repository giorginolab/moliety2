from rdkit import Chem
from utils import mol_to_svg, highlight_by_patterns, IMAGE_SIZE

# Rotatable bond patterns
rotatable_patterns = {
    "DAYLIGHT defn.": Chem.MolFromSmarts("[!$(*#*)&!D1]-!@[!$(*#*)&!D1]"),
    "RDKit defn.":    Chem.MolFromSmarts('[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]')
}

def get_rotatable_bond_indices(mol):
    rot_bond_indices = []
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        if bond.IsInRing():
            continue
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 1 or a2.GetAtomicNum() == 1:
            continue
        if a1.GetDegree() < 2 or a2.GetDegree() < 2:
            continue
        rot_bond_indices.append(bond.GetIdx())
    return rot_bond_indices

def highlight_rotatable_bonds(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    rot_bonds = get_rotatable_bond_indices(mol)
    if not rot_bonds:
        return None
    img = mol_to_svg(mol, IMAGE_SIZE, highlightBonds=rot_bonds, legend="Rotatable Bonds")
    return img

def process_rotatable(smiles: str):
    images = []
    img = highlight_rotatable_bonds(smiles)
    if img:
        images.append((img, "Local algorithm"))
    img = highlight_by_patterns(smiles, rotatable_patterns)
    images.extend(img)
    if images:
        return images, "Rotatable bonds highlighted."
    else:
        return [], "No rotatable bonds found"
