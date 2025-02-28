"""
Module for identifying and visualizing rotatable bonds in molecular structures.
Provides functionality to detect single bonds that can rotate using both
a local algorithm and SMARTS pattern matching approaches.
"""

from rdkit import Chem
from utils import mol_to_svg, highlight_by_patterns, IMAGE_SIZE

# Rotatable bond patterns
rotatable_patterns = {
    "DAYLIGHT defn.": Chem.MolFromSmarts("[!$(*#*)&!D1]-!@[!$(*#*)&!D1]"),
    "RDKit defn.": Chem.MolFromSmarts("[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]"),
}


def get_rotatable_bond_indices(mol):
    """
    Identifies rotatable bonds in a molecule using local structural analysis.

    A bond is considered rotatable if it:
    - Is a single bond
    - Is not in a ring
    - Neither atom is a hydrogen
    - Both atoms have at least 2 neighbors

    Args:
        mol: RDKit molecule object

    Returns:
        list: Indices of rotatable bonds in the molecule
    """
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
    """
    Creates an SVG visualization of a molecule with rotatable bonds highlighted.

    Args:
        smiles: SMILES string representation of the molecule

    Returns:
        str: SVG string of the molecule with rotatable bonds highlighted, or
             None if no rotatable bonds are found or if SMILES is invalid
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    rot_bonds = get_rotatable_bond_indices(mol)
    if not rot_bonds:
        return None
    img = mol_to_svg(
        mol, IMAGE_SIZE, highlightBonds=rot_bonds, legend="Rotatable Bonds"
    )
    return img


def process_rotatable(smiles: str):
    """
    Processes a molecule to identify rotatable bonds using multiple methods
    and generates visualizations for each method.

    Args:
        smiles: SMILES string representation of the molecule

    Returns:
        tuple: (list of (image, caption) tuples, status message string)
               Images show the molecule with rotatable bonds highlighted using
               different detection methods
    """
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
