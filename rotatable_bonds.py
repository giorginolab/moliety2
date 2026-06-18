"""Compatibility wrappers for rotatable bond helpers."""

from rdkit import Chem

from moliety2.features import get_rotatable_bond_indices, rotatable_bonds
from moliety2.rendering import mol_to_svg


def highlight_rotatable_bonds(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    bond_indices = get_rotatable_bond_indices(mol)
    if not bond_indices:
        return None

    return mol_to_svg(mol, highlight_bonds=bond_indices, legend="Rotatable Bonds")


def process_rotatable(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [], "Invalid SMILES."
    return rotatable_bonds(mol, smiles, 7.0, 7.0)
