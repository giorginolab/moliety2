"""Compatibility imports for older callers."""

from rdkit import Chem

from moliety2.rendering import IMAGE_SIZE
from moliety2.rendering import highlight_by_patterns as _highlight_by_patterns
from moliety2.rendering import mol_to_svg as _mol_to_svg


def mol_to_svg(
    mol: Chem.Mol,
    size: tuple[int, int] = IMAGE_SIZE,
    highlightAtoms: list[int] | None = None,
    highlightBonds: list[int] | None = None,
    legend: str = "",
    atomLabels: dict[int, str] | None = None,
) -> str:
    return _mol_to_svg(
        mol,
        size=size,
        highlight_atoms=highlightAtoms,
        highlight_bonds=highlightBonds,
        legend=legend,
        atom_labels=atomLabels,
    )


def highlight_by_patterns(smiles: str, pattern_dict: dict[str, Chem.Mol]) -> list[tuple[str, str]] | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return _highlight_by_patterns(mol, pattern_dict)
