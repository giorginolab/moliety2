"""RDKit rendering and highlighting helpers."""

from __future__ import annotations

import tempfile

from rdkit import Chem
from rdkit.Chem.Draw import MolDraw2DSVG

IMAGE_SIZE = (400, 400)


def mol_to_svg(
    mol: Chem.Mol,
    size: tuple[int, int] = IMAGE_SIZE,
    highlight_atoms: list[int] | None = None,
    highlight_bonds: list[int] | None = None,
    legend: str = "",
    atom_labels: dict[int, str] | None = None,
) -> str:
    """Render an RDKit molecule to a temporary SVG file for Gradio."""
    drawable_mol = Chem.Mol(mol) if atom_labels else mol
    if atom_labels:
        for atom_idx, label in atom_labels.items():
            drawable_mol.GetAtomWithIdx(atom_idx).SetProp("atomNote", label)

    drawer = MolDraw2DSVG(size[0], size[1])
    drawer.DrawMolecule(
        drawable_mol,
        highlightAtoms=highlight_atoms,
        highlightBonds=highlight_bonds,
        legend=legend,
    )
    drawer.FinishDrawing()

    tmp = tempfile.NamedTemporaryFile(suffix=".svg", mode="w", encoding="utf-8", delete=False)
    tmp.write(drawer.GetDrawingText())
    tmp.close()
    return tmp.name


def matching_bond_indices(mol: Chem.Mol, atom_indices: set[int] | tuple[int, ...]) -> list[int]:
    """Return bonds whose two endpoint atoms are in the provided atom set."""
    atoms = set(atom_indices)
    return [
        bond.GetIdx()
        for bond in mol.GetBonds()
        if bond.GetBeginAtomIdx() in atoms and bond.GetEndAtomIdx() in atoms
    ]


def highlight_by_patterns(
    mol: Chem.Mol,
    pattern_dict: dict[str, Chem.Mol | None],
) -> list[tuple[str, str]]:
    """Generate one highlighted image for each matching SMARTS pattern."""
    images: list[tuple[str, str]] = []
    for name, pattern in pattern_dict.items():
        if pattern is None:
            continue

        matches = mol.GetSubstructMatches(pattern)
        if not matches:
            continue

        highlight_atoms: set[int] = set()
        highlight_bonds: set[int] = set()
        for match in matches:
            highlight_atoms.update(match)
            highlight_bonds.update(matching_bond_indices(mol, match))

        images.append(
            (
                mol_to_svg(
                    mol,
                    highlight_atoms=sorted(highlight_atoms),
                    highlight_bonds=sorted(highlight_bonds),
                    legend=name,
                ),
                name,
            )
        )
    return images
