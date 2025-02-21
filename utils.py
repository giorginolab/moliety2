"""
Utility functions for molecular visualization and pattern matching using RDKit.
Provides SVG rendering capabilities and substructure highlighting functionality.
"""

from rdkit import Chem
from rdkit.Chem.Draw import MolDraw2DSVG
import tempfile

IMAGE_SIZE = (400, 400)

def mol_to_svg(mol: Chem.Mol, size: tuple[int, int], highlightAtoms: list[int] = None, 
               highlightBonds: list[int] = None, legend: str = "", atomLabels: dict[int, str] = None) -> str:
    """
    Converts an RDKit molecule to an SVG image file with optional highlighting and labels.

    Args:
        mol: RDKit molecule object to render
        size: Tuple of (width, height) for the output image
        highlightAtoms: List of atom indices to highlight
        highlightBonds: List of bond indices to highlight
        legend: Text to display as image legend
        atomLabels: Dictionary mapping atom indices to label strings

    Returns:
        str: Path to the temporary SVG file containing the rendered molecule
    """
    drawer = MolDraw2DSVG(size[0], size[1])
    opts = drawer.drawOptions()
    
    if atomLabels:
        mol = Chem.Mol(mol)
        for atom_idx, label in atomLabels.items():
            atom = mol.GetAtomWithIdx(atom_idx)
            atom.SetProp('atomNote', label)
    
    drawer.DrawMolecule(mol, highlightAtoms=highlightAtoms, highlightBonds=highlightBonds, legend=legend)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    tmp = tempfile.NamedTemporaryFile(suffix=".svg", delete=False)
    tmp.write(svg.encode("utf-8"))
    tmp.close()
    return tmp.name

def highlight_by_patterns(smiles: str, pattern_dict: dict[str, Chem.Mol]) -> list[tuple[str, str]]:
    """
    Generates visualizations of a molecule with substructures matching SMARTS patterns highlighted.

    Args:
        smiles: SMILES string representation of the molecule
        pattern_dict: Dictionary mapping pattern names to RDKit molecule objects representing SMARTS patterns

    Returns:
        list: List of tuples containing (image_path, pattern_name) for each matching pattern.
              Returns None if the SMILES string is invalid.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    images = []
    for name, pattern in pattern_dict.items():
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            highlight_atoms = set()
            highlight_bonds = set()
            for match in matches:
                highlight_atoms.update(match)
                for bond in mol.GetBonds():
                    a1 = bond.GetBeginAtomIdx()
                    a2 = bond.GetEndAtomIdx()
                    if a1 in match and a2 in match:
                        highlight_bonds.add(bond.GetIdx())
            img = mol_to_svg(mol, IMAGE_SIZE, highlightAtoms=list(highlight_atoms), highlightBonds=list(highlight_bonds), legend=name)
            images.append((img, name))
    return images
