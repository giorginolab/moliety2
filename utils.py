# utils.py
from rdkit import Chem
from rdkit.Chem.Draw import MolDraw2DSVG
import tempfile

IMAGE_SIZE = (400, 400)

def mol_to_svg(mol, size, highlightAtoms=None, highlightBonds=None, legend="", atomLabels=None):
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

def highlight_by_patterns(smiles: str, pattern_dict: dict):
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
