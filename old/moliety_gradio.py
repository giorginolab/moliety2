import gradio as gr
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdFMCS
from io import BytesIO
import base64
import PIL.Image

# Known moieties SMARTS patterns
known_moieties = {
    'amide': '[NX3][CX3](=[OX1])[#6]',
    'ester': '[#6][CX3](=[OX1])[OX2H0][#6]',
    'carboxylic_acid': '[CX3](=O)[OX2H1]',
    'alcohol': '[OX2H]',
    'aromatic_ring': '[c]'
}

def smiles_to_image(smiles_code, highlight_atoms):
    try:
        mol = Chem.MolFromSmiles(smiles_code)
        if mol is None:
            return None
        img = Draw.MolToImage(mol, highlightAtoms=highlight_atoms)
        return img
    except Exception as e:
        return None

def highlight_moieties(smiles_code):
    mol = Chem.MolFromSmiles(smiles_code)
    if mol is None:
        return "Invalid SMILES"
    
    highlight_atoms = []
    for name, smarts in known_moieties.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            highlight_atoms.extend(match)
    
    return smiles_to_image(smiles_code, highlight_atoms)

gui = gr.Interface(
    fn=highlight_moieties,
    inputs=gr.Textbox(label="SMILES Code"),
    outputs=gr.Image(label="Highlighted Moieties"),
    title="Moiety Highlighter",
    description="Enter a SMILES string to highlight known moieties."
)

gui.launch()
