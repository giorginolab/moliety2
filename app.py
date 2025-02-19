import gradio as gr
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

# -----------------------------
# Functional Group Definitions
# -----------------------------
functional_groups = {
    "Hydroxyl": "[OX2H]",
    "Carbonyl": "[CX3]=[OX1]",
    "Amine": "[NX3;H2,H1;!$(NC=O)]",
    "Carboxylic Acid": "C(=O)[OX2H1]",
    "Ester": "C(=O)O[C]",
    "Ether": "[OD2]([#6])[#6]",
    "Halide": "[F,Cl,Br,I]"
}

# Pre-compile the SMARTS queries for functional groups.
compiled_patterns = {name: Chem.MolFromSmarts(smart) for name, smart in functional_groups.items()}

def highlight_functional_groups(smiles: str):
    """
    Given a SMILES string, find all matching functional groups and return
    a list of images (as PIL Images) of the molecule with each group highlighted.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None  # invalid SMILES

    images = []
    # Loop over each functional group
    for group_name, pattern in compiled_patterns.items():
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            # Collect all atoms and bonds that match this group.
            highlight_atoms = set()
            highlight_bonds = set()
            for match in matches:
                highlight_atoms.update(match)
                # For bonds, check each bond in the molecule
                for bond in mol.GetBonds():
                    a1 = bond.GetBeginAtomIdx()
                    a2 = bond.GetEndAtomIdx()
                    if a1 in match and a2 in match:
                        highlight_bonds.add(bond.GetIdx())
            # Draw the molecule with the highlighted atoms and bonds,
            # including the group name as a legend.
            img = Draw.MolToImage(
                mol,
                size=(300, 300),
                highlightAtoms=list(highlight_atoms),
                highlightBonds=list(highlight_bonds),
                legend=group_name
            )
            images.append(img)
    return images

def process_functional_groups(smiles: str):
    """
    Process SMILES and return images highlighting functional groups.
    """
    images = highlight_functional_groups(smiles)
    if images is None or len(images) == 0:
        return [], "No functional groups recognized or invalid SMILES."
    return images, f"Found {len(images)} functional group(s)."

# -----------------------------
# Rotatable Bond Functions
# -----------------------------
def get_rotatable_bond_indices(mol):
    """
    Identify rotatable bonds in the molecule and return their bond indices.
    A simple definition is used:
      - Bond is single and not in a ring.
      - Both atoms are heavy atoms (atomic number > 1).
      - Both atoms have degree > 1 (i.e. non-terminal).
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
    Given a SMILES string, highlight all rotatable bonds in the molecule.
    Returns a PIL Image with rotatable bonds highlighted.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    rot_bonds = get_rotatable_bond_indices(mol)
    if not rot_bonds:
        return None
    img = Draw.MolToImage(
        mol,
        size=(300, 300),
        highlightBonds=rot_bonds,
        legend="Rotatable Bonds"
    )
    return img

def process_rotatable(smiles: str):
    """
    Process SMILES and return an image highlighting the rotatable bonds.
    """
    img = highlight_rotatable_bonds(smiles)
    if img is None:
        return None, "No rotatable bonds recognized or invalid SMILES."
    return img, "Rotatable bonds highlighted."

# -----------------------------
# Combined Processing Function
# -----------------------------
def process_smiles_mode(smiles: str, mode: str):
    """
    Depending on the selected mode, process the SMILES string either to highlight:
      - Functional Groups (returns a list of images), or
      - Rotatable Bonds (returns a single image wrapped in a list).
    """
    if mode == "Functional Groups":
        images, status_msg = process_functional_groups(smiles)
        return images, status_msg
    elif mode == "Rotatable Bonds":
        img, status_msg = process_rotatable(smiles)
        if img is None:
            return [], status_msg
        return [img], status_msg
    else:
        return [], "Invalid mode selected."

# -----------------------------
# Gradio Interface
# -----------------------------
with gr.Blocks() as demo:
    gr.Markdown("# Molecule Feature Highlighter")
    gr.Markdown(
        "Enter a SMILES string below and select a highlighting mode. "
        "You can choose to highlight known functional groups or to display the rotatable bonds in the molecule."
    )
    
    with gr.Row():
        smiles_input = gr.Textbox(
            label="Enter SMILES string", 
            placeholder="e.g. CC(=O)OC1=CC=CC=C1C(=O)O",
            lines=1
        )
        mode_dropdown = gr.Dropdown(
            label="Highlight Mode",
            choices=["Functional Groups", "Rotatable Bonds"],
            value="Functional Groups"
        )
    
    gallery = gr.Gallery(label="Highlighted Features", columns=3, height="auto")
    status = gr.Textbox(label="Status", interactive=False)

    # Submit using both SMILES and selected mode.
    smiles_input.submit(process_smiles_mode, inputs=[smiles_input, mode_dropdown], outputs=[gallery, status])
    mode_dropdown.change(process_smiles_mode, inputs=[smiles_input, mode_dropdown], outputs=[gallery, status])

demo.launch()

