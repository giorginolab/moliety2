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

# -----------------------------
# Interligand Moieties Definitions
# -----------------------------
def load_interligand_moieties():
    moieties = {}
    try:
        with open("SMARTS_InteLigand.txt", "r") as f:
            for line in f:
                line = line.strip()
                # Skip blank lines or comments
                if not line or line.startswith("#"):
                    continue
                # Expect lines of the form "MoietyName: SMARTS"
                if ":" in line:
                    parts = line.split(":", 1)
                    name = parts[0].strip()
                    pattern = parts[1].strip()
                    # Optionally, only keep lines that look like a SMARTS (e.g. start with '[')
                    if pattern and pattern.startswith("["):
                        moieties[name] = pattern
    except Exception as e:
        print("Error loading SMARTS_InteLigand.txt:", e)
    return moieties

interligand_moieties = load_interligand_moieties()
compiled_interligand_patterns = {name: Chem.MolFromSmarts(smart) 
                                 for name, smart in interligand_moieties.items()}

# -----------------------------
# Highlighting Functions
# -----------------------------
def highlight_functional_groups(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    images = []
    for group_name, pattern in compiled_patterns.items():
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
    images = highlight_functional_groups(smiles)
    if images is None or len(images) == 0:
        return [], "No functional groups recognized or invalid SMILES."
    return images, f"Found {len(images)} functional group(s)."

def highlight_interligand_moieties(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    images = []
    for moiety_name, pattern in compiled_interligand_patterns.items():
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
            img = Draw.MolToImage(
                mol,
                size=(300, 300),
                highlightAtoms=list(highlight_atoms),
                highlightBonds=list(highlight_bonds),
                legend=moiety_name
            )
            images.append(img)
    return images

def process_interligand_moieties(smiles: str):
    images = highlight_interligand_moieties(smiles)
    if images is None or len(images) == 0:
        return [], "No interligand moieties recognized or invalid SMILES."
    return images, f"Found {len(images)} interligand moiety(ies)."

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
    img = Draw.MolToImage(
        mol,
        size=(300, 300),
        highlightBonds=rot_bonds,
        legend="Rotatable Bonds"
    )
    return img

def process_rotatable(smiles: str):
    img = highlight_rotatable_bonds(smiles)
    if img is None:
        return None, "No rotatable bonds recognized or invalid SMILES."
    return img, "Rotatable bonds highlighted."

def process_smiles_mode(smiles: str, mode: str):
    if mode == "Functional Groups":
        images, status_msg = process_functional_groups(smiles)
        return images, status_msg
    elif mode == "Rotatable Bonds":
        img, status_msg = process_rotatable(smiles)
        if img is None:
            return [], status_msg
        return [img], status_msg
    elif mode == "Interligand moieties":
        images, status_msg = process_interligand_moieties(smiles)
        return images, status_msg
    else:
        return [], "Invalid mode selected."

# -----------------------------
# Gradio Interface
# -----------------------------
with gr.Blocks() as demo:
    gr.Markdown("# Molecule Feature Highlighter")
    gr.Markdown(
        "Enter a SMILES string below and select a highlighting mode. "
        "You can choose to highlight known functional groups, interligand moieties, or display the rotatable bonds in the molecule."
    )
    
    with gr.Row():
        smiles_input = gr.Textbox(
            label="Enter SMILES string", 
            placeholder="e.g. CC(=O)OC1=CC=CC=C1C(=O)O",
            lines=1
        )
        mode_dropdown = gr.Dropdown(
            label="Highlight Mode",
            choices=["Functional Groups", "Rotatable Bonds", "Interligand moieties"],
            value="Functional Groups"
        )
    
    gallery = gr.Gallery(label="Highlighted Features", columns=3, height="auto")
    status = gr.Textbox(label="Status", interactive=False)

    smiles_input.submit(process_smiles_mode, inputs=[smiles_input, mode_dropdown], outputs=[gallery, status])
    mode_dropdown.change(process_smiles_mode, inputs=[smiles_input, mode_dropdown], outputs=[gallery, status])

demo.launch()
