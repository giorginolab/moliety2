import gradio as gr
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

IMAGE_SIZE=(600,600)

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
compiled_patterns = {name: Chem.MolFromSmarts(smart) 
                     for name, smart in functional_groups.items()}

# -----------------------------
# Interligand Moieties Definitions
# -----------------------------
def load_interligand_moieties():
    moieties = {}
    try:
        with open("SMARTS_InteLigand.txt", "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if ":" in line:
                    parts = line.split(":", 1)
                    name = parts[0].strip()
                    pattern = parts[1].strip()
                    if pattern and pattern.startswith("["):
                        moieties[name] = pattern
    except Exception as e:
        print("Error loading SMARTS_InteLigand.txt:", e)
    return moieties

interligand_moieties = load_interligand_moieties()
compiled_interligand_patterns = {name: Chem.MolFromSmarts(smart)
                                 for name, smart in interligand_moieties.items()}

# -----------------------------
# Generic Highlighting Function
# -----------------------------
def highlight_by_patterns(smiles: str, pattern_dict: dict):
    """
    Generic function to highlight substructures in a molecule given a dictionary
    of name:SMARTS patterns.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    images = []
    for name, pattern in pattern_dict.items():
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
                size=IMAGE_SIZE,
                highlightAtoms=list(highlight_atoms),
                highlightBonds=list(highlight_bonds),
                legend=name
            )
            images.append(img)
    return images

# These two functions simply call the generic one with different pattern dictionaries.
def highlight_functional_groups(smiles: str):
    return highlight_by_patterns(smiles, compiled_patterns)

def highlight_interligand_moieties(smiles: str):
    return highlight_by_patterns(smiles, compiled_interligand_patterns)

def process_functional_groups(smiles: str):
    images = highlight_functional_groups(smiles)
    if images is None or len(images) == 0:
        return [], "No functional groups recognized or invalid SMILES."
    return images, f"Found {len(images)} functional group(s)."

def process_interligand_moieties(smiles: str):
    images = highlight_interligand_moieties(smiles)
    if images is None or len(images) == 0:
        return [], "No interligand moieties recognized or invalid SMILES."
    return images, f"Found {len(images)} interligand moiety(ies)."

# -----------------------------
# Rotatable Bond Functions
# -----------------------------
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
        size=IMAGE_SIZE,
        highlightBonds=rot_bonds,
        legend="Rotatable Bonds"
    )
    return img

def process_rotatable(smiles: str):
    img = highlight_rotatable_bonds(smiles)
    if img is None:
        return None, "No rotatable bonds recognized or invalid SMILES."
    return img, "Rotatable bonds highlighted."

# -----------------------------
# Chiral Center Functions
# -----------------------------
def highlight_chiral_centers(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if not chiral_centers:
        return None
    highlight_atoms = [idx for idx, _ in chiral_centers]
    legend = "Chiral Centers: " + ", ".join(f"{idx} ({ch})" for idx, ch in chiral_centers)
    img = Draw.MolToImage(
        mol,
        size=IMAGE_SIZE,
        highlightAtoms=highlight_atoms,
        legend=legend
    )
    return img

def process_chiral_centers(smiles: str):
    img = highlight_chiral_centers(smiles)
    if img is None:
        return None, "No chiral centers recognized or invalid SMILES."
    return img, "Chiral centers highlighted."

# -----------------------------
# Combined Processing Function
# -----------------------------
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
    elif mode == "Chiral Centers":
        img, status_msg = process_chiral_centers(smiles)
        if img is None:
            return [], status_msg
        return [img], status_msg
    else:
        return [], "Invalid mode selected."

# -----------------------------
# Gradio Interface
# -----------------------------
with gr.Blocks() as demo:
    gr.Markdown("# Moliety: Molecular Feature Highlighter")
    gr.Markdown(
        "Enter a SMILES string below and select a highlighting mode. "
        "You can choose to highlight functional groups, interligand moieties, rotatable bonds, or chiral centers."
    )
    gr.Markdown("www.giorginolab.it")
    
    with gr.Row():
        smiles_input = gr.Textbox(
            label="Enter SMILES string", 
            placeholder="e.g. CC(=O)OC1=CC=CC=C1C(=O)O",
            lines=1
        )
        mode_dropdown = gr.Dropdown(
            label="Highlight Mode",
            choices=["Functional Groups", "Rotatable Bonds", "Interligand moieties", "Chiral Centers"],
            value="Functional Groups"
        )
    
    gallery = gr.Gallery(label="Highlighted Features", columns=3, height="auto")
    status = gr.Textbox(label="Status", interactive=False)

    smiles_input.submit(process_smiles_mode, inputs=[smiles_input, mode_dropdown], outputs=[gallery, status])
    mode_dropdown.change(process_smiles_mode, inputs=[smiles_input, mode_dropdown], outputs=[gallery, status])

demo.launch()
