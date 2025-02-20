import gradio as gr
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from PIL import Image
from rdkit.Chem.Draw import MolDraw2DCairo
from io import BytesIO
import yaml  
from rdkit.Chem.Scaffolds import MurckoScaffold  # New import

IMAGE_SIZE=(800,800)

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
        with open("data/SMARTS_InteLigand.txt", "r") as f:
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
# Rotatable bond Definition 
# https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html
# -----------------------------

rotatable_patterns = {"DAYLIGHT defn.": Chem.MolFromSmarts("[!$(*#*)&!D1]-!@[!$(*#*)&!D1]")}


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
        if pattern is None:  # Skip uncompiled patterns
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
            img = Draw.MolToImage(
                mol,
                size=IMAGE_SIZE,
                highlightAtoms=list(highlight_atoms),
                highlightBonds=list(highlight_bonds),
                legend=name
            )
            images.append((img,name))
    return images


def process_functional_groups(smiles: str):
    images = highlight_by_patterns(smiles, compiled_patterns)
    if images is None or len(images) == 0:
        return [], "No functional groups recognized or invalid SMILES."
    return images, f"Found {len(images)} functional group(s)."

def process_interligand_moieties(smiles: str):
    images = highlight_by_patterns(smiles, compiled_interligand_patterns)
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
    drawer = MolDraw2DCairo(IMAGE_SIZE[0], IMAGE_SIZE[1])
    options = drawer.drawOptions()
    options.addAtomIndices = True  # Show atom number labels
    drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms, legend=legend)
    drawer.FinishDrawing()
    img_data = drawer.GetDrawingText()
    img = Image.open(BytesIO(img_data))
    return img

def process_chiral_centers(smiles: str):
    img = highlight_chiral_centers(smiles)
    if img is None:
        return None, "No chiral centers recognized or invalid SMILES."
    return img, "Chiral centers highlighted."

# -----------------------------
# Potential Stereo Functions
# -----------------------------
def process_potential_stereo(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES."
    # Use RDKit's FindPotentialStereo (assuming available)
    potential_stereocenters = Chem.FindPotentialStereo(mol)
    if not potential_stereocenters:
        return None, "No potential stereo centers found."
    images = []
    for sinfo in potential_stereocenters:
        highlight_atoms = [sinfo.centeredOn] 
        images.append(Draw.MolToImage(
            mol,
            size=IMAGE_SIZE,
            highlightAtoms=highlight_atoms,
            legend=sinfo.type.name
        ))
    return images, f"Found {len(highlight_atoms)} potential stereo center(s)."

# -----------------------------
# DAYLIGHT SMARTS Functions
# -----------------------------
def load_yaml_smarts():
    """
    Load and compile SMARTS from the YAML file.
    """
    try:
        with open("data/daylight_smarts.yml", "r") as f:
            data = yaml.safe_load(f)
    except Exception as e:
        print("Error loading daylight_smarts.yml:", e) 
        return {}
    compiled_yaml = {}
    for group in data.get("groups", []):
        group_name = group.get("name", "Unnamed Group")
        for subgroup in group.get("subgroups", []):
            subgroup_name = subgroup.get("name", "Unnamed Subgroup")
            if "subsubgroups" in subgroup:
                for subsub in subgroup.get("subsubgroups", []):
                    subsub_name = subsub.get("name", "Unnamed Subsubgroup")
                    for rule in subsub.get("rules", []):
                        if "smarts" in rule:
                            key = f"{group_name} > {subgroup_name} > {subsub_name} > {rule.get('name', 'Unnamed Rule')}"
                            compiled_yaml[key] = Chem.MolFromSmarts(rule.get("smarts"))
            elif "rules" in subgroup:
                for rule in subgroup.get("rules", []):
                    if "smarts" in rule:
                        key = f"{group_name}: {subgroup_name} - {rule.get('name', 'Unnamed Rule')}"
                        compiled_yaml[key] = Chem.MolFromSmarts(rule.get("smarts"))
    return compiled_yaml

def process_daylight_smarts_examples(smiles: str):
    """
    Highlight substructures using the SMARTS defined in the YAML file.
    """
    patterns = load_yaml_smarts()
    images = highlight_by_patterns(smiles, patterns)
    if images is None or len(images) == 0:
        return [], "No SMARTS examples recognized or invalid SMILES." 
    return images, f"Found {len(images)} SMARTS matches." 

# -----------------------------
# Scaffold Highlight Function
# -----------------------------
def process_scaffold(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [], "Invalid SMILES."
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None:
        return [], "No scaffold found."
    match = mol.GetSubstructMatch(scaffold)
    if not match:
        return [], "Scaffold not found as substructure."
    img = Draw.MolToImage(
        mol,
        size=IMAGE_SIZE,
        highlightAtoms=list(match),
        legend="Murcko Scaffold"
    )
    return [(img, "Scaffold")], "Scaffold highlighted."

# -----------------------------
# Combined Processing Function
# -----------------------------
def process_smiles_mode(smiles: str, mode: str):
    if mode == "Functional Groups":
        images, status_msg = process_functional_groups(smiles)
        return images, status_msg
    elif mode == "Rotatable Bonds":
        images, status_msg = process_rotatable(smiles)
        return images, status_msg
    elif mode == "Interligand Moieties":
        images, status_msg = process_interligand_moieties(smiles)
        return images, status_msg
    elif mode == "Chiral Centers":
        img, status_msg = process_chiral_centers(smiles)
        if img is None:
            return [], status_msg
        return [img], status_msg
    elif mode == "Potential Stereo":  
        images, status_msg = process_potential_stereo(smiles)
        return images, status_msg
    elif mode == "DAYLIGHT SMARTS Examples":  
        images, status_msg = process_daylight_smarts_examples(smiles)
        return images, status_msg
    elif mode == "Murcko Scaffold Highlight":  # New mode branch
        images, status_msg = process_scaffold(smiles)
        return images, status_msg
    else:
        return [], "Invalid mode selected."

# -----------------------------
# Gradio Interface
# -----------------------------
with gr.Blocks() as demo:
    gr.Markdown("# Moliety: Molecular Feature Highlighter")
    gr.Markdown("www.giorginolab.it")
    gr.Markdown(
        "Boost your impostor syndrome by uploading a molecule in SMILES form and count all the moieties you were supposed to know by heart. <br/>"
        "Enter a SMILES string below and select a highlighting mode. "
        "You can choose to highlight functional groups, interligand moieties, rotatable bonds, chiral centers, or potential stereo centers."
    )
    gr.Markdown("**WARNING: Mostly AI-generated and untested! Use at own risk.**")
    gr.Markdown("Based on SMARTS patterns provided with [OpenBabel](https://github.com/openbabel/openbabel/blob/master/data/SMARTS_InteLigand.txt) and [DAYLIGHT SMARTS examples](https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html).")
    
    with gr.Row():
        smiles_input = gr.Textbox(
            label="Enter SMILES string", 
            placeholder="e.g. CC(=O)OC1=CC=CC=C1C(=O)O",
            lines=3
        )
        mode_dropdown = gr.Dropdown(
            label="Highlight Mode",
            choices=["Functional Groups", "Rotatable Bonds", "Interligand Moieties", "Chiral Centers", 
                     "Potential Stereo", "DAYLIGHT SMARTS Examples", "Murcko Scaffold Highlight"],  # Added new mode
            value="Functional Groups"
        )
    
    # Update gr.Examples component to maintain existing examples and add new ones from the README
    gr.Examples(
        examples=[
            ["CC(=O)Oc1ccccc1C(=O)O", "Functional Groups"],
            ["CC(=O)OC1=CC=CC=C1C(=O)O", "Functional Groups"],
            ["CCOC(=O)C1=CC=CC=C1", "Rotatable Bonds"],
            ["CC(C(=O)O)N", "Chiral Centers"],
            ["CC(C)Cc1ccc(cc1)C(C)C(=O)O", "Functional Groups"],
            ["CC1=C(C=C(C=C1)C(=O)NC2=C3C(=CC(=CC3=C(C=C2)S(=O)(=O)O)S(=O)(=O)O)S(=O)(=O)O)NC(=O)C4=CC(=CC=C4)NC(=O)NC5=CC=CC(=C5)C(=O)NC6=C(C=CC(=C6)C(=O)NC7=C8C(=CC(=CC8=C(C=C7)S(=O)(=O)O)S(=O)(=O)O)S(=O)(=O)O)C", "Functional Groups"]
        ],
        example_labels=["Aspirin", "Aspirin (kekulized)", "Ethylbenzoate", "DL-Alanine", "Ibuprofen", "Suramin"],
        inputs=[smiles_input, mode_dropdown],
        label="Examples"
    )
    
    gallery = gr.Gallery(label="Highlighted Features", columns=3, height="auto")
    status = gr.Textbox(label="Status", interactive=False)

    smiles_input.submit(process_smiles_mode, inputs=[smiles_input, mode_dropdown], outputs=[gallery, status])
    mode_dropdown.change(process_smiles_mode, inputs=[smiles_input, mode_dropdown], outputs=[gallery, status])

demo.launch()
