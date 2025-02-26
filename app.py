import gradio as gr
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold 
from rdkit.Chem import AllChem  # Add this import
from rotatable_bonds import process_rotatable
from dimorphite_dl import dimorphite_dl
import traceback

from file_helpers import load_interligand_moieties, load_yaml_smarts
from utils import mol_to_svg, highlight_by_patterns, IMAGE_SIZE

# Remove mol_to_svg function as it's now in utils.py
# Remove highlight_by_patterns function as it's now in utils.py

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
    "Halide": "[F,Cl,Br,I]",
}
compiled_patterns = {
    name: Chem.MolFromSmarts(smart) for name, smart in functional_groups.items()
}


# -----------------------------
# Interligand Moieties Definitions
# -----------------------------
interligand_moieties = load_interligand_moieties()
compiled_interligand_patterns = {
    name: Chem.MolFromSmarts(smart) for name, smart in interligand_moieties.items()
}



# -----------------------------
# Generic Highlighting Function
# -----------------------------


def process_by_patterns(smiles: str, patterns: dict, not_found_msg: str):
    images = highlight_by_patterns(smiles, patterns)
    if not images:
        return [], not_found_msg
    return images, f"Found {len(images)} match(es)."

# Modified process_functional_groups: removed SMILES validity check
def functional_groups(smiles: str):
    images = highlight_by_patterns(smiles, compiled_patterns)
    if not images:
        return [], "No functional groups recognized."
    return images, f"Found {len(images)} match(es)."

def interligand_moieties(smiles: str):
    return process_by_patterns(smiles, compiled_interligand_patterns, "No interligand moieties recognized or invalid SMILES.")

def daylight_smarts_examples(smiles: str):
    patterns = load_yaml_smarts()
    return process_by_patterns(smiles, patterns, "No SMARTS examples recognized or invalid SMILES.")


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
    
    # Create labels dictionary for chiral centers
    atom_labels = {}
    for idx, chirality in chiral_centers:
        atom_labels[idx] = chirality  # Will show R or S (or ?)
    
    legend = "Chiral Centers: " + ", ".join(
        f"{idx} ({ch})" for idx, ch in chiral_centers
    )
    img = mol_to_svg(mol, IMAGE_SIZE, 
                     highlightAtoms=highlight_atoms, 
                     legend=legend,
                     atomLabels=atom_labels)
    return img


def chiral_centers(smiles: str):
    img = highlight_chiral_centers(smiles)
    if img is None:
        return None, "No chiral centers recognized or invalid SMILES."
    return [(img, "Chiral centers")], "Chiral centers highlighted."


# -----------------------------
# Potential Stereo Functions
# -----------------------------
def stereocenters(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES."
    # Use RDKit's FindPotentialStereo
    potential_stereocenters = Chem.FindPotentialStereo(mol)
    if not potential_stereocenters:
        return None, "No potential stereo centers found."
    
    highlight_atoms = []
    atom_labels = {}
    for sinfo in potential_stereocenters:
        highlight_atoms.append(sinfo.centeredOn)
        atom_labels[sinfo.centeredOn] = sinfo.type.name
    
    # Create single image with all centers highlighted
    svg = mol_to_svg(mol, IMAGE_SIZE, 
                     highlightAtoms=highlight_atoms, 
                     legend="Potential Stereogenic Centers",
                     atomLabels=atom_labels)
    
    return [(svg, "Potential Stereogenic Centers")], f"Found {len(potential_stereocenters)} potential stereogenic center(s)."


# -----------------------------
# Updated Scaffold Highlight Function
# -----------------------------
def scaffold(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [], "Invalid SMILES."
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None:
        return [], "No scaffold found."
    match = mol.GetSubstructMatch(scaffold)
    if not match:
        return [], "Scaffold not found as substructure."
    # Compute bonds where both atoms are in the match
    highlight_bonds = []
    for bond in mol.GetBonds():
        if bond.GetBeginAtomIdx() in match and bond.GetEndAtomIdx() in match:
            highlight_bonds.append(bond.GetIdx())
    # Modified to output SVG
    img = mol_to_svg(mol, IMAGE_SIZE, highlightAtoms=list(match), highlightBonds=highlight_bonds, legend="Murcko Scaffold")
    return [(img, "Murcko Scaffold")], "Scaffold highlighted."


# -----------------------------
# Hybridization State Functions
# -----------------------------
def hybridization(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [], "Invalid SMILES."
    
    # Create atom labels dictionary with hybridization states
    atom_labels = {}
    highlight_atoms = []
    for atom in mol.GetAtoms():
        hyb = atom.GetHybridization()
        # Skip atoms with undefined hybridization
        if hyb != Chem.HybridizationType.UNSPECIFIED:
            atom_labels[atom.GetIdx()] = hyb.name
            highlight_atoms.append(atom.GetIdx())
    
    if not highlight_atoms:
        return [], "No hybridization states to display."
    
    # Generate image with hybridization labels
    img = mol_to_svg(mol, IMAGE_SIZE, 
                     highlightAtoms=highlight_atoms,
                     legend="Hybridization States",
                     atomLabels=atom_labels)
    
    return [(img, "Hybridization States")], "Hybridization states highlighted."


# -----------------------------
# Gasteiger Charges Function
# -----------------------------
def gasteiger_charges(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return [], "Invalid SMILES."
    
    # Add explicit hydrogens to the molecule
    mol = Chem.AddHs(mol)
    
    # Compute Gasteiger charges using AllChem
    AllChem.ComputeGasteigerCharges(mol)
    
    # Create atom labels dictionary with charges
    atom_labels = {}
    highlight_atoms = []
    for atom in mol.GetAtoms():
        charge = atom.GetDoubleProp('_GasteigerCharge')
        atom_labels[atom.GetIdx()] = f"{charge:.3f}"
        highlight_atoms.append(atom.GetIdx())
    
    if not highlight_atoms:
        return [], "Could not compute Gasteiger charges."
    
    # Generate image with charge labels, showing hydrogens
    img = mol_to_svg(mol, IMAGE_SIZE, 
                     highlightAtoms=highlight_atoms,
                     legend="Gasteiger Charges (including H)",
                     atomLabels=atom_labels)
    
    return [(img, "Gasteiger Charges")], "Gasteiger charges computed and displayed (including hydrogens)."


# -----------------------------
# Protonation Function
# -----------------------------
def protonate_ph7(smiles: str):
    """Protonate molecule at pH 7 using dimorphite-dl."""
    try:
        protonated_mols = dimorphite_dl.run_with_mol_list(
            [Chem.MolFromSmiles(smiles)],
            min_ph=7.0,
            max_ph=7.0,
            pka_precision=1.0
        )
        
        if not protonated_mols:
            return [], "No protonation variants found."
        
        images = []
        for i, mol in enumerate(protonated_mols):
            # Generate SVG for each protonated variant
            svg = mol_to_svg(mol, IMAGE_SIZE, 
                           legend=f"Protonated variant {i+1} at pH 7")
            images.append((svg, f"Variant {i+1}"))
        
        return images, f"Found {len(images)} protonation variant(s) at pH 7."
    except Exception as e:
        print(traceback.format_exc())
        return [], f"Error during protonation: {str(e)}"


# -----------------------------
# Combined Processing Function
# -----------------------------
# Modified process_smiles_mode: add SMILES validity check for Functional Groups
def process_smiles_main(smiles: str, mode: str):
    if Chem.MolFromSmiles(smiles) is None:
        return [], "Invalid SMILES."

    if mode == "Functional Groups":
        images, status_msg = functional_groups(smiles)
    elif mode == "Rotatable Bonds":
        images, status_msg = process_rotatable(smiles)  # Simplified call
    elif mode == "Interligand Moieties":
        images, status_msg = interligand_moieties(smiles)
    elif mode == "Chiral Centers":
        images, status_msg = chiral_centers(smiles)
    elif mode == "Potential Stereogenic Centers":
        images, status_msg = stereocenters(smiles)
    elif mode == "DAYLIGHT SMARTS Examples":
        images, status_msg = daylight_smarts_examples(smiles)
    elif mode == "Murcko Scaffold":
        images, status_msg = scaffold(smiles)
    elif mode == "Hybridization":
        images, status_msg = hybridization(smiles)
    elif mode == "Gasteiger Charges":
        images, status_msg = gasteiger_charges(smiles)
    elif mode == "Protonation (pH 7)":  # Add new mode
        images, status_msg = protonate_ph7(smiles)
    else:
        return [], "Invalid mode selected."

    return images, status_msg


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
    gr.Markdown(
        "Based on SMARTS patterns provided with [OpenBabel](https://github.com/openbabel/openbabel/blob/master/data/SMARTS_InteLigand.txt) and [DAYLIGHT SMARTS examples](https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html)."
    )

    with gr.Row():
        smiles_input = gr.Textbox(
            label="Enter SMILES string",
            placeholder="e.g. CC(=O)OC1=CC=CC=C1C(=O)O",
            lines=3,
        )
        mode_dropdown = gr.Dropdown(
            label="Highlight Mode",
            choices=[
                "Functional Groups",
                "Rotatable Bonds",
                "Interligand Moieties",
                "Chiral Centers",
                "Potential Stereogenic Centers",
                "DAYLIGHT SMARTS Examples",
                "Murcko Scaffold",
                "Hybridization",  # Add new mode
                "Gasteiger Charges",  # Add new mode
                "Protonation (pH 7)"  # Add new mode
            ],
            value="Functional Groups",
        )

    # Update gr.Examples component with new examples
    gr.Examples(
        examples=[
            ["CC(=O)Oc1ccccc1C(=O)O", "Functional Groups"],
            ["CC(=O)OC1=CC=CC=C1C(=O)O", "Functional Groups"],
            ["CCOC(=O)C1=CC=CC=C1", "Rotatable Bonds"],
            ["CC(C(=O)O)N", "Chiral Centers"],
            ["CC(C)Cc1ccc(cc1)C(C)C(=O)O", "Functional Groups"],
            [
                "CC1=C(C=C(C=C1)C(=O)NC2=C3C(=CC(=CC3=C(C=C2)S(=O)(=O)O)S(=O)(=O)O)S(=O)(=O)O)NC(=O)C4=CC(=CC=C4)NC(=O)NC5=CC=CC(=C5)C(=O)NC6=C(C=CC(=C6)C(=O)NC7=C8C(=CC(=CC8=C(C=C7)S(=O)(=O)O)S(=O)(=O)O)S(=O)(=O)O)C",
                "Functional Groups",
            ],
            # New examples for Hybridization
            ["C=CC#N", "Hybridization"],  # Shows SP, SP2, and SP3 carbons
            ["C1=CC=CC=C1", "Hybridization"],  # Benzene ring showing SP2
            # Examples for Gasteiger Charges
            ["CCO", "Gasteiger Charges"],  # Simple alcohol showing charge distribution
            ["CC(=O)O", "Gasteiger Charges"],  # Acetic acid showing polar groups
            # Add examples for missing modes
            ["O=C(O)C1N2C(=O)C3C(N=CN3C)C2=O", "Interligand Moieties"],  # Caffeine-like structure
            ["CC(Cl)CC(F)CN", "Potential Stereogenic Centers"],  # Multiple potential stereocenters
            ["c1ccc2c(c1)cccc2", "DAYLIGHT SMARTS Examples"],  # Naphthalene for aromatic patterns
            ["CC1=C(C2=C(C=C1)C=CC=C2)CC(=O)O", "Murcko Scaffold"],  # Naproxen scaffold
            ["CC(=O)O", "Protonation (pH 7)"],  # Acetic acid
            ["NCc1ccccc1", "Protonation (pH 7)"],  # Benzylamine
        ],
        example_labels=[
            "Aspirin",
            "Aspirin (kekulized)",
            "Ethylbenzoate",
            "DL-Alanine",
            "Ibuprofen",
            "Suramin",
            "Acrylonitrile",
            "Benzene",
            "Ethanol",
            "Acetic acid",
            "Caffeine-like",
            "Multiple stereocenters",
            "Naphthalene",
            "Naproxen",
            "Acetic acid protonation",
            "Benzylamine protonation"
        ],
        inputs=[smiles_input, mode_dropdown],
        label="Examples",
    )

    gallery = gr.Gallery(label="Highlighted Features", columns=3, height="auto")
    status = gr.Textbox(label="Status", interactive=False)

    smiles_input.submit(
        process_smiles_main,
        inputs=[smiles_input, mode_dropdown],
        outputs=[gallery, status],
    )
    mode_dropdown.change(
        process_smiles_main,
        inputs=[smiles_input, mode_dropdown],
        outputs=[gallery, status],
    )

demo.launch()
