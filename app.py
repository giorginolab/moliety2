import gradio as gr
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem  # Add this import
from rotatable_bonds import process_rotatable
from dimorphite_dl import dimorphite_dl
import traceback

from file_helpers import load_interligand_moieties, load_smarts_patterns_from_csv
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
    return process_by_patterns(
        smiles,
        compiled_interligand_patterns,
        "No interligand moieties recognized or invalid SMILES.",
    )


def daylight_smarts_examples(smiles: str):
    patterns = load_smarts_patterns_from_csv()
    return process_by_patterns(
        smiles, patterns, "No SMARTS examples recognized or invalid SMILES."
    )


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
    img = mol_to_svg(
        mol,
        IMAGE_SIZE,
        highlightAtoms=highlight_atoms,
        legend=legend,
        atomLabels=atom_labels,
    )
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
    svg = mol_to_svg(
        mol,
        IMAGE_SIZE,
        highlightAtoms=highlight_atoms,
        legend="Potential Stereogenic Centers",
        atomLabels=atom_labels,
    )

    return [
        (svg, "Potential Stereogenic Centers")
    ], f"Found {len(potential_stereocenters)} potential stereogenic center(s)."


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
    img = mol_to_svg(
        mol,
        IMAGE_SIZE,
        highlightAtoms=list(match),
        highlightBonds=highlight_bonds,
        legend="Murcko Scaffold",
    )
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
    img = mol_to_svg(
        mol,
        IMAGE_SIZE,
        highlightAtoms=highlight_atoms,
        legend="Hybridization States",
        atomLabels=atom_labels,
    )

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
        charge = atom.GetDoubleProp("_GasteigerCharge")
        atom_labels[atom.GetIdx()] = f"{charge:.3f}"
        highlight_atoms.append(atom.GetIdx())

    if not highlight_atoms:
        return [], "Could not compute Gasteiger charges."

    # Generate image with charge labels, showing hydrogens
    img = mol_to_svg(
        mol,
        IMAGE_SIZE,
        highlightAtoms=highlight_atoms,
        legend="Gasteiger Charges (including H)",
        atomLabels=atom_labels,
    )

    return [
        (img, "Gasteiger Charges")
    ], "Gasteiger charges computed and displayed (including hydrogens)."


# -----------------------------
# Protonation Function
# -----------------------------
def protonate_ph(smiles: str, min_ph: float, max_ph: float):
    """Protonate molecule at given pH range using dimorphite-dl."""
    try:
        protonated_mols = dimorphite_dl.run_with_mol_list(
            [Chem.MolFromSmiles(smiles)],
            min_ph=min_ph,
            max_ph=max_ph,
            pka_precision=1.0,
        )

        if not protonated_mols:
            return [], "No protonation variants found."

        images = []
        for i, mol in enumerate(protonated_mols):
            # Generate SVG for each protonated variant
            svg = mol_to_svg(
                mol,
                IMAGE_SIZE,
                legend=f"Protonated variant {i + 1} at pH {min_ph}-{max_ph}",
            )
            images.append((svg, f"Variant {i + 1}"))

        return (
            images,
            f"Found {len(images)} protonation variant(s) at pH {min_ph}-{max_ph}.",
        )
    except Exception as e:
        print(traceback.format_exc())
        return [], f"Error during protonation: {str(e)}"


# -----------------------------
# Combined Processing Function
# -----------------------------
# Modified process_smiles_mode: add SMILES validity check for Functional Groups
def process_smiles_main(
    smiles: str, mode: str, min_ph: float = 7.0, max_ph: float = 7.0
):
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
    elif mode == "Protonation":  # Modified mode name
        images, status_msg = protonate_ph(smiles, min_ph, max_ph)
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
        "Based on SMARTS patterns provided with [OpenBabel](https://github.com/openbabel/openbabel/blob/master/data/SMARTS_InteLigand.txt) and [DAYLIGHT SMARTS examples](https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html). Protonation mode supported by Durrant Lab's [dimorphite_dl](https://durrantlab.pitt.edu/dimorphite-dl/) library."
    )

    with gr.Row():
        smiles_input = gr.Textbox(
            label="Enter SMILES string",
            placeholder="e.g. CC(=O)OC1=CC=CC=C1C(=O)O",
            lines=2,
            submit_btn=True,
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
                "Protonation",  # Modified mode name
            ],
            value="Functional Groups",
        )

    # Add pH controls in accordion
    with gr.Accordion("pH Settings", visible=False) as ph_accordion:
        with gr.Row():
            min_ph = gr.Slider(
                minimum=0, maximum=14, value=7.0, step=0.5, label="Minimum pH"
            )
            max_ph = gr.Slider(
                minimum=0, maximum=14, value=7.0, step=0.5, label="Maximum pH"
            )

    # Update visibility of pH controls based on mode
    def update_accordion_visibility(mode):
        return gr.update(visible=(mode == "Protonation"))

    mode_dropdown.change(
        update_accordion_visibility, inputs=[mode_dropdown], outputs=[ph_accordion]
    )

    # Update gr.Examples component with new examples
    gr.Examples(
        examples=[
            "CC(=O)Oc1ccccc1C(=O)O",
            "CC(=O)OC1=CC=CC=C1C(=O)O",
            "CCOC(=O)C1=CC=CC=C1",
            "CC(C(=O)O)N",
            "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
            "CC1=C(C=C(C=C1)C(=O)NC2=C3C(=CC(=CC3=C(C=C2)S(=O)(=O)O)S(=O)(=O)O)S(=O)(=O)O)NC(=O)C4=CC(=CC=C4)NC(=O)NC5=CC=CC(=C5)C(=O)NC6=C(C=CC(=C6)C(=O)NC7=C8C(=CC(=CC8=C(C=C7)S(=O)(=O)O)S(=O)(=O)O)S(=O)(=O)O)C",
            "C=CC#N",
            "C1=CC=CC=C1",
            "CCO",
            "CC(=O)O",
            "CN1C=NC2=C1C(=O)N(C)C(=O)N2C",
            "CC(Cl)CC(F)CN",
            "c1ccc2c(c1)cccc2",
            "CC1=C(C2=C(C=C1)C=CC=C2)CC(=O)O",
            "NCc1ccccc1",
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
            "Caffeine",
            "4-chloro-2-fluoropentan-1-amine",
            "Naphthalene",
            "Naproxen",
            "Benzylamine",
        ],
        inputs=smiles_input,
        examples_per_page=30,
    )

    gallery = gr.Gallery(label="Highlighted Features", columns=3, height="auto")
    status = gr.Textbox(label="Status", interactive=False)

    # Update the submit handlers to include pH values
    smiles_input.submit(
        process_smiles_main,
        inputs=[smiles_input, mode_dropdown, min_ph, max_ph],
        outputs=[gallery, status],
    )
    mode_dropdown.change(
        process_smiles_main,
        inputs=[smiles_input, mode_dropdown, min_ph, max_ph],
        outputs=[gallery, status],
    )
    min_ph.change(
        process_smiles_main,
        inputs=[smiles_input, mode_dropdown, min_ph, max_ph],
        outputs=[gallery, status],
    )
    max_ph.change(
        process_smiles_main,
        inputs=[smiles_input, mode_dropdown, min_ph, max_ph],
        outputs=[gallery, status],
    )

demo.launch()
