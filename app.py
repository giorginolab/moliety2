import gradio as gr
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

# Define a dictionary of functional groups with corresponding SMARTS patterns.
functional_groups = {
    "Hydroxyl": "[OX2H]",
    "Carbonyl": "[CX3]=[OX1]",
    "Amine": "[NX3;H2,H1;!$(NC=O)]",
    "Carboxylic Acid": "C(=O)[OX2H1]",
    "Ester": "C(=O)O[C]",
    "Ether": "[OD2]([#6])[#6]",
    "Halide": "[F,Cl,Br,I]"
}

# Pre-compile the SMARTS queries.
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
            # and include the group name as a legend.
            img = Draw.MolToImage(
                mol,
                size=(300, 300),
                highlightAtoms=list(highlight_atoms),
                highlightBonds=list(highlight_bonds),
                legend=group_name
            )
            images.append(img)
    return images

def process_smiles(smiles):
    """
    Process the input SMILES string and return a tuple:
      - list of PIL Images for each recognized functional group
      - a status message.
    """
    images = highlight_functional_groups(smiles)
    if images is None or len(images) == 0:
        return [], "No functional groups recognized or invalid SMILES."
    return images, f"Found {len(images)} functional group(s)."

# Build the Gradio interface.
with gr.Blocks() as demo:
    gr.Markdown("# Functional Group and Moiety Recognizer")
    gr.Markdown(
        "Enter a SMILES string below. The app uses RDKit to identify common "
        "functional groups and then displays images of the molecule with each group highlighted."
    )
    smiles_input = gr.Textbox(
        label="Enter SMILES string", 
        placeholder="e.g. CC(=O)OC1=CC=CC=C1C(=O)O"
    )
    gallery = gr.Gallery(label="Recognized Functional Groups", columns=3, height="auto")
    status = gr.Textbox(label="Status", interactive=False)

    # When the user submits the SMILES, call process_smiles.
    smiles_input.submit(process_smiles, inputs=smiles_input, outputs=[gallery, status])

demo.launch()

