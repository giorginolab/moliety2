"""Gradio UI assembly for Moliety2."""

from __future__ import annotations

import gradio as gr

from .features import FEATURE_MODES, depict_input_smiles, process_smiles_main

EXAMPLES = [
    ("CC(=O)Oc1ccccc1C(=O)O", "Aspirin"),
    ("CC(=O)OC1=CC=CC=C1C(=O)O", "Aspirin (kekulized)"),
    ("CCOC(=O)C1=CC=CC=C1", "Ethylbenzoate"),
    ("CC(C(=O)O)N", "DL-Alanine"),
    ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "Ibuprofen"),
    (
        "CC1=C(C=C(C=C1)C(=O)NC2=C3C(=CC(=CC3=C(C=C2)S(=O)(=O)O)S(=O)(=O)O)S(=O)(=O)O)NC(=O)C4=CC(=CC=C4)NC(=O)NC5=CC=CC(=C5)C(=O)NC6=C(C=CC(=C6)C(=O)NC7=C8C(=CC(=CC8=C(C=C7)S(=O)(=O)O)S(=O)(=O)O)S(=O)(=O)O)C",
        "Suramin",
    ),
    ("C=CC#N", "Acrylonitrile"),
    ("C1=CC=CC=C1", "Benzene"),
    ("CCO", "Ethanol"),
    ("CC(=O)O", "Acetic acid"),
    ("CN1C=NC2=C1C(=O)N(C)C(=O)N2C", "Caffeine"),
    ("CC(Cl)CC(F)CN", "4-chloro-2-fluoropentan-1-amine"),
    ("c1ccc2c(c1)cccc2", "Naphthalene"),
    ("CC1=C(C2=C(C=C1)C=CC=C2)CC(=O)O", "Naproxen"),
    ("NCc1ccccc1", "Benzylamine"),
]


def update_accordion_visibility(mode: str):
    selected_mode = next((feature_mode for feature_mode in FEATURE_MODES if feature_mode.name == mode), None)
    return gr.update(visible=bool(selected_mode and selected_mode.needs_ph))


def process_smiles_ui(smiles: str, mode: str, min_ph: float, max_ph: float):
    input_image = depict_input_smiles(smiles)
    output_images, status = process_smiles_main(smiles, mode, min_ph, max_ph)
    return input_image, output_images, status


def build_demo() -> gr.Blocks:
    mode_names = [mode.name for mode in FEATURE_MODES]

    with gr.Blocks() as demo:
        gr.Markdown("# Moliety: Molecular Feature Highlighter")
        gr.Markdown("www.giorginolab.it")
        gr.Markdown(
            "Enter a SMILES string below and select a molecular feature mode. "
            "The app highlights matching moieties, stereochemical features, scaffolds, charges, or protonation variants."
        )
        gr.Markdown("**WARNING: Mostly AI-generated and untested! Use at own risk.**")
        gr.Markdown(
            "Based on SMARTS patterns from [OpenBabel](https://github.com/openbabel/openbabel/blob/master/data/SMARTS_InteLigand.txt) "
            "and [DAYLIGHT SMARTS examples](https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html). "
            "Protonation mode uses Durrant Lab's [dimorphite_dl](https://durrantlab.pitt.edu/dimorphite-dl/) library."
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
                choices=mode_names,
                value=mode_names[0],
            )

        with gr.Accordion("pH Settings", visible=False) as ph_accordion:
            with gr.Row():
                min_ph = gr.Slider(
                    minimum=0,
                    maximum=14,
                    value=7.0,
                    step=0.5,
                    label="Minimum pH",
                )
                max_ph = gr.Slider(
                    minimum=0,
                    maximum=14,
                    value=7.0,
                    step=0.5,
                    label="Maximum pH",
                )

        mode_dropdown.change(
            update_accordion_visibility,
            inputs=[mode_dropdown],
            outputs=[ph_accordion],
        )

        gr.Examples(
            examples=[example for example, _label in EXAMPLES],
            example_labels=[label for _example, label in EXAMPLES],
            inputs=smiles_input,
            examples_per_page=30,
        )

        gr.Markdown("## Input")
        input_image = gr.Image(
            label="Input Molecule",
            type="filepath",
            height=420,
            interactive=False,
        )
        gr.Markdown("---\n## Outputs")
        gallery = gr.Gallery(label="Feature Outputs", columns=3, height="auto")
        status = gr.Textbox(label="Status", interactive=False)

        inputs = [smiles_input, mode_dropdown, min_ph, max_ph]
        for component in (smiles_input, mode_dropdown, min_ph, max_ph):
            event = component.submit if component is smiles_input else component.change
            event(process_smiles_ui, inputs=inputs, outputs=[input_image, gallery, status])

    return demo
