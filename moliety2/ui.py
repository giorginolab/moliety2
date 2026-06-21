"""Gradio UI assembly for Moliety2."""

from __future__ import annotations

import gradio as gr

from .features import FEATURE_MODES, depict_input_smiles, feature_mode_source, process_smiles_main

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

APP_CSS = """
:root {
    --moliety-ink: #15202b;
    --moliety-muted: #5b6673;
    --moliety-line: #d8e0e7;
    --moliety-panel: rgba(255, 255, 255, 0.92);
    --moliety-teal: #12806f;
    --moliety-gold: #c98511;
    --moliety-blue: #2f5f98;
}

.gradio-container {
    background:
        linear-gradient(135deg, rgba(18, 128, 111, 0.10), transparent 34%),
        linear-gradient(315deg, rgba(201, 133, 17, 0.14), transparent 30%),
        #f6f8fb;
    color: var(--moliety-ink);
}

.app-shell {
    max-width: 1180px;
    margin: 0 auto;
}

.app-header {
    padding: 30px 0 18px;
    border-bottom: 1px solid var(--moliety-line);
    margin-bottom: 22px;
}

.brand-row {
    display: flex;
    align-items: flex-end;
    justify-content: space-between;
    gap: 18px;
    flex-wrap: wrap;
}

.brand-title {
    color: var(--moliety-ink);
    font-size: clamp(2.1rem, 4vw, 4.15rem);
    line-height: 0.98;
    font-weight: 850;
    margin: 0;
}

.brand-subtitle {
    color: var(--moliety-muted);
    font-size: 1.03rem;
    line-height: 1.55;
    margin: 14px 0 0;
    max-width: 700px;
}

.brand-meta {
    display: flex;
    gap: 8px;
    flex-wrap: wrap;
    justify-content: flex-end;
}

.site-link {
    border: 1px solid var(--moliety-line);
    background: rgba(255, 255, 255, 0.72);
    border-radius: 999px;
    color: #2e475c;
    font-size: 0.82rem;
    font-weight: 700;
    padding: 7px 11px;
    text-decoration: none;
    white-space: nowrap;
}

.site-link:hover {
    border-color: rgba(18, 128, 111, 0.42);
    color: var(--moliety-teal);
}

.notice {
    border-left: 4px solid var(--moliety-gold);
    background: rgba(255, 255, 255, 0.72);
    border-radius: 8px;
    color: #5d430f;
    padding: 10px 13px;
    margin-top: 18px;
    font-size: 0.92rem;
}

.workspace-panel,
.results-panel,
.examples-panel {
    background: var(--moliety-panel);
    border: 1px solid var(--moliety-line);
    border-radius: 8px;
    box-shadow: 0 18px 46px rgba(42, 57, 77, 0.08);
}

.workspace-panel {
    padding: 18px;
}

.control-column {
    min-width: 320px;
}

.preview-column {
    min-width: 320px;
}

.section-label {
    color: var(--moliety-blue);
    font-size: 0.78rem;
    font-weight: 800;
    letter-spacing: 0.07em;
    text-transform: uppercase;
    margin: 0 0 10px;
}

.source-note {
    color: var(--moliety-muted);
    border: 1px solid rgba(47, 95, 152, 0.16);
    background: rgba(47, 95, 152, 0.06);
    border-radius: 8px;
    padding: 10px 12px;
}

.source-note p {
    margin: 0;
}

.run-button button {
    min-height: 46px;
    border-radius: 8px;
    font-weight: 800;
    box-shadow: 0 10px 24px rgba(18, 128, 111, 0.22);
}

.input-preview .image-container,
.input-preview .wrap,
.input-preview img {
    background: #ffffff;
    border-radius: 8px;
}

.examples-panel,
.results-panel {
    margin-top: 18px;
    padding: 16px;
}

.results-panel .gallery {
    border-radius: 8px;
}

.status-box textarea {
    font-weight: 700;
    color: var(--moliety-ink);
}

footer {
    opacity: 0.72;
}

@media (max-width: 760px) {
    .app-header {
        padding-top: 18px;
    }

    .brand-title {
        font-size: 2.35rem;
    }

    .brand-meta {
        justify-content: flex-start;
    }

    .workspace-panel,
    .examples-panel,
    .results-panel {
        padding: 12px;
    }
}
"""


def build_theme() -> gr.Theme:
    return gr.themes.Soft(
        primary_hue="teal",
        secondary_hue="amber",
        neutral_hue="slate",
        radius_size=gr.themes.sizes.radius_sm,
        spacing_size=gr.themes.sizes.spacing_md,
        font=[
            gr.themes.GoogleFont("Inter"),
            "ui-sans-serif",
            "system-ui",
            "sans-serif",
        ],
        font_mono=[
            gr.themes.GoogleFont("JetBrains Mono"),
            "ui-monospace",
            "SFMono-Regular",
            "monospace",
        ],
    )


def update_mode_context(mode: str):
    selected_mode = next((feature_mode for feature_mode in FEATURE_MODES if feature_mode.name == mode), None)
    return gr.update(visible=bool(selected_mode and selected_mode.needs_ph)), feature_mode_source(mode)


def process_smiles_ui(smiles: str, mode: str, min_ph: float, max_ph: float):
    input_image = depict_input_smiles(smiles)
    output_images, status = process_smiles_main(smiles, mode, min_ph, max_ph)
    return input_image, output_images, status


def build_demo() -> gr.Blocks:
    mode_names = [mode.name for mode in FEATURE_MODES]

    with gr.Blocks(title="Moliety", fill_width=True) as demo:
        with gr.Column(elem_classes=["app-shell"]):
            gr.HTML(
                """
                <header class="app-header">
                    <div class="brand-row">
                        <div>
                            <h1 class="brand-title">Moliety</h1>
                            <p class="brand-subtitle">
                                Molecular feature highlighter for SMARTS patterns, stereochemistry,
                                scaffolds, charges, and protonation variants.
                            </p>
                        </div>
                        <div class="brand-meta" aria-label="Application capabilities">
                            <a class="site-link" href="https://www.giorginolab.it" target="_blank" rel="noopener noreferrer">
                                www.giorginolab.it
                            </a>
                        </div>
                    </div>
                    <div class="notice">
                        Mostly AI-generated and untested. Use at own risk.
                    </div>
                </header>
                """
            )

            with gr.Row(elem_classes=["workspace-panel"], equal_height=True):
                with gr.Column(scale=7, elem_classes=["control-column"]):
                    gr.Markdown("Input", elem_classes=["section-label"])
                    smiles_input = gr.Textbox(
                        label="SMILES string",
                        value=EXAMPLES[0][0],
                        placeholder="e.g. CC(=O)OC1=CC=CC=C1C(=O)O",
                        lines=3,
                        submit_btn=True,
                    )
                    with gr.Row():
                        mode_dropdown = gr.Dropdown(
                            label="Highlight mode",
                            choices=mode_names,
                            value=mode_names[0],
                            scale=3,
                        )
                        run_button = gr.Button(
                            "Analyze molecule",
                            variant="primary",
                            scale=1,
                            elem_classes=["run-button"],
                        )
                    source_markdown = gr.Markdown(
                        feature_mode_source(mode_names[0]),
                        elem_classes=["source-note"],
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

                with gr.Column(scale=5, elem_classes=["preview-column"]):
                    gr.Markdown("Molecule", elem_classes=["section-label"])
                    input_image = gr.Image(
                        label="Input molecule",
                        type="filepath",
                        height=430,
                        interactive=False,
                        buttons=["fullscreen"],
                        elem_classes=["input-preview"],
                    )

            mode_dropdown.change(
                update_mode_context,
                inputs=[mode_dropdown],
                outputs=[ph_accordion, source_markdown],
            )

            with gr.Column(elem_classes=["examples-panel"]):
                gr.Markdown("Examples", elem_classes=["section-label"])
                gr.Examples(
                    examples=[example for example, _label in EXAMPLES],
                    example_labels=[label for _example, label in EXAMPLES],
                    inputs=smiles_input,
                    examples_per_page=30,
                )

            with gr.Column(elem_classes=["results-panel"]):
                gr.Markdown("Results", elem_classes=["section-label"])
                gallery = gr.Gallery(
                    label="Feature outputs",
                    columns=3,
                    height="auto",
                    preview=True,
                )
                status = gr.Textbox(
                    label="Status",
                    interactive=False,
                    lines=1,
                    elem_classes=["status-box"],
                )

            inputs = [smiles_input, mode_dropdown, min_ph, max_ph]
            run_button.click(process_smiles_ui, inputs=inputs, outputs=[input_image, gallery, status])
            demo.load(process_smiles_ui, inputs=inputs, outputs=[input_image, gallery, status])
            for component in (smiles_input, mode_dropdown, min_ph, max_ph):
                event = component.submit if component is smiles_input else component.change
                event(process_smiles_ui, inputs=inputs, outputs=[input_image, gallery, status])

    demo.theme = build_theme()
    demo.css = APP_CSS
    return demo
