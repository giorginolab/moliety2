---
title: Moliety2
colorFrom: indigo
colorTo: yellow
sdk: gradio
sdk_version: 6.10.0
app_file: app.py
pinned: false
license: gpl-3.0
short_description: Molecular moieties highlighter
---

# Moliety: Molecular moieties highlighter

www.giorginolab.it

Enter a SMILES string and select a molecular feature mode. Moliety can highlight functional groups, interligand moieties, rotatable bonds, chiral centers, potential stereogenic centers, DAYLIGHT SMARTS examples, Murcko scaffolds, hybridization states, Gasteiger charges, and protonation variants.

**WARNING: Mostly AI-generated and untested! Use at own risk.**

Based on SMARTS patterns provided with [OpenBabel](https://github.com/openbabel/openbabel/blob/master/data/SMARTS_InteLigand.txt) and [DAYLIGHT SMARTS examples](https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html). Protonation mode is supported by Durrant Lab's [dimorphite_dl](https://durrantlab.pitt.edu/dimorphite-dl/) library.

Hosted online at https://huggingface.co/spaces/tonigi/moliety2 .

## Local installation

Uses RDKit and Gradio:

```bash
uv sync
uv run gradio app.py
```

Run tests with:

```bash
uv run python -m unittest
```

## Hugging Face Space dependencies

HF Gradio Spaces install Python packages from `requirements.txt` during build. `uv.lock` is not used by the default Space launcher, so runtime dependencies needed by `app.py` must be listed in `requirements.txt`.
