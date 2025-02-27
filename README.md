---
title: Moliety2
colorFrom: indigo
colorTo: yellow
sdk: gradio
sdk_version: 5.19.0
app_file: app.py
pinned: false
license: gpl-3.0
short_description: Recognize functional groups
---

# Moliety: Molecular Feature Highlighter

www.giorginolab.it

Boost your impostor syndrome by uploading a molecule in SMILES form and count all the moieties you were supposed to know by heart.

Enter a SMILES string and select a highlighting mode. You can choose to highlight functional groups, interligand moieties, rotatable bonds, or chiral centers.

**WARNING: Mostly AI-generated and untested! Use at own risk.**

Based on SMARTS patterns provided with [OpenBabel](https://github.com/openbabel/openbabel/blob/master/data/SMARTS_InteLigand.txt) and [DAYLIGHT SMARTS examples](https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html). Protonation mode supported by Durrant Lab's [dimorphite_dl](https://durrantlab.pitt.edu/dimorphite-dl/) library.

Hosted online at https://huggingface.co/spaces/tonigi/moliety2 .




## Local installation

Uses RDKit and Gradio. Therefore:

```bash
python -m venv env
source env/bin/activate
pip install gradio rdkit
gradio app.py
```