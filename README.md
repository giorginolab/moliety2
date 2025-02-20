---
title: Moliety2
colorFrom: indigo
colorTo: yellow
sdk: gradio
sdk_version: 5.16.1
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

Based on SMARTS patterns provided with [OpenBabel](https://github.com/openbabel/openbabel/blob/master/data/SMARTS_InteLigand.txt) (SMARTS Patterns for Functional Group Classification by Christian Laggner) and DAYLIGHT [SMARTS examples](https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html).

Hosted online at https://huggingface.co/spaces/tonigi/moliety2 .




## Quick tests

#### Ibuprofen
        CC(C)Cc1ccc(cc1)[C@H](C)C(O)=O
#### Aspirin
        CC(=O)Oc1ccccc1C(=O)O
#### Suramin
        CC1=C(C=C(C=C1)C(=O)NC2=C3C(=CC(=CC3=C(C=C2)S(=O)(=O)O)S(=O)(=O)O)S(=O)(=O)O)NC(=O)C4=CC(=CC=C4)NC(=O)NC5=CC=CC(=C5)C(=O)NC6=C(C=CC(=C6)C(=O)NC7=C8C(=CC(=CC8=C(C=C7)S(=O)(=O)O)S(=O)(=O)O)S(=O)(=O)O)C


## Local installation

Uses RDKit and Gradio. Therefore:

```bash
python -m venv env
source env/bin/activate
pip install gradio rdkit
gradio app.py
```