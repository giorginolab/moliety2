---
title: Moliety2
colorFrom: indigo
colorTo: yellow
sdk: gradio
sdk_version: 6.19.0
app_file: app.py
pinned: false
license: gpl-3.0
short_description: Molecular moieties highlighter
---

# Moliety: Molecular moieties highlighter

www.giorginolab.it

Enter a SMILES string and select a molecular feature mode. Moliety can highlight functional groups, interligand moieties, SMARTS-RX moieties, rotatable bonds, chiral centers, potential stereogenic centers, DAYLIGHT SMARTS examples, Murcko scaffolds, hybridization states, Gasteiger charges, and protonation variants.

**WARNING: Mostly AI-generated and untested! Use at own risk.**

Based on SMARTS patterns provided with
 * [OpenBabel](https://github.com/openbabel/openbabel/blob/master/data/SMARTS_InteLigand.txt) 
 * [DAYLIGHT SMARTS examples](https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html). 
 * SMARTS-RX: Kogej, T., Kannas, C., Genheden, S. et al. [SMARTS-RX: a SMARTS-based representation of chemical functions for reactivity analysis](https://doi.org/10.1186/s13321-025-01136-8). J Cheminform 17, 177 (2025).
 
Protonation mode is supported by Durrant Lab's [dimorphite_dl](https://durrantlab.pitt.edu/dimorphite-dl/) library.

Hosted online at https://huggingface.co/spaces/tonigi/moliety2 .

## Local installation

Uses RDKit and Gradio:

```bash
uv sync
uv run gradio app.py
```


## Hugging Face Space dependencies

HF Gradio Spaces install Python packages from `requirements.txt` during build. `uv.lock` is not used by the default Space launcher, so runtime dependencies needed by `app.py` must be listed in `requirements.txt`.
