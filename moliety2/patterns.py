"""SMARTS definitions and cached pattern loading."""

from __future__ import annotations

import csv
import json
from functools import lru_cache
from pathlib import Path

from rdkit import Chem, rdBase

DATA_DIR = Path(__file__).resolve().parent.parent / "data"

FUNCTIONAL_GROUP_SMARTS = {
    "Hydroxyl": "[OX2H]",
    "Carbonyl": "[CX3]=[OX1]",
    "Amine": "[NX3;H2,H1;!$(NC=O)]",
    "Carboxylic Acid": "C(=O)[OX2H1]",
    "Ester": "C(=O)O[C]",
    "Ether": "[OD2]([#6])[#6]",
    "Halide": "[F,Cl,Br,I]",
}

ROTATABLE_BOND_SMARTS = {
    "DAYLIGHT defn.": "[!$(*#*)&!D1]-!@[!$(*#*)&!D1]",
    "RDKit defn.": "[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]",
}


def _compile_smarts(smarts: str, *, quiet: bool = False) -> Chem.Mol | None:
    if not quiet:
        return Chem.MolFromSmarts(smarts)

    with rdBase.BlockLogs():
        return Chem.MolFromSmarts(smarts)


def _compile_dict(patterns: dict[str, str], *, quiet: bool = False) -> dict[str, Chem.Mol | None]:
    return {name: _compile_smarts(smarts, quiet=quiet) for name, smarts in patterns.items()}


@lru_cache(maxsize=1)
def functional_group_patterns() -> dict[str, Chem.Mol | None]:
    return _compile_dict(FUNCTIONAL_GROUP_SMARTS)


@lru_cache(maxsize=1)
def rotatable_patterns() -> dict[str, Chem.Mol | None]:
    return _compile_dict(ROTATABLE_BOND_SMARTS)


@lru_cache(maxsize=1)
def interligand_smarts() -> dict[str, str]:
    moieties: dict[str, str] = {}
    with (DATA_DIR / "SMARTS_InteLigand.txt").open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#") or ":" not in line:
                continue

            name, pattern = (part.strip() for part in line.split(":", 1))
            if pattern and pattern.startswith("["):
                moieties[name] = pattern

    return moieties


@lru_cache(maxsize=1)
def interligand_patterns() -> dict[str, Chem.Mol | None]:
    moieties = interligand_smarts()
    return _compile_dict(moieties, quiet=True)


@lru_cache(maxsize=1)
def daylight_smarts_patterns() -> dict[str, Chem.Mol]:
    patterns: dict[str, Chem.Mol] = {}
    with (DATA_DIR / "smarts_examples_parsed.csv").open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            pattern = row["Pattern"].strip()
            if not pattern:
                continue

            mol = _compile_smarts(pattern, quiet=True)
            if mol is None:
                continue

            name = " > ".join(
                [
                    row["Main Topic"],
                    row["Subtopic"],
                    row["Sub-sub-topic"],
                    row["Rule Name"],
                ]
            )
            patterns[name] = mol

    return patterns


@lru_cache(maxsize=1)
def smartsrx_patterns() -> dict[str, Chem.Mol]:
    patterns: dict[str, Chem.Mol] = {}
    with (DATA_DIR / "smartsrx" / "smartsrx.json").open(encoding="utf-8") as handle:
        rows = json.load(handle)["data"]

    for row in rows:
        pattern = row["smarts"].strip()
        if not pattern:
            continue

        mol = _compile_smarts(pattern, quiet=True)
        if mol is None:
            continue

        name = " > ".join(
            [
                row["category"],
                row["subcategory"],
                row["specific_type"],
            ]
        )
        patterns[name] = mol

    return patterns
