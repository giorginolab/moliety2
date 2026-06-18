"""Molecular feature handlers used by the Gradio app."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import dimorphite_dl
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold

from .patterns import (
    daylight_smarts_patterns,
    functional_group_patterns,
    interligand_patterns,
    rotatable_patterns,
    smartsrx_patterns,
)
from .rendering import IMAGE_SIZE, highlight_by_patterns, matching_bond_indices, mol_to_svg

GalleryResult = tuple[list[tuple[str, str]], str]
FeatureHandler = Callable[[Chem.Mol, str, float, float], GalleryResult]


@dataclass(frozen=True)
class FeatureMode:
    name: str
    handler: FeatureHandler
    needs_ph: bool = False


def _mol_from_smiles(smiles: str) -> Chem.Mol | None:
    with rdBase.BlockLogs():
        return Chem.MolFromSmiles(smiles or "")


def depict_input_smiles(smiles: str) -> str | None:
    mol = _mol_from_smiles(smiles)
    if mol is None:
        return None
    return mol_to_svg(mol, legend="Input molecule")


def _pattern_mode(
    mol: Chem.Mol,
    patterns: dict[str, Chem.Mol | None],
    not_found_message: str,
) -> GalleryResult:
    images = highlight_by_patterns(mol, patterns)
    if not images:
        return [], not_found_message
    return images, f"Found {len(images)} match(es)."


def functional_groups(mol: Chem.Mol, _smiles: str, _min_ph: float, _max_ph: float) -> GalleryResult:
    return _pattern_mode(mol, functional_group_patterns(), "No functional groups recognized.")


def interligand_moieties(mol: Chem.Mol, _smiles: str, _min_ph: float, _max_ph: float) -> GalleryResult:
    return _pattern_mode(
        mol,
        interligand_patterns(),
        "No interligand moieties recognized.",
    )


def daylight_smarts_examples(mol: Chem.Mol, _smiles: str, _min_ph: float, _max_ph: float) -> GalleryResult:
    return _pattern_mode(
        mol,
        daylight_smarts_patterns(),
        "No SMARTS examples recognized.",
    )


def smartsrx_moieties(mol: Chem.Mol, _smiles: str, _min_ph: float, _max_ph: float) -> GalleryResult:
    return _pattern_mode(
        mol,
        smartsrx_patterns(),
        "No SMARTS-RX moieties recognized.",
    )


def get_rotatable_bond_indices(mol: Chem.Mol) -> list[int]:
    rotatable_bonds: list[int] = []
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue
        if bond.IsInRing():
            continue

        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        if begin_atom.GetAtomicNum() == 1 or end_atom.GetAtomicNum() == 1:
            continue
        if begin_atom.GetDegree() < 2 or end_atom.GetDegree() < 2:
            continue

        rotatable_bonds.append(bond.GetIdx())
    return rotatable_bonds


def rotatable_bonds(mol: Chem.Mol, _smiles: str, _min_ph: float, _max_ph: float) -> GalleryResult:
    images: list[tuple[str, str]] = []
    bond_indices = get_rotatable_bond_indices(mol)
    if bond_indices:
        images.append(
            (
                mol_to_svg(
                    mol,
                    highlight_bonds=bond_indices,
                    legend="Rotatable Bonds",
                ),
                "Local algorithm",
            )
        )

    images.extend(highlight_by_patterns(mol, rotatable_patterns()))
    if not images:
        return [], "No rotatable bonds found."
    return images, "Rotatable bonds highlighted."


def chiral_centers(mol: Chem.Mol, _smiles: str, _min_ph: float, _max_ph: float) -> GalleryResult:
    centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    if not centers:
        return [], "No chiral centers recognized."

    atom_labels = {idx: chirality for idx, chirality in centers}
    highlight_atoms = list(atom_labels)
    legend = "Chiral Centers: " + ", ".join(f"{idx} ({ch})" for idx, ch in centers)
    return [
        (
            mol_to_svg(
                mol,
                highlight_atoms=highlight_atoms,
                legend=legend,
                atom_labels=atom_labels,
            ),
            "Chiral centers",
        )
    ], "Chiral centers highlighted."


def stereocenters(mol: Chem.Mol, _smiles: str, _min_ph: float, _max_ph: float) -> GalleryResult:
    centers = Chem.FindPotentialStereo(mol)
    if not centers:
        return [], "No potential stereo centers found."

    highlight_atoms = [center.centeredOn for center in centers]
    atom_labels = {center.centeredOn: center.type.name for center in centers}
    return [
        (
            mol_to_svg(
                mol,
                highlight_atoms=highlight_atoms,
                legend="Potential Stereogenic Centers",
                atom_labels=atom_labels,
            ),
            "Potential Stereogenic Centers",
        )
    ], f"Found {len(centers)} potential stereogenic center(s)."


def scaffold(mol: Chem.Mol, _smiles: str, _min_ph: float, _max_ph: float) -> GalleryResult:
    scaffold_mol = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold_mol is None:
        return [], "No scaffold found."

    match = mol.GetSubstructMatch(scaffold_mol)
    if not match:
        return [], "Scaffold not found as substructure."

    return [
        (
            mol_to_svg(
                mol,
                highlight_atoms=list(match),
                highlight_bonds=matching_bond_indices(mol, match),
                legend="Murcko Scaffold",
            ),
            "Murcko Scaffold",
        )
    ], "Scaffold highlighted."


def hybridization(mol: Chem.Mol, _smiles: str, _min_ph: float, _max_ph: float) -> GalleryResult:
    atom_labels: dict[int, str] = {}
    for atom in mol.GetAtoms():
        hybridization_type = atom.GetHybridization()
        if hybridization_type != Chem.HybridizationType.UNSPECIFIED:
            atom_labels[atom.GetIdx()] = hybridization_type.name

    if not atom_labels:
        return [], "No hybridization states to display."

    return [
        (
            mol_to_svg(
                mol,
                highlight_atoms=list(atom_labels),
                legend="Hybridization States",
                atom_labels=atom_labels,
            ),
            "Hybridization States",
        )
    ], "Hybridization states highlighted."


def gasteiger_charges(mol: Chem.Mol, _smiles: str, _min_ph: float, _max_ph: float) -> GalleryResult:
    mol_with_hydrogens = Chem.AddHs(mol)
    AllChem.ComputeGasteigerCharges(mol_with_hydrogens)

    atom_labels = {
        atom.GetIdx(): f"{atom.GetDoubleProp('_GasteigerCharge'):.3f}"
        for atom in mol_with_hydrogens.GetAtoms()
    }
    if not atom_labels:
        return [], "Could not compute Gasteiger charges."

    return [
        (
            mol_to_svg(
                mol_with_hydrogens,
                highlight_atoms=list(atom_labels),
                legend="Gasteiger Charges (including H)",
                atom_labels=atom_labels,
            ),
            "Gasteiger Charges",
        )
    ], "Gasteiger charges computed and displayed (including hydrogens)."


def protonate_ph(_mol: Chem.Mol, smiles: str, min_ph: float, max_ph: float) -> GalleryResult:
    if min_ph > max_ph:
        return [], "Minimum pH must be less than or equal to maximum pH."

    try:
        protonated_smiles = dimorphite_dl.protonate_smiles(
            smiles,
            ph_min=min_ph,
            ph_max=max_ph,
            precision=1.0,
        )
    except Exception as exc:
        return [], f"Error during protonation: {exc}"

    if not protonated_smiles:
        return [], "No protonation variants found."

    images: list[tuple[str, str]] = []
    for index, variant_smiles in enumerate(protonated_smiles, start=1):
        variant_mol = _mol_from_smiles(variant_smiles)
        if variant_mol is None:
            continue

        images.append(
            (
                mol_to_svg(
                    variant_mol,
                    IMAGE_SIZE,
                    legend=f"Protonated variant {index} at pH {min_ph}-{max_ph}",
                ),
                f"Variant {index}",
            )
        )

    if not images:
        return [], "No valid protonation variants found."
    return images, f"Found {len(images)} protonation variant(s) at pH {min_ph}-{max_ph}."


FEATURE_MODES = [
    FeatureMode("Functional Groups", functional_groups),
    FeatureMode("Rotatable Bonds", rotatable_bonds),
    FeatureMode("Interligand Moieties", interligand_moieties),
    FeatureMode("SMARTS-RX Moieties", smartsrx_moieties),
    FeatureMode("Chiral Centers", chiral_centers),
    FeatureMode("Potential Stereogenic Centers", stereocenters),
    FeatureMode("DAYLIGHT SMARTS Examples", daylight_smarts_examples),
    FeatureMode("Murcko Scaffold", scaffold),
    FeatureMode("Hybridization", hybridization),
    FeatureMode("Gasteiger Charges", gasteiger_charges),
    FeatureMode("Protonation", protonate_ph, needs_ph=True),
]
FEATURE_MODE_BY_NAME = {mode.name: mode for mode in FEATURE_MODES}


def process_smiles_main(
    smiles: str,
    mode: str,
    min_ph: float = 7.0,
    max_ph: float = 7.0,
) -> GalleryResult:
    mol = _mol_from_smiles(smiles)
    if mol is None:
        return [], "Invalid SMILES."

    feature_mode = FEATURE_MODE_BY_NAME.get(mode)
    if feature_mode is None:
        return [], "Invalid mode selected."

    return feature_mode.handler(mol, smiles, min_ph, max_ph)
