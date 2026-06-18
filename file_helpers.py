"""Compatibility wrappers for SMARTS loading helpers."""

from moliety2.patterns import daylight_smarts_patterns, interligand_smarts


def load_interligand_moieties():
    return interligand_smarts()


def load_smarts_patterns_from_csv():
    return daylight_smarts_patterns()
