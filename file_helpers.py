# Description: Helper functions for loading files and data.
from rdkit import Chem


def load_interligand_moieties():
    moieties = {}
    try:
        with open("data/SMARTS_InteLigand.txt", "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if ":" in line:
                    parts = line.split(":", 1)
                    name = parts[0].strip()
                    pattern = parts[1].strip()
                    if pattern and pattern.startswith("["):
                        moieties[name] = pattern
    except Exception as e:
        print("Error loading SMARTS_InteLigand.txt:", e)
    return moieties


def load_smarts_patterns_from_csv():
    """
    Load and compile SMARTS from the CSV file.
    """
    import csv
    compiled_patterns = {}
    try:
        with open("data/smarts_examples_parsed.csv", "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row['Pattern']:  # Skip empty patterns
                    name = f"{row['Main Topic']} > {row['Subtopic']} > {row['Sub-sub-topic']} > {row['Rule Name']}"
                    pattern = row['Pattern'].strip()
                    try:
                        mol = Chem.MolFromSmarts(pattern)
                        if mol:  # Only add if pattern compiles successfully
                            compiled_patterns[name] = mol
                    except:
                        print(f"Failed to compile SMARTS pattern: {pattern}")
    except Exception as e:
        print("Error loading smarts_examples_parsed.csv:", e)
    return compiled_patterns
