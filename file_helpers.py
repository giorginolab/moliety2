import yaml
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

def load_yaml_smarts():
    """
    Load and compile SMARTS from the YAML file.
    """
    try:
        with open("data/daylight_smarts.yml", "r") as f:
            data = yaml.safe_load(f)
    except Exception as e:
        print("Error loading daylight_smarts.yml:", e)
        return {}
    compiled_yaml = {}
    for group in data.get("groups", []):
        group_name = group.get("name", "Unnamed Group")
        for subgroup in group.get("subgroups", []):
            subgroup_name = subgroup.get("name", "Unnamed Subgroup")
            if "subsubgroups" in subgroup:
                for subsub in subgroup.get("subsubgroups", []):
                    subsub_name = subsub.get("name", "Unnamed Subsubgroup")
                    for rule in subsub.get("rules", []):
                        if "smarts" in rule:
                            key = f"{group_name} > {subgroup_name} > {subsub_name} > {rule.get('name', 'Unnamed Rule')}"
                            compiled_yaml[key] = Chem.MolFromSmarts(rule.get("smarts"))
            elif "rules" in subgroup:
                for rule in subgroup.get("rules", []):
                    if "smarts" in rule:
                        key = f"{group_name}: {subgroup_name} - {rule.get('name', 'Unnamed Rule')}"
                        compiled_yaml[key] = Chem.MolFromSmarts(rule.get("smarts"))
    return compiled_yaml
