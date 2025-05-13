import os
import pandas as pd
from snakemake.utils import validate
import glob
import yaml

# Input validation
if not os.path.exists(config["cells_dna"]) and not os.path.exists(config["cells_rna"]):
    raise ValueError("Neither RNA or DNA seedfiles provided. Stopping.")

# Wildcard constraints
wildcard_constraints:
    binsize="\d+"

# Config validation function
def validate_config():
    required_keys = [
        "outdir",
        "path",
        "ref"
    ]
    for key in required_keys:
        if key not in config:
            raise KeyError(f"Required config key '{key}' missing")
            
    # Validate DNA config if DNA analysis enabled
    if "cells_dna" in config:
        required_dna_keys = [
            "binsize",
            "min_mapq",
            "max_insert",
            "min_count",
            "gamma",
            "dim_reduction",
            "normalize_to_panel"
        ]
        for key in required_dna_keys:
            if key not in config["dna"]:
                raise KeyError(f"Required DNA config key 'dna.{key}' missing")
        
        required_ref_keys = [
            "chr_list",
            "genome_bed",
            "blacklist",
            "maptrack"
        ]
        for key in required_ref_keys:
            if key not in config["ref"]:
                raise KeyError(f"Required reference config key 'ref.{key}' missing")
                
    # Validate RNA config if RNA analysis enabled        
    if "cells_rna" in config:
        required_rna_keys = [
            "min_count",
            "min_genes"
        ]
        for key in required_rna_keys:
            if key not in config["rna"]:
                raise KeyError(f"Required RNA config key 'rna.{key}' missing")

validate_config()

# Helper functions
# DNA 
def is_dna_analysis():
    return os.path.exists(config["cells_dna"])

def get_patients_dna():
    if is_dna_analysis():
        return list(dict.fromkeys(cells_dna['patient_id'].dropna()))
    return []

def get_cells_dna(patient_id):
    """Get unique cell IDs for a specific patient."""
    if not is_dna_analysis():
        return []
    cells = cells_dna['cell'][cells_dna.patient_id == patient_id]
    return list(dict.fromkeys(cells))

def get_dna_fq1(wildcards):
    # Check patient cells first
    if cells_dna is not None and wildcards.cell in cells_dna.cell.values:
        return cells_dna[cells_dna.cell == wildcards.cell]["fq1"].values
    # Check normal cells
    elif cells_normals is not None and wildcards.cell in cells_normals.cell.values:
        return cells_normals[cells_normals.cell == wildcards.cell]["fq1"].values
    else:
        raise ValueError(f"Cell {wildcards.cell} not found in patient or normal cells")

def get_dna_fq2(wildcards):
    # Check patient cells first  
    if cells_dna is not None and wildcards.cell in cells_dna.cell.values:
        return cells_dna[cells_dna.cell == wildcards.cell]["fq2"].values
    # Check normal cells
    elif cells_normals is not None and wildcards.cell in cells_normals.cell.values:
        return cells_normals[cells_normals.cell == wildcards.cell]["fq2"].values
    else:
        raise ValueError(f"Cell {wildcards.cell} not found in patient or normal cells")

# RNA
def is_rna_analysis():
    return os.path.exists(config["cells_rna"])

def has_rna_data(patient_id):
    return cells_rna is not None and patient_id in cells_rna['patient_id'].values

def get_patients_rna():
    if is_rna_analysis():
        return list(dict.fromkeys(cells_rna['patient_id'].dropna()))
    return []

def get_cells_rna(patient_id):
    """Get unique cell IDs for a specific patient."""
    if not is_rna_analysis():
        return []
    cells = cells_rna['cell'][cells_rna.patient_id == patient_id]
    return list(dict.fromkeys(cells))

def get_rna_fq1(wildcards):
    return cells_rna[cells_rna.cell==wildcards.cell]["fq1"]

def get_rna_fq2(wildcards):
    return cells_rna[cells_rna.cell==wildcards.cell]["fq2"]

# Loading functions
def load_sample_data():
    data = {"dna": None, "rna": None}
    
    if is_dna_analysis():
        data["dna"] = pd.read_table(
            config["cells_dna"], 
            sep="\t", 
            header=0
        )
        
        # Add patient IDs
        if os.path.exists(config["metadata"]):
            meta = pd.read_table(config["metadata"], header=0)
            meta = meta[['plate_id', 'patient_id']].drop_duplicates()
            data["dna"]['plate_id'] = data["dna"].cell.str[:8]
            data["dna"] = pd.merge(data["dna"], meta, on='plate_id', how='left')
            data["dna"]['patient_id'] = data["dna"]['patient_id'].fillna(
                data["dna"].cell.str[:6]
            )
        else:
            data["dna"]['patient_id'] = data["dna"].cell.str[:6]
            
    if is_rna_analysis():
        data["rna"] = pd.read_table(
            config["cells_rna"],
            sep="\t",
            header=0
        ).set_index("cell", drop=False)
        
        # Add patient IDs similarly for RNA
        if os.path.exists(config["metadata"]):  
            if "meta" not in locals():
                meta = pd.read_table(config["metadata"], header=0)
                meta = meta[['plate_id', 'patient_id']].drop_duplicates()
            data["rna"]['plate_id'] = data["rna"].cell.str[:8]
            data["rna"] = pd.merge(data["rna"], meta, on='plate_id', how='left')
            data["rna"]['patient_id'] = data["rna"]['patient_id'].fillna(
                data["rna"].cell.str[:6]
            )
        else:
            data["rna"]['patient_id'] = data["rna"].cell.str[:6]

    # Filter by target patient if defined
    if "target_patient" in config:
        if data["dna"] is not None:
            data["dna"] = data["dna"][data["dna"]['patient_id'].isin(config["target_patient"])]
        if data["rna"] is not None:
            data["rna"] = data["rna"][data["rna"]['patient_id'].isin(config["target_patient"])]

        # First check: warn if any patient doesn't exist in either DNA or RNA
        found_patients = set()
        if data["dna"] is not None:
            found_patients.update(data["dna"]['patient_id'].unique())
        if data["rna"] is not None:
            found_patients.update(data["rna"]['patient_id'].unique())
            
        if len(found_patients) != len(config["target_patient"]):
            missing_patients = set(config["target_patient"]) - found_patients
            logger.warning(f"No data found for patients: {missing_patients}")

        # Second check: ensure at least one dataset has data
        if (data["dna"] is None or len(data["dna"]) == 0) and (data["rna"] is None or len(data["rna"]) == 0):
            raise ValueError(f"No DNA or RNA data found for patients {config['target_patient']}")

    # Print summaries
    if data["dna"] is not None:
        patient_summary = data["dna"].groupby('patient_id').size()
        print("\n----------------------------------------")
        print("\033[1mDNA samples\033[0m")
        for patient, count in patient_summary.items():
            print(f"  {patient}: {count} cells")
        print(f"Total: {len(data['dna'])} cells")
            
    if data["rna"] is not None:
        patient_summary = data["rna"].groupby('patient_id').size()
        print("----------------------------------------")
        print("\033[1mRNA samples\033[0m")
        for patient, count in patient_summary.items():
            print(f"  {patient}: {count} cells")
        print(f"Total: {len(data['rna'])} cells")
    
    return data

# Load normal cells if specified
def load_normal_data():
    normal_dfs = []
    if is_dna_analysis():
        normal_files = [config["dna"][key] for key in ["normals_bins", "normals_scaling"] 
                       if key in config["dna"] and os.path.exists(config["dna"][key])]
        for file in set(normal_files):
            df = pd.read_table(file, sep="\t", header=0)
            if not {'fq1', 'fq2'}.issubset(df.columns):
                raise ValueError(f"Normal panel file {file} must contain 'fq1' and 'fq2' columns. Need to run find_normal_fq.py first?")
            normal_dfs.append(df)

    if normal_dfs:
        normals = pd.concat(normal_dfs).drop_duplicates()
        print("----------------------------------------")
        print(f"Normal panel: {len(normals)} cells")
        return normals
            
    return None

# Load sample data
sample_data = load_sample_data()
cells_dna = sample_data["dna"]
cells_rna = sample_data["rna"]

# Get normal cell panels
cells_normals = load_normal_data()
# 1. For outlier bins 
if is_dna_analysis() and os.path.exists(config["dna"]["normals_bins"]):
    normals_bins = pd.read_table(config["dna"]["normals_bins"], header=0)
    normals_bins_id = os.path.splitext(os.path.basename(config["dna"]["normals_bins"]))[0]
else:
    normals_bins = pd.DataFrame()
    normals_bins_id = None

# 2. For normal cell scaling
if is_dna_analysis() and os.path.exists(config["dna"]["normals_scaling"]):
    normals_scaling = pd.read_table(config["dna"]["normals_scaling"], header=0)
    normals_scaling_id = os.path.splitext(os.path.basename(config["dna"]["normals_scaling"]))[0]
else:
    normals_scaling = pd.DataFrame()
    normals_scaling_id = None

# Load patient-specific parameters
patient_params_file = config["dna"]["patient_params"] if "patient_params" in config["dna"] else None
patient_params = {}
if os.path.exists(patient_params_file):
    with open(patient_params_file) as f:
        patient_params = yaml.safe_load(f)

# def get_patient_param(patient_id, param_name):
#     if patient_id in patient_params and param_name in patient_params[patient_id]:
#         return patient_params[patient_id][param_name]
#     return config["dna"][param_name]

def get_patient_param(patient_id, param_name):
    """Get parameter value for patient, falling back to config default if not specified."""
    value = patient_params.get(patient_id, {}).get(param_name, config["dna"][param_name])
    return value if isinstance(value, list) else [value]
