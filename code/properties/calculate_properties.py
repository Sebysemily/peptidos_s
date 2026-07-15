#!/usr/bin/env python3
import csv
from pathlib import Path
import warnings

# ==========================================
# BLOQUE DE COMPATIBILIDAD JUPYTER / SNAKEMAKE
# ==========================================
if "snakemake" not in locals():

    class MockSnakemake:
        input = [""]
        output = ["results/plot_prueba.png"]
        wildcards = type("Wildcards", (object,), {"muestra": "prueba"})()
        params = {"threshold": 0.05}

    snakemake = MockSnakemake()
    print("⚠️ Modo Interactivo (Molten/Jupyter) detectado. Usando variables Mock.")
# ==========================================

try:
    from modlamp.descriptors import GlobalDescriptor, PeptideDescriptor
except ImportError:
    import subprocess
    import sys
    print("Modlamp no encontrado. Instalando modlamp sin dependencias (--no-deps)...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "--no-deps", "modlamp==4.3.0", "--quiet"])
    from modlamp.descriptors import GlobalDescriptor, PeptideDescriptor


def read_fasta_records(path):
    records = []
    header = None
    seq_lines = []
    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    sequence = "".join(part.strip() for part in seq_lines)
                    if sequence:
                        records.append((header[1:], sequence))
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
    if header is not None:
        sequence = "".join(part.strip() for part in seq_lines)
        if sequence:
            records.append((header[1:], sequence))
    return records


def calculate_all_properties(sequences):
    results = {}
    
    # --- Global Descriptors ---
    print("Calculating GlobalDescriptors...")
    gd = GlobalDescriptor(sequences)
    
    gd.calculate_charge(ph=7.0, amide=False)
    results['net_charge'] = [row[0] for row in gd.descriptor]
    
    gd.calculate_MW()
    results['molecular_weight'] = [row[0] for row in gd.descriptor]
    
    gd.charge_density()
    results['charge_density'] = [row[0] for row in gd.descriptor]
    
    gd.isoelectric_point()
    results['isoelectric_point'] = [row[0] for row in gd.descriptor]
    
    gd.instability_index()
    results['instability_index'] = [row[0] for row in gd.descriptor]
    
    gd.aromaticity()
    results['aromaticity'] = [row[0] for row in gd.descriptor]
    
    gd.aliphatic_index()
    results['aliphatic_index'] = [row[0] for row in gd.descriptor]
    
    gd.boman_index()
    results['boman_index'] = [row[0] for row in gd.descriptor]
    
    gd.hydrophobic_ratio()
    results['hydrophobic_ratio'] = [row[0] for row in gd.descriptor]
    
    # --- Peptide Descriptors ---
    scales = [
        'AASI', 'ABHPRK', 'argos', 'bulkiness', 'charge_phys', 'charge_acid', 
        'cougar', 'eisenberg', 'Ez', 'flexibility', 'grantham', 'gravy', 
        'hopp-woods', 'ISAECI', 'janin', 'kytedoolittle', 'levitt_alpha', 
        'MSS', 'MSW', 'pepArc', 'pepcats', 'polarity', 'PPCALI', 
        'refractivity', 't_scale', 'TM_tend'
    ]
    
    print("Calculating PeptideDescriptors...")
    # Suppress warnings for scales that crash on moment
    warnings.filterwarnings("ignore")
    
    for scale in scales:
        pd = PeptideDescriptor(sequences, scale)
        
        # Global mean
        pd.calculate_global(window=1000, modality='mean')
        results[f'global_{scale}'] = [row[0] for row in pd.descriptor]
        
        # Max moment
        try:
            pd.calculate_moment(window=1000, angle=100, modality='max')
            results[f'moment_{scale}'] = [row[0] for row in pd.descriptor]
        except Exception:
            # Multi-dimensional scales will crash here
            results[f'moment_{scale}'] = [None] * len(sequences)
            
    return results


def main():
    fastas = snakemake.input
    if isinstance(fastas, str):
        fastas = [fastas]

    all_records = []
    for fasta in fastas:
        all_records.extend(read_fasta_records(fasta))

    if not all_records:
        print("No sequences found.")
        return

    headers = [rec[0] for rec in all_records]
    sequences = [rec[1] for rec in all_records]

    print(f"Calculating properties for {len(sequences)} sequences...")
    properties_dict = calculate_all_properties(sequences)
    
    # Prepare CSV columns
    columns = ["peptide_id", "sequence", "length"] + list(properties_dict.keys())
    
    output_csv = Path(snakemake.output[0])
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    print(f"Writing {len(columns)} columns to {output_csv}...")
    with output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        
        for i in range(len(sequences)):
            row = [headers[i], sequences[i], len(sequences[i])]
            for col in properties_dict.keys():
                val = properties_dict[col][i]
                if val is not None:
                    row.append(round(val, 4))
                else:
                    row.append("NA")
            writer.writerow(row)

    print(f"Saved properties to {output_csv}")


if __name__ == "__main__":
    main()
