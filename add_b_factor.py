
import os
from Bio.PDB import PDBParser
import pandas as pd
from tqdm import tqdm

DATA_FOLDER = "/nfs/scistore20/bronsgrp/jsaleh/data"

def extract_ca_bfactors(pdb_file, chain_id=None):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(os.path.basename(pdb_file), pdb_file)

    b_factors = []
    for model in structure:
        for chain in model:
            if chain_id and chain.id != chain_id:
                continue

            for residue in chain:
                if 'CA' in residue:
                    ca_atom = residue['CA']
                    b_factors.append(ca_atom.get_bfactor())
    return b_factors


def percentage_high_bfactor(b_factors, threshold=75.0):
    if len(b_factors) == 0:
        return 0
    count_above_threshold = sum(b >= threshold for b in b_factors)
    percentage = (count_above_threshold / len(b_factors)) * 100
    return percentage


def analyze_pdb_directory(folder_path, chain_id=None, bfactor_threshold=75.0, percent_cutoff=80.0, af_model = 0):
    pdb_results = []

    for folder in tqdm(os.listdir(folder_path)):

        bmrb_id = folder.split("_")[0] if folder_path.endswith("v3") else folder
        af_file = f"af2_0{af_model}.pdb" if folder_path.endswith("v3") else f"af_0{af_model}.pdb"

        pdb_file_path = os.path.join(folder_path, folder, af_file)
        if pdb_file_path.lower().endswith('.pdb'):
            pdb_path = os.path.join(folder_path, pdb_file_path)
            b_factors = extract_ca_bfactors(pdb_path, chain_id)
            perc_high_bfactor = percentage_high_bfactor(b_factors, bfactor_threshold)

            meets_criteria = perc_high_bfactor >= percent_cutoff

            pdb_results.append({
                'bmrb_id': int(bmrb_id),
                'percentage_high_bfactor': perc_high_bfactor,
                'meets_criteria': meets_criteria
            })

    return pd.DataFrame(pdb_results)


# Example usage:
if __name__ == "__main__":
    folder_path = f'{DATA_FOLDER}/af-bmrb-h-v4'
    chain_of_interest = None  # Set to 'A', 'B', etc., or None for all chains
    try:
        results_df = pd.read_csv(f'{DATA_FOLDER}/bfactor_analysis_v4.csv')
    except FileNotFoundError:
        results_df = analyze_pdb_directory(
            folder_path=folder_path,
            chain_id=chain_of_interest,
            bfactor_threshold=75.0,
            percent_cutoff=80.0
        )
        results_df.to_csv(f'{DATA_FOLDER}/bfactor_analysis_v4.csv')

    # Filter proteins meeting criteria
    proteins_meeting_criteria = results_df[results_df['meets_criteria']]
    print("\nProteins with ≥80% residues having B-factor ≥75:")
    print(proteins_meeting_criteria[['bmrb_id', 'percentage_high_bfactor']])
