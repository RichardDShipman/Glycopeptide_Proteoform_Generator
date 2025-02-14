#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Oct 24 2024

Glycopeptide Proteoform Generator

Author: richarddshipman
"""

import itertools
from pathlib import Path
from collections import defaultdict
import pandas as pd
import argparse
from datetime import datetime
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed

# Set up logging configuration
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

def generate_proteoforms_with_limit(protein_data: dict, limit: int = 100) -> list:
    """
    Generate unique proteoform combinations for a given protein.
    
    Args:
        protein_data (dict): A dictionary where keys are glycosylation sites and values 
                             are lists (of lists) containing glyco options (tuples).
        limit (int): Maximum number of proteoforms to generate.
    
    Returns:
        list: A list of unique proteoform combinations.
    """
    protein_forms = []
    for site, glyco_options in protein_data.items():
        # Each 'glyco_options' is a list-of-lists; use product to get all combinations for this site.
        combinations = list(itertools.product(*glyco_options))
        protein_forms.append(combinations)
    
    # Generate the Cartesian product across all glycosylation sites.
    proteoforms_iter = itertools.product(*protein_forms)
    
    limited_proteoforms = set()
    count = 0
    for proteoform in proteoforms_iter:
        if count >= limit:
            break
        limited_proteoforms.add(proteoform)
        count += 1
    
    return list(limited_proteoforms)

def process_protein(protein: str, protein_data: dict, limit: int, base_output_dir: Path) -> (str, int):
    """
    Process a single protein: generate proteoforms and write details to a text file.
    
    Args:
        protein (str): Protein identifier.
        protein_data (dict): Glycopeptide data for the protein.
        limit (int): Maximum number of proteoforms to generate.
        base_output_dir (Path): Output directory where the protein file is saved.
    
    Returns:
        tuple: (protein, number of proteoforms generated)
    """
    proteoforms = generate_proteoforms_with_limit(protein_data, limit)
    output_file = base_output_dir / f"{protein}_proteoforms.txt"
    try:
        with output_file.open("w") as file:
            for idx, proteoform in enumerate(proteoforms, 1):
                glyco_sites = []
                # Each proteoform is a tuple with one element per site
                for site_combo in proteoform:
                    for glycosylation_site, glycan in site_combo:
                        # Attempt to convert the glycosylation site to int for sorting; if not, keep as is.
                        try:
                            site_int = int(glycosylation_site)
                        except (ValueError, TypeError):
                            site_int = glycosylation_site
                        glyco_sites.append((site_int, glycan))
                glyco_sites.sort(key=lambda x: x[0])
                formatted_sites = ":".join([f"{site}-{glycan}" for site, glycan in glyco_sites])
                file.write(f"{protein}_PF_{idx}, {formatted_sites}\n")
    except Exception as e:
        logging.error(f"Error writing file for protein {protein}: {e}")
    
    return protein, len(proteoforms)

def main(input_file: str, limit: int, protein_col: str, site_col: str, glycan_col: str) -> None:
    """
    Main function to generate glycopeptide proteoforms from an input CSV file.
    
    The function performs the following steps:
      1. Reads the input CSV file into a pandas DataFrame.
      2. Processes rows to create a dictionary grouping glycopeptides by protein and glycosylation site.
      3. Generates proteoform combinations for each protein with a specified limit.
      4. Writes individual proteoform details and counts to output files.
      5. Merges individual protein files into a single CSV.
      6. Creates log and summary files.
    
    Args:
        input_file (str): Path to the input CSV file.
        limit (int): Maximum number of proteoforms per protein.
        protein_col (str): Name of the protein column in the CSV.
        site_col (str): Name of the glycosylation site column.
        glycan_col (str): Name of the glycan column.
    """
    # Read the CSV file
    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        logging.error(f"Error reading input file {input_file}: {e}")
        return

    # Group glycopeptides by protein and glycosylation site
    glycopeptides = defaultdict(lambda: defaultdict(set))
    for _, row in df.iterrows():
        protein = row[protein_col]
        glycosylation_site = row[site_col]
        glycan = row[glycan_col]
        glycopeptides[protein][glycosylation_site].add(glycan)
    
    # Format glycopeptide data for processing
    formatted_glycopeptides = [
        {
            'protein': protein,
            'glycosylation_site': site,
            'glycans': list(glycans)
        }
        for protein, sites in glycopeptides.items()
        for site, glycans in sites.items()
    ]
    
    # Build a dictionary mapping each protein to its glycopeptide data.
    protein_dict = defaultdict(lambda: defaultdict(list))
    for glycopeptide in formatted_glycopeptides:
        protein = glycopeptide["protein"]
        glycosylation_site = glycopeptide["glycosylation_site"]
        glycans = glycopeptide["glycans"]
        # Include an option for no glycosylation (None) plus each glycan.
        glyco_combinations = [(glycosylation_site, None)] + [(glycosylation_site, glycan) for glycan in glycans]
        # Use list(set(...)) to ensure uniqueness.
        protein_dict[protein][glycosylation_site].append(list(set(glyco_combinations)))
    
    # Create an output directory using the base name of the input file.
    input_basename = Path(input_file).stem
    base_output_dir = Path("data") / input_basename
    try:
        base_output_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        logging.error(f"Error creating output directory {base_output_dir}: {e}")
        return
    
    # Process each protein in parallel.
    counts = {}
    with ProcessPoolExecutor() as executor:
        future_to_protein = {
            executor.submit(process_protein, protein, protein_data, limit, base_output_dir): protein
            for protein, protein_data in protein_dict.items()
        }
        for future in as_completed(future_to_protein):
            protein = future_to_protein[future]
            try:
                protein, total_proteoforms = future.result()
                counts[protein] = total_proteoforms
                logging.info(f"{protein}: Total number of proteoforms: {total_proteoforms}")
            except Exception as e:
                logging.error(f"Error processing protein {protein}: {e}")
    
    # Write proteoform counts to a CSV file.
    counts_file_path = base_output_dir / f"00_proteoform_counts_{input_basename}.csv"
    try:
        with counts_file_path.open("w") as counts_file:
            counts_file.write("protein,total_proteoforms\n")
            for protein, total in counts.items():
                counts_file.write(f"{protein},{total}\n")
    except Exception as e:
        logging.error(f"Error writing counts file: {e}")
    
    logging.info("Proteoform counts have been written to the output directory.")
    
    # Merge individual proteoform text files into a single CSV file.
    merged_proteoforms_path = base_output_dir / f"01_merged_proteoforms_{input_basename}.csv"
    try:
        with merged_proteoforms_path.open("w") as merged_file:
            merged_file.write("protein,proteoform_id,glycosylation_sites\n")
            for protein in protein_dict.keys():
                proteoform_file_path = base_output_dir / f"{protein}_proteoforms.txt"
                if proteoform_file_path.exists():
                    with proteoform_file_path.open("r") as file:
                        for line in file:
                            if line.strip():
                                proteoform_id, glycosylation_sites = line.split(', ', 1)
                                merged_file.write(f"{protein},{proteoform_id},{glycosylation_sites.strip()}\n")
    except Exception as e:
        logging.error(f"Error merging proteoform files: {e}")
    
    # Create a log file and summary CSV.
    log_file_path = base_output_dir / f"02_input_log_{Path(input_file).name}.txt"
    try:
        with log_file_path.open("w") as log_file:
            log_file.write("Glycopeptide Proteoform Generator\n")
            log_file.write(f"Run Date and Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            log_file.write("\nInput Parameters:\n")
            log_file.write(f"Input File: {input_file}\n")
            log_file.write(f"Limit: {limit}\n")
            log_file.write(f"Protein Column: {protein_col}\n")
            log_file.write(f"Glycosylation Site Column: {site_col}\n")
            log_file.write(f"Glycan Column: {glycan_col}\n")
            log_file.write("\nSummary:\n")
            total_proteins = len(protein_dict)
            log_file.write(f"Total number of proteins: {total_proteins}\n")
            total_proteoforms = sum(counts.values())
            log_file.write(f"Total number of proteoforms generated: {total_proteoforms}\n")
            
            summary_data = []
            for protein, protein_data in protein_dict.items():
                total_sites = len(protein_data)
                total_glycans = sum(len(glyco_options) for glyco_options in protein_data.values())
                proteoforms = generate_proteoforms_with_limit(protein_data, limit)
                summary_data.append((protein, total_sites, total_glycans, len(proteoforms)))
            summary_df = pd.DataFrame(summary_data, 
                                      columns=["Protein", "TotalGlycosylationSites", "TotalGlycans", "TotalProteoforms"])
            summary_csv_path = base_output_dir / f"03_summary_{Path(input_file).name}.csv"
            summary_df.to_csv(summary_csv_path, index=False)
    except Exception as e:
        logging.error(f"Error writing log/summary files: {e}")
    
    logging.info("Log and summary files have been written to the output directory.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate proteoforms from glycopeptides data.')
    parser.add_argument('-i', '--input', default='human_proteoform_glycosylation_sites_gptwiki.csv', # Test data from GlyGen
                        help='Path to the input CSV file.')
    parser.add_argument('-l', '--limit', type=int, default=10,
                        help='Maximum number of proteoforms per protein.')
    parser.add_argument('-p', '--protein', default='protein', help='Protein column name.')
    parser.add_argument('-s', '--site', default='glycosylation_site', help='Glycosylation site column name.')
    parser.add_argument('-g', '--glycan', default='glycan', help='Glycan column name.')
    
    args = parser.parse_args()
    main(args.input, args.limit, args.protein, args.site, args.glycan)
    