#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Oct 24 2024

Glycopeptide Proteoform Generator

Author: richarddshipman, Review with: Gemini 2.5 Pro and ChatGPT 4o
"""

import itertools
from pathlib import Path
from collections import defaultdict
import pandas as pd
import argparse
from datetime import datetime
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
import traceback # Added for detailed error logging

# Set up logging configuration
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

def generate_proteoforms_with_limit(protein_data: dict, limit: int = 100) -> list:
    """
    Generate unique proteoform combinations for a given protein using a more
    efficient, refactored approach.
    
    Args:
        protein_data (dict): A dictionary where keys are glycosylation sites and values 
                             are lists containing one list of glyco options (tuples).
        limit (int): Maximum number of proteoforms to generate.
    
    Returns:
        list: A list of unique proteoform combinations.
    """
    # The values of protein_data are lists like [[(site1, G1), (site1, G2)], ...].
    # We extract the inner list of options for each site.
    option_lists = [options[0] for options in protein_data.values()]
    
    # Generate the Cartesian product of all glyco options across all sites.
    proteoforms_iter = itertools.product(*option_lists)
    
    # Use itertools.islice for a more direct and efficient way to enforce the limit.
    # This avoids iterating manually with a counter.
    limited_proteoforms = list(itertools.islice(proteoforms_iter, limit))
    
    # Using a set ensures uniqueness, which is good practice although itertools.product
    # already produces unique tuples if the input option lists are unique.
    return list(set(limited_proteoforms))

def process_protein(protein: str, protein_data: dict, limit: int) -> tuple[str, int, list]:
    """
    Process a single protein: generate proteoforms and prepare data for output.
    This function returns data for later aggregation, rather than writing to disk directly.
    
    Args:
        protein (str): Protein identifier.
        protein_data (dict): Glycopeptide data for the protein.
        limit (int): Maximum number of proteoforms to generate.
    
    Returns:
        tuple: (protein, number of proteoforms generated, list of formatted proteoform strings)
    """
    proteoforms_raw = generate_proteoforms_with_limit(protein_data, limit)
    formatted_proteoform_strings = []
    
    for idx, proteoform in enumerate(proteoforms_raw, 1):
        # Flatten the proteoform tuple structure to a list of (site, glycan) pairs.
        # A raw proteoform looks like: ((site1, glycanA), (site2, None))
        current_protein_glyco_states = list(proteoform)
        
        # Sort glycosylation sites for consistent output format.
        # The improved sort key handles non-integer site IDs more robustly.
        def sort_key(x):
            site = x[0]
            try:
                # Attempt to convert site to integer for numerical sorting.
                return int(site)
            except (ValueError, TypeError):
                # If conversion fails, convert to string for alphabetical sorting.
                return str(site)

        current_protein_glyco_states.sort(key=sort_key)
        
        # Format the sorted sites into the desired string format "site-glycan".
        formatted_sites = ":".join([f"{site}-{glycan if glycan is not None else 'None'}" for site, glycan in current_protein_glyco_states])
        formatted_proteoform_strings.append(f"{protein}_PF_{idx},{formatted_sites}")
    
    return protein, len(proteoforms_raw), formatted_proteoform_strings


def main(input_file: str, limit: int, protein_col: str, site_col: str, glycan_col: str) -> None:
    """
    Main function to generate glycopeptide proteoforms from an input CSV file.
    
    Performs the following steps:
      1. Reads the input CSV file into a pandas DataFrame.
      2. Groups glycopeptides by protein and glycosylation site.
      3. Generates proteoform combinations for each protein in parallel.
      4. Writes proteoform counts and merged proteoform details to output files.
      5. Creates a log and a summary file with run details.
    
    Args:
        input_file (str): Path to the input CSV file.
        limit (int): Maximum number of proteoforms per protein.
        protein_col (str): Name of the protein column in the CSV.
        site_col (str): Name of the glycosylation site column.
        glycan_col (str): Name of the glycan column.
    """
    start_time = time.time() # Record start time for performance measurement

    # Read the CSV file
    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        logging.error(f"Error reading input file {input_file}: {e}")
        return

    # Group glycopeptides by protein and site, storing unique glycans in a set.
    # Structure: {protein_id: {glycosylation_site_id: {glycan1, glycan2, ...}}, ...}
    glycopeptides_raw = defaultdict(lambda: defaultdict(set))
    for _, row in df.iterrows():
        protein = row[protein_col]
        glycosylation_site = row[site_col]
        glycan = row[glycan_col]
        glycopeptides_raw[protein][glycosylation_site].add(glycan)
    
    # Prepare the final dictionary for processing.
    # Structure: {protein: {site: [[(site, None), (site, glycan1), ...]]}}
    protein_dict_for_processing = defaultdict(lambda: defaultdict(list))
    for protein, sites_data in glycopeptides_raw.items():
        for glycosylation_site, glycans_set in sites_data.items():
            # For each site, create a list of all possible states: 
            # not glycosylated (None) or one of the observed glycans.
            glyco_options = [(glycosylation_site, None)] + \
                            [(glycosylation_site, glycan) for glycan in glycans_set]
            
            # The `itertools.product` function takes multiple iterables as arguments.
            # To get the Cartesian product of glyco options across different sites,
            # each site's list of options must be wrapped in another list.
            protein_dict_for_processing[protein][glycosylation_site].append(glyco_options)
    
    # Create an output directory named after the input file.
    input_basename = Path(input_file).stem
    base_output_dir = Path("data") / input_basename
    try:
        base_output_dir.mkdir(parents=True, exist_ok=True)
        logging.info(f"Output directory created: {base_output_dir}")
    except Exception as e:
        logging.error(f"Error creating output directory {base_output_dir}: {e}")
        return
    
    # Process each protein in parallel.
    counts = {}
    all_merged_proteoform_lines = []
    proteins_hitting_limit = []

    logging.info("Starting parallel processing of proteins...")
    with ProcessPoolExecutor() as executor:
        future_to_protein = {
            executor.submit(process_protein, protein, protein_data, limit): protein
            for protein, protein_data in protein_dict_for_processing.items()
        }
        for future in as_completed(future_to_protein):
            protein = future_to_protein[future]
            try:
                # Retrieve the result from the parallel process
                processed_protein, total_proteoforms, proteoform_strings = future.result()
                counts[processed_protein] = total_proteoforms
                all_merged_proteoform_lines.extend(proteoform_strings)
                logging.info(f"{processed_protein}: Generated {total_proteoforms} proteoforms.")

                # If the number of proteoforms generated equals the limit,
                # it's likely that more combinations were possible.
                if total_proteoforms == limit and limit > 0:
                    proteins_hitting_limit.append(processed_protein)

            except Exception as e:
                # Log error with a full traceback for better debugging.
                logging.error(f"Error processing protein {protein}: {e}\n{traceback.format_exc()}")
    
    logging.info("Finished parallel processing.")

    # Write proteoform counts to a CSV file.
    counts_file_path = base_output_dir / f"00_proteoform_counts_{input_basename}.csv"
    try:
        with counts_file_path.open("w") as counts_file:
            counts_file.write("protein,total_proteoforms\n")
            for protein, total in sorted(counts.items()):
                counts_file.write(f"{protein},{total}\n")
        logging.info(f"Proteoform counts saved to: {counts_file_path}")
    except Exception as e:
        logging.error(f"Error writing counts file: {e}")
    
    # Merge all proteoform data into a single CSV file.
    merged_proteoforms_path = base_output_dir / f"01_merged_proteoforms_{input_basename}.csv"
    try:
        with merged_proteoforms_path.open("w") as merged_file:
            merged_file.write("protein,proteoform_id,glycosylation_sites\n")
            # Sort lines for consistent output regardless of processing order
            for line in sorted(all_merged_proteoform_lines):
                parts = line.split(',', 1)
                if len(parts) == 2:
                    proteoform_id = parts[0].strip()
                    glycosylation_sites = parts[1].strip()
                    protein_id = proteoform_id.rsplit('_PF_', 1)[0] 
                    merged_file.write(f"{protein_id},{proteoform_id},{glycosylation_sites}\n")
                else:
                    logging.warning(f"Skipping malformed proteoform line: {line}")
        logging.info(f"Merged proteoforms saved to: {merged_proteoforms_path}")
    except Exception as e:
        logging.error(f"Error writing merged proteoforms file: {e}")
    
    # Create a detailed log file and a summary CSV.
    log_file_path = base_output_dir / f"02_input_log_{Path(input_file).name}.txt"
    summary_csv_path = base_output_dir / f"03_summary_{Path(input_file).name}.csv"
    
    try:
        with log_file_path.open("w") as log_file:
            log_file.write("Glycopeptide Proteoform Generator\n")
            log_file.write(f"Run Date and Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            log_file.write("\nInput Parameters:\n")
            log_file.write(f"  Input File: {input_file}\n")
            log_file.write(f"  Limit: {limit}\n")
            log_file.write(f"  Protein Column: {protein_col}\n")
            log_file.write(f"  Glycosylation Site Column: {site_col}\n")
            log_file.write(f"  Glycan Column: {glycan_col}\n")
            log_file.write("\n--- Processing Summary ---\n")
            
            total_proteins = len(protein_dict_for_processing)
            log_file.write(f"Total number of proteins processed: {total_proteins}\n")
            total_proteoforms_generated = sum(counts.values())
            log_file.write(f"Total number of proteoforms generated: {total_proteoforms_generated}\n")

            if proteins_hitting_limit:
                log_file.write(f"Proteins that hit the proteoform limit ({limit}): {len(proteins_hitting_limit)}/{total_proteins}\n")
                log_file.write(f"  IDs: {', '.join(sorted(proteins_hitting_limit))}\n")
            else:
                log_file.write("No proteins hit the specified proteoform limit.\n")

            end_time = time.time()
            runtime_seconds = end_time - start_time
            log_file.write(f"Total script runtime: {runtime_seconds:.2f} seconds\n")

        # Prepare data for summary CSV
        summary_data = []
        for protein, protein_data in sorted(protein_dict_for_processing.items()):
            total_sites = len(protein_data)
            # Sum total unique glycans observed for the protein
            total_glycans_for_protein = sum(len(site_glycans[0]) - 1 for site_glycans in protein_data.values()) # -1 for 'None'
            generated_proteoforms_for_protein = counts.get(protein, 0) 
            summary_data.append((protein, total_sites, total_glycans_for_protein, generated_proteoforms_for_protein))
        
        summary_df = pd.DataFrame(summary_data, 
                                  columns=["Protein", "TotalGlycosylationSites", "TotalGlycans", "TotalProteoformsGenerated"])
        summary_df.to_csv(summary_csv_path, index=False)
        
        logging.info(f"Log file saved to: {log_file_path}")
        logging.info(f"Summary CSV saved to: {summary_csv_path}")

    except Exception as e:
        logging.error(f"Error writing log/summary files: {e}")
    
    logging.info("All output files have been generated successfully.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Generate glycopeptide proteoforms from a CSV file.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input', default='human_proteoform_glycosylation_sites_gptwiki.csv',
                        help='Path to the input CSV file.')
    parser.add_argument('-l', '--limit', type=int, default=1000,
                        help='Maximum number of proteoforms to generate per protein.')
    parser.add_argument('-p', '--protein', default='protein', help='Name of the protein column.')
    parser.add_argument('-s', '--site', default='glycosylation_site', help='Name of the glycosylation site column.')
    parser.add_argument('-g', '--glycan', default='glycan', help='Name of the glycan column.')
    
    args = parser.parse_args()
    main(args.input, args.limit, args.protein, args.site, args.glycan)
