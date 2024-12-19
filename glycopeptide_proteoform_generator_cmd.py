#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Oct 24 2024

@author: richarddshipman
"""

import itertools
import os
import pandas as pd
from collections import defaultdict
import argparse
from datetime import datetime

# Function to generate proteoforms for a given protein with a limit
def generate_proteoforms_with_limit(protein_name, protein_data, limit=100):
    # Generate combinations of glycosylation for each protein
    protein_forms = []
    for protein, glyco_options in protein_data.items():
        combinations = list(itertools.product(*glyco_options))
        protein_forms.append(combinations)
    
    # Generate all possible proteoforms by combining glycopeptides
    proteoforms = itertools.product(*protein_forms)
    
    # Limit the number of proteoforms, ensure uniqueness with set
    limited_proteoforms = set()
    count = 0
    for proteoform in proteoforms:
        if count >= limit:
            break
        limited_proteoforms.add(proteoform)  # Add to set to ensure uniqueness
        count += 1
    
    return list(limited_proteoforms)  # Convert back to list for further processing

# Main function to generate proteoforms from glycopeptides data
def main(input_file, limit, protein_col, site_col, glycan_col):
    """
    Main function to generate glycopeptide proteoforms from an input CSV file.
    Args:
        input_file (str): Path to the input CSV file containing protein, glycosylation_site, and glycan columns.
        limit (int): Limit on the number of proteoforms to generate per protein.
        protein_col (str): Name of the protein column in the input CSV file.
        site_col (str): Name of the glycosylation site column in the input CSV file.
        glycan_col (str): Name of the glycan column in the input CSV file.
    The function performs the following steps:
    1. Reads the input CSV file into a pandas DataFrame.
    2. Processes each row to create a dictionary of glycopeptides grouped by protein and glycosylation site.
    3. Formats the glycopeptides into a list of dictionaries.
    4. Groups glycopeptides by protein and generates combinations of glycosylation sites and glycans.
    5. Creates an output directory named after the input file (excluding extension).
    6. Writes the proteoform counts to a CSV file and individual proteoform details to text files.
    7. Removes duplicate glycosylation sites from the proteoform text files.
    8. Merges all proteoform text files into a single CSV file.
    Outputs:
        - A directory containing:
            - A CSV file with proteoform counts for each protein.
            - Text files with detailed proteoform information for each protein.
            - A merged CSV file with all proteoform details.
            - A log file with details of the input parameters and a report on counts.
    """
    # Read the CSV file
    df = pd.read_csv(input_file)

    # Initialize a dictionary to hold the formatted data
    glycopeptides = defaultdict(lambda: defaultdict(set))

    # Process each row in the DataFrame
    for _, row in df.iterrows():
        protein = row[protein_col]
        glycosylation_site = row[site_col]
        glycan = row[glycan_col]
        
        glycopeptides[protein][glycosylation_site].add(glycan)

    # Convert the defaultdict to the desired format
    formatted_glycopeptides = [
        {
            'protein': protein,
            'glycosylation_site': site,
            'glycans': list(glycans)
        }
        for protein, sites in glycopeptides.items()
        for site, glycans in sites.items()
    ]

    # Group glycopeptides by protein
    protein_dict = defaultdict(lambda: defaultdict(list))

    for glycopeptide in formatted_glycopeptides:
        protein = glycopeptide["protein"]
        glycosylation_site = glycopeptide["glycosylation_site"]
        glycans = glycopeptide["glycans"]
        
        # Add an option for no glycosylation (None) along with the possible glycans
        glyco_combinations = [(glycosylation_site, None)] + [(glycosylation_site, glycan) for glycan in glycans]
        protein_dict[protein][glycosylation_site].append(list(set(glyco_combinations)))

    # Create a directory for the proteoform output named after the input file (excluding extension)
    base_output_dir = os.path.join("data", os.path.splitext(os.path.basename(input_file))[0])
    os.makedirs(base_output_dir, exist_ok=True)

    # Open the CSV file for writing the proteoform counts
    with open(os.path.join(base_output_dir, f'00_proteoform_counts_{input_file}'), 'w') as counts_file:
        counts_file.write('protein,total_proteoforms\n')  # Write the header

        # Process each protein and write results to a text file and CSV
        for protein, protein_data in protein_dict.items():
            proteoforms = generate_proteoforms_with_limit(protein, protein_data, limit)
            total_proteoforms = len(proteoforms)
            
            # Write results to a text file in the output directory
            with open(os.path.join(base_output_dir, f"{protein}_proteoforms.txt"), "w") as file:
                for idx, proteoform in enumerate(proteoforms, 1):
                    file.write(f"{protein}_PF_{idx}, ")
                    for protein_combo in proteoform:
                        for (glycosylation_site, glycan) in protein_combo:
                            file.write(f"{glycosylation_site}-{glycan} ")
                    file.write("\n")
            
            # Write the count for this protein to the CSV file
            counts_file.write(f'{protein},{total_proteoforms}\n')

            # Print summary to console
            print(f"{protein}: Total number of proteoforms: {total_proteoforms}")

    print("Proteoform counts have been written to the output directory.")

    # Check _proteoforms.txt files in the output directory for duplicate glycosylation sites for a proteoform_id, remove duplicates
    for protein in protein_dict.keys():
        proteoform_file_path = os.path.join(base_output_dir, f"{protein}_proteoforms.txt")
        if os.path.exists(proteoform_file_path):
            with open(proteoform_file_path, 'r') as file:
                lines = file.readlines()
            
            with open(proteoform_file_path, 'w') as file:
                for line in lines:
                    if line.strip():  # Skip empty lines
                        parts = line.split(', ')[1:]  # Skip the proteoform identifier
                        glycosylation_pairs = set()
                        unique_parts = []
                        for part in parts[0].split(' '):  # Split the second column by space
                            glyco_pair = tuple(part.split('-'))  # Create a tuple of (site, glycan)
                            if glyco_pair not in glycosylation_pairs:
                                glycosylation_pairs.add(glyco_pair)
                                unique_parts.append(part)
                        if unique_parts:
                            file.write(f"{line.split(', ')[0]}, {' '.join(unique_parts)}")
                        else:
                            print(f"Duplicate glycosylation site found and removed in proteoform: {line.strip()}")

    # Merge all _proteoforms.txt files in the output directory into a single CSV file
    merged_proteoforms_path = os.path.join(base_output_dir, f"01_merged_proteoforms_{input_file}")
    with open(merged_proteoforms_path, 'w') as merged_file:
        merged_file.write('protein,proteoform_id,glycosylation_sites\n')  # Write the header

        for protein in protein_dict.keys():
            proteoform_file_path = os.path.join(base_output_dir, f"{protein}_proteoforms.txt")
            if os.path.exists(proteoform_file_path):
                with open(proteoform_file_path, 'r') as file:
                    lines = file.readlines()
                    for line in lines:
                        if line.strip():  # Skip empty lines
                            proteoform_id, glycosylation_sites = line.split(', ', 1)
                            merged_file.write(f"{protein},{proteoform_id},{glycosylation_sites.strip()}\n")

    # Create a log file with details of the input parameters and a report on counts

    log_file_path = os.path.join(base_output_dir, f"02_input_log_{os.path.basename(input_file)}.txt")
    with open(log_file_path, 'w') as log_file:
        # Write the title and run date/time
        log_file.write("Glycopeptide Proteoform Generator\n")
        log_file.write(f"Run Date and Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        log_file.write("\nInput Parameters:\n")
        log_file.write(f"Input File: {input_file}\n")
        log_file.write(f"Limit: {limit}\n")
        log_file.write(f"Protein Column: {protein_col}\n")
        log_file.write(f"Glycosylation Site Column: {site_col}\n")
        log_file.write(f"Glycan Column: {glycan_col}\n")
        
        # Write the summary section
        log_file.write("\nSummary:\n")
        total_proteins = len(protein_dict)
        log_file.write(f"Total number of proteins: {total_proteins}\n")

        total_proteoforms = sum(len(generate_proteoforms_with_limit(protein, protein_data, limit)) for protein, protein_data in protein_dict.items())
        log_file.write(f"Total number of proteoforms generated: {total_proteoforms}\n")

        log_file.write("\nProtein,TotalGlycosylationSites,TotalGlycans,TotalProteoforms\n")
        summary_data = []
        for protein, protein_data in protein_dict.items():
            total_sites = len(protein_data)
            total_glycans = sum(len(glyco_options) for glyco_options in protein_data.values())
            proteoforms = generate_proteoforms_with_limit(protein, protein_data, limit)
            total_proteoforms = len(proteoforms)
            summary_data.append((protein, total_sites, total_glycans, total_proteoforms))
            log_file.write(f"{protein},{total_sites},{total_glycans},{total_proteoforms}\n")

        # Optionally, save the summary data to a CSV file for easier data extraction
        summary_df = pd.DataFrame(summary_data, columns=["Protein", "TotalGlycosylationSites", "TotalGlycans", "TotalProteoforms"])
        summary_csv_path = os.path.join(base_output_dir, f"03_summary_{os.path.basename(input_file)}.csv")
        summary_df.to_csv(summary_csv_path, index=False)

    print("Log file has been written to the output directory.")

# This line checks if this script is being run as the main program
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate proteoforms from glycopeptides data.')
    parser.add_argument('-i', '--input', default='human_proteoform_glycosylation_sites_gptwiki.csv', help='Path to the input CSV file with glycopeptides data.')
    parser.add_argument('-l', '--limit', type=int, default=10, help='Maximum number of proteoforms to generate per protein.')
    parser.add_argument('-p', '--protein', default='protein', help='Name of the protein column in the input CSV file.')
    parser.add_argument('-s', '--site', default='glycosylation_site', help='Name of the glycosylation site column in the input CSV file.')
    parser.add_argument('-g', '--glycan', default='glycan', help='Name of the glycan column in the input CSV file.')

    args = parser.parse_args()
    main(args.input, args.limit, args.protein, args.site, args.glycan)