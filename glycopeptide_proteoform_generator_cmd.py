import itertools
import os
import pandas as pd
from collections import defaultdict
import argparse

# Function to generate proteoforms for a given protein with a limit
def generate_proteoforms_with_limit(protein_name, protein_data, limit=100):
    # Generate combinations of glycosylation for each peptide
    peptide_forms = []
    for peptide, glyco_options in protein_data.items():
        combinations = list(itertools.product(*glyco_options))
        peptide_forms.append(combinations)
    
    # Generate all possible proteoforms by combining glycopeptides
    proteoforms = itertools.product(*peptide_forms)
    
    # Limit the number of proteoforms
    limited_proteoforms = []
    count = 0
    for proteoform in proteoforms:
        if count >= limit:
            break
        limited_proteoforms.append(proteoform)
        count += 1
    
    return limited_proteoforms

# Main function to generate proteoforms from glycopeptides data
def main(input_file, limit):
    # Read the CSV file
    # csv should have protein, peptide, glycosylation_site, and glycan columns
    # Notion of peptide and glycan defines level of detail for a proteoform
    df = pd.read_csv(input_file)

    # Initialize a dictionary to hold the formatted data
    glycopeptides = defaultdict(lambda: defaultdict(lambda: {'glycosylation_site': None, 'glycans': set()}))

    # Process each row in the DataFrame
    for _, row in df.iterrows():
        protein = row['protein']
        peptide = row['peptide']
        glycosylation_site = row['glycosylation_site']
        glycan = row['glycan']
        
        glycopeptides[protein][peptide]['glycosylation_site'] = glycosylation_site
        glycopeptides[protein][peptide]['glycans'].add(glycan)

    # Convert the defaultdict to the desired format
    formatted_glycopeptides = [
        {
            'protein': protein,
            'peptide': peptide,
            'glycosylation_site': details['glycosylation_site'],
            'glycans': list(details['glycans'])
        }
        for protein, peptides in glycopeptides.items()
        for peptide, details in peptides.items()
    ]

    # Group glycopeptides by protein
    protein_dict = defaultdict(lambda: defaultdict(list))

    for glycopeptide in formatted_glycopeptides:
        protein = glycopeptide["protein"]
        peptide = glycopeptide["peptide"]
        glycosylation_site = glycopeptide["glycosylation_site"]
        glycans = glycopeptide["glycans"]
        
        # Add an option for no glycosylation (None) along with the possible glycans
        protein_dict[protein][peptide].append([(glycosylation_site, None)] + [(glycosylation_site, glycan) for glycan in glycans])
    
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
                    for peptide_combo in proteoform:
                        for (glycosylation_site, glycan) in peptide_combo:
                            file.write(f"{glycosylation_site}-{glycan} ")
                    file.write("\n")
            # Write the count for this protein to the CSV file
            with open(os.path.join(base_output_dir, f'00_proteoform_counts_{input_file}'), 'a') as counts_file:
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
                        glycosylation_sites = set()
                        unique_parts = []
                        for part in parts[0].split(' '):  # Split the second column by space
                            glycosylation_site = part.split('-')[0]  # Check only the number left of '-'
                            if glycosylation_site not in glycosylation_sites:
                                glycosylation_sites.add(glycosylation_site)
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
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate proteoforms from glycopeptides data.')
    parser.add_argument('-i', '--input', required=True, help='Path to the input CSV file with glycopeptides data.')
    parser.add_argument('-l', '--limit', type=int, default=100, help='Maximum number of proteoforms to generate per protein.')

    args = parser.parse_args()
    main(args.input, args.limit)