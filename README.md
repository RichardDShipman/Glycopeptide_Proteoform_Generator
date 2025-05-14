# Glycopeptide Proteoform Generator

Richard Shipman – 25OCT2024

GitHub: https://github.com/RichardDShipman

This repository hosts an optimized Python script, glycopeptide_proteoform_generator_cmd.py, designed for processing glycopeptide data to generate proteoforms by exploring combinations of glycopeptides across the glycosylation sites associated with a given protein. The script reads glycopeptide data from a CSV file, formats the information, and generates proteoforms while allowing users to set a configurable limit on the number of proteoform combinations per protein. This feature enables users to control output complexity based on their requirements and hardware setup. The generated results are stored in both text and CSV file formats within the designated data folder.

## Features

Read and Process Data:
Reads glycopeptide data from a CSV file and processes it to create a structured dictionary of glycopeptides grouped by protein.

Generate Proteoforms:
Generates combinations of glycosylation for peptides and creates all possible proteoforms, with a limit on the number of proteoforms generated per protein to manage large datasets.

Parallel Processing:
Leverages Python’s ProcessPoolExecutor to process each protein concurrently for faster performance on multi-core systems.

Robust File Handling and Logging:
Uses the pathlib module for safe file and path management, along with the built-in logging module for detailed runtime reporting and error handling.

Output Results:
Saves proteoform counts to a CSV file and detailed proteoform information to individual text files for each protein.

Merge and Summarize:
Merges individual proteoform text files into a single CSV file and generates log and summary files containing input parameters and proteoform counts.

##Requirements

	•	Python 3.7+

	•	Pandas

You can install the required Python libraries using pip:

```sh
pip install pandas
```

## Usage

1.	Prepare the Input File

The input CSV file should contain the following columns:

protein: The protein identifier.

glycosylation_site: The site of glycosylation.

glycan: The glycan. (There are no strict format requirements for glycan data.)

Example CSV format: (The human_proteoform_glycosylation_sites_gptwiki.csv file is used as an example.)

```csv
"protein","glycosylation_site","amino_acid","glycan","glycosylation_type","xref_key","xref_id","src_xref_key","src_xref_id","glycopeptide_id","composition","glycan_xref_key","glycan_xref_id","n_sequon","n_sequon_type","start_pos","end_pos","start_aa","end_aa","site_seq"
"O00754-1","692","Asn","G28681TP","N-linked","protein_xref_gptwiki","O00754@N692","protein_xref_gptwiki","O00754@N692","PE001986","HexNAc(2)Hex(3)","glycan_xref_gptwiki","G28681TP","NFS","NXS","692","692","Asn","Asn","N"
"O00754-1","930","Asn","G41247ZX","N-linked","protein_xref_gptwiki","O00754@N930","protein_xref_gptwiki","O00754@N930","PE002013","HexNAc(2)Hex(6)","glycan_xref_gptwiki","G41247ZX","NLS","NXS","930","930","Asn","Asn","N"
"O14672-1","278","Asn","G80920RR","N-linked","protein_xref_gptwiki","O14672@N278","protein_xref_gptwiki","O14672@N278","PE001964","HexNAc(2)Hex(9)","glycan_xref_gptwiki","G80920RR","NTT","NXT","278","278","Asn","Asn","N"
"O43405-1","100","Asn","G27126ED","N-linked","protein_xref_gptwiki","O43405@N100","protein_xref_gptwiki","O43405@N100","PE002197","HexNAc(4)Hex(6)NeuAc(1)","glycan_xref_gptwiki","G27126ED","NYS","NXS","100","100","Asn","Asn","N"
"O43852-1","131","Asn","G25987BV","N-linked","protein_xref_gptwiki","O43852@N131","protein_xref_gptwiki","O43852@N131","PE001250","HexNAc(4)Hex(3)dHex(1)","glycan_xref_gptwiki","G25987BV","NAT","NXT","131","131","Asn","Asn","N"
```

2.	Run the Script

Execute the script using Python. For example, to generate a maximum of 10 proteoforms per protein, run:

```sh
python glycopeptide_proteoform_generator_cmd.py -i human_proteoform_glycosylation_sites_gptwiki.csv -l 10
```

## Arguements

	•	-i: Input CSV file with glycopeptide data.
	•	-l: Limit the number of proteoforms generated for each protein.
	•	-p: Name of the protein column in the input CSV file.
	•	-s: Name of the glycosylation site column in the input CSV file.
	•	-g: Name of the glycan column in the input CSV file.
The script will process the data, generate proteoforms concurrently, and write the results to the data folder along with counts and merged files.

3.	Output Folder and Files

A folder named after the input CSV file (without extension) is created under the data directory. This folder will contain:

protein_proteoforms.txt:

For each protein, a text file with unique proteoform IDs and their corresponding glycosylation site–glycan combinations.

Example:

```csv
P00450-1_PF_1, 138-G10486CT:358-G10486CT:397-G06247RL:762-G27947YN
P00450-1_PF_2, 138-G10486CT:358-G10486CT:397-G06247RL:762-G10486CT
P00450-1_PF_3, 138-G10486CT:358-G10486CT:397-G06247RL:762-G47737VJ

```

#### 00_proteoform_counts_filename.csv:

A CSV file listing each protein and the total number of proteoforms generated.

```csv
protein,total_proteoforms
O00754-1,4
O14672-1,2
O43405-1,2
O43852-1,4

```

#### 01_merged_proteoforms_filename.csv:

A merged CSV file combining proteoform data from all proteins, with columns for protein, proteoform ID, and glycosylation site details.

```csv
protein,proteoform_id,glycosylation_sites
O00754-1,O00754-1_PF_1,692-G28681TP:930-None
O00754-1,O00754-1_PF_2,692-None:930-None
O00754-1,O00754-1_PF_3,692-None:930-G41247ZX
O00754-1,O00754-1_PF_4,692-G28681TP:930-G41247ZX

```


#### 02_input_log_filename.txt:

A log file detailing the input parameters and a summary of the processing run.

```txt
Glycopeptide Proteoform Generator
Run Date and Time: 2025-01-31 16:57:55

Input Parameters:
Input File: human_proteoform_glycosylation_sites_gptwiki.csv
Limit: 10
Protein Column: protein
Glycosylation Site Column: glycosylation_site
Glycan Column: glycan

Summary:
Total number of proteins: 140
Total number of proteoforms generated: 887

```

#### 03_summary_<filename>.csv:
A CSV file containing summary data for each protein, including the number of glycosylation sites, total glycans, and proteoform counts.

```csv
Protein,TotalGlycosylationSites,TotalGlycans,TotalProteoforms
O00754-1,2,2,4
O14672-1,1,1,2
O43405-1,1,1,2
O43852-1,1,1,4

```

#### Customization (Parameters)

Proteoform Limit:

Adjust the -l parameter to control the maximum number of proteoforms generated per protein. Because glycoproteins can generate millions of possible proteoforms, setting an appropriate limit helps manage computational and storage requirements.

Parallel Processing:

The script processes each protein concurrently. This behavior can be adjusted if necessary by modifying the parallel processing section in the source code.

Preparing human_proteoform_glycosylation_sites_gptwiki.csv Data (Glycopeptide Data Example CSV File)

The CSV should include at least the following columns: protein, glycosylation_site, and glycan.

For your convenience, the following columns were renamed:
- uniprotkb_canonical_ac → protein
- glycosylation_site_uniprotkb → glycosylation_site
- saccharide → glycan

This has been pre-processed for the file human_proteoform_glycosylation_sites_gptwiki.csv.

Note: Data sourced from the GlyGen data repository. For more details on the glycoproteomics data, please refer to the URL provided below.

## Test Data: Human Glycosylation Sites [GPTwiki]

Human Glycosylation Sites [GPTwiki] are provided by the Clinical and Translational Glycoscience Research Center (CTGRC) at Georgetown University. The database contains human (taxid:9606) protein glycosylation sites and associated glycan data from the GPTwiki database (GPTwiki Main Page). The proteins listed are part of the GlyGen UniProtKB canonical list (GlyGen Data).
Filename: human_proteoform_glycosylation_sites_gptwiki.csv

## Reference Source Materials

GlyGen: Computational and Informatics Resources for Glycoscience,
Glycobiology, Volume 30, Issue 2, February 2020, Pages 72–73,
https://doi.org/10.1093/glycob/cwz080

## License

This project is licensed under the GPL-3.0 license.
