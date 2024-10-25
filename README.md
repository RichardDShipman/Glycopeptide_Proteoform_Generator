# Glycopeptide Proteoform Generator

Richard Shipman -- 25OCT2024 

GitHub: https://github.com/RichardDShipman

This repository hosts a Python script `glycopeptide_proteoform_generator_cmd.py` designed for processing glycopeptide data to generate proteoforms by exploring combinations of glycopeptides across the glycosylation sites associated with a given protein. The script reads glycopeptide data from a CSV file, formats the information, and generates proteoforms while allowing users to set a configurable limit on the number of proteoforms combinations per protein. This feature enables users to control the output’s complexity based on their requirements and hardware setup. The generated results are stored in both text and CSV file formats within the designated data folder.

## Features

- **Read and Process Data**: Reads glycopeptides data from a CSV file and processes it to create a structured dictionary of glycopeptides grouped by protein.
- **Generate Proteoforms**: Generates combinations of glycosylation for peptides and creates all possible proteoforms, with a limit on the number of proteoforms generated to manage large datasets.
- **Output Results**: Saves proteoform counts to a CSV file and detailed proteoform information to individual text files for each protein.

## Requirements

- Python 3.7+
- Pandas
- itertools

You can install the required Python libraries using pip:

```bash
pip install pandas
```

## Usage

1. **Prepare the Input File**

   The input CSV file should contain the following columns:

   - `protein`: The protein identifier.
   - `glycosylation_site`: The site of glycosylation.
   - `glycan`: The glycan. (Note: There are no requirements for glycan data format.)

   Example CSV format: (`human_proteoform_glycosylation_sites_gptwiki.csv` file used as example, data source below.)


```csv
"protein","glycosylation_site","amino_acid","glycan","glycosylation_type","xref_key","xref_id","src_xref_key","src_xref_id","glycopeptide_id","composition","glycan_xref_key","glycan_xref_id","n_sequon","n_sequon_type","start_pos","end_pos","start_aa","end_aa","site_seq"
"O00754-1","692","Asn","G28681TP","N-linked","protein_xref_gptwiki","O00754@N692","protein_xref_gptwiki","O00754@N692","PE001986","HexNAc(2)Hex(3)","glycan_xref_gptwiki","G28681TP","NFS","NXS","692","692","Asn","Asn","N"
"O00754-1","930","Asn","G41247ZX","N-linked","protein_xref_gptwiki","O00754@N930","protein_xref_gptwiki","O00754@N930","PE002013","HexNAc(2)Hex(6)","glycan_xref_gptwiki","G41247ZX","NLS","NXS","930","930","Asn","Asn","N"
"O14672-1","278","Asn","G80920RR","N-linked","protein_xref_gptwiki","O14672@N278","protein_xref_gptwiki","O14672@N278","PE001964","HexNAc(2)Hex(9)","glycan_xref_gptwiki","G80920RR","NTT","NXT","278","278","Asn","Asn","N"
```


2. **Run the Script**

Execute the script using Python with a limit of 10 on the number of proteoforms generated per protein: 

   ```bash
   python glycopeptide_proteoform_generator_cmd.py -i human_proteoform_glycosylation_sites_gptwiki.csv -l 10
   ```

   - **-i**: input CSV file with glycopeptides.
   - **-l**: limit the number of proteoforms generated for each protein.
   - **<filename>**: glycopeptide data. `human_proteoform_glycosylation_sites_gptwiki.csv` file used as example glycopeptide data.

   The script will process the data, generate proteoforms, and write the results to the `data` folder and `00_proteoform_counts_<filename>.csv` & `01_merged_proteoforms_<filename>.csv` files.

3. **Output Folder and Files**

A folder titled the name of the inputted CSV file containing the the following results:

   - **data**: Contains text files named `<protein>_proteoforms.txt` for each protein, detailing the proteoforms generated and their components.
   - **<protein>_proteoforms.txt**: Unique proteoform ID for each proteoform combination in each `<protein>_proteoforms.txt` text file. Example: 

```CSV
O00754-1_PF_1, 692-G28681TP 930-G41247ZX 
O00754-1_PF_2, 692-None 930-None 
O00754-1_PF_3, 692-None 930-G41247ZX 
O00754-1_PF_4, 692-G28681TP 930-None 
```

Explaination of proteoform format with details used in this script below.

```txt
<protein_id>_PF_<index>, <glycosylation_site>-<glycan> <glycosylation_site>-<glycan> ...
```

   - **<protein_id>**: This represents the unique identifier of the protein (e.g., O00754-1). This ID is sourced from the input data and corresponds to the protein for which the proteoforms are generated.
   - **PF**: This indicates that the entry corresponds to a “Proteoform.”
   - **<index>**: This is a sequential number that identifies the specific proteoform generated for the protein (e.g., 1, 2, 3, etc.).
   - **<glycosylation_site>**: This refers to the specific site on the protein where glycosylation can occur. The site is represented by a number or label (e.g., 692, 930).
   - **<glycan>**: This indicates the type of glycan attached at the corresponding glycosylation site. It can either be a specific identifier for a glycan (e.g., G28681TP, G41247ZX) or None, indicating that no glycan is present at that site.

The following summary files are generated and placed alongside `<protein>_proteoforms.txt` containing total counts and merged proteoform results.

   - **00_proteoform_counts_<filename>.csv**: A CSV file with columns `protein` and `total_proteoforms`, listing the number of proteoforms generated for each protein.

```CSV
protein,total_proteoforms
O00754-1,4
O14672-1,2
O43405-1,2
O43852-1,4
```

   - **01_merged_proteoforms_<filename>.csv**: A CSV file with the columns `protein`, `proteoform_id`, and `glycosylation_sites` containing merged proteoform data.

```CSV
protein,proteoform_id,glycosylation_sites
O00754-1,O00754-1_PF_1,692-G28681TP 930-G41247ZX
O00754-1,O00754-1_PF_2,692-None 930-None
O00754-1,O00754-1_PF_3,692-None 930-G41247ZX
O00754-1,O00754-1_PF_4,692-G28681TP 930-None
O14672-1,O14672-1_PF_1,278-None
O14672-1,O14672-1_PF_2,278-G80920RR
O43405-1,O43405-1_PF_1,100-G27126ED
O43405-1,O43405-1_PF_2,100-None
O43852-1,O43852-1_PF_1,131-G80475RE
```

# Customization (Parameters)

- **Proteoform Limit**: You can adjust the limit parameter `-l` in the generate_proteoforms_with_limit function to control the maximum number of proteoforms generated per protein. Setting an appropriate limit is important, as glycoproteins can generate millions of possible proteoforms, which may quickly exhaust computational resources and storage.

- **NOTE**: If you set the limit above 10 million proteoforms, be prepared for high memory usage and significant processing time. Output text files for large limits can exceed 1GB in size, especially if glycopeptide data involves complex glycosylation patterns across many glycosylation sites.

## Preparing human_proteoform_glycosylation_sites_gptwiki.csv Data (Glycopeptide Data Example CSV File)

**Prepare the Input File**

The CSV should have the following columns: `protein`, `glycosylation_site`, and `glycan`.

Renamed the following columns in the CSV: `uniprotkb_canonical_ac` to `protein`, `glycosylation_site_uniprotkb` to `glycosylation_site`, and `saccharide` to `glycan`.

This step has been done already for the file `human_proteoform_glycosylation_sites_gptwiki.csv`.

**Note**: Data source from GlyGen data repository. Details below on the provided as example glycoproteomics data. Follow URL link for more info.

### Human Glycosylation Sites [GPTwiki]

Human Glycosylation Sites [GPTwiki], provided by the Clinical and Translational Glycoscience Research Center (CTGRC), Georgetown University. The database contains list of human [taxid:9606] proteins with information on glycosylation sites and associated glycans from GPTwiki database [https://edwardslab.bmcb.georgetown.edu/gptwiki/Main_Page]. The listed protein (UniProtKB) accessions are part of the GlyGen UniProtKB canonical list (https://data.glygen.org/GLYDS000001). Filename: `human_proteoform_glycosylation_sites_gptwiki.csv`

# Reference Source Materials 

GlyGen: Computational and Informatics Resources for Glycoscience, Glycobiology, Volume 30, Issue 2, February 2020, Pages 72–73, https://doi.org/10.1093/glycob/cwz080

## License

This project is licensed under the GPL-3.0 license.
