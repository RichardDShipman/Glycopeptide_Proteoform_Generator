# Glycopeptide Proteoform Generator

This repository contains a Python script for processing glycopeptides data to generate text files with proteoforms from combinations of glycopeptides across the glycosylation sites found for a protein. The script reads glycopeptides data from a CSV file, groups and formats the data, and generates proteoforms with a configurable limit on the number of combinations. Results are saved to text files and a CSV file.

## Features

- **Read and Process Data**: Reads glycopeptides data from a CSV file and processes it to create a structured dictionary of glycopeptides grouped by protein.
- **Generate Proteoforms**: Generates combinations of glycosylation for peptides and creates all possible proteoforms, with a limit on the number of proteoforms generated to manage large datasets.
- **Output Results**: Saves proteoform counts to a CSV file and detailed proteoform information to individual text files for each protein.

## Requirements

- Python 3.x
- Pandas
- itertools
- OS

You can install the required Python libraries using pip:

```bash
pip install pandas
```

## Usage

1. **Prepare the Input File**

   Ensure you have a CSV file named `glycopeptides.csv` in the same directory as the script. The CSV file should contain the following columns:

   - `protein`: The protein identifier.
   - `peptide`: The peptide sequence.
   - `glycosylation_site`: The site of glycosylation.
   - `glycan`: The glycan composition.

   Example CSV format:

   ```csv
   protein,peptide,glycosylation_site,glycan
   ACAN,TVYVHAnQTGYPDPSSR,343,N5H5F1S0G0
   ACAN,SnDSGVYR,132,N2H6F0S0G0
   ```

2. **Run the Script**

   Execute the script using Python:

   ```bash
   python glycopeptide_proteoform_generator_cmd.py -i glycopeptides.csv -l 10
   ```

   The script will process the data, generate proteoforms, and write the results to the `data` folder and `00_proteoform_counts.csv` file.

3. **Output Files**

   - **data**: Contains text files named `<protein>_proteoforms.txt` for each protein, detailing the proteoforms generated and their components.
   - **00_proteoform_counts.csv**: Contains a CSV file with columns `protein` and `total_proteoforms`, listing the number of proteoforms generated for each protein.
   - **01_merged_proteoforms_glycopeptides.csv**: Contains a CSV file with merged proteoform data.

## Code Overview

### Import Libraries

```python
import itertools
import os
import pandas as pd
from collections import defaultdict
```

### Read and Process Data

Read the CSV file and process it to group glycopeptides by protein and peptide.

### Generate Proteoforms

The `generate_proteoforms_with_limit` function generates proteoforms with a limit on the total number of combinations.

- **-i**: input glycopeptides file.
- **-l**: limit the number of proteoforms generated.

### Output Results

Write the proteoform counts to `00_proteoform_counts.csv` and detailed proteoforms to individual text files.

## Customization

- **Proteoform Limit**: You can adjust the `limit` parameter in the `generate_proteoforms_with_limit` function to control the maximum number of proteoforms generated per protein.

## Preparing Glycopeptides.csv Data

1. **Prepare the Input File**

   Ensure `glycopeptides.csv` is in the same directory as the script. The CSV should have the following columns: `protein`, `peptide`, `glycosylation_site`, and `glycan`.

   Note: Data source from Glygen data repository. Details below on glycoproteomics data. Extracted list of glycopeptides from csv file.
   
   GLY_001046
   
   Human ccRCC Glycoproteomics (ML Ready)

   The Human Glycosylation Sites (PDC) dataset contains intact glycopeptide abundances, biospecimen and ... 
   
   Filename: human_proteoform_ml_ready_pdc_ccrcc.csv

# Reference Source Materials 

GlyGen: Computational and Informatics Resources for Glycoscience, Glycobiology, Volume 30, Issue 2, February 2020, Pages 72–73, https://doi.org/10.1093/glycob/cwz080

## License

This project is licensed under the GNU GENERAL PUBLIC LICENSE.



