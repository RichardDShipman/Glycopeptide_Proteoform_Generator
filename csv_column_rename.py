import pandas as pd
import argparse

def main():
    """
    Main function to rename columns in a CSV file.

    This function uses argparse to handle command-line arguments for input file path and column names.
    It reads a CSV file, renames specified columns, and saves the updated DataFrame back to the same CSV file.

    Command-line arguments:
    -i, --input: Input CSV file path (default: 'human_proteoform_glycosylation_sites_gptwiki.csv')
    -p, --protein: Column name for protein (default: 'uniprotkb_canonical_ac')
    -s, --site: Column name for glycosylation site (default: 'glycosylation_site_uniprotkb')
    -g, --glycan: Column name for glycan (default: 'saccharide')

    The following columns are renamed:
    - 'uniprotkb_canonical_ac' to 'protein'
    - 'glycosylation_site_uniprotkb' to 'glycosylation_site'
    - 'saccharide' to 'glycan'
    """
    parser = argparse.ArgumentParser(description='Rename columns in a CSV file.')
    parser.add_argument('-i', '--input', type=str, default='human_proteoform_glycosylation_sites_gptwiki.csv', help='Input CSV file path')
    parser.add_argument('-p', '--protein', type=str, default='uniprotkb_canonical_ac', help='Column name for protein')
    parser.add_argument('-s', '--site', type=str, default='glycosylation_site_uniprotkb', help='Column name for glycosylation site')
    parser.add_argument('-g', '--glycan', type=str, default='saccharide', help='Column name for glycan')
    args = parser.parse_args()

    # Load the CSV file
    df = pd.read_csv(args.input)

    # Check if the columns have already been renamed
    if all(col in df.columns for col in ['protein', 'glycosylation_site', 'glycan']):
        print("Columns have already been renamed.")
    else:
        # Rename the columns
        df.rename(columns={
            args.protein: 'protein',
            args.site: 'glycosylation_site',
            args.glycan: 'glycan'
        }, inplace=True)

        # Save the updated DataFrame back to the same CSV file
        df.to_csv(args.input, index=False)
        print("Columns have been renamed and the updated CSV has been saved.")

if __name__ == '__main__':
    main()
