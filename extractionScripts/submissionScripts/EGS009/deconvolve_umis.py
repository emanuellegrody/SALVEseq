import tables
import pandas as pd
import argparse
import numpy as np


def get_umi_sequence(umi_int):
    """Convert numeric UMI to sequence"""
    bases = ['A', 'C', 'G', 'T']
    umi_seq = ''
    # 10x uses 2 bits per base, so we need to decode 2 bits at a time
    for i in range(12):  # UMIs are 12bp long
        base_idx = (umi_int >> (2 * i)) & 3  # Get 2 bits at a time
        umi_seq = bases[base_idx] + umi_seq
    return umi_seq


def extract_10x_molecular_info(molecule_info_path, barcode_file_path, barcode_column=None):
    """Extract UMI information from Cell Ranger molecule_info.h5 file"""

    # Read cell barcodes from CSV and remove "-1" suffix
    barcodes_df = pd.read_csv(barcode_file_path)
    column = barcode_column if barcode_column else barcodes_df.columns[0]
    barcodes_of_interest = set(bc.split('-')[0] for bc in barcodes_df[column])

    # Read molecule info
    with tables.open_file(molecule_info_path, 'r') as f:
        # Read barcode indices and actual barcodes
        barcode_idx = f.root.barcode_idx.read()
        barcodes = f.root.barcodes.read()

        # Convert barcodes to strings
        h5_barcodes = [bc.decode() for bc in barcodes[barcode_idx]]

        # Create DataFrame with molecular info
        df = pd.DataFrame({
            'barcode': h5_barcodes,
            'umi_int': f.root.umi.read(),  # Keep numeric UMI temporarily
            'feature_idx': f.root.feature_idx.read()
        })

        # Convert gene indices to names using features group
        features_ref = f.root.features.id.read()
        df['gene'] = [features_ref[idx].decode() for idx in df['feature_idx']]

    # Filter for barcodes of interest and drop feature_idx column
    df = df[df['barcode'].isin(barcodes_of_interest)].drop('feature_idx', axis=1)

    # Add the "-1" suffix back to the barcodes
    df['barcode'] = df['barcode'] + '-1'

    # Convert numeric UMIs to sequences
    df['umi'] = df['umi_int'].apply(get_umi_sequence)
    df = df.drop('umi_int', axis=1)

    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--molecule_info', required=True, help='Path to molecule_info.h5 file')
    parser.add_argument('-b', '--barcodes', required=True, help='Path to CSV file containing cell barcodes')
    parser.add_argument('-c', '--column', help='Name of column in CSV containing barcodes')
    parser.add_argument('-o', '--output', help='Output CSV path')
    args = parser.parse_args()

    # Process the data
    df = extract_10x_molecular_info(args.molecule_info, args.barcodes, args.column)

    # Save to CSV
    output_path = args.output if args.output else 'molecular_data.csv'
    df.to_csv(output_path, index=False)
    print(f"\nFinal results:")
    print(f"Processed {len(df)} molecules")
    print(f"Found {df['barcode'].nunique()} unique cells")
    print(f"Found {df['gene'].nunique()} unique genes")

if __name__ == "__main__":
    main()