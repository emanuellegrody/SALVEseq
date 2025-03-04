from multiprocessing import Pool
import pandas as pd
import numpy as np
from pathlib import Path
import glob
import gzip
from tqdm import tqdm
from itertools import combinations
import argparse


def hamming_distance(s1: str, s2: str) -> int:
    """Calculate Hamming distance between two strings."""
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def umi_error_correction(df: pd.DataFrame) -> pd.DataFrame:
    """Correct UMI errors within each cell group."""

    def process_group(group):
        if len(group) <= 1:
            return group

        umis = group['UMI'].values
        counts = group['count'].values
        cellids = group['cellID'].values
        targets = group['target'].values
        mask = np.ones(len(umis), dtype=bool)

        for i, j in combinations(range(len(umis)), 2):
            if not (mask[i] and mask[j]):
                continue
            if hamming_distance(umis[i], umis[j]) == 1:
                if counts[i] >= counts[j]:
                    counts[i] += counts[j]
                    mask[j] = False
                else:
                    counts[j] += counts[i]
                    mask[i] = False

        return pd.DataFrame({
            'cellID': cellids[mask],
            'UMI': umis[mask],
            'target': targets[mask],
            'count': counts[mask]
        })

    print("Performing UMI error correction...")
    return df.groupby('cellID', group_keys=False).apply(process_group).reset_index(drop=True)


def find_closest_barcode(barcode: str, whitelist: set) -> str:
    """Find closest barcode in whitelist with Hamming distance 1."""
    if not isinstance(barcode, str):
        return None

    if barcode in whitelist:
        return barcode

    for i in range(len(barcode)):
        for base in 'ACGT':
            if base != barcode[i]:
                candidate = barcode[:i] + base + barcode[i + 1:]
                if candidate in whitelist:
                    return candidate
    return None


def process_barcodes(barcodes: list, whitelist: set) -> dict:
    """Process barcodes in parallel."""
    with Pool(8) as pool:  # Hardcoded to 8 cores
        corrected = list(tqdm(
            pool.starmap(find_closest_barcode,
                         ((b, whitelist) for b in barcodes if pd.notna(b))),
            total=len(barcodes)
        ))

    return dict(zip(barcodes, corrected))


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Process sequencing data with UMI and cell ID correction.')
    parser.add_argument('--reads_directory', required=True,
                        help='Directory containing the input read files')
    parser.add_argument('--output_dir', required=True,
                        help='Directory for output files')

    return parser.parse_args()


def main():
    # Parse command line arguments
    args = parse_arguments()

    # Hardcoded whitelist path
    whitelist_path = "/home/egy2296/packages/cellranger-7.2.0/lib/python/cellranger/barcodes/3M-february-2018.txt.gz"

    # Create output directory if it doesn't exist
    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    # Read whitelist
    print("Reading whitelist...")
    with gzip.open(whitelist_path, 'rt') as f:
        whitelist = set(line.strip() for line in f)

    # Read and combine data files
    print("Reading and combining data files...")
    file_pattern = f"{args.reads_directory}*_shavedReads.txt"
    print(f"Looking for files matching pattern: {file_pattern}")

    matching_files = glob.glob(file_pattern)
    if not matching_files:
        raise ValueError(
            f"No files found matching pattern '{file_pattern}'. Please check the directory path and file naming.")

    dfs = []
    for file_path in matching_files:
        try:
            df = pd.read_csv(file_path)
            dfs.append(df)
        except Exception as e:
            print(f"Error reading file {file_path}:")
            print(e)

    # Combine and process initial data
    print("Processing initial data...")
    data = pd.concat(dfs, ignore_index=True)
    data = (data[['cellID', 'UMI', 'target']]
            .groupby(['cellID', 'UMI', 'target'])
            .size()
            .reset_index(name='count'))

    # Save full dataset
    data.to_csv(f"{args.output_dir}VISER_KLRB1_full.csv", index=False)

    # UMI error correction
    print("Correcting UMI errors...")
    data_umi_corrected = umi_error_correction(data)
    data_umi_corrected.to_csv(f"{args.output_dir}VISER_KLRB1_UMIcorrect.csv", index=False)

    # Cell ID correction
    print("Correcting cell IDs...")
    unique_barcodes = data_umi_corrected['cellID'].unique()
    correction_dict = process_barcodes(unique_barcodes, whitelist)

    # Apply corrections and create correction map
    print("Applying corrections and creating correction map...")
    correction_map = pd.DataFrame({
        'original_cellID': list(correction_dict.keys()),
        'corrected_cellID': list(correction_dict.values())
    })
    correction_map.to_csv(f"{args.output_dir}correction_map.csv", index=False)

    # Apply corrections to the data
    print("Applying corrections and filtering...")
    data_umi_corrected['corrected_cellID'] = data_umi_corrected['cellID'].map(correction_dict)

    # Create detailed output with corrected cellID, UMI, and target
    detailed_data = (data_umi_corrected[
                         (data_umi_corrected['count'] > 1) &
                         (data_umi_corrected['corrected_cellID'].notna()) &
                         (data_umi_corrected['corrected_cellID'].isin(whitelist))
                         ][['corrected_cellID', 'UMI', 'target', 'count']]
                     .rename(columns={'corrected_cellID': 'cellID'}))

    detailed_data.to_csv(f"{args.output_dir}VISER_KLRB1_clean_target.csv", index=False)

    # Create final clean dataset (aggregated counts per cellID)
    clean_data = (detailed_data
                  .groupby('cellID')['count']
                  .sum()
                  .reset_index())

    # Print statistics
    print(f"\nFiltering statistics:")
    print(f"Total unique cellIDs before filtering: {len(unique_barcodes)}")
    print(f"Total unique cellIDs after correction and filtering: {len(clean_data)}")
    print(f"Percentage retained: {(len(clean_data) / len(unique_barcodes) * 100):.2f}%")

    # Save final output
    print("\nSaving final outputs...")
    clean_data.to_csv(f"{args.output_dir}VISER_KLRB1_clean.csv", index=False)

    print("Processing complete!")


if __name__ == "__main__":
    main()