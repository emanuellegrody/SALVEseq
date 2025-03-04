from multiprocessing import Pool
import pandas as pd
import numpy as np
from pathlib import Path
import glob
import gzip
from tqdm import tqdm
from itertools import combinations


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
        cellids = group['cellID'].values  # Store cellIDs
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
            'cellID': cellids[mask],  # Include cellIDs in output
            'UMI': umis[mask],
            'count': counts[mask]
        })

    print("Performing UMI error correction...")
    return df.groupby('cellID', group_keys=False).apply(process_group).reset_index(drop=True)

def find_closest_barcode(barcode: str, whitelist: set) -> str:
    """Find closest barcode in whitelist with Hamming distance 1."""
    if not isinstance(barcode, str):  # Handle non-string input
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


def process_barcodes(barcodes: list, whitelist: set, num_cores: int) -> dict:
    """Process barcodes in parallel."""
    with Pool(num_cores) as pool:
        corrected = list(tqdm(
            pool.starmap(find_closest_barcode,
                         ((b, whitelist) for b in barcodes if pd.notna(b))),
            total=len(barcodes)
        ))

    return dict(zip(barcodes, corrected))


def main():
    # Set paths
    reads_directory = "/projects/b1042/GoyalLab/egrody/20230424_EGS001/VISER/analysis/KLRB1/seqIO/"
    output_dir = "/projects/b1042/GoyalLab/egrody/20230424_EGS001/VISER/analysis/KLRB1/postSeqIO/"
    whitelist_path = "/home/egy2296/packages/cellranger-7.2.0/lib/python/cellranger/barcodes/3M-february-2018.txt.gz"

    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Read whitelist
    print("Reading whitelist...")
    with gzip.open(whitelist_path, 'rt') as f:
        whitelist = set(line.strip() for line in f)

    # Read and combine data files
    print("Reading and combining data files...")
    dfs = []
    for file_path in glob.glob(f"{reads_directory}*_shavedReads.txt"):
        print(f"Processing {Path(file_path).name}")
        df = pd.read_csv(file_path)
        dfs.append(df)

    # Combine and process initial data
    print("Processing initial data...")
    data = pd.concat(dfs, ignore_index=True)
    data = (data[['cellID', 'UMI']]
            .groupby(['cellID', 'UMI'])
            .size()
            .reset_index(name='count'))

    # Save full dataset
    data.to_csv(f"{output_dir}VISER_KLRB1_full.csv", index=False)

    # UMI error correction
    print("Correcting UMI errors...")
    data_umi_corrected = umi_error_correction(data)
    data_umi_corrected.to_csv(f"{output_dir}VISER_KLRB1_UMIcorrect.csv", index=False)

    # Before cell ID correction
    print("\nBefore cell ID correction:")
    print("Checking for AACAAAGAGGTTGGTG in input data:")
    print(data_umi_corrected[data_umi_corrected['cellID'] == "AACAAAGAGGTTGGTG"])

    # Getting unique barcodes
    print("\nGetting unique barcodes...")
    unique_barcodes = data_umi_corrected['cellID'].unique()
    print("AACAAAGAGGTTGGTG in unique_barcodes:", "AACAAAGAGGTTGGTG" in unique_barcodes)

    # Cell ID correction
    print("\nCorrecting cell IDs...")
    correction_dict = process_barcodes(unique_barcodes, whitelist, num_cores=8)

    # After correction
    print("\nAfter correction mapping:")
    if "AACAAAGAGGTTGGTG" in correction_dict:
        print("Mapped to:", correction_dict["AACAAAGAGGTTGGTG"])
    else:
        print("Cell ID not found in correction dictionary")

    # Before applying corrections - only check cellID since corrected_cellID doesn't exist yet
    print("\nBefore applying corrections:")
    relevant_rows = data_umi_corrected[data_umi_corrected['cellID'] == "AACAAAGAGGTTGGTG"]
    print(relevant_rows)

    # Apply corrections
    data_umi_corrected['corrected_cellID'] = data_umi_corrected['cellID'].map(correction_dict)

    # Now we can check both columns since corrected_cellID exists
    print("\nAfter applying corrections but before filtering:")
    relevant_rows = data_umi_corrected[
        (data_umi_corrected['cellID'] == "AACAAAGAGGTTGGTG") |
        (data_umi_corrected['corrected_cellID'] == "AACAAAGAGGTTGGTG")
        ]
    print(relevant_rows)

    # Create final clean dataset
    clean_data = (data_umi_corrected[
                      (data_umi_corrected['count'] > 1) &
                      (data_umi_corrected['corrected_cellID'].notna())
                      ]
                  .groupby('corrected_cellID')
                  .size()
                  .reset_index(name='count')
                  .rename(columns={'corrected_cellID': 'cellID'}))

    # Check final result
    print("\nIn final clean data:")
    print(clean_data[clean_data['cellID'] == "AACAAAGAGGTTGGTG"])

    # Additional check for high counts
    print("\nAll entries with count > 1000:")
    print(clean_data[clean_data['count'] > 1000])

    # Save final output
    print("Saving final output...")
    clean_data.to_csv(f"{output_dir}VISER_KLRB1_clean.csv", index=False)
    # Save correction map
    correction_df = pd.DataFrame(list(correction_dict.items()),
                                 columns=['original', 'corrected'])
    correction_df.to_csv(f"{output_dir}correction_map.csv", index=False)
    print("Processing complete!")


if __name__ == "__main__":
    main()
