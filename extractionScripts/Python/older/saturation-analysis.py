import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from kneed import KneeLocator
import argparse
import os


def parse_arguments():
    parser = argparse.ArgumentParser(description='Perform saturation analysis on scRNA-seq data.')
    parser.add_argument('csv_path', type=str, help='Path to the input CSV file')
    return parser.parse_args()


def subsample_and_analyze(data, depth):
    subsampled = data.sample(n=depth, weights='read', replace=True)
    barcodes_detected = subsampled['cell_barcode'].nunique()
    umis_detected = subsampled.groupby('cell_barcode')['umi'].nunique().median()
    genes_detected = subsampled.groupby('cell_barcode')['umi'].nunique().median()
    return pd.DataFrame({
        'depth': [depth],
        'barcodes_detected': [barcodes_detected],
        'median_umis': [umis_detected],
        'median_genes': [genes_detected]
    })


def find_elbow(x, y):
    kneedle = KneeLocator(x, y, S=1.0, curve='concave', direction='increasing')
    return kneedle.elbow


def main():
    args = parse_arguments()

    # Load the CSV file
    df = pd.read_csv(args.csv_path)

    # Define subsampling depths
    total_reads = df['read'].sum()
    depths = np.logspace(2, np.log10(total_reads), num=20).astype(int)
    depths = np.unique(depths)

    # Perform subsampling and analysis
    results = pd.concat([subsample_and_analyze(df, d) for d in depths])

    # Find saturation points
    barcode_saturation = find_elbow(results['depth'], results['barcodes_detected'])
    umi_saturation = find_elbow(results['depth'], results['median_umis'])
    gene_saturation = find_elbow(results['depth'], results['median_genes'])

    # Plot the results
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(20, 6))

    sns.lineplot(data=results, x='depth', y='barcodes_detected', ax=ax1)
    ax1.axvline(x=barcode_saturation, color='r', linestyle='--', label=f'Saturation at {barcode_saturation:.0f}')
    ax1.set_xscale('log')
    ax1.set_xlabel('Sequencing Depth (Reads)')
    ax1.set_ylabel('Unique Barcodes Detected')
    ax1.set_title('Barcode Saturation')
    ax1.legend()

    sns.lineplot(data=results, x='depth', y='median_umis', ax=ax2)
    ax2.axvline(x=umi_saturation, color='r', linestyle='--', label=f'Saturation at {umi_saturation:.0f}')
    ax2.set_xscale('log')
    ax2.set_xlabel('Sequencing Depth (Reads)')
    ax2.set_ylabel('Median UMIs per Cell')
    ax2.set_title('UMI Saturation')
    ax2.legend()

    sns.lineplot(data=results, x='depth', y='median_genes', ax=ax3)
    ax3.axvline(x=gene_saturation, color='r', linestyle='--', label=f'Saturation at {gene_saturation:.0f}')
    ax3.set_xscale('log')
    ax3.set_xlabel('Sequencing Depth (Reads)')
    ax3.set_ylabel('Median Genes per Cell')
    ax3.set_title('Gene Saturation')
    ax3.legend()

    plt.tight_layout()

    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(args.csv_path)
    os.makedirs(output_dir, exist_ok=True)

    # Save the plot in the same directory as the input CSV
    plot_path = os.path.join(output_dir, 'saturation_plot.png')
    plt.savefig(plot_path)
    print(f"Saturation plot saved to: {plot_path}")

    # Calculate and print statistics
    max_barcodes = results['barcodes_detected'].max()
    max_umis = results['median_umis'].max()
    max_genes = results['median_genes'].max()

    barcode_saturation_percentage = (np.interp(barcode_saturation, results['depth'],
                                               results['barcodes_detected']) / max_barcodes) * 100
    umi_saturation_percentage = (np.interp(umi_saturation, results['depth'], results['median_umis']) / max_umis) * 100
    gene_saturation_percentage = (np.interp(gene_saturation, results['depth'],
                                            results['median_genes']) / max_genes) * 100

    total_cells = df['cell_barcode'].nunique()

    print(
        f"\nPredicted Barcode Saturation Point: {barcode_saturation:.0f} reads ({barcode_saturation_percentage:.1f}% of max)")
    print(f"Predicted UMI Saturation Point: {umi_saturation:.0f} reads ({umi_saturation_percentage:.1f}% of max)")
    print(f"Predicted Gene Saturation Point: {gene_saturation:.0f} reads ({gene_saturation_percentage:.1f}% of max)")
    print(f"\nReads per cell at Barcode Saturation: {barcode_saturation / total_cells:.0f}")
    print(f"Reads per cell at UMI Saturation: {umi_saturation / total_cells:.0f}")
    print(f"Reads per cell at Gene Saturation: {gene_saturation / total_cells:.0f}")


if __name__ == "__main__":
    main()