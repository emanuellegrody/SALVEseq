import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
import argparse
import os


def parse_arguments():
    parser = argparse.ArgumentParser(description='Perform saturation analysis on sparse single-cell data.')
    parser.add_argument('csv_path', type=str, help='Path to the input CSV file')
    parser.add_argument('--output_dir', type=str, default='.', help='Directory to save the output plot')
    return parser.parse_args()


def subsample_and_analyze(data, depth):
    subsampled = data.sample(n=depth, replace=True)
    umi_counts = subsampled.groupby('cell_barcode').size()
    return umi_counts


def plot_umi_distribution(umi_counts, ax, depth):
    kde = gaussian_kde(np.log1p(umi_counts))
    x_range = np.linspace(0, np.log1p(umi_counts.max()), 1000)
    ax.plot(x_range, kde(x_range), label=f'Depth: {depth}')
    ax.set_xlabel('log(UMIs + 1)')
    ax.set_ylabel('Density')
    ax.legend()


def main():
    args = parse_arguments()

    # Load your data
    df = pd.read_csv(args.csv_path)

    # Define subsampling depths
    total_reads = len(df)
    depths = np.logspace(2, np.log10(total_reads), num=10).astype(int)

    # Prepare plots
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(20, 20))

    # UMI distribution plot
    for depth in depths:
        umi_counts = subsample_and_analyze(df, depth)
        plot_umi_distribution(umi_counts, ax1, depth)
    ax1.set_title('UMI Distribution at Different Sequencing Depths')

    # Cumulative UMI curve
    cumulative_umis = [subsample_and_analyze(df, depth).sum() for depth in depths]
    ax2.plot(depths, cumulative_umis)
    ax2.set_xscale('log')
    ax2.set_xlabel('Sequencing Depth')
    ax2.set_ylabel('Cumulative UMIs')
    ax2.set_title('Cumulative UMIs vs Sequencing Depth')

    # High-UMI cell focus
    high_umi_threshold = np.percentile(df.groupby('cell_barcode').size(), 90)  # Top 10% of cells
    high_umi_cells = df.groupby('cell_barcode').size()[lambda x: x > high_umi_threshold].index
    high_umi_umis = [subsample_and_analyze(df[df['cell_barcode'].isin(high_umi_cells)], depth).mean() for depth in
                     depths]
    ax3.plot(depths, high_umi_umis)
    ax3.set_xscale('log')
    ax3.set_xlabel('Sequencing Depth')
    ax3.set_ylabel('Mean UMIs in High-UMI Cells')
    ax3.set_title('Saturation in High-UMI Cells')

    # Fraction of cells above threshold
    umi_threshold = 10  # Adjust based on your data
    fraction_above_threshold = [
        (subsample_and_analyze(df, depth) > umi_threshold).mean() for depth in depths
    ]
    ax4.plot(depths, fraction_above_threshold)
    ax4.set_xscale('log')
    ax4.set_xlabel('Sequencing Depth')
    ax4.set_ylabel(f'Fraction of Cells with > {umi_threshold} UMIs')
    ax4.set_title('Fraction of Cells Above UMI Threshold')

    plt.tight_layout()

    # Save the plot
    output_file = os.path.join(args.output_dir, 'sparse_data_saturation_analysis.png')
    plt.savefig(output_file)
    print(f"Saturation analysis plot saved to: {output_file}")


if __name__ == "__main__":
    main()