#!/usr/bin/env python3
"""
Visualize Cell Masks Colored by Total Spot Count

Creates colored mask images where each cell is colored according to its
n_total_spots value.

Usage:
    python bigfish_visualize.py \
        --mask_dir /path/to/masks \
        --stats_csv /path/to/all_cells_stats.csv \
        --output_dir /path/to/output
"""

# !/usr/bin/env python3

import os
import re
import argparse
import numpy as np
import pandas as pd
import tifffile
from matplotlib import cm
from matplotlib.colors import Normalize


def find_mask_for_fov(mask_dir, fov):
    """
    Find the corresponding mask file for a given FOV.

    Mask files are expected to be named like:
    ROIs2Label_20251119_EGS023_plate_021_XY1_Camera - FarRed_Z07.tif

    The FOV (e.g., '20251119_EGS023_plate_021_XY1') appears between
    'ROIs2Label_' and '_Camera'.

    Parameters
    ----------
    mask_dir : str
        Directory containing mask TIF files
    fov : str
        FOV identifier to match (e.g., '20251119_EGS023_plate_021_XY1')

    Returns
    -------
    str or None
        Path to mask file, or None if not found
    """
    for filename in os.listdir(mask_dir):
        if not filename.startswith('ROIs2Label_'):
            continue

        # Extract the portion between 'ROIs2Label_' and '_Camera'
        # This should exactly match the FOV
        remainder = filename[len('ROIs2Label_'):]

        # Find where '_Camera' starts
        camera_idx = remainder.find('_Camera')
        if camera_idx == -1:
            continue

        extracted_fov = remainder[:camera_idx]

        if extracted_fov == fov:
            return os.path.join(mask_dir, filename)

    return None


def load_cell_mask(mask_path):
    """
    Load a cell segmentation mask (label image) from a TIF file.

    Parameters
    ----------
    mask_path : str
        Path to the mask TIF file

    Returns
    -------
    np.ndarray
        Label image where each cell has a unique integer value, background is 0
    """
    mask = tifffile.imread(mask_path)

    # Ensure it's 2D
    if mask.ndim > 2:
        mask = mask[0] if mask.ndim == 3 else mask[0, 0]
        print("  Warning: Mask was >2D, using first slice")

    # Ensure integer type
    if not np.issubdtype(mask.dtype, np.integer):
        mask = mask.astype(np.int32)

    return mask


def create_spot_count_colormap():
    """
    Create a custom rainbow colormap with bright blue (low) to bright red (high).

    Uses a brighter blue (R50, G70, B200) at the low end to stand out against
    black background, transitioning through cyan, green, yellow to bright red.

    Returns
    -------
    matplotlib.colors.LinearSegmentedColormap
        Custom colormap for mapping normalized values to colors
    """
    from matplotlib.colors import LinearSegmentedColormap

    # Define color stops: bright blue -> cyan -> green -> yellow -> bright red
    # Colors specified as RGB tuples normalized to [0, 1]
    colors = [
        (50 / 255, 70 / 255, 200 / 255),  # Bright blue (0.0) - visible against black
        (0 / 255, 200 / 255, 255 / 255),  # Cyan (0.25)
        (0 / 255, 230 / 255, 0 / 255),  # Green (0.5)
        (255 / 255, 230 / 255, 0 / 255),  # Yellow (0.75)
        (255 / 255, 50 / 255, 0 / 255),  # Bright red (1.0)
    ]

    # Create colormap with evenly spaced color stops
    cmap = LinearSegmentedColormap.from_list('bright_jet', colors, N=256)
    return cmap


def color_mask_by_spots(mask, cell_spot_counts, vmax=None):
    """
    Color a label mask by spot counts using a rainbow colormap.

    Parameters
    ----------
    mask : np.ndarray
        2D label image where each cell has a unique integer ID (background=0)
    cell_spot_counts : dict
        Mapping of cell_id -> n_total_spots
    vmax : float, optional
        Maximum value for color scaling. If None, uses max of spot counts.

    Returns
    -------
    np.ndarray
        RGB image with shape (height, width, 3), dtype uint8
    """
    height, width = mask.shape

    # Initialize output as black (RGB)
    colored = np.zeros((height, width, 3), dtype=np.uint8)

    # Determine color scale
    if cell_spot_counts:
        counts = list(cell_spot_counts.values())
        if vmax is None:
            vmax = max(counts) if counts else 1
    else:
        vmax = 1

    # Ensure vmax is at least 1 to avoid division by zero
    vmax = max(vmax, 1)

    # Create colormap and normalizer
    cmap = create_spot_count_colormap()
    norm = Normalize(vmin=0, vmax=vmax)

    # Get unique cell IDs (excluding background)
    cell_ids = np.unique(mask)
    cell_ids = cell_ids[cell_ids > 0]

    # Color each cell
    for cell_id in cell_ids:
        # Get spot count for this cell (default to 0 if not found)
        spot_count = cell_spot_counts.get(cell_id, 0)

        # Normalize and get color
        normalized_value = norm(spot_count)
        rgba = cmap(normalized_value)

        # Convert to uint8 RGB (ignore alpha)
        rgb = (np.array(rgba[:3]) * 255).astype(np.uint8)

        # Apply color to all pixels belonging to this cell
        cell_pixels = mask == cell_id
        colored[cell_pixels] = rgb

    return colored


def process_fov(fov, mask_dir, stats_df, output_dir, global_vmax=None):
    """
    Process a single FOV: load mask, color by spots, save as TIFF.

    Parameters
    ----------
    fov : str
        FOV identifier
    mask_dir : str
        Directory containing mask files
    stats_df : pd.DataFrame
        DataFrame with cell statistics
    output_dir : str
        Output directory for colored images
    global_vmax : float, optional
        Global maximum for consistent color scaling across FOVs

    Returns
    -------
    dict or None
        Processing result with statistics, or None if mask not found
    """
    # Find mask file
    mask_path = find_mask_for_fov(mask_dir, fov)
    if mask_path is None:
        print("  Mask not found for FOV: {}".format(fov))
        return None

    print("Processing FOV: {}".format(fov))
    print("  Mask: {}".format(os.path.basename(mask_path)))

    # Load mask
    mask = load_cell_mask(mask_path)
    print("  Mask shape: {}, {} cells".format(mask.shape, len(np.unique(mask)) - 1))

    # Get spot counts for this FOV
    fov_stats = stats_df[stats_df['fov'] == fov]
    cell_spot_counts = dict(zip(fov_stats['cell_id'], fov_stats['n_total_spots']))

    print("  Found spot counts for {} cells".format(len(cell_spot_counts)))
    if cell_spot_counts:
        counts = list(cell_spot_counts.values())
        print("  Spot range: {} - {} (mean: {:.1f})".format(
            min(counts), max(counts), np.mean(counts)))

    # Color the mask
    colored = color_mask_by_spots(mask, cell_spot_counts, vmax=global_vmax)

    # Save as TIFF
    output_filename = "{}_colored_by_spots.tif".format(fov)
    output_path = os.path.join(output_dir, output_filename)

    # Save as RGB TIFF using tifffile
    # Transpose to (3, H, W) for proper RGB TIFF format
    tifffile.imwrite(output_path, colored, photometric='rgb')
    print("  Saved: {}".format(output_filename))

    return {
        'fov': fov,
        'n_cells': len(cell_spot_counts),
        'min_spots': min(cell_spot_counts.values()) if cell_spot_counts else 0,
        'max_spots': max(cell_spot_counts.values()) if cell_spot_counts else 0,
        'output_path': output_path
    }


def create_colorbar_image(vmax, output_path, height=400, width=50):
    """
    Create and save a standalone colorbar image for reference.

    Parameters
    ----------
    vmax : float
        Maximum value for the colorbar
    output_path : str
        Path to save the colorbar image
    height : int
        Height of the colorbar in pixels
    width : int
        Width of the colorbar in pixels
    """
    import matplotlib.pyplot as plt
    from matplotlib.colorbar import ColorbarBase
    from matplotlib.colors import Normalize

    fig, ax = plt.subplots(figsize=(1.5, 6))

    cmap = create_spot_count_colormap()
    norm = Normalize(vmin=0, vmax=vmax)

    cb = ColorbarBase(ax, cmap=cmap, norm=norm, orientation='vertical')
    cb.set_label('Number of total spots', fontsize=12)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print("Saved colorbar: {}".format(output_path))


def main():
    parser = argparse.ArgumentParser(
        description='Visualize cell masks colored by total spot count',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
  python visualize_masks_by_spots.py \\
    --mask_dir /path/to/masks \\
    --stats_csv /path/to/all_cells_stats.csv \\
    --output_dir /path/to/output

The output will be TIFF images with cells colored by their n_total_spots value:
  - Blue = 0 spots
  - Red = maximum spots
  - Black = background (no cell mask)
        """
    )

    parser.add_argument(
        '--mask_dir',
        type=str,
        required=True,
        help='Directory containing ROIs2Label_*.tif mask files'
    )

    parser.add_argument(
        '--stats_csv',
        type=str,
        required=True,
        help='Path to all_cells_stats.csv from bigfish pipeline'
    )

    parser.add_argument(
        '--output_dir',
        type=str,
        required=True,
        help='Output directory for colored mask images'
    )

    parser.add_argument(
        '--global_scale',
        action='store_true',
        help='Use global max across all FOVs for consistent color scaling'
    )

    parser.add_argument(
        '--vmax',
        type=float,
        default=None,
        help='Manual maximum value for color scale (overrides auto-detection)'
    )

    args = parser.parse_args()

    # Validate inputs
    if not os.path.isdir(args.mask_dir):
        print("ERROR: Mask directory does not exist: {}".format(args.mask_dir))
        return 1

    if not os.path.isfile(args.stats_csv):
        print("ERROR: Stats CSV does not exist: {}".format(args.stats_csv))
        return 1

    os.makedirs(args.output_dir, exist_ok=True)

    # Load statistics
    print("Loading statistics from: {}".format(args.stats_csv))
    stats_df = pd.read_csv(args.stats_csv)

    required_cols = ['fov', 'cell_id', 'n_total_spots']
    missing_cols = [c for c in required_cols if c not in stats_df.columns]
    if missing_cols:
        print("ERROR: Missing required columns: {}".format(missing_cols))
        return 1

    print("Loaded {} cell records from {} FOVs".format(
        len(stats_df), stats_df['fov'].nunique()))

    # Determine color scale
    if args.vmax is not None:
        global_vmax = args.vmax
        print("Using manual vmax: {}".format(global_vmax))
    elif args.global_scale:
        global_vmax = stats_df['n_total_spots'].max()
        print("Using global vmax: {}".format(global_vmax))
    else:
        global_vmax = None
        print("Using per-FOV scaling")

    # Get unique FOVs
    fovs = stats_df['fov'].unique()
    print("\nProcessing {} FOVs...".format(len(fovs)))
    print("=" * 60)

    results = []
    for fov in fovs:
        result = process_fov(
            fov, args.mask_dir, stats_df, args.output_dir,
            global_vmax
        )
        if result:
            results.append(result)

    # Save colorbar reference
    if results:
        # Use the maximum from processed data for colorbar
        actual_max = max(r['max_spots'] for r in results)
        colorbar_vmax = global_vmax if global_vmax else actual_max

        colorbar_path = os.path.join(args.output_dir, 'colorbar_reference.png')
        create_colorbar_image(colorbar_vmax, colorbar_path)

    # Print summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print("FOVs processed: {}".format(len(results)))

    if results:
        total_cells = sum(r['n_cells'] for r in results)
        overall_max = max(r['max_spots'] for r in results)
        print("Total cells: {}".format(total_cells))
        print("Max spots in any cell: {}".format(overall_max))
        print("\nOutput files saved to: {}".format(args.output_dir))

    return 0


if __name__ == "__main__":
    exit(main())