#!/usr/bin/env python3
"""
big-FISH Multi-FOV Processing Pipeline

Calculates spliceforms per cell.
Processes multiple fields of view (FOVs) with automatic z-slice discovery.
Supports optional Blue (DAPI) channel - when absent, uses cell-only segmentation.

Usage:
    python bigfish_pipeline.py --input_dir /path/to/input --output_dir /path/to/output --channels Red Gold70 FarRed
"""

import os
import re
import sys
import csv
import argparse
import traceback
import numpy as np
from pathlib import Path

# scipy
from scipy import ndimage

# scikit-image
from skimage.filters import threshold_otsu, gaussian, sobel
from skimage.morphology import remove_small_objects, binary_dilation, disk, closing
from skimage.measure import label, regionprops
from skimage.feature import peak_local_max
from skimage.segmentation import watershed

# big-FISH
import bigfish
import bigfish.stack as stack
import bigfish.multistack as multistack
import bigfish.segmentation as segmentation
import bigfish.detection as detection
import bigfish.plot as plot

# matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def discover_z_slices(path_input, fov, channels):
    """
    Auto-discover all z-slices available for a given FOV.

    Parameters
    ----------
    path_input : str
        Directory containing TIF files
    fov : str
        FOV identifier
    channels : list
        List of channel names to validate

    Returns
    -------
    list or None
        Sorted list of z-slice identifiers, or None if channels missing
    """
    z_slices = set()
    z_pattern = re.compile(r'(Z\d+)')

    for filename in os.listdir(path_input):
        if fov in filename and filename.endswith('.tif'):
            z_match = z_pattern.search(filename)
            if z_match:
                z_slices.add(z_match.group(1))

    if not z_slices:
        return None

    return sorted(list(z_slices))


def validate_fov_channels(path_input, fov, channels, opt_prefix="Camera - "):
    """
    Validate that all required channels exist for a given FOV.

    Parameters
    ----------
    path_input : str
        Directory containing TIF files
    fov : str
        FOV identifier
    channels : list
        Required channel names
    opt_prefix : str
        Optical configuration prefix in filename

    Returns
    -------
    tuple
        (bool, list) - (all_present, missing_channels)
    """
    files_in_dir = set(os.listdir(path_input))
    missing_channels = []

    for channel in channels:
        # Check if at least one file exists for this FOV+channel combination
        pattern_str = "{}.*{}{}.*\\.tif$".format(
            re.escape(fov),
            re.escape(opt_prefix),
            re.escape(channel)
        )
        channel_pattern = re.compile(pattern_str)
        channel_found = any(channel_pattern.search(f) for f in files_in_dir)

        if not channel_found:
            missing_channels.append(channel)

    return len(missing_channels) == 0, missing_channels


def discover_fovs(path_input, fov_pattern=r'(\d{8}_[A-Z0-9]+_[a-z]+_slide\d+_\d+)'):
    """
    Auto-discover all unique FOVs in the directory.

    Parameters
    ----------
    path_input : str
        Directory containing TIF files
    fov_pattern : str
        Regex pattern matching FOV identifier format

    Returns
    -------
    list
        Sorted list of unique FOV identifiers
    """
    fovs = set()
    pattern = re.compile(fov_pattern)

    for filename in os.listdir(path_input):
        if filename.endswith('.tif'):
            match = pattern.search(filename)
            if match:
                fovs.add(match.group(1))

    if not fovs:
        raise ValueError("No FOVs found in directory: {}".format(path_input))

    return sorted(list(fovs))


def spots_to_structured(spots):
    """
    Convert spots array to structured array for exact matching.

    Parameters
    ----------
    spots : np.ndarray
        Spots array with shape (n_spots, n_coords)

    Returns
    -------
    np.recarray
        Structured array for efficient set operations
    """
    if spots.shape[1] == 2:
        dtype = [('y', spots.dtype), ('x', spots.dtype)]
    elif spots.shape[1] == 3:
        dtype = [('y', spots.dtype), ('x', spots.dtype), ('z', spots.dtype)]
    else:
        raise ValueError("Unexpected number of columns: {}".format(spots.shape[1]))
    return np.core.records.fromarrays(spots.T, dtype=dtype)


def safe_colocalization(spots_1, spots_2, voxel_size):
    """
    Wrapper for detect_spots_colocalization that handles empty arrays.

    Parameters
    ----------
    spots_1 : np.ndarray
        First spots array (n_spots, 2)
    spots_2 : np.ndarray
        Second spots array (n_spots, 2)
    voxel_size : tuple
        Voxel dimensions for distance calculation

    Returns
    -------
    tuple
        (coloc_1, coloc_2, distances) - colocalized spots from each array
    """
    # Handle empty input arrays
    if spots_1.size == 0 or spots_2.size == 0:
        # Return empty arrays with correct shape
        empty_1 = np.empty((0, spots_1.shape[1] if spots_1.ndim > 1 else 2), dtype=spots_1.dtype)
        empty_2 = np.empty((0, spots_2.shape[1] if spots_2.ndim > 1 else 2), dtype=spots_2.dtype)
        return empty_1, empty_2, np.array([])

    # Handle case with too few spots for reliable colocalization
    min_spots_for_colocalization = 3
    if spots_1.shape[0] < min_spots_for_colocalization or spots_2.shape[0] < min_spots_for_colocalization:
        print("  Warning: Too few spots for colocalization ({} vs {}), returning empty".format(
            spots_1.shape[0], spots_2.shape[0]))
        empty_1 = np.empty((0, spots_1.shape[1] if spots_1.ndim > 1 else 2), dtype=spots_1.dtype)
        empty_2 = np.empty((0, spots_2.shape[1] if spots_2.ndim > 1 else 2), dtype=spots_2.dtype)
        return empty_1, empty_2, np.array([])

    # Call original function with error handling for edge cases
    try:
        return multistack.detect_spots_colocalization(
            spots_1=spots_1,
            spots_2=spots_2,
            voxel_size=voxel_size,
            return_indices=False
        )
    except ValueError as e:
        if "False is not in list" in str(e):
            # Elbow detection failed - likely due to unusual distance distribution
            # Return empty colocalization result
            print("  Warning: Colocalization threshold detection failed, returning empty")
            empty_1 = np.empty((0, spots_1.shape[1] if spots_1.ndim > 1 else 2), dtype=spots_1.dtype)
            empty_2 = np.empty((0, spots_2.shape[1] if spots_2.ndim > 1 else 2), dtype=spots_2.dtype)
            return empty_1, empty_2, np.array([])
        else:
            # Re-raise unexpected errors
            raise


def robust_remove_transcription_site(spots, clusters, nuc_label, image_mip, ndim=2,
                                     ts_radius=15):
    """
    Remove transcription sites with fallback for failed nuclear segmentation (e.g. no DAPI).

    Parameters
    ----------
    spots : np.ndarray
        Spots with cluster assignments, shape (n_spots, 3) for 2D [y, x, cluster_id]
    clusters : np.ndarray
        Cluster information, shape (n_clusters, 4) for 2D [y, x, n_spots, cluster_id]
    nuc_label : np.ndarray
        Nuclear segmentation mask (can be all zeros if segmentation failed)
    image_mip : np.ndarray
        Maximum intensity projection of the channel (used for fallback)
    ndim : int
        Number of spatial dimensions (2 or 3)
    ts_radius : float
        Radius in pixels for the circular TS region around brightest point (default: 15)

    Returns
    -------
    tuple
        (spots_no_ts, foci, ts_spots) - spots after TS removal, foci array, TS spots
    """
    # Check if nuclear segmentation succeeded
    n_nuclei = len(np.unique(nuc_label)) - 1

    if n_nuclei > 0:
        # Nuclear segmentation succeeded - use standard bigfish function
        return multistack.remove_transcription_site(spots, clusters, nuc_label, ndim=ndim)

    # Fallback: nuclear segmentation failed, use intensity-based detection
    print("  Note: Nuclear segmentation failed, using intensity-based TS detection")

    # Handle empty inputs
    if spots.size == 0 or clusters.size == 0:
        empty_spots = np.empty((0, spots.shape[1] if spots.ndim > 1 else 3), dtype=spots.dtype)
        empty_foci = np.empty((0, 2), dtype=np.int64)
        return spots.copy() if spots.size > 0 else empty_spots, empty_foci, empty_spots

    # Find the brightest point in the image (likely the transcription site)
    smoothed = ndimage.gaussian_filter(image_mip.astype(float), sigma=2)
    max_idx = np.unravel_index(np.argmax(smoothed), smoothed.shape)
    ts_center_y, ts_center_x = max_idx

    # Create circular binary mask around the brightest point
    y_coords, x_coords = np.ogrid[:image_mip.shape[0], :image_mip.shape[1]]
    distance_from_center = np.sqrt((y_coords - ts_center_y) ** 2 + (x_coords - ts_center_x) ** 2)
    ts_mask = (distance_from_center <= ts_radius).astype(np.uint8)

    # Use bigfish identify_objects_in_region to find clusters within the TS mask
    cluster_coords = clusters[:, :2]
    clusters_in_ts, clusters_outside_ts = multistack.identify_objects_in_region(
        ts_mask, cluster_coords, ndim=2
    )

    # Get cluster IDs that are in the TS region
    # Match cluster coordinates back to their IDs
    ts_cluster_ids = []
    for ts_cluster_coord in clusters_in_ts:
        # Find matching cluster in original array
        for i in range(clusters.shape[0]):
            if np.allclose(clusters[i, :2], ts_cluster_coord):
                cluster_id = int(clusters[i, -1]) if clusters.shape[1] > 3 else i
                ts_cluster_ids.append(cluster_id)
                break

    # Get the cluster assignment column index in spots array
    cluster_col = 2 if ndim == 2 else 3

    # Separate spots into TS and non-TS based on cluster membership
    if len(ts_cluster_ids) > 0:
        ts_spot_mask = np.isin(spots[:, cluster_col], ts_cluster_ids)
        spots_no_ts = spots[~ts_spot_mask]
        ts_spots = spots[ts_spot_mask]
        foci = np.array([[ts_center_y, ts_center_x]], dtype=np.int64)
    else:
        spots_no_ts = spots.copy()
        ts_spots = np.empty((0, spots.shape[1]), dtype=spots.dtype)
        foci = np.empty((0, 2), dtype=np.int64)

    print("  Intensity-based TS removal: TS center at ({}, {}), radius={}, {} clusters in TS, {} spots removed".format(
        ts_center_y, ts_center_x, ts_radius, len(ts_cluster_ids), ts_spots.shape[0]))

    return spots_no_ts, foci, ts_spots


def load_and_stack_fov(path_input, fov, z_slices, channels, opt_prefix="Camera - "):
    """
    Load and stack images for a single FOV.

    Parameters
    ----------
    path_input : str
        Directory containing TIF files
    fov : str
        FOV identifier
    z_slices : list
        List of z-slice identifiers (e.g., ['Z01', 'Z02', ...])
    channels : list
        List of channel names
    opt_prefix : str
        Optical configuration prefix in filename

    Returns
    -------
    np.ndarray
        Stacked image array with shape [n_channels, n_z, height, width]
    """
    stacks = []

    for channel in channels:
        channel_stack = []
        for z in z_slices:
            # Build filename pattern
            pattern_str = "{}.*{}{}.*{}.*\\.tif$".format(
                re.escape(fov),
                re.escape(opt_prefix),
                re.escape(channel),
                re.escape(z)
            )
            pattern = re.compile(pattern_str)

            # Find matching file
            matching_files = [f for f in os.listdir(path_input) if pattern.search(f)]

            if len(matching_files) == 0:
                raise FileNotFoundError("No file found for {} {} {} in {}".format(
                    fov, channel, z, path_input))
            elif len(matching_files) > 1:
                print("Warning: Multiple files match {} {} {}, using first".format(fov, channel, z))

            # Load image
            filepath = os.path.join(path_input, matching_files[0])
            img = stack.read_image(filepath)
            channel_stack.append(img)

        # Stack z-slices for this channel
        channel_array = np.stack(channel_stack, axis=0)
        stacks.append(channel_array)

    # Stack all channels
    return np.stack(stacks, axis=0)


def segment_cells_by_intensity(cytoplasm_mip, nuc_label=None, min_cell_area=4000,
                               max_cell_area=50000, threshold_method='hybrid',
                               dilation_radius=10, min_circularity=0.4,
                               bright_sigma_threshold=20, background_sigma=50,
                               contrast_percentile=97, smoothing_sigma=3):
    """
    Two-stage hybrid cell segmentation for images with mixed bright and dim cells.

    This approach handles both:
    - Sparse images with few very bright cells (uses MAD-based thresholding)
    - Dense images with many dim cells (uses local contrast enhancement)

    Parameters
    ----------
    cytoplasm_mip : np.ndarray
        Maximum intensity projection of cytoplasmic channel (can be combined signal)
    nuc_label : np.ndarray or None
        Nuclear segmentation labels (used for matching if available)
    min_cell_area : int
        Minimum cell area in pixels (default: 4000)
    max_cell_area : int
        Maximum cell area in pixels (default: 50000)
    threshold_method : str
        'hybrid' (recommended): two-stage bright + contrast detection
        'mad': MAD-based only (for sparse bright cells)
        'contrast': local contrast only (for uniform brightness cells)
        'otsu': legacy Otsu thresholding (not recommended)
    dilation_radius : int
        Pixels to dilate final cell masks (default: 10)
    min_circularity : float
        Minimum circularity (4*pi*area/perimeter^2) to accept a region (default: 0.4)
        Values: 0-1, where 1.0 is a perfect circle
    bright_sigma_threshold : float
        N-sigma threshold for detecting extremely bright cells in Stage 1 (default: 20)
        Higher values = more stringent, only the brightest cells detected
    background_sigma : float
        Gaussian sigma for background estimation in Stage 2 (default: 50)
        Should be larger than cell diameter
    contrast_percentile : float
        Percentile threshold for local contrast detection in Stage 2 (default: 97)
        Lower values = more cells detected
    smoothing_sigma : float
        Gaussian smoothing sigma for noise reduction (default: 3)

    Returns
    -------
    tuple
        (nuc_label, cell_label, cell_intensities) - nuclear labels, cell labels,
        and dict mapping cell_id to mean intensity

    """
    import bigfish.multistack as multistack

    img = cytoplasm_mip.astype(float)
    cell_intensities = {}

    # ========================================================================
    # Stage 1: Detect extremely bright cells using MAD-based threshold
    # ========================================================================
    smoothed = gaussian(img, sigma=smoothing_sigma)
    median_val = np.median(smoothed)
    mad = np.median(np.abs(smoothed - median_val))
    sigma_est = 1.4826 * mad

    bright_thresh = median_val + bright_sigma_threshold * sigma_est
    print("  Stage 1 (bright cells): median={:.1f}, MAD_sigma={:.1f}, threshold={:.1f} ({}-sigma)".format(
        median_val, sigma_est, bright_thresh, bright_sigma_threshold))

    bright_mask = smoothed > bright_thresh
    bright_mask = remove_small_objects(bright_mask, min_size=min_cell_area // 2)
    bright_mask = closing(bright_mask, disk(3))
    bright_labeled = label(bright_mask)

    n_bright = len(np.unique(bright_labeled)) - 1
    print("  Stage 1: {} extremely bright region(s) detected".format(n_bright))

    # ========================================================================
    # Stage 2: Detect remaining cells using local contrast enhancement
    # ========================================================================
    if threshold_method in ['hybrid', 'contrast']:
        # Subtract large-scale background
        background = gaussian(img, sigma=background_sigma)
        local_contrast = img - background
        lc_smooth = gaussian(local_contrast, sigma=smoothing_sigma)

        # Percentile-based threshold on local contrast
        contrast_thresh = np.percentile(lc_smooth, contrast_percentile)
        print("  Stage 2 (local contrast): background_sigma={}, threshold={:.1f} (p{})".format(
            background_sigma, contrast_thresh, contrast_percentile))

        contrast_mask = lc_smooth > contrast_thresh
        contrast_mask = remove_small_objects(contrast_mask, min_size=min_cell_area // 2)
        contrast_mask = closing(contrast_mask, disk(5))
    else:
        contrast_mask = np.zeros_like(img, dtype=bool)
        print("  Stage 2: skipped (threshold_method={})".format(threshold_method))

    # ========================================================================
    # Combine masks and filter by size/circularity
    # ========================================================================
    if threshold_method == 'mad':
        # MAD-only mode: use only bright mask
        combined_mask = bright_mask
    elif threshold_method == 'contrast':
        # Contrast-only mode: use only contrast mask
        combined_mask = contrast_mask
    else:
        # Hybrid mode: combine both
        combined_mask = bright_mask | contrast_mask

    combined_labeled = label(combined_mask)
    props = regionprops(combined_labeled, intensity_image=img)

    # Filter by area and circularity
    cell_label = np.zeros_like(img, dtype=int)
    cell_id = 0

    for prop in props:
        # Area filter
        if prop.area < min_cell_area:
            continue
        if prop.area > max_cell_area:
            print("  Rejected region: area={} > max_cell_area={}".format(prop.area, max_cell_area))
            continue

        # Circularity filter: 4*pi*area / perimeter^2
        circ = (4 * np.pi * prop.area) / (prop.perimeter ** 2) if prop.perimeter > 0 else 0
        if circ < min_circularity:
            continue

        cell_id += 1
        cell_label[combined_labeled == prop.label] = cell_id
        cell_intensities[cell_id] = prop.mean_intensity

    n_cells = cell_id

    # Optional dilation (per-cell to avoid merging)
    if dilation_radius > 0 and n_cells > 0:
        dilated = np.zeros_like(cell_label)
        for i in range(1, cell_label.max() + 1):
            cell_mask = cell_label == i
            cell_dilated = binary_dilation(cell_mask, disk(dilation_radius))
            dilated[cell_dilated & (dilated == 0)] = i
        cell_label = dilated

    # Handle nuclear labels
    if nuc_label is None:
        nuc_label = np.zeros_like(cell_label)
        print("  Hybrid segmentation: {} cells detected".format(n_cells))
    else:
        # Match nuclei to cells if both exist
        if nuc_label.max() > 0 and cell_label.max() > 0:
            nuc_label, cell_label = multistack.match_nuc_cell(
                nuc_label, cell_label, single_nuc=False, cell_alone=True
            )
            n_cells = len(np.unique(cell_label)) - 1
        print("  Hybrid segmentation: {} cells, {} nuclei".format(
            n_cells, len(np.unique(nuc_label)) - 1))

    return nuc_label, cell_label, cell_intensities


# ============================================================================
# DIAGNOSTIC VISUALIZATION FUNCTIONS
# ============================================================================

def save_spot_diagnostic_images(path_output, pol_mip, env_mip, nef_mip,
                                pol_no_ts, env_no_ts, nef_no_ts,
                                nuc_label, cell_label, fov,
                                ts_radius=15):
    """
    Save diagnostic images showing detected spots overlaid on MIPs.
    """
    # Check if nuclear segmentation failed
    nuc_seg_failed = nuc_label.max() == 0

    # Create figure with 3 subplots (one per channel)
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    title_suffix = " [Intensity-based TS]" if nuc_seg_failed else ""
    fig.suptitle('Detected Spots (after TS removal) - {}{}'.format(fov, title_suffix), fontsize=14)

    # Channel names and colors
    channel_info = [
        ('pol (Red)', pol_mip, pol_no_ts, 'cyan'),
        ('env (Gold70)', env_mip, env_no_ts, 'magenta'),
        ('nef (FarRed)', nef_mip, nef_no_ts, 'yellow')
    ]

    for ax, (name, mip, spots, color) in zip(axes, channel_info):
        # Display MIP with auto-contrast
        vmin, vmax = np.percentile(mip, [1, 99.5])
        ax.imshow(mip, cmap='gray', vmin=vmin, vmax=vmax)

        # Overlay spots
        if spots.size > 0:
            # spots are in (y, x) format
            ax.scatter(spots[:, 1], spots[:, 0], c=color, s=8, alpha=0.7,
                       marker='o', linewidths=0.5, edgecolors='white')

        # Add cell boundary contours if cells exist (check for any non-zero values)
        cell_mask = cell_label > 0
        if np.any(cell_mask):
            ax.contour(cell_mask, colors='lime', linewidths=0.5, alpha=0.5)

        # Add nuclear boundary contours if nuclei exist
        nuc_mask = nuc_label > 0
        if np.any(nuc_mask):
            ax.contour(nuc_mask, colors='blue', linewidths=0.5, alpha=0.5)
        elif nuc_seg_failed:
            # Show circular TS region around brightest point when nuc seg failed
            smoothed = ndimage.gaussian_filter(mip.astype(float), sigma=2)
            max_idx = np.unravel_index(np.argmax(smoothed), smoothed.shape)
            ts_center_y, ts_center_x = max_idx

            # Draw circle around TS center
            circle = plt.Circle((ts_center_x, ts_center_y), ts_radius,
                                fill=False, color='red', linewidth=2, linestyle='--')
            ax.add_patch(circle)
            # Mark the center
            ax.plot(ts_center_x, ts_center_y, 'r+', markersize=10, markeredgewidth=2)

        ax.set_title('{}\n{} spots'.format(name, spots.shape[0] if spots.size > 0 else 0))
        ax.axis('off')

    # Add legend
    legend_elements = [
        mpatches.Patch(facecolor='none', edgecolor='lime', label='Cell boundary'),
        mpatches.Patch(facecolor='none', edgecolor='blue', label='Nuclear boundary'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
                   markersize=8, label='Detected spots')
    ]
    if nuc_seg_failed:
        legend_elements.append(
            plt.Line2D([0], [0], marker='o', color='red', markerfacecolor='none',
                       markersize=10, linestyle='--', label='TS region (r={})'.format(ts_radius))
        )
    fig.legend(handles=legend_elements, loc='lower center', ncol=4 if nuc_seg_failed else 3, fontsize=10)

    plt.tight_layout()
    plt.subplots_adjust(bottom=0.12)

    # Save figure
    output_path = os.path.join(path_output, 'spot_detection_diagnostic.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)


# ============================================================================
# PROCESSING FUNCTIONS
# ============================================================================

def process_fov(fov, path_input, path_output_base, channels, pipeline_params, has_dapi=True):
    """
    Process a single FOV through the complete big-FISH pipeline.

    Workflow:
    1. Validate all channels present
    2. Auto-discover z-slices
    3. Load and stack images
    4. Segment nuclei (if DAPI available) and cells
    5. Detect spots in each channel
    6. Cluster spots into foci
    7. Classify splicing patterns (US/SS/MS)
    8. Save results

    Parameters
    ----------
    fov : str
        FOV identifier
    path_input : str
        Input directory containing TIF files
    path_output_base : str
        Base output directory
    channels : list
        List of channel names (3 or 4 channels)
    pipeline_params : dict
        Pipeline parameters
    has_dapi : bool
        Whether Blue/DAPI channel is present

    Returns
    -------
    dict or None
        Processing summary, or None if FOV was skipped
    """
    print("\n" + "=" * 70)
    print("Processing FOV: {}".format(fov))
    print('=' * 70)

    # Validate channels
    all_present, missing = validate_fov_channels(path_input, fov, channels)
    if not all_present:
        print("WARNING: Skipping FOV {} - missing channels: {}".format(fov, missing))
        return None

    # Discover z-slices
    z_slices = discover_z_slices(path_input, fov, channels)
    if z_slices is None or len(z_slices) == 0:
        print("WARNING: No z-slices found for FOV {}".format(fov))
        return None

    print("Found {} z-slices: {}".format(len(z_slices), z_slices))

    # Load images
    print("Loading images...")
    image_stack = load_and_stack_fov(path_input, fov, z_slices, channels)
    print("Image stack shape: {} [channels, z, y, x]".format(image_stack.shape))

    # Separate channels based on configuration
    if has_dapi:
        # 4-channel mode: Blue, Red, Gold70, FarRed
        dapi = image_stack[0]  # Blue
        pol = image_stack[1]  # Red
        env = image_stack[2]  # Gold70
        nef = image_stack[3]  # FarRed
    else:
        # 3-channel mode: Red, Gold70, FarRed (no DAPI)
        dapi = None
        pol = image_stack[0]  # Red
        env = image_stack[1]  # Gold70
        nef = image_stack[2]  # FarRed

    # Create output directory
    path_output = os.path.join(path_output_base, fov)
    os.makedirs(path_output, exist_ok=True)

    # Segmentation
    print("Segmenting nuclei and cells...")
    seg_params = pipeline_params['segmentation']

    # Create MIPs
    nef_mip = stack.maximum_projection(nef)
    pol_mip = stack.maximum_projection(pol)
    env_mip = stack.maximum_projection(env)

    # Get segmentation parameters
    min_cell_area = seg_params.get('min_cell_area', 4000)
    max_cell_area = seg_params.get('max_cell_area', 50000)
    threshold_method = seg_params.get('threshold_method', 'hybrid')
    bright_sigma_threshold = seg_params.get('bright_sigma_threshold', 20)
    background_sigma = seg_params.get('background_sigma', 50)
    contrast_percentile = seg_params.get('contrast_percentile', 97)
    min_circularity = seg_params.get('min_circularity', 0.4)
    smoothing_sigma = seg_params.get('smoothing_sigma', 3)
    dilation_radius = seg_params.get('cell_dilation_radius', 20)

    # Create combined FISH signal for cell detection
    def normalize_image(img):
        img_min, img_max = img.min(), img.max()
        if img_max > img_min:
            return (img - img_min) / (img_max - img_min)
        return np.zeros_like(img, dtype=float)

    combined_fish = (normalize_image(pol_mip.astype(float)) +
                     normalize_image(env_mip.astype(float)) +
                     normalize_image(nef_mip.astype(float)))

    if has_dapi and dapi is not None:
        # Standard workflow with nuclear segmentation
        dapi_mip = stack.maximum_projection(dapi)

        # Use Otsu thresholding instead of fixed threshold for nuclei
        nuc_threshold = seg_params.get('nuc_threshold', None)
        if nuc_threshold is None or nuc_threshold == 'otsu':
            nuc_threshold = int(threshold_otsu(dapi_mip))  # Convert to native int for bigfish
            print("  Nuclear threshold (Otsu): {}".format(nuc_threshold))
        else:
            nuc_threshold = int(nuc_threshold)  # Ensure native int
            print("  Nuclear threshold (manual): {}".format(nuc_threshold))

        nuc_mask = segmentation.thresholding(dapi_mip, threshold=nuc_threshold)
        nuc_mask = segmentation.clean_segmentation(
            nuc_mask,
            small_object_size=seg_params['nuc_small_object_size'],
            fill_holes=True
        )
        nuc_label = segmentation.label_instances(nuc_mask)

        # Cell segmentation using combined FISH signal with hybrid approach
        nuc_label, cell_label, cell_intensities = segment_cells_by_intensity(
            combined_fish,
            nuc_label,
            min_cell_area=min_cell_area,
            max_cell_area=max_cell_area,
            threshold_method=threshold_method,
            dilation_radius=dilation_radius,
            min_circularity=min_circularity,
            bright_sigma_threshold=bright_sigma_threshold,
            background_sigma=background_sigma,
            contrast_percentile=contrast_percentile,
            smoothing_sigma=smoothing_sigma
        )
    else:
        # No DAPI: use cell-only segmentation
        print("No DAPI channel - using cell-only segmentation")
        nuc_label, cell_label, cell_intensities = segment_cells_by_intensity(
            combined_fish,
            nuc_label=None,
            min_cell_area=min_cell_area,
            max_cell_area=max_cell_area,
            threshold_method=threshold_method,
            dilation_radius=dilation_radius,
            min_circularity=min_circularity,
            bright_sigma_threshold=bright_sigma_threshold,
            background_sigma=background_sigma,
            contrast_percentile=contrast_percentile,
            smoothing_sigma=smoothing_sigma
        )

    # Count nuclei and cells
    n_nuclei = len(np.unique(nuc_label)) - 1  # subtract 1 for background
    n_cells = len(np.unique(cell_label)) - 1

    print("Segmented {} nuclei and {} cells".format(n_nuclei, n_cells))

    # Spot detection
    print("Detecting spots...")
    spot_params = pipeline_params['spot_detection']

    spots_pol, threshold_pol = detection.detect_spots(
        images=pol_mip,
        return_threshold=True,
        voxel_size=spot_params['voxel_size'],
        spot_radius=spot_params['spot_radius']
    )
    print("pol: {} spots (threshold: {})".format(spots_pol.shape[0], threshold_pol))

    spots_env, threshold_env = detection.detect_spots(
        images=env_mip,
        return_threshold=True,
        voxel_size=spot_params['voxel_size'],
        spot_radius=spot_params['spot_radius']
    )
    print("env: {} spots (threshold: {})".format(spots_env.shape[0], threshold_env))

    spots_nef, threshold_nef = detection.detect_spots(
        images=nef_mip,
        return_threshold=True,
        voxel_size=spot_params['voxel_size'],
        spot_radius=spot_params['spot_radius']
    )
    print("nef: {} spots (threshold: {})".format(spots_nef.shape[0], threshold_nef))

    # Mask spots by cell segmentation (only if cells were detected)
    if cell_label.max() > 0:
        spots_pol, _ = multistack.identify_objects_in_region(cell_label, spots_pol, ndim=2)

    # Cluster detection
    print("Decomposing dense regions and clustering foci...")
    cluster_params = pipeline_params['clustering']

    # pol channel
    pol_post_decomposition, _, _ = detection.decompose_dense(
        image=pol_mip,
        spots=spots_pol,
        voxel_size=spot_params['voxel_size'],
        spot_radius=spot_params['spot_radius'],
        alpha=cluster_params['alpha'],
        beta=cluster_params['beta'],
        gamma=cluster_params['gamma']
    )
    pol_post_clustering, clusters_pol = detection.detect_clusters(
        spots=pol_post_decomposition,
        voxel_size=spot_params['voxel_size'],
        radius=cluster_params['cluster_radius'],
        nb_min_spots=cluster_params['nb_min_spots']
    )

    # env channel
    env_post_decomposition, _, _ = detection.decompose_dense(
        image=env_mip,
        spots=spots_env,
        voxel_size=spot_params['voxel_size'],
        spot_radius=spot_params['spot_radius'],
        alpha=cluster_params['alpha'],
        beta=cluster_params['beta'],
        gamma=cluster_params['gamma']
    )
    env_post_clustering, clusters_env = detection.detect_clusters(
        spots=env_post_decomposition,
        voxel_size=spot_params['voxel_size'],
        radius=cluster_params['cluster_radius'],
        nb_min_spots=cluster_params['nb_min_spots']
    )

    # nef channel
    nef_post_decomposition, _, _ = detection.decompose_dense(
        image=nef_mip,
        spots=spots_nef,
        voxel_size=spot_params['voxel_size'],
        spot_radius=spot_params['spot_radius'],
        alpha=cluster_params['alpha'],
        beta=cluster_params['beta'],
        gamma=cluster_params['gamma']
    )
    nef_post_clustering, clusters_nef = detection.detect_clusters(
        spots=nef_post_decomposition,
        voxel_size=spot_params['voxel_size'],
        radius=cluster_params['cluster_radius'],
        nb_min_spots=cluster_params['nb_min_spots']
    )

    print("pol: {} clusters".format(clusters_pol.shape[0] if clusters_pol.size > 0 else 0))
    print("env: {} clusters".format(clusters_env.shape[0] if clusters_env.size > 0 else 0))
    print("nef: {} clusters".format(clusters_nef.shape[0] if clusters_nef.size > 0 else 0))

    # Transcription site removal (uses intensity-based fallback when nuc_label is empty)
    print("Removing transcription sites...")
    ts_params = pipeline_params.get('transcription_site', {'ts_radius': 15})
    ts_radius = ts_params.get('ts_radius', 15)

    pol_no_ts, pol_foci, pol_ts = robust_remove_transcription_site(
        pol_post_clustering, clusters_pol, nuc_label, pol_mip, ndim=2, ts_radius=ts_radius
    )
    env_no_ts, env_foci, env_ts = robust_remove_transcription_site(
        env_post_clustering, clusters_env, nuc_label, env_mip, ndim=2, ts_radius=ts_radius
    )
    nef_no_ts, nef_foci, nef_ts = robust_remove_transcription_site(
        nef_post_clustering, clusters_nef, nuc_label, nef_mip, ndim=2, ts_radius=ts_radius
    )

    print("After TS removal: pol={}, env={}, nef={}".format(
        pol_no_ts.shape[0] if pol_no_ts.size > 0 else 0,
        env_no_ts.shape[0] if env_no_ts.size > 0 else 0,
        nef_no_ts.shape[0] if nef_no_ts.size > 0 else 0
    ))

    # Colocalization analysis
    print("Performing colocalization analysis...")
    coloc_params = pipeline_params['colocalization']
    voxel_size = coloc_params['voxel_size']

    # Extract just coordinates (first 2 columns) for colocalization
    pol_coords = pol_no_ts[:, :2] if pol_no_ts.size > 0 else np.empty((0, 2), dtype=np.int64)
    env_coords = env_no_ts[:, :2] if env_no_ts.size > 0 else np.empty((0, 2), dtype=np.int64)
    nef_coords = nef_no_ts[:, :2] if nef_no_ts.size > 0 else np.empty((0, 2), dtype=np.int64)

    # Classify spots: US (unspliced), SS (singly spliced), MS (multiply spliced)
    #   US (unspliced): nef + env + pol (full-length transcript, all three probes bind)
    #   SS (singly spliced): nef + env only (pol intron removed)
    #   MS (multiply spliced): nef only (both introns removed)
    print("Classifying splicing patterns...")

    # Step 1: Find nef spots colocalized with env (these are US + SS combined)
    USSS_nef, USSS_env, _ = safe_colocalization(nef_coords, env_coords, voxel_size)

    # Step 2: From nef+env spots, find those also colocalized with pol (these are US)
    US_nef, US_pol, _ = safe_colocalization(USSS_nef, pol_coords, voxel_size)

    print("Colocalization: nef+env={}, nef+env+pol={}".format(
        USSS_nef.shape[0], US_nef.shape[0]))

    # Step 3: SS = nef+env spots that are NOT also colocalized with pol
    USSS_struct = spots_to_structured(USSS_nef) if USSS_nef.size > 0 else spots_to_structured(
        np.empty((0, 2), dtype=np.int64))
    US_struct = spots_to_structured(US_nef) if US_nef.size > 0 else spots_to_structured(
        np.empty((0, 2), dtype=np.int64))

    SS_mask = ~np.isin(USSS_struct, US_struct)
    SS_spots = USSS_nef[SS_mask] if USSS_nef.size > 0 else np.empty((0, 2), dtype=np.int64)

    # Step 4: MS = nef spots that are NOT colocalized with env
    nef_struct = spots_to_structured(nef_coords) if nef_coords.size > 0 else spots_to_structured(
        np.empty((0, 2), dtype=np.int64))

    MS_mask = ~np.isin(nef_struct, USSS_struct)
    MS_spots = nef_coords[MS_mask] if nef_coords.size > 0 else np.empty((0, 2), dtype=np.int64)

    # US spots are the nef coordinates from the triple colocalization
    US_spots = US_nef if US_nef.size > 0 else np.empty((0, 2), dtype=np.int64)

    print("Classification: US={}, SS={}, MS={}".format(
        US_spots.shape[0], SS_spots.shape[0], MS_spots.shape[0]))

    # Save diagnostic images
    save_spot_diagnostic_images(
        path_output, pol_mip, env_mip, nef_mip,
        pol_coords, env_coords, nef_coords,
        nuc_label, cell_label, fov, ts_radius=ts_radius
    )

    # Per-cell statistics
    cell_results = []
    if cell_label.max() > 0:
        cell_ids = [i for i in np.unique(cell_label) if i > 0]
        for cell_id in cell_ids:
            cell_mask = cell_label == cell_id
            cell_area = np.sum(cell_mask)

            # Nuclear area (if available)
            if nuc_label.max() > 0:
                nuc_mask = nuc_label == cell_id
                nuc_area = np.sum(nuc_mask)
            else:
                nuc_area = 0

            # Count spots in this cell
            def count_spots_in_cell(spots, mask):
                if spots.size == 0:
                    return 0
                in_cell = mask[spots[:, 0].astype(int), spots[:, 1].astype(int)]
                return np.sum(in_cell)

            n_US = count_spots_in_cell(US_spots, cell_mask)
            n_SS = count_spots_in_cell(SS_spots, cell_mask)
            n_MS = count_spots_in_cell(MS_spots, cell_mask)
            n_pol = count_spots_in_cell(pol_coords, cell_mask)
            n_env = count_spots_in_cell(env_coords, cell_mask)
            n_nef = count_spots_in_cell(nef_coords, cell_mask)

            cell_results.append({
                'fov': fov,
                'cell_id': cell_id,
                'cell_area': cell_area,
                'nuc_area': nuc_area,
                'n_US': n_US,
                'n_SS': n_SS,
                'n_MS': n_MS,
                'n_pol': n_pol,
                'n_env': n_env,
                'n_nef': n_nef,
                'n_total_spots': n_US + n_SS + n_MS
            })

    # Save per-FOV CSV
    if cell_results:
        csv_path = os.path.join(path_output, 'cell_stats.csv')
        fieldnames = ['fov', 'cell_id', 'cell_area', 'nuc_area',
                      'n_US', 'n_SS', 'n_MS', 'n_pol', 'n_env', 'n_nef', 'n_total_spots']
        with open(csv_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(cell_results)

    # Return summary
    return {
        'fov': fov,
        'n_nuclei': n_nuclei,
        'n_cells': n_cells,
        'cell_results': cell_results,
        'fov_totals': {
            'n_US': US_spots.shape[0],
            'n_SS': SS_spots.shape[0],
            'n_MS': MS_spots.shape[0]
        }
    }


def run_pipeline(path_input, path_output_base, channels, pipeline_params, fov_list=None, fov_pattern=None):
    """
    Run the pipeline on all FOVs.

    Parameters
    ----------
    path_input : str
        Input directory
    path_output_base : str
        Output directory
    channels : list
        Channel names (3 or 4 channels)
    pipeline_params : dict
        Pipeline parameters
    fov_list : list, optional
        Specific FOVs to process
    fov_pattern : str, optional
        Regex pattern for FOV discovery

    Returns
    -------
    tuple
        (results, errors, skipped)
    """
    # Determine if DAPI is present
    has_dapi = 'Blue' in channels

    # Validate channel count
    expected_count = 4 if has_dapi else 3
    if len(channels) != expected_count:
        if has_dapi:
            raise ValueError("With Blue channel, expected 4 channels, got {}: {}".format(len(channels), channels))
        else:
            raise ValueError("Without Blue channel, expected 3 channels, got {}: {}".format(len(channels), channels))

    # Discover FOVs if not provided
    if fov_list is None:
        if fov_pattern is None:
            fov_pattern = r'(\d{8}_[A-Z0-9]+_[a-z]+_slide\d+_\d+)'
        fov_list = discover_fovs(path_input, fov_pattern)

    print("\n" + "=" * 70)
    print("big-FISH Multi-FOV Pipeline")
    print('=' * 70)
    print("Input directory: {}".format(path_input))
    print("Output directory: {}".format(path_output_base))
    print("Channels: {}".format(channels))
    print("DAPI present: {}".format(has_dapi))
    print("Found {} FOVs to process".format(len(fov_list)))
    print('=' * 70)

    results = []
    errors = []
    skipped = []

    for i, fov in enumerate(fov_list, 1):
        print("\n\nProcessing FOV {}/{}".format(i, len(fov_list)))
        try:
            result = process_fov(fov, path_input, path_output_base, channels, pipeline_params, has_dapi=has_dapi)
            if result is None:
                skipped.append(fov)
            else:
                results.append(result)
        except Exception as e:
            error_msg = "Error processing FOV {}: {}".format(fov, str(e))
            print("\nERROR: {}".format(error_msg))
            traceback.print_exc()
            errors.append({'fov': fov, 'error': str(e)})
            continue

    # Print summary
    print("\n" + "=" * 70)
    print("PIPELINE SUMMARY")
    print("=" * 70)
    print("Total FOVs found: {}".format(len(fov_list)))
    print("Successfully processed: {}".format(len(results)))
    print("Skipped (missing channels): {}".format(len(skipped)))
    print("Errors: {}".format(len(errors)))

    if results:
        # Collect all cell results across FOVs
        all_cell_results = []
        for r in results:
            all_cell_results.extend(r.get('cell_results', []))

        # Save combined cell statistics CSV (append if exists, create if not)
        if all_cell_results:
            combined_csv_path = os.path.join(path_output_base, 'all_cells_stats.csv')
            fieldnames = ['fov', 'cell_id', 'cell_area', 'nuc_area',
                          'n_US', 'n_SS', 'n_MS', 'n_pol', 'n_env', 'n_nef', 'n_total_spots']

            # Check if file exists to determine write mode and header behavior
            file_exists = os.path.isfile(combined_csv_path)

            with open(combined_csv_path, 'a', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                # Only write header if file is new
                if not file_exists:
                    writer.writeheader()
                writer.writerows(all_cell_results)

            if file_exists:
                print("\nAppended {} rows to: {}".format(len(all_cell_results), combined_csv_path))
            else:
                print("\nCreated combined cell statistics: {}".format(combined_csv_path))

        # Print per-cell summary
        print("\n" + "-" * 70)
        print("PER-CELL RESULTS")
        print("-" * 70)
        print("{:<40} {:>6} {:>6} {:>6} {:>6} {:>8}".format(
            "FOV / Cell", "US", "SS", "MS", "Total", "Area"))
        print("-" * 70)

        for r in results:
            fov = r['fov']
            cell_results = r.get('cell_results', [])

            if cell_results:
                for cell in cell_results:
                    cell_name = "{}/cell_{}".format(fov, cell['cell_id'])
                    print("{:<40} {:>6} {:>6} {:>6} {:>6} {:>8}".format(
                        cell_name[:40],
                        cell['n_US'],
                        cell['n_SS'],
                        cell['n_MS'],
                        cell['n_total_spots'],
                        cell['cell_area']
                    ))
            else:
                print("{:<40} {:>6}".format(fov, "no cells"))

        # Print FOV-level summary
        print("\n" + "-" * 70)
        print("FOV SUMMARY")
        print("-" * 70)
        total_cells = sum(r['n_cells'] for r in results)
        total_US = sum(r['fov_totals']['n_US'] for r in results)
        total_SS = sum(r['fov_totals']['n_SS'] for r in results)
        total_MS = sum(r['fov_totals']['n_MS'] for r in results)

        print("Total cells across all FOVs: {}".format(total_cells))
        print("Total spots - US: {}, SS: {}, MS: {}".format(total_US, total_SS, total_MS))

        for r in results:
            fov_totals = r.get('fov_totals', {})
            print("\n{}:".format(r['fov']))
            print("  Cells: {}, Nuclei: {}".format(r['n_cells'], r['n_nuclei']))
            print("  FOV totals - US: {}, SS: {}, MS: {}".format(
                fov_totals.get('n_US', 0),
                fov_totals.get('n_SS', 0),
                fov_totals.get('n_MS', 0)))

    if skipped:
        print("\nSkipped FOVs ({}):".format(len(skipped)))
        for fov in skipped:
            print("  {}".format(fov))

    if errors:
        print("\nErrors encountered ({}):".format(len(errors)))
        for e in errors:
            print("  {}: {}".format(e['fov'], e['error']))

    print("=" * 70)

    return results, errors, skipped


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main entry point for command-line execution"""
    parser = argparse.ArgumentParser(
        description='big-FISH Multi-FOV Processing Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage (with DAPI):
  python bigfish_pipeline.py \\
    --input_dir /path/to/TIFF \\
    --output_dir /path/to/output \\
    --channels Blue Red Gold70 FarRed

Example usage (without DAPI):
  python bigfish_pipeline.py \\
    --input_dir /path/to/TIFF \\
    --output_dir /path/to/output \\
    --channels Red Gold70 FarRed

When DAPI is absent:
  - Nuclear segmentation is skipped
  - Cell segmentation uses Otsu thresholding on FarRed channel (nef)
  - Transcription site removal uses intensity-based fallback
        """
    )

    parser.add_argument(
        '--input_dir',
        type=str,
        required=True,
        help='Input directory containing TIF files'
    )

    parser.add_argument(
        '--output_dir',
        type=str,
        required=True,
        help='Output directory for processed results'
    )

    parser.add_argument(
        '--channels',
        type=str,
        nargs='+',
        default=["Blue", "Red", "Gold70", "FarRed"],
        help='Channel names in order. Use 4 channels with DAPI (Blue Red Gold70 FarRed) or 3 without (Red Gold70 FarRed)'
    )

    parser.add_argument(
        '--fov_pattern',
        type=str,
        default=None,
        help='Regex pattern for FOV discovery (optional)'
    )

    args = parser.parse_args()

    # Validate directories
    if not os.path.isdir(args.input_dir):
        print("ERROR: Input directory does not exist: {}".format(args.input_dir))
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    # Pipeline parameters
    pipeline_params = {
        'segmentation': {
            # Nuclear segmentation (DAPI) - use 'otsu' for adaptive or a number for fixed
            'nuc_threshold': 'otsu',
            'nuc_small_object_size': 2000,
            # Cell segmentation - hybrid two-stage approach
            # Stage 1: Detect extremely bright cells using MAD-based threshold
            # Stage 2: Detect remaining cells using local contrast enhancement
            'threshold_method': 'hybrid',  # 'hybrid', 'mad', 'contrast', or 'otsu'
            'bright_sigma_threshold': 20,  # N-sigma for bright cell detection
            'background_sigma': 50,        # Gaussian sigma for background estimation
            'contrast_percentile': 97,     # Percentile threshold for contrast detection
            'min_circularity': 0.4,        # Minimum circularity (0-1, circle=1)
            'min_cell_area': 4000,          # Minimum cell size in pixels
            'max_cell_area': 50000,        # Maximum cell size in pixels
            'cell_dilation_radius': 20,    # Dilation of final masks
            'smoothing_sigma': 3,          # Gaussian smoothing for noise reduction
        },
        'spot_detection': {
            'voxel_size': (107, 107),
            'spot_radius': (130, 130)
        },
        'clustering': {
            'alpha': 0.7,
            'beta': 1,
            'gamma': 5,
            'cluster_radius': 600,
            'nb_min_spots': 6
        },
        'colocalization': {
            'voxel_size': (107, 107)
        },
        'transcription_site': {
            'ts_radius': 15
        }
    }

    # Run pipeline
    results, errors, skipped = run_pipeline(
        path_input=args.input_dir,
        path_output_base=args.output_dir,
        channels=args.channels,
        pipeline_params=pipeline_params,
        fov_pattern=args.fov_pattern
    )

    # Exit with appropriate code
    if errors:
        sys.exit(1)
    else:
        sys.exit(0)


if __name__ == "__main__":
    main()