#!/usr/bin/env python3

"""
Splice Site Classification and Filtering Tool

Description:
    Classifies donor and acceptor splice sites based on genomic positions and
    filters results based on specified criteria. Designed to work with output
    from bamsort_splice.py.

Input:
    - CSV file: Output from bamsort_splice.py containing splice junction data
    - Donor/Acceptor reference positions: Defined in script or via command line

Output:
    - Filtered CSV with classified splice sites

Filtering:
    - Removes duplicate entries (CB, DonorSite, AcceptorSite, UB)
    - Filters for specific donor sites (default: D4)
    - Adds Target annotation

Usage:
    python3 filter_splice_sites.py input_dir output_dir [options]
"""

import pandas as pd
import numpy as np
import argparse
import os
import sys
from pathlib import Path


def classify_site(position, reference_positions, reference_labels):
    """
    Classify a genomic position to the nearest reference site.
    
    Args:
        position: Genomic position to classify (int or float)
        reference_positions: List of reference positions
        reference_labels: Corresponding labels for reference positions
    
    Returns:
        - Reference label (e.g., "D1", "A2") if within 3bp of reference
        - "mnc" (may not be canonical) if 4-12bp from nearest reference
        - "nc" (not canonical) if >12bp from nearest reference or invalid input
    
    Technical rationale:
        Three-tier classification system:
        1. ≤3bp: Canonical splice site (accounts for minor alignment shifts)
        2. 4-12bp: May not be canonical (potential cryptic/alternative sites)
        3. >12bp: Not canonical (too far from known sites)
        
        This approach is more nuanced than binary classification, allowing
        identification of potential cryptic splice sites while maintaining
        confidence in canonical site assignments.
    """
    # Handle NA/missing values
    if pd.isna(position) or position == "":
        return "nc"
    
    try:
        position = float(position)
    except (ValueError, TypeError):
        return "nc"
    
    # Calculate distances to all reference positions
    distances = np.abs(np.array(reference_positions) - position)
    
    # Filter out any NA distances (shouldn't happen with valid inputs)
    valid_mask = ~np.isnan(distances)
    valid_distances = distances[valid_mask]
    
    if len(valid_distances) == 0:
        return "nc"
    
    # Find minimum distance and corresponding label
    min_distance = np.min(valid_distances)
    min_idx = np.argmin(distances)
    closest_label = reference_labels[min_idx]
    
    # Three-tier classification
    if min_distance <= 3:
        return str(closest_label)
    elif min_distance <= 12:
        return "mnc"
    else:
        return "nc"


def process_sample(input_file, output_file, donors_df, acceptors_df, 
                   target_name="D1", filter_donor="D4"):
    """
    Process a single sample file: classify sites, filter, and deduplicate.
    
    Technical implementation:
        1. Read CSV with error handling for missing/corrupted files
        2. Apply three-tier classification:
           - Canonical (≤3bp from reference)
           - May not be canonical (4-12bp from reference)
           - Not canonical (>12bp from reference)
        3. Use pandas' drop_duplicates for memory-efficient deduplication
        4. Filter in single pass to minimize memory footprint
    
    Args:
        input_file: Path to input CSV from bamsort_splice.py
        output_file: Path for output CSV
        donors_df: DataFrame with 'Donor' positions and 'DonorSite' labels
        acceptors_df: DataFrame with 'Acceptor' positions and 'AcceptorSite' labels
        target_name: Target annotation to add (default: "D1")
        filter_donor: Donor site to filter for (default: "D4")
    
    Returns:
        Number of records written, or None if error
    """
    try:
        # Read input file
        if not os.path.exists(input_file):
            print(f"  WARNING: File not found: {input_file}")
            return None
        
        data = pd.read_csv(input_file)
        
        if len(data) == 0:
            print(f"  WARNING: Empty file: {input_file}")
            return None
        
        # Classify donor sites
        # Using apply here is acceptable since we need custom logic per row
        # Alternative vectorized approaches would require complex numpy operations
        data['DonorSite'] = data['Donor'].apply(
            lambda x: classify_site(x, donors_df['Donor'].tolist(), 
                                   donors_df['DonorSite'].tolist())
        )
        
        # Classify acceptor sites
        data['AcceptorSite'] = data['Acceptor'].apply(
            lambda x: classify_site(x, acceptors_df['Acceptor'].tolist(), 
                                   acceptors_df['AcceptorSite'].tolist())
        )
        
        # Select relevant columns
        # Note: Using select instead of drop to be explicit about output schema
        file_data = data[['CB', 'UB', 'DonorSite', 'AcceptorSite']].copy()
        
        # Remove duplicates
        # keep='first' retains the first occurrence, which preserves original order
        file_data = file_data.drop_duplicates(
            subset=['CB', 'DonorSite', 'AcceptorSite', 'UB'], 
            keep='first'
        )
        
        # Add target annotation
        file_data['Target'] = target_name
        
        # Filter for specific donor site
        if filter_donor:
            file_data = file_data[file_data['DonorSite'] == filter_donor]
        
        # Write output
        file_data.to_csv(output_file, index=False)
        
        return len(file_data)
        
    except Exception as e:
        print(f"  ERROR processing {input_file}: {e}", file=sys.stderr)
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Classify and filter splice sites from bamsort_splice.py output",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Reference splice site positions (default for mac239):
  Acceptors: A1=4658, A2=5141, A3=5746, A4=5818, A5=5959, A7=8249
  Donors: D1=431, D2=4730, D3=5217, D4=6043

Classification rules:
  - Distance ≤3bp: Canonical site (returns site label, e.g., "D4", "A2")
  - Distance 4-12bp: May not be canonical (returns "mnc")
  - Distance >12bp: Not canonical (returns "nc")

Examples:
  # Process all samples with default settings
  python3 filter_splice_sites.py /input/dir /output/dir

  # Filter for different donor site
  python3 filter_splice_sites.py /input/dir /output/dir --filter-donor D2

  # Process specific samples only
  python3 filter_splice_sites.py /input/dir /output/dir --samples D0 D195 W2
        """
    )
    
    parser.add_argument("input_dir", 
                       help="Directory containing splice site CSV files from bamsort_splice.py")
    parser.add_argument("output_dir", 
                       help="Directory for filtered output CSV files")
    parser.add_argument("--samples", nargs='+', 
                       default=["D0", "D195", "Invitro", "Pacute", "W0", "W2", "D83"],
                       help="Sample names to process (default: D0 D195 Invitro Pacute W0 W2 D83)")
    parser.add_argument("--targets", nargs='+',
                       default=["_D1_PL", "_env_PL", "_pol_PL", "_SSenv_PL"],
                       help="Target suffixes to process (default: _D1_PL _env_PL _pol_PL _SSenv_PL)")
    parser.add_argument("--filter-donor", type=str, default="D4",
                       help="Donor site to filter for (default: D4)")
    parser.add_argument("--no-filter", action="store_true",
                       help="Skip donor site filtering (keep all sites)")
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Define reference splice sites
    # These positions are specific to mac239 HIV reference genome
    # Source: HIV-1 splice site annotations from standard references
    acceptors_df = pd.DataFrame({
        'Acceptor': [4658, 5141, 5746, 5818, 5959, 8249],
        'AcceptorSite': ['A1', 'A2', 'A3', 'A4', 'A5', 'A7']
    })
    
    donors_df = pd.DataFrame({
        'Donor': [431, 4730, 5217, 6043],
        'DonorSite': ['D1', 'D2', 'D3', 'D4']
    })
    
    print(f"Input directory: {args.input_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Samples to process: {', '.join(args.samples)}")
    print(f"Targets to process: {', '.join(args.targets)}")
    print(f"Donor site filter: {args.filter_donor if not args.no_filter else 'None (all sites)'}")
    print(f"Classification: ≤3bp=canonical, 4-12bp=mnc, >12bp=nc")
    print("-" * 60)
    
    # Statistics
    total_samples = 0
    processed_samples = 0
    total_records = 0
    
    # Process each sample and target combination
    for sample_name in args.samples:
        for target_name in args.targets:
            total_samples += 1
            
            # Construct input filename
            input_filename = f"{sample_name}{target_name}_splicesites.csv"
            input_file = os.path.join(args.input_dir, input_filename)
            
            # Construct output filename
            # Extract target prefix (e.g., "D1" from "_D1_PL")
            target_prefix = target_name.lstrip('_').split('_')[0]
            output_filename = f"{sample_name}{target_name}_splicesites_filtered.csv"
            output_file = os.path.join(args.output_dir, output_filename)
            
            print(f"\nProcessing: {sample_name}{target_name}")
            
            # Process the file
            filter_donor = None if args.no_filter else args.filter_donor
            num_records = process_sample(
                input_file, output_file, donors_df, acceptors_df,
                target_prefix, filter_donor
            )
            
            if num_records is not None:
                processed_samples += 1
                total_records += num_records
                print(f"  Wrote {num_records} records to {output_filename}")
    
    # Summary
    print("\n" + "=" * 60)
    print("Processing Summary:")
    print(f"Total sample/target combinations: {total_samples}")
    print(f"Successfully processed: {processed_samples}")
    print(f"Failed/skipped: {total_samples - processed_samples}")
    print(f"Total filtered records written: {total_records}")
    print(f"Output directory: {args.output_dir}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"Fatal error: {e}", file=sys.stderr)
        sys.exit(1)
