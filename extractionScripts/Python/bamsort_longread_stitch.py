#!/usr/bin/env python3

"""
UMI Fragment Stitching - Step 2 (with Isoform Classification)

Description:
    Takes the output from bamsort_longread_fragments.py and performs the actual
    reconstruction of original molecules by merging complementary UMI fragments.
    Additionally classifies each molecule into isoforms (US, MS, SS, or any).

Algorithm Overview:
    1. Load candidate pairs from _candidates.csv
    2. Build a graph where UMIs are nodes and COMPLEMENTARY relationships are edges
    3. Find connected components - each component represents one original molecule
    4. For each component, compute the merged alignment span
    5. Classify isoform based on merged span coverage
    6. Output reconstructed molecules with isoform assignments

Isoform Classification Logic (MAC239 coordinates):
    US (Unspliced): Coverage within region 985-5212 (≥50bp)
    MS (Multiply Spliced): Coverage before 6597, gap 6597-8806, coverage after 8806
    SS (Singly Spliced): Coverage before 985, gap 985-5212, coverage 6597-8806
    any: All other patterns

Output Files:
    1. _molecules.csv: Reconstructed molecules with component and isoform assignments
    2. _molecule_chains.csv: Ordered fragment chains within each molecule
    3. _stitching_stats.csv: Statistics on stitching success and isoform distribution

Usage:
    python3 bamsort_longread_stitch.py candidates.csv pairs.csv output_prefix [--tolerance 5]
"""

import csv
import sys
import argparse
from collections import defaultdict

# Isoform classification regions (MAC239 coordinates)
REGION_US_START = 985
REGION_US_END = 5212
REGION_MS_GAP_START = 6597
REGION_MS_GAP_END = 8806
REGION_SS_GAP_START = 985
REGION_SS_GAP_END = 5212
REGION_SS_AFTER_START = 6597
REGION_SS_AFTER_END = 8806

MIN_ALIGNMENT_LENGTH = 50  # Minimum bp for US classification


class UnionFind:
    """
    Union-Find data structure for efficient connected component detection.

    Uses path compression and union by rank for O(alpha(n)) amortized operations,
    where alpha is the inverse Ackermann function (effectively constant).
    """

    def __init__(self):
        self.parent = {}
        self.rank = {}

    def make_set(self, x):
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0

    def find(self, x):
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])  # Path compression
        return self.parent[x]

    def union(self, x, y):
        px, py = self.find(x), self.find(y)
        if px == py:
            return
        # Union by rank
        if self.rank[px] < self.rank[py]:
            px, py = py, px
        self.parent[py] = px
        if self.rank[px] == self.rank[py]:
            self.rank[px] += 1

    def get_components(self):
        """Return dict mapping root -> list of members."""
        components = defaultdict(list)
        for x in self.parent:
            components[self.find(x)].append(x)
        return dict(components)


def has_coverage_in_region(start, end, region_start, region_end,
                           min_length=1, tolerance=5):
    """
    Check if a span overlaps with target region within tolerance.

    Args:
        start, end: Span boundaries
        region_start: Start of target region
        region_end: End of target region
        min_length: Minimum overlap length required (bp)
        tolerance: +/- bp tolerance for boundaries

    Returns:
        bool: True if sufficient overlap exists
    """
    region_start_tol = region_start - tolerance
    region_end_tol = region_end + tolerance

    overlap_start = max(start, region_start_tol)
    overlap_end = min(end, region_end_tol)
    overlap_length = max(0, overlap_end - overlap_start)

    return overlap_length >= min_length


def has_gap_in_region(start, end, gap_start, gap_end, tolerance=5):
    """
    Check if span does NOT cover the specified gap region (within tolerance).

    For a molecule to have a "gap" in a region, its span should not extend
    continuously through that region. Since we only have span boundaries
    (not detailed CIGAR), we check if the span covers the gap region.

    Args:
        start, end: Span boundaries
        gap_start: Start of expected gap region
        gap_end: End of expected gap region
        tolerance: +/- bp tolerance for boundaries

    Returns:
        bool: True if gap exists (span does not fully cover the region)
    """
    gap_start_tol = gap_start + tolerance
    gap_end_tol = gap_end - tolerance

    if gap_start_tol >= gap_end_tol:
        return True

    # If span fully contains the gap region, there is no gap
    if start <= gap_start_tol and end >= gap_end_tol:
        return False

    return True


def classify_isoform_from_span(merged_start, merged_end, tolerance=5):
    """
    Classify molecule into isoform based on merged span.

    Classification priority:
    1. US (Unspliced): ≥50bp coverage within 985-5212
    2. MS (Multiply Spliced): coverage before 6597, gap 6597-8806, coverage after 8806
    3. SS (Singly Spliced): coverage before 985, gap 985-5212, coverage 6597-8806
    4. any: Everything else

    Note: This is a simplified classification based on span boundaries only.
    For precise classification from CIGAR-parsed aligned blocks, use
    bamsort_isoform_longread.py on the original BAM file.

    Args:
        merged_start: Start of merged molecule span
        merged_end: End of merged molecule span
        tolerance: +/- bp tolerance for boundaries

    Returns:
        str: Isoform classification (US, MS, SS, or any)
    """
    if merged_start is None or merged_end is None:
        return "any"

    span_length = merged_end - merged_start

    # Check US: Significant coverage within 985-5212 region
    if has_coverage_in_region(merged_start, merged_end,
                              REGION_US_START, REGION_US_END,
                              min_length=MIN_ALIGNMENT_LENGTH,
                              tolerance=tolerance):
        return "US"

    # Check MS: coverage before 6597, gap 6597-8806, coverage after 8806
    has_before_ms = has_coverage_in_region(merged_start, merged_end,
                                           0, REGION_MS_GAP_START,
                                           min_length=1, tolerance=tolerance)
    has_gap_ms = has_gap_in_region(merged_start, merged_end,
                                   REGION_MS_GAP_START, REGION_MS_GAP_END,
                                   tolerance=tolerance)
    has_after_ms = has_coverage_in_region(merged_start, merged_end,
                                          REGION_MS_GAP_END, float('inf'),
                                          min_length=1, tolerance=tolerance)

    if has_before_ms and has_gap_ms and has_after_ms:
        return "MS"

    # Check SS: coverage before 985, gap 985-5212, coverage 6597-8806
    has_before_ss = has_coverage_in_region(merged_start, merged_end,
                                           0, REGION_SS_GAP_START,
                                           min_length=1, tolerance=tolerance)
    has_gap_ss = has_gap_in_region(merged_start, merged_end,
                                   REGION_SS_GAP_START, REGION_SS_GAP_END,
                                   tolerance=tolerance)
    has_after_ss = has_coverage_in_region(merged_start, merged_end,
                                          REGION_SS_AFTER_START, REGION_SS_AFTER_END,
                                          min_length=1, tolerance=tolerance)

    if has_before_ss and has_gap_ss and has_after_ss:
        return "SS"

    return "any"


def load_candidates(candidates_file):
    """
    Load candidate complementary pairs.

    Returns:
        dict: {cellID: [(umi1, umi2, data_dict), ...]}
    """
    cell_pairs = defaultdict(list)

    with open(candidates_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            cell_id = row['cellID']
            pair_data = {
                'umi1': row['UMI_1'],
                'umi2': row['UMI_2'],
                'chrom': row['chrom'],
                'start_1': int(row['start_1']),
                'end_1': int(row['end_1']),
                'start_2': int(row['start_2']),
                'end_2': int(row['end_2']),
                'aligned_len_1': int(row['aligned_len_1']),
                'aligned_len_2': int(row['aligned_len_2']),
                'overlap_bp': int(row['overlap_bp']),
                'gap_bp': int(row['gap_bp']),
                'combined_span': int(row['combined_span']),
                'span_gain': int(row['span_gain']),
                'strand': row['strand']
            }
            cell_pairs[cell_id].append(pair_data)

    return dict(cell_pairs)


def load_all_umis(pairs_file):
    """
    Load all UMI information from pairs file to get complete UMI data.

    Returns:
        dict: {cellID: {umi: {'start': int, 'end': int, 'aligned_len': int, 'strand': str, 'chrom': str}}}
    """
    cell_umis = defaultdict(dict)

    with open(pairs_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            cell_id = row['cellID']
            chrom = row['chrom']

            # Extract UMI_1 data
            umi1 = row['UMI_1']
            if umi1 not in cell_umis[cell_id]:
                cell_umis[cell_id][umi1] = {
                    'start': int(row['start_1']),
                    'end': int(row['end_1']),
                    'aligned_len': int(row['aligned_len_1']),
                    'strand': row['strand_1'],
                    'chrom': chrom
                }

            # Extract UMI_2 data
            umi2 = row['UMI_2']
            if umi2 not in cell_umis[cell_id]:
                cell_umis[cell_id][umi2] = {
                    'start': int(row['start_2']),
                    'end': int(row['end_2']),
                    'aligned_len': int(row['aligned_len_2']),
                    'strand': row['strand_2'],
                    'chrom': chrom
                }

    return dict(cell_umis)


def build_molecule_components(cell_pairs, cell_umis, tolerance=5):
    """
    Build connected components representing reconstructed molecules.

    For each cell:
    1. Create a graph where UMIs are nodes
    2. Add edges for COMPLEMENTARY pairs
    3. Find connected components
    4. Each component = one reconstructed molecule
    5. Classify isoform based on merged span

    Args:
        cell_pairs: Dict of candidate pairs by cell
        cell_umis: Dict of UMI data by cell
        tolerance: Tolerance for isoform classification boundaries

    Returns:
        dict: {cellID: {molecule_id: {'umis': [...], 'merged_start': int, 'merged_end': int, 'isoform': str, ...}}}
    """
    molecules_by_cell = {}

    for cell_id, umis in cell_umis.items():
        uf = UnionFind()

        # Initialize all UMIs as separate sets
        for umi in umis:
            uf.make_set(umi)

        # Union complementary pairs
        if cell_id in cell_pairs:
            for pair in cell_pairs[cell_id]:
                uf.union(pair['umi1'], pair['umi2'])

        # Get components
        components = uf.get_components()

        # Build molecule data for each component
        molecules = {}
        for mol_idx, (root, member_umis) in enumerate(components.items()):
            # Sort UMIs by genomic position for chain ordering
            sorted_umis = sorted(member_umis, key=lambda u: umis[u]['start'])

            # Compute merged span
            starts = [umis[u]['start'] for u in member_umis]
            ends = [umis[u]['end'] for u in member_umis]
            total_aligned = sum(umis[u]['aligned_len'] for u in member_umis)

            merged_start = min(starts)
            merged_end = max(ends)
            merged_span = merged_end - merged_start

            # Determine strand (should be consistent)
            strands = set(umis[u]['strand'] for u in member_umis)
            strand = strands.pop() if len(strands) == 1 else 'mixed'

            # Determine chromosome (should be consistent for stitched molecules)
            chroms = set(umis[u]['chrom'] for u in member_umis)
            chrom = chroms.pop() if len(chroms) == 1 else 'mixed'

            # Classify isoform based on merged span
            isoform = classify_isoform_from_span(merged_start, merged_end, tolerance)

            molecules[f"mol_{mol_idx}"] = {
                'umis': sorted_umis,
                'num_fragments': len(member_umis),
                'chrom': chrom,
                'merged_start': merged_start,
                'merged_end': merged_end,
                'merged_span': merged_span,
                'total_aligned_bp': total_aligned,
                'strand': strand,
                'is_stitched': len(member_umis) > 1,
                'isoform': isoform
            }

        molecules_by_cell[cell_id] = molecules

    return molecules_by_cell


def compute_chain_details(molecules_by_cell, cell_umis):
    """
    Compute detailed chain information showing how fragments connect.

    For stitched molecules, this shows the order of fragments and
    the gaps/overlaps between consecutive fragments.

    Returns:
        list: Chain records for output
    """
    chains = []

    for cell_id, molecules in molecules_by_cell.items():
        umis = cell_umis[cell_id]

        for mol_id, mol_data in molecules.items():
            if mol_data['num_fragments'] == 1:
                # Single fragment - no chain to report
                umi = mol_data['umis'][0]
                chains.append({
                    'cellID': cell_id,
                    'molecule_id': mol_id,
                    'fragment_order': 1,
                    'umi': umi,
                    'start': umis[umi]['start'],
                    'end': umis[umi]['end'],
                    'aligned_len': umis[umi]['aligned_len'],
                    'gap_to_next': None,
                    'overlap_with_next': None
                })
            else:
                # Multiple fragments - compute chain
                sorted_umis = mol_data['umis']

                for i, umi in enumerate(sorted_umis):
                    gap_to_next = None
                    overlap_with_next = None

                    if i < len(sorted_umis) - 1:
                        next_umi = sorted_umis[i + 1]
                        current_end = umis[umi]['end']
                        next_start = umis[next_umi]['start']

                        if current_end <= next_start:
                            gap_to_next = next_start - current_end
                        else:
                            overlap_with_next = current_end - next_start

                    chains.append({
                        'cellID': cell_id,
                        'molecule_id': mol_id,
                        'fragment_order': i + 1,
                        'umi': umi,
                        'start': umis[umi]['start'],
                        'end': umis[umi]['end'],
                        'aligned_len': umis[umi]['aligned_len'],
                        'gap_to_next': gap_to_next,
                        'overlap_with_next': overlap_with_next
                    })

    return chains


def compute_statistics(molecules_by_cell):
    """Compute stitching and isoform statistics."""
    stats = {
        'total_cells': len(molecules_by_cell),
        'total_molecules': 0,
        'single_fragment_molecules': 0,
        'stitched_molecules': 0,
        'max_fragments_per_molecule': 0,
        'total_umis_stitched': 0,
        'cells_with_stitching': 0,
        'isoform_counts': defaultdict(int),
        'stitched_isoform_counts': defaultdict(int),
        'single_isoform_counts': defaultdict(int)
    }

    fragment_counts = []

    for cell_id, molecules in molecules_by_cell.items():
        cell_has_stitching = False

        for mol_id, mol_data in molecules.items():
            stats['total_molecules'] += 1
            num_frags = mol_data['num_fragments']
            fragment_counts.append(num_frags)

            # Track isoform
            isoform = mol_data['isoform']
            stats['isoform_counts'][isoform] += 1

            if num_frags == 1:
                stats['single_fragment_molecules'] += 1
                stats['single_isoform_counts'][isoform] += 1
            else:
                stats['stitched_molecules'] += 1
                stats['total_umis_stitched'] += num_frags
                stats['stitched_isoform_counts'][isoform] += 1
                cell_has_stitching = True

                if num_frags > stats['max_fragments_per_molecule']:
                    stats['max_fragments_per_molecule'] = num_frags

        if cell_has_stitching:
            stats['cells_with_stitching'] += 1

    # Convert defaultdicts to regular dicts for output
    stats['isoform_counts'] = dict(stats['isoform_counts'])
    stats['stitched_isoform_counts'] = dict(stats['stitched_isoform_counts'])
    stats['single_isoform_counts'] = dict(stats['single_isoform_counts'])

    # Compute distribution
    if fragment_counts:
        stats['mean_fragments_per_molecule'] = sum(fragment_counts) / len(fragment_counts)
        stats['fragment_distribution'] = dict(
            sorted([(k, fragment_counts.count(k)) for k in set(fragment_counts)])
        )

    return stats


def write_outputs(molecules_by_cell, chains, stats, output_prefix):
    """Write all output files."""

    # Molecules file - now includes isoform column
    molecules_file = f"{output_prefix}_molecules.csv"
    with open(molecules_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'cellID', 'molecule_id', 'num_fragments', 'umis',
            'chrom', 'merged_start', 'merged_end', 'merged_span',
            'total_aligned_bp', 'strand', 'is_stitched', 'isoform'
        ])

        for cell_id in sorted(molecules_by_cell.keys()):
            for mol_id, mol_data in molecules_by_cell[cell_id].items():
                writer.writerow([
                    cell_id, mol_id, mol_data['num_fragments'],
                    ';'.join(mol_data['umis']),
                    mol_data['chrom'],
                    mol_data['merged_start'], mol_data['merged_end'],
                    mol_data['merged_span'], mol_data['total_aligned_bp'],
                    mol_data['strand'], mol_data['is_stitched'],
                    mol_data['isoform']
                ])

    # Chains file
    chains_file = f"{output_prefix}_molecule_chains.csv"
    with open(chains_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'cellID', 'molecule_id', 'fragment_order', 'umi',
            'start', 'end', 'aligned_len', 'gap_to_next', 'overlap_with_next'
        ])

        for chain in chains:
            writer.writerow([
                chain['cellID'], chain['molecule_id'], chain['fragment_order'],
                chain['umi'], chain['start'], chain['end'], chain['aligned_len'],
                chain['gap_to_next'] if chain['gap_to_next'] is not None else '',
                chain['overlap_with_next'] if chain['overlap_with_next'] is not None else ''
            ])

    # Stats file - now includes isoform statistics
    stats_file = f"{output_prefix}_stitching_stats.csv"
    with open(stats_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['metric', 'value'])

        # Basic stats
        for key in ['total_cells', 'total_molecules', 'single_fragment_molecules',
                    'stitched_molecules', 'cells_with_stitching',
                    'max_fragments_per_molecule', 'total_umis_stitched']:
            writer.writerow([key, stats[key]])

        if 'mean_fragments_per_molecule' in stats:
            writer.writerow(['mean_fragments_per_molecule',
                             f"{stats['mean_fragments_per_molecule']:.2f}"])

        # Fragment distribution
        if 'fragment_distribution' in stats:
            writer.writerow(['', ''])
            writer.writerow(['fragments_per_molecule', 'count'])
            for k, v in stats['fragment_distribution'].items():
                writer.writerow([k, v])

        # Isoform distribution (all molecules)
        writer.writerow(['', ''])
        writer.writerow(['isoform_distribution_all', ''])
        writer.writerow(['isoform', 'count'])
        for iso in ['US', 'MS', 'SS', 'any']:
            count = stats['isoform_counts'].get(iso, 0)
            writer.writerow([iso, count])

        # Isoform distribution (stitched only)
        writer.writerow(['', ''])
        writer.writerow(['isoform_distribution_stitched', ''])
        writer.writerow(['isoform', 'count'])
        for iso in ['US', 'MS', 'SS', 'any']:
            count = stats['stitched_isoform_counts'].get(iso, 0)
            writer.writerow([iso, count])

        # Isoform distribution (single fragment only)
        writer.writerow(['', ''])
        writer.writerow(['isoform_distribution_single', ''])
        writer.writerow(['isoform', 'count'])
        for iso in ['US', 'MS', 'SS', 'any']:
            count = stats['single_isoform_counts'].get(iso, 0)
            writer.writerow([iso, count])

    return molecules_file, chains_file, stats_file


def main():
    parser = argparse.ArgumentParser(
        description="Stitch UMI fragments into reconstructed molecules with isoform classification",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Algorithm:
  1. Load COMPLEMENTARY pairs from candidates file
  2. Build graph: UMIs = nodes, COMPLEMENTARY relationships = edges
  3. Find connected components (each = one original molecule)
  4. Order fragments within each molecule by genomic position
  5. Classify isoform based on merged span
  6. Output reconstructed molecules with isoform assignments

Isoform Classification (MAC239 coordinates):
  US (Unspliced):       >=50bp coverage within 985-5212
  MS (Multiply Spliced): Coverage before 6597, gap 6597-8806, coverage after 8806
  SS (Singly Spliced):   Coverage before 985, gap 985-5212, coverage 6597-8806
  any:                   All other patterns

Output Files:
  <prefix>_molecules.csv:       Reconstructed molecule assignments with isoform
  <prefix>_molecule_chains.csv: Ordered fragment chains
  <prefix>_stitching_stats.csv: Summary statistics including isoform distribution

Example:
  python3 bamsort_longread_stitch.py results_candidates.csv results_pairs.csv stitched
  python3 bamsort_longread_stitch.py results_candidates.csv results_pairs.csv stitched --tolerance 10
        """
    )

    parser.add_argument("candidates_file", help="Candidates CSV from bamsort_longread_fragments.py")
    parser.add_argument("pairs_file", help="Pairs CSV from bamsort_longread_fragments.py")
    parser.add_argument("output_prefix", help="Prefix for output files")
    parser.add_argument("--tolerance", type=int, default=5,
                        help="Coordinate tolerance in bp for isoform classification (default: 5)")

    args = parser.parse_args()

    print("UMI Fragment Stitching with Isoform Classification")
    print("=" * 60)
    print(f"Candidates file: {args.candidates_file}")
    print(f"Pairs file:      {args.pairs_file}")
    print(f"Output prefix:   {args.output_prefix}")
    print(f"Tolerance:       +/- {args.tolerance}bp")
    print("=" * 60)

    try:
        print("\nLoading candidate pairs...")
        cell_pairs = load_candidates(args.candidates_file)
        print(f"  Loaded pairs from {len(cell_pairs)} cells")

        print("\nLoading UMI data...")
        cell_umis = load_all_umis(args.pairs_file)
        print(f"  Loaded UMI data from {len(cell_umis)} cells")

        print("\nBuilding molecule components and classifying isoforms...")
        molecules_by_cell = build_molecule_components(cell_pairs, cell_umis,
                                                      tolerance=args.tolerance)

        print("\nComputing chain details...")
        chains = compute_chain_details(molecules_by_cell, cell_umis)

        print("\nComputing statistics...")
        stats = compute_statistics(molecules_by_cell)

        print("\nWriting output files...")
        mol_f, chain_f, stats_f = write_outputs(
            molecules_by_cell, chains, stats, args.output_prefix
        )

        print("\n" + "=" * 60)
        print("STITCHING SUMMARY")
        print("=" * 60)
        print(f"Total cells analyzed:          {stats['total_cells']}")
        print(f"Total reconstructed molecules: {stats['total_molecules']}")
        print(f"Single-fragment molecules:     {stats['single_fragment_molecules']}")
        print(f"Stitched molecules:            {stats['stitched_molecules']}")
        print(f"Cells with stitching:          {stats['cells_with_stitching']}")
        print(f"Max fragments per molecule:    {stats['max_fragments_per_molecule']}")
        print(f"Total UMIs stitched:           {stats['total_umis_stitched']}")

        if 'fragment_distribution' in stats:
            print("\nFragment distribution:")
            for k, v in sorted(stats['fragment_distribution'].items()):
                print(f"  {k} fragment(s): {v} molecules")

        print("\n" + "-" * 60)
        print("ISOFORM DISTRIBUTION (all molecules)")
        print("-" * 60)
        total_mol = stats['total_molecules']
        for iso in ['US', 'MS', 'SS', 'any']:
            count = stats['isoform_counts'].get(iso, 0)
            pct = 100 * count / max(1, total_mol)
            print(f"  {iso}: {count} ({pct:.1f}%)")

        if stats['stitched_molecules'] > 0:
            print("\n" + "-" * 60)
            print("ISOFORM DISTRIBUTION (stitched molecules only)")
            print("-" * 60)
            total_stitched = stats['stitched_molecules']
            for iso in ['US', 'MS', 'SS', 'any']:
                count = stats['stitched_isoform_counts'].get(iso, 0)
                pct = 100 * count / max(1, total_stitched)
                print(f"  {iso}: {count} ({pct:.1f}%)")

        print("=" * 60)
        print(f"\nOutput files:")
        print(f"  Molecules: {mol_f}")
        print(f"  Chains:    {chain_f}")
        print(f"  Stats:     {stats_f}")

    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()