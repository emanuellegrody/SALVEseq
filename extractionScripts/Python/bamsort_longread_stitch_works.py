#!/usr/bin/env python3

"""
UMI Fragment Stitching - Step 2

Description:
    Takes the output from umi_fragment_analysis.py and performs the actual
    reconstruction of original molecules by merging complementary UMI fragments.

Algorithm Overview:
    1. Load candidate pairs from _candidates.csv
    2. Build a graph where UMIs are nodes and COMPLEMENTARY relationships are edges
    3. Find connected components - each component represents one original molecule
    4. For each component, compute the merged alignment span
    5. Optionally, output a "pseudo-BAM" or annotation file for visualization

Output Files:
    1. _molecules.csv: Reconstructed molecules with component assignments
    2. _molecule_chains.csv: Ordered fragment chains within each molecule
    3. _stitching_stats.csv: Statistics on stitching success

Usage:
    python3 bamsort_longread_stitch.py candidates.csv pairs.csv output_prefix
"""

import csv
import sys
import argparse
from collections import defaultdict


class UnionFind:
    """
    Union-Find data structure for efficient connected component detection.

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
        dict: {cellID: {umi: {'start': int, 'end': int, 'aligned_len': int, 'strand': str}}}
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


def build_molecule_components(cell_pairs, cell_umis):
    """
    Build connected components representing reconstructed molecules.

    For each cell:
    1. Create a graph where UMIs are nodes
    2. Add edges for COMPLEMENTARY pairs
    3. Find connected components
    4. Each component = one reconstructed molecule

    Returns:
        dict: {cellID: {molecule_id: {'umis': [...], 'merged_start': int, 'merged_end': int, ...}}}
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

            molecules[f"mol_{mol_idx}"] = {
                'umis': sorted_umis,
                'num_fragments': len(member_umis),
                'chrom': chrom,
                'merged_start': merged_start,
                'merged_end': merged_end,
                'merged_span': merged_span,
                'total_aligned_bp': total_aligned,
                'strand': strand,
                'is_stitched': len(member_umis) > 1
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
    """Compute stitching statistics."""
    stats = {
        'total_cells': len(molecules_by_cell),
        'total_molecules': 0,
        'single_fragment_molecules': 0,
        'stitched_molecules': 0,
        'max_fragments_per_molecule': 0,
        'total_umis_stitched': 0,
        'cells_with_stitching': 0
    }
    
    fragment_counts = []
    span_gains = []
    
    for cell_id, molecules in molecules_by_cell.items():
        cell_has_stitching = False
        
        for mol_id, mol_data in molecules.items():
            stats['total_molecules'] += 1
            num_frags = mol_data['num_fragments']
            fragment_counts.append(num_frags)
            
            if num_frags == 1:
                stats['single_fragment_molecules'] += 1
            else:
                stats['stitched_molecules'] += 1
                stats['total_umis_stitched'] += num_frags
                cell_has_stitching = True
                
                if num_frags > stats['max_fragments_per_molecule']:
                    stats['max_fragments_per_molecule'] = num_frags
        
        if cell_has_stitching:
            stats['cells_with_stitching'] += 1
    
    # Compute distribution
    if fragment_counts:
        stats['mean_fragments_per_molecule'] = sum(fragment_counts) / len(fragment_counts)
        stats['fragment_distribution'] = dict(
            sorted([(k, fragment_counts.count(k)) for k in set(fragment_counts)])
        )
    
    return stats


def write_outputs(molecules_by_cell, chains, stats, output_prefix):
    """Write all output files."""
    
    # Molecules file
    molecules_file = f"{output_prefix}_molecules.csv"
    with open(molecules_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([
            'cellID', 'molecule_id', 'num_fragments', 'umis',
            'chrom', 'merged_start', 'merged_end', 'merged_span',
            'total_aligned_bp', 'strand', 'is_stitched'
        ])
        
        for cell_id in sorted(molecules_by_cell.keys()):
            for mol_id, mol_data in molecules_by_cell[cell_id].items():
                writer.writerow([
                    cell_id, mol_id, mol_data['num_fragments'],
                    ';'.join(mol_data['umis']),
                    mol_data['chrom'],
                    mol_data['merged_start'], mol_data['merged_end'],
                    mol_data['merged_span'], mol_data['total_aligned_bp'],
                    mol_data['strand'], mol_data['is_stitched']
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
    
    # Stats file
    stats_file = f"{output_prefix}_stitching_stats.csv"
    with open(stats_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['metric', 'value'])
        for key, value in stats.items():
            if key != 'fragment_distribution':
                writer.writerow([key, value])
        
        # Write fragment distribution
        if 'fragment_distribution' in stats:
            writer.writerow(['', ''])
            writer.writerow(['fragments_per_molecule', 'count'])
            for k, v in stats['fragment_distribution'].items():
                writer.writerow([k, v])
    
    return molecules_file, chains_file, stats_file


def main():
    parser = argparse.ArgumentParser(
        description="Stitch UMI fragments into reconstructed molecules",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Algorithm:
  1. Load COMPLEMENTARY pairs from candidates file
  2. Build graph: UMIs = nodes, COMPLEMENTARY relationships = edges
  3. Find connected components (each = one original molecule)
  4. Order fragments within each molecule by genomic position
  5. Output reconstructed molecules with chain information

Output Files:
  <prefix>_molecules.csv:       Reconstructed molecule assignments
  <prefix>_molecule_chains.csv: Ordered fragment chains
  <prefix>_stitching_stats.csv: Summary statistics

Interpretation:
  - 'is_stitched=True' indicates a molecule reconstructed from multiple UMIs
  - 'num_fragments' shows how many UMIs were merged
  - 'merged_span' is the reconstructed genomic coverage
  - Chain files show fragment order and gaps/overlaps between fragments

Example:
  python3 bamsort_longread_stitch.py results_candidates.csv results_pairs.csv stitched
        """
    )
    
    parser.add_argument("candidates_file", help="Candidates CSV from umi_fragment_analysis.py")
    parser.add_argument("pairs_file", help="Pairs CSV from umi_fragment_analysis.py")
    parser.add_argument("output_prefix", help="Prefix for output files")
    
    args = parser.parse_args()
    
    print("UMI Fragment Stitching")
    print("=" * 60)
    print(f"Candidates file: {args.candidates_file}")
    print(f"Pairs file:      {args.pairs_file}")
    print(f"Output prefix:   {args.output_prefix}")
    print("=" * 60)
    
    try:
        print("\nLoading candidate pairs...")
        cell_pairs = load_candidates(args.candidates_file)
        print(f"  Loaded pairs from {len(cell_pairs)} cells")
        
        print("\nLoading UMI data...")
        cell_umis = load_all_umis(args.pairs_file)
        print(f"  Loaded UMI data from {len(cell_umis)} cells")
        
        print("\nBuilding molecule components...")
        molecules_by_cell = build_molecule_components(cell_pairs, cell_umis)
        
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
