"""
rMATS Delta-PSI Extractor

Parses rMATS-turbo output files to extract PSI and delta-PSI values
for specific target exons. Searches across all five event types
(SE, A3SS, A5SS, MXE, RI) since the splicing pattern may vary by gene.

Coordinates are converted from 1-based (genome browser) to 0-based
(rMATS exonStart_0base) for matching.

Usage:
    python extract_deltaPSI.py --rmats-dir rmats_output/
"""

import os
import csv
import argparse

# Target exons: (gene_symbol, chr, start_1based, end)
# rMATS uses 0-based starts, so we subtract 1 from start when matching.
TARGETS = [
    ("Kras",   "chr6",  145170799, 145170923),
    ("Nf2",    "chr11", 4730577,   4730622),
    ("Slk",    "chr19", 47616086,  47616179),
    ("Lsm14b", "chr2",  179673586, 179673664),
    ("Dgkd",   "chr1",  87869090,  87869221),
]

EVENT_TYPES = ["SE", "A3SS", "A5SS", "MXE", "RI"]

# Columns that hold exon coordinates differ by event type.
# For SE: exonStart_0base, exonEnd (the target/skipped exon)
# For A3SS/A5SS: longExonStart_0base, longExonEnd, shortES, shortEE
# For MXE: 1stExonStart_0base, 1stExonEnd, 2ndExonStart_0base, 2ndExonEnd
# For RI: riExonStart_0base, riExonEnd
# We check all coordinate columns for each event type.
EXON_COORD_COLUMNS = {
    "SE":   [("exonStart_0base", "exonEnd")],
    "A3SS": [("longExonStart_0base", "longExonEnd"),
             ("shortES", "shortEE")],
    "A5SS": [("longExonStart_0base", "longExonEnd"),
             ("shortES", "shortEE")],
    "MXE":  [("1stExonStart_0base", "1stExonEnd"),
             ("2ndExonStart_0base", "2ndExonEnd")],
    "RI":   [("riExonStart_0base", "riExonEnd")],
}


def parse_psi_values(psi_string):
    """Parse comma-separated PSI values, skipping NA entries."""
    values = []
    for v in psi_string.split(","):
        v = v.strip()
        if v and v != "NA":
            values.append(float(v))
    return values


def mean(values):
    if not values:
        return float("nan")
    return sum(values) / len(values)


def search_event_file(filepath, event_type, targets_0based):
    """
    Search a single rMATS output file for target exons.

    Args:
        filepath: path to [EVENT].MATS.JCEC.txt
        event_type: one of SE, A3SS, A5SS, MXE, RI
        targets_0based: list of (gene, chr, start_0based, end)

    Returns:
        list of dicts with match info
    """
    results = []
    coord_pairs = EXON_COORD_COLUMNS[event_type]

    with open(filepath, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            gene_symbol = row.get("geneSymbol", "")
            chrom = row.get("chr", "")

            for target_gene, target_chr, target_start, target_end in targets_0based:
                if gene_symbol != target_gene:
                    continue
                if chrom != target_chr:
                    continue

                # Check all coordinate column pairs for this event type
                for start_col, end_col in coord_pairs:
                    if start_col not in row or end_col not in row:
                        continue
                    try:
                        exon_start = int(row[start_col])
                        exon_end = int(row[end_col])
                    except (ValueError, TypeError):
                        continue

                    if exon_start == target_start and exon_end == target_end:
                        psi_wt = parse_psi_values(row.get("IncLevel1", ""))
                        psi_ko = parse_psi_values(row.get("IncLevel2", ""))

                        results.append({
                            "gene": target_gene,
                            "event_type": event_type,
                            "chr": chrom,
                            "exon_start_0based": exon_start,
                            "exon_end": exon_end,
                            "coord_columns": "{},{}".format(start_col, end_col),
                            "IJC_WT": row.get("IJC_SAMPLE_1", ""),
                            "SJC_WT": row.get("SJC_SAMPLE_1", ""),
                            "IJC_KO": row.get("IJC_SAMPLE_2", ""),
                            "SJC_KO": row.get("SJC_SAMPLE_2", ""),
                            "PSI_WT_per_rep": row.get("IncLevel1", ""),
                            "PSI_KO_per_rep": row.get("IncLevel2", ""),
                            "mean_PSI_WT": mean(psi_wt),
                            "mean_PSI_KO": mean(psi_ko),
                            "delta_PSI": row.get("IncLevelDifference", ""),
                            "PValue": row.get("PValue", ""),
                            "FDR": row.get("FDR", ""),
                        })

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Extract delta-PSI for target exons from rMATS output"
    )
    parser.add_argument(
        "--rmats-dir",
        default="rmats_output/",
        help="Path to rMATS output directory (default: rmats_output/)"
    )
    parser.add_argument(
        "--output",
        default="target_exon_deltaPSI.tsv",
        help="Output TSV file (default: target_exon_deltaPSI.tsv)"
    )
    args = parser.parse_args()

    # Coordinates in TARGETS already match rMATS 0-based starts directly.
    targets_0based = [
        (gene, chrom, start, end)
        for gene, chrom, start, end in TARGETS
    ]

    all_results = []
    found_genes = set()

    for event_type in EVENT_TYPES:
        filepath = os.path.join(args.rmats_dir, "{}.MATS.JCEC.txt".format(event_type))
        if not os.path.exists(filepath):
            print("WARNING: {} not found, skipping".format(filepath))
            continue

        # Count total events in file
        with open(filepath) as f:
            n_events = sum(1 for _ in f) - 1
        print("Searching {} ({} events)".format(filepath, n_events))

        results = search_event_file(filepath, event_type, targets_0based)
        all_results.extend(results)
        for r in results:
            found_genes.add(r["gene"])

    # Report results
    print("\n" + "=" * 70)
    print("RESULTS: Delta-PSI (WT - KO) for target exons")
    print("  b1 = WT (SRR22206371, SRR22206372)")
    print("  b2 = KO (SRR22206369, SRR22206370)")
    print("=" * 70)

    if not all_results:
        print("\nNo target exons found in any event type file.")
        print("Possible causes:")
        print("  - Gene symbols may differ in the GTF (check capitalization)")
        print("  - Coordinates may not match (check genome build)")
        print("  - Exons may not show alternative splicing in this dataset")
        print("\nSearching for gene names regardless of coordinates:")
        for event_type in EVENT_TYPES:
            filepath = os.path.join(args.rmats_dir, "{}.MATS.JCEC.txt".format(event_type))
            if not os.path.exists(filepath):
                continue
            with open(filepath, "r") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    gene = row.get("geneSymbol", "")
                    for target_gene, _, _, _ in TARGETS:
                        if gene == target_gene:
                            print("  Found {} in {} at {}:{}-{}".format(
                                gene, event_type,
                                row.get("chr", ""),
                                row.get(list(row.keys())[5], ""),
                                row.get(list(row.keys())[6], "")
                            ))
    else:
        for r in all_results:
            print("\n{} [{}]  {}:{}-{}".format(
                r["gene"], r["event_type"], r["chr"],
                r["exon_start_0based"], r["exon_end"]
            ))
            print("  Inclusion junction counts WT: {}".format(r["IJC_WT"]))
            print("  Skipping junction counts WT:  {}".format(r["SJC_WT"]))
            print("  Inclusion junction counts KO: {}".format(r["IJC_KO"]))
            print("  Skipping junction counts KO:  {}".format(r["SJC_KO"]))
            print("  PSI WT (per replicate): {}".format(r["PSI_WT_per_rep"]))
            print("  PSI KO (per replicate): {}".format(r["PSI_KO_per_rep"]))
            print("  Mean PSI WT: {:.4f}".format(r["mean_PSI_WT"]))
            print("  Mean PSI KO: {:.4f}".format(r["mean_PSI_KO"]))
            print("  Delta PSI (WT - KO): {}".format(r["delta_PSI"]))
            print("  P-value: {}".format(r["PValue"]))
            print("  FDR: {}".format(r["FDR"]))
            try:
                fdr = float(r["FDR"])
                dpsi = abs(float(r["delta_PSI"]))
                sig = "YES" if fdr < 0.05 and dpsi > 0.1 else "NO"
            except (ValueError, TypeError):
                sig = "N/A"
            print("  Significant (FDR<0.05, |dPSI|>0.1): {}".format(sig))

        # Write TSV
        fieldnames = [
            "gene", "event_type", "chr", "exon_start_0based", "exon_end",
            "IJC_WT", "SJC_WT", "IJC_KO", "SJC_KO",
            "PSI_WT_per_rep", "PSI_KO_per_rep",
            "mean_PSI_WT", "mean_PSI_KO", "delta_PSI",
            "PValue", "FDR"
        ]
        with open(args.output, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                                    extrasaction="ignore")
            writer.writeheader()
            writer.writerows(all_results)
        print("\nResults saved to {}".format(args.output))

    # Check for missing genes
    missing = set(g for g, _, _, _ in TARGETS) - found_genes
    if missing:
        print("\nWARNING: No events found for: {}".format(", ".join(sorted(missing))))
        print("Check if these genes have alternative splicing detected in this dataset,")
        print("or verify gene symbol capitalization matches the GTF.")


if __name__ == "__main__":
    main()