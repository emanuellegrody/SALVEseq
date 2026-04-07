"""
Extract Discordant Isoform Reads

Finds CB:UMI pairs where long read and short read isoform calls disagree,
then extracts matching reads from both BAM files for visual inspection in IGV.

Outputs:
- CSV listing all discordant CB:UMI pairs with both classifications
- Subset BAM from long read with only discordant reads
- Subset BAM from short read with only discordant reads

Usage:
    python3 extract_discordant_reads.py \
        --longread-csv tat_isoforms.csv \
        --shortread-csv GFP_invitro_UMI_read_counts_full.csv \
        --longread-bam tat.genome.dedup.bam \
        --shortread-bam possorted_genome_bam.bam \
        --output-dir discordant_reads/ \
        --sample tat \
        --longread-class MS --shortread-class US
"""

import pysam
import csv
import argparse
import os
from collections import defaultdict


def load_longread_isoforms(csv_path):
    """Load long read isoform CSV → dict of (cellID, UMI) → isoform."""
    result = {}
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = (row['cellID'], row['UMI'])
            result[key] = row['isoform']
    return result


def load_shortread_isoforms(csv_path):
    """
    Load short read UMI_read_counts_full.csv → dict of (cellID, UMI) → isoform.
    Strips -1 suffix from cellID. Consolidates categories to US/MS/SS/any.
    For duplicate CB:UMI, keeps most specific (US > MS > SS > any).
    """
    rank = {'US': 0, 'MS': 1, 'SS': 2}
    result = {}
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            cell = row['cellID'].rstrip('-1') if row['cellID'].endswith('-1') else row['cellID']
            umi = row['UMI']
            cat = row['category']

            # Map to simplified isoform
            if cat == 'US':
                iso = 'US'
            elif cat == 'MS':
                iso = 'MS'
            elif cat == 'SS':
                iso = 'SS'
            else:
                iso = 'any'

            key = (cell, umi)
            if key not in result or rank.get(iso, 3) < rank.get(result[key], 3):
                result[key] = iso
    return result


def find_discordant_pairs(longread_isos, shortread_isos,
                          longread_class=None, shortread_class=None):
    """
    Find CB:UMI pairs present in both with disagreeing isoform calls.
    Optionally filter to specific disagreement pattern.
    """
    shared_keys = set(longread_isos.keys()) & set(shortread_isos.keys())
    discordant = {}

    for key in shared_keys:
        lr = longread_isos[key]
        sr = shortread_isos[key]
        if lr == sr:
            continue
        if longread_class and lr != longread_class:
            continue
        if shortread_class and sr != shortread_class:
            continue
        discordant[key] = (lr, sr)

    return discordant


def extract_reads_from_bam(bam_path, cb_umi_set, output_bam_path,
                           reference="mac239", strip_suffix=False):
    """
    Extract reads matching (CB, UMI) pairs from a BAM file.

    Args:
        bam_path: Input BAM file
        cb_umi_set: Set of (cellID, UMI) tuples to extract
        output_bam_path: Output BAM file path
        reference: Reference to filter to (None for all)
        strip_suffix: If True, strip -1 from CB tag before matching
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    out = pysam.AlignmentFile(output_bam_path, "wb", header=bam.header)

    extracted = 0
    total = 0
    seen = set()

    if reference and reference in [ref['SN'] for ref in bam.header['SQ']]:
        reads_iter = bam.fetch(reference)
    else:
        reads_iter = bam.fetch()

    for read in reads_iter:
        total += 1
        if read.is_unmapped:
            continue
        if read.is_secondary or read.is_supplementary:
            continue

        try:
            cb = read.get_tag("CB")
        except KeyError:
            continue

        # Try UB first, then UR
        umi = None
        for tag in ("UB", "UR"):
            try:
                umi = read.get_tag(tag)
                break
            except KeyError:
                continue
        if umi is None:
            continue

        if strip_suffix:
            cb = cb.rstrip('-1') if cb.endswith('-1') else cb

        key = (cb, umi)
        if key in cb_umi_set and key not in seen:
            out.write(read)
            seen.add(key)
            extracted += 1

    out.close()
    bam.close()

    # Sort and index
    if extracted > 0:
        tmp_path = output_bam_path + ".unsorted.tmp"
        os.rename(output_bam_path, tmp_path)
        pysam.sort("-o", output_bam_path, tmp_path)
        os.remove(tmp_path)
        pysam.index(output_bam_path)

    return extracted, total


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract reads with discordant isoform calls between "
                    "long read and short read"
    )
    parser.add_argument("--longread-csv", required=True,
                        help="Long read isoform CSV (from bamsort_isoform_longread.py)")
    parser.add_argument("--shortread-csv", required=True,
                        help="Short read UMI_read_counts_full.csv")
    parser.add_argument("--longread-bam", required=True,
                        help="Long read BAM file")
    parser.add_argument("--shortread-bam", required=True,
                        help="Short read BAM file (CellRanger possorted)")
    parser.add_argument("--output-dir", required=True,
                        help="Output directory for subset BAMs and CSV")
    parser.add_argument("--sample", required=True,
                        help="Sample name (for output filenames)")
    parser.add_argument("--longread-class", default=None,
                        help="Filter to this long read class (e.g. MS)")
    parser.add_argument("--shortread-class", default=None,
                        help="Filter to this short read class (e.g. US)")
    parser.add_argument("--reference", default="mac239",
                        help="Reference name to filter BAM reads (default: mac239)")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # Load isoform calls
    print("Loading isoform classifications...")
    longread_isos = load_longread_isoforms(args.longread_csv)
    shortread_isos = load_shortread_isoforms(args.shortread_csv)
    print(f"  Long read: {len(longread_isos)} CB:UMI pairs")
    print(f"  Short read: {len(shortread_isos)} CB:UMI pairs")

    # Find discordant pairs
    discordant = find_discordant_pairs(
        longread_isos, shortread_isos,
        longread_class=args.longread_class,
        shortread_class=args.shortread_class
    )

    pattern = ""
    if args.longread_class or args.shortread_class:
        pattern = f" (LR={args.longread_class or 'any'} vs SR={args.shortread_class or 'any'})"
    print(f"  Discordant pairs{pattern}: {len(discordant)}")

    if len(discordant) == 0:
        print("No discordant pairs found.")
        exit(0)

    # Write discordant CSV
    disc_csv = os.path.join(args.output_dir,
                            f"{args.sample}_discordant.csv")
    with open(disc_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['cellID', 'UMI', 'longread_isoform', 'shortread_isoform'])
        for (cb, umi), (lr, sr) in sorted(discordant.items()):
            writer.writerow([cb, umi, lr, sr])
    print(f"  Discordant list: {disc_csv}")

    cb_umi_set = set(discordant.keys())

    # Extract from long read BAM
    print("Extracting from long read BAM...")
    longread_out = os.path.join(args.output_dir,
                                f"{args.sample}_discordant_longread.bam")
    longread_extracted, _ = extract_reads_from_bam(
        args.longread_bam, cb_umi_set, longread_out,
        reference=args.reference, strip_suffix=False
    )
    print(f"  Extracted {longread_extracted} long reads → {longread_out}")

    # Extract from short read BAM (CB tags have -1 suffix)
    print("Extracting from short read BAM...")
    shortread_out = os.path.join(args.output_dir,
                                 f"{args.sample}_discordant_shortread.bam")
    shortread_extracted, _ = extract_reads_from_bam(
        args.shortread_bam, cb_umi_set, shortread_out,
        reference=args.reference, strip_suffix=True
    )
    print(f"  Extracted {shortread_extracted} short reads → {shortread_out}")

    print("Done.")
