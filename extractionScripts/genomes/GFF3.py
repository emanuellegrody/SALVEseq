#!/usr/bin/env python3
"""
Convert rMATS SE.MATS.JC.txt output to the GFF3 format expected by
BRIE2's brie-count.

BRIE's parser (gtf_utils.load_genes) requires a strict hierarchy:
  gene -> mRNA (transcript) -> exon

For each skipping exon event, there are two isoforms:
  - Isoform 1 (inclusion):  upstream_exon -- skipped_exon -- downstream_exon
  - Isoform 2 (exclusion):  upstream_exon -- downstream_exon

rMATS uses 0-based start coordinates; GFF3 uses 1-based.
"""

import csv
import sys

def convert(input_file, output_file):
    with open(input_file, "r") as infile, open(output_file, "w") as out:
        reader = csv.DictReader(infile, delimiter="\t")
        for row in reader:
            event_id = row["ID"]
            chrom = row["chr"]
            strand = row["strand"]
            gene_id = row["GeneID"]

            # 0-based to 1-based conversion for start coordinates
            up_start = int(row["upstreamES"]) + 1
            up_end = int(row["upstreamEE"])
            se_start = int(row["exonStart_0base"]) + 1
            se_end = int(row["exonEnd"])
            dn_start = int(row["downstreamES"]) + 1
            dn_end = int(row["downstreamEE"])

            # Gene spans from earliest start to latest end
            gene_start = min(up_start, se_start, dn_start)
            gene_end = max(up_end, se_end, dn_end)

            # Unique IDs per event to avoid collisions when a gene
            # has multiple SE events
            eid = f"{gene_id}.{event_id}"

            # Gene line
            out.write(f"{chrom}\trMATS\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\t"
                      f"ID={eid};gene_id={eid};gene_name={gene_id}\n")

            # Isoform 1: inclusion (3 exons)
            inc_id = f"{eid}.in"
            inc_start = min(up_start, se_start, dn_start)
            inc_end = max(up_end, se_end, dn_end)
            out.write(f"{chrom}\trMATS\tmRNA\t{inc_start}\t{inc_end}\t.\t{strand}\t.\t"
                      f"ID={inc_id};Parent={eid}\n")
            out.write(f"{chrom}\trMATS\texon\t{up_start}\t{up_end}\t.\t{strand}\t.\t"
                      f"ID={inc_id}.1;Parent={inc_id}\n")
            out.write(f"{chrom}\trMATS\texon\t{se_start}\t{se_end}\t.\t{strand}\t.\t"
                      f"ID={inc_id}.2;Parent={inc_id}\n")
            out.write(f"{chrom}\trMATS\texon\t{dn_start}\t{dn_end}\t.\t{strand}\t.\t"
                      f"ID={inc_id}.3;Parent={inc_id}\n")

            # Isoform 2: exclusion (2 exons, skipped exon absent)
            exc_id = f"{eid}.ex"
            exc_start = min(up_start, dn_start)
            exc_end = max(up_end, dn_end)
            out.write(f"{chrom}\trMATS\tmRNA\t{exc_start}\t{exc_end}\t.\t{strand}\t.\t"
                      f"ID={exc_id};Parent={eid}\n")
            out.write(f"{chrom}\trMATS\texon\t{up_start}\t{up_end}\t.\t{strand}\t.\t"
                      f"ID={exc_id}.1;Parent={exc_id}\n")
            out.write(f"{chrom}\trMATS\texon\t{dn_start}\t{dn_end}\t.\t{strand}\t.\t"
                      f"ID={exc_id}.2;Parent={exc_id}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: python {sys.argv[0]} <SE.MATS.JC.txt> <output.gff3>")
        sys.exit(1)
    convert(sys.argv[1], sys.argv[2])
    print(f"Done. Wrote BRIE-compatible GFF3 to {sys.argv[2]}")