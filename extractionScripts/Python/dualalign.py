import pysam
import pandas as pd

def extract_paired_alignments(bam_file):
    alignments = []

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam:
            if read.is_paired and not read.is_unmapped:
                alignment_info = {
                    'read_name': read.query_name,
                    'read1_chr': read.reference_name if read.is_read1 else None,
                    'read1_pos': read.reference_start if read.is_read1 else None,
                    'read1_cigar': read.cigarstring if read.is_read1 else None,
                    'read2_chr': read.next_reference_name if read.is_read1 else read.reference_name,
                    'read2_pos': read.next_reference_start if read.is_read1 else read.reference_start,
                    'read2_cigar': None if read.is_read1 else read.cigarstring,
                    'mapping_quality': read.mapping_quality,
                    'is_proper_pair': read.is_proper_pair
                }
                alignments.append(alignment_info)

    return pd.DataFrame(alignments)

# Usage
df = extract_paired_alignments('/Volumes/fsmresfiles/Basic_Sciences/CDB/GoyalLab/People/EmmieGrody/Data/EGS/extractedData/EGS004/300cy/counts/Mmul_10_mac239annot_JK85_D13/outs/possorted_genome_bam.bam')