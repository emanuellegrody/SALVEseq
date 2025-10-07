#!/usr/bin/env python3
"""
Simplified XF flag calculator for STAR BAM files
"""

import pysam
import argparse
import sys
from collections import defaultdict
import logging

# XF Flag bit definitions
XF_BITS = {
    'CONF_MAPPED': 1,        # Confidently mapped to transcriptome
    'LOW_SUPPORT_UMI': 2,    # Discarded due to lower read support
    'GENE_DISCORDANT': 4,    # Maps to multiple genes
    'UMI_COUNT': 8,          # Representative molecule for UMI counting
    'CONF_FEATURE': 16,      # Confidently assigned feature
    'FILTERED_TARGET_UMI': 32 # Removed by targeted UMI filtering -- NOT USED
}

def setup_logging(verbose=False):
    """Setup logging configuration"""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=level
    )
    return logging.getLogger(__name__)

def is_confidently_mapped(read, min_mapq=255):
    """Check if read is confidently mapped"""
    if read.is_unmapped:
        return False
    if read.mapping_quality < min_mapq:
        return False
    # Check for unique mapping (NH:i:1 from STAR)
    try:
        if read.get_tag('NH') != 1:
            return False
    except KeyError:
        pass
    return True

def has_multiple_genes(read, gene_tag='GX'):
    """Check if read maps to multiple genes"""
    try:
        gene = read.get_tag(gene_tag)
        return ';' in gene or ',' in gene
    except KeyError:
        return False

def has_valid_feature_assignment(read, gene_tag='GX'):
    """Check if read has valid gene assignment"""
    try:
        gene = read.get_tag(gene_tag)
        return gene and gene.strip() != ''
    except KeyError:
        return False

def get_read_key(read, barcode_tag='CB', umi_tag='UB', gene_tag='GX'):
    """Get unique key for UMI grouping"""
    try:
        barcode = read.get_tag(barcode_tag)
        umi = read.get_tag(umi_tag)
        gene = read.get_tag(gene_tag)
        # Take first gene if multiple
        gene = gene.split(';')[0].split(',')[0]
        return (barcode, umi, gene)
    except KeyError:
        return None

def find_umi_representatives(bam_file, barcode_tag='CB', umi_tag='UB', gene_tag='GX', min_mapq=255):
    """
    Single pass to find UMI representatives
    Returns dict mapping (barcode, umi, gene) -> representative_read_name
    """
    logger = logging.getLogger(__name__)
    logger.info("Finding UMI representatives...")
    
    umi_groups = defaultdict(list)
    
    with pysam.AlignmentFile(bam_file, 'rb') as bam:
        for read in bam:
            # Only consider confidently mapped reads with valid assignments
            if not is_confidently_mapped(read, min_mapq):
                continue
            if not has_valid_feature_assignment(read, gene_tag):
                continue
                
            key = get_read_key(read, barcode_tag, umi_tag, gene_tag)
            if key is None:
                continue
                
            # Store minimal read info for representative selection
            read_info = {
                'name': read.query_name,
                'mapq': read.mapping_quality,
                'alignment_score': read.get_tag('AS') if read.has_tag('AS') else 0
            }
            umi_groups[key].append(read_info)
    
    # Select representatives (highest MAPQ, then alignment score)
    representatives = {}
    for key, reads in umi_groups.items():
        if reads:
            rep = max(reads, key=lambda x: (x['mapq'], x['alignment_score']))
            representatives[key] = rep['name']
    
    logger.info(f"Found {len(representatives)} UMI representatives")
    return representatives

def calculate_xf_flag(read, representatives, barcode_tag='CB', umi_tag='UB', gene_tag='GX', min_mapq=255):
    """Calculate XF flag for a read"""
    xf_flag = 0
    
    # Check basic mapping quality
    is_conf_mapped = is_confidently_mapped(read, min_mapq)
    has_valid_feature = has_valid_feature_assignment(read, gene_tag)
    maps_multiple_genes = has_multiple_genes(read, gene_tag)
    
    # Set CONF_MAPPED bit
    if is_conf_mapped:
        xf_flag |= XF_BITS['CONF_MAPPED']
    
    # Set CONF_FEATURE bit  
    if has_valid_feature:
        xf_flag |= XF_BITS['CONF_FEATURE']
    
    # Set GENE_DISCORDANT bit if maps to multiple genes
    if maps_multiple_genes:
        xf_flag |= XF_BITS['GENE_DISCORDANT']
    
    # Check if this is a UMI representative (only for single-gene assignments)
    if is_conf_mapped and has_valid_feature and not maps_multiple_genes:
        key = get_read_key(read, barcode_tag, umi_tag, gene_tag)
        if key and key in representatives and representatives[key] == read.query_name:
            xf_flag |= XF_BITS['UMI_COUNT']
    
    return xf_flag

def add_xf_flags_to_bam(input_bam, output_bam, barcode_tag='CB', umi_tag='UB', gene_tag='GX', min_mapq=255):
    """
    Add XF flags to BAM file in two passes:
    1. Find UMI representatives
    2. Add XF flags to all reads
    """
    logger = logging.getLogger(__name__)
    
    # First pass: find representatives
    representatives = find_umi_representatives(input_bam, barcode_tag, umi_tag, gene_tag, min_mapq)
    
    # Second pass: add XF flags
    logger.info("Adding XF flags to reads...")
    stats = defaultdict(int)
    
    with pysam.AlignmentFile(input_bam, 'rb') as input_file:
        with pysam.AlignmentFile(output_bam, 'wb', template=input_file) as output_file:
            for read in input_file:
                # Calculate XF flag
                xf_flag = calculate_xf_flag(read, representatives, barcode_tag, umi_tag, gene_tag, min_mapq)
                
                # Add XF tag
                read.set_tag('xf', xf_flag, 'i')
                
                # Write read
                output_file.write(read)
                
                # Update stats
                stats[f'xf_{xf_flag}'] += 1
                stats['total'] += 1
    
    # Print summary stats
    logger.info(f"Processed {stats['total']:,} reads")
    
    xf_counts = [(k, v) for k, v in stats.items() if k.startswith('xf_')]
    xf_counts.sort(key=lambda x: int(x[0].split('_')[1]))
    
    logger.info("XF flag distribution:")
    for xf_key, count in xf_counts:
        xf_value = int(xf_key.split('_')[1])
        percentage = 100.0 * count / stats['total']
        description = get_xf_description(xf_value)
        logger.info(f"  XF={xf_value:2d}: {count:8,} ({percentage:5.1f}%) - {description}")

def get_xf_description(xf_value):
    """Get human-readable description of XF flag"""
    descriptions = {
        0: "Unmapped/low quality/missing tags",
        1: "Mapped only", 
        4: "Multiple genes only",
        5: "Mapped + Multiple genes",
        16: "Feature assigned only",
        17: "Mapped + Feature (duplicate)",
        20: "Feature + Multiple genes",
        21: "Mapped + Feature + Multiple genes", 
        25: "Mapped + Feature + Counted (UMI)",
    }
    return descriptions.get(xf_value, f"Other ({bin(xf_value)})")

def main():
    parser = argparse.ArgumentParser(
        description="Add XF flags to STAR BAM files (simplified version)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  python add_xf_flags_simple.py input.bam output.bam
  
  # With custom tags
  python add_xf_flags_simple.py input.bam output.bam --barcode-tag CR --umi-tag UR --gene-tag GN
  
  # Lower MAPQ threshold
  python add_xf_flags_simple.py input.bam output.bam --min-mapq 10
        """
    )
    
    parser.add_argument('input_bam', help='Input BAM file')
    parser.add_argument('output_bam', help='Output BAM file with XF flags')
    
    parser.add_argument('--barcode-tag', default='CB', 
                       help='BAM tag for cell barcode (default: CB)')
    parser.add_argument('--umi-tag', default='UB',
                       help='BAM tag for UMI (default: UB)')
    parser.add_argument('--gene-tag', default='GX',
                       help='BAM tag for gene assignment (default: GX)')
    parser.add_argument('--min-mapq', type=int, default=255,
                       help='Minimum mapping quality (default: 255)')
    parser.add_argument('--verbose', '-v', action='store_true',
                       help='Verbose logging')
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(args.verbose)
    
    # Validate inputs
    try:
        pysam.AlignmentFile(args.input_bam, 'rb').close()
    except:
        logger.error(f"Cannot read input BAM: {args.input_bam}")
        sys.exit(1)
    
    try:
        # Process BAM file
        logger.info(f"Processing: {args.input_bam} -> {args.output_bam}")
        add_xf_flags_to_bam(
            args.input_bam, 
            args.output_bam,
            args.barcode_tag,
            args.umi_tag, 
            args.gene_tag,
            args.min_mapq
        )
        
        # Index output BAM
        logger.info("Indexing output BAM...")
        pysam.index(args.output_bam)
        
        logger.info("Complete!")
        
    except Exception as e:
        logger.error(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
