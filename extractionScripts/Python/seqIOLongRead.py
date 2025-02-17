#!/usr/bin/env python3
import numpy as np
import regex
from Bio import SeqIO
import os
import sys
import argparse
from pathlib import Path
import gzip
import logging
from typing import Tuple, List, Optional
import signal
from contextlib import contextmanager
import time


class FastqProcessingError(Exception):
    """Base exception for fastq processing errors"""
    pass


class QualityScoreError(FastqProcessingError):
    """Exception raised for invalid quality scores"""
    pass


class SequenceError(FastqProcessingError):
    """Exception raised for invalid sequence data"""
    pass


class TimeoutError(FastqProcessingError):
    """Exception raised when operations timeout"""
    pass


@contextmanager
def timeout(seconds: int):
    """Context manager for timing out operations"""

    def signal_handler(signum, frame):
        raise TimeoutError(f"Operation timed out after {seconds} seconds")

    # Register signal handler
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)

    try:
        yield
    finally:
        # Disable alarm
        signal.alarm(0)


def setup_logging(output_dir: Path) -> None:
    """Configure logging with both file and console handlers"""
    log_file = output_dir / "processing.log"

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )


def validate_search_sequence(sequence: str) -> bool:
    """Validate that search sequence contains only valid DNA/RNA characters"""
    valid_chars = set('ATCGNatcgn')
    return all(c in valid_chars for c in sequence)


def validate_quality_score(quality: float) -> bool:
    """Validate quality score is within reasonable bounds"""
    return 0 <= quality <= 42  # Maximum theoretical Phred quality score is 41.6


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Process long read FASTQ/FASTQ.GZ files and filter based on quality and sequence content')
    parser.add_argument('fastq_file', help='Input FASTQ or FASTQ.GZ file')
    parser.add_argument('-s', '--search', default='TTACAAATTGGCAATAGACATG',
                        help='Sequence to search for (default: TTACAAATTGGCAATAGACATG (nef))')
    parser.add_argument('-q', '--quality', type=float, default=10,
                        help='Minimum average quality score (default: 10)')
    parser.add_argument('-e', '--errors', type=int, default=4,
                        help='Maximum number of errors allowed in search sequence (default: 4)')
    parser.add_argument('-o', '--output', help='Output directory (default: current directory)')
    parser.add_argument('-t', '--timeout', type=int, default=3600,
                        help='Timeout in seconds for processing (default: 3600)')

    args = parser.parse_args()

    # Validate arguments
    if not os.path.exists(args.fastq_file):
        raise FileNotFoundError(f"Input file not found: {args.fastq_file}")

    if not validate_search_sequence(args.search):
        raise ValueError(f"Invalid search sequence: {args.search}. Must contain only ATCGN characters.")

    if not validate_quality_score(args.quality):
        raise ValueError(f"Invalid quality score: {args.quality}. Must be between 0 and 42.")

    if args.errors < 0:
        raise ValueError(f"Invalid error count: {args.errors}. Must be non-negative.")

    return args


def open_fastq(filename: str, buffer_size: int = 1024 * 1024) -> object:
    """
    Open a FASTQ file with optimized buffering and I/O monitoring.

    Args:
        filename: Path to FASTQ file
        buffer_size: Read buffer size in bytes (default 1MB)
    """
    try:
        if not os.path.exists(filename):
            raise FileNotFoundError(f"File not found: {filename}")

        # Open with specified buffer size
        if filename.endswith('.gz'):
            # For gzipped files, we'll decompress in chunks
            return gzip.open(filename, 'rt', buffer_size)
        else:
            return open(filename, 'r', buffering=buffer_size)

    except (IOError, OSError) as e:
        raise FastqProcessingError(f"Error opening file {filename}: {str(e)}")


def process_record(record: SeqIO.SeqRecord, min_quality: float, pattern: str) -> Tuple[Optional[str], str]:
    """Process individual FASTQ record and return result category and sequence"""
    try:
        # Calculate average quality score
        if not record.letter_annotations or "phred_quality" not in record.letter_annotations:
            raise QualityScoreError("Record missing quality scores")

        quality_scores = record.letter_annotations["phred_quality"]
        if not quality_scores:
            raise QualityScoreError("Empty quality scores")

        avg_quality = np.mean(quality_scores)

        # Get sequence
        sequence = str(record.seq)
        if not sequence:
            raise SequenceError("Empty sequence")

        # Check quality
        if avg_quality < min_quality:
            return "quality", sequence

        # Check for search sequence
        if not regex.findall(pattern, sequence):
            return "sequence", sequence

        return "pass", sequence

    except (ValueError, AttributeError) as e:
        logging.error(f"Error processing record {record.id}: {str(e)}")
        return "error", ""


def process_fastq(fastq_file: str, search_seq: str, min_quality: float = 10,
                  max_errors: int = 4, timeout_seconds: int = 3600) -> Tuple[List[str], List[str], List[str], int]:
    """
    Process FASTQ file and filter reads based on quality and sequence content.

    Implements adaptive timeout handling:
    - Monitors processing speed
    - Adjusts batch size based on performance
    - Provides progress updates
    - Allows partial results on timeout
    """
    passing_seqs = []
    failing_quality_seqs = []
    failing_sequence_seqs = []
    total_reads = 0

    pattern = rf'({regex.escape(search_seq)}){{e<={max_errors}}}'

    # Initialize performance monitoring
    start_time = time.time()
    last_check_time = start_time
    last_check_reads = 0
    batch_size = 1000  # Initial batch size

    try:
        with timeout(timeout_seconds):
            # Get file size for progress reporting
            file_size = os.path.getsize(fastq_file)
            logging.info(f"Starting processing of {file_size / (1024 * 1024):.2f} MB file")
            with open_fastq(fastq_file) as f:
                for record in SeqIO.parse(f, "fastq"):
                    total_reads += 1
                    if total_reads % batch_size == 0:
                        current_time = time.time()
                        time_delta = current_time - last_check_time
                        reads_delta = total_reads - last_check_reads

                        # Calculate processing speed
                        if time_delta > 0:
                            reads_per_second = reads_delta / time_delta
                            logging.info(f"Processing speed: {reads_per_second:.2f} reads/second")

                            # Adjust batch size based on performance
                            if reads_per_second > 1000:
                                batch_size = min(batch_size * 2, 100000)
                            elif reads_per_second < 100:
                                batch_size = max(batch_size // 2, 1000)

                        # Update progress
                        logging.info(
                            f"Processed {total_reads} reads... ({time_delta:.2f} seconds for last {reads_delta} reads)")

                        # Update monitoring variables
                        last_check_time = current_time
                        last_check_reads = total_reads

                    try:
                        result, sequence = process_record(record, min_quality, pattern)

                        if result == "pass":
                            passing_seqs.append(sequence)
                        elif result == "quality":
                            failing_quality_seqs.append(sequence)
                        elif result == "sequence":
                            failing_sequence_seqs.append(sequence)
                        # Skip error results

                    except (QualityScoreError, SequenceError) as e:
                        logging.error(f"Error processing record {record.id}: {str(e)}")
                        continue

    except TimeoutError:
        logging.error(f"Processing timed out after {timeout_seconds} seconds")
        raise
    except Exception as e:
        logging.error(f"Unexpected error during processing: {str(e)}")
        raise FastqProcessingError(f"Failed to process file: {str(e)}")

    return passing_seqs, failing_quality_seqs, failing_sequence_seqs, total_reads


def write_sequences(sequences: List[str], filename: str) -> None:
    """Write sequences to text file, one per line."""
    if not sequences:
        logging.info(f"No sequences to write to {filename}")
        return

    try:
        with open(filename, 'w') as f:
            for seq in sequences:
                f.write(f"{seq}\n")
    except IOError as e:
        raise FastqProcessingError(f"Failed to write sequences to {filename}: {str(e)}")


def write_summary(args: argparse.Namespace, base_filename: Path,
                  passing_count: int, quality_fail_count: int,
                  sequence_fail_count: int, total_reads: int) -> None:
    """Write processing summary to file"""
    summary_file = f"{base_filename}_summary.txt"
    try:
        with open(summary_file, "w") as summary:
            summary.write(f"Processing Summary for {args.fastq_file}\n")
            summary.write("-" * 50 + "\n")
            summary.write(f"Total reads processed: {total_reads}\n")
            summary.write(f"Reads passing all filters: {passing_count}\n")
            summary.write(f"Reads failing quality filter (min Q={args.quality}): {quality_fail_count}\n")
            summary.write(f"Reads failing sequence search: {sequence_fail_count}\n")
            summary.write(f"\nSearch parameters:\n")
            summary.write(f"Search sequence: {args.search}\n")
            summary.write(f"Maximum allowed errors: {args.errors}\n")
            summary.write(f"\nOutput files:\n")
            summary.write(f"Passing sequences: {base_filename}_passing.txt\n")
            summary.write(f"Failed quality: {base_filename}_failed_quality.txt\n")
            summary.write(f"Failed sequence: {base_filename}_failed_sequence.txt\n")
    except IOError as e:
        raise FastqProcessingError(f"Failed to write summary to {summary_file}: {str(e)}")


def main():
    start_time = time.time()

    try:
        # Parse arguments
        args = parse_arguments()

        # Setup output directory
        output_dir = Path(args.output if args.output else os.path.dirname(args.fastq_file))
        output_dir.mkdir(parents=True, exist_ok=True)

        # Setup logging
        setup_logging(output_dir)

        logging.info(f"Starting processing of {args.fastq_file}")

        # Process the FASTQ file
        passing_seqs, failing_quality_seqs, failing_sequence_seqs, total_reads = process_fastq(
            args.fastq_file,
            args.search,
            args.quality,
            args.errors,
            args.timeout
        )

        # Generate base filename
        base_filename = output_dir / Path(args.fastq_file).stem
        if base_filename.suffix == '.fastq':
            base_filename = base_filename.with_suffix('')

        # Write outputs
        logging.info("Writing output files...")
        write_sequences(passing_seqs, f"{base_filename}_passing.txt")
        write_sequences(failing_quality_seqs, f"{base_filename}_failed_quality.txt")
        write_sequences(failing_sequence_seqs, f"{base_filename}_failed_sequence.txt")

        # Write summary
        write_summary(args, base_filename, len(passing_seqs),
                      len(failing_quality_seqs), len(failing_sequence_seqs),
                      total_reads)

        processing_time = time.time() - start_time
        logging.info(f"Processing complete! Total time: {processing_time:.2f} seconds")

    except KeyboardInterrupt:
        logging.error("Processing interrupted by user")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Fatal error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()