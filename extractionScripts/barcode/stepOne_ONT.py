#!/usr/bin/env python3
#==================================================================================================
# stepOne from ONT nanopore barcode recovery pipeline
#
# Adapted from stepOne_EG.py (Illumina paired-end) for Oxford Nanopore long reads
# processed with 10X Genomics 3' chemistry.
#
# In ONT 10X reads, the cell barcode (16bp), UMI (12bp), and cDNA are all in a single
# read. The 10X Read1 adapter flanks the barcode/UMI region:
#   [10X adapter][cell barcode 16bp][UMI 12bp][polyT...][cDNA]
# The cDNA portion contains the SALVEseq construct: [stagger][primer][barcode]
#
# Reads may be in either orientation, so both forward and reverse complement are checked.
#
# Since there is no Illumina index, sample demultiplexing is done entirely by
# stagger+primer matching in the cDNA.
#
# Optional but encouraged barcode correction using --whitelist or --filtered-barcodes
#
# Usage:
#   python stepOne_ONT.py \
#       --fastq /path/to/reads.fastq.gz \
#       --staggers /path/to/staggers.txt \
#       --outdir /path/to/output/ \
#       --mode GFP /
#       --whitelist /path/to/whitelist.txt
#==================================================================================================

import argparse
import itertools
import math
import os
import sys
import regex
from gzip import open as gzopen

VALID_BASES = set("ACGTN")

# 10X Read1 adapter — flanks the cell barcode on the 5' side
TENX_ADAPTER = "CTACACGACGCTCTTCCGATCT"
TENX_ADAPTER_LEN = len(TENX_ADAPTER)  # 22bp

# 10X 3' v3 chemistry dimensions
CELL_BARCODE_LEN = 16
UMI_LEN = 12

# Minimum polyT length to confirm correct orientation
MIN_POLYT_LEN = 8

# Complement table for reverse complement (no BioPython dependency)
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq):
    """Return reverse complement of a sequence string."""
    return seq.translate(_COMP)[::-1]


# Pre-compute adapter k-mers for fast screening
KMER_K = 8

def build_kmer_set(seq, k=KMER_K):
    """Build a set of k-mers from a sequence for fast substring pre-screening."""
    return {seq[i:i+k] for i in range(len(seq) - k + 1)}

# Build kmer lists (ordered) for str.find()-based scanning
def build_kmer_list(seq, k=KMER_K):
    """Build an ordered list of unique k-mers for str.find()-based scanning."""
    seen = set()
    kmers = []
    for i in range(len(seq) - k + 1):
        km = seq[i:i+k]
        if km not in seen:
            seen.add(km)
            kmers.append(km)
    return kmers

ADAPTER_KMERS = build_kmer_list(TENX_ADAPTER)
ADAPTER_RC = reverse_complement(TENX_ADAPTER)
ADAPTER_RC_KMERS = build_kmer_list(ADAPTER_RC)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract SALVEseq barcodes from Oxford Nanopore 10X 3'v3 single-cell reads."
    )
    parser.add_argument("--fastq", required=True, help="Path to input FASTQ.gz (ONT reads)")
    parser.add_argument("--staggers", required=True,
                        help="Tab-separated file: sample_name<TAB>stagger_sequence")
    parser.add_argument("--outdir", required=True, help="Root output directory")
    parser.add_argument("--mode", choices=["GFP", "Vpx"], default="GFP",
                        help="Primer mode: GFP (default) or Vpx")
    parser.add_argument("--adapter-mismatches", type=int, default=3,
                        help="Max allowed mismatches when matching 10X adapter (default: 3)")
    parser.add_argument("--primer-mismatches", type=int, default=4,
                        help="Max allowed mismatches when matching primer (default: 4)")
    parser.add_argument("--whitelist", default=None,
                        help="Path to 10X barcode whitelist (e.g. 3M-february-2018.txt.gz). "
                             "Enables cell barcode error correction against the full whitelist.")
    parser.add_argument("--filtered-barcodes", default=None,
                        help="Path to CellRanger barcodes.tsv.gz (from filtered or raw count "
                             "matrix). When provided, only these barcodes are used for "
                             "correction instead of the full whitelist, mimicking the "
                             "BLAZE-filtered approach used by scnanoseq.")
    parser.add_argument("--max-bc-edit-dist", type=int, default=2,
                        help="Max edit distance for cell barcode correction (default: 2)")
    parser.add_argument("--min-posterior-prob", type=float, default=0.975,
                        help="Min posterior probability for barcode correction (default: 0.975)")
    return parser.parse_args()


def load_staggers(path):
    """Parse stagger file. Returns list of (sample_name, stagger_sequence) tuples."""
    samples = []
    seen_names = set()
    with open(path) as f:
        for line_num, line in enumerate(f, 1):
            line = line.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) == 1:
                name = parts[0].strip()
                stagger_seq = ""
            elif len(parts) == 2:
                name = parts[0].strip()
                stagger_seq = parts[1].strip().upper()
            else:
                sys.exit(f"ERROR: staggers.txt line {line_num}: expected 1-2 tab-separated "
                         f"fields (sample_name, [stagger_sequence]), got {len(parts)}")

            if name in seen_names:
                sys.exit(f"ERROR: staggers.txt line {line_num}: duplicate sample name '{name}'")
            seen_names.add(name)

            if stagger_seq and not set(stagger_seq).issubset(VALID_BASES):
                bad = set(stagger_seq) - VALID_BASES
                sys.exit(f"ERROR: staggers.txt line {line_num}: stagger for '{name}' "
                         f"contains invalid characters: {bad}")

            samples.append((name, stagger_seq))

    if not samples:
        sys.exit("ERROR: staggers.txt contains no samples")
    return samples


def get_primer_seq(mode):
    """Return primer sequence based on mode."""
    if mode == "GFP":
        return "GGACGAGCTGTACAAGTAGG"  # 20bp
    elif mode == "Vpx":
        return "CTAGGGGAAGGACATGGGGCAGG"  # 23bp
    else:
        sys.exit(f"Unknown mode: {mode}")


def has_kmer_hit(seq, kmer_list):
    """Fast check: does seq contain any k-mer from the list?
    Uses str.find() (C-level) instead of Python-level character iteration.
    Returns the position of the first hit, or -1 if none found."""
    best = -1
    for km in kmer_list:
        pos = seq.find(km)
        if pos >= 0:
            if best < 0 or pos < best:
                best = pos
            # Early exit: can't find anything earlier than position 0
            if best == 0:
                return 0
    return best


def find_adapter_and_extract(seq, qual, adapter_pattern, adapter_mismatches):
    """Search for the 10X adapter in a sequence and extract cell barcode + UMI.

    Uses a three-tier approach per orientation (forward, then reverse complement):
      1. Exact str.find() for the full adapter (~80% of hits, essentially free)
      2. K-mer pre-screen to locate candidate region, then fuzzy regex in a small window
      3. No full-read fuzzy fallback (too slow on long ONT reads)

    Returns (cell_barcode, cb_qual, umi, cDNA, orientation) or None if adapter not found.
    """
    needed = CELL_BARCODE_LEN + UMI_LEN + MIN_POLYT_LEN
    slen = len(seq)

    # --- Forward orientation ---
    # Tier 1: exact match (C-level str.find, very fast)
    pos = seq.find(TENX_ADAPTER)
    if pos >= 0:
        result = _extract_from_pos(seq, qual, pos + TENX_ADAPTER_LEN, needed, "fwd")
        if result:
            return result

    # Tier 2: k-mer pre-screen + fuzzy regex in small window
    hit_pos = has_kmer_hit(seq, ADAPTER_KMERS)
    if hit_pos >= 0:
        search_start = max(0, hit_pos - TENX_ADAPTER_LEN)
        search_end = min(slen, hit_pos + TENX_ADAPTER_LEN * 2)
        match = adapter_pattern.search(seq, pos=search_start, endpos=search_end)
        if match:
            result = _extract_from_match(seq, qual, match, needed, "fwd")
            if result:
                return result

    # --- Reverse complement orientation ---
    rc_seq = reverse_complement(seq)
    rc_qual = qual[::-1]

    # Tier 1: exact match on RC
    pos = rc_seq.find(TENX_ADAPTER)
    if pos >= 0:
        result = _extract_from_pos(rc_seq, rc_qual, pos + TENX_ADAPTER_LEN, needed, "rev")
        if result:
            return result

    # Tier 2: k-mer pre-screen + fuzzy regex in small window on RC
    hit_pos = has_kmer_hit(rc_seq, ADAPTER_KMERS)
    if hit_pos >= 0:
        search_start = max(0, hit_pos - TENX_ADAPTER_LEN)
        search_end = min(slen, hit_pos + TENX_ADAPTER_LEN * 2)
        match = adapter_pattern.search(rc_seq, pos=search_start, endpos=search_end)
        if match:
            result = _extract_from_match(rc_seq, rc_qual, match, needed, "rev")
            if result:
                return result

    return None


def _extract_from_pos(seq, qual, after_adapter, needed, orientation):
    """Extract cell barcode, CB quality, UMI, and cDNA from a known adapter end position."""
    remaining = len(seq) - after_adapter
    if remaining < needed:
        return None

    cb = seq[after_adapter:after_adapter + CELL_BARCODE_LEN]
    cb_qual = qual[after_adapter:after_adapter + CELL_BARCODE_LEN]
    umi = seq[after_adapter + CELL_BARCODE_LEN:after_adapter + CELL_BARCODE_LEN + UMI_LEN]
    post_umi = seq[after_adapter + CELL_BARCODE_LEN + UMI_LEN:]

    polyt_region = post_umi[:MIN_POLYT_LEN]
    t_count = polyt_region.count("T")
    if t_count < MIN_POLYT_LEN - 2:
        return None

    cDNA_start = 0
    for i, base in enumerate(post_umi):
        if base != "T":
            if i + 1 < len(post_umi) and post_umi[i + 1] == "T":
                continue
            cDNA_start = i
            break
    cDNA = post_umi[max(cDNA_start, MIN_POLYT_LEN):]
    return cb, cb_qual, umi, cDNA, orientation


def _extract_from_match(seq, qual, match, needed, orientation):
    """Extract cell barcode, CB quality, UMI, and cDNA from an adapter match."""
    after_adapter = match.end()
    remaining = len(seq) - after_adapter
    if remaining < needed:
        return None

    cb = seq[after_adapter:after_adapter + CELL_BARCODE_LEN]
    cb_qual = qual[after_adapter:after_adapter + CELL_BARCODE_LEN]
    umi = seq[after_adapter + CELL_BARCODE_LEN:after_adapter + CELL_BARCODE_LEN + UMI_LEN]
    post_umi = seq[after_adapter + CELL_BARCODE_LEN + UMI_LEN:]

    # Verify polyT presence (allows some errors)
    polyt_region = post_umi[:MIN_POLYT_LEN]
    t_count = polyt_region.count("T")
    if t_count < MIN_POLYT_LEN - 2:
        return None

    # Find end of polyT stretch, then cDNA follows
    cDNA_start = 0
    for i, base in enumerate(post_umi):
        if base != "T":
            if i + 1 < len(post_umi) and post_umi[i + 1] == "T":
                continue
            cDNA_start = i
            break
    cDNA = post_umi[max(cDNA_start, MIN_POLYT_LEN):]
    return cb, cb_qual, umi, cDNA, orientation


# ---------------------------------------------------------------------------
# Cell barcode error correction (following scnanoseq / correct_barcodes.py)
# ---------------------------------------------------------------------------

QUAL_OFFSET = 33  # Phred+33 (standard for both Illumina and ONT)


def load_whitelist(path):
    """Load a 10X barcode whitelist into a Python set for O(1) lookup.
    Supports both plain text and gzipped (.gz) files."""
    whitelist = set()
    opener = gzopen if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            bc = line.strip().split()[0]  # take first column only
            if bc:
                whitelist.add(bc)
    return whitelist


def load_filtered_barcodes(path):
    """Load CellRanger barcodes.tsv.gz and strip the '-1' suffix.
    Returns a set of 16bp cell barcodes."""
    barcodes = set()
    opener = gzopen if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        for line in f:
            bc = line.strip()
            if bc:
                # Strip the -1 (or -N) suffix that CellRanger appends
                if "-" in bc:
                    bc = bc.rsplit("-", 1)[0]
                barcodes.add(bc)
    return barcodes


def generate_mutations(barcode, max_edit_dist):
    """Generate all sequences within Hamming distance max_edit_dist of barcode.

    For 16bp with dist=2: 48 + 1080 = 1128 candidates — very fast to enumerate.
    Follows scnanoseq get_mutated_bcs() approach (substitutions only).
    """
    for edit_dist in range(1, max_edit_dist + 1):
        for locs in itertools.combinations(range(len(barcode)), edit_dist):
            bc_list = [[base] for base in barcode]
            for loc in locs:
                orig_base = bc_list[loc][0]
                bc_list[loc] = [b for b in "ACGT" if b != orig_base]
            for combo in itertools.product(*bc_list):
                yield "".join(combo)


def error_probability_at(qual_char):
    """Convert a Phred+33 quality character to the probability of sequencing error."""
    return math.pow(10, -1 * (ord(qual_char) - QUAL_OFFSET) / 10)


def correct_cell_barcode(barcode, barcode_qual, whitelist, bc_counts,
                         max_edit_dist, min_posterior_prob):
    """Correct a cell barcode against the 10X whitelist.

    Following scnanoseq correct_barcodes.py:
      1. If barcode is already in the whitelist, return it as-is.
      2. Generate all Hamming-distance 1..max_edit_dist mutants.
      3. For each mutant in the whitelist, score it using:
         - Product of error probabilities at mismatch positions (from quality scores)
         - Barcode abundance prior (count + 1 pseudo-count, as in scnanoseq)
      4. Pick the candidate with the highest posterior probability.
      5. Only accept if posterior > min_posterior_prob (default 0.975).

    Returns (corrected_barcode, was_corrected) or (None, False) if uncorrectable.
    """
    if barcode in whitelist:
        return barcode, False

    candidates = []
    for mutant in generate_mutations(barcode, max_edit_dist):
        if mutant in whitelist:
            # Calculate likelihood: product of error probs at mismatch positions
            likelihood = 1.0
            for i in range(len(barcode)):
                if barcode[i] != mutant[i]:
                    likelihood *= error_probability_at(barcode_qual[i])
            # Prior: barcode abundance + 1 pseudo-count (as in scnanoseq)
            prior = bc_counts.get(mutant, 0) + 1
            score = likelihood * prior
            if score > 0:
                candidates.append((mutant, score))

    if not candidates:
        return None, False

    # Posterior probability: normalize scores
    total_score = sum(s for _, s in candidates)
    best_bc, best_score = max(candidates, key=lambda x: x[1])
    posterior = best_score / total_score

    if posterior >= min_posterior_prob:
        return best_bc, True
    return None, False


def read_fastq_gz(filepath):
    """FASTQ.gz reader yielding (name, sequence, quality) tuples."""
    with gzopen(filepath, "rt") as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().rstrip("\n")
            f.readline()  # + line
            qual = f.readline().rstrip("\n")
            yield header.rstrip("\n"), seq, qual


class SampleState:
    """Per-sample output file handles, counters, and uniqueness sets."""

    # Max bases from the end of cDNA to search for the RC primer
    PRIMER_SEARCH_TAIL = 150

    def __init__(self, name, stagger_seq, primer_seq, outdir, max_primer_errors=6):
        self.name = name
        self.stagger_seq = stagger_seq
        self.stagger_len = len(stagger_seq)
        self.primer_seq = primer_seq
        self.primer_len = len(primer_seq)
        self.region_len = self.stagger_len + self.primer_len

        # RC of stagger and primer for searching the 3' end of cDNA
        self.primer_rc = reverse_complement(primer_seq)
        self.stagger_rc = reverse_complement(stagger_seq) if stagger_seq else ""

        # --- Forward pattern: stagger + primer at START of cDNA (Illumina-like) ---
        if stagger_seq:
            self.demux_pattern = regex.compile(
                rf'({regex.escape(stagger_seq)}){{e<=2}}({primer_seq}){{e<={max_primer_errors}}}',
                flags=regex.BESTMATCH
            )
        else:
            self.demux_pattern = regex.compile(
                rf'({primer_seq}){{e<={max_primer_errors}}}',
                flags=regex.BESTMATCH
            )

        # --- RC pattern: primer_RC + stagger_RC at END of cDNA (ONT typical) ---
        # In ONT reads: [...transcript...][barcode_RC][primer_RC][stagger_RC]
        if stagger_seq:
            self.demux_pattern_rc = regex.compile(
                rf'({self.primer_rc}){{e<={max_primer_errors}}}({regex.escape(self.stagger_rc)}){{e<=2}}',
                flags=regex.BESTMATCH
            )
        else:
            self.demux_pattern_rc = regex.compile(
                rf'({self.primer_rc}){{e<={max_primer_errors}}}',
                flags=regex.BESTMATCH
            )

        self.sample_dir = os.path.join(outdir, name)
        os.makedirs(self.sample_dir, exist_ok=True)

        # Output file handles
        self.f_joined = open(os.path.join(self.sample_dir, "joinedCBUMIBarcode.txt"), "w")
        self.f_bad_barcode = open(os.path.join(self.sample_dir, "badBarcode.txt"), "w")
        self.f_unique_screened = open(os.path.join(self.sample_dir, "uniqueScreenedReads.txt"), "w")
        self.f_unique_shaved = open(os.path.join(self.sample_dir, "uniqueShavedReads.txt"), "w")

        # Counters
        self.total_by_primer = 0
        self.bad_barcode_count = 0
        self.screened_count = 0
        self.shaved_count = 0

        # Uniqueness sets
        self.unique_screened = set()
        self.unique_shaved = set()

    def close(self):
        for f in (self.f_joined, self.f_bad_barcode,
                  self.f_unique_screened, self.f_unique_shaved):
            f.close()

    def write_summary(self):
        with open(os.path.join(self.sample_dir, "summaryFile.txt"), "w") as summary:
            summary.write(f"sample {self.name}\n")
            summary.write(f"stagger_sequence {self.stagger_seq if self.stagger_seq else '(none)'}\n")
            summary.write(f"stagger_length {self.stagger_len}\n")
            summary.write(f"primer_region_length {self.region_len}\n")
            summary.write(f"total assigned by stagger+primer {self.total_by_primer}\n")
            summary.write(f"total badBarcode reads {self.bad_barcode_count}\n")
            summary.write(f"total screenedReads reads {self.screened_count}\n")
            summary.write(f"total screened Unique Reads reads {len(self.unique_screened)}\n")
            summary.write(f"total shavedReads reads {self.shaved_count}\n")
            summary.write(f"total shaved Unique Reads reads {len(self.unique_shaved)}\n")


def main():
    args = parse_args()

    # --- Load stagger definitions ---
    samples = load_staggers(args.staggers)
    primer_seq = get_primer_seq(args.mode)

    if args.mode == "Vpx":
        overridden = [(name, s) for name, s in samples if s]
        if overridden:
            for name, s in overridden:
                print(f"NOTE: stagger sequence for '{name}' ('{s}') ignored in Vpx mode.")
        samples = [(name, "") for name, _ in samples]

    # --- Build SampleState objects ---
    all_states = []
    for name, stagger_seq in samples:
        state = SampleState(name, stagger_seq, primer_seq, args.outdir,
                            max_primer_errors=args.primer_mismatches)
        all_states.append(state)
    all_states.sort(key=lambda s: s.stagger_len, reverse=True)

    # --- Load whitelist for cell barcode correction ---
    whitelist = None
    bc_counts = {}  # tracks observed barcode abundance for prior
    if args.filtered_barcodes:
        if args.whitelist:
            print("NOTE: --filtered-barcodes provided, ignoring --whitelist")
        print(f"Loading filtered barcodes: {args.filtered_barcodes}...", flush=True)
        whitelist = load_filtered_barcodes(args.filtered_barcodes)
        print(f"  Loaded {len(whitelist):,} cell barcodes (CellRanger-filtered)")
        print(f"  Max CB edit distance: {args.max_bc_edit_dist}")
        print(f"  Min posterior probability: {args.min_posterior_prob}")
    elif args.whitelist:
        print(f"Loading full whitelist: {args.whitelist}...", flush=True)
        whitelist = load_whitelist(args.whitelist)
        print(f"  Loaded {len(whitelist):,} whitelisted barcodes")
        print(f"  Max CB edit distance: {args.max_bc_edit_dist}")
        print(f"  Min posterior probability: {args.min_posterior_prob}")

    # --- Compile adapter pattern with fuzzy matching ---
    # BESTMATCH stops early once optimal match is found, avoiding full backtracking
    adapter_pattern = regex.compile(
        rf'({TENX_ADAPTER}){{e<={args.adapter_mismatches}}}',
        flags=regex.BESTMATCH
    )

    # --- Constants ---
    bad_barcode_pattern = regex.compile(r'AAAA|TTTT|GGGG|CCCC|NN')

    # --- Global counters ---
    total_reads = 0
    no_adapter_reads = 0
    no_primer_reads = 0
    cb_exact_match = 0
    cb_corrected = 0
    cb_uncorrectable = 0

    print(f"Mode: {args.mode} | Primer: {primer_seq} ({len(primer_seq)}bp)")
    print(f"Adapter mismatches allowed: {args.adapter_mismatches}")
    print(f"Primer mismatches allowed: {args.primer_mismatches}")
    print(f"10X adapter: {TENX_ADAPTER} ({TENX_ADAPTER_LEN}bp)")
    print(f"Cell barcode length: {CELL_BARCODE_LEN}bp | UMI length: {UMI_LEN}bp")
    print(f"Samples: {len(all_states)}")
    for state in all_states:
        stag_label = state.stagger_seq if state.stagger_seq else "(none)"
        print(f"  {state.name}: stagger={stag_label} "
              f"({state.stagger_len}bp), region={state.region_len}bp")

    # --- Single-pass processing ---
    print("Processing reads (single pass, extracting barcode/UMI + primer matching)...",
          flush=True)

    try:
        for header, seq, qual in read_fastq_gz(args.fastq):
            total_reads += 1

            # --- Step 1: Find 10X adapter, extract cell barcode + UMI ---
            # Progress reporting
            if total_reads % 100_000 == 0:
                print(f"  Processed {total_reads:,} reads "
                      f"(adapter: {total_reads - no_adapter_reads:,}, "
                      f"no adapter: {no_adapter_reads:,}, "
                      f"no primer: {no_primer_reads:,})...", flush=True)

            result = find_adapter_and_extract(seq, qual, adapter_pattern, args.adapter_mismatches)
            if result is None:
                no_adapter_reads += 1
                continue

            cell_barcode, cb_qual, umi, cDNA, orientation = result

            # --- Step 1b: Cell barcode error correction ---
            if whitelist is not None:
                corrected_cb, was_corrected = correct_cell_barcode(
                    cell_barcode, cb_qual, whitelist, bc_counts,
                    args.max_bc_edit_dist, args.min_posterior_prob)
                if corrected_cb is None:
                    cb_uncorrectable += 1
                    continue
                if was_corrected:
                    cb_corrected += 1
                else:
                    cb_exact_match += 1
                # Update abundance counts for the prior
                bc_counts[corrected_cb] = bc_counts.get(corrected_cb, 0) + 1
                cell_barcode = corrected_cb

            # --- Step 2: Match stagger+primer in cDNA to assign sample ---
            # ONT reads have the construct at the 3' end of cDNA in RC orientation:
            #   [...transcript...][barcode_RC][primer_RC][stagger_RC]
            # Three-tier search: exact str.find (free), then fuzzy regex on tail, then forward.
            matched_state = None
            barcode = None
            cDNA_len = len(cDNA)
            for state in all_states:
                # Tier 1: exact str.find for primer_RC (very fast, ~34% of hits)
                rc_pos = cDNA.rfind(state.primer_rc)
                if rc_pos >= 0:
                    matched_state = state
                    bc_rc_end = rc_pos
                    bc_rc_start = max(0, bc_rc_end - 70)
                    barcode = reverse_complement(cDNA[bc_rc_start:bc_rc_end])
                    break
                # Tier 2: fuzzy regex for primer_RC in last PRIMER_SEARCH_TAIL bp
                tail_start = max(0, cDNA_len - SampleState.PRIMER_SEARCH_TAIL)
                m = state.demux_pattern_rc.search(cDNA, pos=tail_start)
                if m:
                    matched_state = state
                    bc_rc_end = m.start()
                    bc_rc_start = max(0, bc_rc_end - 70)
                    barcode = reverse_complement(cDNA[bc_rc_start:bc_rc_end])
                    break
                # Tier 3: forward pattern at start of cDNA (first 100bp)
                m = state.demux_pattern.search(cDNA, endpos=min(100, cDNA_len))
                if m:
                    matched_state = state
                    barcode_start = m.end()
                    barcode = cDNA[barcode_start:barcode_start + 70]
                    break

            if matched_state is None:
                no_primer_reads += 1
                continue

            matched_state.total_by_primer += 1

            if not barcode:
                continue

            matched_state.f_joined.write(f"{cell_barcode} {umi} {barcode}\n")

            # --- Step 3: Barcode quality filters ---
            if args.mode != "Vpx" and bad_barcode_pattern.search(barcode):
                matched_state.bad_barcode_count += 1
                matched_state.f_bad_barcode.write(f"{cell_barcode} {umi} {barcode}\n")
                continue

            # Passed all filters
            matched_state.screened_count += 1
            screened_key = (cell_barcode, umi, barcode)
            if screened_key not in matched_state.unique_screened:
                matched_state.unique_screened.add(screened_key)
                matched_state.f_unique_screened.write(f"{cell_barcode} {umi} {barcode}\n")

            shaved_key = (cell_barcode, umi, barcode)
            matched_state.shaved_count += 1
            if shaved_key not in matched_state.unique_shaved:
                matched_state.unique_shaved.add(shaved_key)
                matched_state.f_unique_shaved.write(f"{cell_barcode} {umi} {barcode}\n")

    finally:
        for state in all_states:
            state.close()

    # --- Reporting ---
    print(f"Finished processing {total_reads:,} reads.")
    print(f"No adapter found: {no_adapter_reads:,}")
    if whitelist is not None:
        print(f"Cell barcode correction:")
        print(f"  Exact whitelist match: {cb_exact_match:,}")
        print(f"  Corrected: {cb_corrected:,}")
        print(f"  Uncorrectable (dropped): {cb_uncorrectable:,}")
    print(f"No primer match: {no_primer_reads:,}")

    for state in all_states:
        state.write_summary()
        print(f"  {state.name}: "
              f"{state.total_by_primer:,} by stagger+primer, "
              f"{state.screened_count:,} screened, "
              f"{len(state.unique_shaved):,} unique shaved")

    # --- Write global summary ---
    with open(os.path.join(args.outdir, "globalSummary.txt"), "w") as gs:
        gs.write(f"mode {args.mode}\n")
        gs.write(f"primer {primer_seq}\n")
        gs.write(f"primer_length {len(primer_seq)}\n")
        gs.write(f"adapter_mismatches_allowed {args.adapter_mismatches}\n")
        gs.write(f"primer_mismatches_allowed {args.primer_mismatches}\n")
        gs.write(f"total raw reads {total_reads}\n")
        gs.write(f"no adapter found {no_adapter_reads}\n")
        if whitelist is not None:
            gs.write(f"cb_exact_whitelist_match {cb_exact_match}\n")
            gs.write(f"cb_corrected {cb_corrected}\n")
            gs.write(f"cb_uncorrectable {cb_uncorrectable}\n")
        gs.write(f"no primer match {no_primer_reads}\n")
        for state in all_states:
            stag_label = state.stagger_seq if state.stagger_seq else "(none)"
            gs.write(f"sample {state.name} "
                     f"stagger {stag_label} "
                     f"reads_by_primer {state.total_by_primer} "
                     f"screened {state.screened_count} "
                     f"unique_shaved {len(state.unique_shaved)}\n")

    print("Done. Output written to: " + args.outdir)


if __name__ == "__main__":
    main()
