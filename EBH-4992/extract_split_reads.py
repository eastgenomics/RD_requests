#!/usr/bin/env python3
"""
Extract soft-clipped reads from a BAM file at a specified region.
Outputs clipped sequences as FASTA for upload to NCBI BLAST.

Usage:
    python extract_split_reads.py -b sample1.bam -r chr9:132944000-132946000 -o clipped_reads.fa

Requirements:
    pip install pysam
"""

import sys
import re
from collections import defaultdict

try:
    import pysam
except ImportError:
    sys.exit("ERROR: pysam not installed. Run: pip install pysam")

BAM = "141351816-26013R0074-26NGCEN7-9527-M-99347387_markdup.bam"
REGION = "9:132935000-132946000"
OUTPUT = "tsc1_clipped_reads.fa"

MIN_CLIP_LEN = 20
MIN_MAPQ = 20
MIN_SUPPORTING_READS = 2


def get_clipped_sequences(read):
    """Return list of (genomic_position, side, sequence) for soft-clipped ends."""
    clips = []
    cigar = read.cigartuples
    if cigar is None:
        return clips

    SOFT_CLIP = 4

    if cigar[0][0] == SOFT_CLIP and cigar[0][1] >= MIN_CLIP_LEN:
        clips.append((read.reference_start, "5prime", read.query_sequence[:cigar[0][1]]))

    if cigar[-1][0] == SOFT_CLIP and cigar[-1][1] >= MIN_CLIP_LEN:
        clips.append((read.reference_end, "3prime", read.query_sequence[-cigar[-1][1]:]))

    return clips


def extract_split_reads(bam_path, region):
    match = re.match(r"(.+):(\d+)-(\d+)", region)
    if not match:
        sys.exit(f"ERROR: Could not parse region '{region}'. Expected format: chr9:132944000-132946000")
    chrom, start, end = match.group(1), int(match.group(2)), int(match.group(3))

    bam = pysam.AlignmentFile(bam_path, "rb")

    clips_by_position = defaultdict(list)
    total_reads = 0
    clipped_reads = 0

    print(f"Scanning {region} for soft-clipped reads (min clip: {MIN_CLIP_LEN}bp, min MAPQ: {MIN_MAPQ})...")

    for read in bam.fetch(chrom, start, end):
        if read.is_unmapped or read.is_duplicate or read.mapping_quality < MIN_MAPQ:
            continue

        total_reads += 1
        clips = get_clipped_sequences(read)

        if clips:
            clipped_reads += 1
            for pos, side, seq in clips:
                clips_by_position[pos].append((side, seq))

    bam.close()

    print(f"  Total reads: {total_reads:,}")
    print(f"  Reads with soft clips >= {MIN_CLIP_LEN}bp: {clipped_reads:,}")
    print(f"  Unique clip positions: {len(clips_by_position)}")

    return clips_by_position


def summarise_positions(clips_by_position):
    print(f"\n{'Position':<20} {'Count':>6} {'5prime':>8} {'3prime':>8}")
    print("-" * 46)
    for pos, clips in sorted(clips_by_position.items(), key=lambda x: len(x[1]), reverse=True)[:20]:
        n_5 = sum(1 for side, _ in clips if side == "5prime")
        n_3 = sum(1 for side, _ in clips if side == "3prime")
        print(f"{pos:<20} {len(clips):>6} {n_5:>8} {n_3:>8}")


def write_fasta(clips_by_position, output_path):
    sorted_positions = sorted(clips_by_position.items(), key=lambda x: len(x[1]), reverse=True)
    seen_seqs = set()
    written = 0

    with open(output_path, "w") as f:
        for pos, clips in sorted_positions:
            if len(clips) < MIN_SUPPORTING_READS:
                continue
            for i, (side, seq) in enumerate(clips):
                if seq in seen_seqs or len(set(seq)) < 3:
                    continue
                seen_seqs.add(seq)
                f.write(f">pos_{pos}_{side}_{i+1}_len{len(seq)}\n{seq}\n")
                written += 1

    print(f"\nWrote {written} unique clipped sequences to: {output_path}")
    print(f"Upload {output_path} to https://blast.ncbi.nlm.nih.gov/Blast.cgi")
    return written


def main():
    clips_by_position = extract_split_reads(BAM, REGION)

    if not clips_by_position:
        print("No soft-clipped reads found. Try widening the region or lowering MIN_CLIP_LEN.")
        sys.exit(0)

    summarise_positions(clips_by_position)
    written = write_fasta(clips_by_position, OUTPUT)
    if written == 0:
        print("No sequences written. Try lowering MIN_SUPPORTING_READS.")


if __name__ == "__main__":
    main()
