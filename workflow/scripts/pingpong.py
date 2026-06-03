#!/usr/bin/env python3
"""
Ping‑pong signature analysis for piRNAs.

Input: BAM file of small RNA reads mapped to repeat sequences only.
Example line of BAM:
13314-67        16      L1MdTf_III      5457    255     29M     *       0       0       CACCATCCTTAGTTATCAGGGAAATGCAA   IIIIIIIIIIIIIIIIIIIIIIIIIIIII   XA:i:3  MD:Z:2A8A1C15   NM:i:3  XM:i:3

Note: the read name '13314-67' encodes that this sequence was observed 67 times in the original collapsed FASTA file.

Output:
  - CSV with columns: repeat_id, distance, count, sample
    Where 'distance' is the 5' overlap distance between all combinations of
    sense and antisense reads.
  - CSV with columns: repeat_id, strand, nucleotide, count, sample
    Nucleotide frequencies at ping-pong positions for pairs at distance 10:
      sense strand    -> nucleotide at position 10 of the sense read
      antisense strand -> nucleotide at position 1 of the antisense read

According to description in De Fazio et al. 2011 (Nature):
For ping-pong analysis, only reads mapping to repeat elements were considered.
For each pair of sense/antisense overlapping reads, the distance between their
5′ ends was recorded and counts were represented as relative frequencies within
samples for each repeat element.

The steps in this script, replicate the Perl scripts in:
https://github.com/rberrens/SPOCD1-piRNA_directed_DNA_met/tree/master/piRNA_analysis

Essence of the ping-pong analysis:

+ read:          ------>           5′ = 100
- read:     <------                 5′ = 110
Distance:        --

Distance = |110 - 100| = 10 nt

This is calculated for all combinations of sense and antisense reads mapping
to the same repeat element.
(Restricted to distances up to a certain window, e.g. 30 nt.)

"""

import csv
from collections import defaultdict, Counter
import pysam


# Helper functions
# --------------------------------------------------
def load_reads(bam_path):
    bam = pysam.AlignmentFile(bam_path, "rb")
    for read in bam.fetch(until_eof=True):
        # Yield will make one read available at a time,
        # rather than loading all reads into a list (generator)
        yield read


def read_is_sense(read):
    return not read.is_reverse


def read_count_from_header(read):
    """
    Extract the original count from the read name, e.g., 'seq1-12'
    """
    return int(read.query_name.split("-")[1])


# Main ping‑pong computation
# --------------------------------------------------
def compute_ping_pong(bam_path, window=30):
    """
    Compute ping‑pong distances and nucleotide bias at distance 10.

    window: maximum distance to record (e.g. 30 nt).
    Returns:
      results:   dict {repeat_id -> Counter(distance)}
      nt_counts: dict {repeat_id -> {"sense": Counter(nt),
                                     "antisense": Counter(nt)}}
                 Counts are weighted by the number of pairs at distance 10.
                 sense     -> position 10 nucleotide of the sense read
                 antisense -> position 1 nucleotide of the antisense read
    """
    # Lookup table for reverse-complementing a single base.
    RC = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}

    # data structure: te_list[repeat_id]["F" or "R"][five_p] = [(count, nt), ...]
    # F (sense) reads store the nucleotide at position 10 of the read (10A bias).
    # R (antisense) reads store the nucleotide at position 1 of the read (1U bias).
    te_list = defaultdict(lambda: {"F": defaultdict(list), "R": defaultdict(list)})

    for read in load_reads(bam_path):
        te_id = read.reference_name
        count = read_count_from_header(read)
        seq = read.query_sequence or ""

        if read_is_sense(read):
            # Sense (+ strand) read: SEQ is stored 5'→3', so seq[9] is
            # directly position 10 of the piRNA (10A bias).
            five_p = read.reference_start
            nt = seq[9] if len(seq) >= 10 else None
            te_list[te_id]["F"][five_p].append((count, nt))
        else:
            # Antisense (- strand) read: BAM stores SEQ as the reverse
            # complement of the original read, so seq[0] is the complement
            # of the 3' end — NOT position 1. Position 1 (5' end) of the
            # original piRNA is encoded as complement(seq[-1]).
            five_p = read.reference_start + read.query_length
            nt = RC.get(seq[-1]) if seq else None
            te_list[te_id]["R"][five_p].append((count, nt))

    # Compute distances
    # defaultdict of Counter allows safe incrementing,
    # even when keys do not yet exist
    results = defaultdict(Counter)
    nt_counts = defaultdict(lambda: {"sense": Counter(), "antisense": Counter()})

    for te, strands in te_list.items():
        F = strands["F"]
        R = strands["R"]
        for f_pos, f_reads in F.items():
            sum_F = sum(c for c, _ in f_reads)
            for r_pos, r_reads in R.items():
                sum_R = sum(c for c, _ in r_reads)

                # Calculate absolute distance
                distance = abs(f_pos - r_pos)
                if distance <= window:
                    # sum_F and sum_R are sums of counts for reads at these
                    # positions. Each read at f_pos can pair with each read
                    # at r_pos.
                    results[te][distance] += sum_F * sum_R

                    if distance == 10:
                        # Weight each sense read's pos-10 nt by the total
                        # antisense count it pairs with (sum_R).
                        for c_F, nt_F in f_reads:
                            if nt_F is not None:
                                nt_counts[te]["sense"][nt_F] += c_F * sum_R
                        # Weight each antisense read's pos-1 nt by the total
                        # sense count it pairs with (sum_F).
                        for c_R, nt_R in r_reads:
                            if nt_R is not None:
                                nt_counts[te]["antisense"][nt_R] += sum_F * c_R

    return results, nt_counts


# Get inputs/outputs from Snakemake
# --------------------------------------------------
bam = snakemake.input["bam"]
window = snakemake.params["window"]
sample = snakemake.wildcards["sample"]
out = snakemake.output["pingpong"]
out_nt_bias = snakemake.output["nt_bias"]

# Run ping‑pong analysis
res, nt_counts = compute_ping_pong(bam, window=window)

# Write distance CSV
with open(out, "w", newline="") as fh:
    writer = csv.writer(fh)
    writer.writerow(["repeat_id", "distance", "count", "sample"])
    for te, counter in res.items():
        for d, c in sorted(counter.items()):
            writer.writerow([te, d, c, sample])

# Remove "N" from nt_counts before writing nucleotide bias CSV
nt_counts = {
    te: {
        strand: Counter({nt: c for nt, c in counter.items() if nt != "N"})
        for strand, counter in strands.items()
    }
    for te, strands in nt_counts.items()
}

# Write nucleotide bias CSV (ping-pong pairs at distance 10 only)
with open(out_nt_bias, "w", newline="") as fh:
    writer = csv.writer(fh)
    writer.writerow(["repeat_id", "strand", "nucleotide", "count", "sample"])
    for te, strands in nt_counts.items():
        for strand, counter in strands.items():
            for nt, c in sorted(counter.items()):
                writer.writerow([te, strand, nt, c, sample])
