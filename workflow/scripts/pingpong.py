#!/usr/bin/env python3
"""
Ping‑pong signature analysis for piRNAs.

Input: BAM file of small RNA reads mapped to repeat sequences only.
Example line of BAM:
13314-67        16      L1MdTf_III      5457    255     29M     *       0       0       CACCATCCTTAGTTATCAGGGAAATGCAA   IIIIIIIIIIIIIIIIIIIIIIIIIIIII   XA:i:3  MD:Z:2A8A1C15   NM:i:3  XM:i:3

Note: the read name '13314-67' encodes that this sequence was observed 67 times in the original collapsed FASTA file.

Output: CSV with columns: repeat_id, distance, count
Where 'distance' is the 5' overlap distance between all combinations of sense and antisense reads.

According to description in De Fazio et al. 2011 (Nature):
For ping-pong analysis, only reads mapping to repeat elements were considered. For each pair of sense/antisense overlapping reads, the distance between their 5′ ends was recorded and counts were represented as relative frequencies within samples for each repeat element.

The steps in this script, replicate the Perl scripts in:
https://github.com/rberrens/SPOCD1-piRNA_directed_DNA_met/tree/master/piRNA_analysis

Essence of the ping-pong analysis:

+ read:          ------>           5′ = 100
- read:     <------                 5′ = 110
Distance:        --

Distance = |110 - 100| = 10 nt

This is calculated for all combinations of sense and antisense reads mapping to the same repeat element.
(Restricted to distances up to a certain window, e.g. 30 nt.)

"""

import csv
import argparse
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
    Compute ping‑pong distances.

    window: maximum distance to record (e.g. 30 nt).
    Returns: dict {repeat_id -> Counter(distance)}
    """
    # data structure: te_list[repeat_id]["F" or "R"][five_p] = [list of counts]
    te_list = defaultdict(lambda: {"F": defaultdict(list), "R": defaultdict(list)})

    for read in load_reads(bam_path):
        te_id = read.reference_name
        count = read_count_from_header(read)

        # Compute 5' end
        if read_is_sense(read):
            five_p = read.reference_start
        else:
            five_p = read.reference_start + read.query_length

        if read_is_sense(read):
            te_list[te_id]["F"][five_p].append(count)
        else:
            te_list[te_id]["R"][five_p].append(count)

    # Compute distances
    # defaultdict of Counter allows safe incrementing,
    # even when keys do not yet exist
    results = defaultdict(Counter)

    for te, strands in te_list.items():
        F = strands["F"]
        R = strands["R"]
        for f_pos, f_counts in F.items():
            # f_pos: position of 5' end of forward read
            # f_counts: list of counts for reads at this position
            sum_F = sum(f_counts)
            for r_pos, r_counts in R.items():
                # r_pos: position of 5' end of reverse read
                # r_counts: list of counts for reads at this position
                sum_R = sum(r_counts)

                # Calculate absolute distance
                distance = abs(f_pos - r_pos)
                if distance <= window:
                    # sum_F and sum_R are sums of counts for reads at these positions
                    # Each read at f_pos can pair with each read at r_pos
                    results[te][distance] += sum_F * sum_R

    return results


# Command line interface
# --------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Ping-pong signature analysis")
    parser.add_argument("bam", help="Input BAM file of small RNA alignment to repeats")
    parser.add_argument("--window", type=int, default=30, help="Distance window")
    parser.add_argument(
        "--out", default="pingpong.csv", help="Output CSV file (default: pingpong.csv)"
    )
    args = parser.parse_args()

    res = compute_ping_pong(args.bam, window=args.window)

    # Write CSV
    with open(args.out, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["repeat_id", "distance", "count"])
        for te, counter in res.items():
            for d, c in sorted(counter.items()):
                writer.writerow([te, d, c])
