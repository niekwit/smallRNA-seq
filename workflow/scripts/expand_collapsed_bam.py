#!/usr/bin/env python3
"""
Expand a BAM file of collapsed small RNA reads back to per-read alignments.

Input: BAM file where identical sequences were collapsed into a single
alignment record, with the read abundance encoded as a "-<count>" suffix
on the read name (e.g. "455305-1", "2653-13").

Output: BAM file (unsorted, same order as input) where each alignment
record is duplicated <count> times, so downstream tools that count reads
(e.g. deepTools bamCoverage) see the true per-read abundance instead of
one record per unique sequence.
"""

import pysam


def read_count(read_name):
    return int(read_name.rsplit("-", 1)[-1])


in_bam = snakemake.input["bam"]
out_bam = snakemake.output["expanded_bam"]
threads = snakemake.threads

with pysam.AlignmentFile(in_bam, "rb", threads=threads) as read_bam:
    with pysam.AlignmentFile(
        out_bam, "wb", template=read_bam, threads=threads
    ) as out_bam:
        for read in read_bam.fetch(until_eof=True):
            for _ in range(read_count(read.query_name)):
                out_bam.write(read)
