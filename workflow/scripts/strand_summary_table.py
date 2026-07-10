import logging
import os
import re

import pandas as pd

# Load Snakemake variables
# --------------------------------------------------
count_summary_files = snakemake.input["txt"]
pingpong_logs = snakemake.input["pingpong_logs"]
scale_factor_files = snakemake.input["scale_factors"]
summary = snakemake.output["summary"]

# Set up logging
# --------------------------------------------------
logging.basicConfig(
    filename=snakemake.log[0],
    filemode="w",
    format="%(levelname)s:%(asctime)s:%(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
)

# Regex for lines such as:
# INFO:2026-07-10 13:15:46:IAP: FWD reads=10992, REV reads=15784, pairs analysed (distance <= 30)=2884451
TE_LINE_RE = re.compile(
    r"^INFO:\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}:"
    r"(?P<te>.+): FWD reads=(?P<fwd>\d+), REV reads=(?P<rev>\d+), "
    r"pairs analysed \(distance <= \d+\)=(?P<pairs>\d+)$"
)


def sample_from_count_summary(path):
    # results/te_small/{sample}/count_summary.txt
    return os.path.basename(os.path.dirname(path))


def sample_from_pingpong_log(path):
    # logs/pingpong/{sample}_analysis.log
    return os.path.basename(path).replace("_analysis.log", "")


def sample_from_scale_factor(path):
    # results/te_bigwig/{sample}_scale_factor_fwd.txt
    return os.path.basename(path).replace("_scale_factor_fwd.txt", "")


# Total mapped reads per sample (sum across all TEsmall feature types,
# i.e. the sample's full annotated small-RNA library, not just TE-mapped
# reads) from the TEsmall count summary.
# --------------------------------------------------
total_mapped_reads = {}
for count_file in count_summary_files:
    sample = sample_from_count_summary(count_file)
    logging.info(f"Reading total mapped reads for sample {sample} from {count_file}")
    df = pd.read_csv(count_file, sep="\t")
    total_mapped_reads[sample] = df[sample].sum()

# Per-sample scale factor (used to normalise raw FWD/REV read counts).
# Only the forward-strand scale factor is available as input, but its
# magnitude is identical for the reverse strand (stranded_te_bigwig.py
# applies the same TEsmall-derived factor to both strands, only negating
# it for the reverse strand so coverage renders below zero), so its
# absolute value is reused for both FWD and REV normalisation here.
# --------------------------------------------------
scale_factors = {}
for scale_factor_file in scale_factor_files:
    sample = sample_from_scale_factor(scale_factor_file)
    with open(scale_factor_file) as f:
        scale_factors[sample] = abs(float(f.read().strip()))
    logging.info(f"Scale factor for sample {sample}: {scale_factors[sample]}")

# Parse per-TE FWD/REV/pairs counts from each sample's ping-pong analysis log
# --------------------------------------------------
rows = []
for log_file in pingpong_logs:
    sample = sample_from_pingpong_log(log_file)
    scale_factor = scale_factors[sample]

    with open(log_file) as f:
        for line in f:
            match = TE_LINE_RE.match(line.strip())
            if match is None:
                continue

            te = match.group("te")
            fwd_reads = int(match.group("fwd"))
            rev_reads = int(match.group("rev"))
            pairs = int(match.group("pairs"))

            rows.append(
                {
                    "sample": sample,
                    "te": te,
                    "total_mapped_reads": total_mapped_reads[sample],
                    "fwd_reads": fwd_reads,
                    "rev_reads": rev_reads,
                    "fwd_reads_normalised": round(fwd_reads * scale_factor),
                    "rev_reads_normalised": round(rev_reads * scale_factor),
                    "pairs_analysed": pairs,
                }
            )

    logging.info(f"Parsed {len(rows)} TE rows so far after sample {sample}")

# Save table to CSV
# --------------------------------------------------
summary_df = pd.DataFrame(rows).sort_values(["te", "sample"]).reset_index(drop=True)
summary_df.to_csv(summary, index=False)

logging.info(f"Wrote strand summary table with {len(summary_df)} rows to {summary}")
