import os
import logging

import pandas as pd
from snakemake.shell import shell


def tesmall_normalisation(method, count_file):
    """
    Calculate normalisation factor based on TEsmall count summary.

    method: "totalReads" or "miRNAReads"
    count_file: path to TEsmall count summary file (results/te_small/{sample}/count_summary.txt)
    """
    count_column_name = os.path.basename(os.path.dirname(count_file))

    df = pd.read_csv(count_file, sep="\t")
    if method == "totalReads":
        total_reads = df[count_column_name].sum()

        # normalise to per million reads
        norm_factor = 1 / total_reads * 1e6
    elif method == "miRNAReads":
        miRNA_reads = df[df["ftype"] == "miRNA"][count_column_name].sum()

        # normalise to per million miRNA reads
        norm_factor = 1 / miRNA_reads * 1e6
    else:
        raise ValueError(f"Unknown normalisation method: {method}")
    return norm_factor


# Get inputs/outputs/params from Snakemake
# --------------------------------------------------
tesmall_counts = snakemake.input["txt"]
bam = snakemake.input["bam"]
bigwig = snakemake.output["bigwig"]

threads = snakemake.threads
bin_size = snakemake.params["bin_size"]
normalise = snakemake.params["normalise"]
strand = snakemake.params["strand"]

# Set up logging
# --------------------------------------------------
log = snakemake.log[0]
logging.basicConfig(
    filename=log,
    filemode="w",
    format="%(levelname)s:%(asctime)s:%(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
)

# Generate bamCoverage command
# --------------------------------------------------

if strand == "forward":
    strand_flag = "--samFlagExclude 16"
    direction = 1
else:
    strand_flag = "--samFlagInclude 16"
    direction = -1

if normalise in ["RPKM", "CPM", "BPM", "RPGC", "None"]:
    # Use standard normalisation methods
    normalise_flag = f"--normalizeUsing {normalise} --scaleFactor {direction}"
elif normalise == "totalReads":
    # Normalise by total reads
    # Calculate scale factor based on total reads from the TEsmall count summary
    scale_factor = tesmall_normalisation("totalReads", tesmall_counts) * direction
    normalise_flag = f"--normalizeUsing None --scaleFactor {scale_factor}"
elif normalise == "miRNAReads":
    # Normalise by miRNA reads
    # Calculate scale factor based on total miRNA reads from the TEsmall count summary
    scale_factor = tesmall_normalisation("miRNAReads", tesmall_counts) * direction
    normalise_flag = f"--normalizeUsing None --scaleFactor {scale_factor}"
else:
    raise ValueError(f"Unknown normalisation method: {normalise}")

cmd = (
    f"bamCoverage "
    f"-b {bam} "
    f"-o {bigwig} "
    f"--binSize {bin_size} "
    f"{strand_flag} "
    f"{normalise_flag} "
    f"--numberOfProcessors {threads} "
    f"2> {log}"
)

logging.info(f"Running bamCoverage with command: {cmd}")

# Run bamCoverage command
# --------------------------------------------------
shell(cmd)

# Write the scale factor to a file for record-keeping
# --------------------------------------------------
scale_factor_file = snakemake.output["scale_factor"]
with open(scale_factor_file, "w") as f:
    f.write(f"{scale_factor}\n")
