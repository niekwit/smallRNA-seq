import logging

import pyBigWig
from Bio import SeqIO
from snakemake.shell import shell

# Load variables from Snakemake
fwd_bw = snakemake.input["fwd_bw"]
rev_bw = snakemake.input["rev_bw"]
all_bw = fwd_bw + rev_bw
te_fasta = snakemake.input["te_fasta"]
ini_files = snakemake.output["ini"]
pdf_files = snakemake.output["pdf"]
width = snakemake.params["width"]
conditions = snakemake.params["conditions"]
log = snakemake.log["script"]
log_stdout = snakemake.log["stdout"]

# Set up logging
logging.basicConfig(
    format="%(levelname)s:%(asctime)s:%(message)s",
    level=logging.DEBUG,
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.FileHandler(log)],
)

# Get the reference conditions
# this will be plotted first in the pyGenomeTracks output
ref_condition = snakemake.config["reference_condition"]

# Order the conditions so that the reference condition is first
ordered_conditions = [ref_condition] + [c for c in conditions if c != ref_condition]

# Load TE sequences from FASTA file
# and store the length of each TE sequence in a dictionary
logging.info(f"Loading TE sequences from {te_fasta}")
te_seq = {}
for te in SeqIO.parse(te_fasta, "fasta"):
    te_seq[te.id] = len(te.seq)

# Open each bigWig and find maximum value across all conditions/strands
# for each TE sequence --> max_value for pyGenomeTracks
logging.info("Finding maximum values for each TE sequence across all bigWig files")
max_values = {}
for bw in all_bw:
    with pyBigWig.open(bw) as bw_file:
        for te in te_seq.keys():
            # Get the length of the TE sequence
            te_length = te_seq[te]

            # Get the maximum value for this TE sequence in the bigWig file
            # https://www.biostars.org/p/244898/
            max_value = bw_file.stats(te, type="max", start=0, end=te_length)[0]
            if max_value is not None:
                max_value = int(max_value * 1.05)
                if te not in max_values:
                    max_values[te] = max_value
                else:
                    max_values[te] = max(max_values[te], max_value)

# Iterate over each TE sequence:
# 1. Create ini file for pyGenomeTracks
# 2. Run pyGenomeTracks to generate the track image
for te in te_seq.keys():
    # Create ini file for pyGenomeTracks
    logging.info(f"Creating ini file for TE: {te}")
    ini = [x for x in ini_files if te in x][0]
    with open(ini, "w") as f:
        f.write("[x-axis]\n\n")
        for c in conditions:
            for strand in ["fwd", "rev"]:
                f.write(f"[{c} {strand}]\n")
                bw_file = [x for x in all_bw if c in x and strand in x][0]
                f.write(f"file = {bw_file}\n")
                if strand == "fwd":
                    f.write(f"title = {c}\n")
                f.write("color = #2166AC\n" if strand == "fwd" else "color = #dd3b3b\n")
                if strand == "rev":
                    f.write("overlay_previous = yes\n")
                f.write("height = 3\n")
                f.write("type = fill\n")
                max_value = max_values[te]
                f.write(f"min_value = -{max_value}\n")
                f.write(f"max_value = {max_value}\n\n")
                if strand == "rev":
                    f.write("[spacer]\n")
                    f.write("height = 0.3\n\n")

    # Now run pyGenomeTracks to generate the track image
    logging.info(f"Generating track image for TE: {te}")
    pdf = [x for x in pdf_files if te in x][0]
    cmd = f"pyGenomeTracks --tracks {ini} --region {te}:1-{te_seq[te]} -o {pdf} --width {width}"
    logging.debug(f"Running command: {cmd}")
    shell(cmd)

logging.info("All TE tracks generated")
