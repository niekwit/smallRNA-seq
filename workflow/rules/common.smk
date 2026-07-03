import glob
import re
import sys
import pandas as pd
from snakemake.utils import validate
from snakemake.logging import logger


# validate config file
validate(config, schema="../../config/schemas/config.schema.yaml")


def samples():

    fastq = glob.glob("reads/*/*.fastq.gz")

    assert len(fastq) > 0, "No FASTQ files found in reads/ directory!"

    return [os.path.basename(f).split(".")[0] for f in fastq]


def get_dependencies(wildcards):
    """Input function for TEsmall rule"""

    inputs = {
        "fastq": f"reads/{wildcards.sample}/{wildcards.sample}.fastq.gz",
    }

    # When the genome database does not yet exist, LEADER_SAMPLE is set to the
    # first sample so it runs alone and builds the database; all other samples
    # wait for it. When the database already exists LEADER_SAMPLE is None and
    # every sample runs in parallel immediately.
    if LEADER_SAMPLE is not None and wildcards.sample != LEADER_SAMPLE:
        inputs["leader_finished"] = (
            f"results/te_small/{LEADER_SAMPLE}/count_summary.txt"
        )

    return inputs


def comparisons(config):
    reference_condition = config["reference_condition"]

    samples_dict = config["samples"]
    other_conditions = []
    for sample, condition in samples_dict.items():
        if condition != reference_condition:
            other_conditions.append(condition)

    other_conditions = list(set(other_conditions))

    return [f"{cond}_vs_{reference_condition}" for cond in other_conditions]


def conditions(config):
    samples_dict = config["samples"]
    return list(set(samples_dict.values()))


def te(config):

    fasta_file = config["pingpong"]["fasta"]

    # Cleaned up fasta headers are wildcard values
    te_list = []
    with open(fasta_file, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                te_name = line[1:].strip().split()[0]
                te_list.append(te_name)

    return te_list


def mirbase_genome(config):

    genome = config["genome"]

    if re.match("hg", genome):
        return "Homo_sapiens"
    elif re.match("mm", genome):
        return "Mus_musculus"
    else:
        raise ValueError(f"Unsupported genome: {genome}")


def ensembl_ncrna_url(config):
    genome = config["genome"]
    release = config["release"]

    genome_map = {
        "mm39": ("mus_musculus", "GRCm39", "Mus_musculus"),
        "hg38": ("homo_sapiens", "GRCh38", "Homo_sapiens"),
    }

    if genome not in genome_map:
        raise ValueError(f"Unsupported genome for Ensembl ncRNA download: {genome}")

    species, assembly, species_cap = genome_map[genome]
    filename = f"{species_cap}.{assembly}.ncrna.fa.gz"
    return f"https://ftp.ensembl.org/pub/release-{release}/fasta/{species}/ncrna/{filename}"


def check_mirna_correction_settings(config):
    if not config["length_distribution"]["apply_mirna_correction"]:
        return

    # Extract -m INT from cutadapt extra arguments
    cutadapt_extra = config["cutadapt"].get("extra", "")
    m = re.search(r"-m\s+(\d+)", cutadapt_extra)
    if m and int(m.group(1)) > 24:
        logger.warning(
            f"cutadapt -m is set to {m.group(1)}, which will discard reads "
            "shorter than this length. miRNAs are typically ≤24 nt, so miRNA-based "
            "correction may be inaccurate. Consider setting -m to 16 or lower."
        )

    # Check tesmall minlen when TEsmall is enabled
    if config["length_distribution"]["mirna_correction_method"] == "align":
        minlen = config["tesmall"]["minlen"]
        if minlen > 24:
            logger.warning(
                f"tesmall minlen is set to {minlen}, which excludes reads "
                "shorter than this length. miRNAs are typically ≤24 nt, so miRNA-based "
                "correction may be inaccurate. Consider setting minlen to 16 or lower."
            )


def plot_length_distribution_input(wildcards):
    input = {
        "counts": expand(
            "results/length_distribution/{sample}_length_distribution.txt",
            sample=SAMPLES,
        )
    }

    apply_correction = config["length_distribution"]["apply_mirna_correction"]
    method = config["length_distribution"]["mirna_correction_method"]
    if apply_correction and method == "align":
        input["mirna_fasta"] = expand("results/mirna/{sample}.fasta", sample=SAMPLES)

    return input


def length_normalisation_method(config):
    if not config["length_distribution"]["apply_mirna_correction"]:
        return "total_count_normalisation"
    else:
        method = config["length_distribution"]["mirna_correction_method"]
        if method == "align":
            return "mirna_aligned_count_normalisation"
        elif method == "length":
            return "mirna_length_count_normalisation"
        else:
            raise ValueError(f"Unsupported miRNA correction method: {method}")
