import glob
import re
import pandas as pd

# from snakemake.utils import validate


# validate sample sheet and config file
# validate(samples, schema="../../config/schemas/samples.schema.yaml")
# validate(config, schema="../../config/schemas/config.schema.yaml")


def samples():

    fastq = glob.glob("reads/*/*.fastq.gz")

    assert len(fastq) > 0, "No FASTQ files found in reads/ directory!"

    return [os.path.basename(f).split(".")[0] for f in fastq]


def get_dependencies(wildcards):
    inputs = {
        "fastq": f"reads/{wildcards.sample}/{wildcards.sample}.fastq.gz",
        "log": "logs/install_te_small.log",
    }

    # Check if genome data base exists
    # if not, a leader job will download it
    # otherwise run all jobs in parallel
    # fasta = os.path.join(config["tesmall"]["dbfolder"], "genomes", GENOME, "sequence", "genome.fa")
    # bed = os.path.join(config["tesmall"]["dbfolder"], "genomes", GENOME, "annotations", "exon.bed")

    # if not os.path.exists(fasta) or not os.path.exists(bed):
    #    inputs["leader_finished"] = f"results/te_small/{LEADER_SAMPLE}_count_summary.txt"

    # If this job is NOT the leader, make it wait for the leader
    if wildcards.sample != LEADER_SAMPLE:
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
        input["sam"] = expand(
            "results/mirna_correction/{sample}_count.sam", sample=SAMPLES
        )

    return input
