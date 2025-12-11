import glob
import pandas as pd
#from snakemake.utils import validate


# validate sample sheet and config file
#validate(samples, schema="../../config/schemas/samples.schema.yaml")
#validate(config, schema="../../config/schemas/config.schema.yaml")


def samples():

    fastq = glob.glob("reads/*/*.fastq.gz")
    return [os.path.basename(f).split(".")[0] for f in fastq]


def get_dependencies(wildcards):
    inputs = {
        "fastq": f"reads/{wildcards.sample}/{wildcards.sample}.fastq.gz",
        "log": "logs/install_te_small.log"
    }
    
    # Check if genome data base exists
    # if not, a leader job will download it
    # otherwise run all jobs in parallel
    #fasta = os.path.join(config["tesmall"]["dbfolder"], "genomes", GENOME, "sequence", "genome.fa")
    #bed = os.path.join(config["tesmall"]["dbfolder"], "genomes", GENOME, "annotations", "exon.bed")

    #if not os.path.exists(fasta) or not os.path.exists(bed):
    #    inputs["leader_finished"] = f"results/te_small/{LEADER_SAMPLE}_count_summary.txt"

    # If this job is NOT the leader, make it wait for the leader
    if wildcards.sample != LEADER_SAMPLE:
        inputs["leader_finished"] = f"results/te_small/{LEADER_SAMPLE}/count_summary.txt"
        
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