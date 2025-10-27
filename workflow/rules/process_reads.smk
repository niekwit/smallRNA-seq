# Adapter sequence removal
# -----------------------------------------------------
rule trim_adapters:
    input:
        read="reads/{sample}.fastq.gz",
    output:
        trimmed="results/trimmed/{sample}_trimmed.fastq.gz",
    params:
        adapter=config["cutadapt"]["adapter"],
        extra=config["cutadapt"]["extra"]
    threads: 4
    conda:
        "../envs/process_reads.yaml"
    shell:
        "cutadapt "
        "-j {threads} "
        "-m 16 "
        "--trimmed-only "
        "-a {params.adapter} "
        "-o {output.trimmed} "
        "{params.extra} "     
        "{input.read}"



# Get just sequence from reads (every 4th line, starting at line 2)
# -----------------------------------------------------
rule get_sequence:
    input:
        trimmed="results/trimmed/{sample}_trimmed.fastq.gz",
    output:
        seqs="results/seqs/{sample}_seqs.txt",
    conda:
        "../envs/process_reads.yaml"
    shell:
        "zcat {input.trimmed} | sed -n '2~4p' > {output.seqs}"


# Count unique sequences
# -----------------------------------------------------
rule count_unique_sequences:
    input:
        seqs="results/seqs/{sample}_seqs.txt",
    output:
        counts="results/counts/{sample}_unique_counts.txt",
    conda:
        "../envs/process_reads.yaml"
    shell:
        r"""sort {input.seqs} | uniq -c | awk '{{print $2"\t"$1}}' > {output.counts}"""


# Create count fasta file
# Format: >sequence:length:count
# -----------------------------------------------------
rule create_count_fasta:
    input:
        counts="results/counts/{sample}_unique_counts.txt",
    output:
        fasta="results/fasta/{sample}_counts.fasta",
    conda:
        "../envs/process_reads.yaml"
    shell:
        r"""awk '{{print ">" $1 ":" length($1) ":" $2 "\n" $1}}' {input.counts} > {output.fasta}"""

