# Adapter sequence removal
# -----------------------------------------------------
rule trim_adapters:
    input:
        "reads/{sample}/{sample}.fastq.gz",
    output:
        fastq=temp("results/trimmed/{sample}.fastq.gz"),
        qc="results/trimmed/{sample}_qc.txt",
    params:
        adapters=f"-a {config['cutadapt']['adapter']}",
        extra=config["cutadapt"]["extra"],
    log:
        "logs/cutadapt/{sample}.log"
    threads: 4
    wrapper:
        "v7.9.0/bio/cutadapt/se"

# Get just sequence from reads (every 4th line, starting at line 2)
# -----------------------------------------------------
rule get_sequence:
    input:
        trimmed="results/trimmed/{sample}.fastq.gz",
    output:
        seqs="results/seqs/{sample}_seqs.txt",
    log:
        "logs/get_sequence/{sample}.log"
    conda:
        "../envs/process_reads.yaml"
    shell:
        "zcat {input.trimmed} | sed -n '2~4p' > {output.seqs}"


# Generate length counts
# Format: length \t count \t sample
# -----------------------------------------------------
rule length_counts:
    input:
        seqs="results/seqs/{sample}_seqs.txt",
    output:
        txt="results/length_distribution/{sample}_length_distribution.txt",
    log:
        "logs/length_counts/{sample}.log"
    conda:
        "../envs/process_reads.yaml"
    shell:
        "awk '{{print length($0)}}' {input.seqs} | "
        "sort | "
        "uniq -c | "
        r"sed 's/^\s*//;s/\s/\t/g;s/$/\t{wildcards.sample}/' > {output.txt} "


# Plot length distribution of trimmed sequences
# -----------------------------------------------------
rule plot_length_distribution:
    input:
        counts=expand("results/length_distribution/{sample}_length_distribution.txt", sample=SAMPLES),
    output:
        pdf="results/plots/length_distribution.pdf",
        csv="results/plots/length_distribution.csv",
    log:
        "logs/plot_length_distribution.log"
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_length_distribution.R"