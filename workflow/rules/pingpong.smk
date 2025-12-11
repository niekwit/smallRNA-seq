# Build bowtie1 index for ping-pong analysis
# --------------------------------------------
rule bowtie_index:
    input:
        fasta=config["pingpong"]["fasta"]
    output:
        index_files=expand("resources/bowtie1_index/pingpong.{ext}", ext=EXT)
    params:
        index_prefix=lambda wildcards, output: output.index_files[0].replace(".1.ebwt", "")
    threads: 1
    log:
        "logs/pingpong/build_bowtie1_index.log"
    conda:
        "../envs/pingpong.yaml"
    shell:
        "bowtie-build {input.fasta} resources/bowtie1_index/pingpong {log}"


# Collapse unique sequences and count occurrences
# ------------------------------------------------------------
rule collapse_sequences:
    input:
        fastq="results/te_small/{sample}/{sample}.rm_rRNA.fastq",
    output:
        fasta="results/te_small/{sample}/{sample}.rm_rRNA.collapsed.fasta",
    threads: 2
    log:
        "logs/pingpong/{sample}_collapse.log"
    conda:
        "../envs/pingpong.yaml"
    shell:
        "fastx_collapser -i {input.fastq} -o {output.fasta} {log}"


# Align collapsed sequences to bowtie1 index
# Bowtie parameters based on:
# https://github.com/rberrens/SPOCD1-piRNA_directed_DNA_met/blob/master/piRNA_analysis/snakemake/13_map_consensus_L1s_IAPs.py
# ------------------------------------------------------------
rule align:
    input:
        fasta="results/te_small/{sample}/{sample}.rm_rRNA.collapsed.fasta",
        index=expand("resources/bowtie1_index/pingpong.{ext}", ext=EXT)
    output:
        bam="results/pingpong/{sample}.bam"
    params:
        index_prefix="resources/bowtie1_index/pingpong",
        mismatch=config["pingpong"]["mismatch"]
    threads: 4
    log:
        "logs/pingpong/{sample}_align.log"
    conda:
        "../envs/pingpong.yaml"
    shell:
        "bowtie "
        "-f --sam "
        "--threads {threads} "
        "-v {params.mismatch} "
        "-a --best -k 1 "
        "-x {params.index_prefix} "
        "{input.fasta} 2> {log} | "
        "samtools view "
        "-F 4 -bS - | "
        "samtools sort -o {output.bam}"

# Perform ping-pong analysis
# ------------------------------------------------------------
rule pingpong_analysis:
    input:
        bam="results/pingpong/{sample}.bam",
    output:
        pingpong="results/pingpong/{sample}.csv",
    params:
        window=config["pingpong"]["window"]
    threads: 2
    log:
        "logs/pingpong/{sample}_analysis.log"
    conda:
        "../envs/pingpong.yaml"
    script:
        "../scripts/pingpong.py"

# Plot ping-pong results
# ------------------------------------------------------------
rule plot_pingpong:
    input:
        pingpong=expand("results/pingpong/{sample}.csv", sample=SAMPLES)
    output:
        pdf="results/plots/pingpong.pdf",
        csv="results/plots/pingpong.csv"
    params:
    threads: 1
    log:
        "logs/pingpong/plot_pingpong.log"
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_pingpong.R"