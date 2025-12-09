# Collapse unique sequences and count occurrences
# ------------------------------------------------------------
rule collapse_sequences:
    input:
        fastq="results/te_small/{sample}/{sample}.rm_rRNA.fastq",
    output:
        fasta="results/te_small/{sample}/{sample}.rm_rRNA.collapsed.fasta",
        counts="results/pingpong/{sample}_counts.txt",
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
        mismatch=3
    threads: 4
    log:
        "logs/pingpong/{sample}_align.log"
    conda:
        "../envs/pingpong.yaml"
    shell:
        "bowtie"
        "-f --sam "
        "--threads {threads} "
        "-v {params.mismatch} "
        "-a --best -k 1 "
        "-x {params.index_prefix} "
        "{input.fasta} 2> {log} | "
        "samtools view "
        "-F 4 -bS - | "
        "samtools sort -o {output.bam} -"


# Index BAM files
# ------------------------------------------------------------
rule index_bam:
    input:
        bam="results/pingpong/{sample}.bam"
    output:
        bai="results/pingpong/{sample}.bam.bai"
    threads: 1
    log:
        "logs/pingpong/{sample}_index_bam.log"
    conda:
        "../envs/pingpong.yaml"
    shell:
        "samtools index {input.bam} {output.bai} {log}"


# Perform ping-pong analysis
# ------------------------------------------------------------