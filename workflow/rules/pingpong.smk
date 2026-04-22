# Download Ensembl ncRNA FASTA for the configured genome and release
# ------------------------------------------------------------
rule download_ncrna_fasta:
    output:
        fa=f"resources/ncrna/{GENOME}.ncrna.fa.gz",
    params:
        url=NCRNA_URL,
    retries: 3
    log:
        f"logs/ncrna/download.log",
    conda:
        "../envs/pingpong.yaml"
    shell:
        "wget -q -O {output.fa} {params.url} 2> {log}"


# Extract excluded ncRNA biotypes to use as the filter reference.
# Reads mapping to these sequences will be removed before ping-pong
# analysis; lncRNAs are intentionally excluded from this filter so
# that reads from lncRNA-hosted piRNA clusters are retained.
# ------------------------------------------------------------
rule filter_ncrna_fasta:
    input:
        fa=f"resources/ncrna/{GENOME}.ncrna.fa.gz",
    output:
        fa=f"resources/ncrna/{GENOME}.excluded_ncrna.fa",
    log:
        f"logs/ncrna/filter_fasta.log",
    conda:
        "../envs/pingpong.yaml"
    shell:
        "seqkit grep "
        "-r "
        '-p "gene_biotype:(miRNA|misc_RNA|Mt_rRNA|Mt_tRNA|ribozyme|rRNA|scaRNA|scRNA|snoRNA|snRNA|sRNA)" '
        "{input.fa} "
        "-o {output.fa} "
        "2> {log}"


# Build bowtie1 index from the excluded ncRNA sequences
# ------------------------------------------------------------
rule bowtie_index_ncrna:
    input:
        fa=f"resources/ncrna/{GENOME}.excluded_ncrna.fa",
    output:
        index_files=expand(
            f"resources/bowtie1_index/ncrna_{GENOME}.{{ext}}", ext=EXT
        ),
    params:
        index_prefix=f"resources/bowtie1_index/ncrna_{GENOME}",
    threads: 1
    log:
        f"logs/ncrna/bowtie_index.log",
    conda:
        "../envs/pingpong.yaml"
    shell:
        "bowtie-build {input.fa} {params.index_prefix} > {log} 2>&1"


# Filter collapsed reads against excluded ncRNA index.
# Non-aligning reads (piRNA candidates) continue to ping-pong analysis.
# Aligning reads are saved as a BAM for QC/troubleshooting.
# ------------------------------------------------------------
rule filter_ncrna_reads:
    input:
        fasta="results/fasta/{sample}.collapsed.fasta",
        index=expand(
            f"resources/bowtie1_index/ncrna_{GENOME}.{{ext}}", ext=EXT
        ),
    output:
        fasta="results/fasta/{sample}.ncrna_filtered.fasta",
        bam="results/ncrna_filter/{sample}.bam",
    params:
        index_prefix=f"resources/bowtie1_index/ncrna_{GENOME}",
    threads: 4
    log:
        "logs/ncrna/{sample}_filter_reads.log",
    conda:
        "../envs/pingpong.yaml"
    shell:
        "bowtie "
        "-f --sam "
        "--threads {threads} "
        "-v 0 "
        "-a --best -k 1 "
        "--un {output.fasta} "
        "-x {params.index_prefix} "
        "{input.fasta} 2> {log} | "
        "samtools view -F 4 -bS - | "
        "samtools sort -o {output.bam}"


# Build bowtie1 index for ping-pong analysis
# --------------------------------------------
rule bowtie_index:
    input:
        fasta=config["pingpong"]["fasta"],
    output:
        index_files=expand("resources/bowtie1_index/pingpong.{ext}", ext=EXT),
    params:
        index_prefix=lambda wildcards, output: output.index_files[0].replace(
            ".1.ebwt", ""
        ),
    threads: 1
    log:
        "logs/bowtie/pingpong.log",
    conda:
        "../envs/pingpong.yaml"
    shell:
        "bowtie-build {input.fasta} resources/bowtie1_index/pingpong {log}"


# Collapse unique sequences and count occurrences
# ------------------------------------------------------------
rule collapse_sequences:
    input:
        fastq="results/trimmed/{sample}.fastq.gz",
    output:
        fasta="results/fasta/{sample}.collapsed.fasta",
    threads: 2
    log:
        "logs/pingpong/{sample}_collapse.log",
    conda:
        "../envs/pingpong.yaml"
    shell:
        "zcat {input.fastq} | fastx_collapser -o {output.fasta} {log}"


# Align collapsed sequences to bowtie1 index
# Bowtie parameters based on:
# https://github.com/rberrens/SPOCD1-piRNA_directed_DNA_met/blob/master/piRNA_analysis/snakemake/13_map_consensus_L1s_IAPs.py
# ------------------------------------------------------------
rule align:
    input:
        fasta="results/fasta/{sample}.ncrna_filtered.fasta",
        index=expand("resources/bowtie1_index/pingpong.{ext}", ext=EXT),
    output:
        bam="results/pingpong/{sample}.bam",
    params:
        index_prefix=lambda wildcards, input: input.index[0].replace(".1.ebwt", ""),
        mismatch=config["pingpong"]["mismatch"],
    threads: 4
    log:
        "logs/pingpong/{sample}_align.log",
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
        nt_bias="results/pingpong/{sample}_nt_bias.csv",
    params:
        window=config["pingpong"]["window"],
    threads: 2
    log:
        "logs/pingpong/{sample}_analysis.log",
    conda:
        "../envs/pingpong.yaml"
    script:
        "../scripts/pingpong.py"


# Plot ping-pong results
# ------------------------------------------------------------
rule plot_pingpong:
    input:
        pingpong=expand("results/pingpong/{sample}.csv", sample=SAMPLES),
    output:
        pdf="results/plots/pingpong.pdf",
        csv="results/plots/pingpong.csv",
    threads: 1
    log:
        "logs/pingpong/plot_pingpong.log",
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/plot_pingpong.R"


# Plot 1U/10A bias
# ------------------------------------------------------------
rule plot_sequence_bias_pirna:
    input:
        bam=expand("results/pingpong/{sample}.bam", sample=SAMPLES),
        nt_bias=expand("results/pingpong/{sample}_nt_bias.csv", sample=SAMPLES),
    output:
        logo="results/plots/pirna_sequence_bias/{te}/{comparison}_logo.pdf",
        bar="results/plots/pirna_sequence_bias/{te}/{comparison}_bargraph.pdf",
        csv="results/plots/pirna_sequence_bias/{te}/{comparison}_frequencies.csv",
    threads: 2
    log:
        "logs/pingpong/{te}_plot_sequence_bias_{comparison}.log",
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/sequence_bias_pirna.R"


# Plot per-TE positional read coverage (forward/reverse strand)
# ------------------------------------------------------------
rule te_pos_coverage:
    input:
        bam=expand("results/pingpong/{sample}.bam", sample=SAMPLES),
        fasta=config["pingpong"]["fasta"],
    output:
        pdf="results/plots/te_pos_coverage/{te}.pdf",
    threads: 2
    log:
        "logs/pingpong/{te}_te_pos_coverage.log",
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/te_pos_coverage.R"


"""
# Get length distribution of TE-aligned reads
# ------------------------------------------------------------
rule length_distribution_aligned_to_TE:
    input:
        bam=expand("results/pingpong/{sample}.bam", sample=SAMPLES),
    output:
        pdf="results/plots/length_distribution.pdf",
        csv="results/plots/length_distribution.csv"
    threads: 1
    log:
        "logs/pingpong/length_distribution.log"
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/length_distribution.R"
"""
